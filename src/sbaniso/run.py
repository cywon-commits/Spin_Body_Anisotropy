"""One YAML = one run = one HDF5 file.

    python -m sbaniso.run configs/validation/xy_L32.yaml [--threads 8]

YAML keys
    name            free text (default: file stem)
    model           {r: 4.1, G: 1.0}  (rigid couplings)  or  {couplings: [J1, J2, J4]}
    L               linear size (square lattice, periodic)
    temperatures    {Tmin, Tmax, n}  (geometric ladder)  or  {list: [T0, T1, ...]}
    ntherm, nmeas   sweeps of equilibration and of measurement
    clusters        [sigma_phi, sigma_psi, phi] enable flags (default [1, 1, 1]).  The C3 model
                    near the crossing REQUIRES sigma_psi = 1 (CLAUDE.md, MC rule 2).
    seed            integer; all random streams derive from it
    start           random | ordered   (default random)
    checkpoint_every  sweeps between checkpoints (default 1000)
    output          HDF5 path (default runs/<name>.h5, relative to the repository root)

HDF5 layout
    attrs           config (YAML text), git hash, package versions, created/updated, status
    params/         Ts, betas, couplings, L, flags, ntherm, nmeas, seed, obs_names
    timeseries/obs  (R, nmeas, NOBS) per temperature slot, measured sweeps only
    diagnostics/    swap_acc, swap_tries, roundtrips per replica, ncl, deltas, cluster mean sizes,
                    Metropolis acceptance, wall-clock seconds
    final/          S, P per replica, perm (replica at each temperature)
    checkpoint/     the full PTState (incl. RNG states) after the last completed chunk

A run is resumed automatically when the output exists with status 'running': the time series is
truncated to the checkpoint and the run continues bit-identically.  A file with status 'complete'
is left untouched (pass --force to overwrite).
"""
import argparse
import datetime as _dt
import os
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np
import yaml

from .c3model import NOBS, OBS_NAMES, couplings_from_r
from .pt import PTState, advance, geometric_ladder, work_arrays

REPO = Path(__file__).resolve().parents[2]


def git_hash():
    try:
        h = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
                                    stderr=subprocess.DEVNULL).strip()
        dirty = subprocess.call(["git", "diff", "--quiet", "HEAD", "--", "src"], cwd=REPO,
                                stderr=subprocess.DEVNULL)
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def src_digest():
    """SHA-256 over the package sources (identifies the code even when the tree is uncommitted)."""
    import hashlib
    h = hashlib.sha256()
    for f in sorted((REPO / "src" / "sbaniso").glob("*.py")):
        h.update(f.name.encode()); h.update(f.read_bytes())
    return h.hexdigest()[:16]


def load_config(path_or_dict):
    if isinstance(path_or_dict, dict):
        cfg = dict(path_or_dict); text = yaml.safe_dump(cfg, sort_keys=False)
        stem = cfg.get("name", "run")
    else:
        text = Path(path_or_dict).read_text(); cfg = yaml.safe_load(text)
        stem = Path(path_or_dict).stem
    cfg.setdefault("name", stem)
    cfg.setdefault("clusters", [1, 1, 1])
    cfg.setdefault("start", "random")
    cfg.setdefault("checkpoint_every", 1000)
    cfg.setdefault("output", f"runs/{cfg['name']}.h5")
    return cfg, text


def resolve(cfg):
    m = cfg["model"]
    c = np.asarray(m["couplings"], float) if "couplings" in m else couplings_from_r(m["r"], m.get("G", 1.0))
    t = cfg["temperatures"]
    Ts = np.asarray(t["list"], float) if "list" in t else geometric_ladder(t["Tmin"], t["Tmax"], t["n"])
    out = Path(cfg["output"])
    if not out.is_absolute():
        out = REPO / out
    return c, Ts, out


def _versions():
    import numba, scipy
    return f"numpy {np.__version__}, numba {numba.__version__}, scipy {scipy.__version__}, h5py {h5py.__version__}"


def _write_state(g, state):
    for k, v in state.to_dict().items():
        if k in g:
            del g[k]
        g[k] = v


def _write_diagnostics(f, state, secs):
    g = f.require_group("diagnostics")
    tries = np.maximum(state.swap_tries, 1)
    calls = np.maximum(state.acc_calls, 1)
    nmeas_done = max(1, state.it - int(f["params/ntherm"][()]))
    vals = dict(swap_acc=state.swap_acc / tries, swap_tries=state.swap_tries,
                roundtrips_per_replica=state.rt_rep, roundtrips=int(state.rt_rep.sum()),
                ncl=state.ncl, deltas=state.deltas, cluster_mean_size=state.acc_sizes / calls,
                metropolis_acc=state.metro_acc / nmeas_done, seconds=secs)
    for k, v in vals.items():
        if k in g:
            del g[k]
        g[k] = v


def run_config(config, threads=None, force=False, max_sweeps=None, verbose=True):
    """Run (or resume) the run described by `config` (path or dict).  Returns the HDF5 path.
    max_sweeps: stop after this many sweeps in this call (for tests of checkpoint/resume)."""
    cfg, text = load_config(config)
    c, Ts, out = resolve(cfg)
    if threads:
        import numba
        numba.set_num_threads(int(threads))
    L = int(cfg["L"]); R = len(Ts); ntherm = int(cfg["ntherm"]); nmeas = int(cfg["nmeas"])
    flags = np.asarray(cfg["clusters"], np.int64)
    total = ntherm + nmeas
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and not force:
        with h5py.File(out, "r") as f:
            status = f.attrs.get("status", "")
            if f.attrs.get("config", "") != text:
                raise RuntimeError(f"{out} exists with a different config; use --force")
        if status == "complete":
            if verbose:
                print(f"{out}: complete, skipped")
            return out
    else:
        if out.exists():
            out.unlink()
        state = PTState(L, Ts, cfg["seed"], flags, ordered=(cfg["start"] == "ordered"))
        with h5py.File(out, "w") as f:
            f.attrs.update(config=text, name=cfg["name"], git=git_hash(), src_sha256=src_digest(), versions=_versions(),
                           created=_dt.datetime.now().isoformat(timespec="seconds"), status="running")
            p = f.create_group("params")
            p["Ts"] = Ts; p["betas"] = 1 / Ts; p["couplings"] = c; p["L"] = L; p["flags"] = flags
            p["ntherm"] = ntherm; p["nmeas"] = nmeas; p["seed"] = int(cfg["seed"])
            p["obs_names"] = np.array(OBS_NAMES, dtype="S")
            f.create_group("timeseries").create_dataset(
                "obs", shape=(R, 0, NOBS), maxshape=(R, nmeas, NOBS), dtype="f8",
                chunks=(1, min(nmeas, 4096), NOBS))
            _write_state(f.create_group("checkpoint"), state)
            f.attrs["seconds"] = 0.0

    with h5py.File(out, "r+") as f:
        state = PTState.from_dict({k: f["checkpoint"][k][()] for k in PTState.FIELDS})
        ds = f["timeseries/obs"]          # drop rows written after the last checkpoint, if any
        if ds.shape[1] != max(0, state.it - ntherm):
            ds.resize(max(0, state.it - ntherm), axis=1)
        secs = float(f.attrs["seconds"])
        git_now = git_hash()
        if f.attrs["git"] != git_now and verbose:
            print(f"warning: resuming {out} written at {f.attrs['git']} with code at {git_now}")

    work = work_arrays(R, L)
    every = int(cfg["checkpoint_every"])
    budget = total - state.it if max_sweeps is None else min(total - state.it, int(max_sweeps))
    end = state.it + budget
    while state.it < end:
        n = min(every - state.it % every, end - state.it)
        t0 = time.time()
        rows = advance(state, L, 1 / Ts, c, flags, ntherm, n, work)
        secs += time.time() - t0
        with h5py.File(out, "r+") as f:
            ds = f["timeseries/obs"]
            n0 = ds.shape[1]
            if rows.shape[1]:
                ds.resize(n0 + rows.shape[1], axis=1)
                ds[:, n0:, :] = rows
            _write_state(f["checkpoint"], state)
            _write_diagnostics(f, state, secs)
            f.attrs["seconds"] = secs
            f.attrs["updated"] = _dt.datetime.now().isoformat(timespec="seconds")
            if state.it >= total:
                g = f.require_group("final")
                for k, v in dict(S=state.S, P=state.P, perm=state.perm).items():
                    if k in g:
                        del g[k]
                    g[k] = v
                f.attrs["status"] = "complete"
        if verbose:
            print(f"{cfg['name']}: sweep {state.it}/{total}  ({secs:.0f} s)", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("configs", nargs="+")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    for cf in a.configs:
        run_config(cf, threads=a.threads, force=a.force)


if __name__ == "__main__":
    main()
