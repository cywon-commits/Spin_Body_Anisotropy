"""Full spin-body model (task N3): polyspin PT with single-spin anisotropy, measured in the language
of the C3 effective model so that sbaniso.analysis applies unchanged.

    python -m sbaniso.fullmodel configs/n3/zH4_0.80_L16.yaml [--threads 4]

Per body a (spins S_ai, lab axis z), with c = z.S and w = (x + i y).S:
    q30_a = sum_i P3(c)            Ising variable sigma (octupole sign, TR odd)
    q33_a = sum_i w^3              e^{i phi}   (octupole phase)
    q43_a = sum_i c w^3            e^{i psi}   (hexadecapole phase, TR even)
    K_a   = sum_alpha n_alpha^4    axis selection (1/3: C3 along z, 1: S4 along z); body S4 axes
                                   e_alpha ~ S_0 + S_alpha
each normalised by its value for an ideal tetrahedron with a C3 axis along z (|q| = 1).
Stiffness for a uniform twist of the body angle alpha about z: rho_alpha = [<S2> - beta <S1^2>] / N per
bond direction (polyspin.stiffness, exact also with the anisotropy, which is invariant under rotations
about z).  psi = 3 alpha, so rho_psi = rho_alpha / 9; the KT criterion rho_psi = 2T/pi is
rho_alpha = 18 T/pi (docs/TASK.md N3).

Stored columns (same layout as sbaniso.c3model.OBS_NAMES):
    E (total energy, J units), m_sigma = sum q30, cx = S2_x/9, sx = S1_x/3, cy = S2_y/9, sy = S1_y/3,
    Re/Im sum q33, Re/Im sum q43, |F(k_min)|^2 of q30, q33, q43
so that analysis.derived gives U4 (sigma Binder), Ups = rho_psi, xi/L, etc.  Extra: K (mean K_a) and
the mean |q30| per body are stored in timeseries/extra.

YAML: name, model {J, K, zG, zH4, a2, a4}, L, temperatures {Tmin, Tmax, n} | {list}, ntherm, nmeas,
seed, start (random | ordered), checkpoint_every, output, overlap_move (default true), kappa (20).
"""
import argparse
import datetime as _dt
import time
from pathlib import Path

import h5py
import numpy as np
import yaml
from numba import njit, prange

from polyspin.stiffness import _bond_sums

from .c3model import NOBS, OBS_NAMES
from .run import REPO, git_hash, src_digest

NEXTRA = 2


def _norms():
    from polyspin.groundstate import oriented_tetrahedron
    S = oriented_tetrahedron([1, 1, 1], 0.0)
    c = S[:, 2]; w = S[:, 0] + 1j * S[:, 1]
    return (abs(np.sum(0.5 * (5 * c ** 3 - 3 * c))), abs(np.sum(w ** 3)), abs(np.sum(c * w ** 3)))


@njit(cache=True)
def _measure_one(spins, L, nrm, out, extra):
    N, n, _ = spins.shape
    m = 0.0; fr = 0.0; fi = 0.0; sr = 0.0; si = 0.0; kk = 0.0; aq = 0.0
    k = 2.0 * np.pi / L
    ax_s = 0.0 + 0.0j; ay_s = 0.0 + 0.0j; ax_1 = 0.0 + 0.0j; ay_1 = 0.0 + 0.0j
    ax_2 = 0.0 + 0.0j; ay_2 = 0.0 + 0.0j
    for a in range(N):
        q30 = 0.0; q33 = 0.0 + 0.0j; q43 = 0.0 + 0.0j
        for i in range(n):
            c = spins[a, i, 2]
            w = spins[a, i, 0] + 1j * spins[a, i, 1]
            w3 = w * w * w
            q30 += 0.5 * (5.0 * c * c * c - 3.0 * c)
            q33 += w3
            q43 += c * w3
        q30 /= nrm[0]; q33 /= nrm[1]; q43 /= nrm[2]
        Ka = 0.0
        for al in range(1, 4):
            ex = spins[a, 0, 0] + spins[a, al, 0]; ey = spins[a, 0, 1] + spins[a, al, 1]
            ez = spins[a, 0, 2] + spins[a, al, 2]
            r = np.sqrt(ex * ex + ey * ey + ez * ez)
            if r > 0:
                Ka += (ez / r) ** 4
        kk += Ka; aq += abs(q30)
        m += q30; fr += q33.real; fi += q33.imag; sr += q43.real; si += q43.imag
        x = a % L; y = a // L
        ex_ = np.cos(k * x) + 1j * np.sin(k * x); ey_ = np.cos(k * y) + 1j * np.sin(k * y)
        ax_s += q30 * ex_; ay_s += q30 * ey_
        ax_1 += q33 * ex_; ay_1 += q33 * ey_
        ax_2 += q43 * ex_; ay_2 += q43 * ey_
    out[1] = m; out[6] = fr; out[7] = fi; out[8] = sr; out[9] = si
    out[10] = 0.5 * (abs(ax_s) ** 2 + abs(ay_s) ** 2)
    out[11] = 0.5 * (abs(ax_1) ** 2 + abs(ay_1) ** 2)
    out[12] = 0.5 * (abs(ax_2) ** 2 + abs(ay_2) ** 2)
    extra[0] = kk / N; extra[1] = aq / N


@njit(parallel=True, cache=True)
def measure_all(spins_all, rep_at, energies, nbr, p, L, nrm, out, extra):
    """out[t] / extra[t] for the replica at temperature index t."""
    R = rep_at.shape[0]
    zax = np.array([0.0, 0.0, 1.0])
    for t in prange(R):
        r = rep_at[t]
        _measure_one(spins_all[r], L, nrm, out[t], extra[t])
        s1x, s2x = _bond_sums(spins_all[r], nbr, p[3], p[4], p[5], zax, 0)
        s1y, s2y = _bond_sums(spins_all[r], nbr, p[3], p[4], p[5], zax, 1)
        out[t, 0] = energies[r]
        out[t, 2] = s2x / 9.0; out[t, 3] = s1x / 3.0
        out[t, 4] = s2y / 9.0; out[t, 5] = s1y / 3.0


def load_config(path_or_dict):
    if isinstance(path_or_dict, dict):
        cfg = dict(path_or_dict); text = yaml.safe_dump(cfg, sort_keys=False)
    else:
        text = Path(path_or_dict).read_text(); cfg = yaml.safe_load(text)
        cfg.setdefault("name", Path(path_or_dict).stem)
    cfg.setdefault("start", "random"); cfg.setdefault("checkpoint_every", 500)
    cfg.setdefault("overlap_move", True); cfg.setdefault("kappa", 20.0)
    cfg.setdefault("output", f"runs/{cfg['name']}.h5")
    return cfg, text


def _setup(cfg):
    from polyspin.groundstate import phase_of_gs
    from polyspin.lattice import Lattice
    from polyspin.model import Params
    from polyspin.pt import ParallelTempering
    m = cfg["model"]
    lat = Lattice(2, int(cfg["L"]))
    P = Params.from_z(zG=m["zG"], zH4=m["zH4"], z=lat.z, J=m.get("J", 1.0), K=m["K"],
                      a2=m.get("a2", 0.0), a4=m.get("a4", 0.0))
    t = cfg["temperatures"]
    Ts = np.asarray(t["list"], float) if "list" in t else np.geomspace(t["Tmin"], t["Tmax"], t["n"])
    _, _, gs = phase_of_gs(P)
    pt = ParallelTempering(lat, P, Ts, seed=cfg["seed"], init=cfg["start"], ref_state=gs,
                           kappa=float(cfg["kappa"]), use_overlap=bool(cfg["overlap_move"]))
    return lat, P, Ts, pt


def run_config(config, threads=None, verbose=True, max_sweeps=None):
    cfg, text = load_config(config)
    if threads:
        import numba
        numba.set_num_threads(int(threads))
    out_path = Path(cfg["output"]); out_path = out_path if out_path.is_absolute() else REPO / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lat, P, Ts, pt = _setup(cfg)
    R = len(Ts); ntherm = int(cfg["ntherm"]); nmeas = int(cfg["nmeas"]); total = ntherm + nmeas
    nrm = np.array(_norms())
    secs = 0.0
    if out_path.exists():
        with h5py.File(out_path, "r") as f:
            if f.attrs["config"] != text:
                raise RuntimeError(f"{out_path} exists with a different config")
            if f.attrs["status"] == "complete":
                if verbose:
                    print(f"{out_path}: complete, skipped")
                return out_path
            pt.load_state_dict({k: f["checkpoint"][k][()] for k in f["checkpoint"]})
            secs = float(f.attrs["seconds"])
        with h5py.File(out_path, "r+") as f:
            ds = f["timeseries/obs"]; nrows = max(0, pt.sweeps_done - ntherm)
            ds.resize(nrows, axis=1); f["timeseries/extra"].resize(nrows, axis=1)
    else:
        with h5py.File(out_path, "w") as f:
            f.attrs.update(config=text, name=cfg["name"], git=git_hash(), src_sha256=src_digest(),
                           polyspin_git=_polyspin_git(), status="running", seconds=0.0,
                           created=_dt.datetime.now().isoformat(timespec="seconds"))
            g = f.create_group("params")
            g["Ts"] = Ts; g["betas"] = 1 / Ts; g["couplings"] = P.array(); g["L"] = lat.size
            g["ntherm"] = ntherm; g["nmeas"] = nmeas; g["seed"] = int(cfg["seed"])
            g["obs_names"] = np.array(OBS_NAMES, dtype="S"); g["extra_names"] = np.array(["K", "abs_q30"], dtype="S")
            ts = f.create_group("timeseries")
            ts.create_dataset("obs", shape=(R, 0, NOBS), maxshape=(R, nmeas, NOBS), dtype="f8",
                              chunks=(1, min(nmeas, 4096), NOBS))
            ts.create_dataset("extra", shape=(R, 0, NEXTRA), maxshape=(R, nmeas, NEXTRA), dtype="f8",
                              chunks=(1, min(nmeas, 4096), NEXTRA))
            ck = f.create_group("checkpoint")
            for k, v in pt.state_dict().items():
                ck[k] = v
    _check_kernels(pt)
    every = int(cfg["checkpoint_every"])
    obs = np.zeros((R, NOBS)); ext = np.zeros((R, NEXTRA))
    end = total if max_sweeps is None else min(total, pt.sweeps_done + int(max_sweeps))
    while pt.sweeps_done < end:
        n = min(every - pt.sweeps_done % every, end - pt.sweeps_done)
        t0 = time.time(); rows = []; rows_x = []
        for _ in range(n):
            pt.sweep_block(1)
            pt.exchange()
            it = pt.sweeps_done
            if it <= ntherm:
                if it % 100 == 0:
                    pt.tune()
                if it == ntherm:
                    pt.reset_swap_stats(); pt.reset_acceptance()
                    pt.round_trips = 0
            else:
                measure_all(pt.spins, pt.rep_at, pt.energies, pt.nbr, pt.p, lat.size, nrm, obs, ext)
                rows.append(obs.copy()); rows_x.append(ext.copy())
        secs += time.time() - t0
        with h5py.File(out_path, "r+") as f:
            if rows:
                for key, arr in (("obs", np.stack(rows, 1)), ("extra", np.stack(rows_x, 1))):
                    ds = f["timeseries"][key]; n0 = ds.shape[1]
                    ds.resize(n0 + arr.shape[1], axis=1); ds[:, n0:, :] = arr
            ck = f["checkpoint"]
            for k, v in pt.state_dict().items():
                if k in ck:
                    del ck[k]
                ck[k] = v
            dg = f.require_group("diagnostics")
            for k, v in dict(swap_acc=pt.swap_rates(), roundtrips=int(pt.round_trips), seconds=secs,
                             delta=pt.delta, delta_r=pt.delta_r, acceptance=pt.acceptance()).items():
                if k in dg:
                    del dg[k]
                dg[k] = v
            f.attrs["seconds"] = secs
            f.attrs["updated"] = _dt.datetime.now().isoformat(timespec="seconds")
            if pt.sweeps_done >= total:
                f.attrs["status"] = "complete"
        if verbose:
            print(f"{cfg['name']}: sweep {pt.sweeps_done}/{total} ({secs:.0f} s, min swap "
                  f"{pt.swap_rates().min():.2f}, round trips {pt.round_trips})", flush=True)
    return out_path


def _check_kernels(pt):
    """Guard against stale numba caches (a cached caller compiled before the anisotropy was added
    silently ignores it): after one sweep on a copy, the engine's running energies must equal a fresh
    energy_total with the full parameter array."""
    from polyspin.model import energy_total
    from polyspin.pt import run_sweeps
    spins = pt.spins.copy(); states = pt.states.copy(); cnt = np.zeros_like(pt.counters)
    e = np.zeros(pt.R)
    run_sweeps(spins, pt.nbr, pt.p, pt._betas_by_rep(), pt.delta[pt.temp_of], pt.delta_r[pt.temp_of],
               pt.kappa, True, True, pt.use_overlap, False, states, cnt, 1, e)
    ref = np.array([energy_total(spins[r], pt.nbr, pt.p) for r in range(pt.R)])
    if not np.allclose(e, ref, rtol=1e-10, atol=1e-8):
        raise RuntimeError("polyspin kernels disagree with energy_total (stale numba cache?): "
                           "delete polyspin/__pycache__/*.nb[ic] and rerun")


def _polyspin_git():
    import subprocess
    import polyspin
    d = Path(polyspin.__file__).resolve().parents[1]
    try:
        h = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=d, text=True).strip()
        tag = subprocess.run(["git", "describe", "--tags", "--exact-match"], cwd=d, text=True,
                             capture_output=True).stdout.strip()
        return h + (f" ({tag})" if tag else "")
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("configs", nargs="+"); ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    for c in a.configs:
        run_config(c, threads=a.threads)


if __name__ == "__main__":
    main()
