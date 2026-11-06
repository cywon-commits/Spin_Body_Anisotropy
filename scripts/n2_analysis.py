"""N2 analysis following docs/N2_analysis_plan.md (fixed before the production runs).

    python scripts/n2_analysis.py curves    # reweight runs/n2/r*_L*.h5 -> data/n2/curves/
    python scripts/n2_analysis.py decide    # criteria (A)-(D) -> data/n2/n2_summary.json
    python scripts/n2_analysis.py all

Errors: joint per-size block jackknife, as in scripts/n1_analysis.py.  The three interior
temperatures T_q = T_I + q (T_KT - T_I), q = 1/4, 1/2, 3/4, are set from the central values and held
fixed in the jackknife replicas (plan, section 2B).
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts"))
import n1_analysis as n1  # noqa: E402
from sbaniso.analysis import curves, equilibration_check, load_run, tau_table  # noqa: E402

n1.KEYS = n1.KEYS + ("Mpsi2", "Mphi2", "m2")
KEYS = n1.KEYS
RUNS = REPO / "runs" / "n2"
CUR = REPO / "data" / "n2" / "curves"
OUT = REPO / "data" / "n2" / "n2_summary.json"
RANGES = {5.0: (6.75, 8.45), 6.0: (7.0, 9.75)}
QS = (0.25, 0.5, 0.75)


def grid(r):
    a, b = RANGES[r]
    return np.round(np.arange(a, b + 1e-9, 0.005), 6)


def _stage1(path):
    name = Path(path).stem
    run = load_run(path)
    r = float(name.split("_")[0][1:])
    v, e, p, mh = curves(run, grid(r))
    eq = equilibration_check(run); tau = tau_table(run)
    np.savez_compressed(CUR / f"{name}.npz", T=grid(r), L=run.L, Ts=run.Ts,
                        **{f"v_{q}": v[q] for q in KEYS}, **{f"e_{q}": e[q] for q in KEYS},
                        **{f"p_{q}": p[q] for q in KEYS}, hist=np.zeros(1), T_Cmax=np.nan)
    eq["tau_by_T"] = {q: t.tolist() for q, t in tau.items()}
    return name, eq


def stage1():
    CUR.mkdir(parents=True, exist_ok=True)
    files = sorted(str(f) for f in RUNS.glob("r*_L*.h5") if "ordered" not in f.name)
    with Pool(4) as pool:
        res = pool.map(_stage1, files, chunksize=1)
    (CUR / "equilibration.json").write_text(json.dumps(dict(res), indent=1, default=float))
    print(f"stage 1: {len(res)} runs")


def at(c, key, t):
    return float(np.interp(t, c.T, c.v[key]))


def transitions(T, cs, cov):
    pairs = list(zip(cs, cs[1:]))
    tc = np.array([n1.pair_crossing(T, a, b)[0] for a, b in pairs])
    Lbar = np.array([np.sqrt(a.L * b.L) for a, b in pairs], float)
    good = np.isfinite(tc)
    TI = n1.gls_extrapolate(Lbar[good], tc[good], cov[np.ix_(good, good)], n1.X_MAIN)[0]
    TKT, chi = n1.wm_T(T, cs)
    return dict(TI=TI, TKT=TKT, window=TKT - TI, T_cross=tc, WM_chi2dof=chi)


def at_Tq(cs, Tq):
    """Criteria (B), (C) at fixed temperatures Tq."""
    lnL = np.log([c.L for c in cs]); N = np.array([c.L ** 2 for c in cs], float)
    o = {}
    for k, t in enumerate(Tq):
        U = np.array([at(c, "U4", t) for c in cs])
        mg = np.array([at(c, "Ups", t) - 2 * t / np.pi for c in cs])
        chi_psi = N * np.array([at(c, "Mpsi2", t) for c in cs])
        chi_phi = N * np.array([at(c, "Mphi2", t) for c in cs])
        o[f"U4_{k}"] = U
        o[f"U4diff_{k}"] = U[0] - U[-1]
        o[f"U4ratio_{k}"] = U[-1] / U[0]
        o[f"ups_margin_{k}"] = mg
        o[f"ups_margin_min_{k}"] = mg.min()
        o[f"eta_psi_{k}"] = 2 - np.polyfit(lnL, np.log(chi_psi), 1)[0]
        o[f"eta_phi_{k}"] = 2 - np.polyfit(lnL, np.log(chi_phi), 1)[0]
        o[f"eta_gap_{k}"] = o[f"eta_phi_{k}"] - o[f"eta_psi_{k}"]
        o[f"xiL_phi_{k}"] = np.array([at(c, "xiL_phi", t) for c in cs])
        o[f"xiL_psi_{k}"] = np.array([at(c, "xiL_psi", t) for c in cs])
    return o


def jackknife(fn, cs):
    central = fn(cs)
    var = {k: np.zeros_like(np.atleast_1d(np.asarray(v, float))) for k, v in central.items()}
    for i, c in enumerate(cs):
        reps = []
        for b in range(c.nb):
            cs2 = list(cs); cs2[i] = c.replica(b)
            reps.append(fn(cs2))
        for k in central:
            arr = np.array([np.atleast_1d(np.asarray(r[k], float)) for r in reps])
            var[k] += (c.nb - 1) / c.nb * np.nansum((arr - np.nanmean(arr, 0)) ** 2, 0)
    err = {k: (np.sqrt(v) if np.ndim(central[k]) else float(np.sqrt(v[0]))) for k, v in var.items()}
    return central, err


def tolist(x):
    return x.tolist() if hasattr(x, "tolist") else x


def stage2():
    eqs = json.loads((CUR / "equilibration.json").read_text())
    S = dict(plan="docs/N2_analysis_plan.md", per_r={})
    for r in RANGES:
        files = sorted(CUR.glob(f"r{r}_L*.npz"), key=lambda f: int(f.stem.split("_L")[1]))
        cs = [n1.Curves(f) for f in files]; T = cs[0].T
        cov, _ = n1.crossing_cov(T, cs)
        tr, dtr = jackknife(lambda x: transitions(T, x, cov), cs)
        Tq = [tr["TI"] + q * (tr["TKT"] - tr["TI"]) for q in QS]
        q, dq = jackknife(lambda x: at_Tq(x, Tq), cs)
        d = dict(L=[c.L for c in cs], Tq=Tq)
        for k, v in {**tr, **q}.items():
            d[k] = tolist(v); d["d" + k] = tolist({**dtr, **dq}[k])
        # criteria
        A = bool(tr["window"] > 3 * dtr["window"])
        B, C = [], []
        for k in range(3):
            b1 = bool(q[f"U4diff_{k}"] > 3 * dq[f"U4diff_{k}"] and q[f"U4ratio_{k}"] < 0.5)
            mg, dmg = q[f"ups_margin_{k}"], dq[f"ups_margin_{k}"]
            j = int(np.argmin(mg))
            b2 = bool(np.all(mg > 0) and mg[j] > 3 * dmg[j])
            B.append(dict(sigma_disorder=b1, psi_quasi_order=b2))
            ep, dep = q[f"eta_psi_{k}"], dq[f"eta_psi_{k}"]
            c1 = bool(ep >= 0 - 3 * dep and ep <= 0.25 + 3 * dep)
            c2 = bool(q[f"eta_gap_{k}"] > 3 * dq[f"eta_gap_{k}"] and q[f"eta_phi_{k}"] > 1)
            C.append(dict(eta_psi_KT=c1, phi_exponential=c2))
        eqr = {f.stem: eqs[f.stem]["pass"] for f in files}
        d["criteria"] = dict(A_window=A, B=B, C=C, D_equilibration=all(eqr.values()))
        d["equilibration"] = {f.stem: {k: eqs[f.stem][k] for k in ("pass", "quarter_z", "quarter_threshold",
                                                                  "tau_max", "roundtrips", "swap_min")} for f in files}
        d["confirmed"] = bool(A and all(all(x.values()) for x in B + C) and all(eqr.values()))
        S["per_r"][str(r)] = d
        print(f"r={r}: T_I={tr['TI']:.4f}({dtr['TI']:.4f}) T_KT={tr['TKT']:.4f}({dtr['TKT']:.4f}) "
              f"window={tr['window']:.4f}({dtr['window']:.4f}) A={A}  Tq={np.round(Tq, 3).tolist()}")
        for k in range(3):
            print(f"   q={QS[k]}: U4={np.round(q[f'U4_{k}'], 4).tolist()} diff={q[f'U4diff_{k}']:.4f}({dq[f'U4diff_{k}']:.4f}) "
                  f"ratio={q[f'U4ratio_{k}']:.3f} | Ups-2T/pi min={q[f'ups_margin_min_{k}']:.4f}({dq[f'ups_margin_min_{k}']:.4f}) | "
                  f"eta_psi={q[f'eta_psi_{k}']:.4f}({dq[f'eta_psi_{k}']:.4f}) eta_phi={q[f'eta_phi_{k}']:.3f}({dq[f'eta_phi_{k}']:.3f}) "
                  f"| xi_phi/L={np.round(q[f'xiL_phi_{k}'], 3).tolist()} | B={B[k]} C={C[k]}")
        print(f"   confirmed={d['confirmed']}")
    a = load_run(RUNS / "r6.0_L64.h5"); b = load_run(RUNS / "r6.0_L64_ordered.h5")
    S["ordered_vs_random_r6.0_L64"] = dict(n1.compare_runs(a, b), ordered_equilibration=equilibration_check(b)["pass"])
    print("ordered vs random:", S["ordered_vs_random_r6.0_L64"])
    OUT.write_text(json.dumps(S, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("curves", "all"):
        stage1()
    if what in ("decide", "all"):
        stage2()
