"""N3 analysis following docs/N3_analysis_plan.md (fixed before the production runs).

    python scripts/n3_analysis.py curves     # reweight runs/n3/{eff,full}_*.h5 -> data/n3/curves/
    python scripts/n3_analysis.py decide     # -> data/n3/n3_summary.json
    python scripts/n3_analysis.py all

Same estimators and joint per-size jackknife as scripts/n1_analysis.py.  Full-model runs store the
C3-model observables (sbaniso.fullmodel); Ups is rho_alpha / 9 = rho_psi.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import chi2 as chi2dist

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts"))
import n1_analysis as n1  # noqa: E402
from n2_analysis import jackknife  # noqa: E402
from sbaniso.analysis import curves, equilibration_check, load_run, tau_table  # noqa: E402

RUNS = REPO / "runs" / "n3"
CUR = REPO / "data" / "n3" / "curves"
OUT = REPO / "data" / "n3" / "n3_summary.json"
ZH = ("0.72", "0.80", "0.88")


def tgrid(run):
    return np.round(np.arange(np.ceil(run.Ts.min() * 2e4) / 2e4, run.Ts.max() + 1e-12, 0.00025), 6)


def _stage1(path):
    name = Path(path).stem
    run = load_run(path)
    T = tgrid(run)
    v, e, p, mh = curves(run, T)
    eq = equilibration_check(run); tau = tau_table(run)
    k = int(np.nanargmax(v["C"]))
    lw = mh.log_weights(1 / T[k]); w = np.exp(lw - lw.max())
    h, edges = np.histogram(mh.Ef, bins=60, weights=w)
    extra = {}
    with h5py.File(path) as f:
        if "extra" in f["timeseries"]:
            x = f["timeseries/extra"][()]
            extra = dict(K_by_T=x[..., 0].mean(1), absq30_by_T=x[..., 1].mean(1))
    np.savez_compressed(CUR / f"{name}.npz", T=T, L=run.L, Ts=run.Ts,
                        **{f"v_{q}": v[q] for q in n1.KEYS}, **{f"e_{q}": e[q] for q in n1.KEYS},
                        **{f"p_{q}": p[q] for q in n1.KEYS}, hist=h, hist_edges=edges, T_Cmax=T[k], **extra)
    eq["tau_by_T"] = {q: t.tolist() for q, t in tau.items()}
    eq["swap_acc"] = run.diagnostics["swap_acc"].tolist()
    return name, eq


def stage1():
    CUR.mkdir(parents=True, exist_ok=True)
    files = sorted(str(f) for f in RUNS.glob("*.h5") if f.name.startswith(("eff_", "full_")) and "ordered" not in f.name)
    with Pool(4) as pool:
        res = pool.map(_stage1, files, chunksize=1)
    (CUR / "equilibration.json").write_text(json.dumps(dict(res), indent=1, default=float))
    print(f"stage 1: {len(res)} runs")


def load_family(prefix, zh):
    files = sorted(CUR.glob(f"{prefix}_zH4_{zh}_L*.npz"), key=lambda f: int(f.stem.split("_L")[1]))
    return [n1.Curves(f) for f in files], files


def transitions(T, cs, cov):
    pairs = list(zip(cs, cs[1:]))
    tc = np.array([n1.pair_crossing(T, a, b)[0] for a, b in pairs])
    uc = np.array([n1.pair_crossing(T, a, b)[1] for a, b in pairs])
    Lbar = np.array([np.sqrt(a.L * b.L) for a, b in pairs], float)
    good = np.isfinite(tc)
    TI = np.nan
    if good.sum() >= 2:
        cg = cov[np.ix_(good, good)]
        cg = cg + np.eye(len(cg)) * max(1e-12, 1e-9 * np.trace(cg))   # guard against singular cov
        try:
            TI = n1.gls_extrapolate(Lbar[good], tc[good], cg, n1.X_MAIN)[0]
        except np.linalg.LinAlgError:
            TI = np.nan
    TKT, chi = n1.wm_T(T, cs)
    ups_cross = np.array([n1.crossing(T, c.v["Ups"] - 2 * T / np.pi) for c in cs])
    Tc_C = np.array([T[int(np.nanargmax(c.v["C"]))] for c in cs])
    Cmax = np.array([np.nanmax(c.v["C"]) for c in cs])
    U4min = np.array([np.nanmin(c.v["U4"]) for c in cs])
    return dict(TI=TI, TKT=TKT, dT=TI - TKT, T_cross=tc, U_cross=uc, WM_chi2dof=chi, ups_cross=ups_cross,
                Tc_C=Tc_C, Cmax=Cmax, U4min=U4min, Tcross_mean=np.nanmean(tc))


def bimod(c):
    """(dip ratio, barrier dF = ln(peak / dip)); dF only defined for a bimodal histogram (ratio < 0.8)."""
    h = np.convolve(c.hist.astype(float), np.ones(3) / 3, mode="same")
    r = n1.bimodality_hist(c.hist.astype(float))
    return r, (float(-np.log(r)) if r < 0.8 else None)


def family(prefix, zh):
    cs, files = load_family(prefix, zh)
    if not cs:
        return None
    T = cs[0].T
    cov, _ = n1.crossing_cov(T, cs)
    cen, err = jackknife(lambda x: transitions(T, x, cov), cs)
    d = dict(L=[c.L for c in cs])
    for k, v in cen.items():
        d[k] = v.tolist() if hasattr(v, "tolist") else float(v)
        d["d" + k] = err[k].tolist() if hasattr(err[k], "tolist") else float(err[k])
    L = np.array(d["L"], float); C = np.array(cen["Cmax"]); dC = np.array(err["Cmax"])
    fits = {}
    for tag, xv in (("lnL", np.log(L)), ("L2", L ** 2)):
        a, b, _, ch, dof = n1.fit_line(xv, C, dC)
        fits[tag] = dict(a=a, b=b, chi2=ch, dof=dof, p=float(chi2dist.sf(ch, dof)) if dof > 0 else None)
    d["Cmax_fits"] = fits; d["Cmax_over_L2"] = (C / L ** 2).tolist()
    bm = [bimod(c) for c in cs]
    d["bimodality_at_Cmax"] = [b[0] for b in bm]; d["dF_at_Cmax"] = [b[1] for b in bm]
    if prefix == "full":
        d["K_by_T"] = {str(c.L): dict(T=np.load(f)["Ts"].tolist(), K=np.load(f)["K_by_T"].tolist()) for c, f in zip(cs, files)}
    return d


def stage2():
    eqs = json.loads((CUR / "equilibration.json").read_text())
    S = dict(plan="docs/N3_analysis_plan.md", points={})
    for zh in ZH:
        P = {}
        for prefix in ("eff", "full"):
            d = family(prefix, zh)
            if d is not None:
                d["equilibration"] = {k: {q: eqs[k][q] for q in ("pass", "reasons", "quarter_z", "quarter_threshold", "tau_max",
                                                                 "roundtrips", "swap_min")} for k in eqs if k.startswith(f"{prefix}_zH4_{zh}_")}
                P[prefix] = d
        e, f = P["eff"], P["full"]
        # (1) TASK rule
        same_sign = bool(np.sign(e["dT"]) == np.sign(f["dT"]))
        zTI = abs(f["TI"] - e["TI"]) / np.hypot(f["dTI"], e["dTI"])
        zTK = abs(f["TKT"] - e["TKT"]) / np.hypot(f["dTKT"], e["dTKT"])
        # (2) coincidence at the largest L
        t_u = f["Tcross_mean"]; dt_u = f["dTcross_mean"]
        t_y = f["ups_cross"][-1]; dt_y = f["dups_cross"][-1]
        t_c = f["Tc_C"][-1]; dt_c = f["dTc_C"][-1]
        trio = [(t_u, dt_u), (t_y, dt_y), (t_c, dt_c)]
        coincide = all(abs(a[0] - b[0]) <= 3 * np.hypot(a[1], b[1]) + 1e-12 for i, a in enumerate(trio) for b in trio[i + 1:])
        # (3) first order
        dF = f["dF_at_Cmax"]
        fo_hist = bool(f["bimodality_at_Cmax"][-1] < 0.8 and all(x is not None for x in dF) and np.all(np.diff(dF) > 0))
        fo_C = bool(f["Cmax_fits"]["L2"]["p"] is not None and f["Cmax_fits"]["L2"]["p"] >= 0.01
                    and f["Cmax_fits"]["lnL"]["p"] < 0.01)
        fo_U = bool(f["U4min"][-1] < 0 and np.all(np.diff(f["U4min"]) < 0))
        P["decision"] = dict(
            task_rule=dict(same_order=same_sign, z_TI=zTI, z_TKT=zTK,
                           agree=bool(same_sign and zTI <= 3 and zTK <= 3)),
            coincident=dict(T_U4cross=t_u, T_ups=t_y, T_Cmax=t_c, errs=[dt_u, dt_y, dt_c], coincide=bool(coincide)),
            first_order=dict(hist=fo_hist, Cmax_L2=fo_C, U4_negative=fo_U, first_order=bool(fo_hist or fo_C)))
        S["points"][zh] = P
        print(f"zH4={zh}: EFF T_I={e['TI']:.5f}({e['dTI']:.5f}) T_KT={e['TKT']:.5f}({e['dTKT']:.5f}) dT={e['dT']:+.5f}")
        print(f"          FULL T_I={f['TI']:.5f}({f['dTI']:.5f}) T_KT={f['TKT']:.5f}({f['dTKT']:.5f}) dT={f['dT']:+.5f}  "
              f"Tc_C={np.round(f['Tc_C'], 5).tolist()} ups_cross={np.round(f['ups_cross'], 5).tolist()}")
        print(f"          Cmax={np.round(f['Cmax'], 2).tolist()} C/L2={np.round(f['Cmax_over_L2'], 4).tolist()} "
              f"fits lnL p={f['Cmax_fits']['lnL']['p']:.3g} L2 p={f['Cmax_fits']['L2']['p']:.3g} "
              f"bimod={np.round(f['bimodality_at_Cmax'], 2).tolist()} dF={f['dF_at_Cmax']} U4min={np.round(f['U4min'], 3).tolist()}")
        print(f"          decision: {P['decision']}")
    a = load_run(RUNS / "full_zH4_0.80_L32.h5"); b = load_run(RUNS / "full_zH4_0.80_L32_ordered.h5")
    S["ordered_vs_random_zH4_0.80_L32"] = dict(n1.compare_runs(a, b), ordered_equilibration=equilibration_check(b))
    print("ordered vs random:", {k: v for k, v in S["ordered_vs_random_zH4_0.80_L32"].items() if k != "ordered_equilibration"})
    OUT.write_text(json.dumps(S, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else (bool(o) if isinstance(o, np.bool_) else float(o))))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("curves", "all"):
        stage1()
    if what in ("decide", "all"):
        stage2()
