"""N1 analysis following docs/N1_analysis_plan.md (sections 1-5, fixed before computing).

    python scripts/n1_analysis.py curves     # stage 1: reweight every run -> data/n1/curves/*.npz
    python scripts/n1_analysis.py decide     # stage 2: decisions -> data/n1/n1_summary.json
    python scripts/n1_analysis.py supplementary   # labelled, not pre-registered (runs after decide)
    python scripts/n1_analysis.py all

Errors: every derived quantity Q (crossings, extrapolated T_I, T_KT, Delta T, nu, ...) is a function
of the reweighted curves of all sizes at one r.  Its error is a joint block jackknife: one size at a
time, one block at a time, that size's curves are replaced by their leave-one-block-out versions and
Q is recomputed; variances are summed over sizes (different sizes = independent runs).  This keeps
the covariance of T_I and T_KT, which come from the same runs.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2dist

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from sbaniso.analysis import (MultiHistogram, XI_L_KT_XY, _parabolic_min, _wm_scan,  # noqa: E402
                              bimodality, crossing, curves, equilibration_check, load_run,
                              per_temperature, tau_table)

RUNS = REPO / "runs" / "n1"
CUR = REPO / "data" / "n1" / "curves"
OUT = REPO / "data" / "n1" / "n1_summary.json"
TGRID = np.round(np.arange(6.2, 8.0 + 1e-9, 0.005), 6)
RS = (3.6, 3.8, 4.0, 4.1, 4.2, 4.3, 4.4, 4.6, 4.8, 5.0)
KEYS = ("U4", "Ups", "C", "xiL_psi", "xiL_phi", "xiL_sigma", "absm", "E")
X_MAIN = 2.75; X_ALT = 3.0; EPS_WM = 0.0025
U4_WINDOW = (0.25, 0.66)


# ----------------------------------------------------------------------------- stage 1


def _stage1(path):
    name = Path(path).stem
    run = load_run(path)
    v, e, p, mh = curves(run, TGRID)
    eq = equilibration_check(run)
    tau = tau_table(run)
    # energy histogram at the C maximum (reweighted)
    k = int(np.nanargmax(v["C"]))
    lw = mh.log_weights(1 / TGRID[k]); w = np.exp(lw - lw.max())
    h, edges = np.histogram(mh.Ef, bins=60, weights=w)
    np.savez_compressed(CUR / f"{name}.npz", T=TGRID, L=run.L, Ts=run.Ts,
                        **{f"v_{q}": v[q] for q in KEYS}, **{f"e_{q}": e[q] for q in KEYS},
                        **{f"p_{q}": p[q] for q in KEYS},
                        hist=h, hist_edges=edges, T_Cmax=TGRID[k],
                        **{f"tau_{q}": t for q, t in tau.items()})
    eq["tau_by_T"] = {q: t.tolist() for q, t in tau.items()}
    return name, eq


def stage1():
    CUR.mkdir(parents=True, exist_ok=True)
    files = sorted(str(f) for f in RUNS.glob("*.h5") if "ordered" not in f.name)
    files.sort(key=lambda f: -load_run_size(f))
    with Pool(4) as pool:
        res = pool.map(_stage1, files, chunksize=1)
    (CUR / "equilibration.json").write_text(json.dumps(dict(res), indent=1, default=float))
    print(f"stage 1: {len(res)} runs")


def load_run_size(f):
    import h5py
    with h5py.File(f) as h:
        return int(h["params/L"][()]) * len(h["params/Ts"][()])


# ----------------------------------------------------------------------------- stage 2 helpers


class Curves:
    """Reweighted curves of one size: central values v, errors e, pseudo-values p (nb, nT)."""

    def __init__(self, f):
        d = np.load(f)
        self.L = int(d["L"]); self.T = d["T"]
        self.v = {q: d[f"v_{q}"] for q in KEYS}; self.e = {q: d[f"e_{q}"] for q in KEYS}
        self.p = {q: d[f"p_{q}"] for q in KEYS}
        self.hist = d["hist"]; self.T_Cmax = float(d["T_Cmax"])
        self.nb = self.p["U4"].shape[0]

    def replica(self, b):
        """Same object with central values replaced by leave-block-b-out values."""
        c = Curves.__new__(Curves)
        c.L, c.T, c.e, c.nb = self.L, self.T, self.e, self.nb
        c.v = {q: self.p[q][b] for q in KEYS}
        return c


def pair_crossing(T, c1, c2, key="U4", window=U4_WINDOW):
    u1, u2 = c1.v[key], c2.v[key]
    ok = (u1 > window[0]) & (u1 < window[1])
    t = crossing(T, np.where(ok, u1 - u2, np.nan))
    return t, (np.interp(t, T, u1) if np.isfinite(t) else np.nan)


def wm_T(T, cs):
    Ls = [c.L for c in cs]
    prof, _ = _wm_scan(T, Ls, [c.v["Ups"] for c in cs], [c.e["Ups"] for c in cs],
                       np.linspace(-2 * np.log(min(Ls)) + 0.2, 60, 2000))
    return _parabolic_min(T, prof), float(np.nanmin(prof))


def gls_extrapolate(Lbar, y, cov, x):
    """GLS fit y = T_inf + a Lbar^-x; returns (T_inf, a, chi2, dof)."""
    A = np.vstack([np.ones_like(Lbar), Lbar ** -x]).T
    W = np.linalg.pinv(cov)
    M = A.T @ W @ A
    sol = np.linalg.solve(M, A.T @ W @ y)
    res = y - A @ sol
    return sol[0], sol[1], float(res @ W @ res), len(y) - 2


def quantities(T, cs, cov_cross=None):
    """All per-r quantities from a list of Curves (sorted by L).  cov_cross: covariance of the pair
    crossings (needed for the GLS extrapolation; None -> diagonal 1)."""
    out = {}
    pairs = list(zip(cs, cs[1:]))
    tc = np.array([pair_crossing(T, a, b)[0] for a, b in pairs])
    uc = np.array([pair_crossing(T, a, b)[1] for a, b in pairs])
    out["T_cross"] = tc; out["U_cross"] = uc
    Lbar = np.array([np.sqrt(a.L * b.L) for a, b in pairs], float)
    good = np.isfinite(tc)
    cov = cov_cross if cov_cross is not None else np.eye(len(tc))
    cov = cov[np.ix_(good, good)]
    for tag, x in (("TI", X_MAIN), ("TI_x3", X_ALT)):
        t0, a, ch, dof = gls_extrapolate(Lbar[good], tc[good], cov, x) if good.sum() >= 2 else (np.nan,) * 4
        out[tag] = t0
    # x free: grid over x, GLS inside (needs >= 3 points)
    if good.sum() >= 3:
        xs = np.linspace(0.5, 6, 221)
        chis = [gls_extrapolate(Lbar[good], tc[good], cov, x)[2] for x in xs]
        xb = xs[int(np.argmin(chis))]
        out["TI_xfree"] = gls_extrapolate(Lbar[good], tc[good], cov, xb)[0]; out["x_free"] = xb
    else:
        out["TI_xfree"] = np.nan; out["x_free"] = np.nan
    out["TKT"], out["WM_chi2dof"] = wm_T(T, cs)
    out["dT"] = out["TI"] - out["TKT"]
    out["dT_corr"] = out["TI"] - (1 + EPS_WM) * out["TKT"]
    # xi_psi/L crossings (auxiliary)
    out["xi_T"] = np.array([pair_crossing(T, a, b, "xiL_psi", (0.3, 1.5))[0] for a, b in pairs])
    out["xi_val"] = np.array([pair_crossing(T, a, b, "xiL_psi", (0.3, 1.5))[1] for a, b in pairs])
    # nu from |dU4/dT| at T_I
    lnL = np.log([c.L for c in cs]); sl = []
    for c in cs:
        g = np.gradient(c.v["U4"], T)
        sl.append(abs(np.interp(out["TI"], T, g)) if np.isfinite(out["TI"]) else np.nan)
    sl = np.array(sl)
    if np.all(np.isfinite(sl)) and np.all(sl > 0):
        k = np.polyfit(lnL, np.log(sl), 1)[0]
        out["nu"] = 1 / k
    else:
        out["nu"] = np.nan
    out["Cmax"] = np.array([np.nanmax(c.v["C"]) for c in cs])
    return out


def jackknife_all(T, cs, cov_cross=None):
    """Central quantities and joint jackknife errors (see module docstring)."""
    central = quantities(T, cs, cov_cross)
    var = {k: np.zeros_like(np.atleast_1d(np.asarray(v, float))) for k, v in central.items()}
    reps_store = {}
    for i, c in enumerate(cs):
        reps = []
        for b in range(c.nb):
            cs2 = list(cs); cs2[i] = c.replica(b)
            reps.append(quantities(T, cs2, cov_cross))
        for k in central:
            arr = np.array([np.atleast_1d(np.asarray(r[k], float)) for r in reps])
            m = np.nanmean(arr, 0)
            var[k] += (c.nb - 1) / c.nb * np.nansum((arr - m) ** 2, 0)
        reps_store[c.L] = reps
    err = {k: np.sqrt(v) if np.ndim(central[k]) else float(np.sqrt(v[0])) for k, v in var.items()}
    return central, err, reps_store


def crossing_cov(T, cs):
    """Covariance of the pair crossings from the per-size jackknife."""
    pairs = len(cs) - 1
    cov = np.zeros((pairs, pairs))
    base = np.array([pair_crossing(T, a, b)[0] for a, b in zip(cs, cs[1:])])
    for i, c in enumerate(cs):
        arr = []
        for b in range(c.nb):
            cs2 = list(cs); cs2[i] = c.replica(b)
            arr.append([pair_crossing(T, a, bb)[0] for a, bb in zip(cs2, cs2[1:])])
        arr = np.array(arr); d = arr - np.nanmean(arr, 0)
        d = np.nan_to_num(d)
        cov += (c.nb - 1) / c.nb * d.T @ d
    return cov, base


def fit_line(x, y, s):
    """Weighted straight line y = a + b x; returns (a, b, cov, chi2, dof)."""
    A = np.vstack([np.ones_like(x), x]).T / s[:, None]
    sol, *_ = np.linalg.lstsq(A, y / s, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    ch = float((((A @ sol) - y / s) ** 2).sum())
    return sol[0], sol[1], cov, ch, len(x) - 2


def zero_of_line(a, b, cov):
    r0 = -a / b
    g = np.array([-1 / b, a / b ** 2])
    return r0, float(np.sqrt(g @ cov @ g))


def compare_runs(a, b, names=("E", "U4", "Ups", "C", "absm")):
    pa, pb = per_temperature(a), per_temperature(b)
    zs = np.concatenate([(pa[k][0] - pb[k][0]) / np.hypot(pa[k][1], pb[k][1]) for k in names])
    return dict(max_abs_z=float(np.abs(zs).max()), chi2dof=float(np.mean(zs ** 2)), n=len(zs))


# ----------------------------------------------------------------------------- stage 2


def stage2():
    eqs = json.loads((CUR / "equilibration.json").read_text())
    S = dict(plan="docs/N1_analysis_plan.md", x_main=X_MAIN, eps_wm=EPS_WM, per_r={})
    rows = []
    for r in RS:
        files = sorted(CUR.glob(f"r{r}_L*.npz"), key=lambda f: int(f.stem.split("_L")[1]))
        cs = [Curves(f) for f in files]
        T = cs[0].T
        cov, base = crossing_cov(T, cs)
        cen, err, _ = jackknife_all(T, cs, cov)
        Ls = [c.L for c in cs]
        d = dict(L=Ls, pairs=[[a.L, b.L] for a, b in zip(cs, cs[1:])])
        for k, v in cen.items():
            d[k] = v.tolist() if np.ndim(v) else float(v)
            d["d" + k] = err[k].tolist() if np.ndim(err[k]) else float(err[k])
        d["crossing_cov"] = cov.tolist()
        # decision
        if r <= 4.0:
            d["rule"] = "T_I - 1.0025 T_KT > 3 sigma"
            d["z_rule"] = cen["dT_corr"] / err["dT_corr"]; d["pass"] = bool(d["z_rule"] > 3)
        elif r >= 4.6:
            d["rule"] = "T_KT - T_I > 3 sigma"
            d["z_rule"] = -cen["dT"] / err["dT"]; d["pass"] = bool(d["z_rule"] > 3)
        else:
            d["rule"] = None; d["z_rule"] = None; d["pass"] = None
        # specific heat: ln L vs L^2 fits
        C = np.asarray(cen["Cmax"]); dC = np.asarray(err["Cmax"]); Lf = np.array(Ls, float)
        fits = {}
        for tag, xv in (("lnL", np.log(Lf)), ("L2", Lf ** 2)):
            a, b, _, ch, dof = fit_line(xv, C, dC)
            fits[tag] = dict(a=a, b=b, chi2=ch, dof=dof, p=float(chi2dist.sf(ch, dof)) if dof > 0 else None)
        d["Cmax_fits"] = fits
        d["Cmax_over_L2"] = (C / Lf ** 2).tolist()
        d["bimodality"] = {c.L: float(bimodality_hist(c.hist)) for c in cs}
        d["T_Cmax"] = {c.L: c.T_Cmax for c in cs}
        d["equilibration"] = {f.stem: {k: eqs[f.stem][k] for k in ("pass", "quarter_z", "quarter_threshold",
                                                                  "tau_max", "roundtrips", "swap_min")}
                              for f in files}
        d["tau_max_by_L"] = {f.stem: max(max(v) for v in eqs[f.stem]["tau_by_T"].values()) for f in files}
        S["per_r"][str(r)] = d
        rows.append((r, cen["TI"], err["TI"], cen["TKT"], err["TKT"], cen["dT"], err["dT"]))
        print(f"r={r}: T_I={cen['TI']:.4f}({err['TI']:.4f})  T_KT={cen['TKT']:.4f}({err['TKT']:.4f})  "
              f"dT={cen['dT']:+.4f}({err['dT']:.4f})  rule={d['rule']} z={d['z_rule']}  pass={d['pass']}",
              flush=True)
    rr = np.array([x[0] for x in rows]); TI = np.array([x[1] for x in rows]); dTI = np.array([x[2] for x in rows])
    TK = np.array([x[3] for x in rows]); dTK = np.array([x[4] for x in rows])
    dT = np.array([x[5] for x in rows]); ddT = np.array([x[6] for x in rows])

    def crossing_point(sel, TKscale):
        dd = TI[sel] - TKscale * TK[sel]
        a, b, cov, ch, dof = fit_line(rr[sel], dd, ddT[sel])
        r0, dr0 = zero_of_line(a, b, cov)
        aI, bI, cI, _, _ = fit_line(rr[sel], TI[sel], dTI[sel])
        aK, bK, cK, _, _ = fit_line(rr[sel], TKscale * TK[sel], TKscale * dTK[sel])
        Tx = 0.5 * ((aI + bI * r0) + (aK + bK * r0))
        v = np.array([1, r0])
        dTx = 0.5 * np.sqrt(v @ cI @ v + v @ cK @ v + ((bI + bK) * dr0) ** 2)
        return dict(r_x=r0, dr_x=dr0, T_x=Tx, dT_x=dTx, chi2=ch, dof=dof)

    allr = np.ones(len(rr), bool); mid = (rr >= 4.0) & (rr <= 4.6)
    main = crossing_point(allr, 1.0); corr = crossing_point(allr, 1 + EPS_WM)
    S["crossing"] = dict(main=main, wm_corrected=corr, syst_r=abs(corr["r_x"] - main["r_x"]),
                         syst_T=abs(corr["T_x"] - main["T_x"]),
                         window_4p0_4p6=crossing_point(mid, 1.0),
                         prediction=dict(r_x=4.31, T_x=7.07))
    # ordered vs random
    S["ordered_vs_random"] = {}
    for r in (4.2, 4.4):
        a = load_run(RUNS / f"r{r}_L64.h5"); b = load_run(RUNS / f"r{r}_L64_ordered.h5")
        S["ordered_vs_random"][str(r)] = dict(compare_runs(a, b), ordered_equilibration=equilibration_check(b)["pass"])
    OUT.write_text(json.dumps(S, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    c = S["crossing"]
    print(f"r_x = {c['main']['r_x']:.3f} +- {c['main']['dr_x']:.3f} (stat) +- {c['syst_r']:.3f} (syst);  "
          f"T_x = {c['main']['T_x']:.4f} +- {c['main']['dT_x']:.4f} (stat) +- {c['syst_T']:.4f} (syst); "
          f"line chi2/dof = {c['main']['chi2']:.1f}/{c['main']['dof']}")
    print("window 4.0-4.6:", {k: round(v, 4) for k, v in c["window_4p0_4p6"].items()})
    print("ordered vs random:", S["ordered_vs_random"])


def bimodality_hist(h):
    """bimodality() works on samples; this is the same dip ratio on a given (weighted) histogram."""
    h = np.convolve(h, np.ones(3) / 3, mode="same")
    i1 = int(np.argmax(h)); best = 1.0
    for j in range(len(h)):
        if abs(j - i1) < 4:
            continue
        lo, hi = sorted((i1, j))
        if h[j] > 0.2 * h[i1] and h[j] == h[max(0, j - 2):j + 3].max():
            best = min(best, h[lo:hi + 1].min() / h[j])
    return best




# ----------------------------------------------------------------------------- supplementary
# NOT pre-registered (docs/N1_report.md says so): added after the pre-registered linear fit of
# Delta T(r) turned out to be a poor model (chi2/dof ~ 66) and the xi_psi/L pair crossings scattered.


def supplementary():
    S = json.loads(OUT.read_text())
    P = S["per_r"]; rs = np.array(sorted(float(k) for k in P))
    dT = np.array([P[str(r)]["dT"] for r in rs]); ddT = np.array([P[str(r)]["ddT"] for r in rs])
    TI = np.array([P[str(r)]["TI"] for r in rs]); TK = np.array([P[str(r)]["TKT"] for r in rs])
    sup = {}
    m = S["crossing"]["main"]
    sup["linear_scaled_error"] = dict(r_x=m["r_x"], dr_x=m["dr_x"] * np.sqrt(m["chi2"] / m["dof"]),
                                      T_x=m["T_x"], dT_x=m["dT_x"] * np.sqrt(m["chi2"] / m["dof"]))
    # bracketing interpolation between the last positive and first negative Delta T
    def bracket(scale):
        d = TI - scale * TK
        i = int(np.where((d[:-1] > 0) & (d[1:] < 0))[0][0])
        f = d[i] / (d[i] - d[i + 1]); r0 = rs[i] + f * (rs[i + 1] - rs[i])
        dr = (rs[i + 1] - rs[i]) * np.hypot(ddT[i] * d[i + 1], ddT[i + 1] * d[i]) / (d[i] - d[i + 1]) ** 2
        Tx = TI[i] + f * (TI[i + 1] - TI[i])
        return r0, dr, Tx
    r0, dr, Tx = bracket(1.0); r1, _, Tx1 = bracket(1 + EPS_WM)
    sup["bracket_4p3_4p4"] = dict(r_x=r0, dr_x=dr, syst_r=abs(r1 - r0), T_x=Tx, syst_T=abs(Tx1 - Tx))
    # quadratic fit of Delta T over all r
    A = np.vstack([np.ones_like(rs), rs, rs ** 2]).T / ddT[:, None]
    sol, *_ = np.linalg.lstsq(A, dT / ddT, rcond=None)
    ch = float((((A @ sol) - dT / ddT) ** 2).sum())
    roots = np.roots(sol[::-1]); roots = roots[(roots.real > 4.0) & (roots.real < 4.6)].real
    sup["quadratic_all_r"] = dict(r_x=float(roots[0]) if len(roots) else None, chi2=ch, dof=len(rs) - 3)
    # xi_psi/L: value at T_KT for each L, and how distinguishable the sizes are below T_KT
    xi = {}
    for r in rs:
        files = sorted(CUR.glob(f"r{r}_L*.npz"), key=lambda f: int(f.stem.split("_L")[1]))
        cs = [Curves(f) for f in files]; T = cs[0].T; tk = P[str(r)]["TKT"]
        at = {c.L: (float(np.interp(tk, T, c.v["xiL_psi"])), float(np.interp(tk, T, c.e["xiL_psi"]))) for c in cs}
        below = (T > tk - 0.4) & (T < tk)
        zs = []
        for a, b in zip(cs, cs[1:]):
            z = (a.v["xiL_psi"] - b.v["xiL_psi"]) / np.hypot(a.e["xiL_psi"], b.e["xiL_psi"])
            zs.append(float(np.nanmedian(np.abs(z[below]))))
        xi[str(r)] = dict(xiL_at_TKT=at, median_abs_z_pairs_below_TKT=zs)
    sup["xi_psi"] = xi
    S["supplementary_not_preregistered"] = sup
    OUT.write_text(json.dumps(S, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    print(json.dumps({k: v for k, v in sup.items() if k != "xi_psi"}, indent=1, default=float))
    for r, v in xi.items():
        print(r, {L: round(x[0], 3) for L, x in v["xiL_at_TKT"].items()}, "err~", round(max(x[1] for x in v["xiL_at_TKT"].values()), 3),
              "median|z| below T_KT:", np.round(v["median_abs_z_pairs_below_TKT"], 2).tolist())


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("curves", "all"):
        stage1()
    if what in ("decide", "all"):
        stage2()
    if what in ("supplementary", "decide", "all"):
        supplementary()
