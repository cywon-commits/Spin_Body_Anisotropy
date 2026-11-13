"""Coupling-scaling study (docs/N3s_analysis_plan.md, fixed before the production runs).

    python scripts/n3s_analysis.py      -> data/n3s/n3s_summary.json, printed table

Effective model: T_I (Binder, x = 2.75) and T_KT (WM) from reweighted curves with the joint
jackknife (same as N1/N3); T_eff = max(T_I, T_KT).  For s = 1 the N3 runs eff_zH4_0.80_* are used.
Full model (largest L = 24, random and ordered starts): T_Ups = last temperature where
Ups pi/(2T) crosses 1 (linear interpolation between simulated temperatures, no reweighting);
T_C = parabolic maximum of C(T); centre = mean of the two starts, error = half their width.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts"))
import n1_analysis as n1  # noqa: E402
from n2_analysis import jackknife  # noqa: E402
from sbaniso.analysis import curves, equilibration_check, load_run, per_temperature  # noqa: E402

SS = ("1.0", "0.3", "0.1", "0.03")
P = json.loads((REPO / "data" / "n3s" / "params.json").read_text())
CUR = REPO / "data" / "n3s" / "curves"
OUT = REPO / "data" / "n3s" / "n3s_summary.json"


def eff_files(s):
    if s == "1.0":
        return [REPO / f"runs/n3/eff_zH4_0.80_L{L}.h5" for L in (16, 24, 32, 64)]
    return [REPO / f"runs/n3s/eff_s{s}_L{L}.h5" for L in (16, 24, 32, 64)]


def full_files(s, L):
    rnd = REPO / (f"runs/n3/full_zH4_0.80_L{L}.h5" if s == "1.0" else f"runs/n3s/full_s{s}_L{L}_random.h5")
    return rnd, REPO / f"runs/n3s/full_s{s}_L{L}_ordered.h5"


def eff_curves(s):
    CUR.mkdir(parents=True, exist_ok=True)
    out = []
    for f in eff_files(s):
        tgt = CUR / f"eff_s{s}_{f.stem.split('_')[-1]}.npz"
        if not tgt.exists():
            run = load_run(f)
            T = np.linspace(run.Ts.min(), run.Ts.max(), 361)
            v, e, p, _ = curves(run, T)
            np.savez_compressed(tgt, T=T, L=run.L, Ts=run.Ts, **{f"v_{q}": v[q] for q in n1.KEYS},
                                **{f"e_{q}": e[q] for q in n1.KEYS}, **{f"p_{q}": p[q] for q in n1.KEYS},
                                hist=np.zeros(1), T_Cmax=np.nan)
        out.append(n1.Curves(tgt))
    return out


def eff_transitions(T, cs, cov):
    pairs = list(zip(cs, cs[1:]))
    tc = np.array([n1.pair_crossing(T, a, b)[0] for a, b in pairs])
    Lbar = np.array([np.sqrt(a.L * b.L) for a, b in pairs], float)
    g = np.isfinite(tc)
    cg = cov[np.ix_(g, g)] + np.eye(g.sum()) * 1e-9 * max(np.trace(cov[np.ix_(g, g)]), 1e-30)
    TI = n1.gls_extrapolate(Lbar[g], tc[g], cg, n1.X_MAIN)[0]
    TKT, chi = n1.wm_T(T, cs)
    return dict(TI=TI, TKT=TKT, Teff=max(TI, TKT), T_cross=tc, WM_chi2dof=chi)


def last_down_crossing(T, y):
    i = np.where((y[:-1] > 0) & (y[1:] <= 0))[0]
    if not len(i):
        return np.nan
    i = i[-1]
    return float(T[i] + y[i] / (y[i] - y[i + 1]) * (T[i + 1] - T[i]))


def parabolic_max(T, C):
    k = int(np.argmax(C))
    if k == 0 or k == len(T) - 1:
        return float(T[k])
    a, b, c = np.polyfit(T[k - 1:k + 2], C[k - 1:k + 2], 2)
    return float(-b / (2 * a)) if a < 0 else float(T[k])


def full_point(f):
    run = load_run(f); pt = per_temperature(run)
    T = run.Ts
    with h5py.File(f) as h:
        K = h["timeseries/extra"][..., 0].mean(1)
    tu = last_down_crossing(T, pt["Ups"][0] * np.pi / (2 * T) - 1)
    tC = parabolic_max(T, pt["C"][0])
    eq = equilibration_check(run)
    return dict(T_ups=tu, T_C=tC, Ts=T.tolist(), K=K.tolist(),
                eq={k: eq[k] for k in ("pass", "reasons", "tau_max", "roundtrips", "swap_min", "quarter_z")})


def main():
    S = dict(plan="docs/N3s_analysis_plan.md", s={})
    for s in SS:
        cs = eff_curves(s); T = cs[0].T
        cov, _ = n1.crossing_cov(T, cs)
        cen, err = jackknife(lambda x: eff_transitions(T, x, cov), cs)
        d = dict(params=P[s], eff={k: (v.tolist() if hasattr(v, "tolist") else float(v)) for k, v in cen.items()},
                 eff_err={k: (v.tolist() if hasattr(v, "tolist") else float(v)) for k, v in err.items()}, full={})
        for L in (16, 24):
            r, o = full_files(s, L)
            d["full"][str(L)] = dict(random=full_point(r), ordered=full_point(o))
        F = d["full"]["24"]; Teff = cen["Teff"]; dTeff = err["Teff"]
        res = {}
        for key in ("T_ups", "T_C"):
            a, b = F["random"][key], F["ordered"][key]
            c = 0.5 * (a + b); w = abs(a - b)
            res[key] = dict(centre=c, half_width=0.5 * w, width=w, random=a, ordered=b,
                            ratio=c / Teff, ratio_err=float(np.hypot(0.5 * w / Teff, c * dTeff / Teff ** 2)))
        tu = res["T_ups"]["centre"]
        Kat = {}
        for frac in (0.9, 1.0, 1.1):
            vals = [float(np.interp(frac * tu, F[st]["Ts"], F[st]["K"])) for st in ("random", "ordered")]
            Kat[str(frac)] = float(np.mean(vals))
        d["result"] = dict(T_eff=Teff, dT_eff=dTeff, which="T_I" if cen["TI"] >= cen["TKT"] else "T_KT", **res, K_at=Kat,
                           l_Gamma=P[s]["l_Gamma"], Gamma_over_Teff=P[s]["Gamma"] / Teff)
        S["s"][s] = d
        R = d["result"]
        print(f"s={s:5s} l_G={P[s]['l_Gamma']:.2f} Gamma/Teff={R['Gamma_over_Teff']:+.2f} | eff T_I={cen['TI']:.6f} T_KT={cen['TKT']:.6f} "
              f"T_eff={Teff:.6f}({dTeff:.6f}) [{R['which']}] | full T_ups rand/ord {R['T_ups']['random']:.6f}/{R['T_ups']['ordered']:.6f} "
              f"R_ups={R['T_ups']['ratio']:.3f}+-{R['T_ups']['ratio_err']:.3f} | T_C rand/ord {R['T_C']['random']:.6f}/{R['T_C']['ordered']:.6f} "
              f"R_C={R['T_C']['ratio']:.3f}+-{R['T_C']['ratio_err']:.3f} | K(0.9,1,1.1 T_ups)={[round(v, 3) for v in Kat.values()]}")
    R3 = S["s"]["0.03"]["result"]
    ratios = [S["s"][s]["result"]["T_ups"]["ratio"] for s in SS]
    S["decision"] = dict(
        R_ups_s003=R3["T_ups"]["ratio"], K_at_Tups_s003=R3["K_at"]["1.0"],
        represents=bool(abs(R3["T_ups"]["ratio"] - 1) <= 0.05 and R3["K_at"]["1.0"] <= 0.36),
        monotonic_increase=bool(np.all(np.diff(ratios) > 0)), ratios=ratios)
    print("decision:", S["decision"])
    OUT.write_text(json.dumps(S, indent=1, default=float))


if __name__ == "__main__":
    main()
