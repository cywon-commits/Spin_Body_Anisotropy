"""N3d: Debye-Waller corrected effective model vs the full model (docs/N3d_analysis_plan.md).

    python scripts/n3d_analysis.py    -> data/n3d/n3d_summary.json, printed table

T_eff^R = max(T_I, T_KT) of the effective model with the reduced couplings (runs/n3d/eff_*), N1 estimators
with the joint jackknife.  Full model: T_Ups from data/n3s/n3s_summary.json (L = 24 and 16, both starts).
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts"))
import n1_analysis as n1  # noqa: E402
from n2_analysis import jackknife  # noqa: E402
from n3s_analysis import eff_transitions  # noqa: E402
from sbaniso.analysis import curves, load_run  # noqa: E402

SS = ("1.0", "0.3", "0.1", "0.03")
CUR = REPO / "data" / "n3d" / "curves"
OUT = REPO / "data" / "n3d" / "n3d_summary.json"
N3S = json.loads((REPO / "data" / "n3s" / "n3s_summary.json").read_text())["s"]
DW = json.loads((REPO / "data" / "n3s" / "debye_waller.json").read_text())
GAUSS = json.loads((REPO / "data" / "n3s" / "tilt_gauss.json").read_text())


def eff_curves(s):
    CUR.mkdir(parents=True, exist_ok=True)
    out = []
    for L in (16, 24, 32, 64):
        tgt = CUR / f"eff_s{s}_L{L}.npz"
        if not tgt.exists():
            run = load_run(REPO / f"runs/n3d/eff_s{s}_L{L}.h5")
            T = np.linspace(run.Ts.min(), run.Ts.max(), 361)
            v, e, p, _ = curves(run, T)
            np.savez_compressed(tgt, T=T, L=run.L, Ts=run.Ts, **{f"v_{q}": v[q] for q in n1.KEYS},
                                **{f"e_{q}": e[q] for q in n1.KEYS}, **{f"p_{q}": p[q] for q in n1.KEYS},
                                hist=np.zeros(1), T_Cmax=np.nan)
        out.append(n1.Curves(tgt))
    return out


def quantities(T, cs, cov):
    q = eff_transitions(T, cs, cov)
    q["ups_cross_L"] = np.array([n1.crossing(T, c.v["Ups"] - 2 * T / np.pi) for c in cs])
    return q


def main():
    S = dict(plan="docs/N3d_analysis_plan.md", s={})
    for s in SS:
        cs = eff_curves(s); T = cs[0].T
        cov, _ = n1.crossing_cov(T, cs)
        cen, err = jackknife(lambda x: quantities(T, x, cov), cs)
        Te, dTe = float(cen["Teff"]), float(err["Teff"])
        F = N3S[s]["full"]
        a, b = F["24"]["random"]["T_ups"], F["24"]["ordered"]["T_ups"]
        c = 0.5 * (a + b); hw = 0.5 * abs(a - b)
        R = c / Te; dR = float(np.hypot(hw / Te, c * dTe / Te ** 2))
        # same-L comparison (L = 16, 24)
        sameL = {}
        for i, L in enumerate((16, 24)):
            te = float(cen["ups_cross_L"][i]); dte = float(err["ups_cross_L"][i])
            fr, fo = F[str(L)]["random"]["T_ups"], F[str(L)]["ordered"]["T_ups"]
            sameL[str(L)] = dict(eff=te, deff=dte, full_random=fr, full_ordered=fo,
                                 ratio=0.5 * (fr + fo) / te, ratio_random=fr / te, ratio_ordered=fo / te)
        # sensitivity to the K baseline (T = 0 ground state K0 instead of 1/3), first order:
        # scale T_eff^R by the change of the factor of the coupling that dominates it
        th, thb = DW[s]["theta2"], DW[s]["theta2_GSbase"]
        if cen["TI"] >= cen["TKT"]:
            scale = (1 - 6 * thb) / (1 - 6 * th); which = "T_I (J1 factor)"
        else:
            scale = (1 - 5.5 * thb) / (1 - 5.5 * th); which = "T_KT (J4 factor)"
        d = dict(couplings_R=DW[s]["JR"], theta2=th, TI=float(cen["TI"]), dTI=float(err["TI"]),
                 TKT=float(cen["TKT"]), dTKT=float(err["TKT"]), T_eff_R=Te, dT_eff_R=dTe,
                 which="T_I" if cen["TI"] >= cen["TKT"] else "T_KT", WM_chi2dof=float(cen["WM_chi2dof"]),
                 T_ups_full=c, T_ups_full_halfwidth=hw, ratio=R, ratio_err=dR, pass_5pct=bool(abs(R - 1) <= 0.05),
                 ratio_unreduced=N3S[s]["result"]["T_ups"]["ratio"], sameL=sameL,
                 sensitivity_K0=dict(theta2_K0=thb, scale=scale, via=which, ratio_K0=R / scale),
                 gauss=GAUSS[s])
        S["s"][s] = d
        print(f"s={s:5s} th2={th:.4f} | T_I^R={cen['TI']:.6f}({err['TI']:.6f}) T_KT^R={cen['TKT']:.6f}({err['TKT']:.6f}) "
              f"T_eff^R={Te:.6f} [{d['which']}] | T_ups^full={c:.6f}+-{hw:.6f} | R={R:.3f}+-{dR:.3f} (unreduced "
              f"{d['ratio_unreduced']:.3f}) pass={d['pass_5pct']} | sameL16 {sameL['16']['ratio']:.3f} sameL24 "
              f"{sameL['24']['ratio']:.3f} | K0-sens R={R / scale:.3f}")
    S["decision"] = dict(all_within_5pct=all(S["s"][s]["pass_5pct"] for s in SS),
                         ratios=[S["s"][s]["ratio"] for s in SS])
    print("decision:", S["decision"])
    OUT.write_text(json.dumps(S, indent=1, default=float))


if __name__ == "__main__":
    main()
