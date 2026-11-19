"""N3f (docs/N3f_analysis_plan.md): exact Debye-Waller bond averages <d^l_mm(beta_ab)> of the relative C3 tilt.

    python scripts/n3f_analysis.py theta     # theta^2 at T_Ups^full -> data/n3f/debye_waller.json, eff pilot configs
    python scripts/n3f_analysis.py configs   # effective-model production configs from the L = 16 pilots
    python scripts/n3f_analysis.py decide    # -> data/n3f/n3f_summary.json
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts"))
import n1_analysis as n1  # noqa: E402
from n2_analysis import jackknife  # noqa: E402
from n3d_analysis import quantities  # noqa: E402
from n3s_analysis import last_down_crossing  # noqa: E402
from sbaniso.analysis import curves, equilibration_check, load_run, per_temperature  # noqa: E402

SS = ("1.0", "0.03")
P = json.loads((REPO / "data" / "n3s" / "params.json").read_text())
N3S = json.loads((REPO / "data" / "n3s" / "n3s_summary.json").read_text())["s"]
N3D = json.loads((REPO / "data" / "n3d" / "n3d_summary.json").read_text())["s"]
SC = json.loads((REPO / "data" / "n3s" / "selfconsistent.json").read_text())
OUT = REPO / "data" / "n3f"
N3E = json.loads((REPO / "data" / "n3e" / "n3e_summary.json").read_text())["s"]


def full_point(s, start):
    f = REPO / f"runs/n3f/full_s{s}_L24_{start}.h5"
    run = load_run(f); pt = per_temperature(run); T = run.Ts
    with h5py.File(f) as h:
        x = h["timeseries/extra"][()]
        hist = h["checkpoint_beta_hist"][()]
    th2 = np.nanmean(x[..., 2], axis=1); K = x[..., 0].mean(1)
    fac = np.nanmean(x[..., 3:6], axis=1)                # (R, 3): <d3_00>, <d3_33>, <d4_33>
    tt = np.nanmean(x[..., 6], axis=1)
    tu = last_down_crossing(T, pt["Ups"][0] * np.pi / (2 * T) - 1)
    k = int(np.argmin(abs(T - tu)))
    return dict(T_ups=tu, Ts=T.tolist(), theta2=th2.tolist(), K=K.tolist(),
                theta2_at_Tups=float(np.interp(tu, T, th2)), K_at_Tups=float(np.interp(tu, T, K)),
                factors_at_Tups=[float(np.interp(tu, T, fac[:, q])) for q in range(3)],
                Cnn_at_Tups=float(np.interp(tu, T, tt) / np.interp(tu, T, th2)),
                beta_hist_nearest_T=dict(T=float(T[k]), counts=hist[k].tolist()),
                eq={k: equilibration_check(run)[k] for k in ("pass", "tau_max", "roundtrips")})


def theta():
    OUT.mkdir(exist_ok=True)
    D = {}
    for s in SS:
        r, o = full_point(s, "random"), full_point(s, "ordered")
        th2 = 0.5 * (r["theta2_at_Tups"] + o["theta2_at_Tups"])
        fac = [0.5 * (a + b) for a, b in zip(r["factors_at_Tups"], o["factors_at_Tups"])]
        cnn = 0.5 * (r["Cnn_at_Tups"] + o["Cnn_at_Tups"])
        p = P[s]
        JR = [p["J1"] * fac[0], p["J2"] * fac[1], p["J4"] * fac[2]]
        old = N3S[s]["full"]["24"]
        D[s] = dict(factors=fac, C_nn=cnn, factors_per_body_N3e=[1 - 6 * th2, 1 - 1.5 * th2, 1 - 5.5 * th2],
                    theta2=th2, theta2_random=r["theta2_at_Tups"], theta2_ordered=o["theta2_at_Tups"],
                    theta2_Kbased=N3D[s]["theta2"], JR=JR, T_ups_random=r["T_ups"], T_ups_ordered=o["T_ups"],
                    T_ups_previous=[old["random"]["T_ups"], old["ordered"]["T_ups"]],
                    full=dict(random=r, ordered=o))
        print(f"s={s}: factors {np.round(fac, 4).tolist()} (per-body {np.round(D[s]['factors_per_body_N3e'], 4).tolist()}) "
              f"C_nn {cnn:.3f} | theta2 rand/ord {r['theta2_at_Tups']:.4f}/{o['theta2_at_Tups']:.4f} -> {th2:.4f} "
              f"(K-based {N3D[s]['theta2']:.4f}) | T_ups rand/ord {r['T_ups']:.6f}/{o['T_ups']:.6f} "
              f"(previous {old['random']['T_ups']:.6f}/{old['ordered']['T_ups']:.6f}) | JR={np.round(JR, 7).tolist()}")
        Tu = 0.5 * (r["T_ups"] + o["T_ups"])
        (REPO / "configs/n3f/pilot").mkdir(parents=True, exist_ok=True)
        c = dict(name=f"n3f_pilot_eff_s{s}_L16", model=dict(couplings=JR), L=16,
                 temperatures=dict(Tmin=round(0.5 * Tu, 8), Tmax=round(1.4 * Tu, 8), n=24), ntherm=2000,
                 nmeas=8000, clusters=[1, 1, 1], seed=95000 + SS.index(s), start="random",
                 checkpoint_every=5000, output=f"runs/n3f/pilot_eff_s{s}_L16.h5")
        (REPO / f"configs/n3f/pilot/eff_s{s}_L16.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
    (OUT / "debye_waller.json").write_text(json.dumps(D, indent=1, default=float))


def configs():
    D = json.loads((OUT / "debye_waller.json").read_text())
    (REPO / "configs/n3f/eff").mkdir(parents=True, exist_ok=True)
    order = []
    for i, s in enumerate(SS):
        run = load_run(REPO / f"runs/n3f/pilot_eff_s{s}_L16.h5")
        T = np.linspace(run.Ts.min(), run.Ts.max(), 400); v, _, _, _ = curves(run, T)
        mid = 0.5 * (n1.crossing(T, v["Ups"] - 2 * T / np.pi) + n1.crossing(T, v["U4"] - 0.45))
        D[s]["pilot_mid"] = float(mid)
        for j, (L, n) in enumerate(((16, 16), (24, 20), (32, 24), (64, 40))):
            name = f"eff_s{s}_L{L}"
            c = dict(name=f"n3f_{name}", model=dict(couplings=D[s]["JR"]), L=L,
                     temperatures=dict(Tmin=round(float(0.8 * mid), 8), Tmax=round(float(1.2 * mid), 8), n=n), ntherm=2000,
                     nmeas=40000, clusters=[1, 1, 1], seed=96000 + 10 * i + j, start="random",
                     checkpoint_every=2000, output=f"runs/n3f/{name}.h5")
            (REPO / f"configs/n3f/eff/{name}.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
            order.append((-L, name))
        print(f"s={s}: pilot mid {mid:.6f}")
    (REPO / "configs/n3f/eff/ORDER").write_text("".join(f"configs/n3f/eff/{n}.yaml\n" for _, n in sorted(order)))
    (OUT / "debye_waller.json").write_text(json.dumps(D, indent=1, default=float))


def decide():
    D = json.loads((OUT / "debye_waller.json").read_text())
    (OUT / "curves").mkdir(exist_ok=True)
    S = dict(plan="docs/N3f_analysis_plan.md", s={})
    for s in SS:
        cs = []
        for L in (16, 24, 32, 64):
            tgt = OUT / "curves" / f"eff_s{s}_L{L}.npz"
            if not tgt.exists():
                run = load_run(REPO / f"runs/n3f/eff_s{s}_L{L}.h5")
                T = np.linspace(run.Ts.min(), run.Ts.max(), 361)
                v, e, p, _ = curves(run, T)
                np.savez_compressed(tgt, T=T, L=run.L, Ts=run.Ts, **{f"v_{q}": v[q] for q in n1.KEYS},
                                    **{f"e_{q}": e[q] for q in n1.KEYS}, **{f"p_{q}": p[q] for q in n1.KEYS},
                                    hist=np.zeros(1), T_Cmax=np.nan)
            cs.append(n1.Curves(tgt))
        T = cs[0].T
        cov, _ = n1.crossing_cov(T, cs)
        cen, err = jackknife(lambda x: quantities(T, x, cov), cs)
        Te, dTe = float(cen["Teff"]), float(err["Teff"])
        a, b = D[s]["T_ups_random"], D[s]["T_ups_ordered"]
        c = 0.5 * (a + b); hw = 0.5 * abs(a - b)
        R = c / Te; dR = float(np.hypot(hw / Te, c * dTe / Te ** 2))
        te24 = float(cen["ups_cross_L"][1])
        eqs = [equilibration_check(load_run(REPO / f"runs/n3f/eff_s{s}_L{L}.h5"))["pass"] for L in (16, 24, 32, 64)]
        S["s"][s] = dict(theta2=D[s]["theta2"], theta2_Kbased=D[s]["theta2_Kbased"], JR=D[s]["JR"],
                         TI=float(cen["TI"]), TKT=float(cen["TKT"]), T_eff_R=Te, dT_eff_R=dTe,
                         which="T_I" if cen["TI"] >= cen["TKT"] else "T_KT", T_ups_full=c, halfwidth=hw,
                         ratio=R, ratio_err=dR, pass_5pct=bool(abs(R - 1) <= 0.05),
                         sameL24=dict(eff=te24, ratio=c / te24, ratio_random=a / te24, ratio_ordered=b / te24),
                         ratio_N3d=N3D[s]["ratio"], ratio_N3e=N3E[s]["ratio"], factors=D[s]["factors"], C_nn=D[s]["C_nn"], eff_equilibrated=eqs, selfconsistent=SC[s])
        print(f"s={s}: theta2 {D[s]['theta2']:.4f} (K-based {D[s]['theta2_Kbased']:.4f}) | T_I^R={cen['TI']:.6f} "
              f"T_KT^R={cen['TKT']:.6f} T_eff^R={Te:.6f} | T_ups^full={c:.6f}+-{hw:.6f} | R={R:.3f}+-{dR:.3f} "
              f"(N3e {N3E[s]['ratio']:.3f}, N3d {N3D[s]['ratio']:.3f}) pass={abs(R - 1) <= 0.05} | sameL24 {c / te24:.3f} | eff eq {eqs}")
    S["decision"] = dict(both_within_5pct=all(S["s"][s]["pass_5pct"] for s in SS))
    print("decision:", S["decision"])
    (OUT / "n3f_summary.json").write_text(json.dumps(S, indent=1, default=float))


if __name__ == "__main__":
    {"theta": theta, "configs": configs, "decide": decide}[sys.argv[1]]()
