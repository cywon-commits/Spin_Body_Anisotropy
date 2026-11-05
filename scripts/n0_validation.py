"""N0 validation analysis: reads runs/n0/*.h5, writes data/n0/validation_summary.json.

    python scripts/n0_validation.py

Checks (docs/TASK.md, N0; thresholds fixed here before the runs were analysed)
  XY      beta*Ups(L) at beta_KT = 1.1199 vs Hasenbusch (2005) Table 1 (L = 16, 32, 64);
          Ups = 2T/pi crossings T_L; Weber-Minnhagen T_KT over all L: |T_KT - 0.8929| <= 0.01
  Ising   Binder crossings (16,32), (32,64): |T - 2.269| <= 0.01, |U* - 0.611| <= 0.01
  r = 0, J2 = 0   the same two checks in units of the relevant coupling, all cluster moves on
  pilot   L = 32, r = 3.8, 4.1, 4.4: U4, Ups, C, |m|, E per temperature vs data/pilot_v3 (z-scores)
          r = 4.1 random vs ordered start
  every run  quarter-split z, PT round trips, swap acceptance, tau_int
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from sbaniso.analysis import (binder_crossing, curves, load_run, per_temperature,  # noqa: E402
                              quarter_check, tau_table, ups_crossing, weber_minnhagen)

RUNS = REPO / "runs" / "n0"
OUT = REPO / "data" / "n0" / "validation_summary.json"
J1 = 160 / 81; J2 = 128 / 81; J4 = 4 * 128 / 81
BETA_KT = 1.1199; T_KT = 1 / BETA_KT; T_ISING = 2 / np.log(1 + np.sqrt(2)); U_STAR = 0.61069
# Hasenbusch, J. Phys. A 38, 5869 (2005), Table 1: beta*Upsilon at beta = 1.1199
HASENBUSCH = {16: (0.72536, 0.00007), 32: (0.70883, 0.00007), 64: (0.69785, 0.00007),
              128: (0.69001, 0.00007)}


def f(x):
    return None if x is None or not np.isfinite(x) else float(x)


def diagnostics(name, run):
    zq, qd = quarter_check(run)
    tau = tau_table(run)
    d = run.diagnostics
    return dict(name=name, L=run.L, R=len(run.Ts), nmeas=run.obs.shape[1],
                quarter_z_max=f(zq.max()), quarter_threshold=qd["threshold"], roundtrips=int(d["roundtrips"]),
                swap_acc_min=f(d["swap_acc"].min()), swap_acc_max=f(d["swap_acc"].max()),
                tau_max={k: f(v.max()) for k, v in tau.items()},
                seconds=f(d["seconds"]), git=str(run.attrs.get("git")),
                src_sha256=str(run.attrs.get("src_sha256")))


def kt_family(prefix, Ls, J):
    """KT checks for runs '<prefix>_L*' with XY coupling J (temperatures in units of J)."""
    runs = {L: load_run(RUNS / f"{prefix}_L{L}.h5") for L in Ls}
    Tlo = max(r.Ts.min() for r in runs.values()); Thi = min(r.Ts.max() for r in runs.values())
    T = np.linspace(Tlo, Thi, 201)
    curs = {L: curves(runs[L], T) for L in Ls}
    out = dict(J=J, sizes=list(Ls), hasenbusch=[], crossings=[], diagnostics=[])
    t = np.array([T_KT * J])
    for L in Ls:
        v, e, _, _ = curves(runs[L], t, mh=curs[L][3])
        bu = float(v["Ups"][0] / t[0]); dbu = float(e["Ups"][0] / t[0])
        row = dict(L=L, betaUps=bu, err=dbu)
        if L in HASENBUSCH:
            ref, dref = HASENBUSCH[L]
            row.update(ref=ref, z=(bu - ref) / np.hypot(dbu, dref))
        out["hasenbusch"].append(row)
        c = ups_crossing(T, curs[L])
        out["crossings"].append(dict(L=L, T=f(c["T"] / J), dT=f(c["dT"] / J)))
        out["diagnostics"].append(diagnostics(f"{prefix}_L{L}", runs[L]))
    wm = weber_minnhagen(T, [curs[L] for L in Ls], Ls)
    out["weber_minnhagen"] = dict(T=wm["T"] / J, dT=wm["dT"] / J, chi2dof=wm["chi2dof"], c=wm["c"])
    big = [L for L in Ls if L >= 24]
    if len(big) >= 2 and big != list(Ls):
        wm2 = weber_minnhagen(T, [curs[L] for L in big], big)
        out["weber_minnhagen_Lmin24"] = dict(T=wm2["T"] / J, dT=wm2["dT"] / J, chi2dof=wm2["chi2dof"],
                                             c=wm2["c"])
    wmT = out["weber_minnhagen"]["T"]; wmdT = out["weber_minnhagen"]["dT"]
    out["pass"] = dict(
        weber_minnhagen=bool(abs(wmT - T_KT) <= 0.01),
        hasenbusch=bool(all(abs(r["z"]) < 4 for r in out["hasenbusch"] if "z" in r)),
        crossings_decrease=bool(np.all(np.diff([c["T"] for c in out["crossings"]]) < 0)))
    out["T_KT_ref"] = T_KT
    return out


def ising_family(prefix, Ls, J):
    runs = {L: load_run(RUNS / f"{prefix}_L{L}.h5") for L in Ls}
    Tlo = max(r.Ts.min() for r in runs.values()); Thi = min(r.Ts.max() for r in runs.values())
    T = np.linspace(Tlo, Thi, 401)
    curs = {L: curves(runs[L], T) for L in Ls}
    out = dict(J=J, sizes=list(Ls), binder=[], diagnostics=[])
    for a, b in zip(Ls, Ls[1:]):
        c = binder_crossing(T, curs[a], curs[b])
        out["binder"].append(dict(pair=[a, b], T=f(c["T"] / J), dT=f(c["dT"] / J), U=f(c["U"]),
                                  dU=f(c["dU"])))
    for L in Ls:
        out["diagnostics"].append(diagnostics(f"{prefix}_L{L}", runs[L]))
    last = out["binder"][-1]
    out["pass"] = dict(T=bool(abs(last["T"] - T_ISING) <= 0.01), U=bool(abs(last["U"] - 0.611) <= 0.01))
    out["T_ref"] = T_ISING; out["U_ref"] = U_STAR
    return out


def compare(a, b, names=("E", "U4", "Ups", "C", "absm")):
    """Per-temperature z-scores between two runs at the same temperatures."""
    pa, pb = per_temperature(a), per_temperature(b)
    assert np.allclose(a.Ts, b.Ts)
    res = {}
    for k in names:
        z = (pa[k][0] - pb[k][0]) / np.hypot(pa[k][1], pb[k][1])
        res[k] = dict(z=[f(x) for x in z], max_abs_z=f(np.abs(z).max()), chi2dof=f(np.mean(z ** 2)))
    allz = np.concatenate([np.asarray(res[k]["z"], float) for k in names])
    return res, dict(n=len(allz), frac_gt3=f(np.mean(np.abs(allz) > 3)), chi2dof=f(np.mean(allz ** 2)),
                     max_abs_z=f(np.abs(allz).max()))


def main():
    S = {}
    S["xy"] = kt_family("xy", (16, 24, 32, 48, 64), 1.0)
    S["ising"] = ising_family("ising", (16, 32, 64), 1.0)
    S["r0_kt"] = kt_family("r0_kt", (16, 32, 64), J2)
    S["r0_ising"] = ising_family("r0_ising", (16, 32, 64), J1)
    S["j2zero_kt"] = kt_family("j2zero_kt", (16, 32, 64), J4)
    S["j2zero_ising"] = ising_family("j2zero_ising", (16, 32, 64), J1)
    reg = {}
    for r in (3.8, 4.1, 4.4):
        new = load_run(RUNS / f"pilot_r{r}_L32.h5")
        old = load_run(REPO / "data" / "pilot_v3" / f"pt_r{r}_L32.npz")
        per, tot = compare(new, old)
        reg[str(r)] = dict(per_observable=per, total=tot, diagnostics=diagnostics(f"pilot_r{r}_L32", new))
    a = load_run(RUNS / "pilot_r4.1_L32.h5"); b = load_run(RUNS / "pilot_r4.1_L32_ordered.h5")
    per, tot = compare(a, b)
    reg["4.1_random_vs_ordered"] = dict(per_observable=per, total=tot,
                                        diagnostics=diagnostics("pilot_r4.1_L32_ordered", b))
    S["regression"] = reg
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(S, indent=1, default=float))
    # console summary
    for k in ("xy", "r0_kt", "j2zero_kt"):
        w = S[k]["weber_minnhagen"]
        print(f"{k:12s} WM T_KT/J = {w['T']:.4f} +- {w['dT']:.4f} (chi2/dof {w['chi2dof']:.2f})  "
              f"betaUps z = {[round(r.get('z', np.nan), 2) for r in S[k]['hasenbusch']]}  "
              f"T_L = {[round(c['T'], 4) for c in S[k]['crossings']]}  pass {S[k]['pass']}")
    for k in ("ising", "r0_ising", "j2zero_ising"):
        print(f"{k:12s} Binder " + "; ".join(f"{b['pair']}: T/J = {b['T']:.4f}+-{b['dT']:.4f} U* = {b['U']:.4f}+-{b['dU']:.4f}"
                                            for b in S[k]["binder"]) + f"  pass {S[k]['pass']}")
    for k, v in reg.items():
        print(f"regression {k}: {v['total']}")
    diags = [d for k in S if k != "regression" for d in S[k]["diagnostics"]] + \
            [v["diagnostics"] for v in reg.values()]
    print("worst quarter z %.2f, min round trips %d, min swap acc %.3f, max tau(E) %.1f, max tau(|m|) %.1f" % (
        max(d["quarter_z_max"] for d in diags), min(d["roundtrips"] for d in diags),
        min(d["swap_acc_min"] for d in diags), max(d["tau_max"]["E"] for d in diags),
        max(d["tau_max"]["absm"] for d in diags)))


if __name__ == "__main__":
    main()
