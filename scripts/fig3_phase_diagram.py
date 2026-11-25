"""Fig. 3: the Ising and KT lines of the C3 effective model cross at a decoupled tetracritical point
r_x = 4.340(1)(24), T_x = 6.995(1)(7) G, beyond which the order of the transitions is inverted.
(a) T-r phase diagram: N1 (r = 3.6-5.0) and N2 (r = 5.0, 6.0) T_I (Binder crossings, L -> inf) and
    T_KT (Weber-Minnhagen), joined by thin lines (phase boundaries), the decoupling approximation of theory_summary P6-P8 (lines), the four phases.
(b) T_I - T_KT vs r near the crossing; open symbols use T_KT x 1.0025 (Weber-Minnhagen bias bound).
Analytic lines: KT, T = 0.893 (J4 + J2 <ss'>) with the exact Onsager nearest-neighbour correlation of an
Ising model whose T_c is T_x^an = 2.269 (J1 + 0.723 J2); Ising, T = T_x^an (psi correlation held at its
critical value c_psi = 0.723; valid near the crossing only).  Reads data/n1/n1_summary.json, data/n2/n2_summary.json.

    python scripts/fig3_phase_diagram.py   -> figures/fig3_phase_diagram.{pdf,png}
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import ellipk

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from figstyle import OI, WIDTH, panel, plt, save  # noqa: E402

J1, J2 = 160 / 81, 128 / 81
KC = 0.5 * np.log(1 + np.sqrt(2))
TX_AN = 2.269 * (J1 + 0.723 * J2)


def onsager_nn(K):
    """Exact nearest-neighbour correlation <s s'> of the square-lattice Ising model at coupling K."""
    k = 2 * np.sinh(2 * K) / np.cosh(2 * K) ** 2
    return 0.5 / np.tanh(2 * K) * (1 + 2 / np.pi * (2 * np.tanh(2 * K) ** 2 - 1) * ellipk(k ** 2))


def kt_line(r):
    f = lambda T: T - 0.893 * J2 * (r + onsager_nn(KC * TX_AN / T + 1e-12))
    return brentq(f, 3, 15)


N1 = json.loads((REPO / "data" / "n1" / "n1_summary.json").read_text())
N2 = json.loads((REPO / "data" / "n2" / "n2_summary.json").read_text())["per_r"]
P = dict(N1["per_r"]); eps = N1["eps_wm"]
for r in ("6.0",):
    P[r] = N2[r]
rs = sorted(P, key=float); r = np.array([float(x) for x in rs])
g = lambda k: np.array([P[x][k] for x in rs])
TI, dTI, TK, dTK = g("TI"), g("dTI"), g("TKT"), g("dTKT")
n2 = np.array([x in N2 for x in rs])
cr = N1["crossing"]["main"]; syst = N1["crossing"]["syst_r"]

fig, (a, b) = plt.subplots(2, 1, figsize=(WIDTH, 5.0), gridspec_kw=dict(height_ratios=[1.35, 1]))
rr = np.linspace(3.5, 6.15, 200)
a.plot(rr, [kt_line(x) for x in rr], c=OI["blue"], lw=0.8, ls="--", label="KT, decoupling approx.")
a.plot(rr, np.full_like(rr, TX_AN), c=OI["red"], lw=0.8, ls=":",
       label=r"Ising, decoupling approx. ($c_\psi$ = 0.723)")
a.plot(r, TI, "-", c=OI["red"], lw=0.5); a.plot(r, TK, "-", c=OI["blue"], lw=0.5)
a.errorbar(r, TI, dTI, fmt="o", c=OI["red"], label=r"Ising $T_I$")
a.errorbar(r, TK, dTK, fmt="s", c=OI["blue"], label=r"KT $T_{KT}$")
a.plot(r[n2], TI[n2], "o", mfc="none", mec="k", ms=6, mew=0.6)
a.plot(r[n2], TK[n2], "s", mfc="none", mec="k", ms=6, mew=0.6, label="inverted-order runs")
r5 = N2["5.0"]; a.plot([5.0, 5.0], [r5["TI"], r5["TKT"]], c="k", lw=0.6)
a.errorbar([cr["r_x"]], [cr["T_x"]], xerr=[np.hypot(cr["dr_x"], syst)], yerr=[cr["dT_x"]], fmt="*", c="k", ms=7,
           label=r"crossing $(r_\times,T_\times)$")
a.text(4.85, 8.45, "disordered", ha="left", fontsize=6.5)
a.annotate("σ order", xy=(3.72, 6.6), xytext=(3.6, 7.3), fontsize=6, c=OI["red"],
           arrowprops=dict(arrowstyle="-|>", lw=0.5, color=OI["red"], shrinkA=0, shrinkB=0, mutation_scale=5))
a.text(5.6, 7.85, "ψ QLRO", fontsize=6, ha="center", va="center", c=OI["blue"])
a.text(5.0, 6.35, "σ order + ψ QLRO\n(octupolar QLRO)", fontsize=6, ha="center", va="center")
a.set_xlim(3.5, 6.15); a.set_ylim(5.9, 9.2)
a.set_xlabel(r"$r=J_4/J_2$"); a.set_ylabel(r"$T/G$")
a.legend(loc="upper left", fontsize=5.8, handlelength=1.6)

m = r <= 5.0
d = TI[m] - TK[m]; dd = np.hypot(dTI[m], dTK[m])
b.axhline(0, c="k", lw=0.5)
b.errorbar(r[m], d, dd, fmt="o", c=OI["black"], label=r"$T_I-T_{KT}$")
b.plot(r[m], TI[m] - TK[m] * (1 + eps), "o", mfc="none", c=OI["grey"], label=r"$T_I-1.0025\,T_{KT}$")
rl = np.linspace(3.55, 5.05, 100)
b.plot(rl, TX_AN - np.array([kt_line(x) for x in rl]), "--", c=OI["purple"], lw=0.8, label="decoupling approx.")
b.axvspan(cr["r_x"] - syst, cr["r_x"] + syst, color="0.85", lw=0, zorder=0)
b.axvline(cr["r_x"], c="k", lw=0.5, ls=":")
b.set_xlabel(r"$r=J_4/J_2$"); b.set_ylabel(r"$(T_I-T_{KT})/G$")
b.legend(loc="lower left", fontsize=6)
panel(a, "a", x=-0.2, y=1.03); panel(b, "b", x=-0.2, y=1.03)
fig.tight_layout(h_pad=0.5)
save(fig, "fig3_phase_diagram", REPO)
