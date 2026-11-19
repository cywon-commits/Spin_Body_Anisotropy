"""Fig. 2: uniaxial anisotropy selects a body axis through a single cubic coefficient Gamma,
E_eff(n) = -Gamma sum_alpha n_alpha^4, Gamma = (35/9) a4 + 2 a2^2 (3/lam_E - 1/lam_T2)  (theory_summary P2, P3).
(a) Fixed-orientation relaxed single-body energy vs K = sum n_alpha^4 (43 orientations, a4 = -0.02,
    K = 2, zG = 0.08, zH4 = 0.80): a straight line of slope -Gamma (formula) ending at the S4 (K = 1)
    and C3 (K = 1/3) axes.
(b) Selection by a2 alone: S4 along z when 27 kappa1 - 63 kappa3 + 35 kappa4 = 27J + 9zG - 8zH4 > 0
    (kappa2 cancels, no K dependence); the line lies above the T-Sq boundary zH4 = 0.692K + 0.923zG for
    K <= 2, so the C3 phase needs a4 < 0.  Diamonds: the N3 points.
(c) Relative deviation of the numerical Gamma from the formula vs anisotropy amplitude: linear,
    i.e. next order (a4^2, a2^3); the +-a2 average is exact to 0.06 %.
Reads data/n4/gamma_check.json (scripts/n4_gamma_check.py).

    python scripts/fig2_axis_selection.py   -> figures/fig2_axis_selection.{pdf,png}
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from figstyle import OI, WIDTH, panel, plt, save  # noqa: E402

D = json.loads((REPO / "data" / "n4" / "gamma_check.json").read_text())
fig, (a, b, c) = plt.subplots(3, 1, figsize=(WIDTH, 6.6))

K = np.array([o["K"] for o in D["orientations"]]); E = np.array([o["E"] for o in D["orientations"]])
G = D["gamma_formula_a4"]
c0 = np.mean(E + G * K)
a.plot(K, (E - c0) * 1e3, "o", c=OI["blue"], mfc="none", label="relaxed single body")
kk = np.linspace(1 / 3, 1, 10)
a.plot(kk, -G * kk * 1e3, "-", c="k", lw=0.8, label=rf"$-\Gamma K$, $\Gamma$ = {G:.4f} (formula)")
for k0, name, xy in ((1 / 3, r"$C_3\parallel\hat z$ (selected)", (0.42, 26.5)), (1.0, r"$S_4\parallel\hat z$", (0.8, 76))):
    a.annotate(name, xy=(k0, -G * k0 * 1e3), xytext=xy, fontsize=6.5, va="center",
               arrowprops=dict(arrowstyle="-", lw=0.5, shrinkA=1, shrinkB=3))
a.set_xlabel(r"$K=\sum_\alpha n_\alpha^4$"); a.set_ylabel(r"$E-E_0$ ($10^{-3}J$)")
a.legend(loc="upper left")

zg = np.linspace(0, 1.2, 50)
b.fill_between(zg, (27 + 9 * zg) / 8, 4.2, color=OI["sky"], alpha=0.35, lw=0)
b.fill_between(zg, 0, (27 + 9 * zg) / 8, color=OI["orange"], alpha=0.25, lw=0)
b.plot(zg, (27 + 9 * zg) / 8, c="k", lw=0.9)
b.text(0.05, 3.85, r"$a_4=0$: $C_3\parallel\hat z$", fontsize=6.5)
b.text(0.05, 2.75, r"$a_4=0$: $S_4\parallel\hat z$", fontsize=6.5)
b.text(1.18, 3.42, r"$27\kappa_1-63\kappa_3+35\kappa_4=0$", fontsize=6, ha="right", va="bottom")
for Kk, ls in ((0.5, ":"), (2.0, "--")):
    b.plot(zg, 0.692 * Kk + 0.923 * zg, c=OI["red"], lw=0.8, ls=ls, label=f"T–Sq boundary, K = {Kk:g}")
b.plot([0.08] * 3, [0.72, 0.80, 0.88], "D", c=OI["blue"], ms=3, label="N3 points (K = 2)")
b.set_xlim(0, 1.2); b.set_ylim(0, 4.2)
b.set_xlabel(r"$zG/J$"); b.set_ylabel(r"$zH_4/J$")
b.legend(loc="center right", fontsize=5.8)

T = D["table"]
amp4 = [-t["a4"] for t in T if t["a2"] == 0]; dev4 = [t["numeric"] / t["formula"] - 1 for t in T if t["a2"] == 0]
c.plot(amp4, np.array(dev4) * 100, "o-", c=OI["blue"], label=r"$a_4<0$ (amplitude $|a_4|$)")
for t, mk in zip([t for t in T if t["a2"] != 0], ("^", "v")):
    c.plot([abs(t["a2"])], [(t["numeric"] / t["formula"] - 1) * 100], mk, c=OI["red"],
           label=rf"$a_2={t['a2']:+.2f}$ (amplitude $|a_2|$)")
p2 = [t for t in T if t["a2"] != 0]
avg = 0.5 * (p2[0]["numeric"] + p2[1]["numeric"]) / p2[0]["formula"] - 1
c.plot([0.05], [avg * 100], "*", c="k", ms=6, label=r"$\pm a_2$ average")
c.axhline(0, c="k", lw=0.5)
c.set_xlim(0, 0.06); c.set_xlabel("anisotropy amplitude"); c.set_ylabel(r"$\Gamma_{\rm num}/\Gamma_{\rm formula}-1$ (%)")
c.legend(loc="upper left", fontsize=5.8)
for ax, l in zip((a, b, c), "abc"):
    panel(ax, l, x=-0.2, y=1.03)
fig.tight_layout(h_pad=0.5)
save(fig, "fig2_axis_selection", REPO)
