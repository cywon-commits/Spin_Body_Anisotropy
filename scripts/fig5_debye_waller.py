"""Fig. 5: the effective model with Debye-Waller-reduced couplings reproduces the full spin-body
model within 5 % only for l_Gamma < 1; at l_Gamma = 4.5 per-body and relative-tilt prescriptions
bracket the full-model value.
(a) R = T_Ups^full / T_eff vs l_Gamma: bare couplings (N3s), per-body DW factor with the directly
    measured tilt (N3e), exact bond average <d^l_mm(beta_ab)> of the relative tilt (N3f); grey band +-5 %.
(b) Distribution of the relative tilt beta_ab between neighbouring C3 axes at T_Ups^full (random start),
    with the Gaussian of variance <beta^2> = 2<theta^2>(1 - C_nn).
(c) Coupling factors J^R/J for J1 (l = 3, m = 0), J2 (3, 3), J4 (4, 3): per-body vs relative tilt.
Reads data/n3s, data/n3e, data/n3f.

    python scripts/fig5_debye_waller.py   -> figures/fig5_debye_waller.{pdf,png}
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from figstyle import OI, WIDTH, panel, plt, save  # noqa: E402

N3S = json.loads((REPO / "data" / "n3s" / "n3s_summary.json").read_text())["s"]
N3E = json.loads((REPO / "data" / "n3e" / "n3e_summary.json").read_text())["s"]
N3F = json.loads((REPO / "data" / "n3f" / "n3f_summary.json").read_text())["s"]
DWF = json.loads((REPO / "data" / "n3f" / "debye_waller.json").read_text())
SS = ("1.0", "0.3", "0.1", "0.03"); S2 = ("1.0", "0.03")
lG = {s: N3S[s]["params"]["l_Gamma"] for s in SS}

fig, (a, b, c) = plt.subplots(3, 1, figsize=(WIDTH, 6.6), gridspec_kw=dict(height_ratios=[1.15, 1, 0.9]))
a.axhspan(0.95, 1.05, color="0.88", zorder=0, lw=0); a.axhline(1, c="k", lw=0.5)
x = np.array([lG[s] for s in SS])
a.errorbar(x, [N3S[s]["result"]["T_ups"]["ratio"] for s in SS], [N3S[s]["result"]["T_ups"]["ratio_err"] for s in SS],
           fmt="o--", c=OI["grey"], label="bare couplings")
x2 = np.array([lG[s] for s in S2])
a.errorbar(x2 * 0.97, [N3E[s]["ratio"] for s in S2], [N3E[s]["ratio_err"] for s in S2], fmt="s-", c=OI["red"],
           mfc="white", label=r"per-body DW, $1-\frac{1}{2}[l(l+1)-m^2]\langle\theta^2\rangle$")
a.errorbar(x2 * 1.03, [N3F[s]["ratio"] for s in S2], [N3F[s]["ratio_err"] for s in S2], fmt="D-", c=OI["blue"],
           label=r"relative tilt, $\langle d^{\,l}_{mm}(\beta_{ab})\rangle$")
a.set_xscale("log"); a.set_xticks(x); a.set_xticklabels([f"{v:.2f}" for v in x]); a.minorticks_off()
a.set_xlabel(r"$\ell_\Gamma=\sqrt{\rho/|\Gamma|}$"); a.set_ylabel(r"$R=T_\Upsilon^{\rm full}/T_{\rm eff}$")
a.set_ylim(0.55, 1.36); a.legend(loc="upper left", fontsize=6)

edges = np.linspace(0, np.pi, 91); mid = 0.5 * (edges[1:] + edges[:-1]); w = edges[1] - edges[0]
for s, col, mk in zip(S2, (OI["blue"], OI["orange"]), ("-", "--")):
    h = np.array(DWF[s]["full"]["random"]["beta_hist_nearest_T"]["counts"], float); h /= h.sum() * w
    b2 = 2 * DWF[s]["theta2"] * (1 - DWF[s]["C_nn"])
    g = mid / (b2 / 2) * np.exp(-mid ** 2 / b2)          # 2D Gaussian tilt: P(beta) ~ beta exp(-beta^2/<beta^2>)
    deg = np.degrees(mid)
    b.plot(deg, h, mk, c=col, label=rf"$\ell_\Gamma$ = {lG[s]:.2f}: MC")
    b.plot(deg, g, ":", c=col, lw=0.9)
b.plot([], [], ":", c="k", lw=0.9, label=r"Gaussian, $\langle\beta^2\rangle=2\langle\theta^2\rangle(1-C_{nn})$")
b.axvline(35.26, c="0.5", lw=0.5, ls="-.")
b.text(36.5, b.get_ylim()[1] * 0.55, "axis fold-back\n35.3°", fontsize=6, c="0.35")
b.set_xlim(0, 70); b.set_xlabel(r"relative tilt $\beta_{ab}$ (deg)"); b.set_ylabel(r"$P(\beta_{ab})$ (rad$^{-1}$)")
b.legend(loc="upper right", fontsize=5.8)

names = (r"$J_1$" "\n" r"$(3,0)$", r"$J_2$" "\n" r"$(3,3)$", r"$J_4$" "\n" r"$(4,3)$")
xx = np.arange(3); wd = 0.2
for k, (s, hatch) in enumerate(zip(S2, ("", "///"))):
    pb = DWF[s]["factors_per_body_N3e"]; rl = DWF[s]["factors"]
    c.bar(xx + (2 * k - 1.5) * wd, pb, wd, color=OI["red"], alpha=0.85, hatch=hatch, edgecolor="white", lw=0.3,
          label=f"per-body, $\\ell_\\Gamma$ = {lG[s]:.2f}")
    c.bar(xx + (2 * k - 0.5) * wd, rl, wd, color=OI["blue"], alpha=0.85, hatch=hatch, edgecolor="white", lw=0.3,
          label=f"relative, $\\ell_\\Gamma$ = {lG[s]:.2f}")
c.set_xticks(xx); c.set_xticklabels(names); c.set_ylim(0, 1.25); c.set_ylabel(r"$J^R/J$")
c.legend(ncol=2, loc="upper center", fontsize=5.6, columnspacing=0.8, handlelength=1.4)
for ax, l in zip((a, b, c), "abc"):
    panel(ax, l, x=-0.2, y=1.03)
fig.tight_layout(h_pad=0.5)
save(fig, "fig5_debye_waller", REPO)
