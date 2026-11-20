"""Fig. 4: on both sides of the crossing (r = 4.2 and 4.4) the Ising and KT transitions are separate
and each has its own scaling: (a) sigma Binder cumulant vs (T - T_I) L^{1/nu}, nu = 1;
(b) Weber-Minnhagen form at T = T_KT: pi Ups_psi/(2T) vs 1/(2 ln L + C) for each L, with the reference
line of slope 1, pi Ups/(2T) = 1 + 1/(2 ln L + C).  C is fitted at T_KT with the slope held at 1
(weighted least squares); T_KT itself is the N1 simultaneous WM fit, so this panel shows the quality
of that fit (chi^2/dof in the legend), not an independent test;
(c) C_max vs ln L with the fit C_max = a + b ln L of the N1 analysis.
Reads data/n1/curves/*.npz and data/n1/n1_summary.json.

    python scripts/fig4_crossing_collapse.py   -> figures/fig4_crossing_collapse.{pdf,png}
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from figstyle import OI, SIZE_MARKERS, WIDTH, panel, plt, save, size_colors  # noqa: E402

S = json.loads((REPO / "data" / "n1" / "n1_summary.json").read_text())["per_r"]
RS = ("4.2", "4.4")
STY = {"4.2": dict(ls="-"), "4.4": dict(ls="--")}

fig, axs = plt.subplots(3, 1, figsize=(WIDTH, 6.3))
a, b, c = axs
CW = {}
for r in RS:
    d = S[r]; Ls = d["L"]; cols = size_colors(len(Ls))
    for i, L in enumerate(Ls):
        z = np.load(REPO / "data" / "n1" / "curves" / f"r{r}_L{L}.npz")
        x = (z["T"] - d["TI"]) * L
        m = np.abs(x) < 35
        a.plot(x[m], z["v_U4"][m], c=cols[i], label=f"L = {L}" if r == RS[0] else None, **STY[r])
for r, col, mk in zip(RS, (OI["blue"], OI["red"]), ("o", "s")):
    d = S[r]; Ls = np.array(d["L"]); y, dy = [], []
    for L in Ls:
        z = np.load(REPO / "data" / "n1" / "curves" / f"r{r}_L{L}.npz")
        f = np.pi / (2 * z["T"])
        y.append(np.interp(d["TKT"], z["T"], z["v_Ups"] * f)); dy.append(np.interp(d["TKT"], z["T"], z["e_Ups"] * f))
    y, dy = np.array(y), np.array(dy)
    res = lambda C: np.sum(((y - 1 - 1 / (2 * np.log(Ls) + C)) / dy) ** 2)
    C = CW[r] = minimize_scalar(res, bounds=(-2 * np.log(Ls.min()) + 0.2, 20), method="bounded").x
    xv = 1 / (2 * np.log(Ls) + C)
    b.errorbar(xv, y, dy, fmt=mk, c=col, mfc="none" if r == "4.4" else col, ms=3.5,
               label=rf"r = {r}: $C$ = {C:.2f}, WM $\chi^2$/dof = {d['WM_chi2dof']:.1f}")
xl = np.linspace(0, 0.16, 10)
b.plot(xl, 1 + xl, c="k", lw=0.8, label=r"slope 1: $\pi\Upsilon/2T=1+1/(2\ln L+C)$")
for r, col, mk in zip(RS, (OI["blue"], OI["red"]), ("o", "s")):
    d = S[r]; lnL = np.log(d["L"])
    c.errorbar(lnL, d["Cmax"], d["dCmax"], fmt=mk, c=col, label=f"r = {r}")
    f = d["Cmax_fits"]["lnL"]; xx = np.linspace(lnL.min() - 0.1, lnL.max() + 0.1, 10)
    c.plot(xx, f["a"] + f["b"] * xx, c=col, lw=0.8, ls=STY[r]["ls"])
a.plot([], [], c="0.3", ls="-", label="r = 4.2"); a.plot([], [], c="0.3", ls="--", label="r = 4.4")
a.set_xlabel(r"$(T-T_I)\,L$"); a.set_ylabel(r"$U_4$ ($\sigma$)")
a.legend(ncol=2, loc="lower left", fontsize=6)
b.set_xlabel(r"$1/(2\ln L+C)$  (at $T=T_{KT}$)")
b.set_ylabel(r"$\pi\Upsilon_\psi/(2T_{KT})$")
b.legend(loc="upper left", fontsize=5.6)
c.set_xlabel(r"$\ln L$"); c.set_ylabel(r"$C_{\max}$ per site")
c.legend(loc="lower right")
for ax, l in zip(axs, "abc"):
    panel(ax, l, x=-0.2, y=1.02)
fig.tight_layout(h_pad=0.6)
save(fig, "fig4_crossing_collapse", REPO)
