"""Fig. 4: on both sides of the crossing (r = 4.2 and 4.4) the Ising and KT transitions are separate
and each has its own scaling: (a) sigma Binder cumulant vs (T - T_I) L^{1/nu}, nu = 1;
(b) psi helicity modulus divided by the Weber-Minnhagen log factor 1 + 1/(2 ln L + C) vs the KT
variable (T/T_KT - 1)(ln L + C/2)^2, with C = mean over L of 1/(pi Ups/(2T) - 1) - 2 ln L at T_KT;
(c) C_max vs ln L with the fit C_max = a + b ln L of the N1 analysis.
Reads data/n1/curves/*.npz and data/n1/n1_summary.json.

    python scripts/fig4_crossing_collapse.py   -> figures/fig4_crossing_collapse.{pdf,png}
"""
import json
import sys
from pathlib import Path

import numpy as np

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
    cl = []
    for L in Ls:
        z = np.load(REPO / "data" / "n1" / "curves" / f"r{r}_L{L}.npz")
        y0 = np.interp(d["TKT"], z["T"], z["v_Ups"] * np.pi / (2 * z["T"]))
        cl.append(1 / (y0 - 1) - 2 * np.log(L))
    CW[r] = C = float(np.mean(cl))
    for i, L in enumerate(Ls):
        z = np.load(REPO / "data" / "n1" / "curves" / f"r{r}_L{L}.npz")
        T = z["T"]
        lab = f"L = {L}" if r == RS[0] else None
        x = (T - d["TI"]) * L
        m = np.abs(x) < 35
        a.plot(x[m], z["v_U4"][m], c=cols[i], label=lab, **STY[r])
        lam = np.log(L) + C / 2
        xk = (T / d["TKT"] - 1) * lam ** 2
        y = z["v_Ups"] * np.pi / (2 * T) / (1 + 1 / (2 * lam))
        m = (xk > -0.8) & (xk < 0.8)
        b.plot(xk[m], y[m], c=cols[i], **STY[r])
for r, col, mk in zip(RS, (OI["blue"], OI["red"]), ("o", "s")):
    d = S[r]; lnL = np.log(d["L"])
    c.errorbar(lnL, d["Cmax"], d["dCmax"], fmt=mk, c=col, label=f"r = {r}")
    f = d["Cmax_fits"]["lnL"]; xx = np.linspace(lnL.min() - 0.1, lnL.max() + 0.1, 10)
    c.plot(xx, f["a"] + f["b"] * xx, c=col, lw=0.8, ls=STY[r]["ls"])
a.plot([], [], c="0.3", ls="-", label="r = 4.2"); a.plot([], [], c="0.3", ls="--", label="r = 4.4")
a.set_xlabel(r"$(T-T_I)\,L$"); a.set_ylabel(r"$U_4$ ($\sigma$)")
a.legend(ncol=2, loc="lower left", fontsize=6)
b.axhline(1, c="k", lw=0.5, ls=":")
b.set_xlabel(r"$(T/T_{KT}-1)\,\lambda_L^2$,  $\lambda_L=\ln L + C/2$")
b.set_ylabel(r"$\pi\Upsilon_\psi/[2T(1+1/2\lambda_L)]$")
b.text(0.97, 0.95, "\n".join(f"r = {r}: C = {CW[r]:.2f}" for r in RS), transform=b.transAxes, ha="right", va="top", fontsize=6)
c.set_xlabel(r"$\ln L$"); c.set_ylabel(r"$C_{\max}$ per site")
c.legend(loc="lower right")
for ax, l in zip(axs, "abc"):
    panel(ax, l, x=-0.2, y=1.02)
fig.tight_layout(h_pad=0.6)
save(fig, "fig4_crossing_collapse", REPO)
