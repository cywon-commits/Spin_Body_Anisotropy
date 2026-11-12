"""Figure (d): full spin-body model vs the C3 effective model at the same (relaxed) couplings.
Top: Ups pi/(2T) (= rho_alpha pi/(18 T) for the full model); bottom: axis-selection indicator <K>
(full model; 1/3 = C3 axis locked along z, 0.6 = isotropic).  Columns: zH4 = 0.72, 0.80, 0.88.

    python scripts/fig_n3_full_vs_eff.py   -> figures/n3_full_vs_eff.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
S = json.loads((REPO / "data" / "n3" / "n3_summary.json").read_text())
CUR = REPO / "data" / "n3" / "curves"
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.4), sharex=True, sharey="row")
cm = plt.get_cmap("viridis")
for j, zh in enumerate(("0.72", "0.80", "0.88")):
    ax, bx = axs[0, j], axs[1, j]
    for i, L in enumerate((16, 24, 32)):
        c = np.load(CUR / f"full_zH4_{zh}_L{L}.npz"); T = c["T"]
        ax.plot(T, c["v_Ups"] * np.pi / (2 * T), c=cm(i / 3), lw=1, label=f"full L={L}")
        bx.plot(c["Ts"], c["K_by_T"], "o-", c=cm(i / 3), ms=2, lw=0.8)
    for i, L in enumerate((16, 32, 64)):
        c = np.load(CUR / f"eff_zH4_{zh}_L{L}.npz"); T = c["T"]
        ax.plot(T, c["v_Ups"] * np.pi / (2 * T), c="tab:red", lw=0.8, ls=("--", "-.", ":")[i], label=f"eff L={L}")
    e = S["points"][zh]["eff"]
    for a in (ax, bx):
        a.axvline(e["TI"], c="tab:red", lw=0.6); a.axvline(e["TKT"], c="tab:blue", lw=0.6)
    ax.axhline(1, c="k", lw=0.5, ls="--"); bx.axhline(1 / 3, c="k", lw=0.5, ls="--"); bx.axhline(0.6, c="0.5", lw=0.5, ls=":")
    ax.set_title(f"zH4 = {zh}  (r_eff = {dict(zip(('0.72','0.80','0.88'),(3.68,4.05,4.42)))[zh]})", fontsize=8)
    ax.set_xlim(0.066, 0.16); ax.set_ylim(-0.1, 3.0); bx.set_xlabel(r"$T/J$")
axs[0, 0].set_ylabel(r"$\Upsilon_\psi\,\pi/(2T)$"); axs[1, 0].set_ylabel(r"$\langle K\rangle$ (full model)")
axs[0, 0].legend(fontsize=5.5, frameon=False, loc="upper right")
plt.tight_layout()
plt.savefig(REPO / "figures" / "n3_full_vs_eff.png", dpi=200)
