"""Figure: the inverted-order window at r = 5.0 and 6.0 (N2).  Top: sigma Binder U4(L, T);
bottom: Upsilon pi / (2T); shaded: T_I < T < T_KT.  Reads data/n2/curves and n2_summary.json.

    python scripts/fig_n2_window.py   -> figures/n2_window.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
S = json.loads((REPO / "data" / "n2" / "n2_summary.json").read_text())
fig, axs = plt.subplots(2, 2, figsize=(6.8, 4.6), sharex="col")
cm = plt.get_cmap("viridis")
for j, r in enumerate(("5.0", "6.0")):
    d = S["per_r"][r]
    for i, L in enumerate(d["L"]):
        c = np.load(REPO / "data" / "n2" / "curves" / f"r{r}_L{L}.npz")
        T = c["T"]; col = cm(i / (len(d["L"]) - 1))
        axs[0, j].plot(T, c["v_U4"], c=col, lw=1, label=f"L = {L}")
        axs[1, j].plot(T, c["v_Ups"] * np.pi / (2 * T), c=col, lw=1)
    for ax in axs[:, j]:
        ax.axvspan(d["TI"], d["TKT"], color="0.85", zorder=0)
        for t in d["Tq"]:
            ax.axvline(t, c="0.5", lw=0.5, ls=":")
    axs[1, j].axhline(1, c="k", lw=0.6, ls="--")
    axs[0, j].set_title(f"r = {r}:  $T_I$ = {d['TI']:.3f},  $T_{{KT}}$ = {d['TKT']:.3f}", fontsize=8)
    axs[1, j].set_xlabel(r"$T/G$")
axs[0, 0].set_ylabel(r"$U_4$ ($\sigma$)"); axs[1, 0].set_ylabel(r"$\Upsilon\,\pi/(2T)$")
axs[0, 0].legend(fontsize=6, frameon=False)
plt.tight_layout()
(REPO / "figures").mkdir(exist_ok=True)
plt.savefig(REPO / "figures" / "n2_window.png", dpi=200)
