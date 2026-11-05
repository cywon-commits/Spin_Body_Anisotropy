"""Figure (c): specific-heat maximum vs L for every r, from data/n1/n1_summary.json.

    python scripts/fig_n1_cmax.py   -> figures/n1_cmax.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
S = json.loads((REPO / "data" / "n1" / "n1_summary.json").read_text())
fig, ax = plt.subplots(figsize=(3.4, 2.8))
cm = plt.get_cmap("viridis")
rs = sorted(S["per_r"], key=float)
for i, k in enumerate(rs):
    d = S["per_r"][k]; L = np.array(d["L"], float)
    ax.errorbar(L, d["Cmax"], d["dCmax"], fmt="o-", ms=3, lw=0.8, capsize=1.5, c=cm(i / (len(rs) - 1)),
                label=f"r = {k}")
ax.set_xscale("log"); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([32, 48, 64, 96, 128]); ax.set_xticklabels(["32", "48", "64", "96", "128"])
ax.set_xlabel(r"$L$"); ax.set_ylabel(r"$C_{\max}$ per site")
ax.legend(fontsize=5.5, ncol=2, frameon=False)
plt.tight_layout()
(REPO / "figures").mkdir(exist_ok=True)
plt.savefig(REPO / "figures" / "n1_cmax.png", dpi=200)
