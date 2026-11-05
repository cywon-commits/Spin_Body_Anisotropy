"""Figure (a): T-r phase diagram near the crossing from data/n1/n1_summary.json.

    python scripts/fig_n1_phase_diagram.py   -> figures/n1_phase_diagram.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
S = json.loads((REPO / "data" / "n1" / "n1_summary.json").read_text())
P = S["per_r"]; r = np.array(sorted(float(k) for k in P))
g = lambda k: np.array([P[f"{x}"][k] for x in r])
TI, dTI, TK, dTK = g("TI"), g("dTI"), g("TKT"), g("dTKT")
c = S["crossing"]["main"]; syst = S["crossing"]["syst_r"]

fig, ax = plt.subplots(figsize=(3.4, 2.8))
ax.errorbar(r, TI, dTI, fmt="o", c="tab:red", ms=3.5, capsize=2, label=r"Ising $T_I$ (Binder, $L\to\infty$)")
ax.errorbar(r, TK, dTK, fmt="s", c="tab:blue", ms=3.5, capsize=2, label=r"KT $T_{KT}$ (Weber–Minnhagen)")
ax.plot(r, TK * (1 + S["eps_wm"]), "s", c="tab:blue", mfc="none", ms=3.5, alpha=0.5,
        label=r"$T_{KT}\times1.0025$ (fit-bias bound)")
ax.errorbar([c["r_x"]], [c["T_x"]], xerr=[np.hypot(c["dr_x"], syst)], yerr=[c["dT_x"]], fmt="D", c="k",
            ms=4, capsize=2, label=r"crossing $(r_\times, T_\times)$")
ax.plot([4.31], [7.07], "*", c="0.4", ms=9, label="analytic prediction")
ax.set_xlabel(r"$r = J_4/J_2$"); ax.set_ylabel(r"$T/G$")
ax.legend(fontsize=6, loc="upper left", frameon=False)
plt.tight_layout()
(REPO / "figures").mkdir(exist_ok=True)
plt.savefig(REPO / "figures" / "n1_phase_diagram.png", dpi=200)
