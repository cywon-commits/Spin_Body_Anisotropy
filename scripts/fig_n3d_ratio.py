"""Figure: T_Ups^full / T_eff before and after the Debye-Waller reduction of the effective couplings,
vs the confinement length l_Gamma (s = 1, 0.3, 0.1, 0.03).  Error bars: half the random/ordered bracket
of the full model.  Grey band: +-5 % (decision of docs/N3d_analysis_plan.md).

    python scripts/fig_n3d_ratio.py   -> figures/n3d_ratio.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
D = json.loads((REPO / "data" / "n3d" / "n3d_summary.json").read_text())["s"]
S = json.loads((REPO / "data" / "n3s" / "n3s_summary.json").read_text())["s"]
ss = ("1.0", "0.3", "0.1", "0.03")
l = np.array([S[s]["params"]["l_Gamma"] for s in ss])
R0 = np.array([S[s]["result"]["T_ups"]["ratio"] for s in ss]); dR0 = np.array([S[s]["result"]["T_ups"]["ratio_err"] for s in ss])
R1 = np.array([D[s]["ratio"] for s in ss]); dR1 = np.array([D[s]["ratio_err"] for s in ss])
Rk = np.array([D[s]["sensitivity_K0"]["ratio_K0"] for s in ss])
fig, ax = plt.subplots(figsize=(3.6, 2.8))
ax.axhspan(0.95, 1.05, color="0.9", zorder=0); ax.axhline(1, c="k", lw=0.6)
ax.errorbar(l, R0, dR0, fmt="o--", ms=4, capsize=2, c="tab:gray", label="bare couplings")
ax.errorbar(l, R1, dR1, fmt="o-", ms=4, capsize=2, c="tab:blue", label="Debye–Waller reduced")
ax.plot(l, Rk, "x", c="tab:blue", ms=5, label=r"reduced, $K_0$ baseline")
ax.set_xscale("log"); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(l); ax.set_xticklabels([f"{x:.2f}" for x in l])
ax.set_xlabel(r"$\ell_\Gamma=\sqrt{\rho/|\Gamma|}$"); ax.set_ylabel(r"$T_\Upsilon^{full}/T_{eff}$")
ax.legend(fontsize=6.5, frameon=False, loc="lower left")
plt.tight_layout(); plt.savefig(REPO / "figures" / "n3d_ratio.png", dpi=200)
