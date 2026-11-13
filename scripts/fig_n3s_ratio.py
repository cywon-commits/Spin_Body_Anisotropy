"""Figure: full-model transition temperature over the effective-model prediction vs the confinement
length l_Gamma = sqrt(rho/|Gamma|) (coupling scale s = 1, 0.3, 0.1, 0.03; a4 = -0.02, K = 2), and the
axis-selection indicator <K> at T_Ups^full.  Error bars: half the random/ordered bracket.

    python scripts/fig_n3s_ratio.py   -> figures/n3s_ratio.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
S = json.loads((REPO / "data" / "n3s" / "n3s_summary.json").read_text())["s"]
ss = ("1.0", "0.3", "0.1", "0.03")
l = np.array([S[s]["result"]["l_Gamma"] for s in ss])
Ru = np.array([S[s]["result"]["T_ups"]["ratio"] for s in ss]); dRu = np.array([S[s]["result"]["T_ups"]["ratio_err"] for s in ss])
Rc = np.array([S[s]["result"]["T_C"]["ratio"] for s in ss]); dRc = np.array([S[s]["result"]["T_C"]["ratio_err"] for s in ss])
K = np.array([S[s]["result"]["K_at"]["1.0"] for s in ss])
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.8))
a.errorbar(l, Ru, dRu, fmt="o-", ms=4, capsize=2, label=r"$T_\Upsilon^{full}/T_{eff}$")
a.errorbar(l * 1.03, Rc, dRc, fmt="s--", ms=3, capsize=2, mfc="none", label=r"$T_C^{full}/T_{eff}$")
a.axhspan(0.95, 1.05, color="0.9", zorder=0); a.axhline(1, c="k", lw=0.6)
for s, x, y in zip(ss, l, Ru):
    a.annotate(f"s={s}", (x, y), textcoords="offset points", xytext=(4, -10), fontsize=6)
a.set_xscale("log"); a.set_xlabel(r"$\ell_\Gamma=\sqrt{\rho/|\Gamma|}$"); a.set_ylabel("ratio"); a.legend(fontsize=6.5, frameon=False)
b.plot(l, K, "o-", ms=4); b.axhline(1 / 3, c="k", lw=0.6, ls="--"); b.axhline(0.36, c="tab:red", lw=0.6, ls=":")
b.set_xscale("log"); b.set_xlabel(r"$\ell_\Gamma$"); b.set_ylabel(r"$\langle K\rangle$ at $T_\Upsilon^{full}$")
for ax in (a, b):
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([0.77, 1.41, 2.44, 4.46]); ax.set_xticklabels(["0.77", "1.41", "2.44", "4.46"])
plt.tight_layout(); plt.savefig(REPO / "figures" / "n3s_ratio.png", dpi=200)
