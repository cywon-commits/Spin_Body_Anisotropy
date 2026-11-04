import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import i0
from scipy.optimize import brentq
J1, J2 = 160/81, 128/81
r = np.array([3.8, 4.1, 4.4])
TI = np.array([6.807, 6.939, 6.963])          # Binder 48-64
TK64 = np.array([6.557, 6.876, 7.173])        # Upsilon = 2T/pi, L=64
TKWM = np.array([6.435, 6.745, 7.045])          # Weber-Minnhagen
rr = np.linspace(3.6, 4.6, 100)
# first-order (Yosefin-Domany type) lines with exact critical correlations away from crossing
TKT_line = 0.8929 * (rr * J2 + J2 / np.sqrt(2))   # sigma critical-ish: lower bound form
fig, ax = plt.subplots(figsize=(3.4, 2.8))
ax.plot(r, TI, "o-", c="tab:red", label=r"Ising $T_I$ (Binder, $L$=48/64)")
ax.plot(r, TK64, "s--", c="tab:blue", mfc="none", label=r"KT, $\Upsilon=2T/\pi$ at $L$=64")
ax.plot(r, TKWM, "s-", c="tab:blue", label="KT, Weber–Minnhagen fit")
ax.plot([4.31], [7.07], "k*", ms=10, label=r"analytic $(r_\times, T_\times)$")
ax.set_xlabel(r"$r=J_4/J_2$"); ax.set_ylabel(r"$T/G$")
ax.legend(fontsize=6.5, loc="lower right"); ax.set_title("PT pilot (fixed cluster counts): lines cross between r = 4.1 and 4.4", fontsize=8)
plt.tight_layout(); plt.savefig("pilot_crossing.png", dpi=200)
