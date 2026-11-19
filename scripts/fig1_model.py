"""Fig. 1: a tetrahedral spin body under uniaxial anisotropy keeps one body axis along z, and its
in-plane state is an Ising variable (octupole sign sigma) times an XY phase (hexadecapole phase psi).
(a) Spin body (four unit spins from the centre to the vertices) with a C3 axis along z (selected for
    Gamma < 0) and one of the three S4 axes; V = sum_i [a2 P2(z.S_i) + a4 P4(z.S_i)].
(b) sigma = sign of the l = 3, m = 0 octupole q30: apex up or down; time reversal S -> -S inverts the body
    (apex down, base turned by pi), so sigma is TR odd.
(c) Top view: rotation alpha of the body about z; the l = 3, |m| = 3 octupole phase phi = 3 alpha
    (+ pi for sigma = -1) and the l = 4, |m| = 3 hexadecapole phase psi = phi - pi(1 - sigma)/2 = 3 alpha
    are invariant under alpha -> alpha + 2 pi/3.  Inset: multipole weights of the C3 body (derivation note
    Sec. 5) and the couplings they give.  Geometry only; no data file is read.

    python scripts/fig1_model.py   -> figures/fig1_model.{pdf,png}
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from figstyle import OI, WIDTH, panel, plt, save  # noqa: E402

AZ, EL = np.radians(-28), np.radians(16)


def body(alpha=0.25, up=True):
    s = 1 if up else -1
    V = [np.array([0, 0, s])]
    for k in range(3):
        a = alpha + 2 * np.pi * k / 3
        V.append(np.array([np.sqrt(8) / 3 * np.cos(a), np.sqrt(8) / 3 * np.sin(a), -s / 3]))
    return np.array(V)


def proj(P):
    P = np.atleast_2d(P)
    x = P[:, 0] * np.cos(AZ) - P[:, 1] * np.sin(AZ)
    d = P[:, 0] * np.sin(AZ) + P[:, 1] * np.cos(AZ)
    y = P[:, 2] * np.cos(EL) + d * np.sin(EL)
    return np.c_[x, y], d                                  # d > 0: towards the back


def draw_body(ax, V, c0=(0, 0), scale=1.0, spin_col=OI["blue"], lw=0.8):
    P, depth = proj(V); P = P * scale + np.array(c0)
    for i in range(4):
        for j in range(i + 1, 4):
            back = depth[i] > 0.3 and depth[j] > 0.3
            ax.plot(*P[[i, j]].T, c="0.55", lw=0.6, ls=":" if back else "-", zorder=1)
    for i in range(4):
        ax.annotate("", xy=P[i], xytext=c0, arrowprops=dict(arrowstyle="-|>", lw=lw, color=spin_col,
                                                            mutation_scale=7, shrinkA=0, shrinkB=0), zorder=3)
    ax.plot(*c0, "o", c="k", ms=2, zorder=4)
    return P


fig = plt.figure(figsize=(WIDTH, 4.3))
a = fig.add_axes([0.0, 0.47, 0.55, 0.53]); b = fig.add_axes([0.55, 0.47, 0.45, 0.53])
c = fig.add_axes([0.0, 0.0, 1.0, 0.43])
for ax in (a, b, c):
    ax.set_aspect("equal"); ax.axis("off")

V = body()
draw_body(a, V, spin_col=OI["blue"])
z0, z1 = proj(np.array([[0, 0, -1.25], [0, 0, 1.45]]))[0]
a.annotate("", xy=z1, xytext=z0, arrowprops=dict(arrowstyle="-|>", color=OI["red"], lw=1.1, mutation_scale=8), zorder=2)
a.text(z1[0] + 0.07, z1[1] - 0.08, r"$\hat z\parallel C_3$", color=OI["red"], fontsize=7)
m1 = 0.5 * (V[0] + V[1]); m2 = 0.5 * (V[2] + V[3]); u = (m1 - m2) / np.linalg.norm(m1 - m2)
s0, s1 = proj(np.array([-1.25 * u, 1.25 * u]))[0]
a.plot(*np.c_[s0, s1], c=OI["green"], lw=0.9, ls="--", zorder=2)
a.text(s1[0] - 0.05, s1[1] + 0.06, r"$S_4$", color=OI["green"], fontsize=7, ha="right")
a.text(0.0, -1.33, r"$V=\sum_i[a_2P_2(\hat z\cdot\mathbf{S}_i)+a_4P_4(\hat z\cdot\mathbf{S}_i)]$" "\n"
       r"$E_{\rm eff}=-\Gamma\sum_\alpha n_\alpha^4$;  $\Gamma<0$: $C_3\parallel\hat z$", ha="center", va="top", fontsize=5.6)
a.set_xlim(-1.35, 1.35); a.set_ylim(-1.75, 1.6)

for x0, up, lab in ((-0.7, True, r"$\sigma=+1$"), (0.7, False, r"$\sigma=-1$")):
    draw_body(b, body(0.25 if up else 0.25 + np.pi, up), c0=(x0, 0.15), scale=0.55, spin_col=OI["blue"] if up else OI["orange"])
    b.text(x0, -0.75, lab, ha="center", fontsize=7)
b.annotate("", xy=(0.3, 0.15), xytext=(-0.3, 0.15), arrowprops=dict(arrowstyle="<|-|>", lw=0.7, mutation_scale=6))
b.text(0.0, 0.25, r"$\mathcal{T}$", ha="center", fontsize=7)
b.text(0.0, -1.05, r"$\sigma=\mathrm{sgn}\,q_{30}$ (octupole, TR odd)", ha="center", fontsize=6)
b.set_xlim(-1.35, 1.35); b.set_ylim(-1.75, 1.6)

al = 0.35; R = 0.8; cx = -1.05
th = np.linspace(0, 2 * np.pi, 100)
c.plot(cx + R * np.cos(th), R * np.sin(th), c="0.8", lw=0.5)
for k in range(3):
    t = al + 2 * np.pi * k / 3
    c.annotate("", xy=(cx + R * np.cos(t), R * np.sin(t)), xytext=(cx, 0),
               arrowprops=dict(arrowstyle="-|>", color=OI["blue"], lw=0.8, mutation_scale=7, shrinkA=0, shrinkB=0))
c.plot(cx, 0, "o", c="k", ms=3); c.text(cx - 0.08, -0.2, "apex", fontsize=5.5, ha="right")
c.plot([cx, cx + 1.0], [0, 0], c="k", lw=0.5)
tt = np.linspace(0, al, 20); c.plot(cx + 0.45 * np.cos(tt), 0.45 * np.sin(tt), c=OI["red"], lw=0.8)
c.text(cx + 0.5, 0.07, r"$\alpha$", color=OI["red"], fontsize=7)
c.text(cx, -0.98, "top view", ha="center", fontsize=6)
txt = (r"$\varphi=3\alpha+\pi(1-\sigma)/2$  ($q_{33}$, octupole)" "\n"
       r"$\psi=\varphi-\pi(1-\sigma)/2=3\alpha$  ($q_{43}$, hexadecapole)" "\n"
       r"invariant under $\alpha\to\alpha+2\pi/3$")
c.text(-0.1, 0.62, txt, fontsize=5.8, va="top", linespacing=1.5)
rows = [(r"$l=3,\ m=0$", "400/81", r"$J_1\sigma\sigma'$"), (r"$l=3,\ |m|=3$", "320/81", r"$J_2\cos\Delta\varphi$"),
        (r"$l=4,\ |m|=3$", "2240/729", r"$J_4\sigma\sigma'\cos\Delta\varphi$")]
y = -0.1
c.text(-0.1, y, "channel", fontsize=5.6, weight="bold"); c.text(0.85, y, "weight", fontsize=5.6, weight="bold")
c.text(1.55, y, "coupling", fontsize=5.6, weight="bold")
for i, (ch, w, j) in enumerate(rows):
    yy = y - 0.24 * (i + 1)
    c.text(-0.1, yy, ch, fontsize=5.6); c.text(0.85, yy, w, fontsize=5.6); c.text(1.55, yy, j, fontsize=5.6)
c.set_xlim(-1.95, 2.75); c.set_ylim(-1.05, 0.9)
for x, y, l in ((0.01, 0.99, "a"), (0.56, 0.99, "b"), (0.01, 0.435, "c")):
    fig.text(x, y, f"({l})", ha="left", va="top", fontweight="bold")
save(fig, "fig1_model", REPO)
