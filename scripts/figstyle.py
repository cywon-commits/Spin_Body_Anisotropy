"""Common style for the manuscript figures (fig1-fig5): single column 3.4 in, English labels,
Okabe-Ito colour-blind-safe palette, panel labels (a), (b), (c)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

WIDTH = 3.4
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
          blue="#0072B2", red="#D55E00", purple="#CC79A7", grey="#999999")
SIZE_MARKERS = ("o", "s", "^", "D", "v")

plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "lines.markersize": 3.5, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False, "savefig.dpi": 300, "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans", "errorbar.capsize": 1.5,
})


def size_colors(n):
    """Sequential colours for n system sizes (viridis, light end trimmed)."""
    cm = plt.get_cmap("viridis")
    return [cm(0.85 * i / max(n - 1, 1)) for i in range(n)]


def panel(ax, label, x=0.02, y=0.97, **kw):
    ax.text(x, y, f"({label})", transform=ax.transAxes, ha="left", va="top", fontweight="bold", **kw)


def save(fig, name, repo):
    (repo / "figures").mkdir(exist_ok=True)
    fig.savefig(repo / "figures" / f"{name}.pdf")
    fig.savefig(repo / "figures" / f"{name}.png", dpi=300)
