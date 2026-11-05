"""N2 configs (docs/N2_analysis_plan.md, section 1): designed constant-acceptance ladders with the
replica counts of configs/T1_runs.csv, written as explicit temperature lists to configs/n2/.

    python scripts/make_n2_configs.py     (needs runs/n2/pilot_r*_L32.h5)
"""
import csv
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.special import erfc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from design_ladder import heat_curve  # noqa: E402

RANGES = {5.0: (6.75, 8.45), 6.0: (7.0, 9.75)}
GROWTH = 0.21


def ladder_fixed_n(T, C, L0, L, Tmin, Tmax, n):
    """n temperatures with equal predicted acceptance (same construction as design_ladder.design)."""
    g = 1 + GROWTH * np.log(L / L0)
    sigma = lambda t: np.sqrt(L * L * t * t * np.interp(t, T, C) * g)
    bs = np.linspace(1 / Tmax, 1 / Tmin, 4000)
    ds = np.array([sigma(1 / b) for b in bs])
    s = np.concatenate([[0], np.cumsum(0.5 * (ds[1:] + ds[:-1]) * np.diff(bs))])
    Ts = np.sort(1 / np.interp(np.linspace(0, s[-1], n), s, bs))
    acc = [erfc((1 / Ts[k] - 1 / Ts[k + 1]) * 0.5 * (sigma(Ts[k]) + sigma(Ts[k + 1])) / 2) for k in range(n - 1)]
    return Ts, float(min(acc))


def main():
    rows = [r for r in csv.DictReader(open(REPO / "configs" / "T1_runs.csv")) if r["task"] == "N2"]
    out = REPO / "configs" / "n2"; out.mkdir(exist_ok=True)
    heat = {r: heat_curve(REPO / f"runs/n2/pilot_r{r}_L32.h5") for r in RANGES}
    order = []
    for i, row in enumerate(rows):
        r = float(row["r"]); L = int(row["L"]); n = int(row["nrep"])
        L0, T, C = heat[r]
        Ts, amin = ladder_fixed_n(T, C, L0, L, *RANGES[r], n)
        name = f"r{r}_L{L}"
        d = dict(name=f"n2_{name}", model=dict(r=r), L=L,
                 temperatures=dict(list=[round(float(t), 6) for t in Ts]),
                 ntherm=int(row["ntherm"]), nmeas=int(row["nmeas"]), clusters=[1, 1, 1],
                 seed=40000 + i, start="random", checkpoint_every=1000, output=f"runs/n2/{name}.h5")
        (out / f"{name}.yaml").write_text(yaml.safe_dump(d, sort_keys=False))
        order.append((L, r, name)); print(f"{name}: {n} replicas, predicted min acceptance {amin:.3f}")
        if r == 6.0 and L == 64:
            d2 = dict(d, name="n2_r6.0_L64_ordered", seed=40100, start="ordered", output="runs/n2/r6.0_L64_ordered.h5")
            (out / "r6.0_L64_ordered.yaml").write_text(yaml.safe_dump(d2, sort_keys=False))
            order.append((L, r + 0.01, "r6.0_L64_ordered"))
    order.sort(key=lambda x: (-x[0], x[1]))     # largest first: better packing on 4 parallel slots
    (out / "ORDER").write_text("".join(f"configs/n2/{n}.yaml\n" for *_, n in order))


if __name__ == "__main__":
    main()
