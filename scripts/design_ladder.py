"""Design PT temperature ladders from a pilot run's specific heat.

    python scripts/design_ladder.py runs/n2/pilot_r5.0_L32.h5 --Tmin 6.75 --Tmax 8.45 --L 32 48 64 96

For Gaussian energy distributions with standard deviation sigma_E at both neighbours, the swap
acceptance per attempt is A = erfc(dbeta * sigma_E / 2), sigma_E^2 = N T^2 C(T) (C per site).
The ladder is built by stepping dbeta = 2 erfcinv(A_target) / sigma_E(beta) from beta_max down to
beta_min, then spreading the steps uniformly in the stepping variable so that Tmin and Tmax are
hit exactly.  C(T) is the pilot's reweighted specific heat (pilot size L0), scaled for size L by
g(L) = 1 + growth * ln(L / L0) to allow for the logarithmic growth of the specific-heat peak
(conservative: applied at every T).  Target 0.3 leaves margin above the 0.2 requirement
(TASK.md, section 1).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc, erfcinv

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sbaniso.analysis import curves, load_run  # noqa: E402


def heat_curve(pilot):
    run = load_run(pilot)
    T = np.linspace(run.Ts.min(), run.Ts.max(), 400)
    v, _, _, _ = curves(run, T)
    return run.L, T, v["C"]


def design(T, C, L0, L, Tmin, Tmax, target=0.3, growth=0.21):
    g = 1 + growth * np.log(L / L0)
    Cf = lambda t: np.interp(t, T, C) * g
    sigma = lambda t: np.sqrt(L * L * t * t * Cf(t))
    # cumulative "distance" s(beta) = int sigma dbeta / (2 erfcinv(target)); one unit per step
    bs = np.linspace(1 / Tmax, 1 / Tmin, 4000)
    ds = np.array([sigma(1 / b) for b in bs]) / (2 * erfcinv(target))
    s = np.concatenate([[0], np.cumsum(0.5 * (ds[1:] + ds[:-1]) * np.diff(bs))])
    n = int(np.ceil(s[-1])) + 1
    beta = np.interp(np.linspace(0, s[-1], n), s, bs)
    Ts = np.sort(1 / beta)
    acc = [float(erfc((1 / Ts[k] - 1 / Ts[k + 1]) * 0.5 * (sigma(Ts[k]) + sigma(Ts[k + 1])) / 2))
           for k in range(n - 1)]
    return Ts, np.array(acc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pilot"); ap.add_argument("--Tmin", type=float, required=True)
    ap.add_argument("--Tmax", type=float, required=True); ap.add_argument("--L", type=int, nargs="+")
    ap.add_argument("--target", type=float, default=0.3)
    a = ap.parse_args()
    L0, T, C = heat_curve(a.pilot)
    for L in a.L:
        Ts, acc = design(T, C, L0, L, a.Tmin, a.Tmax, a.target)
        print(f"L = {L}: {len(Ts)} replicas, predicted acceptance {acc.min():.3f}-{acc.max():.3f}")


if __name__ == "__main__":
    main()
