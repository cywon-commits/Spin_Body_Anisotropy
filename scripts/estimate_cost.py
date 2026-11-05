"""Cost estimate for a run list (TASK N1: re-estimate before running).

    python scripts/estimate_cost.py configs/T1_runs.csv [--task N1]

Benchmarks one replica-sweep (Metropolis + all three cluster types, ncl tuned as in production)
on one core at each L near the crossing (r = 4.3, T = 7.0), then sums R * (ntherm + nmeas) * t(L)
per row and by priority.
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numba
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sbaniso.c3model import couplings_from_r  # noqa: E402
from sbaniso.pt import PTState, advance  # noqa: E402


def ms_per_replica_sweep(L, nrep=4, ntune=300, nbench=None):
    numba.set_num_threads(1)
    Ts = np.linspace(6.9, 7.1, nrep); c = couplings_from_r(4.3); flags = (1, 1, 1)
    st = PTState(L, Ts, 1, flags)
    advance(st, L, 1 / Ts, c, flags, ntune, ntune)          # equilibrate + tune ncl
    nbench = nbench or max(20, int(200 * (32 / L) ** 2))
    t0 = time.time()
    advance(st, L, 1 / Ts, c, flags, ntune, nbench)
    return 1e3 * (time.time() - t0) / (nrep * nbench)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("--task", default="N1"); ap.add_argument("--cores", type=int, default=8)
    a = ap.parse_args()
    rows = [r for r in csv.DictReader(open(a.csv)) if r["task"] == a.task]
    Ls = sorted({int(r["L"]) for r in rows})
    ms_per_replica_sweep(16, nrep=2, ntune=20, nbench=5)    # JIT warm-up
    t = {L: ms_per_replica_sweep(L) for L in Ls}
    for L in Ls:
        print(f"L = {L:4d}: {t[L]:6.2f} ms per replica-sweep (1 core)")
    tot = {}
    for r in rows:
        L = int(r["L"]); h = int(r["nrep"]) * (int(r["ntherm"]) + int(r["nmeas"])) * t[L] / 3.6e6
        tot[r["priority"]] = tot.get(r["priority"], 0) + h
    for p, h in sorted(tot.items()):
        print(f"priority {p}: {h:6.1f} core-hours  (~{h / a.cores:.1f} h on {a.cores} cores at full efficiency)")
    print(f"total: {sum(tot.values()):.1f} core-hours")


if __name__ == "__main__":
    main()
