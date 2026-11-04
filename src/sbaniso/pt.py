"""Parallel tempering for the C3 effective model (pilot driver, to be extended to YAML/HDF5).

    python -m sbaniso.pt --r 4.1 --L 64 --nrep 26 --Tmin 6.3 --Tmax 8.0 \
                         --ntherm 1000 --nmeas 10000 --seed 1 --out pt_r4.1_L64.npz

Each sweep: every replica does full_step (Metropolis + clusters), all replicas are measured,
then neighbour swaps on alternating even/odd pairs.  Output: per-temperature time series of the
NOBS observables of c3model.measure, swap acceptance, round trips, provenance.
"""
import argparse
import subprocess
import time

import numpy as np
from numba import njit

from .c3model import NOBS, couplings_from_r, full_step, measure, new_state, seed, work_arrays


@njit(cache=True)
def run_pt(S, P, L, betas, c, deltas, ntherm, nmeas, flags, stack, inc, mem, out, swapacc, rt, perm,
           ncl):
    """ncl[k, :] = clusters per sweep at temperature k.  Tuned in the first half of the
    equilibration (target: ~N flipped sites per type and sweep), then frozen."""
    R = betas.shape[0]
    N = L * L
    E = np.zeros(R)
    sizes = np.zeros(3)
    acc_sizes = np.zeros((R, 3))
    acc_calls = np.zeros((R, 3))
    for k in range(R):
        for q in range(3):
            ncl[k, q] = flags[q]
    obs = np.zeros(NOBS)
    direction = np.zeros(R, dtype=np.int64)
    for it in range(ntherm + nmeas):
        for k in range(R):
            rep = perm[k]
            sizes[:] = 0.0
            a = full_step(S[rep], P[rep], L, betas[k], c, deltas[k], stack, inc, mem, ncl[k], sizes)
            if it < ntherm // 2:
                for q in range(3):
                    if flags[q]:
                        acc_sizes[k, q] += sizes[q]
                        acc_calls[k, q] += ncl[k, q]
                        if it > 50:
                            mean = acc_sizes[k, q] / acc_calls[k, q]
                            ncl[k, q] = max(1, int(N / mean + 0.5))
            if it < ntherm:
                if a > 0.55:
                    deltas[k] = min(np.pi, deltas[k] * 1.02)
                elif a < 0.40:
                    deltas[k] *= 0.98
            measure(S[rep], P[rep], L, c, obs)
            E[k] = obs[0]
            if it >= ntherm:
                t = it - ntherm
                for q in range(NOBS):
                    out[k, t, q] = obs[q]
        for k in range(it % 2, R - 1, 2):
            d = (betas[k + 1] - betas[k]) * (E[k + 1] - E[k])
            if d >= 0 or np.random.random() < np.exp(d):
                tmp = perm[k]; perm[k] = perm[k + 1]; perm[k + 1] = tmp
                tmpe = E[k]; E[k] = E[k + 1]; E[k + 1] = tmpe
                if it >= ntherm:
                    swapacc[k] += 1
        if it >= ntherm:
            if direction[perm[0]] == 0:
                direction[perm[0]] = 1
            if direction[perm[R - 1]] == 1:
                direction[perm[R - 1]] = 0
                rt[0] += 1




def geometric_ladder(Tmin, Tmax, R):
    return Tmin * (Tmax / Tmin) ** (np.arange(R) / (R - 1))


def git_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def run(c, L, Ts, ntherm, nmeas, seed_value, flags=(1, 1, 1), ordered=False):
    rng = np.random.default_rng(seed_value)
    seed(int(seed_value) % (2 ** 31))
    R = len(Ts)
    S = np.empty((R, L, L), dtype=np.int64); P = np.empty((R, L, L))
    for k in range(R):
        S[k], P[k] = new_state(L, rng, ordered)
    stack, inc, mem = work_arrays(L)
    out = np.zeros((R, nmeas, NOBS))
    swapacc = np.zeros(R - 1); rt = np.zeros(1, dtype=np.int64)
    perm = np.arange(R); deltas = np.ones(R)
    ncl = np.zeros((R, 3), dtype=np.int64)
    t0 = time.time()
    run_pt(S, P, L, 1.0 / np.asarray(Ts), np.asarray(c, float), deltas, ntherm, nmeas,
           np.array(flags, dtype=np.int64), stack, inc, mem, out, swapacc, rt, perm, ncl)
    return dict(out=out.astype(np.float32), swap=swapacc / nmeas, roundtrips=int(rt[0]), ncl=ncl,
                secs=time.time() - t0, Ts=np.asarray(Ts), c=np.asarray(c), L=L, seed=seed_value,
                ntherm=ntherm, nmeas=nmeas, flags=np.array(flags), git=git_hash())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=float, required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--nrep", type=int, required=True)
    ap.add_argument("--Tmin", type=float, required=True)
    ap.add_argument("--Tmax", type=float, required=True)
    ap.add_argument("--ntherm", type=int, default=1000)
    ap.add_argument("--nmeas", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flags", type=str, default="1,1,1", help="sigma_phi,sigma_psi,phi clusters")
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()
    c = couplings_from_r(a.r)
    Ts = geometric_ladder(a.Tmin, a.Tmax, a.nrep)
    res = run(c, a.L, Ts, a.ntherm, a.nmeas, a.seed, tuple(int(x) for x in a.flags.split(",")))
    out = a.out or f"pt_r{a.r}_L{a.L}.npz"
    np.savez_compressed(out, r=a.r, **res)
    print(f"r={a.r} L={a.L}: {res['secs']:.0f}s, round trips {res['roundtrips']}, "
          f"min swap {res['swap'].min():.2f} -> {out}")


if __name__ == "__main__":
    main()
