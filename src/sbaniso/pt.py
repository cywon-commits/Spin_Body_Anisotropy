"""Parallel tempering engine for the C3 effective model.

The run is advanced in chunks by `pt_chunk` (numba, replicas updated in parallel with prange);
everything needed to continue a run lives in a `PTState` so that a run can be checkpointed and
resumed bit-identically (sbaniso.run handles the YAML -> HDF5 side).

One sweep (iteration `it`):
  1. every temperature slot k updates its replica perm[k]: Metropolis sweep + a FIXED number of
     clusters of each enabled type (ncl[k]); each replica draws from its own xoshiro stream,
  2. every replica is measured (13 observables, c3model.OBS_NAMES),
  3. neighbour swaps on alternating even/odd pairs (master stream).
Tuning (equilibration only):  ncl[k] is re-estimated during the first half of the equilibration
(target: ~N flipped sites per cluster type and sweep) and then FROZEN; the Metropolis step deltas[k]
is tuned during the whole equilibration toward an acceptance in [0.40, 0.55].
Round trips: a replica is labelled 'up' when it visits the lowest temperature and 'down' when it
visits the highest; a full round trip (bottom -> top -> bottom) is counted when an 'up'-labelled
replica that has been at the top returns to the bottom.  Pilot v3 counted bottom -> top passages
(half round trips) as 'round trips'.
"""
import numpy as np
from numba import njit, prange

from .c3model import NOBS, full_step, measure, new_state
from .rng import make_streams, uniform


@njit(cache=True, parallel=True)
def pt_chunk(S, P, L, betas, c, flags, deltas, ncl, acc_sizes, acc_calls, perm, rng, rng_master,
             label, rt_rep, swap_tries, swap_acc, metro_acc, it0, nsweeps, ntherm, stacks, incs,
             mems, out):
    """Advance sweeps it0 .. it0+nsweeps-1.  Measurements of sweeps it >= ntherm are written to
    out[k, it - max(it0, ntherm), :] (out has room for the measured sweeps of this chunk).
    Returns the number of rows written."""
    R = betas.shape[0]
    N = L * L
    E = np.zeros(R)
    obs = np.zeros((R, NOBS))
    sizes = np.zeros((R, 3))
    row = 0
    for it in range(it0, it0 + nsweeps):
        for k in prange(R):
            rep = perm[k]
            for q in range(3):
                sizes[k, q] = 0.0
            a = full_step(S[rep], P[rep], L, betas[k], c, deltas[k], stacks[k], incs[k], mems[k],
                          ncl[k], sizes[k], rng[rep])
            if it < ntherm // 2:
                for q in range(3):
                    if flags[q]:
                        acc_sizes[k, q] += sizes[k, q]
                        acc_calls[k, q] += ncl[k, q]
                        if it > 50:
                            mean = acc_sizes[k, q] / acc_calls[k, q]
                            ncl[k, q] = max(1, int(N / mean + 0.5))
            if it < ntherm:
                if a > 0.55:
                    deltas[k] = min(np.pi, deltas[k] * 1.02)
                elif a < 0.40:
                    deltas[k] *= 0.98
            else:
                metro_acc[k] += a
            measure(S[rep], P[rep], L, c, obs[k])
            E[k] = obs[k, 0]
        if it >= ntherm:
            for k in range(R):
                for q in range(NOBS):
                    out[k, row, q] = obs[k, q]
            row += 1
        for k in range(it % 2, R - 1, 2):
            d = (betas[k + 1] - betas[k]) * (E[k + 1] - E[k])
            ok = d >= 0 or uniform(rng_master) < np.exp(d)
            if it >= ntherm:
                swap_tries[k] += 1
            if ok:
                tmp = perm[k]; perm[k] = perm[k + 1]; perm[k + 1] = tmp
                tmpe = E[k]; E[k] = E[k + 1]; E[k + 1] = tmpe
                if it >= ntherm:
                    swap_acc[k] += 1
        # round trips (counted during measurement only): label 1 = last extreme visited was bottom
        # (k = 0), 2 = it has been at the top since; back at the bottom with label 2 -> one trip
        b = perm[0]; t = perm[R - 1]
        if it >= ntherm and label[b] == 2:
            rt_rep[b] += 1
        label[b] = 1
        if label[t] == 1:
            label[t] = 2
    return row


class PTState:
    """Complete dynamic state of a PT run (everything a checkpoint must hold)."""

    FIELDS = ("S", "P", "deltas", "ncl", "acc_sizes", "acc_calls", "perm", "rng", "rng_master",
              "label", "rt_rep", "swap_tries", "swap_acc", "metro_acc", "it")

    def __init__(self, L, Ts, seed, flags, ordered=False):
        R = len(Ts)
        init = np.random.default_rng(np.random.SeedSequence(int(seed)).spawn(1)[0])
        self.S = np.empty((R, L, L), dtype=np.int64)
        self.P = np.empty((R, L, L))
        for k in range(R):
            self.S[k], self.P[k] = new_state(L, init, ordered)
        self.deltas = np.ones(R)
        self.ncl = np.array([list(flags)] * R, dtype=np.int64)
        self.acc_sizes = np.zeros((R, 3))
        self.acc_calls = np.zeros((R, 3))
        self.perm = np.arange(R, dtype=np.int64)
        streams = make_streams(seed, R + 1)
        self.rng = np.ascontiguousarray(streams[:R])
        self.rng_master = streams[R].copy()
        self.label = np.zeros(R, dtype=np.int64)
        self.rt_rep = np.zeros(R, dtype=np.int64)
        self.swap_tries = np.zeros(R - 1, dtype=np.int64)
        self.swap_acc = np.zeros(R - 1, dtype=np.int64)
        self.metro_acc = np.zeros(R)
        self.it = 0

    def to_dict(self):
        return {f: getattr(self, f) for f in self.FIELDS}

    @classmethod
    def from_dict(cls, d):
        self = cls.__new__(cls)
        for f in cls.FIELDS:
            v = d[f]
            setattr(self, f, int(v) if f == "it" else np.array(v))
        return self


def work_arrays(R, L):
    N = L * L
    return (np.zeros((R, N + 4, 2), dtype=np.int64), np.zeros((R, L, L), dtype=np.bool_),
            np.zeros((R, N + 4, 2), dtype=np.int64))


def advance(state, L, betas, c, flags, ntherm, nsweeps, work=None):
    """Advance `state` by nsweeps sweeps; return the measured rows (R, n_meas_rows, NOBS)."""
    R = len(betas)
    stacks, incs, mems = work if work is not None else work_arrays(R, L)
    nmeas_rows = max(0, state.it + nsweeps - max(state.it, ntherm))
    out = np.zeros((R, nmeas_rows, NOBS))
    rows = pt_chunk(state.S, state.P, L, np.asarray(betas, float), np.asarray(c, float),
                    np.asarray(flags, np.int64), state.deltas, state.ncl, state.acc_sizes,
                    state.acc_calls, state.perm, state.rng, state.rng_master, state.label,
                    state.rt_rep, state.swap_tries, state.swap_acc, state.metro_acc, state.it,
                    nsweeps, ntherm, stacks, incs, mems, out)
    assert rows == nmeas_rows
    state.it += nsweeps
    return out


def geometric_ladder(Tmin, Tmax, R):
    return Tmin * (Tmax / Tmin) ** (np.arange(R) / (R - 1))
