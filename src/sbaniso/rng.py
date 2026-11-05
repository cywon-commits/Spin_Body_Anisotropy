"""Per-replica random streams: xoshiro256** inside numba kernels.

Every replica owns a state row st = uint64[4]; kernels draw only from the row they are given, so
replica updates can run in parallel (prange) and the result does not depend on thread scheduling.
Streams are seeded from numpy's SeedSequence(seed).spawn(n), which gives statistically independent
initial states.  Reference: Blackman & Vigna, "Scrambled linear pseudorandom number generators"
(2021); the same generator family as polyspin (Paper I).
"""
import numpy as np
from numba import njit

_U53 = 1.0 / 9007199254740992.0  # 2^-53


@njit(cache=True, inline="always")
def _rotl(x, k):
    return (x << np.uint64(k)) | (x >> np.uint64(64 - k))


@njit(cache=True)
def next_u64(st):
    s0 = st[0]; s1 = st[1]; s2 = st[2]; s3 = st[3]
    result = _rotl(s1 * np.uint64(5), 7) * np.uint64(9)
    t = s1 << np.uint64(17)
    s2 ^= s0
    s3 ^= s1
    s1 ^= s2
    s0 ^= s3
    s2 ^= t
    s3 = _rotl(s3, 45)
    st[0] = s0; st[1] = s1; st[2] = s2; st[3] = s3
    return result


@njit(cache=True)
def uniform(st):
    """Uniform double in [0, 1) with 53 random bits."""
    return np.float64(next_u64(st) >> np.uint64(11)) * _U53


@njit(cache=True)
def randint(st, n):
    """Integer in [0, n).  Bias <= n * 2^-53 (irrelevant for lattice sizes)."""
    return int(uniform(st) * n)


def make_streams(seed, n):
    """n independent xoshiro256** states (shape (n, 4), uint64) from one integer seed."""
    children = np.random.SeedSequence(int(seed)).spawn(int(n))
    st = np.array([c.generate_state(4, dtype=np.uint64) for c in children], dtype=np.uint64)
    for k in range(n):  # the all-zero state is the only invalid one
        if not st[k].any():
            st[k, 0] = np.uint64(1)
    return st
