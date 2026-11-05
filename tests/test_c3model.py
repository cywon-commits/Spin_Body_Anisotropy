import numpy as np
import pytest

from sbaniso.c3model import (couplings_from_body, couplings_from_r, energy, full_step, local_e,
                             measure, new_state, sweep, wolff_phi, wolff_sigma_phi,
                             wolff_sigma_psi, work_arrays, J1_PER_G, J2_PER_G, J4_PER_H4, NOBS)
from sbaniso.rng import make_streams

RNG = np.random.default_rng(0)


def test_rigid_couplings():
    c = couplings_from_r(4.0)
    assert np.allclose(c, [160 / 81, 128 / 81, 4 * 128 / 81])
    # r = (4/9) zH4/zG
    c2 = couplings_from_body(zG=1.0, zH4=9.0)
    assert np.isclose(c2[2] / c2[1], 4.0)


@pytest.mark.parametrize("r", [0.0, 1.0, 4.3])
def test_local_energy_matches_total(r):
    L = 6; c = couplings_from_r(r); s, p = new_state(L, RNG)
    E0 = energy(s, p, L, c)
    for _ in range(50):
        i, j = RNG.integers(L, size=2)
        sn = -s[i, j] if RNG.random() < 0.5 else s[i, j]
        pn = p[i, j] + RNG.normal()
        dE = local_e(s, p, i, j, L, sn, pn, c) - local_e(s, p, i, j, L, s[i, j], p[i, j], c)
        s2, p2 = s.copy(), p.copy(); s2[i, j] = sn; p2[i, j] = pn
        assert np.isclose(energy(s2, p2, L, c) - E0, dE, atol=1e-10)


def test_measure_energy_and_symmetries():
    L = 8; c = couplings_from_r(2.5); s, p = new_state(L, RNG); o = np.zeros(NOBS)
    measure(s, p, L, c, o)
    assert np.isclose(o[0], energy(s, p, L, c))
    E = energy(s, p, L, c)
    assert np.isclose(energy(s, p + 0.7, L, c), E)                 # U(1)
    assert np.isclose(energy(-s, p + np.pi, L, c), E)              # Z2 at fixed psi
    # psi-frame form of the Hamiltonian
    psi = p - np.pi * (1 - s) / 2
    Epsi = 0.0
    for i in range(L):
        for j in range(L):
            for a, b in (((i + 1) % L, j), (i, (j + 1) % L)):
                ss = s[i, j] * s[a, b]
                Epsi -= c[0] * ss + (c[2] + c[1] * ss) * np.cos(psi[i, j] - psi[a, b])
    assert np.isclose(Epsi, E)


def test_decoupled_limits_are_exact():
    """J4 = 0: energy depends on s only through J1 s s' after phi -> psi?  No: J4 = 0 means
    H = -J1 ss' - J2 cos(dphi): s and phi decouple.  J2 = 0: s and psi decouple."""
    L = 6; s, p = new_state(L, RNG)
    c = np.array([1.3, 0.8, 0.0])
    E = energy(s, p, L, c)
    Eis = -c[0] * sum(s[i, j] * (s[(i + 1) % L, j] + s[i, (j + 1) % L]) for i in range(L) for j in range(L))
    assert np.isclose(E, Eis + energy(np.ones_like(s), p, L, np.array([0, c[1], 0.0])))
    c = np.array([1.3, 0.0, 0.8]); psi = p - np.pi * (1 - s) / 2
    E = energy(s, p, L, c)
    assert np.isclose(E, Eis + energy(np.ones_like(s), psi, L, np.array([0, c[2], 0.0])))


def test_sigma_psi_cluster_preserves_psi():
    L = 10; c = couplings_from_r(4.1); s, p = new_state(L, RNG); stack, inc, mem = work_arrays(L)
    st = make_streams(3, 1)[0]
    for _ in range(20):
        psi0 = np.exp(1j * p) * s
        wolff_sigma_psi(s, p, L, 1 / 7.0, c, stack, inc, mem, st)
        assert np.allclose(np.exp(1j * p) * s, psi0)
        assert not inc.any()


def test_phi_cluster_preserves_sigma_and_sigma_phi_preserves_phi():
    L = 10; c = couplings_from_r(4.1); s, p = new_state(L, RNG); stack, inc, mem = work_arrays(L)
    st = make_streams(4, 1)[0]
    s0 = s.copy(); wolff_phi(s, p, L, 1 / 7.0, c, stack, inc, mem, st); assert (s == s0).all()
    p0 = p.copy(); wolff_sigma_phi(s, p, L, 1 / 7.0, c, stack, inc, mem, st); assert np.allclose(p, p0)


def _series(flags, c, L, T, n, nskip, sd):
    s, p = new_state(L, np.random.default_rng(sd)); stack, inc, mem = work_arrays(L)
    st = make_streams(sd, 1)[0]; o = np.zeros(NOBS); xs = []
    for it in range(n + nskip):
        if flags is None:
            sweep(s, p, L, 1 / T, c, 1.0, st)
        else:
            full_step(s, p, L, 1 / T, c, 1.0, stack, inc, mem, np.array(flags, dtype=np.int64), np.zeros(3), st)
        if it >= nskip:
            measure(s, p, L, c, o); xs.append([o[0] / L ** 2, (o[1] / L ** 2) ** 2, (o[8] ** 2 + o[9] ** 2) / L ** 4])
    xs = np.array(xs); nb = 40
    b = xs[: len(xs) // nb * nb].reshape(nb, -1, 3).mean(1)
    return b.mean(0), b.std(0) / np.sqrt(nb)


@pytest.mark.parametrize("flags", [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (3, 3, 3)])
def test_cluster_updates_sample_boltzmann(flags):
    """Each cluster type (on top of Metropolis) must reproduce Metropolis-only averages of
    E, m_sigma^2 and |M_psi|^2 on a small lattice near the crossing."""
    c = couplings_from_r(4.3); L = 6; T = 7.0
    ref, dref = _series(None, c, L, T, 120000, 5000, 11)
    val, dval = _series(flags, c, L, T, 20000, 1000, 12)
    z = np.abs(val - ref) / np.sqrt(dref ** 2 + dval ** 2)
    assert (z < 4.5).all(), (flags, val, ref, z)


def _exact_ising(L, T):
    import itertools
    st = np.array(list(itertools.product([-1, 1], repeat=L)))
    row = np.array([sum(x[i] * x[(i + 1) % L] for i in range(L)) for x in st])
    inter = st @ st.T; mag = st.sum(1)
    def lnZ(b, h=0.0):
        W = np.exp(b * (row[:, None] + inter) + h * mag[:, None])
        return np.log(np.trace(np.linalg.matrix_power(W, L)))
    b = 1 / T; e = 1e-5; h = 1e-3
    E = -(lnZ(b + e) - lnZ(b - e)) / (2 * e) / L ** 2
    m2 = (lnZ(b, h) + lnZ(b, -h) - 2 * lnZ(b)) / h ** 2 / L ** 4
    return E, m2


@pytest.mark.parametrize("flags", [None, (1, 0, 0), (0, 1, 0), (2, 2, 0)])
def test_exact_ising_6x6(flags):
    """J2 = J4 = 0 is the pure Ising model: compare with the exact 6x6 transfer-matrix result.
    This test fails (by > 50 sigma) for the 'clusters until N sites flipped' stopping rule."""
    L, T = 6, 2.5
    E, m2 = _exact_ising(L, T)
    val, err = _series(flags, np.array([1.0, 0.0, 0.0]), L, T, 60000 if flags is None else 30000, 2000, 5)
    assert abs(val[0] - E) < 4.5 * err[0] + 1e-3, (val, E)
    assert abs(val[1] - m2) < 4.5 * err[1] + 1e-3, (val, m2)
