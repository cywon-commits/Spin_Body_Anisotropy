"""Full-model measurement (sbaniso.fullmodel) on configurations with known answers.  Skipped when
polyspin (Paper I, anisotropy branch) is not installed."""
import numpy as np
import pytest

polyspin = pytest.importorskip("polyspin")
from polyspin.groundstate import oriented_tetrahedron  # noqa: E402
from polyspin.lattice import Lattice  # noqa: E402
from polyspin.model import Params  # noqa: E402
from polyspin.stiffness import twisted_inter_energy  # noqa: E402
from scipy.spatial.transform import Rotation as Rot  # noqa: E402

from sbaniso.c3model import NOBS  # noqa: E402
from sbaniso.fullmodel import NEXTRA, _norms, measure_all  # noqa: E402


def _measure(spins_lat, P, L):
    lat = Lattice(2, L)
    out = np.zeros((1, NOBS)); ext = np.zeros((1, NEXTRA))
    measure_all(spins_lat[None], np.array([0]), np.array([1.234]), lat.nbr, P.array(), L,
                np.array(_norms()), out, ext)
    return out[0], ext[0]


def test_uniform_c3_state():
    L = 6; P = Params.from_z(zG=0.08, zH4=0.8, K=2.0, a4=-0.02)
    S = oriented_tetrahedron([1, 1, 1], 0.4)
    o, x = _measure(np.tile(S, (L * L, 1, 1)), P, L)
    N = L * L
    assert np.isclose(abs(o[1]), N)                              # |q30| = 1 per body
    assert np.isclose(np.hypot(o[6], o[7]), N) and np.isclose(np.hypot(o[8], o[9]), N)
    assert np.isclose(x[0], 1 / 3) and np.isclose(x[1], 1.0)     # K = 1/3 (C3 along z)
    assert np.isclose(o[10], 0, atol=1e-9) and np.isclose(o[0], 1.234)


def test_sigma_and_phase_relations():
    """Flipping the C3 axis (sigma -> -sigma) and rotating by alpha about z: q30 changes sign,
    q33 rotates by 3 alpha (+ pi for the flip), q43 = sigma-consistent: e^{i psi} = sigma e^{i phi}."""
    L = 4; P = Params.from_z(zG=0.08, zH4=0.8, K=2.0)
    S = oriented_tetrahedron([1, 1, 1], 0.0)
    flip = Rot.from_rotvec([np.pi, 0, 0])
    vals = []
    for sig in (1, -1):
        for al in (0.0, 0.3, 1.1):
            B = Rot.from_rotvec([0, 0, al]).apply(S if sig == 1 else flip.apply(S))
            o, _ = _measure(np.tile(B, (L * L, 1, 1)), P, L)
            vals.append((sig, al, o[1] / L ** 2, complex(o[6], o[7]) / L ** 2, complex(o[8], o[9]) / L ** 2))
    q33_0 = vals[0][3]; q43_0 = vals[0][4]
    for sig, al, m, q33, q43 in vals:
        assert np.isclose(np.sign(m), sig)
        assert np.isclose(abs(q33), 1) and np.isclose(abs(q43), 1)
        # psi = 3 alpha up to a sigma-independent constant
        if sig == 1:
            assert np.isclose(q43 / q43_0, np.exp(3j * al)) and np.isclose(q33 / q33_0, np.exp(3j * al))
    # e^{i psi} / (sigma e^{i phi}) is the same constant for both sigma
    r = [v[4] / (v[0] * v[3]) for v in vals]
    assert np.allclose(r, r[0])


def test_stiffness_matches_finite_twist():
    """S2/9 - with S1 = 0 at a uniform state - equals d^2 E/d psi^2 per bond from a finite twist."""
    L = 4; P = Params.from_z(zG=0.08, zH4=0.8, K=2.0, a4=-0.02); lat = Lattice(2, L)
    rng = np.random.default_rng(1)
    S = np.array([Rot.from_rotvec(0.2 * rng.normal(size=3)).apply(oriented_tetrahedron([1, 1, 1], 0.0))
                  for _ in range(L * L)])
    o, _ = _measure(S, P, L)
    h = 1e-4; z = np.array([0.0, 0.0, 1.0]); p = P.array()
    e = [twisted_inter_energy(S, lat.nbr, p, z, d, 0) for d in (-h, 0.0, h)]
    d2 = (e[0] - 2 * e[1] + e[2]) / h ** 2; d1 = (e[2] - e[0]) / (2 * h)
    assert np.isclose(o[2] * 9, d2, rtol=1e-5) and np.isclose(o[3] * 3, d1, rtol=1e-5, atol=1e-9)
