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
from sbaniso.fullmodel import NBETA, NEXTRA, _norms, measure_all  # noqa: E402


def _measure(spins_lat, P, L):
    lat = Lattice(2, L)
    out = np.zeros((1, NOBS)); ext = np.zeros((1, NEXTRA))
    measure_all(spins_lat[None], np.array([0]), np.array([1.234]), lat.nbr, P.array(), L,
                np.array(_norms()), out, ext, True, np.zeros((1, NBETA)))
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
    assert np.allclose(x[2:6], [0.0, 1.0, 1.0, 1.0], atol=1e-8) and abs(x[6]) < 1e-8   # no tilt: factors 1


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


# ----------------------------------------------------------------------------- C3-axis tilt (N3e)
from sbaniso.fullmodel import _V, _VM, bond_factors, c3_axis, c3_tilt2  # noqa: E402


def _gs(s):
    from polyspin.groundstate import find_ground_state
    return find_ground_state(Params.from_z(zG=0.08 * s, zH4=0.80 * s, K=2.0, a4=-0.02))[1]


@pytest.mark.parametrize("s", [1.0, 0.03])
def test_tilt_zero_for_statically_deformed_ground_state(s):
    """The T = 0 body is deformed by the anisotropy (K = sum n^4 gives 0.324-0.340, not 1/3);
    the Kabsch-fit tilt must still be zero."""
    S = _gs(s)
    assert c3_tilt2(S, _V, _VM) < 1e-8


@pytest.mark.parametrize("theta", [0.05, 0.2, 0.4, 0.6])
def test_tilt_reproduces_rigid_rotation(theta):
    """Unambiguous up to ~35 deg (0.615 rad): beyond that another C3 axis is closer to z."""
    rng = np.random.default_rng(7)
    S = _gs(1.0)
    for _ in range(5):
        phi = rng.uniform(0, 2 * np.pi)
        R = Rot.from_rotvec(theta * np.array([np.cos(phi), np.sin(phi), 0.0]))
        assert np.isclose(c3_tilt2(R.apply(S), _V, _VM), np.sin(theta) ** 2, atol=1e-8)


def test_tilt_independent_of_labels_and_spin_about_z():
    rng = np.random.default_rng(3)
    S = Rot.from_rotvec([0.3, -0.1, 0.0]).apply(_gs(0.03))
    ref = c3_tilt2(S, _V, _VM)
    for _ in range(10):
        S2 = Rot.from_rotvec([0, 0, rng.uniform(0, 6)]).apply(S)[rng.permutation(4)]
        assert np.isclose(c3_tilt2(np.ascontiguousarray(S2), _V, _VM), ref, atol=1e-10)


def test_tilt_folds_back_beyond_half_axis_angle():
    """C3 axis lines are arccos(1/3) = 70.53 deg apart and the nearest one to z is reported, so theta^2
    is bounded by its value with z along an S4 axis, 1 - 1/3 = 2/3 (large tilts are folded back)."""
    rng = np.random.default_rng(5)
    S = _gs(1.0)
    for _ in range(50):
        R = Rot.random(random_state=rng)
        th2 = c3_tilt2(np.ascontiguousarray(R.apply(S)), _V, _VM)
        assert th2 <= 2 / 3 + 1e-9



# ----------------------------------------------------------------------------- relative tilt (N3f)
def test_bond_factors_small_beta_expansions():
    """1 - d(beta) = c beta^2 + O(beta^4) with c = 3 (J1), 0.75 (J2), 2.75 (J4), i.e. 1/2[l(l+1) - m^2]/2."""
    for beta in (1e-2, 3e-3):
        f = np.array(bond_factors(np.cos(beta)))
        assert np.allclose((1 - f) / beta ** 2, [3.0, 0.75, 2.75], rtol=1e-3)
    assert np.allclose(bond_factors(1.0), (1.0, 1.0, 1.0))


def test_bond_factors_are_wigner_d():
    """Compare with Wigner small-d computed from the explicit sum (independent implementation)."""
    from math import comb, factorial
    def d(l, m, b):
        tot = 0.0
        for k in range(0, 2 * l + 1):
            if l + m - k < 0 or l - m - k < 0 or k - m + m < 0:
                continue
            num = (-1) ** k * np.sqrt(float(factorial(l + m) * factorial(l - m) * factorial(l + m) * factorial(l - m)))
            den = factorial(l + m - k) * factorial(k) * factorial(l - m - k) * factorial(k)
            tot += num / den * np.cos(b / 2) ** (2 * l - 2 * k) * np.sin(b / 2) ** (2 * k)
        return tot
    for b in (0.1, 0.7, 1.9, 2.8):
        assert np.allclose(bond_factors(np.cos(b)), (d(3, 0, b), d(3, 3, b), d(4, 3, b)), atol=1e-12)


def test_common_tilt_gives_no_reduction_and_axis_points_up():
    """Two bodies tilted together (common rotation) have beta = 0: all bond factors 1, theta_a.theta_b =
    theta^2 (C_nn = 1); a flipped body (sigma -> -sigma) still gives an axis with c_z > 0."""
    S = _gs(1.0)
    R = Rot.from_rotvec([0.25, 0.1, 0.0])
    ca = np.zeros(3); cb = np.zeros(3)
    c3_axis(np.ascontiguousarray(R.apply(S)), _V, _VM, ca)
    c3_axis(np.ascontiguousarray(R.apply(Rot.from_rotvec([0, 0, 1.3]).apply(S))), _V, _VM, cb)
    assert ca[2] > 0 and np.allclose(bond_factors(ca @ cb), (1, 1, 1))
    cf = np.zeros(3)
    c3_axis(np.ascontiguousarray(Rot.from_rotvec([np.pi, 0, 0]).apply(S)), _V, _VM, cf)
    assert cf[2] > 0.999
