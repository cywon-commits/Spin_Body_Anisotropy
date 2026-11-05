"""Unit tests of the analysis tools on problems with known answers."""
import numpy as np
import pytest

from sbaniso.analysis import (MultiHistogram, Run, binder_crossing, crossing, curves,
                              extrapolate_crossings, jackknife, per_temperature, quarter_check,
                              tau_int, weber_minnhagen)
from sbaniso.pt import PTState, advance

from test_c3model import _exact_ising


def test_tau_int_ar1():
    rng = np.random.default_rng(0); rho = 0.8; n = 400000
    x = np.empty(n); x[0] = 0
    e = rng.normal(size=n)
    for t in range(1, n):
        x[t] = rho * x[t - 1] + e[t]
    tau, dtau = tau_int(x)
    assert abs(tau - (1 + rho) / (2 * (1 - rho))) < 4 * dtau + 0.05
    assert abs(tau_int(rng.normal(size=100000))[0] - 0.5) < 0.05


def test_jackknife_of_mean_equals_standard_error():
    rng = np.random.default_rng(1); x = rng.normal(size=(20, 500))
    pseudo = np.array([np.delete(x, b, 0).mean() for b in range(20)])
    _, err = jackknife(pseudo)
    assert np.isclose(err, x.mean(1).std(ddof=1) / np.sqrt(20))


def test_crossing_interpolation():
    T = np.linspace(0, 1, 11)
    assert np.isclose(crossing(T, T - 0.437), 0.437)
    assert np.isnan(crossing(T, T + 1))


@pytest.fixture(scope="module")
def ising_pt():
    """PT on the 6x6 Ising model (J2 = J4 = 0), 6 temperatures in [2.1, 3.0]."""
    L = 6; Ts = np.linspace(2.1, 3.0, 6); c = np.array([1.0, 0.0, 0.0]); flags = (1, 1, 0)
    st = PTState(L, Ts, 21, flags)
    rows = advance(st, L, 1 / Ts, c, flags, 2000, 62000)
    return Run(L, Ts, c, rows)


def test_reweighting_reproduces_exact_6x6_ising(ising_pt):
    run = ising_pt
    mh = MultiHistogram(run.obs[..., 0], run.betas)
    for T in (2.27, 2.5, 2.83):                      # between simulated temperatures
        v, e, _, _ = curves(run, np.array([T]), mh=mh)
        Eex, m2ex = _exact_ising(6, T)
        assert abs(v["E"][0] - Eex) < 4.5 * e["E"][0] + 1e-4, (T, v["E"][0], Eex, e["E"][0])
        assert abs(v["m2"][0] - m2ex) < 4.5 * e["m2"][0] + 1e-4, (T, v["m2"][0], m2ex)


def test_reweighting_matches_direct_averages_at_simulated_T(ising_pt):
    v, e, _, _ = curves(ising_pt, ising_pt.Ts)
    d = per_temperature(ising_pt)
    for k in ("E", "U4", "C", "absm"):
        assert np.all(np.abs(v[k] - d[k][0]) < 3 * d[k][1] + 1e-9), k


def test_quarter_check_passes_on_equilibrium_data(ising_pt):
    z, d = quarter_check(ising_pt)
    assert z.max() < d["threshold"]


def test_weber_minnhagen_finds_synthetic_T_kt():
    """Synthetic Ups(L, T): exactly the WM form at T0 = 0.893, an L-dependence that the form cannot
    absorb elsewhere.  The fit must return T0 within the grid spacing."""
    T = np.linspace(0.85, 0.95, 101); Ls = [16, 32, 64, 128]; T0 = 0.893; c0 = 2.0
    rng = np.random.default_rng(3); curs = []
    for L in Ls:
        y = 1 + 1 / (2 * np.log(L) + c0) + 0.4 * (T0 - T) * np.log(L) ** 1.5
        ups = y * 2 * T / np.pi
        err = np.full_like(T, 2e-4)
        pseudo = ups[None, :] + rng.normal(scale=2e-4 / np.sqrt(19), size=(20, len(T)))
        curs.append(({"Ups": ups}, {"Ups": err}, {"Ups": pseudo}))
    r = weber_minnhagen(T, curs, Ls)
    assert abs(r["T"] - T0) <= 0.0011 and r["chi2dof"] < 1e-6
    assert abs(r["c"] - c0) < 0.1


def test_binder_crossing_on_synthetic_scaling():
    """U4(L, T) = f((T - Tc) L^{1/nu}) with f decreasing: all sizes cross exactly at Tc."""
    Tc = 2.269; T = np.linspace(2.2, 2.35, 301)
    f = lambda x: 0.611 - 0.3 * np.tanh(x)
    mk = lambda L: ({"U4": f((T - Tc) * L)}, {"U4": np.full_like(T, 1e-3)},
                    {"U4": f((T - Tc) * L)[None, :] + np.zeros((20, 1))})
    r = binder_crossing(T, mk(32), mk(64))
    assert abs(r["T"] - Tc) < 1e-3 and abs(r["U"] - 0.611) < 1e-3


def test_extrapolation_recovers_parameters():
    Lm = np.array([40, 56, 80, 112.0]); Tinf = 6.9; x = 1.75
    Tc = Tinf + 0.8 * Lm ** -x
    r = extrapolate_crossings(Lm, Tc, np.full(4, 1e-4), x=x)
    assert abs(r["T_inf"] - Tinf) < 1e-6
