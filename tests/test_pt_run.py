"""Infrastructure tests: RNG, Fourier observables, PT engine determinism, YAML -> HDF5 runs,
checkpoint/resume."""
import h5py
import numba
import numpy as np
import pytest

from sbaniso.c3model import NOBS, OBS_NAMES, couplings_from_r, measure, new_state
from sbaniso.pt import PTState, advance, geometric_ladder
from sbaniso.rng import make_streams, next_u64, uniform
from sbaniso.run import run_config


def _xoshiro_ref(s, n):
    M = (1 << 64) - 1
    rot = lambda x, k: ((x << k) | (x >> (64 - k))) & M
    s = [int(x) for x in s]; out = []
    for _ in range(n):
        out.append(rot((s[1] * 5) & M, 7) * 9 & M)
        t = (s[1] << 17) & M
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rot(s[3], 45)
    return out


def test_xoshiro_matches_reference_and_streams_differ():
    st = make_streams(123, 3)
    ref = _xoshiro_ref(st[0], 10)
    assert [int(next_u64(st[0])) for _ in range(10)] == ref
    assert len({tuple(r) for r in make_streams(5, 50)}) == 50
    u = np.array([uniform(st[1]) for _ in range(100000)])
    assert abs(u.mean() - 0.5) < 0.005 and u.min() >= 0 and u.max() < 1
    v = np.array([uniform(st[2]) for _ in range(100000)])
    assert abs(np.corrcoef(u, v)[0, 1]) < 0.015


def test_fourier_observables_match_numpy():
    L = 8; c = couplings_from_r(4.1)
    s, p = new_state(L, np.random.default_rng(1)); o = np.zeros(NOBS)
    measure(s, p, L, c, o)
    fields = [s.astype(float), np.exp(1j * p), s * np.exp(1j * p)]
    for q, f in zip((10, 11, 12), fields):
        F = np.fft.fft2(f)
        # measure() uses exp(+i k x); fft uses exp(-i k x) -> compare with the -k component
        ref = 0.5 * (abs(F[-1, 0]) ** 2 + abs(F[0, -1]) ** 2)
        assert np.isclose(o[q], ref), (OBS_NAMES[q], o[q], ref)
    assert np.isclose(o[1], s.sum()) and np.isclose(o[8] + 1j * o[9], fields[2].sum())


def _run_engine(threads, seed=3, L=12, n=300, ntherm=100):
    Ts = geometric_ladder(6.5, 7.8, 6); c = couplings_from_r(4.2); flags = (1, 1, 1)
    old = numba.get_num_threads(); numba.set_num_threads(threads)
    try:
        st = PTState(L, Ts, seed, flags)
        rows = advance(st, L, 1 / Ts, c, flags, ntherm, n)
    finally:
        numba.set_num_threads(old)
    return st, rows


def test_engine_independent_of_thread_count():
    st1, r1 = _run_engine(1)
    st4, r4 = _run_engine(min(4, numba.config.NUMBA_NUM_THREADS))
    assert np.array_equal(r1, r4)
    assert np.array_equal(st1.perm, st4.perm) and np.array_equal(st1.rng, st4.rng)


def test_engine_chunks_equal_single_call():
    _, r_all = _run_engine(2, n=300)
    Ts = geometric_ladder(6.5, 7.8, 6); c = couplings_from_r(4.2); flags = (1, 1, 1)
    st = PTState(12, Ts, 3, flags)
    parts = [advance(st, 12, 1 / Ts, c, flags, 100, n) for n in (37, 63, 50, 150)]
    assert np.array_equal(np.concatenate(parts, axis=1), r_all)


def test_engine_measures_correct_energy_and_ncl_frozen():
    Ts = geometric_ladder(6.5, 7.8, 5); c = couplings_from_r(4.2); flags = (1, 1, 1); L = 10
    st = PTState(L, Ts, 9, flags)
    advance(st, L, 1 / Ts, c, flags, 200, 100)          # first half of equilibration: tuning
    ncl = st.ncl.copy()
    rows = advance(st, L, 1 / Ts, c, flags, 200, 300)   # rest: ncl must not change
    assert np.array_equal(st.ncl, ncl)
    # the last row was measured before the last swap step: every row must be the measurement of
    # one replica's final configuration, and every replica must appear once
    meas = np.zeros((len(Ts), NOBS))
    for rep in range(len(Ts)):
        measure(st.S[rep], st.P[rep], L, c, meas[rep])
    match = [[np.allclose(rows[k, -1], meas[rep]) for rep in range(len(Ts))] for k in range(len(Ts))]
    assert all(sum(m) == 1 for m in match) and np.array(match).sum(0).tolist() == [1] * len(Ts)


def _cfg(tmp_path, **kw):
    cfg = dict(name="t", model=dict(r=4.1), L=10, temperatures=dict(Tmin=6.5, Tmax=7.8, n=5),
               ntherm=100, nmeas=400, seed=11, checkpoint_every=70, output=str(tmp_path / "t.h5"))
    cfg.update(kw)
    return cfg


def test_run_writes_complete_hdf5(tmp_path):
    out = run_config(_cfg(tmp_path), verbose=False)
    with h5py.File(out) as f:
        assert f.attrs["status"] == "complete"
        assert f["timeseries/obs"].shape == (5, 400, NOBS)
        assert [x.decode() for x in f["params/obs_names"][()]] == list(OBS_NAMES)
        for k in ("swap_acc", "roundtrips", "ncl", "deltas", "metropolis_acc", "seconds"):
            assert k in f["diagnostics"]
        assert f["final/S"].shape == (5, 10, 10) and int(f["checkpoint/it"][()]) == 500
        assert "seed: 11" in f.attrs["config"] and f.attrs["git"]


def test_resume_is_bit_identical(tmp_path):
    a = run_config(_cfg(tmp_path, output=str(tmp_path / "a.h5")), verbose=False)
    cfg = _cfg(tmp_path, output=str(tmp_path / "b.h5"))
    run_config(cfg, verbose=False, max_sweeps=160)       # stop mid-run (after equilibration)
    with h5py.File(tmp_path / "b.h5", "r+") as f:
        assert f.attrs["status"] == "running"
        n = f["timeseries/obs"].shape[1]
        f["timeseries/obs"].resize(n + 3, axis=1)        # simulate rows written after a crash
    run_config(cfg, verbose=False, max_sweeps=95)
    run_config(cfg, verbose=False)
    with h5py.File(a) as fa, h5py.File(tmp_path / "b.h5") as fb:
        assert np.array_equal(fa["timeseries/obs"][()], fb["timeseries/obs"][()])
        assert np.array_equal(fa["final/S"][()], fb["final/S"][()])


def test_complete_run_is_not_rerun_and_config_change_is_refused(tmp_path):
    out = run_config(_cfg(tmp_path), verbose=False)
    t = out.stat().st_mtime_ns
    run_config(_cfg(tmp_path), verbose=False)
    assert out.stat().st_mtime_ns == t
    with pytest.raises(RuntimeError):
        run_config(_cfg(tmp_path, seed=12), verbose=False)
