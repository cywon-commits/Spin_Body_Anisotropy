"""Analysis of PT runs (HDF5 from sbaniso.run; pilot npz files also accepted).

Statistics
    tau_int          integrated autocorrelation time, Sokal automatic window (c = 6)
    jackknife        block jackknife for any function of block means
Reweighting
    MultiHistogram   Ferrenberg-Swendsen / MBAR free energies over all PT temperatures, then
                     reweighted moments at any beta inside the ladder.  Errors: block jackknife
                     with the free energies held at their full-data values (standard approximation;
                     the same time block is deleted at all temperatures, since PT couples them).
Observables (per temperature or reweighted; N = L^2)
    U4     = 1 - <m^4> / (3 <m^2>^2)                         (Ising Binder, m = sum sigma)
    Ups    = [<cx + cy> - beta <sx^2 + sy^2>] / (2N)        (helicity modulus, phi twist = psi twist)
    C      = beta^2 (<E^2> - <E>^2) / N
    absm   = <|m|> / N
    xiL_X  = sqrt(<|F_X(0)|^2> / <|F_X(k_min)|^2> - 1) / (2 L sin(pi/L)),  X = sigma, phi, psi
Fits and crossings
    crossing         where f_L1(T) - f_L2(T) changes sign (Binder), or f_L(T) = g(T) (Ups = 2T/pi)
    weber_minnhagen  T_KT minimising chi^2/dof of Ups(L) pi/(2T) = 1 + 1/(2 ln L + c) over all L
    extrapolate_crossings  T(L1, L2) = T_inf + a Lm^{-x}, x free or fixed (x = 1/nu + omega)
Diagnostics
    quarter_check    quarter-split averages per temperature (CLAUDE.md, MC rule 3)
    bimodality       energy-histogram dip ratio (1 = unimodal)
"""
from dataclasses import dataclass, field

import h5py
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

NB = 20

# ----------------------------------------------------------------------------- loading


@dataclass
class Run:
    L: int
    Ts: np.ndarray
    couplings: np.ndarray
    obs: np.ndarray                   # (R, n, NOBS) in c3model.OBS_NAMES order (13 columns)
    diagnostics: dict = field(default_factory=dict)
    attrs: dict = field(default_factory=dict)

    @property
    def betas(self):
        return 1.0 / self.Ts

    @property
    def N(self):
        return self.L * self.L


def load_run(path):
    path = str(path)
    if path.endswith(".npz"):  # pilot v3 format: 10 observables, no Fourier columns
        d = np.load(path)
        o = d["out"].astype(np.float64)
        obs = np.concatenate([o, np.full(o.shape[:2] + (3,), np.nan)], axis=2)
        diag = dict(swap_per_sweep=d["swap"], roundtrips_halfway=int(d["roundtrips"]), ncl=d["ncl"])
        return Run(int(d["L"]), d["Ts"].astype(float), d["c"].astype(float), obs, diag,
                   dict(git=str(d["git"]), seed=int(d["seed"])))
    with h5py.File(path, "r") as f:
        diag = {k: f["diagnostics"][k][()] for k in f["diagnostics"]} if "diagnostics" in f else {}
        return Run(int(f["params/L"][()]), f["params/Ts"][()], f["params/couplings"][()],
                   f["timeseries/obs"][()], diag, dict(f.attrs))


# ----------------------------------------------------------------------------- statistics


def tau_int(x, c=6.0):
    """Integrated autocorrelation time (tau = 1/2 + sum rho(t)), Sokal's automatic window
    W = min{W : W >= c tau(W)}.  Returns (tau, error estimate)."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 10 or np.var(x) == 0:
        return 0.5, 0.0
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n] / (np.arange(n, 0, -1))
    rho = acf / acf[0]
    tau = 0.5
    for W in range(1, n):
        tau += rho[W]
        if W >= c * tau:
            break
    tau = max(tau, 0.5)
    return tau, tau * np.sqrt((4 * W + 2) / n)


def block_sums(x, nb=NB, axis=-1):
    """Split along `axis` into nb equal blocks (dropping the remainder) and sum each block."""
    x = np.moveaxis(np.asarray(x), axis, -1)
    n = x.shape[-1] // nb * nb
    return x[..., :n].reshape(x.shape[:-1] + (nb, n // nb)).sum(-1)


def jackknife(pseudo, central=None):
    """pseudo: (nb, ...) leave-one-block-out estimates.  Returns (bias-corrected?) central, error.
    The central value returned is the full-data value if given (no bias correction)."""
    pseudo = np.asarray(pseudo)
    nb = pseudo.shape[0]
    mean = pseudo.mean(0)
    err = np.sqrt((nb - 1) / nb * ((pseudo - mean) ** 2).sum(0))
    return (mean if central is None else central), err


# ----------------------------------------------------------------------------- moments


MOMENT_NAMES = ("one", "E", "E2", "m2", "m4", "absm", "cxy", "sxy2", "Mphi2", "Fk_phi", "Mpsi2",
                "Fk_psi", "Fk_sigma")


def sample_moments(obs):
    """Per-sample quantities whose averages give every observable.  obs: (..., 13)."""
    E = obs[..., 0]; m = obs[..., 1]
    return np.stack([np.ones_like(E), E, E * E, m * m, m ** 4, np.abs(m),
                     obs[..., 2] + obs[..., 4], obs[..., 3] ** 2 + obs[..., 5] ** 2,
                     obs[..., 6] ** 2 + obs[..., 7] ** 2, obs[..., 11],
                     obs[..., 8] ** 2 + obs[..., 9] ** 2, obs[..., 12],
                     obs[..., 10]], axis=-1)


def derived(mom, beta, L):
    """Observables from averaged moments mom[..., len(MOMENT_NAMES)] (already normalised)."""
    N = L * L
    g = {k: mom[..., i] for i, k in enumerate(MOMENT_NAMES)}
    s = 2 * L * np.sin(np.pi / L)

    def xiL(F0, Fk):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.sqrt(np.maximum(F0 / Fk - 1, 0)) / s

    return dict(
        E=g["E"] / N,
        C=beta ** 2 * (g["E2"] - g["E"] ** 2) / N,
        U4=1 - g["m4"] / (3 * g["m2"] ** 2),
        absm=g["absm"] / N,
        m2=g["m2"] / N ** 2,
        Ups=(g["cxy"] - beta * g["sxy2"]) / (2 * N),
        Mphi2=g["Mphi2"] / N ** 2,
        Mpsi2=g["Mpsi2"] / N ** 2,
        xiL_sigma=xiL(g["m2"], g["Fk_sigma"]),
        xiL_phi=xiL(g["Mphi2"], g["Fk_phi"]),
        xiL_psi=xiL(g["Mpsi2"], g["Fk_psi"]),
    )


def per_temperature(run, nb=NB):
    """Observables at the simulated temperatures with block-jackknife errors.
    Returns {name: (value[R], error[R])}."""
    mom = sample_moments(run.obs)                     # (R, n, M)
    bs = block_sums(mom, nb, axis=1)                  # (R, M, nb)
    tot = bs.sum(-1)
    full = derived(tot / tot[:, :1], run.betas, run.L)
    jk = (tot[..., None] - bs)                        # (R, M, nb)
    jk = np.moveaxis(jk / jk[:, :1], -1, 0)           # (nb, R, M)
    pseudo = derived(jk, run.betas[None, :], run.L)
    return {k: jackknife(pseudo[k], full[k]) for k in full}


# ----------------------------------------------------------------------------- reweighting


class MultiHistogram:
    """Multiple-histogram reweighting (Ferrenberg-Swendsen 1989; solved as the MBAR convex
    problem, Shirts & Chodera 2008).  E: (R, n) total energies at inverse temperatures betas."""

    def __init__(self, E, betas, chunk=200000):
        self.E = np.asarray(E, float)
        self.betas = np.asarray(betas, float)
        R, n = self.E.shape
        self.lnn = np.full(R, np.log(n))
        self.Ef = self.E.ravel()
        self.chunk = chunk
        self.f = self._solve()
        self.logD = self._logD(self.f)                 # per sample: ln sum_j n_j e^{f_j - b_j E}

    def _logD(self, f):
        out = np.empty_like(self.Ef)
        for a in range(0, len(self.Ef), self.chunk):
            e = self.Ef[a:a + self.chunk]
            out[a:a + self.chunk] = logsumexp(self.lnn[:, None] + f[:, None] - np.outer(self.betas, e), axis=0)
        return out

    def _objective(self, g):
        f = np.concatenate([[0.0], g])
        val = 0.0; grad = np.zeros_like(f)
        for a in range(0, len(self.Ef), self.chunk):
            e = self.Ef[a:a + self.chunk]
            A = self.lnn[:, None] + f[:, None] - np.outer(self.betas, e)
            lse = logsumexp(A, axis=0)
            val += lse.sum()
            grad += np.exp(A - lse).sum(1)
        n = np.exp(self.lnn)
        val -= (n * f).sum(); grad -= n
        return val / n.sum(), grad[1:] / n.sum()

    def _solve(self):
        # f_k = -ln Z_k (f_0 = 0).  Initial guess from the two one-sided neighbour estimates
        # f_{k+1} - f_k = -ln <e^{-db E}>_k = +ln <e^{+db E}>_{k+1},  db = beta_{k+1} - beta_k
        R = len(self.betas)
        f0 = np.zeros(R)
        for k in range(R - 1):
            db = self.betas[k + 1] - self.betas[k]
            fwd = -(logsumexp(-db * self.E[k]) - np.log(self.E.shape[1]))
            bwd = logsumexp(db * self.E[k + 1]) - np.log(self.E.shape[1])
            f0[k + 1] = f0[k] + 0.5 * (fwd + bwd)
        res = minimize(self._objective, f0[1:], jac=True, method="L-BFGS-B",
                       options=dict(maxiter=5000, gtol=1e-12, ftol=1e-15))
        return np.concatenate([[0.0], res.x])

    def log_weights(self, beta):
        """ln of the (unnormalised) weight of every sample at inverse temperature beta."""
        return -beta * self.Ef - self.logD

    def moments(self, X, betas_out, nb=NB, max_elems=2.5e7):
        """Reweighted averages of the per-sample arrays X (R, n, M) at each beta in betas_out.
        Returns full (B, M) and jackknife pseudo-values (nb, B, M) (normalised by weight sum).
        Block b = the same time window at every temperature (samples beyond nb*(n//nb) dropped)."""
        R, n, M = X.shape
        bs = n // nb
        E = self.E[:, :nb * bs].reshape(R, nb, bs)
        lD = self.logD.reshape(R, n)[:, :nb * bs].reshape(R, nb, bs)
        Xb = X[:, :nb * bs].reshape(R, nb, bs, M)
        bo = np.asarray(betas_out, float)
        B = len(bo)
        # common shift per beta: max log-weight over all samples
        shift = np.full(B, -np.inf)
        for r in range(R):
            lw = -bo[:, None] * E[r].reshape(-1)[None, :] - lD[r].reshape(-1)[None, :]
            shift = np.maximum(shift, lw.max(1))
        Sb = np.zeros((nb, B, M)); Wb = np.zeros((nb, B))
        step = max(1, int(max_elems // max(1, B * bs)))
        for blk in range(nb):
            for r0 in range(0, R, step):
                e = E[r0:r0 + step, blk].reshape(-1); d = lD[r0:r0 + step, blk].reshape(-1)
                w = np.exp(-bo[:, None] * e[None, :] - d[None, :] - shift[:, None])   # (B, s)
                Sb[blk] += w @ Xb[r0:r0 + step, blk].reshape(-1, M)
                Wb[blk] += w.sum(1)
        St = Sb.sum(0); Wt = Wb.sum(0)
        full = St / Wt[:, None]
        pseudo = (St[None] - Sb) / (Wt[None] - Wb)[..., None]
        return full, pseudo


def curves(run, Tgrid, nb=NB, mh=None):
    """Reweighted observables on Tgrid with jackknife pseudo-values.
    Returns (values {name: (B,)}, errors {name: (B,)}, pseudo {name: (nb, B)}, mh)."""
    mh = mh or MultiHistogram(run.obs[..., 0], run.betas)
    X = sample_moments(run.obs)
    bo = 1.0 / np.asarray(Tgrid)
    full, pseudo = mh.moments(X, bo, nb)
    v = derived(full, bo, run.L)
    p = derived(pseudo, bo[None, :], run.L)
    e = {k: jackknife(p[k])[1] for k in v}
    return v, e, p, mh


# ----------------------------------------------------------------------------- crossings and fits


def crossing(T, f, which="last"):
    """T where f changes sign (linear interpolation); nan if none."""
    f = np.asarray(f)
    idx = [i for i in range(len(f) - 1) if np.isfinite(f[i]) and np.isfinite(f[i + 1])
           and np.sign(f[i]) != np.sign(f[i + 1]) and f[i] != 0]
    if not idx:
        return np.nan
    i = idx[-1] if which == "last" else idx[0]
    return T[i] + f[i] / (f[i] - f[i + 1]) * (T[i + 1] - T[i])


def _interp_at(T, y, t):
    return np.interp(t, T, y) if np.isfinite(t) else np.nan


def binder_crossing(T, cur1, cur2, window=(0.25, 0.66), key="U4"):
    """Crossing of `key` (default U4) for two sizes, restricted to values of the smaller L inside
    `window` (for U4 this excludes the trivial low-T crossing at U4 -> 2/3).
    cur = (values, errors, pseudo) from curves().  Returns dict(T, dT, U, dU) (U = value of `key`
    at the crossing) with errors from the two independent jackknifes in quadrature.
    Also used for xi/L crossings (xi_crossing)."""
    def one(u1, u2):
        ok = (u1 > window[0]) & (u1 < window[1])
        t = crossing(T, np.where(ok, u1 - u2, np.nan))
        return t, _interp_at(T, u1, t)
    t0, u0 = one(cur1[0][key], cur2[0][key])
    var_t = 0.0; var_u = 0.0
    for reps in ([one(pu, cur2[0][key]) for pu in cur1[2][key]],       # vary L1, hold L2
                 [one(cur1[0][key], pu) for pu in cur2[2][key]]):      # vary L2, hold L1
        reps = np.array(reps)
        nb = len(reps)
        var_t += (nb - 1) / nb * np.nansum((reps[:, 0] - np.nanmean(reps[:, 0])) ** 2)
        var_u += (nb - 1) / nb * np.nansum((reps[:, 1] - np.nanmean(reps[:, 1])) ** 2)
    return dict(T=t0, dT=np.sqrt(var_t), U=u0, dU=np.sqrt(var_u))


XI_L_KT_XY = 0.7507   # xi_2nd/L of the XY model at T_KT, L -> infinity (Hasenbusch 2005)


def xi_crossing(T, cur1, cur2, key="xiL_psi", window=(0.3, 1.5)):
    """Crossing of xi/L (auxiliary KT indicator, docs/N1_analysis_plan.md): T and xi/L at the
    crossing of two sizes, to be compared with XI_L_KT_XY.  Not a decision rule."""
    r = binder_crossing(T, cur1, cur2, window=window, key=key)
    return dict(T=r["T"], dT=r["dT"], xiL=r["U"], dxiL=r["dU"])


def ups_crossing(T, cur):
    """T where Ups(T) = 2T/pi (the highest such T), with jackknife error."""
    t0 = crossing(T, cur[0]["Ups"] - 2 * T / np.pi)
    reps = np.array([crossing(T, p - 2 * T / np.pi) for p in cur[2]["Ups"]])
    return dict(T=t0, dT=jackknife(reps)[1])


def _wm_scan(T, Ls, ups, errs, cgrid):
    """chi^2/dof profile of the Weber-Minnhagen fit; returns (chi2dof[T], best c[T]).
    c is minimised at each T: grid search, then golden-section refinement in the bracketing cell."""
    lnL = np.log(np.asarray(Ls, float))
    y = np.array([u * np.pi / (2 * T) for u in ups])                  # (nL, nT)
    s = np.array([np.maximum(e, 1e-6) * np.pi / (2 * T) for e in errs])

    def chi(cv):                                                      # cv: (nT,) -> (nT,)
        pred = 1 + 1 / (2 * lnL[:, None] + cv[None, :])
        return (((y - pred) / s) ** 2).sum(0)

    pred = 1 + 1 / (2 * lnL[:, None] + cgrid[None, :])                # (nL, nc)
    full = (((y[:, :, None] - pred[:, None, :]) / s[:, :, None]) ** 2).sum(0)   # (nT, nc)
    j = full.argmin(1)
    a = cgrid[np.maximum(j - 1, 0)]; b = cgrid[np.minimum(j + 1, len(cgrid) - 1)]
    g = (np.sqrt(5) - 1) / 2
    x1 = b - g * (b - a); x2 = a + g * (b - a); f1 = chi(x1); f2 = chi(x2)
    for _ in range(60):
        left = f1 < f2
        b = np.where(left, x2, b); a = np.where(left, a, x1)
        x2n = np.where(left, x1, a + g * (b - a)); x1n = np.where(left, b - g * (b - a), x2)
        x1, x2 = x1n, x2n; f1 = chi(x1); f2 = chi(x2)
    cbest = 0.5 * (a + b)
    dof = max(1, len(Ls) - 1)
    return chi(cbest) / dof, cbest


def weber_minnhagen(T, curs, Ls, cgrid=None):
    """T_KT = argmin_T chi^2/dof of Ups(L, T) pi/(2T) = 1 + 1/(2 ln L + c) (c free at each T),
    all L at once.  curs: list of curves() outputs for Ls (common T grid).  Error: jackknife,
    replacing one size's curve at a time by its pseudo-values, variances added over sizes.
    Returns dict(T, dT, chi2dof, c, profile)."""
    Ls = list(Ls)
    if cgrid is None:
        cgrid = np.linspace(-2 * np.log(min(Ls)) + 0.2, 60, 2000)
    ups = [c[0]["Ups"] for c in curs]; errs = [c[1]["Ups"] for c in curs]
    prof, cbest = _wm_scan(T, Ls, ups, errs, cgrid)
    i = int(np.nanargmin(prof))
    Tbest = _parabolic_min(T, prof)
    var = 0.0
    for a, c in enumerate(curs):
        reps = []
        for pu in c[2]["Ups"]:
            u2 = list(ups); u2[a] = pu
            p2, _ = _wm_scan(T, Ls, u2, errs, cgrid)
            reps.append(_parabolic_min(T, p2))
        reps = np.array(reps); nb = len(reps)
        var += (nb - 1) / nb * ((reps - reps.mean()) ** 2).sum()
    return dict(T=float(Tbest), dT=float(np.sqrt(var)), chi2dof=float(prof[i]), c=float(cbest[i]),
                profile=prof)


def _parabolic_min(T, prof):
    """Minimum of a sampled profile, refined by a parabola through the lowest point and its
    neighbours (grid point itself at the edges)."""
    i = int(np.nanargmin(prof))
    if i == 0 or i == len(T) - 1:
        return float(T[i])
    y0, y1, y2 = prof[i - 1], prof[i], prof[i + 1]
    den = y0 - 2 * y1 + y2
    if den <= 0:
        return float(T[i])
    return float(T[i] + 0.5 * (y0 - y2) / den * (T[i + 1] - T[i]))


def extrapolate_crossings(Lm, Tc, dTc, x=None):
    """Fit T(Lm) = T_inf + a Lm^{-x} (x fixed if given, else free).  Lm = mean size of each pair.
    Returns dict(T_inf, dT_inf, a, x, dx, chi2dof)."""
    from scipy.optimize import curve_fit
    Lm = np.asarray(Lm, float); Tc = np.asarray(Tc, float); dTc = np.asarray(dTc, float)
    if x is not None:
        A = np.vstack([np.ones_like(Lm), Lm ** -x]).T / dTc[:, None]
        sol, *_ = np.linalg.lstsq(A, Tc / dTc, rcond=None)
        cov = np.linalg.inv(A.T @ A)
        chi = (((A @ sol) - Tc / dTc) ** 2).sum()
        return dict(T_inf=sol[0], dT_inf=np.sqrt(cov[0, 0]), a=sol[1], x=x, dx=0.0,
                    chi2dof=chi / max(1, len(Lm) - 2))
    f = lambda L, T0, a, xx: T0 + a * L ** -xx
    p, cov = curve_fit(f, Lm, Tc, p0=[Tc[-1], 0.0, 2.0], sigma=dTc, absolute_sigma=True, maxfev=20000)
    chi = (((f(Lm, *p) - Tc) / dTc) ** 2).sum()
    return dict(T_inf=p[0], dT_inf=np.sqrt(cov[0, 0]), a=p[1], x=p[2], dx=np.sqrt(cov[2, 2]),
                chi2dof=chi / max(1, len(Lm) - 3))


# ----------------------------------------------------------------------------- diagnostics


def quarter_check(run, names=("E", "absm", "Ups", "U4"), nb=NB):
    """Quarter-split test (CLAUDE.md, MC rule 3).  For every temperature and observable, each
    quarter of the series is compared with the mean of the other three quarters.  The error of a
    quarter is taken from the full-run block jackknife (sigma_q = 2 sigma_full; sigma of the other
    three = (2/sqrt(3)) sigma_full), which is far less noisy than errors from a few blocks inside
    each quarter (those follow a heavy-tailed t distribution and gave spurious |z| > 5).
    Returns (z_max per temperature, details incl. the Bonferroni threshold for p = 0.01)."""
    R, n, _ = run.obs.shape
    q = n // 4
    full = per_temperature(run, nb)
    parts = [per_temperature(Run(run.L, run.Ts, run.couplings, run.obs[:, i * q:(i + 1) * q]), nb=5)
             for i in range(4)]
    zs = np.zeros(R); detail = {}
    for name in names:
        v = np.array([x[name][0] for x in parts])                   # (4, R)
        sig = np.maximum(full[name][1], 1e-300)
        z = np.array([(v[i] - np.delete(v, i, 0).mean(0)) / (sig * np.sqrt(4 + 4 / 3))
                      for i in range(4)])
        detail[name] = np.abs(z).max(0)
        zs = np.maximum(zs, detail[name])
    from scipy.stats import norm
    ncomp = R * len(names) * 4
    detail["threshold"] = float(norm.isf(0.005 / ncomp))            # two-sided p = 0.01 overall
    detail["ncomp"] = ncomp
    return zs, detail


def tau_table(run, names=("E", "absm", "m2", "Mpsi2")):
    """tau_int at every temperature for a few observables (sweeps)."""
    obs = run.obs; N = run.N
    series = dict(E=obs[..., 0], absm=np.abs(obs[..., 1]), m2=obs[..., 1] ** 2,
                  Mpsi2=obs[..., 8] ** 2 + obs[..., 9] ** 2, Mphi2=obs[..., 6] ** 2 + obs[..., 7] ** 2,
                  cxy=obs[..., 2] + obs[..., 4])
    return {k: np.array([tau_int(series[k][r])[0] for r in range(len(run.Ts))]) for k in names}


def equilibration_check(run, nb=NB, min_blocks_per_tau=10, min_roundtrips=20, min_swap=0.2):
    """Combined equilibration criterion (CLAUDE.md MC rule 3; docs/N0_report.md, appendix).
    A run passes only if ALL hold:
      quarter      quarter_check z_max below its Bonferroni threshold
      tau          every jackknife block spans >= min_blocks_per_tau tau_int, i.e.
                   tau_int <= n / (nb * min_blocks_per_tau) for E, |m|, m^2, |M_psi|^2 at every T.
                   Without this the block errors absorb a slow drift and the quarter test loses
                   its power (negative control: tau = 893 sweeps, quarter z = 4.22 = threshold)
      roundtrips   full PT round trips >= min_roundtrips
      swap         per-attempt swap acceptance >= min_swap for every neighbour pair
    Returns dict(pass, reasons, quarter_z, quarter_threshold, tau_max, tau_limit, roundtrips,
    swap_min)."""
    z, d = quarter_check(run, nb=nb)
    tau = tau_table(run)
    tau_max = max(float(v.max()) for v in tau.values())
    tau_limit = run.obs.shape[1] / (nb * min_blocks_per_tau)
    diag = run.diagnostics
    rt = int(diag.get("roundtrips", -1)); sw = float(np.min(diag["swap_acc"])) if "swap_acc" in diag else np.nan
    checks = dict(quarter=bool(z.max() < d["threshold"]), tau=bool(tau_max <= tau_limit),
                  roundtrips=bool(rt >= min_roundtrips), swap=bool(np.isfinite(sw) and sw >= min_swap))
    return dict(**{"pass": all(checks.values())}, checks=checks,
                reasons=[k for k, v in checks.items() if not v], quarter_z=float(z.max()),
                quarter_threshold=d["threshold"], tau_max=tau_max, tau_limit=tau_limit,
                roundtrips=rt, swap_min=sw)


def bimodality(E, bins=40):
    """Ratio of the histogram dip between the main peak and the second-highest local maximum to
    that maximum (1 = unimodal, < ~0.8 suggests two peaks)."""
    h, _ = np.histogram(E, bins=bins)
    h = np.convolve(h, np.ones(3) / 3, mode="same")
    i1 = int(np.argmax(h)); best = 1.0
    for j in range(len(h)):
        if abs(j - i1) < 4:
            continue
        lo, hi = sorted((i1, j))
        if h[j] > 0.2 * h[i1] and h[j] == h[max(0, j - 2):j + 3].max():
            best = min(best, h[lo:hi + 1].min() / h[j])
    return best
