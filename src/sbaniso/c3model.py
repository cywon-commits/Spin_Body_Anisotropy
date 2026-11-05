"""C3-axis effective model of tetrahedral spin bodies = general XY-Ising model.

    H = - sum_<ab> [ J1 s_a s_b + (J2 + J4 s_a s_b) cos(phi_a - phi_b) ]          (phi frame)
      = - sum_<ab> [ J1 s_a s_b + (J4 + J2 s_a s_b) cos(psi_a - psi_b) ]          (psi frame)

    psi = phi - pi (1 - s) / 2,   i.e.  exp(i psi) = s exp(i phi).

Variables: s = +-1 (sign of the body octupole q30, odd under time reversal),
phi = 3*alpha + pi(1-s)/2 (full octupole phase), psi = 3*alpha (hexadecapole q43 phase, TR even).
Granato-Kosterlitz-Lee-Nightingale notation:  A = J4, B = J2, C = J1.

Rigid tetrahedron (n = 4) with inter-body couplings G x^3 + H4 x^4 per spin pair:
    J1 = 160/81 G,  J2 = 128/81 G,  J4 = 512/729 H4,   r = J4/J2 = (4/9) zH4/zG.

Couplings are passed to every kernel as c = np.array([J1, J2, J4]).  Square lattice, periodic.
Random numbers: every kernel takes st = uint64[4], the xoshiro256** state of the replica it updates
(sbaniso.rng), so replicas can be updated in parallel and runs are reproducible from one seed.

Updates
    sweep            Metropolis: s flip / (s flip, phi + pi) / phi rotation
    wolff_sigma_phi  s cluster at fixed phi, bond J1 + J4 cos(dphi)   (frustrated when J4 > J1)
    wolff_sigma_psi  s cluster at fixed psi (flip s and phi -> phi + pi), bond J1 + J2 cos(dpsi)
                     -> unfrustrated whenever J1 > J2.  REQUIRED near the crossing (tau ~ 1 sweep
                     instead of 200-1500 at L = 64, see docs/pilot_report.md)
    wolff_phi        Wolff embedding for phi at fixed s, bond J2 + J4 s s'
All cluster moves accept any sign of the bond (a bond joins when K s s' > 0, resp. K u u' > 0).
"""
import numpy as np
from numba import njit

from .rng import randint, uniform

J1_PER_G = 160.0 / 81.0
J2_PER_G = 128.0 / 81.0
J4_PER_H4 = 512.0 / 729.0


def couplings_from_r(r, G=1.0):
    """Rigid-tetrahedron couplings at ratio r = J4/J2 (units of G)."""
    return np.array([J1_PER_G * G, J2_PER_G * G, r * J2_PER_G * G])


def couplings_from_body(zG, zH4, z=4):
    G, H4 = zG / z, zH4 / z
    return np.array([J1_PER_G * G, J2_PER_G * G, J4_PER_H4 * H4])


@njit(cache=True)
def local_e(s, p, i, j, L, si, pi_, c):
    e = 0.0
    for d in range(4):
        a = (i + (d == 0) - (d == 1)) % L
        b = (j + (d == 2) - (d == 3)) % L
        ss = si * s[a, b]
        e -= c[0] * ss + (c[1] + c[2] * ss) * np.cos(pi_ - p[a, b])
    return e


@njit(cache=True)
def energy(s, p, L, c):
    e = 0.0
    for i in range(L):
        for j in range(L):
            for (a, b) in (((i + 1) % L, j), (i, (j + 1) % L)):
                ss = s[i, j] * s[a, b]
                e -= c[0] * ss + (c[1] + c[2] * ss) * np.cos(p[i, j] - p[a, b])
    return e


@njit(cache=True)
def sweep(s, p, L, beta, c, delta, st):
    acc = 0
    for _ in range(L * L):
        i = randint(st, L)
        j = randint(st, L)
        u = uniform(st)
        s_new = s[i, j]
        p_new = p[i, j]
        if u < 0.3:
            s_new = -s_new
        elif u < 0.4:
            s_new = -s_new
            p_new = p_new + np.pi
        else:
            p_new = p_new + delta * (2 * uniform(st) - 1)
        dE = local_e(s, p, i, j, L, s_new, p_new, c) - local_e(s, p, i, j, L, s[i, j], p[i, j], c)
        if dE <= 0 or uniform(st) < np.exp(-beta * dE):
            s[i, j] = s_new
            p[i, j] = p_new
            acc += 1
    return acc / (L * L)


@njit(cache=True)
def wolff_sigma_phi(s, p, L, beta, c, stack, inc, mem, st):
    """s cluster at fixed phi: Ising bonds K = J1 + J4 cos(dphi)."""
    i0 = randint(st, L); j0 = randint(st, L)
    inc[i0, j0] = True
    stack[0, 0] = i0; stack[0, 1] = j0
    n = 1; size = 0
    while n > 0:
        n -= 1
        i = stack[n, 0]; j = stack[n, 1]
        mem[size, 0] = i; mem[size, 1] = j
        size += 1
        for d in range(4):
            a = (i + (d == 0) - (d == 1)) % L
            b = (j + (d == 2) - (d == 3)) % L
            if inc[a, b]:
                continue
            x = (c[0] + c[2] * np.cos(p[i, j] - p[a, b])) * s[i, j] * s[a, b]
            if x > 0 and uniform(st) < 1 - np.exp(-2 * beta * x):
                inc[a, b] = True
                stack[n, 0] = a; stack[n, 1] = b; n += 1
    for k in range(size):
        i = mem[k, 0]; j = mem[k, 1]
        inc[i, j] = False
        s[i, j] = -s[i, j]
    return size


@njit(cache=True)
def wolff_sigma_psi(s, p, L, beta, c, stack, inc, mem, st):
    """s cluster at fixed psi: flip s and phi -> phi + pi.  Ising bonds K = J1 + J2 cos(dpsi)."""
    i0 = randint(st, L); j0 = randint(st, L)
    inc[i0, j0] = True
    stack[0, 0] = i0; stack[0, 1] = j0
    n = 1; size = 0
    while n > 0:
        n -= 1
        i = stack[n, 0]; j = stack[n, 1]
        mem[size, 0] = i; mem[size, 1] = j
        size += 1
        for d in range(4):
            a = (i + (d == 0) - (d == 1)) % L
            b = (j + (d == 2) - (d == 3)) % L
            if inc[a, b]:
                continue
            ss = s[i, j] * s[a, b]
            cos_dpsi = ss * np.cos(p[i, j] - p[a, b])
            x = (c[0] + c[1] * cos_dpsi) * ss
            if x > 0 and uniform(st) < 1 - np.exp(-2 * beta * x):
                inc[a, b] = True
                stack[n, 0] = a; stack[n, 1] = b; n += 1
    for k in range(size):
        i = mem[k, 0]; j = mem[k, 1]
        inc[i, j] = False
        s[i, j] = -s[i, j]
        p[i, j] = p[i, j] + np.pi
    return size


@njit(cache=True)
def wolff_phi(s, p, L, beta, c, stack, inc, mem, st):
    """Wolff embedding (reflection about a random axis) for phi at fixed s; bonds J2 + J4 s s'."""
    ax = 2 * np.pi * uniform(st)
    i0 = randint(st, L); j0 = randint(st, L)
    inc[i0, j0] = True
    stack[0, 0] = i0; stack[0, 1] = j0
    n = 1; size = 0
    while n > 0:
        n -= 1
        i = stack[n, 0]; j = stack[n, 1]
        mem[size, 0] = i; mem[size, 1] = j
        size += 1
        ui = np.cos(p[i, j] - ax)
        for d in range(4):
            a = (i + (d == 0) - (d == 1)) % L
            b = (j + (d == 2) - (d == 3)) % L
            if inc[a, b]:
                continue
            x = (c[1] + c[2] * s[i, j] * s[a, b]) * ui * np.cos(p[a, b] - ax)
            if x > 0 and uniform(st) < 1 - np.exp(-2 * beta * x):
                inc[a, b] = True
                stack[n, 0] = a; stack[n, 1] = b; n += 1
    for k in range(size):
        i = mem[k, 0]; j = mem[k, 1]
        inc[i, j] = False
        p[i, j] = 2 * ax + np.pi - p[i, j]
    return size


@njit(cache=True)
def full_step(s, p, L, beta, c, delta, stack, inc, mem, ncl, sizes, st):
    """One Metropolis sweep, then a FIXED number of clusters of each type:
    ncl = (n_sigma_phi, n_sigma_psi, n_phi).  sizes[k] accumulates the flipped sites (for tuning).

    The number of clusters must not depend on the cluster sizes drawn in this step: the common rule
    'build clusters until ~N sites were flipped' is a state-dependent stopping time and biases the
    sampled distribution (6x6 Ising at T = 2.5: <m^2> = 0.610 instead of the exact 0.5286;
    see tests/test_c3model.py::test_exact_ising_6x6).  Tune ncl during equilibration only."""
    acc = sweep(s, p, L, beta, c, delta, st)
    for k in range(ncl[0]):
        sizes[0] += wolff_sigma_phi(s, p, L, beta, c, stack, inc, mem, st)
    for k in range(ncl[1]):
        sizes[1] += wolff_sigma_psi(s, p, L, beta, c, stack, inc, mem, st)
    for k in range(ncl[2]):
        sizes[2] += wolff_phi(s, p, L, beta, c, stack, inc, mem, st)
    return acc


@njit(cache=True)
def measure(s, p, L, c, out):
    """out[0]=E, [1]=sum s, [2..5]=cx,sx,cy,sy (helicity parts for a twist of phi, identical to a
    twist of psi), [6],[7]=Re,Im sum exp(i phi), [8],[9]=Re,Im sum exp(i psi) = s exp(i phi),
    [10],[11],[12] = |F(k_min)|^2 of s, exp(i phi), exp(i psi), averaged over k = (2pi/L, 0) and
    (0, 2pi/L)  (second-moment correlation length: xi = sqrt(<|F(0)|^2>/<|F(k)|^2> - 1)/(2 sin(pi/L)))."""
    E = 0.0; m = 0.0
    cx = 0.0; sx = 0.0; cy = 0.0; sy = 0.0
    fr = 0.0; fi = 0.0; sr = 0.0; si = 0.0
    # Fourier sums along x (index i) and y (index j): complex amplitudes of each field
    ax_s = 0.0 + 0.0j; ay_s = 0.0 + 0.0j
    ax_f1 = 0.0 + 0.0j; ay_f1 = 0.0 + 0.0j      # exp(i phi) at +k
    ax_f2 = 0.0 + 0.0j; ay_f2 = 0.0 + 0.0j      # exp(i psi) at +k
    w = 2.0 * np.pi / L
    for i in range(L):
        ei = np.cos(w * i) + 1j * np.sin(w * i)
        for j in range(L):
            ej = np.cos(w * j) + 1j * np.sin(w * j)
            a = (i + 1) % L
            b = (j + 1) % L
            ssx = s[i, j] * s[a, j]
            ssy = s[i, j] * s[i, b]
            kx = c[1] + c[2] * ssx
            ky = c[1] + c[2] * ssy
            dx = p[i, j] - p[a, j]
            dy = p[i, j] - p[i, b]
            E -= c[0] * (ssx + ssy) + kx * np.cos(dx) + ky * np.cos(dy)
            cx += kx * np.cos(dx); sx += kx * np.sin(dx)
            cy += ky * np.cos(dy); sy += ky * np.sin(dy)
            m += s[i, j]
            co = np.cos(p[i, j]); sn = np.sin(p[i, j])
            fr += co; fi += sn
            sr += s[i, j] * co; si += s[i, j] * sn
            z1 = co + 1j * sn
            z2 = s[i, j] * z1
            ax_s += s[i, j] * ei; ay_s += s[i, j] * ej
            ax_f1 += z1 * ei; ay_f1 += z1 * ej
            ax_f2 += z2 * ei; ay_f2 += z2 * ej
    out[0] = E; out[1] = m
    out[2] = cx; out[3] = sx; out[4] = cy; out[5] = sy
    out[6] = fr; out[7] = fi; out[8] = sr; out[9] = si
    out[10] = 0.5 * (abs(ax_s) ** 2 + abs(ay_s) ** 2)
    out[11] = 0.5 * (abs(ax_f1) ** 2 + abs(ay_f1) ** 2)
    out[12] = 0.5 * (abs(ax_f2) ** 2 + abs(ay_f2) ** 2)


OBS_NAMES = ("E", "m_sigma", "cx", "sx", "cy", "sy", "Mphi_re", "Mphi_im", "Mpsi_re", "Mpsi_im",
             "Fk_sigma", "Fk_phi", "Fk_psi")
NOBS = len(OBS_NAMES)


def new_state(L, rng=None, ordered=False):
    rng = np.random.default_rng() if rng is None else rng
    if ordered:
        return np.ones((L, L), dtype=np.int64), np.zeros((L, L))
    s = rng.choice(np.array([-1, 1]), size=(L, L)).astype(np.int64)
    p = 2 * np.pi * rng.random((L, L))
    return s, p


def work_arrays(L):
    N = L * L
    return (np.zeros((N + 4, 2), dtype=np.int64), np.zeros((L, L), dtype=np.bool_),
            np.zeros((N + 4, 2), dtype=np.int64))
