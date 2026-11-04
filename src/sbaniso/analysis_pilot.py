"""Analysis of the PT pilot runs (C3 effective model).

Per (r, L): Binder U4 of m_sigma, helicity modulus Upsilon (phi twist = psi twist),
specific heat, energy-histogram bimodality at the C peak.  Errors from 20 blocks.
Ising T_I from Binder crossings of size pairs; KT from Upsilon = 2T/pi per L and from a
Weber-Minnhagen fit  Upsilon*pi/(2T) = 1 + 1/(2 ln L + c)  over all L.
"""
import glob, json
import numpy as np

NB = 20


def blocks(x):
    n = len(x) // NB * NB
    return x[:n].reshape(NB, -1)


def observables(d):
    Ts = d["Ts"]; L = int(d["L"]); N = L * L
    out = d["out"].astype(np.float64)
    res = {k: [] for k in ("U4", "dU4", "Ups", "dUps", "C", "absm")}
    for k, T in enumerate(Ts):
        b = 1.0 / T
        E, m, cx, sx, cy, sy = [out[k, :, i] for i in range(6)]

        def u4(m):
            m2 = np.mean(m ** 2); m4 = np.mean(m ** 4)
            return 1 - m4 / (3 * m2 ** 2)

        def ups(cx, sx, cy, sy):
            return ((np.mean(cx) - b * np.mean(sx ** 2)) + (np.mean(cy) - b * np.mean(sy ** 2))) / (2 * N)

        U = u4(m); Ub = np.array([u4(x) for x in blocks(m)])
        Y = ups(cx, sx, cy, sy)
        Yb = np.array([ups(*x) for x in zip(blocks(cx), blocks(sx), blocks(cy), blocks(sy))])
        res["U4"].append(U); res["dU4"].append(Ub.std() / np.sqrt(NB))
        res["Ups"].append(Y); res["dUps"].append(Yb.std() / np.sqrt(NB))
        res["C"].append(b * b * np.var(E) / N); res["absm"].append(np.mean(np.abs(m)))
    return Ts, {k: np.array(v) for k, v in res.items()}


def crossing(T, f):
    """last sign change of f (interp), return T or nan"""
    idx = [i for i in range(len(f) - 1) if np.sign(f[i]) != np.sign(f[i + 1])]
    if not idx:
        return np.nan
    i = idx[-1]
    return T[i] + f[i] / (f[i] - f[i + 1]) * (T[i + 1] - T[i])


def bimodality(E):
    """ratio of histogram dip to smaller peak (1 = unimodal, < ~0.8 suggests two peaks)"""
    h, _ = np.histogram(E, bins=40)
    h = np.convolve(h, np.ones(3) / 3, mode="same")
    i1 = np.argmax(h)
    best = 1.0
    for j in range(len(h)):
        if abs(j - i1) < 4:
            continue
        lo, hi = sorted((i1, j))
        dip = h[lo:hi + 1].min()
        if h[j] > 0.2 * h[i1] and h[j] == h[max(0, j - 2):j + 3].max():
            best = min(best, dip / h[j])
    return best


def main():
    files = sorted(glob.glob("pt_r*_L*.npz"))
    data = {}
    for f in files:
        d = np.load(f)
        data[(float(d["r"]), int(d["L"]))] = d
    summary = {}
    for r in sorted({k[0] for k in data}):
        Ls = sorted(L for (rr, L) in data if rr == r)
        raw = {L: observables(data[(r, L)]) for L in Ls}
        Tg = np.linspace(6.3, 8.0, 341)
        obs = {L: (Tg, {k: np.interp(Tg, raw[L][0], v) for k, v in raw[L][1].items()}) for L in Ls}
        line = {"r": r, "L": Ls}
        # Ising crossings
        TI = {}
        for a, b in zip(Ls, Ls[1:]):
            Ta, oa = obs[a]; Tb, ob = obs[b]
            ok = (oa["U4"] > 0.25) & (oa["U4"] < 0.645)
            f = np.where(ok, oa["U4"] - ob["U4"], np.nan)
            idx = [i for i in range(len(f) - 1) if ok[i] and ok[i + 1] and np.sign(f[i]) != np.sign(f[i + 1])]
            TI[f"{a}-{b}"] = (Ta[idx[0]] + f[idx[0]] / (f[idx[0]] - f[idx[0] + 1]) * (Ta[1] - Ta[0])) if idx else np.nan
            TI[f"U*{a}-{b}"] = float(oa["U4"][idx[0]]) if idx else np.nan
        line["T_I_cross"] = TI
        # KT per L
        line["T_KT_L"] = {L: crossing(obs[L][0], obs[L][1]["Ups"] - 2 * obs[L][0] / np.pi) for L in Ls}
        # Weber-Minnhagen
        T = obs[Ls[0]][0]
        chi = []
        for k, t in enumerate(T):
            ys = np.array([obs[L][1]["Ups"][k] * np.pi / (2 * t) for L in Ls])
            es = np.array([max(obs[L][1]["dUps"][k], 1e-4) * np.pi / (2 * t) for L in Ls])
            best = np.inf
            for c in np.linspace(-5, 30, 701):
                pred = 1 + 1 / (2 * np.log(np.array(Ls)) + c)
                best = min(best, np.sum(((ys - pred) / es) ** 2))
            chi.append(best / max(1, len(Ls) - 1))
        chi = np.array(chi)
        line["T_KT_WM"] = float(T[np.argmin(chi)])
        line["WM_chi2_min"] = float(chi.min())
        # first-order indicators
        line["minU4"] = {L: float(raw[L][1]["U4"].min()) for L in Ls}
        Lmax = Ls[-1]
        k = int(np.argmax(raw[Lmax][1]["C"]))
        line["Cpeak_T_Lmax"] = float(raw[Lmax][0][k])
        line["Cpeak_L"] = {L: float(raw[L][1]["C"].max()) for L in Ls}
        line["Cpeak_T_L"] = {L: float(raw[L][0][np.argmax(raw[L][1]["C"])]) for L in Ls}
        line["E_bimodality_at_Cpeak_Lmax"] = bimodality(data[(r, Lmax)]["out"][k, :, 0].astype(float))
        line["swap_min"] = {L: float(data[(r, L)]["swap"].min()) for L in Ls}
        line["roundtrips"] = {L: int(data[(r, L)]["roundtrips"]) for L in Ls}
        summary[r] = line
        print(json.dumps(line, default=lambda x: round(float(x), 3)))
    return data, summary


if __name__ == "__main__":
    main()
