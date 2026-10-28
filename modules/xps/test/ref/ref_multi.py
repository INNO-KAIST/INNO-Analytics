"""Independent reference for SPEC_multi.md §1-3 (series processing, PCA/TFA, reference NNLS).
Written from the SPEC equations with numpy/scipy only (no translation of the JS).
    python3 -I test/ref/ref_multi.py      (from modules/xps)  ->  test/ref/ref_multi.json
"""
import json, os
import numpy as np
from scipy import stats
from scipy.optimize import nnls, lsq_linear

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261028)
L = lambda a: np.asarray(a, float).tolist()


def sign_fix(U, V):
    """element of largest |value| of each left singular vector positive (SPEC_multi §2.1)"""
    U = U.copy(); V = V.copy()
    for k in range(U.shape[1]):
        i = np.argmax(np.abs(U[:, k]))
        if U[i, k] < 0:
            U[:, k] *= -1; V[:, k] *= -1
    return U, V


def gl(E, c, w, m=0.3):
    g = np.exp(-4 * np.log(2) * (E - c) ** 2 / w ** 2)
    l = 1 / (1 + 4 * (E - c) ** 2 / w ** 2)
    return (1 - m) * g + m * l


def mixture(nE=160, nL=24, counts=True):
    """3 chemical states (Si 2p-like) with known depth profiles + constant background, Poisson noise"""
    E = np.linspace(96, 108, nE)
    S = np.stack([gl(E, 99.4, 1.0), gl(E, 101.5, 1.3), gl(E, 103.6, 1.5)], 1)
    t = np.linspace(0, 1, nL)
    C = np.stack([1 / (1 + np.exp(-(t - 0.45) / 0.08)), np.exp(-((t - 0.4) / 0.15) ** 2), 1 / (1 + np.exp((t - 0.35) / 0.08))], 1)
    mu = 3000 * S @ C.T + 150 * (1 + 3 * t ** 2)[None, :]     # level-dependent constant background (an independent factor)
    D = rng.poisson(mu).astype(float)
    return E, S, C, mu, D


def preprocess(D, scale, center):
    X = D.copy()
    rowW = colW = mean = None
    if scale == "area":
        colW = 1 / X.sum(0); X = X * colW
    elif scale == "poisson":
        g = X.mean(1); h = X.mean(0)
        g = np.maximum(g, g.max() * 1e-6); h = np.maximum(h, h.max() * 1e-6)
        rowW = 1 / np.sqrt(g); colW = 1 / np.sqrt(h)
        X = X * rowW[:, None] * colW[None, :]
    if center == "mean":
        mean = X.mean(1); X = X - mean[:, None]
    return X, rowW, colW, mean


def malinowski(lam, r, c, alpha=0.05):
    R, C = max(r, c), min(r, c)
    rows = []
    for n in range(1, C):
        tail = lam[n:C].sum()
        RE = np.sqrt(tail / (R * (C - n)))
        IND = RE / (C - n) ** 2
        j = np.arange(n + 1, C + 1)
        F = ((R - j + 1) * (C - j + 1)).sum() / ((R - n + 1) * (C - n + 1)) * lam[n - 1] / tail
        SL = stats.f.sf(F, 1, C - n)
        rows.append(dict(n=n, RE=RE, IND=IND, F=F, SL=SL))
    nIND = min(rows, key=lambda o: o["IND"])["n"]
    nF = 0
    for o in rows:
        if o["SL"] < alpha: nF = o["n"]
        else: break
    return rows, nIND, nF


out = {}

# ---- 1. SVD vs numpy (various shapes) ----
svd_cases = []
for (r, c, rank) in [(120, 10, 10), (8, 40, 8), (300, 25, 3), (60, 60, 60)]:
    A = rng.normal(size=(r, rank)) @ rng.normal(size=(rank, c))
    if rank < min(r, c):
        A += 1e-3 * rng.normal(size=(r, c))
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    U, V = sign_fix(U, Vt.T)
    svd_cases.append(dict(cols=L(A.T), s=L(s), U3=L(U[:, :3].T), V3=L(V[:, :3].T)))
out["svd"] = svd_cases

# ---- 2. PCA on a Poisson mixture, all preprocessing options ----
E, S, C, mu, D = mixture()
out["mix"] = dict(E=L(E), cols=L(D.T), S=L(S.T), C=L(C.T), mu=L(mu.T))
pcs = []
for scale in ["none", "poisson", "area"]:
    for center in ["none", "mean"]:
        X, rowW, colW, mean = preprocess(D, scale, center)
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        U, V = sign_fix(U, Vt.T)
        lam = s ** 2
        rows, nIND, nF = malinowski(lam, X.shape[0], X.shape[1] - (1 if center == "mean" else 0))
        # reconstruction with 3 factors in original units
        Xr = (U[:, :3] * s[:3]) @ V[:, :3].T
        if mean is not None: Xr = Xr + mean[:, None]
        if rowW is not None: Xr = Xr / rowW[:, None]
        if colW is not None: Xr = Xr / colW[None, :]
        pcs.append(dict(scale=scale, center=center, s=L(s), U3=L(U[:, :3].T), mal=rows, nIND=nIND, nF=nF, rec3=L(Xr.T)))
out["pca"] = pcs

# ---- 2.3 noise-limited factor count: weighted (Pearson/Poisson) low-rank fit, alternating weighted LS (MLPCA) ----
def wlra(D, n, cf=1.0, tol=1e-10, maxit=300):
    r, c = D.shape
    if n == 0:
        return float((D ** 2 / (np.maximum(D, 1 / cf) / cf)).sum()), 0
    X, rowW, colW, _ = preprocess(D, "poisson", "none")
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    R0 = (U[:, :n] * s[:n]) @ Vt[:n] / rowW[:, None] / colW[None, :]
    U0, s0, V0t = np.linalg.svd(R0, full_matrices=False)
    U0, V0 = sign_fix(U0, V0t.T)
    A = U0[:, :n] * s0[:n]; B = V0[:, :n].copy()
    rec = R0; chi_old = np.inf
    for it in range(maxit):
        Wt = cf / np.maximum(rec, 1 / cf)
        for j in range(c):
            Aw = A * Wt[:, j:j + 1]
            B[j] = np.linalg.solve(A.T @ Aw, Aw.T @ D[:, j])
        for i in range(r):
            Bw = B * Wt[i][:, None]
            A[i] = np.linalg.solve(B.T @ Bw, Bw.T @ D[i])
        rec = A @ B.T
        chi = float((cf / np.maximum(rec, 1 / cf) * (D - rec) ** 2).sum())
        if abs(chi_old - chi) <= tol * chi: break
        chi_old = chi
    return chi, it + 1

tab = []; nNoise = None
r, c = D.shape
for n in range(0, 6):
    chi, it = wlra(D, n)
    dof = (r - n) * (c - n); lim = dof + 2 * np.sqrt(2 * dof)
    tab.append(dict(n=n, chi2=chi, dof=dof, within=bool(chi <= lim), iter=it))
    if nNoise is None and chi <= lim: nNoise = n
out["noise"] = dict(table=tab, nNoise=nNoise)
U, s, Vt = np.linalg.svd(D, full_matrices=False)

# ---- 2.4 target testing (Malinowski) on the raw (non-centred, unscaled) mixture, n = 4 (3 states + background) ----
n = 4
lam = s ** 2
R_, C_ = max(r, c), min(r, c)
RE = np.sqrt(lam[n:C_].sum() / (R_ * (C_ - n)))
targets = {"state1": S[:, 0], "state2": S[:, 1], "state3": S[:, 2], "flat": np.ones(r),
           "wrong_shifted": gl(E, 100.4, 1.0), "wrong_narrow": gl(E, 101.5, 0.5)}
tt = []
Un, sn = U[:, :n], s[:n]
for name, x in targets.items():
    t = (Un.T @ x) / sn
    xh = Un @ (Un.T @ x)
    AET = np.sqrt(((x - xh) ** 2).sum() / r)
    REP = RE * np.linalg.norm(t)
    RET = np.sqrt(max(AET ** 2 - REP ** 2, 0))
    tt.append(dict(name=name, x=L(x), AET=AET, REP=REP, RET=RET, SPOIL=RET / REP))
out["target"] = dict(n=n, RE=RE, tests=tt)

# ---- 2.5 OPA + ITTFA (deterministic algorithms) ----
Y = D / np.linalg.norm(D, axis=0)
m = Y.mean(1); m /= np.linalg.norm(m)
picks = []
for k in range(3):
    refs = [m] if k == 0 else [Y[:, j] for j in picks]
    best, bj = -1, -1
    for j in range(c):
        if j in picks: continue
        M = np.stack(refs + [Y[:, j]], 1)
        d = np.linalg.det(M.T @ M)
        if d > best: best, bj = d, j
    picks.append(bj)
V3 = Vt[:3, :].T
prof = []
for k in range(3):
    x = np.zeros(c); x[picks[k]] = 1
    for it in range(2000):
        y = V3 @ (V3.T @ x); y[y < 0] = 0; y /= y.max()
        dd = np.abs(y - x).max(); x = y
        if dd < 1e-10: break
    prof.append(x)
out["opa"] = picks

# SIMPLISMA purest energy channels (Windig & Guilment 1991), offset f = 3 %
mu_ = D.mean(1); sd_ = D.std(1); al = 0.03 * mu_.max()
lam_ = np.sqrt(sd_ ** 2 + (mu_ + al) ** 2)
Yn = D / lam_[:, None] / np.sqrt(c)
sp = []
for k in range(4):
    best, bi = -np.inf, -1
    for i in range(r):
        if i in sp: continue
        w = 1.0 if k == 0 else np.linalg.det(np.stack([Yn[q] for q in sp] + [Yn[i]]) @ np.stack([Yn[q] for q in sp] + [Yn[i]]).T)
        pv = w * sd_[i] / (mu_[i] + al)
        if pv > best: best, bi = pv, i
    sp.append(bi)
out["simplisma"] = sp
out["ittfa"] = L(np.array(prof))

# ---- 3. NNLS vs scipy ----
nn = []
for (m_, n_) in [(40, 4), (200, 6), (30, 10)]:
    A = rng.normal(size=(m_, n_)); b = rng.normal(size=m_) + A @ (np.abs(rng.normal(size=n_)) * (rng.random(n_) > 0.4))
    x, _ = nnls(A, b)
    nn.append(dict(cols=L(A.T), b=L(b), x=L(x)))
out["nnls"] = nn

# ---- 3.2 reference fit (weighted, IRLS 2 passes, free background terms) ----
cf = 1.0
refs = [S[:, 0], S[:, 1], S[:, 2]]
fits = []
Em, hw = (E[0] + E[-1]) / 2, (E[-1] - E[0]) / 2
for j in [0, 7, 12, 23]:
    y = D[:, j]
    for nb in [0, 1, 2]:
        cols = refs + ([np.ones(r)] if nb >= 1 else []) + ([(E - Em) / hw] if nb >= 2 else [])
        A = np.stack(cols, 1)
        lb = np.array([0.0] * 3 + [-np.inf] * nb); ub = np.full(A.shape[1], np.inf)
        w = cf / np.maximum(y, 1 / cf)
        for p in range(2):
            sw = np.sqrt(w)
            res = lsq_linear(A * sw[:, None], y * sw, bounds=(lb, ub), method="bvls", tol=1e-14)
            x = res.x
            muf = A @ x
            if p == 0: w = cf / np.maximum(muf, 1 / cf)
        # w is now the pass-2 weight (from pass-1 μ); χ² and covariance use it
        act = [k for k in range(A.shape[1]) if x[k] != 0 or k >= 3]
        Ap = A[:, act]
        cov = np.linalg.inv(Ap.T @ (w[:, None] * Ap))
        sig = np.full(A.shape[1], np.nan); sig[act] = np.sqrt(np.diag(cov))
        chi = float((w * (y - A @ x) ** 2).sum())
        fits.append(dict(level=j, bg=nb, coef=L(x[:3]), bgCoef=L(x[3:]), sig=[None if np.isnan(v) else v for v in sig[:3]], chi2=chi))
out["reffit"] = fits

# ---- 1.2 region integration (linear + Shirley) with Poisson σ ----
def endpoints(E, y, n):
    N = len(E); n = max(1, min(n, N // 4 or 1))
    return y[:n].mean(), y[-n:].mean(), E[:n].mean(), E[-n:].mean(), n

def trapz(E, f):
    return float(((f[1:] + f[:-1]) * np.diff(E)).sum() / 2)

def shirley(E, y, n):
    ylo, yhi, _, _, _ = endpoints(E, y, n)
    B = np.full_like(y, ylo)
    for it in range(200):
        f = y - B
        Q = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) * np.diff(E) / 2)])
        nb = ylo + (yhi - ylo) * Q / Q[-1]
        mx = np.abs(nb - B).max(); B = nb
        if mx < 1e-7 * (abs(yhi - ylo) + 1e-12): break
    return B

def linear_bg(E, y, n):
    ylo, yhi, Elo, Ehi, _ = endpoints(E, y, n)
    return ylo + (yhi - ylo) * (E - Elo) / (Ehi - Elo)

def area_coefs(E, n):
    """∂A/∂y for A = trapz(y − B_linear)"""
    N = len(E); n = max(1, min(n, N // 4 or 1))
    a = np.zeros(N)
    for i in range(N):
        e = np.zeros(N); e[i] = 1
        a[i] = trapz(E, e - linear_bg(E, e, n))
    return a

integ = []
win = (E >= 97) & (E <= 106)
Ew = E[win]
a = area_coefs(Ew, 5)
for kind in ["linear", "shirley"]:
    rows = []
    for j in range(c):
        y = D[win, j]
        B = linear_bg(Ew, y, 5) if kind == "linear" else shirley(Ew, y, 5)
        A_ = trapz(Ew, y - B)
        sA = np.sqrt((a ** 2 * np.maximum(y, 1)).sum())
        rows.append([A_, sA])
    integ.append(dict(bg=kind, lo=97, hi=106, rows=rows))
out["integrate"] = integ

# ---- 1.3 composition per level (SPEC §4: n = A/(S·T·λ) with T=λ=1, n.d. if A≤0 or A<2σ; delta-method σ) ----
def quant(A, sA, Sf):
    n = A / Sf; sn = sA / Sf
    use = (A > 0) & ~((sA > 0) & (A < 2 * sA))
    tot = n[use].sum()
    at = np.where(use, 100 * n / tot, np.nan)
    sat = np.full(len(A), np.nan)
    for i in range(len(A)):
        if not use[i]: continue
        g = np.array([100 * ((tot if i == j else 0) - n[i]) / tot ** 2 for j in range(len(A))])
        sat[i] = np.sqrt((g[use] ** 2 * sn[use] ** 2).sum())
    return at, sat

areas3 = np.stack([3000 * C[:, k] * trapz(E, S[:, k]) + rng.normal(size=c) * 30 for k in range(3)], 0)
sig3 = np.full_like(areas3, 30.0)
Sf = np.array([0.82, 0.82, 1.1])
comp = []
for j in range(c):
    at, sat = quant(areas3[:, j], sig3[:, j], Sf)
    comp.append(dict(at=[None if np.isnan(v) else v for v in at], sAt=[None if np.isnan(v) else v for v in sat]))
out["composition"] = dict(areas=L(areas3), sig=L(sig3), S=L(Sf), levels=comp)

# ---- 1.4 information-depth correction: c = I − λ dI/dz, dI/dz by local quadratic LS (hw points each side) ----
z = np.cumsum(np.r_[0, rng.uniform(0.3, 0.7, 29)])
truth = 1 / (1 + np.exp(-(z - 6) / 0.4))
lamE = 1.6
# measured I(z0) = ∫ c(z) e^{-(z-z0)/λ}/λ dz  (fine numeric integral)
zz = np.linspace(0, 40, 400001)
cz = 1 / (1 + np.exp(-(zz - 6) / 0.4))
I = np.array([np.trapezoid(cz[zz >= z0] * np.exp(-(zz[zz >= z0] - z0) / lamE) / lamE, zz[zz >= z0]) for z0 in z])
def deriv(z, y, hw):
    n = len(z); d = np.zeros(n)
    for i in range(n):
        j0, j1 = max(0, i - hw), min(n - 1, i + hw)
        while j1 - j0 < 2:
            if j0 > 0: j0 -= 1
            elif j1 < n - 1: j1 += 1
            else: break
        p = np.polyfit(z[j0:j1 + 1] - z[i], y[j0:j1 + 1], 2)
        d[i] = p[1]
    return d
X2 = np.stack([I, 1 - I + 0.0], 0)  # second element = complementary signal, λ different
lam2 = np.array([lamE, 2.4])
corr = np.maximum(0, X2 - lam2[:, None] * np.stack([deriv(z, X2[0], 2), deriv(z, X2[1], 2)]))
out["infodepth"] = dict(z=L(z), X=L(X2), lam=L(lam2), hw=2, n=L(corr), truthA=L(truth))

json.dump(out, open(os.path.join(HERE, "ref_multi.json"), "w"))
print("wrote ref_multi.json", {k: (len(v) if hasattr(v, "__len__") else v) for k, v in out.items()})
