"""Independent reference for SPEC_multi.md §5 (spectral image stacks). numpy/scipy only.
    python3 -I test/ref/ref_image.py   (from modules/xps)  ->  test/ref/ref_image.json
"""
import json, os
import numpy as np
from scipy.optimize import lsq_linear
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(65536)
L = lambda a: np.asarray(a, float).ravel().tolist()

nx, ny, nE = 20, 16, 60
E = 280 + 0.15 * np.arange(nE)
def g(c, w): return np.exp(-4 * np.log(2) * (E - c) ** 2 / w ** 2)
r1, r2 = g(284.8, 1.2), g(286.6, 1.4)
X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
a = 300 * (X + 1) / nx
b = 250 * (((X - 9) ** 2 + (Y - 7) ** 2) < 20) + 30
bg = 40 + 2 * Y
mu = a.ravel()[:, None] * r1[None] + b.ravel()[:, None] * r2[None] + bg.ravel()[:, None]     # (np, nE) pixel-major
cf = 2.0
D = rng.poisson(mu * cf) / cf
D = D.astype(np.float32).astype(float)          # the engine stores Float32
npix = nx * ny
out = dict(nx=nx, ny=ny, E=L(E), data=L(D), cfac=cf, r1=L(r1), r2=L(r2), a=L(a), b=L(b))

# window image
lo, hi = 284.0, 286.0
m = (E >= lo) & (E <= hi)
out["window"] = dict(lo=lo, hi=hi, sum=L(D[:, m].sum(1)))

# area maps on [283, 289.5]
m = (E >= 283) & (E <= 289.5); Ew = E[m]; N = m.sum(); n = 3
w = np.zeros(N); h = np.diff(Ew); w[:-1] += h / 2; w[1:] += h / 2
def linbg(y):
    ylo, yhi = y[:n].mean(), y[-n:].mean(); Elo, Ehi = Ew[:n].mean(), Ew[-n:].mean()
    return ylo + (yhi - ylo) * (Ew - Elo) / (Ehi - Elo)
def shirley(y):
    ylo, yhi = y[:n].mean(), y[-n:].mean(); B = np.full(N, ylo)
    for it in range(200):
        f = y - B; Q = np.r_[0, np.cumsum((f[1:] + f[:-1]) * h / 2)]
        nb = ylo + (yhi - ylo) * Q / Q[-1]; mx = np.abs(nb - B).max(); B = nb
        if mx < 1e-7 * (abs(yhi - ylo) + 1e-12): break
    return B
coef = np.array([((np.eye(N)[i] - linbg(np.eye(N)[i])) * w).sum() for i in range(N)])
Am, sAm, Sm = [], [], []
for p in range(npix):
    y = D[p, m]
    Am.append(((y - linbg(y)) * w).sum())
    Sm.append(((y - shirley(y)) * w).sum())
    sAm.append(np.sqrt((coef ** 2 * np.maximum(y, 1 / cf) / cf).sum()))
out["area"] = dict(lo=283, hi=289.5, n=n, linear=L(Am), shirley=L(Sm), sigma=L(sAm))

# second area map (window over peak 2 only) and at% maps
m2 = (E >= 285.8) & (E <= 288.8)
A2 = np.array([np.trapezoid(D[p, m2] - (D[p, m2][:n].mean() + (D[p, m2][-n:].mean() - D[p, m2][:n].mean()) * (E[m2] - E[m2][:n].mean()) / (E[m2][-n:].mean() - E[m2][:n].mean())), E[m2]) for p in range(npix)])
S1, S2 = 1.0, 1.7
n1, n2 = np.array(Am) / S1, A2 / S2
s1 = np.array(sAm) / S1
at1 = np.full(npix, np.nan); at2 = np.full(npix, np.nan); sat1 = np.full(npix, np.nan)
for p in range(npix):
    use1 = Am[p] > 0 and not (sAm[p] > 0 and Am[p] < 2 * sAm[p]); use2 = A2[p] > 0
    tot = (n1[p] if use1 else 0) + (n2[p] if use2 else 0)
    if not tot > 0: continue
    at1[p] = 100 * n1[p] / tot if use1 else 0; at2[p] = 100 * n2[p] / tot if use2 else 0
    if use1:
        g11 = 100 * (tot - n1[p]) / tot ** 2
        sat1[p] = abs(g11) * s1[p]          # area 2 given without σ
out["atomic"] = dict(lo2=285.8, hi2=288.8, A2=L(A2), S=[S1, S2], at1=L(at1), at2=L(at2), sat1=L(sat1))

# ROI spectrum (rect x 3..7, y 2..5)
mask = np.zeros((ny, nx), bool); mask[2:6, 3:8] = True; mk = mask.ravel()
s = D[mk].sum(0)
out["roi"] = dict(rect=dict(x0=3, y0=2, x1=7, y1=5), mean=L(s / mk.sum()), sigma=L(np.sqrt(s / cf) / mk.sum()))

# binning 2x2 (sum)
Db = D.reshape(ny, nx, nE).reshape(ny // 2, 2, nx // 2, 2, nE).sum((1, 3)).reshape(-1, nE)
out["bin2"] = L(Db)

# PCA (energies = variables, pixels = observations)
def prep(scale, center):
    Xm = D.T.copy()       # nE × npix
    rowW = np.ones(nE); colW = np.ones(npix)
    if scale == "poisson":
        gm = Xm.mean(1); hm = Xm.mean(0)
        rowW = 1 / np.sqrt(np.maximum(gm, gm.max() * 1e-6)); colW = 1 / np.sqrt(np.maximum(hm, hm.max() * 1e-6))
    Xm = Xm * rowW[:, None] * colW[None, :]
    mean = Xm.mean(1) if center == "mean" else np.zeros(nE)
    return Xm - mean[:, None], rowW, colW, mean
pcs = []
for scale in ["none", "poisson"]:
    for center in ["none", "mean"]:
        Xm, rowW, colW, mean = prep(scale, center)
        U, sv, Vt = np.linalg.svd(Xm, full_matrices=False)
        for k in range(U.shape[1]):
            i = np.argmax(np.abs(U[:, k]))
            if U[i, k] < 0: U[:, k] *= -1; Vt[k] *= -1
        lam = sv ** 2
        R_, C_ = max(nE, npix - (center == "mean")), min(nE, npix - (center == "mean"))
        IND = [np.sqrt(lam[k:C_].sum() / (R_ * (C_ - k))) / (C_ - k) ** 2 for k in range(1, C_)]
        rec3 = ((U[:, :3] * sv[:3]) @ Vt[:3] + mean[:, None]) / rowW[:, None] / colW[None, :]
        pcs.append(dict(scale=scale, center=center, eig=L(lam[:10]), U3=L(U[:, :3].T), scores3=L(Vt[:3]), IND=L(IND[:10]), nIND=int(np.argmin(IND) + 1), rec3=L(rec3.T)))
out["pca"] = pcs

# noise table on the full stack with the poisson-scaled PCA reconstruction (method "pca")
Xm, rowW, colW, mean = prep("poisson", "none")
U, sv, Vt = np.linalg.svd(Xm, full_matrices=False)
tab = []
for k in range(0, 5):
    rec = ((U[:, :k] * sv[:k]) @ Vt[:k] + mean[:, None]) / rowW[:, None] / colW[None, :]
    v = np.maximum(rec if k > 0 else D.T, 1 / cf) / cf
    chi = (((D.T - rec) ** 2) / v).sum()
    tab.append(dict(n=k, chi2=chi, dof=(nE - k) * (npix - k)))
out["noise"] = tab

# per-pixel NNLS with constant background, unweighted and Poisson IRLS
A = np.stack([r1, r2, np.ones(nE)], 1)
res = {}
for mode in ["none", "poisson"]:
    C = np.zeros((npix, 3)); Sg = np.full((npix, 3), np.nan); RC = np.zeros(npix)
    for p in range(npix):
        y = D[p]
        wv = cf / np.maximum(y, 1 / cf) if mode == "poisson" else np.ones(nE)
        for it in range(2 if mode == "poisson" else 1):
            sw = np.sqrt(wv)
            x = lsq_linear(A * sw[:, None], y * sw, bounds=([0, 0, -np.inf], [np.inf] * 3), method="bvls", tol=1e-14).x
            if mode == "poisson" and it == 0: wv = cf / np.maximum(A @ x, 1 / cf)
        act = [k for k in range(3) if x[k] != 0 or k == 2]
        Ap = A[:, act]; cov = np.linalg.inv(Ap.T @ (wv[:, None] * Ap))
        for i, k in enumerate(act): Sg[p, k] = np.sqrt(cov[i, i])
        C[p] = x; RC[p] = (wv * (y - A @ x) ** 2).sum() / (nE - len(act))
    res[mode] = dict(coef=L(C[:, :2].T), bg=L(C[:, 2]), sig=[None if np.isnan(v) else v for v in Sg[:, :2].T.ravel()], redchi=L(RC))
out["nnls"] = res
def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o): return None
    return o
json.dump(clean(out), open(os.path.join(HERE, "ref_image.json"), "w"), allow_nan=False)
print("wrote ref_image.json")
