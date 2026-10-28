"""Independent reference for SPEC_multi.md §4 (ARXPS). numpy/scipy only, written from the SPEC equations.
    python3 -I test/ref/ref_arxps.py   (from modules/xps)  ->  test/ref/ref_arxps.json
"""
import json, os
import numpy as np
from scipy.optimize import brentq, least_squares, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1976)
L = lambda a: np.asarray(a, float).tolist()
cosd = lambda t: np.cos(np.radians(t))
out = {}

# ---- Seah & Dench (1979) ----
sd = []
for E in [200, 500, 1000, 1400]:
    for a in [0.22, 0.27, 0.31]:
        sd.append(dict(E=E, a=a, element=a * (538 / E ** 2 + 0.41 * np.sqrt(a * E)), inorganic=a * (2170 / E ** 2 + 0.72 * np.sqrt(a * E))))
    sd.append(dict(E=E, rho=1200.0, organic=1000 * (49 / E ** 2 + 0.11 * np.sqrt(E)) / 1200.0))
out["seahDench"] = sd

# ---- single overlayer ----
def R_of(d, th, Lo, Ls, Rinf):
    c = cosd(th)
    return Rinf * (1 - np.exp(-d / (Lo * c))) / np.exp(-d / (Ls * c))

ov = []
for (Lo, Ls, Rinf) in [(3.4, 3.1, 0.9), (2.0, 2.0, 1.3), (1.5, 2.8, 0.5)]:
    for th in [0, 35, 60]:
        for d in [0.3, 1.8, 6.0]:
            R = R_of(d, th, Lo, Ls, Rinf)
            dd = brentq(lambda x: R_of(x, th, Lo, Ls, Rinf) - R, 0, 100, xtol=1e-14)
            ov.append(dict(Lo=Lo, Ls=Ls, Rinf=Rinf, th=th, R=R, d=dd))
out["overlayer"] = ov

# ---- overlayer fit over angles, weighted ----
th = np.array([0, 15, 30, 40, 50, 55, 60, 65.0])
fits = []
for (Lo, Ls, Rinf, d) in [(3.4, 3.1, 0.9, 1.8), (2.6, 2.6, 1.1, 0.7)]:
    Rt = R_of(d, th, Lo, Ls, Rinf)
    sR = 0.02 * Rt + 0.002
    R = Rt + sR * rng.normal(size=len(th))
    for fitR in [False, True]:
        def res(p):
            dd = p[0]; rr = p[1] if fitR else Rinf
            return (R - R_of(dd, th, Lo, Ls, rr)) / sR
        p0 = [1.0, Rinf] if fitR else [1.0]
        sol = least_squares(res, p0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
        # analytic Jacobian of the model for the covariance
        c = cosd(th); dd = sol.x[0]; rr = sol.x[1] if fitR else Rinf
        dR = rr * (np.exp(-dd / (Lo * c)) / (Lo * c) * np.exp(dd / (Ls * c)) + (1 - np.exp(-dd / (Lo * c))) * np.exp(dd / (Ls * c)) / (Ls * c))
        J = [dR / sR] + ([R_of(dd, th, Lo, Ls, rr) / rr / sR] if fitR else [])
        J = np.stack(J, 1)
        C = np.linalg.inv(J.T @ J)
        fits.append(dict(Lo=Lo, Ls=Ls, Rinf=Rinf, fitRinf=fitR, th=L(th), R=L(R), sR=L(sR), d=sol.x[0], sd=np.sqrt(C[0, 0]),
                         Rfit=rr, sRinf=np.sqrt(C[1, 1]) if fitR else None, chi2=float((res(sol.x) ** 2).sum())))
out["fitOverlayer"] = fits

# ---- multilayer: carbon contamination / oxide / substrate ----
model = dict(layers=[dict(name="C", d=0.6), dict(name="ox", d=1.5)],
             species=[dict(name="C1s", layer=0, Iinf=0.8, L=[3.0, 2.8]),
                      dict(name="Si4+", layer=1, Iinf=1.0, L=[3.6, 3.4]),
                      dict(name="Si0", layer=2, Iinf=1.1, L=[3.4, 3.1])])

def ml_I(d, t):
    c = cosd(t); I = []
    for sp in model["species"]:
        m = sp["layer"]; v = sp["Iinf"]
        for j in range(min(m, 2)): v *= np.exp(-d[j] / (sp["L"][j] * c))
        if m < 2: v *= 1 - np.exp(-d[m] / (sp["L"][m] * c))
        I.append(v)
    return np.array(I)

d_true = np.array([0.6, 1.5])
Ith = np.stack([ml_I(d_true, t) for t in th], 1)        # K × A
G = 1000 * (0.5 + 0.5 * cosd(th))
Imeas = Ith * G[None, :]
sI = 0.01 * Imeas
Imeas = Imeas + sI * rng.normal(size=Imeas.shape)

def prof_res(d):
    r = []
    for a, t in enumerate(th):
        M = ml_I(d, t); w = 1 / sI[:, a] ** 2
        g = (w * Imeas[:, a] * M).sum() / (w * M * M).sum()
        r.append((Imeas[:, a] - g * M) / sI[:, a])
    return np.concatenate(r)

sol = least_squares(prof_res, [1.0, 1.0], bounds=([0, 0], [1e3, 1e3]), xtol=1e-15, ftol=1e-15, gtol=1e-15, method="trf")
J = sol.jac; C = np.linalg.inv(J.T @ J)
out["multilayer"] = dict(model=model, th=L(th), I=L(Imeas), sI=L(sI), forward=[L(ml_I(d_true, t)) for t in th],
                         d=L(sol.x), sd=L(np.sqrt(np.diag(C))), chi2=float((sol.fun ** 2).sum()), d_true=L(d_true))

# ---- relative depth plot ----
grz = [len(th) - 1]; nor = [0]
rdp = [float(np.log(Imeas[k, grz].sum() / Imeas[k, nor].sum())) for k in range(3)]
out["rdp"] = dict(rdp=rdp, s=[float(np.sqrt((sI[k, grz] ** 2).sum() / Imeas[k, grz].sum() ** 2 + (sI[k, nor] ** 2).sum() / Imeas[k, nor].sum() ** 2)) for k in range(3)])

# ---- maximum entropy ----
thm = np.array([0, 10, 20, 30, 40, 50, 55, 60, 65, 70.0])
Lm = np.array([2.5, 3.0, 2.8])
J_ = 16
zmax = 3 * Lm.max() * cosd(thm).max()
z = np.arange(J_) * zmax / (J_ - 1)

def kernel(th, Lk, z):
    W = np.zeros((len(Lk), len(th), len(z)))
    for k, l in enumerate(Lk):
        for a, t in enumerate(th):
            lc = l * cosd(t)
            e = np.exp(-z / lc); e2 = np.r_[e[1:], 0]
            W[k, a] = lc * (e - e2)
    return W

W = kernel(thm, Lm, z)
# truth: A-rich top 1 nm, B buried peak at ~2 nm, C substrate (fine-grid truth intensities)
zf = np.linspace(0, 60, 600001)
cA = 0.9 / (1 + np.exp((zf - 1.0) / 0.15))
cB = 0.6 * np.exp(-((zf - 2.2) / 0.6) ** 2)
cC = np.clip(1 - cA - cB, 0, None)
tot = cA + cB + cC
cf = np.stack([cA, cB, cC]) / tot
Itrue = np.array([[np.trapezoid(cf[k] * np.exp(-zf / (Lm[k] * cosd(t))), zf) for t in thm] for k in range(3)])
Gm = 0.7 + 0.3 * cosd(thm)
Im = Itrue * Gm[None, :]
sIm = 0.01 * Im
Im = Im + sIm * rng.normal(size=Im.shape)
a0 = int(np.argmin(thm)); m0 = Im[:, a0] / Im[:, a0].sum()
mdef = np.tile(np.maximum(m0, 1e-6)[:, None], (1, J_)); mdef = mdef / mdef.sum(0)

def softmax(u):
    u = u.reshape(3, J_); e = np.exp(u - u.max(0)); return (e / e.sum(0))

def Q(u, alpha):
    c = softmax(u)
    M = np.einsum("kj,kaj->ka", c, W)
    w = 1 / sIm ** 2
    G = (w * Im * M).sum(0) / (w * M * M).sum(0)
    R = (Im - G[None, :] * M) / sIm
    chi = (R ** 2).sum()
    cc = np.maximum(c, 1e-300)              # SPEC_multi §4.5: c floored at 1e-300 inside the logarithms
    S = (cc - mdef - cc * np.log(cc / mdef)).sum()
    g = np.einsum("ka,a,kaj->kj", -R / sIm, G, W) + alpha * np.log(cc / mdef)
    gu = c * (g - (c * g).sum(0)[None, :])
    return 0.5 * chi - alpha * S, gu.ravel(), chi, S, c

u_test = rng.normal(size=3 * J_) * 0.3
f, gu, chi, S, _ = Q(u_test, 2.0)
out["maxent_obj"] = dict(u=L(u_test), alpha=2.0, f=f, grad=L(gu), chi2=chi, S=S)

def solve(alpha, u0=None):
    u0 = np.log(mdef).ravel() if u0 is None else u0
    r = minimize(lambda u: Q(u, alpha)[:2], u0, jac=True, method="L-BFGS-B", options=dict(maxiter=50000, maxcor=20, ftol=1e-16, gtol=1e-11))
    _, _, chi, S, c = Q(r.x, alpha)
    return r.x, chi, S, c

sols = []
for alpha in [10.0, 30.0, 300.0]:   # α below the discrepancy value (≈4) is ill-posed: minimizers are not unique
    u, chi, S, c = solve(alpha)
    sols.append(dict(alpha=alpha, chi2=chi, S=S, c=L(c)))
# discrepancy α: χ²(α) = K·A − A, continuation from large α (half-decade ladder, warm starts), then bisection in ln α
# warm-started from the larger-α end of the bracket (SPEC_multi §4.5)
N = 3 * len(thm) - len(thm)
prev = None; k = 0
while True:
    a = 1e6 * 10 ** (-k / 2)
    if a < 1e-4 * 0.999: break
    uu, chi, S, c = solve(a, prev[0] if prev else None)
    cur = (uu, chi, a)
    if chi <= N: break
    prev = cur; k += 1
hi, lo = np.log(prev[2]), np.log(cur[2]); uHi = prev[0]; alphaN = cur[2]
for it in range(40):
    mid = 0.5 * (lo + hi); uu, chi, S, c = solve(np.exp(mid), uHi)
    if chi > N: hi = mid; uHi = uu
    else: lo = mid; alphaN = np.exp(mid); cN = c; chiN = chi
    if hi - lo < 1e-4: break
out["maxent"] = dict(th=L(thm), L=L(Lm), J=J_, z=L(z), I=L(Im), sI=L(sIm), W00=L(W[0, 0]), m=L(mdef[:, 0]), solutions=sols, alphaN=float(alphaN), chiN=float(chiN), cN=L(cN), N=N,
                     truth=dict(z=L(zf[::1000][:200]), cA=L(cf[0][::1000][:200]), cB=L(cf[1][::1000][:200])))
json.dump(out, open(os.path.join(HERE, "ref_arxps.json"), "w"))
print("wrote ref_arxps.json; MaxEnt α(χ²=N) =", alphaN, "χ² =", chiN)
