#!/usr/bin/env python3
"""Independent Python reference for modules/xps/SPEC_proc.md (numpy / scipy).

Written from SPEC_proc.md, not from src/xps_proc.js. Third-party routines are used wherever SPEC_proc
says the JS must agree with them (scipy.signal.savgol_filter, scipy.ndimage.gaussian_filter1d,
scipy.interpolate.PchipInterpolator, scipy.special.voigt_profile / erf, scipy.optimize.least_squares,
scipy.integrate.quad). Everything else is re-implemented with a *different* algorithm from the JS
(dense matrices instead of loops, cumulative integrals instead of overlap loops).

Output: proc_ref.json next to this file.     Run:  python3 -I test/ref/proc_ref.py   (~1 min)
"""
import json, math, os, sys
import numpy as np
from scipy import signal, ndimage, interpolate, special, integrate, optimize, fft

HERE = os.path.dirname(os.path.abspath(__file__))
FW2SIG = 1 / (2 * math.sqrt(2 * math.log(2)))
KB = 8.617333262e-5
rng = np.random.default_rng(20261028)
L = lambda a: [float(v) for v in np.asarray(a, float)]
out = {"generator": "proc_ref.py", "spec": "SPEC_proc.md v1.0", "numpy": np.__version__}
import scipy; out["scipy"] = scipy.__version__

# ---------------------------------------------------------------- special functions
xs = np.concatenate([np.linspace(-6, 6, 49), [0.0, 1e-8, 2.4999, 2.5, 2.5001, 5.9]])
out["erf"] = {"x": L(xs), "erf": L(special.erf(xs)), "erfc": L(special.erfc(xs))}
vc = []
for (wG, wL) in [(0.7, 0.4), (1.0, 0.05), (0.2, 1.5), (0.5, 0.0), (0.0, 0.6)]:
    x = np.linspace(-5, 5, 41)
    if wG == 0: v = (2 / (np.pi * wL)) / (1 + 4 * x**2 / wL**2)
    else: v = special.voigt_profile(x, wG * FW2SIG, wL / 2)
    vc.append({"wG": wG, "wL": wL, "x": L(x), "v": L(v)})
out["voigt"] = vc

# ---------------------------------------------------------------- §3 Savitzky–Golay vs scipy
def sig(n):
    t = np.linspace(0, 1, n)
    return 1000 * np.exp(-((t - 0.4) / 0.08) ** 2) + 300 * np.exp(-((t - 0.7) / 0.03) ** 2) + 50 * t + rng.normal(0, 10, n)
sg = []
y50, y101 = sig(50), sig(101)
for y in (y50, y101):
    for w in (5, 7, 11, 21, 31):
        for order in range(0, 6):
            if order >= w: continue
            for deriv in (0, 1, 2, 3):
                for mode in ("interp", "mirror", "nearest", "constant", "wrap"):
                    if mode != "interp" and (order not in (2, 4) or deriv > 2): continue
                    cval = 2.5 if mode == "constant" else 0.0
                    r = signal.savgol_filter(y, w, order, deriv=deriv, delta=0.05, mode=mode, cval=cval)
                    sg.append({"n": len(y), "w": w, "order": order, "deriv": deriv, "delta": 0.05, "mode": mode, "cval": cval, "out": L(r)})
out["sg_y"] = {"50": L(y50), "101": L(y101)}
out["sg"] = sg
# savgol variance propagation: var_out = sum_j w_ij^2 var_j, rows from scipy via unit impulses (linear filter)
vv = np.abs(rng.normal(100, 30, 50))
cases = []
for (w, order, mode) in [(7, 2, "interp"), (11, 3, "mirror"), (9, 4, "nearest")]:
    W = np.array([signal.savgol_filter(np.eye(50)[k], w, order, mode=mode) for k in range(50)]).T
    cases.append({"w": w, "order": order, "mode": mode, "var_out": L((W**2) @ vv)})
out["sg_var"] = {"var": L(vv), "cases": cases}

# ---------------------------------------------------------------- §3 Gaussian vs scipy
ga = []
for n in (40, 20):
    y = sig(n)
    for s in (0.5, 1.3, 3.7, 7.0):
        for mode in ("reflect", "mirror", "nearest", "constant", "wrap"):
            ga.append({"y": L(y), "sigma": s, "mode": mode, "out": L(ndimage.gaussian_filter1d(y, s, mode=mode, truncate=4.0))})
out["gauss"] = ga

# ---------------------------------------------------------------- §2 resampling
xe = np.cumsum(np.concatenate([[0], rng.uniform(0.05, 0.15, 39)])) + 280.0
ye = 500 + 3000 * np.exp(-((xe - 282.0) / 0.5) ** 2) - 400 * np.exp(-((xe - 283.5) / 0.3) ** 2) + rng.normal(0, 20, xe.size)
ye[10:13] = ye[10]                       # flat run -> zero slopes
ve = np.abs(ye) * 0.3 + 5
xn = np.linspace(xe[0], xe[-1], 157)
pch = interpolate.PchipInterpolator(xe, ye)(xn)
# first-order variance by central finite differences of scipy's interpolant (independent of the analytic rows)
J = np.zeros((xn.size, xe.size))
for k in range(xe.size):
    h = 1e-6 * max(1.0, abs(ye[k])); yp = ye.copy(); ym = ye.copy(); yp[k] += h; ym[k] -= h
    J[:, k] = (interpolate.PchipInterpolator(xe, yp)(xn) - interpolate.PchipInterpolator(xe, ym)(xn)) / (2 * h)
out["pchip"] = {"x": L(xe), "y": L(ye), "var": L(ve), "xn": L(xn), "out": L(pch), "var_out": L((J**2) @ ve)}

def edges(x):
    e = np.empty(x.size + 1); e[1:-1] = (x[1:] + x[:-1]) / 2; e[0] = x[0] - (x[1] - x[0]) / 2; e[-1] = x[-1] + (x[-1] - x[-2]) / 2; return e
def rebin_ref(x, y, xn):
    """cumulative-integral formulation: C(e) = integral of piecewise-constant density; y_j = ΔC / covered width"""
    ei, eo = edges(x), edges(xn)
    C = np.concatenate([[0], np.cumsum(y * np.diff(ei))])
    def Cat(e): return np.interp(e, ei, C)
    a, b = np.clip(eo[:-1], ei[0], ei[-1]), np.clip(eo[1:], ei[0], ei[-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        return (Cat(b) - Cat(a)) / (b - a)
xr = np.linspace(xe[0] + 0.2, xe[-1] - 0.2, 23)
W = np.array([rebin_ref(xe, np.eye(xe.size)[k], xr) for k in range(xe.size)]).T
out["rebin"] = {"x": L(xe), "y": L(ye), "var": L(ve), "xn": L(xr), "out": L(rebin_ref(xe, ye, xr)), "var_out": L((W**2) @ ve)}
# partial-cover edges: output grid wider than input
xr2 = np.linspace(xe[0] - 0.1, xe[-1] + 0.1, 15)
out["rebin2"] = {"xn": L(xr2), "out": L(rebin_ref(xe, ye, xr2))}

xu = np.linspace(100, 104, 41); yu = 100 + 50 * np.sin(xu * 3) + rng.normal(0, 2, xu.size)
xq = np.linspace(100, 104, 133)
def lanczos_ref(x, y, xq, a):
    st = x[1] - x[0]; res = []
    for q in xq:
        u = (q - x[0]) / st; i = np.arange(x.size); t = u - i
        w = np.where(np.abs(t) < a, np.sinc(t) * np.sinc(t / a), 0.0)
        res.append(np.dot(w, y) / w.sum())
    return np.array(res)
out["lanczos"] = {"x": L(xu), "y": L(yu), "xn": L(xq), "a3": L(lanczos_ref(xu, yu, xq, 3)), "a5": L(lanczos_ref(xu, yu, xq, 5))}

# ---------------------------------------------------------------- §4 PSF kernels
def k_gauss(fw, st):
    s = fw * FW2SIG; R = max(1, math.ceil(6 * s / st + 0.5)); j = np.arange(-R, R + 1)
    h = special.ndtr((j + 0.5) * st / s) - special.ndtr((j - 0.5) * st / s); return h / h.sum()
def k_voigt(wG, wL, st, tail=1e-3):
    s = wG * FW2SIG; g = wL / 2; fw = 0.5346 * wL + math.sqrt(0.2166 * wL**2 + wG**2)
    X = max(6 * s, 2 * g / (math.pi * tail) if wL > 0 else 0); X = min(X, max(10 * fw, 6 * s))
    R = max(1, math.ceil(X / st + 0.5)); j = np.arange(-R, R + 1)
    h = np.array([integrate.quad(lambda x: special.voigt_profile(x, s, g), (k - 0.5) * st, (k + 0.5) * st, epsabs=1e-15, epsrel=1e-13)[0] for k in j])
    return h / h.sum()
out["kernels"] = [{"type": "gauss", "fwhm": fw, "step": st, "h": L(k_gauss(fw, st))} for fw, st in ((1.0, 0.05), (0.35, 0.1), (0.08, 0.05))] + \
                 [{"type": "voigt", "wG": a, "wL": b, "step": st, "h": L(k_voigt(a, b, st))} for a, b, st in ((0.8, 0.3, 0.05), (0.5, 0.05, 0.1))]

# ---------------------------------------------------------------- §1 noise model
yn = 200 + 5000 * np.exp(-((np.linspace(-3, 3, 121)) / 0.6) ** 2); ynp = rng.poisson(yn).astype(float)
def noise_ref(y):
    n = y.size; lv = np.array([y[max(0, i - 2):min(n, i + 3)].mean() for i in range(n)])
    fl = 1e-3 * np.abs(y).max(); d2 = y[:-2] - 2 * y[1:-1] + y[2:]
    kappa = np.median(d2**2 / (6 * np.maximum(lv[1:-1], fl))) / 0.45493642311957283
    s = 1.482602218505602 * np.median(np.abs(d2 - np.median(d2))); return kappa, s * s / 6
k_, s2_ = noise_ref(ynp)
out["noise"] = {"y": L(ynp), "kappa": float(k_), "s2": float(s2_), "chi2_1_median_check": float(special.chdtri(1, 0.5))}

# ---------------------------------------------------------------- §5 deconvolution (dense matrices)
def Hdense(M, h):
    """reflect boundary built from numpy's 'symmetric' padding of the index vector (independent of the JS map)"""
    R = (h.size - 1) // 2; ext = np.pad(np.arange(M), R, mode="symmetric"); H = np.zeros((M, M))
    for i in range(M):
        for j in range(-R, R + 1):
            H[i, ext[i - j + R]] += h[j + R]
    return H
def passband_ref(y, h, var, eps=0.01):
    """noise-only DCT modes (SPEC_proc §5.4) via scipy.fft.dct(type 2, ortho)"""
    N = y.size; R = (h.size - 1) // 2; j = np.arange(-R, R + 1); hs = 0.5 * (h + h[::-1])
    lam = np.array([np.sum(hs * np.cos(np.pi * m * j / N)) for m in range(N)])
    P = fft.dct(y, type=2, norm="ortho") ** 2 / var.mean()
    hw = max(2, N // 100); Pb = np.array([P[max(0, m - hw):m + hw + 1].mean() for m in range(N)])
    below = np.nonzero(Pb <= 1.5)[0]; mc = below[0] if below.size else N
    perp = (np.abs(lam) < eps) | ((np.arange(N) >= mc) & (np.abs(lam) < 0.5))
    return int(perp.sum()), float(P[perp].sum()), int(mc)
def deconv_ref(y, var, h, method, k_fixed=None, tau=1.0, r0=None, a=0.0, b=None, maxIter=20000, positive=True, reblur=True, rule="projected"):
    if r0 is None: r0 = 1.5 if method == "jansson" else 1.0
    N = y.size; R = 0; yp = y.copy(); M = N; H = Hdense(M, h)
    if b is None: b = 4 * y.max()
    if method == "rl": o = np.maximum(yp, 1e-6 * yp.max()); Ht1 = H.T @ np.ones(M)
    elif method == "jansson": o = np.clip(yp, a, b)
    else: o = np.maximum(yp, 0) if positive else yp.copy()
    def chi2(o): r = (yp - H @ o)[R:R + N]; return float(np.sum(r * r / var))
    hist = [chi2(o)]; k = 0
    if rule == "projected": nP, cP, _ = passband_ref(y, h, var); target = cP + tau * (N - nP)
    else: target = tau * N
    if k_fixed is None and hist[0] <= target: return o[R:R + N], 0, hist
    while True:
        if k >= (k_fixed if k_fixed is not None else maxIter): break
        Ho = H @ o
        if method == "rl":
            ratio = np.where(Ho > 1e-300, yp / np.where(Ho > 1e-300, Ho, 1), 1.0); o = o * (H.T @ ratio) / Ht1
        elif method == "jansson":
            g = H.T @ (yp - Ho) if reblur else (yp - Ho)
            r = np.maximum(0, r0 * (1 - 2 * np.abs(o - (a + b) / 2) / (b - a))); o = np.clip(o + r * g, a, b)
        elif method == "landweber":
            o = o + r0 * (H.T @ (yp - Ho)); o = np.maximum(o, 0) if positive else o
        else:
            o = o + r0 * (yp - Ho); o = np.maximum(o, 0) if positive else o
        k += 1; hist.append(chi2(o))
        if k_fixed is None and (hist[-1] <= target or (k >= 50 and hist[k // 2] - hist[-1] < 1)): break
    return o[R:R + N], k, hist
E = np.round(np.arange(-3, 3 + 1e-9, 0.05), 10)
truth = 150 + 4000 * (special.voigt_profile(E + 0.35, 0.25 * FW2SIG, 0.05) * 0.25 + 0.5 * special.voigt_profile(E - 0.35, 0.25 * FW2SIG, 0.05) * 0.25)
hk = k_gauss(0.8, 0.05)
blur = Hdense(E.size, hk) @ truth
yd = rng.poisson(blur).astype(float); vd = np.maximum(yd, 1)
dc = []
for method, kf, rb, rule in (("jansson", 40, True, "projected"), ("jansson", None, True, "projected"), ("jansson", None, True, "classic"), ("jansson", 25, False, "projected"),
                             ("landweber", 50, True, "projected"), ("landweber", None, True, "projected"), ("vancittert", 8, True, "projected"),
                             ("rl", 60, True, "projected"), ("rl", None, True, "projected"), ("rl", None, True, "classic")):
    o, k, hist = deconv_ref(yd, vd, hk, method, k_fixed=kf, reblur=rb, rule=rule)
    dc.append({"method": method, "iter": kf, "reblur": rb, "rule": rule, "k": k, "o": L(o), "chi2": hist[-1], "hist_head": L(hist[:5])})
nP, cP, mc = passband_ref(yd, hk, vd)
out["passband"] = {"nPerp": nP, "chi2Perp": cP, "mc": mc}
out["deconv"] = {"E": L(E), "y": L(yd), "var": L(vd), "psf_fwhm": 0.8, "cases": dc}
# adjoint check: H^T of the clamp operator (dense) applied to a random vector
rv = rng.normal(0, 1, 60); hs = k_gauss(0.3, 0.05); Hs = Hdense(60, hs)
out["adjoint"] = {"r": L(rv), "fwhm": 0.3, "step": 0.05, "HTr": L(Hs.T @ rv), "Hr": L(Hs @ rv)}

# ---------------------------------------------------------------- §6 Fermi edge (quad, independent of the closed-form split)
def fermi_quad(E, p, T=300.0, axis="BE"):
    s = 1 if axis == "BE" else -1; kT = KB * T; sg = p["wG"] * FW2SIG
    def g(u): return (p["A"] + p["S"] * (u - p["EF"])) * special.expit(s * (u - p["EF"]) / kT)
    res = []
    for e in E:
        lo, hi = e - 12 * sg, e + 12 * sg; pts = [p["EF"]] if lo < p["EF"] < hi else None
        f = lambda u: g(u) * math.exp(-0.5 * ((e - u) / sg) ** 2) / (sg * math.sqrt(2 * math.pi))
        v = integrate.quad(f, lo, hi, points=pts, epsabs=1e-12, epsrel=1e-12, limit=500)[0]
        res.append(p["b0"] + p["b1"] * (e - p["EF"]) + v)
    return np.array(res)
fm = []
for p, T, axis in (({"EF": 0.05, "wG": 0.35, "A": 2000, "S": -300, "b0": 50, "b1": 10}, 300, "BE"),
                   ({"EF": 0.0, "wG": 0.01, "A": 1000, "S": 100, "b0": 5, "b1": 0}, 300, "BE"),
                   ({"EF": 0.2, "wG": 0.12, "A": 800, "S": 40, "b0": 20, "b1": -5}, 30, "BE"),
                   ({"EF": 1486.3, "wG": 0.5, "A": 900, "S": -50, "b0": 30, "b1": 0}, 300, "KE")):
    Ef = np.linspace(p["EF"] - 1.2, p["EF"] + 1.2, 61)
    fm.append({"p": p, "T": T, "axis": axis, "E": L(Ef), "y": L(fermi_quad(Ef, p, T, axis))})
out["fermi_model"] = fm
# fit vs least_squares (same data, same variance)
Ef = np.round(np.arange(-1.0, 1.0 + 1e-9, 0.02), 10); pt = {"EF": 0.05, "wG": 0.35, "A": 2000, "S": -300, "b0": 50, "b1": 0}
yf = rng.poisson(fermi_quad(Ef, pt)).astype(float); vf = np.maximum(yf, 1)
names = ["EF", "wG", "A", "S", "b0", "b1"]
def res_f(q): return (yf - fermi_quad(Ef, dict(zip(names, q)))) / np.sqrt(vf)
ls = optimize.least_squares(res_f, [0.0, 0.3, 1900, 0, 60, 0], method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, x_scale=[0.02, 0.02, 2000, 1000, 50, 50])
Jm = ls.jac; cov = np.linalg.inv(Jm.T @ Jm)
out["fermi_fit"] = {"E": L(Ef), "y": L(yf), "var": L(vf), "p": dict(zip(names, L(ls.x))), "sigma": dict(zip(names, L(np.sqrt(np.diag(cov))))), "chi2": float(2 * ls.cost)}

# ---------------------------------------------------------------- §7 Voigt reference-peak fit
Ep = np.round(np.arange(366.0, 370.0 + 1e-9, 0.025), 10)
def peak(E, q): return q[1] * special.voigt_profile(E - q[0], q[2] * FW2SIG, q[3] / 2) + q[4] + q[5] * (E - q[0])
qt = [368.25, 30000.0, 0.55, 0.30, 400.0, 20.0]
yp_ = rng.poisson(peak(Ep, qt)).astype(float); vp = np.maximum(yp_, 1)
def res_p(q4): q = [q4[0], q4[1], q4[2], 0.30, q4[3], q4[4]]; return (yp_ - peak(Ep, q)) / np.sqrt(vp)
lsp = optimize.least_squares(res_p, [368.2, 28000, 0.5, 380, 0], method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, x_scale=[0.02, 1000, 0.02, 50, 50])
covp = np.linalg.inv(lsp.jac.T @ lsp.jac)
out["peak_fit"] = {"E": L(Ep), "y": L(yp_), "var": L(vp), "wL": 0.30, "p": dict(zip(["E0", "A", "wG", "b0", "b1"], L(lsp.x))),
                   "sigma": dict(zip(["E0", "A", "wG", "b0", "b1"], L(np.sqrt(np.diag(covp))))), "chi2": float(2 * lsp.cost)}

json.dump(out, open(os.path.join(HERE, "proc_ref.json"), "w"))
print("wrote proc_ref.json:", len(sg), "SG cases,", len(ga), "Gaussian cases,", len(dc), "deconvolution cases")
