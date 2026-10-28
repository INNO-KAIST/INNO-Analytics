#!/usr/bin/env python3
"""Independent reference for SPEC.md §12–§15 and §2 spline / §3.3b Simplex / §4 at% Monte Carlo (numpy / scipy only).

Written from SPEC.md, not from xps_core.js. Different numerical routes on purpose:
  - interpolation: scipy.interpolate.PchipInterpolator / numpy.interp (the SPEC names PCHIP = scipy's definition)
  - template broadening: adaptive quadrature (scipy.integrate.quad) of interpolant × Gaussian
  - LCF: scipy.optimize.nnls (linear case) and scipy.optimize.least_squares (shift / broadening free)
  - spline background: scipy.interpolate.CubicSpline(bc_type='natural')
  - Simplex: scipy.optimize.minimize(method='Nelder-Mead') and least_squares (true minimum)
  - collinearity: numpy.linalg.cond / numpy.corrcoef + inv
  - at% Monte Carlo: numpy Generator (independent random stream)
Output: test/ref/ref_casa.json (read by test/test_xps_casa.js)
Run:    python3 -I test/ref/ref_casa.py            (~1 min)
"""
import json, os
import numpy as np
from scipy import integrate, interpolate, optimize, special

HERE = os.path.dirname(os.path.abspath(__file__))
LN2 = np.log(2.0)
F2S = 1.0 / (2.0 * np.sqrt(2.0 * LN2))
rng = np.random.default_rng(20261028)


def G(x, F):
    return (2.0 / F) * np.sqrt(LN2 / np.pi) * np.exp(-4 * LN2 * x * x / (F * F))


def L(x, F):
    return (2.0 / (np.pi * F)) / (1.0 + 4.0 * x * x / (F * F))


def gls(x, F, m):
    return (1 - m) * G(x, F) + m * L(x, F)


def voigt(x, wG, wL):
    return special.voigt_profile(np.asarray(x, float), wG * F2S, wL / 2.0)


# ---- SPEC §2: endpoints, linear background, local quadratic smoothing; §3.1 variance ----
def endpoints(E, y, n=5):
    N = len(E)
    n = max(1, min(n, N // 4 or 1))
    return y[:n].mean(), y[-n:].mean(), E[:n].mean(), E[-n:].mean()


def bg_linear(E, y, n=5):
    ylo, yhi, elo, ehi = endpoints(E, y, n)
    return ylo + (yhi - ylo) / (ehi - elo) * (E - elo)


def smooth(E, y):
    st = (E[-1] - E[0]) / (len(E) - 1)
    hw = max(2, int(np.floor(0.25 / st + 0.5)))
    out = np.empty_like(y)
    for i in range(len(y)):
        a, b = max(0, i - hw), min(len(y) - 1, i + hw)
        t = E[a:b + 1] - E[i]
        if b - a < 3:
            out[i] = y[a:b + 1].mean()
            continue
        out[i] = np.polyfit(t, y[a:b + 1], 2)[-1]
    return out


def var_from(yhat, c):
    c = c if c and c > 0 else 1.0
    vmin = max(1.0 / c, 1e-6 * yhat.max())
    return np.maximum(yhat, vmin) / c


# ---- SPEC §12: template ----
def interp_fn(x, v, kind):
    if kind == 'linear':
        return lambda t: np.where((t >= x[0]) & (t <= x[-1]), np.interp(t, x, v), 0.0)
    p = interpolate.PchipInterpolator(x, v, extrapolate=False)
    return lambda t: np.nan_to_num(p(np.asarray(t, float)), nan=0.0)


def interp_integral(x, v, kind):
    if kind == 'linear':
        return np.trapezoid(v, x)
    return float(interpolate.PchipInterpolator(x, v).integrate(x[0], x[-1]))


def make_template(E, y, shift, eref, kind, n=5):
    z = y - bg_linear(E, y, n)
    x = E + shift - eref
    I = interp_integral(x, z, kind)
    return x, z / I, I


def tpl_conv_quad(x, v, kind, wG, X):
    f = interp_fn(x, v, kind)
    out = []
    for xi in X:
        g = lambda t: float(f(np.array([t]))[0]) * G(xi - t, wG)
        a, b = max(x[0], xi - 8 * wG), min(x[-1], xi + 8 * wG)
        if b <= a:
            out.append(0.0)
            continue
        pts = [p for p in x if a < p < b]
        s = 0.0
        edges = [a] + pts + [b]
        for k in range(len(edges) - 1):
            s += integrate.quad(g, edges[k], edges[k + 1], epsabs=0, epsrel=1e-11, limit=200)[0]
        out.append(s)
    return np.array(out)


def tpl_conv_fine(x, v, kind, wG, X, h=0.0005):
    """second route (used inside least_squares): very fine grid + numpy convolution, linear read-out"""
    f = interp_fn(x, v, kind)
    if wG <= 1e-9:
        return f(X)
    sg = wG * F2S
    g0, g1 = x[0] - 6 * sg, x[-1] + 6 * sg
    t = np.arange(g0, g1 + h, h)
    base = f(t)
    K = int(np.ceil(6 * sg / h))
    kk = np.arange(-K, K + 1) * h
    ker = np.exp(-0.5 * (kk / sg) ** 2)
    ker /= ker.sum()
    conv = np.convolve(base, ker, mode='same')
    return np.interp(X, t, conv, left=0.0, right=0.0)


# ---------------------------------------------------------------- fixtures
out = {}

# 1) PCHIP: values (nodes, between, outside) and exact integral
pc = []
for case in range(4):
    nn = [7, 20, 60, 3][case]
    x = np.sort(rng.uniform(0, 10, nn)) if case != 3 else np.array([0.0, 1.0, 2.5])
    x = np.unique(x)
    v = rng.normal(0, 1, len(x)) if case % 2 else np.abs(np.sin(x)) + 0.1 * rng.random(len(x))
    if case == 2:
        v[10:14] = 0.3   # flat run
        v[20] = v[21]
    q = np.concatenate([x, np.linspace(x[0] - 0.5, x[-1] + 0.5, 157)])
    p = interpolate.PchipInterpolator(x, v, extrapolate=False)
    val = np.nan_to_num(p(q), nan=0.0)
    pc.append({'x': x.tolist(), 'v': v.tolist(), 'q': q.tolist(), 'val': val.tolist(),
               'slopes': p.derivative()(x).tolist(), 'int': float(interpolate.PchipInterpolator(x, v).integrate(x[0], x[-1])),
               'int_lin': float(np.trapezoid(v, x))})
out['pchip'] = pc

# 2) template + broadening (quad)
n = 241
E = 280 + 0.05 * np.arange(n)
truth = 30000 * gls(E - 284.8, 1.1, 0.3) + 2000 + 15 * (E - 280)
yref = rng.poisson(truth).astype(float)
tp = []
for kind in ('pchip', 'linear'):
    x, v, I = make_template(E, yref, 0.37, 285.0, kind)
    X = np.linspace(x[0] - 1, x[-1] + 1, 61)
    br = {}
    for wG in (0.3, 0.8):
        br[str(wG)] = tpl_conv_quad(x, v, kind, wG, X).tolist()
    tp.append({'kind': kind, 'E': E.tolist(), 'y': yref.tolist(), 'shift': 0.37, 'eref': 285.0, 'x': x.tolist(), 'v': v.tolist(),
               'area0': I, 'X': X.tolist(), 'plain': interp_fn(x, v, kind)(X).tolist(), 'broad': br})
out['tpl'] = tp

# 3) interpolation error of sampled analytic shapes (documentation number for SPEC §12.2)
ie = []
for step in (0.025, 0.05, 0.1, 0.2):
    xs = np.arange(-8, 8 + step / 2, step)
    for nm, fn in (('GLS F=1 m=0.3', lambda t: gls(t, 1.0, 0.3)), ('Voigt wG=0.8 wL=0.3', lambda t: voigt(t, 0.8, 0.3))):
        xf = np.linspace(-7.5, 7.5, 30001)
        tv = fn(xf)
        e_p = np.max(np.abs(interpolate.PchipInterpolator(xs, fn(xs))(xf) - tv)) / tv.max()
        e_l = np.max(np.abs(np.interp(xf, xs, fn(xs)) - tv)) / tv.max()
        ie.append({'step': step, 'shape': nm, 'pchip': e_p, 'linear': e_l})
out['interp_err'] = ie

# 4) LCF, linear case: two references, static linear background, 2-pass IRLS (SPEC §3.1) → NNLS
n = 301
E = 280 + 0.05 * np.arange(n)
r1 = 200000 * gls(E - 284.8, 1.0, 0.3) + 500
r2 = 200000 * (0.6 * gls(E - 285.9, 1.4, 0.2) + 0.4 * gls(E - 286.6, 1.2, 0.2)) + 500
y1, y2 = rng.poisson(r1).astype(float), rng.poisson(r2).astype(float)
T = []
for yy in (y1, y2):
    x, v, I = make_template(E, yy, 0.0, 285.0, 'pchip')
    T.append({'x': x.tolist(), 'v': v.tolist()})
tgt = 1000 + 20 * (E - 280) + 6000 * gls(E - 284.8, 1.0, 0.3) + 4000 * (0.6 * gls(E - 285.9, 1.4, 0.2) + 0.4 * gls(E - 286.6, 1.2, 0.2))
y = rng.poisson(tgt).astype(float)
c = 1.0
B = bg_linear(E, y)
cols = np.array([interp_fn(np.array(t['x']), np.array(t['v']), 'pchip')(E - 285.0) for t in T]).T
v1 = var_from(smooth(E, y), c)
a1, _ = optimize.nnls(cols / np.sqrt(v1)[:, None], (y - B) / np.sqrt(v1))
mu1 = B + cols @ a1
v2 = var_from(mu1, c)
a2, _ = optimize.nnls(cols / np.sqrt(v2)[:, None], (y - B) / np.sqrt(v2))
mu2 = B + cols @ a2
chi2 = float(np.sum((y - mu2) ** 2 / v2))
Jw = cols / np.sqrt(v2)[:, None]
cov = np.linalg.inv(Jw.T @ Jw) * chi2 / (n - 2)
Xs = cols / np.sqrt(v2)[:, None]
Xn = Xs / np.linalg.norm(Xs, axis=0)
cond = float(np.linalg.cond(Xn))
vif = np.diag(np.linalg.inv(np.corrcoef(Xs.T))).tolist()
out['lcf_lin'] = {'E': E.tolist(), 'y': y.tolist(), 'c': c, 'tpl': T, 'eref': 285.0, 'areas': a2.tolist(), 'chi2': chi2,
                  'sig': np.sqrt(np.diag(cov)).tolist(), 'v2': v2.tolist(), 'cond': cond, 'vif': vif,
                  'frac': (100 * a2 / a2.sum()).tolist()}

# 5) LCF with shift + broadening free (least_squares, TRF, bounds), 2-pass IRLS
tgt2 = 1000 + 20 * (E - 280) + 6000 * tpl_conv_fine(np.array(T[0]['x']), np.array(T[0]['v']), 'pchip', 0.5, E - 285.12) \
    + 4000 * interp_fn(np.array(T[1]['x']), np.array(T[1]['v']), 'pchip')(E - 284.95)
y = rng.poisson(tgt2).astype(float)
B = bg_linear(E, y)
x1, v1t, x2, v2t = (np.array(T[0]['x']), np.array(T[0]['v']), np.array(T[1]['x']), np.array(T[1]['v']))


def model(p):
    a1_, be1, w1, a2_, be2 = p
    return B + a1_ * tpl_conv_fine(x1, v1t, 'pchip', w1, E - be1) + a2_ * interp_fn(x2, v2t, 'pchip')(E - be2)


lo = [0, 285.0 + 0.1 - 3, 0, 0, 285.0 - 3]
hi = [np.inf, 285.0 + 0.1 + 3, 8, np.inf, 285.0 + 3]
p0 = [5000, 285.1, 0.3, 5000, 285.0]
vv = var_from(smooth(E, y), 1.0)
res = optimize.least_squares(lambda p: (y - model(p)) / np.sqrt(vv), p0, bounds=(lo, hi), x_scale=[1000, 0.05, 0.1, 1000, 0.05], xtol=1e-14, ftol=1e-14, gtol=1e-14)
vv = var_from(model(res.x), 1.0)
res = optimize.least_squares(lambda p: (y - model(p)) / np.sqrt(vv), res.x, bounds=(lo, hi), x_scale=[1000, 0.05, 0.1, 1000, 0.05], xtol=1e-14, ftol=1e-14, gtol=1e-14)
chi2 = float(np.sum(res.fun ** 2))
J = res.jac
covs = np.linalg.inv(J.T @ J) * chi2 / (n - 5)
out['lcf_shift'] = {'E': E.tolist(), 'y': y.tolist(), 'x': res.x.tolist(), 'sig': np.sqrt(np.diag(covs)).tolist(), 'chi2': chi2,
                    'start': p0, 'names': ['a.area', 'a.be', 'a.wG', 'b.area', 'b.be']}

# 6) spline background (natural CubicSpline through knots; outside end knots: straight line with end slope)
sp = []
for case in range(3):
    n = 201
    E = 520 + 0.08 * np.arange(n)
    y = 3000 + 400 * np.sin((E - 520) / 4) + 20000 * gls(E - 531, 1.4, 0.3) + rng.normal(0, 30, n)
    if case == 0:
        bg = {'n': 5, 'knots': [{'e': 524.0}, {'e': 527.5, 'y': 3100.0}, {'e': 535.3}], 'kh': 2}
    elif case == 1:
        bg = {'n': 3, 'ends': False, 'knots': [{'e': 523.0}, {'e': 528.0}, {'e': 534.0}, {'e': 531.0, 'y': 3200.0}], 'kh': 3}
    else:
        bg = {'n': 5, 'knots': []}
    ks = []
    if bg.get('ends', True):
        ylo, yhi, elo, ehi = endpoints(E, y, bg['n'])
        ks += [(elo, ylo), (ehi, yhi)]
    for k in bg['knots']:
        if 'y' in k:
            ks.append((k['e'], k['y']))
        else:
            i = int(np.argmin(np.abs(E - k['e'])))
            a, b = max(0, i - bg['kh']), min(n - 1, i + bg['kh'])
            ks.append((k['e'], y[a:b + 1].mean()))
    ks.sort()
    xs, ys = np.array([k[0] for k in ks]), np.array([k[1] for k in ks])
    if len(xs) >= 3:
        cs = interpolate.CubicSpline(xs, ys, bc_type='natural')
    else:
        cs = interpolate.CubicSpline(xs, ys, bc_type='natural') if len(xs) == 2 else None
    Bv = cs(E)
    lo_, hi_ = E < xs[0], E > xs[-1]
    Bv[lo_] = cs(xs[0]) + cs(xs[0], 1) * (E[lo_] - xs[0])
    Bv[hi_] = cs(xs[-1]) + cs(xs[-1], 1) * (E[hi_] - xs[-1])
    sp.append({'E': E.tolist(), 'y': y.tolist(), 'bg': bg, 'B': Bv.tolist(), 'knots': [list(k) for k in ks]})
out['spline'] = sp

# 7) Simplex: same chi2 (fixed variance v) for GLS×2 + static linear background
n = 241
E = 282 + 0.05 * np.arange(n)
tr = 1500 + 10 * (E - 282) + 8000 * gls(E - 284.8, 1.2, 0.3) + 3000 * gls(E - 286.1, 1.3, 0.3)
y = rng.poisson(tr).astype(float)
B = bg_linear(E, y)
v = var_from(smooth(E, y), 1.0)
start = [284.9, 1.0, 0.25, 7000, 286.0, 1.5, 0.35, 3500]   # be, F, m, area (×2)
blo = [start[0] - 3, 0.05, 0, 0, start[4] - 3, 0.05, 0, 0]
bhi = [start[0] + 3, 10, 1, np.inf, start[4] + 3, 10, 1, np.inf]


def mdl(p):
    return B + p[3] * gls(E - p[0], p[1], p[2]) + p[7] * gls(E - p[4], p[5], p[6])


def chi(p):
    return float(np.sum((y - mdl(p)) ** 2 / v))


nm = optimize.minimize(chi, start, method='Nelder-Mead', bounds=list(zip(blo, bhi)),
                       options={'xatol': 1e-10, 'fatol': 1e-10, 'maxfev': 400000, 'maxiter': 400000, 'adaptive': True})
for _ in range(3):   # restarts at the found point (standard practice; scipy has no built-in restart)
    nm = optimize.minimize(chi, nm.x, method='Nelder-Mead', bounds=list(zip(blo, bhi)),
                           options={'xatol': 1e-10, 'fatol': 1e-10, 'maxfev': 400000, 'maxiter': 400000, 'adaptive': True})
ls = optimize.least_squares(lambda p: (y - mdl(p)) / np.sqrt(v), start, bounds=(blo, bhi), xtol=1e-15, ftol=1e-15, gtol=1e-15,
                            x_scale=[0.05, 0.05, 0.05, 500, 0.05, 0.05, 0.05, 500])
J = ls.jac
c2 = float(np.sum(ls.fun ** 2))
sg = np.sqrt(np.diag(np.linalg.inv(J.T @ J) * c2 / (n - 8)))
out['simplex'] = {'E': E.tolist(), 'y': y.tolist(), 'v': v.tolist(), 'start': start, 'nm_x': nm.x.tolist(), 'nm_chi2': float(nm.fun),
                  'ls_x': ls.x.tolist(), 'ls_chi2': c2, 'sig': sg.tolist()}

# 8) at% Monte Carlo (areas ~ N(A, sA) independent; at% = n_i/Σn)
rows = [{'area': 52000.0, 'sArea': 900.0, 'S': 1.0}, {'area': 61000.0, 'sArea': 2400.0, 'S': 2.93},
        {'area': 9800.0, 'sArea': 1500.0, 'S': 7.0}, {'area': 3100.0, 'sArea': 1200.0, 'S': 1.8}]
A = np.array([r['area'] for r in rows])
sA = np.array([r['sArea'] for r in rows])
S = np.array([r['S'] for r in rows])
M = 400000
draws = (A + sA * rng.standard_normal((M, len(rows)))) / S
at = 100 * draws / draws.sum(axis=1, keepdims=True)
nn = A / S
tot = nn.sum()
delta = []
for i in range(len(rows)):
    g = np.array([100 * ((tot if i == j else 0) - nn[i]) / tot ** 2 for j in range(len(rows))])
    delta.append(float(np.sqrt(np.sum((g * sA / S) ** 2))))
out['quantmc'] = {'rows': rows, 'at': (100 * nn / tot).tolist(), 'mean': at.mean(0).tolist(), 'sd': at.std(0, ddof=1).tolist(),
                  'lo': np.percentile(at, 2.5, axis=0).tolist(), 'hi': np.percentile(at, 97.5, axis=0).tolist(), 'delta': delta, 'M': M}

# 9) CasaXPS region convention (SPEC §2b) on the CasaXPS validation files, Avantage Shirley (§2c), Avantage tail shape (§1.10),
#    Avantage transmission (§4), Tougaard with threshold T0 (§2)
def read_vms(path):
    L = open(path).read().split('\n'); k = L.index('762'); v = [float(L[k + 5 + 2 * i]) for i in range(381)]
    E = 279.08 + 0.05 * np.arange(381)
    return E, np.array(v[::-1])   # file is KE ascending = BE descending


def casa_ends(y, i0, i1, w):
    a, b = i0 + 1, i1 - 1
    return y[max(0, a - w):a + w + 1].mean(), y[b - w:min(len(y), b + w + 1)].mean()


def casa_bg(E, y, i0, i1, w, kind):
    ylo, yhi = casa_ends(y, i0, i1, w)
    R = slice(i0, i1 + 1); Er, yr = E[R], y[R]; a, b = 1, len(Er) - 2
    if kind == 'linear':
        return Er, yr, ylo + (yhi - ylo) * (Er - Er[a]) / (Er[b] - Er[a])
    B = np.full(len(Er), ylo); B[b:] = yhi
    for _ in range(200):
        r = yr - B
        Q = np.concatenate([[0.0], np.cumsum(r[a:b])])          # Q[k] = sum_{j=a}^{a+k-1}
        nb = np.empty_like(B); nb[:a + 1] = ylo; nb[b:] = yhi
        nb[a:b + 1] = ylo + (yhi - ylo) * Q / Q[-1]
        if np.max(np.abs(nb - B)) < 1e-7 * abs(yhi - ylo): B = nb; break
        B = nb
    return Er, yr, B


def casa_stats(Er, yr, B):
    a, b = 1, len(Er) - 2; r = (yr - B)[a:b + 1]
    return float(np.trapezoid(r, Er[a:b + 1])), float(r.max()), float(r.min())


cv = []
for f, kind, w in [('C1s', 'linear', 0), ('C1s', 'linear', 1), ('C1s', 'linear', 3), ('C1s', 'shirley', 1), ('Syn2', 'shirley', 1), ('Syn1', 'linear', 1)]:
    E, y = read_vms(os.path.join(HERE, '..', 'data', 'casaval', f + '.vms'))
    Er, yr, B = casa_bg(E, y, 19, 361, w, kind)
    cv.append({'file': f, 'type': kind, 'aw': w, 'B': B.tolist(), 'stats': casa_stats(Er, yr, B)})
out['casa_conv'] = cv


def endpoints_n(E, y, n):
    n = max(1, min(n, len(E) // 4 or 1)); return y[:n].mean(), y[-n:].mean()


def cumtrapz(E, f):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(E))])


av = []
for f, it, cl in [('C1s', 10, 0.0), ('C1s', 1, 0.0), ('Syn2', 10, 1e-4), ('Syn1', 10, 0.0)]:
    E, y = read_vms(os.path.join(HERE, '..', 'data', 'casaval', f + '.vms'))
    E, y = E[19:362], y[19:362]
    ylo, yhi = endpoints_n(E, y, 5); d = yhi - ylo
    if d <= 0:
        lo_, hi_ = E[:5].mean(), E[-5:].mean(); B = ylo + (yhi - ylo) * (E - lo_) / (hi_ - lo_); fb = 'linear'
    else:
        B = np.full(len(E), ylo); Aold = None; fb = None
        for _ in range(it):
            Q = cumtrapz(E, y - B); B = ylo + d * Q / Q[-1]
            A = cumtrapz(E, y - B)[-1]
            if Aold is not None and abs(A - Aold) <= cl * abs(A): break
            Aold = A
    av.append({'file': f, 'it': it, 'cl': cl, 'E': E.tolist(), 'y': y.tolist(), 'B': B.tolist(), 'fallback': fb})
out['avantage_shirley'] = av


def agt_raw(t, F, m, ET, TM, CT, ch, tu=None):
    g = np.exp(-4 * LN2 * (1 - m) * t * t / F / F) / (1 + 4 * m * t * t / F / F)
    if t <= 0: return g
    u = F if tu == 'fwhm' else ch   # CasaXPS GL(m)T(k): Dx in units of f
    return g + (1 - g) * (TM * CT + (1 - TM) * np.exp(-ET * t / u))


ag = []
for prm in [dict(F=1.2, m=0.3, ET=0.05, TM=0.3, CT=0.1, ch=0.05), dict(F=0.9, m=0.0, ET=0.1, TM=0.0, CT=0.0, ch=0.1), dict(F=1.5, m=0.6, ET=0.02, TM=1.0, CT=0.05, ch=0.05), dict(F=1.2, m=0.3, ET=0.8, TM=0.0, CT=0.0, ch=0.05, tu='fwhm')]:
    lo, hi = -6.0, 12.0
    N = integrate.quad(lambda t: agt_raw(t, **prm), lo, 0, epsabs=0, epsrel=1e-13, limit=400)[0] + integrate.quad(lambda t: agt_raw(t, **prm), 0, hi, epsabs=0, epsrel=1e-13, limit=400)[0]
    X = np.linspace(lo, hi, 181)
    ag.append({'p': prm, 'lo': lo, 'hi': hi, 'x': X.tolist(), 'v': [agt_raw(t, **prm) / N for t in X]})
out['agt'] = ag
coef = [0.12, -0.83, 0.05, -0.004]
tf = []
for KE, PE in [(1200, 50), (500, 20), (1450, 100), (300, 150)]:
    x = np.log10(KE / PE); tf.append([KE, PE, PE * 10 ** sum(c * x ** k for k, c in enumerate(coef))])
out['tf_avantage'] = {'coef': coef, 'rows': tf}
E, y = read_vms(os.path.join(HERE, '..', 'data', 'casaval', 'C1s.vms')); E, y = E[19:362], y[19:362]
ylo, yhi = endpoints_n(E, y, 5); elo, ehi = E[:5].mean(), E[-5:].mean()
tg = []
for C_, D_, u3, T0 in [(1643, 0, False, 0), (1643, 0, False, 3.0), (551, 436, True, 3.0)]:
    K = lambda T: 0.0 if T <= T0 else (T / ((C_ - T * T) ** 2 + D_ * T * T) if u3 else T / (C_ + T * T) ** 2)
    I = np.zeros(len(E))
    for i in range(1, len(E)):
        fa = np.array([K(E[i] - E[j]) * (y[j] - ylo) for j in range(i + 1)])
        I[i] = np.sum(0.5 * (fa[1:] + fa[:-1]) * np.diff(E[:i + 1]))
    Ih = np.interp(ehi, E, I); sc = max(0.0, (yhi - ylo) / Ih)
    tg.append({'C': C_, 'D': D_, 'u3': u3, 'T0': T0, 'B': (ylo + sc * I).tolist()})
out['tougaard_T0'] = {'E': E.tolist(), 'y': y.tolist(), 'cases': tg}

json.dump(out, open(os.path.join(HERE, 'ref_casa.json'), 'w'))
print('ref_casa.json written:', {k: (len(v) if isinstance(v, list) else 'obj') for k, v in out.items()})
print('interp error (max/peak):')
for r in ie:
    print('  step %.3f  %-22s pchip %.2e  linear %.2e' % (r['step'], r['shape'], r['pchip'], r['linear']))
