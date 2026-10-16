#!/usr/bin/env python3
"""Independent reference implementation of modules/xps/SPEC.md (numpy / scipy only).

Written from SPEC.md, NOT from xps_core.js.  Different numerical routes on purpose:
  - line shapes: adaptive quadrature (scipy.integrate.quad) for every normalisation / convolution,
    scipy.special.voigt_profile for Voigt
  - Shirley: plain fixed-point iteration written independently
  - fitting: scipy.optimize.least_squares (TRF, native bounds) — not Levenberg–Marquardt
  - covariance: analytic-free (finite-difference J from least_squares) + same phi definition
Output: test/ref/ref.json (read by test/test_xps.js)
Run:    python3 -I test/ref/ref_xps.py
"""
import json, os, sys
import numpy as np
from scipy import integrate, special, optimize

HERE = os.path.dirname(os.path.abspath(__file__))
LN2 = np.log(2.0)
F2S = 1.0 / (2.0 * np.sqrt(2.0 * LN2))


# ---------------- line shapes (SPEC §1) ----------------
def G(x, F):
    return (2.0 / F) * np.sqrt(LN2 / np.pi) * np.exp(-4 * LN2 * x * x / (F * F))


def L(x, F):
    return (2.0 / (np.pi * F)) / (1.0 + 4.0 * x * x / (F * F))


def gls(x, F, m):
    return (1 - m) * G(x, F) + m * L(x, F)


def glp_raw(x, F, m):
    return np.exp(-4 * LN2 * (1 - m) * x * x / (F * F)) / (1 + 4 * m * x * x / (F * F))


def glp(x, F, m):
    N = 2 * integrate.quad(lambda t: glp_raw(t, F, m), 0, np.inf, limit=400, epsabs=0, epsrel=1e-12)[0]
    return glp_raw(np.asarray(x, float), F, m) / N


def voigt(x, wG, wL):
    return special.voigt_profile(np.asarray(x, float), wG * F2S, wL / 2.0)


def la_raw(t, F, a, b):
    e = a if t >= 0 else b   # SPEC 1.6 v1.1: alpha on the high-BE side
    return (1 + 4 * t * t / (F * F)) ** (-e)


def la(x, F, a, b, wG):
    N = integrate.quad(lambda t: la_raw(t, F, a, b), -np.inf, 0, limit=400, epsrel=1e-12)[0] + \
        integrate.quad(lambda t: la_raw(t, F, a, b), 0, np.inf, limit=400, epsrel=1e-12)[0]
    out = []
    for xi in np.atleast_1d(x):
        if wG <= 1e-9:
            out.append(la_raw(xi, F, a, b) / N)
            continue
        s = wG * F2S
        v = integrate.quad(lambda u: la_raw(xi - u, F, a, b) * np.exp(-0.5 * (u / s) ** 2) / (s * np.sqrt(2 * np.pi)),
                           -8 * s, 8 * s, points=[xi], limit=400, epsrel=1e-11)[0]
        out.append(v / N)
    return np.array(out)


def ds_raw(t, al, g):
    return np.cos(np.pi * al / 2 + (1 - al) * np.arctan(-t / g)) / (t * t + g * g) ** ((1 - al) / 2)


def ds_conv(xi, al, wL, wG):
    g = wL / 2
    if wG <= 1e-9:
        return ds_raw(xi, al, g)
    s = wG * F2S
    return integrate.quad(lambda u: ds_raw(xi - u, al, g) * np.exp(-0.5 * (u / s) ** 2) / (s * np.sqrt(2 * np.pi)),
                          -8 * s, 8 * s, points=[xi], limit=400, epsrel=1e-11)[0]


def ds(x, al, wL, wG, lo, hi):
    N = integrate.quad(lambda t: ds_conv(t, al, wL, wG), lo, hi, points=[0.0], limit=600, epsrel=1e-10)[0]
    return np.array([ds_conv(xi, al, wL, wG) for xi in np.atleast_1d(x)]) / N


def glt(x, F, m, T, tau):
    out = []
    for xi in np.atleast_1d(x):
        conv = integrate.quad(lambda t: np.exp(-t / tau) / tau * gls(xi - t, F, m), 0, np.inf, limit=600, epsrel=1e-11)[0]
        out.append((1 - T) * gls(xi, F, m) + T * conv)
    return np.array(out)


def shape(name, x, p, lo=None, hi=None):
    x = np.asarray(x, float)
    if name == 'GLS':
        return gls(x, p['F'], p['m'])
    if name == 'GLP':
        return glp(x, p['F'], p['m'])
    if name == 'V':
        return voigt(x, p['wG'], p['wL'])
    if name == 'LA':
        return la(x, p['F'], p['alpha'], p['beta'], p['wG'])
    if name == 'DS':
        return ds(x, p['alpha'], p['wL'], p['wG'], lo, hi)
    if name == 'GLT':
        return glt(x, p['F'], p['m'], p['T'], p['tau'])
    raise ValueError(name)


SHAPE_P = {'GLS': ['F', 'm'], 'GLP': ['F', 'm'], 'V': ['wG', 'wL'], 'DS': ['alpha', 'wL', 'wG'],
           'LA': ['F', 'alpha', 'beta', 'wG'], 'GLT': ['F', 'm', 'T', 'tau']}


def comp_eval(c, E, emin, emax):
    p = {k: c[k] for k in SHAPE_P[c['shape']]}
    if c.get('doublet'):
        d = c['doublet']
        A1, A2 = c['area'] / (1 + d['r']), c['area'] * d['r'] / (1 + d['r'])
        p2 = dict(p)
        for k in ('F', 'wL'):   # SPEC 1.9 v1.1: lifetime width only
            if k in p2:
                p2[k] *= d['q']
        e2 = c['be'] + d['dE']
        return A1 * shape(c['shape'], E - c['be'], p, emin - c['be'], emax - c['be']) + \
            A2 * shape(c['shape'], E - e2, p2, emin - e2, emax - e2)
    return c['area'] * shape(c['shape'], E - c['be'], p, emin - c['be'], emax - c['be'])


# ---------------- backgrounds (SPEC §2) ----------------
def cumtrapz(E, f):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(E))])


def ends(E, y, n):
    n = max(1, min(n, len(E) // 4 or 1))
    return y[:n].mean(), y[-n:].mean(), E[:n].mean(), E[-n:].mean()


def smooth_q(E, y, hw):
    out = np.empty_like(y)
    for i in range(len(y)):
        j0, j1 = max(0, i - hw), min(len(y) - 1, i + hw)
        t = E[j0:j1 + 1] - E[i]
        if j1 - j0 < 3:
            out[i] = y[j0:j1 + 1].mean()
            continue
        out[i] = np.polyval(np.polyfit(t, y[j0:j1 + 1], 2), 0.0)
    return out


def hw_of(E):
    st = (E[-1] - E[0]) / (len(E) - 1)
    return max(2, int(round(0.25 / st)))


def bg_linear(E, y, n):
    ylo, yhi, Elo, Ehi = ends(E, y, n)
    return ylo + (yhi - ylo) * (E - Elo) / (Ehi - Elo)


def bg_shirley(E, y, n, smart=False):
    ylo, yhi, _, _ = ends(E, y, n)
    d = yhi - ylo
    B = np.full_like(y, ylo)
    ys = smooth_q(E, y, hw_of(E)) if smart else None
    for _ in range(200):
        Q = cumtrapz(E, y - B)
        nb = ylo + d * Q / Q[-1]
        if smart:
            nb = np.minimum(nb, ys)
        done = np.max(np.abs(nb - B)) < 1e-7 * (abs(d) + 1e-12)
        B = nb
        if done:
            break
    return B


def bg_tougaard(E, y, n, C=1643.0, D=0.0, u3=False):
    ylo, yhi, Elo, Ehi = ends(E, y, n)
    K = (lambda T: T / ((C - T * T) ** 2 + D * T * T)) if u3 else (lambda T: T / (C + T * T) ** 2)
    I = np.zeros_like(y)
    for i in range(1, len(E)):
        T = E[i] - E[:i + 1]
        f = np.where(T > 0, K(np.maximum(T, 0)), 0.0) * (y[:i + 1] - ylo)
        I[i] = np.trapezoid(f, E[:i + 1])
    Ih = np.interp(Ehi, E, I)
    s = (yhi - ylo) / Ih
    return ylo + max(0.0, s) * I   # SPEC 2 v1.1: s >= 0


def bg_static(E, y, bg):
    t, n = bg['type'], bg.get('n', 5)
    if t == 'linear':
        return bg_linear(E, y, n)
    if t == 'shirley':
        return bg_shirley(E, y, n)
    if t == 'smart':
        return bg_shirley(E, y, n, True)
    if t == 'tougaard':
        return bg_tougaard(E, y, n, bg.get('C', 1643), bg.get('D', 0), bg.get('u3', False))
    if t == 'none':
        return np.zeros_like(y)
    return None


# ---------------- fit (SPEC §3) with scipy least_squares ----------------
DEF = {'F': (0.05, 10), 'm': (0, 1), 'wG': (0, 8), 'wL': (0.005, 8), 'alpha': (0, 0.45), 'beta': (0.55, 20),
       'T': (0, 1), 'tau': (0.02, 20), 'dE': (0.05, 40), 'r': (0.01, 5), 'q': (0.3, 3)}


def bounds_of(c, k):
    if k == 'be':
        return c['be'] - 3, c['be'] + 3
    if k == 'area':
        return 0.0, np.inf
    if c['shape'] == 'LA' and k == 'alpha':
        return 0.55, 20
    return DEF[k]


def names(c):
    a = ['be', 'area'] + SHAPE_P[c['shape']]
    if c.get('doublet'):
        a += ['dE', 'r', 'q']
    return a


def getp(c, k):
    return c['doublet'][k] if k in ('dE', 'r', 'q') else c[k]


def setp(c, k, v):
    if k in ('dE', 'r', 'q'):
        c['doublet'][k] = v
    else:
        c[k] = v


def fit(prob):
    E, y = np.asarray(prob['E']), np.asarray(prob['y'])
    comps = json.loads(json.dumps(prob['comps']))
    cons = prob.get('cons', {})
    ids = {c['id']: i for i, c in enumerate(comps)}
    free, links = [], []
    for ci, c in enumerate(comps):
        for k in names(c):
            cs = cons.get(c['id'], {}).get(k, {'mode': 'free'})
            if cs['mode'] == 'fixed':
                continue
            if cs['mode'] == 'link':
                links.append((ci, k, ids[cs['ref']], cs.get('k', 1.0), cs.get('b', 0.0)))
                continue
            lo, hi = bounds_of(c, k)
            free.append((ci, k, cs.get('lo', lo), cs.get('hi', hi)))
    bgt = prob['bg']['type']
    if bgt == 'active':
        free.append((-1, 'c0', -np.inf, np.inf))
        free.append((-1, 'k', 0.0, np.inf))
    bgp = dict(prob['bg'])
    Bst = bg_static(E, y, prob['bg'])
    emin, emax = E[0], E[-1]

    def put(v):
        for (ci, k, lo, hi), x in zip(free, v):
            if ci < 0:
                bgp[k] = x
            else:
                setp(comps[ci], k, x)
        for _ in range(4):
            for ci, k, ri, mul, add in links:
                setp(comps[ci], k, getp(comps[ri], k) * mul + add)

    def model(v):
        put(v)
        P = sum(comp_eval(c, E, emin, emax) for c in comps)
        B = Bst if Bst is not None else bgp['c0'] + bgp['k'] * cumtrapz(E, P)
        return B + P

    ymax = np.max(np.abs(y))
    cpc = prob.get('c') or 1.0

    def var(yh):
        vmin = max(1 / cpc, 1e-6 * float(np.max(yh)))
        return np.maximum(yh, vmin) / cpc

    v0 = np.array([prob['bg'][k] if ci < 0 else getp(comps[ci], k) for ci, k, lo, hi in free], float)
    lo = np.array([f[2] for f in free]); hi = np.array([f[3] for f in free])
    with np.errstate(invalid='ignore'):
        eps_lo = np.where(np.isfinite(lo), lo + 1e-9 * (np.abs(lo) + 1), -np.inf)
        eps_hi = np.where(np.isfinite(hi), hi - 1e-9 * (np.abs(hi) + 1), np.inf)
    w = var(smooth_q(E, y, hw_of(E)))
    # initial areas: weighted NNLS with the non-linear parameters fixed (own implementation via scipy.nnls)
    put(v0)
    Bi = Bst if Bst is not None else bg_shirley(E, y, prob['bg'].get('n', 5))
    area_idx = [j for j, f in enumerate(free) if f[1] == 'area']
    cols = []
    for j in area_idx:
        ci = free[j][0]
        col = np.zeros_like(y)
        grp = [ci] + [l[0] for l in links if l[2] == ci and l[1] == 'area']
        for g in grp:
            cc = json.loads(json.dumps(comps[g])); cc['area'] = 1.0 if g == ci else [l[3] for l in links if l[0] == g and l[1] == 'area'][0]
            col += comp_eval(cc, E, emin, emax)
        cols.append(col / np.sqrt(w))
    if cols:
        a, _ = optimize.nnls(np.array(cols).T, (y - Bi) / np.sqrt(w))
        for j, aj in zip(area_idx, a):
            v0[j] = max(aj, 1e-9 * (abs(y[0]) + 1))
    if bgt == 'active':
        put(v0)
        P0 = sum(comp_eval(c, E, emin, emax) for c in comps)
        ylo, yhi, _, _ = ends(E, y, prob['bg'].get('n', 5))
        Q = cumtrapz(E, P0)
        v0[-2] = ylo; v0[-1] = max(0.0, (yhi - ylo) / Q[-1])
    v0 = np.clip(v0, eps_lo, eps_hi)
    sol = None
    for it in range(2):
        sw = 1 / np.sqrt(w)
        sol = optimize.least_squares(lambda v: (y - model(v)) * sw, v0, bounds=(lo, hi), method='trf',
                                     x_scale='jac', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        v0 = sol.x
        mu = model(v0)
        if it == 0:
            w = var(mu)
    mu = model(v0)
    r = (y - mu) / np.sqrt(w)
    chi2 = float(np.sum(r * r))
    # covariance: central-difference Jacobian at the solution, same phi = chi2/(N-p)
    J = np.empty((len(E), len(v0)))
    for j in range(len(v0)):
        h = max(abs(v0[j]) * 1e-5, 1e-4 if free[j][1] in ('be', 'dE') else 1e-6)
        a, b = v0.copy(), v0.copy(); a[j] -= h; b[j] += h
        J[:, j] = (model(b) - model(a)) / (2 * h) / np.sqrt(w)
    put(v0)
    phi = chi2 / (len(E) - len(v0))
    Cv = np.linalg.inv(J.T @ J) * phi
    sig = np.sqrt(np.diag(Cv))
    return {'free': [[f[0], f[1]] for f in free], 'x': v0.tolist(), 'sig': sig.tolist(), 'chi2': chi2, 'phi': phi,
            'comps': comps, 'bg': {k: bgp.get(k) for k in ('c0', 'k')}, 'nfev': int(sol.nfev), 'status': int(sol.status)}


# ---------------- TPP-2M (SPEC §4) ----------------
def tpp2m(E, rho, Nv, M, Eg):
    U = Nv * rho / M
    Ep = 28.8 * np.sqrt(U)
    beta = -0.10 + 0.944 / np.sqrt(Ep ** 2 + Eg ** 2) + 0.069 * rho ** 0.1
    gamma = 0.191 / np.sqrt(rho)
    C = 1.97 - 0.91 * U
    D = 53.4 - 20.8 * U
    return E / (Ep ** 2 * (beta * np.log(gamma * E) - C / E + D / E ** 2))


# ---------------- cases ----------------
def main():
    out = {'shapes': [], 'bgs': [], 'fits': [], 'tpp2m': [], 'wofz': []}
    for xr in (0.0, 0.3, 1.0, 2.5, 5.0, 12.0, 40.0, -3.0):
        for yi in (0.0, 1e-4, 0.05, 0.4, 1.0, 3.0, 10.0, 60.0):
            w = special.wofz(complex(xr, yi))
            out['wofz'].append([xr, yi, w.real, w.imag])
    for xr in (0.0, 0.5, 1.0, 3.0, 8.0, 30.0):
        out['wofz'].append(['erfcx', xr, float(special.erfcx(xr)), 0.0])
    xs = np.array([-4, -2.5, -1.3, -0.7, -0.31, -0.05, 0, 0.05, 0.3, 0.6, 1.1, 1.9, 3.2, 5.0])
    cases = [('GLS', {'F': 1.2, 'm': 0.3}), ('GLS', {'F': 0.8, 'm': 1.0}),
             ('GLP', {'F': 1.2, 'm': 0.3}), ('GLP', {'F': 1.0, 'm': 0.0}), ('GLP', {'F': 1.0, 'm': 0.9}), ('GLP', {'F': 1.5, 'm': 1.0}),
             ('V', {'wG': 1.0, 'wL': 0.3}), ('V', {'wG': 0.3, 'wL': 1.2}), ('V', {'wG': 0.8, 'wL': 0.05}),
             ('LA', {'F': 1.0, 'alpha': 1.2, 'beta': 2.5, 'wG': 0.5}), ('LA', {'F': 0.6, 'alpha': 0.8, 'beta': 1.4, 'wG': 0.9}),
             ('LA', {'F': 1.0, 'alpha': 1.0, 'beta': 3.0, 'wG': 0.0}),
             ('GLT', {'F': 1.0, 'm': 0.3, 'T': 0.4, 'tau': 1.5}), ('GLT', {'F': 1.4, 'm': 0.1, 'T': 0.2, 'tau': 0.4}),
             ('DS', {'alpha': 0.1, 'wL': 0.3, 'wG': 0.6}), ('DS', {'alpha': 0.25, 'wL': 0.5, 'wG': 0.4}), ('DS', {'alpha': 0.05, 'wL': 0.2, 'wG': 0.0})]
    for s, p in cases:
        lo, hi = -6.0, 9.0
        v = shape(s, xs, p, lo, hi)
        out['shapes'].append({'shape': s, 'p': p, 'lo': lo, 'hi': hi, 'x': xs.tolist(), 'v': np.atleast_1d(v).tolist()})
        print('shape', s, p, 'max', float(np.max(v)))

    # backgrounds on a reproducible synthetic step spectrum (and on a real region, in the JS test)
    rs = np.random.default_rng(11)
    E = np.round(np.arange(279.0, 298.0001, 0.05), 6)
    P = 9000 * gls(E - 284.8, 1.1, 0.3) + 2500 * gls(E - 286.3, 1.3, 0.3)
    yb = 2000 + 0.09 * cumtrapz(E, P) + P
    y = rs.poisson(yb * 0.4) / 0.4
    for bg in ({'type': 'linear', 'n': 5}, {'type': 'shirley', 'n': 5}, {'type': 'smart', 'n': 5}, {'type': 'shirley', 'n': 1},
               {'type': 'tougaard', 'n': 5, 'C': 1643, 'D': 0, 'u3': False}, {'type': 'tougaard', 'n': 5, 'C': 551, 'D': 436, 'u3': True}):
        B = bg_static(E, y, bg)
        out['bgs'].append({'bg': bg, 'E': E.tolist(), 'y': y.tolist(), 'B': B.tolist()})

    # fits: synthetic data with known truth (data are generated here and shipped in ref.json)
    def synth(Ex, comps, bgfun, c, seed):
        mu = sum(comp_eval(cc, Ex, Ex[0], Ex[-1]) for cc in comps)
        mu = bgfun(Ex, mu)
        r = np.random.default_rng(seed)
        return (r.poisson(mu * c) / c).astype(float)

    def C(id, shape_, be, **k):
        d = {'id': id, 'shape': shape_, 'be': be, 'area': k.pop('area', 1.0), 'doublet': k.pop('doublet', None)}
        d.update(k)
        return d

    E1 = np.round(np.arange(279.0, 298.0001, 0.05), 6)
    t1 = [C('a', 'GLS', 284.8, F=1.1, m=0.3, area=9000), C('b', 'GLS', 286.3, F=1.3, m=0.3, area=2500), C('c', 'GLS', 288.8, F=1.2, m=0.3, area=1500)]
    s1 = [C('a', 'GLS', 284.6, F=1.0, m=0.3), C('b', 'GLS', 286.0, F=1.0, m=0.3), C('c', 'GLS', 289.0, F=1.0, m=0.3)]
    y1 = synth(E1, t1, lambda E, P: 2000 + 0.08 * cumtrapz(E, P) + P, 0.3, 1)
    fits = [('GLS3_shirley', E1, y1, 0.3, {'type': 'shirley', 'n': 5}, s1, {'b': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}, 'c': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}}),
            ('GLS3_smart', E1, y1, 0.3, {'type': 'smart', 'n': 5}, s1, {'b': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}, 'c': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}}),
            ('GLS3_active', E1, y1, 0.3, {'type': 'active', 'n': 5, 'c0': 2000.0, 'k': 0.07}, s1, {'b': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}, 'c': {'m': {'mode': 'link', 'ref': 'a', 'k': 1}}})]
    E2 = np.round(np.arange(771.0, 812.0001, 0.05), 6)
    t2 = [C('m', 'GLP', 778.2, F=1.0, m=0.4, area=60000, doublet={'dE': 15.05, 'r': 0.5, 'q': 1.3}),
          C('o', 'GLP', 780.3, F=2.4, m=0.3, area=50000, doublet={'dE': 15.7, 'r': 0.5, 'q': 1.0})]
    s2 = [C('m', 'GLP', 778.0, F=1.2, m=0.3, doublet={'dE': 15.05, 'r': 0.5, 'q': 1.0}),
          C('o', 'GLP', 780.6, F=2.0, m=0.3, doublet={'dE': 15.7, 'r': 0.5, 'q': 1.0})]
    y2 = synth(E2, t2, lambda E, P: 40000 + 0.25 * cumtrapz(E, P) + P, 0.25, 2)
    fits.append(('GLPdoublet_smart', E2, y2, 0.25, {'type': 'smart', 'n': 6}, s2,
                 {'m': {'dE': {'mode': 'fixed'}, 'r': {'mode': 'fixed'}}, 'o': {'dE': {'mode': 'fixed'}, 'r': {'mode': 'fixed'}, 'q': {'mode': 'fixed'}}}))
    E3 = np.round(np.arange(274.0, 292.0001, 0.05), 6)
    t3 = [C('v', 'V', 282.0, wG=0.7, wL=0.35, area=20000), C('w', 'V', 283.4, wG=0.9, wL=0.35, area=6000)]
    s3 = [C('v', 'V', 281.8, wG=0.9, wL=0.2), C('w', 'V', 283.7, wG=0.9, wL=0.2)]
    y3 = synth(E3, t3, lambda E, P: 3000 + 0.05 * cumtrapz(E, P) + P, 0.5, 3)
    fits.append(('Voigt2_tougaard', E3, y3, 0.5, {'type': 'tougaard', 'n': 5, 'C': 1643, 'D': 0, 'u3': False}, s3, {'w': {'wL': {'mode': 'link', 'ref': 'v', 'k': 1}}}))
    E4 = np.round(np.arange(772.0, 786.0001, 0.05), 6)
    t4 = [C('d', 'DS', 778.1, alpha=0.12, wL=0.45, wG=0.55, area=40000)]
    s4 = [C('d', 'DS', 778.0, alpha=0.05, wL=0.3, wG=0.7)]
    y4 = synth(E4, t4, lambda E, P: 30000 + 20 * (E - E[0]) + P, 0.25, 4)
    fits.append(('DS_linear', E4, y4, 0.25, {'type': 'linear', 'n': 5}, s4, {}))
    t5 = [C('l', 'LA', 778.1, F=0.9, alpha=1.0, beta=1.8, wG=0.5, area=40000)]
    s5 = [C('l', 'LA', 778.0, F=1.0, alpha=1.0, beta=1.4, wG=0.6)]
    y5 = synth(E4, t5, lambda E, P: 30000 + 0.15 * cumtrapz(E, P) + P, 0.25, 5)
    fits.append(('LA_shirley', E4, y5, 0.25, {'type': 'shirley', 'n': 5}, s5, {}))
    t6 = [C('t', 'GLT', 778.1, F=1.0, m=0.2, T=0.35, tau=1.2, area=40000)]
    s6 = [C('t', 'GLT', 778.0, F=1.2, m=0.3, T=0.2, tau=1.0)]
    y6 = synth(E4, t6, lambda E, P: 30000 + 0.1 * cumtrapz(E, P) + P, 0.25, 6)
    fits.append(('GLT_shirley', E4, y6, 0.25, {'type': 'shirley', 'n': 5}, s6, {}))

    for name, Ex, yx, c, bg, st, cons in fits:
        prob = {'E': Ex.tolist(), 'y': yx.tolist(), 'c': c, 'bg': bg, 'comps': st, 'cons': cons}
        r = fit(prob)
        out['fits'].append({'name': name, 'prob': prob, 'ref': r})
        print('fit', name, 'chi2', round(r['chi2'], 3), 'phi', round(r['phi'], 4), 'nfev', r['nfev'], 'status', r['status'])

    for mat, E_ in (({'rho': 2.33, 'Nv': 4, 'M': 28.086, 'Eg': 1.1}, [200, 500, 1000, 1486.6 - 99]),
                    ({'rho': 19.3, 'Nv': 11, 'M': 196.97, 'Eg': 0}, [300, 1000, 1402]),
                    ({'rho': 6.0, 'Nv': 24, 'M': 150, 'Eg': 3}, [400, 1200])):
        for e in E_:
            out['tpp2m'].append({'mat': mat, 'E': e, 'lam': float(tpp2m(e, **mat))})
    json.dump(out, open(os.path.join(HERE, 'ref.json'), 'w'))
    print('wrote ref.json')


if __name__ == '__main__':
    main()
