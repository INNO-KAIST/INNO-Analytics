#!/usr/bin/env python3
"""Independent reference (numpy / scipy / mpmath) for the CasaXPS Cookbook 2026 line shapes (SPEC_shapes2.md).
Written from the Cookbook formulas (pp. 65-68) in the Cookbook's own coordinate x_c (unit FWHM, x_c < 0 = high BE)
and mapped to INNO's x = E - E0 by x_c = -x / F (the mapping is the documented assumption under test).
Methods differ from the JS engine on purpose: adaptive quad for convolutions/normalisation (no closed forms),
FFT convolution + cumulative Simpson on a very fine grid for the ST/SB/DS-Shirley cumulative integrals,
mpmath for one high-precision check of the TLA normalisation. No CasaXPS code is used.
Writes test/ref/ref_shapes2.json.   Run: python3 -I test/ref/ref_shapes2.py   (~2-4 min)"""
import os, json, math
import numpy as np
from scipy import integrate, special
from scipy.signal import fftconvolve
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
LN2 = math.log(2); RHO = 1 / 144.87
C4 = 4 * LN2
QK = dict(limit=800, epsabs=1e-15, epsrel=1e-13)


def lor(t): return 1.0 / (1.0 + 4.0 * t * t)


def gnorm(w): return w * math.sqrt(math.pi / C4)          # integral of exp(-4ln2 (t/w)^2)


def gcdf(z, w):  # CDF of the unit-area Gaussian of FWHM w
    return 0.5 * special.erfc(-2 * math.sqrt(LN2) * z / w)


def quad(f, a, b, pts=None):
    if pts:
        pts = sorted(p for p in pts if a < p < b)
    if math.isinf(a) or math.isinf(b):
        if pts:
            s, lo = 0.0, a
            for p in pts + [b]:
                s += integrate.quad(f, lo, p, **QK)[0]; lo = p
            return s
        return integrate.quad(f, a, b, **QK)[0]
    return integrate.quad(f, a, b, points=pts or None, **QK)[0]


# ---------------- Cookbook profiles (x_c coordinate, unit FWHM) ----------------
def prof(case):
    s, q = case['shape'], case['p']
    if s == 'TLA':
        a, mu = q['alpha'], q['mu']
        return lambda t: lor(t) ** a * (math.pi / 2 - math.atan(2 * t) + math.pi / mu)
    if s in ('GLN', 'SGLN', 'GLP', 'GLS'):
        m = q['m']
        if s in ('GLN', 'GLP'):
            base = lambda t: math.exp(-(1 - m) * C4 * t * t) / (1 + 4 * m * t * t)
        elif s == 'SGLN':
            base = lambda t: m * lor(t) + (1 - m) * math.exp(-C4 * t * t)
        else:   # GLS: area mix of unit-area G and L, both FWHM 1
            base = lambda t: (1 - m) * math.exp(-C4 * t * t) / gnorm(1) + m * lor(t) / (math.pi / 2)
        tail = q.get('tail')
        if tail:
            j = 1.0 if tail == 'T' else q['tj']; kap = q['tk'] if tail == 'T' else 13.8 / q['tk']
            return lambda t: base(t) + ((1 - base(t)) * j * math.exp(kap * t) if t < 0 else 0.0)
        return base
    if s == 'LA':   # lg: alpha for x_c <= 0, beta for x_c > 0
        a, b = q['alpha'], q['beta']
        return lambda t: lor(t) ** (a if t <= 0 else b)
    if s == 'V':     # Lorentzian of FWHM 1 (scale = wL); Gaussian handled as G1
        return lor
    if s == 'DS':
        a = q['alpha']
        return lambda t: math.cos(math.pi * a / 2 + (1 - a) * math.atan(2 * t)) / (1 + 4 * t * t) ** ((1 - a) / 2)
    raise ValueError(s)


def scale(case):   # INNO width that defines the Cookbook unit
    q = case['p']; return q['wL'] if case['shape'] in ('V', 'DS') else q['F']


def g1w(case):     # base Gaussian FWHM in x_c units
    q = case['p']; F = scale(case)
    return (q.get('wG') or 0.0) / F


def kinks(case):
    s = case['shape']; k = [0.0]
    if s == 'DS' and case['p'].get('dsl'): k.append(-case['p']['dsl'])
    return k


# ---------------- quad-based: value of (f (x) G)(x_c) ----------------
def conv_at(f, xc, w, brk):
    if w <= 0: return f(xc)
    a, b = xc - 7.5 * w, xc + 7.5 * w
    return quad(lambda t: f(t) * math.exp(-C4 * ((xc - t) / w) ** 2), a, b, [p for p in brk]) / gnorm(w)


def area_R(f, brk):
    return quad(f, -math.inf, 0.0, [p for p in brk if p < 0]) + quad(f, 0.0, math.inf, [p for p in brk if p > 0])


def win_int(f, w, a, b, brk):   # integral over x_c in [a, b] of (f (x) G)
    if w <= 0: return quad(f, a, b, brk)
    return quad(lambda t: f(t) * (gcdf(b - t, w) - gcdf(a - t, w)), a - 8 * w, b + 8 * w, brk + [a, b])


# ---------------- fine-grid machinery (ST, SB, DS-Shirley) ----------------
H = 1 / 1024.0


def gkern(w):
    K = int(math.ceil(8 * w / H)); t = np.arange(-K, K + 1) * H
    return np.exp(-C4 * (t / w) ** 2) / gnorm(w) * H, K


def gconv_valid(v, w):
    if w <= 0: return v, 0
    k, K = gkern(w); return fftconvolve(v, k, mode='valid'), K


def fine_profile(case, f, L):
    """Returns (t, Fv, Cv, A, tails) on t in [-L, L]: Fv = base (x) G1, Cv = int_t^inf Fv, A = int_R base."""
    w1 = g1w(case); K1 = int(math.ceil(8 * w1 / H)) if w1 > 0 else 0
    n = int(round(L / H)); tt = np.arange(-n - K1, n + K1 + 1) * H
    vf = np.vectorize(f, otypes=[float])
    b = vf(tt)
    Fv, _ = gconv_valid(b, w1); t = np.arange(-n, n + 1) * H
    A = area_R(f, kinks(case))
    # mass of Fv beyond +-L, exactly: int f(s) * P(G1 beyond)
    if w1 > 0:
        mR = quad(lambda s: f(s) * (1 - gcdf(L - s, w1)), L - 9 * w1, math.inf, [L]) + 0.0
        mL = quad(lambda s: f(s) * gcdf(-L - s, w1), -math.inf, -L + 9 * w1, [-L])
    else:
        mR = quad(f, L, math.inf); mL = quad(f, -math.inf, -L)
    cum = integrate.cumulative_simpson(Fv[::-1], dx=H, initial=0.0)[::-1]   # int_t^L
    Cv = cum + mR
    return t, Fv, Cv, A, mL, mR


def st_case(case):
    q = case['p']; f = prof(case); mu = q['stMu']; F = scale(case); w2 = (q.get('stG') or 0.0) / F
    L = 240.0
    t, Fv, Cv, A, mL, mR = fine_profile(case, f, L)
    S = Fv * (Cv + A / mu)
    # integral of S over R: interior + exact tails (int F C = A mL - mL^2/2 on the left, mR^2/2 on the right)
    N = integrate.simpson(S, dx=H) + (A * mL - mL * mL / 2 + A / mu * mL) + (mR * mR / 2 + A / mu * mR)
    if w2 > 0:
        Sv, K2 = gconv_valid(S, w2); t = t[K2:len(t) - K2]
    else:
        Sv = S
    cs = CubicSpline(t, Sv / N)
    return (lambda x: float(cs(-x / F)) / F), cs, F, N, A, (A * A * (0.5 + 1 / mu))


def sb_case(case, lo, hi):
    f = prof(case); F = scale(case)
    t, Fv, Cv, A, mL, mR = fine_profile(case, f, 120.0)
    cs = CubicSpline(t, Cv)
    W = F * cs.integrate(-hi / F, -lo / F)
    return (lambda x: float(cs(-x / F)) / W), F


def dsq_case(case, lo, hi):
    q = case['p']; F = q['wL']; a = q['alpha']; l = q['dsl']; w = g1w(case); f = prof(case)
    R = max(-lo / F, 0) + 12 * w + 40.0
    K = int(math.ceil(8 * w / H)) if w > 0 else 0
    n0 = int(round((max(l, hi / F) + 12 * w + 2) / H)); n1 = int(round(R / H))
    tt = np.arange(-n0 - K, n1 + K + 1) * H
    P = np.vectorize(f, otypes=[float])(tt)
    tailR = quad(f, tt[-1], math.inf)
    Sx = integrate.cumulative_simpson(P[::-1], dx=H, initial=0.0)[::-1] + tailR       # S[P](t) = int_t^inf P
    Pl = f(-l); Sl = quad(f, -l, math.inf, [0.0])                                      # independent check of the grid S at -l
    k = Pl / Sl
    Q = np.where(tt >= -l, P - k * Sx, 0.0)
    if w > 0:
        Qv, _ = gconv_valid(Q, w); t = tt[K:len(tt) - K]
    else:
        Qv, t = Q, tt
    cs = CubicSpline(t, Qv) if w > 0 else None
    if w > 0:
        W = F * cs.integrate(-hi / F, -lo / F); val = lambda x: float(cs(-x / F)) / W
    else:
        S_at = lambda u: quad(f, u, math.inf, [0.0])
        Qf = lambda u: 0.0 if u < -l else f(u) - k * S_at(u)
        W = F * quad(Qf, -hi / F, -lo / F, [-l, 0.0]); val = lambda x: Qf(-x / F) / W
    return val, {'k': k, 'S_l_grid': float(np.interp(-l, tt, Sx)), 'S_l_quad': Sl}


def tp_val(case):
    q = case['p']; F = q['F']; ta, tb = q['ta'], q['tb']; wl = (q.get('wL') or 0.0) / F
    R = lambda t: ta + (tb - ta) * (t + 0.5) if abs(t) <= 0.5 else 0.0
    A = quad(R, -0.5, 0.5)
    if wl <= 0:
        return (lambda x: R(-x / F) / (A * F)), A
    Lz = lambda u: (2 / (math.pi * wl)) / (1 + 4 * u * u / (wl * wl))
    return (lambda x: quad(lambda t: R(t) * Lz(-x / F - t), -0.5, 0.5, [-x / F]) / (A * F)), A


# ---------------- FWHM ----------------
def fwhm(val, x0s, R):
    xs = np.linspace(-R, R, 801); v = np.array([val(x) for x in xs]); im = int(np.argmax(v))
    a, b = xs[max(0, im - 1)], xs[min(len(xs) - 1, im + 1)]
    r = minimize_scalar(lambda x: -val(x), bounds=(a, b), method='bounded', options={'xatol': 1e-12})
    xm, ym = r.x, -r.fun; h = ym / 2
    i = im
    while i > 0 and v[i] >= h: i -= 1
    j = im
    while j < len(xs) - 1 and v[j] >= h: j += 1
    A = brentq(lambda x: val(x) - h, xs[i], xm, xtol=1e-14)
    B = brentq(lambda x: val(x) - h, xm, xs[j], xtol=1e-14)
    return B - A, xm, ym


# ---------------- cases (INNO parameters; 'n' = CasaXPS integer, wG/wL derived with RHO) ----------------
def C(shape, p, **k):
    d = {'shape': shape, 'p': dict(p)}; d.update(k); return d


cases = [
    C('TLA', {'F': 1.0, 'alpha': 1.0, 'mu': 2.0, 'wG': 0.0}),
    C('TLA', {'F': 0.8, 'alpha': 1.5, 'mu': 0.5, 'wG': 0.4}),
    C('TLA', {'F': 1.2, 'alpha': 0.8, 'mu': 5.0, 'wG': 0.9}),
    C('TLA', {'F': 1.0, 'alpha': 2.5, 'mu': 1.0}, n=150),
    C('GLN', {'F': 1.2, 'm': 0.3}, n=50),
    C('GLN', {'F': 1.0, 'm': 0.5, 'wG': 0.6, 'tail': 'T', 'tk': 2.0}),
    C('GLN', {'F': 1.0, 'm': 0.3, 'wG': 0.3, 'tail': 'PHI', 'tj': 1.0, 'tk': 2.0}),
    C('GLN', {'F': 1.0, 'm': 0.9, 'wG': 0.0, 'tail': 'PHI', 'tj': 2.5, 'tk': 5.0}),
    C('GLN', {'F': 1.2, 'm': 0.3, 'wG': 0.0, 'tail': 'T', 'tk': 2.0}),
    C('SGLN', {'F': 1.2, 'm': 0.3}, n=50),
    C('SGLN', {'F': 1.0, 'm': 0.5, 'wG': 0.5, 'tail': 'T', 'tk': 1.5}),
    C('SGLN', {'F': 1.0, 'm': 0.2, 'wG': 0.0, 'tail': 'PHI', 'tj': 0.5, 'tk': 3.0}),
    C('GLS', {'F': 1.0, 'm': 0.3, 'mod': 'ST', 'stMu': 1.7, 'stG': 0.0}),
    C('GLP', {'F': 1.2, 'm': 0.5, 'mod': 'ST', 'stMu': 1.0, 'stG': 0.4}),
    C('LA', {'F': 0.8, 'alpha': 1.0, 'beta': 1.0, 'mod': 'ST', 'stMu': 1.0, 'stG': 0.0}, n=143),
    C('LA', {'F': 1.0, 'alpha': 2.0, 'beta': 3.0, 'wG': 0.3, 'mod': 'ST', 'stMu': 0.5, 'stG': 0.3}),
    C('V', {'wG': 0.6, 'wL': 0.3, 'mod': 'ST', 'stMu': 2.0, 'stG': 0.0}),
    C('TLA', {'F': 1.0, 'alpha': 1.5, 'mu': 2.0, 'wG': 0.3, 'mod': 'ST', 'stMu': 1.0, 'stG': 0.2}),
    C('GLN', {'F': 1.0, 'm': 0.3, 'tail': 'T', 'tk': 2.0, 'mod': 'ST', 'stMu': 1.0, 'stG': 0.0}, n=50),
    C('LA', {'F': 1.0, 'alpha': 1.0, 'beta': 1.0, 'mod': 'SB'}, n=43, win=[-12.0, 8.0]),
    C('GLS', {'F': 1.2, 'm': 0.3, 'mod': 'SB'}, win=[-6.0, 10.0]),
    C('V', {'wG': 0.5, 'wL': 0.4, 'mod': 'SB'}, win=[-5.0, 15.0]),
    C('DS', {'alpha': 0.1, 'wL': 0.4}, n=50, win=[-8.0, 12.0]),
    C('DS', {'alpha': 0.3, 'wL': 0.3, 'wG': 0.6}, win=[-6.0, 14.0]),
    C('DS', {'alpha': 0.2, 'wL': 0.5, 'wG': 0.0, 'dsl': 40.0}, win=[-10.0, 15.0]),
    C('DS', {'alpha': 0.3, 'wL': 0.2, 'dsl': 40.0}, n=600, win=[-8.0, 12.0]),
    C('DS', {'alpha': 0.25, 'wL': 0.3, 'wG': 0.5, 'dsl': 10.0}, win=[-8.0, 12.0]),
    C('TP', {'F': 3.0, 'ta': 1.0, 'tb': 3.0, 'wL': 0.5}),
    C('TP', {'F': 10.0, 'ta': 1.0, 'tb': 1.0, 'wL': 1.0}),
    C('TP', {'F': 4.0, 'ta': 1.0, 'tb': 15.0, 'wL': 0.0}),
]

out = []
for ci, c in enumerate(cases):
    q = c['p']
    if 'n' in c:   # CasaXPS integer -> width with the LA rule (assumed for these shapes)
        if c['shape'] == 'TP': q['wL'] = c['n'] * RHO * q['F']
        elif c['shape'] == 'DS': q['wG'] = c['n'] * RHO * q['wL']
        else: q['wG'] = c['n'] * RHO * q['F']
    s = c['shape']; F = scale(c) if s != 'TP' else q['F']
    rec = {'shape': s, 'p': q, 'n': c.get('n'), 'win': c.get('win')}
    W = None
    if s == 'TP':
        val, A = tp_val(c); rec['A_c'] = A
    elif q.get('mod') == 'ST':
        val, cs, F, N, A, Nan = st_case(c); rec['N_numeric'] = N; rec['N_identity'] = Nan
    elif q.get('mod') == 'SB':
        val, F = sb_case(c, *c['win'])
    elif s == 'DS' and q.get('dsl'):
        val, info = dsq_case(c, *c['win']); rec['dsq'] = info
    else:
        f = prof(c); w = g1w(c); brk = kinks(c)
        if c.get('win'):
            lo, hi = c['win']; W = F * win_int(f, w, -hi / F, -lo / F, brk)
            val = (lambda f, w, brk, W, F: (lambda x: conv_at(f, -x / F, w, brk) / W))(f, w, brk, W, F)
        else:
            A = area_R(f, brk); rec['A_c'] = A
            val = (lambda f, w, brk, A, F: (lambda x: conv_at(f, -x / F, w, brk) / (A * F)))(f, w, brk, A, F)
            # area inside [-R, R] (INNO x), exactly: (1/A) int f(s) [Phi_g(Rc - s) - Phi_g(-Rc - s)] ds
            Rc = 40.0; rec['R_in'] = Rc * F
            rec['area_in'] = (win_int(f, w, -Rc, Rc, brk)) / A
    # sample points: dense around the peak, some far out (INNO x, eV)
    wid = (q.get('F') or q.get('wL') or 1.0) + 2 * (q.get('wG') or 0) + (q.get('stG') or 0)
    xs = list(np.linspace(-5 * wid, 5 * wid, 47) + 0.0137 * wid)
    xs += [-12 * wid + 0.031, 9.7 * wid, 15.3 * wid]
    if s == 'TP' and not q.get('wL'):
        xs = [x for x in xs if abs(abs(x) - q['F'] / 2) > 1e-3]
    if c.get('win'):
        lo, hi = c['win']; xs = [x for x in xs if lo <= x <= hi]
    rec['x'] = [float(x) for x in xs]; rec['v'] = [float(val(x)) for x in xs]
    if q.get('mod') == 'ST':   # ST area inside [-R, R] from the spline (INNO x)
        Rc = 40.0; rec['R_in'] = Rc * F; rec['area_in'] = float(cs.integrate(-Rc, Rc))
    if s == 'TP':
        Rc = 40.0 if q.get('wL') else 0.6; rec['R_in'] = Rc * q['F']
        rec['area_in'] = quad(lambda x: val(x), -Rc * q['F'], Rc * q['F'], [-q['F'] / 2, q['F'] / 2]) if q.get('wL') else 1.0
    if q.get('mod') != 'SB' and not (s == 'TP' and not q.get('wL')):
        fw, xm, ym = fwhm(val, None, 8 * wid); rec['fwhm'] = fw; rec['xmax'] = xm; rec['ymax'] = ym
    out.append(rec)
    print(ci, s, q.get('mod') or '', json.dumps({k: v for k, v in q.items()}), 'fwhm', rec.get('fwhm'), 'area_in', rec.get('area_in'), flush=True)

# one high-precision check of the TLA normalisation (mpmath, 30 digits) — independent of the JS closed form
mp.mp.dps = 30
tla_mp = []
for a, mu in [(1.0, 2.0), (1.5, 0.5), (0.8, 5.0), (2.5, 1.0), (0.55, 3.0)]:
    # A: symmetric part via the Beta function; x_c > 0 half: int_0^inf l^a atan(1/(2t)) dt (decays like t^(-2a-1)) + (pi/mu) B/4
    a_ = mp.mpf(a); B = mp.beta(mp.mpf(1) / 2, a_ - mp.mpf(1) / 2) / 2          # int_R (1+4t^2)^-a dt
    A = (mp.pi / 2 + mp.pi / mu) * B
    Lft = mp.quad(lambda t: (1 / (1 + 4 * t * t)) ** a_ * mp.atan(1 / (2 * t)), [0, 1, 10, 100, 1e4, mp.inf]) + (mp.pi / mu) * B / 2
    tla_mp.append({'alpha': a, 'mu': mu, 'A': float(A), 'left_frac': float(Lft / A)})
json.dump({'spec': 'SPEC_shapes2.md', 'rho': RHO, 'cases': out, 'tla_mp': tla_mp}, open(os.path.join(HERE, 'ref_shapes2.json'), 'w'))
print('wrote ref_shapes2.json', len(out), 'cases')
