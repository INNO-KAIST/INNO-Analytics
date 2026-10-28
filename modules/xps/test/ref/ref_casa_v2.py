#!/usr/bin/env python3
"""Independent reference (numpy/scipy) for the CasaXPS 2.3.27 V2 readings (SPEC_casa_v2.md).
Simulates the CasaXPS fit as pinned from the readings: region Start/End channels i0 = 19, i1 = 361 (BE 280.03 / 297.13),
Linear background with endpoints = mean of y over (i0+1)±1 and (i1−1)±1 (Av. Width 1) anchored at E[i0+1], E[i1−1],
unweighted least squares over channels i0+1 … i1−1, 'RMS' = sqrt(Σ r²), DoF = (i1 − i0 − 1) − p.
Line shapes are written here from their published definitions (no CasaXPS code). Writes test/ref/ref_casa_v2.json.
Run: python3 -I test/ref/ref_casa_v2.py"""
import os, json
import numpy as np
from scipy import integrate, special
from scipy.optimize import least_squares, brentq

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data', 'casaval')
LN2 = np.log(2); S2 = 2 * np.sqrt(2 * LN2); I0, I1 = 19, 361; RHO = 1 / 144.87
E = 279.08 + 0.05 * np.arange(381)


def load(name):
    L = open(os.path.join(D, name + '.vms')).read().split('\n'); k = L.index('762')
    return np.array([float(s) for s in L[k + 5:k + 5 + 762:2]])[::-1]   # KE order → BE ascending


def linbg(y, w=1):
    a, b = I0 + 1, I1 - 1; lo = y[a - w:a + w + 1].mean(); hi = y[b - w:b + w + 1].mean()
    return lo + (hi - lo) * (E - E[a]) / (E[b] - E[a])


def unit_area(f):
    return integrate.quad(f, -np.inf, 0, limit=800)[0] + integrate.quad(f, 0, np.inf, limit=800)[0]


def fwhm_of(f, R=60.0):   # true FWHM of a single-peaked f (max near 0)
    t = np.linspace(-3, 3, 60001); v = f(t); im = int(np.argmax(v)); h = v[im] / 2; x0 = t[im]
    a = brentq(lambda s: f(np.array([s]))[0] - h, x0 - R, x0); b = brentq(lambda s: f(np.array([s]))[0] - h, x0, x0 + R)
    return b - a


def lorpow(t, a, b):   # a: high-BE side (t >= 0)
    l = 1 / (1 + 4 * t * t); return np.where(t >= 0, l ** a, l ** b)


def shape_fn(kind, q):
    """returns g(t) with F = 1 (unnormalised) — width scales as g((x−pos)/F)"""
    if kind == 'LA': return lambda t: lorpow(t, q['a'], q['b'])
    if kind == 'SGL':   # height mix, one F
        h = q['p'] / 100; return lambda t: (1 - h) * np.exp(-4 * LN2 * t * t) + h / (1 + 4 * t * t)
    if kind == 'GLT':   # GL(m)T(k): tail on the high-BE side, length unit F
        m, k = q['m'], q['k']
        def g(t):
            gl = np.exp(-4 * LN2 * (1 - m) * t * t) / (1 + 4 * m * t * t)
            return np.where(t > 0, gl + (1 - gl) * np.exp(-k * np.abs(t)), gl)
        return g
    raise ValueError(kind)


def fit(name, model, p0):
    y = load(name); B = linbg(y); ii = np.arange(I0 + 1, I1); x = E[ii]; d = y[ii]; b = B[ii]
    r = least_squares(lambda p: b + model(x, *p) - d, p0, x_scale='jac', xtol=1e-15, ftol=1e-15, gtol=1e-15)
    res = d - b - model(x, *r.x)
    return r.x, float(np.sqrt((res ** 2).sum())), len(ii) - len(p0)


out = {'source': __doc__.split('\n')[0], 'rho': RHO, 'fits': []}
CASES = [
    ('V2_LA_1_1_F1', 'LA', {'a': 1, 'b': 1}, 'LA(1,1,0)', dict(area=14573.8, fwhm=0.974061, pos=285.0004, rms=864.459, dof=338)),
    ('V2_LA_2_2_F1', 'LA', {'a': 2, 'b': 2}, 'LA(2,2,0)', dict(area=14995.1, fwhm=0.643126, pos=285.0000, rms=18.793, dof=338)),
    ('V2_LA_153_243_F2', 'LA', {'a': 1.53, 'b': 2.43}, 'LA(1.53,2.43,0)', dict(area=14979.3, fwhm=1.15857, pos=284.9999, rms=60.9276, dof=338, asym=0.2275)),
    ('Shape_LA153_243', 'LA', {'a': 1.53, 'b': 2.43}, 'LA(1.53,2.43,0)', dict(area=14997.0, fwhm=0.580056, pos=285.0000, rms=15.6509, dof=338, asym=0.2275)),
    ('V2_SGLh30_F1.2', 'SGL', {'p': 30}, 'SGL(30)', dict(area=14888.6, fwhm=1.19223, pos=285.0002, rms=442.323, dof=338)),
    ('V2_SGLprinted30_F1.2', 'SGL', {'p': 30}, 'SGL(30)', dict(area=13038.1, fwhm=1.4036, pos=285.0008, rms=1404.47, dof=338)),
    ('Shape_GLS30', 'SGL', {'p': 30}, 'SGL(30)', dict(area=15279.4, fwhm=1.17737, pos=285.0001, rms=882.717, dof=338)),
    ('Shape_GLS30', 'SGL', {'p': 19.5}, 'SGL(19.5)', dict(area=14771.0, fwhm=1.20112, pos=285.0001, rms=249.968, dof=338)),
    ('V2_GL30T2_F1.2', 'GLT', {'m': 0.3, 'k': 2}, 'GL(30)T(2)', dict(area=14999.8, fwhm=1.18051, pos=285.0000, rms=6.80617, dof=337, asym=0.2138)),
]
for name, kind, q, casa, rd in CASES:
    g = shape_fn(kind, q); A1 = unit_area(g)
    model = lambda x, pos, area, F, g=g, A1=A1: area * g((x - pos) / F) / (A1 * F)
    p, rms, dof = fit(name, model, [285.01, 14000.0, 1.0 if kind == 'LA' else 1.2])
    F = float(p[2]); out['fits'].append(dict(file=name, kind=kind, q=q, casa=casa, casa_read=rd, pos=float(p[0]), area=float(p[1]), F=F,
                                            fwhm_true=F * fwhm_of(g), rms=rms, dof=dof))
SCAN = {100: (15820.8, 0.787418, 2409.41), 200: (15001.9, 0.832929, 554.538), 250: (14722.9, 0.845688, 376.266), 300: (14505.1, 0.855121, 725.425),
        400: (14190.6, 0.867488, 1349.96), 600: (13824.0, 0.88082, 2092.17), 1000: (13490.2, 0.892121, 2766.1)}
out['scan'] = []
for n, (A, fw, R) in SCAN.items():
    rho = n * RHO
    model = lambda x, pos, area, wL, rho=rho: area * special.voigt_profile(x - pos, rho * wL / S2, wL / 2)
    p, rms, dof = fit('V2_Voigt_wG0.6_wL0.4', model, [285.01, 14000.0, 0.4])
    wL = float(p[2]); f = lambda t, wL=wL, rho=rho: special.voigt_profile(t, rho * wL / S2, wL / 2)
    out['scan'].append(dict(n=n, casa_read=dict(area=A, fwhm=fw, rms=R, pos=285.0001), pos=float(p[0]), area=float(p[1]), F=wL, wG=rho * wL,
                            fwhm_true=fwhm_of(f), rms=rms, dof=dof))
json.dump(out, open(os.path.join(HERE, 'ref_casa_v2.json'), 'w'), indent=1)
for r in out['fits']:
    c = r['casa_read']; print(f"{r['file']:22s} {r['casa']:16s} area {r['area']:9.1f}/{c['area']:9.1f} ({r['area']/c['area']-1:+.1e})  fwhm {r['fwhm_true']:.6f}/{c['fwhm']:.6f}  "
                              f"pos {r['pos']:.4f}/{c['pos']:.4f}  rms {r['rms']:9.3f}/{c['rms']:9.3f}  dof {r['dof']}/{c['dof']}")
for r in out['scan']:
    c = r['casa_read']; print(f"LA(1,1,{r['n']:4d}) area {r['area']:9.1f}/{c['area']:9.1f} ({r['area']/c['area']-1:+.1e})  fwhm {r['fwhm_true']:.6f}/{c['fwhm']:.6f} ({r['fwhm_true']/c['fwhm']-1:+.1e})  rms {r['rms']:8.1f}/{c['rms']:8.1f}")
