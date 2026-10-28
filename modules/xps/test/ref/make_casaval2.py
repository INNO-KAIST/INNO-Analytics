#!/usr/bin/env python3
"""Second set of synthetic (noise-free) spectra for CasaXPS readings (SPEC §13). Same grid and VAMAS layout as test/data/casaval/Shape_GLP30.vms
(BE 279.08–298.08 eV, 0.05 eV, flat 1000 CPS, area 15000, E0 = 285 eV). Writes test/data/casaval/V2_*.vms (+ copies to argv[1] if given).
Run: python3 -I test/ref/make_casaval2.py [/mnt/user-data/outputs/casaval]"""
import os, sys, json, shutil
import numpy as np
from scipy import integrate, special

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data', 'casaval')
LN2 = np.log(2); E = 279.08 + 0.05 * np.arange(381); x = E - 285.0
tmpl = open(os.path.join(D, 'Shape_GLP30.vms')).read().split('\n')
k = tmpl.index('762')   # number of ordinate values (intensity, transmission pairs)
head, tail = tmpl[:k + 5], ['end of experiment', '']


def lorpow(t, F, a, b):
    l = 1 / (1 + 4 * t * t / F / F)
    return np.where(t >= 0, l ** a, l ** b)   # a: high-BE side (SPEC §1.6)


def unit(f):
    A = integrate.quad(f, -np.inf, 0, limit=500)[0] + integrate.quad(f, 0, np.inf, limit=500)[0]
    return lambda t: f(t) / A


cases = {
    'V2_LA_1_1_F1': unit(lambda t: lorpow(t, 1.0, 1, 1)),
    'V2_LA_2_2_F1': unit(lambda t: lorpow(t, 1.0, 2, 2)),
    'V2_LA_153_243_F2': unit(lambda t: lorpow(t, 2.0, 1.53, 2.43)),
    'V2_SGLh30_F1.2': unit(lambda t: 0.7 * np.exp(-4 * LN2 * t * t / 1.44) + 0.3 / (1 + 4 * t * t / 1.44)),
    # SGL exactly as printed in the CasaXPS 2.3.15 manual p.301: (1−m)·G(x:e,f,1) + m·L(x:e,f,m), m = 0.3, f = 1.2
    'V2_SGLprinted30_F1.2': unit(lambda t: 0.7 * np.exp(-4 * LN2 * 0.99 * t * t / 1.44) + 0.3 / (1 + 4 * 0.3 * t * t / 1.44)),
    # CasaXPS GL(30)T(2) per p.301 with the tail on the high-BE side: GL + (1−GL)·exp(−k|x−e|/f), f = 1.2, k = 2
    'V2_GL30T2_F1.2': unit(lambda t: (lambda g: np.where(t > 0, g + (1 - g) * np.exp(-2 * np.abs(t) / 1.2), g))(np.exp(-4 * LN2 * 0.7 * t * t / 1.44) / (1 + 4 * 0.3 * t * t / 1.44))),
}
for wG in (0.3, 0.6, 1.0, 1.5):
    cases[f'V2_Voigt_wG{wG}_wL0.4'] = (lambda wG: (lambda t: special.voigt_profile(t, wG / (2 * np.sqrt(2 * LN2)), 0.2)))(wG)
truth = {}
for name, f in cases.items():
    y = 1000 + 15000 * np.array([float(f(t)) for t in x])
    ke_order = y[::-1]   # first channel = lowest KE = highest BE
    body = []
    for v in ke_order:
        body += [f'{v:.6f}', '1']
    open(os.path.join(D, name + '.vms'), 'w').write('\n'.join(head + body + tail))
    truth[name] = {'E0': 285.0, 'area': 15000, 'bg': 1000}
json.dump(truth, open(os.path.join(D, 'V2_truth.json'), 'w'), indent=1)
if len(sys.argv) > 1:
    for name in list(cases) + ['V2_truth']:
        ext = '.json' if name == 'V2_truth' else '.vms'
        shutil.copy(os.path.join(D, name + ext), sys.argv[1])
print('written', list(cases))
