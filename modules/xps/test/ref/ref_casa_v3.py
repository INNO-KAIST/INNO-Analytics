#!/usr/bin/env python3
"""Independent reference (numpy/scipy) for the CasaXPS 2.3.27 Round-3 readings (casa_readings_v3.md, SPEC_casa_v3.md).
Simulates the CasaXPS fit (casa3_lib: region 280.03/297.13, Linear bg Av.W 1, unweighted LSQ over i0+1…i1−1, RMS = √Σr²)
with the Cookbook-2026 shape definitions as implemented in INNO (written independently here; no CasaXPS code), and derives
the CasaXPS DoF rule, the 'fwhm' number / asymmetry index of asymmetric shapes, and the DS / ST position offsets.
Writes test/ref/ref_casa_v3.json.  Run: python3 -I test/ref/ref_casa_v3.py   (~1 min)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from casa3_lib import *                                   # noqa: E402,F403
from scipy.optimize import minimize_scalar, brentq        # noqa: E402

READ = {   # Round 3 (+ V2 rows reused for the width/asym analysis): area, fwhm, pos, rms, dof, asym
    'K_GL30n100_F1_2': ('GL(30,100)', 14999.5, 1.19986, 285.0000, 2.77937, 151, 0),
    'K_SGL30n100_F1_2': ('SGL(30,100)', 14896.6, 1.19441, 285.0001, 359.721, 338, 0),
    'K_TLA12_3_0_F1': ('TLA(1.2,3,0)', 14898.2, 0.911717, 285.0004, 290.636, 338, 0.1524),
    'K_TLA12_3_100_F1': ('TLA(1.2,3,100)', 14946.4, 0.9551, 285.0003, 180.058, 338, 0.0817),
    'K_DS01_100_F1': ('DS(0.1,100)', 20811.2, 0.907537, 285.0523, 1570.19, 338, 0.1178),
    'K_ST17_LA11200_F1': ('ST(1.7)LA(1,1,200)', 14908.6, 0.951278, 285.0988, 322.28, 338, 0.0071),
    'K_LA2_3_0_F1': ('LA(2,3,0)', 14996.9, 0.892084, 285.0000, 22.9486, 338, 0.1936),
    'K_LA1_1_250_F1': ('LA(1,1,250)', 14878.1, 0.994132, 285.0001, 417.314, 338, 0.0),
    'Shape_GLP30': ('GL(30)', 14999.8, 1.15343, 285.0000, 8.31557, 156, 0),
    'Shape_LA153_243': ('LA(1.53,2.43,0)', 14997.0, 0.580056, 285.0000, 15.6509, 338, 0.2275),
    'V2_LA_153_243_F2': ('LA(1.53,2.43,0)', 14979.3, 1.15857, 284.9999, 60.9276, 338, 0.2275),
    'V2_GL30T2_F1.2': ('GL(30)T(2)', 14999.8, 1.18051, 285.0000, 6.80617, 337, 0.2138),
}


def glt(m, k):   # GL(m)T(k): tail on the high-BE side, length unit F (SPEC_casa_v2 §2)
    g0 = gl(m); return lambda t: np.where(t > 0, g0(t) + (1 - g0(t)) * np.exp(-k * np.abs(t)), g0(t))


SHAPE = {   # INNO definition of each Casa string (F = Lorentzian/base width; n → Gaussian FWHM n·F/144.87)
    'K_GL30n100_F1_2': lambda: Grid(gl(.3), s=100 * RHO, T=60, dt=0.001),
    'K_SGL30n100_F1_2': lambda: Grid(sgl(.3), s=100 * RHO, T=100, dt=0.001),
    'K_TLA12_3_0_F1': lambda: Grid(tla(1.2, 3), T=100, dt=0.001),
    'K_TLA12_3_100_F1': lambda: Grid(tla(1.2, 3), s=100 * RHO, T=100, dt=0.001),
    'K_DS01_100_F1': lambda: Grid(ds(0.1), s=100 * RHO, T=300, dt=0.0025, A=1.0),   # A = 1: amplitude only (area rule open)
    'K_ST17_LA11200_F1': lambda: Grid(lor, s=200 * RHO, st_mu=1.7, T=100, dt=0.001),
    'K_LA2_3_0_F1': lambda: Grid(la(2, 3), T=100, dt=0.001),
    'K_LA1_1_250_F1': lambda: Grid(lor, s=250 * RHO, T=100, dt=0.001),
    'Shape_GLP30': lambda: Grid(gl(.3), T=60, dt=0.001),
    'Shape_LA153_243': lambda: Grid(la(1.53, 2.43), T=100, dt=0.001),
    'V2_LA_153_243_F2': lambda: Grid(la(1.53, 2.43), T=100, dt=0.001),
    'V2_GL30T2_F1.2': lambda: Grid(glt(.3, 2), T=60, dt=0.001),
}
F0 = {'K_GL30n100_F1_2': 1.0, 'K_SGL30n100_F1_2': 0.95, 'K_TLA12_3_0_F1': 1.14, 'K_TLA12_3_100_F1': 0.81, 'K_DS01_100_F1': 0.64,
      'K_ST17_LA11200_F1': 0.5, 'K_LA2_3_0_F1': 1.73, 'K_LA1_1_250_F1': 0.43, 'Shape_GLP30': 1.2, 'Shape_LA153_243': 1.0,
      'V2_LA_153_243_F2': 2.0, 'V2_GL30T2_F1.2': 1.2}


def peak(g):
    """refined maximum x̂ and half-maximum points a (low BE) < x̂ < b (high BE), t units"""
    t, v = g.t, g.v; i = int(np.argmax(v))
    xm = minimize_scalar(lambda u: -float(g.f(u)), bracket=(t[i - 1], t[i], t[i + 1]), tol=1e-12).x; vm = float(g.f(xm))
    a = brentq(lambda u: g.f(u) - vm / 2, xm - 40, xm); b = brentq(lambda u: g.f(u) - vm / 2, xm, xm + 40)
    return xm, a, b


def casa_dof(name, p=3):
    y = load(name); B = linbg(y); d = (y - B)[II]
    return int(np.count_nonzero(np.abs(d) > 1e-9)) - p


out = {'source': __doc__.split('\n')[0], 'fits': [], 'shape_only': []}
for name, (casa, A, W, P, R, dof, asy) in READ.items():
    g = SHAPE[name](); r = fit_grid(name, g, F0=F0[name]); F = r['F']; xm, a, b = peak(g)
    L, H = xm - a, b - xm
    row = dict(file=name, casa=casa, read=dict(area=A, fwhm=W, pos=P, rms=R, dof=dof, asym=asy),
               pos=r['pos'], area=r['area'], F=F, rms=r['rms'], dof_rule=casa_dof(name),
               true_fwhm=F * (b - a), w_low=2 * F * L, asym_low=1 - L / H, xmax_eV=r['pos'] + F * xm)
    # implied Casa maximum offset δ (t units) from the Casa fwhm number: W = 2F(L + δ)
    dlt = W / (2 * F) - L; row['delta_from_W'] = dlt; row['asym_with_delta'] = 1 - (L + dlt) / (H - dlt) if asy else 0
    if name == 'K_DS01_100_F1':
        al = 0.1; row['ds_shift_pred'] = F * abs(1 / np.tan(np.pi / (2 - al)))
        row['ds_casa_pos_pred'] = r['pos'] + row['ds_shift_pred']
        row['ds_area_over_ampF'] = A / r['area']     # area param here = amplitude × F (shape not normalised)
        row['area'] = None
    out['fits'].append(row)

# shape-only asym index for the two extra readings on V2_GL30T2_F1.2 (GL(30)T(1) 0.4151, GL(0)T(2) 0.1973): δ implied
for lbl, base, asy in [('GL(30)T(1)', glt(.3, 1), 0.4151), ('GL(0)T(2)', glt(0, 2), 0.1973), ('GL(30)T(2)', glt(.3, 2), 0.2138)]:
    g = Grid(base, T=60, dt=0.001); xm, a, b = peak(g); L, H = xm - a, b - xm
    q = 1 - asy; dlt = (q * H - L) / (1 + q)
    out['shape_only'].append(dict(shape=lbl, asym_read=asy, asym_low=1 - L / H, delta_from_asym=dlt, w_low=2 * L, w_with_delta=2 * (L + dlt)))

json.dump(out, open(os.path.join(HERE, 'ref_casa_v3.json'), 'w'), indent=1, default=float)
print(f"{'file':18s} {'casa':18s} {'area sim/Casa':>22s} {'pos sim/Casa':>18s} {'RMS sim/Casa':>20s} DoF rule/Casa | fwhm: true  2(x̂−a)  Casa | asym 1−L/H  Casa")
for r in out['fits']:
    c = r['read']; ar = f"{r['area']:9.1f}/{c['area']:9.1f}" if r['area'] else f"{'(DS)':>9s}/{c['area']:9.1f}"
    print(f"{r['file']:18s} {r['casa']:18s} {ar:>22s} {r['pos']:9.4f}/{c['pos']:.4f} {r['rms']:9.3f}/{c['rms']:9.3f} {r['dof_rule']:4d}/{c['dof']:4d} |"
          f" {r['true_fwhm']:.6f} {r['w_low']:.6f} {c['fwhm']:.6f} | {r['asym_low']:.4f} {c['asym']:.4f}  δ {r['delta_from_W']:+.5f} asymδ {r['asym_with_delta']:.4f}")
    if 'ds_shift_pred' in r:
        print(f"   DS: Casa pos predicted pos + wL·|cot(π/(2−α))| = {r['ds_casa_pos_pred']:.4f} (read {c['pos']}); Casa area / (amplitude·wL) = {r['ds_area_over_ampF']:.5f} (π = {np.pi:.5f})")
for s in out['shape_only']:
    print(f"   {s['shape']:11s} asym read {s['asym_read']:.4f}  1−L/H {s['asym_low']:.4f}  implied δ {s['delta_from_asym']:.5f} F  → Casa fwhm ≈ {s['w_with_delta']:.5f} F")
