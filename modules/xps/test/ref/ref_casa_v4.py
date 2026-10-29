#!/usr/bin/env python3
"""Independent reference (numpy/scipy) for the CasaXPS 2.3.27 Round-4 readings (casa_readings_v3.md 'Round 4', SPEC_casa_v4.md).
Re-uses casa3_lib (CasaXPS fit simulator: Linear bg Av.W 1, unweighted LSQ over channels i0+1…i1−1, RMS = √Σr²; Cookbook-2026
shapes written independently, no CasaXPS code).  Settles / tests:
  A. DS(α,n): area normalisation  Casa Area = amp·wL·∫_{−600}^{+600} P(t,α) dt   (P = Cookbook peak form, t in wL units)
              position            Casa pos  = INNO origin + wL·|cot(π/(2−α))|
              window dependence   (refit of K_DS01_100_F1 on 280.03–290.03; channel i1 = 220 reproduces DoF 197)
  B. width / asymmetry number for every read shape:  x̂ = vertex of the parabola through the shape at the grid nodes
     k·d (d = F/40, nodes aligned to the base origin) around the largest node;  half level = vertex value / 2;
     l = x̂ − a, h = b − x̂  (a < x̂ < b half-maximum points; a on the low-BE side)
     → fwhm = 2·min(l, h),  asym = 1 − l/h   (sign: negative when the low-BE side is wider)
  C. ST(μ): position candidates (none fits), DoF rule.
Writes test/ref/ref_casa_v4.json.  Run: python3 -I test/ref/ref_casa_v4.py   (~2 min)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from casa3_lib import *                                   # noqa: E402,F403
from scipy.optimize import minimize_scalar, brentq        # noqa: E402

GRID_D = 1 / 40          # CasaXPS maximum search: node spacing in units of the base width F (DS: wL)
DS_T = 600.0             # CasaXPS DS normalisation range ±600·wL

# file: (Casa string, area, fwhm, pos, rms, dof, asym, round)
READ = {
    'K_DS0_100_F1': ('DS(0,100)', 15180.3, 0.986552, 285.0002, 693.123, 338, 0.0, 4),
    'K_DS02_100_F1': ('DS(0.2,100)', 27119.4, 0.795054, 285.0927, 2522.3, 338, 0.2550, 4),
    'K_DS01_0_F1': ('DS(0.1,0)', 21288.6, 0.868703, 285.0689, 1649.87, 338, 0.1603, 4),
    'K_DS01_100_F1': ('DS(0.1,100)', 20811.2, 0.907537, 285.0523, 1570.19, 338, 0.1178, 3),
    'K_ST17_LA110_F1': ('ST(1.7)LA(1,1,0)', 14734.8, 0.891263, 285.0737, 596.667, 338, 0.1529, 4),
    'K_ST3_LA11200_F1': ('ST(3)LA(1,1,200)', 14928.0, 0.911423, 285.1231, 279.251, 338, 0.0440, 4),
    'K_ST17_LA11200_F1': ('ST(1.7)LA(1,1,200)', 14908.6, 0.951278, 285.0988, 322.28, 338, 0.0071, 3),
    'K_ST17_GL30_F1_2': ('ST(1.7)GL(30)', 14998.0, 1.11634, 285.1485, 37.373, 161, -0.0131, 4),
    'K_LA243_153_0_F1': ('LA(2.43,1.53,0)', 14904.2, 0.866617, 285.0002, 320.787, 338, -0.2944, 4),
    'K_LA2_3_30_F1': ('LA(2,3,30)', 14996.5, 0.930135, 285.0000, 19.4075, 338, 0.1296, 4),
    'K_LA2_3_0_F1': ('LA(2,3,0)', 14996.9, 0.892084, 285.0000, 22.9486, 338, 0.1936, 3),
    'K_TLA12_3_0_F1': ('TLA(1.2,3,0)', 14898.2, 0.911717, 285.0004, 290.636, 338, 0.1524, 3),
    'K_TLA12_3_100_F1': ('TLA(1.2,3,100)', 14946.4, 0.9551, 285.0003, 180.058, 338, 0.0817, 3),
    'Shape_LA153_243': ('LA(1.53,2.43,0)', 14997.0, 0.580056, 285.0000, 15.6509, 338, 0.2275, 1),
    'V2_LA_153_243_F2': ('LA(1.53,2.43,0)', 14979.3, 1.15857, 284.9999, 60.9276, 338, 0.2275, 2),
    'V2_GL30T2_F1.2': ('GL(30)T(2)', 14999.8, 1.18051, 285.0000, 6.80617, 337, 0.2138, 2),
}
WINDOW = ('K_DS01_100_F1', 'DS(0.1,100) region 280.03–290.03', 19946.5, 0.876718, 285.0489, 2202.23, 197, 0.1178)


def glt(m, k):
    g0 = gl(m); return lambda t: np.where(t > 0, g0(t) + (1 - g0(t)) * np.exp(-k * np.abs(t)), g0(t))


SHAPE = {   # (grid factory, kind, alpha for DS, F0)
    'K_DS0_100_F1': (lambda: Grid(ds(0.0), s=100 * RHO, T=300, dt=0.0025, A=1.0), 'DS', 0.0, 0.73),
    'K_DS02_100_F1': (lambda: Grid(ds(0.2), s=100 * RHO, T=300, dt=0.0025, A=1.0), 'DS', 0.2, 0.57),
    'K_DS01_0_F1': (lambda: Grid(ds(0.1), T=300, dt=0.0025, A=1.0), 'DS', 0.1, 0.88),
    'K_DS01_100_F1': (lambda: Grid(ds(0.1), s=100 * RHO, T=300, dt=0.0025, A=1.0), 'DS', 0.1, 0.64),
    'K_ST17_LA110_F1': (lambda: Grid(lor, st_mu=1.7, T=200, dt=0.001), 'ST', None, 1.0),
    'K_ST3_LA11200_F1': (lambda: Grid(lor, s=200 * RHO, st_mu=3, T=200, dt=0.001), 'ST', None, 0.5),
    'K_ST17_LA11200_F1': (lambda: Grid(lor, s=200 * RHO, st_mu=1.7, T=200, dt=0.001), 'ST', None, 0.5),
    'K_ST17_GL30_F1_2': (lambda: Grid(gl(.3), st_mu=1.7, T=60, dt=0.001), 'ST', None, 1.25),
    'K_LA243_153_0_F1': (lambda: Grid(la(2.43, 1.53), T=100, dt=0.001), 'LA', None, 1.5),
    'K_LA2_3_30_F1': (lambda: Grid(la(2, 3), s=30 * RHO, T=100, dt=0.001), 'LA', None, 1.59),
    'K_LA2_3_0_F1': (lambda: Grid(la(2, 3), T=100, dt=0.001), 'LA', None, 1.73),
    'K_TLA12_3_0_F1': (lambda: Grid(tla(1.2, 3), T=100, dt=0.001), 'TLA', None, 1.14),
    'K_TLA12_3_100_F1': (lambda: Grid(tla(1.2, 3), s=100 * RHO, T=100, dt=0.001), 'TLA', None, 0.81),
    'Shape_LA153_243': (lambda: Grid(la(1.53, 2.43), T=100, dt=0.001), 'LA', None, 1.0),
    'V2_LA_153_243_F2': (lambda: Grid(la(1.53, 2.43), T=100, dt=0.001), 'LA', None, 2.0),
    'V2_GL30T2_F1.2': (lambda: Grid(glt(.3, 2), T=60, dt=0.001), 'GLT', None, 1.2),
}


def P(t, al): return np.cos(np.pi * al / 2 + (1 - al) * np.arctan(-2 * t)) / (1 + 4 * t * t) ** ((1 - al) / 2)


def ds_int(al, a, b):   # ∫_a^b P(t,α) dt  (t in wL units)
    pts = [x for x in (-50, -5, 0, 5, 50) if a < x < b]; s = [a] + pts + [b]
    return sum(integrate.quad(lambda u: P(u, al), s[i], s[i + 1], limit=4000, epsrel=1e-13)[0] for i in range(len(s) - 1))


def true_peak(g):
    t, v = g.t, g.v; i = int(np.argmax(v))
    xm = minimize_scalar(lambda u: -float(g.f(u)), bracket=(t[i - 1], t[i], t[i + 1]), tol=1e-12).x; vm = float(g.f(xm))
    a = brentq(lambda u: g.f(u) - vm / 2, xm - 40, xm); b = brentq(lambda u: g.f(u) - vm / 2, xm, xm + 40)
    return xm, a, b


def casa_width(g, d=GRID_D):
    """rule B: parabola vertex on nodes k·d around the largest node; half = vertex value / 2; fwhm = 2 min(l,h), asym = 1 − l/h (t units)"""
    xm = true_peak(g)[0]; k = int(round(xm / d)); f = lambda u: float(g.f(u))
    while f((k + 1) * d) > f(k * d): k += 1
    while f((k - 1) * d) > f(k * d): k -= 1
    ym, y0, yp = f((k - 1) * d), f(k * d), f((k + 1) * d); den = ym - 2 * y0 + yp
    xs = k * d + 0.5 * d * (ym - yp) / den; ys = y0 - 0.125 * (ym - yp) ** 2 / den
    a = brentq(lambda u: f(u) - ys / 2, xs - 40, xs); b = brentq(lambda u: f(u) - ys / 2, xs, xs + 40)
    l, h = xs - a, b - xs
    return dict(xhat=xs, lo=a, hi=b, fwhm=2 * min(l, h), asym=1 - l / h, side='low' if l <= h else 'high')


def fit_window(name, g, i0, i1, p0):
    y = load(name); a, b = i0 + 1, i1 - 1
    lo = y[a - 1:a + 2].mean(); hi = y[b - 1:b + 2].mean(); B = lo + (hi - lo) * (E - E[a]) / (E[b] - E[a])
    ch = np.arange(a, b + 1); x = E[ch]; d = y[ch] - B[ch]
    r = least_squares(lambda p: p[1] * g.unit(x, p[0], p[2]) - d, p0, x_scale='jac', xtol=1e-14, ftol=1e-14, gtol=1e-14)
    res = d - r.x[1] * g.unit(x, r.x[0], r.x[2])
    return dict(pos=r.x[0], area=r.x[1], F=r.x[2], rms=float(np.sqrt((res ** 2).sum())), dof=int(np.count_nonzero(np.abs(d) > 1e-9)) - 3)


out = {'source': __doc__.split('\n')[0], 'rules': dict(grid_d=GRID_D, ds_T=DS_T), 'fits': [], 'ds': [], 'st': []}
print(f"{'file':18s} {'casa':19s} {'area sim/Casa':>21s} {'pos sim/Casa':>19s} {'RMS sim/Casa':>19s} DoF | fwhm rule/Casa (rel) | asym rule/Casa | side")
for name, (casa, A, W, Pp, R, dof, asy, rnd) in READ.items():
    mk, kind, al, F0 = SHAPE[name]; g = mk(); r = fit_grid(name, g, F0=F0); F = r['F']
    y = load(name); B = linbg(y); dof_rule = int(np.count_nonzero(np.abs((y - B)[II]) > 1e-9)) - 3
    cw = casa_width(g); xm, a, b = true_peak(g)
    shift = F * abs(1 / np.tan(np.pi / (2 - al))) if kind == 'DS' and al > 0 else 0.0
    row = dict(file=name, casa=casa, round=rnd, kind=kind, read=dict(area=A, fwhm=W, pos=Pp, rms=R, dof=dof, asym=asy),
               pos=r['pos'], F=F, ampF=r['area'], rms=r['rms'], dof_rule=dof_rule,
               fwhm_rule=F * cw['fwhm'], asym_rule=cw['asym'], side=cw['side'], xhat_minus_max=cw['xhat'] - xm,
               true_fwhm=F * (b - a), pos_pred=r['pos'] + shift)
    if kind == 'DS':
        I = ds_int(al, -DS_T, DS_T); row.update(K_read=A / r['area'], K_rule=I, area_pred=r['area'] * I, posShift_rule=shift, posShift_read=Pp - r['pos'])
        out['ds'].append(dict(file=name, alpha=al, n=casa, K_read=A / r['area'], K_rule=I, ratio=A / (r['area'] * I),
                              shift_read=Pp - r['pos'], shift_rule=shift, wL=F))
    if kind == 'ST':
        off = (Pp - r['pos']) / F; t, v = g.t, g.v; dt = t[1] - t[0]; c = np.cumsum(v) * dt; c /= c[-1]
        out['st'].append(dict(file=name, casa=casa, F=F, off_read_F=off, off_read_eV=Pp - r['pos'], max_F=xm, halfmid_F=(a + b) / 2,
                              median_F=float(t[np.searchsorted(c, 0.5)]), d_max_eV=(Pp - r['pos']) - F * xm))
    out['fits'].append(row)
    ap = row.get('area_pred', r['area'])
    print(f"{name:18s} {casa:19s} {ap:9.1f}/{A:9.1f} {row['pos_pred']:9.5f}/{Pp:.4f} {r['rms']:9.2f}/{R:8.2f} {dof_rule:3d}/{dof:3d} |"
          f" {row['fwhm_rule']:.6f}/{W:.6f} ({row['fwhm_rule'] / W - 1:+.1e}) | {cw['asym']:+.4f}/{asy:+.4f} | {cw['side']}")

# window refit (DS(0.1,100), region End 290.03): CasaXPS DoF 197 = channel i1 = 220 (BE 290.08); i1 = 219 gives 196
nm, lab, A, W, Pp, R, dof, asy = WINDOW; g = SHAPE[nm][0](); al = 0.1
for i1 in (219, 220):
    r = fit_window(nm, g, 19, i1, [285.003, 6500, 0.65]); I = ds_int(al, -DS_T, DS_T); xm, a, b = true_peak(g)
    row = dict(file=nm, casa=lab, i0=19, i1=i1, pos=r['pos'], F=r['F'], ampF=r['area'], rms=r['rms'], dof_rule=r['dof'],
               read=dict(area=A, fwhm=W, pos=Pp, rms=R, dof=dof, asym=asy), K_read=A / r['area'], K_rule=I, area_pred=r['area'] * I,
               pos_pred=r['pos'] + r['F'] * abs(1 / np.tan(np.pi / (2 - al))), fwhm_rule=r['F'] * casa_width(g)['fwhm'])
    out['fits'].append(row)
    print(f"window i1={i1} (E {E[i1]:.2f}): area {row['area_pred']:.1f}/{A}  pos {row['pos_pred']:.5f}/{Pp}  RMS {r['rms']:.2f}/{R}  DoF {r['dof']}/{dof}  "
          f"K {row['K_read']:.5f} (default window {next(d['K_read'] for d in out['ds'] if d['file']==nm):.5f})  fwhm {row['fwhm_rule']:.6f}/{W}")

print('\nDS: Casa Area/(amp·wL) vs ∫_{±600}P;  position shift read vs wL·|cot(π/(2−α))|')
for d in out['ds']:
    print(f"  α {d['alpha']:.1f} {d['n']:12s} K read {d['K_read']:.5f} rule {d['K_rule']:.5f} ratio {d['ratio']:.5f} | shift read {d['shift_read']:.5f} rule {d['shift_rule']:.5f} (Δ {d['shift_read'] - d['shift_rule']:+.5f} eV)")
print('\nST position (F units): read offset vs candidate features')
for s in out['st']:
    print(f"  {s['casa']:19s} off {s['off_read_F']:.5f} | max {s['max_F']:.5f} half-mid {s['halfmid_F']:.5f} median {s['median_F']:.5f} | pos − max {s['d_max_eV']:+.4f} eV")
json.dump(out, open(os.path.join(HERE, 'ref_casa_v4.json'), 'w'), indent=1, default=float)

# ======================= Round 5 (casa_readings_v3.md 'Round 5', truth K3_truth.json) =======================
# ST: CasaXPS step = exclusive rectangle cumulative on the F/40 grid ≈ Φ(t − d/2) (step delayed by half a cell towards high BE).
ST_SHIFT = GRID_D / 2


class STG:
    """S(μ)F with Φ(t − sh): base (optionally ⊗G, FWHM s) on a fine grid, cumulative by trapezoid then shifted; normalised to its own numerical area
    (the d/2 delay changes ∫S by −sh·∫F² — CasaXPS's areas follow the numerical area)."""
    def __init__(self, base, mu, sh=ST_SHIFT, s=0.0, T=60.0, dt=0.0005):
        b = Grid(base, s=s, T=T, dt=dt); t, v = b.t, b.v; A = b.A
        left = integrate.quad(lambda u: float(base(np.array([u]))[0]), -np.inf, -T, limit=500)[0]
        cum = left + np.concatenate([[0], np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)])
        Ph = np.interp(t - sh, t, cum, left=0.0); S = v * (Ph + A / mu)
        self.t, self.v, self.A = t, S, float(np.trapezoid(S, t)); self.f = interpolate.CubicSpline(t, S)

    def unit(self, x, pos, F): return self.f((x - pos) / F) / (self.A * F)


READ5 = {   # file: (string, area, fwhm, pos, rms, dof, asym, base, n, mu, F0, pos0)
    'K_ST1_GL0_F1_2': ('ST(1)GL(0)', 14998.2, 1.16736, 285.1241, 28.8428, 144, -0.0560, gl(0), 0, 1, 1.25),
    'K_ST17_GL0_F1_2': ('ST(1.7)GL(0)', 14999.1, 1.1715, 285.1595, 41.0479, 148, -0.0483, gl(0), 0, 1.7, 1.29),
    'K_ST3_GL0_F1_2': ('ST(3)GL(0)', 15000.0, 1.18428, 285.2315, 50.1108, 153, -0.0256, gl(0), 0, 3, 1.34),
    'K_ST10_GL0_F1_2': ('ST(10)GL(0)', 15002.9, 1.1855, 285.2791, 41.3143, 161, 0.0226, gl(0), 0, 10, 1.42),
    'K_ST17_GL30_F0_6': ('ST(1.7)GL(30)', 14999.0, 0.59604, 285.0793, 49.9146, 86, -0.0131, gl(.3), 0, 1.7, 0.67),
    'K_ST17_GL30_F1_2': ('ST(1.7)GL(30)', 14998.0, 1.11634, 285.1485, 37.373, 161, -0.0131, gl(.3), 0, 1.7, 1.25),
    'K_ST17_GL30_F2_4': ('ST(1.7)GL(30)', 14995.5, 2.38363, 285.3171, 26.4777, 337, -0.0131, gl(.3), 0, 1.7, 2.67),
    'K_ST17_LA110_F1': ('ST(1.7)LA(1,1,0)', 14734.8, 0.891263, 285.0737, 596.667, 338, 0.1529, lor, 0, 1.7, 1.0),
    'K_ST17_LA11200_F1': ('ST(1.7)LA(1,1,200)', 14908.6, 0.951278, 285.0988, 322.28, 338, 0.0071, lor, 200, 1.7, 0.5),
    'K_ST3_LA11200_F1': ('ST(3)LA(1,1,200)', 14928.0, 0.911423, 285.1231, 279.251, 338, 0.0440, lor, 200, 3, 0.5),
}
out['st5'] = []
print('\nRound 5 / all ST: Cookbook S(μ)F with the step delayed by d/2 = F/80 (sh), vs sh = 0')
for name, (casa, A, W, Pp, R, dof, asy, base, n, mu, F0) in READ5.items():
    for sh in (0.0, ST_SHIFT):
        T = 200.0 if base is lor else 40.0
        g = STG(base, mu, sh=sh, s=n * RHO, T=T, dt=0.001 if base is lor else 0.0005)
        r = fit_grid(name, g, F0=F0); F = r['F']; cw = casa_width(g); xm = true_peak(g)[0]
        y = load(name); B = linbg(y); dr = int(np.count_nonzero(np.abs((y - B)[II]) > 1e-9)) - 3
        row = dict(file=name, casa=casa, sh=sh, read=dict(area=A, fwhm=W, pos=Pp, rms=R, dof=dof, asym=asy), pos=r['pos'], F=F, area=r['area'], rms=r['rms'],
                   dof_rule=dr, fwhm_rule=F * cw['fwhm'], asym_rule=cw['asym'], xhat_F=cw['xhat'], off_F=(Pp - r['pos']) / F, off_minus_xhat_F=(Pp - r['pos']) / F - cw['xhat'])
        out['st5'].append(row)
        print(f"  {name:18s} sh {sh:.4f} area {r['area']:8.1f}/{A} RMS {r['rms']:7.2f}/{R} DoF {dr}/{dof} fwhm {row['fwhm_rule']:.6f}/{W} ({row['fwhm_rule'] / W - 1:+.1e})"
              f" asym {cw['asym']:+.4f}/{asy:+.4f} | off {row['off_F']:.5f}F x̂ {cw['xhat']:.5f}F Δ {row['off_minus_xhat_F']:+.4f}F")

# DS(0.3,0): truth origin 284.7508 (generated with the v3 casaToComp shift)
g = Grid(ds(0.3), T=300, dt=0.001, A=1.0); r = fit_grid('K_DS03_0_F1', g, F0=0.88, pos0=284.75); F = r['F']; I = ds_int(0.3, -DS_T, DS_T); cw = casa_width(g)
sh = F * abs(1 / np.tan(np.pi / 1.7))
out['ds5'] = dict(file='K_DS03_0_F1', read=dict(area=37178.9, fwhm=0.814936, pos=284.9485, rms=3185.66, dof=338, asym=0.4834), pos=r['pos'], wL=F, ampF=r['area'],
                  rms=r['rms'], pos_pred=r['pos'] + sh, K_read=37178.9 / r['area'], K_rule=I, fwhm_rule=F * cw['fwhm'], asym_rule=cw['asym'])
print(f"\nDS(0.3,0): pos {r['pos'] + sh:.5f}/284.9485 (INNO origin {r['pos']:.5f} + wL|cot| {sh:.5f}); K {37178.9 / r['area']:.5f} rule {I:.5f} ratio {37178.9 / r['area'] / I:.5f};"
      f" RMS {r['rms']:.2f}/3185.66; fwhm {F * cw['fwhm']:.6f}/0.814936; asym {cw['asym']:.4f}/0.4834")
json.dump(out, open(os.path.join(HERE, 'ref_casa_v4.json'), 'w'), indent=1, default=float)

# ======================= Round 6 + position rule (all 16 ST readings) =======================
# CasaXPS tabulates S(μ)F on nodes k·D, D = (base FWHM)/40, step lagged by D/2, and reports pos = base origin + D·argmax_k S(kD).
READ5.update({
    'K_ST13_GL0_F1_2': ('ST(1.3)GL(0)', 14998.5, 1.16827, 285.1573, 34.8693, 146, -0.0543, gl(0), 0, 1.3, 1.27),
    'K_ST2_GL0_F1_2': ('ST(2)GL(0)', 14999.5, 1.17442, 285.1935, 44.0331, 149, -0.0429, gl(0), 0, 2, 1.30),
    'K_ST25_GL0_F1_2': ('ST(2.5)GL(0)', 14999.8, 1.17948, 285.1960, 48.4436, 151, -0.0339, gl(0), 0, 2.5, 1.32),
    'K_ST4_GL0_F1_2': ('ST(4)GL(0)', 15001.7, 1.19204, 285.2351, 50.5211, 155, -0.0121, gl(0), 0, 4, 1.36),
    'K_ST6_GL0_F1_2': ('ST(6)GL(0)', 15002.6, 1.19598, 285.2742, 48.7939, 158, 0.0053, gl(0), 0, 6, 1.39),
})
out['st6'] = []
print('\nST position rule: pos = origin + D·(largest node), D = base FWHM/40; shape lag D/2')
for name, (casa, A, W, Pp, R, dof, asy, base, n, mu, F0) in READ5.items():
    lor_b = base is lor; T = 200.0 if lor_b else 40.0; dt = 0.001 if lor_b else 0.0005
    bg = Grid(base, s=n * RHO, T=T, dt=dt); xb, ab, bb = true_peak(bg); K = bb - ab; D = K / 40
    g = STG(base, mu, sh=D / 2, s=n * RHO, T=T, dt=dt); r = fit_grid(name, g, F0=F0); F = r['F']
    f = lambda u: float(g.f(u)); k = int(round(true_peak(g)[0] / D))
    while f((k + 1) * D) > f(k * D): k += 1
    while f((k - 1) * D) > f(k * D): k -= 1
    cw = casa_width(g, d=D)
    row = dict(file=name, casa=casa, K=K, D=D, k=k, pos=r['pos'], F=F, area=r['area'], rms=r['rms'], pos_pred=r['pos'] + F * k * D,
               fwhm_rule=F * cw['fwhm'], asym_rule=cw['asym'], read=dict(area=A, fwhm=W, pos=Pp, rms=R, dof=dof, asym=asy))
    out['st6'].append(row)
    print(f"  {name:18s} K {K:.5f} k {k:2d} | pos {row['pos_pred']:.5f}/{Pp} ({(row['pos_pred'] - Pp) * 1e3:+.2f} meV) | area {r['area']:.1f}/{A} RMS {r['rms']:.1f}/{R}"
          f" | fwhm {row['fwhm_rule'] / W - 1:+.1e} asym {cw['asym']:+.4f}/{asy:+.4f}")
json.dump(out, open(os.path.join(HERE, 'ref_casa_v4.json'), 'w'), indent=1, default=float)
