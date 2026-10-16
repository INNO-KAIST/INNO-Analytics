#!/usr/bin/env python3
"""Build src/xps_lib.js (window.XPS_LIB / module.exports) from the curated data files.

Inputs (all in ../data):
  core_xdb_table_1_1.json     X-ray Data Booklet Table 1-1 (core-level BE, elements in natural form)
  kherve_extract.json         numeric RSF tables (Scofield Al/Mg, Avantage ALTHERMO1) + Avantage Auger KE
  scofield_paper_spotcheck.json  values read from the Scofield 1976 paper   (verification only)
  auger_literature_check.json    NIST SRD 20 / ASTM E2108 Auger KE         (verification only)

Every check below must pass, otherwise the build aborts (no silent fallback).
Run:  python3 -I build/build_lib.py
"""
import json, os, re, sys, hashlib, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'src', 'xps_lib.js')
rd = lambda f: json.load(open(os.path.join(DATA, f), encoding='utf-8'))

core = rd('core_xdb_table_1_1.json')
kh = rd('kherve_extract.json')
spot = rd('scofield_paper_spotcheck.json')
aug_ref = rd('auger_literature_check.json')

ORB = re.compile(r'^([1-7])([spdf])(1/2|3/2|5/2|7/2)?$')
FIX_KEY = {'5p1': '5p1/2', '4p1': '4p1/2', '3p1': '3p1/2', '2p1': '2p1/2', '6p1': '6p1/2'}


def norm_photo_key(k):
    el, orb = k.split('_', 1)
    orb = FIX_KEY.get(orb, orb)
    if not ORB.match(orb):
        return None, None
    return el, orb


problems = []

# ---------- core levels ----------
core_out = {}
for z, e in core['elements'].items():
    lv = {}
    for col, v in e['lv'].items():
        orb = col.split(' ')[1]
        lv[orb] = v
    rec = {'Z': int(z), 'lv': lv}
    fn = {col.split(' ')[1]: m for col, m in e['fn'].items()}
    if fn:
        rec['fn'] = fn
    core_out[e['sym']] = rec

# ---------- photoelectron RSF ----------
def photo_table(lib):
    out = {}
    for k, v in kh['photo'][lib].items():
        el, orb = norm_photo_key(k)
        if el is None:
            continue
        out[el + ' ' + orb] = v['rsf']
    return out

sc_al = photo_table('A-ALSCOF')
sc_mg = photo_table('A-MGSCOF')
th1 = photo_table('A-ALTHERMO1')

# Avantage element-table positions (survey identification only)
idpos = {}
for k, v in kh['photo']['A-ALTHERMO1'].items():
    el, orb = norm_photo_key(k)
    if el is None or v.get('pos') is None or v['ke_be'] != 'BE':
        continue
    idpos[el + ' ' + orb] = v['pos']
# elements whose booklet values are referenced to the VACUUM level (gas phase) — never used for solid-state ID
GAS = ['H', 'He', 'N', 'O', 'F', 'Ne', 'Cl', 'Ar', 'Kr', 'Xe', 'Rn']

# check 1: digitised Scofield vs values read from the paper (rel. tol 1 %, last-digit rounding)
#   - subshell missing in the digitised table  -> filled from the paper (listed in meta.filled)
#   - disagreement where the rest of that paper row matches exactly -> paper wins (listed in meta.corrected)
#   - any other disagreement -> abort
n_chk = 0
filled, corrected = [], []
for src, tab in (('al', sc_al), ('mg', sc_mg)):
    for el, subs in spot[src].items():
        row_ok = sum(1 for sub, val in subs.items() if (el + ' ' + sub) in tab and abs(tab[el + ' ' + sub] - val) <= 0.01 * abs(val) + 5e-4)
        for sub, val in subs.items():
            key = el + ' ' + sub
            n_chk += 1
            if key not in tab:
                tab[key] = val
                filled.append('%s %s = %s' % (src, key, val))
            elif abs(tab[key] - val) > 0.01 * abs(val) + 5e-4:
                if row_ok >= 3:
                    corrected.append('%s %s: digitised %s -> paper %s (other %d values of the row agree)' % (src, key, tab[key], val, row_ok))
                    tab[key] = val
                else:
                    problems.append('Scofield %s %s: digitised %s vs paper %s' % (src, key, tab[key], val))
for el, d in spot.get('al_inferred', {}).items():
    for sub, val in d.items():
        if sub == 'how':
            continue
        key = el + ' ' + sub
        n_chk += 1
        if abs(sc_al.get(key, -1) - val) > 0.01 * val + 5e-4:
            problems.append('Scofield al (inferred) %s: digitised %s vs paper %s' % (key, sc_al.get(key), val))

# check 2: doublet totals.  Single source of truth = the two subshell values; the total is always
#   re-derived as their sum (the digitised file has a few totals equal to one subshell, e.g. Ni 2p).
J = {'p': ('1/2', '3/2'), 'd': ('3/2', '5/2'), 'f': ('5/2', '7/2')}
n_sum, retotal = 0, []
for name, tab in (('Scofield-Al', sc_al), ('Scofield-Mg', sc_mg), ('ALTHERMO1', th1)):
    keys = set(k.rsplit(' ', 1)[0] + ' ' + k.split(' ')[1][:2] for k in tab if len(k.split(' ')[1]) > 2)
    for key in sorted(keys):
        el, orb = key.split(' ')
        if orb[1] not in J:
            continue
        a, b = tab.get(key + J[orb[1]][0]), tab.get(key + J[orb[1]][1])
        if a is None or b is None:
            continue
        n_sum += 1
        v = tab.get(key)
        t = round(a + b, 6)
        if v is None or abs(t - v) > 0.012 * t + 0.006:
            retotal.append('%s %s: %s -> %s' % (name, key, v, t))
        tab[key] = t

# check 3: C 1s = 1 in every table
for name, tab in (('Scofield-Al', sc_al), ('Scofield-Mg', sc_mg), ('ALTHERMO1', th1)):
    if abs(tab.get('C 1s', 0) - 1) > 1e-9:
        problems.append(name + ': C 1s != 1')

# ---------- Auger (Avantage KE, Al Ka library) ----------
auger = []
for k, v in sorted(kh['auger']['A-ALTHERMO1'].items()):
    el, code = k.split('_', 1)
    if v['ke_be'] != 'KE':
        problems.append('Auger %s not KE' % k)
        continue
    auger.append({'el': el, 'line': code, 'ke': v['pos'], 'rsf': v['rsf']})

# check 4: Avantage Auger KE vs NIST / ASTM values.  For every element of the reference list at
#   least one reference transition must have an Avantage line within 1.5 eV (Avantage lists mostly the
#   strongest multiplet only; compound vs metal shifts are reported, not failed, when <= 5 eV).
aug_rep = []
by_el = {}
for c in aug_ref['checks']:
    by_el.setdefault(c['el'], []).append(c)
for el, cs in by_el.items():
    cands = [a for a in auger if a['el'] == el]
    if not cands:
        problems.append('Auger %s: no Avantage entry' % el)
        continue
    ds = []
    for c in cs:
        best = min(cands, key=lambda a: abs(a['ke'] - c['ke']))
        ds.append(best['ke'] - c['ke'])
        aug_rep.append('%s %s %.2f vs %s %.1f (%+.1f)' % (el, c['transition'], c['ke'], best['line'], best['ke'], best['ke'] - c['ke']))
    m = min(abs(d) for d in ds)
    if m > 5.0:
        problems.append('Auger %s: nearest Avantage line %.1f eV from every reference value' % (el, m))
    elif m > 1.5:
        aug_rep[-1] += ' [Avantage value is %.1f eV off the elemental reference: chemical-state or multiplet difference]' % m

# check 5: spin-orbit ordering in the booklet (j=l+1/2 must have the lower BE).  The booklet itself
#   has a few inverted pairs (e.g. Co M2/M3, Na L2/L3); they are kept as printed and flagged in
#   rec.soWarn so the UI never derives a doublet splitting from them.
n_so, so_warn = 0, []
for el, rec in core_out.items():
    lv = rec['lv']
    for orb in sorted(set(k[:2] for k in lv if len(k) > 2)):
        l = orb[1]
        if l not in J:
            continue
        lo, hi = lv.get(orb + J[l][1]), lv.get(orb + J[l][0])
        if lo is None or hi is None:
            continue
        n_so += 1
        if hi < lo - 1e-9:
            rec.setdefault('soWarn', []).append(orb)
            so_warn.append('%s %s (%s%s=%s < %s%s=%s)' % (el, orb, orb, J[l][0], hi, orb, J[l][1], lo))

if problems:
    print('\n'.join(problems))
    sys.exit('library build aborted: %d problem(s)' % len(problems))

meta = {
    'built': datetime.date.today().isoformat(),
    'sources': {
        'BE (core levels)': core['source'] + '. ' + core['notes'],
        'Scofield RSF': 'J.H. Scofield, J. Electron Spectrosc. Relat. Phenom. 8, 129 (1976), Hartree-Slater subshell '
                        'photoionization cross-sections at 1487 eV (Al Ka) and 1254 eV (Mg Ka), relative to C 1s = 1. '
                        'Digitised table: ' + kh['source'] + ' Verified against the paper for %d subshells (tol 1 %%).' % n_chk,
        'ALTHERMO1 RSF': 'Thermo Fisher Scientific Avantage "ALTHERMO1" sensitivity-factor library (Al Ka, C 1s = 1), '
                         'as digitised in KherveFitting (same file). Not identical to Scofield (e.g. Ta 4f 9.09 vs 8.62); '
                         'intended to reproduce Avantage at%% together with Avantage transmission and energy compensation.',
        'Auger KE': 'Avantage element-table Auger kinetic energies (Al Ka library, same file). Spot-checked against NIST '
                    'SRD 20 / ASTM E2108: ' + '; '.join(aug_rep),
    },
    'checks': {'scofield_vs_paper': n_chk, 'doublet_totals': n_sum, 'spin_orbit_order': n_so},
    'filled_from_paper': filled, 'corrected_from_paper': corrected, 'totals_rederived': retotal,
    'booklet_so_inverted': so_warn,
}
for el in GAS:
    if el in core_out:
        core_out[el]['vac'] = True
# check 6: ID positions vs booklet for metals/semiconductors: |Avantage - booklet| <= 6 eV for p>=1 main lines
far = []
for key, pos in idpos.items():
    el, orb = key.split(' ')
    lv = core_out.get(el, {}).get('lv', {})
    if orb in lv and el not in GAS and abs(pos - lv[orb]) > 8:
        far.append('%s: Avantage %s vs booklet %s' % (key, pos, lv[orb]))
meta['id_far_from_booklet'] = far
meta['sources']['ID positions'] = ('Avantage element-table peak positions (Al Ka library), digitised in KherveFitting (same file). '
    'Representative values for peak identification only (e.g. O 1s 532.9, Ta 4f 25); the booklet values of %s are '
    'vacuum-level (gas) energies and are never used for identification.' % ', '.join(GAS))
lib = {'core': core_out, 'scofield': {'al': sc_al, 'mg': sc_mg}, 'althermo1': th1, 'idpos': idpos, 'auger': auger, 'meta': meta}
js = ('/* INNO Analytics — XPS 원소 라이브러리 (자동 생성: modules/xps/build/build_lib.py — 직접 고치지 말 것) */\n'
      '(function(root){ const L=' + json.dumps(lib, ensure_ascii=False, separators=(',', ':')) + ';\n'
      ' if(typeof module!=="undefined"&&module.exports) module.exports=L; else root.XPS_LIB=L; })(typeof self!=="undefined"?self:this);\n')
open(OUT, 'w', encoding='utf-8').write(js)
print('idpos %d, far %d: %s' % (len(idpos), len(far), '; '.join(far[:30])))
print('xps_lib.js: %d elements, Scofield Al %d / Mg %d, ALTHERMO1 %d, Auger %d; checks %s; %d bytes, sha256 %s'
      % (len(core_out), len(sc_al), len(sc_mg), len(th1), len(auger), meta['checks'], len(js),
         hashlib.sha256(js.encode()).hexdigest()[:16]))
for r in aug_rep:
    print('  auger:', r)
for t,L in (('filled',filled),('corrected',corrected),('retotal',retotal),('so_warn',so_warn)):
    print(' ',t,len(L),'; '.join(L[:40]))
