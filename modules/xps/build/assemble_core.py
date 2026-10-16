#!/usr/bin/env python3
"""Consolidate the raw X-ray Data Booklet Table 1-1 transcriptions (data/raw_fetch/*.json)
into data/core_xdb_table_1_1.json.

Pass-1 files give every element (values + footnote markers); pass-2 files re-transcribe the
rows flagged by the xraydb comparison as printed cell strings.  Every pass-2 cell is compared
with pass 1 (value, blank/non-blank, markers); any disagreement aborts.
Run:  python3 build/assemble_core.py
"""
import json, glob, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RAW = os.path.join(DATA, 'raw_fetch')

COLS = ["K 1s", "L1 2s", "L2 2p1/2", "L3 2p3/2", "M1 3s", "M2 3p1/2", "M3 3p3/2", "M4 3d3/2", "M5 3d5/2",
        "N1 4s", "N2 4p1/2", "N3 4p3/2", "N4 4d3/2", "N5 4d5/2", "N6 4f5/2", "N7 4f7/2",
        "O1 5s", "O2 5p1/2", "O3 5p3/2", "O4 5d3/2", "O5 5d5/2", "P1 6s", "P2 6p1/2", "P3 6p3/2"]
MARK = '*†ab'
# Known glyph typos on the web page (value is unambiguous); documented in README.
GLYPH_FIX = {(48, "N5 4d5/2"): ("l0.7", 10.7)}


def norm_marks(m):
    if m is None:
        return ''
    if isinstance(m, list):
        m = ''.join(m)
    return ''.join(c for c in MARK if c in m)  # canonical order *, †, a, b


def num(s):
    return float(s) if '.' in str(s) else int(s)


def main():
    p1, p2 = {}, {}
    for f in sorted(glob.glob(os.path.join(RAW, '*.json'))):
        doc = json.load(open(f, encoding='utf-8'))
        if doc['kind'] == 'pass1':
            for z, row in doc['rows'].items():
                p1[int(z)] = row
        else:
            for sym, cells in doc['rows'].items():
                p2[sym] = cells
    missing = sorted(set(range(1, 93)) - set(p1))
    if missing:
        sys.exit('missing Z: %s' % missing)

    elements, problems, glyph_notes = {}, [], []
    for z in range(1, 93):
        row = p1[z]
        sym = row['sym']
        # markers: list form ["K 1s*", ...] or dict form {col: "*"}
        mk = {}
        m = row.get('markers', {})
        if isinstance(m, list):
            for s in m:
                for c in COLS:
                    if s.startswith(c):
                        mk[c] = norm_marks(s[len(c):])
        else:
            for c, v in m.items():
                mk[c] = norm_marks(v)
        lv, fn = {}, {}
        for c in COLS:
            if c not in row:
                continue
            raw = row[c]
            if (z, c) in GLYPH_FIX:
                assert str(raw) == GLYPH_FIX[(z, c)][0], (z, c, raw)
                val = GLYPH_FIX[(z, c)][1]
                glyph_notes.append('%s %s printed "%s" on the web page; read as %s' % (sym, c, raw, val))
            else:
                val = num(raw)
            lv[c] = val
            if mk.get(c):
                fn[c] = mk[c]
        # pass-2 comparison
        if sym in p2:
            cells = p2[sym]
            for c, cell in zip(COLS, cells):
                if cell in ('BLANK', 'DASH'):
                    if c in lv:
                        problems.append('%s %s: pass1=%s pass2=%s' % (sym, c, lv[c], cell))
                    continue
                v2 = re.sub('[' + MARK + ']', '', cell)
                m2 = norm_marks(cell[len(v2):] if cell.startswith(v2) else re.sub('[0-9.]', '', cell))
                if c not in lv or abs(float(v2) - float(lv[c])) > 1e-9:
                    problems.append('%s %s: pass1=%s pass2=%s' % (sym, c, lv.get(c), cell))
                if m2 != fn.get(c, ''):
                    problems.append('%s %s marker: pass1="%s" pass2="%s"' % (sym, c, fn.get(c, ''), m2))
        elements[str(z)] = {'sym': sym, 'lv': lv, 'fn': fn, 'pass2_verified': sym in p2}

    if problems:
        print('\n'.join(problems))
        sys.exit('pass1/pass2 disagreement: resolve before building')

    out = {
        'source': 'G.P. Williams, "Electron Binding Energies", in X-ray Data Booklet (Lawrence Berkeley National '
                  'Laboratory, Center for X-ray Optics), Section 1.1, Table 1-1 "Electron binding energies, in '
                  'electron volts, for the elements in their natural forms". https://xdb.lbl.gov/Section1/Table_1-1a.htm , '
                  'https://xdb.lbl.gov/Section1/Table_1-1b.htm , https://xdb.lbl.gov/Section1/Table_1-1c.htm ; '
                  'introduction https://xdb.lbl.gov/Section1/Sec_1-1.html',
        'notes': 'Transcribed 2026-10-15 via WebFetch (raw transcriptions in data/raw_fetch/). Cross-checked cell by cell '
                 'against xraydb 4.5.8 (Elam/Ravel/Sieber edge table) and by a second targeted transcription of every row '
                 'that disagreed with xraydb. Columns are X-ray level (spectroscopic) names. Blank and dash cells of the '
                 'booklet are omitted. Reference levels (booklet text): vacuum level for the rare gases and H2, N2, O2, F2, '
                 'Cl2; Fermi level for metals; top of the valence band for semiconductors.',
        'footnotes': {
            '*': 'Values taken from Ref. 2 (M. Cardona and L. Ley, Eds., Photoemission in Solids I: General Principles, Springer-Verlag, Berlin, 1978), with additional corrections.',
            '†': 'Values taken from Ref. 3 (J.C. Fuggle and N. Martensson, "Core-level binding energies in metals", J. Electron Spectrosc. Relat. Phenom. 21, 275 (1980)).',
            'a': 'One-particle approximation not valid owing to short core-hole lifetime.',
            'b': 'Value derived from Ref. 1.',
            'unmarked': 'Values taken from Ref. 1 (J.A. Bearden and A.F. Burr, "Reevaluation of X-ray atomic energy levels", Rev. Mod. Phys. 39, 125 (1967)).'
        },
        'glyph_fixes': glyph_notes,
        'columns': COLS,
        'elements': elements,
    }
    json.dump(out, open(os.path.join(DATA, 'core_xdb_table_1_1.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    n = sum(len(e['lv']) for e in elements.values())
    print('core: %d elements, %d values, pass-2 verified rows: %d' % (len(elements), n, len(p2)))


if __name__ == '__main__':
    main()
