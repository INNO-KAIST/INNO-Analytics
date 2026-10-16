#!/usr/bin/env python3
"""Extract numeric facts (sensitivity factors, Auger kinetic energies) from the KherveFitting
library file into data/kherve_extract.json.  Only numbers and line labels are taken; no code.

Usage:
  git clone --depth 1 https://github.com/gkerherve/KherveFitting.git /tmp/kf
  python3 -I build/extract_kherve.py /tmp/kf/KherveFitting_library.json
"""
import json, os, sys, hashlib, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'data', 'kherve_extract.json')
LIBS = ['A-ALSCOF', 'A-MGSCOF', 'A-ALTHERMO1']
REF_AUGER_LIBS = ['C-Any', 'ScientaC1s-Ga']   # generic element Auger positions, used only for cross-checks


def main(path):
    raw = open(path, 'rb').read()
    d = json.loads(raw)
    commit = None
    try:
        commit = subprocess.check_output(['git', '-C', os.path.dirname(path), 'log', '-1', '--format=%H'],
                                         text=True).strip()
    except Exception:
        pass
    photo = {lib: {} for lib in LIBS}
    auger = {lib: {} for lib in LIBS}
    ref = {lib: {} for lib in REF_AUGER_LIBS}
    for key, entry in sorted(d.items()):
        for lib in LIBS:
            if lib in entry:
                e = entry[lib]
                if e['auger'] == 0:
                    photo[lib][key] = {'rsf': e['rsf'], 'ke_be': e['ke_be'], 'pos': e.get('position')}
                else:
                    auger[lib][key] = {'pos': e['position'], 'ke_be': e['ke_be'], 'rsf': e['rsf']}
        for lib in REF_AUGER_LIBS:
            if lib in entry and entry[lib]['auger'] == 1:
                ref[lib][key] = {'pos': entry[lib]['position'], 'ke_be': entry[lib]['ke_be']}
    out = {
        'source': 'Numeric values digitised in KherveFitting (G. Kerherve), file KherveFitting_library.json, '
                  'https://github.com/gkerherve/KherveFitting (commit %s, sha256 %s). Library A-ALSCOF / A-MGSCOF = '
                  'J.H. Scofield, "Hartree-Slater subshell photoionization cross-sections at 1254 and 1487 eV", '
                  'J. Electron Spectrosc. Relat. Phenom. 8, 129-137 (1976), doi:10.1016/0368-2048(76)80015-1, relative '
                  'to C 1s; A-ALTHERMO1 = Thermo Fisher Scientific Avantage "ALTHERMO1" sensitivity-factor library '
                  '(Al Ka). Auger kinetic energies in the A-* libraries are those of the Avantage element tables.'
                  % (commit, hashlib.sha256(raw).hexdigest()),
        'notes': 'Only numeric facts are extracted (GPL code is not used). Keys are "Sym_orbital" as in the source; '
                 'photoelectron "pos" = Avantage element-table position, used ONLY for survey peak identification '
                 '(coarse representative values, e.g. O 1s 532.9, Ta 4f 25); core-level energies for display and doublet '
                 'splittings come from the X-ray Data Booklet. C-Any / ScientaC1s-Ga Auger positions have no stated provenance and are used only as a '
                 'secondary consistency reference, never as output values.',
        'photo': photo, 'auger': auger, 'auger_reference_only': ref,
    }
    json.dump(out, open(OUT, 'w'), indent=1, sort_keys=False)
    print('wrote', OUT, {l: len(photo[l]) for l in LIBS}, {l: len(auger[l]) for l in LIBS})


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/kf/KherveFitting_library.json')
