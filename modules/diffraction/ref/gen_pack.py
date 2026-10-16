"""Diffraction 내장 결정상 팩: ref/pack/*.cif (COD 원본) → src/diff_pack.js · test/ref_pack.json · ref/pack/SOURCES.md
· 기록(PhaseRec v1)은 gemmi.cif 로 CIF 를 직접 읽어 만든다 (diff_core.js 의 CIF 해석과 독립).
· 기준값: gemmi(격자 d, 대칭 연산 소멸칙) · pymatgen CifParser(점유율 포함 자리 펼치기) × gemmi c4322 전자 산란 인자.
사용: python3 -I ref/gen_pack.py"""
import os, re, json, math, hashlib, warnings
from fractions import Fraction
import numpy as np, gemmi
from pymatgen.io.cif import CifParser
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..'); PACK = os.path.join(HERE, 'pack')
VERSION = '2026-10-17'; GMAX = 10.0  # 1/nm
WY = 'Wyckoff 1963 (AMCSD) 표준 참조 구조'
# short, COD id, 이름, 표시 화학식, CIF SHA-256 앞 16자리(COD 에서 받은 원문), 고른 이유, 기록 note, B 덮어쓰기
P = [
 ('Au', '9008463', 'Gold', 'Au', 'c7c350025f5f5887', WY + ', a = 4.07825 Å (현대 값 4.0782 Å 와 일치)', 'Calibration standard (fcc).', None),
 ('Al', '9008460', 'Aluminium', 'Al', '6a5e880d63896bbc', WY + ', a = 4.04958 Å', 'Calibration standard (fcc).', None),
 ('Si', '9011998', 'Silicon', 'Si', 'a3dfc7a2fcd58df8', '298.15 K 다중반사 정밀 격자상수 측정(a = 5.430941 Å), DOI 있음', 'Calibration standard (diamond, origin choice 1).', None),
 ('Cu', '9008468', 'Copper', 'Cu', '0804d4ddcc3748a4', WY + ', a = 3.61496 Å', 'Calibration standard (fcc).', None),
 ('Ni', '9008476', 'Nickel', 'Ni', '45d7bfec776bf7a4', WY + ', a = 3.52387 Å', 'Calibration standard (fcc).', None),
 ('Pt', '9008480', 'Platinum', 'Pt', '0373322030e51bc9', WY + ', a = 3.9231 Å', 'Calibration standard (fcc).', None),
 ('Ag', '9008459', 'Silver', 'Ag', 'e19ea35d6106a90d', WY + ', a = 4.0862 Å', 'Calibration standard (fcc).', None),
 ('aNaYF4', '1517676', 'α-NaYF4 (cubic)', 'NaYF4', '6649962009b99ef6', 'COD 의 화학량론 α-NaYF4 유일 항목(Fm-3m, Na/Y 4a 자리 0.5/0.5 무질서, F 8c). 1517677 은 Na0.39Y0.61F2.2 비화학량론 결함 구조라 제외', 'Disordered Na/Y on 4a (0.5/0.5), Z = 2. Old (1964) lattice constant, no displacement parameters.', None),
 ('bNaYF4', '1517672', 'β-NaYF4 (hexagonal)', 'Na1.5Y1.5F6', 'ec7b33f744b23fcc', '고압 연구(Grzechnik 2002)의 상압 P-6 구조 — 1a Y, 1f Y/Na 0.5/0.5, 2h Na 0.5 (CIF 점유율 그대로). 같은 논문의 1517673–75 는 고압', 'P-6, Z = 1 (Na1.5Y1.5F6 = 1.5 NaYF4). Ambient-pressure block of a high-pressure study; partial occupancies as published.', None),
 ('CaF2', '1000043', 'Fluorite CaF2', 'CaF2', 'ed4fc6d25246923c', '중성자 회절 정밀화(R = 0.008), 화학량론 CaF2. 단 a = 5.462(30) Å 로 불확도가 크다(흔한 값 5.463 Å 와 0.02 % 차)', 'Fluorite structure.', None),
 ('LiYF4', '2001129', 'LiYF4 (scheelite)', 'LiYF4', '9eca30c100ec66ca', '상온 단결정 정밀화(Acta C, 레이저 호스트), I41/a 원점 2', 'I41/a origin choice 2. CIF Uiso values (0.10-0.20 A^2) are implausible for RT and were not used (B = 0).', 0.0),
 ('Y2O3', '1009013', 'Y2O3 (bixbyite, C-type)', 'Y2O3', 'f8173f5a23fa50f0', 'COD 의 순수 Y2O3 Ia-3 중 위치 매개변수 정밀도가 가장 높음(R = 0.011). 1513300(a = 10.604 Å, 295 K)은 Uovl = 0.00058 Å² 로 비정상적이고 좌표 불확도가 큼', 'Ia-3. Y2O3 phase of a Zr-doped two-phase Rietveld study (formula given as pure Y2O3); a = 10.6056 A is 0.015 % above the usual 10.604 A. No displacement parameters.', None),
 ('aFe', '9008536', 'α-Fe (ferrite, bcc)', 'Fe', '58d627740a3a7eae', WY + ', 298 K, a = 2.8665 Å', 'bcc iron.', None),
 ('Fe3O4', '9007706', 'Magnetite Fe3O4', 'Fe3O4', 'c9cd18c4fe7c7cf7', '단결정 상온 정밀화(Fleet 1984, 시료 633 Fe3.005O4), a = 8.3969 Å 가 표준 값과 일치. 1981 년 항목(9007644)은 a = 8.3941 Å', 'Fd-3m origin choice 2. Only anisotropic U given in the CIF; B = 0 (isotropic B taken only from Uiso/Biso).', None),
 ('TiO2-rutile', '9015662', 'Rutile TiO2', 'TiO2', 'cc2a7010508f50f6', '중성자 분말 상온 정밀화(Howard 1991), 루틸·아나타제 같은 논문', 'P42/mnm. Only anisotropic U given; B = 0.', None),
 ('TiO2-anatase', '9015929', 'Anatase TiO2', 'TiO2', '9df468431a17fa79', '중성자 분말 상온 정밀화(Howard 1991)', 'I41/amd origin choice 1. Only anisotropic U given; B = 0.', None),
 ('ZnO', '2300112', 'Zincite ZnO (wurtzite)', 'ZnO', 'bf25e35b06c55d9f', '293 K 단결정 정밀화(R = 0.0149, Uiso 있음) — 고압 연구의 상압(블록 zno0) 측정, a·c 가 흔한 값(3.2498, 5.2066 Å)과 0.02 % 이내', 'P63mc, ambient-pressure block of a high-pressure study.', None),
 ('Al2O3', '9007496', 'Corundum α-Al2O3', 'Al2O3', 'fbc95a4d3d0800ce', '정밀 구조인자 연구(Lewis 1982), a = 4.7602, c = 12.9933 Å', 'R-3c, hexagonal axes. Only anisotropic U given; B = 0.', None),
 ('NaCl', '9008678', 'Halite NaCl', 'NaCl', 'f6719e52179a4908', WY + ', a = 5.64056 Å', 'Rock-salt structure.', None),
 ('MgO', '9008671', 'Periclase MgO', 'MgO', 'bf0b4bae1964dab2', WY + ', a = 4.2112 Å', 'Rock-salt structure.', None),
 ('graphite-2H', '9011577', 'Graphite-2H', 'C', '35d1661c3131acc9', '천연 결정 중성자 회절(Trucano & Chen 1975), 2H 적층 P63/mmc', 'P63/mmc (2H stacking). Only anisotropic U given; B = 0.', None),
]
NOT_FOUND = [('α-NaLuF4', 'Na·Lu·F 원소 검색(3–6 원소) 결과 NaLuF4 구조 없음'), ('β-NaLuF4', '같음'), ('NaGdF4 (육방)', 'Na·Gd·F 원소 검색(3–6 원소) 결과 NaGdF4 구조 없음')]

def num(s):
    s = gemmi.cif.as_string(s) if s is not None else ''
    s = re.sub(r'\(.*?\)', '', s).strip()
    return None if s in ('', '.', '?') else float(s)

def element(sym):
    m = re.match(r'([A-Z][a-z]?)', sym or '')
    if not m: return None
    for cand in (m.group(1), m.group(1)[0]):
        try:
            if gemmi.Element(cand).atomic_number > 0: return cand
        except Exception: pass
    return None

def fmt_op(op):
    # diff_core.js 의 DC.formatSymop 과 같은 꼴: 축 순서 x,y,z, 평행이동은 끝에 n/d (0 ≤ t < 1)
    out = []
    for i in range(3):
        s = ''
        for j, ax in enumerate('xyz'):
            r = op.rot[i][j] // gemmi.Op.DEN
            if r: s += ('+' if s else '') + ax if r > 0 else '-' + ax
        t = Fraction(op.tran[i] % gemmi.Op.DEN, gemmi.Op.DEN)
        if t: s += ('+' if s else '') + f'{t.numerator}/{t.denominator}'
        out.append(s or '0')
    return ','.join(out)

def citation(b):
    au = [gemmi.cif.as_string(v) for v in b.find_values('_publ_author_name')]
    au = (', '.join(au[:3]) + (' et al.' if len(au) > 3 else '')) if au else ''
    g = lambda k: gemmi.cif.as_string(b.find_value(k)) if b.find_value(k) else ''
    pp = g('_journal_page_first') + ('-' + g('_journal_page_last') if g('_journal_page_last') else '')
    return f"{au} ({g('_journal_year')}) {g('_journal_name_full')} {g('_journal_volume')}, {pp}".strip()

recs, ref, rows_md = [], {}, []
for short, cod, name, formula, sha, why, note, Bover in P:
    path = os.path.join(PACK, f'{short}_{cod}.cif')
    raw = open(path, 'rb').read(); h = hashlib.sha256(raw).hexdigest()
    assert h.startswith(sha), f'{short}: CIF SHA-256 {h[:16]} != {sha} (COD 원문과 다름)'
    b = gemmi.cif.read(path).sole_block()
    assert b.find_value('_cod_database_code') and gemmi.cif.as_string(b.find_value('_cod_database_code')) == cod
    a, bb, c = (num(b.find_value(k)) for k in ('_cell_length_a', '_cell_length_b', '_cell_length_c'))
    al, be, ga = (num(b.find_value(k)) or 90.0 for k in ('_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'))
    trip = []
    for tag in ('_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz'):
        col = b.find_values(tag)
        if len(col): trip = [gemmi.cif.as_string(v) for v in col]; break
    assert trip, f'{short}: 대칭 연산 없음'
    ops = [gemmi.Op(t.replace(' ', '')) for t in trip]
    gops = gemmi.GroupOps(ops)
    sgn = int(num(b.find_value('_space_group_IT_number') or b.find_value('_symmetry_Int_Tables_number')))
    hm = gemmi.cif.as_string(b.find_value('_space_group_name_H-M_alt') or b.find_value('_symmetry_space_group_name_H-M'))
    sgf = gemmi.find_spacegroup_by_ops(gops)
    assert sgf and sgf.number == sgn, f'{short}: 연산의 공간군 {sgf and sgf.number} != CIF #{sgn}'
    t = b.find('_atom_site_', ['label', '?type_symbol', 'fract_x', 'fract_y', 'fract_z', '?occupancy', '?U_iso_or_equiv', '?B_iso_or_equiv'])
    assert len(t), f'{short}: 원자 자리 없음'
    atoms = []
    for r in t:
        lab = gemmi.cif.as_string(r[0]); el = element(gemmi.cif.as_string(r[1]) if r.has(1) else lab) or element(lab)
        occ = num(r[5]) if r.has(5) else None; U = num(r[6]) if r.has(6) else None; B = num(r[7]) if r.has(7) else None
        Bv = B if B is not None else (8 * math.pi ** 2 * U if U is not None else 0.0)
        if Bover is not None: Bv = Bover
        atoms.append(dict(label=lab, el=el, x=num(r[2]), y=num(r[3]), z=num(r[4]), occ=1.0 if occ is None else occ, B=round(Bv, 6)))
    doi = gemmi.cif.as_string(b.find_value('_journal_paper_doi')) if b.find_value('_journal_paper_doi') else ''
    rec = dict(id=f'pack:{short}', name=name, formula=formula,
               cell=dict(a=round(a / 10, 10), b=round(bb / 10, 10), c=round(c / 10, 10), al=al, be=be, ga=ga),
               sg=dict(number=sgn, symbol=hm), symops=[fmt_op(o) for o in ops], atoms=atoms,
               source=dict(kind='COD', codId=cod, doi=doi, citation=citation(b)), verified=True, note=note)
    recs.append(rec)
    # ---- 기준값 ----
    uc = gemmi.UnitCell(a, bb, c, al, be, ga)
    st = CifParser(path, site_tolerance=1e-3, frac_tolerance=0).parse_structures(primitive=False)[0]
    Bel = {}
    for at in atoms: Bel.setdefault(at['el'], set()).add(at['B'])
    assert all(len(v) == 1 for v in Bel.values()), f'{short}: 같은 원소에 다른 B — 원소별 대응 불가'
    Bel = {k: next(iter(v)) for k, v in Bel.items()}
    sites = [([(sp.symbol, occ) for sp, occ in s.species.items()], np.array(s.frac_coords)) for s in st]
    ff = {el: gemmi.Element(el).c4322 for el in Bel}
    H = [math.ceil(GMAX * x / 10) for x in (a, bb, c)]
    rows = []
    for hh in range(-H[0], H[0] + 1):
        for k in range(-H[1], H[1] + 1):
            for l in range(-H[2], H[2] + 1):
                if hh == k == l == 0: continue
                dA = uc.calculate_d([hh, k, l]); g = 10 / dA
                if g > GMAX: continue
                s2 = 1 / (4 * dA * dA); hv = np.array([hh, k, l])
                F = 0j
                for sp, xyz in sites:
                    f = sum(o * ff[e].calculate_sf(s2) * math.exp(-Bel[e] * s2) for e, o in sp)
                    F += f * np.exp(2j * np.pi * hv.dot(xyz))
                rows.append([hh, k, l, float(f'{dA / 10:.12g}'), bool(gops.is_systematically_absent([hh, k, l])), float(f'{abs(F) ** 2:.9g}')])
    Z = num(b.find_value('_cell_formula_units_Z'))
    comp = {e.symbol: float(v) for e, v in st.composition.element_composition.items()}  # 산화수 없는 원소 기호로
    fsum = gemmi.cif.as_string(b.find_value('_chemical_formula_sum'))
    rho = num(b.find_value('_exptl_crystal_density_diffrn')) if b.find_value('_exptl_crystal_density_diffrn') else None
    ref[short] = dict(cif=f'{short}_{cod}.cif', sg=sgn, hm=hm, Z=Z, formula_sum=fsum, comp_cell=comp, density_cif=rho, density_pmg=float(st.density),
                      cellA=[a, bb, c, al, be, ga], nops=len(ops), rows=rows)
    print(f'{short:13s} #{sgn:3d} {hm:14s} ops {len(ops):3d}  atoms {len(atoms)}  pmg sites {len(st):3d}  hkl {len(rows):5d} allowed {sum(1 for r in rows if not r[4]):5d}  Z {Z:g}  rho {rho} / {st.density:.3f}')
    a_s = (f'a = {a / 10:.6g}' if a == bb == c else f'a = {a / 10:.6g}, c = {c / 10:.6g}' if a == bb else f'{a / 10:.6g} / {bb / 10:.6g} / {c / 10:.6g}')
    rows_md.append(f'| pack:{short} | {name} | {formula} (Z = {Z:g}) | {hm} (#{sgn}) | {a_s} | [{cod}](https://www.crystallography.net/cod/{cod}.html) | {doi or "— (CIF 에 없음)"} | {citation(b)} | {why} |')

# ---- src/diff_pack.js ----
js = ['/* INNO Analytics — Diffraction 내장 결정상 팩 (PhaseRec v1). ref/gen_pack.py 로 생성 — 손으로 고치지 말 것.',
      '   출처: COD (ref/pack/*.cif 원문, 출처표 ref/pack/SOURCES.md). 길이 nm·각도 도, B Å² (Uiso×8π² 또는 Biso, 없으면 0). */',
      '(function(){', '"use strict";', f'const P={{version:"{VERSION}",phases:[']
js += [json.dumps(r, ensure_ascii=False, separators=(',', ':')) + ',' for r in recs]
js[-1] = js[-1][:-1]
js += [']};', 'if(typeof module!=="undefined"&&module.exports) module.exports=P;', '(typeof self!=="undefined"?self:globalThis).DIFF_PACK=P;', '})();', '']
open(os.path.join(ROOT, 'src', 'diff_pack.js'), 'w').write('\n'.join(js))
json.dump(dict(version=VERSION, gmax=GMAX, phases=ref), open(os.path.join(ROOT, 'test', 'ref_pack.json'), 'w'), separators=(',', ':'))

# ---- SOURCES.md ----
md = [f'# 내장 결정상 팩 출처표 (DIFF_PACK {VERSION})', '',
      '`ref/gen_pack.py` 가 이 폴더의 CIF 에서 `src/diff_pack.js` 와 함께 생성한다. CIF 는 COD 에서 받은 원문 그대로이며, 생성기가 SHA-256 앞 16 자리를 확인한다.',
      '격자상수는 nm (CIF Å ÷ 10). B(Å²) 는 `_atom_site_U_iso_or_equiv`×8π² 또는 `_atom_site_B_iso_or_equiv`, 둘 다 없으면 0 (비등방 U 만 있는 CIF 는 B = 0).', '',
      '| id | 이름 | 화학식 | 공간군 (CIF 설정) | 격자상수 (nm) | COD | DOI | 인용 | 고른 이유 |', '|---|---|---|---|---|---|---|---|---|'] + rows_md + [
      '', '## COD 에 없는 상', '', '| 상 | 결과 |', '|---|---|'] + [f'| {n} | 없음 — {w} (2026-10-17 COD 검색) |' for n, w in NOT_FOUND] + [
      '', '## 확인할 점', '',
      '- LiYF4 (2001129): CIF 의 Uiso(Y 0.100, Li 0.20, F 0.170 Å²)는 상온의 흔한 값(0.005–0.015 Å²)보다 열 배가량 커서 입력 오류로 보고 쓰지 않았다 (B = 0). 좌표·격자는 그대로.',
      '- Fe3O4·루틸·아나타제·커런덤·흑연: CIF 에 비등방 U 만 있어 B = 0. 루틸·아나타제(9015662, 9015929)의 비등방 U 열은 이름(U11 U22 U33 U12 U13 U23)과 값 배열(U11 U12 U13 U22 U23 U33 로 읽어야 정방정 자리 대칭에 맞음)이 어긋나 있다 — 쓰지 않았으므로 팩에는 영향 없음.',
      '- α-NaYF4 (1517676, Roy & Roy 1964): 오래된 격자상수(a = 5.47 Å, 유효숫자 3)이고 변위 매개변수 없음. 실험 시료의 α 상 격자상수는 조성·크기에 따라 다를 수 있다.',
      '- β-NaYF4 (1517672): CIF 의 0.6667/0.3333 좌표를 2/3·1/3 로 바꾸지 않고 그대로 썼다. Na2 (2h) z = 0.689(5) 점유 0.5.',
      '- Y2O3 (1009013): Zr 도핑 2 상 시료의 Y2O3 상 (화학식은 순수 Y2O3 로 기록). 측정 온도 미기재(상온으로 봄).',
      '- Wyckoff (1963) 원소·암염 항목은 AMCSD 경유이며 DOI·측정 온도가 없는 것이 많다 (α-Fe 는 298 K 명시).', '',
      '## 사용 조건', '',
      '- COD 는 수록 자료를 공개 영역(public domain)으로 둔다. CIF 머리말 문구는 두 가지다: (1) "All data on this site have been placed in the public domain by the contributors." — 1000043, 1009013, 1517672, 1517676, 9008460; (2) AMCSD·IUCr 제공 항목(90xxxxx, 2001129, 2300112) — "The file may be used within the scientific community so long as proper attribution is given to the journal article from which the data were obtained." 그래서 각 상의 `source.citation`·`doi` 로 원 논문을 밝힌다.', '']
open(os.path.join(PACK, 'SOURCES.md'), 'w').write('\n'.join(md))
print(f'{len(recs)} phases → src/diff_pack.js, test/ref_pack.json, ref/pack/SOURCES.md')
