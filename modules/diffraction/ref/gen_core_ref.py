"""Diffraction core 기준 자료: 시험 구조(CIF) + gemmi(소멸칙·전자 구조인자 c4322)·pymatgen(TEMCalculator d·정대축 반사) 결과.
사용: python3 -I ref/gen_core_ref.py   → test/cif/*.cif, test/ref_core.json"""
import os, json, numpy as np, gemmi, warnings
warnings.filterwarnings('ignore')
from pymatgen.core import Structure, Lattice
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.diffraction.tem import TEMCalculator
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', 'test')
S = {
 'Au_Fm-3m': Structure.from_spacegroup('Fm-3m', Lattice.cubic(4.0782), ['Au'], [[0, 0, 0]]),
 'Si_Fd-3m': Structure.from_spacegroup('Fd-3m', Lattice.cubic(5.431), ['Si'], [[0, 0, 0]]),
 'NaCl_Fm-3m': Structure.from_spacegroup('Fm-3m', Lattice.cubic(5.64), ['Na', 'Cl'], [[0, 0, 0], [0.5, 0.5, 0.5]]),
 'TiO2_rutile': Structure.from_spacegroup('P4_2/mnm', Lattice.tetragonal(4.594, 2.959), ['Ti', 'O'], [[0, 0, 0], [0.3048, 0.3048, 0]]),
 'bNaYF4_P-6': Structure.from_spacegroup('P-6', Lattice.hexagonal(5.96, 3.53), ['Na', 'Y', 'F', 'F'], [[2/3, 1/3, 0.5], [0, 0, 0], [0.64, 0.73, 0], [0.24, 0.76, 0.5]]),
 'LiYF4_I41a': Structure.from_spacegroup('I4_1/a', Lattice.tetragonal(5.16, 10.74), ['Li', 'Y', 'F'], [[0, 0.25, 0.625], [0, 0.25, 0.125], [0.2186, 0.0857, 0.4565]]),
 'Fe3O4_Fd-3m': Structure.from_spacegroup('Fd-3m', Lattice.cubic(8.396), ['Fe', 'Fe', 'O'], [[0.125, 0.125, 0.125], [0.5, 0.5, 0.5], [0.2549, 0.2549, 0.2549]]),
 'ZrO2_P21c': Structure.from_spacegroup('P2_1/c', Lattice.monoclinic(5.1505, 5.2116, 5.3173, 99.23), ['Zr', 'O', 'O'], [[0.2754, 0.0395, 0.2083], [0.0700, 0.3317, 0.3447], [0.4496, 0.7569, 0.4792]]),
 'Y2O3_Ia-3': Structure.from_spacegroup('Ia-3', Lattice.cubic(10.604), ['Y', 'Y', 'O'], [[0.25, 0.25, 0.25], [0.9677, 0, 0.25], [0.3907, 0.1520, 0.3804]]),
 'Al2O3_R-3c': Structure.from_spacegroup('R-3c', Lattice.hexagonal(4.759, 12.991), ['Al', 'O'], [[0, 0, 0.3523], [0.3064, 0, 0.25]]),
}
ref = {}
for name, st in S.items():
    p = os.path.join(OUT, 'cif', name + '.cif'); CifWriter(st, symprec=1e-3).write_file(p)
    ss = gemmi.read_small_structure(p); ss.setup_cell_images()
    doc = gemmi.cif.read(p); blk = doc.sole_block()
    trip = list(blk.find_loop('_symmetry_equiv_pos_as_xyz'))
    gops = gemmi.GroupOps([gemmi.Op(t.strip("'\"")) for t in trip])
    sg = gemmi.find_spacegroup_by_ops(gops)
    cell = ss.cell
    # 구조인자 정답: 단위세포 원자 = pymatgen 이 CIF 를 다시 읽어 대칭으로 펼친 자리(독립 구현), 산란 인자 = gemmi c4322 (같은 IT C 4.3.2.2 표)
    from pymatgen.io.cif import CifParser
    cifst = CifParser(p).parse_structures(primitive=False)[0]
    # 원자 B = 0 (CifWriter 는 Uiso 를 쓰지 않음)
    gmax = 6.0  # 1/nm
    rows = []
    hmax = [int(np.ceil(gmax * cell.a / 10)), int(np.ceil(gmax * cell.b / 10)), int(np.ceil(gmax * cell.c / 10))]
    for h in range(-hmax[0], hmax[0] + 1):
        for k in range(-hmax[1], hmax[1] + 1):
            for l in range(-hmax[2], hmax[2] + 1):
                if h == k == l == 0: continue
                d = cell.calculate_d([h, k, l]) / 10   # nm
                if 1 / d > gmax: continue
                absent = gops.is_systematically_absent([h, k, l])
                s2 = 1 / (4 * (d * 10) ** 2)
                F = sum(gemmi.Element(site.specie.symbol).c4322.calculate_sf(s2) * np.exp(2j * np.pi * np.dot([h, k, l], site.frac_coords)) for site in cifst)
                rows.append([h, k, l, d, bool(absent), abs(F) ** 2])
    # pymatgen 정대축 반사 (d, 위치) — [001] 과 [011]
    tem = TEMCalculator(voltage=200, beam_direction=(0, 0, 1))
    zones = {}
    for zb in [(0, 0, 1), (0, 1, 1), (1, 1, 1)]:
        tem = TEMCalculator(voltage=200, beam_direction=zb)
        pts = tem.zone_axis_filter(tem.generate_points(-4, 4), laue_zone=0)
        pts = [tuple(int(x) for x in q) for q in pts if any(q)]
        dsp = tem.get_interplanar_spacings(st, pts)
        zones[','.join(map(str, zb))] = sorted([[list(k), v / 10] for k, v in dsp.items() if v / 10 > 1 / gmax], key=lambda x: (-x[1], x[0]))[:40]
    ref[name] = dict(cif=name + '.cif', sg=sg.number, hm=sg.hm, rows=rows, zones=zones)
    print(name, sg.hm, sg.number, len(rows), 'refl', sum(1 for r in rows if not r[4]), 'allowed')
json.dump(ref, open(os.path.join(OUT, 'ref_core.json'), 'w'))
