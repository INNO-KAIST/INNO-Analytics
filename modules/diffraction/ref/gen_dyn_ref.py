"""동역학 회절(diff_dyn.js) 독립 기준 자료.
  (1) 독립 Bloch: gemmi(CIF·단위세포 원자·c4322 전자 산란 인자) + numpy/scipy — 빔 목록을 따로 세고, 흡수는 섭동 없이
      scipy.linalg.expm(2πi t H) 로 정확히 (고유값 풀이도 쓰지 않음).
  (2) 다중 슬라이스: abTEM 1.0 퍼텐셜 (IT C 4.3.2.2 계수 + DW 로 b_i+B, B = 0.5 Å², projection='finite' — 'infinite' 는 Si [011] 초격자에서 평균
      퍼텐셜이 13 % 작게 나와 쓰지 않음). 'zolz' = 한 반복의 투영 퍼텐셜을 8 등분 얇은 슬라이스로 numpy FFT 전파 (0 차 Laue 영역만),
      'abtem' = abTEM 기본 3D 슬라이스 ≈1 Å (HOLZ 포함). B = 0 이면 Au·Fe 는 8 Å⁻¹ 밖으로 수 % 가 흩어져 어느 방법도 수렴하지 않으므로 B = 0.5.
  사용: python3 -I ref/gen_dyn_ref.py [bloch|ms|all]   → test/ref_dyn.json
"""
import os, sys, json, time, numpy as np, gemmi, scipy.linalg as sla
from scipy import constants as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', 'test', 'ref_dyn.json')
PACK = os.path.join(HERE, 'pack')
CIF = dict(Si='Si_9011998.cif', Au='Au_9008463.cif', aFe='aFe_9008536.cif', ZnO='ZnO_2300112.cif', rutile='TiO2-rutile_9015662.cif',
           bNaYF4='bNaYF4_1517672.cif', Y2O3='Y2O3_1009013.cif', LiYF4='LiYF4_2001129.cif')

def lam_gamma(kV):
    E = kV * 1e3 * C.e; m = C.m_e; c = C.c
    lam = C.h / np.sqrt(2 * m * E * (1 + E / (2 * m * c * c)))
    return lam * 1e10, 1 + E / (m * c * c)        # Å, γ

def load(name, dw=True):
    ss = gemmi.read_small_structure(os.path.join(PACK, CIF[name]))
    sites = [(s.element, np.array(s.fract.tolist()), s.occ, (8 * np.pi ** 2 * s.u_iso) if dw else 0.0) for s in ss.get_all_unit_cell_sites()]
    sg = gemmi.find_spacegroup_by_name(ss.spacegroup_hm)
    cen = [np.array(c) / 24 for c in sg.operations().cen_ops]
    return ss.cell, sites, cen

def beams(cell, cen, uvw, gB_nm):
    gB = gB_nm / 10; out = []
    H = [int(np.ceil(gB * x)) + 1 for x in (cell.a, cell.b, cell.c)]
    for h in range(-H[0], H[0] + 1):
        for k in range(-H[1], H[1] + 1):
            for l in range(-H[2], H[2] + 1):
                if h * uvw[0] + k * uvw[1] + l * uvw[2]: continue
                g2 = cell.calculate_1_d2([h, k, l]) if (h or k or l) else 0.0
                if g2 > gB * gB * (1 + 1e-12): continue
                if any(abs(np.dot([h, k, l], c) - round(np.dot([h, k, l], c))) > 1e-6 for c in cen): continue
                out.append(((h, k, l), g2))
    out.sort(key=lambda x: (x[1], x[0]))
    return [b[0] for b in out]

def Ug(cell, sites, gam, h):
    s2 = (cell.calculate_1_d2(list(h)) if any(h) else 0.0) / 4
    F = sum(occ * el.c4322.calculate_sf(s2) * np.exp(-B * s2) * np.exp(2j * np.pi * np.dot(h, x)) for el, x, occ, B in sites)
    return gam * np.conj(F) / (np.pi * cell.volume)   # Å⁻²

def bloch(name, uvw, kV, gB, thick_nm, kappa=0.0, kt_nm=(0, 0, 0), dw=True):
    cell, sites, cen = load(name, dw)
    lam, gam = lam_gamma(kV); K0 = 1 / lam
    B = beams(cell, cen, uvw, gB); N = len(B)
    M = np.array(cell.frac.mat)              # frac = M · cart  →  g_cart = Mᵀ h
    gc = np.array([M.T @ np.array(h, float) for h in B])
    cache = {}
    def U(h):
        if h not in cache: cache[h] = Ug(cell, sites, gam, h)
        return cache[h]
    U0 = U((0, 0, 0)).real; K = np.sqrt(K0 * K0 + U0)
    kt = np.array(kt_nm) / 10
    s = np.array([-(2 * np.dot(kt, g) + np.dot(g, g)) / (2 * K) for g in gc])
    Hm = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j: Hm[i, i] = s[i] + 1j * kappa * U0 / (2 * K)
            else:
                d = tuple(B[i][q] - B[j][q] for q in range(3)); Hm[i, j] = (1 + 1j * kappa) * U(d) / (2 * K)
    e0 = np.zeros(N, complex); e0[B.index((0, 0, 0))] = 1
    I = []
    for t in thick_nm:
        psi = sla.expm(2j * np.pi * (t * 10) * Hm) @ e0; I.append(np.abs(psi) ** 2)
    # 고유값 검사 (흡수 없을 때 Hermite)
    herm = float(np.max(np.abs(Hm - Hm.conj().T))) if kappa == 0 else None
    return dict(phase=name, zone=list(uvw), kV=kV, gB=gB, kappa=kappa, kt=list(kt_nm), dw=dw, thickness=list(thick_nm), N=N,
                beams=[list(b) for b in B], I=np.array(I).T.tolist(), lam=lam / 10, gamma=gam, K=K * 10, U0=U0 * 100, hermErr=herm)

# ---------------- abTEM ----------------
SUPER = {  # 행 = 새 축 (관례 격자 벡터 정수 조합), 셋째 = 빔 방향, det > 0
    ('Si', (0, 1, 1)): [[1, 0, 0], [0, 1, -1], [0, 1, 1]],
    ('Si', (0, 0, 1)): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    ('Au', (0, 0, 1)): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    ('aFe', (1, 1, 1)): [[1, -1, 0], [1, 1, -2], [1, 1, 1]],
    ('ZnO', (1, 1, 0)): [[-1, 1, 0], [0, 0, 1], [1, 1, 0]],
    ('rutile', (0, 0, 1)): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
}

def supercell(name, uvw):
    from ase import Atoms
    cell, sites, cen = load(name, dw=False)
    Pm = np.array(SUPER[(name, tuple(uvw))]); assert round(np.linalg.det(Pm)) > 0
    A = np.array(cell.orth.mat).T            # 행 = a, b, c (Å, 직교 좌표)
    An = Pm @ A                              # 새 축
    n = round(abs(np.linalg.det(Pm)))
    # 새 셀 안 원자: 관례 셀을 넉넉히 반복하고 새 분율 좌표 [0,1) 인 것
    Ainv = np.linalg.inv(An); pos = []; sym = []
    R = int(np.ceil(np.abs(Pm).sum(axis=0).max())) + 1
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            for k in range(-R, R + 1):
                for el, x, occ, B in sites:
                    c = (x + np.array([i, j, k])) @ A; f = c @ Ainv
                    if np.all(f > -1e-7) and np.all(f < 1 - 1e-7):
                        pos.append(f); sym.append(el.name)
    assert len(pos) == n * len(sites), (len(pos), n, len(sites))
    L = np.linalg.norm(An, axis=1)
    for a in range(3):
        for b in range(a + 1, 3): assert abs(np.dot(An[a], An[b])) < 1e-8, '직교 아님'
    cart = np.array(pos) * L
    at = Atoms(symbols=sym, positions=cart, cell=np.diag(L), pbc=True)
    U = An / L[:, None]                      # 단위 새 축 (행)
    return at, cell, U, L

def peng_params(sites, B):
    """IT C 4.3.2.2 (gemmi c4322, diff_scat.js 와 같은 표) 5-가우스 계수에 Debye–Waller 를 b_i + B 로 넣은 abTEM Peng 매개변수"""
    out = {}
    for el, x, occ, _ in sites:
        c = el.c4322; out[el.name] = [list(c.a), [bb + B for bb in c.b]]
    return out

def multislice(name, uvw, kV, thick_nm, mode, B=0.5, sampling=0.04, tilt=(0, 0), gmax_nm=12.0, nsub=8):
    """mode 'zolz': abTEM 으로 한 반복 길이의 투영 퍼텐셜 V_p 를 만들고 (finite, 이웃 꼬리 포함), 그것을 nsub 등분한 얇은 슬라이스로
       직접 (numpy FFT) 전파 → 0 차 Laue 영역만의 동역학 (Δz → 0 이면 Bloch ZOLZ 와 같아야 함).
       mode 'abtem': abTEM 기본 다중 슬라이스 (≈1 Å 3D 슬라이스, HOLZ 포함). 둘 다 B 로 번진 정적 퍼텐셜 (흡수 없음)."""
    import abtem
    from abtem.parametrizations import PengParametrization
    from abtem.core.energy import energy2sigma, energy2wavelength
    abtem.config.set({'device': 'cpu'})
    at, cell, U, L = supercell(name, uvw)
    cell0, sites, cen = load(name, dw=False)
    par = PengParametrization(parameters=peng_params(sites, B))
    c = L[2]; M = np.array(cell0.frac.mat); Bm = beams(cell0, cen, uvw, gmax_nm)
    out = []
    lam = energy2wavelength(kV * 1e3); sig = energy2sigma(kV * 1e3)
    if mode == 'zolz':
        pot = abtem.Potential(at, sampling=sampling, slice_thickness=c, parametrization=par, projection='finite')
        Vp = pot.build().compute().array.sum(0)            # V·Å
        nx, ny = Vp.shape; kx = np.fft.fftfreq(nx, L[0] / nx); ky = np.fft.fftfreq(ny, L[1] / ny)
        K2 = kx[:, None] ** 2 + ky[None, :] ** 2
        kmax = min(np.abs(kx).max(), np.abs(ky).max()) * 2 / 3
        tx, ty = np.array(tilt) / 1000 / lam               # 기울임: 입사 파의 접선 파수
        dz = c / nsub
        prop = np.exp(-1j * np.pi * lam * dz * ((kx[:, None] + tx) ** 2 + (ky[None, :] + ty) ** 2)) * (K2 <= kmax ** 2)
        trans = np.exp(1j * sig * Vp / nsub)
        psi = np.ones((nx, ny), complex); z = 0; tdone = []
        targets = sorted(set(max(1, int(round(t * 10 / c))) for t in thick_nm)); res = {}
        for n in range(1, targets[-1] + 1):
            for _ in range(nsub):
                psi = np.fft.ifft2(np.fft.fft2(psi * trans) * prop)
            if n in targets: res[n] = np.fft.fft2(psi) / psi.size
        for t in thick_nm:
            n = max(1, int(round(t * 10 / c))); out.append((n * c / 10, res[n]))
    else:
        for t in thick_nm:
            nrep = max(1, int(round(t * 10 / c)))
            pot = abtem.Potential(at * (1, 1, nrep), sampling=sampling, slice_thickness=c / int(np.ceil(c / 1.0)), parametrization=par, projection='finite')
            ex = abtem.PlaneWave(energy=kV * 1e3, tilt=tilt).multislice(pot).compute().array
            out.append((nrep * c / 10, np.fft.fft2(ex) / ex.size))
    rows = []; tot = []
    for tt, Psi in out:
        I = []
        for h in Bm:
            g = M.T @ np.array(h, float); gx, gy = g @ U[0], g @ U[1]
            ix, iy = gx * L[0], gy * L[1]
            assert abs(ix - round(ix)) < 1e-6 and abs(iy - round(iy)) < 1e-6, (h, ix, iy)
            I.append(float(abs(Psi[int(round(ix)) % Psi.shape[0], int(round(iy)) % Psi.shape[1]]) ** 2))
        rows.append(I); tot.append(float(np.sum(np.abs(Psi) ** 2)))
    return dict(phase=name, zone=list(uvw), kV=kV, mode=mode, B=B, nsub=nsub, sampling=sampling, tilt=list(tilt), repeat=c / 10, beams=[list(b) for b in Bm],
                thickness=[o[0] for o in out], I=np.array(rows).T.tolist(), total=tot)

if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    ref = json.load(open(OUT)) if os.path.exists(OUT) else {}
    if what in ('bloch', 'all'):
        T = [10, 30, 60]; L = []
        t0 = time.time()
        L.append(bloch('Si', (0, 1, 1), 200, 20, T))
        L.append(bloch('Si', (0, 1, 1), 200, 20, T, kappa=0.1))
        L.append(bloch('Si', (0, 0, 1), 200, 20, T))
        L.append(bloch('Au', (0, 0, 1), 200, 30, T))
        L.append(bloch('Au', (0, 0, 1), 200, 30, T, kappa=0.1))
        L.append(bloch('aFe', (1, 1, 1), 200, 30, T))
        L.append(bloch('ZnO', (1, 1, 0), 200, 20, T))
        L.append(bloch('ZnO', (1, 1, 0), 200, 20, T, kappa=0.1))
        L.append(bloch('rutile', (0, 0, 1), 200, 15, T, kappa=0.05))
        L.append(bloch('bNaYF4', (0, 0, 1), 200, 14, T, kappa=0.1))
        L.append(bloch('Y2O3', (0, 0, 1), 200, 8, T))
        L.append(bloch('Si', (0, 1, 1), 200, 20, T, kt_nm=(0.9, -0.4, 0.4)))     # 기울임 (kt ⟂ [011])
        L.append(bloch('Au', (0, 0, 1), 300, 25, [5, 20, 45], kappa=0.08))
        ref['bloch'] = L
        print('bloch', len(L), 'cases', round(time.time() - t0, 1), 's', [c['N'] for c in L])
    if what in ('ms', 'all'):
        T = [10, 30, 60]; L = []; t0 = time.time()
        for nm, z in [('Si', (0, 1, 1)), ('Si', (0, 0, 1)), ('Au', (0, 0, 1)), ('aFe', (1, 1, 1)), ('ZnO', (1, 1, 0)), ('rutile', (0, 0, 1))]:
            for mode in ('zolz', 'abtem'):
                L.append(multislice(nm, z, 200, T, mode)); print(nm, z, mode, round(time.time() - t0, 1), 's', 'Σ|ψ|²', [round(x, 4) for x in L[-1]['total']], flush=True)
        L.append(multislice('Au', (0, 0, 1), 200, T, 'zolz', nsub=16)); print('Au nsub16', flush=True)          # 슬라이스 수렴 확인
        L.append(multislice('Au', (0, 0, 1), 200, T, 'zolz', sampling=0.025)); print('Au 0.025', flush=True)    # 표본화 수렴 확인
        L.append(multislice('rutile', (0, 0, 1), 200, [30], 'zolz', tilt=(5.0, 1.5)))                            # GM: 기울이면 (100) 이 생김
        ref['ms'] = L
    json.dump(ref, open(OUT, 'w'))
    print('→', OUT)
