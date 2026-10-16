"""실측 Olympus/Evident OIR 의 정답 — 공개 시료 (ref/real/oir, CC BY 4.0, SOURCE.txt).
세 독립 판독을 대조한다: oirfile(C. Gohlke, 2026.9.6) · Bio-Formats 7.3.1 (ref/bf/BFDump.java, 블랙박스) · 같은 출처의 .ome.tif 내보내기(있을 때, tifffile 로 읽음 —
내보내기는 Bio-Formats 5.2.0-DEV 가 만든 것이라 7.3.1 과 완전히 독립은 아님).
정답은 화소가 서로 같은 것만. 다르면 truth.json 의 'decisions' 에 어느 쪽을 따르는지와 까닭을 적는다.
화소 정답은 크기가 커서(≈ 460 MB) 면마다 uint16 LE 바이트의 SHA-1 로 저장한다 (값 비교와 같은 엄밀함).
앱 배치: 면 = c + nc·(z + nz·t'), t' = λ + nλ·t (λ 를 T 에 잇는다).
사용: python3 -I ref/gen_oir_real.py <OIR 폴더> <bioformats_package.jar> <BFDump 컴파일 폴더>"""
import sys, os, json, glob, subprocess, tempfile, hashlib, re, shutil
import numpy as np, oirfile, tifffile

D, JAR, OUT = sys.argv[1:4]
T = os.path.join(D, 'truth'); os.makedirs(T, exist_ok=True)


def bf(path):
    td = tempfile.mkdtemp()
    r = subprocess.run(['java', '-Xmx8g', '-cp', JAR + ':' + OUT, 'BFDump', path, td, '0'], capture_output=True, text=True, env=dict(os.environ, JAVA_TOOL_OPTIONS=''))
    S = None
    for l in r.stdout.splitlines():
        if l.startswith('S 0 '): t = l.split(); S = dict(zip(t[2::2], t[3::2]))
    if S is None:
        shutil.rmtree(td); return None, ([l for l in r.stderr.splitlines() if 'Exception' in l] or ['?'])[0][:160]
    X, Y, Z, C, Tt = (int(S[k]) for k in ['X', 'Y', 'Z', 'C', 'T'])
    S['arr'] = np.fromfile(os.path.join(td, 's0.f32'), '<f4').reshape(Tt, Z, C, Y, X); shutil.rmtree(td)
    return S, None


def canon(o):
    """oirfile → (T, L, Z, C, Y, X)"""
    a = o.asarray(); dims = list(o.dims)
    if 'S' in dims: dims[dims.index('S')] = 'C'
    for ax in 'TLZC':
        if ax not in dims: a = a[np.newaxis]; dims.insert(0, ax)
    return np.transpose(a, [dims.index(x) for x in 'TLZCYX'])


def sha(p): return hashlib.sha1(np.ascontiguousarray(p, '<u2').tobytes()).hexdigest()


def rel(a, b): return a is not None and b is not None and abs(a - b) <= 1e-6 * max(abs(a), abs(b))


truth, table = {}, []
for p in sorted(glob.glob(os.path.join(D, '*.oir'))):
    fn = os.path.basename(p); rec = {'file': fn}; dec = []
    with oirfile.OirFile(p) as o:
        sizes = dict(o.sizes); sc = dict(o.coord_scales)
        guids = o._get_ordered_channel_guids() if o._pixel_map else []
        names = [next((c.name for c in o.channels if c.id == g), g) for g in guids]
        a = canon(o) if sizes else None
        ref = o.reference.asarray() if o.reference is not None else None
        ax = {k: dict(v) for k, v in o._axis_info.items()}
    S, err = bf(p)
    if a is None:
        rec.update(error='지도(Map)', oirfile='크기 없음 (화소 블록 0 — FRAMEPROPERTIES 만, 45858×30169 지도 설정)', bf=err or '읽음',
                   decisions=['화소 블록이 없는 지도(Map) 파일 — oirfile 은 빈 배열, Bio-Formats 는 오류 → 앱은 분명한 오류 문구'])
        truth[fn] = rec; table.append((fn, '-', 'X', 'err' if err else '?', '-', '오류(화소 없음)')); print(fn, '화소 없음', err); continue
    Tn, L, Z, C, Y, X = a.shape
    app = a.reshape(Tn * L, Z, C, Y, X)  # t' = t·L + λ
    # Bio-Formats: λ 를 C 로 (C = nc·nλ). 같은 화소인지 λ·채널 배치를 맞춰 비교
    bf_eq = None
    if S is not None:
        b = S['arr']
        if L > 1: cand = a.transpose(0, 2, 1, 3, 4, 5).reshape(Tn, Z, L * C, Y, X)
        else: cand = a.reshape(Tn, Z, C, Y, X)
        bf_eq = b.shape == cand.shape and np.array_equal(b, cand.astype(np.float32))
        if not bf_eq and b.shape == cand.shape: rec['bf_diff_planes'] = int(sum(not np.array_equal(b[t, z, c], cand[t, z, c]) for t in range(Tn) for z in range(Z) for c in range(b.shape[2])))
    ome = os.path.join(D, fn[:-4] + '.ome.tif'); ome_eq = None
    if os.path.exists(ome):
        with tifffile.TiffFile(ome) as t:
            s = t.series[0]; oa = s.asarray(); om = t.ome_metadata
        pa = re.search(r'<Pixels [^>]*>', om).group(0); g = lambda k: (re.search(r"\b" + k + r'="([^"]*)"', pa) or [None, None])[1]
        oax = s.axes; full = oa
        for axn in 'TZC':
            if axn not in oax: full = full[np.newaxis]; oax = axn + oax
        full = np.transpose(full, [oax.index(x) for x in 'TZCYX'])
        ome_eq = full.shape == app.shape and np.array_equal(full, app)
        rec['ome'] = dict(axes=s.axes, shape=list(s.shape), px=float(g('PhysicalSizeX')), pz=float(g('PhysicalSizeZ')) if g('PhysicalSizeZ') else None, sizeZ=int(g('SizeZ')), sizeT=int(g('SizeT')), creator=(re.search(r'Creator="([^"]*)"', om) or [None, ''])[1], equal=bool(ome_eq))
    o_px = sc.get('X'); b_px = float(S['px']) if S and S['px'] != '-' else None
    o_pz = sc.get('Z') if Z > 1 else None; b_pz = float(S['pz']) if S and S['pz'] != '-' and Z > 1 else None
    rec.update(w=X, h=Y, nc=C, nz=Z, nt=Tn * L, nt0=Tn, nl=L, bit=16, px=o_px, chNames=names, oirfile_sizes=sizes, oirfile_scales={k: float(v) for k, v in sc.items()},
               bf=(err or {k: S[k] for k in ['X', 'Y', 'Z', 'C', 'T', 'type', 'px', 'pz', 'ch']}), bf_equal=bf_eq, ome_equal=ome_eq,
               axes_oirfile=ax)
    if S and not rel(o_px, b_px): dec.append(f'화소 크기 oirfile {o_px} / Bio-Formats {b_px}')
    if Z > 1:
        rec['zStep'] = o_pz if rel(o_pz, b_pz) else b_pz
        if not rel(o_pz, b_pz): dec.append(f'z step oirfile {o_pz} / Bio-Formats {b_pz} → Bio-Formats (lsmimage length z)')
    else: rec['zStep'] = 0
    if S:
        bch = S['ch'].split('|'); want = names if L == 1 else names
        if L == 1 and bch != [n.replace(' ', '_') for n in names]: dec.append(f'채널 이름 oirfile {names} / Bio-Formats {bch}')
    if L > 1 and S is not None:
        b = S['arr']; rep = b.shape[2] == L * C and all(np.array_equal(b[:, :, k], a[:, 0, :, 0].astype(np.float32)) for k in range(L))
        rec['bf_lambda'] = 'Bio-Formats 는 λ 자리 모두에 λ1 면을 되풀이 (그 λ1 면은 oirfile 과 같음)' if rep else 'Bio-Formats λ 배치 미상'
        # λ 번호 ↔ 파장: FRAMEPROPERTIES 이름(lNNN…)과 LAMBDA 위치가 축 설정 start + (l−1)·step 과 맞는지 (oirfile 과 별개로 파일 자체에서)
        with oirfile.OirFile(p) as o2: fps = o2.xml_metadata.get(oirfile.METADATA.FRAMEPROPERTIES, [])
        ok_l = []
        for x in fps:
            mi = re.search(r'id="l(\d+)', x); mp = re.search(r'<commonframe:axisType>LAMBDA</commonframe:axisType>\s*<commonframe:position>([\d.]+)<', x)
            if mi and mp: ok_l.append(abs(float(mp.group(1)) - (ax['LAMBDA']['start'] + (int(mi.group(1)) - 1) * ax['LAMBDA']['step'])) < 1e-6)
        rec['lambda_check'] = f'면 {len(ok_l)} 개 중 {sum(ok_l)} 개의 λ 위치 = start + (l−1)·step'
    if L > 1: dec.append(f'λ {L}: oirfile = L 축, Bio-Formats = C {L * C} 개 (화소 {"같음" if bf_eq else "다름"}) → 앱은 λ 를 T 에 잇는다 (LIF 와 같음)')
    if 'T' in sizes and Tn > 1: rec['fi'] = float(sc['T']) if 'T' in sc else None
    zax = ax.get('ZSTACK'); rec['axisWarn'] = bool(zax and zax.get('maxSize', 0) > Z and Z > 1)  # oirfile 이 고른 축(마지막 정의) 기준 — 참고용
    if ref is not None:
        r3 = ref.reshape(-1, ref.shape[-2], ref.shape[-1]); rec['ref'] = dict(w=int(r3.shape[2]), h=int(r3.shape[1]), nc=int(r3.shape[0]), sha=[sha(x) for x in r3], src='oirfile 만 (Bio-Formats 는 따로 내지 않음)')
    rec['sha'] = [sha(app[t, z, c]) for t in range(Tn * L) for z in range(Z) for c in range(C)]
    agree = bool(bf_eq) and (ome_eq is not False)
    rec['truth'] = ('oirfile = Bio-Formats' + (' = ome.tif' if ome_eq else '')) if agree else ('oirfile' + (' = ome.tif' if ome_eq else '') + f' (Bio-Formats {"다름" if S else "오류"})')
    if dec: rec['decisions'] = dec
    truth[fn] = rec
    table.append((fn, f'{X}×{Y} C{C} Z{Z} T{Tn}' + (f' λ{L}' if L > 1 else ''), 'ok', 'ok' if bf_eq else ('err' if err else 'X'), '-' if ome_eq is None else ('ok' if ome_eq else 'X'), rec['truth']))
    print(fn, table[-1][1:], dec)
json.dump(truth, open(os.path.join(T, 'truth.json'), 'w'), ensure_ascii=False, indent=1)
print(f'\n{"파일":52s} {"크기":24s} oirfile BF  ome.tif  정답')
for r in table: print(f'{r[0][:52]:52s} {r[1]:24s} {r[2]:7s} {r[3]:4s} {r[4]:8s} {r[5]}')
