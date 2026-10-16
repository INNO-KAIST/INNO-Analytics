"""합성 Olympus/Evident OIR (FluoView FV3000/FV4000) + 정답 — test/io/oir/ 에 쓴다.
공개된 실측 .oir 시료가 없어 합성 파일로만 시험한다. 정답 = oirfile(C. Gohlke) 과 Bio-Formats 7.3.1 이 같은 값을 읽은 경우만.
· 파일 구조(블록 index·UID/PIXEL 짝·METADATA/FRAMEPROPERTIES XML)는 oirfile 이 읽는 대로, 블록을 놓는 차례와 XML 요소 이름은
  Bio-Formats 를 블랙박스로 시험해(어느 바이트를 읽는지·어느 XML 요소를 찾는지 관찰) 두 구현이 모두 읽도록 맞췄다 — Bio-Formats 원본 코드는 보지 않았다.
· 둘이 다르게 읽는 경우는 synth.json 의 'disagree' 에 적고, 앱이 따를 쪽을 'follow' 에 적는다.
· Bio-Formats 7.3.1 블랙박스 관찰 (I/O 위치 추적·XML DOM 질의 기록·공개 API 로만):
  - 16 바이트부터 u32 -1 표지를 찾고 4 바이트를 건너뛴 뒤 블록을 차례로 읽는다 (index 는 쓰지 않음) → 표지를 머리 뒤에 둔다.
  - XML 은 '<?xml' 을 찾아 앞의 u32 를 길이로 읽는다. 바로 읽은 FRAMEPROPERTIES 뒤에는 u32 하나를 더 건너뜀 → 첫 FRAMEPROPERTIES 뒤 4 바이트 틈.
  - 화소 블록 사이에 FRAMEPROPERTIES 가 끼면 다음 UID 를 놓친다 → 둘째 면부터의 FRAMEPROPERTIES 는 화소 뒤에 모은다.
  - lsmimage 에서 commonimage:imageInfo/width·height, commonphase:channel(id·order·name·imageDefinition depth·bitCounts·length x·y·z),
    commonimage:axis(commonparam:axis·maxSize·startPosition·step) 를 찾는다; FRAMEPROPERTIES 는 lsmframe:frameProperties/commonframe:imageDefinition/base:width…
    와 commonframe:general/base:name(없으면 오류). commonimage:acquisition 이 있으면 채널 수가 0 이 되어(필요 요소 미상) 넣지 않았다.
  - 실측 파일의 블록 차례는 알 수 없다 — 위 배치는 두 구현이 모두 읽는 한 가지일 뿐이다.
사용: python3 -I ref/gen_oir_synth.py <bioformats_package.jar> <BFDump 컴파일 폴더>"""
import sys, os, json, struct, subprocess, tempfile, shutil
import numpy as np
import oirfile

JAR, OUT = sys.argv[1:3]
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test', 'io', 'oir'); os.makedirs(D, exist_ok=True)
rng = np.random.default_rng(20261016)

NS = ' '.join(f'xmlns:{p}="http://www.olympus.co.jp/hpf/model/{p}"' for p in
              ['lsmimage', 'commonimage', 'commonphase', 'commonparam', 'base', 'opticalelement', 'fvCommonphase', 'lsmparam', 'fvCommonparam'])
I = lambda *v: struct.pack('<' + 'i' * len(v), *v)
AX = {'Z': 'ZSTACK', 'T': 'TIMELAPSE', 'L': 'LAMBDA'}


def guid(k):
    return f'{0x1a2b3c4d + 0x111 * k:08x}-5e6f-4a1b-9c2d-{0x3e4f5a6b7c8d + k:012x}'


def lsm_xml(sp):
    """lsmimage:imageProperties — 영상 크기·채널·축·화소 크기 (Bio-Formats 와 oirfile 이 모두 찾는 요소)"""
    W, H = sp['W'], sp['H']
    axes = ''
    for a in 'ZTL':
        n = sp.get('maxSize', {}).get(a, sp['n' + a.lower()])
        if n > 1 or a in sp.get('axes_enabled', ''):
            st, step = sp.get(a + '0', 0.0), sp.get(a + 'step', 1.0)
            axes += (f'<commonimage:axis enable="true"><commonparam:axis>{AX[a]}</commonparam:axis><commonparam:maxSize>{n}</commonparam:maxSize>'
                     f'<commonparam:startPosition>{st}</commonparam:startPosition><commonparam:endPosition>{st + (n - 1) * step}</commonparam:endPosition>'
                     f'<commonparam:step>{step}</commonparam:step></commonimage:axis>')
    chs = ''
    for c in range(sp['nc']):
        nm = sp['names'][c]; wl = sp.get('wl', [(None, None)] * sp['nc'])[c]
        zl = f'<commonparam:z>{sp["Zstep"]}</commonparam:z>' if sp['nz'] > 1 and sp.get('Zstep') else ''
        wls = (f'<opticalelement:startWavelength>{wl[0]}</opticalelement:startWavelength><opticalelement:endWavelength>{wl[1]}</opticalelement:endWavelength>' if wl[0] else '')
        chs += (f'<commonphase:channel id="{sp["guids"][c]}" order="{sp["order"][c]}"><commonphase:name>{nm}</commonphase:name>'
                f'<commonphase:imageDefinition><commonphase:depth>{sp["depth"]}</commonphase:depth><commonphase:bitCounts>{sp["bits"]}</commonphase:bitCounts></commonphase:imageDefinition>'
                f'<commonphase:length><commonparam:x>{sp["px"]}</commonparam:x><commonparam:y>{sp["py"]}</commonparam:y>{zl}</commonphase:length>{wls}'
                f'<fvCommonphase:pinholeDiameter>100.0</fvCommonphase:pinholeDiameter></commonphase:channel>')
    return (f'<?xml version="1.0" encoding="UTF-8"?><lsmimage:imageProperties {NS}>'
            f'<commonimage:general><base:creationDateTime>2026-09-17T14:03:21.512+09:00</base:creationDateTime></commonimage:general>'
            f'<commonimage:imageInfo><commonimage:width>{W}</commonimage:width><commonimage:height>{H}</commonimage:height>{axes}'
            f'<commonphase:phase><commonphase:group>{chs}</commonphase:group></commonphase:phase></commonimage:imageInfo>'
            f'<commonimage:length><commonparam:x>{sp["px"]}</commonparam:x><commonparam:y>{sp["py"]}</commonparam:y></commonimage:length>'
            f'</lsmimage:imageProperties>').encode()


def fileinfo_xml():
    return ('<?xml version="1.0" encoding="UTF-8"?><fileinfo:fileInfo xmlns:fileinfo="http://www.olympus.co.jp/hpf/model/fileinfo" '
            'xmlns:base="http://www.olympus.co.jp/hpf/model/base"><base:version>2.1.2.0</base:version><fileinfo:system>FV31S-SW</fileinfo:system></fileinfo:fileInfo>').encode()


def lut_xml(g):
    return (f'<?xml version="1.0" encoding="UTF-8"?><lut:LUT xmlns:lut="http://www.olympus.co.jp/hpf/model/lut"><lut:channelId>{g}</lut:channelId>'
            f'<lut:linear>true</lut:linear><lut:data></lut:data></lut:LUT>').encode()


def frame_xml(sp, pos):
    ax = ''.join(f'<commonframe:axisValue><commonframe:axisType>{t}</commonframe:axisType><commonframe:position>{p}</commonframe:position></commonframe:axisValue>' for t, p in pos)
    return (f'<?xml version="1.0" encoding="UTF-8"?><lsmframe:frameProperties xmlns:lsmframe="http://www.olympus.co.jp/hpf/model/lsmframe" '
            f'xmlns:commonframe="http://www.olympus.co.jp/hpf/model/commonframe" xmlns:base="http://www.olympus.co.jp/hpf/model/base">'
            f'<commonframe:general><base:name>{sp.get("base", "Image")}</base:name></commonframe:general><commonframe:imageDefinition><base:width>{sp["W"]}</base:width>'
            f'<base:height>{sp["H"]}</base:height><base:depth>{sp["depth"]}</base:depth><base:bitCounts>{sp["bits"]}</base:bitCounts>'
            f'<base:colorType>{sp.get("colortype", "GlayScale")}</base:colorType></commonframe:imageDefinition>{ax}</lsmframe:frameProperties>').encode()


def uid(sp, t, l, z, c, b):
    s = (f'l{l + 1:03d}' if sp['nl'] > 1 else '') + (f'z{z + 1:03d}' if sp['nz'] > 1 else '') + (f't{t + 1:03d}' if sp['nt'] > 1 else '')
    return f'{s}_0_1_{sp["guids"][c]}_{b}'.encode()


def spec(**k):
    sp = dict(W=8, H=6, nc=1, nz=1, nt=1, nl=1, depth=2, bits=12, px=0.31, py=0.31, blocks=1, frames=True)
    sp.update(k)
    sp.setdefault('names', [f'CH{c + 1}' for c in range(sp['nc'])])
    sp.setdefault('guids', [guid(c) for c in range(sp['nc'])])
    sp.setdefault('order', [c + 1 for c in range(sp['nc'])])
    return sp


def imgdef_xml(w, h):
    return (f'<?xml version="1.0" encoding="UTF-8"?><base:imageDefinition xmlns:base="http://www.olympus.co.jp/hpf/model/base">'
            f'<base:width>{w}</base:width><base:height>{h}</base:height></base:imageDefinition>').encode()


def xml_block(t, *xs):
    c = b''.join(I(len(x)) + x for x in xs); return I(len(c), t) + c


def pix_pair(u, p):
    return [(I(len(u) + 12, 3, 0, 0, len(u)) + u, b''), (I(len(p), 4) + p, b'')]


def assemble(blocks):
    out = bytearray(48) + I(-1, 0); offs = []  # -1 표지 + u32 (Bio-Formats 는 16 바이트부터 -1 을 찾고 4 바이트 건너뜀)
    for b, gap in blocks:
        offs.append(len(out)); out += b + gap
    idx = len(out); out += I(-1) + b''.join(struct.pack('<Q', o) for o in offs) + struct.pack('<Q', 0)
    out[0:48] = b'OLYMPUSRAWFORMAT' + struct.pack('<IIIIQQ', 12, 0, 1, 2, len(out), idx)
    return bytes(out)


def build(sp, split_t=None):
    """블록 배치: 머리 48 | -1 표지 | 첫 면 FRAMEPROPERTIES + 틈 | [METADATA (meta_early)] | UID PIXEL… | 나머지 FRAMEPROPERTIES | METADATA | block index
    split_t: 이 T 부터의 면은 보조 파일(이름_00001)로 — (본 파일, 보조 파일, 자료)"""
    W, H, nc, nz, nt, nl, dep = (sp[k] for k in ['W', 'H', 'nc', 'nz', 'nt', 'nl', 'depth'])
    dt = {1: '<u1', 2: '<u2', 4: '<f4'}[dep]
    Hd = sp.get('lines', H)  # 선 주사: 면 자료의 줄 수 (FRAMEPROPERTIES height 보다 많음)
    mx = (1 << sp['bits']) - 1 if dep < 4 else 1000.0
    data = (rng.random((nt, nl, nz, nc, Hd, W)) * 1e3).astype(dt) if dep == 4 else rng.integers(0, mx + 1, (nt, nl, nz, nc, Hd, W)).astype(dt)
    blocks, later, comp = [], [], []  # (블록 바이트, 블록 뒤 틈) — FRAMEPROPERTIES 뒤에는 4 바이트 틈 (Bio-Formats 가 XML 다음에 u32 하나를 건너뜀)
    skip = set(sp.get('drop', []))
    meta = (xml_block(0, fileinfo_xml(), lsm_xml(sp), *[lut_xml(g) for g in sp['guids']], *([imgdef_xml(*sp['ref'])] if sp.get('ref') else [])), b'')
    for t in range(nt):
        dst = comp if split_t is not None and t >= split_t else blocks
        for l in range(nl):
            for z in range(nz):
                if sp['frames'] or (t, l, z) == (0, 0, 0):
                    pos = []
                    if nz > 1: pos.append(('ZSTACK', round(sp.get('Z0', 0.0) + z * sp.get('Zstep', 1.0), 6)))
                    if nt > 1: pos.append(('TIMELAPSE', round((sp.get('T0', 0.0) + t * sp.get('Tstep', 1.0)) * 1000, 3)))
                    if nl > 1: pos.append(('LAMBDA', round(sp.get('L0', 0.0) + l * sp.get('Lstep', 1.0), 3)))
                    fb = (xml_block(1, frame_xml(sp, pos)), I(0))
                    if (t, l, z) == (0, 0, 0):
                        blocks.append(fb)
                        if sp.get('meta_early'): blocks.append(meta)
                    else: (comp if dst is comp else later).append(fb)
                for c in range(nc):
                    raw = data[t, l, z, c].tobytes(); nb = sp['blocks']
                    rows = [Hd * k // nb for k in range(nb + 1)]
                    for b in range(nb):
                        if (t, l, z, c, b) in skip: continue
                        dst += pix_pair(uid(sp, t, l, z, c, b), raw[rows[b] * W * dep:rows[b + 1] * W * dep])
    if sp.get('ref'):  # 참조(개관) 영상 — REF_LSM0_<GUID>_<블록>
        rw, rh = sp['ref']; sp['refdata'] = rng.integers(0, 4096, (nc, rh, rw)).astype(dt)
        for c in range(nc): blocks += pix_pair(f'REF_LSM0_{sp["guids"][c]}_0'.encode(), sp['refdata'][c].tobytes())
    blocks += later  # 둘째 면부터의 FRAMEPROPERTIES 는 화소 뒤에 모아 둔다 (화소 블록 사이에 두면 Bio-Formats 가 다음 UID 를 놓침)
    if not sp.get('meta_early'): blocks.append(meta)
    if split_t is not None: return assemble(blocks), assemble(comp), data
    return assemble(blocks), data


def bfdump(fp):
    td = tempfile.mkdtemp()
    r = subprocess.run(['java', '-Xmx2g', '-cp', JAR + ':' + OUT, 'BFDump', fp, td, '0'], capture_output=True, text=True, env=dict(os.environ, JAVA_TOOL_OPTIONS=''))
    S = None
    for l in r.stdout.splitlines():
        if l.startswith('S 0 '): t = l.split(); S = dict(zip(t[2::2], t[3::2]))
    if S is None: return None, ([l for l in r.stderr.strip().splitlines() if not l.startswith('Picked up')] or ['?'])[0][:200]
    X, Y, Z, C, T = (int(S[k]) for k in ['X', 'Y', 'Z', 'C', 'T']); rgb = int(S['rgb'])
    a = np.fromfile(os.path.join(td, 's0.f32'), '<f4'); shutil.rmtree(td)
    S['arr'] = a.reshape(T, Z, C * rgb if rgb > 1 else C, Y, X) if a.size == T * Z * C * Y * X else a
    return S, None


def oir_canon(fp):
    """oirfile → (T, L, Z, C, Y, X) 배열·크기·눈금·채널 이름"""
    with oirfile.OirFile(fp) as o:
        a = o.asarray(); dims = list(o.dims); sizes = dict(o.sizes); sc = dict(o.coord_scales)
        names = [str(x) for x in o.coords.get('C', o.coords.get('S', []))]
        if o.is_linescan: dims[dims.index('T')] = 'Y'
    if 'S' in dims: dims[dims.index('S')] = 'C'
    for ax in 'TLZC':
        if ax not in dims: a = a[np.newaxis]; dims.insert(0, ax)
    a = np.transpose(a, [dims.index(x) for x in 'TLZCYX'])
    return a.astype(np.float32), sizes, sc, names


def to_app(a):
    """(T, L, Z, C, Y, X) → 앱 면 순서 c + nc·(z + nz·t'), t' = λ + nλ·t (λ 를 T 에 잇는다)"""
    T, L, Z, C, Y, X = a.shape
    return np.ascontiguousarray(a.transpose(0, 1, 2, 3, 4, 5).reshape(T * L, Z, C, Y, X)).astype(np.float32)


def close(a, b):
    if a is None or b is None: return a is None and b is None
    return abs(a - b) <= 1e-9 * max(abs(a), abs(b))


# (파일, spec, 설명, 기대) — 기대: Bio-Formats 와 다를 때 따를 쪽과 까닭, 알림·오류 문구 등
CASES = [
    ('gray12.oir', spec(), 'C1 12-bit (16-bit 저장)', {}),
    ('gray8_odd.oir', spec(W=7, H=5, depth=1, bits=8, px=0.124, py=0.124), '8-bit, 홀수 너비 7', {}),
    ('gray16.oir', spec(W=9, H=4, bits=16, px=1.2427, py=1.2427), '16-bit', {}),
    ('c3.oir', spec(nc=3, names=['DAPI', 'Alexa 488', 'Alexa 594'], wl=[(430, 470), (500, 540), (610, 680)], order=[2, 1, 3]), '3 채널 (order 속성이 GUID 순서와 다름)', {}),
    ('z5.oir', spec(nz=5, Z0=12.5, Zstep=0.37), 'Z 5', {}),
    ('t4.oir', spec(nt=4, Tstep=2.5), 'T 4', {}),
    ('c2z3t2.oir', spec(nc=2, nz=3, nt=2, Z0=3.0, Zstep=0.5, Tstep=1.5), 'C2·Z3·T2', {}),
    ('blocks3.oir', spec(W=10, H=9, blocks=3, nz=2, Zstep=2.0), '면마다 PIXEL 블록 3 개', {}),
    ('wide.oir', spec(W=515, H=3, px=0.0621, py=0.0621), '너비 515', {}),
    ('c2z4big.oir', spec(W=64, H=48, nc=2, nz=4, nt=3, Zstep=1.25, Tstep=10), 'C2·Z4·T3 64×48 (Z·T 고르기 시험)', {'choose': {'z': [1, 3], 't': [0, 2]}}),
    ('rgb8.oir', spec(nc=3, depth=1, bits=8, colortype='RGB', names=['R', 'G', 'B'], W=6, H=4), 'RGB 8-bit (colorType RGB, 표본별 GUID)',
     {'follow': 'oirfile', 'why': 'Bio-Formats 는 colorType 을 보지 않아 3 채널 회색으로 읽음 — 화소는 같고, 앱은 oirfile 처럼 RGB 로 연다', 'rgb': True, 'pixels_only': True}),
    ('noframes.oir', spec(nz=3, Zstep=0.8, frames=False), '면마다 FRAMEPROPERTIES 없음 (위치는 첫 면 하나)',
     {'follow': 'Bio-Formats', 'why': 'oirfile 은 면 위치가 하나뿐이면 Z 간격을 버리고 축 설정 step 도 안 씀; Bio-Formats 는 0.8 µm — 앱은 축 설정 step 을 쓴다', 'zStep': 0.8}),
    ('miss_mid.oir', spec(nz=4, Zstep=0.5, drop=[(0, 0, 2, 0, 0)]), 'Z 4 중 가운데(z3) PIXEL 블록 없음',
     {'follow': 'Bio-Formats', 'why': 'oirfile 은 있는 Z 번호만 모아 Z 3 으로 당김(z4 가 z3 자리로); Bio-Formats 는 Z 4, 빠진 면 = 0 — 앱은 번호 자리를 지킨다', 'warn': '0 으로 채웠습니다'}),
    ('miss_last.oir', spec(nz=4, Zstep=0.5, drop=[(0, 0, 3, 0, 0)]), 'Z 4 중 마지막 면 없음 (중단된 획득)', {'warn': '획득 설정보다 면이 적습니다'}),
    ('lambda4.oir', spec(nl=4, L0=500, Lstep=10), 'λ 4',
     {'follow': 'oirfile', 'why': 'Bio-Formats 7.3.1 은 λ 를 채널로 세지만 화소는 모두 λ1 면을 되풀이 — 앱은 oirfile 의 L 축을 T 에 잇는다'}),
    ('c2l3.oir', spec(nc=2, nl=3, L0=500, Lstep=10), 'C2·λ3', {'follow': 'oirfile', 'why': '위와 같음'}),
    ('l3z2.oir', spec(nl=3, nz=2, L0=480, Lstep=20, Zstep=1.0), 'λ3·Z2', {'follow': 'oirfile', 'why': '위와 같음'}),
    ('f32.oir', spec(depth=4, bits=32), '깊이 4 바이트 (float 자료)',
     {'follow': 'oirfile', 'why': 'oirfile 은 깊이 4 = 32-bit float, Bio-Formats 는 uint32 로 읽음 (실측으로 확인 못 함) — 앱은 float 로 봐서 값이 그럴듯하면 float, 대부분 비정규수면 uint32 (p9_u32 참고)', 'warn': '깊이 4 바이트'}),
    ('linescan.oir', spec(W=16, H=1, lines=40, blocks=4, axes_enabled='T', maxSize={'T': 40}), '선 주사 (면 높이 1, 40 줄을 블록 4 개로)',
     {'follow': 'oirfile', 'why': 'Bio-Formats 7.3.1 은 Y 1 줄만 읽음; oirfile 은 PIXEL 자료 길이로 40 줄 (kymograph) — 앱은 oirfile', 'warn': '선 주사'}),
    ('ref.oir', spec(nc=2, ref=(12, 10)), '참조(개관) 영상 REF_ 블록 포함',
     {'follow': 'oirfile', 'why': 'Bio-Formats 7.3.1 은 REF_ 블록을 본 영상 채널 자료로 덮어씀; oirfile 은 따로 둠 — 앱은 본 영상 + 참조 영상 2 개', 'ref': True, 'warn': '참조 영상'}),
]


def write_truth(key, arr):
    arr.astype('<f4').tofile(os.path.join(D, key + '.f32'))


def review_probes(out, table):
    """독립 검토(2026-10-16)에서 나온 경계 사례 — 정답 = 만든 자료를 앱 규칙대로 놓은 것. oirfile·Bio-Formats 결과는 함께 적는다."""
    r = np.random.default_rng(1)
    def frame(sp, pos=()): return (xml_block(1, frame_xml(sp, list(pos))), I(0))
    def meta(sp, lsm=None, extra=()): return (xml_block(0, fileinfo_xml(), lsm if lsm is not None else lsm_xml(sp), *[lut_xml(x) for x in sp['guids']], *extra), b'')
    def run(fn, blocks, app, desc, rec, asm=assemble):
        fp = os.path.join(D, fn); open(fp, 'wb').write(asm(blocks)); write_truth(fn, np.asarray(app, np.float32).reshape(-1))
        a = np.asarray(app, np.float32)
        try:
            with oirfile.OirFile(fp) as o: x = o.asarray(); oe = f'{dict(o.sizes)} {o.dtype}' + (' = 만든 자료' if x.size == a.size and np.array_equal(np.sort(x.astype(np.float32).ravel()), np.sort(a.ravel())) and np.array_equal(x.astype(np.float32).ravel(), a.ravel()) else ' ≠ 정답')
        except Exception as e: oe = f'오류 {type(e).__name__}'
        S, err = bfdump(fp)
        be = f'오류 ({err[:80]})' if err else f'X{S["X"]} Y{S["Y"]} Z{S["Z"]} C{S["C"]} T{S["T"]} {S["type"]} pz {S["pz"]} ch {S["ch"]}' + (' = 정답' if S['arr'].size == a.size and np.array_equal(S['arr'].ravel(), a.ravel()) else ' ≠ 정답')
        rec.update(desc=desc, key=fn, truth=f'만든 자료 (앱 규칙) · oirfile: {oe} · Bio-Formats: {be}'); rec.setdefault('px', 0.31); rec.setdefault('zStep', 0); rec.setdefault('bit', 16)
        out[fn] = rec; table.append((fn, f'{rec["w"]}×{rec["h"]} C{rec["nc"]} Z{rec["nz"]} T{rec["nt"]}', oe[-12:], be[-12:], desc[:60]))
    sp = spec(W=8, H=6, nc=2); d = r.integers(0, 4096, (2, 6, 8)).astype('<u2')
    # p1: 채널끼리 엇갈린 블록 + 파일 전체 블록 번호 (c0: 0,2,4 · c1: 1,3,5) → block index 차례대로 잇기 (oirfile·Bio-Formats 와 같음)
    bl = [frame(sp)]; k = 0
    for b in range(3):
        for c in range(2): bl += pix_pair(f'_0_1_{sp["guids"][c]}_{k}'.encode(), d[c, 2 * b:2 * b + 2].tobytes()); k += 1
    run('p1_globalctr.oir', bl + [meta(sp)], d, '채널끼리 엇갈린 PIXEL 블록, 번호는 파일 전체에서 셈', dict(w=8, h=6, nc=2, nz=1, nt=1, chNames=['CH1', 'CH2']))
    bl = [frame(sp)]
    for b in range(3):
        for c in range(2): bl += pix_pair(f'_0_1_{sp["guids"][c]}_{b + 1}'.encode(), d[c, 2 * b:2 * b + 2].tobytes())
    run('p1b_interleave1.oir', bl + [meta(sp)], d, '채널끼리 엇갈린 블록, 면마다 1 부터 번호', dict(w=8, h=6, nc=2, nz=1, nt=1, chNames=['CH1', 'CH2']))
    sp3 = spec(); d3 = r.integers(0, 4096, (6, 8)).astype('<u2'); one = lambda u='_0_1_': pix_pair(f'{u}{sp3["guids"][0]}_0'.encode(), d3.tobytes())
    lsm = lsm_xml(sp3).replace(b'<fvCommonphase:pinholeDiameter>', b'<opticalelement:dye><opticalelement:name>Alexa Fluor 488</opticalelement:name></opticalelement:dye><fvCommonphase:pinholeDiameter>')
    run('p3_nestedname.oir', [frame(sp3)] + one() + [meta(sp3, lsm=lsm)], d3, '채널 안에 염료 name 이 또 있음 → 바로 아래 name (Bio-Formats 와 같음; oirfile 은 안쪽 마지막 name)', dict(w=8, h=6, nc=1, nz=1, nt=1, chNames=['CH1']))
    # p4: ZSTACK 축 정의 2 개 (step 9.0 · 0.5), length z 없음, 면 위치 하나 → 마지막 정의의 step (oirfile 과 같음)
    sp4 = spec(nz=3, Zstep=0.5); d4 = r.integers(0, 4096, (3, 6, 8)).astype('<u2')
    lsm = lsm_xml(sp4).replace(b'<commonparam:z>0.5</commonparam:z>', b'').replace(b'<commonimage:axis enable="true">', b'<commonimage:axis enable="true"><commonparam:axis>ZSTACK</commonparam:axis><commonparam:maxSize>3</commonparam:maxSize><commonparam:step>9.0</commonparam:step></commonimage:axis><commonimage:axis enable="true">', 1)
    bl = [frame(sp4)]
    for z in range(3): bl += pix_pair(f'z{z + 1:03d}_0_1_{sp4["guids"][0]}_0'.encode(), d4[z].tobytes())
    run('p4_twoaxes.oir', bl + [meta(sp4, lsm=lsm)], d4, 'ZSTACK 축 정의 2 개 → 마지막 것 (oirfile 과 같음; Bio-Formats 는 오류)', dict(w=8, h=6, nc=1, nz=3, nt=1, zStep=0.5, chNames=['CH1']))
    def nomarker(blocks):
        b = bytearray(assemble(blocks)); idx = struct.unpack('<Q', b[40:48])[0]; b[idx:idx + 4] = b'\0\0\0\0'; return bytes(b)
    run('p5_nomarker.oir', [frame(sp3)] + one() + [meta(sp3)], d3, 'block index 에 0xFFFFFFFF 표지 없음 (위치는 맞음)', dict(w=8, h=6, nc=1, nz=1, nt=1, chNames=['CH1']), asm=nomarker)
    run('p6_dupuid.oir', [frame(sp3)] + one() + pix_pair(f'_0_1_{sp3["guids"][0]}_0'.encode(), (d3 + 1).tobytes()) + [meta(sp3)], d3, '같은 UID 가 두 번 → 첫 것 + 알림 (Bio-Formats 는 마지막 것, oirfile 은 둘을 이어 높이 12)', dict(w=8, h=6, nc=1, nz=1, nt=1, warn='되풀이', chNames=['CH1']))
    run('p7_snapuid.oir', [frame(sp3)] + one('_') + [meta(sp3)], d3, 'UID 에 _a_b 없음 (스냅샷형)', dict(w=8, h=6, nc=1, nz=1, nt=1, chNames=['CH1']))
    run('p8_fileinfo_only.oir', [frame(sp3)] + one() + [(xml_block(0, fileinfo_xml()), b'')], d3, 'METADATA 에 fileinfo 만 (lsmimage 없음) → 화소만 + 알림', dict(w=8, h=6, nc=1, nz=1, nt=1, px=None, warn='메타데이터(lsmimage)가 없습니다', nometa=True))
    sp9 = spec(depth=4, bits=32); d9 = r.integers(0, 4096, (6, 8)).astype('<u4')
    run('p9_u32.oir', [frame(sp9)] + pix_pair(f'_0_1_{sp9["guids"][0]}_0'.encode(), d9.tobytes()) + [meta(sp9)], d9, '깊이 4 바이트, uint32 자료 (float 로 보면 비정규수) → uint32 (Bio-Formats 와 같음)', dict(w=8, h=6, nc=1, nz=1, nt=1, bit=32, warn='uint32', chNames=['CH1']))
    sp10 = spec(nc=2, nt=3, Tstep=2.0); d10 = r.integers(0, 4096, (3, 2, 6, 8)).astype('<u2'); bl = []; later = []
    for t in range(3):
        for c in range(2):
            f = frame(sp10, [('TIMELAPSE', 1000 * (5.0 + 2.5 * t))]); (bl if (t, c) == (0, 0) else later).append(f)
            bl += pix_pair(f't{t + 1:03d}_0_1_{sp10["guids"][c]}_0'.encode(), d10[t, c].tobytes())
    run('p10_fp_perchan.oir', bl + later + [meta(sp10)], d10, 'FRAMEPROPERTIES 가 채널마다 → frame interval 2.5 s (면 위치)', dict(w=8, h=6, nc=2, nz=1, nt=3, fi='2.5 s', chNames=['CH1', 'CH2']))
    sp11 = spec(nz=2, Zstep=1.0); d11 = r.integers(0, 4096, (2, 6, 8)).astype('<u2')
    bl = [frame(sp11)] + pix_pair(f'z001_0_1_{sp11["guids"][0]}_0'.encode(), d11[0].tobytes()) + pix_pair(f'z002_0_1_{sp11["guids"][0]}_0'.encode(), d11[1].tobytes() + b'\0' * 16)
    run('p11_extraline.oir', bl + [meta(sp11)], d11, '한 면에 한 줄 남는 자료 → 높이 6 유지 (Bio-Formats 와 같음; oirfile 은 7)', dict(w=8, h=6, nc=1, nz=2, nt=1, zStep=1.0, warn='남는 바이트', chNames=['CH1']))
    spU = spec(guids=[guid(0).upper()])
    run('p12_upper.oir', [frame(spU)] + pix_pair(f'_0_1_{spU["guids"][0]}_0'.encode(), d3.tobytes()) + [meta(spU)], d3, 'GUID 가 대문자 (oirfile 은 못 읽음, Bio-Formats 는 읽음)', dict(w=8, h=6, nc=1, nz=1, nt=1, chNames=['CH1']))
    spL = spec(nc=2, names=['GFP', 'RFP'], order=[2, 1]); dL = r.integers(0, 4096, (2, 6, 8)).astype('<u2'); lsm = lsm_xml(dict(spL, guids=[x.upper() for x in spL['guids']]))
    run('p13_case_mismatch.oir', [frame(spL)] + sum((pix_pair(f'_0_1_{spL["guids"][c]}_0'.encode(), dL[c].tobytes()) for c in range(2)), []) + [meta(spL, lsm=lsm)], dL[::-1],
        'lsmimage GUID 는 대문자, UID 는 소문자 → 대소문자 무시, order 순 RFP·GFP', dict(w=8, h=6, nc=2, nz=1, nt=1, chNames=['RFP', 'GFP']))
    spR = spec(nc=3, depth=1, bits=8, colortype='RGB', names=['B', 'G', 'R'], W=6, H=4); dR = r.integers(0, 256, (3, 4, 6)).astype('u1')
    run('p14_rgb_bgr.oir', [frame(spR)] + sum((pix_pair(f'_0_1_{spR["guids"][c]}_0'.encode(), dR[c].tobytes()) for c in range(3)), []) + [meta(spR)], dR[::-1],
        'RGB, 채널 차례 B·G·R (이름) → R·G·B 로', dict(w=6, h=4, nc=3, nz=1, nt=1, bit=8, rgb=True, warn='RGB 성분'))
    run('p15_refonly.oir', [frame(sp3)] + pix_pair(f'REF_CAMERA0_{sp3["guids"][0]}_0'.encode(), d3.tobytes()) + [meta(sp3)], d3,
        'REF_ 블록만 (본 영상 없음) → 참조 영상만 + 알림', dict(w=8, h=6, nc=1, nz=1, nt=1, warn='참조(개관) 영상만', chNames=['CH1'], refonly=True))
    sp16 = spec(nz=3, nt=2, Zstep=0.5); d16 = np.random.default_rng(2).integers(0, 4096, (2, 3, 6, 8)).astype('<u2'); bl = []; later = []
    for t in range(2):
        for z in range(3):
            f = (xml_block(1, frame_xml(sp16, [('ZSTACK', round(10 + 0.5 * z + 0.013 * t, 4)), ('TIMELAPSE', 1000.0 * t)])), I(0)); (bl if (t, z) == (0, 0) else later).append(f)
            bl += pix_pair(f'z{z + 1:03d}t{t + 1:03d}_0_1_{sp16["guids"][0]}_0'.encode(), d16[t, z].tobytes())
    run('p16_zdrift.oir', bl + later + [(xml_block(0, fileinfo_xml(), lsm_xml(sp16)), b'')], d16, 'T 사이 Z 표류 0.013 µm → 같은 T 안에서 잰 간격 0.5 (oirfile 0.013, Bio-Formats 0.5)',
        dict(w=8, h=6, nc=1, nz=3, nt=2, zStep=0.5, fi='1 s', chNames=['CH1']))
    # POIR 변형: data descriptor·한글 이름 / ZIP64 extra (작은 항목) / 묶음 전체 메모리 한도·취소
    src = open(os.path.join(D, 'gray12.oir'), 'rb').read()
    class Sink(io.RawIOBase):
        def __init__(s): s.b = bytearray()
        def writable(s): return True
        def write(s, x): s.b += x; return len(x)
    k = Sink()
    with zipfile.ZipFile(k, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        with z.open('Image_01/데이터 한글.oir', 'w') as f: f.write(src)
    open(os.path.join(D, 'q1_datadesc.poir'), 'wb').write(bytes(k.b))
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        with z.open('Image_01/a.oir', 'w', force_zip64=True) as f: f.write(src)
    open(os.path.join(D, 'q2_zip64entry.poir'), 'wb').write(bio.getvalue())
    bio = io.BytesIO(); big = open(os.path.join(D, 'c2z4big.oir'), 'rb').read()
    with zipfile.ZipFile(bio, 'w') as z:
        for i in range(3): z.writestr(f'Image_0{i}/c2z4big.oir', big)
    open(os.path.join(D, 'q3_three.poir'), 'wb').write(bio.getvalue())
    for q, inner, key in (('q1_datadesc.poir', 'Image_01/데이터 한글.oir', 'gray12.oir'), ('q2_zip64entry.poir', 'Image_01/a.oir', 'gray12.oir')):
        with oirfile.PoirFile(os.path.join(D, q)) as pf: assert list(pf) == [inner]
        out[q] = dict(desc='POIR ' + ('data descriptor·한글 이름' if q.startswith('q1') else 'ZIP64 extra 가 있는 작은 항목'), images=[key], inner=[inner], truth='oirfile PoirFile = 단일 파일 정답')
    out['q3_three.poir'] = dict(desc='POIR 에 64×48·C2·Z4·T3 OIR 3 개 — 메모리 한도는 묶음 전체', images=['c2z4big.oir'] * 3, inner=[f'Image_0{i}/c2z4big.oir' for i in range(3)], budget=True)


import zipfile, io
if __name__ == '__main__':
    only = sys.argv[3:] or None
    out = {}; table = []
    for fn, sp, desc, exp in CASES:
        if only and fn not in only: continue
        sp['base'] = fn[:-4]; b, data = build(sp); fp = os.path.join(D, fn); open(fp, 'wb').write(b)
        a, sizes, sc, names = oir_canon(fp)
        with oirfile.OirFile(fp) as o: chn = [c.name for c in o.channels]
        # 만든 자료를 oirfile 이 보여야 할 모양으로: 채널 = order 속성 순, 빠진 면 = 0 (가운데) 또는 없음 (끝)
        want = data.astype(np.float32)[:, :, :, np.argsort(sp['order'], kind='stable')]
        for (t, l, z, c, bb) in sp.get('drop', []): want[t, l, z, c] = 0
        if fn == 'miss_last.oir': want = want[:, :, :3]
        oir_ok = a.shape == want.shape and np.array_equal(a, want)
        S, err = bfdump(fp)
        T, L, Z, C, Y, X = a.shape
        app = to_app(a)  # (T·L, Z, C, Y, X)
        bfa = None if err else S['arr']
        bf_pix = bfa is not None and bfa.shape == app.shape and np.array_equal(bfa, app)
        bf_want = bfa is not None and bfa.shape == to_app(want).shape and np.array_equal(bfa, to_app(want))  # Bio-Formats = 만든 자료
        bf_px = None if err or S['px'] == '-' else float(S['px']); bf_pz = None if err or S['pz'] == '-' else float(S['pz'])
        o_px = sc.get('X'); o_pz = sc.get('Z')
        dims_ok = not err and (int(S['X']), int(S['Y']), int(S['Z']), int(S['C']) * int(S['rgb']), int(S['T'])) == (X, Y, Z, C, T * L)
        agree = bool(dims_ok and bf_pix and close(bf_px, o_px) and close(bf_pz, o_pz if Z > 1 else None) if not (Z > 1 and bf_pz is None and o_pz is None) else dims_ok and bf_pix and close(bf_px, o_px))
        # 정답: 두 구현이 같으면 그 값, 다르면 'follow' 쪽 (그 쪽이 만든 자료와 같을 때만)
        rec = dict(desc=desc, w=X, h=Y, nc=C, nz=Z, nt=T * L, bit={1: 8, 2: 16, 4: 32}[sp['depth']], px=o_px if o_px else bf_px, chNames=chn if chn else None)
        rec['zStep'] = exp.get('zStep', o_pz if Z > 1 else 0) or 0
        if 'rgb' in exp: rec['rgb'] = True
        if sp.get('wl'):  # 방출 파장 → 표시 색 (앱 규칙), 채널은 order 순
            lut = lambda w: 'Blue' if w < 470 else 'Cyan' if w < 505 else 'Green' if w < 560 else 'Yellow' if w < 590 else 'Red' if w < 700 else 'Magenta'
            rec['chLut'] = [lut(sum(sp['wl'][i]) / 2) for i in np.argsort(sp['order'], kind='stable')]
        for k in ('warn', 'choose'): 
            if k in exp: rec[k] = exp[k]
        dis = []
        if not agree:
            if err: dis.append(f'Bio-Formats 오류: {err}')
            else:
                if not dims_ok: dis.append(f'크기 BF X{S["X"]} Y{S["Y"]} Z{S["Z"]} C{S["C"]} T{S["T"]} / oirfile X{X} Y{Y} Z{Z} C{C} T{T} L{L}')
                if not bf_pix: dis.append('화소 다름')
                if not close(bf_px, o_px): dis.append(f'화소 크기 BF {bf_px} / oirfile {o_px}')
                if Z > 1 and not close(bf_pz, o_pz): dis.append(f'z step BF {bf_pz} / oirfile {o_pz}')
                if S['type'] != {1: 'uint8', 2: 'uint16', 4: 'float'}[sp['depth']]: dis.append(f'형식 BF {S["type"]}')
        rec['truth'] = 'oirfile = Bio-Formats' if agree else f"{exp.get('follow', '?')} (다름: {exp.get('why', '')})"
        if not agree and exp.get('pixels_only') and bf_pix: rec['truth'] = f"oirfile = Bio-Formats (화소) · {exp.get('why')}"
        if dis: rec['disagree'] = dis
        good = agree or (exp.get('follow') == 'oirfile' and oir_ok) or (exp.get('follow') == 'Bio-Formats' and bf_want)
        if not good: print('!! 정답 없음', fn, dis); continue
        truth = bfa if exp.get('follow') == 'Bio-Formats' else app
        rec['key'] = fn; rec['nz'] = truth.shape[1]; write_truth(fn, truth.reshape(-1))
        if exp.get('ref'):
            with oirfile.OirFile(fp) as o: r = o.reference.asarray().astype(np.float32)
            assert np.array_equal(r, sp['refdata'].astype(np.float32)); write_truth(fn + '#ref', r.reshape(-1)); rec['ref'] = dict(w=r.shape[-1], h=r.shape[-2], nc=r.shape[0] if r.ndim == 3 else 1)
        out[fn] = rec
        table.append((fn, f'{X}×{Y} C{C} Z{Z} T{T}' + (f' L{L}' if L > 1 else ''), 'ok' if oir_ok else 'X', 'ok' if agree else ('pix' if bf_pix else ('BF=자료' if bf_want else ('err' if err else 'X'))), rec['truth'][:70]))
    # 잘린 파일·손상 — 정답 = 만든 자료 (oirfile 은 index 가 없으면 못 읽음)
    if not only:
        sp = spec(nz=4, Zstep=0.5, W=8, H=6); sp['base'] = 'trunc'; b, data = build(sp)
        idx = struct.unpack('<Q', b[40:48])[0]
        open(os.path.join(D, 'noindex.oir'), 'wb').write(b[:idx + 2])  # block index 만 잘림 → 차례 읽기
        S, err = bfdump(os.path.join(D, 'noindex.oir'))
        try: oir_canon(os.path.join(D, 'noindex.oir')); oe = 'oirfile 읽음'
        except Exception as e: oe = f'oirfile 실패 ({type(e).__name__})'
        bf_ok = not err and np.array_equal(S['arr'], to_app(data.astype(np.float32)))
        write_truth('noindex.oir', to_app(data.astype(np.float32)).reshape(-1))
        out['noindex.oir'] = dict(desc='block index 가 잘림 (끝 일부만 남음)', w=8, h=6, nc=1, nz=4, nt=1, bit=16, px=0.31, zStep=0.5, chNames=['CH1'], key='noindex.oir', warn='index',
                                  truth=f'만든 자료 = Bio-Formats({"같음" if bf_ok else "다름"}); {oe}', disagree=[oe])
        table.append(('noindex.oir', '8×6 Z4', oe, 'ok' if bf_ok else 'X', '만든 자료'))
        k = b.find(b'z003_0_1'); cut = k + 47 + 8 + 40  # z3 면 PIXEL 안에서 잘림 → 메타데이터(파일 끝)도 없음 → 분명한 오류
        open(os.path.join(D, 'trunc_pix.oir'), 'wb').write(b[:cut])
        S, err = bfdump(os.path.join(D, 'trunc_pix.oir'))
        try: oir_canon(os.path.join(D, 'trunc_pix.oir')); oe = 'oirfile 읽음'
        except Exception as e: oe = f'oirfile 실패 ({type(e).__name__})'
        tp = data[0, 0, :3, 0].astype(np.float32).reshape(-1).copy(); tp[2 * 48 + (cut - (k + 47 + 8)) // 2:] = 0  # z3 은 앞 20 화소만
        write_truth('trunc_pix.oir', tp)
        out['trunc_pix.oir'] = dict(desc='z3 PIXEL 안에서 잘림 (index·메타데이터 없음 → 화소만)', w=8, h=6, nc=1, nz=3, nt=1, bit=16, px=None, zStep=0, chNames=['CH1'], key='trunc_pix.oir', warn='메타데이터(lsmimage)가 없습니다', nometa=True,
                                    truth=f'만든 자료 (잘린 뒤 = 0); {oe}; Bio-Formats {"오류: " + err if err else "읽음"}', disagree=[oe, f'Bio-Formats: {err or "읽음"}'])
        table.append(('trunc_pix.oir', '8×6 Z3(/4)', oe, 'err' if err else '?', '만든 자료 (화소만)'))
        open(os.path.join(D, 'header_only.oir'), 'wb').write(b[:40])
        out['header_only.oir'] = dict(desc='머리 40 바이트만', error='OIR 머리')
        sp = spec(nc=2, nt=4, Tstep=3.0); sp['base'] = 'split'; main, comp, data = build(sp, split_t=2)
        open(os.path.join(D, 'split.oir'), 'wb').write(main); open(os.path.join(D, 'split_00001'), 'wb').write(comp)
        with oirfile.OirFile(os.path.join(D, 'split.oir'), multifile=False) as o: a1 = o.asarray()
        with oirfile.OirFile(os.path.join(D, 'split.oir'), multifile=True) as o: a2 = o.asarray(); ncomp = len(o.companion_files)
        assert a1.shape[0] == 2 and a2.shape[0] == 4 and ncomp == 1
        td = tempfile.mkdtemp(); shutil.copy(os.path.join(D, 'split.oir'), td)  # 보조 파일 없이 (본 파일만 끌어 놓은 경우)
        S, err = bfdump(os.path.join(td, 'split.oir')); bf_split = not err and np.array_equal(S['arr'], a1.astype(np.float32)[:, None])
        S2, err2 = bfdump(os.path.join(D, 'split.oir'))  # 보조 파일이 옆에 있을 때
        write_truth('split.oir', a1.transpose(0, 1, 2, 3).reshape(-1).astype(np.float32))
        out['split.oir'] = dict(desc='큰 획득이 본 파일 + 보조 파일 split_00001 로 나뉨 (본 파일만 열기)', w=8, h=6, nc=2, nz=1, nt=2, bit=16, px=0.31, zStep=0, chNames=['CH1', 'CH2'], key='split.oir',
                                warn='보조 파일', truth=f'본 파일만: oirfile(multifile=False) T 2 {"= Bio-Formats" if bf_split else "≠ Bio-Formats"}; 보조 파일과 함께: oirfile(multifile=True) T 4, Bio-Formats {"오류 (" + err2 + ")" if err2 else "T " + S2["T"]}',
                                disagree=[f'보조 파일이 옆에 있을 때 Bio-Formats: {err2 or "T " + S2["T"]} / oirfile T 4 (앱은 보조 파일을 읽지 않음)'])
        table.append(('split.oir', '8×6 C2 T2(/4)', 'ok', 'ok' if bf_split else 'X', '본 파일만 = oirfile multifile=False = BF'))
        out['split_00001'] = dict(desc='보조 파일만 (메타데이터 없음)', error='보조 파일')
        # POIR: OIR 두 개를 ZIP (하나는 stored, 하나는 deflate)
        bio = io.BytesIO()
        with zipfile.ZipFile(bio, 'w') as z:
            z.write(os.path.join(D, 'gray12.oir'), 'Image_01/gray12.oir', compress_type=zipfile.ZIP_STORED)
            z.write(os.path.join(D, 'c3.oir'), 'Image_02/c3.oir', compress_type=zipfile.ZIP_DEFLATED)
            z.writestr('Image_01/readme.txt', 'not an oir')
        open(os.path.join(D, 'pair.poir'), 'wb').write(bio.getvalue())
        with oirfile.PoirFile(os.path.join(D, 'pair.poir')) as pf: names = list(pf); arrs = [pf[n].asarray() for n in names]
        assert np.array_equal(arrs[0], np.fromfile(os.path.join(D, 'gray12.oir.f32'), '<f4').reshape(arrs[0].shape))
        out['pair.poir'] = dict(desc='POIR (ZIP): OIR 2 개 (stored·deflate) + txt', images=['gray12.oir', 'c3.oir'], inner=names, truth='oirfile PoirFile 의 각 OIR = 단일 파일 정답')
        table.append(('pair.poir', '2 × OIR', 'ok', '-', 'PoirFile'))
    if not only: review_probes(out, table)
    json.dump(out, open(os.path.join(D, 'synth.json'), 'w'), ensure_ascii=False, indent=1)
    print(f'{"파일":16s} {"크기":22s} oirfile=만든자료  oirfile=BF  정답')
    for r in table: print(f'{r[0]:16s} {r[1]:22s} {r[2]:16s} {r[3]:10s} {r[4]}')
    print('합성 OIR', len(out))
