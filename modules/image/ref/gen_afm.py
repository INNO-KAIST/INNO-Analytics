"""AFM·SPM 형식(GWY·GSF·Nanoscope·ARDF) 시험 자료와 정답.

정답은 앱과 독립된 구현에서 나온다.
 - Gwyddion 2.64 (`gwyddion --convert-to-gwy`) 가 원본을 읽어 GWY 로 쓴 결과를 gwyfile(파이썬, 독립 GWY 파서)로 읽는다.
   Nanoscope·GSF·GWYO·ARDF 모두 Gwyddion 의 가져오기 모듈이 정답의 근거.
 - GWYP 파일은 gwyfile 로 직접 읽은 값이 정답이고, Gwyddion 변환 결과와도 같은지 확인한다.
 - 합성 Nanoscope 는 머리에 적은 원시 정수·눈금으로 손으로 값을 계산(원시/256^Bytes × Z scale × soft scale)하고, Gwyddion 결과와 같아야 한다.
 - 합성 ARDF 는 magic-afm(R. Sheridan) 의 ARDF 파서가 읽은 영상이 정답(첫 줄이 아래 → 앱은 위아래를 뒤집어 보여 줌).
   Gwyddion 결과(뒤집지 않음, 정사각만 맞음)와도 같은 값인지 확인한다.
사용: python3 -I gen_afm.py <합성 출력 폴더> [<실측 폴더>]
"""
import sys, os, json, struct, zlib, subprocess, tempfile, shutil, glob
import numpy as np
import gwyfile
from gwyfile.objects import GwyContainer, GwyDataField, GwySIUnit, GwyObject
import inspect, gwyfile.objects as _go
# Gwyddion 이 meta 문자열을 UTF-8 이 아닌 채로 쓰는 경우가 있다(Nanoscope 머리의 latin1) — gwyfile 이 멈추지 않게 문자열만 'replace' 로 푼다
exec(compile(inspect.getsource(_go.component_from_buffer).replace(".decode('utf-8')", ".decode('utf-8', 'replace')"), _go.__file__, 'exec'), _go.__dict__)

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
REAL = sys.argv[2] if len(sys.argv) > 2 else None
GWYHOME = tempfile.mkdtemp(prefix='gwyhome-')

# Gwyddion 의 가져오기 모듈을 그대로 부르는 작은 변환기 (ref/gwyconv/gwyconv.c). `gwyddion --convert-to-gwy` 는 GTK 를 초기화하지 않아
# 진행 창을 여는 ARDF 모듈이 실패하므로 쓰지 않는다. 빌드: libgwyddion20-dev, libgtk2.0-dev, libfftw3-dev
CONV = os.path.expanduser('~/.cache/gwyconv/gwyconv')
if not os.path.exists(CONV):
    os.makedirs(os.path.dirname(CONV), exist_ok=True)
    flags = subprocess.check_output(['pkg-config', '--cflags', '--libs', 'gwyddion'], text=True).split()
    subprocess.check_call(['gcc', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gwyconv', 'gwyconv.c'), '-o', CONV, '-Wno-deprecated-declarations'] + flags)
def gwy_convert(path):
    """Gwyddion 으로 읽어 GWY 로 쓰고, gwyfile 로 연 GwyContainer 를 돌려준다 (실패하면 None)."""
    out = tempfile.mktemp(suffix='.gwy')
    env = dict(os.environ, HOME=GWYHOME)
    r = subprocess.run(['xvfb-run', '-a', CONV, path, out], env=env, capture_output=True, timeout=300)
    if not os.path.exists(out):
        print('   Gwyddion:', r.stderr.decode(errors='replace').strip().splitlines()[-1:] if r.stderr else '')
        return None
    c = gwyfile.load(out); os.remove(out); return c

def unit(o):
    return o['unitstr'] if o is not None and 'unitstr' in o else ''

def fields(c):
    """Gwyddion 컨테이너 → [(종류, 제목, dict)] (데이터 채널·mask·brick, id 순)"""
    out = []
    ids = sorted(int(k.split('/')[1]) for k in c if k.count('/') == 2 and k.endswith('/data') and k.split('/')[1].isdigit())
    for i in ids:
        d = c[f'/{i}/data']; t = c.get(f'/{i}/data/title', f'Channel {i}')
        out.append(('data', t, dict(w=d['xres'], h=d['yres'], xreal=d.get('xreal', 1.0), yreal=d.get('yreal', 1.0),
                                   uxy=unit(d.get('si_unit_xy')), uz=unit(d.get('si_unit_z')), a=np.asarray(d['data'], np.float64).reshape(d['yres'], d['xres']))))
        if f'/{i}/mask' in c:
            m = c[f'/{i}/mask']
            out.append(('mask', t + ' · mask', dict(w=m['xres'], h=m['yres'], xreal=m.get('xreal', 1.0), yreal=m.get('yreal', 1.0), uxy=unit(m.get('si_unit_xy')), uz='',
                                                  a=(np.asarray(m['data'], np.float64).reshape(m['yres'], m['xres']) > 0.5) * 255.0)))
    bids = sorted(int(k.split('/')[2]) for k in c if k.startswith('/brick/') and k.count('/') == 2)
    for i in bids:
        b = c[f'/brick/{i}']; t = c.get(f'/brick/{i}/title', f'Volume {i}')
        a = np.asarray(b['data'], np.float64).reshape(b['zres'], b['yres'], b['xres'])
        out.append(('brick', t, dict(w=b['xres'], h=b['yres'], nz=b['zres'], xreal=b.get('xreal', 1.0), yreal=b.get('yreal', 1.0), zreal=b.get('zreal', 1.0),
                                    uxy=unit(b.get('si_unit_x')), uz=unit(b.get('si_unit_w')), uzz=unit(b.get('si_unit_z')), a=a)))
    return out

TRUTH = {}
def save(folder, fname, recs, extra=None, flip=False):
    """recs: fields() 형식. flip=True 면 정답을 위아래 뒤집어 저장(앱 표시 방향)."""
    tdir = os.path.join(folder, 'truth'); os.makedirs(tdir, exist_ok=True)
    L = []
    for j, (kind, t, d) in enumerate(recs):
        a = d['a'][..., ::-1, :] if flip else d['a']
        fn = f'{fname}#{j}.f64'; np.ascontiguousarray(a, np.float64).tofile(os.path.join(tdir, fn))
        r = {k: v for k, v in d.items() if k != 'a'}; r.update(kind=kind, title=t, file=fn); L.append(r)
    TRUTH.setdefault(folder, {})[fname] = dict(images=L, **(extra or {}))

def same(a, b, tol=1e-12):
    a = np.asarray(a, np.float64); b = np.asarray(b, np.float64)
    return a.shape == b.shape and np.allclose(a, b, rtol=tol, atol=tol * max(1e-300, np.nanmax(np.abs(b)) if b.size else 0), equal_nan=True)

# =============================================================== 합성 GWY (gwyfile 로 씀)
def df(a, xreal, yreal, uxy, uz, xoff=0.0):
    d = GwyDataField(np.ascontiguousarray(a, np.float64), xreal=xreal, yreal=yreal, xoff=xoff, si_unit_xy=GwySIUnit(unitstr=uxy), si_unit_z=GwySIUnit(unitstr=uz))
    return d
rng = np.random.default_rng(7)
yy, xx = np.mgrid[0:24, 0:32]
h = (20e-9 * np.sin(xx / 5.0) * np.cos(yy / 4.0) + 3e-9 * rng.standard_normal((24, 32)))
c = GwyContainer()
c['/0/data'] = df(h, 3e-6, 2.25e-6, 'm', 'm', xoff=1e-6); c['/0/data/title'] = 'Height 높이'
c['/0/mask'] = df((h > 5e-9).astype(np.float64), 3e-6, 2.25e-6, 'm', '')
meta = GwyContainer(); meta['Scan rate'] = '1.5 Hz'; meta['Comment'] = 'synthetic'; c['/0/meta'] = meta
c['/1/data'] = df(0.5 * np.cos(xx / 3.0) + 0.1 * yy / 24, 3e-6, 2.25e-6, 'm', 'V'); c['/1/data/title'] = 'Phase signal'
c['/2/data'] = df(1e6 * (1 + xx + 10 * yy), 3e-6, 2.25e-6, 'm', 'Pa'); c['/2/data/title'] = 'Modulus'
c['/5/data'] = df(2.0 + 0.0 * xx, 40e-9, 30e-9, 'm', 'N/m'); c['/5/data/title'] = 'Stiffness'
# GwyBrick: data[x + xres*(y + yres*z)]
bz = np.stack([1e-9 * (z + 1) * (xx[:8, :10] + 2 * yy[:8, :10]) for z in range(5)])
brick = GwyObject('GwyBrick', dict(xres=10, yres=8, zres=5, xreal=1e-6, yreal=0.8e-6, zreal=2.5e-6, data=bz.ravel(),
                                   si_unit_x=GwySIUnit(unitstr='m'), si_unit_y=GwySIUnit(unitstr='m'), si_unit_z=GwySIUnit(unitstr='m'), si_unit_w=GwySIUnit(unitstr='m')),
                  typecodes=dict(xres='i', yres='i', zres='i', xreal='d', yreal='d', zreal='d', data='D'))
c['/brick/0'] = brick; c['/brick/0/title'] = 'z-stack volume'
cm = GwyObject('GwyGraphCurveModel', dict(xdata=np.linspace(0, 1, 5), ydata=np.linspace(0, 2, 5), description='line'), typecodes=dict(xdata='D', ydata='D', description='s'))
g = GwyObject('GwyGraphModel', dict(curves=[cm], title='profile graph', top_label='', left_label='', right_label='', bottom_label='',
                                    x_unit=GwySIUnit(unitstr='m'), y_unit=GwySIUnit(unitstr='m')),
              typecodes=dict(curves='O', title='s', top_label='s', left_label='s', right_label='s', bottom_label='s'))
c['/0/graph/graph/1'] = g
p = os.path.join(OUT, 'g_multi.gwy'); c.tofile(p)
chk = gwyfile.load(p); assert same(chk['/0/data'].data, h)
recs = fields(chk)
gv = gwy_convert(p); assert gv is not None, 'Gwyddion 이 합성 GWY 를 못 읽음'
gr = fields(gv)
assert [r[0] for r in gr] == [r[0] for r in recs] and all(same(a[2]['a'], b[2]['a']) for a, b in zip(gr, recs)), 'Gwyddion 과 gwyfile 이 다름'
save(OUT, 'g_multi.gwy', recs, dict(warn=['그래프']))

# 합성 GWYO (Gwyddion 1.x): 맨 위 컨테이너만 옛 직렬화 (int32 GType, 이름, 값), meta 는 '/meta/…' 문자열 열쇠
def gwyo(path, entries):
    body = b''
    for k, v in entries:
        if isinstance(v, GwyObject): body += struct.pack('<i', 80) + k.encode() + b'\0' + v.serialize()
        elif isinstance(v, str): body += struct.pack('<i', 64) + k.encode() + b'\0' + v.encode('latin1') + b'\0'
        elif isinstance(v, bool): body += struct.pack('<i', 20) + k.encode() + b'\0' + bytes([int(v)])
        elif isinstance(v, int): body += struct.pack('<i', 24) + k.encode() + b'\0' + struct.pack('<i', v)
        else: body += struct.pack('<i', 60) + k.encode() + b'\0' + struct.pack('<d', v)
    open(path, 'wb').write(b'GWYO' + b'GwyContainer\0' + struct.pack('<I', len(body)) + body)
ho = 1e-8 * (xx[:6, :9] - yy[:6, :9])
gwyo(os.path.join(OUT, 'g_old.gwy'), [('/0/data', df(ho, 9e-7, 6e-7, 'm', 'm')), ('/meta/Operator', 'Synthetic \xb5m'), ('/0/show', True), ('/0/base/min', 0.5), ('/0/data/visible', 1)])
gv = gwy_convert(os.path.join(OUT, 'g_old.gwy')); gr = fields(gv); assert len(gr) == 1 and same(gr[0][2]['a'], ho), 'GWYO: Gwyddion 결과가 다름'
save(OUT, 'g_old.gwy', gr, dict(warn=['GWYO']))

# =============================================================== 합성 GSF (명세: 머리 + NUL 1–4 개로 4 바이트 경계, float32 LE)
def write_gsf(path, a, hdr):
    s = 'Gwyddion Simple Field 1.0\n' + ''.join(f'{k} = {v}\n' for k, v in hdr.items())
    b = s.encode('utf-8'); pad = 4 - len(b) % 4
    open(path, 'wb').write(b + b'\0' * pad + np.ascontiguousarray(a, '<f4').tobytes())
for nm, extra in [('s_a.gsf', {}), ('s_b.gsf', {'Comment': 'x'}), ('s_c.gsf', {'Comment': 'xy'}), ('s_d.gsf', {'Comment': 'xyz'})]:
    a = (1e-9 * (xx[:10, :13] * 3 - yy[:10, :13])).astype(np.float32).astype(np.float64)
    if nm == 's_b.gsf': a[2, 3] = np.nan
    hdr = dict(XRes=13, YRes=10, XReal=1.3e-6, YReal=1.0e-6, XYUnits='m', ZUnits='m', Title='Topo ' + nm, **extra)
    p = os.path.join(OUT, nm); write_gsf(p, a, hdr)
    gv = gwy_convert(p); gr = fields(gv); d0 = gr[0][2]; ok = ~np.isnan(a)   # Gwyddion 은 NaN 을 mask 로 바꾸고 값을 채운다
    assert same(d0['a'][ok], a[ok]), nm; assert abs(d0['xreal'] - 1.3e-6) < 1e-18 and d0['uz'] == 'm'
    save(OUT, nm, [('data', gr[0][1], dict(d0, a=a))], dict(nan=int((~ok).sum())))

# =============================================================== 합성 Nanoscope (값을 손으로 계산 + Gwyddion 대조)
def write_ns(path, version, chans, scan_size, samps, lines, glob_lines=None, aspect='1:1', glob_aspect='1:1', bpp_file=2, extra_scanner='', img_res=None, magic='\\*File list', mode='Image', scan_text=None, nanoscope_e=False):
    """chans: [dict(name, soft, zscale_v, hard_lsb, line_dir, raw(int array H×W), qbpp, image_data_soft)]"""
    HLEN = 16384
    body = []
    data = b''
    offs = []
    off = HLEN
    for ch in chans:
        raw = ch['raw']; bb = raw.astype('<i4' if bpp_file == 4 else '<i2').tobytes()
        offs.append((off, len(bb))); data += bb; off += len(bb)
    L = [magic, f'\\Version: 0x{version:08x}', '\\Date: 10:00:00 AM Thu Oct 16 2026', '\\Start context: OL', f'\\Data length: {HLEN}', '\\Text: ',
         '\\*Equipment list', '\\Description: synthetic',
         '\\*Scanner list', '\\Scanner type: Synthetic', '\\@Sens. Zsens: V 23.12000 nm/V', '\\@Sens. DeflSens: V 50.00000 nm/V', '\\@Sens. Phase: V 1.000000 ', '\\@Sens. StiffSens: V 1.000000 ']
    if extra_scanner: L.append(extra_scanner)
    L += ['\\*Ciao scan list', f'\\Samps/line: {samps}', f'\\Lines: {glob_lines or lines}', f'\\Aspect Ratio: {glob_aspect}', f'\\Scan Size: {scan_size[0]} {scan_size[2]}', f'\\Operating mode: {mode}']
    for ch, (o, n) in zip(chans, offs):
        L += ['\\*Ciao image list', f'\\Data offset: {o}', f'\\Data length: {n}', f'\\Bytes/pixel: {ch.get("qbpp", 2)}', f'\\Start context: OL2',
              *([f'\\Samps/line: {(img_res or (samps, lines))[0]} {(img_res or (samps, lines))[1]}'] if nanoscope_e else [f'\\Samps/line: {(img_res or (samps, lines))[0]}', f'\\Number of lines: {(img_res or (samps, lines))[1]}']),
              f'\\Aspect Ratio: {aspect}', f'\\Line Direction: {ch["line_dir"]}',
              scan_text or f'\\Scan Size: {scan_size[0]} {scan_size[1]} {scan_size[2]}',
              f'\\@2:Image Data: S [{ch["soft"] or ""}] "{ch["name"]}"']
        if ch.get('zline'):
            L += [ch['zline']]
        elif ch.get('old_z'):
            L += [f'\\Z scale: {ch["old_z"]}']
        elif ch['soft']:
            L += [f'\\@2:Z scale: V [{ch["zsoft"]}] ({ch["hard_lsb"]} V/LSB) {ch["zscale_v"]} V']
        else:
            L += [f'\\@2:Z scale: V ({ch["hard_lsb"]} V/LSB) {ch["zscale_v"]} {ch.get("zunit", "V")}']
    L += ['\\*File list end']
    hb = ('\r\n'.join(L) + '\r\n').encode('latin1'); assert len(hb) < HLEN
    open(path, 'wb').write(hb + b'\x1a' + b'\0' * (HLEN - len(hb) - 1) + data)

def ns_case(nm, version, chans, scan, samps, lines, expect, **kw):
    """expect: [(제목 일부, 값 배율(원시 정수 → SI), dx(m), dy(m), SI 단위)] — 값은 원시 정수를 위아래 뒤집어 배율을 곱한 것"""
    p = os.path.join(OUT, nm); write_ns(p, version, chans, scan, samps, lines, **kw)
    gv = gwy_convert(p); assert gv is not None, nm + ': Gwyddion 실패'
    gr = fields(gv); assert len(gr) == len(chans), (nm, len(gr))
    recs = []
    for ch, (kind, t, d), (q, dx, dy, u) in zip(chans, gr, expect):
        mine = ch['raw'][::-1, :].astype(np.float64) * q
        assert same(d['a'], mine, 1e-9), (nm, t, np.nanmax(np.abs(d['a'] - mine)))
        assert abs(d['xreal'] / d['w'] - dx) < 1e-9 * dx and abs(d['yreal'] / d['h'] - dy) < 1e-9 * dy, (nm, d['xreal'] / d['w'], dx, d['yreal'] / d['h'], dy)
        assert d['uz'] == u, (nm, d['uz'], u)
        recs.append((kind, t, dict(d, a=mine)))
    save(OUT, nm, recs, dict(src='hand formula = Gwyddion'))

R = lambda H, W, lo, hi, seed: np.random.default_rng(seed).integers(lo, hi, size=(H, W))
# v5.31, 16-bit, ~m, Height(nm/V) + Phase(soft 단위 없음 → deg) + soft scale 없음(V)
H_, W_ = 24, 32
ns_case('ns_v531.spm', 0x05310001,
        [dict(name='Height', soft='Height', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='440.0000', line_dir='Trace', raw=R(H_, W_, -30000, 30000, 1)),
         dict(name='Phase', soft='Phase', zsoft='Sens. Phase', hard_lsb='0.0003750000', zscale_v='24.57600', line_dir='Retrace', raw=R(H_, W_, -20000, 20000, 2)),
         dict(name='Potential', soft=None, hard_lsb='0.0003750000', zscale_v='24.57600', zunit='V', line_dir='Trace', raw=R(H_, W_, -32768, 32767, 3))],
        (2, 1.5, '~m'), W_, H_,
        [(440.0 * 23.12e-9 / 65536, 2e-6 / W_, 1.5e-6 / H_, 'm'), (24.576 * 1.0 / 65536, 2e-6 / W_, 1.5e-6 / H_, 'deg'), (24.576 / 65536, 2e-6 / W_, 1.5e-6 / H_, 'V')])
# v9.2 이후: 자료가 32-bit, Bytes/pixel 4 → 원시/2^32
ns_case('ns_v92_32bit.spm', 0x09200201,
        [dict(name='Height Sensor', soft='ZSensor', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='440.0000', line_dir='Retrace', qbpp=4, raw=R(H_, W_, -2**31, 2**31 - 1, 4)),
         dict(name='Deflection Error', soft='DeflectionError', zsoft='Sens. DeflSens', hard_lsb='0.0003750000', zscale_v='24.57600', line_dir='Retrace', qbpp=4, raw=R(H_, W_, -2**30, 2**30, 5))],
        (500, 500, 'nm'), W_, H_,
        [(440.0 * 23.12e-9 / 2**32, 500e-9 / W_, 500e-9 / H_, 'm'), (24.576 * 50e-9 / 2**32, 500e-9 / W_, 500e-9 / H_, 'm')], bpp_file=4)
# 직사각 주사(Aspect Ratio 2:1, 전역 Lines = 영상 줄 수) — Gwyddion: yreal × yres / gyres
ns_case('ns_aspect.spm', 0x08150203,
        [dict(name='Height', soft='Height', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='220.0000', line_dir='Retrace', raw=R(16, 32, -30000, 30000, 6))],
        (4, 4, '~m'), 32, 16,   # 실제 파일처럼 Scan Size 는 두 수가 같고, 직사각은 Aspect Ratio·줄 수로 (Gwyddion: yreal × yres/xres)
        [(220.0 * 23.12e-9 / 65536, 4e-6 / 32, 2e-6 / 16, 'm')], glob_lines=16, aspect='2:1', glob_aspect='2:1')

# 영상 절의 크기가 틀리고 전역(Ciao scan list) 크기만 맞는 파일 (Gwyddion 주석: 0x05300001 에서 관찰) → 전역 크기, 가로 길이도 비율대로
ns_case('ns_globalres.spm', 0x05300001,
        [dict(name='Height', soft='Height', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='440.0000', line_dir='Trace', raw=R(48, 64, -30000, 30000, 8))],
        (3, 3, '~m'), 64, 48,
        [(440.0 * 23.12e-9 / 65536, 3e-6 * (64 / 32) / 64, 3e-6 * (48 / 32) / 48, 'm')], img_res=(32, 32))

# 독립 검토(2026-10-16)가 짚은 경우들 — 모두 Gwyddion 과 손 계산이 같아야 함
Hc, Wc = 16, 24; rr = lambda sd: R(Hc, Wc, -30000, 30000, sd)
# \*Force file list 로 시작해도 Operating mode 가 Force Volume 이면 영상 채널을 읽는다 (Gwyddion: 첫 머리는 처음 종류일 뿐)
ns_case('ns_fv_magic.spm', 0x08150203, [dict(name='Height', soft='Height', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='440.0000', line_dir='Trace', raw=rr(11))],
        (2, 2, '~m'), Wc, Hc, [(440.0 * 23.12e-9 / 65536, 2e-6 / Wc, 2e-6 / Hc, 'm')], magic='\\*Force file list', mode='Force Volume')
# Nanoscope E: "Samps/line: 24 16" (Number of lines 없음), Scan Size 두 수 사이 빈칸 둘, @4:Z scale, 옛 형식 Z scale, Phase soft scale 단위 º/V, soft scale 이 없는 경우(V)
ns_case('ns_misc.spm', 0x05300001,
        [dict(name='Height4', soft='Height', zline='\\@4:Z scale: V [Sens. Zsens] (0.006713867 V/LSB) 440.0000 V', line_dir='Trace', raw=rr(12)),
         dict(name='OldZ', soft='OldZ', old_z='0.5 nm', line_dir='Trace', raw=rr(13)),
         dict(name='PhaseDeg', soft='PhaseX', zline='\\@2:Z scale: V [Sens. PhaseX] (0.0003750000 V/LSB) 24.57600 V', line_dir='Retrace', raw=rr(14)),
         dict(name='Missing', soft='Nope', zline='\\@2:Z scale: V [Sens. Nope] (0.0003750000 V/LSB) 24.57600 V', line_dir='Retrace', raw=rr(15))],
        (3, 3, '~m'), Wc, Hc,
        [(440.0 * 23.12e-9 / 65536, 3e-6 / Wc, 3e-6 / Hc, 'm'), (0.5e-9 / 65536, 3e-6 / Wc, 3e-6 / Hc, 'm'), (24.576 * 18.0 / 65536, 3e-6 / Wc, 3e-6 / Hc, 'deg'), (24.576 / 65536, 3e-6 / Wc, 3e-6 / Hc, '')],
        nanoscope_e=True, scan_text='\\Scan Size: 3  3 ~m', extra_scanner='\\@Sens. PhaseX: V 18.00000 \xba/V')
# 영상만 직사각(Aspect Ratio 2:1), 전역은 1:1 → yreal × yres / gyres (Gwyddion 의 흔한 경우)
ns_case('ns_aspect_local.spm', 0x08150203, [dict(name='Height', soft='Height', zsoft='Sens. Zsens', hard_lsb='0.006713867', zscale_v='440.0000', line_dir='Retrace', raw=R(16, 32, -30000, 30000, 16))],
        (4, 4, '~m'), 32, 16, [(440.0 * 23.12e-9 / 65536, 4e-6 / 32, 4e-6 * 16 / 32 / 16, 'm')], glob_lines=32, aspect='2:1', glob_aspect='1:1')

# =============================================================== 합성 ARDF (magic-afm 정답 + Gwyddion 대조)
def blk(name, payload, size=None, flags=0):
    size = size or 16 + len(payload)
    b = bytearray(struct.pack('<II4sI', 0, size, name.encode('latin1'), flags) + payload); assert len(b) <= size, (name, len(b), size)
    b += b'\0' * (size - len(b)); struct.pack_into('<I', b, 0, zlib.crc32(bytes(b[4:size]))); return bytes(b)
def cs(s, n=32): b = s.encode('latin1'); return b + b'\0' * (n - len(b))
def toc(name, entries, stride, pad=0):
    """entries: [(이름, payload bytes(stride-16))], pad: 끝의 0 채움 항목 수 (실제 파일처럼)"""
    body = b''.join(blk(n, p, stride) for n, p in entries) + b'\0' * (stride * pad)
    n = len(entries) + pad
    return blk(name, struct.pack('<QII', 32 + n * stride, n, stride), 32) + body
def text_blk(idx, s):
    t = s.replace('\n', '\r').encode('latin1'); size = 24 + len(t) + 1; size += (-size) % 4
    return blk('TEXT', struct.pack('<II', idx, len(t)) + t, size)
def ardf(path, images, notes, volm=None, zero_line=None, pad_ftoc=0, two_notes=False):
    """images: [(제목, 단위, 배열 lines×points(디스크 순서), x_step, y_step)]"""
    parts = []   # (bytes) — 위치는 순서대로 계산
    pos = 0
    def place(b):
        nonlocal pos; o = pos; parts.append(b); pos += len(b); return o
    n_items = len(images) + (1 if volm else 0)
    place(blk('ARDF', b'', 16))
    ftoc_at = place(b'')   # 자리만 — 나중에 채움
    ftoc_len = 32 + (n_items + pad_ftoc) * 24; parts[-1] = b'\0' * ftoc_len; pos += ftoc_len
    ttoc_at = pos; place(b'\0' * (32 + 32))
    note_at = place(text_blk(0, notes))
    item_ptrs = []
    for k, (title, un, a, dx, dy) in enumerate(images):
        lines, points = a.shape
        imag_at = pos
        item_ptrs.append(('IMAG', imag_at))
        tsz = 32 + 1 * 24
        place(b'\0' * tsz)                                  # IMAG 목차 (NEXT 1 항목) — 나중에
        ttoc2 = pos; place(b'\0' * (32 + 32 * (2 if two_notes else 1)))   # 영상 메모 TTOC
        idef = blk('IDEF', struct.pack('<IIQQdd', points, lines, 0, 0, dx, dy) + cs('m') + cs('m') + cs(title) + cs(un), 184)
        place(idef)
        stride = 16 + 4 * points
        idats = b''.join((b'\0' * stride) if zero_line == (k, r) else blk('IDAT', np.ascontiguousarray(a[r], '<f4').tobytes(), stride) for r in range(lines))
        ibox = blk('IBOX', struct.pack('<QII', 32 + len(idats), lines, stride), 32)
        place(ibox + idats); place(blk('GAMI', b'', 16))
        inote = place(text_blk(0, f'ImageNote: {title}\nChannel: {k}\n'))   # IDEF 는 TTOC 바로 뒤여야 하므로 메모는 영상 뒤에
        inote2 = place(text_blk(1, f'SecondNote: {title}\n')) if two_notes else None
        item_ptrs[-1] = ('IMAG', imag_at, ttoc2, inote, inote2)
    if volm:
        P_, L_, npts = volm
        vol_at = pos; item_ptrs.append(('VOLM', vol_at))
        place(b'\0' * (32 + 24))                            # VOLM 목차 (NSET)
        vtt = pos; place(b'\0' * (32 + 32))
        place(blk('VDEF', struct.pack('<II', P_, L_) + b'\0' * 24 + struct.pack('<ddd', 1e-7, 1e-7, 1e-4) + cs('m') + cs('m') + cs('s') + cs('Ext;Ret;') + struct.pack('<Q', 2), 208))
        for nm in ['Raw', 'Defl', 'ZSnsr']: place(blk('VCHN', cs(nm) + cs('m'), 80))
        xd = 'Ext;Ret;'.encode('latin1'); place(blk('XDEF', struct.pack('<II', 0, len(xd)) + xd, 96))
        vtoc_at = pos; place(b'\0' * (32 + L_ * 40)); place(blk('MLOV', b'', 16))
        line_ptr = []; prev = 0
        for r in range(L_):
            line_ptr.append(pos)
            for q in range(P_):
                vs_at = pos
                vset_size = 48
                vnam = blk('VNAM', struct.pack('<IIII', r * P_ + q, r, q, 8) + b'Curve\0\0\0', 40)
                vdats = b''
                n_ = npts + 2 * ((r * P_ + q) % 3)   # 곡선 길이가 점마다 다름 (force map 처럼 → 목차 간격이 고르지 않음)
                for ci in range(3):
                    f = np.linspace(0, 1, n_, dtype='<f4') * (ci + 1) * 1e-7
                    vdats += blk('VDAT', struct.pack('<10I', r * P_ + q, r, q, n_, ci, 0, n_ // 2, n_, 0, 0) + f.tobytes(), 56 + 4 * n_)
                xdat = blk('XDAT', b'', 16)
                nxt = vs_at + vset_size + len(vnam) + len(vdats) + len(xdat)
                place(blk('VSET', struct.pack('<IIIIQQ', r * P_ + q, r, q, 2, prev, nxt), vset_size) + vnam + vdats + xdat)
                prev = vs_at
        vnote = place(text_blk(0, 'VolumeNote: fmap\n'))
        vt = toc('VTOC', [('VOFF', struct.pack('<IIQQ', r * P_, r, 0, line_ptr[r])) for r in range(L_)], 40)
        volm_fill = (vol_at, vtt, vnote, vtoc_at, vt)
    buf = bytearray(b''.join(parts))
    def put(at, b): buf[at:at + len(b)] = b
    put(16, toc('FTOC', [(e[0], struct.pack('<Q', e[1])) for e in item_ptrs], 24, pad_ftoc))
    put(ttoc_at, toc('TTOC', [('TOFF', struct.pack('<QQ', 0, note_at))], 32))
    imags = [e for e in item_ptrs if e[0] == 'IMAG']
    for k, e in enumerate(imags):
        nxt = imags[k + 1][1] if k + 1 < len(imags) else (item_ptrs[-1][1] if volm else 0)
        put(e[1], toc('IMAG', [('NEXT', struct.pack('<Q', nxt))], 24))
        put(e[2], toc('TTOC', [('TOFF', struct.pack('<QQ', 0, e[3]))] + ([('TOFF', struct.pack('<QQ', 1, e[4]))] if e[4] else []), 32))
    if volm:
        vol_at, vtt, vnote, vtoc_at, vt = volm_fill
        put(vol_at, toc('VOLM', [('NSET', struct.pack('<Q', P_ * L_))], 24))
        put(vtt, toc('TTOC', [('TOFF', struct.pack('<QQ', 0, vnote))], 32))
        put(vtoc_at, vt)
    open(path, 'wb').write(bytes(buf))

import magic_afm.data_readers as MR
def ardf_case(nm, images, volm, square, **kw):
    notes = 'ImagingMode: AC Mode\nFastScanSize: {}\nSlowScanSize: {}\nSpringConstant: 2.1\nInvOLS: 5e-08\nScanRate: 0.5\n'.format(
        images[0][3] * images[0][2].shape[1], images[0][4] * images[0][2].shape[0])
    p = os.path.join(OUT, nm); ardf(p, images, notes, volm, **kw)
    zl = kw.get('zero_line')
    if zl: images[zl[0]][2][zl[1], :] = 0     # 빈 줄: 디스크 값은 0 (magic-afm·Gwyddion 은 0 으로 읽음, 앱은 NaN 으로 표시)
    data = open(p, 'rb').read()
    f = MR.ARDFFile.parse(data)     # magic-afm 이 파일 전체를 받아들여야 한다 (CRC·구조)
    recs = []
    for (title, un, a, dx, dy) in images:
        im = f.images[title]; arr = im.get_image()
        assert arr.shape == a.shape and np.array_equal(arr, a.astype(np.float32)), nm
        (sx, ux), (sy, uy) = im.step_info
        assert ux == 'm' and uy == 'm' and im.units == un
        recs.append(('data', title, dict(w=a.shape[1], h=a.shape[0], xreal=sx * a.shape[1], yreal=sy * a.shape[0], uxy='m', uz=un, a=arr.astype(np.float64))))
    gv = gwy_convert(p); assert gv is not None, nm + ': Gwyddion 실패'
    gr = [r for r in fields(gv) if r[0] == 'data']; assert len(gr) == len(images), (nm, len(gr))
    if square:
        for (k, t, d), (_, _, mine) in zip(gr, recs):
            assert same(d['a'], mine['a']), (nm, t)   # Gwyddion: 디스크 순서 그대로
    if zl: recs[zl[0]][2]['a'][zl[1], :] = np.nan
    save(OUT, nm, recs, dict(src='magic-afm (첫 줄이 아래 → 뒤집음)' + (' = Gwyddion' if square else ''), warn=(['VOLM'] if volm else [])), flip=True)

AY, AX = np.mgrid[0:16, 0:16]
ardf_case('a_fmap.ardf', [('Height', 'm', (5e-9 * np.sin(AX / 3.0) + 1e-9 * AY).astype(np.float32), 1e-7, 1e-7),
                          ('Adhesion', 'N', (1e-9 * (1 + AX + AY * 0.5)).astype(np.float32), 1e-7, 1e-7),
                          ('YoungsModulus', 'Pa', (2e6 * (1 + AY) + 1e5 * AX).astype(np.float32), 1e-7, 1e-7)], (16, 16, 12), True)
BY, BX = np.mgrid[0:12, 0:20]
ardf_case('a_rect.ardf', [('Height', 'm', (1e-8 * BX - 2e-9 * BY).astype(np.float32), 5e-8, 1e-7)], (20, 12, 8), False)
# 목차 끝 0 채움 항목, 영상 메모 TOFF 둘, 비어 있는(0) IDAT 줄 하나 — 독립 검토가 짚은 경우
ardf_case('a_pad.ardf', [('Height', 'm', (3e-9 * AX + 1e-9 * AY).astype(np.float32), 1e-7, 1e-7), ('Phase', 'deg', (AX - AY).astype(np.float32) * 1.5, 1e-7, 1e-7)], (16, 16, 6), True, zero_line=(0, 5), pad_ftoc=2, two_notes=True)

# =============================================================== 실측
if REAL:
    for p in sorted(glob.glob(os.path.join(REAL, '*'))):
        fn = os.path.basename(p)
        if os.path.isdir(p): continue
        head = open(p, 'rb').read(32)
        if head.startswith(b'\\*Force file list'):
            TRUTH.setdefault(REAL, {})[fn] = dict(images=[], expect_error='힘 곡선'); print(fn, '→ 힘 곡선 (오류 기대)'); continue
        gv = gwy_convert(p)
        if gv is None:
            print(fn, 'Gwyddion 도 못 읽음 — 건너뜀'); continue
        gr = fields(gv)
        if head.startswith(b'GWYP'):
            direct = fields(gwyfile.load(p))
            assert len(direct) == len(gr) and all(same(a[2]['a'], b[2]['a']) for a, b in zip(direct, gr)), fn + ': gwyfile 과 Gwyddion 이 다름'
            gr = direct
        save(REAL, fn, gr, dict(src='Gwyddion 2.64' + (' = gwyfile' if head.startswith(b'GWYP') else '')))
        print(fn, '→', len(gr), '영상', [r[1] for r in gr][:8])

for folder, t in TRUTH.items():
    json.dump(t, open(os.path.join(folder, 'truth', 'truth.json'), 'w'), ensure_ascii=False, indent=0)
    print(folder, '정답', len(t), '파일')
shutil.rmtree(GWYHOME, ignore_errors=True)