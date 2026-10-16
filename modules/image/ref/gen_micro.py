"""현미경 형식(nd2·czi) 읽기 시험 파일과 정답을 만든다.
 - nd2: Nikon ND2 v3(chunk map + CLx Lite Variant) 를 형식대로 직접 쓰고, **독립 구현인 Python `nd2` 패키지**가 같은 배열·보정을 읽는지 먼저 확인한다.
        정답 = nd2.ND2File 이 돌려준 배열(우리가 쓴 배열이 아니라) → 앱 순서 c + nc·(z + nz·t) 로 다시 놓은 것.
 - czi: Zeiss 의 libCZI(pylibCZIrw) 로 쓰고, 정답 = 독립 구현 czifile(Gohlke) 이 읽은 배열.
사용: python3 -I gen_micro.py <출력 폴더>   (출력 폴더의 truth.json 에 항목을 더한다)"""
import sys, os, json, struct, zlib, warnings
import numpy as np

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
tpath = os.path.join(out, 'truth.json')
truth = json.load(open(tpath)) if os.path.exists(tpath) else {}
for k in [k for k in truth if k.split('#')[0].endswith(('.nd2', '.czi'))]: del truth[k]
r = np.random.default_rng(11)

def save_truth(key, planes, meta):
    a = np.asarray(planes, np.float32); a.tofile(os.path.join(out, key + '.truth'))
    truth[key] = {'shape': list(a.shape), **meta}

# ============================================================ ND2 (v3) 쓰기 ============================================================
U32, I32, U64, F64, BOOL, STR, BYTES = 'u32', 'i32', 'u64', 'f64', 'bool', 'str', 'bytes'
def lv_item(name, v):
    nm = (name + '\0').encode('utf-16-le'); nl = len(name) + 1
    if isinstance(v, tuple):
        t, x = v
        code, body = {U32: (3, lambda x: struct.pack('<I', x)), I32: (2, lambda x: struct.pack('<i', x)), U64: (5, lambda x: struct.pack('<Q', x)),
                      F64: (6, lambda x: struct.pack('<d', x)), BOOL: (1, lambda x: struct.pack('<B', 1 if x else 0)),
                      STR: (8, lambda x: (x + '\0').encode('utf-16-le')), BYTES: (9, lambda x: struct.pack('<Q', len(x)) + bytes(x))}[t]
        return bytes([code, nl]) + nm + body(x)
    if isinstance(v, bool): return lv_item(name, (BOOL, v))
    if isinstance(v, int): return lv_item(name, (U32, v))
    if isinstance(v, float): return lv_item(name, (F64, v))
    if isinstance(v, str): return lv_item(name, (STR, v))
    if isinstance(v, (bytes, bytearray)): return lv_item(name, (BYTES, v))
    if isinstance(v, (dict, list)):
        items = list(v.items()) if isinstance(v, dict) else [('', x) for x in v]
        head = bytes([11, nl]) + nm
        kids = []; offs = []; pos = len(head) + 12
        for k, x in items:
            b = lv_item(k, x); offs.append(pos); kids.append(b); pos += len(b)
        body = b''.join(kids)
        return head + struct.pack('<IQ', len(items), len(head) + 12 + len(body)) + body + b''.join(struct.pack('<Q', o) for o in offs)
    raise TypeError(type(v))

def lv(top):
    return b''.join(lv_item(k, x) for k, x in top.items())

def spec(wl): return {'uiCount': 1, 'pPoint': {'Point0': {'dTValue': 1.0, 'dWavelength': wl, 'eType': (U32, 1)}}}
MAGIC = 0x0ABECEDA
def chunk(name, data, name_len=None):
    nb = name.encode()
    nl = name_len or max(len(nb), 32)
    return struct.pack('<IIQ', MAGIC, nl, len(data)) + nb.ljust(nl, b'\0') + data

def write_nd2(path, frames, *, H, W, nch, comps_per_ch=1, bpc=16, sig=None, cal=None, loops=(), names=None, compress=False,
              width_bytes=None, missing=(), version='Ver3.0', ne_periods=None, xy_valid=None):
    """frames: list of (H, W, nch*comps) arrays in sequence order (outer loop first)."""
    ncomp = nch * comps_per_ch; bypc = bpc // 8; wb = width_bytes or W * ncomp * bypc
    blobs = []; cmap = {}
    head = struct.pack('<IIQ', MAGIC, 32, 64) + b'ND2 FILE SIGNATURE CHUNK NAME01!' + version.encode().ljust(64, b'\0')
    blobs.append(head); pos = len(head)
    def add(name, data, name_len=None):
        nonlocal pos
        b = chunk(name, data, name_len)
        cmap[name] = (pos, len(data)); blobs.append(b); pos += len(b)
    dt = {8: np.uint8, 16: np.uint16, 32: np.float32}[bpc]
    for i, f in enumerate(frames):
        if i in missing: continue
        a = np.ascontiguousarray(np.asarray(f, dt).reshape(H, W * ncomp))
        rows = a.tobytes() if wb == W * ncomp * bypc else b''.join(row.tobytes().ljust(wb, b'\xee') for row in a)
        stamp = struct.pack('<d', 123.0 + i)
        if compress: add(f'ImageDataSeq|{i}!', stamp + zlib.compress(rows), 4072)
        else: add(f'ImageDataSeq|{i}!', stamp + rows, 4072)
    attrs = {'SLxImageAttributes': {'uiWidth': W, 'uiWidthBytes': wb, 'uiHeight': H, 'uiComp': ncomp, 'uiBpcInMemory': bpc,
             'uiBpcSignificant': sig or bpc, 'uiSequenceCount': len(frames), 'uiTileWidth': W, 'uiTileHeight': H,
             'eCompression': (U32, 0 if compress else 2), 'dCompressionParam': 0.0, 'ePixelType': (U32, 1 if bpc == 32 else 0), 'uiVirtualComponents': ncomp}}
    planes = {f'a{c}': {'uiCompCount': comps_per_ch, 'uiSampleIndex': c, 'sDescription': (names[c] if names else f'Ch{c}'),
                        'uiColor': [0xFF00, 0xFF00FF, 0xFF0000, 0xFF][c % 4], 'pFluorescentProbe': {'m_ExcitationSpectrum': spec([405.0, 488.0, 561.0, 640.0][c % 4]), 'm_EmissionSpectrum': spec([450.0, 525.0, 600.0, 700.0][c % 4])}, 'pFilterPath': {'m_pFilter': {}}} for c in range(nch)}
    meta = {'SLxPictureMetadata': {'dCalibration': float(cal or 0.0), 'bCalibrated': bool(cal), 'dObjectiveMag': 60.0, 'wsObjectiveName': (STR, 'Plan Apo λ 60x Oil'),
            'dObjectiveNA': 1.4, 'sPicturePlanes': {'uiCount': nch, 'uiCompCount': ncomp, 'sPlaneNew': planes, 'sSampleSetting': {}}}}
    add('ImageAttributesLV!', lv(attrs))
    add('ImageMetadataSeqLV|0!', lv(meta))
    def exp_level(ls):
        if not ls: return None
        t, n, extra = ls[0]
        if t == 'T': d = {'eType': 1, 'uLoopPars': {'uiCount': n, 'dStart': 0.0, 'dPeriod': extra.get('period', 1000.0), 'dDuration': 0.0}}
        elif t == 'Z': d = {'eType': 4, 'uLoopPars': {'uiCount': n, 'dZLow': 0.0, 'dZHigh': extra['step'] * (n - 1), 'dZStep': extra['step'], 'dZHome': 0.0, 'bZInverted': False, 'iType': (I32, 0)}}
        elif t == 'P':
            pts = {f'i{k:010}': {'dPosX': 100.0 * k, 'dPosY': 0.0, 'dPosZ': 0.0, 'dPosName': f'pos{k}'} for k in range(n)}
            d = {'eType': 2, 'uLoopPars': {'Points': pts, 'bUseZ': False}, 'pItemValid': bytes(xy_valid or [1] * n)}
        elif t == 'N':
            pp = {f'i{k:010}': {'uiCount': c, 'dStart': 0.0, 'dPeriod': 100.0, 'dDuration': 0.0} for k, c in enumerate(ne_periods)}
            d = {'eType': 8, 'uLoopPars': {'uiCount': len(ne_periods), 'pPeriod': pp, 'pPeriodValid': bytes([1] * len(ne_periods))}}
        nxt = exp_level(ls[1:])
        if nxt: d['ppNextLevelEx'] = [nxt]
        return d
    if loops: add('ImageMetadataLV!', lv({'SLxExperiment': exp_level(list(loops))}))
    add('ImageTextInfoLV!', lv({'SLxImageTextInfo': {'TextInfoItem_5': 'INNO synthetic', 'TextInfoItem_13': 'Plan Apo 60x'}}))
    # chunk map (마지막 40 바이트 = 서명 + filemap 위치)
    cm = b''.join(n.encode() + struct.pack('<QQ', o, s) for n, (o, s) in cmap.items())
    loc = pos
    cm += b'ND2 CHUNK MAP SIGNATURE 0000001!' + struct.pack('<Q', loc)
    blobs.append(chunk('ND2 FILEMAP SIGNATURE NAME 0001!', cm))
    open(path, 'wb').write(b''.join(blobs))

def to_app(arr, sizes, rgb=False):
    """nd2 패키지 배열(차원 순서 sizes) → 앱 plane 목록 (위치 하나), 반환 (planes, nc, nz, nt)."""
    keys = list(sizes)
    a = np.asarray(arr)
    for k in ['T', 'Z', 'C']:
        if k not in keys: a = np.expand_dims(a, 0); keys.insert(0, k)
    if rgb and 'S' not in keys: raise ValueError
    order = ['T', 'Z', 'C', 'Y', 'X'] + (['S'] if rgb else [])
    a = np.transpose(a, [keys.index(k) for k in order])
    nt, nz, nc = a.shape[:3]
    if rgb:  # RGB: 앱은 R·G·B plane 3 장 (c 자리), 3 개 × nz × nt
        a = np.moveaxis(a, -1, 2)  # T Z S C Y X?  (C=1)
        a = a.reshape(nt, nz, 3 * nc, a.shape[-2], a.shape[-1])
    return a.reshape(-1, a.shape[-2], a.shape[-1]), nc, nz, nt

import nd2
def nd2_case(name, frames, check=None, positions=1, **kw):
    p = os.path.join(out, name); write_nd2(p, frames, **kw)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with nd2.ND2File(p) as f:
            arr = f.asarray(); sizes = dict(f.sizes); vs = f.voxel_size(); rgb = f.is_rgb
            chn = [c.channel.name for c in (f.metadata.channels or [])]
    cal = kw.get('cal')
    for pi in range(positions):
        a = arr; sz = dict(sizes)
        if 'P' in sz:
            ax = list(sz).index('P'); a = np.take(arr, pi, axis=ax); sz.pop('P')
        planes, nc, nz, nt = to_app(a, sz, rgb)
        meta = {'nc': nc, 'nz': nz, 'nt': nt, 'fmt': 'nd2', 'names': chn}
        if rgb and kw.get('bpc', 16) == 8: meta['rgb'] = True; meta.pop('nc')
        elif rgb: meta['nc'] = 3; meta['names'] = ['R', 'G', 'B']   # 16-bit RGB → R·G·B 3 채널 (8-bit 로 자르지 않음)
        if cal: meta.update(pw=vs.x, unit='µm')
        if vs.z != 1 and nz > 1: meta['zstep'] = vs.z
        save_truth(name + (f'#{pi}' if positions > 1 else ''), planes, meta)
    if check: check(arr, sizes)
    print(f'{name}: nd2 패키지 sizes={sizes} voxel={tuple(vs)}')

H, W = 24, 33
def frames_u16(n, nc, sig=16): return [r.integers(0, 2 ** sig, (H, W, nc), dtype=np.uint16) for _ in range(n)]
fr = frames_u16(1, 1, 12)
nd2_case('n_single12.nd2', fr, H=H, W=W, nch=1, sig=12, cal=0.1083, names=['DAPI'],
         check=lambda a, s: np.testing.assert_array_equal(a, fr[0][..., 0]))
fr = frames_u16(3 * 4, 2)  # T 바깥 3, Z 안쪽 4, 채널 2
nd2_case('n_t3z4c2.nd2', fr, H=H, W=W, nch=2, cal=0.25, loops=[('T', 3, {}), ('Z', 4, {'step': 0.5})], names=['FITC', 'TRITC'],
         check=lambda a, s: np.testing.assert_array_equal(a[1, 2, 1], fr[1 * 4 + 2][..., 1]))
fr = frames_u16(4 * 3, 2)  # Z 바깥, T 안쪽 (드문 순서)
nd2_case('n_z4t3c2.nd2', fr, H=H, W=W, nch=2, loops=[('Z', 4, {'step': 0.3}), ('T', 3, {})],
         check=lambda a, s: np.testing.assert_array_equal(a[1, 2, 0], fr[1 * 3 + 2][..., 0]))
fr = frames_u16(3, 2)
nd2_case('n_lossless.nd2', fr, H=H, W=W, nch=2, compress=True, loops=[('T', 3, {})], cal=0.5)
fr = frames_u16(1, 1)
nd2_case('n_pad.nd2', fr, H=H, W=W, nch=1, width_bytes=W * 2 + 6)
fr = [r.integers(0, 256, (H, W, 3), dtype=np.uint8)]
nd2_case('n_rgb.nd2', fr, H=H, W=W, nch=1, comps_per_ch=3, bpc=8, cal=1.25)
fr = [r.integers(0, 4096, (H, W, 3), dtype=np.uint16)]
nd2_case('n_rgb16.nd2', fr, H=H, W=W, nch=1, comps_per_ch=3, bpc=16, sig=12)
fr = [r.normal(0, 100, (H, W, 1)).astype(np.float32) for _ in range(2)]
nd2_case('n_f32.nd2', fr, H=H, W=W, nch=1, bpc=32, loops=[('T', 2, {})])
fr = frames_u16(2 * 3, 1)  # 위치 2 × Z 3
nd2_case('n_p2z3.nd2', fr, H=H, W=W, nch=1, loops=[('P', 2, {}), ('Z', 3, {'step': 1.0})], positions=2, cal=0.65)
fr = frames_u16(3, 1)
nd2_case('n_missing.nd2', fr, H=H, W=W, nch=1, loops=[('T', 3, {})], missing={1})
fr = frames_u16(8 * 4, 1)  # 위치 8 × Z 4 → 앱에서 고르기 창
nd2_case('n_p8z4.nd2', fr, H=H, W=W, nch=1, loops=[('P', 8, {}), ('Z', 4, {'step': 0.25})], positions=8, cal=0.2)
fr = frames_u16(5, 1)
nd2_case('n_netime.nd2', fr, H=H, W=W, nch=1, loops=[('N', 5, {})], ne_periods=[2, 3])

# ============================================================ CZI ============================================================
from pylibCZIrw import czi as pyczi
import czifile
def czi_case(name, writes, ptype_note, compression=None, scale=None, names=None, scenes=False, expect=None):
    p = os.path.join(out, name)
    if os.path.exists(p): os.remove(p)
    kw = {'compression_options': compression} if compression else {}
    with pyczi.create_czi(p, exist_ok=True, **kw) as w:
        for data, plane, loc, scene in writes:
            w.write(data=data, plane=plane, location=loc, scene=scene)
        mk = {}
        if scale: mk.update(scale_x=scale * 1e-6, scale_y=scale * 1e-6)
        if names: mk['channel_names'] = {i: n for i, n in enumerate(names)}
        if mk: w.write_metadata(**mk)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        xa = czifile.imread(p, asxarray=True, squeeze=False)
    return p, ''.join(str(d)[0] for d in xa.dims), np.asarray(xa)

def czi_truth(name, planes, meta): save_truth(name, planes, {'fmt': 'czi', **meta})

# 1) Gray16, C2 Z3 T2, 비압축, 보정
H2, W2 = 20, 28
g = {(c, z, t): r.integers(0, 65536, (H2, W2), dtype=np.uint16) for c in range(2) for z in range(3) for t in range(2)}
p, axes, arr = czi_case('c_c2z3t2.czi', [(g[k][..., None], {'C': k[0], 'Z': k[1], 'T': k[2]}, (0, 0), 0) for k in g], 'Gray16', scale=0.1, names=['DAPI', 'GFP'])
sq = arr.squeeze(); ax = [a for a, n in zip(axes, arr.shape) if n > 1]
print('c_c2z3t2.czi czifile axes', axes, arr.shape, '→', ax)
T = {a: i for i, a in enumerate(ax)}; sqT = np.transpose(sq, [T['T'], T['Z'], T['C'], T['Y'], T['X']])
for k in g: assert np.array_equal(sqT[k[2], k[1], k[0]], g[k])
czi_truth('c_c2z3t2.czi', sqT.reshape(-1, H2, W2), {'nc': 2, 'nz': 3, 'nt': 2, 'pw': 0.1, 'unit': 'µm', 'names': ['DAPI', 'GFP']})
# 2) zstd1 (HiLoByteUnpack) Gray16
g2 = [r.integers(0, 65536, (H2, W2), dtype=np.uint16) for _ in range(2)]
p, axes, arr = czi_case('c_zstd1.czi', [(g2[c][..., None], {'C': c}, (0, 0), 0) for c in range(2)], 'Gray16', compression='zstd1:ExplicitLevel=2;PreProcess=HiLoByteUnpack')
sq = arr.squeeze(); assert np.array_equal(sq[0], g2[0]) and np.array_equal(sq[1], g2[1])
czi_truth('c_zstd1.czi', sq, {'nc': 2})
# 3) zstd0 Gray8
g3 = r.integers(0, 256, (H2, W2), dtype=np.uint8)
p, axes, arr = czi_case('c_zstd0_8.czi', [(g3[..., None], {}, (0, 0), 0)], 'Gray8', compression='zstd0:ExplicitLevel=1')
assert np.array_equal(arr.squeeze(), g3); czi_truth('c_zstd0_8.czi', arr.squeeze()[None], {})
# 4) Float32
g4 = r.normal(0, 3, (H2, W2)).astype(np.float32)
p, axes, arr = czi_case('c_f32.czi', [(g4[..., None], {}, (0, 0), 0)], 'Gray32Float', scale=0.05)
assert np.array_equal(arr.squeeze(), g4); czi_truth('c_f32.czi', arr.squeeze()[None], {'pw': 0.05, 'unit': 'µm'})
# 5) Bgr24 → 앱 RGB
g5 = r.integers(0, 256, (H2, W2, 3), dtype=np.uint8)  # B,G,R
p, axes, arr = czi_case('c_bgr24.czi', [(g5, {}, (0, 0), 0)], 'Bgr24')
a5 = arr.squeeze()  # libCZI 는 B,G,R 순서로 저장, czifile 은 R,G,B 로 돌려준다
assert np.array_equal(a5, g5[..., ::-1])
czi_truth('c_bgr24.czi', np.stack([a5[..., 0], a5[..., 1], a5[..., 2]]), {'rgb': True})
# 6) Mosaic: 타일 2×2 (겹침 4 px), 나중 타일이 위 (libCZI 합성 규칙 = M 순서)
tiles = []; tw, th = 16, 12; ov = 4
canvas = np.zeros((2 * th - ov, 2 * tw - ov), np.uint16)
for j in range(2):
    for i in range(2):
        t = r.integers(1, 65536, (th, tw), dtype=np.uint16); x0, y0 = i * (tw - ov) + 10, j * (th - ov) + 20
        tiles.append((t[..., None], {}, (x0, y0), 0)); canvas[y0 - 20:y0 - 20 + th, x0 - 10:x0 - 10 + tw] = t
p, axes, arr = czi_case('c_mosaic.czi', tiles, 'Gray16')
with pyczi.open_czi(p) as rd:
    bb = rd.total_bounding_rectangle; lib = rd.read(roi=(bb.x, bb.y, bb.w, bb.h))[..., 0]
assert np.array_equal(lib, canvas), 'libCZI 합성이 M 순서와 다름'
czi_truth('c_mosaic.czi', lib[None], {'mosaic': True})
# 7) Scene 2 (서로 다른 위치) → 영상 2 개
s0 = r.integers(0, 65536, (H2, W2), dtype=np.uint16); s1 = r.integers(0, 65536, (H2 + 4, W2 - 2), dtype=np.uint16)
p, axes, arr = czi_case('c_scenes.czi', [(s0[..., None], {}, (0, 0), 0), (s1[..., None], {}, (500, 300), 1)], 'Gray16', scale=0.2)
with pyczi.open_czi(p) as rd:
    sb = rd.scenes_bounding_rectangle; ims = [rd.read(scene=k)[..., 0] for k in sorted(sb)]
assert np.array_equal(ims[0], s0) and np.array_equal(ims[1], s1)
czi_truth('c_scenes.czi#0', s0[None], {'pw': 0.2, 'unit': 'µm'}); czi_truth('c_scenes.czi#1', s1[None], {'pw': 0.2, 'unit': 'µm'})

# 8) RGB 채널(Bgr24) + 회색 채널(Gray8) 섞임 → 회색 영상 1 + RGB 영상 1 (정답 = libCZI)
m_rgb = r.integers(0, 256, (H2, W2, 3), dtype=np.uint8); m_g = r.integers(0, 256, (H2, W2, 1), dtype=np.uint8)
p, axes, arr = czi_case('c_mix_rgb_gray.czi', [(m_rgb, {'C': 0}, (0, 0), 0), (m_g, {'C': 1}, (0, 0), 0)], 'mix', names=['BF', 'DAPI'])
with pyczi.open_czi(p) as rd:
    g1 = rd.read(plane={'C': 1})[..., 0]; b0 = rd.read(plane={'C': 0})
assert np.array_equal(g1, m_g[..., 0]) and np.array_equal(b0, m_rgb)
czi_truth('c_mix_rgb_gray.czi#0', g1[None], {'nc': 1, 'names': ['DAPI']})
czi_truth('c_mix_rgb_gray.czi#1', np.stack([b0[..., 2], b0[..., 1], b0[..., 0]]), {'rgb': True})
json.dump(truth, open(tpath, 'w'), ensure_ascii=False, indent=0)
print('정답', sum(1 for k in truth if k.split('#')[0].endswith(('.nd2', '.czi'))), '항목')

# ============================================================ Register 시험 stack (UI) ============================================================
# 장 k: 내용이 화면에서 반시계 θk 회전 + 이동 (dxk, dyk) — 정답을 아는 rigid drift. 기대 결과 = pystackreg 'previous' 의 transform_stack.
import tifffile
from scipy import ndimage
from pystackreg import StackReg
Hr, Wr = 96, 128; yy, xx = np.mgrid[0:Hr + 64, 0:Wr + 64].astype(float)
base = np.zeros((Hr + 64, Wr + 64))
for _ in range(60):
    cy, cx, s, a = r.uniform(0, Hr + 64), r.uniform(0, Wr + 64), r.uniform(2, 6), r.uniform(50, 200)
    base += a * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * s * s))
ang = [0.0, 1.0, 2.5, 3.0, 4.5, 6.0]; dxs = [0.0, 1.3, 2.9, 4.1, 5.2, 7.0]; dys = [0.0, -0.8, -1.5, -2.9, -3.4, -4.6]
frames = []
c0x, c0y = (Wr + 64 - 1) / 2, (Hr + 64 - 1) / 2
for th, dx, dy in zip(ang, dxs, dys):
    t = np.deg2rad(th)
    # 반시계(화면, y 아래) θ 회전: 움직인 장의 화소 p 는 기준의 R(θ)(p − c − d) + c 에서 온다 (y 아래 좌표에서 R(−θ) 가 화면 반시계)
    X = xx - c0x - dx; Y = yy - c0y - dy
    xs = np.cos(t) * X - np.sin(t) * Y + c0x; ys = np.sin(t) * X + np.cos(t) * Y + c0y
    f = ndimage.map_coordinates(base, [ys, xs], order=3, mode='constant')[32:32 + Hr, 32:32 + Wr]
    frames.append(f.astype(np.float32))
stk = np.stack(frames)
tifffile.imwrite(os.path.join(out, 'reg_rigid.tif'), stk, photometric='minisblack')
sr = StackReg(StackReg.RIGID_BODY); tm = sr.register_stack(stk.astype(np.float64), reference='previous'); exp = sr.transform_stack(stk.astype(np.float64), tmats=tm)
exp.astype(np.float32).tofile(os.path.join(out, 'reg_rigid.expected.f32'))
json.dump({'angle': ang, 'dx': dxs, 'dy': dys, 'shape': list(stk.shape)}, open(os.path.join(out, 'reg_rigid.json'), 'w'))
print('reg_rigid.tif: 6 장, 정답 각도', ang)
