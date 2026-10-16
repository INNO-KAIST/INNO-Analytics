"""실측 nd2·czi 파일의 정답 — 독립 구현이 읽은 값.
 - nd2: Python `nd2` 패키지(NIS-Elements SDK 와 대조된 구현) 의 ND2File.asarray / voxel_size / 채널 이름
 - czi: Zeiss 공식 libCZI(pylibCZIrw) 의 read(plane, scene) — mosaic 합성·scene 경계 포함, Scaling·채널 이름은 libCZI 메타데이터
실측 파일 출처는 fetch_real.sh (GitHub 공개 시험 자료). 사용: python3 -I gen_real.py <파일 폴더>"""
import sys, os, json, glob, warnings
import numpy as np
D = sys.argv[1]; warnings.simplefilter('ignore')
out = os.path.join(D, 'truth'); os.makedirs(out, exist_ok=True)
truth = {}
def save(key, planes, meta):
    a = np.ascontiguousarray(np.asarray(planes, np.float32)); a.tofile(os.path.join(out, key.replace('/', '_') + '.truth'))
    truth[key] = {'shape': list(a.shape), **meta}

# ---------------- nd2 ----------------
import nd2
CHOOSE = {'Argo_sphere_water.nd2': {'z': list(range(0, 121, 8))}}   # 1 GB z-stack: 8 장마다 (앱의 고르기 경로 시험)
for p in sorted(glob.glob(os.path.join(D, '*.nd2'))):
    fn = os.path.basename(p)
    try:
        f = nd2.ND2File(p)
    except Exception as e:
        print(fn, 'nd2 패키지도 못 읽음', e); continue
    with f:
        if f.version < (2, 0):
            truth[fn] = {'expect_error': 'v1.0'}; print(fn, 'v1 (JPEG 2000) → 앱은 읽지 않음(오류 기대)'); continue
        sizes = dict(f.sizes); vs = f.voxel_size(); rgb = f.is_rgb
        names = [c.channel.name for c in (f.metadata.channels or [])]
        calib = bool(f._rdr._cached_raw_metadata().get('bCalibrated')) if hasattr(f._rdr, '_cached_raw_metadata') else True
        ch = CHOOSE.get(fn, {})
        if 'z' in ch:
            # 큰 파일: 필요한 프레임만 read_frame 으로 (Z loop 만 있는 경우)
            assert list(sizes)[0] == 'Z' and len(f.experiment) == 1
            arr = np.stack([f.read_frame(z) for z in ch['z']]); sizes2 = {'Z': len(ch['z']), 'Y': sizes['Y'], 'X': sizes['X']}
        else:
            arr = f.asarray(); sizes2 = sizes
        keys = list(sizes2)
        P = sizes2.get('P', 1)
        for pi in range(P):
            a = np.take(arr, pi, axis=keys.index('P')) if 'P' in keys else arr
            k2 = [k for k in keys if k != 'P']
            for k in ['T', 'Z', 'C']:
                if k not in k2: a = np.expand_dims(a, 0); k2.insert(0, k)
            if rgb:
                a = np.transpose(a, [k2.index(x) for x in ['C', 'T', 'Z', 'S', 'Y', 'X']])   # 채널마다 RGB 영상
                imgs = [a[c].reshape(-1, a.shape[-2], a.shape[-1]) for c in range(a.shape[0])]
            else:
                a = np.transpose(a, [k2.index(x) for x in ['T', 'Z', 'C', 'Y', 'X']]); imgs = [a.reshape(-1, a.shape[-2], a.shape[-1])]
            nt, nz = a.shape[1:3] if rgb else a.shape[:2]
            for ci, im in enumerate(imgs):
                idx = pi * len(imgs) + ci
                meta = {'fmt': 'nd2', 'nz': int(nz), 'nt': int(nt)}
                if rgb: meta['rgb'] = True
                else: meta['nc'] = int(a.shape[2]); meta['names'] = names if len(names) == a.shape[2] else None
                if calib and vs.x != 1: meta.update(pw=vs.x, unit='µm')
                if 'Z' in sizes and vs.z != 1:
                    zz = ch.get('z'); meta['zstep'] = vs.z * ((zz[1] - zz[0]) if zz else 1)
                if ch: meta['choose'] = ch
                save(f'{fn}#{idx}', im, meta)
        print(fn, sizes, 'voxel', tuple(round(v, 5) for v in vs), 'rgb' if rgb else '', '→', P * (len(imgs)), '영상')

# ---------------- czi ----------------
from pylibCZIrw import czi as pyczi
def chan_meta(md):
    try:
        img = md['ImageDocument']['Metadata']['Information']['Image']
        chs = img['Dimensions']['Channels']['Channel']
        chs = chs if isinstance(chs, list) else [chs]
        return [c.get('@Name') or c.get('@Id') or '' for c in chs]
    except Exception:
        return None
def scaling(md):
    out = {}
    try:
        it = md['ImageDocument']['Metadata']['Scaling']['Items']['Distance']
        for d in (it if isinstance(it, list) else [it]):
            v = float(d.get('Value') or 0)
            if v > 0: out[d['@Id']] = v
    except Exception:
        pass
    return out
for p in sorted(glob.glob(os.path.join(D, '*.czi'))):
    fn = os.path.basename(p)
    try:
        rd_cm = pyczi.open_czi(p); rd = rd_cm.__enter__()
    except Exception as e:
        print(fn, 'libCZI 도 못 읽음', e); continue
    try:
        tb = rd.total_bounding_box_no_pyramid; pts = rd.pixel_types
        sc = scaling(rd.metadata); names = chan_meta(rd.metadata)
        dims = {k: v for k, v in tb.items() if k not in 'XY'}
        nc = dims.get('C', (0, 1))[1]; nz = dims.get('Z', (0, 1))[1]; nt = dims.get('T', (0, 1))[1]
        extra = [k for k, v in dims.items() if k not in ('C', 'Z', 'T', 'S', 'M') and v[1] > 1]
        scenes = sorted(rd.scenes_bounding_rectangle_no_pyramid) or [None]
        rgbtype = {'Bgr24', 'Bgr48', 'Bgr96Float'} & set(pts.values())
        import itertools
        combos = list(itertools.product(*[range(dims[k][0], dims[k][0] + dims[k][1]) for k in extra])) or [()]
        idx = 0; n_img = 0
        for s in scenes:
            bb = rd.scenes_bounding_rectangle_no_pyramid[s] if s is not None else rd.total_bounding_rectangle_no_pyramid
            for combo in combos:
                planes = {}
                for c in range(nc):
                    for z in range(nz):
                        for t in range(nt):
                            pl = {'C': dims.get('C', (0, 1))[0] + c, 'Z': dims.get('Z', (0, 1))[0] + z, 'T': dims.get('T', (0, 1))[0] + t}
                            for k, v in zip(extra, combo): pl[k] = v
                            a = rd.read(roi=(bb.x, bb.y, bb.w, bb.h), plane=pl, scene=s)
                            planes[(c, z, t)] = a
                meta = {'fmt': 'czi', 'nz': nz, 'nt': nt, 'bbox': [bb.x, bb.y, bb.w, bb.h]}
                if 'X' in sc: meta.update(pw=sc['X'] * 1e6, unit='µm')
                if 'Z' in sc and nz > 1: meta['zstep'] = sc['Z'] * 1e6
                if rgbtype:
                    for c in range(nc):
                        L = []
                        for t in range(nt):
                            for z in range(nz):
                                a = planes[(c, z, t)]
                                L += [a[..., 2], a[..., 1], a[..., 0]]   # libCZI BGR → R,G,B
                        m2 = dict(meta)
                        if pts[0] == 'Bgr24': m2['rgb'] = True
                        else: m2['nc'] = 3
                        save(f'{fn}#{idx}', np.stack(L), m2); idx += 1
                else:
                    L = [planes[(c, z, t)][..., 0] for t in range(nt) for z in range(nz) for c in range(nc)]
                    m2 = dict(meta, nc=nc)
                    if names and len(names) == nc and all(names): m2['names'] = names
                    save(f'{fn}#{idx}', np.stack(L), m2); idx += 1
        print(fn, {k: v[1] for k, v in dims.items()}, 'scenes', len(scenes), 'extra', extra, pts[0] if pts else '', '→', idx, '영상')
    except Exception as e:
        # libCZI 가 못 읽는 형식(Bgr96Float 등) → 다른 독립 구현 czifile(Gohlke) 로 (단일 plane 만)
        try:
            import czifile
            a = np.asarray(czifile.imread(p))
            assert a.ndim == 3 and a.shape[-1] == 3, a.shape   # (Y, X, RGB) — czifile 은 R,G,B 순서
            sc2 = scaling(rd.metadata); meta = {'fmt': 'czi', 'nz': 1, 'nt': 1, 'nc': 3, 'truth_by': 'czifile'}
            if 'X' in sc2: meta.update(pw=sc2['X'] * 1e6, unit='µm')
            save(f'{fn}#0', np.stack([a[..., 0], a[..., 1], a[..., 2]]), meta); print(fn, 'libCZI 불가', repr(e)[:80], '→ czifile 정답', a.shape, a.dtype)
        except Exception as e2:
            print(fn, '정답 생성 실패', repr(e)[:200], repr(e2)[:200])
    finally:
        rd_cm.__exit__(None, None, None)
json.dump(truth, open(os.path.join(out, 'truth.json'), 'w'), ensure_ascii=False, indent=0)
print('정답', len(truth))
