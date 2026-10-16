"""합성 Leica LIF (경계 사례) + 정답 — test/io/lif/ 에 쓴다.
정답: Bio-Formats 7.3.1 (ref/bf/BFDump.java) 이 읽은 series 별 화소 (BF 가 읽는 파일만). 나머지는 기대 동작(건너뜀·알림)을 synth.json 에.
사용: python3 -I ref/gen_lif_synth.py <bioformats_package.jar> <BFDump 컴파일 폴더>"""
import sys, os, json, struct, subprocess, tempfile
import numpy as np
JAR, OUT = sys.argv[1:3]
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test', 'io', 'lif'); os.makedirs(D, exist_ok=True)
rng = np.random.default_rng(20261016)
def u16(s): return s.encode('utf-16le')
def write(fn, images, version=2, tail=b'', cut=0, drop=()):
    """images: dict(name, dims=[(DimID, n, BytesInc, Length, Unit)], chs=[(Resolution, BytesInc, ChannelTag, LUTName[, DataType])], data=bytes, path=[folders], extra='')"""
    els, blocks = [], [('MemBlock_0', b'')]
    for k, im in enumerate(images):
        bid = f'MemBlock_{k + 1}'
        dims = ''.join(f'<DimensionDescription DimID="{d[0]}" NumberOfElements="{d[1]}" Origin="0" Length="{d[3]}" Unit="{d[4]}" BitInc="0" BytesInc="{d[2]}"/>' for d in im['dims'])
        chs = ''.join(f'<ChannelDescription DataType="{c[4] if len(c) > 4 else 0}" ChannelTag="{c[2]}" Resolution="{c[0]}" NameOfMeasuredQuantity="" Min="0" Max="255" Unit="" LUTName="{c[3]}" IsLUTInverted="0" BytesInc="{c[1]}" BitInc="0">{c[5] if len(c) > 5 else ""}</ChannelDescription>' for c in im['chs'])
        e = (f'<Element Name="{im["name"]}" Visibility="1" CopyOption="1" UniqueID="u{k}"><Data><Image TextDescription="">{im.get("extra", "")}<ImageDescription><Channels>{chs}</Channels>'
             f'<Dimensions>{dims}</Dimensions></ImageDescription></Image></Data><Memory Size="{len(im["data"])}" MemoryBlockID="{bid}"/><Children/></Element>')
        for f in reversed(im.get('path', [])):
            fb = f'MemBlock_f{k}{f}'
            e = f'<Element Name="{f}" Visibility="1" CopyOption="1" UniqueID="f{k}{f}"><Data><Experiment/></Data><Memory Size="0" MemoryBlockID="{fb}"/><Children>{e}</Children></Element>'
            blocks.append((fb, b''))
        els.append(e)
        if k not in drop: blocks.append((bid, im['data']))
    xml = f'<LMSDataContainerHeader Version="{version}"><Element Name="root" Visibility="1" CopyOption="1" UniqueID="r"><Data><Experiment/></Data><Memory Size="0" MemoryBlockID="MemBlock_0"/><Children>{"".join(els)}</Children></Element></LMSDataContainerHeader>'
    out = bytearray(); x = u16(xml); out += struct.pack('<IIBI', 0x70, len(x) + 5, 0x2A, len(xml)) + x
    for bid, data in blocks:
        out += struct.pack('<IIBQBI' if version == 2 else '<IIBIBI', 0x70, 0, 0x2A, len(data), 0x2A, len(bid)) + u16(bid) + data
    out += tail
    if cut: out = out[:-cut]
    open(os.path.join(D, fn), 'wb').write(bytes(out))
def gray(name, W, H, nz=1, nt=1, nc=1, bps=1, path=(), extra='', chx=None):
    n = W * H * nz * nt * nc; a = rng.integers(0, 250 if bps == 1 else 60000, n).astype('<u1' if bps == 1 else '<u2')
    # 채널은 면 단위로 나뉘어 (Leica 형식: X·Y 다음 채널, 그다음 Z·T)
    dims = [(1, W, bps, (W - 1) * 0.2e-6, 'm'), (2, H, bps * W, (H - 1) * 0.2e-6, 'm')]
    if nz > 1: dims.append((3, nz, bps * nc * W * H, (nz - 1) * 1e-6, 'm'))
    if nt > 1: dims.append((4, nt, bps * nc * W * H * nz, (nt - 1) * 2.0, 's'))
    chs = [(8 if bps == 1 else 16, bps * W * H * c, 0, ['Red', 'Green', 'Blue', 'Gray'][c % 4]) + ((chx[c],) if chx else ()) for c in range(nc)]
    if chx: chs = [(c[0], c[1], c[2], c[3], 0, c[4]) for c in chs]
    return dict(name=name, dims=dims, chs=chs, data=a.tobytes(), path=list(path), extra=extra)
spec = {}
# 1) RGB 8-bit (BGR 저장, ChannelTag 3·2·1), 줄 끝 채움 (W=5 → 15 바이트를 16 으로)
W, H = 5, 4; row = 16; a = rng.integers(0, 256, row * H).astype('<u1')
write('rgb8.lif', [dict(name='rgb8', dims=[(1, W, 3, 0.8e-6, 'm'), (2, H, row, 0.6e-6, 'm')], chs=[(8, 0, 3, 'Blue'), (8, 1, 2, 'Green'), (8, 2, 1, 'Red')], data=a.tobytes())])
spec['rgb8.lif'] = dict(bf=True, rgb=True)
# 2) RGB 16-bit
W, H = 6, 3; a = rng.integers(0, 65535, W * H * 3).astype('<u2')
write('rgb16.lif', [dict(name='rgb16', dims=[(1, W, 6, 1e-6, 'm'), (2, H, 6 * W, 1e-6, 'm')], chs=[(12, 0, 3, 'Blue'), (12, 2, 2, 'Green'), (12, 4, 1, 'Red')], data=a.tobytes())])
spec['rgb16.lif'] = dict(bf=True, chNames=['R', 'G', 'B'])
# 3) 가운데 series 의 memory block 없음 → 나머지 둘 + 알림
write('missing.lif', [gray('good1', 7, 5), gray('gone', 7, 5), gray('good2', 7, 5)], drop=(1,))
spec['missing.lif'] = dict(images=2, warn='gone')
# 4) 채널마다 자료형이 다른 series → 건너뜀
mixed = gray('mixed', 6, 4, nc=2); mixed['chs'] = [(8, 0, 0, 'Red'), (16, 1, 0, 'Green')]
write('mixed.lif', [gray('ok1', 6, 4), mixed, gray('ok2', 6, 4)])
spec['mixed.lif'] = dict(images=2, warn='mixed')
# 5) 끝이 잘린 파일 → 마지막 series 의 잘린 부분 = 0 (Bio-Formats 와 같음)
t1, t2 = gray('t1', 9, 7, nz=3), gray('t2', 9, 7, nz=3); write('trunc.lif', [t1, t2], cut=10)
# Bio-Formats 7.3.1 은 이 파일에서 멈춤(initFile seek) → 정답 = 만든 자료, 잘린 마지막 10 바이트 = 0
for k, im in enumerate([t1, t2]):
    a = np.frombuffer(im['data'], '<u1').astype(np.float32).copy()
    if k == 1: a[-10:] = 0
    a.tofile(os.path.join(D, f'trunc.lif#{k}.f32'))
spec['trunc.lif'] = dict(series=[dict(key=f'trunc.lif#{k}', w=9, h=7, nz=3, nc=1, nt=1, rgb=1, ch='') for k in range(2)], warn='잘려')
# 6) 같은 잎 이름 → 경로로 구분
write('dup.lif', [gray('Position001', 4, 3, path=['Mark_and_Find_001']), gray('Position001', 4, 3, path=['Mark_and_Find_002'])])
spec['dup.lif'] = dict(names=['dup.lif · Mark_and_Find_001/Position001', 'dup.lif · Mark_and_Find_002/Position001'])
# 7) XT (세로 = 시간) → 화소 높이 = 너비 + 알림
W, H = 8, 6; a = rng.integers(0, 255, W * H).astype('<u1')
write('xt.lif', [dict(name='xt', dims=[(1, W, 1, 7 * 0.5e-6, 'm'), (4, H, W, 5 * 0.25, 's')], chs=[(8, 0, 0, 'Gray')], data=a.tobytes())])
spec['xt.lif'] = dict(pw=0.5, ph=0.5, warn='세로 축이 공간이 아닙니다')
# 8) 크기 다른 series (개관 + 시간 연속) → Z·T 를 가장 큰 값으로 묻고, series 마다 그 안의 것만
write('mixsize.lif', [gray('overview', 8, 4), gray('timelapse', 8, 4, nt=20)])
spec['mixsize.lif'] = dict(bf=True, choose=dict(p=[0, 1], t=[0, 1, 2]), nt=[1, 3])
# 9) 차원 값이 맞지 않음 (NumberOfElements 0, BytesInc 0) → 건너뜀
z0 = gray('z0', 4, 3, nz=2); z0['dims'][2] = (3, 0, 12, 0, 'm'); noinc = gray('noinc', 4, 3); noinc['dims'][0] = (1, 4, 0, 1e-6, 'm')
write('baddims.lif', [gray('fine', 4, 3), z0, noinc])
spec['baddims.lif'] = dict(images=1, warn='NumberOfElements')
# 10) 채널 이름 = ChannelProperty DyeName
props = ['<ChannelProperty><Key>DyeName</Key><Value>Leica/DAPI</Value></ChannelProperty>', '<ChannelProperty><Key>DyeName</Key><Value>Alexa 488</Value></ChannelProperty>']
write('dyes.lif', [gray('dyes', 5, 4, nc=2, chx=props)])
spec['dyes.lif'] = dict(bf=True, chNames=['DAPI', 'Alexa 488'])
# 11) 잘못된 문자 참조 (&#x110000;) 가 있는 속성
im = gray('bad&#x110000;name', 4, 3); write('badref.lif', [im])
spec['badref.lif'] = dict(images=1)

def bf(fp):
    td = tempfile.mkdtemp(); r = subprocess.run(['java', '-Xmx2g', '-cp', JAR + ':' + OUT, 'BFDump', fp, td, '0,1,2,3,4,5'], capture_output=True, text=True, env=dict(os.environ, JAVA_TOOL_OPTIONS=''))
    S = {}
    for l in r.stdout.splitlines():
        if l.startswith('S '): t = l.split(); S[int(t[1])] = dict(zip(t[2::2], t[3::2]))
    return S, td, r
for fn, sp in spec.items():
    if not sp.get('bf'): continue
    S, td, r = bf(os.path.join(D, fn)); assert S, (fn, r.stderr[-500:])
    sp['series'] = []
    for s, d in sorted(S.items()):
        X, Y, Z, C, T, rgb = (int(d[k]) for k in ['X', 'Y', 'Z', 'C', 'T', 'rgb'])
        a = np.fromfile(os.path.join(td, f's{s}.f32'), '<f4'); key = f'{fn}#{s}'; a.tofile(os.path.join(D, key + '.f32'))
        sp['series'].append(dict(key=key, w=X, h=Y, nz=Z, nc=C, nt=T, rgb=rgb, ch=d.get('ch', '')))
    print(fn, [(x['w'], x['h'], x['nc'], x['nz'], x['nt'], x['rgb'], x['ch']) for x in sp['series']])
json.dump(spec, open(os.path.join(D, 'synth.json'), 'w'), ensure_ascii=False, indent=1)
print('합성 LIF', len(spec))
