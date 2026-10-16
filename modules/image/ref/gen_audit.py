"""전수 검수(2026-10-16)에서 찾은 파일 읽기 결함의 회귀 시험 파일 + 정답. 사용: python3 -I gen_audit.py <test/io/audit>"""
import sys, os, json, zlib, struct, numpy as np, tifffile, png
out = sys.argv[1]; os.makedirs(out, exist_ok=True); T = {}
r = np.random.default_rng(3); H, W = 21, 34
def save(name, arr, **meta):
    np.asarray(arr, np.float32).tofile(os.path.join(out, name + '.truth')); T[name] = {'shape': list(np.asarray(arr).shape), **meta}
# 1) ImageJ 한 채널 + 표시 범위
a = r.integers(0, 4000, (H, W), dtype=np.uint16); tifffile.imwrite(os.path.join(out, 'ij_minmax.tif'), a, imagej=True, metadata={'min': 100.0, 'max': 3000.0})
save('ij_minmax.tif', a, display=[100, 3000])
# 2) OME 두 series (C=2 각각, 화소 크기 다름)
s1 = r.integers(0, 256, (2, H, W), dtype=np.uint8); s2 = r.integers(0, 256, (2, H, W), dtype=np.uint8)
with tifffile.TiffWriter(os.path.join(out, 'ome_2series.ome.tif'), ome=True) as tw:
    tw.write(s1, metadata={'axes': 'CYX', 'PhysicalSizeX': 0.2, 'PhysicalSizeY': 0.2}); tw.write(s2, metadata={'axes': 'CYX', 'PhysicalSizeX': 0.4, 'PhysicalSizeY': 0.4})
save('ome_2series.ome.tif#0', s1, nc=2, pw=0.2, unit='µm'); save('ome_2series.ome.tif#1', s2, nc=2, pw=0.4, unit='µm')
# 3) OME PhysicalSizeZ
z = r.integers(0, 256, (3, H, W), dtype=np.uint8); tifffile.imwrite(os.path.join(out, 'ome_z.ome.tif'), z, metadata={'axes': 'ZYX', 'PhysicalSizeX': 0.1, 'PhysicalSizeZ': 0.5})
save('ome_z.ome.tif', z, nz=3, zstep=0.5)
# 4) NPY (H,W,3) uint16 / float32 → 3 채널
n16 = r.integers(0, 65536, (H, W, 3), dtype=np.uint16); np.save(os.path.join(out, 'rgb_u2.npy'), n16); save('rgb_u2.npy', np.moveaxis(n16, 2, 0), nc=3)
# 5) 텍스트: ';' + ',' 소수점, NaN/inf, 줄 끝 구분자
m = np.round(r.normal(0, 50, (6, 5)), 3)
open(os.path.join(out, 'semi_comma.txt'), 'w').write('x;y;z;u;v\n' + '\n'.join(';'.join(f'{v:.3f}'.replace('.', ',') for v in row) for row in m) + '\n'); save('semi_comma.txt', m)
m2 = m.copy(); m2[1, 2] = np.nan; m2[3, 0] = np.inf; m2[4, 4] = -np.inf
open(os.path.join(out, 'naninf.csv'), 'w').write('\n'.join(','.join({True: 'NaN'}.get(np.isnan(v), 'inf' if v == np.inf else '-inf' if v == -np.inf else f'{v:.3f}') for v in row) + ',' for row in m2) + '\n'); save('naninf.csv', m2)
# 6) ImageJ signed 16-bit (calibration function 0, c0=-32768, c1=1)
sv = r.integers(-2000, 2000, (H, W)).astype(np.int32); raw = (sv + 32768).astype(np.uint16)
tifffile.imwrite(os.path.join(out, 'ij_signed16.tif'), raw, description='ImageJ=1.54u\ncf=0\nc0=-32768.0\nc1=1.0\nvunit=Gray Value\nmin=0\nmax=65535\n', metadata=None)
save('ij_signed16.tif', sv)
# 7) 4-bit 팔레트 → RGB
idx = r.integers(0, 16, (H, W), dtype=np.uint8); cm = np.zeros((3, 16), np.uint16)
for k in range(16): cm[:, k] = [(k * 4099) % 65536, (k * 12007) % 65536, (65535 - k * 4000) % 65536]
def tiff_pal4(path, idx, cm):   # 4-bit 팔레트 TIFF 를 직접 씀 (tifffile 은 4-bit 팔레트 쓰기를 받지 않음)
    h, w = idx.shape; rb = (w + 1) // 2; data = bytearray(h * rb)
    for y in range(h):
        for x in range(w):
            data[y * rb + x // 2] |= int(idx[y, x]) << (4 if x % 2 == 0 else 0)
    cmb = struct.pack('<%dH' % cm.size, *cm.reshape(-1).tolist()); off_data = 8; off_cm = off_data + len(data); off_ifd = off_cm + len(cmb)
    tags = [(256, 3, 1, w), (257, 3, 1, h), (258, 3, 1, 4), (259, 3, 1, 1), (262, 3, 1, 3), (273, 4, 1, off_data), (277, 3, 1, 1), (278, 3, 1, h), (279, 4, 1, len(data)), (320, 3, cm.size, off_cm)]
    ifd = struct.pack('<H', len(tags)) + b''.join(struct.pack('<HHI', c, t, n) + (struct.pack('<HH', v, 0) if t == 3 and n == 1 else struct.pack('<I', v)) for c, t, n, v in tags) + struct.pack('<I', 0)
    open(path, 'wb').write(b'II*\x00' + struct.pack('<I', off_ifd) + bytes(data) + cmb + ifd)
tiff_pal4(os.path.join(out, 'pal4.tif'), idx, cm)
rgb = np.stack([cm[c][idx] >> 8 for c in range(3)]); save('pal4.tif', rgb, rgb=True)
# 8) signed int64
i64 = r.integers(-3_000_000, 3_000_000, (H, W)).astype(np.int64); tifffile.imwrite(os.path.join(out, 'i64.tif'), i64); save('i64.tif', i64)
# 9) ImageJ RGB 스택(이어진 자료)
rs = r.integers(0, 256, (3, H, W, 3), dtype=np.uint8); tifffile.imwrite(os.path.join(out, 'ij_rgbstack.tif'), rs, imagej=True, photometric='rgb')
save('ij_rgbstack.tif', np.concatenate([np.moveaxis(f, 2, 0) for f in rs]), rgb=True, nz=3)
# 10) 잘린 TIFF·PNG → 분명한 오류
b = open(os.path.join(out, 'ij_minmax.tif'), 'rb').read(); big = r.integers(0, 65536, (200, 300), dtype=np.uint16); tifffile.imwrite(os.path.join(out, '_full.tif'), big)
bb = open(os.path.join(out, '_full.tif'), 'rb').read(); open(os.path.join(out, 'trunc.tif'), 'wb').write(bb[:len(bb) - 1500]); os.remove(os.path.join(out, '_full.tif')); T['trunc.tif'] = {'expect_error': '잘렸거나'}
with open(os.path.join(out, '_p.png'), 'wb') as f: png.Writer(width=W, height=H, greyscale=True, bitdepth=16).write(f, r.integers(0, 65536, (H, W)).tolist())
pb = open(os.path.join(out, '_p.png'), 'rb').read(); os.remove(os.path.join(out, '_p.png'))
# IDAT 안의 zlib 을 짧은(그러나 올바른) 자료로 바꿔 넣음 → 행 수 부족
i = pb.index(b'IDAT') - 4; ln = struct.unpack('>I', pb[i:i + 4])[0]; rawz = zlib.decompress(pb[i + 8:i + 8 + ln]); z2 = zlib.compress(rawz[:len(rawz) // 2])
chunk = struct.pack('>I', len(z2)) + b'IDAT' + z2 + struct.pack('>I', zlib.crc32(b'IDAT' + z2) & 0xffffffff)
open(os.path.join(out, 'short.png'), 'wb').write(pb[:i] + chunk + pb[i + 12 + ln:]); T['short.png'] = {'expect_error': '짧습니다'}
json.dump(T, open(os.path.join(out, 'truth.json'), 'w'), ensure_ascii=False, indent=0); print(len(T), '항목')
