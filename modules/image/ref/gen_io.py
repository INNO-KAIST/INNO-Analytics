"""파일 읽기 시험용 정답 파일을 만든다 (tifffile·pypng·numpy·직접 쓴 IBW/SXM). 정답 = 같은 난수 배열.
사용: python3 -I gen_io.py <출력 폴더>"""
import sys, os, json, struct, numpy as np, tifffile, png
out=sys.argv[1]; os.makedirs(out,exist_ok=True)
r=np.random.default_rng(5); H,W=37,53
a8=r.integers(0,256,(H,W),dtype=np.uint8); a16=r.integers(0,65536,(H,W),dtype=np.uint16); a32=r.normal(0,10,(H,W)).astype(np.float32)
rgb=r.integers(0,256,(H,W,3),dtype=np.uint8); stk=r.integers(0,65536,(4,H,W),dtype=np.uint16)
truth={}
def T(name,arr,**kw):
    tifffile.imwrite(os.path.join(out,name),arr,**kw); truth[name]=None
def save_truth(name,arr,meta=None):
    np.asarray(arr,np.float32).tofile(os.path.join(out,name+'.truth')); truth[name]={'shape':list(np.asarray(arr).shape),**(meta or {})}
for nm,arr,kw in [('t8.tif',a8,{}),('t16.tif',a16,{}),('t32.tif',a32,{}),('t16be.tif',a16,{'byteorder':'>'}),
                  ('t8_lzw.tif',a8,{'compression':'lzw'}),('t16_lzw_pred.tif',a16,{'compression':'lzw','predictor':True}),
                  ('t32_deflate_fpred.tif',a32,{'compression':'zlib','predictor':3}),('t16_packbits.tif',a16,{'compression':'packbits'}),
                  ('t8_tiled.tif',a8,{'tile':(16,16)}),('t32_tiled_deflate.tif',a32,{'tile':(16,32),'compression':'zlib'}),
                  ('t16_big.tif',a16,{'bigtiff':True}),('trgb.tif',rgb,{'photometric':'rgb'}),('trgb_planar.tif',np.moveaxis(rgb,2,0),{'photometric':'rgb','planarconfig':'separate'}),
                  ('trgb_lzw.tif',rgb,{'photometric':'rgb','compression':'lzw','predictor':True}),
                  ('tstack16.tif',stk,{'photometric':'minisblack'}),('t8_strips.tif',a8,{'rowsperstrip':5})]:
    tifffile.imwrite(os.path.join(out,nm),arr,**kw); save_truth(nm,arr)
# 1-bit
b1=(a8>127)
tifffile.imwrite(os.path.join(out,'t1bit.tif'),np.packbits(b1,axis=1),shape=(H,W),dtype=np.uint8,bitspersample=1,photometric='minisblack',compression=None) if False else None
# ImageJ hyperstack + 보정 (tifffile imagej=True)
hs=r.integers(0,65536,(3,2,2,H,W),dtype=np.uint16)  # T Z C Y X
tifffile.imwrite(os.path.join(out,'tij_hyper.tif'),hs,imagej=True,resolution=(1/0.125,1/0.125),metadata={'axes':'TZCYX','unit':'um','spacing':0.5})
save_truth('tij_hyper.tif',hs,{'nc':2,'nz':2,'nt':3,'pw':0.125,'unit':'µm'})
# OME-TIFF
om=r.integers(0,256,(2,H,W),dtype=np.uint8)
tifffile.imwrite(os.path.join(out,'tome.ome.tif'),om,metadata={'axes':'CYX','PhysicalSizeX':0.25,'PhysicalSizeXUnit':'µm','PhysicalSizeY':0.25,'PhysicalSizeYUnit':'µm'})
save_truth('tome.ome.tif',om,{'nc':2,'pw':0.25,'unit':'µm'})
# OME DimensionOrder XYZCT (C 바깥, Z 안쪽) — 앱 순서 c+nc·z 로 다시 놓아야 함
om2=r.integers(0,65536,(2,3,H,W),dtype=np.uint16)
tifffile.imwrite(os.path.join(out,'tome_czyx.ome.tif'),om2,metadata={'axes':'CZYX'})
save_truth('tome_czyx.ome.tif',om2.transpose(1,0,2,3),{'nc':2,'nz':3})
# FEI 태그 (34682) + 데이터 막대
fei="[User]\nDate=10/16/2026\n[System]\nSystemType=Helios G4\n[Beam]\nHV=5000\n[Scan]\nPixelWidth=2.5e-09\nPixelHeight=2.5e-09\n[Image]\nResolutionX=53\nResolutionY=30\n[PrivateFei]\nDatabarHeight=7\n"
tifffile.imwrite(os.path.join(out,'tfei.tif'),a8,extratags=[(34682,'s',0,fei,True)]); save_truth('tfei.tif',a8,{'pw':2.5,'unit':'nm','databar':7})
zeiss="0\r\nAP_PIXEL_SIZE\r\nImage Pixel Size = 3.906 nm\r\nAP_MAG\r\nMag = 25.00 K X\r\n"
tifffile.imwrite(os.path.join(out,'tzeiss.tif'),a8,extratags=[(34118,'s',0,zeiss,True)]); save_truth('tzeiss.tif',a8,{'pw':3.906,'unit':'nm'})
# Hitachi 옆 파일
tifffile.imwrite(os.path.join(out,'thitachi.tif'),a8); open(os.path.join(out,'thitachi.txt'),'w').write("[SemImageFile]\nInstructName=SU8230\nPixelSize=1.984375\nMagnification=50000\n"); save_truth('thitachi.tif',a8,{'pw':1.984375,'unit':'nm'})
# TIFF 해상도 태그만 (적용하지 않고 후보)
tifffile.imwrite(os.path.join(out,'tres.tif'),a8,resolution=(2000,2000),resolutionunit='CENTIMETER'); save_truth('tres.tif',a8,{'cand':5.0})
# PNG
def P(name,arr,**kw):
    with open(os.path.join(out,name),'wb') as f: png.Writer(**kw).write(f,arr)
P('p8.png',a8.tolist(),width=W,height=H,greyscale=True,bitdepth=8); save_truth('p8.png',a8)
P('p16.png',a16.tolist(),width=W,height=H,greyscale=True,bitdepth=16); save_truth('p16.png',a16)
P('p16i.png',a16.tolist(),width=W,height=H,greyscale=True,bitdepth=16,interlace=True); save_truth('p16i.png',a16)
P('prgb.png',rgb.reshape(H,W*3).tolist(),width=W,height=H,greyscale=False,bitdepth=8); save_truth('prgb.png',rgb)
a4=(a8>>4)
P('p4.png',a4.tolist(),width=W,height=H,greyscale=True,bitdepth=4); save_truth('p4.png',np.round(a4.astype(float)*255/15))
# NPY
np.save(os.path.join(out,'n32.npy'),a32); save_truth('n32.npy',a32)
np.save(os.path.join(out,'n16stk.npy'),stk); save_truth('n16stk.npy',stk)
# 텍스트
np.savetxt(os.path.join(out,'txt.csv'),a32[:10,:12],delimiter=',',fmt='%.6f'); save_truth('txt.csv',np.round(a32[:10,:12],6))
# IBW v5 (Igor): 3 층 (Height, Amplitude, Phase), x=dim0, y=dim1
nx,ny,nl=W,H,3; data=r.normal(0,1e-9,(nl,ny,nx)).astype(np.float32)   # 우리 배열: [l][y_igor][x]
note=b"ScanSize: 1e-06\rImagingMode: AC Mode\r"
labels=[b'',b'',b'Height'.ljust(32,b'\0')[:32]+b'HeightRetrace'.ljust(32,b'\0')+b'AmplitudeRetrace'.ljust(32,b'\0')+b'PhaseRetrace'.ljust(32,b'\0'),b'']
labels[2]=b'\0'*32+b'HeightRetrace'.ljust(32,b'\0')+b'AmplitudeRetrace'.ljust(32,b'\0')+b'PhaseRetrace'.ljust(32,b'\0')
dimEU=[b'',b'',b'',b'']
bh=struct.pack('<hhiiii4i4iiii',5,0,320+nx*ny*nl*4,0,len(note),0,*[len(x) for x in dimEU],*[len(x) for x in labels],0,0,0)
assert len(bh)==64
wh=bytearray(320)
struct.pack_into('<i',wh,12,nx*ny*nl); struct.pack_into('<h',wh,16,2)
wh[28:28+8]=b'afmwave\0'  # WaveHeader5.bname @28
struct.pack_into('<4i',wh,68,nx,ny,nl,0)
struct.pack_into('<4d',wh,84,1e-6/nx,1e-6/ny,1,1)
wh[152:156]=b'm\0\0\0'; wh[156:160]=b'm\0\0\0'
igor=np.zeros(nx*ny*nl,np.float32)
for l in range(nl):
    for y in range(ny): igor[l*nx*ny+y*nx:l*nx*ny+y*nx+nx]=data[l,y]
with open(os.path.join(out,'afm.ibw'),'wb') as f: f.write(bh+bytes(wh)+igor.tobytes()+note+b''.join(dimEU)+b''.join(labels))
for l in range(nl): save_truth(f'afm.ibw#{l}',data[l][::-1],{'pw':1e-6/nx*1e9,'unit':'nm'})
# SXM (Nanonis)
sx,sy=W,H; ch=r.normal(0,1,(2,2,sy,sx)).astype('>f4')  # [채널][fwd/bwd][y][x] (저장 순서 그대로)
hdr=(":NANONIS_VERSION:\n2\n:SCAN_PIXELS:\n       %d       %d\n:SCAN_RANGE:\n           2.000000E-7           1.400000E-7\n:SCAN_DIR:\nup\n:BIAS:\n0.5\n"
     ":DATA_INFO:\n\tChannel\tName\tUnit\tDirection\tCalibration\tOffset\n\t14\tZ\tm\tboth\t1.0\t0.0\n\t0\tCurrent\tA\tboth\t1.0\t0.0\n\n:SCANIT_END:\n\n\n")%(sx,sy)
with open(os.path.join(out,'scan.sxm'),'wb') as f: f.write(hdr.encode()+b'\x1a\x04'+ch.tobytes())
k=0
for c in range(2):
    for d in range(2):
        img=ch[c,d].astype(np.float32)[::-1]   # up: 첫 줄이 아래
        if d==1: img=img[:,::-1]
        save_truth(f'scan.sxm#{k}',img,{'pw':2e-7/sx*1e9,'unit':'nm'}); k+=1
json.dump(truth,open(os.path.join(out,'truth.json'),'w'),indent=0)
print('ok',len(truth))
