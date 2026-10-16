"""실측 회절 파일 → float32 원자료 (시험용, ref/real 안에만; 공개 저장소에서 제외). 사용: python3 -I ref/gen_real.py
DM3 는 rosettasciio (독립 리더) 로, RGB TIFF 는 Pillow 로 밝기(Rec.601) — 앱의 리더와 다른 구현."""
import os, json, numpy as np
from PIL import Image
import rsciio.digitalmicrograph as dm
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'real'); meta={}
for f in sorted(os.listdir(D)):
    p=os.path.join(D,f)
    if f.endswith('.dm3'):
        d=dm.file_reader(p)[0]; a=np.asarray(d['data'],dtype=np.float32); ax=d['axes'][0]; meta[f]=dict(w=a.shape[1],h=a.shape[0],scale=ax['scale'],unit=str(ax['units']))
    elif f.lower().endswith(('.tif','.tiff')):
        im=np.asarray(Image.open(p)).astype(np.float64); a=(0.299*im[...,0]+0.587*im[...,1]+0.114*im[...,2]).astype(np.float32) if im.ndim==3 else im.astype(np.float32); meta[f]=dict(w=a.shape[1],h=a.shape[0])
    else: continue
    a.tofile(p+'.f32'); print(f,meta[f])
json.dump(meta,open(os.path.join(D,'real.json'),'w'),indent=1)
