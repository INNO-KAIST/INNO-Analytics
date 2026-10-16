"""ImageJ 1.54u(소스에서 컴파일, github.com/imagej/ImageJ) 로 기준 결과를 만든다.
사용: python3 -I gen_ref.py <ij.jar> <작업 폴더>
- 정답을 아는 합성 영상(난수 고정)을 만든다 → ImageJ 매크로로 처리 → 결과를 .f32 + manifest.json 으로.
- node test/test_img.js 가 같은 입력을 IMG_CORE 로 처리해 비교한다."""
import sys, os, json, subprocess, numpy as np, tifffile
jar, work = sys.argv[1], sys.argv[2]
IN=os.path.join(work,'in'); OUT=os.path.join(work,'out'); os.makedirs(IN,exist_ok=True); os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(20261016)
H,W=180,240
yy,xx=np.mgrid[0:H,0:W]
def blobs(n,rmin,rmax,amp,seed):
    r=np.random.default_rng(seed); img=np.zeros((H,W))
    for _ in range(n):
        cx,cy=r.uniform(0,W),r.uniform(0,H); rad=r.uniform(rmin,rmax); a=r.uniform(0.5,1)*amp
        img+=a*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*rad**2))
    return img
base=blobs(40,3,9,180,1)+20+15*np.sin(xx/37.0)+0.08*yy   # 배경 기울기 + 입자
g8=np.clip(np.round(base+rng.normal(0,6,(H,W))),0,255).astype(np.uint8)
g16=np.clip(np.round(base*180+300+rng.normal(0,900,(H,W))),0,65535).astype(np.uint16)
g32=(base/50.0-1.3+rng.normal(0,0.1,(H,W))).astype(np.float32)
# 이진: 원·타원·붙은 쌍·구멍·가장자리
b=np.zeros((H,W),bool)
def disk(cx,cy,r): return (xx-cx)**2+(yy-cy)**2<=r*r
def ell(cx,cy,a,c,t):
    ct,st=np.cos(t),np.sin(t); X=(xx-cx)*ct+(yy-cy)*st; Y=-(xx-cx)*st+(yy-cy)*ct; return (X/a)**2+(Y/c)**2<=1
b|=disk(30,30,12); b|=disk(52,34,11)            # 붙은 쌍
b|=disk(120,40,16)&~disk(120,40,6)               # 구멍
b|=disk(120,40,3)                                 # 구멍 안 입자
b|=ell(190,50,22,9,0.6)
b|=disk(5,120,9)                                  # 왼쪽 가장자리
b|=(np.abs(xx-80)<=1)&(yy>90)&(yy<160)            # 가는 막대
b|=disk(150,130,14)|disk(170,140,13)|disk(160,158,12)  # 세 개 붙음
b|=ell(220,120,6,25,0.0)
for i in range(12): b|=disk(20+i*9,170,2.2)        # 작은 점
rr=np.random.default_rng(7)
for i in range(25): b[rr.integers(0,H),rr.integers(0,W)]=True   # 단일 화소
bin8=(b*255).astype(np.uint8)
# 문턱 시험용 히스토그램 모양
def mk(kind,seed):
    r=np.random.default_rng(seed); n=H*W
    if kind=='bimodal': v=np.where(r.random(n)<0.3,r.normal(70,12,n),r.normal(170,20,n))
    elif kind=='skew': v=r.gamma(2.0,25.0,n)+10
    elif kind=='unimodal': v=r.normal(120,25,n)
    elif kind=='trimodal': v=np.choose(r.integers(0,3,n),[r.normal(40,8,n),r.normal(110,10,n),r.normal(200,12,n)])
    elif kind=='small': v=np.where(r.random(n)<0.04,r.normal(200,15,n),r.normal(60,10,n))
    return v.reshape(H,W)
thr_imgs={}
for k,kind in enumerate(['bimodal','skew','unimodal','trimodal','small']):
    v=mk(kind,100+k); thr_imgs[kind+'8']=np.clip(np.round(v),0,255).astype(np.uint8)
    thr_imgs[kind+'16']=np.clip(np.round(v*97+1234),0,65535).astype(np.uint16)
    thr_imgs[kind+'32']=(v/37.0-2.1).astype(np.float32)
imgs={'g8':g8,'g16':g16,'g32':g32,'bin8':bin8}; imgs.update(thr_imgs)
# 스택 (Z Project)
stk=np.stack([np.clip(np.round(base*0.8+rng.normal(0,10,(H,W))+k*3),0,255).astype(np.uint8) for k in range(7)])
for k,v in imgs.items():
    tifffile.imwrite(os.path.join(IN,k+'.tif'),v)
    v.astype(np.float32).tofile(os.path.join(IN,k+'.f32'))
tifffile.imwrite(os.path.join(IN,'stk.tif'),stk); stk.astype(np.float32).tofile(os.path.join(IN,'stk.f32'))
meta={k:{'w':int(v.shape[1]),'h':int(v.shape[0]),'bit':{np.uint8:8,np.uint16:16,np.float32:32}[v.dtype.type]} for k,v in imgs.items()}
meta['stk']={'w':W,'h':H,'bit':8,'n':7}

# ---------- 매크로 ----------
cases=[]   # (이름, 입력, 매크로 본문, 종류: img|thr|csv|txt)
def C(name,src,body,kind='img'): cases.append((name,src,body,kind))
for s in ['g8','g16','g32']:
    for sg in [0.8,2,5,12]: C(f'gauss_{s}_{sg}',s,f'run("Gaussian Blur...", "sigma={sg}");')
    for r in [1,2,3.5,6]: C(f'median_{s}_{r}',s,f'run("Median...", "radius={r}");')
    for f in ['Mean','Minimum','Maximum','Variance']:
        for r in [2,4.5]: C(f'{f.lower()}_{s}_{r}',s,f'run("{f}...", "radius={r}");')
    C(f'outB_{s}',s,'run("Remove Outliers...", "radius=2 threshold=%s which=Bright");'%({'g8':20,'g16':3000,'g32':0.3}[s]))
    C(f'outD_{s}',s,'run("Remove Outliers...", "radius=2 threshold=%s which=Dark");'%({'g8':20,'g16':3000,'g32':0.3}[s]))
    C(f'tophat_{s}',s,'run("Top Hat...", "radius=5");'); C(f'tophatL_{s}',s,'run("Top Hat...", "radius=5 light");')
    C(f'despeckle_{s}',s,'run("Despeckle");')
    for f in ['Smooth','Sharpen','Find Edges']: C(f'{f.replace(" ","").lower()}_{s}',s,f'run("{f}");')
    C(f'unsharp_{s}',s,'run("Unsharp Mask...", "radius=2 mask=0.60");')
    for r in [5,20,50]: C(f'rb_{s}_{r}',s,f'run("Subtract Background...", "rolling={r}");')
    C(f'rbL_{s}',s,'run("Subtract Background...", "rolling=15 light");')
    C(f'rbNS_{s}',s,'run("Subtract Background...", "rolling=8 disable");')
    C(f'rbC_{s}',s,'run("Subtract Background...", "rolling=25 create");')
    C(f'conv_{s}',s,'run("Convolve...", "text1=[1 2 1\\n2 4 2\\n1 2 1\\n] normalize");')
    C(f'bp_{s}',s,'run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 autoscale saturate");')
    C(f'bpH_{s}',s,'run("Bandpass Filter...", "filter_large=25 filter_small=2 suppress=Horizontal tolerance=5 autoscale saturate");')
    C(f'bpN_{s}',s,'run("Bandpass Filter...", "filter_large=30 filter_small=4 suppress=None tolerance=5");')
    for op,arg in [('Log',''),('Exp',''),('Square Root',''),('Square',''),('Gamma...','value=0.6'),('Multiply...','value=1.3'),('Add...','value=25'),('Subtract...','value=17.6'),('Divide...','value=2.5'),('Min...','value=60'),('Max...','value=150')]:
        nm=op.replace('...','').replace(' ','').lower(); C(f'math_{nm}_{s}',s,f'run("{op}"'+(f', "{arg}"' if arg else '')+');')
    if s!='g32': C(f'enh_{s}',s,'run("Enhance Contrast...", "saturated=0.35 normalize");'); C(f'eq_{s}',s,'run("Enhance Contrast...", "saturated=0.35 equalize");'); C(f'eqC_{s}',s,'setKeyDown("alt"); run("Enhance Contrast...", "saturated=0.35 equalize"); setKeyDown("none");')
    else: C(f'enh_{s}',s,'run("Enhance Contrast...", "saturated=0.35 normalize");')
    C(f'inv_{s}',s,'run("Invert");')
    for prom in [10,30] if s=='g8' else ([2000,6000] if s=='g16' else [0.2,0.6]):
        for out,on in [('Count','count'),('List','list'),('[Single Points]','single'),('[Maxima Within Tolerance]','within'),('[Segmented Particles]','seg')]:
            C(f'max_{s}_{prom}_{on}',s,f'run("Find Maxima...", "prominence={prom} output={out}");',kind='csv' if on in('count','list') else 'imgnew')
        C(f'max_{s}_{prom}_listX',s,f'run("Find Maxima...", "prominence={prom} exclude output=List");',kind='csv')
        C(f'max_{s}_{prom}_listS',s,f'run("Find Maxima...", "prominence={prom} strict output=List");',kind='csv')
        C(f'max_{s}_{prom}_listL',s,f'run("Find Maxima...", "prominence={prom} light output=List");',kind='csv')
    # 형 변환
    if s!='g8': C(f'to8_{s}',s,'resetMinAndMax(); run("8-bit");')
    if s=='g32': C(f'to16_{s}',s,'resetMinAndMax(); run("16-bit");')
# 문턱: 모든 방법 × dark/light
methods=["Default","Huang","Intermodes","IsoData","Li","MaxEntropy","Mean","MinError","Minimum","Moments","Otsu","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"]
for s in list(thr_imgs)+['g8','g16','g32']:
    body=''.join(f'setAutoThreshold("{m} dark"); getThreshold(l,u); print(f,"{m}|dark|"+d2s(l,9)+"|"+d2s(u,9)); setAutoThreshold("{m}"); getThreshold(l,u); print(f,"{m}|light|"+d2s(l,9)+"|"+d2s(u,9));' for m in methods)
    C(f'thr_{s}',s,body,kind='txt')
# 이진
bopt='setOption("BlackBackground", true); '
for op in ['Erode','Dilate','Open','Close']:
    for it,cnt in [(1,1),(2,3),(3,5)]:
        C(f'bin_{op.lower()}_{it}_{cnt}','bin8',bopt+f'run("Options...", "iterations={it} count={cnt} black do=[{op}]");')
# Outline: ImageJ 1.54u 는 Black background 에서 snapshot 사본을 읽어 바깥 테두리를 만든다(확인함).
# 문서 정의(물체의 안쪽 1 화소 테두리)는 흰 배경 규칙에서 나오므로 그 결과를 기준으로 쓴다.
C('bin_outlineW','bin8','run("Invert"); setOption("BlackBackground", false); run("Outline"); run("Invert");')
for op in ['Fill Holes','Skeletonize','Watershed','Ultimate Points']:
    C('bin_'+op.replace(' ','').lower(),'bin8',bopt+f'run("{op}");')
C('bin_edm8','bin8',bopt+'run("Options...", "iterations=1 count=1 black edm=Overwrite"); run("Distance Map");')
C('bin_edm32','bin8',bopt+'run("Options...", "iterations=1 count=1 black edm=32-bit"); run("Distance Map");',kind='imgnew')
C('bin_g8_mask','g8','setAutoThreshold("Otsu dark"); '+bopt+'run("Convert to Mask");')
# 입자 분석
meas='run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret\'s integrated median skewness kurtosis redirect=None decimal=9");'
for nm,opt in [('pa','size=0-Infinity'),('paX','size=0-Infinity exclude'),('paH','size=0-Infinity include'),('paS','size=20-400'),('paC','size=0-Infinity circularity=0.60-1.00'),('paXH','size=5-Infinity exclude include')]:
    C(nm,'bin8',bopt+meas+f'run("Analyze Particles...", "{opt} display clear");',kind='csv')
C('paCal','bin8',bopt+meas+'run("Properties...", "channels=1 slices=1 frames=1 pixel_width=2.5 pixel_height=2.5 voxel_depth=1 unit=nm"); run("Analyze Particles...", "size=0-Infinity display clear");',kind='csv')
C('paRed','bin8',bopt+'open("'+os.path.join(IN,'g8.tif')+'"); rename("gray"); selectImage(1); run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret\'s integrated median skewness kurtosis redirect=gray decimal=9"); run("Analyze Particles...", "size=0-Infinity display clear");',kind='csv')
# ROI 측정
for nm,mk in [('roiRect','makeRectangle(31,22,57,41);'),('roiOval','makeOval(40,30,61,45);'),('roiOval2','makeOval(20,20,50,50);'),('roiPoly','makePolygon(30,20,120,35,140,110,60,150,20,90);'),('roiFree','makeSelection("freehand", newArray(30,40,55,70,90,95,85,60,40,28), newArray(30,22,25,20,35,60,90,100,85,60));'),
                  ('roiOvalOut','makeOval(-10,5,30,20);'),('roiRectOut','makeRectangle(-10,5,30,20);'),('roiOvalOut2','makeOval(220,150,40,50);')]:
    for s in ['g8','g16','g32']:
        C(f'{nm}_{s}',s,meas+mk+'run("Measure");',kind='csv')
# 계산기·투영
C('calc_add','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Add create","g8.tif","bin8.tif");',kind='imgnew')
C('calc_sub32','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Subtract create 32-bit","g8.tif","bin8.tif");',kind='imgnew')
C('calc_diff','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Difference create","g8.tif","bin8.tif");',kind='imgnew')
C('calc_avg','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Average create","g8.tif","bin8.tif");',kind='imgnew')
C('calc_div','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Divide create","g8.tif","bin8.tif");',kind='imgnew')
C('calc_mul','g8','open("'+os.path.join(IN,'bin8.tif')+'"); imageCalculator("Multiply create","g8.tif","bin8.tif");',kind='imgnew')
for pj,nm in [('[Average Intensity]','avg'),('[Max Intensity]','max'),('[Min Intensity]','min'),('[Sum Slices]','sum'),('[Standard Deviation]','sd'),('Median','median')]:
    C('zp_'+nm,'stk',f'run("Z Project...", "projection={pj}");',kind='imgnew')

lines=['setBatchMode(true);']
for name,src,body,kind in cases:
    lines.append(f'open("{os.path.join(IN,src+".tif")}");')
    if kind=='txt':
        lines.append(f'f=File.open("{os.path.join(OUT,name+".txt")}");'); lines.append(body); lines.append('File.close(f);')
    else:
        lines.append(body)
        if kind=='csv': lines.append(f'if (nResults>0) {{ saveAs("Results","{os.path.join(OUT,name+".csv")}"); run("Clear Results"); }} else {{ File.saveString("","{os.path.join(OUT,name+".csv")}"); }}')
        else: lines.append(f'saveAs("Tiff","{os.path.join(OUT,name+".tif")}");')
    lines.append('run("Close All");')
lines.append('print("DONE"); eval("script","System.exit(0);");')
mac=os.path.join(work,'ref.ijm'); open(mac,'w').write('\n'.join(lines))
env=dict(os.environ); env['JAVA_TOOL_OPTIONS']=''
r=subprocess.run(['xvfb-run','-a','java','-Xmx3g','-cp',jar,'ij.ImageJ','-batch',mac],capture_output=True,text=True,env=env,timeout=3000)
print(r.stdout[-2000:], r.stderr[-3000:])
man={'meta':meta,'cases':[]}
for name,src,body,kind in cases:
    e={'name':name,'src':src,'kind':kind}
    if kind in('img','imgnew'):
        p=os.path.join(OUT,name+'.tif')
        if os.path.exists(p):
            a=tifffile.imread(p); e['bit']={np.uint8:8,np.uint16:16,np.float32:32,np.int16:16}[a.dtype.type]; e['w']=int(a.shape[-1]); e['h']=int(a.shape[-2])
            a.astype(np.float32).tofile(os.path.join(OUT,name+'.f32'))
        else: e['missing']=True
    elif kind=='csv':
        p=os.path.join(OUT,name+'.csv'); e['rows']=[]
        if os.path.exists(p):
            import csv
            with open(p) as fh:
                rd=csv.DictReader(fh)
                for row in rd: e['rows'].append({k.strip():(float(v) if v not in('','NaN') else None) for k,v in row.items() if k and k.strip()!=' '})
    elif kind=='txt':
        p=os.path.join(OUT,name+'.txt'); e['thr']={}
        if os.path.exists(p):
            for ln in open(p):
                q=ln.strip().split('|')
                if len(q)==4: e['thr'][q[0]+'|'+q[1]]=[float(q[2]),float(q[3])]
    man['cases'].append(e)
json.dump(man,open(os.path.join(work,'manifest.json'),'w'))
print('cases',len(cases),'missing',sum(1 for c in man['cases'] if c.get('missing')))
