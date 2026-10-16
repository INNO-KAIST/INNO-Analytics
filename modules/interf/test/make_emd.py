"""합성 Velox .emd — 조성·강도 참값을 아는 EDS 점 스펙트럼 (간섭 검사 라인 추천 시험용)"""
import json,math,numpy as np,h5py,sys
s=open('/home/claude/ed/base/analytics_base.html').read()
g=lambda k:json.loads(s[s.index(k)+len(k):s.index('\n',s.index(k))].rstrip().rstrip(';'))
XRL=g('const XRL = '); MAC=g('const MACDB = '); YT=json.load(open('../yield_table.json'))
NS={'K':2,'L3':4,'M5':6}
Si=MAC['mac']['Si']
def eff(E):
    m=math.exp(np.interp(math.log(E*1000),np.log(Si[0]),np.log(Si[1])))*2.33; return math.exp(-m*1e-5)*(1-math.exp(-m*0.045))
def fy(el,f,kV=200):
    Ec,wb,lvl,_=YT[el][f]; return NS.get(lvl,2)*math.log(kV*1000/Ec)/Ec*wb
def fwhm(E,fmn=129.2):
    s2=(fmn/2.355)**2-0.12*3.86*5895; return 2.355*math.sqrt(s2+0.12*3.86*E*1000)/1000
comp=json.loads(sys.argv[1]) if len(sys.argv)>1 else {"Na":1,"F":4,"Lu":0.6,"Yb":0.38,"Tm":0.02,"Cu":0.05,"C":0.5,"O":0.1}
nch,disp,off=4096,5.0,-250.0; E=(off+disp*np.arange(nch))/1000
lines=[]
for el,c in comp.items():
    for f,ls in XRL[el]['lines'].items():
        if f not in YT[el]: continue
        for nm,e,rel in ls: lines.append((el,f,nm,e,c*fy(el,f)*rel*eff(e)))
mx=max(l[4] for l in lines); SC=3e4/mx
y=np.zeros(nch)
for el,f,nm,e,I in lines:
    sg=fwhm(e)/2.355; y+=SC*I*np.exp(-0.5*((E-e)/sg)**2)/(sg*math.sqrt(2*math.pi))*disp/1000
cont=np.where(E>0.12, 4000*(20-np.clip(E,0.12,20))/np.clip(E,0.12,None)/20, 0)*np.array([eff(max(e,0.12)) for e in E])
rng=np.random.default_rng(7); cnt=rng.poisson(np.clip(y+cont,0,None)); cnt[E<0.12]=0
truth={f"{el} {f}":SC*c*fy(el,f) for el,c in comp.items() for f in YT[el] }
json.dump({"comp":comp,"scale":SC,"truth_principal_counts":{k:v for k,v in truth.items()}},open('truth.json','w'),indent=1)
md={"Detectors":{"Detector-1":{"DetectorType":"AnalyticalDetector","DetectorName":"SuperXG21","Dispersion":str(disp),"OffsetEnergy":str(off),"BeginEnergy":"120","ElectronicsNoise":"16.75"}},
    "BinaryResult":{"Detector":"SuperXG21"},"Optics":{"AccelerationVoltage":"200000"}}
b=np.frombuffer(json.dumps(md).encode(),dtype=np.uint8)
with h5py.File('synthetic_eds.emd','w') as f:
    uid='a1b2c3d4e5f6'; f.create_dataset(f'/Data/Spectrum/{uid}/Data',data=cnt.reshape(nch,1).astype(np.uint32))
    m=np.zeros((60000,1),np.uint8); m[:len(b),0]=b; f.create_dataset(f'/Data/Spectrum/{uid}/Metadata',data=m)
    f.create_dataset('/Version',data=np.frombuffer(json.dumps({"version":"0.2","format":"Velox"}).encode(),dtype=np.uint8).reshape(-1,1))
print('ok',cnt.sum(), {k:round(v) for k,v in truth.items()})
