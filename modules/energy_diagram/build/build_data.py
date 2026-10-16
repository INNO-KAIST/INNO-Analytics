"""ED_DATA 생성기 — research/levels 의 JSON 들을 앱 내장용으로 압축한다.
출처는 각 JSON 의 provenance 를 그대로 따른다 (README 참조)."""
import json, re, sys, os
D=os.path.join(os.path.dirname(__file__),'..','data')
ln3=json.load(open(f'{D}/ln3_laf3.json'))['levels']
cur=json.load(open(f'{D}/curated_transitions.json'))['transitions']
hp=json.load(open(f'{D}/host_phonons.json'))
ts=json.load(open(f'{D}/ts_grids.json'))

NF={'Ce':1,'Pr':2,'Nd':3,'Pm':4,'Sm':5,'Eu':6,'Gd':7,'Tb':8,'Dy':9,'Ho':10,'Er':11,'Tm':12,'Yb':13}
# UCNP 관례 표기 (Carnall 표기와 다른 경우) — 표시만 바꾸고 원래 표기는 note 로 남긴다
RELABEL={('Er','2G9/2(1)'):('2H9/2','Carnall 1989 표기 ²G₉/₂ (²G·²H·⁴F 강혼합). UCNP 문헌 관례는 ²H₉/₂ (≈410 nm 방출 준위).'),
         ('Nd','2H9/2(1)'):('2H9/2',None),('Nd','2H11/2(1)'):('2H11/2',None),('Nd','2G7/2(1)'):('2G7/2',None),
         ('Nd','2G9/2(1)'):('2G9/2',None),('Nd','2D3/2(1)'):('2D3/2',None),('Nd','2D5/2(1)'):('2D5/2',None),
         ('Er','2H11/2(1)'):('2H11/2',None),('Er','2D5/2(1)'):('2D5/2',None),('Ho','5F2(1)'):('5F2',None)}
SRC={'exp':'e','calc':'c'}
out={}
for ion,L in ln3.items():
    lv=[]
    for x in L:
        E=x['E_cm1']
        if E>52000: continue
        if x['source']=='exp' and x.get('exp_stark_range'): rg=x['exp_stark_range']
        elif x.get('carnall_calc_stark_range'): rg=x['carnall_calc_stark_range']
        else: rg=x.get('calc_stark_range') or [E,E]
        src=x['source']
        m=x.get('method','')
        code='e' if src=='exp' else ('c' if 'Carnall 1989 calculated' in m else ('f' if 'free-ion' in m or 'free ion' in m else 'q'))
        tu=x['term_unique']; lab=re.sub(r'\(\d+\)$','',tu); note=x.get('note') or ''
        if (ion,tu) in RELABEL:
            lab,n2=RELABEL[(ion,tu)]
            if n2: note=(n2+' '+note).strip()
        comp=', '.join(x.get('composition',[])[:3])
        lv.append([tu,lab,round(E),round(rg[0]),round(rg[1]),code,x.get('dominant_pct') or 0,comp,note])
    lv.sort(key=lambda r:r[2])
    out[ion]={'N':NF[ion],'lv':lv}
# 전이 목록: 준위 id 로 해석
trans=[]
for t in cur:
    ion=t['ion'].replace('3+','')
    if ion not in out: continue
    ids={}
    for r in out[ion]['lv']:                                  # 에너지 오름차순: 같은 표기면 가장 낮은 준위
        ids.setdefault(r[0],r[0]); ids.setdefault(re.sub(r'\(\d+\)$','',r[0]),r[0])
    for r in out[ion]['lv']: ids.setdefault(r[1],r[0])
    ups=[u for u in t['upper'].split('+')]
    upid=[ids.get(u) for u in ups]
    if None not in upid:                                      # 열적 결합 쌍은 낮은 준위를 대표로
        E={r[0]:r[2] for r in out[ion]['lv']}; upid.sort(key=lambda i:E[i])
    lo=ids.get(t['lower'])
    if None in upid or lo is None:
        print('skip',t['ion'],t['upper'],t['lower'],file=sys.stderr); continue
    kind={'em':'em','abs':'abs','abs/em':'both','em/abs':'both','ir':'em'}.get(t['type'],'em')
    src=t.get('typical_source') or ''
    trans.append([ion,upid,lo,kind,t.get('typical_nm'),src[:140]])
hosts=[]
for h in hp['hosts']:
    hosts.append({'host':h['host'],'w':h.get('E_ph_cm1'),'src':h.get('source'),'primary':h.get('primary'),'note':h.get('note'),'approx':bool(h.get('approx'))})
# TS grids -> int (E/B*100)
g={}
KEEP={'d3':['4A2#0','4T2#0','2E#0','2T1#0','4T1#0','2T2#0','4T1#1'],
      'd5':['6A1#0','4T1#0','4T2#0','4A1#0','4E#0','4T2#1','4E#1','4T1#1','2T2#0','2A2#0','2T1#0'],
      'd8':['3A2#0','3T2#0','1E#0','3T1#0','1A1#0','1T2#0','3T1#1']}
for k,G in ts['grids'].items():
    labs=[l for l in G['labels'] if l in KEEP[k]]
    g[k]={'labels':labs,'E':[[[int(round(v*100)) for v in row] for row in G['E'][G['labels'].index(l)]] for l in labs]}
data={'ln3':out,'trans':trans,'hosts':hosts,'ts':{'DQ':ts['DQ'],'CB':ts['CB'],'g':g},
      'refs':{
        'ln3':'W. T. Carnall, G. L. Goodman, K. Rajnak, R. S. Rana, J. Chem. Phys. 90, 3443 (1989) — LaF3:Ln³⁺ Stark 준위 표(ANL-88-8). qlanth(zia-lab) 디지털화본 + 같은 파라미터로 재계산해 다중항 귀속.',
        'ts':'Y. Tanabe, S. Sugano, J. Phys. Soc. Jpn. 9, 753 (1954); 정확한 dⁿ 대각화(팔면체, 스핀-궤도 무시)로 생성한 격자에서 보간.',
        'dorenbos':'P. Dorenbos, J. Lumin. 104, 239 (2003); D(2+)=0.64·D(3+)−0.233 eV.',
        'nist':'NIST Atomic Spectra Database (Kramida et al.), 자유 이온 준위.',
        'gap':'다중포논 이완은 준위 간격을 최대 포논 에너지로 나눈 포논 수 p 로 판정 (energy-gap law 경험칙: van Dijk & Schuurmans, J. Chem. Phys. 78, 5317 (1983)).'}}
js='/* ED_DATA — 자동 생성(build_data.py). 직접 수정 금지. */\nwindow.ED_DATA='+json.dumps(data,ensure_ascii=False,separators=(',',':'))+';\n'
open(os.path.join(os.path.dirname(__file__),'..','src','ed_data.js'),'w').write(js)
print('bytes',len(js.encode()), 'levels',sum(len(v['lv']) for v in out.values()),'trans',len(trans))
