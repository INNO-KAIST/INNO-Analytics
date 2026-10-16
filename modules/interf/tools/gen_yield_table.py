"""interf_suggest.js 에 박아 넣는 계열별 yield 표를 만든다 (xraydb 4.5.8 = Elam, Ravel & Sieber 2002).
el → fam → [흡수단 E eV, 형광수율 × branching(주선), 초기준위, 주선]. 주선 = XRL 에서 rel 이 가장 큰 선.
사용: python3 gen_yield_table.py <analytics html> > yield_table.json"""
import json,sys,xraydb as x
s=open(sys.argv[1],encoding='utf-8').read(); i=s.index('const XRL = ')+12
XRL=json.loads(s[i:s.index('\n',i)].rstrip().rstrip(';')); out={}
for el,r in XRL.items():
    xl=x.xray_lines(el); d={}
    for fam,lines in r['lines'].items():
        nm=max(lines,key=lambda q:q[2])[0]; L=xl.get(nm)
        if L is None: continue
        lvl=L.initial_level.split(',')[0]; ed=x.xray_edge(el,lvl)
        if ed is None: continue
        d[fam]=[round(ed.energy,1),round(ed.fyield*L.intensity,6),lvl,nm]
    out[el]=d
print(json.dumps(out,separators=(',',':')))
