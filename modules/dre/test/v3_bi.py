import numpy as np, json, subprocess
from pa3 import all_steady, rhs
e={"k31":1830,"sG":4.55e-24,"sE":5.2e-20,"w2nr":2.34,"w3nr":3.29,"kQ22":25500,"kQ23":63.7,"b32":0.066}
c=8.0
p=dict(sG=e['sG']*1e-4,sE=e['sE']*1e-4,w2r=162.60,w2nr=e['w2nr'],w3r=636.01,w3nr=e['w3nr'],b32=e['b32'],s31=e['k31']*c*c,Q22=e['kQ22']*c**3/(c*c+4.3**2),Q23=e['kQ23']*c**3/(c*c+4.3**2))
I=np.logspace(np.log10(8e4),np.log10(2.2e5),25)
res=[]
for x in I:
    r=all_steady(p,x*1e4,ngrid=20000)
    # 안정성: 축약 Jacobian 고유값
    st=[]
    for n2,n3 in r:
        h=1e-7; f0=rhs(np.array([1-n2-n3,n2,n3]),p,x*1e4)[1:]
        J=np.zeros((2,2))
        for j in range(2):
            d=np.zeros(2); d[j]=h*max(abs([n2,n3][j]),1e-12)
            y=np.array([n2,n3])+d; f1=rhs(np.array([1-y[0]-y[1],y[0],y[1]]),p,x*1e4)[1:]
            J[:,j]=(f1-f0)/d[j]
        st.append(bool(np.all(np.linalg.eigvals(J).real<0)))
    res.append({'I':x,'roots':[[a,b] for a,b in r],'stable':st})
json.dump({'e':e,'c':c,'res':res},open('v3_ref.json','w'))
for q in res: print(round(q['I']), len(q['roots']), q['stable'], [round(b,5) for a,b in q['roots']])
