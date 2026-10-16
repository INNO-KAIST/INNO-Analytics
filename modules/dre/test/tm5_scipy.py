"""Tm 5/6준위 — 방정식을 손으로 적은 독립 구현 (scipy). JS 엔진 결과(tm5_ref_js.json)와 대조."""
import json,numpy as np
from scipy.integrate import solve_ivp
R=json.load(open('tm5_ref_js.json'))
A=lambda v: [v[str(i)] for i in range(len(v))] if isinstance(v,dict) else v
for L in R:
    R[L]['up']=[A(v) for v in R[L]['up']]; R[L]['kin']['n']=[A(v) for v in R[L]['kin']['n']]; EPH=6.62607015e-34*2.99792458e8/1064e-9
def rhs(t,n,p,phi,six):
    n1,n2,n3,n4,n5=n[:5]; n6=n[5] if six else 0.0
    N=p['NRE']*p['c']/100
    W2=p['w2r']+p['w2nr']; CR=p['cCR']*N; U23=p['cUC23']*N; U22=p['cUC22']*N
    gsa=p['sG']*phi*n1; e1=p['sE1']*phi*n2; e2=p['sE2']*phi*n4
    cr=CR*n4*n1; u23=U23*n2*n2; u22=U22*n2*n2
    w35=p['w35']; w6=p['w6']
    d=np.zeros(6 if six else 5)
    # 3H6
    d[0]=-gsa+W2*n2+w35*p['b51']*n3+(1-p['b32'])*p['w3r']*n4+0.30*w6*n5-cr+u23+u22
    # 3F4
    d[1]=-e1-W2*n2+w35*(1-p['b51'])*n3+p['b32']*p['w3r']*n4+0.10*w6*n5+2*cr-2*u23-2*u22
    # 3H5
    d[2]=gsa-w35*n3+p['w3nr']*n4+0.45*w6*n5+u22
    # 3H4
    d[3]=e1-e2-(p['w3r']+p['w3nr'])*n4+0.15*w6*n5-cr+u23
    # 1G4
    d[4]=e2-w6*n5
    if six:
        c2=p['cCR2']*N*n4*n5; w7=p['w7']
        d[3]+=-c2+0.15*w7*n6; d[1]+=c2+0.46*w7*n6; d[4]+=-c2+0.01*w7*n6; d[0]+=0.36*w7*n6; d[2]+=0.02*w7*n6
        d[5]=c2-w7*n6
    return d
worst=0
for L in ['5','6']:
    r=R[L]; p=r['params']; six=L=='6'; nv=6 if six else 5
    # 올림 스윕: 이전 해에서 장시간 적분 (LSODA)
    y=np.zeros(nv); y[0]=1; errs=[]
    for I,nj in zip(r['I'],r['up']):
        phi=I/EPH; s=solve_ivp(rhs,[0,50],y,args=(p,phi,six),method='LSODA',rtol=1e-10,atol=1e-14); y=s.y[:,-1]
        errs.append(max(abs(y[k]-nj[k])/max(abs(nj[k]),1e-9) for k in range(nv) if nj[k]>1e-8))
    k=r['kin']; s=solve_ivp(rhs,[0,k['t'][-1]],np.eye(nv)[0],args=(p,k['I']/EPH,six),method='LSODA',rtol=1e-10,atol=1e-15,t_eval=k['t'])
    ek=max(abs(s.y[j,i]-k['n'][i][j])/max(abs(k['n'][i][j]),1e-12) for i in range(len(k['t'])) for j in range(nv) if k['n'][i][j]>1e-10)
    print(f'{L}준위: 정상상태 최대 상대차 {max(errs):.2e}, 동역학 최대 상대차 {ek:.2e}'); worst=max(worst,max(errs),ek)
print('PASS' if worst<1e-4 else 'FAIL', worst)
