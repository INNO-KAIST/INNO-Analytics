import numpy as np, json, warnings, time
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
warnings.filterwarnings('ignore')
from fitcore import *
from pa3 import rhs
best=np.load('best_S3.npy'); chi_best=float(np.sum(r_log_hier(best)**2))
# 노이즈 수준: 로그 잔차 표준편차 추정(자유도 보정) → Δχ² 임계를 정규화
w=(y_obs-bg)>3*np.sqrt(bg); dof=int(w.sum())-len(FREE)-1; s2=chi_best/dof
thr=3.84*s2   # 95% 점별 신뢰 (χ²₁) — 잔차 분산으로 스케일
prof={}
t0=time.time()
for j,k in enumerate(FREE):
    grid=best[j]+np.linspace(-1.5,1.5,13); vals=[]
    for g in grid:
        idx=[i for i in range(len(FREE)) if i!=j]
        def f(sub): th=best.copy(); th[idx]=sub; th[j]=g; return r_log_hier(th)
        bestc=np.inf
        for x0 in [best[idx], best[idx]+np.random.default_rng(int(g*1000)%999).normal(0,0.5,len(idx))]:
            r=least_squares(f,x0,method='trf',max_nfev=300); bestc=min(bestc,2*r.cost)
        vals.append(bestc)
    vals=np.array(vals); ok=grid[vals-chi_best<=thr]
    prof[k]={'grid':(grid-th_true[j]).tolist(),'dchi':((vals-chi_best)/s2).tolist(),'ci_log10_width':float(ok.max()-ok.min()) if len(ok) else None,
             'open_low': bool(vals[0]-chi_best<=thr),'open_high':bool(vals[-1]-chi_best<=thr)}
    print(k,'CI width %.2f dec'%prof[k]['ci_log10_width'],'open',prof[k]['open_low'],prof[k]['open_high'],'%.0fs'%(time.time()-t0))
# 정상상태는 같게 맞는 두 매개변수 집합(참값 vs 최적)의 동역학 비교: 95% 상승 시간
def params_of(th): p=dict(TRUE); p.update({k:10**v for k,v in zip(FREE,th)}); return p
def rise95(p,Ival):
    sol=solve_ivp(lambda t,n: rhs(n,p,Ival),[0,3.0],[1,0,0],method='Radau',rtol=1e-8,atol=1e-15,dense_output=True)
    t=np.logspace(-6,np.log10(3.0),4000); y=sol.sol(t)[2]; return float(t[np.argmax(y>=0.95*y[-1])])
Ps=np.array([2,5,7,10,20])*1e7
rt={'I_kWcm2':(Ps/1e7).tolist(),'true':[rise95(params_of(th_true),P) for P in Ps],'best':[rise95(params_of(best),P) for P in Ps]}
print('rise95 true',np.round(rt['true'],4)); print('rise95 best',np.round(rt['best'],4))
# FIM 고유값 스펙트럼 (sloppiness)
eps=1e-3; J=np.array([(r_log_hier(best+eps*np.eye(len(FREE))[i])-r_log_hier(best-eps*np.eye(len(FREE))[i]))/(2*eps) for i in range(len(FREE))]).T
ev=np.linalg.eigvalsh(J.T@J); print('FIM eig (log10):',np.round(np.log10(np.maximum(ev,1e-30)),1))
vecs=np.linalg.eigh(J.T@J)[1]; stiff=vecs[:,-1]; sloppy=vecs[:,0]
print('stiffest dir:',dict(zip(FREE,np.round(stiff,2)))); print('sloppiest dir:',dict(zip(FREE,np.round(sloppy,2))))
json.dump({'profile':prof,'rise':rt,'fim_log10_eig':np.log10(np.maximum(ev,1e-30)).tolist(),'stiff':dict(zip(FREE,stiff.tolist())),'sloppy':dict(zip(FREE,sloppy.tolist())),'chi_best':chi_best,'s2':s2},open('profile.json','w'),indent=1)
