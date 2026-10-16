"""합성 데이터(참값 앎)로 DRE 피팅 전략 비교 — 지역 최소 갇힘 빈도와 식별성."""
import numpy as np, json, time, warnings
from scipy.optimize import least_squares
from scipy.stats import qmc
from pa3 import params_lee
from fwd import forward
warnings.filterwarnings('ignore')
rng=np.random.default_rng(7)
TRUE=params_lee()
FREE=['sG','sE','w2nr','w3nr','s31','Q22','Q23']       # w2r, w3r, b32 고정(문헌·JO) → 절대 시간 스케일 고정
I=np.logspace(np.log10(0.3e7),np.log10(60e7),45)
S_true=3.0e6; bg=40.0                                    # 검출 배율, 암전류 바닥(counts)
y0=S_true*forward(TRUE,I)[0]
y_obs=np.maximum(y0+bg+rng.normal(0,np.sqrt(y0+bg))*1.0,1.0)*np.exp(rng.normal(0,0.05,len(I)))   # 샷 노이즈 + 5% 곱 노이즈
th_true=np.log10([TRUE[k] for k in FREE])
LB=th_true-2.0; UB=th_true+2.0                          # 탐색 상자: 참값 ±2 decade
def model(th,branch='low'):
    p=dict(TRUE); p.update({k:10**v for k,v in zip(FREE,th)}); m,nr=forward(p,I,branch); return m
def r_log_hier(th):
    m=model(th)
    if not np.all(np.isfinite(m)) or np.any(m<=0): return np.full(len(I),1e3)
    # 계층 최적화: log y = log(S m + bg) — 배경은 알려진 암전류로 고정, S 는 1차원 닫힌해가 없으므로 log 오프셋 근사 후 1회 뉴턴
    ly=np.log(np.maximum(y_obs-bg,1.0)); lm=np.log(m)
    w=(y_obs-bg)>3*np.sqrt(bg)                             # 배경 3σ 이하 점은 하한 정보로만(여기선 제외)
    logS=np.mean((ly-lm)[w])
    r=np.where(w,ly-lm-logS,0.0)
    return r
def r_log_joint(x):
    th=x[:-1]; logS=x[-1]; m=model(th)
    if not np.all(np.isfinite(m)) or np.any(m<=0): return np.full(len(I),1e3)
    ly=np.log(np.maximum(y_obs-bg,1.0)); w=(y_obs-bg)>3*np.sqrt(bg)
    return np.where(w,ly-np.log(m)-logS,0.0)
def r_lin_joint(x):
    th=np.log10(np.maximum(x[:-1],1e-300)); S=x[-1]; m=model(th)
    if not np.all(np.isfinite(m)): return np.full(len(I),1e9)
    return (S*m+bg-y_obs)/np.maximum(y_obs,1)**0.5
N=60
lhs=qmc.LatinHypercube(d=len(FREE),seed=3).random(N); starts=LB+lhs*(UB-LB)
out={}
def run(name,fun,x0s,**kw):
    res=[];t=time.time()
    for x0 in x0s:
        try: r=least_squares(fun,x0,**kw); res.append((2*r.cost,r.x,r.nfev))
        except Exception as e: res.append((np.inf,x0,0))
    out[name]={'costs':sorted(float(c) for c,_,_ in res),'time':time.time()-t,'nfev_med':float(np.median([n for _,_,n in res]))}
    return res
# S1: 선형 매개변수 + 선형 잔차 + S 동시 (흔한 방식)
x0s=[np.r_[10**s,S_true*10**rng.uniform(-1,1)] for s in starts]
R1=run('S1_linear',r_lin_joint,x0s,method='trf',x_scale='jac',max_nfev=400)
# S2: log 매개변수 + log 잔차 + S 동시
x0s=[np.r_[s,np.log(S_true)+rng.uniform(-2,2)] for s in starts]
R2=run('S2_log_joint',r_log_joint,x0s,method='trf',max_nfev=400)
# S3: log 매개변수 + log 잔차 + S 계층 소거
R3=run('S3_log_hier',r_log_hier,list(starts),method='trf',bounds=(LB-1,UB+1),max_nfev=400)

chi_true=float(np.sum(r_log_hier(th_true)**2))
for k,v in out.items():
    c=np.array(v['costs']); cmin=min(chi_true,min(out['S3_log_hier']['costs']))
    v['success_frac']=float(np.mean(c<=cmin*1.05+1e-3)) if k!='S1_linear' else None
    print(k,'time %.1fs'%v['time'],'nfev_med',v['nfev_med'],'best %.4g'%c.min(),'succ',v['success_frac'])
# S1 성공 판정은 log 잔차로 재평가
c1=[float(np.sum(r_log_hier(np.log10(np.maximum(x[:-1],1e-300)))**2)) for _,x,_ in R1]
cmin=min(chi_true,min(out['S3_log_hier']['costs']))
out['S1_linear']['success_frac']=float(np.mean(np.array(c1)<=cmin*1.05+1e-3))
print('S1 succ (log χ² 재평가)',out['S1_linear']['success_frac'],' chi_true',chi_true,' cmin',cmin)
# 최적해 파라미터 복원
bestx=min(R3,key=lambda r:r[0])[1]
out['recovery_S3']={k:float(v) for k,v in zip(FREE,bestx-th_true)}
print('recovery (log10 best − true):',{k:round(v,3) for k,v in out['recovery_S3'].items()})
json.dump({'chi_true':chi_true,**out},open('fitstudy.json','w'),indent=1,default=float)
np.save('best_S3.npy',bestx)
