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
