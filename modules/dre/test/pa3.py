"""Lee et al. (2021) 3-level Tm PA 모델 — 공개 코드 구조를 그대로 옮긴 정상상태/동역학 기준 구현 (연구용 검증 스크립트)."""
import numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.integrate import solve_ivp
h=6.62607015e-34; c0=2.99792458e8; lam=1064e-9; hnu=h*c0/lam
def params_lee(ct=8.0,f2=0.25,f3=0.5):
    return dict(sG=0.005*0.12e-25, sE=2.0*3.2e-25, w2r=162.60, w2nr=162.60*f2+0.00239, w3r=636.01, w3nr=6.02+162.60*f3,
                b32=0.144, s31=0.1*1.6*ct**2*1000, Q22=0.08*320*ct**3/(ct**2+4.3**2), Q23=0.1*90*ct**3/(ct**2+4.3**2))
def rhs(n,p,I):
    n1,n2,n3=n; phi=I/hnu
    w2=p['w2r']+p['w2nr']; w3r=p['w3r']; w3nr=p['w3nr']
    d1= -p['sG']*phi*n1 + w2*n2 + (1-p['b32'])*w3r*n3 + (p['Q22']+p['Q23'])*n2*n2 - p['s31']*n1*n3
    d2=  p['sG']*phi*n1 - p['sE']*phi*n2 - w2*n2 + (p['b32']*w3r+w3nr)*n3 - (p['Q22']+2*p['Q23'])*n2*n2 + 2*p['s31']*n1*n3
    d3=  p['sE']*phi*n2 - (w3r+w3nr)*n3 + p['Q23']*n2*n2 - p['s31']*n1*n3
    return np.array([d1,d2,d3])
def steady(p,I,guess):
    f=lambda x: rhs(np.array([1-x[0]-x[1],x[0],x[1]]),p,I)[1:]
    x,info,ier,_=fsolve(f,guess,full_output=True,xtol=1e-13)
    return x,ier==1
def curve(p,I_arr,direction=1):
    """continuation(자연 파라미터): 앞 점의 해를 다음 점의 초기값으로 — up/down 방향별"""
    idx=range(len(I_arr)) if direction>0 else range(len(I_arr)-1,-1,-1)
    out=np.zeros((len(I_arr),2)); g=np.array([1e-9,1e-12]) if direction>0 else np.array([0.3,0.05])
    for i in idx:
        x,ok=steady(p,I_arr[i],g); out[i]=x; g=x
    return out

def all_steady(p,I,ngrid=4000):
    """모든 정상해: d3=0 을 n3 에 대한 2차식으로 풀어(두 근 모두) n1=1-n2-n3 를 대입 → d2(n2)=0 의 근을 두 가지(branch)에서 모두 찾음.
    판별식이 음수인 구간은 해가 없으므로 건너뛴다(이전 판본은 0 으로 잘라 가짜 근을 만들었음)."""
    from scipy.optimize import brentq
    phi=I/hnu; w3=p['w3r']+p['w3nr']
    def n3_of(n2,sgn):
        a=p['s31']; b=-(w3+p['s31']*(1-n2)); cc=p['sE']*phi*n2+p['Q23']*n2*n2
        disc=b*b-4*a*cc
        if disc<0: return None
        return (-b+sgn*np.sqrt(disc))/(2*a)
    def g(n2,sgn):
        n3=n3_of(n2,sgn)
        if n3 is None: return None
        n1=1-n2-n3
        if n3<0 or n1<-1e-15: return None
        return rhs(np.array([n1,n2,n3]),p,I)[1]
    xs=np.logspace(-15,np.log10(0.99999),ngrid)
    roots=[]
    for sgn in (-1,1):
        gs=[g(x,sgn) for x in xs]
        for i in range(len(xs)-1):
            if gs[i] is None or gs[i+1] is None: continue
            if np.sign(gs[i])!=np.sign(gs[i+1]):
                r=brentq(lambda z:g(z,sgn),xs[i],xs[i+1],xtol=1e-18,rtol=1e-14); roots.append((r,n3_of(r,sgn)))
    roots.sort()
    return roots
