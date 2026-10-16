"""정상상태 n3(I) — (수정판 2026-09-29) 연속 스윕 + Newton/선탐색, 실패 시 암시적 Euler 의사시간.
이전 판본은 d3=0 을 2차식 작은 근으로만 소거하고 판별식을 0 으로 잘라, 높은 파워에서 가짜 근을 돌려주는 결함이 있었다."""
import numpy as np
from pa3 import hnu
def _fJ(p,phi,y):
    n2,n3=y; n1=1-n2-n3
    w2=p['w2r']+p['w2nr']; w3r=p['w3r']; w3nr=p['w3nr']; b=p['b32']; s=p['s31']; Q22=p['Q22']; Q23=p['Q23']; G=p['sG']*phi; E=p['sE']*phi
    f2=G*n1-E*n2-w2*n2+(b*w3r+w3nr)*n3-(Q22+2*Q23)*n2*n2+2*s*n1*n3
    f3=E*n2-(w3r+w3nr)*n3+Q23*n2*n2-s*n1*n3
    # d/dn2, d/dn3 (n1=1-n2-n3)
    J=np.array([[-G-E-w2-2*(Q22+2*Q23)*n2-2*s*n3, -G+(b*w3r+w3nr)+2*s*(n1-n3)],
                [E+2*Q23*n2+s*n3, -(w3r+w3nr)-s*(n1-n3)]])
    mag=np.array([G*n1+E*n2+w2*n2+(b*w3r+w3nr)*n3+(Q22+2*Q23)*n2*n2+2*s*n1*n3, E*n2+(w3r+w3nr)*n3+Q23*n2*n2+s*n1*n3])
    return np.array([f2,f3]),J,mag
def _valid(y): return y[0]>=0 and y[1]>=0 and y[0]+y[1]<=1
def _newton(p,phi,y,maxit=40,tol=1e-11):
    f,J,mag=_fJ(p,phi,y); sc=mag+1e-300; mf=np.sum((f/sc)**2)
    for it in range(maxit):
        if np.max(np.abs(f)/(mag+1e-300))<tol: return y,True
        try: d=np.linalg.solve(J,-f)
        except np.linalg.LinAlgError: return y,False
        t=1.0; acc=False
        for _ in range(30):
            yn=y+t*d
            if _valid(yn):
                fn,Jn,mn=_fJ(p,phi,yn); m2=np.sum((fn/sc)**2)
                if m2<mf*(1-1e-4*t) or np.max(np.abs(fn)/(mn+1e-300))<tol: y,f,J,mag,mf=yn,fn,Jn,mn,m2; acc=True; break
            t*=0.5
        if not acc: return y,False
    return y, np.max(np.abs(f)/(mag+1e-300))<tol
def steady(p,phi,y0):
    y,ok=_newton(p,phi,np.array(y0,float))
    if ok: return y,True
    y=np.array(y0,float); rmax=max(p['sE']*phi,p['sG']*phi,p['w3r']+p['w3nr'],p['s31'],1.0); dt=1.0/rmax
    for step in range(600):
        z=y.copy(); okst=False
        for it in range(25):
            f,J,_=_fJ(p,phi,z); A=np.eye(2)-dt*J; r=y-z+dt*f
            d=np.linalg.solve(A,r); t=1.0
            for k in range(2):
                if z[k]+d[k]<0: t=min(t,0.9*z[k]/(-d[k]))
            z=np.maximum(z+t*d,0)
            if np.max(np.abs(d)/(z+1e-18))<1e-9: okst=_valid(z); break
        if not okst: dt*=0.25; continue
        y=z; dt*=2
        if step%10==9:
            yy,ok=_newton(p,phi,y,20)
            if ok: return yy,True
    return _newton(p,phi,y,60)
def forward(p,I_arr,branch='low'):
    order=np.argsort(I_arr)
    if branch!='low': order=order[::-1]
    out=np.empty(len(I_arr)); nroots=np.ones(len(I_arr),int)
    y=np.array([0.0,0.0]) if branch=='low' else None
    if y is None:
        y,_=steady(p,I_arr[order[0]]/hnu,[0,0])
    for k in order:
        y,ok=steady(p,I_arr[k]/hnu,y); out[k]=y[1] if ok else np.nan
    return out,nroots
