/* ============================================================
   INNO Analytics — XPS 계산 엔진 (DOM 없음: 브라우저·Web Worker·node 공용)
   수식의 정의는 modules/xps/SPEC.md. 이 파일은 그 구현이다.
   ============================================================ */
(function(root){
"use strict";
const LN2=Math.LN2, PI=Math.PI, SQPI=Math.sqrt(Math.PI), S2=Math.SQRT2;
const FW2SIG=1/(2*Math.sqrt(2*LN2));     // FWHM → σ

/* ============================================================
   특수함수
   ============================================================ */
/* Faddeeva w(z), Im z ≥ 0 — Weideman (SIAM J. Numer. Anal. 31, 1497, 1994), N = 32 */
const WD=(function(){
  const N=32, M=2*N, M2=2*M, L=Math.sqrt(N/Math.SQRT2);
  const f=new Float64Array(M2); // [0, f(k=-M+1..M-1)]
  for(let k=-M+1;k<=M-1;k++){ const th=k*PI/M, t=L*Math.tan(th/2); f[k+M]=Math.exp(-t*t)*(L*L+t*t); }
  f[0]=0;
  const g=new Float64Array(M2); for(let i=0;i<M2;i++) g[i]=f[(i+M2/2)%M2];   // fftshift
  const a=new Float64Array(N+1);
  for(let j=1;j<=N;j++){ let re=0; for(let m=0;m<M2;m++) re+=g[m]*Math.cos(2*PI*j*m/M2); a[j]=re/M2; }
  return {N,L,a};
})();
function faddeeva(x,y){   // 반환 [Re w, Im w], y ≥ 0
  const L=WD.L, a=WD.a, N=WD.N;
  // iz = -y + i x ; L - iz = (L + y) - i x ; L + iz = (L - y) + i x
  const dr=L+y, di=-x;                      // L - i z
  const nr=L-y, ni=x;                       // L + i z
  const den=dr*dr+di*di;
  const Zr=(nr*dr+ni*di)/den, Zi=(ni*dr-nr*di)/den;    // Z = (L+iz)/(L-iz)
  let pr=a[N], pi=0;
  for(let j=N-1;j>=1;j--){ const tr=pr*Zr-pi*Zi+a[j], ti=pr*Zi+pi*Zr; pr=tr; pi=ti; }
  // 1/(L-iz)
  const ir=dr/den, ii=-di/den;
  const i2r=ir*ir-ii*ii, i2i=2*ir*ii;      // 1/(L-iz)^2
  const wr=2*(pr*i2r-pi*i2i)+ir/SQPI, wi=2*(pr*i2i+pi*i2r)+ii/SQPI;
  return [wr,wi];
}
function erfcx(x){   // exp(x²)·erfc(x)
  if(x>=0) return faddeeva(0,x)[0];
  return 2*Math.exp(x*x)-faddeeva(0,-x)[0];
}
const LG=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
function lgamma(x){
  if(x<0.5) return Math.log(PI/Math.abs(Math.sin(PI*x)))-lgamma(1-x);
  x-=1; let a=LG[0]; const t=x+7.5; for(let i=1;i<9;i++) a+=LG[i]/(x+i);
  return 0.5*Math.log(2*PI)+(x+0.5)*Math.log(t)-t+Math.log(a);
}

/* ============================================================
   1) Lineshape (넓이 1, SPEC §1)
   ============================================================ */
const SHAPES={
  GLS:{name:"GL 합 (pseudo-Voigt · Avantage Sum)",p:["F","m"]},
  GLP:{name:"GL 곱 (Avantage Product · Casa GL)",p:["F","m"]},
  V:{name:"Voigt (G⊗L 합성곱)",p:["wG","wL"]},
  DS:{name:"Doniach–Šunjić ⊗ G (금속 비대칭)",p:["alpha","wL","wG"]},
  LA:{name:"LA(α,β) ⊗ G (Casa LA 형태)",p:["F","alpha","beta","wG"]},
  GLT:{name:"GL 합 + 지수 꼬리",p:["F","m","T","tau"]},
};
const PDEF={   // 기본 한계·이름·단위
  be:{lab:"BE",u:"eV"}, area:{lab:"넓이",u:"CPS·eV",lo:0,hi:Infinity},
  F:{lab:"FWHM",u:"eV",lo:0.05,hi:10}, m:{lab:"L 비율",u:"",lo:0,hi:1},
  wG:{lab:"G FWHM",u:"eV",lo:0,hi:8}, wL:{lab:"L FWHM",u:"eV",lo:0.005,hi:8},
  alpha:{lab:"α",u:"",lo:0,hi:0.45}, beta:{lab:"β",u:"",lo:0.55,hi:20},
  T:{lab:"꼬리 비율",u:"",lo:0,hi:1}, tau:{lab:"꼬리 길이",u:"eV",lo:0.02,hi:20},
  dE:{lab:"Δ(SO)",u:"eV",lo:0.05,hi:40}, r:{lab:"넓이 비",u:"",lo:0.01,hi:5}, q:{lab:"폭 배율",u:"",lo:0.3,hi:3},
  c0:{lab:"배경 c₀",u:"CPS",lo:-Infinity,hi:Infinity}, k:{lab:"Shirley k",u:"1/eV",lo:0,hi:Infinity}, s:{lab:"배경 기울기",u:"CPS/eV",lo:-Infinity,hi:Infinity},
};
const LA_ALPHA_LO=0.55;
function gauss(x,F){ return (2/F)*Math.sqrt(LN2/PI)*Math.exp(-4*LN2*x*x/(F*F)); }
function lorentz(x,F){ return (2/(PI*F))/(1+4*x*x/(F*F)); }
function glpNorm(F,m){ const a=4*LN2*(1-m), b=4*m; if(b<1e-12) return F*Math.sqrt(PI/a); if(a<1e-14) return F*PI/Math.sqrt(b); return F*(PI/Math.sqrt(b))*erfcx(Math.sqrt(a/b)); }
function voigt(x,wG,wL){
  if(wG<1e-9) return lorentz(x,Math.max(wL,1e-12));
  const sg=wG*FW2SIG;
  if(wL<1e-12) return gauss(x,wG);
  const w=faddeeva(x/(sg*S2),(wL/2)/(sg*S2));
  return w[0]/(sg*Math.sqrt(2*PI));
}
function laNorm(F,a,b){ return (F/4)*SQPI*(Math.exp(lgamma(a-0.5)-lgamma(a))+Math.exp(lgamma(b-0.5)-lgamma(b))); }
function dsBase(x,alpha,g){ return Math.cos(PI*alpha/2+(1-alpha)*Math.atan(-x/g))/Math.pow(x*x+g*g,(1-alpha)/2); }
/* DS 원시함수: d/dx = dsBase.  α>0: −r^α sin(αθ − πα/2)/α, θ = atan(−x/γ);  α→0: π/2 + atan(x/γ) */
function dsPrim(x,alpha,g){ const th=Math.atan(-x/g); if(alpha<1e-9) return PI/2-th; const r=Math.sqrt(x*x+g*g); return -Math.pow(r,alpha)*Math.sin(alpha*th-PI*alpha/2)/alpha; }
function glsVal(x,F,m){ return (1-m)*gauss(x,F)+m*lorentz(x,F); }

/* 격자에서 Gaussian 합성곱 (핵 ±5σ, 이산 정규화) */
function gconv(v,h,wG){
  const sg=wG*FW2SIG; const K=Math.ceil(5*sg/h); if(K<1) return v;
  const ker=new Float64Array(2*K+1); let s=0; for(let k=-K;k<=K;k++){ const t=k*h/sg; ker[k+K]=Math.exp(-0.5*t*t); s+=ker[k+K]; } for(let k=0;k<ker.length;k++) ker[k]/=s;
  const n=v.length, out=new Float64Array(n);
  for(let i=0;i<n;i++){ let a=0; const j0=Math.max(0,i-K), j1=Math.min(n-1,i+K); for(let j=j0;j<=j1;j++) a+=v[j]*ker[i-j+K]; out[i]=a; }
  return out;
}
/* 높은 BE 쪽 지수 꼬리 합성곱: out(x) = ∫₀^∞ (1/τ)e^{−t/τ} v(x−t) dt  (SPEC 1.8)
   out_i = e^{−λ}·out_{i−1} + ∫₀^h (1/τ)e^{−t/τ} v(x_i − t) dt, 마디 사이 v 는 (v_i, v_{i−1}, v_{i−2}) 2차 보간 → 정확 적분 (λ = h/τ) */
function econv(v,h,tau){
  const n=v.length, out=new Float64Array(n), lam=h/tau, d=Math.exp(-lam);
  let I0,I1,I2;
  if(lam<1e-2){ const l2=lam*lam, l3=l2*lam, l4=l3*lam; I0=lam-l2/2+l3/6-l4/24; I1=lam/2-l2/3+l3/8-l4/30; I2=lam/3-l2/4+l3/10-l4/36; }
  else { I0=-Math.expm1(-lam); I1=(I0-lam*d)/lam; I2=2*I1/lam-d; }
  const a0=(I2-3*I1+2*I0)/2, a1=-I2+2*I1, a2=(I2-I1)/2;   // v_i, v_{i−1}, v_{i−2} 계수
  const b0=I0-I1, b1=I1;                                   // 첫 구간: 1차
  out[0]=0;   // 왼쪽 끝 앞은 0 (격자 여유가 충분)
  if(n>1) out[1]=d*out[0]+b0*v[1]+b1*v[0];
  for(let i=2;i<n;i++) out[i]=d*out[i-1]+a0*v[i]+a1*v[i-1]+a2*v[i-2];
  return out;
}
/* 상대좌표 x(배열)에서 넓이 1 프로파일. ctx={lo,hi}: 피팅 영역(상대좌표, DS 정규화·격자 범위) */
function shapeUnit(shape,x,p,ctx){
  const n=x.length, out=new Float64Array(n);
  if(shape==="GLS"){ const F=p.F,m=p.m; for(let i=0;i<n;i++) out[i]=glsVal(x[i],F,m); return out; }
  if(shape==="GLP"){ const F=p.F,m=p.m, N=glpNorm(F,m), a=4*LN2*(1-m)/(F*F), b=4*m/(F*F); for(let i=0;i<n;i++){ const t=x[i]*x[i]; out[i]=Math.exp(-a*t)/(1+b*t)/N; } return out; }
  if(shape==="V"){ for(let i=0;i<n;i++) out[i]=voigt(x[i],p.wG,p.wL); return out; }
  if(shape==="LA"&&!(p.wG>1e-6)){ const F=p.F, a=p.alpha, b=p.beta, N=laNorm(F,a,b); for(let i=0;i<n;i++){ const t=1+4*x[i]*x[i]/(F*F); out[i]=Math.pow(t,-(x[i]>=0?a:b))/N; } return out; }
  if(shape==="DS"&&!(p.wG>1e-6)){   // 해석: 값과 영역 적분(원시함수) 모두 닫힌 꼴
    const g=p.wL/2, al=p.alpha; let xmin=Infinity,xmax=-Infinity; for(let i=0;i<n;i++){ if(x[i]<xmin) xmin=x[i]; if(x[i]>xmax) xmax=x[i]; }
    const lo=ctx&&ctx.lo!=null?ctx.lo:xmin, hi=ctx&&ctx.hi!=null?ctx.hi:xmax; const N=dsPrim(hi,al,g)-dsPrim(lo,al,g);
    for(let i=0;i<n;i++) out[i]=dsBase(x[i],al,g)/(N>0?N:1); return out; }
  // 격자형
  return gridShape(shape,x,p,ctx);
}
/* 격자 간격: 0.04·2^-k 로 양자화(폭 매개변수가 조금 바뀌어도 격자 마디가 그대로라 목적함수가 매끄럽다) */
function gridH(shape,p){
  let w;
  if(shape==="DS") w=Math.min(p.wL,p.wG>1e-6?p.wG:Infinity);
  else if(shape==="LA") w=Math.min(p.F,p.wG>1e-6?p.wG:Infinity);
  else w=Math.min(p.F,p.tau);
  let h=0.04; while(h>w/16&&h>0.0005) h/=2;
  return h;
}
/* Catmull–Rom (C¹ 3차) 보간: 격자 v[k] (x = (k0+k)·h) → 임의 점 */
function crAt(v,k0,h,x){
  const M=v.length, t=x/h-k0; let k=Math.floor(t);
  if(k<1){ return k<0?v[0]:v[0]+(v[1]-v[0])*(t-k); }
  if(k>=M-2){ return k>=M-1?v[M-1]:v[M-2]+(v[M-1]-v[M-2])*(t-k); }
  const f=t-k, p0=v[k-1], p1=v[k], p2=v[k+1], p3=v[k+2], f2=f*f, f3=f2*f;
  return 0.5*((2*p1)+(-p0+p2)*f+(2*p0-5*p1+4*p2-p3)*f2+(-p0+3*p1-3*p2+p3)*f3);
}
function gridShape(shape,x,p,ctx){
  const n=x.length; const out=new Float64Array(n);
  if(!n) return out;
  let xmin=x[0],xmax=x[0]; for(let i=1;i<n;i++){ if(x[i]<xmin) xmin=x[i]; if(x[i]>xmax) xmax=x[i]; }
  const lo=ctx&&ctx.lo!=null?Math.min(ctx.lo,xmin):xmin, hi=ctx&&ctx.hi!=null?Math.max(ctx.hi,xmax):xmax;
  let h=gridH(shape,p);
  const sg=(p.wG||0)*FW2SIG;
  let big=shape==="DS"?Math.max(p.wL,p.wG||0):shape==="LA"?Math.max(p.F,p.wG||0):Math.max(p.F,p.tau);
  let P=8*big+5*sg; if(shape==="GLT") P=Math.max(P,12*p.tau+8*p.F);
  while((hi-lo+2*P)/h>400000) h*=2;
  // 마디는 상대좌표 x = k·h (k 정수) — 성분 중심 0 이 항상 마디
  const k0=Math.floor((lo-P)/h), k1=Math.ceil((hi+P)/h), M=k1-k0+1;
  const gx=new Float64Array(M); for(let i=0;i<M;i++) gx[i]=(k0+i)*h;
  let v=new Float64Array(M);
  if(shape==="DS"){ const g=p.wL/2, al=p.alpha; for(let i=0;i<M;i++) v[i]=dsBase(gx[i],al,g); if(p.wG>1e-6) v=gconv(v,h,p.wG); }
  else if(shape==="LA"){ const F=p.F, a=p.alpha, b=p.beta, N=laNorm(F,a,b); for(let i=0;i<M;i++){ const t=1+4*gx[i]*gx[i]/(F*F); v[i]=Math.pow(t,-(gx[i]>=0?a:b))/N; } if(p.wG>1e-6) v=gconv(v,h,p.wG); }   // α: 높은 BE 쪽 (SPEC 1.6 v1.1)
  else if(shape==="GLT"){ const F=p.F,m=p.m,T=p.T; const base=new Float64Array(M); for(let i=0;i<M;i++) base[i]=glsVal(gx[i],F,m); const tail=econv(base,h,p.tau); for(let i=0;i<M;i++) v[i]=(1-T)*base[i]+T*tail[i]; }
  // DS: 영역 안 적분으로 정규화
  let norm=1;
  if(shape==="DS"){
    const a=ctx&&ctx.lo!=null?ctx.lo:xmin, b=ctx&&ctx.hi!=null?ctx.hi:xmax;
    norm=gridIntegral(v,k0,h,a,b); if(!(norm>0)) norm=1;
  }
  for(let i=0;i<n;i++) out[i]=crAt(v,k0,h,x[i])/norm;
  return out;
}
/* [a,b] 적분: 안쪽 마디는 Simpson(가능하면)·사다리꼴, 양 끝 조각은 보간값으로 */
function gridIntegral(v,k0,h,a,b){
  const ka=Math.ceil(a/h-k0), kb=Math.floor(b/h-k0);
  if(kb<=ka) return crAt(v,k0,h,(a+b)/2)*(b-a);
  const xa=(k0+ka)*h, xb=(k0+kb)*h;
  let s=0; const m=kb-ka;
  if(m>=2){ const even=m%2===0; const e=even?kb:kb-1; for(let k=ka;k<e;k+=2) s+=h/3*(v[k]+4*v[k+1]+v[k+2]); if(!even) s+=0.5*h*(v[kb-1]+v[kb]); }
  else s=0.5*h*(v[ka]+v[kb]);
  // 끝 조각: 3점 Simpson (보간 중점)
  const seg=(u,w)=>{ if(w-u<=0) return 0; return (w-u)/6*(crAt(v,k0,h,u)+4*crAt(v,k0,h,(u+w)/2)+crAt(v,k0,h,w)); };
  return s+seg(a,xa)+seg(xb,b);
}
/* 실제 FWHM (수치): 넓이 1 프로파일의 반치폭. 단일 봉우리 가정 */
function trueFWHM(shape,p){
  const w=shape==="V"?(0.5346*p.wL+Math.sqrt(0.2166*p.wL*p.wL+p.wG*p.wG)):shape==="DS"?Math.max(p.wL,0.01)+(p.wG||0):shape==="LA"?p.F+(p.wG||0):p.F;
  const R=6*w+(shape==="GLT"?6*p.tau:0), N=2401; const x=new Float64Array(N); for(let i=0;i<N;i++) x[i]=-R+2*R*i/(N-1);
  const v=shapeUnit(shape,x,p,{lo:-R,hi:R});
  let im=0; for(let i=1;i<N;i++) if(v[i]>v[im]) im=i;
  const half=v[im]/2; let a=null,b=null;
  for(let i=im;i>0;i--) if(v[i-1]<half){ a=x[i-1]+(half-v[i-1])/(v[i]-v[i-1])*(x[i]-x[i-1]); break; }
  for(let i=im;i<N-1;i++) if(v[i+1]<half){ b=x[i]+(v[i]-half)/(v[i]-v[i+1])*(x[i+1]-x[i]); break; }
  return {fwhm:a!=null&&b!=null?b-a:null,xmax:x[im],ymax:v[im]};
}

/* ============================================================
   2) 성분 → 세기
   comp = {id, shape, be, area, F,m,wG,wL,alpha,beta,T,tau, doublet:{dE,r,q}|null}
   ============================================================ */
function shapeParams(c,scale){
  const q=scale||1, p={};
  SHAPES[c.shape].p.forEach(k=>p[k]=c[k]);
  if(q!==1){ ["F","wL"].forEach(k=>{ if(p[k]!=null) p[k]*=q; }); }   // 수명 폭만 (SPEC 1.9 v1.1) — wG 는 공통
  if(c.shape==="LA"){ p.alpha=Math.max(LA_ALPHA_LO,p.alpha); p.beta=Math.max(LA_ALPHA_LO,p.beta); }
  return p;
}
/* E: 절대 BE 배열, ctx={emin,emax} → Float64Array */
function compEval(c,E,ctx,parts){
  const n=E.length, out=new Float64Array(n);
  const one=(e0,A,p)=>{ const x=new Float64Array(n); for(let i=0;i<n;i++) x[i]=E[i]-e0; const u=shapeUnit(c.shape,x,p,ctx?{lo:ctx.emin-e0,hi:ctx.emax-e0}:null); const r=new Float64Array(n); for(let i=0;i<n;i++){ r[i]=A*u[i]; out[i]+=r[i]; } return r; };
  if(c.doublet){
    const d=c.doublet, A1=c.area/(1+d.r), A2=c.area*d.r/(1+d.r);
    const a=one(c.be,A1,shapeParams(c,1)), b=one(c.be+d.dE,A2,shapeParams(c,d.q));
    if(parts){ parts.push(a,b); }
  } else { const a=one(c.be,c.area,shapeParams(c,1)); if(parts) parts.push(a); }
  return out;
}

/* ============================================================
   3) Background (SPEC §2)
   ============================================================ */
function cumtrapz(E,f){ const n=E.length, Q=new Float64Array(n); for(let i=1;i<n;i++) Q[i]=Q[i-1]+0.5*(f[i]+f[i-1])*(E[i]-E[i-1]); return Q; }
function endpoints(E,y,n){
  const N=E.length; n=Math.max(1,Math.min(n||5,Math.floor(N/4)||1));
  let a=0,ea=0,b=0,eb=0; for(let i=0;i<n;i++){ a+=y[i]; ea+=E[i]; b+=y[N-1-i]; eb+=E[N-1-i]; }
  return {ylo:a/n,yhi:b/n,Elo:ea/n,Ehi:eb/n,n};
}
/* 국소 2차 다항식 평활 (SPEC §2 smart) */
function smoothQ(E,y,hw){
  const n=y.length, out=new Float64Array(n);
  for(let i=0;i<n;i++){
    const j0=Math.max(0,i-hw), j1=Math.min(n-1,i+hw), x0=E[i];
    let S0=0,S1=0,S2=0,S3=0,S4=0,T0=0,T1=0,T2=0;
    for(let j=j0;j<=j1;j++){ const t=E[j]-x0, t2=t*t; S0++; S1+=t; S2+=t2; S3+=t2*t; S4+=t2*t2; T0+=y[j]; T1+=y[j]*t; T2+=y[j]*t2; }
    // [S0 S1 S2; S1 S2 S3; S2 S3 S4] c = [T0 T1 T2] → c0
    const det=S0*(S2*S4-S3*S3)-S1*(S1*S4-S3*S2)+S2*(S1*S3-S2*S2);
    if(Math.abs(det)<1e-300||j1-j0<3){ out[i]=T0/S0; continue; }
    out[i]=(T0*(S2*S4-S3*S3)-S1*(T1*S4-S3*T2)+S2*(T1*S3-S2*T2))/det;
  }
  return out;
}
function smoothHW(E){ const st=E.length>1?(E[E.length-1]-E[0])/(E.length-1):0.05; return Math.max(2,Math.round(0.25/st)); }
function bgLinear(E,y,n){ const p=endpoints(E,y,n); const sl=(p.yhi-p.ylo)/((p.Ehi-p.Elo)||1); return Float64Array.from(E,e=>p.ylo+sl*(e-p.Elo)); }
function bgShirley(E,y,n,smart){
  const N=E.length, p=endpoints(E,y,n); const d=p.yhi-p.ylo;
  let B=new Float64Array(N).fill(p.ylo);
  if(Math.abs(d)<1e-300&&!smart) return B;
  const ys=smart?smoothQ(E,y,smoothHW(E)):null;
  const tol=1e-7*(Math.abs(d)+1e-12); let it=0, conv=false;
  for(;it<200;it++){
    const f=new Float64Array(N); for(let i=0;i<N;i++) f[i]=y[i]-B[i];
    const Q=cumtrapz(E,f), Qt=Q[N-1];
    const nb=new Float64Array(N);
    for(let i=0;i<N;i++) nb[i]=Math.abs(Qt)>1e-300?p.ylo+d*Q[i]/Qt:p.ylo;
    if(smart) for(let i=0;i<N;i++) if(nb[i]>ys[i]) nb[i]=ys[i];
    let mx=0; for(let i=0;i<N;i++) mx=Math.max(mx,Math.abs(nb[i]-B[i]));
    B=nb; if(mx<tol){ conv=true; break; }
  }
  B.iter=it+1; B.converged=conv;
  return B;
}
function bgTougaard(E,y,n,o){
  const N=E.length, p=endpoints(E,y,n); o=o||{};
  const C=o.C!=null?o.C:1643, D=o.D||0, U3=!!o.u3;
  const K=T=>T<=0?0:(U3?T/((C-T*T)*(C-T*T)+D*T*T):T/((C+T*T)*(C+T*T)));
  const I=new Float64Array(N);
  for(let i=1;i<N;i++){ let s=0; for(let j=0;j<i;j++){ const fa=K(E[i]-E[j])*(y[j]-p.ylo), fb=K(E[i]-E[j+1])*(y[j+1]-p.ylo); s+=0.5*(fa+fb)*(E[j+1]-E[j]); } I[i]=s; }
  // I(E_hi): 선형 보간
  let Ih=I[N-1]; for(let i=1;i<N;i++) if(E[i]>=p.Ehi){ const f=(p.Ehi-E[i-1])/(E[i]-E[i-1]); Ih=I[i-1]+f*(I[i]-I[i-1]); break; }
  let s=Math.abs(Ih)>1e-300?(p.yhi-p.ylo)/Ih:0; const neg=s<0; if(neg) s=0;   // 음의 손실 단면적 금지 (SPEC 2 v1.1)
  const B=Float64Array.from(I,v=>p.ylo+s*v); B.scale=s; B.negClamped=neg; return B;
}
function bgStatic(E,y,bg){
  const t=bg.type, n=bg.n||5;
  if(t==="linear") return bgLinear(E,y,n);
  if(t==="shirley") return bgShirley(E,y,n,false);
  if(t==="smart") return bgShirley(E,y,n,true);
  if(t==="tougaard") return bgTougaard(E,y,n,{C:bg.C,D:bg.D,u3:bg.u3});
  if(t==="none") return new Float64Array(E.length);
  return null;   // active
}
const BG_TYPES={linear:"Linear",shirley:"Shirley (반복)",smart:"Smart (Shirley·데이터 초과 금지)",tougaard:"Tougaard (U2/U3)",active:"Shirley — 피팅 안에서 함께",none:"없음"};

/* ============================================================
   4) 문제 구성 · 매개변수 벡터 · 제약/연결
   prob = {E, y, c, bg:{type,n,C,D,u3, c0,k,s, cons:{c0,k,s}}, comps:[...], cons:{compId:{param:{mode,lo,hi,ref,k,b}}}}
   ============================================================ */
function paramNames(c){ const a=["be","area"].concat(SHAPES[c.shape].p); if(c.doublet) a.push("dE","r","q"); return a; }
function getP(c,k){ return (k==="dE"||k==="r"||k==="q")?c.doublet[k]:c[k]; }
function setP(c,k,v){ if(k==="dE"||k==="r"||k==="q") c.doublet[k]=v; else c[k]=v; }
function defBounds(c,k){
  const d=PDEF[k]||{}; let lo=d.lo!=null?d.lo:-Infinity, hi=d.hi!=null?d.hi:Infinity;
  if(k==="be"){ lo=c.be-3; hi=c.be+3; }
  if(c.shape==="LA"&&(k==="alpha")){ lo=LA_ALPHA_LO; hi=20; }
  return [lo,hi];
}
function clone(o){ return JSON.parse(JSON.stringify(o)); }
/* 자유 매개변수 목록과 연결 해석 */
function buildMap(prob){
  const comps=prob.comps, byId={}; comps.forEach((c,i)=>byId[c.id]=i);
  const free=[], links=[], fixed=[];
  comps.forEach((c,ci)=>{
    paramNames(c).forEach(k=>{
      const cs=((prob.cons||{})[c.id]||{})[k]||{mode:"free"};
      if(cs.mode==="fixed") { fixed.push({ci,k}); return; }
      if(cs.mode==="link"&&cs.ref!=null&&byId[cs.ref]!=null&&byId[cs.ref]!==ci){ links.push({ci,k,ri:byId[cs.ref],rk:cs.rk||k,mul:cs.k!=null?+cs.k:1,add:cs.b!=null?+cs.b:0}); return; }
      let [lo,hi]=defBounds(c,k); if(cs.lo!=null&&isFinite(cs.lo)) lo=+cs.lo; if(cs.hi!=null&&isFinite(cs.hi)) hi=+cs.hi;
      if(k==="area") lo=Math.max(0,lo);
      free.push({ci,k,lo,hi,id:c.id});
    });
  });
  // 배경 매개변수 (active)
  if(prob.bg&&prob.bg.type==="active"){
    const bc=prob.bg.cons||{};
    ["c0","k","s"].forEach(k=>{ const cs=bc[k]||{mode:k==="s"?"fixed":"free"}; if(cs.mode==="fixed") return; const d=PDEF[k];
      // k 상한: 계단 높이 k·∫P 가 데이터 세기 범위의 1.5 배를 넘지 못하게(SPEC 2 v1.1) — 넓이→0·k→∞ 퇴화 방향을 막는다
      const hi=cs.hi!=null?cs.hi:(k==="k"&&prob.bg._khi>0?prob.bg._khi:d.hi);
      free.push({ci:-1,k,lo:cs.lo!=null?cs.lo:d.lo,hi,id:"bg"}); });
  }
  // 연결 사슬: 연결 대상이 또 연결이면 끝까지 따라간다(순환 금지)
  const linkOf=(ci,k)=>links.find(l=>l.ci===ci&&l.k===k);
  links.forEach(l=>{ let g=0, cur=l; const seen=new Set([l.ci+"|"+l.k]); while(true){ const nx=linkOf(cur.ri,cur.rk); if(!nx) break; const key=nx.ci+"|"+nx.k; if(seen.has(key)||g++>20) throw new Error("연결이 순환합니다: "+comps[l.ci].name+" "+l.k); seen.add(key); cur=nx; } });
  return {free,links,fixed};
}
function applyLinks(prob,map){
  // 위상 순서대로 (사슬 길이만큼 반복)
  for(let pass=0;pass<8;pass++){ map.links.forEach(l=>{ const ref=prob.comps[l.ri]; const v=getP(ref,l.rk)*l.mul+l.add; setP(prob.comps[l.ci],l.k,v); }); }
}
function getVec(prob,map){ return map.free.map(f=>f.ci<0?prob.bg[f.k]:getP(prob.comps[f.ci],f.k)); }
function setVec(prob,map,v){ map.free.forEach((f,i)=>{ if(f.ci<0) prob.bg[f.k]=v[i]; else setP(prob.comps[f.ci],f.k,v[i]); }); applyLinks(prob,map); }
/* 무제약 변환 (MINUIT) */
function toU(p,lo,hi){ const L=isFinite(lo), H=isFinite(hi);
  if(L&&H){ const t=Math.max(-1,Math.min(1,2*(p-lo)/(hi-lo)-1)); return Math.asin(t); }
  if(L){ const t=p-lo+1; return Math.sqrt(Math.max(0,t*t-1)); }
  if(H){ const t=hi-p+1; return Math.sqrt(Math.max(0,t*t-1)); }
  return p; }
function fromU(u,lo,hi){ const L=isFinite(lo), H=isFinite(hi);
  if(L&&H) return lo+(hi-lo)*(Math.sin(u)+1)/2;
  if(L) return lo-1+Math.sqrt(u*u+1);
  if(H) return hi+1-Math.sqrt(u*u+1);
  return u; }
function dPdU(u,lo,hi){ const L=isFinite(lo), H=isFinite(hi);
  if(L&&H) return (hi-lo)*Math.cos(u)/2;
  if(L) return u/Math.sqrt(u*u+1);
  if(H) return -u/Math.sqrt(u*u+1);
  return 1; }

/* ============================================================
   5) 모형 평가 (성분 캐시: 한 매개변수를 흔들 때 영향받는 성분만 다시 계산)
   ============================================================ */
function makeEvaluator(prob,map){
  const E=prob.E, n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const Bst=prob.bg.type==="active"?null:(prob._Bst||(prob._Bst=bgStatic(E,prob.y,prob.bg)));
  const nc=prob.comps.length;
  // 자유 매개변수 i 가 바꾸는 성분 목록
  const dep=map.free.map(f=>{ if(f.ci<0) return []; const s=new Set([f.ci]); let grow=true; while(grow){ grow=false; map.links.forEach(l=>{ if(s.has(l.ri)&&!s.has(l.ci)){ s.add(l.ci); grow=true; } }); } return [...s]; });
  let cache=null;
  function compAll(){ cache=prob.comps.map(c=>compEval(c,E,ctx)); }
  function total(cs){
    const P=new Float64Array(n); cs.forEach(v=>{ for(let i=0;i<n;i++) P[i]+=v[i]; });
    let B;
    if(Bst) B=Bst;
    else { const Q=cumtrapz(E,P); const c0=prob.bg.c0||0, k=prob.bg.k||0, s=prob.bg.s||0; B=new Float64Array(n); for(let i=0;i<n;i++) B[i]=c0+k*Q[i]+s*(E[i]-E[0]); }
    const mu=new Float64Array(n); for(let i=0;i<n;i++) mu[i]=B[i]+P[i];
    return {mu,B,P};
  }
  return {
    ctx, Bst, dep,
    full(){ compAll(); return total(cache); },
    /* 자유 매개변수 j 만 바꾼 상태의 μ (캐시 사용, 캐시는 바꾸지 않음) */
    partial(j){ const cs=cache.slice(); dep[j].forEach(ci=>{ cs[ci]=compEval(prob.comps[ci],E,ctx); }); return total(cs).mu; },
    comps(){ return cache; },
  };
}

/* ============================================================
   6) 선형대수 · NNLS · 난수
   ============================================================ */
function cholSolve(A,b,n){
  const L=new Float64Array(n*n);
  for(let i=0;i<n;i++){ for(let j=0;j<=i;j++){ let s=A[i*n+j]; for(let k=0;k<j;k++) s-=L[i*n+k]*L[j*n+k]; if(i===j){ if(!(s>0)) return null; L[i*n+i]=Math.sqrt(s); } else L[i*n+j]=s/L[j*n+j]; } }
  const z=new Float64Array(n); for(let i=0;i<n;i++){ let s=b[i]; for(let k=0;k<i;k++) s-=L[i*n+k]*z[k]; z[i]=s/L[i*n+i]; }
  const x=new Float64Array(n); for(let i=n-1;i>=0;i--){ let s=z[i]; for(let k=i+1;k<n;k++) s-=L[k*n+i]*x[k]; x[i]=s/L[i*n+i]; }
  return x;
}
function invSym(A,n){   // 대칭 양정치 역행렬 (Cholesky). 실패 시 null
  const out=new Float64Array(n*n);
  for(let j=0;j<n;j++){ const e=new Float64Array(n); e[j]=1; const x=cholSolve(A,e,n); if(!x) return null; for(let i=0;i<n;i++) out[i*n+j]=x[i]; }
  return out;
}
/* Lawson–Hanson NNLS: min ||A x − b||, x ≥ 0. A: 열 배열(cols[j] = Float64Array(m)) */
function nnls(cols,b){
  const m=b.length, n=cols.length; const x=new Float64Array(n); const P=new Array(n).fill(false);
  const dot=(u,v)=>{ let s=0; for(let i=0;i<m;i++) s+=u[i]*v[i]; return s; };
  const resid=()=>{ const r=Float64Array.from(b); for(let j=0;j<n;j++) if(x[j]) for(let i=0;i<m;i++) r[i]-=cols[j][i]*x[j]; return r; };
  const lsP=idx=>{ const k=idx.length; const A=new Float64Array(k*k), g=new Float64Array(k); for(let a=0;a<k;a++){ g[a]=dot(cols[idx[a]],b); for(let c=0;c<=a;c++){ const v=dot(cols[idx[a]],cols[idx[c]]); A[a*k+c]=v; A[c*k+a]=v; } }
    let s=null, jit=0; while(!s&&jit<6){ if(jit){ for(let a=0;a<k;a++) A[a*k+a]*=1+1e-10*Math.pow(10,jit); } s=cholSolve(A,g,k); jit++; } return s||new Float64Array(k); };
  for(let it=0;it<3*n+10;it++){
    const r=resid(); const w=cols.map((c,j)=>P[j]?-Infinity:dot(c,r));
    let jm=-1, wm=1e-12*Math.max(1,Math.sqrt(dot(b,b))); for(let j=0;j<n;j++) if(!P[j]&&w[j]>wm){ wm=w[j]; jm=j; }
    if(jm<0) break; P[jm]=true;
    for(let inner=0;inner<3*n+10;inner++){
      const idx=[]; for(let j=0;j<n;j++) if(P[j]) idx.push(j);
      const z=lsP(idx);
      if(z.every(v=>v>0)){ idx.forEach((j,a)=>x[j]=z[a]); break; }
      let al=Infinity; idx.forEach((j,a)=>{ if(z[a]<=0){ const t=x[j]/(x[j]-z[a]); if(t<al) al=t; } });
      if(!isFinite(al)) al=0;
      idx.forEach((j,a)=>{ x[j]=x[j]+al*(z[a]-x[j]); if(x[j]<=1e-14*Math.abs(z[a]||1)){ x[j]=0; P[j]=false; } });
    }
  }
  return x;
}
function rng(seed){ let a=(seed>>>0)||1; const u=()=>{ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; let sp=null; const g=()=>{ if(sp!=null){ const v=sp; sp=null; return v; } let u1=0; while(u1<1e-300) u1=u(); const u2=u(); const r=Math.sqrt(-2*Math.log(u1)); sp=r*Math.sin(2*PI*u2); return r*Math.cos(2*PI*u2); }; return {u,g}; }

/* ============================================================
   7) 피팅 (SPEC §3)
   ============================================================ */
/* SPEC 3.1: v = max(ŷ, v_min)/c, v_min = max(1/c, 1e-6·max ŷ) */
function varianceFrom(yhat,c){ const cc=c>0?c:1; let m=0; for(const t of yhat) if(t>m) m=t; const vmin=Math.max(1/cc,1e-6*m); return Float64Array.from(yhat,v=>Math.max(v,vmin)/cc); }
/* 넓이 초기값: 비선형 매개변수 고정, 가중 NNLS */
function initAreas(prob,map,v){
  const E=prob.E, n=E.length, ctx={emin:E[0],emax:E[n-1]};
  // active 의 넓이 초기화는 직선 배경으로 — 정적 Shirley 는 기울어진 배경에서 수렴하지 않을 수 있다(고정점이 없음)
  const Bst=prob.bg.type==="active"?bgLinear(E,prob.y,prob.bg.n||5):(prob._Bst||(prob._Bst=bgStatic(E,prob.y,prob.bg)));
  // 넓이 묶음: 자유 넓이마다 한 열(그 넓이에 ×k 로 연결된 성분 포함)
  const groups=[]; const freeArea=map.free.map((f,i)=>f.k==="area"&&f.ci>=0?i:-1).filter(i=>i>=0);
  if(!freeArea.length) return;
  const sw=Float64Array.from(v,x=>1/Math.sqrt(x));
  const cols=freeArea.map(fi=>{ const f=map.free[fi]; const col=new Float64Array(n);
    const add=(ci,mul)=>{ const c=prob.comps[ci]; const save=c.area; c.area=mul; const u=compEval(c,E,ctx); c.area=save; for(let i=0;i<n;i++) col[i]+=u[i]; map.links.forEach(l=>{ if(l.ri===ci&&l.k==="area"&&l.rk==="area") add(l.ci,mul*l.mul); }); };
    add(f.ci,1); for(let i=0;i<n;i++) col[i]*=sw[i]; return col; });
  // 고정 넓이 성분은 빼고 맞춘다
  const rhs=new Float64Array(n); for(let i=0;i<n;i++) rhs[i]=prob.y[i]-Bst[i];
  prob.comps.forEach((c,ci)=>{ if(freeArea.some(fi=>map.free[fi].ci===ci)) return; if(map.links.some(l=>l.ci===ci&&l.k==="area")) return; const u=compEval(c,E,ctx); for(let i=0;i<n;i++) rhs[i]-=u[i]; });
  for(let i=0;i<n;i++) rhs[i]*=sw[i];
  const a=nnls(cols,rhs);
  freeArea.forEach((fi,j)=>{ const f=map.free[fi]; const val=Math.max(a[j],1e-9*(Math.abs(prob.y[0])+1)); setP(prob.comps[f.ci],"area",Math.min(val,f.hi)); });
  applyLinks(prob,map);
  if(prob.bg.type==="active"){
    const p=endpoints(E,prob.y,prob.bg.n||5); const P=new Float64Array(n); prob.comps.forEach(c=>{ const u=compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=u[i]; });
    const Q=cumtrapz(E,P); if(prob.bg.c0==null||prob.bg._auto!==false) prob.bg.c0=p.ylo; if(prob.bg.k==null||prob.bg._auto!==false) prob.bg.k=Q[n-1]>0?Math.max(0,(p.yhi-p.ylo)/Q[n-1]):0; if(prob.bg.s==null) prob.bg.s=0;
    let ymn=Infinity, ymx=-Infinity; for(const t of prob.y){ if(t<ymn) ymn=t; if(t>ymx) ymx=t; }
    if(Q[n-1]>0) prob.bg._khi=Math.max(1.5*(ymx-ymn)/Q[n-1],2*prob.bg.k);
  }
}
/* LM 한 번. prob 은 제자리에서 바뀐다. v: 분산 배열 */
function lmFit(prob,map,v,opt){
  opt=opt||{}; const maxit=opt.maxit||500; const E=prob.E, y=prob.y, n=E.length;
  const np=map.free.length; const ev=makeEvaluator(prob,map);
  const sw=Float64Array.from(v,x=>1/Math.sqrt(x));
  const lo=map.free.map(f=>f.lo), hi=map.free.map(f=>f.hi);
  let p=getVec(prob,map).map((x,i)=>Math.min(hi[i],Math.max(lo[i],x)));
  let u=p.map((x,i)=>toU(x,lo[i],hi[i]));
  const chi=mu=>{ let s=0; for(let i=0;i<n;i++){ const r=(y[i]-mu[i])*sw[i]; s+=r*r; } return s; };
  setVec(prob,map,p); let cur=ev.full(); let X=chi(cur.mu);
  if(!np) return {chi2:X,iter:0,conv:true,mu:cur.mu,B:cur.B};
  let lam=1e-3, nu=2, calm=0, it=0, conv=false;
  const scale=map.free.map(f=>f.k==="be"||f.k==="dE"?1e-4:f.k==="area"?0:1e-5);
  for(;it<maxit;it++){
    // Jacobian (u 공간): ∂r/∂u = −sw·∂μ/∂p · dp/du
    const J=new Array(np);
    for(let j=0;j<np;j++){
      const pj=p[j]; let h=Math.max(Math.abs(pj)*1e-6,scale[j]||1e-7); if(map.free[j].k==="area") h=Math.max(Math.abs(pj)*1e-6,1e-9);
      // 경계 안쪽으로 차분
      let a=pj-h, b=pj+h; if(b>hi[j]) { b=pj; } if(a<lo[j]) { a=pj; } if(b===a){ J[j]=new Float64Array(n); continue; }
      const setj=val=>{ const q=p.slice(); q[j]=val; setVec(prob,map,q); };
      setj(b); const mb=ev.partial(j); setj(a); const ma=ev.partial(j); setVec(prob,map,p);
      const d=dPdU(u[j],lo[j],hi[j])/(b-a); const col=new Float64Array(n); for(let i=0;i<n;i++) col[i]=-(mb[i]-ma[i])*d*sw[i]; J[j]=col;
    }
    const r=new Float64Array(n); for(let i=0;i<n;i++) r[i]=(y[i]-cur.mu[i])*sw[i];
    const A=new Float64Array(np*np), g=new Float64Array(np);
    for(let a=0;a<np;a++){ const Ja=J[a]; let s=0; for(let i=0;i<n;i++) s+=Ja[i]*r[i]; g[a]=-s; for(let b=0;b<=a;b++){ const Jb=J[b]; let t=0; for(let i=0;i<n;i++) t+=Ja[i]*Jb[i]; A[a*np+b]=t; A[b*np+a]=t; } }
    const D=new Float64Array(np); for(let a=0;a<np;a++) D[a]=Math.max(A[a*np+a],1e-30);
    let accepted=false, tries=0;
    while(!accepted&&tries++<30){
      const M=A.slice(); for(let a=0;a<np;a++) M[a*np+a]+=lam*D[a];
      let dlt=cholSolve(M,g,np);
      if(!dlt){ lam*=nu; nu*=2; continue; }
      const un=u.map((x,i)=>x+dlt[i]); const pn=un.map((x,i)=>fromU(x,lo[i],hi[i]));
      setVec(prob,map,pn); const tr=ev.full(); const Xn=chi(tr.mu);
      let pred=0; for(let a=0;a<np;a++) pred+=dlt[a]*(lam*D[a]*dlt[a]+g[a]);
      const rho=pred>0?(X-Xn)/pred:-1;
      if(rho>0&&isFinite(Xn)){
        const rel=(X-Xn)/Math.max(X,1e-300);
        u=un; p=pn; cur=tr; X=Xn; accepted=true;
        lam*=Math.max(1/3,1-Math.pow(2*rho-1,3)); nu=2;
        if(rel<1e-10) calm++; else calm=0;
      } else { setVec(prob,map,p); lam*=nu; nu*=2; if(lam>1e16) break; }
    }
    if(!accepted){ setVec(prob,map,p); cur=ev.full(); conv=true; break; }
    if(calm>=3){ conv=true; break; }
  }
  setVec(prob,map,p); cur=ev.full();
  return {chi2:X,iter:it,conv,mu:cur.mu,B:cur.B,comps:ev.comps()};
}
/* 한 번의 완전한 피팅: (넓이 초기화) → LM(평활 분산) → LM(모형 분산) */
function fitOnce(prob,opt){
  opt=opt||{};
  let map=buildMap(prob); applyLinks(prob,map);
  const E=prob.E, y=prob.y; let ymax=0; for(const t of y) ymax=Math.max(ymax,Math.abs(t));
  const ys=smoothQ(E,y,smoothHW(E));
  let v=varianceFrom(ys,prob.c);
  if(opt.initAreas!==false){ initAreas(prob,map,v); if(prob.bg.type==="active"){ map=buildMap(prob); applyLinks(prob,map); } }   // k 상한 반영
  let r=lmFit(prob,map,v,opt); const chi2p1=r.chi2;
  if(opt.irls!==false){ v=varianceFrom(r.mu,prob.c); r=lmFit(prob,map,v,opt); }
  r.v=v; r.map=map; r.chi2p1=chi2p1;   // 1차 χ²(데이터로 정한 분산) — 시작점끼리 비교하는 공통 목적함수
  return r;
}
/* 다중 시작: 비선형 free 매개변수를 Latin hypercube 로 */
function startPoints(prob,N,seed){
  const map=buildMap(prob); const R=rng(seed||1);
  const idx=map.free.map((f,i)=>i).filter(i=>map.free[i].k!=="area"&&map.free[i].ci>=0);
  const pts=[];
  const perm=idx.map(()=>{ const a=[...Array(N).keys()]; for(let i=N-1;i>0;i--){ const j=Math.floor(R.u()*(i+1)); [a[i],a[j]]=[a[j],a[i]]; } return a; });
  for(let s=0;s<N;s++){
    const o={};
    idx.forEach((fi,q)=>{ const f=map.free[fi]; let lo=f.lo, hi=f.hi; const c=prob.comps[f.ci]; const cur=getP(c,f.k);
      if(!isFinite(lo)||!isFinite(hi)){ const w=Math.max(Math.abs(cur)*0.5,0.5); lo=isFinite(lo)?lo:cur-w; hi=isFinite(hi)?hi:cur+w; }
      const t=(perm[q][s]+R.u())/N; o[f.ci+"|"+f.k]=lo+(hi-lo)*t; });
    pts.push(s===0?null:o);    // 첫 시작은 현재 값 그대로
  }
  return pts;
}
function applyStart(prob,o){ if(!o) return; Object.entries(o).forEach(([key,v])=>{ const [ci,k]=key.split("|"); setP(prob.comps[+ci],k,v); }); }
function stripProb(prob){ const o=Object.assign({},prob); delete o._Bst; return o; }
function fitMulti(prob,opt){
  opt=opt||{}; const N=Math.max(1,opt.starts||1); const pts=startPoints(prob,N,opt.seed||1);
  let best=null; const all=[];
  // 시작점 비교는 1차 χ²(모든 시작점에 같은 가중)로 — 2차 분산 v=μ₁ 은 시작점마다 달라 2차 χ² 는 서로 비교할 수 없다 (SPEC 3.3)
  pts.forEach((pt,s)=>{ const P=clone(stripProb(prob)); P.E=prob.E; P.y=prob.y; applyStart(P,pt); try{ const r=fitOnce(P,opt); all.push(r.chi2p1); if(!best||r.chi2p1<best.r.chi2p1) best={P,r}; }catch(e){ all.push(Infinity); } if(opt.progress) opt.progress(s+1,N,best&&best.r.chi2p1); });
  if(!best) throw new Error("모든 시작점에서 피팅 실패");
  const near=all.filter(x=>x-best.r.chi2p1<1).length;
  best.r.multi={N,chi2:all,near,frac:near/N};
  return best;
}

/* ============================================================
   8) 불확도 · 적합도 (SPEC §3.4–3.5)
   ============================================================ */
function jacNatural(prob,map,v){
  const E=prob.E, n=E.length, ev=makeEvaluator(prob,map); ev.full();
  const p=getVec(prob,map), np=p.length, sw=Float64Array.from(v,x=>1/Math.sqrt(x));
  const J=[]; const atB=[];
  map.free.forEach((f,j)=>{
    const rng_=f.hi-f.lo; const near=isFinite(rng_)?Math.min(p[j]-f.lo,f.hi-p[j])<1e-6*rng_:(isFinite(f.lo)&&p[j]-f.lo<1e-9*(Math.abs(p[j])+1));
    atB.push(!!near);
    let h=Math.max(Math.abs(p[j])*1e-5,f.k==="be"||f.k==="dE"?1e-4:1e-6); if(f.k==="area") h=Math.max(Math.abs(p[j])*1e-5,1e-9);
    let a=p[j]-h, b=p[j]+h; if(b>f.hi) b=p[j]; if(a<f.lo) a=p[j];
    const q=p.slice(); q[j]=b; setVec(prob,map,q); const mb=ev.partial(j); q[j]=a; setVec(prob,map,q); const ma=ev.partial(j); setVec(prob,map,p);
    const col=new Float64Array(n); for(let i=0;i<n;i++) col[i]=(mb[i]-ma[i])/(b-a)*sw[i]; J.push(col);
  });
  return {J,atB,p};
}
function covariance(prob,res){
  const map=res.map||buildMap(prob); const {J,atB,p}=jacNatural(prob,map,res.v);
  const n=prob.E.length, np=p.length; const act=[]; for(let j=0;j<np;j++) if(!atB[j]) act.push(j);
  const k=act.length; const A=new Float64Array(k*k);
  for(let a=0;a<k;a++) for(let b=0;b<=a;b++){ let s=0; const Ja=J[act[a]], Jb=J[act[b]]; for(let i=0;i<n;i++) s+=Ja[i]*Jb[i]; A[a*k+b]=s; A[b*k+a]=s; }
  // 척도 맞춘 역행렬(조건수 완화)
  const d=new Float64Array(k); for(let a=0;a<k;a++) d[a]=Math.sqrt(Math.max(A[a*k+a],1e-300));
  const As=new Float64Array(k*k); for(let a=0;a<k;a++) for(let b=0;b<k;b++) As[a*k+b]=A[a*k+b]/(d[a]*d[b]);
  let inv=invSym(As,k); let singular=false;
  if(!inv){ singular=true; for(let a=0;a<k;a++) As[a*k+a]*=1+1e-8; inv=invSym(As,k); }
  const nfree=np; const dof=Math.max(1,n-k);
  const phi=res.chi2/dof;
  const C=new Float64Array(np*np).fill(NaN);
  if(inv) for(let a=0;a<k;a++) for(let b=0;b<k;b++) C[act[a]*np+act[b]]=inv[a*k+b]/(d[a]*d[b])*phi;
  const sig=p.map((x,j)=>atB[j]?null:Math.sqrt(Math.max(0,C[j*np+j])));
  const corr=new Float64Array(np*np).fill(NaN); for(let a=0;a<np;a++) for(let b=0;b<np;b++){ const s=sig[a]*sig[b]; if(s>0) corr[a*np+b]=C[a*np+b]/s; }
  return {C,sig,corr,atB,phi,dof,p,np,singular,names:map.free.map(f=>({id:f.id,k:f.k,ci:f.ci}))};
}
function fitStats(prob,res){
  const y=prob.y, mu=res.mu, v=res.v, n=y.length; const p=(res.map||buildMap(prob)).free.length;
  let chi2=0, ss=0, ybar=0; for(let i=0;i<n;i++) ybar+=y[i]; ybar/=n;
  const r=new Float64Array(n); for(let i=0;i<n;i++){ r[i]=(y[i]-mu[i])/Math.sqrt(v[i]); chi2+=r[i]*r[i]; ss+=(y[i]-ybar)*(y[i]-ybar); }
  let sr=0; for(let i=0;i<n;i++) sr+=(y[i]-mu[i])*(y[i]-mu[i]);
  let a1=0, a0=0, dw=0; for(let i=0;i<n;i++){ a0+=r[i]*r[i]; if(i) { a1+=r[i]*r[i-1]; dw+=(r[i]-r[i-1])*(r[i]-r[i-1]); } }
  const rho1=a0>0?a1/a0:0;
  // runs 시험
  let runs=1, np_=0, nn=0; for(let i=0;i<n;i++){ if(r[i]>=0) np_++; else nn++; if(i&&((r[i]>=0)!==(r[i-1]>=0))) runs++; }
  const mR=1+2*np_*nn/n, vR=2*np_*nn*(2*np_*nn-n)/(n*n*(n-1)); const zR=vR>0?(runs-mR)/Math.sqrt(vR):0;
  const dof=Math.max(1,n-p);
  // 모형 비교용 χ²_d: 데이터만으로 정한 분산(평활 데이터) — 모형마다 같은 가중 (SPEC 3.5 v1.1)
  let ymax=0; for(const t of y) ymax=Math.max(ymax,Math.abs(t));
  const vd=varianceFrom(smoothQ(prob.E,y,smoothHW(prob.E)),prob.c); let chi2d=0; for(let i=0;i<n;i++){ const d=y[i]-mu[i]; chi2d+=d*d/vd[i]; }
  const AICc=n*Math.log(chi2d/n)+2*p+(n-p-1>0?2*p*(p+1)/(n-p-1):Infinity), BIC=n*Math.log(chi2d/n)+p*Math.log(n);
  return {N:n,p,dof,chi2,chi2d,chi2nu:chi2/dof,R2:ss>0?1-sr/ss:null,rmsNorm:Math.sqrt(chi2/n),DW:a0>0?dw/a0:null,rho1,ar1Inflate:rho1>0?Math.sqrt((1+rho1)/(1-Math.min(rho1,0.95))):1,runs,runsZ:zR,AICc,BIC,resid:r};
}

/* 성분별 결과(보이는 값) + 델타 방법 불확도 */
function compSummary(prob,res,cov){
  const map=res.map||buildMap(prob); const E=prob.E, n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const p0=getVec(prob,map);
  const quantities=(P)=>{ // P: 현재 상태의 prob
    const out=[]; let tot=0;
    P.comps.forEach(c=>{ tot+=c.area; });
    P.comps.forEach(c=>{
      const sp=shapeParams(c,1); const tf=trueFWHM(c.shape,sp);
      const rec={be:c.be,area:c.area,areaPct:tot>0?100*c.area/tot:null,fwhm:tf.fwhm,height:(c.doublet?c.area/(1+c.doublet.r):c.area)*tf.ymax,beMax:c.be+tf.xmax};
      if(c.doublet){ rec.area1=c.area/(1+c.doublet.r); rec.area2=c.area*c.doublet.r/(1+c.doublet.r); rec.be2=c.be+c.doublet.dE; }
      out.push(rec);
    });
    return out;
  };
  const base=quantities(prob);
  // 영역 안 넓이(배경 위, 데이터 격자 사다리꼴)
  const comps=prob.comps.map(c=>compEval(c,E,ctx));
  base.forEach((b,i)=>{ const Q=cumtrapz(E,comps[i]); b.areaIn=Q[n-1]; });
  if(!cov) return base;
  // 수치 기울기
  const keys=["be","area","areaPct","fwhm","height","area1","area2","be2","beMax"];
  const np=p0.length; const grads=base.map(()=>({}));
  for(let j=0;j<np;j++){
    if(cov.atB[j]||!(cov.sig[j]>0)) continue;
    const f=map.free[j]; const h=Math.max(cov.sig[j]*1e-3,1e-9*(Math.abs(p0[j])+1));
    const q=p0.slice(); q[j]=p0[j]+h; setVec(prob,map,q); const up=quantities(prob); q[j]=p0[j]-h; setVec(prob,map,q); const dn=quantities(prob); setVec(prob,map,p0);
    base.forEach((b,i)=>keys.forEach(k=>{ if(b[k]==null||up[i][k]==null||dn[i][k]==null) return; (grads[i][k]=grads[i][k]||new Float64Array(np))[j]=(up[i][k]-dn[i][k])/(2*h); }));
  }
  base.forEach((b,i)=>{ b.sig={}; keys.forEach(k=>{ const g=grads[i][k]; if(!g){ if(b[k]!=null) b.sig[k]=0; return; } let s=0; for(let a=0;a<np;a++){ if(!g[a]) continue; for(let c=0;c<np;c++){ if(!g[c]) continue; const v=cov.C[a*np+c]; if(isFinite(v)) s+=g[a]*v*g[c]; } } b.sig[k]=Math.sqrt(Math.max(0,s)); }); });
  return base;
}

/* Monte Carlo (SPEC §3.4): 최적 μ + AR(1) 잡음 → 재피팅 */
function monteCarlo(prob,res,M,opt){
  opt=opt||{}; const R=rng(opt.seed||12345); const st=fitStats(prob,res);
  const rho=Math.max(0,Math.min(0.95,st.rho1)); const phi=st.chi2nu; const n=prob.E.length;
  const base=clone(stripProb(prob)); base.E=prob.E;
  const samples=[]; let fails=0;
  for(let m=0;m<M;m++){
    const y=new Float64Array(n); let e=0;
    for(let i=0;i<n;i++){ const z=R.g(); e=i?rho*e+Math.sqrt(1-rho*rho)*z:z; y[i]=res.mu[i]+e*Math.sqrt(phi*res.v[i]); }
    const P=clone(base); P.E=prob.E; P.y=y;
    try{ const r=fitOnce(P,{initAreas:false,irls:opt.irls!==false,maxit:opt.maxit||200}); samples.push(compSummary(P,r,null).map(q=>({be:q.be,area:q.area,areaPct:q.areaPct,fwhm:q.fwhm,height:q.height}))); }
    catch(e){ fails++; }
    if(opt.progress) opt.progress(m+1,M);
  }
  return {samples,fails,rho,phi};
}
function mcSummarize(samples){
  if(!samples.length) return null; const nc=samples[0].length; const keys=Object.keys(samples[0][0]);
  const pct=(a,q)=>{ const s=a.slice().sort((x,y)=>x-y); const t=(s.length-1)*q, k=Math.floor(t); return s[k]+(s[Math.min(k+1,s.length-1)]-s[k])*(t-k); };
  return [...Array(nc).keys()].map(i=>{ const o={}; keys.forEach(k=>{ const a=samples.map(s=>s[i][k]).filter(v=>v!=null&&isFinite(v)); if(!a.length) return; const m=a.reduce((x,y)=>x+y,0)/a.length; const sd=Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/Math.max(1,a.length-1)); o[k]={mean:m,sd,lo:pct(a,0.025),hi:pct(a,0.975),n:a.length}; }); return o; });
}

/* ============================================================
   9) 잡음 추정 (SPEC §6)
   ============================================================ */
function noiseEstimate(E,y){
  const n=y.length; if(n<8) return null; const ys=smoothQ(E,y,smoothHW(E));
  const q=[]; for(let i=1;i<n-1;i++){ const d=(y[i]-(y[i-1]+y[i+1])/2)/Math.sqrt(1.5); if(ys[i]>0) q.push(d*d/ys[i]); }
  q.sort((a,b)=>a-b); const med=q[Math.floor(q.length/2)];
  return {kappa:med/0.4549};
}

/* ============================================================
   10) 정량 (SPEC §4)
   ============================================================ */
function tfValue(tf,KE,PE){
  if(!tf||tf.model==="none") return 1;
  if(tf.model==="power"){ const pe=PE||tf.pe||50; return Math.pow(pe/KE,tf.x!=null?tf.x:1); }
  if(tf.model==="thermo"){ const pe=PE||tf.pe; if(!(pe>0)) return NaN; const t=Math.log(KE/pe); let s=0; for(let k=tf.coef.length-1;k>=0;k--) s=s*t+tf.coef[k]; return Math.exp(s); }
  if(tf.model==="table"){ const T=tf.table; if(!T||T.length<2) return NaN; if(KE<=T[0][0]) return T[0][1]; if(KE>=T[T.length-1][0]) return T[T.length-1][1];   // 표 밖은 끝값 (SPEC 4)
    const lk=Math.log(KE); let i=1; while(i<T.length-1&&Math.log(T[i][0])<lk) i++; const a=T[i-1], b=T[i]; const f=(lk-Math.log(a[0]))/(Math.log(b[0])-Math.log(a[0])); return Math.exp(Math.log(a[1])+f*(Math.log(b[1])-Math.log(a[1]))); }
  return 1;
}
function tpp2m(E,mat){
  const rho=mat.rho, Nv=mat.Nv, M=mat.M, Eg=mat.Eg||0;
  const U=Nv*rho/M, Ep=28.8*Math.sqrt(U);
  const beta=-0.10+0.944/Math.sqrt(Ep*Ep+Eg*Eg)+0.069*Math.pow(rho,0.1), gamma=0.191/Math.sqrt(rho), C=1.97-0.91*U, D=53.4-20.8*U;
  return E/(Ep*Ep*(beta*Math.log(gamma*E)-C/E+D/(E*E)));   // Å
}
function imfpValue(o,KE){
  if(!o||o.model==="none") return 1;
  if(o.model==="power") return Math.pow(KE,o.n!=null?o.n:0.6);
  if(o.model==="tpp2m") return tpp2m(KE,o.mat);
  return 1;
}
/* rows=[{label, el, area, sArea(σ), S, BE, hv, PE, tf}], opts={imfp} */
function quantify(rows,opts){
  opts=opts||{};
  const out=rows.map(r=>{ const KE=r.hv-r.BE; const T=tfValue(r.tf,KE,r.PE); const L=imfpValue(opts.imfp,KE); const den=r.S*T*L; return Object.assign({},r,{KE,T,L,den,nrm:r.area/den,sNrm:r.sArea!=null?r.sArea/den:null}); });
  const tot=out.reduce((s,r)=>s+(isFinite(r.nrm)?r.nrm:0),0);
  out.forEach((r,i)=>{ r.at=100*r.nrm/tot;
    // 델타: ∂at_i/∂n_j = 100(δij·tot − n_i)/tot²
    let v=0; out.forEach((q,j)=>{ if(q.sNrm==null) return; const g=100*((i===j?tot:0)-r.nrm)/(tot*tot); v+=g*g*q.sNrm*q.sNrm; }); r.sAt=Math.sqrt(v); });
  return out;
}

/* ============================================================
   10b) 화학량론 비교 · Auger parameter (SPEC §10)
   ============================================================ */
/* 화학식 해석: "NaLuF4", "Na(Lu0.92Tm0.08)F4", "Na3LuCl6", "[Ba0.9Sr0.1]TiO3" → {원소: 개수}.
   기호는 syms(주기율표 목록)로 검사 — 두 글자 기호를 먼저, 없으면 한 글자. 괄호 뒤 숫자는 괄호 안 전체에 곱함. */
function parseFormula(s,syms){
  const str=String(s==null?"":s).replace(/\s+/g,"");
  if(!str) throw new Error("화학식이 비었습니다");
  const set=new Set(syms||[]); let i=0;
  const num=()=>{ const m=/^(\d+(\.\d*)?|\.\d+)/.exec(str.slice(i)); if(!m) return 1; i+=m[0].length; return parseFloat(m[0]); };
  const sym=()=>{ const two=str.slice(i,i+2), one=str[i];
    if(/^[A-Z][a-z]$/.test(two)&&(!set.size||set.has(two))) return two;
    if(/^[A-Z]$/.test(one)&&(!set.size||set.has(one))) return one;
    if(/^[A-Z]/.test(one)) throw new Error(`'${/^[A-Z][a-z]?/.exec(str.slice(i))[0]}' 는 원소 기호가 아닙니다`);
    throw new Error(`${i+1} 번째 글자 '${str[i]}' 를 읽을 수 없습니다 (원소 기호·숫자·괄호만)`); };
  const group=depth=>{ const out={};
    while(i<str.length){ const ch=str[i];
      if(ch===")"||ch==="]"){ if(!depth) throw new Error("닫는 괄호가 남습니다"); break; }
      if(ch==="("||ch==="["){ const close=ch==="("?")":"]"; i++; const g=group(depth+1); if(str[i]!==close) throw new Error("괄호가 닫히지 않았습니다"); i++; const k=num(); for(const e in g) out[e]=(out[e]||0)+g[e]*k; continue; }
      const e=sym(); i+=e.length; const k=num(); out[e]=(out[e]||0)+k; }
    return out; };
  const r=group(0); if(i!==str.length) throw new Error("괄호 짝이 맞지 않습니다");
  const ks=Object.keys(r); if(!ks.length) throw new Error("원소가 없습니다");
  ks.forEach(e=>{ if(!(r[e]>0)) throw new Error(`'${e}' 의 개수가 0 입니다`); });
  return r;
}
/* 측정 ↔ 화학식 비교. rows = [{el, nrm, sNrm}] (quantify 결과의 A/(S·T·λ) 와 그 σ; 원소당 한 행).
   화학식에 있고 측정된 원소끼리만 다시 정규화: x_i = n_i / Σ n_j, 기대 e_i = ν_i / Σ ν_j (ν = 화학식 개수).
   σ(x_i): 넓이끼리 독립이라 보고 델타 방법 — ∂x_i/∂n_j = (δ_ij·W − n_i)/W².
   기준 원소 비: (n_i/n_ref)/(ν_i/ν_ref), 상대 σ = √((σ_i/n_i)² + (σ_ref/n_ref)²). 감도계수·T·λ 의 계통 오차는 들어 있지 않다. */
function stoichCompare(rows,formula,refEl){
  const byEl={}, dup=[];
  (rows||[]).forEach(r=>{ if(!(isFinite(r.nrm)&&r.nrm>0)) return; if(byEl[r.el]) dup.push(r.el); else byEl[r.el]=r; });
  const fEls=Object.keys(formula||{});
  const inF=fEls.filter(e=>byEl[e]), missing=fEls.filter(e=>!byEl[e]), outside=Object.keys(byEl).filter(e=>!(e in formula));
  const W=inF.reduce((s,e)=>s+byEl[e].nrm,0), N=inF.reduce((s,e)=>s+formula[e],0);
  const sq=v=>v*v;
  const ref=refEl&&inF.includes(refEl)?refEl:null;
  const items=inF.map(e=>{ const w=byEl[e].nrm; const x=w/W; let v=0, noSig=false;
    inF.forEach(j=>{ const sj=byEl[j].sNrm; if(sj==null){ noSig=true; return; } const g=((j===e?W:0)-w)/(W*W); v+=g*g*sj*sj; });
    const ex=formula[e]/N, sx=Math.sqrt(v);
    const it={el:e,x,sx:noSig?null:sx,ex,R:x/ex,sR:noSig?null:sx/ex};
    if(ref){ const a=byEl[e], b=byEl[ref]; const rat=a.nrm/b.nrm, exp=formula[e]/formula[ref];
      it.ratio=rat; it.ratioExp=exp; it.ratioRel=rat/exp;
      it.sRatioRel=e===ref?0:(a.sNrm==null||b.sNrm==null?null:Math.sqrt(sq(a.sNrm/a.nrm)+sq(b.sNrm/b.nrm))*rat/exp); }
    return it; });
  return {items,missing,outside,dup:[...new Set(dup)],ref,W,N};
}
/* 수정 Auger parameter α′ = BE(광전자선) + KE(Auger 선) = BE_p + hν − BE_A,겉보기.
   BE 눈금의 스펙트럼에서 Auger 봉우리의 겉보기 BE 를 쓰면 대전 이동 δ 가 두 항에서 같은 부호로 들어와 빠진다
   (BE_p = BE_p⁰ + δ, BE_A,겉보기 = hν − KE_A⁰ + δ) → α′ 는 대전 보정과 무관. 두 스펙트럼은 같은 광원(hν)이어야 한다. */
function augerParam(beP,beAapp,hv,sP,sA){
  const a=beP+hv-beAapp; const s=sP!=null&&sA!=null?Math.sqrt(sP*sP+sA*sA):null; return {alpha:a,s,keA:hv-beAapp};
}
/* 구간 [lo, hi] 안 최댓값 위치: 국소 2차 평활(반창 ≈ 0.25 eV) 뒤 가장 큰 점을 세 점 포물선으로 다듬음. 끝점이면 null(구간 밖에 더 큰 값). */
function peakMaxPos(E,y,lo,hi,hwE){
  const n=E.length; if(n<5) return null; const st=Math.abs((E[n-1]-E[0])/(n-1))||1;
  const ys=smoothQ(E,y,Math.max(2,Math.round((hwE||0.25)/st)));
  let im=-1; for(let i=0;i<n;i++){ if(E[i]<lo||E[i]>hi) continue; if(im<0||ys[i]>ys[im]) im=i; }
  if(im<=0||im>=n-1) return null;
  const e1=Math.min(E[im-1],E[im+1]), e2=Math.max(E[im-1],E[im+1]); if(e1<lo||e2>hi) return null;   // 이웃이 구간 밖 = 최댓값이 구간 끝(오름·내림차순 모두)
  const a=ys[im-1], b=ys[im], c=ys[im+1]; const d=a-2*b+c; const off=d<0?Math.max(-0.5,Math.min(0.5,0.5*(a-c)/d)):0;
  return {be:E[im]+off*(E[im+1]-E[im]), y:b, i:im};
}

/* ============================================================
   11) 서베이: 봉우리 찾기 · 원소 판별
   ============================================================ */
function findPeaks(E,y,opt){
  opt=opt||{}; const n=y.length; if(n<9) return [];
  const st=(E[n-1]-E[0])/(n-1);
  const hw=Math.max(2,Math.round((opt.width||1.5)/st));
  const ys=smoothQ(E,y,hw);
  // 음의 2차 미분(평활) — 봉우리에서 양수
  const d2=new Float64Array(n); for(let i=1;i<n-1;i++) d2[i]=-(ys[i+1]-2*ys[i]+ys[i-1])/(st*st);
  const nz=noiseEstimate(E,y); const kap=nz?nz.kappa:1;
  const pk=[];
  for(let i=2;i<n-2;i++){
    if(!(ys[i]>=ys[i-1]&&ys[i]>=ys[i+1]&&d2[i]>0)) continue;
    // 높이 = 지형 prominence(scipy.signal.peak_prominences 와 같은 정의): 양쪽으로 이 봉우리보다 높은 점(또는 창 끝)까지 가며
    // 그 구간 최솟값을 구하고, 둘 중 큰 값을 기준선으로. 큰 봉우리 꼬리 위의 잡음 굴곡은 prominence 가 작아 걸러진다.
    const w=Math.max(hw*3,Math.round((opt.window||4)/st)); const a=Math.max(0,i-w), b=Math.min(n-1,i+w);
    let la=ys[i], lb=ys[i];
    for(let j=i-1;j>=a;j--){ if(ys[j]>ys[i]) break; if(ys[j]<la) la=ys[j]; }
    for(let j=i+1;j<=b;j++){ if(ys[j]>ys[i]) break; if(ys[j]<lb) lb=ys[j]; }
    const base=Math.max(la,lb); const h=ys[i]-base; if(!(h>0)) continue;
    // 평활(2차 다항식, N=2hw+1 점) 뒤 한 점의 잡음 σ ≈ 1.5·σ/√N, prominence 는 두 평활값의 차라 ×√2
    const sd=1.5*Math.SQRT2*Math.sqrt(kap*Math.max(ys[i],1)/(2*hw+1));
    const snr=h/sd; if(snr<(opt.snr||5)) continue;
    pk.push({i,be:E[i],h,snr,y:ys[i]});
  }
  // 가까운 봉우리 병합(1 eV 안 큰 것만)
  pk.sort((a,b)=>b.h-a.h); const keep=[]; pk.forEach(p=>{ if(!keep.some(q=>Math.abs(q.be-p.be)<(opt.merge||1.2))) keep.push(p); });
  return keep.sort((a,b)=>a.be-b.be);
}
/* 서베이 원소 판별 — SPEC §7  (가장 큰 미설명 봉우리부터 설명하는 탐욕적 배정, 수작업 판별 순서와 같음)
   1) 원소마다 광전자 선 목록: 위치 = Avantage 원소표(idpos) → 없으면 Booklet(진공 기준 기체 원소 제외), 가중 w = σ/σmax(Scofield).
      같은 궤도 doublet 이 허용폭 tol 보다 가까우면 하나로 합침. Auger 선(Avantage KE)은 위치만.
   2) 미설명 봉우리 중 가장 큰 것 P 에 대해, 센 선(w ≥ 0.3)이 P 와 맞는 원소마다 δ = P − 선 위치(δ ∈ [δmin, δmax])로
      score = Σ_맞은 선 w(1 − ½d/tol)·(이미 설명된 봉우리면 0.3) − 0.6 Σ_빠진 선(w>0.2) w,  × exp(−(δ/5)²)
      최고 점수 원소가 accept 문턱(0.5)을 넘으면 채택하고 그 원소의 맞은 봉우리(광전자·Auger)를 '설명됨'으로 표시.
      없으면 P 를 '미배정'으로 두고 다음 봉우리로.
   3) 채택 원소 뒤에 나머지 후보(최고 점수 순)를 붙여 돌려준다. accepted=true 가 기본 체크. */
function elementLines(el,rec,lib,hv,tol,src){
  const sc=(lib.scofield&&lib.scofield[src])||{}, idp=lib.idpos||{};
  const lines=[];
  Object.entries(rec.lv||{}).forEach(([sub,bk])=>{ const s=sc[el+" "+sub]; if(s==null) return;
    const pos=idp[el+" "+sub]!=null?idp[el+" "+sub]:(rec.vac?null:bk); if(pos==null||!(pos>0)||pos>hv-30) return; lines.push({name:el+" "+sub,be:pos,s}); });
  lines.sort((a,b)=>a.be-b.be); const merged=[];
  lines.forEach(l=>{ const orb=l.name.split(" ")[1].slice(0,2); const m=merged.find(q=>q.orb===orb&&Math.abs(q.be-l.be)<tol); if(m){ m.be=(m.be*m.s+l.be*l.s)/(m.s+l.s); m.s+=l.s; m.name=el+" "+orb; } else merged.push(Object.assign({orb},l)); });
  const sMax=merged.length?Math.max(...merged.map(l=>l.s)):0; merged.forEach(l=>l.w=sMax>0?l.s/sMax:0);
  const aug=(lib.auger||[]).filter(a=>a.el===el).map(a=>({name:el+" "+a.line+" (Auger)",be:hv-a.ke})).filter(a=>a.be>0&&a.be<hv-30);
  return {pe:merged,aug};
}
/* 기본 제외: XPS 선이 없거나(H, He) 실험실 시료에 사실상 없는 방사성 원소 — 화면에서 바꿀 수 있음 */
const EXCLUDE_DEFAULT=["H","He","Tc","Pm","Po","At","Rn","Fr","Ra","Ac","Pa","Np","Pu","Am","Cm","Bk","Cf","Es","Fm"];
function identify(peaks,lib,o){
  o=o||{}; const hv=o.hv||1486.6, src=o.source||(hv<1300?"mg":"al"), shift=o.shift||0, tol=o.tol||2.0;
  const dmin=o.dmin!=null?o.dmin:-4, dmax=o.dmax!=null?o.dmax:6, thr=o.accept!=null?o.accept:0.5, trig=o.trigPow!=null?o.trigPow:1;
  const excl=o.exclude||EXCLUDE_DEFAULT;
  const obs=peaks.map((p,j)=>({j,be:p.be-shift,snr:p.snr,h:p.h}));
  const EL={}; Object.entries(lib.core||{}).forEach(([el,rec])=>{ if(excl.includes(el)) return; const L=elementLines(el,rec,lib,hv,tol,src); if(L.pe.length) EL[el]=L; });
  const explained=new Set();
  const evalAt=(el,d)=>{ const L=EL[el]; const used=new Set(); let score=0; const found=[], missing=[];
    L.pe.forEach(l=>{ let best=null; obs.forEach(p=>{ if(used.has(p.j)) return; const dd=Math.abs(p.be-(l.be+d)); if(dd<=tol&&(!best||dd<best.dd)) best={p,dd}; });
      if(best){ used.add(best.p.j); score+=l.w*(1-0.5*best.dd/tol)*(explained.has(best.p.j)?0.3:1); found.push({line:l.name,be:l.be,obs:best.p.be,snr:best.p.snr,d:best.dd,j:best.p.j,w:l.w}); }
      else if(l.w>0.2){ score-=0.6*l.w; missing.push(l.name); } });
    L.aug.forEach(a=>{ let best=null; obs.forEach(p=>{ if(used.has(p.j)) return; const dd=Math.abs(p.be-a.be); if(dd<=tol+2&&(!best||dd<best.dd)) best={p,dd}; }); if(best){ used.add(best.p.j); found.push({line:a.name,be:a.be,obs:best.p.be,snr:best.p.snr,d:best.dd,j:best.p.j,auger:true}); } });
    return {score:score*Math.exp(-Math.pow(d/5,2)),d,found,missing}; };
  // P 를 설명하는 선의 세기 w_P 를 곱한다: 큰 봉우리는 어떤 원소의 '센 선'일 가능성이 크다(약한 선의 우연 일치 억제)
  const bestFor=(el,P)=>{ let best=null; EL[el].pe.forEach(l=>{ if(l.w<0.3) return; const d=P.be-l.be; if(d<dmin-tol||d>dmax+tol) return; const dc=Math.max(dmin,Math.min(dmax,d)); const r=evalAt(el,dc); const fP=r.found.find(f=>f.j===P.j); if(!fP) return; r.score*=Math.pow(fP.w||0,trig); if(!best||r.score>best.score) best=r; }); return best; };
  const accepted=[], unassigned=[];
  const order=obs.slice().sort((a,b)=>b.h-a.h);
  for(const P of order){
    if(explained.has(P.j)) continue;
    let win=null; Object.keys(EL).forEach(el=>{ if(accepted.some(c=>c.el===el)) return; const r=bestFor(el,P); if(r&&r.found.some(f=>f.j===P.j)&&(!win||r.score>win.r.score)) win={el,r}; });
    if(win&&win.r.score>thr){ accepted.push({el:win.el,score:win.r.score,shift:+win.r.d.toFixed(2),found:win.r.found,missing:win.r.missing,accepted:true,by:P.be}); win.r.found.forEach(f=>explained.add(f.j)); }
    else unassigned.push(P.be);
  }
  // 나머지 후보: 모든 봉우리 기준 최고 점수(이미 설명된 봉우리는 0.3 배)
  const rest=[];
  Object.keys(EL).forEach(el=>{ if(accepted.some(c=>c.el===el)) return; let best=null; obs.forEach(P=>{ const r=bestFor(el,P); if(r&&r.found.some(f=>f.j===P.j)&&(!best||r.score>best.score)) best=r; }); if(best&&best.score>0) rest.push({el,score:best.score,shift:+best.d.toFixed(2),found:best.found,missing:best.missing,accepted:false,overlap:best.found.filter(f=>!f.auger).every(f=>explained.has(f.j))}); });
  rest.sort((a,b)=>b.score-a.score);
  // 채택되지 않았지만 점수가 문턱을 넘고 자기만의(미설명) 광전자 봉우리가 있는 후보 = '가능'(겹친 봉우리 뒤에 숨은 원소일 수 있음)
  rest.forEach(c=>{ c.possible=c.score>thr&&!c.overlap; });
  const out=accepted.concat(rest);
  out.forEach(c=>{ c.nPE=c.found.filter(f=>!f.auger).length; c.found.forEach(f=>{ delete f.j; delete f.w; }); if(c.overlap==null) c.overlap=false; });
  out.unassigned=unassigned;
  return out;
}
/* 기본 채택(화면의 체크 기본값) = 탐욕적 배정에서 채택된 원소 */
function acceptDefault(c){ return !!c.accepted; }

/* ============================================================
   12) 원소에 독립적인 모형 출발점 (SPEC §9) — 어떤 원소·궤도든 같은 규칙
   ============================================================ */
const ORB_JJ={p:["3/2","1/2"],d:["5/2","3/2"],f:["7/2","5/2"]}, ORB_RR={p:0.5,d:2/3,f:0.75};
/* doublet 정보: Δ 후보 두 출처(Avantage 원소표 · X-ray Data Booklet), 넓이 비 (2j+1), 분해 가능 여부 */
/* 다중항(multiplet)이 지배하는 준위: 열린 껍질(3d 전이금속의 3d, 희토류의 4f)과 내각 구멍의 교환 상호작용이
   스핀–궤도 갈라짐과 비슷하거나 커서 '두 봉우리(doublet)' 모형이 물리적으로 맞지 않는 준위.
   - Ti–Cu 의 3s·3p (예: Mn 3s 교환 갈라짐, Co·Ni 3p 다중항; 닫힌 껍질 상태인 Ti⁴⁺ 3d⁰·Cu⁺/Cu⁰ 3d¹⁰ 는 예외)
   - Ce–Yb 의 4s·4p·4d·5s·5p (예: Gd·Tb·Tm 4d 넓은 다중항; La 4f⁰·Lu 4f¹⁴·Yb²⁺ 4f¹⁴ 는 예외)
   이런 준위는 화학 상태 분해를 하지 않고 넓은 성분들로 '모양(envelope)'만 맞춰 넓이(정량)를 쓴다. 원소·궤도 규칙이며 특정 데이터에 맞춘 것이 아니다. */
const MULT_TM=["Ti","V","Cr","Mn","Fe","Co","Ni","Cu"], MULT_RE=["Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"];
function multipletInfo(el,orb){
  if(MULT_TM.includes(el)&&(orb==="3p"||orb==="3s"))
    return {kind:"3d",why:`${el} ${orb}: 열린 3d 껍질과의 교환 상호작용(다중항)이 스핀–궤도 갈라짐과 비슷 — doublet 모형이 맞지 않음`,closed:el==="Ti"?"Ti⁴⁺ (3d⁰)":el==="Cu"?"Cu⁺·Cu⁰ (3d¹⁰)":null};
  if(MULT_RE.includes(el)&&["4s","4p","4d","5s","5p"].includes(orb))
    return {kind:"4f",why:`${el} ${orb}: 열린 4f 껍질과의 교환 상호작용(다중항)으로 넓게 갈라짐 — doublet 모형이 맞지 않음`,closed:el==="Yb"?"Yb²⁺ (4f¹⁴)":el==="Ce"?"Ce⁴⁺ (형식상 4f⁰, 실제로는 4f¹L 종상태가 섞임)":null};
  return null;
}
function doubletInfo(lib,el,orb){
  if(!el||!orb) return null; const l=orb.slice(-1), J=ORB_JJ[l]; if(!J) return null;
  const rec=(lib.core||{})[el]||{}, lv=rec.lv||{}, idp=lib.idpos||{};
  const bk=(lv[orb+J[0]]>0&&lv[orb+J[1]]>0&&!(rec.soWarn||[]).includes(orb))?+(lv[orb+J[1]]-lv[orb+J[0]]).toFixed(3):null;
  const av=(idp[el+" "+orb+J[0]]!=null&&idp[el+" "+orb+J[1]]!=null&&idp[el+" "+orb+J[1]]>idp[el+" "+orb+J[0]])?+(idp[el+" "+orb+J[1]]-idp[el+" "+orb+J[0]]).toFixed(3):null;
  const bkOK=bk!=null&&bk>0;
  let dE=null, src="";
  if(av!=null){ dE=av; src="Avantage 원소표 "+el+" "+orb+J[1]+" − "+orb+J[0]; }
  else if(bkOK){ dE=bk; src="X-ray Data Booklet "+el+" "+orb+J[1]+" − "+orb+J[0]; }
  return {dE,r:ORB_RR[l],J,src:dE!=null?src:"라이브러리에 쓸 수 있는 값 없음",alt:{avantage:av,booklet:bkOK?bk:null},
    spread:(av!=null&&bkOK)?Math.abs(av-bk):null, resolvable:dE!=null&&dE>=DOUBLET_RESOLVE, multiplet:multipletInfo(el,orb)};
}
const DOUBLET_RESOLVE=1.2;   // eV — 이 이상이면 두 봉우리가 데이터에서 갈라져 Δ 를 데이터로 정할 수 있다
/* Δ 제약: 라이브러리 Δ ± max(0.3 eV, 5 %, 두 출처 차 + 0.1) 안에서 자유. 라이브러리 Δ 의 정밀도가 0.1–0.5 eV 라
   고정하면 분해되지 않는 doublet 에서도 계통 편향·가짜 성분이 생긴다(합성 시험: Rb 3d 라이브러리 1.0 vs 참 1.07 → 20 회 중 3 회 가짜 성분).
   여러 화학 상태의 doublet 은 첫 성분의 Δ 에 연결(매개변수 1 개만 늘어남). */
function dEConstraint(info){
  if(!info||info.dE==null) return {mode:"free",lo:0.05,hi:40};
  const w=Math.max(0.3,0.05*info.dE,(info.spread||0)+0.1); return {mode:"free",lo:+Math.max(0.05,info.dE-w).toFixed(3),hi:+(info.dE+w).toFixed(3)};
}
/* 출발 모형 제안: 봉우리 찾기 → (짝이 보이거나 짝 위치가 영역 안이면) doublet → 넓이는 높이에서.
   o = {bg, el, orb, lib, shape, F}. 반환 {comps, cons, info, peaks} */
function suggestModel(E,y,o){
  o=o||{}; const n=E.length; const bg=o.bg||{type:"smart",n:5};
  const B=bg.type==="active"||bg.type==="none"?bgShirley(E,y,bg.n||5,true):bgStatic(E,y,bg);
  const z=Float64Array.from(y,(v,i)=>v-B[i]);
  const ix=be=>{ let i=Math.round((be-E[0])/((E[n-1]-E[0])/(n-1))); return Math.max(0,Math.min(n-1,i)); };
  let pk=findPeaks(E,y,{width:o.width||0.5,snr:o.snr||4,merge:o.merge||0.9}).filter(q=>q.be>E[3]&&q.be<E[n-4]);
  pk=pk.map(q=>({be:q.be,h:z[ix(q.be)]})).filter(q=>q.h>0);
  pk=pk.sort((a,b)=>b.h-a.h).slice(0,o.maxComps||8);
  // 가장 큰 봉우리의 2 % 미만은 출발 모형에서 뺀다(작은 실제 성분은 잔차 보강이 근거 있을 때 다시 더함)
  if(pk.length){ const hm=pk[0].h; pk=pk.filter(q=>q.h>=(o.minRel!=null?o.minRel:0.02)*hm); }
  const info=doubletInfo(o.lib||{},o.el,o.orb); const mult=multipletInfo(o.el,o.orb);
  const comps=[], cons={}; if(!pk.length) return {comps,cons,info,mult,peaks:[]};
  const asc=pk.slice().sort((a,b)=>a.be-b.be), used=new Set();
  // 다중항 준위: doublet 을 쓰지 않고 넓은 성분(F 출발 2 eV, 0.3–8 eV)으로 모양만
  const F=o.F||(mult?2.0:1.2), shape=o.shape||"GLS", Fhi=mult?8:4;
  asc.forEach(q=>{
    if(used.has(q)) return; used.add(q);
    let partner=null;
    if(info&&info.dE&&!mult){ const tol=Math.max(0.7,0.08*info.dE); partner=asc.find(t=>t!==q&&!used.has(t)&&Math.abs(t.be-(q.be+info.dE))<tol&&t.h<q.h*1.3); }
    // 이 봉우리가 '보이지 않는' 주선의 짝(높은 BE 쪽)일 수 있다: 겹친 다중항에서 다른 상태의 주선이 앞 doublet 의 짝과
    // 합쳐져 봉우리로 따로 안 잡히는 경우(예: Au 4f 두 상태 84.0·86.5 → 86.5 와 87.67 이 한 덩어리, 90.17 만 따로 보임).
    // 짝이 없고, Δ 아래 자리에 이 봉우리 이상의 신호가 있고(주선은 짝보다 크다: 높이 비 1/r > 1), 그 자리에 이미 성분이 없으면 주선을 그 자리에 둔다.
    // 판단은 '이미 둔 성분(짝 포함)이 그 자리에서 설명하는 몫을 뺀 나머지'로 — 앞 doublet 의 짝·꼬리를 새 주선으로 착각하지 않게(검토에서 찾은 반례:
    // Au 4f 84 eV + 91 eV 의 작은 단일 봉우리 → 87.3 eV 에 가짜 doublet). 주선 자리가 기존 성분 중심·짝 위치와 0.5 eV 안이면 두지 않는다.
    let beMain=q.be, hMain=q.h;
    if(!partner&&info&&info.dE&&!mult&&comps.length){ const mb=q.be-info.dE;
      if(mb>E[0]+0.5&&!comps.some(c=>Math.abs(c.be-mb)<0.5||(c.doublet&&Math.abs(c.be+c.doublet.dE-mb)<0.5))){
        const ctx={emin:E[0],emax:E[n-1]}, Em=Float64Array.of(mb); let pred=0; comps.forEach(c=>{ pred+=compEval(c,Em,ctx)[0]; });
        const res=z[ix(mb)]-pred; if(res>=q.h){ beMain=mb; hMain=res; } } }
    const c=newComp(shape,beMain,{F,el:o.el||"",line:o.orb||""});
    const cs=cons[c.id]={};
    if(info&&info.dE&&!mult&&(partner||beMain!==q.be||q.be+info.dE<E[n-1]-0.5)){
      c.doublet={dE:info.dE,r:info.r,q:1}; cs.r={mode:"fixed"}; cs.q={mode:"fixed"};
      if(partner) used.add(partner);
    }
    cs.F={mode:"free",lo:0.3,hi:Fhi}; cs.be={mode:"free",lo:E[0],hi:E[n-1]};
    const hh=trueFWHM(c.shape,shapeParams(c,1)).ymax; c.area=hMain/(hh||1)*(c.doublet?1+c.doublet.r:1);
    comps.push(c);
  });
  // Δ 는 가장 큰 doublet 이 기준(범위 안 자유), 나머지 doublet 은 거기에 연결. 모양 매개변수(m)도 가장 큰 성분에 연결 — 흔한 출발점
  const dA=anchorOf(comps.filter(c=>c.doublet)); comps.filter(c=>c.doublet).forEach(c=>{ cons[c.id].dE=c===dA?dEConstraint(info):{mode:"link",ref:dA.id,k:1,b:0}; });
  const mA=anchorOf(comps); if(mA&&mA.m!=null) comps.forEach(c=>{ if(c!==mA&&c.m!=null) cons[c.id].m={mode:"link",ref:mA.id,k:1,b:0}; });
  return {comps,cons,info,mult,peaks:pk};
}

/* 연결 기준(anchor) = 넓이가 가장 큰 성분 — 작은 성분을 기준으로 삼으면 그 성분이 죽을 때 전체가 퇴화한다 */
function anchorOf(cs){ let b=null; (cs||[]).forEach(c=>{ if(!b||c.area>b.area) b=c; }); return b; }
/* 성분 빼기 — 이 성분을 연결 기준으로 쓰던 매개변수는 남은 성분 중 가장 큰 것을 새 기준으로:
   새 기준은 뺀 성분의 제약(예: Δ 의 자유 범위)을 물려받고, 나머지는 새 기준에 연결 */
function removeComp(P,id){
  const X=P.comps.find(c=>c.id===id); if(!X) return; const xc=P.cons[id]||{};
  P.comps=P.comps.filter(c=>c.id!==id); delete P.cons[id];
  const byK={}; for(const cid in P.cons) for(const k in P.cons[cid]){ const q=P.cons[cid][k]; if(q&&q.mode==="link"&&q.ref===id) (byK[k]=byK[k]||[]).push(cid); }
  for(const k in byK){ const users=byK[k].map(cid=>P.comps.find(c=>c.id===cid)).filter(Boolean); const Y=anchorOf(users); if(!Y) continue;
    const q=P.cons[Y.id][k]; const own=xc[k]?clone(xc[k]):null;
    // 새 기준의 값: 연결식 그대로 (k·X + b) — 값은 이미 그렇게 들어 있으므로 제약만 바꾼다
    if(own&&own.mode!=="link") P.cons[Y.id][k]=own; else delete P.cons[Y.id][k];
    // 나머지는 새 기준에 연결. 연결식이 k=1, b=0(자동 제안이 만드는 꼴)이면 정확히 같은 관계, 아니면 근사 — 출발값이라 피팅이 다시 정한다
    users.forEach(u=>{ if(u!==Y) P.cons[u.id][k].ref=Y.id; }); }
}

/* 잔차로 성분 더하기 (SPEC §9): 극대가 아닌 어깨(예: C 1s 의 C–O)는 봉우리 찾기로 안 보이므로
   맞춘 뒤 정규화 잔차를 성분 폭 정도로 평활해 가장 큰 양의 덩어리(z)에 성분을 더하고, BIC(데이터 가중 χ²_d)가
   10 이상 줄 때만 받아들인다(강한 증거). 원소·궤도와 무관한 규칙. 넓이 몫 < 5 % 인 성분은 빼 보고 ΔBIC ≤ 2 면 뺀다. 반환: {prob, res, stats, added, removed, log} */
function refineByResidual(prob,o){
  o=o||{}; const maxAdd=o.maxAdd!=null?o.maxAdd:4, zMin=o.zMin||4, dBIC=o.dBIC||10;
  const E=prob.E, n=E.length, st=(E[n-1]-E[0])/(n-1);
  // 성분을 더할지는 유연한 배경(피팅 안 Shirley + 기울기)으로 판단한다 — 정적 배경의 모양 오차가 잔차 덩어리로 남아
  // '가짜 성분'으로 채택되는 것을 막기 위해(합성 시험: 정적 배경으로 판단하면 참 1 성분 스펙트럼 143 개 중 106 개에 성분이 더해짐)
  let cur=clone(stripProb(prob)); cur.E=E; cur.y=prob.y; const userBg=clone(prob.bg);
  if(o.decisionBg!==false) cur.bg={type:"active",n:(prob.bg&&prob.bg.n)||5,cons:{s:{mode:"free"}}};
  let r=fitOnce(cur,{}); let S=fitStats(cur,r); const log=[`시작 ${cur.comps.length} 성분 BIC ${S.BIC.toFixed(1)}`]; let added=0, removed=0;
  // 죽은 성분 덜기: 넓이 몫이 작거나(< 5 %) 중심이 창 밖으로 나간 성분을 빼 보고 BIC 가 나빠지지 않으면(ΔBIC ≤ 2) 뺀다.
  // 자동 제안이 겹친 다중항의 위쪽 짝(예: Au 4f5/2 + 다른 상태 4f7/2)을 별도 성분으로 잡으면 피팅이 그 성분을 넓이 0 으로 보내는데,
  // 그대로 두면 퇴화(공분산 특이)와 '성분 수' 오판이 생기므로 여기서 정리한다. 다른 성분이 연결(link)로 참조하던 성분이면 연결 기준을 남은 성분 중 가장 큰 것으로 옮긴다(removeComp).
  // (실측 Ba 4d: 0.7 % 높이의 작은 봉우리가 첫 doublet = Δ·m 연결 기준이 되어 넓이 0 으로 창 밖에 남던 결함)
  const prune=()=>{ for(let pass=0;pass<4&&cur.comps.length>1;pass++){
      const tot=cur.comps.reduce((s,c)=>s+Math.max(0,c.area),0)||1; const refd=new Set();
      for(const id in cur.cons) for(const k in cur.cons[id]){ const q=cur.cons[id][k]; if(q&&q.mode==="link"&&q.ref) refd.add(q.ref); }
      const out=c=>c.be<E[0]||c.be>E[n-1];
      const order=cur.comps.map((c,i)=>i).filter(i=>cur.comps[i].area/tot<0.05||out(cur.comps[i])).sort((a,b)=>cur.comps[a].area-cur.comps[b].area);
      let done=false;
      for(const i of order){ const cand=clone(stripProb(cur)); cand.E=E; cand.y=prob.y; removeComp(cand,cand.comps[i].id);
        let r2,S2; try{ r2=fitOnce(cand,{initAreas:false}); S2=fitStats(cand,r2); }catch(e){ continue; }
        if(S2.BIC<=S.BIC+2){ log.push(`${cur.comps[i].be.toFixed(2)} eV 성분(넓이 ${(100*cur.comps[i].area/tot).toFixed(1)}%) 뺌 → ΔBIC ${(S.BIC-S2.BIC).toFixed(1)}`); cur=cand; r=r2; S=S2; removed++; done=true; break; } }
      if(!done) break; } };
  // 한 걸음마다 후보 수를 모두 맞춰 보고 BIC 가 가장 많이 주는 것 하나를 받는다(≥ 10 일 때만):
  //  (가) 잔차 덩어리(z ≥ 4)에 성분 더하기
  //  (나) doublet 넓이 비 r 풀기 / (다) 폭 배율 q 풀기 — 통계 비 (2j+1) 는 얕은 준위(예: Ba 4d)에서 어긋나는 일이 흔하다(광이온화 단면적의 j 의존).
  //      모든 화학 상태에 같은 값(가장 큰 doublet 에 연결). 분해되지 않는 doublet(Δ < 1.2 eV)은 폭과 퇴화하므로 풀지 않음.
  // (나)(다)를 성분 더하기와 '경쟁'시키는 까닭: 비가 틀린 채 두면 1/2 자리 잔차를 가짜 성분이 메우고(합성 양성 대조),
  // 반대로 비를 먼저 풀면 빠진 화학 상태를 비가 흡수할 수 있다 — 어느 쪽이 데이터를 더 잘 설명하는지로 고른다.
  const relaxed=[];
  const relaxCand=(k,lo,hi)=>{ if(o.relax===false) return null; const dA=anchorOf(cur.comps.filter(c=>c.doublet&&c.doublet.dE>=DOUBLET_RESOLVE)); if(!dA) return null;
    const cc=(cur.cons[dA.id]||{})[k]; if(cc&&cc.mode!=="fixed") return null;
    const cand=clone(stripProb(cur)); cand.E=E; cand.y=prob.y; const v0=dA.doublet[k];
    cand.comps.filter(c=>c.doublet&&c.doublet.dE>=DOUBLET_RESOLVE).forEach(c=>{ cand.cons[c.id]=cand.cons[c.id]||{};
      cand.cons[c.id][k]=c.id===dA.id?{mode:"free",lo:+(v0*lo).toFixed(4),hi:+(v0*hi).toFixed(4)}:{mode:"link",ref:dA.id,k:1,b:0}; c.doublet[k]=v0; });
    const r2=fitOnce(cand,{initAreas:false}), S2=fitStats(cand,r2); const nv=cand.comps.find(c=>c.id===dA.id).doublet[k];
    return {cand,r2,S2,kind:k,txt:`${k==="r"?"넓이 비 r":"폭 배율 q"} 풀기 (${v0.toFixed(3)} → ${nv.toFixed(3)})`}; };
  const addCand=()=>{ const w=Math.max(2,Math.round(0.6/st)); const z=new Float64Array(n);
    for(let i=0;i<n;i++){ let sum=0,c=0; for(let j=Math.max(0,i-w);j<=Math.min(n-1,i+w);j++){ sum+=S.resid[j]; c++; } z[i]=sum/Math.sqrt(c); }
    // ρ₁ 이 크면 이웃 잔차가 상관 → 평활 합의 분산이 커지므로 그만큼 나눈다
    const inf=Math.sqrt((1+Math.max(0,S.rho1))/(1-Math.min(0.95,Math.max(0,S.rho1)))); let im=-1, zm=0;
    for(let i=w;i<n-w;i++){ const zz=z[i]/inf; if(zz>zm){ zm=zz; im=i; } }
    if(im<0||zm<zMin){ log.push(`잔차 최대 z ${zm.toFixed(1)} < ${zMin}`); return null; }
    const be=E[im]; const cand=clone(stripProb(cur)); cand.E=E; cand.y=prob.y;
    const ref=anchorOf(cand.comps); const c=newComp(ref?ref.shape:"GLS",be,{F:1.2,el:ref?ref.el:"",line:ref?ref.line:""});
    if(ref&&ref.m!=null) c.m=ref.m;
    const Fr=ref&&(cand.cons[ref.id]||{}).F; // 기준 성분의 폭 한계를 물려받음(다중항 envelope 은 0.3–8 eV)
    const cs=cand.cons[c.id]={F:Fr&&Fr.mode==="free"?{mode:"free",lo:Fr.lo,hi:Fr.hi}:{mode:"free",lo:0.3,hi:4},be:{mode:"free",lo:E[0],hi:E[n-1]}};
    if(ref&&ref.m!=null) cs.m={mode:"link",ref:ref.id,k:1,b:0};
    const d0=anchorOf(cand.comps.filter(q=>q.doublet));
    if(d0&&be+d0.doublet.dE<E[n-1]-0.5){ c.doublet={dE:d0.doublet.dE,r:d0.doublet.r,q:d0.doublet.q};
      const d0c=cand.cons[d0.id]||{}; ["dE","r","q"].forEach(k=>{ const dc=d0c[k]; cs[k]=!dc||dc.mode==="fixed"?{mode:"fixed"}:{mode:"link",ref:d0.id,k:1,b:0}; }); }
    // 넓이 출발: 잔차 덩어리 높이
    const hh=trueFWHM(c.shape,shapeParams(c,1)).ymax; c.area=Math.max(1e-9,(prob.y[im]-r.mu[im]))/(hh||1)*(c.doublet?1+c.doublet.r:1);
    cand.comps.push(c);
    const r2=fitOnce(cand,{initAreas:false}), S2=fitStats(cand,r2);
    return {cand,r2,S2,kind:"add",txt:`${be.toFixed(2)} eV 성분 추가 (z ${zm.toFixed(1)})`}; };
  prune();
  for(let it=0;it<maxAdd+2;it++){
    const cands=[]; const tryC=f=>{ try{ const x=f(); if(x) cands.push(x); }catch(e){ log.push("후보 피팅 실패: "+e.message); } };
    if(added<maxAdd) tryC(addCand); tryC(()=>relaxCand("r",0.75,1.33)); tryC(()=>relaxCand("q",0.6,1.6));
    if(!cands.length) break;
    cands.forEach(x=>x.dB=S.BIC-x.S2.BIC); cands.sort((a,b)=>b.dB-a.dB); const best=cands[0];
    log.push(cands.map(x=>`${x.txt} ΔBIC ${x.dB.toFixed(1)}`).join(" · "));
    if(best.dB<dBIC){ log.push(`가장 나은 수도 ΔBIC < ${dBIC} → 멈춤`); break; }
    log.push(`→ 채택: ${best.txt}`); cur=best.cand; r=best.r2; S=best.S2; if(best.kind==="add") added++; else relaxed.push(best.kind);
  }
  if(added) prune();
  // 사용자 배경으로 되돌려 다시 맞춘다(성분 구성만 가져감)
  if(o.decisionBg!==false){ const fin=clone(stripProb(cur)); fin.E=E; fin.y=prob.y; fin.bg=userBg; try{ r=fitOnce(fin,{}); S=fitStats(fin,r); cur=fin; }catch(e){ log.push("사용자 배경 재피팅 실패: "+e.message); } }
  return {prob:cur,res:r,stats:S,added,removed,relaxed,log};
}

/* Worker 용: 잔차 보강 결과를 직렬화 가능한 형태로 */
function autoJob(prob,opt){ const rf=refineByResidual(prob,opt||{}); return {comps:rf.prob.comps,cons:rf.prob.cons,bg:rf.prob.bg,added:rf.added,removed:rf.removed,relaxed:rf.relaxed,log:rf.log}; }

/* ============================================================
   12) 모형 민감도 (SPEC §11): 같은 성분 구성을 '합리적인 분석 선택'의 조합으로 다시 맞춰
   관심 값(BE·넓이 %·FWHM)이 선택에 따라 얼마나 흔들리는지 — 분석자의 주관을 숫자로(모형 불확도).
   요인: 배경 종류 · 끝점 평균 점 수(정적 배경만) · 낮은/높은 BE 창 끝 ±dw.
   ============================================================ */
const SENS_BG_LABEL={smart:"Smart",shirley:"Shirley",active:"Shirley 피팅 안(+기울기)",tougaard:"Tougaard",linear:"Linear"};
function sensBg(base,type,n){
  // 지금 배경과 같은 종류면 사용자의 제약(예: 기울기 고정)을 그대로 — 0 번 변형이 사용자가 실제로 맞춘 모형이어야 한다.
  // 다른 종류로 바꿀 때만 그 종류의 기본값(피팅 안 Shirley 는 기울기 자유)
  const same=(base||{}).type===type; const b=clone(base||{}); b.type=type; b.n=n!=null?n:(b.n||5);
  if(type==="active"){ b.c0=null; b.k=null; if(!same){ b.s=0; b.cons={s:{mode:"free"}}; } }
  if(type==="tougaard"&&(base||{}).type!=="tougaard"){ b.C=1643; b.D=0; b.u3=false; }
  return b;
}
/* 변형 목록. model = {lo, hi, bg}; spec = {be}; o = {bgs:[type…], ns:[…], dw, mode:"oat"|"full"}.
   0 번 = 지금 설정(기준). 반환 [{lo,hi,bg,f:{bg,n,lo,hi}}] — f 는 요인 수준 이름(주효과 계산용). */
function sensVariants(model,spec,o){
  o=o||{}; const be=spec.be, E0=be[0], E1=be[be.length-1];
  const dw=o.dw!=null?o.dw:0.5; const cur=model.bg||{type:"smart",n:5}; const curN=cur.n||5;
  const bgTypes=[cur.type].concat((o.bgs||["smart","shirley","active"]).filter(t=>t!==cur.type));
  const nsStatic=[...new Set([curN].concat(o.ns||[3,10]))];
  const bgLv=[]; bgTypes.forEach(t=>{ if(t==="active"||t==="none") bgLv.push({t,n:null}); else (t===cur.type?nsStatic:[curN]).forEach(n=>bgLv.push({t,n})); });
  // 끝점 점 수는 지금 배경 종류에서만 바꾼다(다른 배경은 지금 점 수) — 요인이 섞이지 않게
  const loLv=[0,-dw,dw], hiLv=[0,-dw,dw];
  // 창이 스펙트럼 밖으로 나가면 스펙트럼 끝으로 자르고, 자른 결과가 이미 있는 변형과 같으면 넣지 않는다(기준이 두 번 세어지지 않게)
  const r3=x=>Math.round(x*1000)/1000;
  const mk=(b,dl,dh)=>{ const lo=Math.max(E0,model.lo+dl), hi=Math.min(E1,model.hi+dh); return {lo,hi,bg:sensBg(cur,b.t,b.n),f:{bg:b.t,n:b.n,lo:r3(lo-model.lo),hi:r3(hi-model.hi),nRef:curN}}; };
  const out=[]; const key=v=>[v.f.bg,v.f.n,r3(v.lo),r3(v.hi)].join("|"); const seen=new Set(); const push=v=>{ const k=key(v); if(seen.has(k)) return; seen.add(k); out.push(v); };
  push(mk(bgLv[0],0,0));
  if(o.mode==="full"){ bgLv.forEach(b=>loLv.forEach(dl=>hiLv.forEach(dh=>push(mk(b,dl,dh))))); }
  else { bgLv.forEach(b=>push(mk(b,0,0))); loLv.forEach(dl=>push(mk(bgLv[0],dl,0))); hiLv.forEach(dh=>push(mk(bgLv[0],0,dh))); }
  return out;
}
/* 변형들을 실제로 맞춘다(Worker 안에서). spec={be,y,c}, model={lo,hi,bg,comps,cons}.
   각 변형은 지금 성분 값에서 출발(넓이는 NNLS 로 다시). 성분 중심이 창 밖이면 건너뜀. */
function sensJob(spec,model,variants){
  return variants.map((v,j)=>{ const i=v.idx!=null?v.idx:j;   // 여러 Worker 로 나눌 때 원래 번호
    try{
      // 지금 창 안에 있던 성분 중심이 변형 창에서 잘려 나가면 건너뜀(창 밖에 있던 성분 — 꼬리·가장자리 — 은 사용자 모형 그대로 둔다)
      if(model.comps.some(c=>c.be>=model.lo&&c.be<=model.hi&&(c.be<v.lo||c.be>v.hi))) return {i,ok:false,err:"성분 중심이 창 밖으로 잘림"};
      const P=makeProblem(spec,{lo:v.lo,hi:v.hi,bg:v.bg,comps:model.comps,cons:model.cons});
      const r=fitOnce(P,{}); const st=fitStats(P,r); const tot=P.comps.reduce((s,c)=>s+c.area,0);
      return {i,ok:true,conv:r.conv,chi2nu:st.chi2nu,comps:P.comps.map(c=>{ const tf=trueFWHM(c.shape,shapeParams(c,1)); return {id:c.id,be:c.be,area:c.area,areaPct:tot>0?100*c.area/tot:null,fwhm:tf.fwhm,dE:c.doublet?c.doublet.dE:null}; }),tot};
    }catch(e){ return {i,ok:false,err:String(e&&e.message||e)}; }
  });
}
/* 요약: 성분마다 각 값의 기준·범위·평균·sd(변형들 같은 무게), 요인별 주효과(수준별 평균의 범위). 기준 = 0 번 변형. */
function sensSummary(variants,res){
  const ok=res.filter(r=>r&&r.ok); if(!ok.length) return null; const base=res.find(r=>r&&r.i===0&&r.ok)||null;
  const ids=(base||ok[0]).comps.map(c=>c.id); const Q=["be","areaPct","area","fwhm"];
  const st=a=>{ const n=a.length; const m=a.reduce((s,x)=>s+x,0)/n; const sd=n>1?Math.sqrt(a.reduce((s,x)=>s+(x-m)*(x-m),0)/(n-1)):0; return {n,mean:m,sd,min:Math.min(...a),max:Math.max(...a)}; };
  const comps=ids.map(id=>{ const o={id};
    Q.forEach(q=>{ const vals=ok.map(r=>{ const c=r.comps.find(x=>x.id===id); return c&&c[q]!=null&&isFinite(c[q])?c[q]:null; }).filter(v=>v!=null);
      if(!vals.length) return; const s=st(vals); const bc=base&&base.comps.find(x=>x.id===id); s.base=bc?bc[q]:null;
      // 주효과: 요인마다 수준별 평균 → 그 범위
      // 요인 fk 의 효과는 '다른 요인이 기준 수준인 변형들'끼리만 비교(배경 종류가 바뀌면 끝점 점 수도 함께 바뀌는 혼동을 피함:
      // 배경 효과는 끝점 점 수를 따지지 않고, 끝점 효과는 같은 배경 종류 안에서만)
      const b0=(variants[0]||{}).f||{};
      // 배경 효과: 끝점 점 수가 '다른 배경에 준 값(nRef)' 이거나 없는(피팅 안 배경) 변형만 — 지금 배경의 n=3·10 변형이 섞이지 않게
      const same=(f,fk)=>["bg","n","lo","hi"].every(k=>k===fk||(fk==="bg"&&k==="n"&&(f.n==null||f.n===f.nRef))||f[k]===b0[k]);
      const eff={}; ["bg","n","lo","hi"].forEach(fk=>{ const g={}; ok.forEach(r=>{ const v=variants[r.i]; const c=r.comps.find(x=>x.id===id); if(!v||!c||c[q]==null||!same(v.f,fk)) return; const lv=String(v.f[fk]); (g[lv]=g[lv]||[]).push(c[q]); });
        const ms=Object.values(g).map(a=>a.reduce((s,x)=>s+x,0)/a.length); eff[fk]=ms.length>1?Math.max(...ms)-Math.min(...ms):0; });
      s.eff=eff; o[q]=s; });
    return o; });
  return {n:variants.length,nOk:ok.length,nFail:variants.length-ok.length,base:!!base,comps,chi2nu:st(ok.map(r=>r.chi2nu))};
}

/* ============================================================
   보조: 문제 만들기 (스펙트럼 + 영역 + 배경 + 성분)
   ============================================================ */
function sliceRegion(be,y,lo,hi){ const i0=be.findIndex(e=>e>=lo); let i1=-1; for(let i=be.length-1;i>=0;i--) if(be[i]<=hi){ i1=i; break; } if(i0<0||i1<i0) return null; return {E:be.slice(i0,i1+1),y:y.slice(i0,i1+1),i0,i1}; }
function makeProblem(spec,st){
  const r=sliceRegion(spec.be,spec.y,st.lo,st.hi); if(!r||r.E.length<8) throw new Error("영역 안 데이터 점이 8개 미만입니다");
  return {E:Float64Array.from(r.E),y:Float64Array.from(r.y),c:spec.c||null,bg:clone(st.bg),comps:clone(st.comps),cons:clone(st.cons||{}),i0:r.i0};
}
function newComp(shape,be,o){
  o=o||{}; const c={id:o.id||("c"+Math.random().toString(36).slice(2,8)),name:o.name||"",shape,be,area:o.area||1,el:o.el||"",line:o.line||"",doublet:null,color:o.color||null};
  const F=o.F||1.2;
  if(shape==="GLS"||shape==="GLP"){ c.F=F; c.m=o.m!=null?o.m:0.3; }
  if(shape==="V"){ c.wG=o.wG!=null?o.wG:F*0.85; c.wL=o.wL!=null?o.wL:F*0.25; }
  if(shape==="DS"){ c.alpha=o.alpha!=null?o.alpha:0.1; c.wL=o.wL!=null?o.wL:0.3; c.wG=o.wG!=null?o.wG:F*0.7; }
  if(shape==="LA"){ c.F=o.F||F*0.8; c.alpha=o.alpha!=null?o.alpha:1; c.beta=o.beta!=null?o.beta:1.5; c.wG=o.wG!=null?o.wG:F*0.4; }
  if(shape==="GLT"){ c.F=F; c.m=o.m!=null?o.m:0.3; c.T=o.T!=null?o.T:0.2; c.tau=o.tau!=null?o.tau:1.0; }
  if(o.doublet) c.doublet={dE:o.doublet.dE,r:o.doublet.r,q:o.doublet.q!=null?o.doublet.q:1};
  return c;
}
/* 모양 바꾸기: 현재 폭을 최대한 보존 */
function convertShape(c,shape){
  if(c.shape===shape) return c;
  const tf=trueFWHM(c.shape,shapeParams(c,1)).fwhm||1.2;
  const n=newComp(shape,c.be,{F:tf,id:c.id,name:c.name,el:c.el,line:c.line,color:c.color}); n.area=c.area; if(c.doublet) n.doublet=Object.assign({},c.doublet);
  if((shape==="GLS"||shape==="GLP"||shape==="GLT")&&c.m!=null) n.m=c.m;
  return n;
}

/* 넓이 합(선택 성분)의 σ — 공분산과 연결을 반영 */
function areaSumSigma(prob,map,cov,idxs){
  const p0=getVec(prob,map), np=p0.length; const g=new Float64Array(np);
  const sum=()=>idxs.reduce((s,i)=>s+prob.comps[i].area,0);
  for(let j=0;j<np;j++){ if(cov.atB[j]||!(cov.sig[j]>0)) continue; const h=Math.max(cov.sig[j]*1e-3,1e-12*(Math.abs(p0[j])+1)); const q=p0.slice(); q[j]+=h; setVec(prob,map,q); const a=sum(); q[j]=p0[j]-h; setVec(prob,map,q); const b=sum(); setVec(prob,map,p0); g[j]=(a-b)/(2*h); }
  let s=0; for(let a=0;a<np;a++) if(g[a]) for(let b=0;b<np;b++) if(g[b]){ const v=cov.C[a*np+b]; if(isFinite(v)) s+=g[a]*v*g[b]; }
  return Math.sqrt(Math.max(0,s));
}
/* 피팅 결과 → 직렬화 가능한 요약 (Worker 와 주 스레드 공용) */
function finalize(P,r){
  const cov=covariance(P,r); const st=fitStats(P,r); const summ=compSummary(P,r,cov);
  const map=r.map; const E=P.E, n=E.length, ctx={emin:E[0],emax:E[n-1]};
  return {comps:P.comps,bg:P.bg,chi2:r.chi2,iter:r.iter,conv:r.conv,mu:Array.from(r.mu),B:Array.from(r.B),v:Array.from(r.v),
    stats:Object.assign({},st,{resid:Array.from(st.resid)}),
    cov:{sig:cov.sig,corr:Array.from(cov.corr),C:Array.from(cov.C),atB:cov.atB,phi:cov.phi,dof:cov.dof,singular:cov.singular,names:cov.names,np:cov.np},
    summ,multi:r.multi||null,free:map.free.map(f=>({id:f.id,k:f.k,lo:f.lo,hi:f.hi})),
    sumSig:{all:areaSumSigma(P,map,cov,P.comps.map((c,i)=>i))}};
}
/* 한 묶음의 시작점에서 피팅 → 최선 해 (Worker 작업 단위) */
function fitJob(prob,starts,opt){
  let best=null; const chis=[];
  (starts&&starts.length?starts:[null]).forEach(pt=>{ const P=clone(stripProb(prob)); P.E=prob.E; P.y=prob.y; applyStart(P,pt); try{ const r=fitOnce(P,opt); chis.push(r.chi2p1); if(!best||r.chi2p1<best.r.chi2p1) best={P,r}; }catch(e){ chis.push(Infinity); } });
  if(!best) throw new Error("모든 시작점에서 피팅 실패");
  return {comps:best.P.comps,bg:best.P.bg,chi2:best.r.chi2p1,chis};   // chi2 = 1차 χ² (시작점 비교용)
}
/* 최종: 주어진 매개변수에서 (넓이 재초기화 없이) 다시 다듬고 요약 */
function finalJob(prob,opt){
  const P=clone(stripProb(prob)); P.E=prob.E; P.y=prob.y;
  const r=fitOnce(P,Object.assign({},opt,{initAreas:false}));
  return finalize(P,r);
}
/* MC 작업 단위: prob(최적 상태), mu, v, 잡음 매개변수 */
function mcJob(prob,res,M,seed,opt){
  const P=clone(stripProb(prob)); P.E=prob.E; P.y=prob.y;
  const r={mu:Float64Array.from(res.mu),v:Float64Array.from(res.v),chi2:res.chi2,map:buildMap(P)};
  return monteCarlo(P,r,M,Object.assign({},opt,{seed}));
}

const XPS_CORE={areaSumSigma,finalize,fitJob,finalJob,mcJob,SHAPES,PDEF,BG_TYPES,faddeeva,erfcx,lgamma,gauss,lorentz,voigt,glpNorm,laNorm,shapeUnit,trueFWHM,compEval,shapeParams,
  cumtrapz,endpoints,smoothQ,smoothHW,bgLinear,bgShirley,bgTougaard,bgStatic,
  paramNames,getP,setP,defBounds,buildMap,applyLinks,getVec,setVec,toU,fromU,makeEvaluator,
  nnls,rng,fitOnce,fitMulti,startPoints,applyStart,lmFit,covariance,fitStats,compSummary,monteCarlo,mcSummarize,noiseEstimate,
  tfValue,tpp2m,imfpValue,quantify,parseFormula,stoichCompare,augerParam,peakMaxPos,findPeaks,identify,acceptDefault,elementLines,EXCLUDE_DEFAULT,doubletInfo,multipletInfo,dEConstraint,suggestModel,refineByResidual,autoJob,sensVariants,sensJob,sensSummary,sensBg,SENS_BG_LABEL,removeComp,anchorOf,DOUBLET_RESOLVE,sliceRegion,makeProblem,newComp,convertShape,clone,stripProb};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_CORE; else root.XPS_CORE=XPS_CORE;
})(typeof self!=="undefined"?self:this);
