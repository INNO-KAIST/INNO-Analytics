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
  GLS:{name:"GL sum (pseudo-Voigt, area-weighted)",p:["F","m"]},
  GLP:{name:"GL product (Avantage Product · Casa GL)",p:["F","m"]},
  V:{name:"Voigt (G⊗L convolution)",p:["wG","wL"]},
  DS:{name:"Doniach–Šunjić ⊗ G (asymmetric, metals)",p:["alpha","wL","wG"]},
  LA:{name:"LA(α,β) ⊗ G (CasaXPS LA style)",p:["F","alpha","beta","wG"]},
  GLT:{name:"GL sum + exponential tail",p:["F","m","T","tau"]},
  MEAS:{name:"Measured line shape (reference spectrum)",p:["wG"]},
  AGT:{name:"GL product + Avantage tail (ET, TM, CT)",p:["F","m","ET","TM","CT"]},
  TLA:{name:"TLA — asym. Lorentzian ⊗ G (Casa)",p:["F","alpha","mu","wG"]},   // SPEC_shapes2 §1 — CasaXPS Cookbook 2026 pp. 65–68
  GLN:{name:"GL(m,n) — GL product ⊗ G, tail T/PHI (Casa)",p:["F","m","wG"]},
  SGLN:{name:"SGL(m,n) — GL height sum ⊗ G, tail T/PHI (Casa)",p:["F","m","wG"]},
  TP:{name:"TP — trapezium ⊗ L, background term (Casa)",p:["F","ta","tb","wL"]},   // SPEC §1.10 — Avantage tail function   // SPEC §12 — 모양 = 측정 template(c.tpl), wG = 추가 Gaussian 넓힘
};
const PDEF={   // 기본 한계·이름·단위
  be:{lab:"BE",u:"eV"}, area:{lab:"Area",u:"CPS·eV",lo:0,hi:Infinity},
  F:{lab:"FWHM",u:"eV",lo:0.05,hi:10}, m:{lab:"L frac.",u:"",lo:0,hi:1},
  wG:{lab:"G FWHM",u:"eV",lo:0,hi:8}, wL:{lab:"L FWHM",u:"eV",lo:0.005,hi:8},
  alpha:{lab:"α",u:"",lo:0,hi:0.45}, beta:{lab:"β",u:"",lo:0.55,hi:20},
  T:{lab:"Tail frac.",u:"",lo:0,hi:1}, tau:{lab:"Tail len.",u:"eV",lo:0.02,hi:20},
  ET:{lab:"ET (exp. tail)",u:"/ch",lo:0.001,hi:1}, TM:{lab:"TM (tail mix)",u:"",lo:0,hi:1}, CT:{lab:"CT (const. tail)",u:"",lo:0,hi:1},
  dE:{lab:"Δ(SO)",u:"eV",lo:0.05,hi:40}, r:{lab:"Area ratio",u:"",lo:0.01,hi:5}, q:{lab:"Width fac.",u:"",lo:0.3,hi:3},
  mu:{lab:"μ (asym.)",u:"",lo:0.05,hi:1000}, stMu:{lab:"ST μ",u:"",lo:0.02,hi:1000}, stG:{lab:"ST G FWHM",u:"eV",lo:0,hi:8},
  tk:{lab:"Tail k",u:"",lo:0.05,hi:100}, tj:{lab:"Tail j",u:"",lo:0,hi:20}, ta:{lab:"TP α (hi BE)",u:"",lo:0,hi:1000}, tb:{lab:"TP β (lo BE)",u:"",lo:0,hi:1000},
  c0:{lab:"Bkg c₀",u:"CPS",lo:-Infinity,hi:Infinity}, k:{lab:"Shirley k",u:"1/eV",lo:0,hi:Infinity}, s:{lab:"Bkg slope",u:"CPS/eV",lo:-Infinity,hi:Infinity},
};
const LA_ALPHA_LO=0.55;
/* CasaXPS LA(α,β,n): Gaussian FWHM = n·CASA_LA_RHO·(Lorentzian F) — CasaXPS 2.3.27 LA(1,1,n) n-scan on V2_Voigt_wG0.6_wL0.4 (n = 100…1000,
   fwhm 상대 ≤ 1.6e-4 재현; SPEC_casa_v2 §2). 1/144.87 ± 0.3 %. α = β = 1 에서만 측정됨 */
const CASA_LA_RHO=1/144.87;
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
/* 계산량 세기(결과와 무관, 결정적): 평가 점 수 × 모양별 대략의 점당 비용 + 격자 점 수 — §3.3b 의 현직·갈래를 같은 계산량씩 번갈아 돌리려고 */
const WORK={n:0,kc:0.03}, WORK_PT={GLP:0.5,V:2,LA:12,GLT:4,AGT:3,TP:1.3};
function shapeUnit(shape,x,p,ctx){
  const n=x.length, out=new Float64Array(n); WORK.n+=n*(shape==="LA"&&!(p.wG>1e-6)?1:WORK_PT[shape]||1);
  if(shape==="MEAS") return tplShape(p.tpl,x,p.wG||0);
  if(p.mod||shape==="TLA"||shape==="GLN"||shape==="SGLN"||shape==="TP"||(shape==="DS"&&(p.dsl>0||p.wG>1e-6))) return s2Eval(shape,x,p,ctx);   // SPEC_shapes2
  if(shape==="AGT") return aglShape(x,p,ctx);
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
  if(p.mod||shape==="TLA"||shape==="GLN"||shape==="SGLN"||shape==="TP"||shape==="DS") return s2FWHM(shape,p);   // DS: 정밀 교차점 (v1.9, 이전 선형 보간 오차 ~3e-5)
  const w=shape==="V"?(0.5346*p.wL+Math.sqrt(0.2166*p.wL*p.wL+p.wG*p.wG)):shape==="DS"?Math.max(p.wL,0.01)+(p.wG||0):shape==="LA"?p.F+(p.wG||0):p.F;
  let R=6*w+(shape==="GLT"?6*p.tau:0), N=2401;
  if(shape==="MEAS"){ const tc=tplCache(p.tpl); R=Math.max(Math.abs(tc.x[0]),Math.abs(tc.x[tc.x.length-1]))+5*(p.wG||0)*FW2SIG+0.05; N=4001; } const x=new Float64Array(N); for(let i=0;i<N;i++) x[i]=-R+2*R*i/(N-1);
  const v=shapeUnit(shape,x,p,{lo:-R,hi:R});
  let im=0; for(let i=1;i<N;i++) if(v[i]>v[im]) im=i;
  const half=v[im]/2; let a=null,b=null;
  for(let i=im;i>0;i--) if(v[i-1]<half){ a=x[i-1]+(half-v[i-1])/(v[i]-v[i-1])*(x[i]-x[i-1]); break; }
  for(let i=im;i<N-1;i++) if(v[i+1]<half){ b=x[i]+(v[i]-half)/(v[i]-v[i+1])*(x[i+1]-x[i]); break; }
  return {fwhm:a!=null&&b!=null?b-a:null,xmax:x[im],ymax:v[im],lo:a,hi:b};   // lo/hi: 반높이 점 (낮은/높은 BE 쪽, 원점 기준) — casaWidthAsym
}

/* ============================================================
   2) 성분 → 세기
   comp = {id, shape, be, area, F,m,wG,wL,alpha,beta,T,tau, doublet:{dE,r,q}|null}
   ============================================================ */
/* 새 모양·수식어 (SPEC_shapes2): ST·SB 를 붙일 수 있는 바탕, 꼬리 T/PHI 는 GLN·SGLN 에만, DS 의 Shirley 자름 l = c.dsl */
const S2_MODBASE={GLS:1,GLP:1,LA:1,V:1,TLA:1,GLN:1,SGLN:1};
function compMod(c){ return (c.mod==="ST"||c.mod==="SB")&&S2_MODBASE[c.shape]?c.mod:null; }
function compTail(c){ return (c.shape==="GLN"||c.shape==="SGLN")&&(c.tail==="T"||c.tail==="PHI")?c.tail:null; }
function shapePList(c){ const a=SHAPES[c.shape].p.slice(); const t=compTail(c); if(t==="T") a.push("tk"); if(t==="PHI") a.push("tj","tk"); if(compMod(c)==="ST") a.push("stMu","stG"); return a; }
/* CasaXPS 정수 n → 폭 (SPEC_casa_v2 §2.3 의 LA 규칙 FWHM = n·F/144.87 을 다른 모양에도 '가정'): 어느 매개변수가 n 으로 정해지는가 */
function casaNKey(c){ return c.casaN==null?null:c.shape==="TP"?"wL":(c.shape==="LA"||c.shape==="TLA"||c.shape==="GLN"||c.shape==="SGLN"||c.shape==="DS")?"wG":null; }
function derivedParam(c,k){ return (k===casaNKey(c))||(k==="stG"&&c.stN!=null&&compMod(c)==="ST"); }
function defaultCons(c,k){ if(c.shape==="TP"&&(k==="ta"||k==="wL")) return {mode:"fixed"}; if(k==="stG"&&compMod(c)==="ST"&&!(c.stG>0)) return {mode:"fixed"}; return null; }   // TP: 넓이와 겹쳐 α/β 비만 정해짐, wL 은 풀면 0 으로 가며 피팅이 늘어짐(SPEC_shapes2 §3); ST(μ) = γ 0
function shapeParams(c,scale){
  const q=scale||1, p={};
  shapePList(c).forEach(k=>p[k]=c[k]);
  if(q!==1){ ["F","wL"].forEach(k=>{ if(p[k]!=null) p[k]*=q; }); }   // 수명 폭만 (SPEC 1.9 v1.1) — wG 는 공통
  if(c.shape==="LA"){ p.alpha=Math.max(LA_ALPHA_LO,p.alpha); p.beta=Math.max(LA_ALPHA_LO,p.beta); if(c.casaN!=null) p.wG=Math.max(0,+c.casaN)*CASA_LA_RHO*p.F; }
  if(c.casaN!=null&&c.shape!=="LA"){ const k=casaNKey(c); if(k) p[k]=Math.max(0,+c.casaN)*CASA_LA_RHO*(c.shape==="DS"?p.wL:p.F); }
  const md=compMod(c); if(md){ p.mod=md; if(md==="ST"&&c.stN!=null) p.stG=Math.max(0,+c.stN)*CASA_LA_RHO*(c.shape==="V"?p.wL:p.F); if(md==="ST"&&c.stCasaK>0) p.stLag=c.stCasaK*(c.shape==="V"?p.wL:p.F)/80; if(md==="SB") p.sbArg=c.sbArg||0; }
  const tl=compTail(c); if(tl) p.tail=tl;
  if(c.shape==="DS"&&c.dsl>0) p.dsl=+c.dsl;
  if(c.shape==="MEAS") p.tpl=c.tpl;
  if(c.shape==="AGT"){ p.ch=c.ch; if(c.tu) p.tu=c.tu; }
  return p;
}
/* E: 절대 BE 배열, ctx={emin,emax} → Float64Array */
/* 넓이 1 모양 캐시 (v1.9): 넓이만 바뀌는 평가(LM 의 넓이 열, 넓이 초기화)는 모양을 다시 계산하지 않는다. 비싼(격자) 모양만, 성분마다 최근 4 개.
   키 = 모양 매개변수 + 위치 + 창 + 같은 E 배열 → 값은 그대로(비트 단위 같음) */
const UCACHE=typeof WeakMap!=="undefined"?new WeakMap():null, EID=typeof WeakMap!=="undefined"?new WeakMap():null; let EIDN=0;
function costlyShape(shape,p){ return shape==="TLA"||shape==="GLN"||shape==="SGLN"||shape==="TP"||!!p.mod||shape==="DS"||(shape==="LA"&&p.wG>1e-6)||shape==="GLT"; }
function unitCached(c,E,e0,p,ctx){ const n=E.length;
  const calc=()=>{ const x=new Float64Array(n); for(let i=0;i<n;i++) x[i]=E[i]-e0; return shapeUnit(c.shape,x,p,ctx?{lo:ctx.emin-e0,hi:ctx.emax-e0}:null); };
  if(!UCACHE||!costlyShape(c.shape,p)||typeof E!=="object") return calc();
  let id=EID.get(E); if(id==null){ id=++EIDN; EID.set(E,id); }
  const key=c.shape+"|"+id+"|"+e0+"|"+(ctx?ctx.emin+","+ctx.emax:"")+"|"+JSON.stringify(p);
  let L=UCACHE.get(c); if(!L){ L=[]; UCACHE.set(c,L); } for(const q of L) if(q.k===key) return q.u;
  const u=calc(); L.unshift({k:key,u}); if(L.length>4) L.pop(); return u; }
function compEval(c,E,ctx,parts){
  const n=E.length, out=new Float64Array(n);
  const one=(e0,A,p)=>{ const u=unitCached(c,E,e0,p,ctx); const r=new Float64Array(n); for(let i=0;i<n;i++){ r[i]=A*u[i]; out[i]+=r[i]; } return r; };
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
  const T0=o.T0>0?+o.T0:0;   // CasaXPS U(x; B,C,D,T0): 0 for x <= T0 (Fairley et al. 2023 SI Eq. 9) — 기본 0 은 전과 같음
  const K=T=>T<=T0?0:(U3?T/((C-T*T)*(C-T*T)+D*T*T):T/((C+T*T)*(C+T*T)));
  const I=new Float64Array(N);
  for(let i=1;i<N;i++){ let s=0; for(let j=0;j<i;j++){ const fa=K(E[i]-E[j])*(y[j]-p.ylo), fb=K(E[i]-E[j+1])*(y[j+1]-p.ylo); s+=0.5*(fa+fb)*(E[j+1]-E[j]); } I[i]=s; }
  // I(E_hi): 선형 보간
  let Ih=I[N-1]; for(let i=1;i<N;i++) if(E[i]>=p.Ehi){ const f=(p.Ehi-E[i-1])/(E[i]-E[i-1]); Ih=I[i-1]+f*(I[i]-I[i-1]); break; }
  let s=Math.abs(Ih)>1e-300?(p.yhi-p.ylo)/Ih:0; const neg=s<0; if(neg) s=0;   // 음의 손실 단면적 금지 (SPEC 2 v1.1)
  const B=Float64Array.from(I,v=>p.ylo+s*v); B.scale=s; B.negClamped=neg; return B;
}
/* CasaXPS 호환 영역 규약 (SPEC §2b, CasaXPS 2.3.27 실측으로 정함): 끝점 = 창의 두 번째·끝에서 두 번째 채널(i0+1, i1−1)을 중심으로 ±aw 채널 평균,
   직선은 그 두 채널 에너지에 고정. Shirley 누적은 배타적 직사각 합 Q_i = Σ_{j<i}(y_j − B_j) (j = i0+1 부터). */
function endpointsCasa(E,y,aw,ends){
  const N=E.length, w=Math.max(0,Math.round(aw||0)), a=Math.min(1,N-1), b=Math.max(0,N-2);
  if(ends&&isFinite(ends.ylo)&&isFinite(ends.yhi)) return {ylo:ends.ylo,yhi:ends.yhi,Elo:E[a],Ehi:E[b],ia:a,ib:b,n:2*w+1};   // 창 밖 채널까지 평균(전체 스펙트럼에서 미리 계산)
  const av=c=>{ let s=0,k=0; for(let j=Math.max(0,c-w);j<=Math.min(N-1,c+w);j++){ s+=y[j]; k++; } return s/k; };
  return {ylo:av(a),yhi:av(b),Elo:E[a],Ehi:E[b],ia:a,ib:b,n:2*w+1};
}
function bgShirleyCasa(E,y,p){
  const N=E.length, a=p.ia, b=p.ib, d=p.yhi-p.ylo; let B=new Float64Array(N).fill(p.ylo); for(let i=b;i<N;i++) B[i]=p.yhi;
  const tol=1e-7*(Math.abs(d)+1e-12); let it=0, conv=false;
  for(;it<200;it++){ const nb=new Float64Array(N); let Q=0; const Qs=new Float64Array(N);
    for(let i=a;i<=b;i++){ Qs[i]=Q; Q+=y[i]-B[i]; }
    const Qt=Qs[b]; let mx=0;
    for(let i=0;i<N;i++){ const q=i<=a?0:i>=b?Qt:Qs[i]; nb[i]=Math.abs(Qt)>1e-300?p.ylo+d*q/Qt:p.ylo; mx=Math.max(mx,Math.abs(nb[i]-B[i])); }
    B=nb; if(mx<tol){ conv=true; break; } }
  B.iter=it+1; B.converged=conv; return B;
}
/* Avantage 규약 Shirley v1.8 (SPEC §2c; Thermo Avantage Data System manual, Shirley Background): 반복 수 상한 it(기본 10, 1 = 단순 한 번),
   넓이 변화가 cl(상대) 이하이면 멈춤; 높은 BE(낮은 KE) 끝이 낮은 BE 끝보다 낮거나 같으면 Linear 로 바꿈. 끝점·적분은 INNO §2 그대로.
   v1.9.2: avantageV19(bg) 가 거짓이면(기본) 이것을 쓴다 — 아래 v1.9.2 규약은 bg.avgEV · bg.ah · bg.rule "v19" 로 켠다. */
function bgShirleyAvantageV18(E,y,bg){
  const N=E.length, p=endpoints(E,y,bg.n||5), d=p.yhi-p.ylo, itMax=Math.max(1,Math.round(bg.it!=null?bg.it:10)), cl=bg.cl!=null?+bg.cl:0;
  if(!(d>0)){ const B=bgLinear(E,y,bg.n||5); B.fallback="linear"; return B; }
  let B=new Float64Array(N).fill(p.ylo), Aold=null, it=0, conv=false;
  for(;it<itMax;it++){ const Q=cumtrapz(E,Float64Array.from(y,(v,i)=>v-B[i])), Qt=Q[N-1]; B=Float64Array.from(Q,q=>Math.abs(Qt)>1e-300?p.ylo+d*q/Qt:p.ylo);
    const A=cumtrapz(E,Float64Array.from(y,(v,i)=>v-B[i]))[N-1]; if(Aold!=null&&Math.abs(A-Aold)<=cl*Math.abs(A)){ conv=true; it++; break; } Aold=A; }
  B.iter=it; B.converged=conv||it>=itMax; return B;
}
/* ---- Avantage 규약 v1.9.2 (test/ref/avantage_a2_patch.md — K-Alpha 실측 2 시료, Smart 12 · Shirley 12 영역의 'Backgnd' 곡선에서 역으로 정함; Avantage 코드는 쓰지 않음) ----
   h  = 'Background Average at Start & End' 반폭(채널) = round(avgEV/(2·step)); avgEV 가 없으면 bg.ah, 그것도 없으면 4 (v1.8 Kratos 자료의 9 점).
        a2: 0.50 eV · 0.05 eV → h = 5 (11 점).
   ①  끝값 = 끝 채널 ±h 평균, 창 밖 채널 포함(bg._pad), 스펙트럼 끝에서 자름 (24/24 영역, 인쇄 자리수 안).
   ②  평활 S: ±hs 이동평균, 창 안 채널만 쓰고 낮은 BE 끝 채널은 뺀다 (Smart hs = h, Shirley hs = 2h).
        Smart 의 접점(배경 = 데이터) 값이 이 S 와 24 영역 모두 인쇄 자리수 안에서 같다.
   ③  모양: 고정점 u = p·u₁ + (1−p)·G(q·u₁ + (1−q)·u), G(v) = Δ·Q[max(0, s − v)]/Q_N (포함 누적합, 낮은 BE 끝에서),
        s = S − 낮은 끝값, Δ = 높은 끝값 − 낮은 끝값, u₁ = G(0). p = 0.07, q = 0.25 는 실측 맞춤 상수 (bg.avP, bg.avQ 로 바꿀 수 있음).
        '반복들의 평균'으로 계산하는 Shirley 의 근사로 보인다(첫 반복 무게 ≈ 1/14). Avantage 의 정확한 반복 규칙은 미확인.
   ④  높은 BE 끝값이 낮은 BE 끝값보다 크지 않으면 직선(설명서).                                                              */
/* v1.9.2 규약을 쓰는 조건: bg.rule === "v19", 또는 bg.avgEV / bg.ah 가 주어짐. 그 밖(기존 프로젝트·UI 의 it/cl, SPEC §2c 기준 구현 시험)은 v1.8 그대로; bg.rule === "v18" 은 항상 v1.8. */
function avantageV19(bg){ return !!bg&&bg.rule!=="v18"&&(bg.rule==="v19"||bg.avgEV>0||(bg.ah!=null&&isFinite(bg.ah))); }
function avantageHalf(bg,E){
  if(bg&&bg.ah!=null&&isFinite(bg.ah)) return Math.max(0,Math.round(+bg.ah));
  if(bg&&bg.avgEV>0&&E&&E.length>1){ const st=Math.abs(E[E.length-1]-E[0])/(E.length-1); if(st>0) return Math.max(0,Math.round(bg.avgEV/(2*st)+1e-9)); }
  return 4;
}
/* 오름차순(낮은 BE → 높은 BE) 배열로 돌려 계산하고 원래 순서로 되돌린다. pad: {lo:[창 앞 채널], hi:[창 뒤 채널]} (배열 순서 기준) */
function avantageOriented(E,y,bg,fn){
  const N=E.length, desc=N>1&&E[0]>E[N-1];
  const pl=bg&&bg._pad&&bg._pad.lo?Array.from(bg._pad.lo):[], ph=bg&&bg._pad&&bg._pad.hi?Array.from(bg._pad.hi):[];
  if(!desc) return fn(Array.from(y),pl,ph);
  const B=fn(Array.from(y).reverse(),ph.slice().reverse(),pl.slice().reverse()); const R=Float64Array.from(B).reverse();
  for(const k of ["fallback","iter","converged","above","aboveFrac","maxAbove","half","rule","ends"]) if(B[k]!==undefined) R[k]=B[k];
  if(B.splits) R.splits=B.splits.map(i=>N-1-i).sort((a,b)=>a-b);
  if(B.smoothed) R.smoothed=Float64Array.from(B.smoothed).reverse();
  return R;
}
function avantageEnds(y,pl,ph,h){
  const N=y.length, Y=pl.concat(y,ph), o=pl.length;
  const av=c=>{ let s=0,k=0; for(let j=Math.max(0,c-h);j<=Math.min(Y.length-1,c+h);j++){ s+=Y[j]; k++; } return s/k; };
  return {a:av(o),b:av(o+N-1)};
}
function avantageSmooth(y,hs){   // 창 안, 낮은 BE 끝(0번) 채널 제외
  const N=y.length, S=new Float64Array(N), c=new Float64Array(N+1); for(let i=0;i<N;i++) c[i+1]=c[i]+y[i];
  for(let i=0;i<N;i++){ const lo=Math.max(1,i-hs), hi=Math.min(N-1,i+hs); S[i]=hi>=lo?(c[hi+1]-c[lo])/(hi-lo+1):y[i]; }
  return S;
}
function avantagePQ(s,D,p,q){
  const n=s.length-1, lin=()=>Float64Array.from(s,(_,t)=>D*t/(n||1));
  if(n<1) return lin();
  const G=v=>{ const C=new Float64Array(n+1); for(let t=1;t<=n;t++) C[t]=C[t-1]+Math.max(0,s[t]-v[t]); if(!(C[n]>0)) return lin(); const k=D/C[n]; for(let t=0;t<=n;t++) C[t]*=k; return C; };
  const u1=G(new Float64Array(n+1)), u=Float64Array.from(u1), w=new Float64Array(n+1); let it=0, conv=false;
  for(;it<5000;it++){ for(let t=0;t<=n;t++) w[t]=q*u1[t]+(1-q)*u[t]; const g=G(w); let m=0;
    for(let t=0;t<=n;t++){ const nu=0.5*u[t]+0.5*(p*u1[t]+(1-p)*g[t]); m=Math.max(m,Math.abs(nu-u[t])); u[t]=nu; }
    if(m<=1e-10*Math.abs(D)+1e-12){ conv=true; break; } }
  u.iter=it; u.converged=conv; return u;
}
function bgShirleyAvantage(E,y,bg){
  bg=bg||{}; if(!avantageV19(bg)) return bgShirleyAvantageV18(E,y,bg);
  const N=E.length; if(N<3) return bgLinear(E,y,1);
  const h=avantageHalf(bg,E), p=bg.avP!=null?+bg.avP:0.07, q=bg.avQ!=null?+bg.avQ:0.25, hs=bg.avHs!=null?Math.round(bg.avHs):2*h;
  return avantageOriented(E,y,bg,(Y,pl,ph)=>{
    const {a,b}=avantageEnds(Y,pl,ph,h), S=avantageSmooth(Y,hs), D=b-a;
    let B;
    if(!(D>0)){ B=Float64Array.from(Y,(_,t)=>a+D*t/(N-1)); B.fallback="linear"; }
    else { const u=avantagePQ(Array.from(S,v=>v-a),D,p,q); B=Float64Array.from(u,v=>a+v); B.iter=u.iter; B.converged=u.converged; }
    // 경고용: 평활 데이터 위로 지나는 채널 수 · 그 구간 최대 초과
    let above=0, maxAbove=0; for(let t=1;t<N-1;t++){ const d=B[t]-S[t]; if(d>0){ above++; if(d>maxAbove) maxAbove=d; } }
    B.above=above; B.aboveFrac=above/Math.max(1,N-2); B.maxAbove=maxAbove; B.half=h; B.rule="v19"; B.smoothed=S; B.ends=[a,b];
    return B; });
}
/* Avantage 꼬리 함수 성분 (SPEC §1.10; manual, Tail Function — after Sherwood, Briggs & Seah, Practical Surface Analysis, App. 3):
   GL = exp(−4ln2(1−M)x²/F²)/(1+4M x²/F²) (F = 2β), 높은 BE 쪽(x>0) h = GL + (1−GL)·T, T = TM·CT + (1−TM)·exp(−ET·Dx), Dx = x/ch (채널 수).
   낮은 BE 쪽 h = GL. 상수 꼬리는 ℝ 넓이가 무한이라 넓이는 DS 처럼 피팅 창 안 적분(Simpson, [lo,0]·[0,hi] 따로, 간격 ≤ min(F/40, ch/4)) */
function aglRaw(x,p){ const F=p.F,m=p.m,a=4*LN2*(1-m)/(F*F),b=4*m/(F*F),t=x*x,g=Math.exp(-a*t)/(1+b*t); if(x<=0) return g;
  const ch=p.tu==="fwhm"?F:(p.ch>0?p.ch:0.05), T=(p.TM||0)*(p.CT||0)+(1-(p.TM||0))*Math.exp(-(p.ET!=null?p.ET:0.05)*x/ch); return g+(1-g)*T; }   // tu = "fwhm": CasaXPS GL(m)T(k) — 꼬리 길이 단위 = f (SPEC §1.10)
function aglShape(x,p,ctx){ const n=x.length, out=new Float64Array(n); let xmin=Infinity,xmax=-Infinity; for(let i=0;i<n;i++){ if(x[i]<xmin) xmin=x[i]; if(x[i]>xmax) xmax=x[i]; }
  const lo=ctx&&ctx.lo!=null?ctx.lo:xmin, hi=ctx&&ctx.hi!=null?ctx.hi:xmax, h=Math.min(p.F/40,(p.ch>0?p.ch:0.05)/4);
  const simp=(a,b)=>{ if(!(b>a)) return 0; let k=Math.max(2,Math.ceil((b-a)/h)); if(k%2) k++; const hh=(b-a)/k; let s=aglRaw(a,p)+aglRaw(b,p); for(let i=1;i<k;i++) s+=(i%2?4:2)*aglRaw(a+i*hh,p); return s*hh/3; };
  const N=lo<0&&hi>0?simp(lo,0)+simp(0,hi):simp(lo,hi);
  for(let i=0;i<n;i++) out[i]=aglRaw(x[i],p)/(N>0?N:1); return out; }
/* ============================================================
   1b) CasaXPS Cookbook 2026 선모양 (SPEC_shapes2.md; Cookbook pp. 65–68 의 식만 사용 — CasaXPS 코드 없음)
   TLA(α,μ,w) · GL(m,n)/SGL(m,n) [+T(k) | PHI(j,k)] · TP(α,β,n) · DS(α,n[,l]) · 수식어 ST(μ,γ)·SB(0) (GLS·GLP·LA·V·TLA·GLN·SGLN 에)
   좌표: INNO x = E − E0 (x > 0 = 높은 BE). Cookbook 의 x_c = −x/F (FWHM 단위). LA 의 α(x_c ≤ 0), DS 꼬리, E(j,k)(x_c < 0), ∫_{x_c}^∞ 가
   모두 x_c < 0 = 높은 BE 쪽 → 꼬리·Shirley 계단은 높은 BE 쪽 (CasaXPS LA(1.53,2.43) 실측 대조와 같은 방향, SPEC_casa_v2 §2).
   ============================================================ */
const S2_NEG=0.0015;   // Gaussian σ < 0.0015·(최소 폭) 이면 합성곱 생략 (상대 영향 ≤ σ²·|f″|/2f ≈ 1e-5)
const GKX=[0.991455371120812639,0.949107912342758525,0.864864423359769073,0.741531185599394440,0.586087235467691130,0.405845151377397167,0.207784955007898468,0];
const GKW=[0.022935322010529225,0.063092092629978553,0.104790010322250184,0.140653259715525919,0.169004726639267903,0.190350578064785410,0.204432940075298892,0.209482141084727828];
const GGW=[0.129484966168869693,0.279705391489276668,0.381830050505118945,0.417959183673469388];
function gk15(f,a,b){ const c=(a+b)/2, h=(b-a)/2, fc=f(c); let rk=fc*GKW[7], rg=fc*GGW[3];
  for(let j=0;j<7;j++){ const d=h*GKX[j], s=f(c-d)+f(c+d); rk+=GKW[j]*s; if(j&1) rg+=GGW[(j-1)>>1]*s; } return [rk*h,Math.abs((rk-rg)*h)]; }
/* 적응 Gauss–Kronrod (7–15): 오차 합 ≤ tol 이 될 때까지 오차가 가장 큰 구간을 반으로 */
function quadA(f,a,b,tol,brk){
  if(!(b>a)) return 0; tol=tol||1e-13;
  const pts=[a].concat((brk||[]).filter(t=>t>a&&t<b).sort((u,v)=>u-v),[b]); const L=[];
  for(let i=0;i+1<pts.length;i++){ const r=gk15(f,pts[i],pts[i+1]); L.push([pts[i],pts[i+1],r[0],r[1]]); }
  for(let it=0;it<3000;it++){ let e=0,im=0,S=0; for(let i=0;i<L.length;i++){ e+=L[i][3]; S+=L[i][2]; if(L[i][3]>L[im][3]) im=i; }
    if(e<=tol*Math.max(1,Math.abs(S))) break; const q=L[im], m=(q[0]+q[1])/2, r1=gk15(f,q[0],m), r2=gk15(f,m,q[1]); L.splice(im,1,[q[0],m,r1[0],r1[1]],[m,q[1],r2[0],r2[1]]); }
  let S=0; for(const q of L) S+=q[2]; return S; }
function quadInf(f,a,tol){ return quadA(t=>{ const u=1-t; return f(a+t/u)/(u*u); },0,1,tol); }   // ∫_a^∞
const GL4T=[0.0694318442029737124,0.330009478207571868,0.669990521792428132,0.930568155797026288], GL4W=[0.173927422568726929,0.326072577431273071,0.326072577431273071,0.173927422568726929];
function erfcS(z){ return z>=0?erfcx(z)*Math.exp(-z*z):2-erfcx(-z)*Math.exp(-z*z); }
/* TLA 상수: I1 = ∫_0^∞(1+4u²)^−α du = (√π/4)Γ(α−½)/Γ(α);  J = ∫_0^{π/2} φ sin^{2α−2}φ dφ  (SPEC_shapes2 §1) */
const TLAJ=new Map();
function tlaJ(a){ const k=+a.toPrecision(15); let v=TLAJ.get(k); if(v==null){ v=quadA(ph=>ph*Math.pow(Math.sin(ph),2*a-2),0,PI/2,1e-14); if(TLAJ.size>256) TLAJ.clear(); TLAJ.set(k,v); } return v; }
function tlaI1(a){ return (SQPI/4)*Math.exp(lgamma(a-0.5)-lgamma(a)); }
/* GL(m)·SGL(m) 꼬리 정규화: ∫_0^∞ raw(u)·e^{−κu} du */
function glnRaw(form,m){ if(form==="SGL"){ const c=4*LN2; return u=>(1-m)*Math.exp(-c*u*u)+m/(1+4*u*u); } const a=4*LN2*(1-m), b=4*m; return u=>Math.exp(-a*u*u)/(1+b*u*u); }
const GLNC=new Map();
function glnConst(form,m,tail,j,kk){   // → {R0: ∫raw (u), Z: ∫h (u), kap}
  const key=form+"|"+m+"|"+tail+"|"+j+"|"+kk; let v=GLNC.get(key); if(v) return v;
  const R0=form==="SGL"?(1-m)*Math.sqrt(PI/(4*LN2))+m*PI/2:glpNorm(1,m);
  let Z=R0, kap=0, jj=0;
  if(tail==="T"||tail==="PHI"){ kap=tail==="T"?kk:13.8/kk; jj=tail==="T"?1:j; if(kap>0&&jj!==0){ const raw=glnRaw(form,m); const I=quadInf(u=>raw(u)*Math.exp(-kap*u),0,1e-14); Z=R0+jj*(1/kap-I); } }
  v={R0,Z,kap,j:jj}; if(GLNC.size>256) GLNC.clear(); GLNC.set(key,v); return v; }
/* 바탕 모양 (넓이 1; DS 는 정규화 안 함): f(x), left = ∫_{−∞}^0 f, cdf(선택), s0 = 가장 좁은 구조의 폭, sg = 바탕 Gaussian σ (eV) */
function s2Base(sh,p){
  const sgw=w=>w>0?w*FW2SIG:0;
  if(sh==="GLS"){ const F=p.F,m=p.m, z=2*Math.sqrt(LN2)/F; return {f:x=>glsVal(x,F,m),left:0.5,s0:F,sg:0,cdf:x=>(1-m)*0.5*erfcS(-z*x)+m*(0.5+Math.atan(2*x/F)/PI)}; }
  if(sh==="GLP"){ const F=p.F,m=p.m,N=glpNorm(F,m),a=4*LN2*(1-m)/(F*F),b=4*m/(F*F); return {f:x=>Math.exp(-a*x*x)/(1+b*x*x)/N,left:0.5,s0:F,sg:0}; }
  if(sh==="LA"){ const F=p.F,a=Math.max(LA_ALPHA_LO,p.alpha),b=Math.max(LA_ALPHA_LO,p.beta),N=laNorm(F,a,b); return {f:x=>Math.pow(1+4*x*x/(F*F),-(x>=0?a:b))/N,left:(F/4)*SQPI*Math.exp(lgamma(b-0.5)-lgamma(b))/N,s0:F,sg:sgw(p.wG),kink:a!==b}; }
  if(sh==="V"){ const wG=p.wG,wL=p.wL; return {f:x=>voigt(x,wG,wL),left:0.5,s0:Math.max(1e-3,0.5346*wL+Math.sqrt(0.2166*wL*wL+wG*wG)),sg:0}; }
  if(sh==="TLA"){ const F=p.F,a=Math.max(LA_ALPHA_LO,p.alpha),mu=p.mu, I1=tlaI1(a), J=tlaJ(a), c=PI/2+PI/mu, Z=2*c*I1*F;
    return {f:x=>{ const u=x/F; return Math.pow(1+4*u*u,-a)*(PI/2+Math.atan(2*u)+PI/mu)/Z; },left:F*((PI/mu)*I1+0.5*J)/Z,s0:F,sg:sgw(p.wG)}; }
  if(sh==="GLN"||sh==="SGLN"){ const F=p.F,m=p.m,form=sh==="SGLN"?"SGL":"GL", tail=p.tail==="T"||p.tail==="PHI"?p.tail:"", K=glnConst(form,m,tail,p.tj,p.tk), raw=glnRaw(form,m), Z=K.Z*F, kap=K.kap, jj=K.j;
    const f=tail&&kap>0?x=>{ const u=x/F, r=raw(u); return (u>0?r+(1-r)*jj*Math.exp(-kap*u):r)/Z; }:x=>raw(x/F)/Z;
    return {f,left:0.5*K.R0*F/Z,s0:tail&&kap>0?Math.min(F,4*F/kap):F,sg:sgw(p.wG),kink:!!tail}; }
  if(sh==="DS"){ const g=p.wL/2, al=p.alpha, l=p.dsl>0?p.dsl:0;
    if(l>0){ const xc=l*p.wL, kq=dsBase(xc,al,g)/dsPrim(xc,al,g); return {f:x=>x>xc?0:dsBase(x,al,g)-kq*dsPrim(x,al,g),s0:p.wL,sg:sgw(p.wG),xc,win:true,kink:true}; }
    return {f:x=>dsBase(x,al,g),s0:p.wL,sg:sgw(p.wG),win:true,prim:x=>dsPrim(x,al,g)}; }
  throw new Error("s2Base: "+sh); }
/* 6점 Lagrange (마디 k−2…k+3) 보간 · 칸 적분 */
function lag6(v,k0,h,x){ const M=v.length, t=x/h-k0; let b=Math.floor(t)-2; if(b<0) b=0; if(b>M-6) b=M-6; const f=t-b-2;
  const a=f+2,c=f+1,d=f-1,e=f-2,g=f-3;
  return -c*f*d*e*g/120*v[b]+a*f*d*e*g/24*v[b+1]-a*c*d*e*g/12*v[b+2]+a*c*f*e*g/12*v[b+3]-a*c*f*d*g/24*v[b+4]+a*c*f*d*e/120*v[b+5]; }
function gridInt6(v,k0,h,a,b){   // ∫_a^b (6점 보간 함수) — 온 칸은 (11,−93,802,802,−93,11)/1440, 끝 조각은 GL4
  const M=v.length, ta=a/h-k0, tb=b/h-k0; const part=(u,w)=>{ if(!(w>u)) return 0; let s=0; for(let i=0;i<4;i++) s+=GL4W[i]*lag6(v,k0,h,(k0+u+(w-u)*GL4T[i])*h); return s*(w-u)*h; };
  const c0=Math.ceil(ta), c1=Math.floor(tb); if(c1<=c0) return part(ta,tb);
  let s=0; for(let c=c0;c<c1;c++){ if(c-2<0||c+3>M-1){ s+=part(c,c+1); continue; } s+=h*(11*(v[c-2]+v[c+3])-93*(v[c-1]+v[c+2])+802*(v[c]+v[c+1]))/1440; }
  return s+part(ta,c0)+part(c1,tb); }
/* 격자 Gaussian 합성곱 (사다리꼴 = 매끄러운 피적분에 지수 수렴; 핵 ±8σ, 해석 무게) */
function s2Gconv(v,k0,h,sg,r,xa,xb){   // 출력은 r 칸마다 (x = j·r·h), [xa, xb] 안만 → [값, j0].  핵 ±7σ (잘린 무게 < 3e-12)
  r=r||1; const K=Math.ceil(7*sg/h), w=new Float64Array(K+1), c=h/(sg*Math.sqrt(2*PI)); for(let k=0;k<=K;k++){ const t=k*h/sg; w[k]=c*Math.exp(-0.5*t*t); }
  const n=v.length; let j0=Math.ceil(k0/r), j1=Math.floor((k0+n-1)/r);
  if(xa!=null){ j0=Math.max(j0,Math.floor(xa/(r*h))); j1=Math.min(j1,Math.ceil(xb/(r*h))); }
  const m=Math.max(0,j1-j0+1), out=new Float64Array(m); WORK.n+=m*K*WORK.kc;
  for(let jj=0;jj<m;jj++){ const i=(j0+jj)*r-k0; let s=v[i]*w[0];
    if(i>=K&&i+K<n){ for(let k=1;k<=K;k++) s+=(v[i-k]+v[i+k])*w[k]; }
    else { const a=Math.min(K,i), b=Math.min(K,n-1-i); for(let k=1;k<=a;k++) s+=v[i-k]*w[k]; for(let k=1;k<=b;k++) s+=v[i+k]*w[k]; }
    out[jj]=s; }
  return [out,j0]; }
/* TP(α,β,n): 사다리꼴 R (폭 F, 높은 BE 끝 = α, 낮은 BE 끝 = β) ⊗ Lorentzian(wL) — 닫힌 꼴, 넓이 1 = R/((α+β)F/2) */
function tpEval(x,p){ const n=x.length, out=new Float64Array(n), F=p.F, ta=p.ta, tb=p.tb, w=p.wL||0, Z=0.5*(ta+tb)*F; if(!(Z>0)) return out;
  const c0=0.5*(ta+tb), c1=-(tb-ta)/F, a=-F/2, b=F/2;
  if(!(w>S2_NEG*F)){ for(let i=0;i<n;i++){ const t=x[i]; out[i]=t>=a&&t<=b?(c0+c1*t)/Z:0; } return out; }
  for(let i=0;i<n;i++){ const t=x[i], ua=t-a, ub=t-b;
    const dl=(Math.atan(2*ua/w)-Math.atan(2*ub/w))/PI, dm=(w/(4*PI))*Math.log((w*w+4*ua*ua)/(w*w+4*ub*ub));
    out[i]=((c0+c1*t)*dl-c1*dm)/Z; }
  return out; }
/* 새 모양의 값 (상대좌표 x). 넓이 1 (ℝ), 단 DS·SB 는 피팅 창 [ctx.lo, ctx.hi] 안 넓이 1 (SPEC_shapes2 §2) */
function s2Eval(sh,x,p,ctx){
  const n=x.length; let out=new Float64Array(n); if(!n) return out;
  if(sh==="TP") return tpEval(x,p);
  const B=s2Base(sh,p), mod=p.mod==="ST"||p.mod==="SB"?p.mod:null, win=!!B.win||mod==="SB";
  let xmin=Infinity,xmax=-Infinity; for(let i=0;i<n;i++){ if(x[i]<xmin) xmin=x[i]; if(x[i]>xmax) xmax=x[i]; }
  const lo=ctx&&ctx.lo!=null?ctx.lo:xmin, hi=ctx&&ctx.hi!=null?ctx.hi:xmax;
  const neg=S2_NEG*B.s0, sg1=B.sg>neg?B.sg:0, sg2=mod==="ST"&&p.stG>0&&p.stG*FW2SIG>neg?p.stG*FW2SIG:0;
  const muI=mod==="ST"?1/Math.max(1e-9,p.stMu):0, L=mod==="ST"&&p.stLag>0?p.stLag:0; let stN=1/(0.5+muI);   // L: CasaXPS step lag (base FWHM)/80 (SPEC_casa_v4 §4)
  // 1) 합성곱·누적 없음: 해석 값 (창 정규화는 적응 구적)
  if(!sg1&&!mod){ for(let i=0;i<n;i++) out[i]=B.f(x[i]);
    if(win){ const N=B.prim&&!B.xc?B.prim(hi)-B.prim(lo):quadA(B.f,lo,hi,1e-12,[0,B.xc]); const iN=N>0?1/N:1; for(let i=0;i<n;i++) out[i]*=iN; }
    return out; }
  // 2) 격자 (SPEC_shapes2 §3): 입력 간격 h = 0.04·2^−k ≤ min(s0/16 (꺾임 있으면 /32), σ/2) — 양자화(폭이 조금 바뀌어도 마디 고정).
  //    출력(보간·창 적분) 간격 r·h ≤ max(s0/32, σ_eff/8): 넓은 Gaussian 이면 r 칸마다만 합성곱을 계산. 여유 8(σ1+σ2)
  const qz=t=>{ let q=0.04; while(q>t&&q>1e-7) q/=2; return q; };
  const seff=Math.sqrt(sg1*sg1+sg2*sg2), smooth=sg1&&!(mod==="ST"&&!sg1)&&seff>=B.s0/2;
  let h=qz(Math.min(B.s0/(B.kink||!smooth?32:16),sg1?sg1/2:Infinity,sg2?sg2/2:Infinity));
  const P=8*(sg1+sg2)+16*h; let a=Math.min(lo,xmin)-P, b=Math.max(hi,xmax)+P; if(mod){ a=Math.min(a,-h); b=Math.max(b,h); }
  while((b-a)/h>400000) h*=2;
  let r=smooth?Math.max(1,Math.floor(Math.min(seff/8,sg2?sg2/2:Infinity)/h)):1;
  let k0=Math.floor(a/h); const k1=Math.ceil(b/h), M=k1-k0+1, gx=new Float64Array(M); for(let i=0;i<M;i++) gx[i]=(k0+i)*h;
  let v=new Float64Array(M); for(let i=0;i<M;i++) v[i]=B.f(gx[i]); WORK.n+=M*(mod?2:1)*(sh==="V"?2.5:1);
  let Phi=null;
  if(mod){ Phi=new Float64Array(M);
    if(B.cdf){ for(let i=0;i<M;i++) Phi[i]=B.cdf(gx[i]); }
    else { const i0=-k0, cell=(u,w)=>{ let s=0; for(let q=0;q<4;q++) s+=GL4W[q]*B.f(u+(w-u)*GL4T[q]); return s*(w-u); };
      Phi[i0]=B.left; for(let i=i0+1;i<M;i++) Phi[i]=Phi[i-1]+cell(gx[i-1],gx[i]); for(let i=i0-1;i>=0;i--) Phi[i]=Phi[i+1]-cell(gx[i],gx[i+1]); }
    if(L){ const P0=Phi.slice(); let dz=0; for(let i=0;i<M;i++){ Phi[i]=B.cdf?B.cdf(gx[i]-L):lag6(P0,k0,h,gx[i]-L); dz+=v[i]*(Phi[i]-P0[i]); } stN=1/(0.5+muI+dz*h); } }
  const Pf=Phi, kf=k0, hf=h;   // 고운 격자의 Φ (직접 경로용)
  const xa0=Math.min(lo,xmin), xb0=Math.max(hi,xmax);   // 출력이 필요한 범위 (+ 다음 단계 여유, 보간 여유 8 칸)
  if(sg1){ const r1=sg2?1:r, m1=8*sg2+8*h*Math.max(r1,r); let t=s2Gconv(v,k0,h,sg1,r1,xa0-m1,xb0+m1); v=t[0]; if(Phi){ t=s2Gconv(Phi,k0,h,sg1,r1,xa0-m1,xb0+m1); Phi=t[0]; } k0=t[1]; h*=r1; if(r1>1) r=1; }
  let direct=false;
  if(mod==="ST"){ if(!sg1&&!sg2) direct=true; else for(let i=0;i<v.length;i++) v[i]=v[i]*(Phi[i]+muI)*stN; }
  else if(mod==="SB"){ if(!sg1) direct=true; else v=Phi; }
  if(sg2){ const t=s2Gconv(v,k0,h,sg2,r,xa0-8*h*r,xb0+8*h*r); v=t[0]; k0=t[1]; h*=r; }
  if(direct){   // 바탕 값은 해석, Φ 만 보간 (꺾임이 있는 바탕도 정확)
    for(let i=0;i<n;i++){ const ph=B.cdf?B.cdf(x[i]-L):lag6(Pf,kf,hf,x[i]); out[i]=mod==="SB"?ph:B.f(x[i])*(ph+muI)*stN; } }
  else for(let i=0;i<n;i++) out[i]=lag6(v,k0,h,x[i]);
  if(win){ const N=direct?gridInt6(Pf,kf,hf,lo,hi):gridInt6(v,k0,h,lo,hi), iN=N>0?1/N:1; for(let i=0;i<n;i++) out[i]*=iN; }
  return out; }
/* 새 모양의 실제 FWHM: 거친 표본(2401) → 최댓값·반높이 교차를 좁은 구간에서 다시 표본(보간 오차 ≤ 1e-7 F) */
const S2FW=new Map();
function s2FWHM(sh,p){
  if(p.mod==="SB") return {fwhm:null,xmax:0,ymax:0};
  const key=sh+JSON.stringify(p); const hit=S2FW.get(key); if(hit) return hit;
  let w, R;
  if(sh==="TP"){ w=p.F; R=0.6*p.F+8*(p.wL||0)+0.05*p.F; }
  else { const B=s2Base(sh,p); w=B.s0+2.3548*(B.sg+(p.mod==="ST"?(p.stG||0)*FW2SIG:0)); if(sh==="V") w=B.s0; if((sh==="GLN"||sh==="SGLN")&&p.F) w=Math.max(w,p.F); if(sh==="DS") w=Math.max(p.wL,0.01)+2.3548*B.sg; R=8*w; }
  const ev=xs=>shapeUnit(sh,xs,p,{lo:-R,hi:R});
  const N=2401, x=new Float64Array(N); for(let i=0;i<N;i++) x[i]=-R+2*R*i/(N-1); const v=ev(x), dx=2*R/(N-1);
  let im=0; for(let i=1;i<N;i++) if(v[i]>v[im]) im=i;
  const K=64, fine=(u,w)=>{ const xs=new Float64Array(K+1); for(let i=0;i<=K;i++) xs[i]=u+(w-u)*i/K; return {xs,vs:ev(xs)}; };
  const fm=fine(x[Math.max(0,im-1)],x[Math.min(N-1,im+1)]); let jm=0; for(let i=1;i<=K;i++) if(fm.vs[i]>fm.vs[jm]) jm=i;
  let xm=fm.xs[jm], ym=fm.vs[jm]; if(jm>0&&jm<K){ const y0=fm.vs[jm-1],y1=fm.vs[jm],y2=fm.vs[jm+1], dd=y0-2*y1+y2; if(dd<0){ const s=0.5*(y0-y2)/dd; xm+=s*(fm.xs[1]-fm.xs[0]); ym=y1-0.25*(y0-y2)*s; } }
  const half=ym/2; let A=null,Bx=null;
  const cross=(i0,i1)=>{ const f=fine(x[i0],x[i1]); for(let i=0;i<K;i++){ const a=f.vs[i]-half, b=f.vs[i+1]-half; if(a===0) return f.xs[i]; if(a*b<0) return f.xs[i]+a/(a-b)*(f.xs[i+1]-f.xs[i]); } return null; };   // 간격 dx/64 → 선형 보간 오차 ~1e-9·폭
  for(let i=im;i>0;i--) if(v[i-1]<half){ A=cross(i-1,i); break; }
  for(let i=im;i<N-1;i++) if(v[i+1]<half){ Bx=cross(i,i+1); break; }
  const res={fwhm:A!=null&&Bx!=null?Bx-A:null,xmax:xm,ymax:ym,lo:A,hi:Bx}; if(S2FW.size>512) S2FW.clear(); S2FW.set(key,res); return res; }
/* Avantage Voigt 표기(SPEC §13b): M = fL/(fL+fG) 로 가정(설명서의 비 식은 그림이라 읽지 못함 — 대조 대기), fV = Olivero–Longbothum 0.5346 fL + √(0.2166 fL² + fG²) */
function avantageVoigt(wG,wL){ const fV=0.5346*wL+Math.sqrt(0.2166*wL*wL+wG*wG); return {fV,M:wL+wG>0?wL/(wL+wG):0,fG:wG}; }
function voigtFromAvantage(fV,M){ if(M>=1) return {wG:0,wL:fV}; const r=M/(1-M), fG=fV/(0.5346*r+Math.sqrt(0.2166*r*r+1)); return {wG:fG,wL:r*fG}; }
/* 전체 스펙트럼(be, y)에서 창 [i0, i1] 의 Casa 끝점 평균 — 평균 창(±aw)이 피팅 창 밖 채널을 포함할 수 있어서 */
function casaEnds(be,y,i0,i1,aw){ const N=be.length, w=Math.max(0,Math.round(aw||0)); const av=c=>{ let s=0,k=0; for(let j=Math.max(0,c-w);j<=Math.min(N-1,c+w);j++){ s+=y[j]; k++; } return s/k; };
  return {ylo:av(Math.min(i0+1,N-1)),yhi:av(Math.max(i1-1,0))}; }
/* CasaXPS 영역 표의 Area · Max Height · Min Height 와 같은 정의: 채널 i0+1 … i1−1 에서 (y − B) 의 사다리꼴 적분·최대·최소 */
function casaRegionStats(E,y,B){ const N=E.length, a=Math.min(1,N-1), b=Math.max(0,N-2); let t=0,mx=-Infinity,mn=Infinity;
  for(let i=a;i<=b;i++){ const r=y[i]-B[i]; if(i<b) t+=0.5*(r+y[i+1]-B[i+1])*(E[i+1]-E[i]); if(r>mx) mx=r; if(r<mn) mn=r; } return {area:t,maxH:mx,minH:mn}; }
/* Avantage 규약 Smart v1.8 (SPEC_real_validation §A — Kratos 계 실측 Avantage Peak Table + 'Backgnd.' 열 26 영역에서 역으로 정함; Avantage 코드는 쓰지 않음).
   v1.9.2: avantageV19(bg) 가 거짓이면(기본) 이것을 쓴다.
   ① 평활 S = 9 점 중심 이동평균(창 밖 4 채널까지 — bg._pad = {lo:[창 아래 채널], hi:[창 위 채널]} 의 가장 가까운 4 점, 없으면 창 안에서 자름).
   ② 끝점 = S(창 첫 채널), S(창 끝 채널).
   ③ 구간 [i, j] (i 낮은 BE): S_j > S_i 이면 Shirley(수렴, 피적분 S − B, 낮은 BE 쪽부터 포함 누적합), 아니면 직선.
   ④ 구간 안(끝 제외)에서 B − S 가 가장 큰 채널 k 가 0 보다 크면 그 채널(값 S_k)을 새 끝점으로 [i,k]·[k,j] 를 다시 계산(재귀). */
function bgSmartAvantageV18(E,y,bg){
  const N=E.length; if(N<3) return bgLinear(E,y,1);
  const pl=bg&&bg._pad&&bg._pad.lo?Array.from(bg._pad.lo).slice(-4):[], ph=bg&&bg._pad&&bg._pad.hi?Array.from(bg._pad.hi).slice(0,4):[];
  const Y=pl.concat(Array.from(y),ph), o=pl.length, S=new Float64Array(N);
  for(let i=0;i<N;i++){ const c=i+o; let s=0,k=0; for(let j=Math.max(0,c-4);j<=Math.min(Y.length-1,c+4);j++){ s+=Y[j]; k++; } S[i]=s/k; }
  const B=new Float64Array(N);
  const seg=(i,j)=>{ const a=S[i], b=S[j], n=j-i;
    if(!(b>a)){ for(let t=i;t<=j;t++) B[t]=a+(b-a)*(t-i)/(n||1); return; }
    for(let t=i;t<=j;t++) B[t]=a;
    for(let it=0;it<300;it++){ let q=0, Qt=0; for(let t=i+1;t<=j;t++) Qt+=S[t]-B[t]; if(!(Math.abs(Qt)>1e-300)){ for(let t=i;t<=j;t++) B[t]=a+(b-a)*(t-i)/(n||1); return; }
      let mx=0; const nb=new Float64Array(n+1); nb[0]=a; for(let t=i+1;t<=j;t++){ q+=S[t]-B[t]; nb[t-i]=a+(b-a)*q/Qt; }
      for(let t=i;t<=j;t++){ mx=Math.max(mx,Math.abs(nb[t-i]-B[t])); B[t]=nb[t-i]; } if(mx<1e-9*(Math.abs(b-a)+1)) break; } };
  const rec=(i,j,dp)=>{ seg(i,j); if(j-i<2||dp>60) return; let k=-1, ex=0; for(let t=i+1;t<j;t++){ const d=B[t]-S[t]; if(d>ex){ ex=d; k=t; } } if(k<0) return; rec(i,k,dp+1); rec(k,j,dp+1); };
  rec(0,N-1,0); B.smoothed=S; B.rule="v18"; return B;
}
/* Avantage 규약 Smart v1.9.2 (test/ref/avantage_a2_patch.md; 위 bgShirleyAvantage 의 ①–③ 과 같은 h · 모양):
   ①  끝값 = 끝 채널 ±h 평균(창 밖 포함).   ②  S = ±h 이동평균(창 안, 낮은 BE 끝 채널 제외) — 제약과 피적분 모두.
   ③  구간 [i, j] (i 낮은 BE): 높은 BE 끝값 > 낮은 BE 끝값이면 bgShirleyAvantage ③ 의 모양(Δ = 구간 끝값 차), 아니면 직선.
   ④  구간 안 B − S 최대 채널 k (> 0) 에서 나눔(값 = S_k), 재귀. B.splits = 접점 채널.                                   */
function bgSmartAvantage(E,y,bg){
  bg=bg||{}; if(!avantageV19(bg)) return bgSmartAvantageV18(E,y,bg);
  const N=E.length; if(N<3) return bgLinear(E,y,1);
  const h=avantageHalf(bg,E), p=bg.avP!=null?+bg.avP:0.07, q=bg.avQ!=null?+bg.avQ:0.25, hs=bg.avHs!=null?Math.round(bg.avHs):h;
  return avantageOriented(E,y,bg,(Y,pl,ph)=>{
    const {a,b}=avantageEnds(Y,pl,ph,h), S=avantageSmooth(Y,hs), V=Float64Array.from(S); V[0]=a; V[N-1]=b;
    const B=new Float64Array(N), splits=[];
    const seg=(i,j)=>{ const va=V[i], vb=V[j], n=j-i;
      if(!(vb>va)||n<2){ for(let t=i;t<=j;t++) B[t]=va+(vb-va)*(t-i)/(n||1); return; }
      const s=new Array(n+1); for(let t=0;t<=n;t++) s[t]=S[i+t]-va; const u=avantagePQ(s,vb-va,p,q); for(let t=0;t<=n;t++) B[i+t]=va+u[t]; };
    const rec=(i,j,dp)=>{ seg(i,j); if(j-i<2||dp>200) return; let k=-1, ex=0; for(let t=i+1;t<j;t++){ const d=B[t]-S[t]; if(d>ex){ ex=d; k=t; } } if(k<0) return; splits.push(k); rec(i,k,dp+1); rec(k,j,dp+1); };
    rec(0,N-1,0);
    B.smoothed=S; B.splits=splits.sort((x,z)=>x-z); B.half=h; B.rule="v19"; B.ends=[a,b];
    return B; });
}
/* 전체 스펙트럼에서 창 [i0, i1] 밖 h 채널씩 (Avantage 끝값 평균이 창 밖 채널까지 씀; 기본 h = 4) */
function avantagePad(y,i0,i1,h){ h=h!=null?Math.max(0,Math.round(h)):4; const N=y.length, lo=[], hi=[]; for(let j=Math.max(0,i0-h);j<i0;j++) lo.push(y[j]); for(let j=i1+1;j<=Math.min(N-1,i1+h);j++) hi.push(y[j]); return {lo,hi}; }
function bgStatic(E,y,bg){
  const t=bg.type, n=bg.n||5;
  if(bg.conv==="avantage"&&t==="shirley") return bgShirleyAvantage(E,y,bg);
  if(bg.conv==="avantage"&&t==="smart") return bgSmartAvantage(E,y,bg);
  if(bg.conv==="casa"&&(t==="linear"||t==="shirley")){ const p=endpointsCasa(E,y,bg.aw!=null?bg.aw:1,bg._ends);
    if(t==="shirley") return bgShirleyCasa(E,y,p);
    const sl=(p.yhi-p.ylo)/((p.Ehi-p.Elo)||1); return Float64Array.from(E,e=>p.ylo+sl*(e-p.Elo)); }
  if(t==="linear") return bgLinear(E,y,n);
  if(t==="shirley") return bgShirley(E,y,n,false);
  if(t==="smart") return bgShirley(E,y,n,true);
  if(t==="tougaard") return bgTougaard(E,y,n,{C:bg.C,D:bg.D,u3:bg.u3,T0:bg.T0});
  if(t==="spline") return bgSpline(E,y,bg);
  if(t==="none") return new Float64Array(E.length);
  return null;   // active
}
const BG_TYPES={linear:"Linear",shirley:"Shirley",smart:"Smart",tougaard:"Tougaard (U2/U3)",active:"Shirley (in-fit)",spline:"Spline (knots)",none:"None"};

/* ============================================================
   4) 문제 구성 · 매개변수 벡터 · 제약/연결
   prob = {E, y, c, bg:{type,n,C,D,u3, c0,k,s, cons:{c0,k,s}}, comps:[...], cons:{compId:{param:{mode,lo,hi,ref,k,b}}}}
   ============================================================ */
function paramNames(c){ const a=["be","area"].concat(shapePList(c)); if(c.doublet) a.push("dE","r","q"); return a; }
function getP(c,k){ return (k==="dE"||k==="r"||k==="q")?c.doublet[k]:c[k]; }
function setP(c,k,v){ if(k==="dE"||k==="r"||k==="q") c.doublet[k]=v; else c[k]=v; }
function defBounds(c,k){
  const d=PDEF[k]||{}; let lo=d.lo!=null?d.lo:-Infinity, hi=d.hi!=null?d.hi:Infinity;
  if(k==="be"){ lo=c.be-3; hi=c.be+3; }
  if((c.shape==="LA"||c.shape==="TLA")&&(k==="alpha")){ lo=LA_ALPHA_LO; hi=20; }
  if(c.shape==="TP"&&k==="F"){ lo=0.1; hi=200; }   // 배경 성분은 창 전체만큼 넓을 수 있다
  return [lo,hi];
}
function clone(o){ return JSON.parse(JSON.stringify(o)); }
/* 자유 매개변수 목록과 연결 해석 */
function buildMap(prob){
  const comps=prob.comps, byId={}; comps.forEach((c,i)=>byId[c.id]=i);
  const free=[], links=[], fixed=[];
  comps.forEach((c,ci)=>{
    paramNames(c).forEach(k=>{
      if(derivedParam(c,k)){ fixed.push({ci,k}); return; }   // CasaXPS LA(α,β,n): wG = n·CASA_LA_RHO·F (SPEC_casa_v2 §2); 새 모양도 같은 규칙(가정, SPEC_shapes2)
      const cs=((prob.cons||{})[c.id]||{})[k]||defaultCons(c,k)||{mode:"free"};
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
/* 일정 가중(배경을 뺀 자료 등 세기가 0·음수 근처라 Poisson 분산 모형이 맞지 않을 때): v_i = 상수(2차 차분 잡음의 중앙값 추정).
   절대 척도는 어차피 φ = χ²/(N−p) 로 다시 정해진다 (SPEC §3.1). */
function constVar(y){ const n=y.length, q=[]; for(let i=1;i<n-1;i++){ const d=y[i]-(y[i-1]+y[i+1])/2; q.push(d*d/1.5); } q.sort((a,b)=>a-b);
  let v=q.length?q[Math.floor(q.length/2)]/0.4549:1; if(!(v>0)){ let m=0; for(const t of y) m=Math.max(m,Math.abs(t)); v=Math.max(1e-12,(1e-3*m)**2); } return Float64Array.from(y,()=>v); }
/* CasaXPS 피팅 가중 (SPEC_casa_v2 §1, CasaXPS 2.3.27 V2 실측): 가중 없는 최소제곱, 창 [i0, i1] 의 양 끝 채널은 빼고 i0+1 … i1−1 만 (DoF = N − 2 − p) */
function casaVar(y){ const n=y.length; return Float64Array.from(y,(v,i)=>i===0||i===n-1?1e30:1); }
function dataVar(prob,ys){ return prob.weight==="const"?constVar(prob.y):prob.weight==="casa"?casaVar(prob.y):varianceFrom(ys,prob.c); }
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
/* 멈춤(stall) 규칙 (SPEC §3.3b, v1.9.2) — 옛 LM 보다 최종 χ² 가 나빠지지 않게 하면서 퇴화 피팅의 기어감을 줄인다.
   갇힘: 받아들인 걸음의 λ 가 W 번 연속 1 … 1e6 이고 그 W 번 중 rej 번 이상이 먼저 거절된 걸음을 겪음.
   현직(incumbent) = 옛 LM 그대로(비트 단위). 갇힌 뒤 grace 번이 지나도 끝나지 않고 W 걸음에 ≥ desc·φ̂ 씩 내려가는 중이면
   시도 갈래(branch) 하나를 현직과 번갈아 같은 계산량씩(현직 계산량 × share, 0.25 에서 시작해 shareRamp 걸음마다 1 씩 늘어 최대 share; 합계는
   갈래 전까지 현직이 쓴 것의 budget 배까지) 푼다: rescue = 한계에 붙었거나(bound) 모형이 보지 못하는(flat) 매개변수를 1차 시작값으로 되돌려
   고정 → 나머지 → 모두; 되돌릴 것이 없으면 kick = 같은 점에서 λ = 1e-3. 갈래 안은 v1.9.1 규칙(kick → stalled/converged, 수렴 기준 calmFast).
   갈래는 **잡음 바닥**(χ²/ν ≤ 1 + zFloor·√(2/ν))에 닿았을 때만 결과가 될 수 있다:
     · 현직이 끝까지 가서 갈래보다 높으면 갈래를 채택하고, 2차(IRLS)를 두 끝점(갈래·현직 = 옛 1차 끝점)에서 풀어 최종 χ² 가 낮은 쪽을 고른다;
     · 현직이 아직 가는 중이면 minInc 번 더 준 뒤 최근 두 창(W 걸음)이 모두 ≥ desc·φ̂ 씩 내려가고, 남은 감소의 외삽(감소 비 ≤ qGeo 면 등비, 아니면
       maxit 까지 직선)으로도 갈래보다 (1+cut) 배 넘게 높을 때만 현직을 버린다(평탄하거나 빨라지는 현직은 외삽할 수 없어 끝까지 둔다).
   2차의 현직 쪽(옛 2차)도 같은 규칙으로 갈래 쪽 최종 χ² 와 겨룬다. 갈래가 끝나기 전에 현직이 끝나면 갈래는 버린다(바닥 가까이면 끝까지 풀어 봄).
   adoptAny: true 면 바닥이 아니어도 더 낮은 갈래를 채택(χ² 가 크게 좋아지는 피팅이 많지만 시간이 듦). 갇히지 않거나 grace 안에 끝나는 피팅은 옛 LM 과 같다. */
const LM_STALL={W:8,lam:1,lamHi:1e6,rej:5,tau:0.3,conv:1e-6,grace:8,desc:0.3,tries:1,share0:0.25,shareRamp:8,share:2,budget:1.5,calmFast:1e-7,zFloor:1,minInc:8,qGeo:0.5,cut:0.01,adoptAny:false};   // 창·λ 범위·창 안 거절 수·멈춤 τ·수렴 ε_c · 갈래 전 기다림·갈래 수·현직을 버리는 여유 (SPEC §3.3b)
const FIT_HOOK={fit:null,lm:null,br:null,stall:true};   // 시험용 기록 고리 (test/fitlog_hook.js; 화면에서는 비어 있음)
/* 대칭 행렬 고유값·고유벡터 (Jacobi). 반환 {w: 고유값, V: 열 = 고유벡터(행 우선 k×k)} */
function eigSym(A,k){
  const a=Float64Array.from(A), V=new Float64Array(k*k); for(let i=0;i<k;i++) V[i*k+i]=1;
  for(let sweep=0;sweep<100;sweep++){ let off=0, dg=0; for(let p=0;p<k;p++){ dg+=a[p*k+p]*a[p*k+p]; for(let q=p+1;q<k;q++) off+=a[p*k+q]*a[p*k+q]; } if(off<=1e-30*dg) break;
    for(let p=0;p<k;p++) for(let q=p+1;q<k;q++){ const apq=a[p*k+q]; if(Math.abs(apq)<1e-300) continue;
      const th=(a[q*k+q]-a[p*k+p])/(2*apq), t=Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1)), cs=1/Math.sqrt(t*t+1), sn=t*cs;
      for(let r=0;r<k;r++){ const arp=a[r*k+p], arq=a[r*k+q]; a[r*k+p]=cs*arp-sn*arq; a[r*k+q]=sn*arp+cs*arq; }
      for(let r=0;r<k;r++){ const apr=a[p*k+r], aqr=a[q*k+r]; a[p*k+r]=cs*apr-sn*aqr; a[q*k+r]=sn*apr+cs*aqr; }
      for(let r=0;r<k;r++){ const vrp=V[r*k+p], vrq=V[r*k+q]; V[r*k+p]=cs*vrp-sn*vrq; V[r*k+q]=sn*vrp+cs*vrq; } } }
  return {w:Array.from({length:k},(_,i)=>a[i*k+i]),V};
}
/* 멈춘 피팅에서 잘 정해지지 않은 매개변수 고르기 (SPEC §3.3b): 마지막 Jacobian 의 JᵀJ (자연 매개변수 공간)
   - bound: 한계까지 남은 거리 < 1e-3·(hi − lo) (한쪽 한계면 < 1e-3·(|p| + 1))
   - flat : 전 범위(유한 한계) 또는 ±(|p|+0.1) 를 움직여도 Δχ² < 0.01·φ̂ (모형이 그 매개변수를 거의 보지 못함)
   - degenerate: 상관 행렬 R = D⁻¹ JᵀJ D⁻¹ 의 고유값 < 1e-9 인 방향에 성분 |v_j|² > 0.2 (다른 매개변수와 구별 안 됨) */
function illDetermined(map,p,A,u,phi){
  const np=p.length, out=[]; const dn=new Float64Array(np);
  for(let j=0;j<np;j++){ const f=map.free[j], d=dPdU(u[j],f.lo,f.hi); dn[j]=Math.abs(d)>1e-300?1/d:0; }
  const An=new Float64Array(np*np); for(let a=0;a<np;a++) for(let b=0;b<np;b++) An[a*np+b]=A[a*np+b]*dn[a]*dn[b];
  const why=p.map(()=>[]);
  map.free.forEach((f,j)=>{ const L=isFinite(f.lo), H=isFinite(f.hi);
    const dist=Math.min(L?p[j]-f.lo:Infinity,H?f.hi-p[j]:Infinity), ref=L&&H?f.hi-f.lo:Math.abs(p[j])+1;
    if(dist<1e-3*ref) why[j].push("bound");
    const span=L&&H?f.hi-f.lo:Math.abs(p[j])+0.1; if(!(An[j*np+j]*span*span>=0.01*phi)) why[j].push("flat"); });
  const act=[]; for(let j=0;j<np;j++) if(!why[j].length&&An[j*np+j]>0) act.push(j);
  const k=act.length;
  if(k>1){ const d=act.map(j=>Math.sqrt(An[j*np+j])); const R=new Float64Array(k*k); for(let a=0;a<k;a++) for(let b=0;b<k;b++) R[a*k+b]=An[act[a]*np+act[b]]/(d[a]*d[b]);
    const e=eigSym(R,k); for(let m=0;m<k;m++) if(e.w[m]<1e-9) for(let a=0;a<k;a++){ const v=e.V[a*k+m]; if(v*v>0.2&&why[act[a]].indexOf("degenerate")<0) why[act[a]].push("degenerate"); } }
  map.free.forEach((f,j)=>{ if(why[j].length) out.push({id:f.id,ci:f.ci,k:f.k,j,p:p[j],lo:f.lo,hi:f.hi,why:why[j]}); });
  return out;
}
/* LM 한 걸음씩 (SPEC §3.3). makeLM(prob,map,v,opt) → L: L.step() 이 한 반복을 하고 끝나면 L.done, 결과는 L.out().
   step 전에 prob 의 매개변수가 L 의 것이어야 한다(L.step 은 끝날 때 prob 를 L 의 상태로 둔다 — 여러 LM 을 번갈아 돌리는 쪽이 되돌림).
   정지 이유(stop): "converged" (상대 χ² 변화 < 1e-10 이 3 번 연속 · 어떤 감쇠로도 χ² 가 줄지 않음 · kick 뒤 감소 < ε_c·φ̂), "stalled" (§3.3b), "maxIter".
   opt._fast: v1.9.1 의 kick/멈춤 규칙으로 빨리 끝냄(시도 갈래 안, §3.3b). 아니면 옛 LM 그대로이고 갇힘(trapped)만 알려 준다 */
function makeLM(prob,map,v,opt){
  opt=opt||{}; const maxit=opt.maxit||500; const E=prob.E, y=prob.y, n=E.length;
  const np=map.free.length; const ev=makeEvaluator(prob,map);
  const sw=Float64Array.from(v,x=>1/Math.sqrt(x));
  const lo=map.free.map(f=>f.lo), hi=map.free.map(f=>f.hi);
  let p=getVec(prob,map).map((x,i)=>Math.min(hi[i],Math.max(lo[i],x)));
  let u=p.map((x,i)=>toU(x,lo[i],hi[i]));
  const chi=mu=>{ let s=0; for(let i=0;i<n;i++){ const r=(y[i]-mu[i])*sw[i]; s+=r*r; } return s; };
  setVec(prob,map,p); let cur=ev.full(); let X=chi(cur.mu);
  let nEff=0; for(let i=0;i<n;i++) if(sw[i]>1e-10) nEff++;
  const L={done:false,it:0,X,p,u,np,nEff,trapped:false,XH:[X],lastA:null};
  if(!np){ L.done=true; L.out=()=>({chi2:X,iter:0,conv:true,stop:"converged",mu:cur.mu,B:cur.B}); return L; }
  let lam=1e-3, nu=2, calm=0, it=0, conv=false, stop=null;
  const stall=opt.stall!==false&&FIT_HOOK.stall!==false, ST=opt.stallCfg?Object.assign({},LM_STALL,opt.stallCfg):LM_STALL; const hX=FIT_HOOK.lm?[X]:null, hL=hX?[]:null, hT=hX?[]:null;
  const fast=!!opt._fast, calmTol=fast?ST.calmFast:1e-10;   // 갈래 안은 조금 느슨하게 끝냄(채택되면 2차가 다시 다듬음)
  let damp=0, kicks=0, Xkick=X, kickNext=false, lastA=null, gainKick=null, probeLam=null; const rejW=[];
  const scale=map.free.map(f=>f.k==="be"||f.k==="dE"?1e-4:f.k==="area"?0:1e-5);
  function finish(){ L.done=true; if(!stop) stop="maxIter"; setVec(prob,map,p); cur=ev.full(); if(hX) FIT_HOOK.lm({X:hX,lam:hL,tries:hT,n:nEff,np,stop,kicks}); }
  function body(){   // 옛 for 루프의 몸통 한 번. true = 끝(break)
    if(kickNext){ lam=1e-3; nu=2; kickNext=false; }
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
    lastA=A;
    const D=new Float64Array(np); for(let a=0;a<np;a++) D[a]=Math.max(A[a*np+a],1e-30);
    let accepted=false, tries=0, lamUsed=lam;
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
        u=un; p=pn; cur=tr; X=Xn; accepted=true; lamUsed=lam;
        lam*=Math.max(1/3,1-Math.pow(2*rho-1,3)); nu=2;
        if(rel<calmTol) calm++; else calm=0; if(opt.onIter) opt.onIter(it,X);   // 진행 표시 (Worker → 화면)
        if(hX){ hX.push(X); hL.push(lamUsed); hT.push(tries); }
      } else { setVec(prob,map,p); lam*=nu; nu*=2; if(lam>1e16) break; }
    }
    if(!accepted){ setVec(prob,map,p); cur=ev.full(); conv=true; stop="converged"; return true; }
    L.XH.push(X);
    if(calm>=3){ conv=true; stop="converged"; return true; }
    // 멈춤 규칙 (SPEC §3.3b): 갇힘 판정
    if(stall){
      if(probeLam==null) probeLam=lamUsed;   // kick(또는 시작) 뒤 첫 걸음의 λ: < 1 이면 Gauss–Newton 에 가까운 걸음이 받아들여짐
      damp=lamUsed>=ST.lam&&lamUsed<=ST.lamHi?damp+1:0;   // λ > 1e6: 걸음이 GN 의 1e-6 미만 = 수치 잡음 한계(수렴 막바지) — 1e-10 기준에 맡김
      rejW.push(tries>1?1:0); if(rejW.length>ST.W) rejW.shift();
      L.trapped=damp>=ST.W&&rejW.reduce((a,b)=>a+b,0)>=ST.rej;   // 감쇠 지배 + 걸음이 자꾸 넘침(거절) = 갇힘. 거절 없이 λ 가 줄어 가는 고원은 기다린다
      if(fast&&L.trapped){   // 갈래 안 (v1.9.1): kick 한 번 → 다시 갇히면 kick 뒤 감소 Δ 로 수렴/멈춤
        const phi=X/Math.max(1,nEff-np), gain=Xkick-X;
        const probed=kicks>0||opt.afterStall;
        if(probed&&gain<ST.conv*phi&&probeLam<ST.lam){ conv=true; stop="converged"; gainKick=gain/phi; return true; }   // GN 걸음이 통했는데도 더 줄 것이 없음 = 최솟값
        if(probed&&gain<ST.tau*phi){ stop="stalled"; gainKick=gain/phi; return true; }
        kicks++; kickNext=true; Xkick=X; damp=0; rejW.length=0; probeLam=null; }
    }
    return false;
  }
  L.step=function(){ if(L.done) return; if(it>=maxit){ finish(); return; }
    const end=body(); L.X=X; L.p=p; L.u=u; L.lastA=lastA;
    if(end){ finish(); return; } it++; L.it=it; if(it>=maxit) finish(); };
  L.out=function(){ const o={chi2:X,iter:it,conv,stop,kicks,mu:cur.mu,B:cur.B,comps:ev.comps()};
    if((stop==="stalled"||stop==="maxIter")&&lastA){ const phi=X/Math.max(1,nEff-np); o.stall={gain:gainKick,phi,flags:illDetermined(map,p,lastA,u,phi)}; }
    return o; };
  L.abandon=function(){ L.done=true; stop="abandoned"; setVec(prob,map,p); cur=ev.full(); };
  return L;
}
/* LM 한 번 (SPEC §3.3, §3.3b). prob 은 제자리에서 바뀐다. v: 분산 배열.
   현직(incumbent) = 옛 LM 그대로. 갇힌 뒤 조건이 맞으면 시도 갈래 하나를 현직과 한 걸음씩 번갈아 푼다(같은 분산 → χ² 로 바로 비교).
   opt._pStart: rescue 가 되돌릴 1차 시작값. opt._beat {X, nit}: 2차 두 끝점 중 현직 쪽 — 갈래 쪽 최종 χ² 와만 겨룸 */
function lmFit(prob,map,v,opt){
  opt=opt||{}; const wStart=WORK.n, L=makeLM(prob,map,v,opt);
  const stall=opt.stall!==false&&FIT_HOOK.stall!==false, ST=opt.stallCfg?Object.assign({},LM_STALL,opt.stallCfg):LM_STALL;
  if(L.done||opt._fast||!stall){ while(!L.done) L.step(); return L.out(); }
  const maxit=opt.maxit||500, np=L.np, nEff=L.nEff, lo=map.free.map(f=>f.lo), hi=map.free.map(f=>f.hi);
  const pS=opt._pStart&&opt._pStart.length===np?opt._pStart:L.p.slice();
  let tFire=-1, br=null, nBr=0, cand=null, adopted=null, cutBreak=false, itx=0;
  if(opt._beat){ cand={X:opt._beat.X,at:0,nit:opt._beat.nit,virt:true}; nBr=ST.tries; }
  /* 시도 갈래 시작 (현직의 지금 점에서): rescue = 한계에 붙었거나(bound) 모형이 보지 못하는(flat) 매개변수를 1차 시작값으로 되돌려 고정 →
     나머지 풂 → 모두 풀어 다듬기; 되돌릴 것이 없으면 kick = 같은 점에서 λ = 1e-3. 갈래 안은 v1.9.1 규칙(_fast) */
  function startBranch(k){
    const p=L.p, phi=L.X/Math.max(1,nEff-np), fl=illDetermined(map,p,L.lastA,L.u,phi), q=p.slice(), reset=[];
    fl.forEach(f=>{ if(f.why.indexOf("bound")<0&&f.why.indexOf("flat")<0) return; const j=f.j, ref=isFinite(lo[j])&&isFinite(hi[j])?hi[j]-lo[j]:Math.abs(p[j])+1;
      if(Math.abs(pS[j]-p[j])>1e-3*ref){ reset.push(j); q[j]=pS[j]; } });   // 시작값이 지금과 다를 때만 되돌릴 수 있음
    const o2=Object.assign({},opt,{_fast:true,_pStart:null,_beat:null,afterStall:false,onIter:null}), nuB=Math.max(1,nEff-np);
    const B={kind:reset.length?"rescue":"kick",at:k,Xinc:L.X,nit:0,phase:reset.length?1:2,done:false,wInc:0,wBr:0,wMax:ST.budget*(WORK.n-wStart),
      reset:reset.map(j=>({id:map.free[j].id,ci:map.free[j].ci,k:map.free[j].k,from:p[j],to:pS[j]})),
      flags:fl.map(f=>({id:f.id,ci:f.ci,k:f.k,p:f.p,lo:f.lo,hi:f.hi,why:f.why}))};
    setVec(prob,map,q); B.lm=makeLM(prob,reset.length?{free:map.free.filter((f,j)=>reset.indexOf(j)<0),links:map.links,fixed:map.fixed}:map,v,o2);
    B.q=getVec(prob,map); setVec(prob,map,L.p);
    B.step=function(){ setVec(prob,map,B.q); B.lm.step(); B.nit++;
      if(B.lm.done){ if(B.phase===1){ B.phase=2; B.lm=makeLM(prob,map,v,o2); } else { B.done=true; B.r=B.lm.out(); B.X=B.r.chi2; } }
      else if(B.phase===2){ const H=B.lm.XH, Lh=H.length-1, W=ST.W;   // 다듬기 단계가 잦아드는데(등비 외삽) 잡음 바닥 근처(3z)에 못 미칠 것 → 쓸모없음(채택 조건이 바닥) → 버림
        if(Lh>=2*W){ const ga=H[Lh-2*W]-H[Lh-W], gb=H[Lh-W]-H[Lh], q=ga>0?gb/ga:Infinity; if(q<=ST.qGeo&&H[Lh]-gb*q/(1-q)>nuB*(1+3*ST.zFloor*Math.sqrt(2/nuB))) B.dead=true; } }
      B.q=getVec(prob,map); setVec(prob,map,L.p); };
    return B;
  }
  while(!L.done){
    const w0=WORK.n; L.step(); if(L.done) break;
    const k=L.it-1, XH=L.XH, X=L.X;   // 방금 끝난 반복의 번호
    if(L.trapped&&tFire<0) tFire=k;
    if(!br&&!cand&&nBr<ST.tries&&tFire>=0&&k-tFire>=ST.grace&&XH.length>ST.W&&XH[XH.length-1-ST.W]-X>=ST.desc*X/Math.max(1,nEff-np)){ br=startBranch(k); nBr++;  }   // 갇힌 뒤 grace 번 지나도 W 걸음에 ≥ desc·φ̂ 씩 기어 내려가는 중
    if(br){ br.wInc+=(WORK.n-w0)*Math.min(ST.share,ST.share0+(k-br.at)/ST.shareRamp); while(!br.done&&!br.dead&&br.wBr<=br.wInc&&br.wBr<=br.wMax){ const w1=WORK.n; br.step(); itx++; br.wBr+=WORK.n-w1; }
      if(br.done){ cand=br; br=null; } else if(br.wBr>br.wMax||br.dead) br=null; }   // 현직 한 걸음마다 갈래는 share 배 계산량까지; 갈래 전체는 현직이 갈래 전까지 쓴 것의 budget 배까지(넘으면 버림)
    if(cand){
      if(!(cand.X<X)) cand=null;   // 현직이 앞섬 → 갈래 버림(현직은 옛 경로 그대로)
      else if(k-cand.at>=ST.minInc){ const W=ST.W, Lx=XH.length-1;
        if(Lx>=2*W){ const ga=XH[Lx-2*W]-XH[Lx-W], gb=XH[Lx-W]-XH[Lx], lin=gb/W*(maxit-1-k), q=ga>0?gb/ga:Infinity, gEst=q<=ST.qGeo?Math.min(gb*q/(1-q),lin):lin;
          const nu=Math.max(1,nEff-np), thr=ST.desc*X/nu; if(cand.X<=nu*(1+ST.zFloor*Math.sqrt(2/nu))&&ga>=thr&&gb>=thr&&X-gEst>cand.X*(1+ST.cut)){ adopted=cand; cutBreak=true; break; } } }
      // 현직을 버리는 것은 꾸준히 내려가는 중(최근 두 창 W 걸음마다 감소 ≥ desc·φ̂)이고 남은 감소를 외삽해도(등비, 빨라지면 직선, maxit 까지)
      // 갈래보다 (1+cut) 배 넘게 높을 때만. 평탄한 현직은 곧 수렴하거나 갑자기 빠져나올 수 있어(외삽 불가) 끝까지 둔다
    }
  }
  const nuF=Math.max(1,nEff-np), atFloor=(x,z)=>x<=nuF*(1+z*ST.zFloor*Math.sqrt(2/nuF));   // 잡음 바닥: χ²ν ≤ 1 + z·√(2/ν)
  if(br&&!adopted&&br.lm.X<L.X&&(ST.adoptAny||atFloor(br.lm.X,3))){ while(!br.done){ br.step(); itx++; } cand=br; br=null; }   // 현직이 먼저 끝났는데 갈래가 이미 더 낮고 바닥 가까이(3z 안) → 갈래를 끝까지 (아니면 버림: 현직 = 옛 결과)
  if(!adopted&&cand&&!cand.virt&&cand.X<L.X&&(ST.adoptAny||atFloor(cand.X,1))) adopted=cand;   // 현직이 끝까지 갔는데 갈래가 더 낮고 잡음 바닥에 있음(같은 분산)
  if(FIT_HOOK.br&&(cand||adopted)&&!(adopted||cand).virt){ const c=adopted||cand; FIT_HOOK.br({tFire,at:c.at,kind:c.kind,nit:c.nit,candX:c.X,Xend:L.X,adopt:!!adopted,cut:cutBreak,XH:L.XH.slice(Math.max(0,tFire-8))}); }
  if(adopted&&adopted.virt){ L.abandon(); const o=L.out(); o.tries=0; return o; }   // 상대(갈래 쪽 최종)보다 낮아질 수 없음 → 그만
  if(!adopted){ if(!L.done){ L.abandon(); } const o=L.out(); o.iter+=itx; o.tries=nBr; return o; }
  // 갈래 채택: 현직 끝점(끝까지 갔으면 옛 LM 의 끝점)은 fitOnce 가 2차를 거기서도 풀도록 넘김
  let inc=null; if(cutBreak) L.abandon(); else { const lo_=L.out(); inc={p:L.p.slice(),X:L.X,mu:lo_.mu,stop:lo_.stop}; }   // 현직 끝점(= 옛 LM 의 1차 끝점)
  setVec(prob,map,adopted.q); const a=adopted.r, ev=makeEvaluator(prob,map); const cu=ev.full();
  const out={chi2:a.chi2,iter:L.it+itx,conv:a.conv,stop:a.stop,kicks:a.kicks||0,mu:cu.mu,B:cu.B,comps:ev.comps(),tries:nBr,inc,
    branch:{kind:adopted.kind,reset:adopted.reset,flags:adopted.flags,Xinc:adopted.Xinc,X:adopted.X,at:adopted.at,cut:cutBreak}};
  if(a.stall) out.stall=a.stall; return out;
}
/* 한 번의 완전한 피팅: (넓이 초기화) → LM(평활 분산) → LM(모형 분산) */
function fitOnce(prob,opt){
  opt=opt||{};
  let map=buildMap(prob); applyLinks(prob,map);
  const E=prob.E, y=prob.y; let ymax=0; for(const t of y) ymax=Math.max(ymax,Math.abs(t));
  const ys=smoothQ(E,y,smoothHW(E));
  let v=dataVar(prob,ys);
  if(opt.initAreas!==false){ initAreas(prob,map,v); if(prob.bg.type==="active"){ map=buildMap(prob); applyLinks(prob,map); } }   // k 상한 반영
  const fitter=(opt.optimizer||prob.optimizer)==="simplex"?nmFit:lmFit;   // SPEC §3.3 (Simplex: §14)
  const o1=Object.assign({},opt,{_pStart:getVec(prob,map)});   // rescue 가 되돌릴 시작값 (SPEC §3.3b)
  let r=fitter(prob,map,v,o1); const it1=r.iter;
  const r1=r; let chi2p1=r1.chi2, pick=null;
  if(opt.irls!==false&&prob.weight!=="const"&&prob.weight!=="casa"){
    const o2=r1.stop==="stalled"?Object.assign({},o1,{afterStall:true}):o1;   // 1차가 멈췄으면 2차의 시작(λ = 1e-3)을 Gauss–Newton 시도로 셈 (SPEC §3.3b)
    if(r1.inc){   // 1차가 갈래를 골랐음 → 2차를 두 끝점(갈래·현직)에서 풀어 최종 χ² 가 낮은 쪽 (SPEC §3.3b)
      const vA=varianceFrom(r1.mu,prob.c), rA=fitter(prob,map,vA,o2), pA=getVec(prob,map);
      // 현직 끝점 = 옛 1차 끝점 → 옛 2차 그대로(갈래 쪽 최종이 잡음 바닥에 있고 이쪽이 그보다 낮아질 수 없다고 판단될 때만 그만)
      setVec(prob,map,r1.inc.p); const vB=varianceFrom(r1.inc.mu,prob.c), rB=fitter(prob,map,vB,Object.assign({},o1,{_beat:{X:rA.chi2,nit:rA.iter}}));
      pick={branch:rA.chi2,incumbent:rB.chi2};
      if(rB.chi2<rA.chi2){ r=rB; v=vB; chi2p1=r1.inc.X; rB.iter+=rA.iter; pick.use="incumbent"; }
      else { setVec(prob,map,pA); r=rA; v=vA; rA.iter+=rB.iter; pick.use="branch"; }
    } else { v=varianceFrom(r1.mu,prob.c); r=fitter(prob,map,v,o2); }
  }
  // 정지 이유: 두 번 중 하나라도 멈췄으면(stalled) 그것을, 아니면 마지막 것을 (SPEC §3.3b)
  const viaInc=!!(pick&&pick.use==="incumbent");
  if(r!==r1){ r.stop1=viaInc?r1.inc.stop:r1.stop; r.kicks=(r.kicks||0)+(viaInc?0:(r1.kicks||0));
    if(!viaInc&&r1.stop==="stalled"&&r.stop!=="stalled"){ r.stop="stalled"; if(!r.stall) r.stall=r1.stall; }
    if(!viaInc&&r1.branch&&!r.branch) r.branch=r1.branch; }
  r.trig=(r1.tries||0)+(r!==r1?(r.tries||0):0); if(pick) r.pick=pick;   // trig: 갈래를 푼 횟수(0 이면 옛 LM 과 같은 경로)
  r.v=v; r.map=map; r.chi2p1=chi2p1; r.iter1=it1;   // 1차 χ²(데이터로 정한 분산) — 시작점끼리 비교하는 공통 목적함수
  if(FIT_HOOK.fit) FIT_HOOK.fit(prob,r,opt);   // 시험용 기록 고리 (test/fitlog_hook.js; 화면에서는 비어 있음)
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
  const nfree=np; const dof=Math.max(1,n-k-(prob.weight==="casa"?2:0));
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
  // CasaXPS 'RMS'·DoF (weight = "casa"; SPEC_casa_v2 §1, SPEC_casa_v3 §2 실측): √Σ(y − μ)² over channels i0+1 … i1−1;
  // DoF = #{i0+1 … i1−1 : y_i ≠ B_i} − p (B = region background; channels where the data equal the background exactly are not counted:
  // GL(30,100) 154−3 = 151, GL(30) 159−3 = 156, GL(30)T(2) 340−3 = 337). Active (in-fit) background: all channels counted.
  let casa=null; if(prob.weight==="casa"){ const Bc=prob.bg&&prob.bg.type!=="active"?(prob._Bst||bgStatic(prob.E,y,prob.bg)):null; let s2=0, ne=0;
    for(let i=1;i<n-1;i++){ const d=y[i]-mu[i]; s2+=d*d; if(!Bc||Math.abs(y[i]-Bc[i])>1e-12*Math.max(1,Math.abs(y[i]))) ne++; }
    casa={rms:Math.sqrt(s2),dof:Math.max(1,ne-p),n:ne}; }
  // 모형 비교용 χ²_d: 데이터만으로 정한 분산(평활 데이터) — 모형마다 같은 가중 (SPEC 3.5 v1.1)
  let ymax=0; for(const t of y) ymax=Math.max(ymax,Math.abs(t));
  const vd=dataVar(prob,smoothQ(prob.E,y,smoothHW(prob.E))); let chi2d=0; for(let i=0;i<n;i++){ const d=y[i]-mu[i]; chi2d+=d*d/vd[i]; }
  const AICc=n*Math.log(chi2d/n)+2*p+(n-p-1>0?2*p*(p+1)/(n-p-1):Infinity), BIC=n*Math.log(chi2d/n)+p*Math.log(n);
  return {N:n,p,dof,chi2,chi2d,chi2nu:chi2/dof,R2:ss>0?1-sr/ss:null,rmsNorm:Math.sqrt(chi2/n),DW:a0>0?dw/a0:null,rho1,ar1Inflate:rho1>0?Math.sqrt((1+rho1)/(1-Math.min(rho1,0.95))):1,runs,runsZ:zR,AICc,BIC,resid:r,casa};
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
  // 어떤 자유 매개변수에도 의존하지 않는 값(고정·고정에 연결)은 σ = null (화면 '—')
  base.forEach((b,i)=>{ b.sig={}; keys.forEach(k=>{ const g=grads[i][k]; if(!g||!g.some(x=>x!==0)){ if(b[k]!=null) b.sig[k]=null; return; } let s=0; for(let a=0;a<np;a++){ if(!g[a]) continue; for(let c=0;c<np;c++){ if(!g[c]) continue; const v=cov.C[a*np+c]; if(isFinite(v)) s+=g[a]*v*g[c]; } } b.sig[k]=Math.sqrt(Math.max(0,s)); }); });
  return base;
}

/* Monte Carlo (SPEC §3.4): 최적 μ + AR(1) 잡음 → 재피팅 */
function monteCarlo(prob,res,M,opt){
  opt=opt||{}; const R=rng(opt.seed||12345); const st=fitStats(prob,res);
  const rho=Math.max(0,Math.min(0.95,st.rho1)); const phi=st.chi2nu; const n=prob.E.length;
  const base=clone(stripProb(prob)); base.E=prob.E;
  const samples=[]; let fails=0;
  // 기준 스펙트럼 잡음(SPEC §12.5): 측정 선모양 성분의 template 을 원자료에 잡음을 더해 다시 만든다(배경 빼기·넓이 정규화 포함)
  const refNoise=opt.refNoise!==false&&base.comps.some(c=>c.shape==="MEAS"&&c.tpl&&c.tpl.raw);
  for(let m=0;m<M;m++){
    const y=new Float64Array(n); let e=0;
    for(let i=0;i<n;i++){ const z=R.g(); e=i?rho*e+Math.sqrt(1-rho*rho)*z:z; y[i]=res.mu[i]+e*Math.sqrt(phi*res.v[i]); }
    const P=clone(base); P.E=prob.E; P.y=y;
    if(refNoise){ const pt={}; P.comps.forEach(c=>{ if(c.shape!=="MEAS"||!c.tpl||!c.tpl.raw) return; const id=c.tpl.id; if(!pt[id]) pt[id]=perturbTpl(c.tpl,R); c.tpl=pt[id]; }); }
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
  if(tf.model==="avantage"){ const pe=PE||tf.pe; if(!(pe>0)) return NaN; const t=Math.log10(KE/pe); let s=0; for(let k=tf.coef.length-1;k>=0;k--) s=s*t+tf.coef[k]; return pe*Math.pow(10,s); }   // manual, Transmission Function: log10[A/(PE·XSF)] = Σ aₖ (log10(KE/PE))ᵏ
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
  if(o.model==="power") return Math.pow(KE,o.n!=null?o.n:0.7);
  if(o.model==="tpp2m") return tpp2m(KE,o.mat);
  return 1;
}
/* rows=[{label, el, area, sArea(σ), S, BE, hv, PE, tf}], opts={imfp} */
/* 넓이 ≤ 0 이거나 2σ 미만(σ 를 알 때)이면 '검출 안 됨'(nd): at% 를 내지 않고 Σ 에서 뺀다 — 음의 넓이가 다른 원소의 at% 를 부풀리지 않게 (SPEC §4) */
function quantify(rows,opts){
  opts=opts||{};
  const out=rows.map(r=>{ const KE=r.hv-r.BE; const T=tfValue(r.tf,KE,r.PE); const L=imfpValue(opts.imfp,KE); const den=r.S*T*L;
    const nd=!(r.area>0)||(r.sArea!=null&&r.sArea>0&&r.area<2*r.sArea);
    return Object.assign({},r,{KE,T,L,den,nrm:r.area/den,sNrm:r.sArea!=null?r.sArea/den:null,nd}); });
  const use=r=>!r.nd&&isFinite(r.nrm);
  const tot=out.reduce((s,r)=>s+(use(r)?r.nrm:0),0);
  out.forEach((r,i)=>{ if(!use(r)||!(tot>0)){ r.at=null; r.sAt=null; return; } r.at=100*r.nrm/tot;
    // 델타: ∂at_i/∂n_j = 100(δij·tot − n_i)/tot²
    let v=0; out.forEach((q,j)=>{ if(!use(q)||q.sNrm==null) return; const g=100*((i===j?tot:0)-r.nrm)/(tot*tot); v+=g*g*q.sNrm*q.sNrm; }); r.sAt=Math.sqrt(v); });
  return out;
}
/* 템플릿·모형을 BE 축에서 d 만큼 옮김(성분 BE + 'free' BE 한계). link 의 b(간격)는 상대값이라 그대로. 원본은 건드리지 않고 복사본을 돌려준다.
   저장(+Δ: 측정 → 보정 BE) 뒤 불러오기(−Δ)는 어떤 Δ 에서도 항등이어야 한다. */
function shiftModelBE(comps,cons,d){
  const cs=clone(comps||[]), ks=clone(cons||{}); d=+d||0;
  cs.forEach(c=>{ c.be+=d; });
  Object.values(ks).forEach(o=>{ const b=o&&o.be; if(b&&b.mode!=="link"){ if(b.lo!=null&&isFinite(b.lo)) b.lo+=d; if(b.hi!=null&&isFinite(b.hi)) b.hi+=d; } });
  return {comps:cs,cons:ks};
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
  /* v1.9.5: Booklet 먼저 (원자 준위 편집값; Bearden & Burr 1967, Fuggle & Mårtensson 1980 등 원측정 — xraydb 와 칸별 대조됨).
     Avantage 원소표 위치는 화학 상태 표지라 정수로 반올림된 것이 많아(예: Ce 3d 18 vs 18.6, Xe 3d 15 vs 12.6, Ge 2p 26.5 vs 31.1)
     둘 다 있을 때 357 개 중 78 개가 ≥ 0.3 eV 달랐다 → Booklet 이 없거나 soWarn 일 때만 Avantage 를 씀. */
  /* 단, Δ < 2 eV 에서는 Booklet 의 0.1 eV 반올림이 상대 오차로 커서(Si 2p 0.4 vs 고분해능 문헌 ≈0.6, S 2p 1.1 vs ≈1.2) 둘 다 정밀하지 않다 —
     이 범위는 v1.9.4 까지와 같이 Avantage 먼저(시험으로 검증된 동작 유지), Δ 제약 범위(±max(0.3, 5 %, 두 출처 차+0.1))가 차이를 흡수. */
  if(bkOK&&(bk>=2||av==null)){ dE=bk; src="X-ray Data Booklet "+el+" "+orb+J[1]+" − "+orb+J[0]; }
  else if(av!=null){ dE=av; src="Avantage 원소표 "+el+" "+orb+J[1]+" − "+orb+J[0]+(bkOK?" (Δ < 2 eV: Booklet 반올림 정밀도 부족)":" (Booklet 값 없음)"); }
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
const SENS_BG_LABEL={smart:"Smart",shirley:"Shirley",active:"Shirley (in-fit, +slope)",tougaard:"Tougaard",linear:"Linear"};
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
  // Spline 배경(SPEC §11 v1.7): 사용자 매듭 하나씩 ±dk 옮김(자동 y 는 새 자리에서 다시, 끌어 둔 y 는 그대로) — 다른 요인은 기준 수준(두 방식 모두 하나씩)
  if(cur.type==="spline"&&(cur.knots||[]).length){ const dk=o.dk!=null?o.dk:dw;
    cur.knots.forEach((k,ki)=>[-dk,dk].forEach(dd=>{ const v=mk(bgLv[0],0,0); v.bg.knots[ki].e=k.e+dd; v.f.kn=ki+(dd<0?"-":"+"); out.push(v); })); }
  return out;
}
/* 변형들을 실제로 맞춘다(Worker 안에서). spec={be,y,c}, model={lo,hi,bg,comps,cons}.
   각 변형은 지금 성분 값에서 출발(넓이는 NNLS 로 다시). 성분 중심이 창 밖이면 건너뜀. */
function sensJob(spec,model,variants){
  return variants.map((v,j)=>{ const i=v.idx!=null?v.idx:j;   // 여러 Worker 로 나눌 때 원래 번호
    try{
      // 지금 창 안에 있던 성분 중심이 변형 창에서 잘려 나가면 건너뜀(창 밖에 있던 성분 — 꼬리·가장자리 — 은 사용자 모형 그대로 둔다)
      if(model.comps.some(c=>c.be>=model.lo&&c.be<=model.hi&&(c.be<v.lo||c.be>v.hi))) return {i,ok:false,err:"성분 중심이 창 밖으로 잘림"};
      const P=makeProblem(spec,{lo:v.lo,hi:v.hi,bg:v.bg,comps:model.comps,cons:model.cons,weight:model.weight,optimizer:model.optimizer});
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
      const same=(f,fk)=>(fk==="kn"||(f.kn||"0")===(b0.kn||"0"))&&["bg","n","lo","hi"].every(k=>k===fk||(fk==="bg"&&k==="n"&&(f.n==null||f.n===f.nRef))||f[k]===b0[k]);
      const eff={}; ["bg","n","lo","hi"].concat(variants.some(v=>v&&v.f&&v.f.kn)?["kn"]:[]).forEach(fk=>{ const g={}; ok.forEach(r=>{ const v=variants[r.i]; const c=r.comps.find(x=>x.id===id); if(!v||!c||c[q]==null||!same(v.f,fk)) return; const lv=String(v.f[fk]); (g[lv]=g[lv]||[]).push(c[q]); });
        const ms=Object.values(g).map(a=>a.reduce((s,x)=>s+x,0)/a.length); eff[fk]=ms.length>1?Math.max(...ms)-Math.min(...ms):0; });
      s.eff=eff; o[q]=s; });
    return o; });
  return {n:variants.length,nOk:ok.length,nFail:variants.length-ok.length,base:!!base,comps,chi2nu:st(ok.map(r=>r.chi2nu))};
}

/* ============================================================
   13) 측정 선모양 · 기준 스펙트럼 선형 결합 (SPEC §12)
   CasaXPS 'line shapes defined from spectra' · Avantage 'NLLS fitting to reference levels' 와 같은 목적 — 식은 SPEC 에서 정의
   ============================================================ */
/* PCHIP 기울기 — scipy.interpolate.PchipInterpolator 와 같은 정의(Fritsch–Butland 가중 조화평균, 끝은 3점 한쪽식 + 모양 보존) */
function pchipSlopes(x,y){
  const n=x.length, m=new Float64Array(n); if(n<2) return m;
  const h=new Float64Array(n-1), d=new Float64Array(n-1);
  for(let k=0;k<n-1;k++){ h[k]=x[k+1]-x[k]; d[k]=(y[k+1]-y[k])/h[k]; }
  if(n===2){ m[0]=m[1]=d[0]; return m; }
  for(let k=1;k<n-1;k++){ const a=d[k-1], b=d[k];
    if(Math.sign(a)!==Math.sign(b)||a===0||b===0){ m[k]=0; continue; }
    const w1=2*h[k]+h[k-1], w2=h[k]+2*h[k-1]; m[k]=(w1+w2)/(w1/a+w2/b); }
  const edge=(h0,h1,m0,m1)=>{ let t=((2*h0+h1)*m0-h0*m1)/(h0+h1); if(Math.sign(t)!==Math.sign(m0)) t=0; else if(Math.sign(m0)!==Math.sign(m1)&&Math.abs(t)>3*Math.abs(m0)) t=3*m0; return t; };
  m[0]=edge(h[0],h[1],d[0],d[1]); m[n-1]=edge(h[n-2],h[n-3],d[n-2],d[n-3]);
  return m;
}
/* 보간값: 마디 밖은 0 (template 은 배경을 뺀 봉우리). m=null 이면 선형 */
function interpAt(x,y,m,t){
  const n=x.length; if(!(t>=x[0]&&t<=x[n-1])) return 0;
  let lo=0, hi=n-1; while(hi-lo>1){ const mid=(lo+hi)>>1; if(x[mid]<=t) lo=mid; else hi=mid; }
  const h=x[hi]-x[lo], s=(t-x[lo])/h;
  if(!m) return y[lo]+(y[hi]-y[lo])*s;
  const s2=s*s, s3=s2*s;
  return (2*s3-3*s2+1)*y[lo]+(s3-2*s2+s)*h*m[lo]+(3*s2-2*s3)*y[hi]+(s3-s2)*h*m[hi];
}
/* 보간 함수의 정확한 적분: Hermite 3차 조각 ∫ = h(y_k+y_{k+1})/2 + h²(m_k − m_{k+1})/12, 선형 = 사다리꼴 */
function interpIntegral(x,y,kind){
  const n=x.length, m=kind==="linear"?null:pchipSlopes(x,y); let s=0;
  for(let k=0;k<n-1;k++){ const h=x[k+1]-x[k]; s+=0.5*h*(y[k]+y[k+1]); if(m) s+=h*h*(m[k]-m[k+1])/12; }
  return s;
}
/* template = {kind:"tpl", id, name, interp, eref(보정 BE), x(= 보정 BE − eref, 오름차순), v(넓이 1), area0, raw:{E,y,bg|B,shift,noise,kappa}}
   E 는 측정(보정 전) BE 오름차순, o.shift = 그 스펙트럼의 대전 보정 Δ. eref 를 주지 않으면 배경 뺀 신호의 (평활) 최댓값 위치. */
function makeTemplate(E,y,o){
  o=o||{}; const n=E.length; if(n<5) throw new Error("기준 스펙트럼 영역의 점이 5개 미만입니다");
  const bg=o.B?null:clone(o.bg||{type:"shirley",n:5});
  let B=o.B?Float64Array.from(o.B):bgStatic(E,y,bg); if(!B) B=bgShirley(E,y,(bg&&bg.n)||5,true);
  const z=Float64Array.from(y,(t,i)=>t-B[i]); const sh=+o.shift||0;
  let eref=o.eref;
  if(eref==null||!isFinite(eref)){ const zs=smoothQ(E,z,smoothHW(E)); let im=0; for(let i=1;i<n;i++) if(zs[i]>zs[im]) im=i;
    let e=E[im]; if(im>0&&im<n-1){ const a=zs[im-1],b=zs[im],c=zs[im+1], dd=a-2*b+c; if(dd<0){ const off=Math.max(-0.5,Math.min(0.5,0.5*(a-c)/dd)); e=E[im]+off*(off>0?E[im+1]-E[im]:E[im]-E[im-1]); } }
    eref=e+sh; }
  const interp=o.interp==="linear"?"linear":"pchip";
  const x=Array.from(E,e=>e+sh-eref);
  const I=interpIntegral(x,Array.from(z),interp);
  if(!(I>0)) throw new Error("배경을 뺀 기준 스펙트럼의 넓이가 0 이하입니다 — 영역·배경을 확인하세요");
  const v=Array.from(z,t=>t/I);
  let raw=null;
  if(o.keepRaw!==false){ const noise=o.noise==="const"?"const":"poisson"; const nz=noise==="poisson"?noiseEstimate(E,y):null;
    raw={E:Array.from(E),y:Array.from(y),bg,B:o.B?Array.from(o.B):null,shift:sh,noise,kappa:nz&&nz.kappa>0?nz.kappa:null}; }
  return {kind:"tpl",id:o.id||("t"+Math.random().toString(36).slice(2,8)),name:o.name||"",interp,eref,x,v,area0:I,raw};
}
/* 원자료에 잡음을 더해 template 을 다시 만든다(Monte Carlo). 분산: poisson 이면 κ·max(ŷ, 1e-6·max ŷ) (κ: SPEC §6 2차 차분 추정), const 면 constVar.
   eref 는 고정(성분 위치의 뜻이 바뀌지 않게) */
function perturbTpl(t,R){
  const r=t.raw; if(!r) return t; const E=Float64Array.from(r.E), y0=Float64Array.from(r.y), n=E.length;
  let va;
  if(r.noise==="const") va=constVar(y0);
  else { const ys=smoothQ(E,y0,smoothHW(E)); let mx=0; for(const q of ys) mx=Math.max(mx,q); const k=r.kappa>0?r.kappa:1; va=Float64Array.from(ys,q=>k*Math.max(q,1e-6*mx)); }
  const y=new Float64Array(n); for(let i=0;i<n;i++) y[i]=y0[i]+Math.sqrt(va[i])*R.g();
  return makeTemplate(E,y,{bg:r.bg,B:r.B,shift:r.shift,eref:t.eref,interp:t.interp,id:t.id,name:t.name,keepRaw:false});
}
/* 캐시: 기울기 · 넓힘 격자 (template 객체마다; JSON 복사본은 새 객체라 새로 만든다) */
const TPLC=typeof WeakMap!=="undefined"?new WeakMap():null;
function tplCache(t){
  if(!t||!t.x) throw new Error("측정 선모양의 template 이 없습니다");
  let c=TPLC&&TPLC.get(t); if(c) return c;
  const x=Float64Array.from(t.x), v=Float64Array.from(t.v); let st=Infinity; for(let i=1;i<x.length;i++) st=Math.min(st,x[i]-x[i-1]);
  c={x,v,m:t.interp==="linear"?null:pchipSlopes(x,v),step:st,conv:[]}; if(TPLC) TPLC.set(t,c); return c;
}
/* 넓이 1 측정 선모양 (상대좌표 X = E − be). wG > 0: 균일 격자(h = 0.04·2^−k ≤ min(wG/16, 점 간격/4))에서 보간값 ⊗ 이산 정규화 Gaussian(±5σ), 격자 → 점은 Catmull–Rom */
function tplShape(t,X,wG){
  const c=tplCache(t), n=X.length, out=new Float64Array(n);
  if(!(wG>1e-6)){ for(let i=0;i<n;i++) out[i]=interpAt(c.x,c.v,c.m,X[i]); return out; }
  let g=null; for(const q of c.conv) if(q.wG===wG){ g=q; break; }
  if(!g){ let h=0.04; const lim=Math.min(wG/16,c.step/4); while(h>lim&&h>0.0005) h/=2;
    const P=5*wG*FW2SIG+h, a=c.x[0]-P, b=c.x[c.x.length-1]+P; while((b-a)/h>400000) h*=2;
    const k0=Math.floor(a/h), k1=Math.ceil(b/h), M=k1-k0+1; const base=new Float64Array(M);
    for(let i=0;i<M;i++) base[i]=interpAt(c.x,c.v,c.m,(k0+i)*h);
    g={wG,h,k0,M,v:gconv(base,h,wG)}; c.conv.push(g); if(c.conv.length>6) c.conv.shift(); }
  const xa=g.k0*g.h, xb=(g.k0+g.M-1)*g.h;
  for(let i=0;i<n;i++){ const xi=X[i]; out[i]=xi<xa||xi>xb?0:crAt(g.v,g.k0,g.h,xi); }
  return out;
}
/* 대칭 행렬 고윳값 (Jacobi 회전, 작은 행렬용) */
function symEig(A,k){
  const a=Float64Array.from(A);
  for(let sweep=0;sweep<100;sweep++){ let off=0; for(let p=0;p<k;p++) for(let q=p+1;q<k;q++) off+=a[p*k+q]*a[p*k+q]; if(off<1e-30) break;
    for(let p=0;p<k;p++) for(let q=p+1;q<k;q++){ const apq=a[p*k+q]; if(Math.abs(apq)<1e-300) continue;
      const th=(a[q*k+q]-a[p*k+p])/(2*apq), t=Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1)), cs=1/Math.sqrt(t*t+1), sn=t*cs;
      for(let r=0;r<k;r++){ const arp=a[r*k+p], arq=a[r*k+q]; a[r*k+p]=cs*arp-sn*arq; a[r*k+q]=sn*arp+cs*arq; }
      for(let r=0;r<k;r++){ const apr=a[p*k+r], aqr=a[q*k+r]; a[p*k+r]=cs*apr-sn*aqr; a[q*k+r]=sn*apr+cs*aqr; } } }
  return Array.from({length:k},(_,i)=>a[i*k+i]).sort((x,y)=>x-y);
}
/* 공선성 진단 (SPEC §12.4): 측정 선모양 성분(넓이 1, 지금 위치·넓힘)을 데이터 점에서 계산해 1/√v 가중한 열 X_j.
   조건수 κ = √(λmax/λmin), λ = 열을 단위 길이로 맞춘 XᵀX 의 고윳값(Belsley 척도 조건수). VIF_j = (R⁻¹)_jj, R = 열의 상관행렬(평균을 뺀 열).
   경고: κ > 30 또는 VIF > 10 */
function lcfDiag(prob,v){
  const E=prob.E, n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const idx=prob.comps.map((c,i)=>c.shape==="MEAS"&&c.tpl?i:-1).filter(i=>i>=0); if(idx.length<2) return null;
  const X=idx.map(i=>{ const c=Object.assign({},prob.comps[i],{area:1,doublet:null}); const u=compEval(c,E,ctx); if(v) for(let k=0;k<n;k++) u[k]/=Math.sqrt(v[k]); return u; });
  const k=X.length, dot=(a,b)=>{ let s=0; for(let i=0;i<n;i++) s+=a[i]*b[i]; return s; };
  const nr=X.map(c=>Math.sqrt(dot(c,c))); const G=new Float64Array(k*k);
  for(let a=0;a<k;a++) for(let b=0;b<k;b++) G[a*k+b]=nr[a]>0&&nr[b]>0?dot(X[a],X[b])/(nr[a]*nr[b]):(a===b?1:0);
  const ev=symEig(G,k); const cond=ev[0]>0?Math.sqrt(ev[k-1]/ev[0]):Infinity;
  const Xc=X.map(c=>{ let m=0; for(let i=0;i<n;i++) m+=c[i]; m/=n; return Float64Array.from(c,t=>t-m); }); const nc=Xc.map(c=>Math.sqrt(dot(c,c)));
  const Rm=new Float64Array(k*k); for(let a=0;a<k;a++) for(let b=0;b<k;b++) Rm[a*k+b]=nc[a]>0&&nc[b]>0?dot(Xc[a],Xc[b])/(nc[a]*nc[b]):(a===b?1:0);
  const Ri=invSym(Rm,k); const vif=Ri?Array.from({length:k},(_,i)=>Ri[i*k+i]):Array(k).fill(Infinity);
  return {cond,vif,ids:idx.map(i=>prob.comps[i].id),warn:!(cond<=30)||vif.some(x=>!(x<=10))};
}

/* ============================================================
   14) Spline 배경 (SPEC §2 spline) — CasaXPS 'Cubic Spline' 과 같은 목적
   매듭: 끝점 둘(endpoints 평균, bg.ends !== false) + 사용자 매듭 {e, y?}. y 가 없으면 가장 가까운 점 ±kh(기본 2) 점 평균.
   자연 3차 spline(양 끝 2차 도함수 0), 끝 매듭 밖은 끝 기울기로 직선 연장.
   ============================================================ */
function splineKnots(E,y,bg){
  const N=E.length, out=[]; bg=bg||{};
  if(bg.ends!==false){ const p=endpoints(E,y,bg.n||5); out.push({e:p.Elo,y:p.ylo,end:true},{e:p.Ehi,y:p.yhi,end:true}); }
  const hw=bg.kh!=null?Math.max(0,Math.round(bg.kh)):2;
  (bg.knots||[]).forEach(k=>{ if(!k||!isFinite(k.e)) return; let yv=k.y;
    if(yv==null||!isFinite(yv)){ let i=0, d=Infinity; for(let j=0;j<N;j++){ const dd=Math.abs(E[j]-k.e); if(dd<d){ d=dd; i=j; } }
      const a=Math.max(0,i-hw), b=Math.min(N-1,i+hw); let s=0; for(let j=a;j<=b;j++) s+=y[j]; yv=s/(b-a+1); }
    out.push({e:+k.e,y:+yv,auto:k.y==null}); });
  out.sort((a,b)=>a.e-b.e);
  const u=[]; out.forEach(k=>{ const L=u[u.length-1]; if(L&&Math.abs(L.e-k.e)<1e-9){ if(L.end) u[u.length-1]=k; } else u.push(k); });
  return u;
}
/* 자연 3차 spline 의 마디 2차 도함수 M (M_0 = M_{K−1} = 0, 삼중대각 Thomas) */
function naturalSplineM(xs,ys){
  const K=xs.length, M=new Float64Array(K); if(K<3) return M;
  const a=new Float64Array(K), b=new Float64Array(K), c=new Float64Array(K), d=new Float64Array(K);
  for(let i=1;i<K-1;i++){ const h0=xs[i]-xs[i-1], h1=xs[i+1]-xs[i]; a[i]=h0; b[i]=2*(h0+h1); c[i]=h1; d[i]=6*((ys[i+1]-ys[i])/h1-(ys[i]-ys[i-1])/h0); }
  for(let i=2;i<K-1;i++){ const w=a[i]/b[i-1]; b[i]-=w*c[i-1]; d[i]-=w*d[i-1]; }
  for(let i=K-2;i>=1;i--) M[i]=(d[i]-(i<K-2?c[i]*M[i+1]:0))/b[i];
  return M;
}
function splineAt(xs,ys,M,t){
  const K=xs.length; if(!K) return 0; if(K===1) return ys[0];
  if(t<=xs[0]){ const h=xs[1]-xs[0]; const s=(ys[1]-ys[0])/h-h*(2*M[0]+M[1])/6; return ys[0]+s*(t-xs[0]); }
  if(t>=xs[K-1]){ const h=xs[K-1]-xs[K-2]; const s=(ys[K-1]-ys[K-2])/h+h*(M[K-2]+2*M[K-1])/6; return ys[K-1]+s*(t-xs[K-1]); }
  let lo=0, hi=K-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(xs[m]<=t) lo=m; else hi=m; }
  const h=xs[hi]-xs[lo], A=(xs[hi]-t)/h, Bq=(t-xs[lo])/h;
  return A*ys[lo]+Bq*ys[hi]+((A*A*A-A)*M[lo]+(Bq*Bq*Bq-Bq)*M[hi])*h*h/6;
}
function bgSpline(E,y,bg){
  const ks=splineKnots(E,y,bg), xs=ks.map(k=>k.e), ys=ks.map(k=>k.y), M=naturalSplineM(xs,ys);
  const B=Float64Array.from(E,e=>splineAt(xs,ys,M,e)); B.knots=ks; return B;
}

/* ============================================================
   15) Simplex (Nelder–Mead) 최적화 (SPEC §14) — CasaXPS 의 Marquardt/Simplex 선택과 같은 목적
   변수: MINUIT 변환 u 를 매개변수마다 척도 s_j 로 나눈 z (s_j = 전형적 걸음 dp_j 에 해당하는 u 변화). 계수 ρ=1, χ=2, γ=½, σ=½ (Nelder & Mead 1965; Lagarias et al. 1998).
   정지: max|z_i − z_0| ≤ xtol 그리고 max|f_i − f_0| ≤ ftol·max(1,|f_0|). 수렴 뒤 같은 점에서 새 simplex 로 다시 시작해 개선이 없을 때까지(최대 6 회).
   ============================================================ */
function nmFit(prob,map,v,opt){
  opt=opt||{}; const E=prob.E, y=prob.y, n=E.length, np=map.free.length; const ev=makeEvaluator(prob,map);
  const sw=Float64Array.from(v,x=>1/Math.sqrt(x));
  const lo=map.free.map(f=>f.lo), hi=map.free.map(f=>f.hi);
  const p0=getVec(prob,map).map((x,i)=>Math.min(hi[i],Math.max(lo[i],x)));
  const u0=p0.map((x,i)=>toU(x,lo[i],hi[i]));
  const chi=mu=>{ let s=0; for(let i=0;i<n;i++){ const r=(y[i]-mu[i])*sw[i]; s+=r*r; } return s; };
  setVec(prob,map,p0); let cur=ev.full();
  if(!np) return {chi2:chi(cur.mu),iter:0,conv:true,stop:"converged",mu:cur.mu,B:cur.B,comps:ev.comps()};
  let ymax=0; for(const t of y) ymax=Math.max(ymax,Math.abs(t));
  const s=map.free.map((f,j)=>{ const p=p0[j]; let dp;
    if(f.k==="be"||f.k==="dE") dp=0.1; else if(Math.abs(p)>1e-8) dp=0.1*Math.abs(p); else if(isFinite(lo[j])&&isFinite(hi[j])) dp=0.05*(hi[j]-lo[j]); else dp=1e-3*Math.max(1,ymax);
    let pu=Math.min(hi[j],p+dp); if(pu===p) pu=Math.max(lo[j],p-dp); let d=Math.abs(toU(pu,lo[j],hi[j])-u0[j]);
    if(!(d>1e-9)) d=Math.abs(toU(Math.max(lo[j],p-dp),lo[j],hi[j])-u0[j]); return d>1e-9?d:1e-3; });
  let nfev=0; const maxfev=opt.maxfev||Math.max(4000,1500*(np+1));
  const F=z=>{ nfev++; const p=z.map((zi,j)=>fromU(u0[j]+s[j]*zi,lo[j],hi[j])); setVec(prob,map,p); const X=chi(ev.full().mu); return isFinite(X)?X:Infinity; };
  const xtol=opt.xtol||1e-7, ftol=opt.ftol||1e-12;
  let best=new Array(np).fill(0), fbest=F(best), conv=false;
  for(let rs=0;rs<6&&nfev<maxfev;rs++){
    const step=rs===0?1:0.2;
    let sim=[best.slice()], fs=[fbest];
    for(let j=0;j<np;j++){ const z=best.slice(); z[j]+=step; sim.push(z); fs.push(F(z)); }
    let ok=false;
    while(nfev<maxfev){
      const ord=fs.map((f,i)=>i).sort((a,b)=>fs[a]-fs[b]); sim=ord.map(i=>sim[i]); fs=ord.map(i=>fs[i]);
      let dx=0, df=0; for(let i=1;i<=np;i++){ df=Math.max(df,Math.abs(fs[i]-fs[0])); for(let j=0;j<np;j++) dx=Math.max(dx,Math.abs(sim[i][j]-sim[0][j])); }
      if(dx<=xtol&&df<=ftol*Math.max(1,Math.abs(fs[0]))){ ok=true; break; }
      const xb=new Array(np).fill(0); for(let i=0;i<np;i++) for(let j=0;j<np;j++) xb[j]+=sim[i][j]/np;
      const W=sim[np], at=c=>xb.map((b,j)=>b+c*(b-W[j]));
      const xr=at(1), fr=F(xr);
      if(fr<fs[0]){ const xe=at(2), fe=F(xe); if(fe<fr){ sim[np]=xe; fs[np]=fe; } else { sim[np]=xr; fs[np]=fr; } continue; }
      if(fr<fs[np-1]){ sim[np]=xr; fs[np]=fr; continue; }
      let shrink=false;
      if(fr<fs[np]){ const xc=at(0.5), fc=F(xc); if(fc<=fr){ sim[np]=xc; fs[np]=fc; } else shrink=true; }
      else { const xcc=at(-0.5), fcc=F(xcc); if(fcc<fs[np]){ sim[np]=xcc; fs[np]=fcc; } else shrink=true; }
      if(shrink){ for(let i=1;i<=np;i++){ sim[i]=sim[i].map((z,j)=>sim[0][j]+0.5*(z-sim[0][j])); fs[i]=F(sim[i]); } }
    }
    const imp=fbest-fs[0]; best=sim[0].slice(); fbest=fs[0];
    if(ok&&rs>0&&imp<=ftol*Math.max(1,Math.abs(fbest))*10){ conv=true; break; }
    if(!ok) break;
  }
  const p=best.map((zi,j)=>fromU(u0[j]+s[j]*zi,lo[j],hi[j])); setVec(prob,map,p); cur=ev.full();
  return {chi2:chi(cur.mu),iter:nfev,conv,stop:nfev>=maxfev?"maxIter":"converged",mu:cur.mu,B:cur.B,comps:ev.comps(),optimizer:"simplex"};
}

/* ============================================================
   16) 여러 스펙트럼에 같은 모형 (SPEC §15) · 정량 Monte Carlo (SPEC §4)
   ============================================================ */
/* 모형(보정 BE 기준) → 다른 스펙트럼: 성분 BE·free BE 한계·spline 매듭을 Δ_src − Δ_dst 만큼, 창도 같은 보정 BE 로(스펙트럼 끝에서 자름) */
function propagateModel(model,shSrc,shDst,beRange){
  const d=(+shSrc||0)-(+shDst||0); const m=shiftModelBE(model.comps,model.cons,d);
  const bg=clone(model.bg||{type:"smart",n:5}); if(bg.knots) bg.knots.forEach(k=>{ k.e+=d; });
  if(bg.type==="active"){ bg.c0=null; bg.k=null; }
  let lo=model.lo+d, hi=model.hi+d; if(beRange){ lo=Math.max(beRange[0],lo); hi=Math.min(beRange[1],hi); }
  const out={lo,hi,bg,comps:m.comps,cons:m.cons,weight:model.weight||"poisson"}; if(model.optimizer==="simplex") out.optimizer="simplex"; return out;
}
/* 한 스펙트럼 피팅 작업(Worker): 넓이는 NNLS 로 다시 시작 */
function batchJob(spec,model,opt){ const P=makeProblem(spec,model); const r=fitOnce(P,opt||{}); const out=finalize(P,r); out.lo=model.lo; out.hi=model.hi; return out; }
/* at% Monte Carlo: 넓이를 N(A, σ_A) 에서 M 번 뽑아 at% 분포 → 평균·sd·2.5/97.5 백분위. 검출 안 됨(nd)은 기준값으로 정한 집합을 고정 */
function quantMC(rows,opts,M,seed){
  const base=quantify(rows,opts); const R=rng(seed||7); M=M||2000; const use=base.map(r=>!r.nd&&isFinite(r.nrm)&&r.den>0);
  const A=base.map(()=>new Float64Array(M)); let ok=0;
  for(let m=0;m<M;m++){ let tot=0; const nn=base.map((r,i)=>{ const g=R.g(); if(!use[i]) return 0; const a=r.area+(r.sArea>0?r.sArea*g:0); const v=a/r.den; tot+=v; return v; });
    if(!(tot>0)) continue; base.forEach((r,i)=>{ if(use[i]) A[i][ok]=100*nn[i]/tot; }); ok++; }
  const pct=(s,q)=>{ const t=(s.length-1)*q, k=Math.floor(t); return s[k]+(s[Math.min(k+1,s.length-1)]-s[k])*(t-k); };
  return base.map((r,i)=>{ if(!use[i]||!ok) return null; const a=Array.from(A[i].subarray(0,ok)).sort((x,y)=>x-y); const mu=a.reduce((s,x)=>s+x,0)/ok;
    const sd=Math.sqrt(a.reduce((s,x)=>s+(x-mu)*(x-mu),0)/Math.max(1,ok-1)); return {mean:mu,sd,lo:pct(a,0.025),hi:pct(a,0.975),n:ok}; });
}

/* ============================================================
   보조: 문제 만들기 (스펙트럼 + 영역 + 배경 + 성분)
   ============================================================ */
function sliceRegion(be,y,lo,hi){ const i0=be.findIndex(e=>e>=lo); let i1=-1; for(let i=be.length-1;i>=0;i--) if(be[i]<=hi){ i1=i; break; } if(i0<0||i1<i0) return null; return {E:be.slice(i0,i1+1),y:y.slice(i0,i1+1),i0,i1}; }
function makeProblem(spec,st){
  const r=sliceRegion(spec.be,spec.y,st.lo,st.hi); if(!r||r.E.length<8) throw new Error("영역 안 데이터 점이 8개 미만입니다");
  const P={E:Float64Array.from(r.E),y:Float64Array.from(r.y),c:spec.c||null,bg:clone(st.bg),comps:clone(st.comps),cons:clone(st.cons||{}),i0:r.i0,weight:st.weight==="const"?"const":st.weight==="casa"?"casa":"poisson"};
  if(st.optimizer==="simplex") P.optimizer="simplex";
  if(P.bg&&P.bg.conv==="casa") P.bg._ends=casaEnds(spec.be,spec.y,r.i0,r.i1,P.bg.aw!=null?P.bg.aw:1);
  if(P.bg&&P.bg.conv==="avantage"&&(P.bg.type==="smart"||P.bg.type==="shirley")) P.bg._pad=avantagePad(spec.y,r.i0,r.i1,avantageHalf(P.bg,r.E));   // v1.9.2: 끝값 ±h 평균에 창 밖 h 채널
  return P;
}
function newComp(shape,be,o){
  o=o||{}; const c={id:o.id||("c"+Math.random().toString(36).slice(2,8)),name:o.name||"",shape,be,area:o.area||1,el:o.el||"",line:o.line||"",doublet:null,color:o.color||null};
  const F=o.F||1.2;
  if(shape==="GLS"||shape==="GLP"){ c.F=F; c.m=o.m!=null?o.m:0.3; }
  if(shape==="V"){ c.wG=o.wG!=null?o.wG:F*0.85; c.wL=o.wL!=null?o.wL:F*0.25; }
  if(shape==="DS"){ c.alpha=o.alpha!=null?o.alpha:0.1; c.wL=o.wL!=null?o.wL:0.3; c.wG=o.wG!=null?o.wG:F*0.7; }
  if(shape==="LA"){ c.F=o.F||F*0.8; c.alpha=o.alpha!=null?o.alpha:1; c.beta=o.beta!=null?o.beta:1.5; c.wG=o.wG!=null?o.wG:F*0.4; }
  if(shape==="GLT"){ c.F=F; c.m=o.m!=null?o.m:0.3; c.T=o.T!=null?o.T:0.2; c.tau=o.tau!=null?o.tau:1.0; }
  if(shape==="MEAS"){ c.tpl=o.tpl||null; c.wG=o.wG!=null?o.wG:0; }
  if(shape==="AGT"){ c.F=F; c.m=o.m!=null?o.m:0.3; c.ET=o.ET!=null?o.ET:0.05; c.TM=o.TM!=null?o.TM:0; c.CT=o.CT!=null?o.CT:0; c.ch=o.ch!=null?o.ch:0.05; if(o.tu) c.tu=o.tu; }
  if(shape==="TLA"){ c.F=o.F||F*0.8; c.alpha=o.alpha!=null?o.alpha:1; c.mu=o.mu!=null?o.mu:2; c.wG=o.wG!=null?o.wG:F*0.4; }
  if(shape==="GLN"||shape==="SGLN"){ c.F=F; c.m=o.m!=null?o.m:0.3; c.wG=o.wG!=null?o.wG:0; if(o.tail==="T"||o.tail==="PHI"){ c.tail=o.tail; c.tk=o.tk!=null?o.tk:(o.tail==="T"?2:2); if(o.tail==="PHI") c.tj=o.tj!=null?o.tj:1; } }
  if(shape==="TP"){ c.F=o.F||F; c.ta=o.ta!=null?o.ta:1; c.tb=o.tb!=null?o.tb:1; c.wL=o.wL!=null?o.wL:0; }
  if(o.casaN!=null) c.casaN=o.casaN;
  if(shape==="DS"&&o.dsl>0) c.dsl=o.dsl;
  if((o.mod==="ST"||o.mod==="SB")&&S2_MODBASE[shape]){ c.mod=o.mod; if(o.mod==="ST"){ c.stMu=o.stMu!=null?o.stMu:1; c.stG=o.stG!=null?o.stG:0; if(o.stN!=null) c.stN=o.stN; } else c.sbArg=o.sbArg||0; }
  if(o.doublet) c.doublet={dE:o.doublet.dE,r:o.doublet.r,q:o.doublet.q!=null?o.doublet.q:1};
  return c;
}
/* 모양 바꾸기: 현재 폭을 최대한 보존 */
function convertShape(c,shape){
  if(c.shape===shape) return c;
  if(shape==="MEAS") throw new Error("측정 선모양은 기준 스펙트럼에서만 만들 수 있습니다");
  const tf=trueFWHM(c.shape,shapeParams(c,1)).fwhm||1.2;
  const n=newComp(shape,c.be,{F:tf,id:c.id,name:c.name,el:c.el,line:c.line,color:c.color}); n.area=c.area; if(c.doublet) n.doublet=Object.assign({},c.doublet);
  if((shape==="GLS"||shape==="GLP"||shape==="GLT"||shape==="GLN")&&c.m!=null) n.m=c.shape==="SGLN"?sglToGLS(100*c.m):c.m;   // SGLN 의 m = 높이 비
  if(shape==="SGLN"&&c.m!=null) n.m=c.shape==="SGLN"?c.m:c.shape==="GLS"?glsToSGL(c.m)/100:c.m;
  if((shape==="GLN"||shape==="SGLN")&&compTail(c)){ n.tail=c.tail; n.tk=c.tk; if(c.tj!=null) n.tj=c.tj; }
  if((shape==="LA"||shape==="TLA")&&(c.shape==="LA"||c.shape==="TLA")){ n.alpha=c.alpha; }
  if(compMod(c)&&S2_MODBASE[shape]){ n.mod=c.mod; ["stMu","stG","stN","sbArg"].forEach(k=>{ if(c[k]!=null) n[k]=c[k]; }); }
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
  return {comps:P.comps,bg:P.bg,chi2:r.chi2,iter:r.iter,conv:r.conv,stop:r.stop||null,stall:r.stall||null,branch:r.branch||null,mu:Array.from(r.mu),B:Array.from(r.B),v:Array.from(r.v),
    stats:Object.assign({},st,{resid:Array.from(st.resid),stop:r.stop||null,stall:r.stall||null}),   // stop: "converged" | "stalled" | "maxIter" (SPEC §3.3b)
    cov:{sig:cov.sig,corr:Array.from(cov.corr),C:Array.from(cov.C),atB:cov.atB,phi:cov.phi,dof:cov.dof,singular:cov.singular,names:cov.names,np:cov.np},
    summ,multi:r.multi||null,free:map.free.map(f=>({id:f.id,k:f.k,lo:f.lo,hi:f.hi})),
    sumSig:{all:areaSumSigma(P,map,cov,P.comps.map((c,i)=>i))},
    lcf:lcfDiag(P,r.v),optimizer:P.optimizer==="simplex"?"simplex":"lm"};
}
/* 한 묶음의 시작점에서 피팅 → 최선 해 (Worker 작업 단위) */
function fitJob(prob,starts,opt){
  let best=null; const chis=[];
  (starts&&starts.length?starts:[null]).forEach(pt=>{ const P=clone(stripProb(prob)); P.E=prob.E; P.y=prob.y; applyStart(P,pt); try{ const r=fitOnce(P,opt); chis.push(r.chi2p1); if(!best||r.chi2p1<best.r.chi2p1) best={P,r}; }catch(e){ chis.push(Infinity); } });
  if(!best) throw new Error("모든 시작점에서 피팅 실패");
  return {comps:best.P.comps,bg:best.P.bg,chi2:best.r.chi2p1,chis,iter:best.r.iter,stop:best.r.stop||null,stall:best.r.stall||null,branch:best.r.branch||null};   // chi2 = 1차 χ² (시작점 비교용)
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


/* ============================================================
   16) CasaXPS 표기 ↔ INNO 성분·배경 (SPEC_casa_v2.md — CasaXPS 2.3.27 V2 실측으로 정함; CasaXPS 코드는 쓰지 않음)
   ============================================================ */
const CASA_AL=PI/2, CASA_AG=Math.sqrt(PI/(4*LN2));   // ∫ 1/(1+4t²) dt, ∫ exp(−4ln2 t²) dt (F = 1)
/* SGL(p) = (1 − p/100)·exp(−4ln2 x²/F²) + (p/100)/(1 + 4x²/F²) (높이 혼합, 같은 F) ≡ GLS(F, m) 넓이 혼합 — 같은 함수, m ↔ p 만 다름 */
function sglToGLS(p){ const h=Math.min(1,Math.max(0,p/100)); return h*CASA_AL/(h*CASA_AL+(1-h)*CASA_AG); }
function glsToSGL(m){ m=Math.min(1,Math.max(0,m)); const a=m/CASA_AL, b=(1-m)/CASA_AG; return a+b>0?100*a/(a+b):0; }
/* F = 1 일 때의 실제 FWHM (척도 불변 모양의 폭 인자) */
function casaWidthFactor(shape,o){ const c=newComp(shape,0,Object.assign({F:1},o)); c.F=1; if(o.wG!=null) c.wG=o.wG; if(o.casaN!=null) c.casaN=o.casaN; return trueFWHM(shape,shapeParams(c,1)).fwhm; }
/* CasaXPS 'fwhm' 숫자와 비대칭 지수 (SPEC_casa_v4 §3, CasaXPS 2.3.27 Round 1–4 실측 14 모양): CasaXPS 의 최댓값 x̂ = 성분 원점에 맞춘 격자 k·d (d = F/40,
   DS 는 wL/40) 에서 가장 큰 마디와 이웃 두 마디를 지나는 포물선의 꼭짓점, 반높이 = 꼭짓점 값/2, a < x̂ < b 반높이 점 (a 낮은 BE 쪽),
   l = x̂ − a, h = b − x̂ → fwhm = 2·min(l, h), asym = 1 − l/h (낮은 BE 쪽이 넓으면 음수). 꺾인 꼭대기(LA α≠β n 0, GL/SGL T·PHI)의 'δ' 는
   이 포물선 꼭짓점이 꼬리 쪽으로 밀리는 것(LA: δ ≈ d(β−α)/(2(α+β))). 대칭: 실제 FWHM, 0. ST: 같은 규칙, 위치·폭 미정 (exact:false). */
const CASA_GRID=1/40, CASA_DS_T=600, CWA=new Map();
/* 격자 FWHM(trueFWHM 일반 경로: 간격 0.005·폭, 선형 보간)의 최댓값·반높이 점을 다듬음 — 비대칭 지수가 L/H 에 민감해서 (LA(1.53,2.43): H 오차 1e-4 → asym 4e-4).
   같은 창(ctx)에서 값을 매기므로 창 정규화 상수는 상쇄된다 */
function refineHalfMax(shape,p,w){
  if(w.fwhm==null||w.lo==null||w.hi==null) return w; const W=w.fwhm, ctx={lo:-60*W,hi:60*W};
  const ev=x=>shapeUnit(shape,Float64Array.of(x),p,ctx)[0], gr=(Math.sqrt(5)-1)/2;
  let a=w.xmax-0.02*W, b=w.xmax+0.02*W, x1=b-gr*(b-a), x2=a+gr*(b-a), f1=ev(x1), f2=ev(x2);
  for(let k=0;k<60&&b-a>1e-12*W;k++){ if(f1>=f2){ b=x2; x2=x1; f2=f1; x1=b-gr*(b-a); f1=ev(x1); } else { a=x1; x1=x2; f1=f2; x2=a+gr*(b-a); f2=ev(x2); } }
  let xm=(a+b)/2, ym=ev(xm); const y0=ev(w.xmax); if(y0>ym){ xm=w.xmax; ym=y0; }
  const half=ym/2, cross=(u,v)=>{ let fu=ev(u)-half, fv=ev(v)-half; if(!(fu*fv<0)) return null; for(let k=0;k<60&&Math.abs(v-u)>1e-12*W;k++){ const m=(u+v)/2, fm=ev(m)-half; if(fu*fm<=0){ v=m; fv=fm; } else { u=m; fu=fm; } } return (u+v)/2; };
  const lo=cross(w.lo-0.02*W,Math.min(w.lo+0.02*W,xm)), hi=cross(Math.max(w.hi-0.02*W,xm),w.hi+0.02*W);
  return lo!=null&&hi!=null?{fwhm:hi-lo,xmax:xm,ymax:ym,lo,hi}:w; }
function casaWidthAsym(c){ const key=c.shape+"|"+JSON.stringify(shapeParams(c,1)); const hit=CWA.get(key); if(hit) return Object.assign({},hit);
  const res=casaWidthAsym0(c); if(CWA.size>256) CWA.clear(); CWA.set(key,res); return Object.assign({},res); }
function casaWidthAsym0(c){
  const md=compMod(c), sp=shapeParams(c,1), Fb=c.shape==="V"?sp.wL:sp.F;
  // ST: CasaXPS tabulates S(μ)F on nodes k·D, D = (base FWHM)/40, with the step lagged by D/2, and puts its position on the largest node (SPEC_casa_v4 §4)
  const D=md==="ST"?(c.stCasaK>0?c.stCasaK:casaStK(c))*Fb/40:0; if(md==="ST") sp.stLag=D/2;
  const w0=trueFWHM(c.shape,sp);
  if(w0.fwhm==null||w0.lo==null||w0.hi==null) return {fwhm:null,asym:null,exact:false,lo:null,hi:null,xmax:null,node:null,side:null,note:"no half-maximum points (CasaXPS fwhm not defined here)"};
  const W=w0.fwhm, F=(c.shape==="DS"?sp.wL:sp.F)||W, d=md==="ST"?D:CASA_GRID*F, ctx={lo:-60*W,hi:60*W};   // 같은 ctx → 창 정규화 상수는 상쇄
  const ev=xs=>shapeUnit(c.shape,Float64Array.from(xs),sp,ctx);
  let k=Math.round(w0.xmax/d), v=ev([(k-1)*d,k*d,(k+1)*d]);
  for(let it=0;it<16;it++){ if(v[2]>v[1]) k++; else if(v[0]>v[1]) k--; else break; v=ev([(k-1)*d,k*d,(k+1)*d]); }
  const den=v[0]-2*v[1]+v[2], xs=den<0?(k+0.5*(v[0]-v[2])/den)*d:k*d, half=(den<0?v[1]-0.125*(v[0]-v[2])*(v[0]-v[2])/den:v[1])/2;
  const cross=(u,w)=>{ for(let pass=0;pass<3;pass++){ const K=64, g=new Float64Array(K+1); for(let i=0;i<=K;i++) g[i]=u+(w-u)*i/K; const f=ev(g); let j=-1;
      for(let i=0;i<K;i++) if((f[i]-half)*(f[i+1]-half)<=0){ j=i; break; } if(j<0) return null; if(pass===2) return g[j]+(f[j]-half)/(f[j]-f[j+1])*(g[j+1]-g[j]); u=g[j]; w=g[j+1]; } return null; };
  const lo=cross(w0.lo-0.03*W,Math.min(w0.lo+0.03*W,xs)), hi=cross(Math.max(w0.hi-0.03*W,xs),w0.hi+0.03*W);
  if(lo==null||hi==null) return {fwhm:null,asym:null,exact:false,lo:null,hi:null,xmax:null,node:null,side:null,note:"half-maximum points not found"};
  const l=xs-lo, h=hi-xs, sym=Math.abs(h-l)<=1e-7*(h+l);
  const stConv=md==="ST"&&(sp.wG>1e-6||sp.stG>1e-9);
  const note=md==="ST"?"ST: CasaXPS shape = S(μ)F with the step lagged by (base FWHM)/80; position = INNO origin + largest node k·(base FWHM)/40 (measured, 16 readings ≤ 0.2 meV); fwhm = 2·min(l,h), asym = 1 − l/h"+(stConv?" (Gaussian-broadened base: −0.5…+0.3 % — approximate)":" (measured ≤ 4e-4 / 1.1e-3)")
    :sym?"CasaXPS fwhm = true FWHM"
    :"CasaXPS fwhm = 2·min(l, h), asym = 1 − l/h; maximum = parabola vertex on nodes k·F/40 (measured: 14 shapes, fwhm ≤ 5.3e-4, asym ≤ 6e-4)";
  return {fwhm:sym?w0.fwhm:2*Math.min(l,h),asym:sym?0:1-l/h,side:l<=h?"low":"high",exact:!stConv,lo,hi,xmax:xs,node:k*d,note};
}
/* ST: base FWHM / base width (F, or wL for V) — scale-free, fixed by the base's shape parameters (m, α, β, n …) */
function casaStK(c){ const b=clone(c); b.mod=null; delete b.stCasaK; const sk=b.shape==="V"?"wL":"F"; b[sk]=1; const w=trueFWHM(b.shape,shapeParams(b,1)).fwhm; return w>0?w:1; }
/* CasaXPS DS area = amplitude × ∫ P over E0 ± 600·wL (SPEC_casa_v4 §2, α 0 / 0.1 / 0.2, n 0 / 100, two windows: ≤ 4e-4); INNO DS area = area inside [emin, emax]
   → CasaXPS area = INNO area × factor. Raw-DS integrals (closed form); the n-Gaussian changes the window integral by ≤ 5e-5. Cut variant (dsl): null (not read). */
function dsCasaAreaFactor(c,emin,emax){ if(!c||c.shape!=="DS"||c.dsl>0||!(emax>emin)) return null; const sp=shapeParams(c,1), g=sp.wL/2, al=c.alpha, T=CASA_DS_T*sp.wL;
  const N600=dsPrim(T,al,g)-dsPrim(-T,al,g), Nw=dsPrim(emax-c.be,al,g)-dsPrim(emin-c.be,al,g); return Nw>0?N600/Nw:null; }
/* INNO 성분 → CasaXPS 표기 + CasaXPS 'fwhm' 칸에 들어갈 수. exact = 실측으로 같은 수가 확인된 경우만 */
function compToCasa0(c){
  const sp=shapeParams(c,1), r=(v,d)=>+(+v).toFixed(d);
  if(c.shape==="GLP") return {notation:`GL(${r(100*c.m,2)})`,fwhm:trueFWHM("GLP",sp).fwhm,exact:true,note:"CasaXPS GL(p) = GLP, m = p/100; CasaXPS fwhm = true FWHM (−1.7e-4 measured)"};
  if(c.shape==="GLS") return {notation:`SGL(${r(glsToSGL(c.m),2)})`,fwhm:c.F,exact:true,note:"CasaXPS SGL(p) = height mix with one F = GLS with m = sglToGLS(p); fwhm = F (= true FWHM)"};
  if(c.shape==="LA"){ const sym=Math.abs(c.alpha-c.beta)<1e-9; const n=c.casaN!=null?+c.casaN:(sp.wG>1e-9?sp.wG/(CASA_LA_RHO*sp.F):0), ni=Math.round(n);
    const not=`LA(${r(c.alpha,3)},${r(c.beta,3)},${ni})`;
    if(!sym){ const cw=casaWidthAsym(c); return {notation:not,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,n,exact:Math.abs(n-ni)<0.05,note:"asymmetric LA: shape, tail side and n rule measured (LA(1.53,2.43,0), LA(2.43,1.53,0), LA(2,3,0), LA(2,3,30)); "+cw.note}; }
    return {notation:not,fwhm:trueFWHM("LA",sp).fwhm,n,exact:Math.abs(n-ni)<0.05&&(c.alpha===1||ni===0),note:"symmetric LA: CasaXPS fwhm = true FWHM; Gaussian FWHM = n·F/144.87 (measured for α = β = 1 only)"}; }
  if(c.shape==="AGT"&&c.tu==="fwhm"&&!(c.TM>0)){ const cw=casaWidthAsym(c); return {notation:`GL(${r(100*c.m,2)})T(${r(c.ET,4)})`,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,exact:true,note:"same shape as CasaXPS GL(p)T(k) (area/position measured); "+cw.note+" (GL(30)T(2), F 1.2: 1.18051 read)"}; }
  return {notation:null,fwhm:null,exact:false,note:"no CasaXPS line shape with the same definition"};
}
/* CasaXPS 표기('GL(30)', 'SGL(30)', 'LA(1.53,2.43,0)', 'LA(1,250)', 'LA(50)', 'GL(30)T(2)') + CasaXPS fwhm → INNO 성분 */
function casaToComp0(str,o){
  o=o||{}; const s=String(str||"").replace(/\s+/g,""), be=o.be!=null?+o.be:0, area=o.area!=null?+o.area:1, fw=o.fwhm>0?+o.fwhm:1.2; let m;
  const mk=(shape,q)=>{ const c=newComp(shape,be,Object.assign({area},q)); if(o.id) c.id=o.id; return c; };
  if((m=s.match(/^GL\(([\d.]+)\)T\(([\d.]+)\)$/i))){ const c=mk("AGT",{F:1,m:+m[1]/100,ET:+m[2],TM:0,CT:0,tu:"fwhm"}); const u=o.innoUnits?trueFWHM("AGT",shapeParams(c,1)).fwhm:casaWidthAsym(c).fwhm; c.F=u>0?fw/u:fw;
    return {comp:c,exact:false,note:"F from the CasaXPS fwhm number via an empirical peak offset (SPEC_casa_v3 §3: CasaXPS fwhm ≈ 0.984·F for GL(30)T(2)) — refit the width"}; }
  if((m=s.match(/^SGL\(([\d.]+)\)$/i))){ const c=mk("GLS",{F:fw,m:sglToGLS(+m[1])}); return {comp:c,exact:true,note:"SGL(p) = GLS with m = "+c.m.toFixed(4)}; }
  if((m=s.match(/^GL\(([\d.]+)\)$/i))){ const mm=+m[1]/100, c=mk("GLP",{F:1,m:mm}); c.F=fw/casaWidthFactor("GLP",{m:mm}); return {comp:c,exact:true,note:"GL(p) = GLP; F from the true FWHM"}; }
  let al,bt,n;
  if((m=s.match(/^LA\(([\d.]+)\)$/i))){ al=bt=1; n=1401-14.01*(+m[1]); }   // Fairley et al. 2023 SI Eq. 6
  else if((m=s.match(/^LA\(([\d.]+),([\d.]+)\)$/i))){ if(+m[2]>=20&&Number.isInteger(+m[2])){ al=bt=+m[1]; n=+m[2]; } else { al=+m[1]; bt=+m[2]; n=0; } }   // 'LA(1,1283)' (CasaXPS 2.3.25 review) = LA(α=β, n); 'LA(1.53,2.43)' = LA(α,β,0) — 해석(미확인)
  else if((m=s.match(/^LA\(([\d.]+),([\d.]+),([\d.]+)\)$/i))){ al=+m[1]; bt=+m[2]; n=+m[3]; }
  if(al!=null){ const c=mk("LA",{F:1,alpha:al,beta:bt,wG:0}); c.casaN=n; const sym=Math.abs(al-bt)<1e-9;
    let cw=null; if(sym||o.innoUnits) c.F=fw/casaWidthFactor("LA",{alpha:al,beta:bt,casaN:n,wG:0})||fw; else { const u0=clone(c); u0.F=1; cw=casaWidthAsym(u0); c.F=cw.fwhm>0?fw/cw.fwhm:fw; }
    return {comp:c,exact:sym?(al===1||n===0):!o.innoUnits&&!!(cw&&cw.exact),note:sym?"CasaXPS fwhm = true FWHM; Gaussian FWHM = n·F/144.87 (wG fixed by n)":o.innoUnits?"asymmetric LA: width kept (true FWHM)":"asymmetric LA: F from the CasaXPS fwhm number (SPEC_casa_v4 §3), "+cw.note}; }
  return null;
}
/* ---- CasaXPS Cookbook 2026 표기 (SPEC_shapes2 §4): [ST(μ[,γ]) | SB(a)] 바탕 [T(k) | PHI(j,k)]
   바탕: GL(m) · GL(m,n) · SGL(m) · SGL(m,n) · LA(…) · TLA(α,μ[,n]) · DS(α,n[,l]) · TP(α,β[,n]). n(정수) → 폭 = n·F/144.87 (LA 실측 규칙을 가정) ---- */
const CASA_SHAPES_HELP="GL(m[,n]), SGL(m[,n]) with optional T(k) or PHI(j,k); LA(m)|LA(α,n)|LA(α,β,n); TLA(α,μ[,n]); DS(α,n[,l]); TP(α,β[,n]); prefixes ST(μ[,γ]) and SB(0)";
function casaTokens(s){ const out=[], re=/([A-Za-z]+)\(([^()]*)\)/y; let i=0; while(i<s.length){ re.lastIndex=i; const m=re.exec(s); if(!m) return null; const a=m[2]===""?[]:m[2].split(",").map(Number); if(a.some(v=>!isFinite(v))) return null; out.push({n:m[1].toUpperCase(),a}); i=re.lastIndex; } return out.length?out:null; }
/* 단위 폭(F 또는 DS 의 wL = 1)에서의 실제 FWHM — 척도 불변 모양 (n 규칙 폭 포함) */
function unitFWHM(c){ const u=clone(c); u.mod=null; const sk=u.shape==="DS"?"wL":"F"; u[sk]=1; return trueFWHM(u.shape,shapeParams(u,1)).fwhm; }
function casaToComp(str,o){
  o=o||{}; const raw=String(str||""), s=raw.replace(/\s+/g,"");
  const fail=msg=>{ if(o.strict) throw new Error(msg||`CasaXPS line shape '${raw}' not recognised — supported: ${CASA_SHAPES_HELP}`); return null; };
  const T=casaTokens(s); if(!T) return fail();
  let mod=null; if(T[0].n==="ST"||T[0].n==="SB"){ mod=T.shift(); if(!T.length) return fail(`'${raw}': ${mod.n} needs a line shape after it (e.g. ${mod.n==="ST"?"ST(1)LA(1,143)":"SB(0)LA(1,43)"})`); }
  if(mod&&(mod.n==="ST"?!(mod.a.length===1||mod.a.length===2)||!(mod.a[0]>0)||(mod.a[1]!=null&&!(mod.a[1]>=0)):mod.a.length>1)) return fail(`'${raw}': ${mod.n}(${mod.a.join(",")}) — ST(μ) or ST(μ,γ) with μ > 0, γ ≥ 0; SB(0)`);
  const b=T[0], tl=T[1]||null; if(T.length>2) return fail();
  if(tl&&!((tl.n==="T"&&tl.a.length===1&&tl.a[0]>0)||(tl.n==="PHI"&&tl.a.length===2&&tl.a[1]>0&&tl.a[0]>=0))) return fail(`'${raw}': tail must be T(k), k > 0, or PHI(j,k), k > 0`);
  if(tl&&!(b.n==="GL"||b.n==="SGL")) return fail(`'${raw}': T(k)/PHI(j,k) tails apply to GL or SGL only`);
  let r=null, notes=[];
  const be=o.be!=null?+o.be:0, area=o.area!=null?+o.area:1, fw=o.fwhm>0?+o.fwhm:1.2;
  const mk=(shape,q)=>{ const c=newComp(shape,be,Object.assign({area},q)); if(o.id) c.id=o.id; return c; };
  const bad=()=>fail(`'${raw}': wrong number or range of arguments for ${b.n}(${b.a.join(",")}) — supported: ${CASA_SHAPES_HELP}`);
  if((b.n==="GL"||b.n==="SGL")&&(b.a.length===2||tl)){ const p=b.a[0], n=b.a.length===2?b.a[1]:0; if(!(p>=0&&p<=100&&n>=0)||b.a.length>2) return bad();
    const q={F:1,m:p/100,casaN:n}; if(tl){ q.tail=tl.n; if(tl.n==="T") q.tk=tl.a[0]; else { q.tj=tl.a[0]; q.tk=tl.a[1]; } }
    r={comp:mk(b.n==="GL"?"GLN":"SGLN",q),exact:true,note:(n>0?"Gaussian FWHM = n·F/144.87 (measured: LA(1,1,n), GL(30,100), SGL(30,100))":"n = 0")+(tl?"; tail per Cookbook p. 67 (shape measured, GL(30)T(2))":"; CasaXPS fwhm = true FWHM")}; }
  else if(b.n==="TLA"){ const [al,mu,n]=[b.a[0],b.a[1],b.a[2]!=null?b.a[2]:0]; if(!(b.a.length===2||b.a.length===3)||!(al>0.5&&mu>0&&n>=0)) return fail(`'${raw}': TLA(α,μ[,n]) needs α > 0.5 (normalisable), μ > 0, n ≥ 0`);
    r={comp:mk("TLA",{F:1,alpha:al,mu,wG:0,casaN:n}),exact:true,note:"TLA(α,μ,n) notation accepted by CasaXPS 2.3.27; shape, tail side (high BE, +π/μ) and n rule (Gaussian FWHM = n·F/144.87) measured (TLA(1.2,3,0), TLA(1.2,3,100))"}; }
  else if(b.n==="DS"){ const [al,n,l]=b.a; if(!(b.a.length===2||b.a.length===3)||!(al>=0&&al<1&&n>=0)||(l!=null&&!(l>0))) return fail(`'${raw}': DS(α,n[,l]) needs 0 ≤ α < 1, n ≥ 0, l > 0`);
    r={comp:mk("DS",{alpha:al,wL:1,wG:0,casaN:n,dsl:l}),exact:false,note:"DS(α,n): shape, tail side and n rule (Gaussian FWHM = n·wL/144.87) measured (α 0–0.2, n 0, 100)"+(l?"; Shirley-cut variant Q(x,α,n,l), Cookbook pp. 66–67 (not read)":"")}; }
  else if(b.n==="TP"){ const [ta,tb,n]=[b.a[0],b.a[1],b.a[2]!=null?b.a[2]:0]; if(!(b.a.length===2||b.a.length===3)||!(ta>=0&&tb>=0&&ta+tb>0&&n>=0)) return fail(`'${raw}': TP(α,β[,n]) needs α, β ≥ 0 (not both 0), n ≥ 0`);
    r={comp:mk("TP",{F:fw,ta,tb,wL:0,casaN:n}),exact:false,note:"TP width F = CasaXPS fwhm (assumed); Lorentzian FWHM = n·F/144.87 assumed"}; r.keepF=true; }
  if(!r){ if(tl) return bad(); const r0=casaToComp0(b.n+"("+b.a.join(",")+")",Object.assign({},o,{strict:false})); if(!r0) return fail(); r=r0; r.keepF=true; if(mod) r.comp=clone(r.comp); }
  const c=r.comp, ex0=r.exact;   // ex0: 바탕 모양의 exact (ST 는 위치를 정한 뒤 다시 정함)
  // 수식어를 먼저 붙인다 — ST 의 폭은 수식어 붙인 모양으로 (SPEC_casa_v3 §5; 바탕 폭으로 F 를 정하면 ST(1.7)LA(1,1,200) 에서 4.5 % 작음)
  if(mod){ if(!S2_MODBASE[c.shape]) return fail(`'${raw}': ${mod.n} applies to GL, SGL, LA, TLA (here: ${c.shape})`);
    c.mod=mod.n; if(mod.n==="ST"){ c.stMu=mod.a[0]; c.stN=mod.a[1]!=null?mod.a[1]:0; c.stG=0; if(!o.innoUnits) c.stCasaK=casaStK(c); } else c.sbArg=mod.a[0]!=null?mod.a[0]:0;
    r.exact=false; r.note=(mod.n==="ST"?"ST(μ,γ): S(μ)F = N·F·[∫_x^∞F + A/μ] ⊗ G(γ), Cookbook p. 66 — shape measured (ST(1.7)LA(1,1,200): S after ⊗G, step on the high-BE side); γ → FWHM γ·F/144.87 assumed; CasaXPS draws it with the step lagged by (base FWHM)/80 and positions it on the largest node of a (base FWHM)/40 grid (SPEC_casa_v4 §4); ":"SB(0)F = ∫_x^∞F (Shirley step), Cookbook p. 68 — area = area inside the fit window; ")+r.note; }
  // 폭: CasaXPS fwhm 숫자 → INNO F (또는 DS wL) — casaWidthAsym 을 거꾸로 (척도 불변: 단위 폭에서의 CasaXPS fwhm 으로 나눔). o.innoUnits: fw 는 실제 FWHM (화면에서 모양만 바꿀 때)
  if(!r.keepF||(mod&&mod.n==="ST")){ const sk=c.shape==="DS"?"wL":"F", u0=clone(c); u0[sk]=1;
    const cw=o.innoUnits?null:casaWidthAsym(u0), u=cw&&cw.fwhm>0?cw.fwhm:(trueFWHM(u0.shape,shapeParams(u0,1)).fwhm||unitFWHM(u0));
    c[sk]=u>0?fw/u:fw; if(cw&&!cw.exact&&!mod&&!r.keepF) r.note+="; "+cw.note+" — refit the width"; }
  // DS 위치: CasaXPS 위치 = INNO 원점 + wL·|cot(π/(2−α))| (높은 BE 쪽; DS(0.1,100) 판독 1 개, SPEC_casa_v3 §4)
  if(c.shape==="DS"&&c.alpha>0&&!o.innoUnits){ c.be=be+c.wL/Math.tan(PI/(2-c.alpha)); r.posShift=be-c.be; r.note+="; position: INNO origin = CasaXPS position − wL·|cot(π/(2−α))| (measured α 0, 0.1, 0.2; ≤ 1.6e-4 eV)"+(c.dsl>0?"; same shift assumed for the Shirley-cut variant":""); }
  else if(mod&&mod.n==="ST"&&!o.innoUnits){ const cw=casaWidthAsym(c); if(cw.node!=null){ c.be=be-cw.node; r.posShift=cw.node; r.exact=!!ex0&&cw.exact; r.note+="; position: INNO origin = CasaXPS position − "+cw.node.toFixed(5)+" eV (largest node of CasaXPS's ST table, step (base FWHM)/40)"; } else r.posShift=null; }
  else r.posShift=0;
  // DS area: CasaXPS Area (amplitude × ∫ over ±600·wL) → INNO area inside the fit window [o.emin, o.emax] (SPEC_casa_v4 §2)
  if(c.shape==="DS"&&!o.innoUnits){ const k=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;
    if(k){ c.area=area/k; r.areaFactor=k; r.exact=Number.isInteger(+c.casaN)&&!(c.dsl>0); r.note+="; area: INNO window area = CasaXPS area / "+k.toFixed(5)+" (CasaXPS DS area = amplitude × ∫ over ±600·wL)"; }
    else r.note+="; area kept as given — pass the fit window {emin, emax} to convert the CasaXPS DS area (CasaXPS area = amplitude × ∫ over ±600·wL)"; }
  delete r.keepF; return r; }
function compToCasa(c,o){
  o=o||{};
  const md=compMod(c), sp=shapeParams(c,1), r=(v,d)=>+(+v).toFixed(d);
  const nOf=(w,F)=>c.casaN!=null?+c.casaN:(w>1e-9&&F>0?w/(CASA_LA_RHO*F):0);
  const nTxt=n=>{ const ni=Math.round(n); return {t:String(ni),ok:Math.abs(n-ni)<0.05}; };
  let base=null;
  if(c.shape==="GLN"||c.shape==="SGLN"){ const n=nTxt(nOf(sp.wG,sp.F)), t=compTail(c), cw=t?casaWidthAsym(c):null;
    base={notation:`${c.shape==="GLN"?"GL":"SGL"}(${r(100*c.m,2)}${t&&n.t==="0"?"":","+n.t})${t==="T"?`T(${r(c.tk,4)})`:t==="PHI"?`PHI(${r(c.tj,4)},${r(c.tk,4)})`:""}`,fwhm:t?cw.fwhm:s2FWHM(c.shape,sp).fwhm,asym:t?cw.asym:0,fwhmExact:true,exact:n.ok,
      note:(c.shape==="SGLN"?"SGL(m,n): m = height fraction × 100":"GL(m,n): GL product")+"; Gaussian FWHM = n·F/144.87 (measured: LA(1,1,n), GL(30,100), SGL(30,100))"+(t?"; tail per Cookbook p. 67 (shape measured, GL(30)T(2)); "+cw.note:"; fwhm = true FWHM")+(n.ok?"":" — n not an integer")}; }
  else if(c.shape==="TLA"){ const n=nTxt(nOf(sp.wG,sp.F)), cw=casaWidthAsym(c); base={notation:`TLA(${r(c.alpha,3)},${r(c.mu,3)},${n.t})`,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,exact:n.ok&&cw.exact,note:"TLA(α,μ,n), Cookbook p. 66 — notation, tail side and n rule measured (TLA(1.2,3,0), TLA(1.2,3,100)); "+cw.note}; }
  else if(c.shape==="DS"){ const n=nTxt(nOf(sp.wG,sp.wL)), cw=casaWidthAsym(c), sh=c.alpha>0?sp.wL*Math.abs(1/Math.tan(PI/(2-c.alpha))):0, ak=o.emin!=null&&o.emax!=null?dsCasaAreaFactor(c,+o.emin,+o.emax):null;
    base={notation:`DS(${r(c.alpha,3)},${n.t}${c.dsl>0?","+r(c.dsl,3):""})`,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:sh,areaFactor:ak,area:ak?c.area*ak:null,exact:!!ak&&n.ok,note:"DS(α,n"+(c.dsl>0?",l":"")+"), Cookbook p. 66: Lorentzian FWHM = wL; shape, tail side and n rule measured (α 0, 0.1, 0.2; n 0, 100); "+cw.note+"; CasaXPS position = INNO origin + wL·|cot(π/(2−α))| (measured ≤ 1.6e-4 eV)"+(c.dsl>0?" — assumed for the cut variant":"")+"; CasaXPS area = amplitude × ∫ over ±600·wL"+(ak?" = INNO window area × "+ak.toFixed(5):" — pass the fit window {emin, emax} to convert (INNO area = area inside the fit window)")}; }
  else if(c.shape==="TP"){ const n=nTxt(nOf(sp.wL,sp.F)); base={notation:`TP(${r(c.ta,3)},${r(c.tb,3)},${n.t})`,fwhm:c.F,exact:false,note:"TP(α,β,n), Cookbook p. 68: trapezium of width F (α at the high-BE end) ⊗ Lorentzian FWHM n·F/144.87 (assumed)"}; }
  else base=compToCasa0(md?Object.assign({},c,{mod:null}):c);
  const dflt=o=>Object.assign({asym:o.fwhm!=null?0:null,posShift:0,fwhmExact:!!o.exact},o);   // 기본: 대칭(asym 0), 위치 같음
  if(!md||!base.notation) return dflt(md&&!base.notation?{notation:null,fwhm:null,exact:false,note:"no CasaXPS notation for this base shape with "+md}:base);
  if(md==="ST"){ const g=nTxt(c.stN!=null?+c.stN:(c.stG>1e-9?c.stG/(CASA_LA_RHO*(c.shape==="V"?sp.wL:sp.F)):0)), cw=casaWidthAsym(c);
    return {notation:`ST(${r(c.stMu,4)}${g.t!=="0"?","+g.t:""})${base.notation}`,fwhm:cw.fwhm,asym:cw.asym,fwhmExact:cw.exact,posShift:cw.node,exact:cw.exact&&base.exact,note:"ST(μ,γ) modifier, Cookbook p. 66 (shape measured, ST(1.7)LA(1,1,200)); γ → Gaussian FWHM γ·F/144.87 assumed; "+cw.note+"; "+base.note}; }
  return {notation:`SB(${r(c.sbArg||0,3)})${base.notation}`,fwhm:null,asym:null,fwhmExact:false,posShift:0,exact:false,note:"SB(0) Shirley profile ∫_x^∞F, Cookbook p. 68; area = area inside the fit window; "+base.note}; }
/* CasaXPS 배경 이름 + Cross Section 칸("B, C, D, T0") → INNO bg. U 2/U 3 Tougaard 에서 B = 0 (U 3 기본값 "0, 551, 436, 3") 이면 CasaXPS 는 Linear 를 그린다(실측) */
function casaBg(name,cross,aw){
  const nm=String(name||"").trim(), w=aw!=null?+aw:1;
  if(/^linear$/i.test(nm)) return {bg:{type:"linear",conv:"casa",aw:w},exact:true,note:"CasaXPS Linear (SPEC §2b)"};
  if(/^shirley$/i.test(nm)) return {bg:{type:"shirley",conv:"casa",aw:w},exact:true,note:"CasaXPS Shirley (SPEC §2b)"};
  const t=nm.match(/^U\s*([234])\s*Tougaard$/i); if(!t) return null;   // U 4 (CasaXPS Cookbook 2026 p. 70: U(x; B,C,D,T0) = Bx/((C − x²)² + Dx²), x > T0) = 같은 핵에 T0 — INNO u3 핵 + T0
  const u3=t[1]!=="2"; let v=String(cross==null?"":cross).split(/[,;\s]+/).filter(x=>x!=="").map(Number);
  if(!v.length) v=u3?[0,551,436,3]:[2866,1643,1,0];   // CasaXPS 2.3.27 이 채우는 값 (U 3 실측; U 2 는 앞선 판독의 칸)
  const [B,C,D,T0]=[v[0],v[1],v[2]!=null?v[2]:0,v[3]!=null?v[3]:0], tg={name:"U "+t[1]+" Tougaard",B,C,D,T0};
  if(!(B>0)) return {bg:{type:"linear",conv:"casa",aw:w,casaTougaard:tg},exact:false,note:"B = 0: treated as Linear — unverified (CasaXPS region table did not recompute Tougaard from typed cells in our check)"};
  return {bg:{type:"tougaard",u3,C:isFinite(C)?C:(u3?551:1643),D:u3?(isFinite(D)?D:0):0,T0:T0>0?T0:0,casaTougaard:tg},exact:false,note:"explicit B: INNO scales B so the background meets the high-BE end; CasaXPS uses the given B — numbers differ (not pinned)"+(!u3&&D?"; U 2 ignores D here":"")};
}
/* CasaXPS VAMAS 주석 모형(XPS_IO.casaCommentParse) → INNO 피팅 상태 {lo, hi, bg, comps, cons, weight:"casa", notes} (v1.8).
   에너지: 측정 BE = hν − KE − shift (CasaXPS 주석의 KE 는 Calib 이동이 들어간 눈금; shift 는 따로 돌려줌 — 시료 대전 보정으로 쓰면 됨). 성분 모양·폭은 casaToComp, 배경은 casaBg.
   연결: Area link j, k → 넓이 = k × 성분 j 넓이; MFWHM link j, k → 같은 모양이면 F = k × F_j (CasaXPS fwhm 비 = F 비), 모양이 다르면 연결하지 않고 알림;
   Position link 는 오프셋 표기를 아직 확인하지 못해(실측 파일은 모두 −1) 연결하지 않고 알림. 한계(lo, hi)는 넓이·위치는 그대로, 폭은 F 로 환산. */
function casaModelToINNO(cm,hv,o){
  if(!cm||!cm.regions||!cm.regions.length) return null; hv=hv>0?+hv:1486.6; o=o||{};
  const R=cm.regions[0], notes=[], sh=+cm.shift||0; const be=ke=>hv-ke-sh;   // CasaXPS 는 영역·성분 에너지를 보정(Calib)된 눈금으로 저장한다(실측: 영역 끝 채널 = hν − KE − shift)
  let lo=Math.min(be(R.keStart),be(R.keEnd)), hi=Math.max(be(R.keStart),be(R.keEnd));
  // 실측(SPEC_real_validation §B): 배경 끝점 = 영역 안 첫·끝 채널 중심 평균, 주석의 Av.Width 2 → ±1 채널. INNO casa 규약(§2b)은 창의 두 번째·끝에서 두 번째 채널을 쓰므로
  // 창을 한 채널씩 넓힌다(step 을 알 때). 맞추는 채널(casaVar: 창 끝 두 채널 제외)도 CasaXPS 와 같아진다.
  const st=Math.abs(o.step)||0; if(st>0){ lo-=st+1e-6; hi+=st+1e-6; } else notes.push("channel step unknown — region not widened by one channel (background end points differ from CasaXPS by one channel)");
  const aw=R.aw!=null&&isFinite(R.aw)?R.aw/2:1;
  const g=casaBg(R.bg,R.cross&&R.cross.length?R.cross.join(","):null,aw); let bg;
  if(g){ bg=g.bg; if(!g.exact) notes.push("background: "+g.note); } else { bg={type:"shirley",conv:"casa",aw}; notes.push(`background '${R.bg}' has no INNO equivalent — Shirley (CasaXPS convention) used`); }
  if((R.startOffset||0)!==0||(R.endOffset||0)!==0) notes.push("region Start/End Offset ≠ 0 is not applied (INNO has no offset term)");
  const comps=[], cons={}, ids=cm.comps.map((c,i)=>"casa"+i), wf=[], psh=[], af=[];   // psh: CasaXPS 위치 − INNO 원점 (DS: wL·|cot|, ST: 가장 큰 마디); af: CasaXPS 넓이 / INNO 넓이 (DS, SPEC_casa_v4 §2; 그 밖 1)
  cm.comps.forEach((c,i)=>{ const r=casaToComp(c.shape,{be:c.pos?be(c.pos.v):0,area:c.area?c.area.v:1,fwhm:c.fwhm?c.fwhm.v:1.2,id:ids[i],emin:lo,emax:hi});   // DS 넓이: INNO 창 [lo, hi] 안 넓이로 환산
    if(!r){ notes.push(`${c.name}: line shape '${c.shape}' not supported — skipped`); comps.push(null); wf.push(null); af.push(1); psh.push(0); return; }
    r.comp.name=c.name; comps.push(r.comp); if(!r.exact) notes.push(`${c.name}: ${r.note}`); psh.push(r.posShift); af.push(r.areaFactor>0?r.areaFactor:1);
    const W=r.comp.shape==="DS"?r.comp.wL:r.comp.F; const k=c.fwhm&&c.fwhm.v>0&&W>0?W/c.fwhm.v:1; wf.push(k); });
  cm.comps.forEach((c,i)=>{ const C=comps[i]; if(!C) return; const q=cons[C.id]={};
    // 넓이: CasaXPS 넓이 = af × INNO 넓이 → 한계는 ÷ af, 연결 A_i = k·A_j (CasaXPS) → INNO k·af_j/af_i (DS 와 다른 성분 사이에서 인자가 실림)
    if(c.area){ if(c.area.link>=0&&comps[c.area.link]) q.area={mode:"link",ref:ids[c.area.link],k:c.area.k*af[c.area.link]/af[i],b:0}; else q.area={mode:"free",lo:Math.max(0,c.area.lo)/af[i],hi:c.area.hi/af[i]}; }
    if(c.fwhm&&(C.F!=null||C.shape==="DS")){ const j=c.fwhm.link, wk=C.shape==="DS"?"wL":"F";
      if(j>=0&&comps[j]){ const P=comps[j]; if(P.shape===C.shape&&Math.abs((C.m||0)-(P.m||0))<1e-12&&C.casaN===P.casaN) q[wk]={mode:"link",ref:ids[j],k:c.fwhm.k,b:0}; else notes.push(`${c.name}: FWHM link to a different line shape not reproduced (fit independently)`); }
      else q[wk]={mode:"free",lo:c.fwhm.lo*wf[i],hi:c.fwhm.hi*wf[i]}; }
    // 위치: CasaXPS 위치 한계·오프셋은 CasaXPS 위치에 걸린다 → INNO 원점 = CasaXPS 위치 − posShift (DS: wL·|cot(π/(2−α))|; ST: 표의 가장 큰 마디, SPEC_casa_v4 §4; 못 정하면 0 으로 두고 알림)
    if(c.pos){ const s0=psh[i]||0;
      if(c.pos.link>=0){ const sj=comps[c.pos.link]?psh[c.pos.link]:0; notes.push(`${c.name}: position link not reproduced (CasaXPS offset notation not verified)`+((psh[i]==null||sj==null||Math.abs((sj||0)-s0)>1e-9)?`; note: CasaXPS offsets apply to CasaXPS positions, and this component's CasaXPS-position − INNO-origin shift (${psh[i]==null?"not pinned":s0.toFixed(4)+" eV"}) differs from its partner's (${sj==null?"not pinned":(sj||0).toFixed(4)+" eV"}) — an INNO offset = CasaXPS offset − (shift difference)`:"")); }
      q.be={mode:"free",lo:be(c.pos.hi)-s0,hi:be(c.pos.lo)-s0}; if(psh[i]==null) notes.push(`${c.name}: ST position offset not pinned — position bounds taken as CasaXPS positions`); }
    // Cookbook 2026 모양(SPEC_shapes2): CasaXPS 처럼 표기 안의 모양 매개변수는 고정 (LA·GL·SGL 은 이전 동작 그대로)
    if(compMod(C)||["TLA","GLN","SGLN","TP","DS"].includes(C.shape)) paramNames(C).forEach(k=>{ if(["be","area","F"].includes(k)||(C.shape==="DS"&&k==="wL")||q[k]||derivedParam(C,k)) return; q[k]={mode:"fixed"}; }); });
  return {lo,hi,bg,comps:comps.filter(Boolean),cons,weight:"casa",shift:cm.shift||0,notes,region:R.name};
}
const XPS_CORE={_hook:FIT_HOOK,_work:WORK,LM_STALL,areaSumSigma,finalize,fitJob,finalJob,mcJob,SHAPES,PDEF,BG_TYPES,faddeeva,erfcx,lgamma,gauss,lorentz,voigt,glpNorm,laNorm,shapeUnit,trueFWHM,compEval,shapeParams,
  cumtrapz,endpoints,smoothQ,smoothHW,bgLinear,bgShirley,bgTougaard,bgStatic,
  paramNames,getP,setP,defBounds,buildMap,applyLinks,getVec,setVec,toU,fromU,makeEvaluator,
  nnls,rng,fitOnce,fitMulti,startPoints,applyStart,lmFit,covariance,fitStats,compSummary,monteCarlo,mcSummarize,noiseEstimate,
  tfValue,tpp2m,imfpValue,quantify,quantMC,endpointsCasa,bgShirleyCasa,casaEnds,bgShirleyAvantage,bgSmartAvantage,avantagePad,bgShirleyAvantageV18,bgSmartAvantageV18,avantageV19,avantageHalf,avantagePQ,avantageSmooth,aglRaw,aglShape,avantageVoigt,voigtFromAvantage,casaRegionStats,pchipSlopes,interpAt,interpIntegral,makeTemplate,perturbTpl,tplShape,tplCache,lcfDiag,symEig,splineKnots,naturalSplineM,splineAt,bgSpline,nmFit,propagateModel,batchJob,shiftModelBE,constVar,parseFormula,stoichCompare,augerParam,peakMaxPos,findPeaks,identify,acceptDefault,elementLines,EXCLUDE_DEFAULT,doubletInfo,multipletInfo,dEConstraint,suggestModel,refineByResidual,autoJob,sensVariants,sensJob,sensSummary,sensBg,SENS_BG_LABEL,removeComp,anchorOf,DOUBLET_RESOLVE,sliceRegion,makeProblem,newComp,convertShape,clone,stripProb,CASA_LA_RHO,casaVar,sglToGLS,glsToSGL,casaWidthFactor,casaWidthAsym,refineHalfMax,CASA_GRID,CASA_DS_T,dsCasaAreaFactor,casaStK,compToCasa,casaToComp,casaBg,casaModelToINNO,
  s2Base,s2Eval,s2FWHM,tlaI1,tlaJ,tpEval,quadA,quadInf,lag6,gridInt6,compMod,compTail,shapePList,casaNKey,derivedParam,defaultCons,casaTokens,unitFWHM,compToCasa0,casaToComp0,CASA_SHAPES_HELP,S2_MODBASE};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_CORE; else root.XPS_CORE=XPS_CORE;
})(typeof self!=="undefined"?self:this);
