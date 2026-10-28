/* ============================================================
   INNO Analytics — XPS spectral pre-processing engine (no DOM: browser / Web Worker / node)
   Definitions: modules/xps/SPEC_proc.md. This file implements that document; where they differ,
   the document is fixed first and the code follows.
     §1 noise model            §2 resampling (rebin / PCHIP / Lanczos / linear)
     §3 smoothing & derivatives (Savitzky–Golay = scipy.signal.savgol_filter, Gaussian = gaussian_filter1d)
     §4 PSF kernels            §5 deconvolution (Jansson, Van Cittert, Richardson–Lucy) + discrepancy stop
     §6 Fermi-edge fit         §7 PSF from a reference peak (Voigt fit)
   Published methods only; no third-party code or data.
   ============================================================ */
(function(root){
"use strict";
const LN2=Math.LN2, PI=Math.PI, SQ2=Math.SQRT2, SQPI=Math.sqrt(PI), SQ2PI=Math.sqrt(2*PI);
const FW2SIG=1/(2*Math.sqrt(2*LN2));          // Gaussian FWHM -> sigma
const KB_EV=8.617333262e-5;                    // Boltzmann constant, eV/K (CODATA 2018, exact)
const Z975=1.959963984540054;                  // two-sided 95 % normal quantile
const F64=a=>Float64Array.from(a);

/* ============================================================
   Special functions
   ============================================================ */
function erfcCF(x){ // x >= 2.5 : erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x+ (1/2)/(x+ 1/(x+ (3/2)/(x+ ...))))
  let f=x; for(let n=60;n>=1;n--) f=x+(n/2)/f;
  return Math.exp(-x*x)/(SQPI*f);
}
function erf(x){
  const ax=Math.abs(x); let v;
  if(ax<2.5){ const x2=ax*ax; let t=ax, s=ax; for(let n=1;n<200;n++){ t*=-x2/n; const d=t/(2*n+1); s+=d; if(Math.abs(d)<1e-17*Math.abs(s)) break; } v=2/SQPI*s; }
  else v=1-erfcCF(ax);
  return x<0?-v:v;
}
function erfc(x){   // x >= 0.5: erfc(x) = exp(-x²)·Re w(ix) (Weideman), keeps full relative accuracy in the tail
  if(x>=0.5) return Math.exp(-x*x)*faddeeva(0,x)[0]; if(x<=-0.5) return 2-Math.exp(-x*x)*faddeeva(0,-x)[0]; return 1-erf(x); }
const Phi=z=>0.5*erfc(-z/SQ2);                 // standard normal CDF
const phi=z=>Math.exp(-0.5*z*z)/SQ2PI;
/* Faddeeva w(z), Im z >= 0 — Weideman, SIAM J. Numer. Anal. 31, 1497 (1994), rational series with N = 40 */
const WD=(function(){
  const N=40, M=2*N, M2=2*M, L=Math.sqrt(N/SQ2);
  const f=new Float64Array(M2);
  for(let k=-M+1;k<=M-1;k++){ const th=k*PI/M, t=L*Math.tan(th/2); f[k+M]=Math.exp(-t*t)*(L*L+t*t); }
  f[0]=0;
  const g=new Float64Array(M2); for(let i=0;i<M2;i++) g[i]=f[(i+M2/2)%M2];
  const a=new Float64Array(N+1);
  for(let j=1;j<=N;j++){ let re=0; for(let m=0;m<M2;m++) re+=g[m]*Math.cos(2*PI*j*m/M2); a[j]=re/M2; }
  return {N,L,a};
})();
function faddeeva(x,y){
  const L=WD.L, a=WD.a, N=WD.N;
  const dr=L+y, di=-x, nr=L-y, ni=x, den=dr*dr+di*di;
  const Zr=(nr*dr+ni*di)/den, Zi=(ni*dr-nr*di)/den;
  let pr=a[N], pi=0;
  for(let j=N-1;j>=1;j--){ const tr=pr*Zr-pi*Zi+a[j], ti=pr*Zi+pi*Zr; pr=tr; pi=ti; }
  const ir=dr/den, ii=-di/den, i2r=ir*ir-ii*ii, i2i=2*ir*ii;
  return [2*(pr*i2r-pi*i2i)+ir/SQPI, 2*(pr*i2i+pi*i2r)+ii/SQPI];
}
/* area-1 profiles; FWHM arguments */
function gauss(x,F){ const s=F*FW2SIG; return Math.exp(-0.5*x*x/(s*s))/(s*SQ2PI); }
function lorentz(x,F){ return (2/(PI*F))/(1+4*x*x/(F*F)); }
function voigt(x,wG,wL){
  if(!(wL>0)) return gauss(x,wG); if(!(wG>0)) return lorentz(x,wL);
  const s=wG*FW2SIG, k=1/(s*SQ2); return faddeeva(x*k,(wL/2)*k)[0]/(s*SQ2PI);
}
const LG=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
function lgamma(x){
  if(x<0.5) return Math.log(PI/Math.abs(Math.sin(PI*x)))-lgamma(1-x);
  x-=1; let a=LG[0]; const t=x+7.5; for(let i=1;i<9;i++) a+=LG[i]/(x+i);
  return 0.5*Math.log(2*PI)+(x+0.5)*Math.log(t)-t+Math.log(a);
}
/* Gauss–Legendre 8-point nodes/weights on [-1,1] */
const GL8X=[-0.9602898564975363,-0.7966664774136267,-0.5255324099163290,-0.1834346424956498,0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363];
const GL8W=[0.1012285362903763,0.2223810344533745,0.3137066458778873,0.3626837833783620,0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763];

/* ============================================================
   RNG (tests / Monte-Carlo): xoshiro128** + Box–Muller + Poisson (Knuth / PTRS, Hörmann 1993)
   ============================================================ */
function rng(seed){
  let s0=(seed>>>0)||1, s1=0x9E3779B9, s2=0x243F6A88, s3=0xB7E15162;
  const rotl=(x,k)=>(x<<k)|(x>>>(32-k));
  const next=()=>{ const r=Math.imul(rotl(Math.imul(s1,5),7),9); const t=s1<<9; s2^=s0; s3^=s1; s1^=s2; s0^=s3; s2^=t; s3=rotl(s3,11); return (r>>>0)/4294967296; };
  for(let i=0;i<16;i++) next();
  let spare=null;
  const u=()=>{ let v; do{ v=next(); }while(v===0); return v; };
  const normal=()=>{ if(spare!==null){ const v=spare; spare=null; return v; } const a=u(), b=u(); const r=Math.sqrt(-2*Math.log(a)); spare=r*Math.sin(2*PI*b); return r*Math.cos(2*PI*b); };
  const poisson=lam=>{
    if(!(lam>0)) return 0;
    if(lam<30){ const L=Math.exp(-lam); let k=0, p=1; do{ k++; p*=next(); }while(p>L); return k-1; }
    const slam=Math.sqrt(lam), loglam=Math.log(lam), b=0.931+2.53*slam, a=-0.059+0.02483*b, invalpha=1.1239+1.1328/(b-3.4), vr=0.9277-3.6224/(b-2);
    for(;;){ const U=next()-0.5, V=next(), us=0.5-Math.abs(U), k=Math.floor((2*a/us+b)*U+lam+0.43);
      if(us>=0.07&&V<=vr) return k;
      if(k<0||(us<0.013&&V>us)) continue;
      if(Math.log(V)+Math.log(invalpha)-Math.log(a/(us*us)+b) <= -lam+k*loglam-lgamma(k+1)) return k; }
  };
  return {next,normal,poisson};
}

/* ============================================================
   Small utilities
   ============================================================ */
function gridInfo(E){
  const n=E.length; if(n<2) return {uniform:false,step:NaN,ascending:false};
  const st=(E[n-1]-E[0])/(n-1); let dev=0, asc=true;
  for(let i=1;i<n;i++){ const d=E[i]-E[i-1]; if(!(d>0)) asc=false; dev=Math.max(dev,Math.abs(d-st)); }
  return {uniform:asc&&dev<=1e-6*Math.abs(st)+1e-12, step:st, ascending:asc, maxDev:dev};
}
function median(a){ const b=Float64Array.from(a).sort(); const n=b.length; return n%2?b[(n-1)/2]:0.5*(b[n/2-1]+b[n/2]); }
function trapz(x,y){ let s=0; for(let i=1;i<y.length;i++) s+=0.5*(y[i]+y[i-1])*(x[i]-x[i-1]); return s; }
/* FWHM of a single-maximum curve by linear interpolation at half maximum (above `base`) */
function fwhm(x,y,base){
  base=base||0; let im=0; for(let i=1;i<y.length;i++) if(y[i]>y[im]) im=i;
  const h=(y[im]-base)/2+base; let l=NaN, r=NaN;
  for(let i=im;i>0;i--) if(y[i-1]<h){ l=x[i-1]+(h-y[i-1])/(y[i]-y[i-1])*(x[i]-x[i-1]); break; }
  for(let i=im;i<y.length-1;i++) if(y[i+1]<h){ r=x[i]+(y[i]-h)/(y[i]-y[i+1])*(x[i+1]-x[i]); break; }
  return {fwhm:r-l,left:l,right:r,imax:im,xmax:x[im],ymax:y[im]};
}
/* boundary index mapping (scipy.ndimage conventions) — returns -1 for 'constant' outside */
function mapIndex(k,n,mode){
  if(k>=0&&k<n) return k;
  if(mode==="constant") return -1;
  if(mode==="nearest") return k<0?0:n-1;
  if(mode==="wrap") return ((k%n)+n)%n;
  if(n===1) return 0;
  if(mode==="mirror"){ const p=2*(n-1); let m=((k%p)+p)%p; return m<n?m:p-m; }       // d c b | a b c d | c b a
  if(mode==="reflect"){ const p=2*n; let m=((k%p)+p)%p; return m<n?m:p-1-m; }        // d c b a | a b c d | d c b a
  throw new Error("unknown boundary mode "+mode);
}

/* ============================================================
   §1 Noise model
   var_i = kappa * max(level_i, floor)  ("poisson-like"), or var_i = s2 ("constant").
   kappa = median_i[(D2 y_i)^2 / (6 level_i)] / 0.454936...  (median of chi^2_1)
   level_i = 5-point moving mean (edges: available points), D2 y_i = y_{i-1} - 2 y_i + y_{i+1}
   s2 = (1.482602 * MAD(D2 y))^2 / 6
   ============================================================ */
const MED_CHI2_1=0.45493642311957283;
function localLevel(y){ const n=y.length, o=new Float64Array(n); for(let i=0;i<n;i++){ let s=0,c=0; for(let j=Math.max(0,i-2);j<=Math.min(n-1,i+2);j++){ s+=y[j]; c++; } o[i]=s/c; } return o; }
function noiseModel(y,opt){
  opt=opt||{}; const n=y.length; if(n<5) throw new Error("noiseModel: need >= 5 points");
  const lv=localLevel(y); let ymax=-Infinity; for(let i=0;i<n;i++) ymax=Math.max(ymax,Math.abs(y[i]));
  const floor=opt.floor!=null?opt.floor:1e-3*ymax;
  const d2=new Float64Array(n-2); for(let i=1;i<n-1;i++) d2[i-1]=y[i-1]-2*y[i]+y[i+1];
  const q=new Float64Array(n-2); for(let i=1;i<n-1;i++) q[i-1]=d2[i-1]*d2[i-1]/(6*Math.max(lv[i],floor));
  const kappa=median(q)/MED_CHI2_1;
  const md=median(d2); const ad=d2.map(v=>Math.abs(v-md)); const s=1.482602218505602*median(ad); const s2=s*s/6;
  let minY=Infinity; for(let i=0;i<n;i++) minY=Math.min(minY,y[i]);
  const kind=opt.kind||(minY>0?"poisson-like":"constant");
  const v=new Float64Array(n);
  for(let i=0;i<n;i++) v[i]=kind==="constant"?s2:kappa*Math.max(lv[i],floor);
  return {kind,kappa,s2,var:v,floor};
}
function resolveVar(y,noise){
  const n=y.length;
  if(noise&&noise.length===n) return {var:F64(noise),kind:"given"};
  if(typeof noise==="number") return {var:new Float64Array(n).fill(noise),kind:"constant-given"};
  if(noise==="poisson"){ const v=new Float64Array(n); for(let i=0;i<n;i++) v[i]=Math.max(y[i],1); return {var:v,kind:"poisson"}; }
  const m=noiseModel(y); return {var:m.var,kind:"auto-"+m.kind,kappa:m.kappa,s2:m.s2};
}

/* ============================================================
   §2 Resampling. Each method returns sparse weight rows W (y_out = W y), so the variance of independent
   inputs propagates exactly as var_out_j = sum_i W_ji^2 var_i (PCHIP: first-order, see SPEC_proc §2.3).
   ============================================================ */
function edgesOf(x){ const n=x.length, e=new Float64Array(n+1);
  if(n===1){ e[0]=x[0]-0.5; e[1]=x[0]+0.5; return e; }
  e[0]=x[0]-(x[1]-x[0])/2; for(let i=1;i<n;i++) e[i]=(x[i-1]+x[i])/2; e[n]=x[n-1]+(x[n-1]-x[n-2])/2; return e; }
/* §2.1 flux-conserving rebin (Carnall 2017): y_j = sum_i y_i |B_i ∩ B'_j| / sum_i |B_i ∩ B'_j| */
function rebinRows(E,En){
  const ei=edgesOf(E), eo=edgesOf(En), n=E.length, rows=[], cover=new Float64Array(En.length);
  let i0=0;
  for(let j=0;j<En.length;j++){
    const a=eo[j], b=eo[j+1]; const idx=[], w=[]; let tot=0;
    while(i0<n&&ei[i0+1]<=a) i0++;
    for(let i=i0;i<n&&ei[i]<b;i++){ const ov=Math.min(b,ei[i+1])-Math.max(a,ei[i]); if(ov>0){ idx.push(i); w.push(ov); tot+=ov; } }
    for(let k=0;k<w.length;k++) w[k]/=tot;
    rows.push(tot>0?{idx,w}:null); cover[j]=tot/(b-a);
  }
  return {rows,cover};
}
function linearRows(E,En){
  const n=E.length, rows=[];
  for(let j=0;j<En.length;j++){ const x=En[j];
    if(x<E[0]||x>E[n-1]){ rows.push(null); continue; }
    let k=bsearch(E,x); if(k>=n-1) k=n-2; const t=(x-E[k])/(E[k+1]-E[k]);
    rows.push({idx:[k,k+1],w:[1-t,t]}); }
  return {rows};
}
function bsearch(E,x){ let lo=0, hi=E.length-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(E[m]<=x) lo=m; else hi=m; } return lo; }
/* §2.3 PCHIP exactly as scipy.interpolate.PchipInterpolator (Fritsch–Carlson 1980 / Fritsch–Butland 1984)
   with analytic linearisation: the weight rows reproduce the value exactly (Euler: d is 1-homogeneous in y). */
function pchipDeriv(E,y){
  const n=E.length, h=new Float64Array(n-1), m=new Float64Array(n-1);
  for(let k=0;k<n-1;k++){ h[k]=E[k+1]-E[k]; m[k]=(y[k+1]-y[k])/h[k]; }
  // d_k as sparse linear form in the slopes m: dm[k] = [[slopeIndex, coef], ...]
  const d=new Float64Array(n), dm=new Array(n);
  if(n===2){ d[0]=d[1]=m[0]; dm[0]=[[0,1]]; dm[1]=[[0,1]]; return {d,dm,h,m}; }
  for(let k=1;k<n-1;k++){
    const m0=m[k-1], m1=m[k];
    if(Math.sign(m0)!==Math.sign(m1)||m0===0||m1===0){ d[k]=0; dm[k]=[]; continue; }
    const w1=2*h[k]+h[k-1], w2=h[k]+2*h[k-1], D=w1/m0+w2/m1;
    d[k]=(w1+w2)/D; dm[k]=[[k-1,(w1+w2)*w1/(m0*m0)/(D*D)],[k,(w1+w2)*w2/(m1*m1)/(D*D)]];
  }
  const edge=(h0,h1,m0,m1,i0,i1)=>{
    let dd=((2*h0+h1)*m0-h0*m1)/(h0+h1); let f=[[i0,(2*h0+h1)/(h0+h1)],[i1,-h0/(h0+h1)]];
    if(Math.sign(dd)!==Math.sign(m0)){ dd=0; f=[]; }
    else if(Math.sign(m0)!==Math.sign(m1)&&Math.abs(dd)>3*Math.abs(m0)){ dd=3*m0; f=[[i0,3]]; }
    return [dd,f];
  };
  [d[0],dm[0]]=edge(h[0],h[1],m[0],m[1],0,1);
  [d[n-1],dm[n-1]]=edge(h[n-2],h[n-3],m[n-2],m[n-3],n-2,n-3);
  return {d,dm,h,m};
}
function pchipRows(E,y,En){
  const n=E.length, P=pchipDeriv(E,y), rows=[], val=new Float64Array(En.length);
  for(let j=0;j<En.length;j++){ const x=En[j];
    if(x<E[0]||x>E[n-1]){ rows.push(null); val[j]=NaN; continue; }
    let k=bsearch(E,x); if(k>=n-1) k=n-2;
    const h=P.h[k], t=(x-E[k])/h, t2=t*t, t3=t2*t;
    const h00=2*t3-3*t2+1, h10=t3-2*t2+t, h01=-2*t3+3*t2, h11=t3-t2;
    const acc=new Map(); const add=(i,w)=>acc.set(i,(acc.get(i)||0)+w);
    add(k,h00); add(k+1,h01);
    [[k,h10*h],[k+1,h11*h]].forEach(([q,c])=>P.dm[q].forEach(([s,a])=>{ add(s+1,c*a/P.h[s]); add(s,-c*a/P.h[s]); }));
    val[j]=h00*y[k]+h10*h*P.d[k]+h01*y[k+1]+h11*h*P.d[k+1];
    const idx=[...acc.keys()]; rows.push({idx,w:idx.map(i=>acc.get(i))});
  }
  return {rows,val};
}
/* §2.4 Lanczos-windowed sinc (Duchon 1979), uniform input only, weights renormalised to sum 1 */
function lanczosRows(E,En,a){
  a=a||3; const g=gridInfo(E); if(!g.uniform) throw new Error("lanczos: input grid must be uniform");
  const n=E.length, st=g.step, rows=[];
  const L=t=>{ if(t===0) return 1; if(Math.abs(t)>=a) return 0; const p=PI*t; return a*Math.sin(p)*Math.sin(p/a)/(p*p); };
  for(let j=0;j<En.length;j++){ const u=(En[j]-E[0])/st;
    if(u<-1e-9||u>n-1+1e-9){ rows.push(null); continue; }
    const k0=Math.floor(u); const idx=[], w=[]; let s=0;
    for(let i=k0-a+1;i<=k0+a;i++){ if(i<0||i>=n) continue; const v=L(u-i); if(v!==0){ idx.push(i); w.push(v); s+=v; } }
    rows.push({idx,w:w.map(v=>v/s)}); }
  return {rows};
}
function applyRows(rows,y,v){
  const m=rows.length, out=new Float64Array(m), ov=v?new Float64Array(m):null;
  for(let j=0;j<m;j++){ const r=rows[j]; if(!r||!r.idx.length){ out[j]=NaN; if(ov) ov[j]=NaN; continue; }
    let s=0, q=0; for(let k=0;k<r.idx.length;k++){ s+=r.w[k]*y[r.idx[k]]; if(v) q+=r.w[k]*r.w[k]*v[r.idx[k]]; } out[j]=s; if(ov) ov[j]=q; }
  return {y:out,var:ov};
}
function uniformGrid(E0,E1,step){ const n=Math.floor((E1-E0)/step+1e-9)+1; const g=new Float64Array(n); for(let i=0;i<n;i++) g[i]=E0+i*step; return g; }
/* Gaussian line sampled at step Δ: relative spectral amplitude at the Nyquist frequency */
function aliasLevel(fw,step){ const s=fw*FW2SIG; return Math.exp(-PI*PI*s*s/(2*step*step)); }
function suggestStep(fw,tol){ tol=tol||1e-3; const s=fw*FW2SIG; return PI*s/Math.sqrt(2*Math.log(1/tol)); }
function resample(E,y,opt){
  opt=opt||{}; E=F64(E); y=F64(y); const gi=gridInfo(E);
  if(!gi.ascending) throw new Error("resample: E must be strictly ascending");
  const En=opt.Enew?F64(opt.Enew):uniformGrid(E[0],E[E.length-1],opt.step||gi.step);
  const stepOut=En.length>1?(En[En.length-1]-En[0])/(En.length-1):gi.step;
  let method=opt.method||"auto";
  if(method==="auto") method=stepOut>=gi.step*(1-1e-9)?"rebin":"pchip";
  const v=opt.var?F64(opt.var):null; let R, extra={};
  if(method==="rebin"){ R=rebinRows(E,En); extra.cover=R.cover; }
  else if(method==="pchip") R=pchipRows(E,y,En);
  else if(method==="lanczos") R=lanczosRows(E,En,opt.a||3);
  else if(method==="linear") R=linearRows(E,En);
  else throw new Error("resample: unknown method "+method);
  const r=applyRows(R.rows,y,v);
  const warnings=[];
  if(method==="rebin"&&stepOut<gi.step*0.999) warnings.push("rebin onto a finer grid gives a staircase; use pchip or lanczos for up-sampling");
  if(extra.cover){ let pc=0; for(let j=0;j<En.length;j++) if(extra.cover[j]<0.999) pc++; if(pc) warnings.push(pc+" output bin(s) only partly covered by input (edges)"); }
  if(method!=="rebin"&&opt.fwhm&&stepOut>opt.fwhm/2) warnings.push("output step > FWHM/2: point interpolation aliases; use rebin");
  const finite=[]; for(let j=0;j<En.length;j++) if(isFinite(r.y[j])) finite.push(j);
  const areaIn=trapz(E,y); const Ef=F64(finite.map(j=>En[j])), yf=F64(finite.map(j=>r.y[j]));
  return {E:En,y:r.y,var:r.var,method,rows:opt.keepRows?R.rows:undefined,cover:extra.cover,
    areaIn,areaOut:trapz(Ef,yf),stepIn:gi.step,stepOut,warnings};
}

/* ============================================================
   §3 Smoothing / derivatives
   Savitzky–Golay: identical to scipy.signal.savgol_filter(x, window, order, deriv, delta, mode, cval)
   (odd window). Coefficients by least squares through a QR of the scaled Vandermonde matrix.
   ============================================================ */
function polyProjector(w,order,c){   // rows of B = R^{-1} Q^T for u_j = (j-c)/s, j = 0..w-1 ; returns {B, s}
  const s=Math.max(1,(w-1)/2), P=order+1;
  const Q=[]; for(let k=0;k<P;k++){ const col=new Float64Array(w); for(let j=0;j<w;j++) col[j]=Math.pow((j-c)/s,k); Q.push(col); }
  const R=[]; for(let k=0;k<P;k++) R.push(new Float64Array(P));
  for(let k=0;k<P;k++){ for(let pass=0;pass<2;pass++) for(let i=0;i<k;i++){ let d=0; for(let j=0;j<w;j++) d+=Q[i][j]*Q[k][j]; R[i][k]+=d; for(let j=0;j<w;j++) Q[k][j]-=d*Q[i][j]; }
    let nn=0; for(let j=0;j<w;j++) nn+=Q[k][j]*Q[k][j]; nn=Math.sqrt(nn); R[k][k]=nn; for(let j=0;j<w;j++) Q[k][j]/=nn; }
  // B = R^{-1} Q^T  (P x w), back substitution
  const B=[]; for(let k=0;k<P;k++) B.push(new Float64Array(w));
  for(let k=P-1;k>=0;k--){ for(let j=0;j<w;j++){ let v=Q[k][j]; for(let i=k+1;i<P;i++) v-=R[k][i]*B[i][j]; B[k][j]=v/R[k][k]; } }
  return {B,s};
}
/* weights evaluating the deriv-th derivative (w.r.t. sample index, then / delta^deriv) of the LS polynomial at position u0 (index) */
function polyEvalWeights(PP,order,deriv,c,pos,delta){
  const w=PP.B[0].length, out=new Float64Array(w); if(deriv>order) return out;
  const u=(pos-c)/PP.s;
  for(let m=deriv;m<=order;m++){ let f=1; for(let q=0;q<deriv;q++) f*=(m-q); const coef=f*Math.pow(u,m-deriv);
    for(let j=0;j<w;j++) out[j]+=coef*PP.B[m][j]; }
  const sc=1/(Math.pow(PP.s,deriv)*Math.pow(delta,deriv)); for(let j=0;j<w;j++) out[j]*=sc;
  return out;
}
function savgolCoeffs(window,order,deriv,delta){
  deriv=deriv||0; delta=delta||1; checkSG(window,order);
  const c=(window-1)/2; return polyEvalWeights(polyProjector(window,order,c),order,deriv,c,c,delta);   // 'dot' order: weight j multiplies y[i-c+j]
}
function checkSG(window,order){ if(!(window%2===1&&window>=1)) throw new Error("savgol: window must be a positive odd integer"); if(!(order>=0&&order<window)) throw new Error("savgol: order must be < window"); }
function savgolRows(n,window,order,opt){
  opt=opt||{}; const deriv=opt.deriv||0, delta=opt.delta||1, mode=opt.mode||"interp";
  checkSG(window,order); const half=(window-1)/2;
  if(mode==="interp"&&window>n) throw new Error("savgol interp: window > number of points");
  const cf=savgolCoeffs(window,order,deriv,delta), rows=new Array(n);
  const isInterior=i=>i-half>=0&&i+half<n;
  for(let i=0;i<n;i++){
    if(mode==="interp"&&!isInterior(i)){ continue; }
    const acc=new Map(); for(let j=0;j<window;j++){ const k=mapIndex(i-half+j,n,mode==="interp"?"constant":mode); if(k<0) continue; acc.set(k,(acc.get(k)||0)+cf[j]); }
    const idx=[...acc.keys()]; rows[i]={idx,w:idx.map(k=>acc.get(k))};
  }
  if(mode==="interp"){
    const PP=polyProjector(window,order,(window-1)/2);
    for(let i=0;i<half;i++){ const wv=polyEvalWeights(PP,order,deriv,(window-1)/2,i,delta); rows[i]={idx:[...Array(window).keys()],w:Array.from(wv)}; }
    for(let i=n-half;i<n;i++){ const wv=polyEvalWeights(PP,order,deriv,(window-1)/2,i-(n-window),delta); rows[i]={idx:[...Array(window).keys()].map(j=>j+n-window),w:Array.from(wv)}; }
  }
  return rows;
}
function savgol(y,window,order,opt){
  opt=opt||{}; y=F64(y); const rows=savgolRows(y.length,window,order,opt);
  const cval=opt.cval||0;
  const r=applyRows(rows,y,opt.var?F64(opt.var):null);
  if(opt.mode==="constant"&&cval!==0){ const cf=savgolCoeffs(window,order,opt.deriv||0,opt.delta||1), half=(window-1)/2, n=y.length;
    for(let i=0;i<n;i++){ let s=0; for(let j=0;j<window;j++){ const k=i-half+j; if(k<0||k>=n) s+=cf[j]*cval; } r.y[i]+=s; } }
  return opt.var?r:r.y;
}
/* Gaussian smoothing = scipy.ndimage.gaussian_filter1d(y, sigma, mode, truncate) (sigma in points) */
function gaussRows(n,sigma,opt){
  opt=opt||{}; const mode=opt.mode||"reflect", trunc=opt.truncate!=null?opt.truncate:4;
  const R=Math.floor(trunc*sigma+0.5), k=new Float64Array(2*R+1); let s=0;
  for(let j=-R;j<=R;j++){ k[j+R]=Math.exp(-0.5*j*j/(sigma*sigma)); s+=k[j+R]; } for(let j=0;j<k.length;j++) k[j]/=s;
  const rows=new Array(n);
  for(let i=0;i<n;i++){ const acc=new Map(); for(let j=-R;j<=R;j++){ const q=mapIndex(i+j,n,mode); if(q<0) continue; acc.set(q,(acc.get(q)||0)+k[j+R]); }
    const idx=[...acc.keys()]; rows[i]={idx,w:idx.map(q=>acc.get(q))}; }
  return {rows,kernel:k,R};
}
function gaussSmooth(y,sigma,opt){
  opt=opt||{}; y=F64(y); if(!(sigma>0)) return opt.var?{y,var:F64(opt.var)}:y;
  const G=gaussRows(y.length,sigma,opt); const r=applyRows(G.rows,y,opt.var?F64(opt.var):null);
  if(opt.mode==="constant"&&opt.cval){ const n=y.length; for(let i=0;i<n;i++){ let s=0; for(let j=-G.R;j<=G.R;j++){ const q=i+j; if(q<0||q>=n) s+=G.kernel[j+G.R]*opt.cval; } r.y[i]+=s; } }
  return opt.var?r:r.y;
}
/* eV-based convenience wrappers (uniform grid required) */
function smooth(E,y,opt){
  opt=opt||{}; const g=gridInfo(E); if(!g.uniform) throw new Error("smooth: resample to a uniform grid first");
  const method=opt.method||"sg";
  if(method==="sg"){ const w=opt.window||Math.max(3,2*Math.round((opt.widthEV||5*g.step)/g.step/2)+1); return savgol(y,w,opt.order!=null?opt.order:2,{mode:opt.mode||"interp",var:opt.var}); }
  if(method==="gauss"){ const sig=(opt.fwhmEV!=null?opt.fwhmEV*FW2SIG/g.step:opt.sigma); return gaussSmooth(y,sig,{mode:opt.mode||"reflect",var:opt.var}); }
  throw new Error("smooth: unknown method");
}
function derivative(E,y,opt){
  opt=opt||{}; const g=gridInfo(E); if(!g.uniform) throw new Error("derivative: resample to a uniform grid first");
  return savgol(y,opt.window||7,opt.order!=null?opt.order:2,{deriv:opt.deriv||1,delta:g.step,mode:opt.mode||"interp",var:opt.var});
}
/* noise-variance reduction factor of a linear filter row (white noise): sum w^2 */
function noiseFactor(rows,i){ const r=rows[i]; let s=0; for(const w of r.w) s+=w*w; return s; }

/* ============================================================
   §4 PSF kernels on a uniform grid of step Δ: h_j = (1/Δ)∫_{(j-1/2)Δ}^{(j+1/2)Δ} psf(x) dx · Δ, j=-R..R, Σh = 1
   (bin-averaged — exact for samples that represent channel averages). Convention: h_j spreads intensity
   from E to E + jΔ.
   ============================================================ */
function binIntegral(f,a,b,minWidth){ const m=Math.max(1,Math.ceil(4*(b-a)/minWidth)); const hw=(b-a)/(2*m); let s=0;
  for(let p=0;p<m;p++){ const c=a+(2*p+1)*hw; for(let q=0;q<8;q++) s+=GL8W[q]*f(c+hw*GL8X[q]); } return s*hw; }
function psfKernel(psf,step,opt){
  opt=opt||{}; if(!psf) throw new Error("psf required");
  const type=psf.type||(psf.h?"array":psf.E?"measured":psf.wL!=null?"voigt":"gauss");
  let h, R, fw, info={};
  if(type==="gauss"){
    const s=psf.fwhm*FW2SIG; R=Math.max(1,Math.ceil((opt.nsig||6)*s/step+0.5)); h=new Float64Array(2*R+1);
    for(let j=-R;j<=R;j++){ const a=((j-0.5)*step)/s, b=((j+0.5)*step)/s;
      h[j+R]= b<=0 ? 0.5*(erfc(-b/SQ2)-erfc(-a/SQ2)) : a>=0 ? 0.5*(erfc(a/SQ2)-erfc(b/SQ2)) : 1-0.5*(erfc(b/SQ2)+erfc(-a/SQ2)); }
    fw=psf.fwhm;
  } else if(type==="voigt"){
    const wG=psf.wG||0, wL=psf.wL||0, s=wG*FW2SIG, g=wL/2;
    fw=0.5346*wL+Math.sqrt(0.2166*wL*wL+wG*wG);
    const tailTol=opt.tailTol||1e-3;
    let X=Math.max((opt.nsig||6)*s, wL>0?2*g/(PI*tailTol):0); X=Math.min(X,opt.maxRadius||Math.max(10*fw,(opt.nsig||6)*s));
    R=Math.max(1,Math.ceil(X/step+0.5)); h=new Float64Array(2*R+1);
    const mw=Math.min(wG>0?wG:Infinity,wL>0?wL:Infinity);
    for(let j=-R;j<=R;j++) h[j+R]=binIntegral(x=>voigt(x,wG,wL),(j-0.5)*step,(j+0.5)*step,mw);
  } else if(type==="measured"){
    const K=psfFromSamples(psf.E,psf.y,step,psf); h=K.h; R=K.R; fw=K.fwhm; info=K.info;
  } else if(type==="array"){
    h=F64(psf.h); if(h.length%2!==1) throw new Error("psf array must have odd length"); R=(h.length-1)/2;
    const x=new Float64Array(h.length); for(let j=0;j<h.length;j++) x[j]=(j-R)*step; fw=fwhm(x,h).fwhm;
  } else throw new Error("unknown psf type "+type);
  let sum=0; for(let j=0;j<h.length;j++) sum+=h[j];
  for(let j=0;j<h.length;j++) h[j]/=sum;
  const x=new Float64Array(h.length); for(let j=0;j<h.length;j++) x[j]=(j-R)*step;
  return {h,R,step,type,fwhm:fw,fwhmSampled:fwhm(x,h).fwhm,mass:sum,info};
}
/* measured PSF: linear baseline through the mean of nb end points, clip < 0, centre at centroid,
   flux-conserving rebin onto offsets jΔ */
function psfFromSamples(E,y,step,o){
  o=o||{}; E=F64(E); y=F64(y); const n=E.length, nb=o.nb||3;
  let x0=0,y0=0,x1=0,y1=0; for(let i=0;i<nb;i++){ x0+=E[i]/nb; y0+=y[i]/nb; x1+=E[n-1-i]/nb; y1+=y[n-1-i]/nb; }
  const z=new Float64Array(n); for(let i=0;i<n;i++) z[i]=Math.max(0,y[i]-(y0+(y1-y0)*(E[i]-x0)/(x1-x0)));
  let a=0,m=0; for(let i=0;i<n;i++){ a+=z[i]; m+=z[i]*E[i]; } const c=m/a;
  const xr=E.map(v=>v-c); const R=Math.ceil(Math.max(-xr[0],xr[n-1])/step);
  const off=new Float64Array(2*R+1); for(let j=-R;j<=R;j++) off[j+R]=j*step;
  const rb=rebinRows(xr,off); const h=applyRows(rb.rows,z,null).y; for(let j=0;j<h.length;j++) if(!isFinite(h[j])) h[j]=0;
  return {h,R,fwhm:fwhm(xr,z).fwhm,info:{centroid:c,baseline:[y0,y1]}};
}

/* ============================================================
   §5 Deconvolution with half-sample-symmetric ("reflect") boundaries:
   (H o)_i = Σ_j h_j o_ρ(i-j),  ρ = scipy 'reflect' index map (d c b a | a b c d | d c b a).
   For a symmetric kernel H is symmetric (Toeplitz + Hankel, diagonalised by the DCT-II) with eigenvalues
   equal to the kernel's DTFT, so Van Cittert/Jansson iterations are stable and Σ(H o) = Σ o exactly.
   Adjoint (needed by RL for asymmetric measured kernels): (Hᵀ r)_k = Σ_i Σ_j h_j r_i [ρ(i-j)=k].
   ============================================================ */
function convRefl(o,h,R,out){
  const M=o.length;
  for(let i=0;i<M;i++){ let s=0;
    if(i-R>=0&&i+R<M){ for(let j=-R;j<=R;j++) s+=h[j+R]*o[i-j]; }
    else for(let j=-R;j<=R;j++) s+=h[j+R]*o[mapIndex(i-j,M,"reflect")];
    out[i]=s; }
  return out;
}
function convReflT(r,h,R,out){
  const M=r.length; out.fill(0);
  for(let i=0;i<M;i++){ const ri=r[i]; if(ri===0) continue;
    if(i-R>=0&&i+R<M){ for(let j=-R;j<=R;j++) out[i-j]+=h[j+R]*ri; }
    else for(let j=-R;j<=R;j++) out[mapIndex(i-j,M,"reflect")]+=h[j+R]*ri; }
  return out;
}
/* §5.4 noise-only modes. For a symmetric kernel with reflect boundaries the eigenvectors of H are the DCT-II
   vectors, λ_m = Σ_j h_j cos(π m j / N). Let P_m = Y_m² / mean(var) (Y = orthonormal DCT-II of y; E P_m = 1 for
   pure noise) and P̄_m its centred running mean over w = max(5, 2⌊N/100⌋+1) modes. The noise-only set is
   ⊥ = {m : |λ_m| < ε₀} ∪ {m ≥ m_c : |λ_m| < ½}, m_c = first m with P̄_m ≤ 1.5 (Picard-type split).
   Noise in ⊥ cannot be fitted by H o, so the discrepancy target is χ²_k ≤ χ²_⊥ + τ (N − N_⊥), χ²_⊥ = Σ_{m∈⊥} P_m
   ("projected discrepancy"). Asymmetric kernels use the symmetrised kernel. */
function passband(y,K,v,eps){
  const N=y.length, h=K.h, R=K.R; eps=eps!=null?eps:0.01;
  let vm=0; for(let i=0;i<N;i++) vm+=v[i]; vm/=N;
  const lam=new Float64Array(N), P=new Float64Array(N);
  for(let m=0;m<N;m++){ let l=0; for(let j=-R;j<=R;j++) l+=0.5*(h[j+R]+h[R-j])*Math.cos(PI*m*j/N); lam[m]=l;
    let Y=0; const w=PI*m/N; for(let i=0;i<N;i++) Y+=y[i]*Math.cos(w*(i+0.5)); Y*=m===0?Math.sqrt(1/N):Math.sqrt(2/N); P[m]=Y*Y/vm; }
  const hw=Math.max(2,Math.floor(N/100)); let mc=N;
  for(let m=0;m<N;m++){ let s=0,c=0; for(let j=Math.max(0,m-hw);j<=Math.min(N-1,m+hw);j++){ s+=P[j]; c++; } if(s/c<=1.5){ mc=m; break; } }
  let nPerp=0, chiPerp=0;
  for(let m=0;m<N;m++) if(Math.abs(lam[m])<eps||(m>=mc&&Math.abs(lam[m])<0.5)){ nPerp++; chiPerp+=P[m]; }
  return {lam,P,mc,nPerp,chiPerp,eps,meanVar:vm};
}
/* core iteration — fixed or discrepancy-stopped. y: data (N), v: its variance */
function iterate(ypad,K,v,o){
  const N=ypad.length;
  const M=ypad.length, h=K.h, R=K.R, method=o.method||"jansson";
  const ob=new Float64Array(M), Ho=new Float64Array(M), tmp=new Float64Array(M), Ht1=new Float64Array(M);
  let a=o.a, b=o.b; const r0=o.r0!=null?o.r0:(method==="jansson"?1.5:1);
  if(method==="rl"){ let ym=0; for(let i=0;i<M;i++) ym=Math.max(ym,ypad[i]); const fl=1e-6*ym;
    for(let i=0;i<M;i++) ob[i]=Math.max(ypad[i],fl); tmp.fill(1); convReflT(tmp,h,R,Ht1); }
  else if(method==="jansson"){ for(let i=0;i<M;i++) ob[i]=Math.min(b,Math.max(a,ypad[i])); }
  else { for(let i=0;i<M;i++) ob[i]=o.positive?Math.max(0,ypad[i]):ypad[i]; }
  const chi2=()=>{ convRefl(ob,h,R,Ho); let c=0; for(let i=0;i<N;i++){ const d=ypad[i]-Ho[i]; c+=d*d/v[i]; } return c; };
  const hist=[chi2()]; let k=0, stop="fixed";
  const target=o.target!=null?o.target:(o.tau!=null?o.tau:1)*N, fixed=o.iter!=null, maxIter=fixed?o.iter:(o.maxIter||20000);
  let rising=0;
  if(!fixed&&hist[0]<=target) stop="discrepancy";
  else for(;;){
    if(k>=maxIter){ stop=fixed?"fixed":"maxIter"; break; }
    if(method==="rl"){ for(let i=0;i<M;i++) tmp[i]=Ho[i]>1e-300?ypad[i]/Ho[i]:1; convReflT(tmp,h,R,Ho); for(let i=0;i<M;i++) ob[i]*=Ho[i]/Ht1[i]; }
    else if(method==="jansson"){ const mid=(a+b)/2, wd=b-a;
      if(o.reblur!==false){ for(let i=0;i<M;i++) tmp[i]=ypad[i]-Ho[i]; convReflT(tmp,h,R,Ho); }   // Kawata–Ichioka reblurring: Hᵀ(y − H o)
      else for(let i=0;i<M;i++) Ho[i]=ypad[i]-Ho[i];
      for(let i=0;i<M;i++){ const r=Math.max(0,r0*(1-2*Math.abs(ob[i]-mid)/wd)); let q=ob[i]+r*Ho[i]; ob[i]=q<a?a:q>b?b:q; } }
    else if(method==="landweber"){ for(let i=0;i<M;i++) tmp[i]=ypad[i]-Ho[i]; convReflT(tmp,h,R,Ho);
      for(let i=0;i<M;i++){ let q=ob[i]+r0*Ho[i]; if(o.positive&&q<0) q=0; ob[i]=q; } }
    else if(method==="vancittert"){ for(let i=0;i<M;i++){ let q=ob[i]+r0*(ypad[i]-Ho[i]); if(o.positive&&q<0) q=0; ob[i]=q; } }
    else throw new Error("deconvolve: unknown method "+method);
    k++; const c=chi2(); hist.push(c);
    if(!fixed){ if(c<=target){ stop="discrepancy"; break; }
      if(c>=hist[k-1]){ if(++rising>=20){ stop="stalled"; break; } } else rising=0;
      if(k>=50&&hist[k>>1]-c<1){ stop="plateau"; break; } }
  }
  return {o:ob,Ho:F64(Ho),k,chi2:hist[hist.length-1],hist,stop};
}
function deconvolve(E,y,opt){
  opt=opt||{}; E=F64(E); y=F64(y); const N=y.length, gi=gridInfo(E), warnings=[];
  if(!gi.uniform) throw new Error("deconvolve: E must be ascending and uniform — resample() first");
  const step=gi.step, method=opt.method||"jansson";
  const K=opt.kernel||psfKernel(opt.psf,step,opt.psfOpt);
  const nv=resolveVar(y,opt.noise!=null?opt.noise:opt.var);
  let ymax=-Infinity, ymin=Infinity; for(let i=0;i<N;i++){ ymax=Math.max(ymax,y[i]); ymin=Math.min(ymin,y[i]); }
  const bounds=opt.bounds||{}; const a=bounds.a!=null?bounds.a:0, b=bounds.b!=null?bounds.b:(opt.bFactor||4)*ymax;
  if(ymin<0&&(method==="rl"||method==="jansson")) warnings.push("data contain negative values; they are clipped by the positivity constraint");
  if(K.fwhm<2*step) warnings.push("PSF FWHM < 2 grid steps: the PSF is under-sampled; deconvolution gain is not meaningful");
  if(K.fwhm>(E[N-1]-E[0])/4) warnings.push("PSF FWHM > 1/4 of the region: edge effects dominate");
  const io={method,a,b,r0:opt.r0,tau:opt.tau,maxIter:opt.maxIter,iter:opt.iter,positive:opt.positive!==false,reblur:opt.reblur};
  if(method==="vancittert"||(method==="jansson"&&opt.reblur===false)) warnings.push("iteration without reblurring amplifies noise in frequency bands where the PSF transfer is ~0 (oversampled data); prefer reblur (default) or rl");
  const tau=opt.tau!=null?opt.tau:1, rule=opt.stopRule||"projected";
  let pb=null; if(rule==="projected"){ pb=passband(y,K,nv.var,opt.passTol); io.target=pb.chiPerp+tau*(N-pb.nPerp); } else io.target=tau*N;
  const it=iterate(y,K,nv.var,io);
  const o=it.o, model=it.Ho, resid=new Float64Array(N);
  for(let i=0;i<N;i++) resid[i]=(y[i]-model[i])/Math.sqrt(nv.var[i]);
  if(it.stop==="maxIter") warnings.push("discrepancy level not reached in "+it.k+" iterations (PSF too wide for this noise, or noise under-estimated)");
  if(it.stop==="plateau") warnings.push("discrepancy level not reached; stopped where doubling the iteration count lowered χ² by < 1 (k = "+it.k+")");
  if(it.stop==="stalled") warnings.push("χ² stopped decreasing before reaching the discrepancy level");
  if(it.k===0&&it.stop==="discrepancy") warnings.push("data already consistent with the PSF-broadened model at the noise level (k = 0): nothing to restore");
  if(method==="jansson"){ let mx=0; for(let i=0;i<N;i++) mx=Math.max(mx,o[i]); if(mx>0.95*b) warnings.push("upper Jansson bound active (max o = "+mx.toPrecision(4)+" ≥ 0.95 b); raise bounds.b / bFactor"); }
  const areaIn=trapz(E,y), areaOut=trapz(E,o);
  const res={E,o,model,resid,k:it.k,chi2:it.chi2,chi2Hist:it.hist,N,target:io.target,stopRule:rule,passband:pb?{nPerp:pb.nPerp,chi2Perp:pb.chiPerp,mc:pb.mc,eps:pb.eps}:null,stop:it.stop,method,
    var:nv.var,noiseKind:nv.kind,bounds:{a,b},kernel:{h:K.h,R:K.R,fwhm:K.fwhm,fwhmSampled:K.fwhmSampled,type:K.type},
    edgePoints:K.R,areaIn,areaOut,areaRatio:areaOut/areaIn,warnings};
  // effective resolution: point-source response of the same k-iteration restoration at the data's level
  if(opt.gainTest!==false&&it.k>0){
    const T=Math.max(4*K.R,20), L=2*T+1, ped=median(y); let hm=0; for(const v of K.h) hm=Math.max(hm,v);
    const A0=(ymax-ped)/hm; const s=new Float64Array(L).fill(ped); for(let j=-K.R;j<=K.R;j++) s[T+j]+=A0*K.h[j+K.R];
    const vt=new Float64Array(L).fill(1);
    const r=iterate(s,K,vt,Object.assign({},io,{iter:it.k}));
    const ot=r.o, xt=new Float64Array(L); for(let i=0;i<L;i++) xt[i]=(i-T)*step;
    const fe=fwhm(xt,ot,ped).fwhm;
    res.resolution={fwhmPSF:K.fwhmSampled,fwhmEff:fe,gain:K.fwhmSampled/fe,limitedBySampling:fe<1.5*step};
  } else res.resolution={fwhmPSF:K.fwhmSampled,fwhmEff:K.fwhmSampled,gain:1};
  // noise amplification: Monte-Carlo re-deconvolution of model + N(0,var) with the same k
  const M=opt.mc!=null?opt.mc:6;
  if(M>=2&&it.k>0){
    const g=rng(opt.seed||12345), acc=new Float64Array(N), acc2=new Float64Array(N);
    for(let m=0;m<M;m++){ const ys=new Float64Array(N); for(let i=0;i<N;i++) ys[i]=model[i]+Math.sqrt(nv.var[i])*g.normal();
      const r=iterate(ys,K,nv.var,Object.assign({},io,{iter:it.k})); for(let i=0;i<N;i++){ const q=r.o[i]; acc[i]+=q; acc2[i]+=q*q; } }
    const sd=new Float64Array(N); let so=0, sv=0, ro=0, ry=0;
    for(let i=0;i<N;i++){ const mu=acc[i]/M; sd[i]=Math.sqrt(Math.max(0,(acc2[i]-M*mu*mu)/(M-1))); if(i>=K.R&&i<N-K.R){ so+=sd[i]*sd[i]; sv+=nv.var[i]; } }
    for(let i=1;i<N-1;i++){ const d=o[i-1]-2*o[i]+o[i+1], e=y[i-1]-2*y[i]+y[i+1]; ro+=d*d; ry+=e*e; }
    res.noise={sigma:sd,gain:Math.sqrt(so/sv),roughness:Math.sqrt(ro/ry),mc:M};
    const lim=opt.maxNoiseGain||3;
    if(res.noise.gain>lim) warnings.push("noise amplification ×"+res.noise.gain.toFixed(2)+" exceeds "+lim+": do not interpret fine structure; use fewer iterations or a narrower PSF");
  }
  return res;
}

/* ============================================================
   Levenberg–Marquardt (Marquardt 1963) with diagonal scaling, box bounds by clamping
   model(p) -> Float64Array(N); jac(p, f) -> array of P columns (optional; central differences otherwise)
   ============================================================ */
function cholSolve(A,b){ const n=b.length, L=A.map(r=>Float64Array.from(r));
  for(let j=0;j<n;j++){ let s=L[j][j]; for(let k=0;k<j;k++) s-=L[j][k]*L[j][k]; if(!(s>0)) return null; L[j][j]=Math.sqrt(s);
    for(let i=j+1;i<n;i++){ let t=L[i][j]; for(let k=0;k<j;k++) t-=L[i][k]*L[j][k]; L[i][j]=t/L[j][j]; } }
  const z=new Float64Array(n); for(let i=0;i<n;i++){ let s=b[i]; for(let k=0;k<i;k++) s-=L[i][k]*z[k]; z[i]=s/L[i][i]; }
  const x=new Float64Array(n); for(let i=n-1;i>=0;i--){ let s=z[i]; for(let k=i+1;k<n;k++) s-=L[k][i]*x[k]; x[i]=s/L[i][i]; } return x; }
function invSPD(A){ const n=A.length, out=[]; for(let i=0;i<n;i++) out.push(new Float64Array(n));
  for(let c=0;c<n;c++){ const e=new Float64Array(n); e[c]=1; const x=cholSolve(A,e); if(!x) return null; for(let r=0;r<n;r++) out[r][c]=x[r]; } return out; }
function lm(model,jac,p0,y,w,opt){
  opt=opt||{}; const Pn=p0.length, N=y.length, free=[]; for(let j=0;j<Pn;j++) if(!(opt.fix&&opt.fix[j])) free.push(j);
  const lo=opt.lo||[], hi=opt.hi||[], scale=opt.scale||new Array(Pn).fill(1);
  const clampP=p=>{ for(const j of free){ if(lo[j]!=null&&p[j]<lo[j]) p[j]=lo[j]; if(hi[j]!=null&&p[j]>hi[j]) p[j]=hi[j]; } return p; };
  const chi=f=>{ let c=0; for(let i=0;i<N;i++){ const d=y[i]-f[i]; c+=w[i]*d*d; } return c; };
  const J=(p,f)=>{ if(jac) return jac(p,f);
    const cols=[]; for(let j=0;j<Pn;j++){ if(opt.fix&&opt.fix[j]){ cols.push(null); continue; } const hh=1e-6*(Math.abs(p[j])+scale[j]);
      const pp=p.slice(), pm=p.slice(); pp[j]+=hh; pm[j]-=hh; const fp=model(pp), fm=model(pm); const c=new Float64Array(N); for(let i=0;i<N;i++) c[i]=(fp[i]-fm[i])/(2*hh); cols.push(c); } return cols; };
  let p=clampP(p0.slice()), f=model(p), c=chi(f), lam=1e-3, it=0, small=0;
  const normal=(p,f)=>{ const cols=J(p,f), Fn=free.length, A=[], g=new Float64Array(Fn);
    for(let a=0;a<Fn;a++){ A.push(new Float64Array(Fn)); const ca=cols[free[a]]; let s=0; for(let i=0;i<N;i++) s+=w[i]*ca[i]*(y[i]-f[i]); g[a]=s;
      for(let b=0;b<=a;b++){ const cb=cols[free[b]]; let t=0; for(let i=0;i<N;i++) t+=w[i]*ca[i]*cb[i]; A[a][b]=t; } }
    for(let a=0;a<Fn;a++) for(let b=a+1;b<Fn;b++) A[a][b]=A[b][a]; return {A,g}; };
  let NE=normal(p,f);
  const maxIter=opt.maxIter||300;
  for(;it<maxIter;it++){
    let accepted=false;
    for(let tries=0;tries<40;tries++){
      const Fn=free.length, Ad=NE.A.map((r,a)=>{ const q=Float64Array.from(r); q[a]+=lam*Math.max(r[a],1e-300); return q; });
      const d=cholSolve(Ad,NE.g); if(!d){ lam*=10; continue; }
      const pn=p.slice(); for(let a=0;a<Fn;a++) pn[free[a]]+=d[a]; clampP(pn);
      const fn=model(pn), cn=chi(fn);
      if(isFinite(cn)&&cn<=c){ let mx=0; for(let a=0;a<Fn;a++) mx=Math.max(mx,Math.abs(pn[free[a]]-p[free[a]])/(Math.abs(p[free[a]])+scale[free[a]]));
        const rel=(c-cn)/Math.max(c,1e-300); p=pn; f=fn; c=cn; lam=Math.max(lam/3,1e-12); accepted=true;
        if(rel<1e-13&&mx<1e-10) small++; else small=0; break; }
      lam*=4;
    }
    if(!accepted||small>=2) break;
    NE=normal(p,f);
  }
  NE=normal(p,f); const cov=invSPD(NE.A);
  return {p,f,chi2:c,iter:it,free,cov,dof:N-free.length};
}

/* ============================================================
   §6 Fermi edge:  y(E) = b0 + b1 (E-EF) + ∫ (A + S (u-EF)) f(u) G_σ(E-u) du
   f(u) = 1/(1+exp(-s (u-EF)/kT)),  s = +1 for BE axis (occupied side at E > EF), -1 for KE.
   Evaluated as closed-form step part + Gauss–Legendre integral of the localised remainder f - H.
   ============================================================ */
function fermiBasis(E,EF,wG,kT,s){
  const N=E.length, sig=Math.max(wG,1e-6)*FW2SIG, MA=new Float64Array(N), MS=new Float64Array(N);
  const V=36*kT, wp=Math.min(kT,sig), np=Math.ceil(V/wp), hw=V/(2*np);
  // remainder nodes v>0 (and mirror v<0): e(v) = -s*sign(v)*L(|v|), L(t)=1/(1+exp(t/kT))
  const nv=np*8, vn=new Float64Array(nv), wn=new Float64Array(nv), Ln=new Float64Array(nv);
  for(let p=0;p<np;p++) for(let q=0;q<8;q++){ const k=p*8+q, v=(2*p+1)*hw+hw*GL8X[q]; vn[k]=v; wn[k]=hw*GL8W[q]; Ln[k]=1/(1+Math.exp(v/kT)); }
  const inv=1/(sig*SQ2PI), i2=1/(2*sig*sig);
  for(let i=0;i<N;i++){
    const d=E[i]-EF, z=d/sig; let a, b;
    if(s>0){ const P=Phi(z); a=P; b=d*P+sig*phi(z); } else { const P=Phi(-z); a=P; b=d*P-sig*phi(z); }
    if(Math.abs(d)<V+9*sig){ let ra=0, rb=0;
      for(let k=0;k<nv;k++){ const v=vn[k], gp=Math.exp(-(d-v)*(d-v)*i2), gm=Math.exp(-(d+v)*(d+v)*i2), L=Ln[k]*wn[k];
        // v>0: e = -s L ; v<0 (−v): e = +s L
        ra+=L*(-gp+gm); rb+=L*(-v*gp-v*gm); }
      a+=s*ra*inv; b+=s*rb*inv; }
    MA[i]=a; MS[i]=b; }
  return {MA,MS};
}
function fermiModel(E,p,opt){
  opt=opt||{}; E=F64(E); const kT=KB_EV*(p.T!=null?p.T:(opt.T!=null?opt.T:300)), s=(opt.axis||"BE")==="KE"?-1:1;
  const B=fermiBasis(E,p.EF,p.wG,kT,s), out=new Float64Array(E.length);
  for(let i=0;i<E.length;i++) out[i]=(p.b0||0)+(p.b1||0)*(E[i]-p.EF)+p.A*B.MA[i]+(p.S||0)*B.MS[i];
  return out;
}
const FERMI_P=["EF","wG","A","S","b0","b1"];
function fermiGuess(E,y,s){
  const N=E.length, gi=gridInfo(E); let ys=y;
  if(gi.uniform&&N>=9) ys=savgol(y,Math.min(N-(N%2===0?1:0),2*Math.max(2,Math.round(N/40))+1),2,{mode:"interp"});
  let d=0, im=1; for(let i=1;i<N-1;i++){ const g=s*(ys[i+1]-ys[i-1]); if(g>d){ d=g; im=i; } }
  const EF=E[im]; const k=Math.max(2,Math.round(N/10));
  let lo=0, hi=0; for(let i=0;i<k;i++){ lo+=y[i]/k; hi+=y[N-1-i]/k; }
  const occ=s>0?hi:lo, emp=s>0?lo:hi; const A=occ-emp;
  const slope=d/(E[1]-E[0])/2; const wG=Math.max(Math.abs(A)/Math.max(slope,1e-300)/1.0645*0.94,2*Math.abs(E[1]-E[0]));
  return {EF,wG:Math.min(wG,(E[N-1]-E[0])/3),A,S:0,b0:emp,b1:0};
}
function fermiFit(E,y,opt){
  opt=opt||{}; E=F64(E); y=F64(y); const N=E.length, T=opt.T!=null?opt.T:300, kT=KB_EV*T, s=(opt.axis||"BE")==="KE"?-1:1;
  const g0=Object.assign(fermiGuess(E,y,s),opt.p0||{});
  const p0=FERMI_P.map(k=>g0[k]); const fix=FERMI_P.map(k=>!!(opt.fix&&opt.fix[k]));
  const span=E[N-1]-E[0], step=Math.abs(E[1]-E[0]);
  const lo=[E[0],Math.max(1e-4,step/20),null,null,null,null], hi=[E[N-1],span,null,null,null,null];
  const yr=Math.max(1,Math.abs(g0.A));
  const scale=[step,step,yr,yr/span,yr,yr/span];
  const model=p=>{ const B=fermiBasis(E,p[0],p[1],kT,s), f=new Float64Array(N); for(let i=0;i<N;i++) f[i]=p[4]+p[5]*(E[i]-p[0])+p[2]*B.MA[i]+p[3]*B.MS[i]; return f; };
  const jac=p=>{ const B=fermiBasis(E,p[0],p[1],kT,s), cols=new Array(6);
    cols[2]=B.MA; cols[3]=B.MS; cols[4]=new Float64Array(N).fill(1); cols[5]=E.map(v=>v-p[0]);
    [0,1].forEach(j=>{ const hh=1e-6*(Math.abs(p[j])+scale[j]), pp=p.slice(), pm=p.slice(); pp[j]+=hh; pm[j]-=hh; const fp=model(pp), fm=model(pm), c=new Float64Array(N); for(let i=0;i<N;i++) c[i]=(fp[i]-fm[i])/(2*hh); cols[j]=c; });
    return cols; };
  const noise=opt.noise!=null?opt.noise:(opt.var||"auto");
  let w, scaleCov=false, r;
  if(noise&&noise.length===N){ w=F64(noise).map(v=>1/v); r=lm(model,jac,p0,y,w,{fix,lo,hi,scale}); }
  else { // Poisson-shaped weights: pass 1 with data, pass 2 with model (Poisson ML approx)
    w=y.map(v=>1/Math.max(v,1)); r=lm(model,jac,p0,y,w,{fix,lo,hi,scale});
    w=r.f.map(v=>1/Math.max(v,1)); r=lm(model,jac,r.p,y,w,{fix,lo,hi,scale});
    scaleCov=noise!=="poisson"; }
  const red=r.chi2/r.dof, fac=scaleCov?red:1;
  const p={}, sg={}, ci={}; FERMI_P.forEach((k,j)=>{ p[k]=r.p[j]; const a=r.free.indexOf(j); sg[k]=a>=0&&r.cov?Math.sqrt(r.cov[a][a]*fac):0; ci[k]=[p[k]-Z975*sg[k],p[k]+Z975*sg[k]]; });
  const cov=r.cov?r.cov.map(row=>Array.from(row).map(v=>v*fac)):null;
  const warnings=[];
  if(p.wG<2*step) warnings.push("fitted Gaussian width < 2 grid steps: resolution not resolved by the sampling");
  if(p.EF-E[0]<2*p.wG||E[N-1]-p.EF<2*p.wG) warnings.push("edge within 2·FWHM of the region boundary: widen the region");
  if(scaleCov&&red>3) warnings.push("reduced χ² = "+red.toFixed(2)+": model (linear DOS) may be inadequate");
  if(!r.cov) warnings.push("singular normal matrix: uncertainties not available");
  return {p,sigma:sg,ci95:ci,cov,names:FERMI_P.filter((k,j)=>!fix[j]),chi2:r.chi2,dof:r.dof,redchi2:red,covScaled:scaleCov,
    model:r.f,T,kT,axis:s>0?"BE":"KE",iter:r.iter,
    EF:{value:p.EF,sigma:sg.EF,ci95:ci.EF},resolution:{fwhm:p.wG,sigma:sg.wG,ci95:ci.wG},psf:{type:"gauss",fwhm:p.wG},warnings};
}

/* ============================================================
   §7 PSF from a reference core line (e.g. Ag 3d5/2): y = A V(E-E0; wG, wL) + b0 + b1 (E-E0)
   The instrumental PSF is the Gaussian part wG; wL (core-hole lifetime) must be supplied or fitted.
   ============================================================ */
const PEAK_P=["E0","A","wG","wL","b0","b1"];
function peakModel(E,p){ const f=new Float64Array(E.length); for(let i=0;i<E.length;i++) f[i]=p[1]*voigt(E[i]-p[0],p[2],p[3])+p[4]+p[5]*(E[i]-p[0]); return f; }
function peakPSF(E,y,opt){
  opt=opt||{}; E=F64(E); y=F64(y); const N=E.length;
  const k=Math.max(2,Math.round(N/12)); let lo=0, hi=0; for(let i=0;i<k;i++){ lo+=y[i]/k; hi+=y[N-1-i]/k; }
  const b0=(lo+hi)/2, b1=(hi-lo)/(E[N-1]-E[0]); let im=0; for(let i=1;i<N;i++) if(y[i]>y[im]) im=i;
  const fw=fwhm(E,y,b0).fwhm||((E[N-1]-E[0])/5), wL0=opt.wL!=null?opt.wL:0.3*fw;
  const wG0=Math.max(1e-3,Math.sqrt(Math.max(fw*fw-wL0*fw,fw*fw*0.25)));
  const A0=(y[im]-b0)/voigt(0,wG0,wL0);
  const p0=[E[im],A0,wG0,wL0,b0,b1].map((v,j)=>opt.p0&&opt.p0[PEAK_P[j]]!=null?opt.p0[PEAK_P[j]]:v);
  const fix=[false,false,false,!opt.fitWL,false,false];
  const step=Math.abs(E[1]-E[0]), yr=Math.max(1,y[im]-b0);
  const lo2=[E[0],0,1e-4,0,null,null], hi2=[E[N-1],null,E[N-1]-E[0],E[N-1]-E[0],null,null];
  const scale=[step,yr*fw,step,step,yr,yr];
  const noise=opt.noise!=null?opt.noise:"auto"; let w, scaleCov=true, r;
  const model=p=>peakModel(E,p);
  if(noise&&noise.length===N){ w=F64(noise).map(v=>1/v); scaleCov=false; r=lm(model,null,p0,y,w,{fix,lo:lo2,hi:hi2,scale}); }
  else { w=y.map(v=>1/Math.max(v,1)); r=lm(model,null,p0,y,w,{fix,lo:lo2,hi:hi2,scale}); w=r.f.map(v=>1/Math.max(v,1)); r=lm(model,null,r.p,y,w,{fix,lo:lo2,hi:hi2,scale}); scaleCov=noise!=="poisson"; }
  const red=r.chi2/r.dof, fac=scaleCov?red:1, p={}, sg={};
  PEAK_P.forEach((k,j)=>{ p[k]=r.p[j]; const a=r.free.indexOf(j); sg[k]=a>=0&&r.cov?Math.sqrt(r.cov[a][a]*fac):0; });
  const warnings=[]; if(opt.fitWL) warnings.push("wG and wL fitted together are strongly correlated; prefer a fixed literature wL");
  warnings.push("Voigt model ignores the metallic (Doniach–Šunjić) asymmetry of the reference line: wG is an upper estimate if the line is asymmetric");
  return {p,sigma:sg,chi2:r.chi2,dof:r.dof,redchi2:red,model:r.f,psf:{type:"gauss",fwhm:p.wG},psfVoigt:{type:"voigt",wG:p.wG,wL:p.wL},warnings};
}

const XPS_PROC={
  version:"1.0",
  // utilities
  erf,erfc,Phi,faddeeva,gauss,lorentz,voigt,lgamma,rng,gridInfo,fwhm,trapz,median,mapIndex,KB_EV,
  // §1
  noiseModel,resolveVar,
  // §2
  resample,rebinRows,pchipRows,pchipDeriv,lanczosRows,linearRows,applyRows,uniformGrid,edgesOf,aliasLevel,suggestStep,
  // §3
  savgol,savgolCoeffs,savgolRows,gaussSmooth,gaussRows,smooth,derivative,noiseFactor,
  // §4–5
  psfKernel,psfFromSamples,passband,convRefl,convReflT,iterate,deconvolve,
  // fitting
  lm,fermiModel,fermiBasis,fermiFit,fermiGuess,peakModel,peakPSF
};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_PROC; else root.XPS_PROC=XPS_PROC;
})(typeof self!=="undefined"?self:this);
