/* ============================================================
   INNO Analytics — XPS angle-resolved (ARXPS) engine
   DOM-free (browser, Web Worker, node). Load after xps_multi.js (uses its small linear algebra).
   Equations, citations and validation: modules/xps/SPEC_multi.md §4.
   Angles θ are emission (take-off) angles measured FROM THE SURFACE NORMAL, in degrees. Lengths in nm.
   ============================================================ */
(function(root){
"use strict";
function getMulti(){ if(root&&root.XPS_MULTI) return root.XPS_MULTI; if(typeof require==="function"){ try{ return require("./xps_multi.js"); }catch(e){} } throw new Error("XPS_MULTI not loaded"); }
const DEG=Math.PI/180;
const cosd=t=>Math.cos(t*DEG);

/* ============================================================
   1) Attenuation lengths (SPEC_multi §4.1)
   ============================================================ */
/* mean atomic spacing a (nm) from molar mass M (g/mol), density ρ (g/cm³) and atoms per formula unit n: a³ = M/(ρ N_A n) */
function atomSize(M,rho,n){ return 1e7*Math.cbrt(M/(rho*6.02214076e23*(n||1))); }
/* Seah & Dench, Surf. Interface Anal. 1 (1979) 2, attenuation length (Briggs & Seah, Practical Surface Analysis Vol. 1):
   elements:            λ_m = 538/E² + 0.41 (aE)^½      [monolayers]   λ = a·λ_m (nm)
   inorganic compounds: λ_m = 2170/E² + 0.72 (aE)^½
   organic compounds:   λ_d = 49/E² + 0.11 E^½          [mg m⁻²]       λ = 1000 λ_d / ρ (ρ in kg m⁻³) nm   (E in eV, a in nm) */
function seahDench(E,o){
  o=o||{}; const kind=o.kind||"element";
  if(kind==="organic"){ const rho=o.rhoKgM3||(o.rho?o.rho*1000:null); if(!(rho>0)) throw new Error("organic: density needed"); return 1000*(49/(E*E)+0.11*Math.sqrt(E))/rho; }
  const a=o.a; if(!(a>0)) throw new Error("atom size a (nm) needed");
  const lm=kind==="inorganic"?2170/(E*E)+0.72*Math.sqrt(a*E):538/(E*E)+0.41*Math.sqrt(a*E);
  return a*lm;
}
/* practical effective attenuation length from IMFP and single-scattering albedo ω (Jablonski & Powell; NIST SRD 82 user guide eq. E.3):
   L = λ_in (1 − 0.735 ω) */
function ealFromImfp(lam,omega){ return lam*(1-0.735*(omega||0)); }

/* ============================================================
   2) Single overlayer (SPEC_multi §4.2)
   I_o = I_o∞ (1 − exp(−d/(L_o cosθ))),  I_s = I_s∞ exp(−d/(L_s cosθ)),  R = I_o/I_s,  R∞ = I_o∞/I_s∞
   (Hill, Royce, Fadley, Wagner, Grunthaner, Chem. Phys. Lett. 44 (1976) 225)
   ============================================================ */
function overlayerR(d,th,o){ const c=cosd(th); return o.Rinf*(1-Math.exp(-d/(o.Lo*c)))*Math.exp(d/(o.Ls*c)); }
function overlayerdR(d,th,o){ const c=cosd(th), eo=Math.exp(-d/(o.Lo*c)), es=Math.exp(d/(o.Ls*c)); return o.Rinf*(eo/(o.Lo*c)*es+(1-eo)*es/(o.Ls*c)); }
/* exact solution for one angle (the equation the Thickogram (Cumpson, Surf. Interface Anal. 29 (2000) 403) solves graphically) */
function overlayerThickness(R,th,o){
  if(!(R>0)) return 0; const c=cosd(th);
  if(Math.abs(o.Lo-o.Ls)<=1e-12*o.Lo) return o.Lo*c*Math.log(1+R/o.Rinf);
  let lo=0, hi=o.Lo*c; while(overlayerR(hi,th,o)<R&&hi<1e6) hi*=2;
  for(let i=0;i<200;i++){ const m=0.5*(lo+hi); if(overlayerR(m,th,o)<R) lo=m; else hi=m; if(hi-lo<=1e-15*hi) break; }
  return 0.5*(lo+hi);
}
/* weighted least-squares fit of d (and optionally R∞) to R(θ) at several angles. sR: σ of each ratio (optional; unit weights if absent).
   Returns d, σ_d (from (JᵀWJ)⁻¹), σ_d·sqrt(φ), χ², per-angle exact d. Gauss–Newton with step halving. */
function fitOverlayer(th,R,sR,o){
  o=Object.assign({},o); const n=th.length; const w=Float64Array.from({length:n},(_,i)=>sR&&sR[i]>0?1/(sR[i]*sR[i]):1);
  const perAngle=Float64Array.from(th,(t,i)=>overlayerThickness(R[i],t,o));
  let d=Math.max(1e-6,perAngle.reduce((a,b)=>a+b,0)/n), Ri=o.Rinf; const fitR=!!o.fitRinf;
  const chiOf=(dd,rr)=>{ const q=Object.assign({},o,{Rinf:rr}); let s=0; for(let i=0;i<n;i++){ const e=R[i]-overlayerR(dd,th[i],q); s+=w[i]*e*e; } return s; };
  let chi=chiOf(d,Ri), it=0;
  for(;it<200;it++){
    const q=Object.assign({},o,{Rinf:Ri}); const p=fitR?2:1; const JTJ=new Float64Array(p*p), JTr=new Float64Array(p);
    for(let i=0;i<n;i++){ const m=overlayerR(d,th[i],q), r=R[i]-m; const J=[overlayerdR(d,th[i],q)]; if(fitR) J.push(m/Ri);
      for(let a=0;a<p;a++){ JTr[a]+=w[i]*J[a]*r; for(let b=0;b<p;b++) JTJ[a*p+b]+=w[i]*J[a]*J[b]; } }
    const M=getMulti(); const st=M.solveGE(JTJ,JTr,p); if(!st) break;
    let lam=1, nd, nr, nc; do{ nd=Math.max(0,d+lam*st[0]); nr=fitR?Math.max(1e-12,Ri+lam*st[1]):Ri; nc=chiOf(nd,nr); lam/=2; }while(nc>chi&&lam>1e-10);
    const rel=Math.abs(chi-nc)/Math.max(chi,1e-300); d=nd; Ri=nr; const imp=chi-nc; chi=nc; if(Math.abs(st[0])<1e-12*Math.max(d,1e-12)||(imp>=0&&rel<1e-15)) break;
  }
  const q=Object.assign({},o,{Rinf:Ri}); const p=fitR?2:1; const JTJ=new Float64Array(p*p);
  for(let i=0;i<n;i++){ const m=overlayerR(d,th[i],q); const J=[overlayerdR(d,th[i],q)]; if(fitR) J.push(m/Ri); for(let a=0;a<p;a++) for(let b=0;b<p;b++) JTJ[a*p+b]+=w[i]*J[a]*J[b]; }
  const C=getMulti().invSPD(JTJ,p); const dof=Math.max(n-p,1), phi=chi/dof;
  const sd=C?Math.sqrt(C[0]):NaN, sRi=fitR&&C?Math.sqrt(C[3]):null;
  return {d,sd,sdScaled:sd*Math.sqrt(Math.max(phi,1)),Rinf:Ri,sRinf:sRi,chi2:chi,dof,phi,perAngle,iter:it,
    warn:th.some(t=>t>65)?["angles > 65° from normal: elastic scattering makes the straight-line model increasingly biased"]:[]};
}
/* Hill plot: y = ln(1 + R/R∞) vs x = 1/cosθ, slope = d/L (L_o = L_s = L), weighted fit through the origin */
function hillPlot(th,R,sR,o){
  const x=th.map(t=>1/cosd(t)), y=R.map(r=>Math.log(1+r/o.Rinf)); const w=th.map((_,i)=>{ if(!(sR&&sR[i]>0)) return 1; const sy=sR[i]/(o.Rinf+R[i]); return 1/(sy*sy); });
  let sxy=0, sxx=0; for(let i=0;i<x.length;i++){ sxy+=w[i]*x[i]*y[i]; sxx+=w[i]*x[i]*x[i]; } const slope=sxy/sxx;
  return {x,y,slope,d:slope*o.Lo,sd:o.Lo/Math.sqrt(sxx)};
}

/* ============================================================
   3) Multilayer forward model and fit (up to 3 layers on a substrate) — SPEC_multi §4.3
   model = {layers:[{name, d}], species:[{name, layer, Iinf, L:[L through layer 0, 1, ...]}]}; layer index = number of layers → substrate.
   I_k(θ) = I_k∞ · [1 − exp(−d_m/(L_k,m cosθ))]^(m<N) · Π_{j<m} exp(−d_j/(L_k,j cosθ))
   ============================================================ */
function multilayerI(model,d,th){
  const N=model.layers.length, c=cosd(th);
  return model.species.map(sp=>{ const m=sp.layer; let v=sp.Iinf; for(let j=0;j<Math.min(m,N);j++) v*=Math.exp(-d[j]/(sp.L[j]*c)); if(m<N) v*=1-Math.exp(-d[m]/(sp.L[m]*c)); return v; });
}
/* profiled residuals: per-angle scale G_a = Σ_k I M/σ² / Σ_k M²/σ² removes any common angle-dependent factor (area, transmission). */
function profiledResid(model,d,th,I,sI,useScale){
  const K=model.species.length, A=th.length, r=new Float64Array(K*A), G=new Float64Array(A);
  for(let a=0;a<A;a++){ const M=multilayerI(model,d,th[a]); let num=0, den=0; for(let k=0;k<K;k++){ const w=1/(sI[k][a]*sI[k][a]); num+=w*I[k][a]*M[k]; den+=w*M[k]*M[k]; }
    G[a]=useScale?(den>0?num/den:0):1; for(let k=0;k<K;k++) r[a*K+k]=(I[k][a]-G[a]*M[k])/sI[k][a]; }
  return {r,G};
}
/* box-constrained Levenberg–Marquardt with forward-difference Jacobian; params p within [lo, hi] */
function lmBox(resFn,p0,lo,hi,o){
  o=o||{}; const M=getMulti(); let p=Float64Array.from(p0); const n=p.length; let r=resFn(p); let f=M.dot(r,r); let mu=1e-3, it=0, J=null;
  const jac=(p,r)=>{ const J=[]; for(let j=0;j<n;j++){ const h=1e-7*Math.max(Math.abs(p[j]),1e-3); const q=Float64Array.from(p); q[j]=Math.min(hi[j],q[j]+h); let hh=q[j]-p[j]; if(hh<=0){ q[j]=p[j]-h; hh=-h; } const rr=resFn(q); J.push(Float64Array.from(rr,(v,i)=>(v-r[i])/hh)); } return J; };
  for(;it<(o.maxIter||300);it++){
    J=jac(p,r); const A=new Float64Array(n*n), g=new Float64Array(n); for(let a=0;a<n;a++){ g[a]=-M.dot(J[a],r); for(let b=0;b<n;b++) A[a*n+b]=M.dot(J[a],J[b]); }
    let accepted=false;
    for(let t=0;t<30;t++){ const B=Float64Array.from(A); for(let a=0;a<n;a++) B[a*n+a]+=mu*Math.max(A[a*n+a],1e-12); const st=M.solveGE(B,g,n); if(!st) { mu*=10; continue; }
      const q=Float64Array.from(p,(v,j)=>Math.min(hi[j],Math.max(lo[j],v+st[j]))); const rq=resFn(q); const fq=M.dot(rq,rq);
      if(fq<f){ const rel=(f-fq)/Math.max(f,1e-300); p=q; r=rq; f=fq; mu=Math.max(mu/3,1e-12); accepted=true; if(rel<1e-14) it=1e9; break; } mu*=4; }
    if(!accepted) break;
  }
  J=jac(p,r); const A=new Float64Array(n*n); for(let a=0;a<n;a++) for(let b=0;b<n;b++) A[a*n+b]=M.dot(J[a],J[b]);
  return {p,r,chi2:f,JTJ:A,iter:it};
}
/* fit layer thicknesses. I[k][a], sI[k][a] (σ; default 1 % of I), o.fix = [bool per layer], o.d0 initial, o.scale (default true: per-angle factor profiled) */
function fitMultilayer(model,th,I,sI,o){
  o=o||{}; const N=model.layers.length; const K=model.species.length, A=th.length;
  sI=sI||I.map(row=>row.map(v=>0.01*Math.abs(v)||1)); const useScale=o.scale!==false;
  const free=[]; for(let j=0;j<N;j++) if(!(o.fix&&o.fix[j])) free.push(j);
  const d0=Float64Array.from(model.layers.map((l,j)=>o.d0?o.d0[j]:(l.d!=null?l.d:1)));
  const full=p=>{ const d=Float64Array.from(d0); free.forEach((j,i)=>d[j]=p[i]); return d; };
  const res=p=>profiledResid(model,full(p),th,I,sI,useScale).r;
  const p0=Float64Array.from(free.map(j=>d0[j]));
  const fit=lmBox(res,p0,free.map(()=>0),free.map(()=>1e3),o);
  const inv=getMulti().invSPD(fit.JTJ,free.length); const dof=Math.max(K*A-free.length-(useScale?A:0),1), phi=fit.chi2/dof;
  const d=full(fit.p); const sd=new Float64Array(N).fill(NaN); if(inv) free.forEach((j,i)=>sd[j]=Math.sqrt(inv[i*free.length+i]));
  const pr=profiledResid(model,d,th,I,sI,useScale);
  return {d,sd,sdScaled:Float64Array.from(sd,v=>v*Math.sqrt(Math.max(phi,1))),chi2:fit.chi2,dof,phi,G:pr.G,iter:fit.iter,
    model:th.map(t=>multilayerI(model,d,t))};
}

/* ============================================================
   4) Relative depth plot (model-independent ordering) — SPEC_multi §4.4
   RDP_k = ln( Σ_{a∈grazing} I_ka / Σ_{a∈normal} I_ka )   (larger = nearer the surface)
   public description: Thermo Fisher application note "Angle Resolved XPS" (doc. D16069); concept: Cumpson, J. Electron Spectrosc.
   Relat. Phenom. 73 (1995) 25. A common angle factor shifts every RDP_k equally → only differences/ordering are meaningful.
   ============================================================ */
function relativeDepth(th,I,sI,o){
  o=o||{}; const A=th.length; const order=th.map((t,a)=>a).sort((a,b)=>th[a]-th[b]);
  const ng=o.nSet||1; const normal=o.normal||order.slice(0,ng), grazing=o.grazing||order.slice(A-ng);
  const rows=I.map((row,k)=>{ let g=0,n=0,vg=0,vn=0; grazing.forEach(a=>{ g+=row[a]; if(sI) vg+=sI[k][a]*sI[k][a]; }); normal.forEach(a=>{ n+=row[a]; if(sI) vn+=sI[k][a]*sI[k][a]; });
    return {k,rdp:Math.log(g/n),s:sI?Math.sqrt(vg/(g*g)+vn/(n*n)):null}; });
  const mean=rows.reduce((s,r)=>s+r.rdp,0)/rows.length; rows.forEach(r=>r.rel=r.rdp-mean);
  return {rows,order:rows.slice().sort((a,b)=>b.rdp-a.rdp).map(r=>r.k),normal,grazing};
}

/* ============================================================
   5) Maximum-entropy depth profile — SPEC_multi §4.5
   Livesey & Smith, J. Electron Spectrosc. Relat. Phenom. 67 (1994) 439; Smith & Livesey, Surf. Interface Anal. 19 (1992) 175.
   Layers j = 0..J−1 with tops z_j (last layer semi-infinite), composition c_kj (Σ_k c_kj = 1 via softmax of u_kj),
   M_ka = Σ_j c_kj W_kaj,  W_kaj = L_k cosθ_a [exp(−z_j/(L_k cosθ_a)) − exp(−z_{j+1}/(L_k cosθ_a))].
   χ² = Σ_ka ((I_ka − G_a M_ka)/σ_ka)² with G_a profiled; S = Σ_kj [c − m − c ln(c/m)];  minimize Q = χ²/2 − α S.
   ============================================================ */
function kernel(th,L,z){
  const J=z.length; return L.map(Lk=>th.map(t=>{ const lc=Lk*cosd(t); const w=new Float64Array(J); for(let j=0;j<J;j++){ const a=Math.exp(-z[j]/lc), b=j<J-1?Math.exp(-z[j+1]/lc):0; w[j]=lc*(a-b); } return w; }));
}
function softmax(u,K,J){ const c=new Float64Array(K*J); for(let j=0;j<J;j++){ let mx=-Infinity; for(let k=0;k<K;k++) mx=Math.max(mx,u[k*J+j]); let s=0; for(let k=0;k<K;k++){ const e=Math.exp(u[k*J+j]-mx); c[k*J+j]=e; s+=e; } for(let k=0;k<K;k++) c[k*J+j]/=s; } return c; }
function meObjective(u,P,alpha){
  const {K,J,A,W,I,sI,m}=P; const c=softmax(u,K,J);
  const M=new Float64Array(K*A); for(let k=0;k<K;k++) for(let a=0;a<A;a++){ let s=0; const w=W[k][a]; for(let j=0;j<J;j++) s+=c[k*J+j]*w[j]; M[k*A+a]=s; }
  const G=new Float64Array(A); let chi=0; const R=new Float64Array(K*A);
  for(let a=0;a<A;a++){ let num=0, den=0; for(let k=0;k<K;k++){ const wt=1/(sI[k][a]*sI[k][a]); num+=wt*I[k][a]*M[k*A+a]; den+=wt*M[k*A+a]*M[k*A+a]; } G[a]=num/den;
    for(let k=0;k<K;k++){ const r=(I[k][a]-G[a]*M[k*A+a])/sI[k][a]; R[k*A+a]=r; chi+=r*r; } }
  let S=0; for(let i=0;i<K*J;i++){ const ci=Math.max(c[i],1e-300); S+=ci-m[i]-ci*Math.log(ci/m[i]); }
  const g=new Float64Array(K*J); // dQ/dc
  for(let k=0;k<K;k++) for(let j=0;j<J;j++){ let s=0; for(let a=0;a<A;a++) s+=-R[k*A+a]*G[a]*W[k][a][j]/sI[k][a]; g[k*J+j]=s+alpha*Math.log(Math.max(c[k*J+j],1e-300)/m[k*J+j]); }
  const gu=new Float64Array(K*J); for(let j=0;j<J;j++){ let t=0; for(let k=0;k<K;k++) t+=c[k*J+j]*g[k*J+j]; for(let k=0;k<K;k++) gu[k*J+j]=c[k*J+j]*(g[k*J+j]-t); }
  return {f:0.5*chi-alpha*S,grad:gu,chi2:chi,S,c,G,M};
}
/* L-BFGS (Nocedal & Wright, Numerical Optimization, alg. 7.4/7.5) with backtracking Armijo line search */
function lbfgs(fg,x0,o){
  o=o||{}; const n=x0.length, mem=o.m||12; let x=Float64Array.from(x0); let F=fg(x); const S=[],Y=[],RHO=[]; let it=0, small=0;
  for(;it<(o.maxIter||5000);it++){
    const g=F.grad; let gn=0; for(let i=0;i<n;i++) gn=Math.max(gn,Math.abs(g[i])); if(gn<(o.gtol||1e-9)) break;
    const q=Float64Array.from(g), al=[]; for(let i=S.length-1;i>=0;i--){ let a=0; for(let t=0;t<n;t++) a+=S[i][t]*q[t]; a*=RHO[i]; al[i]=a; for(let t=0;t<n;t++) q[t]-=a*Y[i][t]; }
    let gam=1; if(S.length){ const l=S.length-1; let sy=0, yy=0; for(let t=0;t<n;t++){ sy+=S[l][t]*Y[l][t]; yy+=Y[l][t]*Y[l][t]; } gam=sy/yy; } else { gam=1/Math.max(gn,1); }
    for(let t=0;t<n;t++) q[t]*=gam;
    for(let i=0;i<S.length;i++){ let b=0; for(let t=0;t<n;t++) b+=Y[i][t]*q[t]; b*=RHO[i]; for(let t=0;t<n;t++) q[t]+=S[i][t]*(al[i]-b); }
    let gd=0; for(let t=0;t<n;t++) gd-=g[t]*q[t]; if(gd>=0){ S.length=0; Y.length=0; RHO.length=0; for(let t=0;t<n;t++) q[t]=g[t]/Math.max(gn,1); gd=0; for(let t=0;t<n;t++) gd-=g[t]*q[t]; }
    let step=1, Fn=null, xn=null;
    for(let ls=0;ls<60;ls++){ xn=Float64Array.from(x,(v,t)=>v-step*q[t]); Fn=fg(xn); if(Fn.f<=F.f+1e-4*step*gd) break; step*=0.5; }
    if(!(Fn.f<=F.f)) break;
    const s=Float64Array.from(xn,(v,t)=>v-x[t]), y=Float64Array.from(Fn.grad,(v,t)=>v-g[t]); let sy=0; for(let t=0;t<n;t++) sy+=s[t]*y[t];
    if(sy>1e-300){ S.push(s); Y.push(y); RHO.push(1/sy); if(S.length>mem){ S.shift(); Y.shift(); RHO.shift(); } }
    const rel=(F.f-Fn.f)/Math.max(Math.abs(F.f),1e-300); x=xn; F=Fn; if(rel<(o.ftol||1e-15)){ if(++small>=5) break; } else small=0;
  }
  return {x,F,iter:it};
}
/* setup: th (deg), L[k] (nm, attenuation length of species k photoelectrons in the film), I[k][a] intensities ÷ sensitivity, sI σ.
   o = {J (layers, default 20), zmax (default 3·max L·cos θ_min), m (default model: 'normal' = apparent composition at smallest θ, or array K), alpha|'auto'} */
function maxentSetup(th,L,I,sI,o){
  o=o||{}; const K=L.length, A=th.length, J=o.J||20; const cmin=Math.min(...th.map(cosd)), cmax=Math.max(...th.map(cosd));
  const zmax=o.zmax||3*Math.max(...L)*cmax; const z=Float64Array.from({length:J},(_,j)=>j*zmax/(J-1));
  sI=sI||I.map(r=>r.map(v=>0.01*Math.abs(v)));
  let m0;
  if(Array.isArray(o.m)) m0=o.m;
  else { const a0=th.indexOf(Math.min(...th)); const tot=I.reduce((s,r)=>s+r[a0],0); m0=I.map(r=>r[a0]/tot); }
  const m=new Float64Array(K*J); for(let k=0;k<K;k++) for(let j=0;j<J;j++) m[k*J+j]=Math.max(m0[k],1e-6);
  for(let j=0;j<J;j++){ let s=0; for(let k=0;k<K;k++) s+=m[k*J+j]; for(let k=0;k<K;k++) m[k*J+j]/=s; }
  return {K,J,A,th,L,I,sI,z,W:kernel(th,L,z),m,cmin};
}
function maxentSolve(P,alpha,u0){
  const u=u0?Float64Array.from(u0):Float64Array.from(P.m,v=>Math.log(v));
  const r=lbfgs(x=>meObjective(x,P,alpha),u,{maxIter:20000,gtol:1e-10});
  return {alpha,u:r.x,c:r.F.c,chi2:r.F.chi2,S:r.F.S,Q:r.F.f,G:r.F.G,M:r.F.M,iter:r.iter};
}
/* α by the discrepancy ('historic MaxEnt') criterion χ²(α) = N, N = K·A − A (data minus profiled per-angle scales).
   Continuation (SPEC_multi §4.5): α_k = α_max·10^(−k/2) from α_max = 1e6 down to α_min = 1e−4, each solve warm-started from the previous
   (large α = unique, near the default model); the first α_k with χ² ≤ N brackets the root with α_{k−1}; then 40 bisection steps in ln α,
   always warm-started from the larger-α end of the bracket (deterministic path). Small-α solutions are ill-posed — never used as starts. */
function maxentProfile(th,L,I,sI,o){
  o=o||{}; const P=maxentSetup(th,L,I,sI,o); const N=o.target||(P.K*P.A-P.A); const curve=[]; let sol;
  if(typeof o.alpha==="number") sol=maxentSolve(P,o.alpha);
  else {
    const aMax=o.alphaMax||1e6, aMin=o.alphaMin||1e-4; let prev=null, cur=null;
    for(let k=0;;k++){ const a=aMax*Math.pow(10,-k/2); if(a<aMin*0.999) break; cur=maxentSolve(P,a,prev?prev.u:null); curve.push({alpha:a,chi2:cur.chi2,S:cur.S}); if(cur.chi2<=N) break; prev=cur; }
    if(!prev&&cur.chi2<=N) sol=cur;                       // the default-like solution already fits within noise
    else if(cur.chi2>N){ sol=cur; sol.notReached=true; }  // noise level never reached (model error or underestimated σ)
    else { let hi=Math.log(prev.alpha), lo=Math.log(cur.alpha), sHi=prev; sol=cur;
      for(let i=0;i<40;i++){ const mid=0.5*(lo+hi); const s=maxentSolve(P,Math.exp(mid),sHi.u); curve.push({alpha:s.alpha,chi2:s.chi2,S:s.S}); if(s.chi2>N){ hi=mid; sHi=s; } else { lo=mid; sol=s; } if(hi-lo<1e-4) break; } }
  }
  const K=P.K, J=P.J; const c=Array.from({length:K},(_,k)=>Float64Array.from({length:J},(_,j)=>sol.c[k*J+j]));
  const dz=P.z[1]-P.z[0];
  return {z:P.z,dz,c,alpha:sol.alpha,chi2:sol.chi2,N,S:sol.S,G:sol.G,iter:sol.iter,curve,notReached:!!sol.notReached,setup:P,
    model:Array.from({length:K},(_,k)=>Float64Array.from({length:P.A},(_,a)=>sol.G[a]*sol.M[k*P.A+a])),
    warn:["ARXPS carries only a few independent depth parameters (Cumpson 1995): the MaxEnt profile is the smoothest one consistent with the data, not a unique profile",
          "below ≈ "+(3*Math.max(...L)*P.cmin).toFixed(2)+"–"+(3*Math.max(...L)).toFixed(2)+" nm the profile reverts to the default model (no information)",
          "straight-line attenuation (no elastic scattering); angles > 65° bias the result"]};
}

const XPS_ARXPS={version:"1.0",atomSize,seahDench,ealFromImfp,overlayerR,overlayerdR,overlayerThickness,fitOverlayer,hillPlot,
  multilayerI,profiledResid,lmBox,fitMultilayer,relativeDepth,kernel,softmax,meObjective,lbfgs,maxentSetup,maxentSolve,maxentProfile};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_ARXPS; else root.XPS_ARXPS=XPS_ARXPS;
})(typeof self!=="undefined"?self:this);
