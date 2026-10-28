/* ============================================================
   INNO Analytics — XPS series engine: depth profiles + multivariate (PCA / TFA / MCR / reference NNLS)
   DOM-free: browser main thread, Web Worker (load after xps_core.js) and node.
   Equations, citations and validation: modules/xps/SPEC_multi.md (this file implements it; if they disagree, fix SPEC first).
   Data orientation everywhere: a data matrix D is r energies x c spectra, stored as an ARRAY OF COLUMNS
   (c Float64Array of length r) — one column per spectrum (etch level, angle, pixel).
   ============================================================ */
(function(root){
"use strict";
function getCore(){ if(root&&root.XPS_CORE) return root.XPS_CORE; if(typeof require==="function"){ try{ return require("./xps_core.js"); }catch(e){} } return null; }

/* ============================================================
   0) Small dense linear algebra (no dependencies)
   ============================================================ */
const EPS=2.220446049250313e-16;
function dot(a,b){ let s=0; for(let i=0;i<a.length;i++) s+=a[i]*b[i]; return s; }
function norm2(a){ return Math.sqrt(dot(a,a)); }
/* lgamma (Lanczos g=7, n=9) — for the F distribution only */
const LGC=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
function lgamma(x){ if(x<0.5) return Math.log(Math.PI/Math.abs(Math.sin(Math.PI*x)))-lgamma(1-x); x-=1; let a=LGC[0]; const t=x+7.5; for(let i=1;i<9;i++) a+=LGC[i]/(x+i); return 0.5*Math.log(2*Math.PI)+(x+0.5)*Math.log(t)-t+Math.log(a); }
/* regularized incomplete beta I_x(a,b) — continued fraction (modified Lentz) */
function betacf(a,b,x){ const FP=1e-300; let qab=a+b, qap=a+1, qam=a-1, c=1, d=1-qab*x/qap; if(Math.abs(d)<FP) d=FP; d=1/d; let h=d;
  for(let m=1;m<=300;m++){ const m2=2*m; let aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=1+aa*d; if(Math.abs(d)<FP) d=FP; c=1+aa/c; if(Math.abs(c)<FP) c=FP; d=1/d; h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=1+aa*d; if(Math.abs(d)<FP) d=FP; c=1+aa/c; if(Math.abs(c)<FP) c=FP; d=1/d; const del=d*c; h*=del; if(Math.abs(del-1)<1e-15) break; }
  return h; }
function betai(a,b,x){ if(x<=0) return 0; if(x>=1) return 1; const bt=Math.exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*Math.log(x)+b*Math.log(1-x));
  return x<(a+1)/(a+b+2)?bt*betacf(a,b,x)/a:1-bt*betacf(b,a,1-x)/b; }
/* P(F > f) for F(d1,d2) */
function fSurvival(f,d1,d2){ if(!(f>0)) return 1; return betai(d2/2,d1/2,d2/(d2+d1*f)); }

/* Cholesky solve of SPD system (row-major n x n). Returns null if not SPD. */
function cholesky(A,n){ const L=new Float64Array(n*n); for(let i=0;i<n;i++) for(let j=0;j<=i;j++){ let s=A[i*n+j]; for(let k=0;k<j;k++) s-=L[i*n+k]*L[j*n+k];
  if(i===j){ if(!(s>0)) return null; L[i*n+i]=Math.sqrt(s); } else L[i*n+j]=s/L[j*n+j]; } return L; }
function cholSolveL(L,b,n){ const y=new Float64Array(n); for(let i=0;i<n;i++){ let s=b[i]; for(let k=0;k<i;k++) s-=L[i*n+k]*y[k]; y[i]=s/L[i*n+i]; }
  for(let i=n-1;i>=0;i--){ let s=y[i]; for(let k=i+1;k<n;k++) s-=L[k*n+i]*y[k]; y[i]=s/L[i*n+i]; } return y; }
function solveSPD(A,b,n){ const L=cholesky(A,n); return L?cholSolveL(L,b,n):null; }
function invSPD(A,n){ const L=cholesky(A,n); if(!L) return null; const inv=new Float64Array(n*n); for(let j=0;j<n;j++){ const e=new Float64Array(n); e[j]=1; const x=cholSolveL(L,e,n); for(let i=0;i<n;i++) inv[i*n+j]=x[i]; } return inv; }
/* general solve (Gaussian elimination, partial pivot) */
function solveGE(A0,b0,n){ const A=Float64Array.from(A0), b=Float64Array.from(b0);
  for(let k=0;k<n;k++){ let p=k, mx=Math.abs(A[k*n+k]); for(let i=k+1;i<n;i++) if(Math.abs(A[i*n+k])>mx){ mx=Math.abs(A[i*n+k]); p=i; } if(mx<1e-300) return null;
    if(p!==k){ for(let j=0;j<n;j++){ const t=A[k*n+j]; A[k*n+j]=A[p*n+j]; A[p*n+j]=t; } const t=b[k]; b[k]=b[p]; b[p]=t; }
    for(let i=k+1;i<n;i++){ const f=A[i*n+k]/A[k*n+k]; if(!f) continue; for(let j=k;j<n;j++) A[i*n+j]-=f*A[k*n+j]; b[i]-=f*b[k]; } }
  const x=new Float64Array(n); for(let i=n-1;i>=0;i--){ let s=b[i]; for(let j=i+1;j<n;j++) s-=A[i*n+j]*x[j]; x[i]=s/A[i*n+i]; } return x; }
function detGram(cols){ const k=cols.length, G=new Float64Array(k*k); for(let a=0;a<k;a++) for(let b=0;b<=a;b++){ const v=dot(cols[a],cols[b]); G[a*k+b]=v; G[b*k+a]=v; }
  // determinant via elimination
  let det=1; const A=G; for(let c=0;c<k;c++){ let p=c, mx=Math.abs(A[c*k+c]); for(let i=c+1;i<k;i++) if(Math.abs(A[i*k+c])>mx){ mx=Math.abs(A[i*k+c]); p=i; } if(mx===0) return 0;
    if(p!==c){ for(let j=0;j<k;j++){ const t=A[c*k+j]; A[c*k+j]=A[p*k+j]; A[p*k+j]=t; } det=-det; } det*=A[c*k+c]; for(let i=c+1;i<k;i++){ const f=A[i*k+c]/A[c*k+c]; for(let j=c;j<k;j++) A[i*k+j]-=f*A[c*k+j]; } } return det; }

/* Symmetric eigen-decomposition, cyclic Jacobi. A row-major n x n (not modified). Returns {w (desc), V: array of n eigenvectors} */
function symEig(A0,n){
  const A=Float64Array.from(A0), V=new Float64Array(n*n); for(let i=0;i<n;i++) V[i*n+i]=1;
  for(let sweep=0;sweep<100;sweep++){
    let off=0, tot=0; for(let i=0;i<n;i++) for(let j=0;j<n;j++){ const v=A[i*n+j]*A[i*n+j]; tot+=v; if(i!==j) off+=v; }
    if(off<=1e-30*tot||off===0) break;
    for(let p=0;p<n-1;p++) for(let q=p+1;q<n;q++){
      const apq=A[p*n+q]; if(Math.abs(apq)<1e-300) continue;
      const app=A[p*n+p], aqq=A[q*n+q]; const th=(aqq-app)/(2*apq); const t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)); const c=1/Math.sqrt(t*t+1), s=t*c;
      for(let k=0;k<n;k++){ const akp=A[k*n+p], akq=A[k*n+q]; A[k*n+p]=c*akp-s*akq; A[k*n+q]=s*akp+c*akq; }
      for(let k=0;k<n;k++){ const apk=A[p*n+k], aqk=A[q*n+k]; A[p*n+k]=c*apk-s*aqk; A[q*n+k]=s*apk+c*aqk; }
      for(let k=0;k<n;k++){ const vkp=V[k*n+p], vkq=V[k*n+q]; V[k*n+p]=c*vkp-s*vkq; V[k*n+q]=s*vkp+c*vkq; }
    }
  }
  const idx=Array.from({length:n},(_,i)=>i).sort((a,b)=>A[b*n+b]-A[a*n+a]);
  return {w:Float64Array.from(idx,i=>A[i*n+i]),V:idx.map(i=>Float64Array.from({length:n},(_,k)=>V[k*n+i]))};
}

/* sign convention shared with the Python reference: the element of largest |value| of each left vector is positive */
function fixSign(u,v){ let im=0; for(let i=1;i<u.length;i++) if(Math.abs(u[i])>Math.abs(u[im])) im=i; if(u[im]<0){ for(let i=0;i<u.length;i++) u[i]=-u[i]; if(v) for(let i=0;i<v.length;i++) v[i]=-v[i]; } }

/* Thin SVD of D (array of c columns, each length r) by one-sided Jacobi (Hestenes) on the side with fewer vectors.
   Returns {U:[k vectors len r], s:Float64Array(k), V:[k vectors len c]}, k = min(r,c), s descending. */
function svdCols(cols){
  const c=cols.length, r=cols[0].length;
  if(c>r){ // work on D^T: its columns are the r rows of D
    const rows=Array.from({length:r},(_,i)=>Float64Array.from({length:c},(_,j)=>cols[j][i]));
    const t=svdCols(rows); t.U.forEach((u,k)=>fixSign(t.V[k],u));   // D = V' S U'^T
    const out={U:t.V,s:t.s,V:t.U}; out.U.forEach((u,k)=>fixSign(u,out.V[k])); return out;
  }
  const A=cols.map(v=>Float64Array.from(v)), n=c; const V=Array.from({length:n},(_,j)=>{ const e=new Float64Array(n); e[j]=1; return e; });
  for(let sweep=0;sweep<80;sweep++){
    let rot=0;
    for(let p=0;p<n-1;p++) for(let q=p+1;q<n;q++){
      const ap=A[p], aq=A[q]; let al=0,be=0,ga=0; for(let i=0;i<r;i++){ al+=ap[i]*ap[i]; be+=aq[i]*aq[i]; ga+=ap[i]*aq[i]; }
      if(Math.abs(ga)<=1e-15*Math.sqrt(al*be)||ga===0) continue; rot++;
      const ze=(be-al)/(2*ga); const t=(ze>=0?1:-1)/(Math.abs(ze)+Math.sqrt(1+ze*ze)); const cs=1/Math.sqrt(1+t*t), sn=cs*t;
      for(let i=0;i<r;i++){ const x=ap[i], y=aq[i]; ap[i]=cs*x-sn*y; aq[i]=sn*x+cs*y; }
      const vp=V[p], vq=V[q]; for(let i=0;i<n;i++){ const x=vp[i], y=vq[i]; vp[i]=cs*x-sn*y; vq[i]=sn*x+cs*y; }
    }
    if(!rot) break;
  }
  const sv=A.map(norm2); const idx=sv.map((_,i)=>i).sort((a,b)=>sv[b]-sv[a]);
  const U=[], s=new Float64Array(n), Vo=[];
  idx.forEach((j,k)=>{ s[k]=sv[j]; const u=Float64Array.from(A[j],x=>sv[j]>0?x/sv[j]:0); const v=Float64Array.from(V[j]); fixSign(u,v); U.push(u); Vo.push(v); });
  return {U,s,V:Vo};
}
/* Large-matrix PCA via the r x r Gram matrix (r = number of energies). Used for image stacks (r ~ 100, c ~ 65536).
   Only the leading k vectors are returned; small singular values lose relative accuracy (squared condition) — documented. */
function svdGram(cols,k){
  const c=cols.length, r=cols[0].length; const G=new Float64Array(r*r);
  for(let j=0;j<c;j++){ const x=cols[j]; for(let a=0;a<r;a++){ const xa=x[a]; if(!xa) continue; const o=a*r; for(let b=0;b<=a;b++) G[o+b]+=xa*x[b]; } }
  for(let a=0;a<r;a++) for(let b=0;b<a;b++) G[b*r+a]=G[a*r+b];
  const E=symEig(G,r); k=Math.min(k||r,r,c);
  const s=new Float64Array(Math.min(r,c)); for(let i=0;i<s.length;i++) s[i]=Math.sqrt(Math.max(E.w[i],0));
  const U=[], V=[];
  for(let i=0;i<k;i++){ const u=Float64Array.from(E.V[i]); const v=new Float64Array(c); for(let j=0;j<c;j++) v[j]=s[i]>0?dot(cols[j],u)/s[i]:0; fixSign(u,v); U.push(u); V.push(v); }
  return {U,s,V,gram:true};
}

/* ============================================================
   1) NNLS (Lawson & Hanson 1974, ch. 23) in Gram form: min ||A x - b||, x >= 0, given G = AᵀA, h = Aᵀb.
   free[]: indices that are unconstrained in sign (handled by splitting into ± columns by the caller is NOT needed:
   here they are kept permanently in the passive set and may take any sign).
   ============================================================ */
function nnlsGram(G,h,n,free){
  const x=new Float64Array(n), P=new Uint8Array(n), F=new Uint8Array(n); if(free) free.forEach(j=>{ F[j]=1; });
  let gmax=0; for(let i=0;i<n*n;i++) gmax=Math.max(gmax,Math.abs(G[i])); let hmax=0; for(let i=0;i<n;i++) hmax=Math.max(hmax,Math.abs(h[i]));
  const tol=10*EPS*Math.max(gmax,1e-300)*n*Math.max(1,hmax/Math.max(gmax,1e-300));
  const solveP=idx=>{ const k=idx.length, A=new Float64Array(k*k), b=new Float64Array(k); for(let a=0;a<k;a++){ b[a]=h[idx[a]]; for(let c=0;c<k;c++) A[a*k+c]=G[idx[a]*n+idx[c]]; }
    return solveSPD(A,b,k)||solveGE(A,b,k)||new Float64Array(k); };
  // unconstrained (free) variables first
  if(F.some(v=>v)){ for(let j=0;j<n;j++) if(F[j]) P[j]=1; const idx=[]; for(let j=0;j<n;j++) if(P[j]) idx.push(j); const z=solveP(idx); idx.forEach((j,a)=>x[j]=z[a]); }
  const grad=()=>{ const w=new Float64Array(n); for(let i=0;i<n;i++){ let s=h[i]; for(let j=0;j<n;j++) if(x[j]) s-=G[i*n+j]*x[j]; w[i]=s; } return w; };
  for(let outer=0;outer<3*n+30;outer++){
    const w=grad(); let jm=-1, wm=tol; for(let j=0;j<n;j++) if(!P[j]&&w[j]>wm){ wm=w[j]; jm=j; }
    if(jm<0) break; P[jm]=1;
    for(let inner=0;inner<3*n+30;inner++){
      const idx=[]; for(let j=0;j<n;j++) if(P[j]) idx.push(j);
      const z=solveP(idx);
      let okAll=true; idx.forEach((j,a)=>{ if(!F[j]&&z[a]<=0) okAll=false; });
      if(okAll){ idx.forEach((j,a)=>x[j]=z[a]); break; }
      let al=Infinity; idx.forEach((j,a)=>{ if(!F[j]&&z[a]<=0){ const t=x[j]/(x[j]-z[a]); if(t<al) al=t; } }); if(!isFinite(al)) al=0;
      idx.forEach((j,a)=>{ x[j]+=al*(z[a]-x[j]); });
      idx.forEach(j=>{ if(!F[j]&&x[j]<=1e-14*Math.max(1,Math.abs(x[j]))){ x[j]=0; P[j]=0; } });
    }
  }
  return x;
}
/* convenience: NNLS from columns (A = cols as columns), with optional row weights wts (multiplies squared residuals) */
function nnlsCols(cols,b,wts,free){
  const n=cols.length, m=b.length; const G=new Float64Array(n*n), h=new Float64Array(n);
  for(let a=0;a<n;a++){ const ca=cols[a]; let s=0; for(let i=0;i<m;i++) s+=(wts?wts[i]:1)*ca[i]*b[i]; h[a]=s;
    for(let c=0;c<=a;c++){ const cc=cols[c]; let t=0; for(let i=0;i<m;i++) t+=(wts?wts[i]:1)*ca[i]*cc[i]; G[a*n+c]=t; G[c*n+a]=t; } }
  return {x:nnlsGram(G,h,n,free),G,h};
}

/* ============================================================
   2) PCA / abstract factor analysis (SPEC_multi §2)
   ============================================================ */
/* preprocessing: scale = "none" | "poisson" (Keenan & Kotula 2004) | "area" (each spectrum / its sum); center = "none" | "mean" */
function preprocess(cols,opt){
  opt=opt||{}; const c=cols.length, r=cols[0].length; const scale=opt.scale||"none", center=opt.center||"none";
  let X=cols.map(v=>Float64Array.from(v)); let rowW=null, colW=null, mean=null;
  if(scale==="area"){ colW=new Float64Array(c); X.forEach((v,j)=>{ let s=0; for(let i=0;i<r;i++) s+=v[i]; colW[j]=s>0?1/s:1; for(let i=0;i<r;i++) v[i]*=colW[j]; }); }
  else if(scale==="poisson"){ // D~ = G^-1/2 D H^-1/2, G = mean spectrum (row means), H = mean of each spectrum (column means)
    const g=new Float64Array(r), hh=new Float64Array(c); X.forEach((v,j)=>{ let s=0; for(let i=0;i<r;i++){ g[i]+=v[i]; s+=v[i]; } hh[j]=s/r; }); for(let i=0;i<r;i++) g[i]/=c;
    const gpos=Math.max(...g)*1e-6||1e-300, hpos=Math.max(...hh)*1e-6||1e-300;
    rowW=Float64Array.from(g,v=>1/Math.sqrt(Math.max(v,gpos))); colW=Float64Array.from(hh,v=>1/Math.sqrt(Math.max(v,hpos)));
    X.forEach((v,j)=>{ for(let i=0;i<r;i++) v[i]*=rowW[i]*colW[j]; });
  } else if(scale!=="none") throw new Error("unknown scale "+scale);
  if(center==="mean"){ mean=new Float64Array(r); X.forEach(v=>{ for(let i=0;i<r;i++) mean[i]+=v[i]/c; }); X.forEach(v=>{ for(let i=0;i<r;i++) v[i]-=mean[i]; }); }
  return {X,rowW,colW,mean,scale,center};
}
function undo(pp,v,j){ // processed-space column j -> original units
  const r=v.length, o=new Float64Array(r); for(let i=0;i<r;i++){ let x=v[i]+(pp.mean?pp.mean[i]:0); if(pp.rowW) x/=pp.rowW[i]; if(pp.colW) x/=pp.colW[j]; o[i]=x; } return o; }

/* Malinowski statistics on eigenvalues λ (= s²) of an r x c matrix:
   RE(n) = sqrt( Σ_{j>n} λ_j / (R (C − n)) ), IND(n) = RE(n)/(C − n)², R = max(r,c), C = min(r,c)  [Malinowski, Anal. Chem. 49 (1977) 612]
   F(1, C−n) = [Σ_{j=n+1}^{C} (R−j+1)(C−j+1)] / [(R−n+1)(C−n+1)] · λ_n / Σ_{j=n+1}^{C} λ_j  [Malinowski, J. Chemometrics 3 (1988) 49] */
function malinowski(lam,r,c,alpha){
  const R=Math.max(r,c), C=Math.min(r,c); alpha=alpha||0.05; const out=[];
  for(let n=1;n<C;n++){
    let tail=0; for(let j=n;j<C;j++) tail+=lam[j];          // λ index j (0-based) = factor j+1
    const RE=Math.sqrt(tail/(R*(C-n))), IND=RE/((C-n)*(C-n));
    let num=0; for(let j=n+1;j<=C;j++) num+=(R-j+1)*(C-j+1);
    const F=tail>0?num/((R-n+1)*(C-n+1))*lam[n-1]/tail:Infinity; const SL=fSurvival(F,1,C-n);
    const REV=lam[n-1]/((R-n+1)*(C-n+1));
    out.push({n,RE,IND,F,SL,REV});
  }
  let nIND=1, best=Infinity; out.forEach(o=>{ if(o.IND<best){ best=o.IND; nIND=o.n; } });
  let nF=0; for(const o of out){ if(o.SL<alpha) nF=o.n; else break; }
  return {table:out,nIND,nF,alpha};
}
/* pca(cols, {scale, center, nKeep, method:"auto"|"jacobi"|"gram", alpha}) */
function pca(cols,opt){
  opt=opt||{}; const c=cols.length, r=cols[0].length; if(c<2) throw new Error("PCA needs at least 2 spectra");
  const pp=preprocess(cols,opt);
  const method=opt.method&&opt.method!=="auto"?opt.method:((Math.min(r,c)<=30||r*c<=2e5)?"jacobi":"gram");
  const S=method==="jacobi"?svdCols(pp.X):svdGram(pp.X,opt.nKeep||Math.min(r,c,20));
  const lam=Float64Array.from(S.s,v=>v*v); const tot=lam.reduce((a,b)=>a+b,0);
  const varPct=Float64Array.from(lam,v=>tot>0?100*v/tot:0);
  const mal=malinowski(lam,r,pp.mean?c-1:c,opt.alpha);   // mean-centring removes one dimension
  return {r,c,pp,method,U:S.U,s:S.s,V:S.V,eig:lam,varPct,malinowski:mal,warn:pp.mean?["Malinowski RE/IND/F were derived for non-centred data; with mean-centring c−1 is used (approximation)"]:[]};
}
/* reconstruction with n factors, in original units */
function reconstruct(P,n){
  const r=P.r, c=P.c, out=[]; n=Math.min(n,P.U.length);
  for(let j=0;j<c;j++){ const v=new Float64Array(r); for(let k=0;k<n;k++){ const f=P.s[k]*P.V[k][j]; const u=P.U[k]; for(let i=0;i<r;i++) v[i]+=f*u[i]; } out.push(undo(P.pp,v,j)); }
  return out;
}
/* Weighted low-rank approximation for Poisson data (maximum-likelihood PCA, Wentzell et al., J. Chemometrics 11 (1997) 339):
   min Σ_ij w_ij (D_ij − (A Bᵀ)_ij)², w_ij = cfac / max(D̂_ij, 1/cfac) (Pearson weights from the current model, refreshed every sweep).
   Alternating weighted least squares; start = rank-n SVD of the Poisson-scaled data (SPEC_multi §2.3). Returns D̂ (columns) and χ². */
function wlra(cols,n,o){
  o=o||{}; const c=cols.length, r=cols[0].length, cf=o.cfac>0?o.cfac:1;
  if(n===0){ let chi=0; cols.forEach(d=>{ for(let i=0;i<r;i++){ const v=Math.max(d[i],1/cf)/cf; chi+=d[i]*d[i]/v; } }); return {rec:cols.map(()=>new Float64Array(r)),chi2:chi,iter:0}; }
  const P0=o.pca0||pca(cols,{scale:"poisson"}); const R0=reconstruct(P0,n);
  // factor start: A = G^½·U_n S_n, B = H^½·V_n (the rank-n Poisson-scaled SVD mapped back; A Bᵀ = R0). Every ALS product depends on A only
  // through its column span, so any factorization of R0 gives the same iterates (the Python reference factorizes R0 by SVD).
  const rw=P0.pp.rowW, cw=P0.pp.colW; let A=Array.from({length:n},(_,k)=>Float64Array.from(P0.U[k],(v,i)=>v*P0.s[k]/(rw?rw[i]:1))); let B=Array.from({length:n},(_,k)=>Float64Array.from(P0.V[k],(v,j)=>v/(cw?cw[j]:1)));
  let rec=R0, chiOld=Infinity, chi=0, it=0; const G=new Float64Array(n*n), h=new Float64Array(n);
  const W=new Float64Array(r*c);   // W[j*r+i]
  for(;it<(o.maxIter||300);it++){
    for(let j=0;j<c;j++){ const m=rec[j], o2=j*r; for(let i=0;i<r;i++){ const v=m[i]; W[o2+i]=cf/(v>1/cf?v:1/cf); } }
    // B step: each column j
    for(let j=0;j<c;j++){ G.fill(0); h.fill(0); const d=cols[j], o2=j*r; for(let i=0;i<r;i++){ const wi=W[o2+i]; for(let a=0;a<n;a++){ const ta=wi*A[a][i]; h[a]+=ta*d[i]; for(let b=0;b<=a;b++) G[a*n+b]+=ta*A[b][i]; } }
      for(let a=0;a<n;a++) for(let b=0;b<a;b++) G[b*n+a]=G[a*n+b]; const x=solveSPD(G,h,n)||solveGE(G,h,n); if(x) for(let a=0;a<n;a++) B[a][j]=x[a]; }
    // A step: each row i
    for(let i=0;i<r;i++){ G.fill(0); h.fill(0); for(let j=0;j<c;j++){ const wi=W[j*r+i], d=cols[j][i]; for(let a=0;a<n;a++){ const ta=wi*B[a][j]; h[a]+=ta*d; for(let b=0;b<=a;b++) G[a*n+b]+=ta*B[b][j]; } }
      for(let a=0;a<n;a++) for(let b=0;b<a;b++) G[b*n+a]=G[a*n+b]; const x=solveSPD(G,h,n)||solveGE(G,h,n); if(x) for(let a=0;a<n;a++) A[a][i]=x[a]; }
    rec=new Array(c); for(let j=0;j<c;j++){ const v=new Float64Array(r); for(let a=0;a<n;a++){ const ba=B[a][j], Aa=A[a]; for(let i=0;i<r;i++) v[i]+=Aa[i]*ba; } rec[j]=v; }
    chi=0; for(let j=0;j<c;j++){ const d=cols[j], m=rec[j]; for(let i=0;i<r;i++){ const e=d[i]-m[i], v=m[i]; chi+=cf/(v>1/cf?v:1/cf)*e*e; } }
    if(Math.abs(chiOld-chi)<=(o.tol||1e-10)*chi) break; chiOld=chi;
  }
  return {rec,chi2:chi,iter:it+1};
}
/* Noise-limited number of factors (SPEC_multi §2.3). cols in CPS (or counts), cfac = counts per unit (dwell×scans; 1 for counts).
   χ²_n = Σ (D − D̂_n)² / v, v_ij = max(D̂_n,ij, 1/cfac)/cfac (Pearson, Poisson variance from the reconstruction; n = 0 uses the data);  dof_n = (r − n)(c' − n), c' = c (−1 if centred).
   n_noise = smallest n with χ²_n ≤ dof_n + z·sqrt(2 dof_n) (z = 2). */
function noiseFactors(cols,P,cfac,o){
  o=o||{}; const z=o.z!=null?o.z:2; const r=P.r, c=P.c, cc=P.pp.mean&&o.method==="pca"?c-1:c; cfac=cfac>0?cfac:1; const method=o.method||"mlpca";
  const nMax=Math.min(o.nMax||Math.min(r,cc)-1,P.U.length); const tab=[]; let nNoise=null;
  const P0=method==="mlpca"?(P.pp.scale==="poisson"&&!P.pp.mean&&P.U.length>=Math.min(nMax,Math.min(r,c))?P:pca(cols,{scale:"poisson",nKeep:nMax+1})):null;
  for(let n=0;n<=nMax;n++){
    if(method==="mlpca"){ const W=wlra(cols,n,{cfac,pca0:P0}); const dof=(r-n)*(c-n), lim=dof+z*Math.sqrt(2*dof);
      tab.push({n,chi2:W.chi2,dof,ratio:W.chi2/dof,limit:lim,within:W.chi2<=lim,iter:W.iter}); if(nNoise==null&&W.chi2<=lim){ nNoise=n; if(!o.full) break; } continue; }
    const R=n>0?reconstruct(P,n):cols.map((v,j)=>P.pp.mean?undo(P.pp,new Float64Array(r),j):new Float64Array(r));
    let chi=0; for(let j=0;j<c;j++){ const d=cols[j], m=R[j]; for(let i=0;i<r;i++){ const v=Math.max(n>0?m[i]:d[i],1/cfac)/cfac; const e=d[i]-m[i]; chi+=e*e/v; } }
    const dof=(r-n)*(cc-n); const lim=dof+z*Math.sqrt(2*dof);
    tab.push({n,chi2:chi,dof,ratio:chi/dof,limit:lim,within:chi<=lim});
    if(nNoise==null&&chi<=lim) nNoise=n;
  }
  return {table:tab,nNoise,z,method};
}
/* Target testing (Malinowski, Anal. Chim. Acta 103 (1978) 339; "Factor Analysis in Chemistry", 3rd ed., 2002, ch. 4).
   Target x (length r, original units) is row-scaled like the data; t = S⁻¹Uᵀx (transformation vector), x̂ = U Uᵀ x,
   AET = sqrt(Σ(x − x̂)²/r), REP = RE(n)·‖t‖, RET = sqrt(max(AET² − REP², 0)), SPOIL = RET/REP. */
function targetTest(P,n,x){
  if(P.pp.mean) throw new Error("target testing needs a non-centred PCA (center:'none')");
  const r=P.r; const xs=Float64Array.from(x,(v,i)=>P.pp.rowW?v*P.pp.rowW[i]:v);
  const t=new Float64Array(n), xh=new Float64Array(r);
  for(let k=0;k<n;k++){ const p=dot(P.U[k],xs); t[k]=p/P.s[k]; const u=P.U[k]; for(let i=0;i<r;i++) xh[i]+=p*u[i]; }
  let ss=0; for(let i=0;i<r;i++){ const e=xs[i]-xh[i]; ss+=e*e; }
  const AET=Math.sqrt(ss/r); const RE=n<Math.min(P.r,P.c)?P.malinowski.table[n-1].RE:0; const REP=RE*norm2(t);
  const RET=Math.sqrt(Math.max(AET*AET-REP*REP,0)); const SPOIL=REP>0?RET/REP:Infinity;
  const verdict=SPOIL<3?"acceptable":SPOIL<=6?"moderately acceptable":"rejected";
  const xhat=Float64Array.from(xh,(v,i)=>P.pp.rowW?v/P.pp.rowW[i]:v);
  // cosine similarity of target with its projection (scale-free)
  const cos=dot(xs,xh)/(norm2(xs)*norm2(xh)||1);
  return {t,xhat,AET,REP,RET,SPOIL,verdict,cos};
}

/* Orthogonal projection approach (Cuesta Sánchez et al., Anal. Chem. 68 (1996) 79): picks the n most mutually dissimilar spectra */
function opa(cols,n){
  const c=cols.length, r=cols[0].length; const Y=cols.map(v=>{ const s=norm2(v)||1; return Float64Array.from(v,x=>x/s); });
  const mean=new Float64Array(r); Y.forEach(v=>{ for(let i=0;i<r;i++) mean[i]+=v[i]/c; }); { const s=norm2(mean)||1; for(let i=0;i<r;i++) mean[i]/=s; }
  const picks=[];
  for(let k=0;k<n;k++){ const refs=k===0?[mean]:picks.map(j=>Y[j]); let best=-1, bj=-1;
    for(let j=0;j<c;j++){ if(picks.includes(j)) continue; const d=detGram(refs.concat([Y[j]])); if(d>best){ best=d; bj=j; } }
    picks.push(bj); }
  return picks;
}
/* SIMPLISMA purest variables in the ENERGY domain (Windig & Guilment, Anal. Chem. 63 (1991) 1425):
   for channel i over the c spectra: μ_i, σ_i (population), offset α = f·max μ; λ_i = sqrt(σ_i² + (μ_i + α)²); y_ij = D_ij/λ_i;
   weight w_i = det(Yᵀ Y / c) of [selected y's, y_i] (w = 1 for the first); purity p_i = w_i σ_i/(μ_i + α). Returns n channel indices. */
function simplisma(cols,n,f){
  f=f!=null?f:0.03; const c=cols.length, r=cols[0].length; const mu=new Float64Array(r), sd=new Float64Array(r);
  for(let i=0;i<r;i++){ let s=0; for(let j=0;j<c;j++) s+=cols[j][i]; mu[i]=s/c; let v=0; for(let j=0;j<c;j++){ const d=cols[j][i]-mu[i]; v+=d*d; } sd[i]=Math.sqrt(v/c); }
  const al=f*Math.max(...mu); const Y=Array.from({length:r},(_,i)=>{ const l=Math.sqrt(sd[i]*sd[i]+(mu[i]+al)*(mu[i]+al)); return Float64Array.from({length:c},(_,j)=>cols[j][i]/l/Math.sqrt(c)); });
  const picks=[];
  for(let k=0;k<n;k++){ let best=-Infinity, bi=-1;
    for(let i=0;i<r;i++){ if(picks.includes(i)) continue; const w=k===0?1:detGram(picks.map(q=>Y[q]).concat([Y[i]])); const p=w*sd[i]/(mu[i]+al); if(p>best){ best=p; bi=i; } }
    picks.push(bi); }
  return picks;
}
/* ITTFA (Gemperline, J. Chem. Inf. Comput. Sci. 24 (1984) 206; Vandeginste et al., Anal. Chim. Acta 173 (1985) 253):
   needle targets in the level (concentration) domain at the OPA-selected levels, iterate x ← V_n V_nᵀ x, x[x<0] = 0, x /= max x. */
function ittfa(cols,P,n,o){
  o=o||{}; const c=P.c; const picks=o.needles||opa(cols,n); const Vn=P.V.slice(0,n); const prof=[];
  for(let k=0;k<n;k++){
    let x=new Float64Array(c); x[picks[k]]=1; let it=0;
    for(;it<(o.maxIter||2000);it++){
      const y=new Float64Array(c); Vn.forEach(v=>{ const p=dot(v,x); for(let j=0;j<c;j++) y[j]+=p*v[j]; });
      for(let j=0;j<c;j++) if(y[j]<0) y[j]=0; const mx=Math.max(...y)||1; for(let j=0;j<c;j++) y[j]/=mx;
      let d=0; for(let j=0;j<c;j++) d=Math.max(d,Math.abs(y[j]-x[j])); x=y; if(d<(o.tol||1e-10)) break;
    }
    prof.push({C:x,iter:it,needle:picks[k]});
  }
  // spectra by NNLS of each energy row on the concentration profiles
  const S=spectraFromC(cols,prof.map(p=>p.C));
  return {C:prof.map(p=>p.C),S,needles:picks,iter:prof.map(p=>p.iter)};
}
function spectraFromC(cols,C){ // D (r x c) ≈ S Cᵀ ; for each energy row i: min ||D_i· − C s_i||, s_i ≥ 0
  const n=C.length, c=cols.length, r=cols[0].length; const G=new Float64Array(n*n); for(let a=0;a<n;a++) for(let b=0;b<=a;b++){ const v=dot(C[a],C[b]); G[a*n+b]=v; G[b*n+a]=v; }
  const S=Array.from({length:n},()=>new Float64Array(r)); const h=new Float64Array(n);
  for(let i=0;i<r;i++){ for(let a=0;a<n;a++){ let s=0; const ca=C[a]; for(let j=0;j<c;j++) s+=ca[j]*cols[j][i]; h[a]=s; } const x=nnlsGram(G,h,n); for(let a=0;a<n;a++) S[a][i]=x[a]; }
  return S;
}
function concFromS(cols,S,wts){ // for each spectrum j: min ||d_j − S c_j||, c_j ≥ 0
  const n=S.length, r=S[0].length; const G=new Float64Array(n*n); for(let a=0;a<n;a++) for(let b=0;b<=a;b++){ let v=0; for(let i=0;i<r;i++) v+=S[a][i]*S[b][i]; G[a*n+b]=v; G[b*n+a]=v; }
  const C=Array.from({length:n},()=>new Float64Array(cols.length)); const h=new Float64Array(n);
  cols.forEach((d,j)=>{ for(let a=0;a<n;a++) h[a]=dot(S[a],d); const x=nnlsGram(G,h,n); for(let a=0;a<n;a++) C[a][j]=x[a]; });
  return C;
}
/* MCR-ALS with non-negativity in both modes (Tauler, Chemom. Intell. Lab. Syst. 30 (1995) 133).
   init: "simplisma" (default; profiles at the purest energy channels), "opa" (spectra of the most dissimilar levels), "ittfa", or n initial spectra.
   Each spectrum is normalized to unit area after every S step (scale moved into C). Optional closure: Σ_k C_kj = 1 (applied after the C step). */
function mcrAls(cols,n,o){
  o=o||{}; const r=cols[0].length, c=cols.length;
  let S, Cinit=null;
  if(Array.isArray(o.init)) S=o.init.map(v=>Float64Array.from(v));
  else if(o.init==="ittfa"){ const P=o.pca||pca(cols,{}); S=ittfa(cols,P,n).S; }
  else if(o.init==="opa"){ const pk=opa(cols,n); S=pk.map(j=>Float64Array.from(cols[j])); }
  else { // default: SIMPLISMA purest energy channels → initial concentration profiles → first S step
    const pk=simplisma(cols,n,o.offset); Cinit=pk.map(i=>Float64Array.from(cols,d=>Math.max(d[i],0))); S=spectraFromC(cols,Cinit); }
  const normS=()=>S.forEach(v=>{ let s=0; for(let i=0;i<r;i++) s+=v[i]; if(s>0) for(let i=0;i<r;i++) v[i]/=s; });
  normS(); let C=null, lof0=Infinity, it=0; let ss=0; cols.forEach(d=>{ for(let i=0;i<r;i++) ss+=d[i]*d[i]; });
  const lofOf=()=>{ let e=0; for(let j=0;j<c;j++){ const d=cols[j]; for(let i=0;i<r;i++){ let m=0; for(let a=0;a<n;a++) m+=S[a][i]*C[a][j]; e+=(d[i]-m)*(d[i]-m); } } return 100*Math.sqrt(e/ss); };
  let lof=Infinity;
  for(;it<(o.maxIter||1000);it++){
    C=concFromS(cols,S);
    if(o.closure){ for(let j=0;j<c;j++){ let s=0; for(let a=0;a<n;a++) s+=C[a][j]; if(s>0) for(let a=0;a<n;a++) C[a][j]/=s; } }
    S=spectraFromC(cols,C); normS();
    C=concFromS(cols,S); lof=lofOf();
    if(Math.abs(lof0-lof)<=(o.tol||1e-10)*Math.max(lof,1e-12)) break; lof0=lof;
  }
  return {S,C,lof,iter:it+1,warn:["MCR-ALS solutions are subject to rotational ambiguity: with overlapping spectra and no selective levels the recovered S and C are one member of a family of feasible solutions"]};
}

/* ============================================================
   3) Weighted NNLS fit of spectra to reference spectra (SPEC_multi §3)
   ============================================================ */
function shiftCopy(E,y,d){ // y(E − d) by linear interpolation, constant extrapolation
  const n=E.length, o=new Float64Array(n); for(let i=0;i<n;i++){ const e=E[i]-d; if(e<=E[0]){ o[i]=y[0]; continue; } if(e>=E[n-1]){ o[i]=y[n-1]; continue; }
    let lo=0, hi=n-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(E[m]<=e) lo=m; else hi=m; } const f=(e-E[lo])/(E[hi]-E[lo]); o[i]=y[lo]+f*(y[hi]-y[lo]); } return o; }
/* fitReferences(cols, refs, {E, weight:"poisson"|"none", cfac, bg:0|1|2, maxShift, shiftStep, irls:2})
   model: y ≈ Σ_k a_k ref_k(E − δ) + Σ_{m<bg} b_m P_m(E)   (a_k ≥ 0, b free sign; P_0 = 1, P_1 = (E − Ē)/(range/2))
   weights 1/v: v_i = max(y_i, 1/cfac)/cfac on pass 1, then v_i = max(μ_i, 1/cfac)/cfac (IRLS, 2 passes total).
   σ(a) = sqrt(diag((AᵀWA)_PP⁻¹)) on the passive (non-zero) set; σ·sqrt(φ) also reported (φ = χ²/(N − p)). */
function fitOne(E,y,refs,o,shift){
  const m=y.length, K=refs.length, nb=o.bg||0; const cf=o.cfac>0?o.cfac:1;
  const R=shift?refs.map(r=>shiftCopy(E,r,shift)):refs;
  const Em=(E[0]+E[m-1])/2, hw=(E[m-1]-E[0])/2||1; const B=[]; if(nb>=1) B.push(new Float64Array(m).fill(1)); if(nb>=2) B.push(Float64Array.from(E,e=>(e-Em)/hw));
  const cols=R.concat(B); const free=B.map((_,i)=>K+i);
  // w = 1/v with v = max(y, 1/cf)/cf  →  w = cf/max(y, 1/cf)
  let w=null; if(o.weight==="poisson") w=Float64Array.from(y,v=>cf/Math.max(v,1/cf));
  let res=null;
  const passes=o.weight==="poisson"?(o.irls||2):1;
  for(let p=0;p<passes;p++){
    res=nnlsCols(cols,y,w,free);
    if(p<passes-1){ const mu=new Float64Array(m); cols.forEach((cc,k)=>{ const a=res.x[k]; if(a) for(let i=0;i<m;i++) mu[i]+=a*cc[i]; }); for(let i=0;i<m;i++) w[i]=cf/Math.max(mu[i],1/cf); }
  }
  const x=res.x, n=cols.length; const mu=new Float64Array(m); cols.forEach((cc,k)=>{ if(x[k]) for(let i=0;i<m;i++) mu[i]+=x[k]*cc[i]; });
  let chi=0; for(let i=0;i<m;i++){ const e=y[i]-mu[i]; chi+=(w?w[i]:1)*e*e; }
  const act=[]; for(let k=0;k<n;k++) if(x[k]!==0||free.includes(k)) act.push(k);
  const p=act.length; const Gp=new Float64Array(p*p); act.forEach((a,i)=>act.forEach((b,j)=>{ Gp[i*p+j]=res.G[a*n+b]; }));
  const inv=invSPD(Gp,p); const sig=new Float64Array(n).fill(NaN); if(inv) act.forEach((a,i)=>{ sig[a]=Math.sqrt(inv[i*p+i]); });
  const dof=Math.max(m-p,1), phi=chi/dof;
  return {coef:x.slice(0,K),bgCoef:x.slice(K),sig:sig.slice(0,K),sigBg:sig.slice(K),chi2:chi,dof,phi,mu,active:act};
}
function fitReferences(cols,refs,o){
  o=o||{}; const E=o.E||Float64Array.from({length:cols[0].length},(_,i)=>i); const out=[];
  const areas=refs.map(r=>{ let s=0; for(let i=1;i<E.length;i++) s+=0.5*(r[i]+r[i-1])*(E[i]-E[i-1]); return s; });
  cols.forEach(y=>{
    let best=null;
    if(o.maxShift>0){ const st=o.shiftStep||Math.abs(E[1]-E[0])/2; const grid=[]; for(let d=-o.maxShift;d<=o.maxShift+1e-12;d+=st) grid.push(d);
      const ch=grid.map(d=>fitOne(E,y,refs,o,d).chi2); let k=0; ch.forEach((v,i)=>{ if(v<ch[k]) k=i; });
      let d=grid[k]; if(k>0&&k<grid.length-1){ const a=ch[k-1], b=ch[k], c=ch[k+1]; const den=a-2*b+c; if(den>0) d=grid[k]+st*0.5*(a-c)/den; }
      best=fitOne(E,y,refs,o,d); best.shift=d; }
    else { best=fitOne(E,y,refs,o,0); best.shift=0; }
    const amt=Float64Array.from(best.coef,(a,k)=>a*areas[k]); const tot=amt.reduce((s,v)=>s+v,0);
    best.fraction=Float64Array.from(amt,v=>tot>0?100*v/tot:NaN);
    // delta-method σ of fractions (coefficients treated as independent — ignores covariance, documented)
    best.sFraction=Float64Array.from(amt,(v,i)=>{ if(!(tot>0)) return NaN; let s2=0; amt.forEach((q,j)=>{ const g=100*((i===j?tot:0)-v)/(tot*tot); const sj=best.sig[j]*areas[j]; if(isFinite(sj)) s2+=g*g*sj*sj; }); return Math.sqrt(s2); });
    out.push(best);
  });
  return out;
}

/* ============================================================
   4) Series (depth profile / angle series / linescan) processing (SPEC_multi §1)
   ============================================================ */
/* group spectra (from XPS_IO, VAMAS blocks or VGD multi-spectra) into series by region.
   coordinate: meta.vamas.expVals[iVar] if present, else meta.vamas.polarTakeOff (angle), else index. */
function seriesFromSpectra(spectra,o){
  o=o||{}; const groups=new Map();
  spectra.forEach((sp,idx)=>{ const key=(sp.sample||"")+"|"+sp.region; if(!groups.has(key)) groups.set(key,[]); groups.get(key).push({sp,idx}); });
  const out=[];
  groups.forEach((list,key)=>{
    const vm=list[0].sp.meta&&list[0].sp.meta.vamas; let label="index", unit="", get=(s,i)=>i;
    if(o.coord==="angle"||(!o.coord&&vm&&vm.expVars&&!vm.expVars.length&&new Set(list.map(q=>q.sp.meta.vamas&&q.sp.meta.vamas.takeOffPolar)).size>1)){ label="take-off polar angle"; unit="deg"; get=s=>s.meta.vamas.takeOffPolar; }
    else if(vm&&vm.expVars&&vm.expVars.length){ const iv=o.iVar||0; label=vm.expVars[iv].label; unit=vm.expVars[iv].unit; get=(s,i)=>{ const v=parseFloat(s.meta.vamas.expVals[iv]); return isFinite(v)?v:i; }; }
    const pts=list.map((q,i)=>({sp:q.sp,x:get(q.sp,i)})).sort((a,b)=>a.x-b.x);
    out.push({key,region:list[0].sp.region,sample:list[0].sp.sample||"",coordLabel:label,coordUnit:unit,coord:Float64Array.from(pts.map(p=>p.x)),spectra:pts.map(p=>p.sp)});
  });
  return out;
}
/* resample all spectra of a series to a common energy grid (the overlap of all ranges, step = finest step) */
function commonGrid(spectra){
  const lo=Math.max(...spectra.map(s=>s.be[0])), hi=Math.min(...spectra.map(s=>s.be[s.be.length-1]));
  if(!(hi>lo)) throw new Error("spectra of the series do not overlap in energy");
  const same=spectra.every(s=>s.be.length===spectra[0].be.length&&Math.abs(s.be[0]-spectra[0].be[0])<1e-9&&Math.abs(s.be[s.be.length-1]-spectra[0].be[spectra[0].be.length-1])<1e-9);
  if(same) return {E:Float64Array.from(spectra[0].be),Y:spectra.map(s=>Float64Array.from(s.y)),resampled:false};
  const st=Math.min(...spectra.map(s=>s.step||((s.be[s.be.length-1]-s.be[0])/(s.be.length-1)))); const n=Math.floor((hi-lo)/st+1e-9)+1;
  const E=Float64Array.from({length:n},(_,i)=>lo+i*st); return {E,Y:spectra.map(s=>interp(s.be,s.y,E)),resampled:true};
}
function interp(x,y,xs){ const n=x.length, o=new Float64Array(xs.length); for(let k=0;k<xs.length;k++){ const e=xs[k]; if(e<=x[0]){ o[k]=y[0]; continue; } if(e>=x[n-1]){ o[k]=y[n-1]; continue; }
  let lo=0, hi=n-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(x[m]<=e) lo=m; else hi=m; } o[k]=y[lo]+(e-x[lo])/(x[hi]-x[lo])*(y[hi]-y[lo]); } return o; }

/* Linearized area coefficients a_i such that A ≈ Σ a_i y_i for a linear background through endpoint means (exact for linear; used as the
   σ approximation for Shirley/smart/Tougaard — SPEC_multi §1.2). */
function areaCoefs(E,nEnd){
  const N=E.length; const n=Math.max(1,Math.min(nEnd||5,Math.floor(N/4)||1)); const w=new Float64Array(N);
  for(let i=1;i<N;i++){ const h=E[i]-E[i-1]; w[i-1]+=h/2; w[i]+=h/2; }
  let Elo=0, Ehi=0; for(let i=0;i<n;i++){ Elo+=E[i]/n; Ehi+=E[N-1-i]/n; }
  let s0=0, s1=0; for(let k=0;k<N;k++){ const t=(E[k]-Elo)/((Ehi-Elo)||1); s0+=w[k]*(1-t); s1+=w[k]*t; }
  const a=Float64Array.from(w); for(let i=0;i<n;i++){ a[i]-=s0/n; a[N-1-i]-=s1/n; } return a;
}
/* integrate a region at every level: background of type bg (XPS_CORE.bgStatic) on [lo,hi], trapezoid area of (y − B).
   cfac = counts per CPS unit for the Poisson σ (σ² = Σ a_i² max(y_i,1/cfac)/cfac); if unknown, σ from the χ²-free estimate is not given (null). */
function integrateSeries(E,Y,o){
  o=o||{}; const C=getCore(); const bg=o.bg||{type:"linear",n:5};
  let i0=0, i1=E.length-1; if(o.lo!=null) while(i0<E.length&&E[i0]<o.lo) i0++; if(o.hi!=null) while(i1>=0&&E[i1]>o.hi) i1--;
  if(i1-i0<7) throw new Error("fewer than 8 points in the integration window");
  const Ew=E.slice(i0,i1+1); const a=areaCoefs(Ew,bg.n||5); const cf=o.cfac>0?o.cfac:null;
  return Y.map(y0=>{ const y=y0.slice(i0,i1+1); let B;
    if(bg.type==="linear"||!C){ if(bg.type!=="linear"&&!C) throw new Error("XPS_CORE not loaded (needed for "+bg.type+")"); B=C?C.bgLinear(Ew,y,bg.n||5):null; }
    else B=C.bgStatic(Ew,y,bg);
    if(!B) throw new Error("background type "+bg.type+" not usable here (use a static background)");
    let A=0; for(let i=1;i<Ew.length;i++) A+=0.5*((y[i]-B[i])+(y[i-1]-B[i-1]))*(Ew[i]-Ew[i-1]);
    let s=null; if(cf){ let v=0; for(let i=0;i<y.length;i++) v+=a[i]*a[i]*Math.max(y[i],1/cf)/cf; s=Math.sqrt(v); }
    return {area:A,sArea:s,converged:B.converged!==false};
  });
}
/* apply a peak model at every level: fitFn(E, y, levelIndex) → {components:[{name, area, sArea}]}  (plug the core fitter in here) */
function modelSeries(E,Y,fitFn){ return Y.map((y,j)=>{ const r=fitFn(E,y,j)||{}; return {components:(r.components||[]).map(c=>({name:c.name,area:c.area,sArea:c.sArea!=null?c.sArea:null})),raw:r}; }); }

/* composition vs level: rows = [{label, el, S, BE, hv, PE, tf}] (one per quantified region/component),
   areas[k][j] = {area, sArea} of row k at level j. Uses XPS_CORE.quantify (SPEC §4) per level, imfp option passed through. */
function compositionSeries(rows,areas,opts){
  const C=getCore(); if(!C) throw new Error("XPS_CORE not loaded");
  const nl=areas[0].length; const out=[];
  for(let j=0;j<nl;j++){ const q=C.quantify(rows.map((r,k)=>Object.assign({},r,{area:areas[k][j].area,sArea:areas[k][j].sArea})),opts||{});
    out.push(q.map(x=>({label:x.label,at:x.at,sAt:x.sAt,nrm:x.nrm,sNrm:x.sNrm,nd:x.nd}))); }
  return out;   // out[level][row]
}
/* subtract a static background (XPS_CORE.bgStatic: linear | shirley | smart | tougaard) from every spectrum on [lo, hi] — recommended before
   MCR-ALS / target testing: a level-dependent background is otherwise an extra, poorly separable factor (SPEC_multi §2.6). */
function subtractBackground(E,cols,o){
  o=o||{}; const C=getCore(); if(!C) throw new Error("XPS_CORE not loaded"); let i0=0, i1=E.length-1; if(o.lo!=null) while(i0<E.length&&E[i0]<o.lo) i0++; if(o.hi!=null) while(i1>=0&&E[i1]>o.hi) i1--;
  const Ew=E.slice(i0,i1+1); const bg=o.bg||{type:"shirley",n:5};
  return {E:Ew,cols:cols.map(y0=>{ const y=y0.slice(i0,i1+1); const B=C.bgStatic(Ew,y,bg); return Float64Array.from(y,(v,i)=>v-B[i]); }),i0,i1};
}
/* local quadratic least-squares derivative on a non-uniform grid (hw points each side, asymmetric at the ends) */
function derivQ(z,y,hw){
  const n=z.length, d=new Float64Array(n);
  for(let i=0;i<n;i++){ let j0=Math.max(0,i-hw), j1=Math.min(n-1,i+hw); while(j1-j0<2){ if(j0>0) j0--; else if(j1<n-1) j1++; else break; }
    const deg=j1-j0>=2?2:1; let S=[0,0,0,0,0], T=[0,0,0];
    for(let j=j0;j<=j1;j++){ const t=z[j]-z[i]; let p=1; for(let k=0;k<5;k++){ S[k]+=p; p*=t; } T[0]+=y[j]; T[1]+=y[j]*t; T[2]+=y[j]*t*t; }
    if(deg===2){ const A=[S[0],S[1],S[2],S[1],S[2],S[3],S[2],S[3],S[4]]; const x=solveGE(A,T,3); d[i]=x?x[1]:0; }
    else { const det=S[0]*S[2]-S[1]*S[1]; d[i]=det?(S[0]*T[1]-S[1]*T[0])/det:0; } }
  return d;
}
/* Information-depth (escape-depth) correction of a sputter depth profile (SPEC_multi §1.4).
   Measured signal with exponential information-depth response (Hofmann MRI model, Surf. Interface Anal. 21 (1994) 673, term g_λ):
     I(z0) = ∫_{z0}^{∞} c(z) exp(−(z − z0)/λ_e)/λ_e dz,   λ_e = L·cosθ
   ⇒ exact inverse  c(z0) = I(z0) − λ_e dI/dz0  (derived here by differentiating; noise-amplifying, so the derivative is a local quadratic fit).
   X[k] = normalized intensities (n = A/(S·T·λ)) of element k vs depth z (nm); lamE[k] in nm. Returns corrected n (clipped ≥ 0) and at%. */
function infoDepthCorrect(z,X,lamE,o){
  o=o||{}; const hw=o.hw||2; const K=X.length, n=z.length;
  const Cn=X.map((x,k)=>{ const d=derivQ(z,x,hw); return Float64Array.from(x,(v,i)=>Math.max(0,v-lamE[k]*d[i])); });
  const at=Cn.map(()=>new Float64Array(n)); for(let i=0;i<n;i++){ let s=0; for(let k=0;k<K;k++) s+=Cn[k][i]; for(let k=0;k<K;k++) at[k][i]=s>0?100*Cn[k][i]/s:NaN; }
  return {n:Cn,at,warn:["valid for a constant information depth and no sputter-induced broadening (MRI w = σ = 0); amplifies noise — inspect before use"]};
}
/* depth scale: z_j = Σ rate·Δt (constant rate) or user array */
function depthScale(t,rate){ return Float64Array.from(t,v=>v*rate); }

/* ============================================================
   5) Worker entry (optional): self.onmessage → XPS_MULTI.run(job)
   ============================================================ */
function run(job){
  const f=job.fn; const a=job.args||[];
  const fns={pca,reconstruct,wlra,noiseFactors,subtractBackground,targetTest,ittfa,mcrAls,fitReferences,integrateSeries,compositionSeries,infoDepthCorrect};
  if(!fns[f]) throw new Error("unknown function "+f); return fns[f].apply(null,a);
}

const XPS_MULTI={version:"1.0",dot,norm2,lgamma,betai,fSurvival,cholesky,solveSPD,invSPD,solveGE,symEig,svdCols,svdGram,fixSign,nnlsGram,nnlsCols,
  preprocess,undo,malinowski,pca,reconstruct,wlra,noiseFactors,subtractBackground,targetTest,opa,simplisma,ittfa,spectraFromC,concFromS,mcrAls,
  shiftCopy,interp,fitReferences,fitOne,seriesFromSpectra,commonGrid,areaCoefs,integrateSeries,modelSeries,compositionSeries,derivQ,infoDepthCorrect,depthScale,run};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_MULTI; else root.XPS_MULTI=XPS_MULTI;
})(typeof self!=="undefined"?self:this);
