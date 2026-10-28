/* ============================================================
   INNO Analytics — XPS imaging / mapping engine (spectral image stacks)
   DOM-free (browser, Web Worker, node). Load after xps_multi.js. SPEC: modules/xps/SPEC_multi.md §5.
   Stack = {nx, ny, nE, E: Float64Array(nE) (BE ascending), data: Float32Array(nx·ny·nE) PIXEL-MAJOR (data[p·nE + e], p = y·nx + x),
            cfac: counts per data unit (1 if data are counts; dwell·scans if CPS; null if unknown), x?: Float64Array(nx), y?: Float64Array(ny)}
   All per-pixel loops avoid allocation so 256×256×100 runs in a Worker in ~1 s per operation.
   ============================================================ */
(function(root){
"use strict";
function getMulti(){ if(root&&root.XPS_MULTI) return root.XPS_MULTI; if(typeof require==="function"){ try{ return require("./xps_multi.js"); }catch(e){} } throw new Error("XPS_MULTI not loaded"); }

/* ---------- construction ---------- */
function makeStack(o){
  const nx=o.nx, ny=o.ny, nE=o.E.length, np=nx*ny; let E=Float64Array.from(o.E), src=o.data, flip=false;
  if(src.length!==np*nE) throw new Error("data length "+src.length+" ≠ nx·ny·nE = "+np*nE);
  if(nE>1&&E[0]>E[nE-1]){ E=E.reverse(); flip=true; }
  const data=new Float32Array(np*nE);
  if(o.layout==="energy"){ for(let e=0;e<nE;e++){ const ee=flip?nE-1-e:e; const off=e*np; for(let p=0;p<np;p++) data[p*nE+ee]=src[off+p]; } }
  else { if(flip){ for(let p=0;p<np;p++) for(let e=0;e<nE;e++) data[p*nE+nE-1-e]=src[p*nE+e]; } else data.set(src); }
  return {nx,ny,nE,E,data,cfac:o.cfac!=null?o.cfac:null,x:o.x?Float64Array.from(o.x):null,y:o.y?Float64Array.from(o.y):null,name:o.name||""};
}
/* from per-position spectra (e.g. VAMAS MAP blocks: meta.vamas.x, meta.vamas.y). Positions are snapped to the sorted unique x / y values. */
function fromSpectra(spectra,o){
  o=o||{}; const xs=[...new Set(spectra.map(s=>+s.meta.vamas.x))].sort((a,b)=>a-b), ys=[...new Set(spectra.map(s=>+s.meta.vamas.y))].sort((a,b)=>a-b);
  const nx=xs.length, ny=ys.length, M=getMulti(); const g=M.commonGrid(spectra); const nE=g.E.length; const data=new Float32Array(nx*ny*nE); const filled=new Uint8Array(nx*ny);
  spectra.forEach((s,i)=>{ const ix=xs.indexOf(+s.meta.vamas.x), iy=ys.indexOf(+s.meta.vamas.y); const p=iy*nx+ix; filled[p]=1; data.set(g.Y[i],p*nE); });
  const st=makeStack({nx,ny,E:g.E,data,cfac:spectra[0].c||null,x:xs,y:ys,name:spectra[0].region});
  st.missing=filled.length-filled.reduce((a,b)=>a+b,0); return st;
}
/* from a series of images at energies: images = [{E, data: (nx·ny) array}] */
function fromImages(images,nx,ny,o){
  const nE=images.length, np=nx*ny; const ord=images.map((im,i)=>i).sort((a,b)=>images[a].E-images[b].E);
  const data=new Float32Array(np*nE); ord.forEach((i,e)=>{ const d=images[i].data; for(let p=0;p<np;p++) data[p*nE+e]=d[p]; });
  return makeStack(Object.assign({nx,ny,E:ord.map(i=>images[i].E),data},o||{}));
}
function idxRange(E,lo,hi){ let i0=0, i1=E.length-1; if(lo!=null) while(i0<E.length&&E[i0]<lo) i0++; if(hi!=null) while(i1>=0&&E[i1]>hi) i1--; return [i0,i1]; }

/* ---------- simple maps ---------- */
/* image at an energy window: sum of channels with lo ≤ E ≤ hi (mode "sum" | "mean") */
function windowImage(st,lo,hi,mode){
  const [i0,i1]=idxRange(st.E,lo,hi); const np=st.nx*st.ny, nE=st.nE, d=st.data, out=new Float64Array(np);
  for(let p=0;p<np;p++){ let s=0; const o=p*nE; for(let e=i0;e<=i1;e++) s+=d[o+e]; out[p]=mode==="mean"?s/(i1-i0+1):s; } return out;
}
/* background-subtracted peak area map on [lo, hi]: bg "none" | "linear" (endpoint means of n points) | "shirley" (SPEC §2, iterative)
   σ² = Σ a_i² max(y_i, 1/cfac)/cfac with a = linearized area coefficients (exact for linear/none; approximation for Shirley — SPEC_multi §1.2) */
function areaMap(st,o){
  o=o||{}; const M=getMulti(); const [i0,i1]=idxRange(st.E,o.lo,o.hi); const N=i1-i0+1; if(N<4) throw new Error("window has < 4 points");
  const E=st.E.subarray(i0,i1+1), bg=o.bg||"linear", nEnd=Math.max(1,Math.min(o.n||3,Math.floor(N/4)||1));
  const w=new Float64Array(N); for(let i=1;i<N;i++){ const h=E[i]-E[i-1]; w[i-1]+=h/2; w[i]+=h/2; }
  const a=bg==="none"?w:M.areaCoefs(E,nEnd); const np=st.nx*st.ny, nE=st.nE, d=st.data, cf=st.cfac>0?st.cfac:null;
  const area=new Float64Array(np), sArea=cf?new Float64Array(np):null, conv=bg==="shirley"?new Uint8Array(np):null;
  const y=new Float64Array(N), B=new Float64Array(N), Q=new Float64Array(N);
  for(let p=0;p<np;p++){
    const o0=p*nE+i0; for(let i=0;i<N;i++) y[i]=d[o0+i];
    if(bg==="shirley"){ let ylo=0,yhi=0; for(let i=0;i<nEnd;i++){ ylo+=y[i]; yhi+=y[N-1-i]; } ylo/=nEnd; yhi/=nEnd; const dd=yhi-ylo;
      for(let i=0;i<N;i++) B[i]=ylo; let ok=Math.abs(dd)<1e-300; const tol=1e-7*(Math.abs(dd)+1e-12);
      for(let it=0;it<200&&!ok;it++){ Q[0]=0; for(let i=1;i<N;i++) Q[i]=Q[i-1]+0.5*((y[i]-B[i])+(y[i-1]-B[i-1]))*(E[i]-E[i-1]); const Qt=Q[N-1]; let mx=0;
        for(let i=0;i<N;i++){ const nb=Math.abs(Qt)>1e-300?ylo+dd*Q[i]/Qt:ylo; const df=Math.abs(nb-B[i]); if(df>mx) mx=df; B[i]=nb; } if(mx<tol) ok=true; }
      let A=0; for(let i=1;i<N;i++) A+=0.5*((y[i]-B[i])+(y[i-1]-B[i-1]))*(E[i]-E[i-1]); area[p]=A; conv[p]=ok?1:0; }
    else { let A=0; for(let i=0;i<N;i++) A+=a[i]*y[i]; area[p]=A; }
    if(cf){ let v=0; for(let i=0;i<N;i++) v+=a[i]*a[i]*Math.max(y[i],1/cf)/cf; sArea[p]=Math.sqrt(v); }
  }
  return {area,sArea,converged:conv,bg,window:[E[0],E[N-1]]};
}
/* at% maps: items = [{label, area: map, sArea?: map, S, T?, L?}] → n = A/(S·T·L); at% = 100 n/Σn; pixels with Σn ≤ minTotal → NaN.
   n.d. rule as SPEC §4 (A ≤ 0 or A < 2σ ⇒ excluded from the sum at that pixel). σ(at%) by the delta method (independent areas). */
function atomicMaps(items,o){
  o=o||{}; const np=items[0].area.length, K=items.length; const den=items.map(it=>it.S*(it.T||1)*(it.L||1));
  const at=items.map(()=>new Float64Array(np)), sAt=items.map(()=>new Float64Array(np)); const n=new Float64Array(K), sn=new Float64Array(K), use=new Uint8Array(K);
  for(let p=0;p<np;p++){ let tot=0;
    for(let k=0;k<K;k++){ const A=items[k].area[p], s=items[k].sArea?items[k].sArea[p]:null; n[k]=A/den[k]; sn[k]=s!=null?s/den[k]:0; use[k]=(A>0&&!(s!=null&&s>0&&A<2*s))?1:0; if(use[k]) tot+=n[k]; }
    if(!(tot>(o.minTotal||0))){ for(let k=0;k<K;k++){ at[k][p]=NaN; sAt[k][p]=NaN; } continue; }
    for(let k=0;k<K;k++){ if(!use[k]){ at[k][p]=0; sAt[k][p]=NaN; continue; } at[k][p]=100*n[k]/tot; let v=0; for(let j=0;j<K;j++){ if(!use[j]) continue; const g=100*((k===j?tot:0)-n[k])/(tot*tot); v+=g*g*sn[j]*sn[j]; } sAt[k][p]=Math.sqrt(v); } }
  return {at,sAt,labels:items.map(it=>it.label)};
}
function ratioMap(a,b){ return Float64Array.from(a,(v,i)=>b[i]!==0?v/b[i]:NaN); }
/* ROI average spectrum: mask (Uint8Array np) or rect {x0,y0,x1,y1} (inclusive). σ of the mean from Poisson: sqrt(Σ y/cfac)/n_pix */
function roiSpectrum(st,roi){
  const np=st.nx*st.ny, nE=st.nE, d=st.data; let mask=roi&&roi.length===np?roi:null; if(!mask){ mask=new Uint8Array(np); const r=roi||{x0:0,y0:0,x1:st.nx-1,y1:st.ny-1}; for(let y=r.y0;y<=r.y1;y++) for(let x=r.x0;x<=r.x1;x++) mask[y*st.nx+x]=1; }
  const s=new Float64Array(nE); let n=0; for(let p=0;p<np;p++){ if(!mask[p]) continue; n++; const o=p*nE; for(let e=0;e<nE;e++) s[e]+=d[o+e]; }
  const mean=Float64Array.from(s,v=>n?v/n:NaN); const cf=st.cfac>0?st.cfac:null;
  return {E:Float64Array.from(st.E),y:mean,sum:s,n,sigma:cf?Float64Array.from(s,v=>Math.sqrt(Math.max(v,0)/cf)/n):null};
}
/* spatial binning b×b (sum keeps Poisson statistics; mean divides by b², cfac multiplied accordingly) */
function binStack(st,b,mode){
  const nx=Math.floor(st.nx/b), ny=Math.floor(st.ny/b), nE=st.nE, d=st.data, out=new Float32Array(nx*ny*nE), acc=new Float64Array(nE);
  for(let Y=0;Y<ny;Y++) for(let X=0;X<nx;X++){ acc.fill(0); for(let dy=0;dy<b;dy++) for(let dx=0;dx<b;dx++){ const o=((Y*b+dy)*st.nx+(X*b+dx))*nE; for(let e=0;e<nE;e++) acc[e]+=d[o+e]; }
    const o=(Y*nx+X)*nE; for(let e=0;e<nE;e++) out[o+e]=mode==="mean"?acc[e]/(b*b):acc[e]; }
  return {nx,ny,nE,E:Float64Array.from(st.E),data:out,cfac:st.cfac>0?(mode==="mean"?st.cfac*b*b:st.cfac):null,x:null,y:null,name:st.name+" bin"+b};
}

/* per-pixel static background subtraction on [lo, hi] (bg "linear" | "shirley", SPEC §2 endpoint means of n points) → new stack on the window.
   Use before nnlsMap / pcaStack when the background is a Shirley step (a constant/slope term cannot represent it). */
function subtractBackgroundStack(st,o){
  o=o||{}; const [i0,i1]=idxRange(st.E,o.lo,o.hi); const N=i1-i0+1, nEnd=Math.max(1,Math.min(o.n||3,Math.floor(N/4)||1)); const E=st.E.subarray(i0,i1+1);
  const np=st.nx*st.ny, nE=st.nE, d=st.data, out=new Float32Array(np*N); const y=new Float64Array(N), B=new Float64Array(N), Q=new Float64Array(N); const shir=(o.bg||"shirley")==="shirley";
  let El=0, Eh=0; for(let i=0;i<nEnd;i++){ El+=E[i]/nEnd; Eh+=E[N-1-i]/nEnd; }
  for(let p=0;p<np;p++){ const o0=p*nE+i0; for(let i=0;i<N;i++) y[i]=d[o0+i]; let ylo=0,yhi=0; for(let i=0;i<nEnd;i++){ ylo+=y[i]; yhi+=y[N-1-i]; } ylo/=nEnd; yhi/=nEnd; const dd=yhi-ylo;
    if(shir){ for(let i=0;i<N;i++) B[i]=ylo; let ok=Math.abs(dd)<1e-300; const tol=1e-7*(Math.abs(dd)+1e-12);
      for(let it=0;it<200&&!ok;it++){ Q[0]=0; for(let i=1;i<N;i++) Q[i]=Q[i-1]+0.5*((y[i]-B[i])+(y[i-1]-B[i-1]))*(E[i]-E[i-1]); const Qt=Q[N-1]; let mx=0;
        for(let i=0;i<N;i++){ const nb=Math.abs(Qt)>1e-300?ylo+dd*Q[i]/Qt:ylo; const df=Math.abs(nb-B[i]); if(df>mx) mx=df; B[i]=nb; } if(mx<tol) ok=true; } }
    else for(let i=0;i<N;i++) B[i]=ylo+dd*(E[i]-El)/((Eh-El)||1);
    const oo=p*N; for(let i=0;i<N;i++) out[oo+i]=y[i]-B[i]; }
  return Object.assign({},st,{nE:N,E:Float64Array.from(E),data:out,bgSubtracted:true,name:st.name+" −bg"});
}
/* ---------- PCA on the stack (energies = variables, pixels = observations) ----------
   same preprocessing definitions as XPS_MULTI.preprocess ("none" | "poisson" | "area"; centre "none"|"mean"), computed on the fly,
   r×r Gram matrix → Jacobi eigen-decomposition; scores = Uᵀx/s for the leading nKeep factors. */
function pcaStack(st,o){
  o=o||{}; const M=getMulti(); const np=st.nx*st.ny, r=st.nE, d=st.data; const scale=o.scale||"none", center=o.center||"none"; const k=Math.min(o.nKeep||8,r);
  const rowW=new Float64Array(r).fill(1), colW=new Float64Array(np).fill(1);
  if(scale==="poisson"||scale==="area"){ const g=new Float64Array(r); for(let p=0;p<np;p++){ let s=0; const off=p*r; for(let e=0;e<r;e++){ g[e]+=d[off+e]; s+=d[off+e]; } colW[p]=s; }
    if(scale==="poisson"){ let gm=0; for(let e=0;e<r;e++){ g[e]/=np; gm=Math.max(gm,g[e]); } let hm=0; for(let p=0;p<np;p++){ colW[p]/=r; hm=Math.max(hm,colW[p]); }
      for(let e=0;e<r;e++) rowW[e]=1/Math.sqrt(Math.max(g[e],gm*1e-6||1e-300)); for(let p=0;p<np;p++) colW[p]=1/Math.sqrt(Math.max(colW[p],hm*1e-6||1e-300)); }
    else for(let p=0;p<np;p++) colW[p]=colW[p]>0?1/colW[p]:1; }
  const mean=new Float64Array(r); if(center==="mean"){ for(let p=0;p<np;p++){ const off=p*r, cw=colW[p]; for(let e=0;e<r;e++) mean[e]+=d[off+e]*rowW[e]*cw; } for(let e=0;e<r;e++) mean[e]/=np; }
  const G=new Float64Array(r*r), x=new Float64Array(r);
  for(let p=0;p<np;p++){ const off=p*r, cw=colW[p]; for(let e=0;e<r;e++) x[e]=d[off+e]*rowW[e]*cw-mean[e];
    for(let a=0;a<r;a++){ const xa=x[a]; if(xa===0) continue; const o2=a*r; for(let b=0;b<=a;b++) G[o2+b]+=xa*x[b]; } }
  for(let a=0;a<r;a++) for(let b=0;b<a;b++) G[b*r+a]=G[a*r+b];
  const Ei=M.symEig(G,r); const eig=Float64Array.from(Ei.w,v=>Math.max(v,0)); const s=Float64Array.from(eig,Math.sqrt);
  const U=Ei.V.slice(0,k).map(v=>Float64Array.from(v)); const scores=U.map(()=>new Float64Array(np));
  for(let p=0;p<np;p++){ const off=p*r, cw=colW[p]; for(let e=0;e<r;e++) x[e]=d[off+e]*rowW[e]*cw-mean[e]; for(let q=0;q<k;q++){ let t=0; const u=U[q]; for(let e=0;e<r;e++) t+=u[e]*x[e]; scores[q][p]=s[q]>0?t/s[q]:0; } }
  U.forEach((u,q)=>M.fixSign(u,scores[q]));
  const tot=eig.reduce((a,b)=>a+b,0);
  const mal=M.malinowski(eig,r,center==="mean"?np-1:np,o.alpha);
  return {U,s,eig,varPct:Float64Array.from(eig,v=>100*v/tot),scores,rowW,colW,mean,scale,center,malinowski:mal,nx:st.nx,ny:st.ny,
    warn:["small eigenvalues from the Gram matrix have reduced relative accuracy (~1e-8 of the largest); factor counts rely on the leading ones"]};
}
/* PCA noise reduction: reconstruct the stack from n factors (original units) */
function pcaReconstruct(st,P,n){
  const np=st.nx*st.ny, r=st.nE, out=new Float32Array(np*r); n=Math.min(n,P.U.length);
  for(let p=0;p<np;p++){ const off=p*r, cw=P.colW[p]; for(let e=0;e<r;e++){ let v=P.mean[e]; for(let q=0;q<n;q++) v+=P.s[q]*P.scores[q][p]*P.U[q][e]; out[off+e]=v/(P.rowW[e]*cw); } }
  return Object.assign({},st,{data:out,name:st.name+" PCA"+n});
}
/* Noise test on the stack (SPEC_multi §5.4). method "mlpca" (default): weighted low-rank fit (XPS_MULTI.noiseFactors) on an evenly strided
   subset of nSub pixels (default 2048) — statistically correct for Poisson data; method "pca": full stack, Pearson χ² of the reconstruction of P
   (fast; biased when the intensity dynamic range is large, because unweighted factors leak high-variance noise into low-count channels). */
function noiseFactorsStack(st,P,o){
  o=o||{}; if((o.method||"mlpca")==="mlpca"){ const M=getMulti(); const np=st.nx*st.ny, r=st.nE, nSub=Math.min(o.nSub||2048,np), stp=np/nSub; const cols=[];
    for(let k=0;k<nSub;k++){ const p=Math.floor(k*stp); cols.push(Float64Array.from(st.data.subarray(p*r,(p+1)*r))); }
    const Ps=M.pca(cols,{scale:"poisson"}); const res=M.noiseFactors(cols,Ps,st.cfac>0?st.cfac:1,{nMax:o.nMax||Math.min(10,r-1),z:o.z}); res.nSub=nSub; return res; }
  const z=o.z!=null?o.z:2, cf=st.cfac>0?st.cfac:1; const np=st.nx*st.ny, r=st.nE, d=st.data; const nMax=Math.min(o.nMax||P.U.length,P.U.length);
  const cc=P.center==="mean"?np-1:np; const tab=[]; let nNoise=null; const rec=new Float64Array(r);
  for(let n=0;n<=nMax;n++){ let chi=0;
    for(let p=0;p<np;p++){ const off=p*r, cw=P.colW[p]; for(let e=0;e<r;e++){ let v=P.mean[e]; for(let q=0;q<n;q++) v+=P.s[q]*P.scores[q][p]*P.U[q][e]; rec[e]=v/(P.rowW[e]*cw); }
      for(let e=0;e<r;e++){ const y=d[off+e], v=Math.max(n>0?rec[e]:y,1/cf)/cf, df=y-rec[e]; chi+=df*df/v; } }
    const dof=(r-n)*(cc-n), lim=dof+z*Math.sqrt(2*dof); tab.push({n,chi2:chi,dof,ratio:chi/dof,limit:lim,within:chi<=lim}); if(nNoise==null&&chi<=lim){ nNoise=n; if(!o.full) break; } }
  return {table:tab,nNoise,z,method:"pca"};
}

/* ---------- per-pixel NNLS to reference spectra → chemical-state maps ----------
   refs: K arrays on st.E;  bg: 0 | 1 (constant) | 2 (constant + slope), sign-free;  weight: "none" | "poisson" (IRLS, 2 passes: v = max(y,1/cf)/cf then max(μ,1/cf)/cf)
   outputs coef maps (K), σ maps (from (AᵀWA)_PP⁻¹, NOT scaled by φ), χ²/dof map. */
function nnlsMap(st,refs,o){
  o=o||{}; const M=getMulti(); const np=st.nx*st.ny, r=st.nE, d=st.data, K=refs.length, nb=o.bg||0, n=K+nb; const cf=st.cfac>0?st.cfac:1;
  const E=st.E, Em=(E[0]+E[r-1])/2, hw=(E[r-1]-E[0])/2||1; const cols=refs.map(v=>Float64Array.from(v)); if(nb>=1) cols.push(new Float64Array(r).fill(1)); if(nb>=2) cols.push(Float64Array.from(E,e=>(e-Em)/hw));
  const free=[]; for(let i=K;i<n;i++) free.push(i);
  const coef=cols.map(()=>new Float64Array(np)), sig=cols.map(()=>new Float64Array(np)), redchi=new Float64Array(np);
  const G0=new Float64Array(n*n); for(let a=0;a<n;a++) for(let b=0;b<=a;b++){ const v=M.dot(cols[a],cols[b]); G0[a*n+b]=v; G0[b*n+a]=v; }
  const y=new Float64Array(r), w=new Float64Array(r), h=new Float64Array(n), G=new Float64Array(n*n), mu=new Float64Array(r);
  const pois=o.weight==="poisson"; const passes=pois?2:1;
  for(let p=0;p<np;p++){
    const off=p*r; for(let e=0;e<r;e++) y[e]=d[off+e];
    if(pois) for(let e=0;e<r;e++) w[e]=cf/Math.max(y[e],1/cf);
    let x=null;
    for(let pass=0;pass<passes;pass++){
      if(pois){ for(let a=0;a<n;a++){ const ca=cols[a]; let s=0; for(let e=0;e<r;e++) s+=w[e]*ca[e]*y[e]; h[a]=s; for(let b=0;b<=a;b++){ const cb=cols[b]; let t=0; for(let e=0;e<r;e++) t+=w[e]*ca[e]*cb[e]; G[a*n+b]=t; G[b*n+a]=t; } } }
      else { G.set(G0); for(let a=0;a<n;a++){ const ca=cols[a]; let s=0; for(let e=0;e<r;e++) s+=ca[e]*y[e]; h[a]=s; } }
      x=M.nnlsGram(G,h,n,free);
      mu.fill(0); for(let a=0;a<n;a++){ const xa=x[a]; if(xa){ const ca=cols[a]; for(let e=0;e<r;e++) mu[e]+=xa*ca[e]; } }
      if(pois&&pass<passes-1) for(let e=0;e<r;e++) w[e]=cf/Math.max(mu[e],1/cf);
    }
    let chi=0; for(let e=0;e<r;e++){ const df=y[e]-mu[e]; chi+=(pois?w[e]:1)*df*df; }
    const act=[]; for(let a=0;a<n;a++) if(x[a]!==0||a>=K) act.push(a); const q=act.length; const Gp=new Float64Array(q*q); for(let i=0;i<q;i++) for(let j=0;j<q;j++) Gp[i*q+j]=G[act[i]*n+act[j]];
    const inv=M.invSPD(Gp,q); for(let a=0;a<n;a++){ coef[a][p]=x[a]; sig[a][p]=NaN; } if(inv) act.forEach((a,i)=>{ sig[a][p]=Math.sqrt(inv[i*q+i]); });
    redchi[p]=chi/Math.max(r-q,1);
  }
  const areas=refs.map(v=>{ let s=0; for(let i=1;i<r;i++) s+=0.5*(v[i]+v[i-1])*(E[i]-E[i-1]); return s; });
  const frac=refs.map(()=>new Float64Array(np)); for(let p=0;p<np;p++){ let t=0; for(let k=0;k<K;k++) t+=coef[k][p]*areas[k]; for(let k=0;k<K;k++) frac[k][p]=t>0?100*coef[k][p]*areas[k]/t:NaN; }
  return {coef:coef.slice(0,K),sig:sig.slice(0,K),bgCoef:coef.slice(K),redchi,fraction:frac,refAreas:areas,
    note:pois?"σ assumes Poisson variance with the given cfac; compare redchi to 1":"unweighted: σ are in units of the (unknown) noise σ — multiply by sqrt(redchi)"};
}

/* ---------- Worker entry ---------- */
function run(job){ const fns={makeStack,fromImages,windowImage,areaMap,subtractBackgroundStack,atomicMaps,ratioMap,roiSpectrum,binStack,pcaStack,pcaReconstruct,noiseFactorsStack,nnlsMap}; const f=fns[job.fn]; if(!f) throw new Error("unknown function "+job.fn); return f.apply(null,job.args||[]); }

const XPS_IMAGE={version:"1.0",makeStack,fromSpectra,fromImages,idxRange,windowImage,areaMap,subtractBackgroundStack,atomicMaps,ratioMap,roiSpectrum,binStack,pcaStack,pcaReconstruct,noiseFactorsStack,nnlsMap,run};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_IMAGE; else root.XPS_IMAGE=XPS_IMAGE;
})(typeof self!=="undefined"?self:this);
