/* XPS imaging engine tests — node test/test_image.js (from modules/xps)
   A) engine ↔ independent numpy/scipy reference (test/ref/ref_image.py → ref_image.json)
   B) synthetic truth: chemical-state maps (NNLS σ coverage), at% maps, factor count, PCA denoising, VAMAS MAP file end-to-end
   C) performance: 256×256×100 stack, every operation timed (budget: a few seconds each in a Worker) */
"use strict";
const path=require("path"), fs=require("fs");
const IM=require("../src/xps_image.js"), M=require("../src/xps_multi.js");
const {rng,reporter,maxRel,F64,gl,stat,trapz}=require("./synth_util.js");
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref/ref_image.json"),"utf8"));
const R=reporter("XPS_IMAGE"), ok=R.ok; const QUICK=!!process.env.QUICK;

const st=IM.makeStack({nx:REF.nx,ny:REF.ny,E:REF.E,data:REF.data,cfac:REF.cfac});
console.log("== A1 window image, area maps (linear / Shirley) and Poisson σ vs reference");
{ const w=IM.windowImage(st,REF.window.lo,REF.window.hi); ok(maxRel(w,REF.window.sum)<1e-12,"window sum image",maxRel(w,REF.window.sum).toExponential(1));
  const a=IM.areaMap(st,{lo:REF.area.lo,hi:REF.area.hi,bg:"linear",n:REF.area.n}); const sc=Math.max(...REF.area.linear.map(Math.abs));
  ok(maxRel(a.area,REF.area.linear,sc)<1e-11&&maxRel(a.sArea,REF.area.sigma)<1e-11,"linear-background area and σ",`${maxRel(a.area,REF.area.linear,sc).toExponential(1)}, σ ${maxRel(a.sArea,REF.area.sigma).toExponential(1)}`);
  const s=IM.areaMap(st,{lo:REF.area.lo,hi:REF.area.hi,bg:"shirley",n:REF.area.n}); ok(maxRel(s.area,REF.area.shirley,sc)<1e-9,"Shirley area",maxRel(s.area,REF.area.shirley,sc).toExponential(1));
  const A2=IM.areaMap(st,{lo:REF.atomic.lo2,hi:REF.atomic.hi2,bg:"linear",n:3}); ok(maxRel(A2.area,REF.atomic.A2,Math.max(...REF.atomic.A2))<1e-11,"second area map");
  const at=IM.atomicMaps([{label:"a",area:a.area,sArea:a.sArea,S:REF.atomic.S[0]},{label:"b",area:A2.area,S:REF.atomic.S[1]}]);
  let e=0; at.at[0].forEach((v,p)=>{ const r=REF.atomic.at1[p]; if(r==null||!isFinite(r)){ if(isFinite(v)) e=Infinity; return; } e=Math.max(e,Math.abs(v-r),Math.abs(at.at[1][p]-REF.atomic.at2[p])); const sr=REF.atomic.sat1[p]; if(sr!=null&&isFinite(sr)) e=Math.max(e,Math.abs(at.sAt[0][p]-sr)); });
  ok(e<1e-10,"at% maps, σ(at%), n.d. rule",`max abs ${e.toExponential(1)}`); }

console.log("== A2 ROI spectrum, 2×2 binning");
{ const r=IM.roiSpectrum(st,REF.roi.rect); ok(maxRel(r.y,REF.roi.mean)<1e-12&&maxRel(r.sigma,REF.roi.sigma)<1e-12,"ROI mean spectrum and σ of the mean");
  const b=IM.binStack(st,2); ok(maxRel(b.data,REF.bin2,1)<1e-9,"bin 2×2 (sums)"); }

console.log("== A3 PCA on the stack (Gram + Jacobi) vs numpy SVD: eigenvalues, loadings, scores, IND, reconstruction");
REF.pca.forEach(p=>{ const P=IM.pcaStack(st,{scale:p.scale,center:p.center,nKeep:3}); const l0=p.eig[0];
  const ee=maxRel(P.eig.slice(0,10),p.eig,l0); let eu=0, es=0; for(let k=0;k<3;k++){ eu=Math.max(eu,maxRel(P.U[k],p.U3[k],1)); es=Math.max(es,maxRel(P.scores[k],p.scores3[k],1)); }
  const ind=P.malinowski.table.slice(0,10).map(o=>o.IND); const ei=maxRel(ind,p.IND);
  const rec=IM.pcaReconstruct(st,P,3); const er=maxRel(rec.data,p.rec3,Math.max(...p.rec3));
  ok(ee<1e-11&&eu<1e-7&&es<1e-7&&ei<1e-5&&er<1e-6&&P.malinowski.nIND===p.nIND,`scale=${p.scale} center=${p.center}`,`λ ${ee.toExponential(1)}, U ${eu.toExponential(1)}, scores ${es.toExponential(1)}, IND ${ei.toExponential(1)}, rec ${er.toExponential(1)}, nIND ${P.malinowski.nIND}/${p.nIND}`); });

console.log("== A4 noise table (method 'pca', Poisson-scaled reconstruction) vs reference");
{ const P=IM.pcaStack(st,{scale:"poisson",nKeep:5}); const nf=IM.noiseFactorsStack(st,P,{method:"pca",nMax:4,full:true}); let e=0; REF.noise.forEach((o,i)=>e=Math.max(e,Math.abs(nf.table[i].chi2-o.chi2)/o.chi2));
  ok(e<1e-7,"χ²_n, n = 0…4",`max rel ${e.toExponential(1)}; χ²/dof ${nf.table.map(t=>t.ratio.toFixed(3)).join(" ")}`); }

console.log("== A5 per-pixel NNLS (constant background, unweighted and Poisson IRLS) vs scipy lsq_linear(bvls)");
for(const mode of ["none","poisson"]){ const r=IM.nnlsMap(st,[F64(REF.r1),F64(REF.r2)],{bg:1,weight:mode}); const q=REF.nnls[mode]; const np=REF.nx*REF.ny;
  let ec=0, es=0, er=0; for(let k=0;k<2;k++) for(let p=0;p<np;p++){ ec=Math.max(ec,Math.abs(r.coef[k][p]-q.coef[k*np+p])/Math.max(1,Math.abs(q.coef[k*np+p]))); const sr=q.sig[k*np+p]; if(sr!=null) es=Math.max(es,Math.abs(r.sig[k][p]-sr)/sr); }
  for(let p=0;p<np;p++) er=Math.max(er,Math.abs(r.redchi[p]-q.redchi[p])/q.redchi[p]);
  ok(ec<1e-8&&es<1e-8&&er<1e-8,`weight=${mode}`,`coef ${ec.toExponential(1)}, σ ${es.toExponential(1)}, χ²ν ${er.toExponential(1)}`); }

/* ================= B) synthetic truth ================= */
function synthStack(nx,ny,nE,seed,cf){ const G=rng(seed); const E=F64(Array.from({length:nE},(_,i)=>281+i*(10/nE))); const r1=gl(E,284.8,1.1,0.2), r2=gl(E,286.5,1.3,0.2);
  const A1=trapz(E,r1), A2=trapz(E,r2); const n1=F64(Array.from(r1,v=>v/A1)), n2=F64(Array.from(r2,v=>v/A2));     // unit-area references
  const data=new Float32Array(nx*ny*nE), a=new Float64Array(nx*ny), b=new Float64Array(nx*ny), bg=new Float64Array(nx*ny), mu=new Float64Array(nx*ny*nE);
  for(let y=0;y<ny;y++) for(let x=0;x<nx;x++){ const p=y*nx+x; a[p]=600*(0.2+0.8*x/(nx-1)); b[p]=500*(((x-nx/2)**2+(y-ny/2)**2)<(nx/4)**2?1:0.1); bg[p]=30+40*y/(ny-1);
    for(let e=0;e<nE;e++){ const m=a[p]*n1[e]+b[p]*n2[e]+bg[p]; mu[p*nE+e]=m; data[p*nE+e]=G.pois(m*cf)/cf; } }
  return {st:IM.makeStack({nx,ny,E,data,cfac:cf}),a,b,bg,mu,n1,n2,E}; }

console.log("== B1 chemical-state maps: Poisson-weighted NNLS coverage of ±1.96σ, bias");
{ const S=synthStack(48,40,80,3,1); const r=IM.nnlsMap(S.st,[S.n1,S.n2],{bg:1,weight:"poisson"}); const z=[[],[]]; let ins=[0,0];
  for(let p=0;p<S.a.length;p++){ [S.a[p],S.b[p]].forEach((t,k)=>{ const zz=(r.coef[k][p]-t)/r.sig[k][p]; z[k].push(zz); if(Math.abs(zz)<1.96) ins[k]++; }); }
  [0,1].forEach(k=>{ const s=stat(z[k]); ok(ins[k]/z[k].length>0.93&&ins[k]/z[k].length<0.97&&Math.abs(s.m)<0.1,`state ${k+1}`,`coverage ${(100*ins[k]/z[k].length).toFixed(1)} % of ${z[k].length}, mean z ${s.m.toFixed(3)}, sd z ${s.sd.toFixed(3)}`); });
  const rc=stat(Array.from(r.redchi)); ok(Math.abs(rc.m-1)<0.05,"mean reduced χ² ≈ 1 (Poisson model adequate)",rc.m.toFixed(3));
  // at% maps from the coefficients (S = 1, 1.5) vs truth
  const at=IM.atomicMaps([{label:"1",area:r.coef[0],sArea:r.sig[0],S:1},{label:"2",area:r.coef[1],sArea:r.sig[1],S:1.5}]); const zz=[]; let emax=0;
  for(let p=0;p<S.a.length;p++){ const n1=S.a[p], n2=S.b[p]/1.5; const t=100*n1/(n1+n2); if(isFinite(at.at[0][p])&&at.sAt[0][p]>0){ zz.push((at.at[0][p]-t)/at.sAt[0][p]); emax=Math.max(emax,Math.abs(at.at[0][p]-t)); } }
  const sz=stat(zz); ok(Math.abs(sz.m)<0.15&&sz.sd>0.85&&sz.sd<1.15,"at% map: z-scores of (at% − truth)/σ(at%)",`mean ${sz.m.toFixed(3)}, sd ${sz.sd.toFixed(3)}, max |error| ${emax.toFixed(2)} at.%`); }

console.log("== B2 number of factors (2 states + spatially varying background ⇒ 3) and PCA noise reduction");
{ const S=synthStack(64,64,80,8,1); const P=IM.pcaStack(S.st,{scale:"poisson",nKeep:8}); const nf=IM.noiseFactorsStack(S.st,P,{});
  ok(nf.nNoise===3,"noise criterion (MLPCA on 2048-pixel subset)",`n = ${nf.nNoise}; χ²/dof ${nf.table.map(t=>t.ratio.toFixed(3)).join(" ")}; IND → ${P.malinowski.nIND}`);
  const rec=IM.pcaReconstruct(S.st,P,3); let e0=0, e1=0; for(let i=0;i<S.mu.length;i++){ e0+=(S.st.data[i]-S.mu[i])**2; e1+=(rec.data[i]-S.mu[i])**2; }
  ok(Math.sqrt(e0/e1)>4,"PCA(3) reconstruction reduces RMS error vs truth",`raw ${Math.sqrt(e0/S.mu.length).toFixed(2)} → ${Math.sqrt(e1/S.mu.length).toFixed(2)} counts (×${Math.sqrt(e0/e1).toFixed(1)})`); }

console.log("== B3 VAMAS MAP file end-to-end (8×6 positions, C 1s, Shirley background)");
(async()=>{
  const IO=require("../src/xps_io.js"); const dd=path.join(__dirname,"data/synthetic"); const f=path.join(dd,"map_c1s.vms");
  if(fs.existsSync(f)){ const T=JSON.parse(fs.readFileSync(path.join(dd,"map_truth.json"),"utf8")); const rd=await IO.readFile("map_c1s.vms",new Uint8Array(fs.readFileSync(f)));
    const S=IM.fromSpectra(rd.spectra); ok(rd.vamas.mode==="MAP"&&rd.vamas.nX===8&&rd.vamas.nY===6&&S.nx===8&&S.ny===6&&S.missing===0,"MAP header (positions 48, 8×6) and stack assembly",`fov ${JSON.stringify(rd.spectra[0].meta.vamas.fov)}`);
    const sb=IM.subtractBackgroundStack(S,{lo:281,hi:293,bg:"shirley",n:5}); const ref=(c,w)=>{ const v=gl(sb.E,c,w,0); const a=trapz(sb.E,v); return F64(Array.from(v,x=>x/a)); };
    const r=IM.nnlsMap(sb,[ref(T.peaks.A[0],T.peaks.A[1]),ref(T.peaks.B[0],T.peaks.B[1])],{bg:0,weight:"none"});
    let eA=0, eB=0; for(let p=0;p<48;p++){ eA=Math.max(eA,Math.abs(r.coef[0][p]-T.areaA[p])/T.areaA[p]); eB=Math.max(eB,Math.abs(r.coef[1][p]-T.areaB[p])/T.areaB[p]); }
    ok(eA<0.08&&eB<0.15,"state maps from the file vs generator areas",`max rel error A ${(100*eA).toFixed(1)} %, B ${(100*eB).toFixed(1)} % (cps; peak B as low as 450)`);
  } else ok(false,"synthetic MAP VAMAS present (run test/ref/make_vamas.py)");

  console.log("== C performance: 256×256×100 (single thread; a Worker runs the same code)");
  { const S=synthStack(256,256,100,1,1); const T={}; const t=(k,f)=>{ const t0=Date.now(); const v=f(); T[k]=Date.now()-t0; return v; };
    t("window",()=>IM.windowImage(S.st,284,286)); t("area linear",()=>IM.areaMap(S.st,{lo:282,hi:289,bg:"linear"})); t("area shirley",()=>IM.areaMap(S.st,{lo:282,hi:289,bg:"shirley"}));
    t("roi",()=>IM.roiSpectrum(S.st,{x0:10,y0:10,x1:120,y1:200})); t("bin 2",()=>IM.binStack(S.st,2));
    const P=t("pca (poisson, 8)",()=>IM.pcaStack(S.st,{scale:"poisson",nKeep:8})); t("noise test (mlpca subset)",()=>IM.noiseFactorsStack(S.st,P,{}));
    t("pca reconstruct 3",()=>IM.pcaReconstruct(S.st,P,3)); t("nnls (none)",()=>IM.nnlsMap(S.st,[S.n1,S.n2],{bg:1})); t("nnls (poisson IRLS)",()=>IM.nnlsMap(S.st,[S.n1,S.n2],{bg:1,weight:"poisson"}));
    t("bg subtract shirley",()=>IM.subtractBackgroundStack(S.st,{lo:282,hi:289}));
    const worst=Math.max(...Object.values(T)); ok(worst<5000,"every operation < 5 s",Object.entries(T).map(([k,v])=>k+" "+v+" ms").join(", ")); }
  R.done();
})();
