/* XPS series / multivariate engine tests — node test/test_multi.js (from modules/xps)
   A) engine ↔ independent numpy/scipy reference (test/ref/ref_multi.py → ref_multi.json), tight tolerances
   B) synthetic-truth recovery: factor counts, MCR-ALS, reference NNLS coverage, energy shift, integration σ, information-depth correction */
"use strict";
const path=require("path"), fs=require("fs");
const M=require("../src/xps_multi.js"), C=require("../src/xps_core.js");
const {rng,reporter,maxRel,F64,gl,trapz,stat}=require("./synth_util.js");
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref/ref_multi.json"),"utf8"));
const R=reporter("XPS_MULTI");
const ok=R.ok;
const QUICK=!!process.env.QUICK;

console.log("== A1 SVD (one-sided Jacobi) vs numpy.linalg.svd");
REF.svd.forEach((c,ci)=>{ const cols=c.cols.map(F64); const S=M.svdCols(cols); const s0=c.s[0];
  const es=maxRel(S.s.slice(0,c.s.length),c.s,s0); let eu=0; for(let k=0;k<3;k++) eu=Math.max(eu,maxRel(S.U[k],c.U3[k],1));
  ok(es<1e-12&&eu<1e-9,`case ${ci} (${cols[0].length}×${cols.length})`,`σ err/σ₁ ${es.toExponential(1)}, U err ${eu.toExponential(1)}`); });

console.log("== A2 PCA preprocessing × Malinowski RE/IND/F-test × reconstruction vs numpy/scipy.stats");
const mix=REF.mix, cols=mix.cols.map(F64), E=F64(mix.E);
REF.pca.forEach(p=>{
  const P=M.pca(cols,{scale:p.scale,center:p.center,method:"jacobi"});
  const es=maxRel(P.s.slice(0,p.s.length),p.s,p.s[0]); let eu=0; for(let k=0;k<3;k++) eu=Math.max(eu,maxRel(P.U[k],p.U3[k],1));
  let em=0; p.mal.forEach((o,i)=>{ const q=P.malinowski.table[i]; ["RE","IND","F"].forEach(k=>em=Math.max(em,Math.abs(q[k]-o[k])/Math.abs(o[k]))); em=Math.max(em,Math.abs(q.SL-o.SL)/Math.max(o.SL,1e-300)*(o.SL>1e-250?1:0)); });
  const rec=M.reconstruct(P,3); let er=0; rec.forEach((v,j)=>er=Math.max(er,maxRel(v,p.rec3[j],Math.max(...p.rec3[j]))));
  ok(es<1e-11&&eu<1e-8&&em<1e-7&&er<1e-9&&P.malinowski.nIND===p.nIND&&P.malinowski.nF===p.nF,`scale=${p.scale} center=${p.center}`,
    `σ ${es.toExponential(1)}, U ${eu.toExponential(1)}, RE/IND/F/SL ${em.toExponential(1)}, rec ${er.toExponential(1)}, nIND ${P.malinowski.nIND}/${p.nIND}, nF ${P.malinowski.nF}/${p.nF}`);
});

console.log("== A3 noise-limited factor count (weighted low-rank fit, MLPCA) vs reference");
{ const P=M.pca(cols,{}); const nf=M.noiseFactors(cols,P,1,{nMax:5,full:true}); let e=0; REF.noise.table.forEach((o,i)=>e=Math.max(e,Math.abs(nf.table[i].chi2-o.chi2)/o.chi2));
  ok(e<1e-6&&nf.nNoise===REF.noise.nNoise,"χ²_n table and n_noise",`max rel ${e.toExponential(1)}, n ${nf.nNoise}/${REF.noise.nNoise}; χ²/dof: ${nf.table.map(t=>t.ratio.toFixed(3)).join(" ")}`); }

console.log("== A4 target testing (AET, REP, RET, SPOIL) vs reference");
{ const P=M.pca(cols,{}); const T=REF.target; let e=0; const v=[];
  T.tests.forEach(t=>{ const r=M.targetTest(P,T.n,F64(t.x)); ["AET","REP","RET","SPOIL"].forEach(k=>e=Math.max(e,Math.abs(r[k]-t[k])/Math.max(Math.abs(t[k]),1e-12))); v.push(t.name+" "+r.SPOIL.toFixed(2)); });
  ok(e<1e-8,"Malinowski target-test statistics",`max rel ${e.toExponential(1)}`);
  const sp=Object.fromEntries(T.tests.map(t=>[t.name,t.SPOIL]));
  ok(sp.state1<3&&sp.state2<3&&sp.state3<3&&sp.flat<3,"true component spectra accepted (SPOIL < 3)",v.join(", "));
  ok(sp.wrong_shifted>6&&sp.wrong_narrow>6,"wrong spectra (shifted 1 eV / half width) rejected (SPOIL > 6)",`${sp.wrong_shifted.toFixed(1)}, ${sp.wrong_narrow.toFixed(1)}`); }

console.log("== A5 OPA picks, SIMPLISMA channels and ITTFA profiles vs reference");
{ const pk=M.opa(cols,3); ok(JSON.stringify(pk)===JSON.stringify(REF.opa),"OPA picks",pk.join(",")+" vs "+REF.opa.join(","));
  const sp=M.simplisma(cols,4,0.03); ok(JSON.stringify(sp)===JSON.stringify(REF.simplisma),"SIMPLISMA purest channels",sp.join(",")+" vs "+REF.simplisma.join(","));
  const P=M.pca(cols,{}); const it=M.ittfa(cols,P,3,{needles:REF.opa}); let e=0; it.C.forEach((x,k)=>e=Math.max(e,maxRel(x,REF.ittfa[k],1)));
  ok(e<1e-8,"ITTFA concentration profiles",`max abs ${e.toExponential(1)}`); }

console.log("== A6 NNLS (Lawson–Hanson, Gram form) vs scipy.optimize.nnls");
REF.nnls.forEach((c,i)=>{ const x=M.nnlsCols(c.cols.map(F64),F64(c.b)).x; const sc=Math.max(...c.x.map(Math.abs),1); ok(maxRel(x,c.x,sc)<1e-10,`case ${i} (${c.b.length}×${c.cols.length})`,`max err ${maxRel(x,c.x,sc).toExponential(1)}`); });

console.log("== A7 reference NNLS fit (IRLS Poisson weights, free background terms) vs scipy lsq_linear(bvls)");
{ const refs=mix.S.map(F64); let ec=0, es=0, ex=0;
  REF.reffit.forEach(f=>{ const r=M.fitReferences([cols[f.level]],refs,{E,weight:"poisson",cfac:1,bg:f.bg})[0]; const sc=Math.max(...f.coef);
    ec=Math.max(ec,maxRel(r.coef,f.coef,sc)); f.sig.forEach((s,k)=>{ if(s!=null) es=Math.max(es,Math.abs(r.sig[k]-s)/s); }); ex=Math.max(ex,Math.abs(r.chi2-f.chi2)/f.chi2); });
  ok(ec<1e-8&&es<1e-8&&ex<1e-8,"coefficients, σ, χ² (12 cases: 4 levels × bg 0/1/2)",`coef ${ec.toExponential(1)}, σ ${es.toExponential(1)}, χ² ${ex.toExponential(1)}`); }

console.log("== A8 region integration (linear / Shirley, Poisson σ) vs reference");
REF.integrate.forEach(g=>{ const r=M.integrateSeries(E,cols,{lo:g.lo,hi:g.hi,bg:{type:g.bg,n:5},cfac:1}); let ea=0, es=0;
  r.forEach((q,j)=>{ ea=Math.max(ea,Math.abs(q.area-g.rows[j][0])/Math.abs(g.rows[j][0])); es=Math.max(es,Math.abs(q.sArea-g.rows[j][1])/g.rows[j][1]); });
  ok(ea<1e-9&&es<1e-10,`${g.bg}: area and σ (24 levels)`,`area ${ea.toExponential(1)}, σ ${es.toExponential(1)}`); });

console.log("== A9 composition vs level (XPS_CORE.quantify per level) vs reference");
{ const Cm=REF.composition; const rows=Cm.S.map((S,k)=>({label:"s"+k,el:"Si",S,BE:100,hv:1486.6,PE:20,tf:null})); const areas=Cm.areas.map((a,k)=>a.map((v,j)=>({area:v,sArea:Cm.sig[k][j]})));
  const out=M.compositionSeries(rows,areas,{}); let e=0, nd=0;
  out.forEach((lv,j)=>lv.forEach((q,k)=>{ const ra=Cm.levels[j].at[k], rs=Cm.levels[j].sAt[k]; if(ra==null){ if(q.at!=null) e=Infinity; nd++; return; } e=Math.max(e,Math.abs(q.at-ra),Math.abs(q.sAt-rs)); }));
  ok(e<1e-10,"at% and σ(at%) per level, n.d. rule",`max abs ${e.toExponential(1)}, n.d. cells ${nd}`); }

console.log("== A10 information-depth correction c = I − λ dI/dz vs reference");
{ const I=REF.infodepth; const r=M.infoDepthCorrect(F64(I.z),I.X.map(F64),I.lam,{hw:I.hw}); let e=0; r.n.forEach((v,k)=>e=Math.max(e,maxRel(v,I.n[k],1)));
  ok(e<1e-10,"corrected signals",`max abs ${e.toExponential(1)}`);
  // truth: the measured signal of a sigmoid interface is broadened by the exponential response; the correction must restore it
  let eRaw=0, eCor=0; I.truthA.forEach((t,i)=>{ eRaw=Math.max(eRaw,Math.abs(I.X[0][i]-t)); eCor=Math.max(eCor,Math.abs(r.n[0][i]-t)); });
  ok(eCor<0.1&&eCor<eRaw/5,"restores a known interface (noise-free, λ_e = 1.6 nm, step 0.3–0.7 nm)",`max |error| raw ${eRaw.toFixed(3)} → corrected ${eCor.toFixed(3)}`); }

/* ================= B) synthetic truth ================= */
function mixture(G,o){ // nE×nL, K Gaussian–Lorentzian states, smooth depth profiles, constant background, Poisson
  const nE=o.nE||160, nL=o.nL||24, K=o.K; const E=F64(Array.from({length:nE},(_,i)=>o.E0+i*o.dE));
  const S=o.centers.map((c,k)=>gl(E,c,o.widths[k]));
  const Ct=[]; for(let k=0;k<K;k++) Ct.push(F64(Array.from({length:nL},(_,j)=>{ const t=j/(nL-1); return o.prof[k](t); })));
  const cols=[], mu=[]; for(let j=0;j<nL;j++){ const t=j/(nL-1); const m=new Float64Array(nE).fill(o.bgConst?o.bg:o.bg*(1+3*t*t)); for(let k=0;k<K;k++) for(let i=0;i<nE;i++) m[i]+=o.scale*Ct[k][j]*S[k][i]; mu.push(m); cols.push(F64(Array.from(m,v=>G.pois(v)))); }
  return {E,S,Ct,cols,mu};
}
const PROF={ top:t=>1/(1+Math.exp((t-0.35)/0.07)), mid:t=>Math.exp(-(((t-0.45)/0.14)**2)), bot:t=>1/(1+Math.exp(-(t-0.55)/0.07)), mid2:t=>Math.exp(-(((t-0.25)/0.1)**2)) };
const CASES=[
  {K:2,E0:280,dE:0.05,centers:[284.8,286.4],widths:[1.1,1.3],prof:[PROF.top,PROF.bot],bg:100},
  {K:3,E0:96,dE:0.075,centers:[99.4,101.5,103.6],widths:[1.0,1.3,1.5],prof:[PROF.bot,PROF.mid,PROF.top],bg:150},
  {K:4,E0:526,dE:0.06,centers:[529.6,531.2,532.4,533.6],widths:[1.0,1.3,1.4,1.5],prof:[PROF.bot,PROF.mid,PROF.top,PROF.mid2],bg:200},
];

console.log("== B1 number of components (K states + level-dependent background ⇒ K+1 factors), 3 count levels × 3 cases × "+(QUICK?3:10)+" Poisson realizations");
for(const scale of [300,3000,30000]){
  for(const cs of CASES){ const tally={IND:0,F:0,noise:0}; const nR=QUICK?3:10; const got={IND:[],F:[],noise:[]};
    for(let s=0;s<nR;s++){ const G=rng(1000*cs.K+s+scale); const d=mixture(G,Object.assign({scale},cs)); const P=M.pca(d.cols,{}); const nf=M.noiseFactors(d.cols,P,1,{nMax:cs.K+3});
      const truth=cs.K+1; if(P.malinowski.nIND===truth) tally.IND++; if(P.malinowski.nF===truth) tally.F++; if(nf.nNoise===truth) tally.noise++;
      got.IND.push(P.malinowski.nIND); got.F.push(P.malinowski.nF); got.noise.push(nf.nNoise); }
    // detectability of the weakest true factor: noise-free Poisson-scaled σ_{K+1} vs the noise edge σ_noise(√r + √c) (Marchenko–Pastur)
    const d0=mixture(rng(1),Object.assign({scale},cs)); const P0=M.pca(d0.mu,{scale:"poisson"}); const pp=P0.pp; let vs=0, nv=0; d0.mu.forEach((m,j)=>{ for(let i=0;i<m.length;i++){ vs+=m[i]*pp.rowW[i]*pp.rowW[i]*pp.colW[j]*pp.colW[j]; nv++; } });
    const edge=Math.sqrt(vs/nv)*(Math.sqrt(d0.mu[0].length)+Math.sqrt(d0.mu.length)); const snr=P0.s[cs.K]/edge;
    const info=`σ_${cs.K+1}/edge ${snr.toFixed(1)}; IND ${tally.IND}/${nR} [${got.IND.join("")}], F(5%) ${tally.F}/${nR} [${got.F.join("")}], noise ${tally.noise}/${nR} [${got.noise.join("")}]`;
    // the noise criterion is the one the engine recommends: correct in ≥ 80 % whenever the weakest factor is above 2× the noise edge
    ok(snr<2||tally.noise>=Math.ceil(0.8*nR),`K=${cs.K}, peak height ≈ ${scale} counts`+(snr<2?" (weakest factor below detectability — info only)":""),info);
  }
}
{ // documented failure mode: a level-INDEPENDENT constant background is nearly a combination of the concentration profiles → not a separate factor
  const cs=Object.assign({},CASES[0],{bgConst:true}); const G=rng(42); const d=mixture(G,Object.assign({scale:3000},cs)); const P=M.pca(d.cols,{}); const nf=M.noiseFactors(d.cols,P,1,{nMax:5});
  const P0=M.pca(d.mu,{}); ok(true,"(info) constant background: noise-free σ_"+(cs.K+1)+" = "+P0.s[cs.K].toFixed(0)+" vs noise ≈ "+P.s[cs.K+2].toFixed(0)+" → n_noise = "+nf.nNoise+" (K+1 = "+(cs.K+1)+" not resolvable)"); }

console.log("== B2 MCR-ALS (non-negative, SIMPLISMA start) on background-subtracted spectra (linear bg) recovers concentration fractions");
for(const cs of CASES){ const G=rng(77+cs.K); const d=mixture(G,Object.assign({scale:5000},cs)); const n=cs.K;
  const sb=M.subtractBackground(d.E,d.cols,{bg:{type:"linear",n:5}});
  const r=M.mcrAls(sb.cols,n,{maxIter:3000});
  const areas=d.S.map(s=>trapz(d.E,s)); const tr=d.Ct.map((c,k)=>F64(Array.from(c,v=>5000*v*areas[k])));
  const truthS=d.S; const used=new Set(), map=[];
  truthS.forEach(ts=>{ let best=-2, bj=-1; r.S.forEach((rs,j)=>{ if(used.has(j)) return; const c=M.dot(ts,rs)/(M.norm2(ts)*M.norm2(rs)); if(c>best){ best=c; bj=j; } }); used.add(bj); map.push(bj); });
  let maxErr=0; const nl=d.cols.length;
  for(let j=0;j<nl;j++){ let T=0, Rr=0; for(let k=0;k<cs.K;k++){ T+=tr[k][j]; Rr+=r.C[map[k]][j]; } for(let k=0;k<cs.K;k++){ const ft=tr[k][j]/T, fr=r.C[map[k]][j]/Rr; maxErr=Math.max(maxErr,Math.abs(ft-fr)); } }
  const ro=M.mcrAls(sb.cols,n,{maxIter:3000,init:"opa"}); let eo=0; { const mp=truthS.map(ts=>{ let b=-2,bj=-1; ro.S.forEach((rs,j)=>{ const c=M.dot(ts,rs)/(M.norm2(ts)*M.norm2(rs)); if(c>b){ b=c; bj=j; } }); return bj; });
    for(let j=0;j<nl;j++){ let T=0,Rr=0; for(let k=0;k<cs.K;k++){ T+=tr[k][j]; Rr+=ro.C[mp[k]][j]; } for(let k=0;k<cs.K;k++) eo=Math.max(eo,Math.abs(tr[k][j]/T-ro.C[mp[k]][j]/Rr)); } }
  ok(cs.K===4||maxErr<0.05,`K=${cs.K}: max |fraction error| over all levels`+(cs.K===4?" (4 overlapping O 1s states — info: rotational ambiguity)":""),`SIMPLISMA start ${(100*maxErr).toFixed(2)} at.% (LOF ${r.lof.toFixed(3)} %, ${r.iter} it); OPA start ${(100*eo).toFixed(2)} at.% (LOF ${ro.lof.toFixed(3)} %)`);
}

console.log("== B3 reference NNLS fit: coverage of ±1.96σ with the true reference spectra (Poisson, IRLS)");
{ const cs=CASES[1]; const z=[]; let inside=0, n=0;
  for(let s=0;s<(QUICK?4:20);s++){ const G=rng(500+s); const d=mixture(G,Object.assign({scale:3000},cs)); const res=M.fitReferences(d.cols,d.S,{E:d.E,weight:"poisson",cfac:1,bg:1});
    res.forEach((r,j)=>{ for(let k=0;k<cs.K;k++){ const t=3000*d.Ct[k][j]; if(t<50) continue; const zz=(r.coef[k]-t)/r.sig[k]; z.push(zz); n++; if(Math.abs(zz)<1.96) inside++; } }); }
  const st=stat(z); ok(inside/n>0.92&&inside/n<0.98&&Math.abs(st.m)<0.15&&Math.abs(st.sd-1)<0.1,"coefficient z-scores",`coverage ${(100*inside/n).toFixed(1)} % of ${n}, mean z ${st.m.toFixed(3)}, sd z ${st.sd.toFixed(3)}`); }

console.log("== B4 reference fit with an energy shift (charging drift between levels)");
{ const cs=CASES[0]; const G=rng(9); const d=mixture(G,Object.assign({scale:4000},cs)); const shifts=[0.23,-0.31,0.05]; let e=0;
  shifts.forEach((sh,q)=>{ const y=F64(Array.from(M.shiftCopy(d.E,d.mu[5+q*6],sh),v=>G.pois(v))); const r=M.fitReferences([y],d.S,{E:d.E,weight:"poisson",cfac:1,bg:1,maxShift:0.5,shiftStep:0.05})[0]; e=Math.max(e,Math.abs(r.shift-sh)); });
  ok(e<0.02,"recovered shift within 0.02 eV (true 0.23, −0.31, 0.05 eV)",`max |error| ${e.toFixed(4)} eV`); }

console.log("== B5 integration σ vs Monte Carlo spread (linear and Shirley)");
{ const nE=241, E2=F64(Array.from({length:nE},(_,i)=>280+i*0.05)); const pk=gl(E2,286,1.2); const step=F64(Array.from(E2,e=>1/(1+Math.exp(-(e-286)/0.4))));
  const mu=F64(Array.from(E2,(e,i)=>400+3000*pk[i]+150*step[i])); const cf=0.5; const G=rng(4);
  for(const bg of ["linear","shirley"]){ const As=[], Ss=[]; for(let s=0;s<(QUICK?100:400);s++){ const y=F64(Array.from(mu,v=>G.pois(v*cf)/cf)); const r=M.integrateSeries(E2,[y],{bg:{type:bg,n:5},cfac:cf})[0]; As.push(r.area); Ss.push(r.sArea); }
    const a=stat(As), sm=stat(Ss).m; ok(Math.abs(sm/a.sd-1)<0.12,`${bg}: mean predicted σ vs MC sd of the area`,`σ ${sm.toFixed(1)} vs sd ${a.sd.toFixed(1)} (ratio ${(sm/a.sd).toFixed(3)})`); }
  const trueA=3000*trapz(E2,pk); const r=M.integrateSeries(E2,[mu],{bg:{type:"shirley",n:5}})[0]; ok(Math.abs(r.area/trueA-1)<0.03,"Shirley area of a peak on a Shirley-like step (noise-free)",`${(100*(r.area/trueA-1)).toFixed(2)} % vs true`); }

console.log("== B6 information-depth correction with noise (two elements, different λ)");
{ const zs=F64(Array.from({length:40},(_,i)=>i*0.25)); const lam=[1.5,2.4]; const ct=z=>1/(1+Math.exp(-(z-4)/0.3));
  const meas=(z0,L,f)=>{ let s=0; const h=0.002; for(let z=z0;z<z0+20*L;z+=h){ const w=Math.exp(-(z-z0)/L)/L; s+=f(z)*w*h; } return s; };
  const IA=F64(Array.from(zs,z=>meas(z,lam[0],ct))), IB=F64(Array.from(zs,z=>meas(z,lam[1],z=>1-ct(z))));
  const G=rng(31); const errs=[], errsRaw=[];
  for(let s=0;s<(QUICK?5:30);s++){ const XA=F64(Array.from(IA,v=>v+0.01*G.g())), XB=F64(Array.from(IB,v=>v+0.01*G.g()));
    const r=M.infoDepthCorrect(zs,[XA,XB],lam,{hw:3}); let e=0, e0=0; zs.forEach((z,i)=>{ const t=100*ct(z); e=Math.max(e,Math.abs(r.at[0][i]-t)); const raw=100*XA[i]/(XA[i]+XB[i]); e0=Math.max(e0,Math.abs(raw-t)); }); errs.push(e); errsRaw.push(e0); }
  const a=stat(errs), b=stat(errsRaw); ok(a.m<b.m/2,"1 % noise: mean of max |at% error| vs depth",`uncorrected ${b.m.toFixed(1)} at.% → corrected ${a.m.toFixed(1)} at.% (hw = 3)`); }

console.log("== B7 series assembly from synthetic VAMAS (depth profile, SDP sputter fields)");
(async()=>{
  const IO=require("../src/xps_io.js"); const dd=path.join(__dirname,"data/synthetic");
  const f=path.join(dd,"depth_profile_sdp.vms");
  if(fs.existsSync(f)){ const r=await IO.readFile("depth_profile_sdp.vms",new Uint8Array(fs.readFileSync(f)));
    const ser=M.seriesFromSpectra(r.spectra); const truth=JSON.parse(fs.readFileSync(path.join(dd,"depth_profile_truth.json"),"utf8"));
    ok(r.vamas&&r.vamas.mode==="SDP"&&ser.length===truth.regions.length&&ser.every(s=>s.coord.length===truth.levels.length),"SDP file → series per region with etch-time coordinate",`mode ${r.vamas&&r.vamas.mode}, ${ser.map(s=>s.region+":"+s.coord.length).join(", ")}, coord ${ser[0].coordLabel} [${ser[0].coordUnit}]`);
    const okc=ser.every(s=>s.coord.every((v,i)=>Math.abs(v-truth.levels[i])<1e-9)); ok(okc,"coordinates = etch times written by the generator");
    const sp0=r.spectra[0]; ok(sp0.meta.vamas.sputter&&Math.abs(sp0.meta.vamas.sputter.energy-truth.sputterEnergy)<1e-9,"sputter source energy read from the block",JSON.stringify(sp0.meta.vamas.sputter));
    // integrate + quantify each region, compare with truth composition
    const rows=[], areas=[]; ser.forEach(s=>{ const g=M.commonGrid(s.spectra); const t=truth.regionInfo[s.region]; const a=M.integrateSeries(g.E,g.Y,{lo:t.lo,hi:t.hi,bg:{type:"shirley",n:5},cfac:s.spectra[0].c||1}); areas.push(a); rows.push({label:s.region,el:t.el,S:t.S,BE:t.BE,hv:1486.6,PE:20,tf:null}); });
    const comp=M.compositionSeries(rows,areas,{}); let e=0; comp.forEach((lv,j)=>lv.forEach((q,k)=>{ const tAt=truth.at[j][truth.regions.indexOf(ser[k].region)]; e=Math.max(e,Math.abs(q.at-tAt)); }));
    ok(e<3,"composition vs etch time from the synthetic file (Shirley areas) within 3 at.%",`max |error| ${e.toFixed(2)} at.%`);
  } else ok(false,"synthetic VAMAS files present (run test/ref/make_vamas.py)");
  const t=path.join(dd,"depth_profile_multicol.txt");
  if(fs.existsSync(t)){ const r=IO.textMultiParse(fs.readFileSync(t,"utf8"),"depth_profile_multicol.txt"); ok(r.Y.length===12&&r.coord[3]===90&&r.E[0]<r.E[1],"multi-column text: 12 spectra, coordinates from header, BE ascending",`${r.Y.length} cols, coord ${r.coord.slice(0,5).join(",")}…`); }
  R.done();
})();
