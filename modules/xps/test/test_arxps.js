/* ARXPS engine tests — node test/test_arxps.js (from modules/xps)
   A) engine ↔ independent numpy/scipy reference (test/ref/ref_arxps.py → ref_arxps.json)
   B) synthetic-truth recovery: overlayer (incl. end-to-end from a synthetic VAMAS angle series), multilayer, RDP ordering, MaxEnt profiles + failure modes */
"use strict";
const path=require("path"), fs=require("fs");
const AR=require("../src/xps_arxps.js"), M=require("../src/xps_multi.js");
const {rng,reporter,maxRel,F64,gl,stat}=require("./synth_util.js");
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref/ref_arxps.json"),"utf8"));
const R=reporter("XPS_ARXPS"), ok=R.ok; const QUICK=!!process.env.QUICK;
const cosd=t=>Math.cos(t*Math.PI/180);

console.log("== A1 Seah–Dench attenuation lengths");
{ let e=0; REF.seahDench.forEach(o=>{ if(o.organic!=null) e=Math.max(e,Math.abs(AR.seahDench(o.E,{kind:"organic",rhoKgM3:o.rho})-o.organic)/o.organic);
  else { e=Math.max(e,Math.abs(AR.seahDench(o.E,{a:o.a})-o.element)/o.element,Math.abs(AR.seahDench(o.E,{a:o.a,kind:"inorganic"})-o.inorganic)/o.inorganic); } });
  ok(e<1e-14,"element / inorganic / organic forms",`max rel ${e.toExponential(1)}`);
  const a=AR.atomSize(60.08,2.2,3); ok(Math.abs(a-0.24728)<1e-4,"SiO2 mean atomic spacing a = (M/(ρ N_A n))^⅓ = 0.2473 nm (hand calculation)",a.toFixed(5)+" nm"); }

console.log("== A2 single-overlayer thickness (exact solution) vs scipy brentq");
{ let e=0; REF.overlayer.forEach(o=>{ const d=AR.overlayerThickness(o.R,o.th,{Lo:o.Lo,Ls:o.Ls,Rinf:o.Rinf}); e=Math.max(e,Math.abs(d-o.d)/o.d); }); ok(e<1e-10,"27 cases (L_o ≠ L_s and L_o = L_s closed form)",`max rel ${e.toExponential(1)}`); }

console.log("== A3 overlayer fit over angles vs scipy least_squares");
REF.fitOverlayer.forEach((f,i)=>{ const r=AR.fitOverlayer(f.th,f.R,f.sR,{Lo:f.Lo,Ls:f.Ls,Rinf:f.Rinf,fitRinf:f.fitRinf});
  const e=Math.max(Math.abs(r.d-f.d)/f.d,Math.abs(r.sd-f.sd)/f.sd,Math.abs(r.chi2-f.chi2)/f.chi2,f.fitRinf?Math.abs(r.Rinf-f.Rfit)/f.Rfit:0,f.fitRinf?Math.abs(r.sRinf-f.sRinf)/f.sRinf:0);
  ok(e<1e-7,`case ${i} (fit R∞: ${f.fitRinf})`,`d ${r.d.toFixed(5)} ± ${r.sd.toFixed(5)} nm vs ${f.d.toFixed(5)} ± ${f.sd.toFixed(5)}; max rel ${e.toExponential(1)}`); });

console.log("== A4 multilayer forward model and profiled fit vs reference");
{ const ML=REF.multilayer; let e=0; ML.th.forEach((t,a)=>{ const v=AR.multilayerI(ML.model,ML.d_true,t); e=Math.max(e,maxRel(v,ML.forward[a])); }); ok(e<1e-14,"forward intensities (3 species, 2 layers + substrate)",`max rel ${e.toExponential(1)}`);
  const r=AR.fitMultilayer(ML.model,ML.th,ML.I,ML.sI,{d0:[1,1]}); const ed=maxRel(r.d,ML.d), es=maxRel(r.sd,ML.sd), ec=Math.abs(r.chi2-ML.chi2)/ML.chi2;
  ok(ed<1e-6&&es<2e-3&&ec<1e-8,"fitted thicknesses, σ, χ² (per-angle scale profiled)",`d [${Array.from(r.d).map(v=>v.toFixed(4))}] vs [${ML.d.map(v=>v.toFixed(4))}] (rel ${ed.toExponential(1)}), σ rel ${es.toExponential(1)}, χ² rel ${ec.toExponential(1)}`); }

console.log("== A5 relative depth plot vs reference");
{ const ML=REF.multilayer; const r=AR.relativeDepth(ML.th,ML.I,ML.sI); let e=0; r.rows.forEach((q,k)=>{ e=Math.max(e,Math.abs(q.rdp-REF.rdp.rdp[k]),Math.abs(q.s-REF.rdp.s[k])); });
  ok(e<1e-12,"RDP and σ",`max abs ${e.toExponential(1)}`); ok(JSON.stringify(r.order)===JSON.stringify([0,1,2]),"ordering C 1s (contamination) > Si⁴⁺ (oxide) > Si⁰ (substrate)",r.rows.map(q=>q.rdp.toFixed(3)).join(", ")); }

console.log("== A6 MaxEnt: kernel, objective/gradient, fixed-α solutions, discrepancy α vs scipy L-BFGS-B");
{ const X=REF.maxent; const P=AR.maxentSetup(X.th,X.L,X.I,X.sI,{J:X.J});
  ok(maxRel(P.z,X.z,1)<1e-14&&maxRel(P.W[0][0],X.W00,1)<1e-14&&maxRel(Array.from({length:3},(_,k)=>P.m[k*X.J]),X.m,1)<1e-14,"layer grid, kernel W, default model");
  const O=REF.maxent_obj; const q=AR.meObjective(F64(O.u),P,O.alpha); ok(Math.abs(q.f-O.f)/Math.abs(O.f)<1e-12&&maxRel(q.grad,O.grad,Math.max(...O.grad.map(Math.abs)))<1e-11,"objective Q and analytic gradient",`Q ${q.f.toFixed(6)} vs ${O.f.toFixed(6)}`);
  X.solutions.forEach(s=>{ const r=AR.maxentSolve(P,s.alpha); let e=0; for(let k=0;k<3;k++) for(let j=0;j<X.J;j++) e=Math.max(e,Math.abs(r.c[k*X.J+j]-s.c[k][j]));
    ok(e<2e-4&&Math.abs(r.chi2-s.chi2)/s.chi2<1e-4,`α = ${s.alpha}: composition profile`,`max |Δc| ${e.toExponential(1)}, χ² ${r.chi2.toFixed(4)} vs ${s.chi2.toFixed(4)}`); });
  const r=AR.maxentProfile(X.th,X.L,X.I,X.sI,{J:X.J}); let e=0; for(let k=0;k<3;k++) for(let j=0;j<X.J;j++) e=Math.max(e,Math.abs(r.c[k][j]-X.cN[k][j]));
  ok(Math.abs(r.alpha-X.alphaN)/X.alphaN<2e-3&&e<1e-3,"discrepancy α (χ² = N) and its profile",`α ${r.alpha.toFixed(4)} vs ${X.alphaN.toFixed(4)}, χ² ${r.chi2.toFixed(3)} (N = ${r.N}), max |Δc| ${e.toExponential(1)}`); }

/* ================= B) synthetic truth ================= */
console.log("== B1 overlayer thickness: bias and ±1.96σ coverage over noisy realizations");
for(const cfg of [{d:0.5,Lo:3.4,Ls:3.1,Rinf:0.9},{d:1.8,Lo:3.4,Ls:3.1,Rinf:0.9},{d:5,Lo:2.6,Ls:2.9,Rinf:1.2}]){
  const th=[0,15,30,40,50,55,60,65]; const G=rng(Math.round(cfg.d*100)); const z=[], ds=[]; let inside=0; const nR=QUICK?60:300;
  for(let s=0;s<nR;s++){ // Poisson counts on both peaks, N0 counts for a thick (I∞) signal
    const N0=2e5; const Io=th.map(t=>G.pois(N0*cfg.Rinf*(1-Math.exp(-cfg.d/(cfg.Lo*cosd(t)))))), Is=th.map(t=>G.pois(N0*Math.exp(-cfg.d/(cfg.Ls*cosd(t)))));
    const Rr=th.map((t,a)=>Io[a]/Is[a]), sR=th.map((t,a)=>Rr[a]*Math.sqrt(1/Math.max(Io[a],1)+1/Math.max(Is[a],1)));
    const r=AR.fitOverlayer(th,Rr,sR,cfg); ds.push(r.d); const zz=(r.d-cfg.d)/r.sd; z.push(zz); if(Math.abs(zz)<1.96) inside++; }
  const sd=stat(ds), sz=stat(z); ok(Math.abs(sd.m-cfg.d)<3*sd.sd/Math.sqrt(nR)+1e-3*cfg.d&&Math.abs(inside/nR-0.95)<Math.max(0.03,3*Math.sqrt(0.95*0.05/nR)),`d = ${cfg.d} nm (L_o ${cfg.Lo}, L_s ${cfg.Ls})`,`mean ${sd.m.toFixed(4)} ± sd ${sd.sd.toFixed(4)} nm, coverage ${(100*inside/nR).toFixed(1)} %, sd(z) ${sz.sd.toFixed(3)}`); }

console.log("== B2 end-to-end: synthetic VAMAS angle series → reference fit of Si⁰/Si⁴⁺ → R(θ) → thickness");
(async()=>{
  const IO=require("../src/xps_io.js"); const dd=path.join(__dirname,"data/synthetic"); const f=path.join(dd,"arxps_angles.vms");
  if(fs.existsSync(f)){ const T=JSON.parse(fs.readFileSync(path.join(dd,"arxps_truth.json"),"utf8"));
    const rd=await IO.readFile("arxps_angles.vms",new Uint8Array(fs.readFileSync(f))); const ser=M.seriesFromSpectra(rd.spectra)[0];
    ok(ser.coordLabel==="Angle"&&ser.coord.length===T.angles.length&&ser.spectra.every((s,i)=>s.meta.vamas.takeOffPolar===T.angles[i]),"series by 'Angle' variable; take-off polar angle read per block",`${ser.coord.join(", ")} deg`);
    const g=M.commonGrid(ser.spectra); const sb=M.subtractBackground(g.E,g.Y,{lo:95,hi:108.5,bg:{type:"shirley",n:5}});
    const refs=[gl(sb.E,T.peaks.Si0[0],T.peaks.Si0[1],0),gl(sb.E,T.peaks.Si4[0],T.peaks.Si4[1],0)].map(v=>{ const a=v.reduce((s,x)=>s+x,0)*(sb.E[1]-sb.E[0]); return F64(Array.from(v,x=>x/a)); });
    const fits=M.fitReferences(sb.cols,refs,{E:sb.E,weight:"none",bg:0});
    const Rr=fits.map(q=>q.coef[1]/q.coef[0]); const sR=fits.map(q=>Math.sqrt(q.phi)*Math.hypot(q.sig[1]/q.coef[0],q.coef[1]*q.sig[0]/(q.coef[0]*q.coef[0])));
    const r=AR.fitOverlayer(T.angles,Rr,sR,{Lo:T.Lo,Ls:T.Ls,Rinf:T.Rinf});
    ok(Math.abs(r.d-T.d)<Math.max(3*r.sdScaled,0.03),"thickness from the file",`d = ${r.d.toFixed(4)} ± ${r.sdScaled.toFixed(4)} nm (true ${T.d}); χ²ν ${r.phi.toFixed(2)}; per-angle ${Array.from(r.perAngle).map(v=>v.toFixed(3)).join(" ")}`);
    const h=AR.hillPlot(T.angles,Rr,sR,{Lo:T.Lo,Rinf:T.Rinf}); ok(true,"(info) Hill plot slope·L_o (assumes L_o = L_s)",`${h.d.toFixed(4)} nm — biased because L_s ≠ L_o here (use the exact fit)`);
  } else ok(false,"synthetic ARXPS VAMAS present (run test/ref/make_vamas.py)");

  console.log("== B3 multilayer (contamination / oxide / substrate): bias and coverage with an unknown angle factor");
  { const ML=REF.multilayer; const th=ML.th; const dT=[0.6,1.5]; const G=rng(77); const z=[[],[]]; let inside=[0,0]; const nR=QUICK?30:150;
    for(let s=0;s<nR;s++){ const I=[], sI=[]; ML.model.species.forEach((sp,k)=>{ I.push([]); sI.push([]); });
      th.forEach((t,a)=>{ const v=AR.multilayerI(ML.model,dT,t); const Ga=1e5*(0.4+0.6*cosd(t)); v.forEach((x,k)=>{ const c=G.pois(x*Ga); I[k].push(c); sI[k].push(Math.sqrt(Math.max(c,1))); }); });
      const r=AR.fitMultilayer(ML.model,th,I,sI,{d0:[1,1]}); for(let j=0;j<2;j++){ const zz=(r.d[j]-dT[j])/r.sd[j]; z[j].push(zz); if(Math.abs(zz)<1.96) inside[j]++; } }
    for(let j=0;j<2;j++){ const st=stat(z[j]); ok(inside[j]/nR>0.9&&Math.abs(st.m)<0.25,`layer ${j} (true ${dT[j]} nm)`,`coverage ${(100*inside[j]/nR).toFixed(1)} %, mean z ${st.m.toFixed(3)}, sd z ${st.sd.toFixed(3)}`); } }

  console.log("== B4 MaxEnt recovery of simple profiles (1 % noise, discrepancy α) and failure modes");
  const th=[0,10,20,30,40,50,55,60,65,70];
  const synth=(L,cfun,G,noise)=>{ // cfun(z) → fractions array; fine-grid intensities, unknown angle factor
    const I=L.map(()=>[]), sI=L.map(()=>[]); const h=0.002;
    th.forEach(t=>{ const Ga=0.6+0.4*cosd(t); const acc=L.map(()=>0); for(let z=0;z<40;z+=h){ const c=cfun(z+h/2); L.forEach((Lk,k)=>{ acc[k]+=c[k]*Math.exp(-(z+h/2)/(Lk*cosd(t)))*h; }); }
      L.forEach((Lk,k)=>{ const v=acc[k]*Ga; const s=noise*v; I[k].push(v+s*G.g()); sI[k].push(s); }); });
    return {I,sI}; };
  const eqT=(r,k,zEnd)=>{ let s=0; for(let j=0;j<r.z.length-1;j++){ if(r.z[j]>=zEnd) break; s+=r.c[k][j]*(Math.min(r.z[j+1],zEnd)-r.z[j]); } return s; };
  { // (a) uniform 30/70 → stays uniform
    const G=rng(5); const L=[2.5,3.0]; const d=synth(L,()=>[0.3,0.7],G,0.01); const r=AR.maxentProfile(th,L,d.I,d.sI,{J:20}); let e=0; r.c[0].forEach(v=>e=Math.max(e,Math.abs(v-0.3)));
    ok(e<0.05,"uniform composition stays uniform",`max |c_A − 0.30| ${e.toFixed(3)} (α ${r.alpha.toExponential(2)}, χ² ${r.chi2.toFixed(1)}/${r.N})`); }
  { // (b) overlayer A of t on B: equivalent thickness of A, depth of the 50 % crossing
    for(const t of [0.8,1.5,3.0]){ const res=[], cr=[]; const nR=QUICK?3:10;
      for(let s=0;s<nR;s++){ const G=rng(100+s+t*10); const L=[2.5,3.0]; const d=synth(L,z=>z<t?[1,0]:[0,1],G,0.01); const r=AR.maxentProfile(th,L,d.I,d.sI,{J:24});
        res.push(eqT(r,0,r.z[r.z.length-1])); let zc=NaN; for(let j=0;j<r.z.length-1;j++){ const c0=r.c[0][j], c1=r.c[0][j+1]; if(c0>=0.5&&c1<0.5){ zc=r.z[j]+(c0-0.5)/(c0-c1)*(r.z[j+1]-r.z[j]); break; } } cr.push(zc); }
      const a=stat(res), b=stat(cr.filter(isFinite)); const tol=t<2?0.15*t+0.05:0.35*t;
      ok(Math.abs(a.m-t)<tol,`A overlayer t = ${t} nm on B (L 2.5/3.0 nm)`,`equivalent thickness ${a.m.toFixed(3)} ± ${a.sd.toFixed(3)} nm; 50 % crossing ${b.m.toFixed(3)} ± ${b.sd.toFixed(3)} nm (cells ${(r=>r)(24)} × ${(3*3.0/23).toFixed(2)} nm)`); } }
  { // (c) buried layer (truth from the Python reference): centroid and amount of B within 6 nm
    const X=REF.maxent; const r=AR.maxentProfile(X.th,X.L,X.I,X.sI,{J:X.J}); const zt=X.truth.z, cB=X.truth.cB; let am=0, cm=0; for(let i=1;i<zt.length;i++){ const h=zt[i]-zt[i-1]; am+=cB[i]*h; cm+=cB[i]*zt[i]*h; }
    let ar=0, crr=0; for(let j=0;j<r.z.length-1;j++){ const h=r.z[j+1]-r.z[j], zc=r.z[j]+h/2; ar+=r.c[1][j]*h; crr+=r.c[1][j]*zc*h; }
    ok(Math.abs(ar/am-1)<0.3&&Math.abs(crr/ar-cm/am)<0.6,"buried B layer (Gaussian, centre 2.2 nm, σ 0.42 nm) under an A-rich top",`amount ${ar.toFixed(3)} vs ${am.toFixed(3)} nm; centroid ${(crr/ar).toFixed(2)} vs ${(cm/am).toFixed(2)} nm`); }
  { // (d) failure mode: layer deeper than the information depth
    const G=rng(9); const L=[2.5,3.0]; const d=synth(L,z=>z>7&&z<8?[1,0]:[0,1],G,0.01); const r=AR.maxentProfile(th,L,d.I,d.sI,{J:20}); const tot=eqT(r,0,r.z[r.z.length-1]);
    ok(true,"(info) failure mode: 1 nm A layer buried at 7–8 nm (> 2.3 L)",`recovered A amount ${tot.toFixed(3)} nm (true 1.0); the profile is the smooth default — buried layers below ≈ 2–3 L cannot be located`); }
  { // (e) failure mode: underestimated σ → χ² = N not reachable
    const G=rng(11); const L=[2.5,3.0]; const d=synth(L,z=>z<1.5?[1,0]:[0,1],G,0.01); const sI=d.sI.map(r=>r.map(v=>v/20)); const r=AR.maxentProfile(th,L,d.I,sI,{J:20});
    ok(r.notReached||r.alpha<1e-3,"(guard) σ underestimated 20× → flagged 'noise level not reached' or α at the lower bound",`notReached ${r.notReached}, α ${r.alpha.toExponential(2)}, χ² ${r.chi2.toFixed(1)} vs N ${r.N}`); }
  R.done();
})();
