/* Real-data case S3 (CeO₂ sample, K-Alpha VGD, PE 100) for the fit guidance (SPEC §17.5)
   node test/test_real_ce.js            (from modules/xps, about 15 s)
   Data: test/data/real/ce/ (private, not in the public source backup): "C1s Scan.VGD", "Ce3d Scan.VGD", "O1s Scan.VGD", "Si2p Scan.VGD", "Survey.VGD".
   If the folder is absent every item is skipped and the script exits 0.
   The data are a validation case only — nothing in src/ is tuned to them.
   1 charge correction from C 1s (maximum → 284.8 eV)
   2 naive fit = what Auto suggest + Fit does (suggestModel + refineByResidual, Smart background, whole Ce 3d spectrum) → fit check
   3 recipe (Mixed, 10 lines) + staged fit → Ce³⁺ fraction, u‴ fraction, χ²ν, bound hits, fit check; CeO₂-only recipe for comparison */
"use strict";
const path=require("path"), fs=require("fs");
const ROOT=path.join(__dirname,"data","real","ce");
if(!fs.existsSync(path.join(ROOT,"Ce3d Scan.VGD"))){ console.log("test/data/real/ce absent — real-data case S3 skipped (skip)"); process.exit(0); }
const IO=require("../src/xps_io.js"), C=require("../src/xps_core.js"), G=require("../src/xps_guide.js"), L=require("../src/xps_lib.js");
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info!=null&&info!==""?"  ("+info+")":"")); }
const info=(name,v)=>console.log("  info "+name+(v!=null?"  ("+v+")":""));
const rd=async f=>(await IO.readFile(f,new Uint8Array(fs.readFileSync(path.join(ROOT,f))))).spectra[0];
const nb=fin=>fin.cov.atB.filter((b,j)=>b&&fin.free[j].k!=="area").length;
const show=d=>d.map(x=>`[${x.sev}] ${x.label}`).join(" | ");
(async()=>{
  const t0=Date.now();
  const c1=await rd("C1s Scan.VGD"), ce=await rd("Ce3d Scan.VGD"), o1=await rd("O1s Scan.VGD");
  console.log("== 1 charge correction");
  const pk=C.peakMaxPos(c1.be,c1.y,280,292,0.25); const sh=+(284.8-pk.be).toFixed(3);
  ok(ce.reg.el==="Ce"&&ce.reg.orb==="3d"&&ce.meta.passEnergy===100,"Ce 3d region read (PE 100)",`${ce.n} points ${ce.be[0]}–${ce.be[ce.n-1]} eV, c = ${ce.c}`);
  ok(pk.be>283&&pk.be<284,"C 1s maximum (raw) → shift",`${pk.be.toFixed(2)} eV → Δ = +${sh.toFixed(2)} eV`);
  const spec={be:ce.be,y:ce.y,c:ce.c}, sp={be0:ce.be[0],be1:ce.be[ce.n-1]};
  console.log("== 2 naive: Auto suggest + Fit (whole spectrum, Smart)");
  const lo=ce.be[0], hi=ce.be[ce.n-1]; const R0=C.sliceRegion(ce.be,ce.y,lo,hi);
  const sg=C.suggestModel(Float64Array.from(R0.E),Float64Array.from(R0.y),{bg:{type:"smart",n:5},el:"Ce",orb:"3d",lib:L});
  let PN=C.makeProblem(spec,{lo,hi,bg:{type:"smart",n:5},comps:sg.comps,cons:sg.cons}); const rf=C.refineByResidual(PN,{});
  PN=C.makeProblem(spec,{lo,hi,bg:{type:"smart",n:5},comps:rf.prob.comps,cons:rf.prob.cons}); const rN=C.fitOnce(PN,{}); const fN=C.finalize(PN,rN);
  const dN=G.diagnose({E:PN.E,y:PN.y,fin:fN,el:"Ce",orb:"3d",shift:sh,bgType:"smart",spec:sp});
  info("naive",`${PN.comps.length} lines · χ²ν ${fN.stats.chi2nu.toFixed(2)} · bound hits ${nb(fN)} · ${fN.stop}`); info("naive fit check",show(dN));
  ok(dN.some(x=>x.id==="complex"&&x.act&&x.act.type==="recipe")&&dN.some(x=>x.sev==="warn"),"fit check on the naive fit: warnings + pointer to the Ce 3d recipe");
  console.log("== 3 recipe + staged fit");
  const out={};
  for(const opt of ["mix","ce4"]){ const t1=Date.now();
    const M=G.recipeModel("ce3d",{shift:sh,be:ce.be,y:ce.y,opt}); const P=C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons});
    const s=G.stagedFit(P,{}); const fin=C.finalize(s.P,s.r); const dv=G.recipeDerived(fin,{id:"ce3d"});
    const d=G.diagnose({E:s.P.E,y:s.P.y,fin,el:"Ce",orb:"3d",shift:sh,recipe:G.recipeById("ce3d"),bgType:M.bg.type,spec:sp});
    const v3=s.P.comps.find(c=>c.rk==="v3"); const u3be=v3.be+v3.doublet.dE;
    out[opt]={fin,dv,d,align:M.recipe.align,u3be,staged:s.staged};
    info(`recipe ${M.recipe.label}`,`alignment ${M.recipe.align>=0?"+":""}${M.recipe.align.toFixed(2)} eV · Ce³⁺ ${dv[0].v==null?"—":dv[0].v.toFixed(1)+" ± "+dv[0].sig.toFixed(1)} % · u‴ ${dv[1].v.toFixed(2)} ± ${dv[1].sig.toFixed(2)} % · χ²ν ${fin.stats.chi2nu.toFixed(2)} · bound hits ${nb(fin)} · kept ${s.staged.used} (χ²_d staged ${s.staged.chi2dStaged.toFixed(0)} / normal ${s.staged.chi2dNormal.toFixed(0)}) · u‴ raw ${u3be.toFixed(2)} eV · ${Date.now()-t1} ms`);
    info(`fit check (${opt})`,show(d));
    ok(s.staged.chi2dStaged>=0&&fin.stats.chi2d<=s.staged.chi2dNormal*(1+1e-12),`${opt}: staged result not worse than the normal fit (χ²_d)`); }
  const m=out.mix;
  ok(m.dv.every(x=>x.v!=null&&isFinite(x.v)&&x.sig!=null),"Mixed recipe: Ce³⁺ and u‴ fractions with σ",m.dv.map(x=>`${x.label} ${x.v.toFixed(1)}`).join(" · "));
  ok(nb(m.fin)<=nb(fN),"bound hits recipe ≤ naive",`${nb(m.fin)} vs ${nb(fN)}`);
  ok(Math.abs(m.align)>0.5,"alignment flags a charge-reference mismatch: Ce 3d lines sit > 0.5 eV from the C 1s-corrected literature positions",`${m.align.toFixed(2)} eV; u‴ raw ${m.u3be.toFixed(2)} eV vs literature 916.7–917.0 eV`);
  ok(m.dv[1].v>8&&m.dv[1].v<15,"u‴ fraction in the Ce⁴⁺-dominated range (8–15 %)",`${m.dv[1].v.toFixed(2)} %`);
  // O 1s: low-BE lattice-oxide shoulder raw position (CeO₂ lattice O 1s ≈ 529.7 eV, Morgan 2023) — independent check of the charging offset
  const R1=C.sliceRegion(o1.be,o1.y,527,531); const ys=C.smoothQ(Float64Array.from(R1.E),Float64Array.from(R1.y),3); let best=null;
  for(let i=2;i<ys.length-2;i++){ const d2=ys[i+1]-2*ys[i]+ys[i-1]; if(!best||d2<best.d2) best={be:R1.E[i],d2}; }
  info("O 1s lattice shoulder (raw, max curvature)",`${best.be.toFixed(2)} eV`);
  console.log(`\ntest_real_ce: ${pass} passed, ${fail} failed (${((Date.now()-t0)/1000).toFixed(0)} s)`);
  if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; }
})().catch(e=>{ console.error(e); process.exitCode=1; });
