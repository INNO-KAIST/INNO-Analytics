/* Real-data comparison for Auto build (SPEC §18.6) — numbers only, samples anonymised (K1/K2: CasaXPS-curated Ba 4d · C 1s · O 1s;
   A1/A2: K-Alpha VGD mixed-metal oxide (C 1s · O 1s · Ba 4d), no curated deconvolution — compared with the manual session described in SPEC §18;
   S3: CeO₂ (Ce 3d recipe as in test_real_ce)).
   node test/test_real_autobuild.js     (from modules/xps, about 2–4 min). Data: test/data/real/{kara,a2,ce} (private) — absent → SKIPPED.
   For every spectrum: Auto build from Auto suggest and (where a recipe applies) from the recipe; reported: lines (raw BE), established or the reasons,
   groups, χ²ν, and the agreement with the curated model (same line count; every curated line within 0.3 eV of an Auto build line).
   Checks are invariants only (runs complete, input unchanged, caps hold in the final model, a not-established result always carries reasons);
   the agreement itself is reported, not asserted — the data are a validation case and nothing in src/ is tuned to them. */
"use strict";
const path=require("path"), fs=require("fs");
const R0=path.join(__dirname,"data","real");
const A_DIRS=(()=>{ try{ return fs.readdirSync(path.join(R0,"a2")).sort().map(d=>"a2/"+d+"/Group").filter(d=>fs.existsSync(path.join(R0,d,"O1s Scan.VGD"))); }catch(e){ return []; } })();   // 시료 폴더 이름은 저장소에 적지 않음 (비공개 자료)
const PRES=(()=>{ try{ return JSON.parse(fs.readFileSync(path.join(R0,"present.json"),"utf8")); }catch(e){ return ["O","C"]; } })();   // 시료 원소 목록: 비공개 data/real/present.json (gitignore)
const have=["kara/casa_240213/Integrated.vms",...A_DIRS.map(d=>d+"/O1s Scan.VGD"),"ce/Ce3d Scan.VGD"].filter(f=>fs.existsSync(path.join(R0,f)));
if(!have.length){ console.log("\nSKIPPED (data absent): test_real_autobuild — test/data/real not found (private data, not in the public backup); 0 checks run"); process.exit(0); }
const IO=require("../src/xps_io.js"), C=require("../src/xps_core.js"), G=require("../src/xps_guide.js"), A=require("../src/xps_autobuild.js"), L=require("../src/xps_lib.js");
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info!=null&&info!==""?"  ("+info+")":"")); }
const info=(name,v)=>console.log("  info "+name+(v!=null?"  ("+v+")":""));
const rd=async f=>(await IO.readFile(path.basename(f),new Uint8Array(fs.readFileSync(path.join(R0,f))))).spectra;
const t0=Date.now();
/* lines of a model as raw BE (a doublet = its two lines) */
const linesOf=comps=>[].concat(...comps.map(c=>c.doublet?[c.be,c.be+c.doublet.dE]:[c.be])).sort((a,b)=>a-b);
function agree(auto,cur){ const near=cur.map(b=>Math.min(...auto.map(a=>Math.abs(a-b)))); return {sameCount:auto.length===cur.length,within:near.filter(d=>d<=0.3).length,n:cur.length,maxD:Math.max(...near)}; }
function runAB(name,P,ctx,first){ const before=JSON.stringify(P.comps); const t=Date.now(); const R=A.autoBuild(P,ctx,{first}); const F=R.final;
  const fw=R.S.P.comps.map(c=>F.lines.find(l=>l.id===c.id).fwhm); let capOK=true; R.S.P.comps.forEach((a,i)=>R.S.P.comps.forEach((b,j)=>{ if(j<=i) return; const qa=(R.S.P.cons[a.id]||{}).be, qb=(R.S.P.cons[b.id]||{}).be; if(qa&&qa.mode==="link"||qb&&qb.mode==="link") return; if(a.recipe&&a.recipe===b.recipe) return; if(Math.abs(a.be-b.be)<0.5*(fw[i]+fw[j])/2-1e-9) capOK=false; }));
  ok(JSON.stringify(P.comps)===before&&capOK&&(F.established||F.reasons.length>0),`${name}: runs, input unchanged, separation cap holds, result has a status`,`${((Date.now()-t)/1000).toFixed(1)} s`);
  info(`${name}`,`${R.steps.length} steps [${R.steps.map(s=>s.kind).join(",")}] → ${R.S.P.comps.length} comps, χ²ν ${F.chi2nu.toFixed(3)}, ρ₁ ${F.rho1.toFixed(2)}, same min ${F.multi?F.multi.near+"/"+F.multi.N:"—"} · ${F.established?"ESTABLISHED":"not established: "+F.reasons.join("; ")}`);
  info(`${name} lines`,F.lines.map(l=>`${l.be.toFixed(2)}${l.separate?"":"*"}${l.added?"+":""}${l.weak?"w":""} (${l.separate?l.areaPct.toFixed(1)+"%":"grp"})`).join(" · ")+(F.groups.length?` · groups: ${F.groups.map(g=>g.names.length+" lines "+(100*g.sig/g.area).toFixed(1)+"%").join(", ")}`:""));
  return R; }
(async()=>{
/* ---- K1/K2: CasaXPS-curated models ---- */
if(fs.existsSync(path.join(R0,"kara/casa_240213/Integrated.vms"))){
  console.log("== K CasaXPS-curated Ba 4d · C 1s · O 1s (2 samples)");
  const V=await rd("kara/casa_240213/Integrated.vms"); const withM=V.filter(s=>s.meta.casaModel&&s.meta.casaModel.regions.length);
  const smp=b=>/Bare/.test(b)?"K1":"K2"; const shifts={};
  for(const key of ["C1s","O1s","Ba4d"]) for(const sp of withM.filter(s=>s.meta.blockId.includes(key))){
    const cm=sp.meta.casaModel, M=C.casaModelToINNO(cm,sp.hv,{step:sp.step}); const spec={be:sp.be,y:sp.y,c:sp.c}; const nm=`${smp(sp.meta.blockId)} ${sp.reg.el} ${sp.reg.orb}`;
    const cur=linesOf(M.comps); info(`${nm} curated (CasaXPS)`,`${cur.length} lines: ${cur.map(b=>b.toFixed(2)).join(" · ")} (raw BE, window ${M.lo.toFixed(1)}–${M.hi.toFixed(1)})`);
    // (a) from Auto suggest on the curated window
    const R_=C.sliceRegion(sp.be,sp.y,M.lo,M.hi); const sg=C.suggestModel(Float64Array.from(R_.E),Float64Array.from(R_.y),{bg:{type:"smart",n:5},el:sp.reg.el,orb:sp.reg.orb,lib:L});
    let P=C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:{type:"smart",n:5},comps:sg.comps,cons:sg.cons}); const rf=C.refineByResidual(P,{}); P=C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:{type:"smart",n:5},comps:rf.prob.comps,cons:rf.prob.cons});
    const Ra=runAB(nm+" from Auto suggest",P,{el:sp.reg.el,orb:sp.reg.orb,present:PRES},true); const ga=agree(linesOf(Ra.S.P.comps),cur);
    info(`${nm} vs curated (Auto suggest start)`,`count ${linesOf(Ra.S.P.comps).length} vs ${cur.length} · curated lines within 0.3 eV: ${ga.within}/${ga.n} (max ${ga.maxD.toFixed(2)} eV)`);
    // (b) from the recipe (C 1s → its Δ → O 1s, Ba 4d)
    const k=smp(sp.meta.blockId); let Pr=null, ctx={el:sp.reg.el,orb:sp.reg.orb,present:PRES};
    if(key==="C1s"){ const Mc=G.recipeModel("c1s",{shift:0,be:sp.be,y:sp.y,present:ctx.present}); const P0=C.makeProblem(spec,{lo:Mc.lo,hi:Mc.hi,bg:Mc.bg,comps:Mc.comps,cons:Mc.cons}); const f=G.c1sFit(P0,{}); shifts[k]=f.info.delta; Pr=C.makeProblem(spec,{lo:Mc.lo,hi:Mc.hi,bg:f.P.bg,comps:f.P.comps,cons:f.P.cons}); }
    else if(key==="O1s"){ const sh=shifts[k]||0; const b={shift:sh,present:ctx.present,flags:{lat2:false}}; const Mo=G.recipeModel("o1s",Object.assign({be:sp.be,y:sp.y},b)); const P0=C.makeProblem(spec,{lo:Mo.lo,hi:Mo.hi,bg:Mo.bg,comps:Mo.comps,cons:Mo.cons}); const f=G.o1sFit(P0,{build:b}); Pr=C.makeProblem(spec,{lo:Mo.lo,hi:Mo.hi,bg:f.P.bg,comps:f.P.comps,cons:f.P.cons}); ctx.shift=sh; }
    else { const sh=shifts[k]||0; const Mb=G.recipeModel("ba4d",{shift:sh,be:sp.be,y:sp.y,species:G.parDefaults(G.recipeById("ba4d")),flags:G.parFlags(G.recipeById("ba4d"),{})}); const P0=C.makeProblem(spec,{lo:Mb.lo,hi:Mb.hi,bg:Mb.bg,comps:Mb.comps,cons:Mb.cons}); const f=G.stagedFit(P0,{}); Pr=C.makeProblem(spec,{lo:Mb.lo,hi:Mb.hi,bg:f.P.bg,comps:f.P.comps,cons:f.P.cons}); ctx.shift=sh; }
    const Rb=runAB(nm+" from the recipe",Pr,ctx,false); const gb=agree(linesOf(Rb.S.P.comps),cur);
    info(`${nm} vs curated (recipe start)`,`count ${linesOf(Rb.S.P.comps).length} vs ${cur.length} · curated lines within 0.3 eV: ${gb.within}/${gb.n} (max ${gb.maxD.toFixed(2)} eV)`); }
} else info("K data absent — skipped");
/* ---- A1/A2: VGD mixed-metal oxide, recipes C 1s → O 1s, and from Auto suggest ---- */
for(const [k,dir] of A_DIRS.map((d,i)=>["A"+(i+1),d])){ if(!fs.existsSync(path.join(R0,dir,"O1s Scan.VGD"))){ info(k+" absent — skipped"); continue; }
  console.log(`== ${k} VGD mixed-metal oxide: C 1s recipe → Δ → O 1s recipe → Auto build; and from Auto suggest`);
  const c1=(await rd(dir+"/C1s Scan.VGD"))[0], o1=(await rd(dir+"/O1s Scan.VGD"))[0]; const pres=PRES;
  const Mc=G.recipeModel("c1s",{shift:0,be:c1.be,y:c1.y,present:pres}); const sc={be:c1.be,y:c1.y,c:c1.c}; const fc=G.c1sFit(C.makeProblem(sc,{lo:Mc.lo,hi:Mc.hi,bg:Mc.bg,comps:Mc.comps,cons:Mc.cons}),{});
  const sh=fc.info.delta; info(`${k} C 1s recipe`,`C–C/C–H ${fc.info.refBE.toFixed(2)} eV raw → Δ ${sh.toFixed(2)} eV, ${fc.P.comps.length} lines`);
  runAB(`${k} C 1s from the recipe`,C.makeProblem(sc,{lo:Mc.lo,hi:Mc.hi,bg:fc.P.bg,comps:fc.P.comps,cons:fc.P.cons}),{el:"C",orb:"1s",present:pres,shift:0},false);
  const so={be:o1.be,y:o1.y,c:o1.c}; const b={shift:sh,present:pres,flags:{lat2:false}}; const Mo=G.recipeModel("o1s",Object.assign({be:o1.be,y:o1.y},b));
  const fo=G.o1sFit(C.makeProblem(so,{lo:Mo.lo,hi:Mo.hi,bg:Mo.bg,comps:Mo.comps,cons:Mo.cons}),{build:b}); const sto=C.fitStats(fo.P,fo.r);
  info(`${k} O 1s recipe fit`,`${fo.P.comps.length} lines, χ²ν ${sto.chi2nu.toFixed(3)}, ρ₁ ${sto.rho1.toFixed(2)}; low lattice: ${fo.info.latlo?fo.info.latlo.verdict+" (ΔAICc "+(fo.info.latlo.dAICc||0).toFixed(1)+")":"—"}`);
  const Ro=runAB(`${k} O 1s from the recipe`,C.makeProblem(so,{lo:Mo.lo,hi:Mo.hi,bg:fo.P.bg,comps:fo.P.comps,cons:fo.P.cons}),{el:"O",orb:"1s",present:pres,shift:sh},false);
  const lo=Ro.S.P.comps.filter(c=>c.ab).map(c=>(c.be+sh).toFixed(2)); info(`${k} O 1s lines added by Auto build (corrected BE)`,lo.join(", ")||"none");
  const R_=C.sliceRegion(o1.be,o1.y,Mo.lo,Mo.hi); const sg=C.suggestModel(Float64Array.from(R_.E),Float64Array.from(R_.y),{bg:{type:"smart",n:5},el:"O",orb:"1s",lib:L});
  let P=C.makeProblem(so,{lo:Mo.lo,hi:Mo.hi,bg:{type:"smart",n:5},comps:sg.comps,cons:sg.cons}); const rf=C.refineByResidual(P,{}); P=C.makeProblem(so,{lo:Mo.lo,hi:Mo.hi,bg:{type:"smart",n:5},comps:rf.prob.comps,cons:rf.prob.cons});
  runAB(`${k} O 1s from Auto suggest`,P,{el:"O",orb:"1s",present:pres,shift:sh},true);
  const ba=(await rd(dir+"/Ba4d Scan.VGD"))[0]; const sb={be:ba.be,y:ba.y,c:ba.c}; const Rb=G.recipeById("ba4d"); const Mb=G.recipeModel("ba4d",{shift:sh,be:ba.be,y:ba.y,species:G.parDefaults(Rb),flags:G.parFlags(Rb,{})});
  const fb=G.stagedFit(C.makeProblem(sb,{lo:Mb.lo,hi:Mb.hi,bg:Mb.bg,comps:Mb.comps,cons:Mb.cons}),{}); runAB(`${k} Ba 4d from the recipe`,C.makeProblem(sb,{lo:Mb.lo,hi:Mb.hi,bg:fb.P.bg,comps:fb.P.comps,cons:fb.P.cons}),{el:"Ba",orb:"4d",present:pres,shift:sh},false); }
/* ---- S3: Ce 3d recipe (Mixed) as in test_real_ce ---- */
if(fs.existsSync(path.join(R0,"ce/Ce3d Scan.VGD"))){ console.log("== S3 CeO₂: Ce 3d recipe (Mixed) + staged fit → Auto build");
  const c1=(await rd("ce/C1s Scan.VGD"))[0], ce=(await rd("ce/Ce3d Scan.VGD"))[0]; const pk=C.peakMaxPos(c1.be,c1.y,280,292,0.25); const sh=+(284.8-pk.be).toFixed(3);
  const M=G.recipeModel("ce3d",{shift:sh,be:ce.be,y:ce.y,opt:"mix"}); const spec={be:ce.be,y:ce.y,c:ce.c}; const s=G.stagedFit(C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons}),{});
  const R=runAB("S3 Ce 3d from the recipe",C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:s.P.bg,comps:s.P.comps,cons:s.P.cons}),{el:"Ce",orb:"3d",present:["Ce","O","C","Si"],shift:sh},false);
  info("S3 Ce 3d lines added by Auto build (corrected BE)",R.S.P.comps.filter(c=>c.ab).map(c=>(c.be+sh).toFixed(2)).join(", ")||"none"); }
console.log(`\ntest_real_autobuild: ${pass} passed, ${fail} failed, 0 skipped (${((Date.now()-t0)/1000).toFixed(0)} s)`);
if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; }
})().catch(e=>{ console.error(e); process.exitCode=1; });
