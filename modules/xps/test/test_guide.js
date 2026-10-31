/* Fit guidance (SPEC §17, v1.9.3): fit check · region recipes · staged fit
   node test/test_guide.js        (from modules/xps; about 2–3 min; QUICK=1 → fewer synthetic Ce spectra)
   1 recipe data (Ce 3d) · model builder anchors to the charge correction · alignment recovers a known offset
   2 noise statistics from second differences (white · MA(1))
   3 fit check — positive controls: each rule fires on a synthetic spectrum built to trigger it
   4 fit check — false positives on well-posed fits (true model family): all_elements cases, in-fit Shirley + slope (= truth)
     + the same peaks on a flat baseline with static Shirley/Smart (reported: warnings only where χ²ν shows a real misfit)
   5 one-click actions (applyAction) and that the finding disappears after the action + refit
   6 staged fit: never worse than the normal fit (χ²_d), normal path inside = fitOnce bit for bit
   7 synthetic Ce 3d (truth from the recipe, Ce³⁺ fraction 0 · 0.2 · 0.5, Poisson noise, Shirley-type background + slope):
     recipe + staged fit recovers the Ce³⁺ and u‴ fractions; bound hits vs a naive suggestModel fit
   8 v1.9.6 transition-metal 2p3/2 envelope recipes (Biesinger 2011/2017): numbers = journal transcripts (test/ref/biesinger), shapes,
     PE columns, constraint structure, Fe₃O₄ 1:2, synthetic 2–3-species mixtures (Ni Co Fe Mn Cr Cu) recovered, σ, Cu shake-up, fit check */
"use strict";
const C=require("../src/xps_core.js"), G=require("../src/xps_guide.js"), L=require("../src/xps_lib.js"), U=require("./synth_util.js"), AE=require("./all_elements.js");
const QUICK=!!process.env.QUICK;
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info!=null&&info!==""?"  ("+info+")":"")); }
const ids=d=>d.map(x=>x.id), warnIds=d=>d.filter(x=>x.sev!=="info").map(x=>x.id);
const t0=Date.now();

/* ---------- helpers ---------- */
function spectrum(lines,o){ o=o||{}; const rg=U.rng(o.seed||5); const lo=o.lo!=null?o.lo:278, hi=o.hi!=null?o.hi:296, step=o.step||0.05;
  const n=Math.round((hi-lo)/step)+1, E=Float64Array.from({length:n},(_,i)=>lo+i*step), ctx={emin:E[0],emax:E[n-1]};
  const P=new Float64Array(n); lines.forEach(([be,A,F])=>{ const v=C.compEval(C.newComp("GLS",be,{F:F||1.2,m:0.3,area:A}),E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P), cf=o.cf||1; const y=Float64Array.from(P,(p,i)=>rg.pois((2000+(o.slope||0)*(E[i]-lo)+(o.k!=null?o.k:0.05)*Q[i]+p)*cf)/cf);
  return {be:E,y,c:cf}; }
function fitDiag(spec,model,o){ o=o||{}; const P=C.makeProblem(spec,model); const r=C.fitOnce(P,o.fit||{}); const fin=C.finalize(P,r);
  const d=G.diagnose({E:P.E,y:P.y,fin,el:o.el,orb:o.orb,bgType:P.bg.type,spec:{be0:spec.be[0],be1:spec.be[spec.be.length-1]},recipe:o.recipe,dEref:o.dEref});
  return {P,r,fin,d}; }
const comp=(be,A,F,id)=>C.newComp("GLS",be,{F:F||1.2,m:0.3,area:A||10000,id});
const ACT={type:"active",n:5,cons:{s:{mode:"free"}}};

/* ---------- 1 recipe ---------- */
console.log("== 1 Ce 3d recipe: data, anchoring to the charge correction, alignment");
{ const R=G.recipeById("ce3d");
  ok(R&&R.lines.length===5&&R.lines.filter(l=>l.state==="Ce3+").length===2&&R.options.map(o=>o.keys.length*2).join()==="6,4,10","5 spin–orbit pairs (Ce⁴⁺ 3 · Ce³⁺ 2), options 6 · 4 · 10 lines",R.options.map(o=>o.label).join(" · "));
  // literature values reachable from the start ranges (SPEC §17.3 table): v-line positions and ΔE of every accessible source lie inside start ± beTol / the ΔE range
  const lit={v0:[880.2,880.6,880.7],v:[882.2,882.4,882.5],v1:[884.5,884.7,885.6],v2:[888.7,888.8,889.2],v3:[898.1,898.2,898.3]}, dEs=[18.3,18.4,18.5,18.6,18.7];
  const inR=R.lines.every(l=>lit[l.key].every(b=>Math.abs(b-l.be)<=R.beTol+1e-9)), dOK=dEs.every(d=>d>=Math.min(R.dE.lo,(R.dEs&&R.dEs["Ce3+"]||R.dE).lo)&&d<=R.dE.hi);
  ok(inR&&dOK&&Math.abs(R.r-2/3)<1e-12,"literature positions within start ± 0.6 eV, ΔE 18.3–18.7 eV inside the free range, u/v area 2/3",`ΔE ${R.dE.lo}–${R.dE.hi}`);
  const M=G.recipeModel("ce3d",{shift:1.3,opt:"mix",align:false}); const c0=M.comps.find(c=>c.rk==="v3");
  ok(Math.abs(c0.be-(898.2-1.3))<1e-9&&Math.abs(M.lo-(875-1.3))<1e-9&&Math.abs(M.hi-(925-1.3))<1e-9,"anchors to the charge correction: raw BE = literature − Δ (positions and window)",`v‴ ${c0.be.toFixed(2)} · window ${M.lo.toFixed(1)}–${M.hi.toFixed(1)}`);
  const anc=M.comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); const ancS={}; M.comps.forEach(c=>{ const x=ancS[c.rstate]; if(!x||c.area>x.area) ancS[c.rstate]=c; });
  const links=M.comps.every(c=>{ const q=M.cons[c.id].dE, aS=ancS[c.rstate], d=G.dEof?null:null; return c===aS?(q.mode==="free"&&Math.abs(c.doublet.dE-(c.rstate==="Ce3+"?18.2:18.6))<1e-9):(q.mode==="link"&&q.ref===aS.id); })&&M.comps.filter(c=>c!==anc).every(c=>M.cons[c.id].m.mode==="link");
  ok(links&&M.comps.every(c=>c.doublet&&c.doublet.r===R.r&&M.cons[c.id].r.mode==="fixed"&&M.cons[c.id].q.mode==="fixed"&&M.cons[c.id].be.hi-M.cons[c.id].be.lo>1.19&&M.cons[c.id].be.hi-M.cons[c.id].be.lo<1.21),"pairs: ΔE per oxidation state (Ce⁴⁺ 18.6, Ce³⁺ 18.2 eV — Bêche 2008; linked within a state), area ratio and width factor fixed, BE free ± 0.6 eV");
  ok(M.bg.type==="active"&&M.bg.cons.s.mode==="free","background = in-fit Shirley + slope (SPEC §17.3)");
  // alignment: synthetic recipe spectrum shifted by a known amount relative to the charge-corrected positions
  const rg=U.rng(3); const n=601, E=Float64Array.from({length:n},(_,i)=>870+i*0.1), ctx={emin:E[0],emax:E[n-1]};
  const devs=[]; for(const off of [-2.2,-0.4,0.3,1.4,2.6]){ const P=new Float64Array(n); R.lines.forEach(l=>{ if(l.state!=="Ce4+") return; const c=C.newComp("GLS",l.be+off,{F:l.F,m:0.3,area:40000}); c.doublet={dE:18.5,r:2/3,q:1}; const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
    const y=Float64Array.from(P,p=>rg.pois((3000+p)*2)/2); const A=G.recipeModel("ce3d",{shift:0,be:E,y,opt:"ce4"}); devs.push(A.recipe.align-off); }
  ok(devs.every(d=>Math.abs(d)<=0.11),"alignment (NNLS match within ±3 eV) recovers a known offset of the lines",devs.map(d=>d.toFixed(2)).join(" "));
  const Mc=G.recipeModel("ce3d",{shift:0,be:Float64Array.from({length:301},(_,i)=>880+i*0.1),y:new Float64Array(301).fill(1000),opt:"ce4",align:false});
  ok(Mc.recipe.clipped.includes("low")&&Math.abs(Mc.lo-880)<1e-9,"window clipped to the spectrum and reported",Mc.recipe.clipped.join());
  ok(G.complexInfo("Ce","3d")&&G.complexInfo("Ce","3d").recipes[0]==="ce3d"&&G.complexInfo("Ni","2p").kind==="satellite"&&G.complexInfo("Co","3p").kind==="multiplet"&&!G.complexInfo("C","1s")&&!G.complexInfo("Ag","3d"),"complex-region rule: Ce 3d (recipe) · Ni 2p (satellites) · Co 3p (multiplet) · not C 1s / Ag 3d");
}

/* ---------- 2 noise ---------- */
console.log("== 2 noise statistics from second differences (MA(1))");
{ const rg=U.rng(11); const n=4000, wn=Float64Array.from({length:n},()=>rg.g()); const a=G.noiseMA1(wn,null);
  ok(a.rho<0.04&&Math.abs(a.s2-1)<0.06,"white noise: ρ ≈ 0, σ² ≈ 1",`ρ ${a.rho.toFixed(3)} σ² ${a.s2.toFixed(3)}`);
  const th=0.6, e=Float64Array.from({length:n+1},()=>rg.g()); const ma=Float64Array.from({length:n},(_,i)=>(e[i+1]+th*e[i])/Math.sqrt(1+th*th)); const b=G.noiseMA1(ma,null), rt=th/(1+th*th);
  ok(Math.abs(b.rho-rt)<0.04&&Math.abs(b.lrv-(1+2*rt))<0.12,"MA(1) noise (lag-1 ρ 0.44, like VGD): ρ and long-run variance recovered",`ρ ${b.rho.toFixed(3)}/${rt.toFixed(3)} · LRV ${b.lrv.toFixed(3)}/${(1+2*rt).toFixed(3)}`);
  const sm=Float64Array.from({length:n},(_,i)=>1e4*Math.exp(-Math.pow((i-2000)/300,2))+rg.g()); const s=G.noiseMA1(sm,null);
  ok(s.rho<0.05&&Math.abs(s.s2-1)<0.08,"smooth signal (peak 300 channels wide) does not leak into the noise estimate",`ρ ${s.rho.toFixed(3)} σ² ${s.s2.toFixed(3)}`);
}

/* ---------- 3 positive controls ---------- */
console.log("== 3 fit check — every rule fires on a spectrum built to trigger it");
const POS={};
{ const sp1=spectrum([[285,30000]]);
  let x=fitDiag(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(285,20000,1.0,"a")],cons:{a:{F:{mode:"free",lo:0.3,hi:0.9}}}});
  const b=x.d.find(f=>f.id==="bound"&&/FWHM/.test(f.label)); POS.bound=!!b; ok(!!b&&b.sev==="warn"&&/upper bound/.test(b.label),"R1 bound: FWHM range below the true width → 'FWHM at upper bound'",b&&b.label);
  x=fitDiag(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(284.75,15000,1.2,"a"),comp(285.25,15000,1.2,"b")],cons:{a:{F:{mode:"free",lo:0.3,hi:4},be:{mode:"free",lo:284,hi:285.5}},b:{F:{mode:"free",lo:0.3,hi:4},be:{mode:"free",lo:284.5,hi:286}}}});
  POS.corr=ids(x.d).includes("corr"); POS.merge=ids(x.d).includes("merge");
  ok(POS.corr,"R2 correlation: one peak fitted by two free lines → |r| > 0.95 pair",(x.d.find(f=>f.id==="corr")||{}).label);
  ok(POS.merge||ids(x.d).includes("ill")||ids(x.d).includes("dead"),"R4 merge (or ill-determined/dead): the two lines describe one peak",warnIds(x.d).join(","));
  x=fitDiag(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(285,30000,1.2,"a"),comp(290,300,1.2,"b")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}},b:{F:{mode:"fixed"},be:{mode:"free",lo:289,hi:291}}}});
  POS.ill=x.d.some(f=>["ill","dead","neg"].includes(f.id)&&f.ids.includes("b")); ok(POS.ill,"R3 negligible/ill-determined: a line where there is no peak",x.d.filter(f=>f.ids.includes("b")).map(f=>f.label).join(" | "));
  x=fitDiag(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(284.5,10000,2.0,"a"),comp(285.6,10000,0.8,"b")],cons:{}},{fit:{maxit:2}});
  POS.stop=ids(x.d).includes("stop"); ok(POS.stop,"R5 stop: iteration limit → 'Fit did not converge'",x.fin.stop);
  const sp2=spectrum([[285,30000],[289,30000]]);
  x=fitDiag(sp2,{lo:279,hi:288.6,bg:{type:"shirley",n:5},comps:[comp(285,30000,1.2,"a")],cons:{}});
  const ed=x.d.find(f=>f.id==="edge"); POS.edge=!!ed; ok(!!ed&&/high-BE/.test(ed.label)&&ed.act&&ed.act.type==="window","R6 edge: window ends on the flank of the next peak → 'Peak cut at the high-BE end' + widen action",ed&&ed.label+" → "+(ed.act&&ed.act.label));
  x=fitDiag(sp2,{lo:285.4,hi:292,bg:{type:"linear",n:5},comps:[comp(289,30000,1.2,"a")],cons:{}});
  POS.bgabove=ids(x.d).includes("bgabove"); ok(POS.bgabove,"R7 background above data: linear background from an endpoint on a peak flank",(x.d.find(f=>f.id==="bgabove")||{}).label);
  const sp3=spectrum([[285,30000],[287.2,6000]]);
  x=fitDiag(sp3,{lo:278,hi:296,bg:ACT,comps:[comp(285,30000,1.2,"a")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}}}});
  const rs=x.d.find(f=>f.id==="resid"); POS.resid=!!rs; ok(!!rs&&rs.act&&rs.act.type==="add"&&Math.abs(rs.act.be-287.2)<0.5,"R8 residual: a 20 % line left out → 'Missing intensity near 287 eV' + add action",rs&&rs.label);
  x=fitDiag(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(285,30000,1.2,"a")],cons:{}},{el:"Ce",orb:"3d"});
  const cx=x.d.find(f=>f.id==="complex"); POS.complex=!!cx; ok(!!cx&&cx.sev==="info"&&cx.act&&cx.act.type==="recipe","R9 complex region (Ce 3d) → info + 'Apply recipe'",cx&&cx.label);
}

/* ---------- 4 false positives ---------- */
console.log("== 4 fit check — false positives on well-posed fits");
const FP={n:0,warned:0,byRule:{},nS:0,warnedS:0,chiS:[]};
{ const R=C.rng(77); const CS=AE.cases().filter(cs=>!C.multipletInfo(cs.el,cs.orb));
  for(const cs of CS) for(const two of [false,true]){ const S=AE.synth(R,cs,{two});
    const comps=S.truth.map((t,i)=>{ const c=C.clone(t); c.id="c"+i; c.be+=0.2*R.g(); c.F*=1+0.15*R.g(); c.area*=0.7; return c; });
    const cons={}; comps.forEach((c,i)=>{ const q=cons[c.id]={F:{mode:"free",lo:0.3,hi:4}}; if(i) q.m={mode:"link",ref:"c0",k:1,b:0};
      if(c.doublet){ q.r={mode:"fixed"}; q.q={mode:"fixed"}; c.doublet.dE=S.info.dE; q.dE=i?{mode:"link",ref:"c0",k:1,b:0}:C.dEConstraint(S.info); } });
    // (a) true background family: in-fit Shirley + slope
    { const P={E:S.E,y:S.y,c:S.cf,bg:{type:"active",n:5,cons:{s:{mode:"free"}}},comps:C.clone(comps),cons:C.clone(cons)};
      const r=C.fitOnce(P,{}); const fin=C.finalize(P,r); const d=G.diagnose({E:P.E,y:P.y,fin,el:cs.el,orb:cs.orb,bgType:"active",dEref:S.info&&S.info.dE,spec:{be0:S.E[0],be1:S.E[S.E.length-1]}});
      FP.n++; const w=d.filter(x=>x.sev!=="info"); if(w.length){ FP.warned++; console.log(`     ${cs.el} ${cs.orb} ${two?"two":"one"}: ${w.map(x=>x.label).join(" | ")}`); } w.forEach(x=>FP.byRule[x.id]=(FP.byRule[x.id]||0)+1); }
    // (b) flat baseline + Shirley step (no slope) fitted with static Shirley / Smart (every third case, both)
    if(FP.n%3===0){ const n=S.E.length, ctx={emin:S.E[0],emax:S.E[n-1]}; const P0=new Float64Array(n); S.truth.forEach(c=>{ const v=C.compEval(c,S.E,ctx); for(let i=0;i<n;i++) P0[i]+=v[i]; }); const Q=C.cumtrapz(S.E,P0);
      const y=Float64Array.from(P0,(p,i)=>{ const m=(3000+0.06*Q[i]+p)*S.cf; return Math.max(0,Math.round(m+Math.sqrt(m)*R.g()))/S.cf; });
      for(const bt of ["smart","shirley"]){ const P={E:S.E,y,c:S.cf,bg:{type:bt,n:5},comps:C.clone(comps),cons:C.clone(cons)}; const r=C.fitOnce(P,{}); const fin=C.finalize(P,r);
        const d=G.diagnose({E:P.E,y:P.y,fin,el:cs.el,orb:cs.orb,bgType:bt,dEref:S.info&&S.info.dE,spec:{be0:S.E[0],be1:S.E[n-1]}}); FP.nS++; if(d.some(x=>x.sev!=="info")){ FP.warnedS++; FP.chiS.push(fin.stats.chi2nu); console.log(`     [${bt}] ${cs.el} ${cs.orb} ${two?"two":"one"} χ²ν ${fin.stats.chi2nu.toFixed(2)}: ${d.filter(x=>x.sev!=="info").map(x=>x.label).join(" | ")}`); } } } }
  ok(FP.warned===0,`true model family (in-fit Shirley + slope): ${FP.n} fits, warnings on ${FP.warned}`,JSON.stringify(FP.byRule));
  const minChi=FP.chiS.length?Math.min(...FP.chiS):null;
  ok(FP.warnedS<=0.1*FP.nS&&(minChi==null||minChi>2),`static Shirley/Smart on a flat baseline: ${FP.nS} fits, warnings on ${FP.warnedS} — all with χ²ν > 2 (static Shirley failing on wide windows: a real misfit)`,minChi==null?"none":`min χ²ν of warned ${minChi.toFixed(2)}`);
}

/* ---------- 5 actions ---------- */
console.log("== 5 one-click actions");
{ const sp1=spectrum([[285,30000]]);
  // merge → drop the smaller line → refit → finding gone
  let m={lo:278,hi:296,bg:ACT,comps:[comp(284.9,20000,1.2,"a"),comp(285.1,8000,1.2,"b")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}},b:{F:{mode:"link",ref:"a",k:1,b:0},m:{mode:"link",ref:"a",k:1,b:0}}}};
  let x=fitDiag(sp1,m); const f1=x.d.find(f=>["merge","ill","dead"].includes(f.id)&&f.act&&f.act.type==="drop");
  if(f1){ const mm={lo:m.lo,hi:m.hi,bg:m.bg,comps:x.fin.comps,cons:x.P.cons}; G.applyAction(mm,f1.act); const y2=fitDiag(sp1,mm); ok(mm.comps.length===1&&!y2.d.some(f=>f.sev!=="info"),"drop: removes the line (links re-anchored) → refit has no warnings",f1.label); }
  else ok(false,"drop: merge finding expected",ids(x.d).join());
  // bound F → link to the anchor
  const sp2=spectrum([[285,30000],[287,15000]]);
  m={lo:278,hi:296,bg:ACT,comps:[comp(285,30000,1.2,"a"),comp(287,15000,1.0,"b")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}},b:{F:{mode:"free",lo:0.3,hi:0.8}}}};
  x=fitDiag(sp2,m); const fb=x.d.find(f=>f.id==="bound"&&f.act&&f.act.type==="link");
  if(fb){ const mm={lo:m.lo,hi:m.hi,bg:m.bg,comps:x.fin.comps,cons:x.P.cons}; G.applyAction(mm,fb.act); const y2=fitDiag(sp2,mm); ok(mm.cons.b.F.mode==="link"&&!y2.d.some(f=>f.id==="bound"&&/FWHM/.test(f.label)),"link: width at bound linked to the anchor → bound finding gone after refit",fb.label); }
  else ok(false,"link: bound finding with link action expected",ids(x.d).join());
  // edge → widen window
  const sp3=spectrum([[285,30000],[289,30000]]); m={lo:279,hi:288.6,bg:{type:"shirley",n:5},comps:[comp(285,30000,1.2,"a")],cons:{}};
  x=fitDiag(sp3,m); const fe=x.d.find(f=>f.id==="edge"); const mm={lo:m.lo,hi:m.hi,bg:m.bg,comps:x.fin.comps,cons:{}};
  G.applyAction(mm,fe.act,{be0:sp3.be[0],be1:sp3.be[sp3.be.length-1]}); ok(mm.hi>m.hi+1.9&&mm.hi<=sp3.be[sp3.be.length-1],"window: widens the cut side within the spectrum",`${m.hi} → ${mm.hi.toFixed(2)}`);
  // residual → add line → refit → residual finding gone
  const sp4=spectrum([[285,30000],[287.2,6000]]); m={lo:278,hi:296,bg:ACT,comps:[comp(285,30000,1.2,"a")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}}}};
  x=fitDiag(sp4,m); const fr=x.d.find(f=>f.id==="resid"); const ma={lo:m.lo,hi:m.hi,bg:m.bg,comps:x.fin.comps,cons:x.P.cons}; G.applyAction(ma,fr.act); const y4=fitDiag(sp4,ma);
  ok(ma.comps.length===2&&!y4.d.some(f=>f.id==="resid")&&Math.abs(y4.fin.comps[1].be-287.2)<0.1,"add: line at the residual blob → refit explains it",`new line ${y4.fin.comps[1].be.toFixed(2)} eV`);
  // fix / widen / offset / bg
  const mf={lo:0,hi:1,bg:{type:"shirley",n:5},comps:[comp(285,1,1.2,"a"),comp(287,1,1.2,"b")],cons:{a:{be:{mode:"free",lo:284.5,hi:285.5}}}};
  G.applyAction(mf,{type:"fix",id:"a",k:"F",v:1.5}); G.applyAction(mf,{type:"widen",id:"a",k:"be",by:0.3}); G.applyAction(mf,{type:"offset",id:"b",ref:"a"}); G.applyAction(mf,{type:"bg",bg:"smart"});
  ok(mf.comps[0].F===1.5&&mf.cons.a.F.mode==="fixed"&&Math.abs(mf.cons.a.be.lo-284.2)<1e-9&&Math.abs(mf.cons.a.be.hi-285.8)<1e-9&&mf.cons.b.be.mode==="link"&&Math.abs(mf.cons.b.be.b-2)<1e-9&&mf.bg.type==="smart","fix · widen · fix separation · background actions edit the model as stated");
}

/* ---------- 6 staged fit ---------- */
console.log("== 6 staged fit — never worse than the normal fit");
{ const probs=[]; const sp1=spectrum([[285,30000],[286.3,12000],[288.6,4000]]);
  probs.push(C.makeProblem(sp1,{lo:278,hi:296,bg:ACT,comps:[comp(285.3,20000,1.6,"a"),comp(286.0,20000,1.0,"b"),comp(288.2,2000,1.4,"c")],cons:{a:{F:{mode:"free",lo:0.3,hi:4}},b:{F:{mode:"free",lo:0.3,hi:4}},c:{F:{mode:"free",lo:0.3,hi:4}}}}));
  probs.push(C.makeProblem(sp1,{lo:278,hi:296,bg:{type:"smart",n:5},comps:[comp(284.6,5000,2.5,"a"),comp(286.9,5000,0.6,"b"),comp(289.2,1000,0.8,"c")],cons:{}}));
  const R=G.recipeById("ce3d"); const rg=U.rng(21); const n=601, E=Float64Array.from({length:n},(_,i)=>870+i*0.1), ctx={emin:E[0],emax:E[n-1]}; const P=new Float64Array(n);
  R.lines.forEach((l,i)=>{ const c=C.newComp("GLS",l.be+0.2,{F:l.F,m:0.3,area:[8000,40000,20000,50000,40000][i]}); c.doublet={dE:18.45,r:2/3,q:1}; const v=C.compEval(c,E,ctx); for(let j=0;j<n;j++) P[j]+=v[j]; });
  const Q=C.cumtrapz(E,P); const y=Float64Array.from(P,(p,i)=>rg.pois((2500+0.1*Q[i]+p)*2)/2); const sce={be:E,y,c:2};
  const M=G.recipeModel("ce3d",{shift:0,be:E,y,opt:"mix"}); probs.push(C.makeProblem(sce,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons}));
  const M2=G.recipeModel("ce3d",{shift:0,be:E,y,opt:"mix",align:false}); M2.comps.forEach(c=>c.be-=0.5); probs.push(C.makeProblem(sce,{lo:M2.lo,hi:M2.hi,bg:{type:"smart",n:5},comps:M2.comps,cons:M2.cons}));
  let worse=0, same=0, ident=0; const rows=[];
  probs.forEach((p,i)=>{ const s=G.stagedFit(p,{}); const st=C.fitStats(s.P,s.r); const ref=C.clone(C.stripProb(p)); ref.E=p.E; ref.y=p.y; const rr=C.fitOnce(ref,{}); const sr=C.fitStats(ref,rr);
    if(st.chi2d>sr.chi2d*(1+1e-12)) worse++; if(s.staged.used==="normal"&&rr.chi2===s.r.chi2) ident++; if(s.staged.chi2dNormal===sr.chi2d) same++;
    rows.push(`#${i} ${s.staged.used} χ²_d staged ${s.staged.chi2dStaged.toFixed(1)} · normal ${s.staged.chi2dNormal.toFixed(1)}`); });
  ok(worse===0,"kept result χ²_d ≤ normal fit χ²_d on every problem",rows.join(" | "));
  ok(same===probs.length,"normal path inside stagedFit = fitOnce from the same start (χ²_d bit-identical)",`${same}/${probs.length}; normal kept and identical ${ident}`);
  const j=G.stagedJob(probs[2],{}); ok(j.staged&&j.cov&&j.summ&&j.stats&&["staged","normal"].includes(j.staged.used)&&j.staged.stages.length===3,"stagedJob returns the finalize() summary + staged report (Worker)",`${j.staged.used}, stages χ² ${j.staged.stages.map(s=>s.chi2.toFixed(0)).join(" → ")}`);
}

/* ---------- 7 synthetic Ce 3d ---------- */
console.log("== 7 synthetic Ce 3d — recipe + staged fit vs naive Auto suggest");
const SYN={};
{ const R=G.recipeById("ce3d"); const P4={v:0.296,v2:0.405,v3:0.299}, P3={v0:0.288,v1:0.712};   // area shares within each state: one accessible 10-line decomposition (SPEC §17.3)
  function synth(f,seed){ const rg=U.rng(seed); const step=0.1, lo=870, hi=930, n=Math.round((hi-lo)/step)+1; const E=Float64Array.from({length:n},(_,i)=>lo+i*step), ctx={emin:E[0],emax:E[n-1]};
    const off=(rg.u()-0.5)*0.6, dE=18.5+(rg.u()-0.5)*0.4, A=150000;
    const truth=R.lines.map(l=>{ const fr=l.state==="Ce3+"?f*P3[l.key]:(1-f)*P4[l.key]; const c=C.newComp("GLS",l.be+off+(rg.u()-0.5)*0.3,{F:l.F*(1+(rg.u()-0.5)*0.2),m:0.3,area:A*fr}); c.doublet={dE,r:2/3,q:1}; c.rk=l.key; return c; }).filter(c=>c.area>0);
    const P=new Float64Array(n); truth.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; }); const Q=C.cumtrapz(E,P);
    const y=Float64Array.from(P,(p,i)=>rg.pois((2500+8*(E[i]-lo)+0.10*Q[i]+p)*2)/2);
    return {E,y,c:2,f,u3:truth.filter(c=>c.rk==="v3").reduce((s,c)=>s+c.area*0.4,0)/A}; }
  const nS=QUICK?2:4; const res=[];
  for(const f of [0,0.2,0.5]) for(let s=1;s<=nS;s++){ const S=synth(f,1000*f+s+7), spec={be:S.E,y:S.y,c:S.c};
    const M=G.recipeModel("ce3d",{shift:0,be:S.E,y:S.y,opt:"mix"}); const P=C.makeProblem(spec,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons});
    const st=G.stagedFit(P,{}); const fin=C.finalize(st.P,st.r); const dv=G.recipeDerived(fin,{id:"ce3d"});
    const nbR=fin.cov.atB.filter((b,j)=>b&&fin.free[j].k!=="area").length;
    // naive: what Auto suggest does (suggestModel + refineByResidual) on the same window, Smart background, then Fit
    const R0=C.sliceRegion(S.E,S.y,875,925); const sg=C.suggestModel(Float64Array.from(R0.E),Float64Array.from(R0.y),{bg:{type:"smart",n:5},el:"Ce",orb:"3d",lib:L});
    let PN=C.makeProblem(spec,{lo:875,hi:925,bg:{type:"smart",n:5},comps:sg.comps,cons:sg.cons}); const rf=C.refineByResidual(PN,{}); PN=C.makeProblem(spec,{lo:875,hi:925,bg:{type:"smart",n:5},comps:rf.prob.comps,cons:rf.prob.cons});
    const rN=C.fitOnce(PN,{}); const fN=C.finalize(PN,rN); const nbN=fN.cov.atB.filter((b,j)=>b&&fN.free[j].k!=="area").length;
    const dN=G.diagnose({E:PN.E,y:PN.y,fin:fN,el:"Ce",orb:"3d",bgType:"smart"});
    res.push({f,s,d:dv[0].v-100*f,sig:dv[0].sig,du:dv[1].v-100*S.u3,chi:fin.stats.chi2nu,nbR,nbN,nN:PN.comps.length,chiN:fN.stats.chi2nu,used:st.staged.used,naiveRecipe:dN.some(x=>x.id==="complex"&&x.act&&x.act.type==="recipe")});
    console.log(`     f ${f} #${s}: Ce³⁺ ${dv[0].v.toFixed(1)} ± ${dv[0].sig.toFixed(1)} % (true ${100*f}) · u‴ ${dv[1].v.toFixed(2)} % (true ${(100*S.u3).toFixed(2)}) · χ²ν ${fin.stats.chi2nu.toFixed(2)} · ${st.staged.used} · bound hits ${nbR} | naive ${PN.comps.length} lines χ²ν ${fN.stats.chi2nu.toFixed(2)} bound hits ${nbN}`); }
  const mx=Math.max(...res.map(r=>Math.abs(r.d))), mb=res.reduce((s,r)=>s+r.d,0)/res.length, mu=Math.max(...res.map(r=>Math.abs(r.du)));
  SYN.maxD=mx; SYN.meanD=mb; SYN.maxDu=mu; SYN.bR=res.reduce((s,r)=>s+r.nbR,0); SYN.bN=res.reduce((s,r)=>s+r.nbN,0); SYN.n=res.length;
  ok(mx<=4,"Ce³⁺ fraction recovered within ±4 % (absolute) for true 0 · 20 · 50 %",`max |Δ| ${mx.toFixed(2)} %, mean Δ ${mb.toFixed(2)} %`);
  ok(mu<=0.3,"u‴ fraction recovered within ±0.3 %",`max |Δ| ${mu.toFixed(2)} %`);
  ok(SYN.bR<SYN.bN,"bound hits (non-area parameters) recipe + staged < naive Auto suggest",`recipe ${SYN.bR} vs naive ${SYN.bN} over ${res.length} spectra`);
  ok(res.every(r=>r.naiveRecipe),"fit check on every naive fit points to the recipe");
}

/* ---------- 8 transition-metal 2p3/2 envelope recipes (v1.9.6, SPEC §17.7) ---------- */
console.log("== 8 TM 2p3/2 envelope recipes (Biesinger 2011/2017): data = journal transcripts, constraints, synthetic recovery, fit check");
const TM={};
{ const fs=require("fs"), path=require("path"), RD=path.join(__dirname,"ref","biesinger");
  const rows=f=>fs.readFileSync(path.join(RD,f),"utf8").split("\n").filter(l=>l.trim()&&!l.startsWith("#")&&!/^Compound\|/.test(l)).map(l=>l.split("|"));
  const num=x=>x===""||x==null?null:+x;
  /* parse each journal table into {name: [{be|d, pct, F10, F20}]} (column layouts as in the transcript headers) */
  const T={};
  rows("b2011_ni.txt").forEach(r=>{ const pk=[]; pk.push({be:num(r[1]),pct:num(r[2]),F10:num(r[3]),F20:num(r[4])}); for(let i=5;i+4<r.length;i+=5){ if(r[i]==="") break; pk.push({be:num(r[i]),d:num(r[i+1]),pct:num(r[i+2]),F10:num(r[i+3]),F20:num(r[i+4])}); } T["ni|"+r[0]]=pk; });
  rows("b2011_co.txt").forEach(r=>{ const pk=[]; pk.push({be:num(r[1]),pct:num(r[2]),F10:num(r[3]),F20:num(r[4])}); for(let i=5;i+4<r.length;i+=5){ if(r[i]==="") break; pk.push({be:num(r[i]),d:num(r[i+1]),pct:num(r[i+2]),F10:num(r[i+3]),F20:num(r[i+4])}); } T["co|"+r[0]]=pk; });
  rows("b2011_fe.txt").forEach(r=>{ const pk=[]; pk.push({be:num(r[1]),F10:num(r[2]),pct:num(r[3]),F20:null}); for(let i=4;i+3<r.length;i+=4){ if(r[i]==="") break; pk.push({be:num(r[i]),d:num(r[i+1]),F10:num(r[i+2]),pct:num(r[i+3]),F20:null}); } T["fe|"+r[0]]=pk; });
  { const L=fs.readFileSync(path.join(RD,"b2011_cr_mn.txt"),"utf8").split("\n"); let tb="cr";
    L.forEach(l=>{ if(/Table 3/.test(l)) tb="mn"; if(!l.trim()||l.startsWith("#")||/^Compound\|/.test(l)) return; const r=l.split("|"); const nP=tb==="cr"?5:6, F10=num(r[1+3*nP-1]), F20=num(r[1+3*nP]);
      const pk=[{be:num(r[1]),pct:num(r[2]),F10,F20}]; for(let k=1;k<nP;k++){ const i=3+3*(k-1); if(r[i]===""||r[i]==null) break; pk.push({be:num(r[i]),d:num(r[i+1]),pct:num(r[i+2]),F10,F20}); } T[tb+"|"+r[0]]=pk; }); }
  rows("b2017_cu.txt").forEach(r=>{ const pk=[{be:num(r[1]),pct:num(r[2]),F10:null,F20:num(r[3]),gl:num(r[4])}]; for(let i=5;i+3<r.length;i+=4) pk.push({be:num(r[i]),pct:num(r[i+1]),d:num(r[i+2]),F10:null,F20:num(r[i+3])}); T["cu|"+r[0]]=pk; });
  const MAP={ni2p3:{Ni0:"Ni Metal (new line shape)",NiO:"NiO",NiOH2:"Ni(OH)2",NiOOH:"Gamma NiOOH"},co2p3:{Co0:"Co(0)",CoO:"CoO",CoOH2:"Co(OH)2",CoOOH:"CoOOH",Co3O4:"Co3O4"},
    fe2p3:{Fe0:"Fe(0)",FeO:"FeO",Fe2O3:"Ave. Fe2O3",FeOOH:"Ave. FeOOH",Fe3O4_2:"Fe3O4 2+",Fe3O4_3:"Fe3O4 3+"},mn2p3:{Mn0:"Mn(0)",MnO:"Mn(II) MnO",Mn2O3:"Mn(III) Mn2O3",MnOOH:"Mn(III) Manganite (MnOOH)",MnO2:"Mn(IV) MnO2"},
    cr2p3:{Cr0:"Cr(0)",Cr2O3:"Cr(III) Oxide",CrOH3:"Cr(III) hydroxide",Cr6:"Cr(VI) oxide"},cu2p3:{Cu0:"Cu(0)",Cu2O:"Cu2O",CuO:"CuO",CuOH2:"Cu(OH)2"}};
  const ENV=G.RECIPES.filter(r=>r.kind==="env"); let nVal=0, bad=[]; const incons=[], stepDev={n:0,max:0,at:""};
  ENV.forEach(R=>{ R.species.forEach(s=>(s.sets||[s]).forEach(st=>{ const key=R.id.slice(0,2)+"|"+MAP[R.id][st.key]; const t=T[key]; if(!t){ bad.push("missing "+key); return; }
      if(t.length!==st.peaks.length) bad.push(`${key}: ${st.peaks.length} vs ${t.length} peaks`);
      st.peaks.forEach((r,j)=>{ const q=t[j]; if(!q) return; const pos=(R.abs||j===0)?q.be:q.d; const cmp=[[r[0],pos],[r[1],q.pct],[r[2],q.F10],[r[3],q.F20]];
        cmp.forEach(([a,b])=>{ nVal++; if(!(a===b||(a==null&&b==null))) bad.push(`${key} P${j+1}: ${a} vs ${b}`); });
        if(j>0&&q.be!=null&&t[j-1].be!=null&&q.d!=null){ const dd=Math.abs(+(q.be-t[j-1].be).toFixed(3)-q.d); stepDev.n++; if(dd>stepDev.max){ stepDev.max=+dd.toFixed(3); stepDev.at=`${key} P${j+1}`; } if(dd>0.1+1e-9) incons.push(`${key} P${j+1}: BE step ${+(q.be-t[j-1].be).toFixed(3)} vs printed spacing ${q.d} (beyond 0.1 eV rounding)`); } });
      const sum=t.reduce((a,q)=>a+q.pct,0); if(!st.key.startsWith("Fe3O4")&&Math.abs(sum-100)>0.35) incons.push(`${key}: % sum ${+sum.toFixed(1)}`); })); });
  { const py=T["mn|Mn(IV) Pyrolusite (MnO2)"]; const sm=py.reduce((a,q)=>a+q.pct,0); if(Math.abs(sm-100)>0.35) incons.push(`mn|Mn(IV) Pyrolusite (MnO2): % sum ${+sm.toFixed(1)} (not used)`); }
  const f34=T["fe|Fe3O4 2+"].reduce((a,q)=>a+q.pct,0), f33=T["fe|Fe3O4 3+"].reduce((a,q)=>a+q.pct,0); incons.push(`fe|Fe3O4: Fe²⁺ : Fe³⁺ = ${+f34.toFixed(1)} : ${+f33.toFixed(1)} (stoichiometric 33.3 : 66.7)`);
  [["CuO",1.89],["Cu(OH)2",1.57]].forEach(([n,r])=>{ const t=T["cu|"+n]; const sat=t.filter(q=>q.be>=938).reduce((a,q)=>a+q.pct,0), main=t.filter(q=>q.be<938).reduce((a,q)=>a+q.pct,0); incons.push(`cu|${n}: Table 4 main/shake-up ${(main/sat).toFixed(2)} vs quoted A1s/Bs ${r}`); });
  incons.unshift(`BE steps vs printed spacings: max |Δ| ${stepDev.max} eV (${stepDev.at}) over ${stepDev.n} — within the 0.1 eV rounding of the BEs`); TM.incons=incons; TM.stepDev=stepDev; incons.forEach(x=>console.log("     table note: "+x));
  ok(!bad.length&&nVal>300,`recipe numbers = journal transcripts (test/ref/biesinger, ${nVal} values: BE/spacing · % · FWHM PE10 · PE20)`,bad.slice(0,4).join(" | "));
  // shapes and positions
  const sh=(id,k)=>{ const M=G.recipeModel(id,{species:[k],pe:20}); return M.comps; };
  const ni=sh("ni2p3","Ni0")[0], fe=sh("fe2p3","Fe0")[0], cu0=sh("cu2p3","Cu0")[0], cu1=sh("cu2p3","Cu2O")[0], no=sh("ni2p3","NiO")[0];
  ok(ni.shape==="LA"&&ni.alpha===1.1&&ni.beta===2.2&&ni.casaN===10&&fe.shape==="LA"&&fe.alpha===1.2&&fe.beta===4.8&&fe.casaN===3&&cu0.shape==="GLP"&&Math.abs(cu0.m-0.9)<1e-12&&cu1.shape==="GLP"&&Math.abs(cu1.m-0.8)<1e-12&&no.shape==="GLP"&&Math.abs(no.m-0.3)<1e-12,
    "line shapes: Ni metal LA(1.1,2.2,10) · Fe metal LA(1.2,4.8,3) · Cu(0) GL(90) · Cu₂O GL(80) · others GL(30) (casaToComp)");
  const cw=c=>{ const r=C.compToCasa(c); return r.fwhm; };
  ok(Math.abs(cw(ni)-0.95)<2e-3&&Math.abs(cw(no)-1.02)<1e-3&&Math.abs(cw(cu0)-0.83)<1e-3,"FWHMs as CasaXPS 'fwhm' numbers (asymmetric LA: Casa width rule) — Ni metal 0.95 · NiO P1 1.02 · Cu(0) 0.83 (PE 20)",`${cw(ni).toFixed(4)} · ${cw(no).toFixed(4)} · ${cw(cu0).toFixed(4)}`);
  const M10=G.recipeModel("ni2p3",{species:["NiO","NiOOH"],pe:10}), M20=G.recipeModel("ni2p3",{species:["NiO","NiOOH"],pe:20});
  const w=(M,k)=>cw(M.comps.find(c=>c.rk===k));
  ok(Math.abs(w(M10,"NiO_1")-2.04)<1e-3&&Math.abs(w(M20,"NiO_1")-3.25)<1e-3&&Math.abs(w(M20,"NiOOH_1")-6.50)<1e-3&&M20.recipe.pe==="10/20","pass energy picks the FWHM column (NiO P2: PE 10 2.04 · PE 20 3.25); a peak with one column uses it (γ-NiOOH, PE 10 only)",`PE used ${M20.recipe.pe}`);
  const pos=G.recipeModel("ni2p3",{species:["NiOH2"],pe:20}).comps.map(c=>+c.be.toFixed(2)).join(" ");
  ok(pos==="854.9 855.67 857.69 860.48 861.48 866.44","positions = first BE + published spacings (Ni(OH)₂; BEs significant to 0.1 eV, spacings to 0.01 — Biesinger)",pos);
  // constraint structure
  let structOK=true; const info=[];
  ENV.forEach(R=>{ const all=R.species.map(s=>s.key); const M=G.recipeModel(R.id,{species:all,pe:20,shift:0.4}); const P={comps:M.comps,cons:M.cons,bg:M.bg};
    const map=C.buildMap(P); const sets=new Set(M.comps.map(c=>c.rset)); const fr=map.free.filter(f=>f.ci>=0);
    const per={}; fr.forEach(f=>{ const c=M.comps[f.ci]; (per[c.rset]=per[c.rset]||[]).push(f.k); });
    const okSets=[...sets].every(k=>{ const v=(per[k]||[]).sort().join(); return v==="area,be"||(R.stoich&&k==="Fe3O4_3"&&v===""); });
    C.applyLinks(P,map); const rel=M.comps.every(c=>{ const A=M.comps.find(q=>q.rset===c.rset&&q.ranc); const S=envSet(R,c.rset); const j=+c.rk.split("_").pop(); return Math.abs((c.be-A.be)-(S.be[j]-S.be[S.jA]))<1e-6&&Math.abs(c.area/A.area-S.pct[j]/S.pct[S.jA])<1e-9; });
    const win=M.comps.filter(c=>c.ranc&&!(R.stoich&&c.rset==="Fe3O4_3")).every(c=>{ const q=M.cons[c.id].be; return Math.abs(q.hi-q.lo-2*R.beTol)<1e-6&&Math.abs((q.lo+q.hi)/2-c.be)<1e-6; });
    const fixedW=fr.every(f=>f.k!=="F");
    const Mw=G.recipeModel(R.id,{species:all,pe:20,wscale:true}); const mw=C.buildMap({comps:Mw.comps,cons:Mw.cons,bg:Mw.bg}); const nF=mw.free.filter(f=>f.k==="F").length;
    if(!(okSets&&rel&&win&&fixedW&&nF===1&&Math.abs(M.lo-(R.window[0]-0.4))<1e-9)) { structOK=false; info.push(R.id); } info.push(`${R.el} ${sets.size} sets/${M.comps.length} peaks/${fr.length} free`); });
  function envSet(R,key){ const st=G.envSets(R,R.species.map(s=>s.key),20).find(s=>s.key===key); const jA=st.peaks.reduce((a,p,j)=>p.pct>st.peaks[a].pct?j:a,0); return {be:st.peaks.map(p=>p.be),pct:st.peaks.map(p=>p.pct),jA}; }
  ok(structOK,"per set exactly one free position (± tolerance, centred on literature − Δ) and one free area; spacings and area ratios held by links; widths fixed, one free width with 'FWHM×'; window anchored to Δ",info.join(" · "));
  { const M=G.recipeModel("fe2p3",{species:["Fe3O4"],stoich:true}); const P={comps:M.comps,cons:M.cons,bg:M.bg}; const map=C.buildMap(P); M.comps.find(c=>c.rset==="Fe3O4_2"&&c.ranc).area=777; C.applyLinks(P,map);
    const t2=M.comps.filter(c=>c.rset==="Fe3O4_2").reduce((a,c)=>a+c.area,0), t3=M.comps.filter(c=>c.rset==="Fe3O4_3").reduce((a,c)=>a+c.area,0);
    const a2=M.comps.find(c=>c.rset==="Fe3O4_2"&&c.ranc), a3=M.comps.find(c=>c.rset==="Fe3O4_3"&&c.ranc);
    const nfc=map.free.filter(q=>q.ci>=0).length; ok(Math.abs(t3/t2-2)<1e-9&&Math.abs(a3.be-a2.be-1.8)<1e-9&&nfc===2,"Fe₃O₄ option 1:2: Fe³⁺ set total = 2 × Fe²⁺ set total, separation 1.8 eV, one position + one area for the compound",`ratio ${(t3/t2).toFixed(6)}, free ${nfc}`); }
  ok(["Ni","Co","Fe","Mn","Cr","Cu"].every(el=>{ const ci=G.complexInfo(el,"2p"); return ci&&ci.kind==="satellite"&&ci.recipes[0]===el.toLowerCase()+"2p3"; }),"complexInfo: Cr–Cu 2p → satellite region with the envelope recipe attached (Fit check R9 points to it)");

  /* synthetic mixtures built from the recipe itself (per-set offsets ±0.1 eV, Shirley c0 + k∫P, Poisson), fitted with the recipe (defaults ∪ truth species) + staged fit */
  const MIX={ni2p3:[{Ni0:.5,NiO:.3,NiOH2:.2},{NiO:.6,NiOH2:.4},{Ni0:.6,NiOOH:.4}],
    co2p3:[{Co0:.5,Co3O4:.3,CoOH2:.2},{Co0:.6,CoO:.4},{Co3O4:.7,CoOOH:.3}],
    fe2p3:[{Fe0:.5,Fe2O3:.5},{Fe0:.4,FeO:.3,Fe2O3:.3},{Fe0:.5,Fe3O4:.5}],
    mn2p3:[{Mn0:.5,MnO:.5},{MnO:.6,MnO2:.4},{Mn0:.4,Mn2O3:.3,MnO2:.3}],
    cr2p3:[{Cr0:.5,Cr2O3:.5},{Cr0:.4,Cr2O3:.3,CrOH3:.3},{Cr2O3:.7,Cr6:.3}],
    cu2p3:[{Cu0:.5,CuO:.5},{Cu2O:.5,CuOH2:.5},{Cu0:.4,CuO:.3,CuOH2:.3}]};
  const TOL=2.5, res=[]; TM.syn={};
  function synthTM(R,mix,seed,o){ o=o||{}; const keys=Object.keys(mix), rg=U.rng(seed); const T0=G.recipeModel(R.id,{species:keys,pe:20,stoich:true});
    const lo=R.window[0], hi=R.window[1], step=0.05, n=Math.round((hi-lo)/step)+1, E=Float64Array.from({length:n},(_,i)=>lo+i*step), ctx={emin:E[0],emax:E[n-1]};
    const off={}; T0.comps.forEach(c=>{ if(!(c.rset in off)) off[c.rset]=(rg.u()-0.5)*0.2; }); if(off.Fe3O4_3!=null) off.Fe3O4_3=off.Fe3O4_2;
    const tot={}; T0.comps.forEach(c=>{ tot[c.rsp]=(tot[c.rsp]||0)+c.area; }); const truth=T0.comps.map(c=>{ const d=C.clone(c); d.be+=off[c.rset]; d.area=60000*mix[c.rsp]*c.area/tot[c.rsp]; return d; });
    const P=new Float64Array(n); truth.forEach(d=>{ const v=C.compEval(d,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; }); const Q=C.cumtrapz(E,P); const y=Float64Array.from(P,(p,i)=>rg.pois(3000+0.08*Q[i]+p));
    return {E,y,truth}; }
  const ids=Object.keys(MIX);
  ids.forEach(id=>{ const R=G.recipeById(id); TM.syn[id]=[]; MIX[id].slice(0,QUICK?2:3).forEach((mix,mi)=>{ const S=synthTM(R,mix,100+mi); const keys=Object.keys(mix);
      const sp=[...new Set(G.envDefaults(R).concat(keys))]; const M=G.recipeModel(id,{shift:0,be:S.E,y:S.y,species:sp,pe:20,stoich:true});
      const P=C.makeProblem({be:S.E,y:S.y,c:1},M); const st=G.stagedFit(P,{}); const fin=C.finalize(st.P,st.r); const dv=G.recipeDerived(fin,{id,cons:st.P.cons});
      const dvM=G.recipeDerived(fin,{id});   // build-time link metadata instead of cons → same σ
      const spd=dv.filter(d=>d.id.startsWith("sp_")); const dev=spd.map(d=>d.v-100*(mix[d.id.slice(3)]||0)); const mx=Math.max(...dev.map(Math.abs));
      const z=spd.map((d,i)=>(mix[d.id.slice(3)]||0)>0&&d.sig>0?dev[i]/d.sig:NaN);   /* present species only: an absent one sits at area 0 (a bound) where σ is not defined */ const sigSame=dv.every((d,i)=>Math.abs((d.sig||0)-(dvM[i].sig||0))<1e-9*(1+(d.sig||0)));
      let shk=null; if(id==="cu2p3"){ const tsum=S.truth.reduce((a,c)=>a+c.area,0), rOf={CuO:1.89,CuOH2:1.57}; const exp=100*S.truth.reduce((a,c)=>a+(c.rsat?(1+rOf[c.rsp])*c.area:0),0)/tsum; const d=dv.find(x=>x.id==="CuShk"); shk={v:d.v,exp,tr:100*(mix.CuO||0)+100*(mix.CuOH2||0)}; }
      const diag=G.diagnose({E:P.E,y:P.y,fin,el:R.el,orb:"2p",bgType:P.bg.type,recipe:R});
      const r={id,mi,mx,chi:fin.stats.chi2nu,z,sigSame,shk,merge:diag.filter(x=>x.id==="merge").length,perPeak:diag.filter(x=>["ill","dead","neg"].includes(x.id)&&x.act&&x.act.type!=="dropSet").length,used:st.staged.used};
      res.push(r); TM.syn[id].push({mix,fit:Object.fromEntries(spd.map(d=>[d.id.slice(3),+d.v.toFixed(2)])),maxDev:+mx.toFixed(2),chi2nu:+r.chi.toFixed(2)});
      console.log(`     ${R.el} #${mi} ${keys.map(k=>k+" "+100*mix[k]).join(" + ")}: ${spd.map((d,i)=>`${d.label} ${d.v.toFixed(1)}±${d.sig.toFixed(1)}`).join(" · ")} | max |Δ| ${mx.toFixed(2)} % · χ²ν ${r.chi.toFixed(2)} · ${r.used}`+(shk?` · shake-up Cu(II) ${shk.v.toFixed(1)} % (formula on truth ${shk.exp.toFixed(1)}, species truth ${shk.tr})`:"")); }); });
  const MX=Math.max(...res.map(r=>r.mx)); TM.maxDev=MX; TM.chi=[Math.min(...res.map(r=>r.chi)),Math.max(...res.map(r=>r.chi))];
  TM.byEl=Object.fromEntries(ids.map(id=>[id,Math.max(...res.filter(r=>r.id===id).map(r=>r.mx))]));
  ok(MX<=TOL,`species fractions recovered within ±${TOL} % (absolute) on ${res.length} synthetic mixtures of 2–3 species (6 elements)`,ids.map(id=>id.slice(0,2)+" "+TM.byEl[id].toFixed(2)).join(" · "));
  ok(res.every(r=>r.chi<1.3),"every fit reaches the noise floor (χ²ν < 1.3; true model family)",`χ²ν ${TM.chi.map(x=>x.toFixed(2)).join("–")}`);
  const zs=res.flatMap(r=>r.z).filter(isFinite), zr=Math.sqrt(zs.reduce((a,b)=>a+b*b,0)/zs.length); TM.zrms=zr;
  ok(zr<3,"covariance σ of the fractions is of the right size (rms of Δ/σ over the present species < 3)",`rms z ${zr.toFixed(2)} over ${zs.length}`);
  ok(res.every(r=>r.sigSame),"fraction σ from the cons links = σ from the build-time link metadata (delta method through linked areas)");
  const cu=res.filter(r=>r.shk); TM.shk=cu.map(r=>r.shk);
  ok(cu.every(r=>Math.abs(r.shk.v-r.shk.exp)<=TOL),"Cu(II) shake-up estimate (Biesinger 2017 Eq. 2) = the formula on the true areas within ±2.5 % (it reads high vs the species truth by the Table-4 ratio mismatch)",cu.map(r=>`${r.shk.v.toFixed(1)}/${r.shk.exp.toFixed(1)}/${r.shk.tr}`).join(" · "));
  ok(res.every(r=>r.merge===0&&r.perPeak===0),"fit check on envelope fits: no line-merge findings, area findings once per species (not per peak)");
  // drop an absent species via the fit check action → links re-anchored exactly → refit keeps the fractions
  { const R=G.recipeById("ni2p3"), mix={Ni0:.5,NiO:.5}; const S=synthTM(R,mix,321); const M=G.recipeModel("ni2p3",{shift:0,be:S.E,y:S.y,species:["Ni0","NiO","NiOH2","NiOOH"],pe:20,wscale:true});
    let P=C.makeProblem({be:S.E,y:S.y,c:1},M); let r=C.fitOnce(P,{}); let fin=C.finalize(P,r); const d=G.diagnose({E:P.E,y:P.y,fin,el:"Ni",orb:"2p",bgType:P.bg.type,recipe:R});
    const f=d.find(x=>x.act&&x.act.type==="dropSet"); let msg="no dropSet finding";
    if(f){ const mm={lo:M.lo,hi:M.hi,bg:M.bg,comps:fin.comps,cons:P.cons}; const nb=mm.comps.length; G.applyAction(mm,f.act); const map=C.buildMap({comps:mm.comps,cons:mm.cons,bg:mm.bg});
      const P2=C.makeProblem({be:S.E,y:S.y,c:1},mm); const r2=C.fitOnce(P2,{}); const fin2=C.finalize(P2,r2); const dv=G.recipeDerived(fin2,{id:"ni2p3",cons:P2.cons});
      const ni0=dv.find(x=>x.id==="sp_Ni0"), nio=dv.find(x=>x.id==="sp_NiO"); const nF=map.free.filter(q=>q.k==="F").length;
      msg=`${f.label}: ${nb} → ${mm.comps.length} peaks, free F ${nF}, Ni ${ni0.v.toFixed(1)} · NiO ${nio.v.toFixed(1)} %, χ²ν ${fin2.stats.chi2nu.toFixed(2)}`;
      ok(mm.comps.length<nb&&nF===1&&Math.abs(ni0.v-50)<TOL&&Math.abs(nio.v-50)<TOL&&fin2.stats.chi2nu<1.3,"'Remove <species>' (dropSet) removes the whole envelope, re-anchors links (common FWHM scale kept) → refit unchanged",msg); }
    else ok(false,"'Remove <species>' (dropSet) expected on an absent species",d.map(x=>x.label).join(" | ")); }
  // naive fit of a Ni 2p3/2 envelope → fit check points to the recipe
  { const R=G.recipeById("ni2p3"); const S=synthTM(R,{Ni0:.5,NiO:.5},55); const sg=C.suggestModel(S.E,S.y,{bg:{type:"smart",n:5},el:"Ni",orb:"2p",lib:L});
    const P=C.makeProblem({be:S.E,y:S.y,c:1},{lo:S.E[0],hi:S.E[S.E.length-1],bg:{type:"smart",n:5},comps:sg.comps,cons:sg.cons}); const r=C.fitOnce(P,{}); const fin=C.finalize(P,r);
    const d=G.diagnose({E:P.E,y:P.y,fin,el:"Ni",orb:"2p",bgType:"smart"}); const cx=d.find(x=>x.id==="complex");
    ok(cx&&cx.act&&cx.act.type==="recipe"&&cx.act.recipe==="ni2p3"&&/envelope/.test(cx.tip),"naive Ni 2p fit → Fit check 'Apply recipe' with the envelope recipe",cx&&cx.label); }
}
console.log(`\ntest_guide: ${pass} passed, ${fail} failed (${((Date.now()-t0)/1000).toFixed(0)} s)`);
if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; }
module.exports={FP,SYN,POS,TM};
