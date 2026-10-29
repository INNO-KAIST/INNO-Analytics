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
     recipe + staged fit recovers the Ce³⁺ and u‴ fractions; bound hits vs a naive suggestModel fit */
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
  const inR=R.lines.every(l=>lit[l.key].every(b=>Math.abs(b-l.be)<=R.beTol+1e-9)), dOK=dEs.every(d=>d>=R.dE.lo&&d<=R.dE.hi);
  ok(inR&&dOK&&Math.abs(R.r-2/3)<1e-12,"literature positions within start ± 0.6 eV, ΔE 18.3–18.7 eV inside the free range, u/v area 2/3",`ΔE ${R.dE.lo}–${R.dE.hi}`);
  const M=G.recipeModel("ce3d",{shift:1.3,opt:"mix",align:false}); const c0=M.comps.find(c=>c.rk==="v3");
  ok(Math.abs(c0.be-(898.2-1.3))<1e-9&&Math.abs(M.lo-(875-1.3))<1e-9&&Math.abs(M.hi-(925-1.3))<1e-9,"anchors to the charge correction: raw BE = literature − Δ (positions and window)",`v‴ ${c0.be.toFixed(2)} · window ${M.lo.toFixed(1)}–${M.hi.toFixed(1)}`);
  const anc=M.comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); const links=M.comps.filter(c=>c!==anc).every(c=>M.cons[c.id].dE.mode==="link"&&M.cons[c.id].dE.ref===anc.id&&M.cons[c.id].m.mode==="link");
  ok(links&&M.comps.every(c=>c.doublet&&c.doublet.r===R.r&&M.cons[c.id].r.mode==="fixed"&&M.cons[c.id].q.mode==="fixed"&&M.cons[c.id].be.hi-M.cons[c.id].be.lo>1.19&&M.cons[c.id].be.hi-M.cons[c.id].be.lo<1.21),"pairs: shared ΔE (linked), area ratio and width factor fixed, BE free ± 0.6 eV");
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
console.log(`\ntest_guide: ${pass} passed, ${fail} failed (${((Date.now()-t0)/1000).toFixed(0)} s)`);
if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; }
module.exports={FP,SYN,POS};
