/* v1.10 Auto build (SPEC §18): rules one by one, honest uncertainty, recipes, and the synthetic benchmark
   node test/test_autobuild.js     (from modules/xps; QUICK=1 → benchmark subset 24 + 4 spectra, about 1–2 min; full: 400 + 40, about 15–25 min)
   1 statistics: N_eff = N(1 − ρ)/(1 + ρ), ρ = residual ρ₁ of the better model capped at 0.5, static background counts 2 parameters
   2 acceptance: add needs ΔAICc_eff and ΔBIC_eff ≥ 10; simplification kept when ΔAICc_eff > −2; automatic lines dropped by the add criterion (symmetric); recipe lines only when |ΔAICc| < 2; NaN never passes
   3 rules on synthetic spectra: decision background first · add at the residual run (and none on noise) · merge cap · min separation · shared L fraction ·
     tied widths (equal yes, 1.0 vs 2.0 no) · slope on ('slope off' disabled) · release a non-literature FWHM bound · drop a line on noise · recipe lines not tied/released
   4 final analysis: |r| > 0.9 group → sum with σ, members 'not separately determined'; χ²ν > 2 → not established; static background → not established;
     local multi-start uses χ²_d; hidden-neighbour bound covers a missed neighbour; candidates only from cited tables (no naming); input never mutated; log
   5 user-session-like O 1s (lattice + missed low-BE lattice line + surface O, VGD-like smoothing): recipe start → Auto build adds the ≈ 528.6 eV line
   6 benchmark (SPEC §18.4 acceptance criteria A–C on confirmation set 4, seeds 700001… / 800001…; QUICK: first 24 + 4, no false line in-family, runs complete);
     v1.11.1: Fit check R16 (structured residual, SPEC §17.12) on the finals — ≤ 2 % on finals with every true line and no false one */
"use strict";
const C=require("../src/xps_core.js"), G=require("../src/xps_guide.js"), A=require("../src/xps_autobuild.js"), U=require("./synth_util.js"), B=require("./autobuild_bench.js"), REP=require("./autobuild_report.js");
const QUICK=!!process.env.QUICK;
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info!=null&&info!==""?"  ("+info+")":"")); }
const t0=Date.now();
function spec(lines,o){ o=o||{}; const rg=U.rng(o.seed||5); const lo=o.lo!=null?o.lo:280, hi=o.hi!=null?o.hi:296, step=o.step||0.05;
  const n=Math.round((hi-lo)/step)+1, E=Float64Array.from({length:n},(_,i)=>lo+i*step), ctx={emin:E[0],emax:E[n-1]};
  const P=new Float64Array(n); lines.forEach(([be,Ar,F,m])=>{ const v=C.compEval(C.newComp("GLS",be,{F:F||1.2,m:m!=null?m:0.3,area:Ar}),E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P); const mu=Float64Array.from(P,(p,i)=>2000+(o.slope||0)*(E[i]-lo)+0.05*Q[i]+p); let y;
  if(o.smooth){ const x=Float64Array.from(mu,v=>rg.pois(v)); y=Float64Array.from(x,(v,i)=>i===0||i===n-1?v:x[i-1]/6+2*v/3+x[i+1]/6); } else y=Float64Array.from(mu,v=>rg.pois(v));
  return {be:E,y,c:1}; }
const ACT=(slope)=>({type:"active",n:5,c0:null,k:null,s:0,cons:{s:{mode:slope?"free":"fixed"}}});
function prob(sp,comps,cons,bg){ return C.makeProblem(sp,{lo:sp.be[0],hi:sp.be[sp.be.length-1],bg:bg||ACT(false),comps,cons:cons||{}}); }
const gl=(be,A_,F,o)=>{ const c=C.newComp("GLS",be,Object.assign({F:F||1.2,m:0.3},o||{})); c.area=A_||10000; return c; };
const freeCons=(cs,o)=>{ const q={}; cs.forEach(c=>{ q[c.id]={be:{mode:"free",lo:c.be-1.5,hi:c.be+1.5},F:{mode:"free",lo:(o&&o.Flo)||0.3,hi:(o&&o.Fhi)||4}}; if(o&&o.mfix) q[c.id].m={mode:"fixed"}; }); return q; };
const S0=P=>A.startState(P,{first:true});
const cands=(P,ctx)=>{ const S=S0(P); return {S,cd:A.candidates(S,ctx||{},new Set())}; };

console.log("== 1 statistics");
{ const st={N:400,p:5,chi2d:380,rho1:0.3,chi2nu:1}; const S={st,P:{bg:{type:"active"}}};
  const e=A.icEff(S,0.4); const Ne=400*0.6/1.4; ok(Math.abs(e.Neff-Ne)<1e-9&&Math.abs(e.AICc-(Ne*Math.log(380/400)+10+60/(Ne-6)))<1e-9&&Math.abs(e.BIC-(Ne*Math.log(380/400)+5*Math.log(Ne)))<1e-9,"N_eff = N(1 − ρ)/(1 + ρ) in AICc and BIC (Zięba 2010 Eq. 15)",Ne.toFixed(2));
  const Ss={st,P:{bg:{type:"smart"}}}; ok(A.icEff(Ss,0).p===7,"static background counts its two end levels (p + 2)");
  const a={st:Object.assign({},st,{rho1:0.9,chi2d:900}),P:{bg:{type:"active"}}}, b={st:Object.assign({},st,{rho1:0.2,chi2d:400}),P:{bg:{type:"active"}}}, c={st:Object.assign({},st,{rho1:0.8,chi2d:390}),P:{bg:{type:"active"}}};
  ok(Math.abs(A.cmpIC(a,b).rho-0.2)<1e-12&&Math.abs(A.cmpIC(b,a).rho-0.2)<1e-12&&Math.abs(A.cmpIC(a,c).rho-0.5)<1e-12,"ρ = ρ₁ of the better-fitting model (the misfit's ρ₁ 0.9 is not used), capped at 0.5"); }

console.log("== 2 acceptance rules");
{ const d=(x,y,dp)=>({dAICc:x,dBIC:y,dp:dp||1});
  const P=m=>A.passes(m);
  ok(P({cls:"add",kind:"add",d:d(10,10)})&&!P({cls:"add",kind:"add",d:d(30,9.9)})&&!P({cls:"add",kind:"add",d:d(9.9,30)}),"add: both ΔAICc_eff and ΔBIC_eff ≥ 10");
  ok(P({cls:"simp",kind:"tieF",d:d(-1.9,-5,-1)})&&!P({cls:"simp",kind:"tieF",d:d(-2.1,5,-1)}),"simplification kept when ΔAICc_eff > −2");
  ok(P({cls:"simp",kind:"drop",info:{recipe:false},d:d(-9,-30,-3)})&&!P({cls:"simp",kind:"drop",info:{recipe:false},d:d(-10,-10,-3)})&&P({cls:"simp",kind:"drop",info:{recipe:false},d:d(-12,-9,-3)}),"automatic line dropped unless removal costs ≥ 10 in both (symmetric with add)");
  ok(P({cls:"simp",kind:"drop",info:{recipe:true},d:d(-1.5,-9,-3)})&&!P({cls:"simp",kind:"drop",info:{recipe:true},d:d(-2.5,9,-3)}),"recipe line dropped only when |ΔAICc| < 2");
  ok(!P({cls:"add",kind:"add",d:d(NaN,20)})&&!P({cls:"simp",kind:"tieF",d:d(NaN,1)})&&!P({cls:"add",kind:"add",d:d(20,20),capFail:"x"}),"NaN never passes; a cap failure never passes"); }

console.log("== 3 rules on synthetic spectra");
{ // 3a decision background first
  const sp=spec([[285,30000,1.2]],{seed:11,slope:-40}); const c=gl(285.1,25000); const S=S0(prob(sp,[c],freeCons([c]),{type:"smart",n:5})); const s=A.stepFrom(S,{},new Set());
  ok(s.best&&s.best.kind==="bgAct"&&s.best.S2.P.bg.type==="active"&&s.best.S2.P.bg.cons.s.mode==="free","static background → first step: in-fit Shirley + slope (decision background)",s.best&&s.best.label); }
{ // 3b add at the residual run; none on noise
  const sp=spec([[285,30000,1.2],[287.2,6000,1.2]],{seed:12}); const c=gl(285,30000); const R=A.autoBuild(prob(sp,[c],freeCons([c])),{},{final:false,multi:false});
  const st=R.steps.find(x=>x.kind==="add"); const nc=R.S.P.comps.find(q=>q.ab);
  ok(st&&nc&&Math.abs(nc.be-287.2)<0.1&&nc.name===nc.be.toFixed(1)+" eV line"&&st.ev.dAICc>=10&&st.ev.dBIC>=10&&R.S.P.comps.length===2,"missing 20 % line → 'Add a line' (both criteria ≥ 10), final line at 287.2 eV named by its fitted BE, two lines",R.steps.map(x=>x.label.slice(0,40)).join(" | ")+(nc?" → "+nc.name:""));
  let adds=0; for(let k=0;k<(QUICK?6:16);k++){ const sp2=spec([[285,30000,1.2]],{seed:200+k}); const c2=gl(285,30000); const R=A.autoBuild(prob(sp2,[c2],freeCons([c2]),ACT(false)),{},{final:false,multi:false}); if(R.S.P.comps.length>1) adds++; }
  ok(adds===0,"one-line spectra (Poisson): no line added",`${adds} of ${QUICK?6:16}`); }
{ // 3c merge cap and min separation
  const sp=spec([[285,30000,1.2]],{seed:13}); const a=gl(285.0,15000), b=gl(285.15,15000); const cs={}; cs[a.id]={be:{mode:"free",lo:284,hi:286}}; cs[b.id]={be:{mode:"free",lo:284,hi:286}};
  const S=S0(prob(sp,[a,b],cs)); const s=A.stepFrom(S,{},new Set()); ok(s.best&&(s.best.kind==="merge"||s.best.kind==="dead"||s.best.kind==="drop")&&s.best.S2.P.comps.length===1,"two lines on one peak → merge/remove leaves one line",s.best&&s.best.label);
  const sp2=spec([[285,30000,1.2],[285.45,9000,1.2]],{seed:14}); const c=gl(285.1,39000); const {cd}=cands(prob(sp2,[c],freeCons([c])));
  ok(cd.list.filter(m=>m.kind==="add").every(m=>Math.abs(m.info.be-c.be)>=0.5*1.2-0.15),"no add candidate closer than 0.5 × FWHM to an existing line",cd.list.filter(m=>m.kind==="add").map(m=>m.info.be.toFixed(2)).join(" ")); }
{ // 3d shared L fraction; tied widths
  const sp=spec([[285,30000,1.2,0.3],[287.5,15000,1.2,0.3]],{seed:15}); const a=gl(285,30000), b=gl(287.5,15000); const q=freeCons([a,b]);
  const {cd}=cands(prob(sp,[a,b],q)); const sm=cd.list.find(m=>m.kind==="shareM"), ta=cd.list.find(m=>m.key==="tieF|all");
  ok(sm&&A.passes(sm)&&sm.S2.P.cons[b.id].m.mode==="link","free per-line L fractions, common truth → 'one shared L fraction' passes (ΔAICc_eff > −2)",sm&&sm.d.dAICc.toFixed(2));
  ok(ta&&A.passes(ta)&&ta.S2.P.cons[b.id].F.k===1,"equal widths → 'tie all FWHMs' (k = 1) passes",ta&&ta.d.dAICc.toFixed(2));
  const sp2=spec([[285,30000,1.0],[288,15000,2.0]],{seed:16}); const c=gl(285,30000,1.0), d_=gl(288,15000,2.0); const {cd:cd2}=cands(prob(sp2,[c,d_],freeCons([c,d_])));
  const t2=cd2.list.find(m=>m.key==="tieF|all"); ok(t2&&!A.passes(t2),"1.0 vs 2.0 eV widths → tying is rejected",t2&&t2.d.dAICc.toFixed(1)); }
{ // 3e slope: 'slope off' is disabled (SPEC §18.5 — held-out set 1: involved in 5 of 8 false lines); 'slope on' passes on sloped truth
  const sp=spec([[285,30000,1.2]],{seed:17}); const c=gl(285,30000); const {cd}=cands(prob(sp,[c],freeCons([c]),ACT(true)));
  ok(A.AB_TH.slopeOff===false&&!cd.list.some(m=>m.kind==="bgOff"),"'slope off' simplification is disabled (no candidate)");
  const sp2=spec([[285,30000,1.2]],{seed:18,slope:-150}); const c2=gl(285,30000); const {cd:cd2}=cands(prob(sp2,[c2],freeCons([c2]),ACT(false))); const on=cd2.list.find(m=>m.kind==="bgOn");
  ok(on&&A.passes(on),"sloped truth: 'slope on' passes (≥ 10 both)",on&&on.d.dAICc.toFixed(1)); }
{ // 3f release a non-literature FWHM bound
  const sp=spec([[285,30000,1.5]],{seed:19}); const c=gl(285,30000,0.9); const S=S0(prob(sp,[c],freeCons([c],{Flo:0.3,Fhi:0.9}))); const s=A.stepFrom(S,{},new Set());
  ok(s.best&&s.best.kind==="release"&&s.best.S2.P.cons[c.id].F.hi>1.3,"FWHM at a user upper bound below the truth → release (× 1.5) passes",s.best&&s.best.label); }
{ // 3g drop a line on noise (Auto suggest-like extra line)
  const sp=spec([[285,30000,1.2]],{seed:20}); const a=gl(285,30000), b=gl(290.5,300); const q=freeCons([a,b]); const S=S0(prob(sp,[a,b],q)); const s=A.stepFrom(S,{},new Set());
  ok(s.best&&(s.best.kind==="drop"||s.best.kind==="dead")&&s.best.S2.P.comps.length===1,"an extra line on flat background is dropped (cost < 10)",s.best&&s.best.label); }
{ // 3h recipe lines: not tied, not released, no shared-m change; C 1s
  // truth in the recipe's own shape (CasaXPS GL(30) product); a GL-sum truth at these counts makes the tail mismatch look like a line (SPEC §18.5)
  const sp=(()=>{ const rg=U.rng(21); const E=Float64Array.from({length:361},(_,i)=>278+i*0.05), n=E.length, cx={emin:E[0],emax:E[n-1]}; const P_=new Float64Array(n);
    [[284.8,30000],[286.3,7000],[288.8,3000]].forEach(([b,a])=>{ const v=C.compEval(C.casaToComp("GL(30)",{be:b,fwhm:1.3,area:a}).comp,E,cx); for(let i=0;i<n;i++) P_[i]+=v[i]; });
    const Q=C.cumtrapz(E,P_); return {be:E,y:Float64Array.from(P_,(p,i)=>rg.pois(2000+0.05*Q[i]+p)),c:1}; })(); const M=G.recipeModel("c1s",{shift:0,be:sp.be,y:sp.y,present:[]});
  const P=C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons}); const f=G.c1sFit(P,{}); const P2=C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:f.P.bg,comps:f.P.comps,cons:f.P.cons});
  const S=A.startState(P2,{}); const cd=A.candidates(S,{el:"C",orb:"1s"},new Set()); const ks=cd.list.map(m=>m.kind);
  ok(!ks.includes("tieF")&&!ks.includes("shareM")&&!cd.list.some(m=>m.kind==="release"),"recipe model: no tie / shared-m / release candidates on recipe lines (literature constraints kept)",ks.join(","));
  const R=A.autoBuild(P2,{el:"C",orb:"1s"},{first:false}); ok(R.S.P.comps.filter(c=>c.recipe).length>=3&&R.S.P.comps.filter(c=>c.ab).length===0,"C 1s recipe start (three lines in truth): Auto build keeps the recipe lines and adds none",R.S.P.comps.map(c=>c.name).join(", ")); }

console.log("== 4 final analysis");
{ // 4a correlated pair → group with sum
  const sp=spec([[285,20000,1.2],[285.75,20000,1.2]],{seed:31}); const a=gl(285,20000), b=gl(285.75,20000); const S=S0(prob(sp,[a,b],freeCons([a,b])));
  const F=A.finalAnalysis(S,{},{starts:8}); const g=F.groups[0];
  ok(g&&g.ids.length===2&&F.lines.every(l=>!l.separate)&&g.sig>0&&Math.abs(g.area-40000)<5*g.sig,"two lines 0.6 FWHM apart (|r| > 0.9): one group, sum ± σ, members 'not separately determined', sum within 5σ of 40 000",g&&`${g.area.toFixed(0)} ± ${g.sig.toFixed(0)} r ${g.pairs.map(q=>q.r.toFixed(3)).join(",")}`); }
{ // 4b wrong model χ²ν > 2 → not established
  const sp=spec([[285,30000,1.2],[288,20000,1.2]],{seed:32}); const a=gl(285,30000); const q={}; q[a.id]={be:{mode:"free",lo:284,hi:286}}; const S=S0(prob(sp,[a],q));
  const F=A.finalAnalysis(S,{},{starts:8}); ok(!F.established&&F.reasons.some(r=>/χ²ν/.test(r)),"missing 40 % line (χ²ν ≫ 2) → model not established",F.reasons.join("; ")); }
{ // 4c static background → not established
  const sp=spec([[285,30000,1.2]],{seed:33}); const a=gl(285,30000); const S=S0(prob(sp,[a],freeCons([a]),{type:"smart",n:5})); const F=A.finalAnalysis(S,{},{starts:6});
  ok(!F.established&&F.reasons.some(r=>/static/.test(r)),"static background kept → not established (σ excludes the background error)"); }
{ // 4d multi-start: local starts, χ²_d — clean single line reaches the same minimum
  const sp=spec([[285,30000,1.2]],{seed:34,slope:-80}); const a=gl(285.1,25000); const S=S0(prob(sp,[a],freeCons([a]),ACT(true))); const ms=A.multiStart(S,{starts:12});
  ok(ms.frac>=0.75&&Math.abs(ms.cur-S.st.chi2d)<1e-9,"clean one-line spectrum: local multi-start reaches the same χ²_d minimum from ≥ 3/4 of 12 starts",`${ms.near}/${ms.N}`); }
{ // 4e hidden-neighbour bound covers a missed neighbour (sub-threshold 0.7 FWHM neighbour, 12 %)
  let cov=0, n=0, strict=0; for(let k=0;k<(QUICK?4:10);k++){ const sp=spec([[285,30000,1.2],[285.85,3600,1.2]],{seed:400+k}); const a=gl(285.2,33000); const S=S0(prob(sp,[a],freeCons([a]))); const F=A.finalAnalysis(S,{},{multi:false,support:false});
    const L=F.lines[0]; const tb=(285*30000+285.85*3600)/33600; n++; if(Math.abs(L.be-285)<=3*L.beU) cov++; if(Math.abs(L.be-285)<=3*L.beSig) strict++; }
  ok(cov===n,"missed 12 % neighbour at 0.7 FWHM: the true main-line position lies within 3u (u = σ + hidden-line bound/3) in every case",`3u ${cov}/${n}, 3σ only ${strict}/${n}`); }
{ // 4f candidates: cited tables only, filtered by the elements present, never used as names
  const c1=A.candidateSpecies("O","1s",531.5,["Co","O","C"]); const c2=A.candidateSpecies("O","1s",529.5,["Co","O"]); const c3=A.candidateSpecies("O","1s",529.5,["Ti"]);
  ok(c1.some(x=>/OH/.test(x.name)&&/Graf/.test(x.src))&&c2.some(x=>/lattice O \(Co\)/.test(x.name)&&/Graf/.test(x.src))&&!c3.some(x=>/\(Co\)/.test(x.name))&&A.candidateSpecies("Zz","1s",100,[]).length===0,"candidates: recipe table entries with source, lattice O only for metals present, nothing for an unknown region",c1.map(x=>x.name).join(", ")); }
{ // 4g the input problem is never changed; log
  const sp=spec([[285,30000,1.2],[287.2,6000,1.2]],{seed:35}); const c=gl(285,30000); const P=prob(sp,[c],freeCons([c]),{type:"smart",n:5}); const before=JSON.stringify({c:P.comps,q:P.cons,b:P.bg});
  const R=A.autoBuild(P,{},{starts:6}); ok(JSON.stringify({c:P.comps,q:P.cons,b:P.bg})===before&&R.steps.length>=2,"autoBuild works on copies: the input model is unchanged",`${R.steps.length} steps`);
  const md=A.logText({region:"C 1s",steps:R.steps.map(s=>Object.assign({},s,{decision:"accepted"})),stop:R.stop,final:R.final});
  ok(/\| # \| decision \|/.test(md)&&/Stopped:/.test(md)&&/## Result:/.test(md)&&R.steps.every(s=>md.includes(s.label.slice(0,20))),"step log (Markdown): every step, the stop reason and the result"); }

console.log("== 5 user-session-like O 1s (recipe start; missed low-BE lattice line)");
{ const rg=U.rng(77); const E=Float64Array.from({length:241},(_,i)=>525.0+i*0.05), n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const L=[[529.3,30000,1.2],[528.4,6000,1.2],[531.5,16000,1.8],[532.1,4000,1.8],[532.6,3000,1.8]]; const P=new Float64Array(n); L.forEach(([b,a,F])=>{ const c=C.casaToComp("GL(30)",{be:b,fwhm:F,area:a}).comp; const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P); const x=Float64Array.from(P,(p,i)=>rg.pois(1500+0.04*Q[i]+p)); const y=Float64Array.from(x,(v,i)=>i===0||i===n-1?v:x[i-1]/6+2*v/3+x[i+1]/6);
  const sp={be:E,y,c:1}; const M=G.recipeModel("o1s",{shift:0,be:E,y,present:["Ba","Co","Ta"],flags:{lat2:false}}); const P0=C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons});
  const f=G.o1sFit(P0,{build:{shift:0,present:["Ba","Co","Ta"],flags:{lat2:false}}}); const P1=C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:f.P.bg,comps:f.P.comps,cons:f.P.cons});
  const Pc=C.clone(C.stripProb(P1)); Pc.E=P1.E; Pc.y=P1.y; const st1=C.fitStats(Pc,C.fitOnce(Pc,{initAreas:false}));
  const R=A.autoBuild(P1,{el:"O",orb:"1s",present:["Ba","Co","Ta"]},{first:false}); const add=R.S.P.comps.find(c=>c.ab);
  ok(add&&Math.abs(add.be-528.4)<0.3&&R.final.chi2nu<st1.chi2nu,"recipe start without 'Low lattice': Auto build adds a line near 528.4 eV (truth), χ²ν falls",`${add?add.name:"none"} · χ²ν ${st1.chi2nu.toFixed(3)} → ${R.final.chi2nu.toFixed(3)} · ${R.final.established?"established":"not established: "+R.final.reasons.join("; ")}`);
  const sp_=R.final.species.find(s=>add&&s.id===add.id); ok(sp_&&sp_.cands.some(c=>/lattice/i.test(c.name))&&/eV line$/.test(add.name),"the added line gets candidates only (lattice O entries with sources), not a name",sp_&&sp_.cands.map(c=>c.name).join(", ")); }

console.log("== 6 benchmark (SPEC §18.4)");
{ const nG=QUICK?24:400, nR=QUICK?4:40; const R=[]; const tb=Date.now();
  // seeds = confirmation set 3 (SPEC §18.5: set 1 = 100001/200001 pre-registered and failed (A) → 'slope off' disabled; set 2 = 300001/400001; set 3 = 500001/600001 after the FWHM-at-bound check; set 4 = 700001/800001 after restricting it to the upper bound — used here)
  const G0=700001, R0=800001;
  for(let s=G0;s<G0+nG;s++) R.push(B.runCase(B.gen(s)));
  for(let s=R0;s<R0+nR;s++) R.push(B.runRecipeCase(B.genRecipe(s)));
  const errs=R.filter(r=>r.err); ok(!errs.length,`benchmark: ${R.length} spectra ran without error`,errs.slice(0,3).map(r=>r.seed+": "+r.err).join("; "));
  const sym=R.filter(r=>!r.meta.asym&&!r.meta.recipe), rec=R.filter(r=>r.meta.recipe); const S=REP.summarize(sym), SR=REP.summarize(rec);
  console.log(REP.report(R).split("\n").slice(0,QUICK?6:200).join("\n"));
  if(QUICK){ ok(S.fa===0,"QUICK subset: no false line in the in-family spectra",`${S.fa}/${S.n}`); ok(SR.fa===0,"QUICK subset: no false line added on recipe starts",`${SR.fa}/${SR.n}`); }
  else {
    const bins=[["sepBin",["single","0.3-0.5","0.5-1","1-1.5","1.5-3"]],["noise",["poisson","ar1","smooth"]],["H",[300,3000,30000]],["K",[1,2,3,4,5,6]],["slope",[true,false]]];
    let worst={v:0,b:""}; bins.forEach(([k,vs])=>vs.forEach(v=>{ const a=sym.filter(r=>r.meta[k]===v); if(!a.length) return; const q=REP.summarize(a); const f=q.fa/q.n; if(f>worst.v) worst={v:f,b:`${k} ${v} (${q.fa}/${q.n})`}; }));
    ok(S.fa/S.n<=0.02&&worst.v<=0.05&&SR.fa/SR.n<=0.02,"(A) false-add ≤ 2 % in-family, ≤ 5 % in every bin, ≤ 2 % on recipe starts",`${S.fa}/${S.n}; worst bin ${worst.b||"none"}; recipe ${SR.fa}/${SR.n}`);
    const all=REP.summarize(R.filter(r=>!r.meta.asym)); ok(all.posOK/all.estPos>=0.95,"(B) established cases: all positions within 3u in ≥ 95 %",`${all.posOK}/${all.estPos} (strict 3σ ${all.posStrict}/${all.estPos})`);
    ok(all.unresOK/Math.max(1,all.unres)>=0.95,"(C) unresolvable truth pairs handled honestly in ≥ 95 %",`${all.unresOK}/${all.unres}`); }
  // v1.11.1 Fit check R16 (구조가 있는 잔차, SPEC §17.12) 를 Auto build 끝 해에: 참 선을 모두 찾고 가짜 선이 없는 끝 해에서의 경고율 ≤ 2 % ·
  // '확립됨'인데 R16 이 걸린 것은 무엇을 놓쳤나 (실측 사각지대 — 끝 판정은 χ²ν > 2 또는 ρ₁ > 0.5)
  { const sy=R.filter(r=>r.fc&&!r.meta.asym), comp=sy.filter(r=>r.ev.missed===0&&r.ev.falseAdd===0), fA=comp.filter(r=>r.fc.struct);
    const est=sy.filter(r=>r.ev.established&&r.fc.struct), mis=sy.filter(r=>r.ev.missed>0), asy=R.filter(r=>r.fc&&r.meta.asym);
    console.log(`     R16 on Auto build finals: complete ${fA.length}/${comp.length}`+(fA.length?` (${fA.map(r=>`${r.seed} χ²ν ${r.fc.chi2nu.toFixed(2)} ρ₁ ${r.fc.rho1.toFixed(2)}${r.ev.established?"":" not established"}`).join(", ")})`:"")+` · established with R16 ${est.length} (missed lines ${est.map(r=>r.ev.missed).join(",")||"—"}) · finals that missed a line ${mis.filter(r=>r.fc.struct).length}/${mis.length} · asymmetric (LA, outside the model set) ${asy.filter(r=>r.fc.struct).length}/${asy.length}`);
    ok(fA.length<=0.02*comp.length,"Fit check R16 (structured residual) on Auto build finals that found every true line and no false one: ≤ 2 %",`${fA.length}/${comp.length}`); }
  console.log(`  benchmark time ${((Date.now()-tb)/1000).toFixed(0)} s`); }

console.log(`\ntest_autobuild: ${pass} passed, ${fail} failed (${((Date.now()-t0)/1000).toFixed(0)} s)`);
if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; }
