/* Auto build 합성 벤치마크 (SPEC §18.4–18.5): 참값을 아는 스펙트럼 생성 · Auto suggest → Auto build · 참값 대조.
   node test/autobuild_bench.js <seed0> <count> [out.jsonl]   — 한 줄에 한 스펙트럼의 결과(JSON). test_autobuild.js 가 이 모듈을 씀.
   설계(벤치마크 전에 SPEC §18.4 에 적음): 선 1–6 · 이웃 간격 0.3–3 FWHM (구간 0.3–0.5 · 0.5–1 · 1–1.5 · 1.5–3) · 넓이 비 0.03–1 (로그 균등) ·
   GL 합 L 비율 0–0.6(스펙트럼 공통) · 선마다 FWHM ±10 % · 10 % 는 가장 큰 선이 비대칭(LA, 모형 집합 밖) · 배경 = 상수 + Shirley 계단(k·∫P) ± 기울기 ·
   봉우리 높이 300 · 3000 · 30000 count · 잡음: Poisson 60 % · AR(1)(ρ 0.3–0.5) 20 % · Poisson 뒤 [1/6, 2/3, 1/6] 평활(VGD 처럼 χ²ν ≈ 0.5, ρ₁ ≈ 0.44) 20 %. */
"use strict";
const path=require("path");
const C=require("../src/xps_core.js"), A=require("../src/xps_autobuild.js"), U=require("./synth_util.js");
const SEPBIN=[[0.3,0.5],[0.5,1],[1,1.5],[1.5,3]];
function gen(seed,o){ o=o||{}; const R=U.rng(seed*7919+13); const u=(a,b)=>a+(b-a)*R.u(); const pick=a=>a[Math.floor(R.u()*a.length)];
  const K=o.K||pick([1,2,2,3,3,3,4,4,5,6]); const F0=u(0.9,1.6), m=u(0,0.6), step=pick([0.05,0.05,0.1]);
  const seps=[]; for(let i=1;i<K;i++){ const b=pick(SEPBIN); seps.push(u(b[0],b[1])); }
  const asym=o.asym!=null?o.asym:R.u()<0.1;
  const H=pick([300,3000,30000]); const nk=R.u(); const noise=o.noise||(nk<0.6?"poisson":nk<0.8?"ar1":"smooth"); const rho=noise==="ar1"?u(0.3,0.5):0;
  const lines=[]; let be=0; const main=Math.floor(R.u()*K);
  for(let i=0;i<K;i++){ if(i) be+=seps[i-1]*F0; lines.push({be,F:F0*u(0.9,1.1),ratio:i===main?1:Math.exp(u(Math.log(0.03),0))}); }
  const span=be, lo0=-2.5*F0-u(0.5,1.5), hi0=span+2.5*F0+u(0.5,1.5); const off=500+Math.round(u(0,100)); const lo=off+lo0, hi=off+hi0;
  const n=Math.round((hi-lo)/step)+1; const E=Float64Array.from({length:n},(_,i)=>lo+i*step); const ctx={emin:E[0],emax:E[n-1]};
  const comps=lines.map((l,i)=>{ const c=asym&&i===main?C.newComp("LA",off+l.be,{F:l.F*0.8,alpha:u(0.8,1.2),beta:u(2,4),wG:0.35*l.F}):C.newComp("GLS",off+l.be,{F:l.F,m}); c.area=1; return c; });
  // 넓이: 가장 큰 선의 높이가 H 가 되게
  const prof=comps.map(c=>C.compEval(c,E,ctx)); const hMain=Math.max(...prof[main]);
  comps.forEach((c,i)=>{ c.area=lines[i].ratio*H/hMain; });
  const P=new Float64Array(n); comps.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P), Qn=Q[n-1]||1; const Pmax=Math.max(...P);
  const base=u(0.3,2)*Pmax, step_=u(0.1,0.5)*Pmax, slope=R.u()<0.5?u(-0.3,0.3)*Pmax/(hi-lo):0;
  const mu=Float64Array.from(P,(p,i)=>base+step_*Q[i]/Qn+slope*(E[i]-lo)+p);
  let y;
  if(noise==="poisson") y=Float64Array.from(mu,v=>R.pois(v));
  else if(noise==="ar1"){ let e=0; y=Float64Array.from(mu,(v,i)=>{ const z=R.g(); e=i?rho*e+Math.sqrt(1-rho*rho)*z:z; return v+Math.sqrt(v)*e; }); }
  else { const x=Float64Array.from(mu,v=>R.pois(v)); y=Float64Array.from(x,(v,i)=>i===0||i===n-1?v:x[i-1]/6+2*v/3+x[i+1]/6); }
  const truth=comps.map(c=>({be:c.be,area:c.area,F:C.trueFWHM(c.shape,C.shapeParams(c,1)).fwhm,shape:c.shape}));
  const minSep=K>1?Math.min(...seps):null;
  return {seed,spec:{be:E,y,c:1},truth,meta:{K,F0,m,step,H,noise,rho,asym,slope:slope!==0,minSep,sepBin:minSep==null?"single":minSep<0.5?"0.3-0.5":minSep<1?"0.5-1":minSep<1.5?"1-1.5":"1.5-3",minRatio:Math.min(...lines.map(l=>l.ratio))}};
}
/* 시작 모형 = 화면 Auto suggest (suggestModel + refineByResidual, 사용자 배경 = Smart — 새 스펙트럼의 기본) */
function autoSuggest(spec){ const E=spec.be, y=spec.y; const bg={type:"smart",n:5};
  const s=C.suggestModel(E,y,{bg,shape:"GLS"}); if(!s.comps.length) return null;
  const prob=C.makeProblem(spec,{lo:E[0],hi:E[E.length-1],bg,comps:s.comps,cons:s.cons});
  const rf=C.refineByResidual(prob,{}); const P=rf.prob; return C.makeProblem(spec,{lo:E[0],hi:E[E.length-1],bg:P.bg,comps:P.comps,cons:P.cons}); }
/* 참값 대조 */
function evaluate(g,res,o){ o=o||{}; const T=g.truth, F0=g.meta.F0, gate=0.5*F0; const F=res.final; const L=F?F.lines:[];
  // 해상 불가 묶음: 간격 < 0.5 × 평균 FWHM 인 참 선들의 사슬
  const grp=[]; T.forEach((t,i)=>{ if(i&&Math.abs(t.be-T[i-1].be)<0.5*(t.F+T[i-1].F)/2) grp[grp.length-1].push(i); else grp.push([i]); });
  // 일대일 짝짓기 (가까운 것부터)
  const pairs=[]; L.forEach((l,j)=>T.forEach((t,i)=>{ const d=Math.abs(l.be-t.be); if(d<=gate) pairs.push({i,j,d}); })); pairs.sort((a,b)=>a.d-b.d);
  const mt=new Map(), mf=new Map(); pairs.forEach(p=>{ if(mt.has(p.i)||mf.has(p.j)) return; mt.set(p.i,p.j); mf.set(p.j,p.i); });
  const falseAdd=L.length-mf.size, missed=T.length-mt.size;
  const resolvable=grp.length;
  // 위치·넓이 3σ (확립된 경우의 따로 정해진 선): 묶음 안에서 맞춘 선이 하나면 묶음의 중심·합과 비교
  let posOK=true, areaOK=true, nPos=0; const zs=[];
  L.forEach((l,j)=>{ if(!mf.has(j)||!l.separate) return; const i=mf.get(j); const gi=grp.find(gg=>gg.includes(i)); const nIn=gi.filter(k=>mt.has(k)).length;
    let tb=T[i].be, ta=T[i].area; if(nIn===1&&gi.length>1){ ta=gi.reduce((s,k)=>s+T[k].area,0); tb=gi.reduce((s,k)=>s+T[k].area*T[k].be,0)/ta; }
    const blend=T.some((t,k)=>k!==i&&!mt.has(k)&&Math.abs(t.be-T[i].be)<1.5*F0&&!gi.includes(k));
    nPos++; const zb=l.beSig>0?Math.abs(l.be-tb)/l.beSig:Infinity, za=l.areaSig>0?Math.abs(l.area-ta)/l.areaSig:Infinity; zs.push({zb,za,zbU:l.beU>0?Math.abs(l.be-tb)/l.beU:Infinity,zaU:l.areaU>0?Math.abs(l.area-ta)/l.areaU:Infinity,res:gi.length===1,broad:!!l.broad,blend,fw:l.fwhm,fref:l.fwRef}); if(!(zb<=3)) posOK=false; if(!(za<=3)) areaOK=false; });
  // 해상 불가 참 묶음을 정직하게 다뤘나: 맞춘 선 ≤ 1, 또는 그 선들이 한 '따로 정해지지 않음' 묶음, 또는 '확립 안 됨'
  const unres=grp.filter(gg=>gg.length>1).map(gg=>{ const js=gg.filter(k=>mt.has(k)).map(k=>mt.get(k)); const ids=js.map(j=>L[j].id);
    const flagged=js.length<=1||(F&&F.groups.some(G_=>ids.every(id=>G_.ids.includes(id))))||(F&&!F.established); return {n:gg.length,fitted:js.length,ok:!!flagged}; });
  // 놓친 선: 간격·넓이 비별
  const miss=T.map((t,i)=>{ const nb=T.map((q,k)=>k===i?Infinity:Math.abs(q.be-t.be)/F0); const mainA=Math.max(...T.map(q=>q.area)); return {missed:!mt.has(i),sep:Math.min(...nb),ratio:t.area/mainA}; });
  return {K:T.length,nFit:L.length,resolvable,countOK:L.length===T.length,countOKres:L.length===resolvable,falseAdd,missed,established:!!(F&&F.established),reasons:F?F.reasons:[],posOK:nPos?posOK:null,areaOK:nPos?areaOK:null,nPos,zs,unres,miss}; }
function runCase(g,o){ o=o||{}; const t0=Date.now(); const P0=autoSuggest(g.spec); const tS=Date.now()-t0; if(!P0) return {seed:g.seed,meta:g.meta,err:"no start model"};
  const nStart=P0.comps.length; let R; try{ R=A.autoBuild(P0,{},{starts:o.starts}); }catch(e){ return {seed:g.seed,meta:g.meta,err:String(e.message||e)}; }
  const ev=evaluate(g,R); const evS=evaluate(g,{final:{lines:P0.comps.map(c=>({be:c.be,id:c.id,separate:false})),groups:[],established:false}});
  return {seed:g.seed,meta:g.meta,nStart,startFalseAdd:evS.falseAdd,startMissed:evS.missed,steps:R.steps.map(s=>s.kind+":"+s.label.slice(0,50)),stop:R.stop,ev,ms:Date.now()-t0,msSuggest:tS}; }
module.exports={gen,autoSuggest,evaluate,runCase,SEPBIN};

/* 레시피 시작 경우 (SPEC §18.4): C 1s (C–C/C–H + C–O ± C=O ± O–C=O, GL(30), 같은 폭) · O 1s (격자 ± 낮은 격자 528.6 + 표면 O) — 참 위치는 레시피 창 안.
   시작 = 레시피 모형 + 레시피 피팅(c1sFit: 단계 피팅 · 약한 선 F 검정 / o1sFit), 그다음 Auto build. 가짜 선 = Auto build 가 더한 선 중 참에 없는 것 */
const Gd=require("../src/xps_guide.js");
function genRecipe(seed){ const R=U.rng(seed*104729+7); const u=(a,b)=>a+(b-a)*R.u(); const pick=a=>a[Math.floor(R.u()*a.length)];
  const kind=seed%2?"c1s":"o1s"; const H=pick([3000,30000]); const noise=pick(["poisson","poisson","smooth"]); const step=0.05; const lines=[];
  if(kind==="c1s"){ const cc=284.8+u(-0.4,0.4), F=u(1.0,1.5); lines.push([cc,1,F]); lines.push([cc+u(1.3,1.7),u(0.1,0.4),F]); if(R.u()<0.5) lines.push([cc+u(2.9,3.3),u(0.03,0.1),F]); if(R.u()<0.7) lines.push([cc+u(3.9,4.3),u(0.05,0.2),F]); }
  else { const lb=u(529.0,530.0), F=u(1.0,1.4), Fs=u(1.6,2.0); lines.push([lb,1,F]); if(R.u()<0.4) lines.push([lb-u(0.8,1.2),u(0.1,0.3),F]); lines.push([u(531.3,531.7),u(0.3,0.8),Fs]); lines.push([u(531.95,532.25),u(0.1,0.2),Fs]); if(R.u()<0.5) lines.push([u(532.5,532.9),u(0.1,0.2),Fs]); if(R.u()<0.5) lines.push([u(533.3,533.7),u(0.05,0.15),Fs]); }
  lines.sort((a,b)=>a[0]-b[0]);
  const lo=kind==="c1s"?lines[0][0]-6:524.5, hi=kind==="c1s"?lines[0][0]+10:537.5; const n=Math.round((hi-lo)/step)+1; const E=Float64Array.from({length:n},(_,i)=>lo+i*step), ctx={emin:E[0],emax:E[n-1]};
  const comps=lines.map(([b,r,F])=>{ const c=C.casaToComp("GL(30)",{be:b,fwhm:F,area:1}).comp; return c; }); const prof=comps.map(c=>C.compEval(c,E,ctx)); const im=lines.findIndex(l=>l[1]===1); const hM=Math.max(...prof[im]);
  comps.forEach((c,i)=>{ c.area=lines[i][1]*H/hM; }); const P=new Float64Array(n); comps.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P), Qn=Q[n-1]; const Pm=Math.max(...P); const base=u(0.3,1)*Pm, st_=u(0.1,0.4)*Pm, sl=R.u()<0.5?u(-0.2,0.2)*Pm/(hi-lo):0;
  const mu=Float64Array.from(P,(p,i)=>base+st_*Q[i]/Qn+sl*(E[i]-lo)+p); let y;
  if(noise==="poisson") y=Float64Array.from(mu,v=>R.pois(v)); else { const x=Float64Array.from(mu,v=>R.pois(v)); y=Float64Array.from(x,(v,i)=>i===0||i===n-1?v:x[i-1]/6+2*v/3+x[i+1]/6); }
  const truth=comps.map(c=>({be:c.be,area:c.area,F:C.trueFWHM(c.shape,C.shapeParams(c,1)).fwhm,shape:c.shape}));
  const F0=truth.reduce((s,t)=>s+t.F,0)/truth.length;
  return {seed,spec:{be:E,y,c:1},truth,meta:{recipe:kind,K:lines.length,F0,H,noise,step,slope:sl!==0,latlo:kind==="o1s"&&lines.length>1&&lines[0][0]<lines[1][0]-0.5&&lines[1][0]<531,asym:false,sepBin:"recipe-"+kind}}; }
function recipeStart(g){ const sp=g.spec; const R=Gd.recipeById(g.meta.recipe);
  const M=Gd.recipeModel(g.meta.recipe,{shift:0,be:sp.be,y:sp.y,present:[]}); const P=C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:M.bg,comps:M.comps,cons:M.cons});
  const s=g.meta.recipe==="c1s"?Gd.c1sFit(P,{}):Gd.o1sFit(P,{build:{shift:0,present:[]}}); return C.makeProblem(sp,{lo:M.lo,hi:M.hi,bg:s.P.bg,comps:s.P.comps,cons:s.P.cons}); }
function runRecipeCase(g,o){ const t0=Date.now(); let P0; try{ P0=recipeStart(g); }catch(e){ return {seed:g.seed,meta:g.meta,err:"recipe: "+e.message}; }
  let R; try{ R=A.autoBuild(P0,{el:g.meta.recipe==="c1s"?"C":"O",orb:"1s",present:[]},{first:false,starts:o&&o.starts}); }catch(e){ return {seed:g.seed,meta:g.meta,err:String(e.message||e)}; }
  const ev=evaluate(g,R); const added=R.final.lines.filter(l=>l.added); const T=g.truth, gate=0.5*g.meta.F0;
  ev.addedLines=added.length; ev.falseAdd=added.filter(l=>!T.some(t=>Math.abs(t.be-l.be)<=gate)).length;   // 레시피 경우: Auto build 가 더한 선만 셈
  ev.recipeExtra=R.final.lines.filter(l=>l.recipe&&!T.some(t=>Math.abs(t.be-l.be)<=gate)).length;
  ev.latloFound=g.meta.latlo?R.final.lines.some(l=>Math.abs(l.be-T[0].be)<=gate):null;
  return {seed:g.seed,meta:g.meta,nStart:P0.comps.length,startFalseAdd:0,steps:R.steps.map(s=>s.kind+":"+s.label.slice(0,50)),stop:R.stop,ev,ms:Date.now()-t0}; }
module.exports.genRecipe=genRecipe; module.exports.runRecipeCase=runRecipeCase;
if(require.main===module){ const s0=+process.argv[2]||1, cnt=+process.argv[3]||10, out=process.argv[4]; const fs=require("fs"); if(out) fs.writeFileSync(out,"");
  const rec=process.argv[5]==="recipe";
  for(let s=s0;s<s0+cnt;s++){ const g=rec?genRecipe(s):gen(s); const r=rec?runRecipeCase(g):runCase(g); const line=JSON.stringify(r); if(out) fs.appendFileSync(out,line+"\n");
    console.log(`${s} K${g.meta.K} ${g.meta.sepBin} H${g.meta.H} ${g.meta.noise}${g.meta.asym?" asym":""} start ${r.nStart} → ${r.ev?r.ev.nFit:"-"} fa ${r.ev?r.ev.falseAdd:"-"} miss ${r.ev?r.ev.missed:"-"} est ${r.ev?r.ev.established:"-"} pos ${r.ev?r.ev.posOK:"-"} ${r.ms} ms ${r.err||""}`); } }
