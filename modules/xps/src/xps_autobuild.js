/* INNO Analytics — XPS Auto build (v1.10, SPEC §18). 피팅 모형 만들기 고리(잔차 → 선 더하기/빼기 → 제약 → 비교 → 불확도)를 자동으로 *제안*한다.
   원칙 (SPEC §18.1): 제안만 하고 바꾸지 않음(화면이 받아들일 때만 적용, 한 걸음 = 되돌리기 하나) · 원소·화합물과 무관한 규칙만 ·
   화학 이름을 붙이지 않음(새 선은 BE 로; 후보 화학종은 앱 안의 출처 있는 표에서만, '후보'로) · 정직한 불확도(|r| > 0.9 묶음은 합만,
   같은 최솟값 < 1/3 또는 χ²ν > 2 이면 '모형 확립 안 됨') · 레시피 존중(레시피 선의 문헌 범위·제약은 풀지 않음).
   모든 비교는 데이터 분산 χ²_d(§3.5)와 AR(1) 유효 점 수 N_eff = N(1 − ρ₁)/(1 + ρ₁) (Zięba 2010 식 15)로 계산한 AICc·BIC.
   이 파일은 xps_core.js 의 피팅을 바꾸지 않는다 — Auto build 를 쓰지 않으면 모든 결과가 그대로(비트 동일). */
(function(root){
"use strict";
function getCore(){ if(root&&root.XPS_CORE) return root.XPS_CORE; if(typeof require==="function"){ try{ return require("./xps_core.js"); }catch(e){} } return null; }
function getGuide(){ if(root&&root.XPS_GUIDE) return root.XPS_GUIDE; if(typeof require==="function"){ try{ return require("./xps_guide.js"); }catch(e){} } return null; }
let _C=null, _G=null; const CORE=()=>_C||(_C=getCore()), GD=()=>_G||(_G=getGuide());

/* 문턱 (SPEC §18.2 — 벤치마크 전에 정함; 바꾸면 §18.5 를 다시 잼) */
const AB_TH={
  T:10,          // 자유도를 늘리는 걸음(선 더하기 · 한계 풀기 · 배경 기울기 · 정적 → 피팅 안 배경): ΔAICc_eff ≥ 10 그리고 ΔBIC_eff ≥ 10
                 //   (Burnham & Anderson 2002: ΔAIC > 10 이면 나쁜 쪽 지지 거의 없음; Kass & Raftery 1995: ΔBIC > 10 매우 강함)
  simp:2,        // 단순화(공유 L 비율 · 폭 묶기 · 기울기 끄기 · 레시피 선 빼기): 새 모형의 AICc_eff 가 2 넘게 나빠지지 않으면 — |ΔAICc| < 2 이면 매개변수가 적은 쪽
  zRun:4,        // 잔차 덩어리 문턱(데이터 2차 차분의 잡음·이웃 상관으로 정규화, 모형과 무관 — SPEC §17.6)
  runHalf:0.5,   // 덩어리 반폭 = 0.5 × 중앙 FWHM
  minSep:0.5,    // 두 선의 간격 ≥ 0.5 × 평균 FWHM (R4 와 같음; 그보다 가까운 두 Gaussian 의 합은 폭이 조금 넓은 한 선과 4 차 항까지만 다름)
  floorRel:0.7,  // 해상도 바닥 = 0.7 × 가장 좁은 '잘 정해진' 선의 FWHM (σ_F/F < 10 %, 한계 아님), 그리고 ≥ 3 채널
  floorSig:0.10,
  fhiRel:1.5,    // 새 선 FWHM 상한 = 1.5 × 가장 넓은 선
  maxLinesMax:10,// 선 수 상한 = min(10, ⌊창 폭 / 중앙 FWHM⌋)
  maxSteps:12,
  corrGroup:0.9, // |r| > 0.9 인 두 선 → '따로 정해지지 않음', 합만 보고
  broadRel:1.25, // 넓은 선 표지: FWHM > 1.25 × 가장 좁은 잘 정해진 선 (그리고 3σ 넘게) — 참 선끼리 폭 차이(벤치마크 ±10 %)보다 크게
  slopeOff:false,// '기울기 끄기' 단순화 — 꺼짐(§18.5): 검증 묶음 1 에서 가짜 선 8 경우 중 5 에 관여(창 끝 선과 기울기가 서로 바꿔 쓰여, 선이 있으면 기울기가
                 // |ΔAICc| < 2 로 빠지고 그 뒤로 그 선은 뺄 수 없음 — 예 seed 100009 의 609.1 eV 선, 100324 의 FWHM 4 eV(한계) 선)
  chi2nuMax:2, sameMin:1/3, starts:12, rhoMax:0.95, boundRel:1e-3
};
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const med=a=>{ const s=a.slice().sort((x,y)=>x-y); const n=s.length; return n?(n%2?s[n>>1]:(s[n/2-1]+s[n/2])/2):NaN; };
const fx=(v,d)=>v==null||!isFinite(v)?"—":(+v).toFixed(d);
const sg=(v,d)=>v==null||!isFinite(v)?"—":(v>0?"+":"")+(+v).toFixed(d);

/* ---------- 피팅·통계 ---------- */
function clonePs(P){ const C=CORE(); const Q=C.clone(C.stripProb(P)); Q.E=P.E; Q.y=P.y; return Q; }
function fitS(P,o){ const C=CORE(); o=o||{}; try{ const r=C.fitOnce(P,{initAreas:o.init!==false}); const st=C.fitStats(P,r); if(!isFinite(st.chi2d)) return null; return {P,r,st}; }catch(e){ return null; } }
const isStatic=bg=>!!bg&&bg.type!=="active"&&bg.type!=="none";
/* 정적 배경의 두 끝값은 자료에서 정한 수(피팅 안 Shirley 의 c0·k 에 해당) → 매개변수 2 개로 셈 (SPEC §18.2) */
const pEff=S=>S.st.p+(isStatic(S.P.bg)?2:0);
function icEff(S,rho){ const N=S.st.N, Ne=N*(1-rho)/(1+rho), p=pEff(S), L=Ne*Math.log(S.st.chi2d/N);
  return {Neff:Ne,p,AICc:L+2*p+(Ne-p-1>0?2*p*(p+1)/(Ne-p-1):Infinity),BIC:L+p*Math.log(Ne)}; }
/* 두 모형 비교: ρ = 두 잔차 ρ₁ 중 작은 쪽(잡음에 가까운 쪽 — 빠진 선의 구조가 ρ₁ 를 키움), [0, 0.95]. 부호: 양수 = b(후보)가 좋음 */
/* 비교에 쓰는 ρ (SPEC §18.2): 두 모형 중 χ²_d 가 낮은(더 잘 맞는) 쪽의 잔차 ρ₁, [0, 0.5]. 나쁜 쪽 ρ₁ 은 빠진 선의 구조가 키우므로(0.9 넘음) 쓰지 않는다.
   상한 0.5: 실측 VGD 잡음 0.31–0.47(§17.6), 평활 [1/6, 2/3, 1/6] 0.44 — 잡음이 더 상관이면 N_eff 를 크게 보므로, 끝 모형의 ρ₁ > 0.5 는
   '확립 안 됨' 이유로 보고(§18.3). (자료만의 추정 — 2차 차분 MA(1) 은 큰 세기에서 봉우리 곡률이 새고, 4차 차분은 ρ 에 거의 둔감(−0.80 → −0.76) — 은 쓰지 않음) */
const RHO_CAP=0.5;
function rhoNoise(P){ const C=CORE(), G=GD(); const ys=C.smoothQ(P.E,P.y,C.smoothHW(P.E)); const vv=Float64Array.from(ys,x=>Math.max(Math.abs(x),1e-9)); return G.noiseMA1(P.y,vv).rho; }
function cmpIC(a,b){ const bt=b.st.chi2d<a.st.chi2d?b:a; const rho=clamp(bt.st.rho1,0,RHO_CAP); const A=icEff(a,rho), B=icEff(b,rho);
  return {rho,Neff:A.Neff,dAICc:A.AICc-B.AICc,dBIC:A.BIC-B.BIC,dChi2nu:b.st.chi2nu-a.st.chi2nu,dp:B.p-A.p,chi2nu:b.st.chi2nu}; }

/* ---------- 선 정보 · 한계 · 상한 ---------- */
function fwhmOf(c){ const C=CORE(); try{ return C.trueFWHM(c.shape,C.shapeParams(c,1)).fwhm||c.F||1; }catch(e){ return c.F||1; } }
function heightPerArea(c){ const C=CORE(); try{ return C.trueFWHM(c.shape,C.shapeParams(c,1)).ymax||1; }catch(e){ return 1; } }
const consOf=(P,id,k)=>((P.cons||{})[id]||{})[k]||null;
const ownFree=(P,c,k)=>{ const q=consOf(P,c.id,k); return !q||q.mode==="free"; };
function withCov(S){ if(!S.cov){ try{ S.cov=CORE().covariance(S.P,S.r); }catch(e){ S.cov=null; } } return S.cov; }
/* 한계에 닿은 비선형 매개변수 (R1 과 같은 기준) */
function boundHits(S){ const C=CORE(); const map=S.r.map||C.buildMap(S.P), v=C.getVec(S.P,map), out=[];
  map.free.forEach((f,j)=>{ if(f.k==="area") return; const L=isFinite(f.lo), H=isFinite(f.hi); if(!L||!H) return; const rg=f.hi-f.lo, p=v[j];
    if(Math.min(p-f.lo,f.hi-p)<=AB_TH.boundRel*rg) out.push({id:f.id,k:f.k,ci:f.ci,side:p-f.lo<=f.hi-p?"lower":"upper",v:p,lo:f.lo,hi:f.hi}); }); return out; }
function caps(S){ const C=CORE(); const P=S.P, E=P.E, n=E.length, step=(E[n-1]-E[0])/(n-1); const fw=P.comps.map(fwhmOf); const cov=withCov(S);
  const hits=new Set(boundHits(S).map(h=>h.id+"|"+h.k)); const good=[];
  if(cov) cov.names.forEach((q,j)=>{ if(q.k!=="F"||q.ci<0) return; const c=P.comps[q.ci]; const s=cov.sig[j]; if(s>0&&isFinite(s)&&s/c.F<AB_TH.floorSig&&!hits.has(q.id+"|F")) good.push(fw[q.ci]); });
  const base=good.length?Math.min(...good):(fw.length?Math.min(...fw):1);
  const floor=Math.max(3*step,AB_TH.floorRel*base), Fmed=fw.length?med(fw):1, W=E[n-1]-E[0];
  return {fw,floor,Fmed,Fmax:fw.length?Math.max(...fw):Fmed,nMax:Math.max(P.comps.length,Math.min(AB_TH.maxLinesMax,Math.floor(W/Math.max(Fmed,1e-6)))),step,W,floorFrom:good.length?"narrowest well-determined line":"narrowest line"}; }
const sepOK=(be1,F1,be2,F2)=>Math.abs(be1-be2)>=AB_TH.minSep*(F1+F2)/2;

/* ---------- 잔차 덩어리 (모형과 무관한 잡음 척도, SPEC §17.6) ---------- */
function resRuns(S,K){ const G=GD(); const P=S.P, E=P.E, y=P.y, n=E.length, mu=S.r.mu, v=S.r.v; const cp=caps(S);
  const h=Math.max(2,Math.round(AB_TH.runHalf*cp.Fmed/cp.step)); const nz=G.noiseMA1(y,v); const r=new Float64Array(n);
  for(let i=0;i<n;i++) r[i]=(y[i]-mu[i])/Math.sqrt(v[i]);
  const z=new Float64Array(n).fill(0); for(let i=h;i<n-h;i++){ let s=0; for(let j=i-h;j<=i+h;j++) s+=r[j]; z[i]=s/Math.sqrt(nz.lrv*(2*h+1)); }
  const pk=[]; for(let i=h;i<n-h;i++){ if(z[i]>0&&z[i]>=z[i-1]&&z[i]>=z[i+1]) pk.push({i,be:E[i],z:z[i],w:(2*h+1)*cp.step,amp:y[i]-mu[i]}); }
  pk.sort((a,b)=>b.z-a.z); const out=[]; pk.forEach(q=>{ if(out.length<(K||2)&&out.every(o=>Math.abs(o.be-q.be)>=cp.Fmed)) out.push(q); });
  let neg=null; for(let i=h;i<n-h;i++) if(!neg||z[i]<neg.z) neg={be:E[i],z:z[i]};
  return {runs:out,neg,rhoNoise:nz.rho,half:h*cp.step}; }

/* ---------- 모형 조작 (복사본에서만) ---------- */
function removeIds(P,ids){ const C=CORE(), G=GD(); const set=new Set(ids); if(P.comps.some(c=>set.has(c.id)&&c.recipe)) G.relinkWithout(P,set); else ids.forEach(id=>C.removeComp(P,id)); }
const nameOf=(c,sh)=>c.name||((c.be+(sh||0)).toFixed(1)+" eV line");

/* ---------- 후보 걸음 ---------- */
function candidates(S,ctx,excl){
  const C=CORE(), G=GD(); const P=S.P, E=P.E, n=E.length, sh=ctx.shift||0; const cp=caps(S); const cov=withCov(S); const out=[];
  const tot=P.comps.reduce((s,c)=>s+Math.max(0,c.area),0)||1; const nm=c=>nameOf(c,sh);
  const push=(m)=>{ if(excl&&excl.has(m.key)) return; out.push(m); };
  const run=(kind,key,label,cls,build,o)=>{ if(excl&&excl.has(key)) return null; const Q=clonePs(P); let info=null; try{ info=build(Q); }catch(e){ return null; } if(info===false) return null;
    const T=fitS(Q,{init:!(o&&o.noInit)}); if(!T) return null; const d=cmpIC(S,T); const m={kind,key,label,cls,S2:T,d,info:info||{}}; push(m); return m; };
  const hits=boundHits(S);
  // (1) 정리: 넓이 0 또는 중심이 창 밖
  P.comps.forEach((c,i)=>{ if(P.comps.length<2) return; const dead=c.area<=1e-6*tot||hits.some(h=>h.id===c.id&&h.k==="area")||c.be<E[0]||c.be>E[n-1]||(cov&&cov.names.some((q,j)=>q.id===c.id&&q.k==="area"&&cov.atB[j]));
    if(dead) run("dead","dead|"+c.id,`Remove ${nm(c)} (area at 0 or centre outside the window)`,"cleanup",Q=>{ removeIds(Q,[c.id]); return {ids:[c.id]}; }); });
  // (2) 상한: 너무 가까운 두 선 (자유 위치끼리, 같은 레시피의 선끼리는 문헌 범위가 정함)
  for(let i=0;i<P.comps.length;i++) for(let j=i+1;j<P.comps.length;j++){ const a=P.comps[i], b=P.comps[j];
    if(!ownFree(P,a,"be")||!ownFree(P,b,"be")) continue; if(a.recipe&&a.recipe===b.recipe) continue; if(sepOK(a.be,cp.fw[i],b.be,cp.fw[j])) continue;
    const sm=a.area<=b.area?a:b; const key="merge|"+[a.id,b.id].sort().join("|");
    const m=run("merge",key,`Lines ${nm(a)} and ${nm(b)} merge (ΔBE ${fx(Math.abs(a.be-b.be),2)} eV < ${AB_TH.minSep} × FWHM) — remove ${nm(sm)}`,"cap",Q=>{ removeIds(Q,[sm.id]); return {ids:[sm.id]}; });
    if(m&&!(-m.d.dAICc>=AB_TH.T&&-m.d.dBIC>=AB_TH.T)) continue;   // 빼도 나빠지지 않음 → 빼기
    if(m) out.splice(out.indexOf(m),1);
    if(excl&&excl.has(key+"|b")) continue;   // 한 쌍에 'Bound positions' 는 한 번만(검증 묶음 1 seed 100094: 같은 걸음 9 번) — 그래도 가까우면 끝 판정이 '확립 안 됨'
    run("bound",key+"|b",`Lines ${nm(a)} and ${nm(b)} merge but both are needed — bound positions ≥ ${AB_TH.minSep} × FWHM apart`,"cap",Q=>{ const msg=G.applyAction(Q,{type:"bound",ids:[a.id,b.id],sep:+(AB_TH.minSep*(cp.fw[i]+cp.fw[j])/2).toFixed(3)}); if(msg==null) return false; return {ids:[a.id,b.id]}; },{noInit:true}); }
  // (3) 단순화
  // 3a 선 빼기: 자동/사용자 선은 '더하기와 같은 기준'(빼면 ΔAICc_eff 와 ΔBIC_eff 가 모두 ≥ 10 나빠질 때만 남김), 레시피 선은 |ΔAICc| < 2 일 때만 뺌
  if(P.comps.length>1) P.comps.forEach(c=>{ run("drop","drop|"+c.id,`Remove ${nm(c)}`,"simp",Q=>{ removeIds(Q,[c.id]); return {ids:[c.id],recipe:!!c.recipe}; }); });
  // 3b 공유 L 비율 (자유 m 이 둘 이상인 GL 선)
  { const gl=P.comps.filter(c=>(c.shape==="GLS"||c.shape==="GLP")&&c.m!=null&&!C.derivedParam(c,"m")&&!c.recipe); const fr=gl.filter(c=>ownFree(P,c,"m"));
    const anc=fr.reduce((a,c)=>!a||c.area>a.area?c:a,null);
    if(fr.length>=2) run("shareM","shareM","One shared L fraction for the GL lines","simp",Q=>{ fr.forEach(c=>{ if(c===anc) return; const q=Q.cons[c.id]=Q.cons[c.id]||{}; q.m={mode:"link",ref:anc.id,k:1,b:0}; }); Q.comps.forEach(c=>{ if(fr.some(x=>x.id===c.id)) c.m=anc.m; }); return {ids:fr.map(c=>c.id)}; }); }
  // 3c 폭 묶기 (같은 폭: k = 1 — 지금 비로 묶으면 추정값을 상수로 바꿔 매개변수 수를 속임). 전체 하나 / 통계적으로 다르지 않은 무리
  { const fl=P.comps.map((c,i)=>({c,i})).filter(o=>(o.c.shape==="GLS"||o.c.shape==="GLP")&&!o.c.recipe&&!C.derivedParam(o.c,"F")&&ownFree(P,o.c,"F"));
    const sF=o=>{ if(!cov) return null; const j=cov.names.findIndex(q=>q.id===o.c.id&&q.k==="F"); return j>=0&&cov.sig[j]>0?cov.sig[j]:null; };
    const tie=(grp,Q)=>{ const anc=grp.reduce((a,o)=>!a||o.c.area>a.c.area?o:a,null); const lo=Math.min(...grp.map(o=>(consOf(P,o.c.id,"F")||{}).lo??0.3)), hi=Math.max(...grp.map(o=>(consOf(P,o.c.id,"F")||{}).hi??4));
      const qa=Q.cons[anc.c.id]=Q.cons[anc.c.id]||{}; qa.F={mode:"free",lo,hi}; grp.forEach(o=>{ if(o===anc) return; const q=Q.cons[o.c.id]=Q.cons[o.c.id]||{}; q.F={mode:"link",ref:anc.c.id,k:1,b:0}; Q.comps[o.i].F=anc.c.F; }); };
    if(fl.length>=2){ run("tieF","tieF|all","Tie all FWHMs (equal widths)","simp",Q=>{ tie(fl,Q); return {ids:fl.map(o=>o.c.id)}; });
      const srt=fl.slice().sort((a,b)=>a.c.F-b.c.F); const cl=[[srt[0]]]; for(let t=1;t<srt.length;t++){ const g=cl[cl.length-1], p=g[g.length-1], q=srt[t]; const sa=sF(p), sb=sF(q); if(sa!=null&&sb!=null&&Math.abs(q.c.F-p.c.F)<2*Math.hypot(sa,sb)) g.push(q); else cl.push([q]); }
      const multi=cl.filter(g=>g.length>=2); if(multi.length&&cl.length>1) run("tieF","tieF|cl",`Tie FWHMs within ${multi.length} group(s) of similar width`,"simp",Q=>{ multi.forEach(g=>tie(g,Q)); return {groups:multi.map(g=>g.map(o=>o.c.id))}; }); } }
  // 3d 배경 기울기 끄기
  if(AB_TH.slopeOff&&P.bg.type==="active"&&P.bg.cons&&P.bg.cons.s&&P.bg.cons.s.mode!=="fixed") run("bgOff","bg|slopeOff","Background: in-fit Shirley without slope","simp",Q=>{ Q.bg.cons=Object.assign({},Q.bg.cons,{s:{mode:"fixed"}}); Q.bg.s=0; return {}; },{noInit:true});
  // (4) 자유도를 늘리는 걸음
  // 4a 선 더하기: 가장 큰 양의 잔차 덩어리 (상위 2 개, 중앙 FWHM 넘게 떨어진)
  const rr=resRuns(S,2);
  if(P.comps.length<cp.nMax) rr.runs.forEach(q=>{ if(q.z<AB_TH.zRun) return; const be=q.be; if(P.comps.some((c,i)=>!sepOK(be,cp.Fmed,c.be,cp.fw[i]))) return;
    const key="add|"+(be+sh).toFixed(1);
    const m=run("add",key,`Add a line at ${fx(be+sh,2)} eV (residual run z ${fx(q.z,1)})`,"add",Q=>{
      const anc=Q.comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); const shp=anc&&(anc.shape==="GLS"||anc.shape==="GLP")?anc.shape:"GLS";
      const F0=clamp(cp.Fmed,cp.floor,Math.max(cp.floor*1.5,AB_TH.fhiRel*cp.Fmax)); const c=C.newComp(shp,be,{F:F0,el:anc?anc.el:"",line:anc?anc.line:""});
      c.name=(be+sh).toFixed(1)+" eV line"; c.ab=1; const cs=Q.cons[c.id]={};
      if(anc&&anc.m!=null&&c.m!=null){ c.m=anc.m; const qm=consOf(Q,anc.id,"m"); if(qm&&qm.mode==="fixed") cs.m={mode:"fixed"}; else cs.m={mode:"link",ref:qm&&qm.mode==="link"?qm.ref:anc.id,k:1,b:0}; }
      cs.F={mode:"free",lo:+cp.floor.toFixed(4),hi:+Math.max(cp.floor*1.5,AB_TH.fhiRel*cp.Fmax).toFixed(4)};
      const wb=Math.max(1,cp.Fmed); cs.be={mode:"free",lo:+Math.max(E[0],be-wb).toFixed(4),hi:+Math.min(E[n-1],be+wb).toFixed(4)};
      const d0=Q.comps.filter(x=>x.doublet).reduce((a,x)=>!a||x.area>a.area?x:a,null);
      if(d0&&be+d0.doublet.dE<E[n-1]-0.5){ c.doublet={dE:d0.doublet.dE,r:d0.doublet.r,q:d0.doublet.q}; ["dE","r","q"].forEach(k=>{ const dc=consOf(Q,d0.id,k); cs[k]=!dc||dc.mode==="fixed"?{mode:"fixed"}:dc.mode==="link"?{mode:"link",ref:dc.ref,k:1,b:0}:{mode:"link",ref:d0.id,k:1,b:0}; }); }
      c.area=Math.max(1e-9,q.amp)/heightPerArea(c)*(c.doublet?1+c.doublet.r:1); Q.comps.push(c); return {be,run:{be:be+sh,z:q.z,w:q.w},id:c.id}; });
    if(!m) return; const T=m.S2, nc=T.P.comps.find(c=>c.id===m.info.id); const hb=boundHits(T); const tot2=T.P.comps.reduce((s,c)=>s+Math.max(0,c.area),0)||1;
    const fw2=T.P.comps.map(fwhmOf), k=T.P.comps.indexOf(nc);
    // 새 선의 BE·넓이·FWHM 위 한계에 닿으면 그 선은 자료가 아니라 한계가 정한 것(끝 판정과 같은 기준; FWHM 은 검증 묶음 2 seed 300159 뒤 추가 — 창 끝 배경을 먹는 넓은 선)
    // (v1.10 확인 4: 아래 한계(해상도 바닥)는 빼지 않음 — 세션형 O 1s 의 참 528.4 eV 선이 레시피 C=O 를 뺀 뒤 바닥에서 시작해 막혔음; 그 경우는 끝 판정이 '한계' 로 막음)
    const bad=!nc||nc.area<=1e-6*tot2||hb.some(h=>h.id===nc.id&&(h.k==="be"||h.k==="area"||(h.k==="F"&&h.side==="upper")))||T.P.comps.some((c,i)=>i!==k&&!sepOK(nc.be,fw2[k],c.be,fw2[i]))||T.P.comps.some(c=>c.area<=1e-6*tot2);
    if(bad){ m.capFail="the added line merged with another line, reached a range end (BE or FWHM) or killed a line"; } });
  // 4b 한계 풀기 (레시피가 아닌 범위만: 레시피의 문헌 범위는 물리 근거)
  hits.forEach(h=>{ if(h.ci<0) return; const c=P.comps[h.ci]; if(!c||c.recipe) return;
    if(h.k==="F"){ const lo=h.side==="lower"?Math.max(cp.floor,h.lo/1.5):h.lo, hi=h.side==="upper"?Math.min(cp.W/2,h.hi*1.5):h.hi; if(hi-h.hi<1e-6&&h.lo-lo<1e-6) return;
      run("release","release|"+c.id+"|F",`Release ${nm(c)} FWHM ${h.side} bound (${fx(h.side==="upper"?h.hi:h.lo,2)} → ${fx(h.side==="upper"?hi:lo,2)} eV)`,"add",Q=>{ const q=Q.cons[c.id]=Q.cons[c.id]||{}; q.F={mode:"free",lo:+lo.toFixed(4),hi:+hi.toFixed(4)}; return {id:c.id,k:"F"}; },{noInit:true}); }
    else if(h.k==="be"){ const bv=h.side==="lower"?h.lo:h.hi; if(Math.abs(bv-E[0])<=1.5*cp.step||Math.abs(bv-E[n-1])<=1.5*cp.step) return;
      const lo=Math.max(E[0],h.lo-0.3), hi=Math.min(E[n-1],h.hi+0.3);
      run("release","release|"+c.id+"|be",`Widen ${nm(c)} BE range ±0.3 eV`,"add",Q=>{ const q=Q.cons[c.id]=Q.cons[c.id]||{}; q.be={mode:"free",lo:+lo.toFixed(4),hi:+hi.toFixed(4)}; return {id:c.id,k:"be"}; },{noInit:true}); } });
  // 4c 배경: 기울기 켜기 / 정적 → 피팅 안 Shirley (± 기울기)
  if(P.bg.type==="active"&&!(P.bg.cons&&P.bg.cons.s&&P.bg.cons.s.mode!=="fixed")) run("bgOn","bg|slopeOn","Background: in-fit Shirley + slope","add",Q=>{ Q.bg.cons=Object.assign({},Q.bg.cons,{s:{mode:"free"}}); return {}; },{noInit:true});
  // 정적 배경 → 피팅 안 Shirley + 기울기: IC 문턱이 아니라 §9 의 측정(정적 배경으로 판단하면 참 1 성분 143 개 중 106 개에 가짜 성분)과 이 벤치마크 개발 중 측정
  // (정적 Smart 를 둔 단일 선: 위치 4–9σ · 넓이 4–8σ 어긋남 — 정적 끝값의 오차는 σ 에 들어가지 않음)으로 정한 '판단 배경' — 첫 제안. 버리면 요약에 경고
  if(isStatic(P.bg)&&ctx.staticOK!==true) run("bgAct","bg|activeS","Background: in-fit Shirley + slope (instead of static "+P.bg.type+") — decision background","start",Q=>{ Q.bg={type:"active",n:P.bg.n||5,c0:null,k:null,s:0,cons:{s:{mode:"free"}}}; return {}; });
  return {list:out,runs:rr,caps:cp,hits}; }

/* 받아들일 수 있는가 (SPEC §18.2) */
function passes(m){ const d=m.d, T=AB_TH.T;
  if(m.capFail) return false; if(!isFinite(d.dAICc)||!isFinite(d.dBIC)) return m.cls==="cleanup"||m.cls==="start";
  if(m.cls==="cleanup"||m.cls==="cap"||m.cls==="start") return true;
  if(m.kind==="drop") return m.info.recipe?d.dAICc>-AB_TH.simp:!(-d.dAICc>=T&&-d.dBIC>=T);
  if(m.cls==="simp") return d.dAICc>-AB_TH.simp;
  return d.dAICc>=T&&d.dBIC>=T; }

/* 한 걸음: 정리 → 상한 → 단순화(가장 좋은 AICc) → 더하기(가장 큰 ΔAICc) */
function stepFrom(S,ctx,excl){
  const cd=candidates(S,ctx,excl); const L=cd.list; const ok=L.filter(passes);
  const pick=cls=>{ const a=ok.filter(m=>m.cls===cls); if(!a.length) return null; return a.reduce((b,m)=>!b||m.d.dAICc>b.d.dAICc?m:b,null); };
  // 더하기 무리 안에서는 배경 걸음이 먼저(배경 모양 오차를 선이 메우지 않게 — 개발 중 seed 1006: 정적 Smart + 기울기 참값에서 창 끝 585.7 eV 에 가짜 선)
  const bgFirst=()=>{ const a=ok.filter(m=>(m.cls==="add"||m.cls==="simp")&&/^bg/.test(m.kind)); return a.length?a.reduce((b,m)=>!b||m.d.dAICc>b.d.dAICc?m:b,null):null; };
  const best=pick("start")||pick("cleanup")||pick("cap")||bgFirst()||pick("simp")||pick("add");
  const rej=L.filter(m=>!passes(m)).map(m=>({kind:m.kind,key:m.key,label:m.label,dAICc:m.d.dAICc,dBIC:m.d.dBIC,cap:m.capFail||null}));
  return {best,cands:L.length,rejected:rej,runs:cd.runs,caps:cd.caps,hits:cd.hits}; }

/* 걸음 기록(직렬화) */
function stepRec(S,m,ctx,meta){ const sh=ctx.shift||0;
  m.S2.P.comps.forEach(c=>{ if(c.ab&&/^-?\d+\.\d eV line$/.test(c.name||"")) c.name=(c.be+sh).toFixed(1)+" eV line"; });   // 더한 선의 이름 = 맞춘 BE (처음 놓은 자리 아님)
  const hb=boundHits(m.S2).map(h=>`${nameOf(m.S2.P.comps[h.ci]||{},sh)} ${h.k} at ${h.side} bound`);
  return {kind:m.kind,key:m.key,cls:m.cls,label:m.label,
    ev:{chi2nuFrom:S.st.chi2nu,chi2nuTo:m.S2.st.chi2nu,dChi2nu:m.d.dChi2nu,dAICc:m.d.dAICc,dBIC:m.d.dBIC,dp:m.d.dp,rho:m.d.rho,Neff:m.d.Neff,N:S.st.N,
      run:m.info&&m.info.run?m.info.run:(meta.runs&&meta.runs.runs[0]?{be:meta.runs.runs[0].be+sh,z:meta.runs.runs[0].z,w:meta.runs.runs[0].w}:null),bounds:hb,sameMin:m.sameMin||null},
    why:whyText(m),model:{comps:m.S2.P.comps,cons:m.S2.P.cons,bg:m.S2.P.bg},newId:m.info&&m.info.id||null}; }
function whyText(m){ const d=m.d, T=AB_TH.T;
  if(m.cls==="start") return "Decisions are made with an in-fit background: a static background's end levels are fixed from the data, their error is not in any σ, and residual structure from a wrong static shape is taken for missing lines (SPEC §9: 106 of 143 one-line spectra got a false line). The slope is removed again if the data do not need it.";
  if(m.cls==="cleanup") return "A line with zero area or outside the window carries no information.";
  if(m.kind==="merge") return `The two lines are closer than ${AB_TH.minSep} × their mean FWHM (cap); removing the smaller one costs ΔAICc_eff ${fx(-d.dAICc,1)}, ΔBIC_eff ${fx(-d.dBIC,1)} (< ${T}).`;
  if(m.kind==="bound") return `Removing either line costs ≥ ${T} in both criteria, so both stay with BE ranges held ≥ ${AB_TH.minSep} × FWHM apart (cap).`;
  if(m.kind==="multi") return "Multi-start found a lower minimum of the same model.";
  if(m.kind==="drop") return m.info.recipe?`Recipe line: the model without it is not worse (ΔAICc_eff ${sg(d.dAICc,1)} > −${AB_TH.simp}).`:`Without this line ΔAICc_eff ${sg(d.dAICc,1)}, ΔBIC_eff ${sg(d.dBIC,1)}: it does not pass the add criterion (both ≥ ${T}), so it is not supported.`;
  if(m.cls==="simp") return `Fewer parameters (${d.dp}) and not worse: ΔAICc_eff ${sg(d.dAICc,1)} > −${AB_TH.simp}.`;
  return `ΔAICc_eff ${sg(d.dAICc,1)} and ΔBIC_eff ${sg(d.dBIC,1)} both ≥ ${T} (N_eff ${fx(d.Neff,0)}, ρ₁ ${fx(d.rho,2)}).`; }

/* 처음 상태: 지금 모형을 맞춤 (first: 단계 피팅 — 보통 피팅보다 나빠지지 않음, §17.4) */
function startState(prob,o){ const C=CORE(), G=GD(); const P=clonePs(prob);
  if(o&&o.first){ const s=G.stagedFit(P,{}); const st=C.fitStats(s.P,s.r); return {P:s.P,r:s.r,st}; }
  const T=fitS(P,{init:false}); if(!T) throw new Error("the current model could not be fitted"); return T; }

/* 한 걸음 제안 (화면: 하나씩). args = {prob, ctx:{shift, el, orb, present, hv}, excluded:[keys], first} */
function autoStep(prob,ctx,o){ o=o||{}; ctx=ctx||{}; const excl=new Set(o.excluded||[]);
  const S=startState(prob,o); const s=stepFrom(S,ctx,excl);
  const cur={chi2nu:S.st.chi2nu,N:S.st.N,p:S.st.p,rho1:S.st.rho1,nLines:S.P.comps.length};
  if(!s.best) return {step:null,stop:stopWhy(s),cur,rejected:s.rejected,curModel:{comps:S.P.comps,cons:S.P.cons,bg:S.P.bg}};
  return {step:stepRec(S,s.best,ctx,s),cur,rejected:s.rejected}; }
function stopWhy(s){ const capd=s.rejected.filter(x=>x.cap); return capd.length?`no candidate move passes (best addition violated a cap: ${capd[0].cap})`:`no candidate move passes (${s.cands} candidates evaluated)`; }

/* 끝까지 (받아들인다고 보고) — 벤치마크와 'Accept remaining' */
function autoBuild(prob,ctx,o){ o=o||{}; ctx=ctx||{}; const excl=new Set(o.excluded||[]); const t0=Date.now();
  let S=startState(prob,{first:o.first!==false}); const steps=[]; let stop=null, msDone=false; const maxS=o.maxSteps||AB_TH.maxSteps;
  for(let it=0;;it++){
    if(steps.length>=maxS){ stop=`maximum number of steps (${maxS}) reached`; break; }
    const s=stepFrom(S,ctx,excl);
    if(!s.best){ if(o.multi!==false&&!msDone){ msDone=true; const ad=multiAdopt(S,o); if(ad){ steps.push(stepRec(S,ad,ctx,s)); S=ad.S2; continue; } } stop=stopWhy(s); break; }
    steps.push(stepRec(S,s.best,ctx,s)); if(s.best.kind==="bound") excl.add(s.best.key); S=s.best.S2; S.cov=null; }
  const fin=o.final===false?null:finalAnalysis(S,ctx,o);
  return {steps,stop,S,final:fin,ms:Date.now()-t0}; }

/* ---------- 다중 시작 (유일성) ---------- */
/* 시작점: 비선형 자유 매개변수를 Latin hypercube 로 — 범위 = 사용자 한계 ∩ 이웃(BE ± max(FWHM, 0.5 eV), 폭 × 0.5–2, 나머지는 한계 또는 ± 50 %).
   창 전체에서 시작하면 최적화기의 수렴 영역을 재지 해의 유일성(같은 깊이의 다른 해)을 재지 않음(개발 중: 깨끗한 단일 선 H 30000 이 3/12) */
function localStarts(S,N,seed){ const C=CORE(); const P=S.P, map=S.r.map||C.buildMap(P), R=C.rng(seed||1); const v=C.getVec(P,map);
  const idx=map.free.map((f,i)=>i).filter(i=>map.free[i].k!=="area"&&map.free[i].ci>=0);
  const rngOf=i=>{ const f=map.free[i], c=P.comps[f.ci], x=v[i]; let a=f.lo, b=f.hi;
    if(f.k==="be"){ const w=Math.max(fwhmOf(c),0.5); a=Math.max(a,x-w); b=Math.min(b,x+w); }
    else if(["F","wG","wL"].includes(f.k)){ a=Math.max(a,x/2); b=Math.min(b,x*2); }
    else { if(!isFinite(a)) a=x-0.5*Math.abs(x||1); if(!isFinite(b)) b=x+0.5*Math.abs(x||1); }
    return [a,b]; };
  const perm=idx.map(()=>{ const a=[...Array(N).keys()]; for(let i=N-1;i>0;i--){ const j=Math.floor(R.u()*(i+1)); [a[i],a[j]]=[a[j],a[i]]; } return a; });
  const pts=[]; for(let s=0;s<N;s++){ if(s===0){ pts.push(null); continue; } const o={}; idx.forEach((fi,q)=>{ const [a,b]=rngOf(fi); const t=(perm[q][s]+R.u())/N; o[map.free[fi].ci+"|"+map.free[fi].k]=a+(b-a)*t; }); pts.push(o); }
  return pts; }
function multiStart(S,o){ const C=CORE(); const N=(o&&o.starts)||AB_TH.starts; const pts=localStarts(S,N,(o&&o.seed)||1); const res=[];
  // 시작점끼리는 끝난 해의 χ²_d(데이터 분산, §3.5 — 모든 해에 같은 가중)로 비교. 1차 χ²(chi2p1)는 쓰지 않음: 먼 시작에서 1차가 일찍 끝나
  // 같은 최종 해(χ² 188)인데 chi2p1 이 188 · 1310 · 13 594 로 갈렸다(개발 중 seed 1021) — 같은 최솟값을 '다르다'고 셈
  pts.forEach(pt=>{ const Q=clonePs(S.P); C.applyStart(Q,pt); try{ const r=C.fitOnce(Q,{}); const st=C.fitStats(Q,r); res.push({x:isFinite(st.chi2d)?st.chi2d:Infinity,Q,r,st}); }catch(e){ res.push({x:Infinity}); } });
  // v1.11: 같은 최솟값 계산은 core.sameMinimum (Fit 의 Multi-start 와 같은 함수); st.chi2d = core.chi2dOf = startChi2 와 같은 값
  const sm=C.sameMinimum(res.map(q=>q.x),S.st.N,(S.r.map||C.buildMap(S.P)).free.length); const best=res[sm.bi<0?0:sm.bi];
  return {N,near:sm.near,frac:sm.near/N,thr:sm.thr,chi2nu:sm.nu,best,cur:S.st.chi2d}; }
function multiAdopt(S,o){ const ms=multiStart(S,o); S.ms=ms; if(!(ms.best&&ms.best.Q&&ms.best.x<ms.cur-ms.thr)) return null; const C=CORE();
  const T={P:ms.best.Q,r:ms.best.r,st:ms.best.st}; T.ms=null; const d=cmpIC(S,T);
  return {kind:"multi",key:"multi",cls:"cap",label:`Adopt the lower minimum found by multi-start (χ²_d ${fx(ms.cur,1)} → ${fx(ms.best.x,1)})`,S2:T,d,info:{},sameMin:{near:ms.near,N:ms.N}}; }


/* ---------- 숨은 이웃 선의 치우침 한계 (SPEC §18.3) ----------
   모형이 볼 수 없는 선(검출 문턱 아래)이 맞춘 선의 위치·넓이를 얼마나 끌 수 있나 — 선 i 옆 ±{0.5, 0.75, 1, 1.25, 1.5} × FWHM 에 같은 모양의 선 하나를 가정:
   그 선의 (가중) 모양 s 를 모형 Jacobian J 의 공간에 사영 → 남는 부분 s⊥ 만 χ² 를 줄일 수 있으므로, 더하기 규칙(ΔAICc_eff, ΔBIC_eff ≥ 10; 선 하나 = 매개변수 3)을
   통과하지 못하는 가장 큰 넓이 a = √(Δχ²_문턱 / |s⊥|²) (단 ≤ 선 i 의 넓이). 나머지 J·c 는 매개변수가 흡수 → 치우침 a·c (선형). 최댓값을 '숨은 선 한계'로 보고하고
   위치·넓이 불확도 u = σ + 한계/3 (3u = 3σ + 최악 치우침). 개발 중 측정: 확립된 경우의 위치 실패 27 중 24 가 1.5 FWHM 안의 놓친 참 선 — σ 만으로는 덮지 못함 */
function solveSym(A,b,m){ const M=Array.from({length:m},(_,i)=>Array.from({length:m+1},(_,j)=>j<m?A[i*m+j]:b[i])); const sc=Array.from({length:m},(_,i)=>Math.sqrt(Math.max(A[i*m+i],1e-300)));
  for(let i=0;i<m;i++){ for(let j=0;j<m;j++) M[i][j]/=sc[i]*sc[j]; M[i][m]/=sc[i]; M[i][i]+=1e-10; }
  for(let c=0;c<m;c++){ let pv=c; for(let r=c+1;r<m;r++) if(Math.abs(M[r][c])>Math.abs(M[pv][c])) pv=r; [M[c],M[pv]]=[M[pv],M[c]]; const d=M[c][c]; if(Math.abs(d)<1e-14) return null;
    for(let r=0;r<m;r++){ if(r===c) continue; const f=M[r][c]/d; if(!f) continue; for(let k=c;k<=m;k++) M[r][k]-=f*M[c][k]; } }
  return Array.from({length:m},(_,i)=>M[i][m]/M[i][i]/sc[i]); }
function hiddenBound(S,rho){ const C=CORE(); const P=S.P, E=P.E, n=E.length, map=S.r.map||C.buildMap(P), v=S.r.v; const w=Float64Array.from(v,x=>1/Math.sqrt(x));
  const p0=C.getVec(P,map), np=p0.length; const ev=C.makeEvaluator(P,map); ev.full(); const J=[];
  map.free.forEach((f,j)=>{ let h=Math.max(Math.abs(p0[j])*1e-5,f.k==="be"||f.k==="dE"?1e-4:1e-6); if(f.k==="area") h=Math.max(Math.abs(p0[j])*1e-5,1e-9);
    let a=p0[j]-h, b=p0[j]+h; if(b>f.hi) b=p0[j]; if(a<f.lo) a=p0[j]; const q=p0.slice(); q[j]=b; C.setVec(P,map,q); const mb=ev.partial(j); q[j]=a; C.setVec(P,map,q); const ma=ev.partial(j); C.setVec(P,map,p0);
    const col=new Float64Array(n); if(b>a) for(let i=0;i<n;i++) col[i]=(mb[i]-ma[i])/(b-a)*w[i]; J.push(col); });
  const A=new Float64Array(np*np); for(let a=0;a<np;a++) for(let b=0;b<=a;b++){ let t=0; for(let i=0;i<n;i++) t+=J[a][i]*J[b][i]; A[a*np+b]=t; A[b*np+a]=t; }
  const st=S.st, Neff=st.N*(1-rho)/(1+rho), fN=st.N/Neff, phi=st.chi2nu; const d2=Math.max(AB_TH.T+2*3,AB_TH.T+3*Math.log(Neff))*phi*fN;
  const ctx={emin:E[0],emax:E[n-1]}; const out=P.comps.map(()=>null);
  P.comps.forEach((c,i)=>{ const jb=map.free.findIndex(f=>f.ci===i&&f.k==="be"), ja=map.free.findIndex(f=>f.ci===i&&f.k==="area"); if(jb<0&&ja<0) return; const F=fwhmOf(c); let best={be:0,area:0,d:null,a:0};
    [-1.5,-1.25,-1,-0.75,-0.5,0.5,0.75,1,1.25,1.5].forEach(k=>{ const be=c.be+k*F; if(be<E[0]||be>E[n-1]) return; const q=C.clone(c); q.be=be; q.area=1; let s; try{ s=C.compEval(q,E,ctx); }catch(e){ return; }
      const sw=Float64Array.from(s,(x,i)=>x*w[i]); const rhs=new Float64Array(np); for(let a=0;a<np;a++){ let t=0; for(let i=0;i<n;i++) t+=J[a][i]*sw[i]; rhs[a]=t; }
      const cf=solveSym(A,rhs,np); if(!cf) return; let r2=0; for(let i=0;i<n;i++){ let m=0; for(let a=0;a<np;a++) m+=J[a][i]*cf[a]; const d=sw[i]-m; r2+=d*d; }
      const ad=Math.min(Math.sqrt(d2/Math.max(r2,1e-300)),Math.max(0,c.area)); const bb=jb>=0?Math.abs(ad*cf[jb]):0, ba=ja>=0?Math.abs(ad*cf[ja]):0;
      if(bb>best.be) best.be=bb; if(ba>best.area) best.area=ba; if(bb>=best.be&&bb>0){ best.d=k; best.a=ad; } });
    out[i]={be:jb>=0?best.be:null,area:ja>=0?best.area:null,d:best.d,aRel:c.area>0?best.a/c.area:null}; });
  return out; }
/* ---------- 끝 분석: 정직한 불확도 (SPEC §18.3) ---------- */
function finalAnalysis(S,ctx,o){ const C=CORE(); o=o||{}; const sh=(ctx&&ctx.shift)||0; const P=S.P, r=S.r;
  P.comps.forEach(c=>{ if(c.ab&&/^-?\d+\.\d eV line$/.test(c.name||"")) c.name=(c.be+sh).toFixed(1)+" eV line"; });
  const fin=C.finalize(P,r); const st=fin.stats; const rho=st.rho1; const infl=rho>0?Math.sqrt((1+rho)/(1-Math.min(rho,AB_TH.rhoMax))):1;
  const ms=o.multi===false?null:(S.ms||multiStart(S,o)); const cov=fin.cov; const np=cov.np;
  // 묶음: 서로 다른 두 선의 BE·넓이·폭 사이 |r| > 0.9
  const par=P.comps.map((c,i)=>i); const fnd=i=>par[i]===i?i:(par[i]=fnd(par[i])); const pairs=[];
  for(let a=0;a<np;a++) for(let b=a+1;b<np;b++){ const qa=cov.names[a], qb=cov.names[b]; if(qa.ci<0||qb.ci<0||qa.ci===qb.ci) continue; if(!["be","area","F"].includes(qa.k)||!["be","area","F"].includes(qb.k)) continue;
    const v=cov.corr[a*np+b]; if(isFinite(v)&&Math.abs(v)>AB_TH.corrGroup){ pairs.push({a:qa.ci,b:qb.ci,ka:qa.k,kb:qb.k,r:v}); par[fnd(qa.ci)]=fnd(qb.ci); } }
  const gm={}; P.comps.forEach((c,i)=>{ const g=fnd(i); (gm[g]=gm[g]||[]).push(i); });
  const map=r.map||C.buildMap(P); const groups=Object.values(gm).filter(g=>g.length>1).map(g=>{ let s=null; try{ s=C.areaSumSigma(P,map,{atB:cov.atB,sig:cov.sig,C:Float64Array.from(cov.C)},g); }catch(e){}
    const A=g.reduce((t,i)=>t+P.comps[i].area,0); return {idx:g,ids:g.map(i=>P.comps[i].id),names:g.map(i=>nameOf(P.comps[i],sh)),area:A,sig:s!=null?s*infl:null,
      be:g.reduce((t,i)=>t+P.comps[i].area*P.comps[i].be,0)/(A||1)+sh,pairs:pairs.filter(q=>g.includes(q.a)).map(q=>({a:nameOf(P.comps[q.a],sh),b:nameOf(P.comps[q.b],sh),ka:q.ka,kb:q.kb,r:q.r}))}; });
  /* 넓은 선 (SPEC §18.3): 가장 좁은 '잘 정해진' 선(σ_FWHM/FWHM < 10 %, 한계 아님)보다 FWHM 이 broadRel 배 넘게 넓고 그 차이가 3σ 넘음 →
     '풀리지 않은 선이 들어 있을 수 있음' — 위치·넓이는 그 덩어리(envelope)의 값으로만 보고. 개발 중 측정: 확립된 경우의 위치 실패 27 중 24 가
     1.5 FWHM 안의 놓친 참 선(분해 한계 아래) — 모형이 볼 수 없는 것이므로 주장 범위를 줄임 */
  const hbS=new Set(boundHits(S).map(h=>h.ci+"|"+h.k)); const fwS=P.comps.map((c,i)=>{ const q=fin.summ[i]||{}; return {w:q.fwhm,s:q.sig&&q.sig.fwhm}; });
  const wd=fwS.filter((x,i)=>x.w>0&&x.s>0&&x.s/x.w<AB_TH.floorSig&&!hbS.has(i+"|F")).map(x=>x.w); const Fref=wd.length?Math.min(...wd):null;
  const broad=i=>{ const x=fwS[i]; return Fref!=null&&x.w>AB_TH.broadRel*Fref&&(x.s==null||x.w-Fref>3*x.s*infl); };
  const inG=new Set([].concat(...groups.map(g=>g.idx))); const tot=P.comps.reduce((t,c)=>t+Math.max(0,c.area),0)||1; const hb=boundHits(S);
  let hid=null; try{ hid=hiddenBound(S,clamp(rho,0,RHO_CAP)); }catch(e){ hid=null; }
  const lines=P.comps.map((c,i)=>{ const q=fin.summ[i]||{}, s=q.sig||{}; const hb_=hid&&hid[i]; const bS=s.be!=null?s.be*infl:null, aS=s.area!=null?s.area*infl:null;
    return {id:c.id,name:nameOf(c,sh),be:c.be+sh,beSig:bS,area:c.area,areaSig:aS,areaPct:100*Math.max(0,c.area)/tot,fwhm:q.fwhm,
    hidden:hb_,beU:bS!=null?bS+(hb_&&hb_.be!=null?hb_.be/3:0):null,areaU:aS!=null?aS+(hb_&&hb_.area!=null?hb_.area/3:0):null,
    recipe:c.recipe||null,added:!!c.ab,separate:!inG.has(i),broad:broad(i),fwRef:Fref,bounds:hb.filter(h=>h.ci===i).map(h=>h.k+" "+h.side)}; });
  // 선마다 지지(빼면 얼마나 나빠지나) — 보고용
  if(o.support!==false&&P.comps.length>1) lines.forEach((L,i)=>{ const Q=clonePs(P); try{ removeIds(Q,[L.id]); const T=fitS(Q,{}); if(T){ const d=cmpIC(T,S); L.support={dAICc:d.dAICc,dBIC:d.dBIC}; L.weak=!(d.dAICc>=AB_TH.T&&d.dBIC>=AB_TH.T); } }catch(e){} });
  const reasons=[]; if(fin.conv===false||fin.stop==="stalled"||fin.stop==="maxIter") reasons.push(`fit ${fin.stop||"not converged"}`);
  if(st.chi2nu>AB_TH.chi2nuMax) reasons.push(`χ²ν ${fx(st.chi2nu,2)} > ${AB_TH.chi2nuMax} (the model does not describe the data)`);
  // v1.10 검증 뒤 조임(§18.5): 위치·폭이 한계에 닿으면 그 값은 자료가 아니라 한계가 정함(R1) → 확립 안 됨
  const hbF=hb.filter(h=>h.k==="be"||h.k==="F"); if(hbF.length) reasons.push(`${hbF.map(h=>`${nameOf(P.comps[h.ci],sh)} ${h.k==="F"?"FWHM":"BE"} at its ${h.side} bound`).join(", ")} — set by the bound, not by the data`);
  { const fw=P.comps.map(fwhmOf); const bad=[]; for(let i=0;i<P.comps.length;i++) for(let j=i+1;j<P.comps.length;j++){ const a=P.comps[i], b=P.comps[j]; if(!ownFree(P,a,"be")||!ownFree(P,b,"be")||(a.recipe&&a.recipe===b.recipe)) continue; if(!sepOK(a.be,fw[i],b.be,fw[j])) bad.push(`${nameOf(a,sh)} / ${nameOf(b,sh)}`); }
    if(bad.length) reasons.push(`lines closer than ${AB_TH.minSep} × FWHM: ${bad.join(", ")}`); }
  if(isStatic(P.bg)) reasons.push(`static ${P.bg.type} background: its end levels are not fitted, so every σ excludes the background error (Auto build's in-fit background step was not taken)`);
  if(rho>RHO_CAP) reasons.push(`residual ρ₁ ${fx(rho,2)} > ${RHO_CAP}: neighbouring residuals more correlated than the comparison assumed (N_eff overestimated) or structure left in the residual`);
  if(ms&&ms.frac<AB_TH.sameMin) reasons.push(`same minimum in only ${ms.near}/${ms.N} starts (< 1/3; threshold Δχ² < max(1, χ²ν) = ${fx(ms.thr,2)})`);
  const species=(ctx&&ctx.el)?lines.filter(L=>!L.recipe).map(L=>({id:L.id,be:L.be,cands:candidateSpecies(ctx.el,ctx.orb,L.be,ctx.present,ctx.hv)})):[];
  return {established:!reasons.length,reasons,chi2nu:st.chi2nu,rho1:rho,infl,N:st.N,p:st.p,multi:ms?{N:ms.N,near:ms.near,frac:ms.frac,thr:ms.thr}:null,lines,groups,species,fin}; }

/* ---------- 후보 화학종 (앱 안의 출처 있는 표만, 시료에 있는 원소로 거름; 이름은 사용자가 붙임) ---------- */
function candidateSpecies(el,orb,be,present,hv){ const G=GD(); if(!G||!el) return []; const out=[]; const P=present instanceof Set?present:new Set(present||[]);
  const add=(name,lo,hi,src)=>{ if(be>=lo&&be<=hi) out.push({name,lo:+lo.toFixed(2),hi:+hi.toFixed(2),src}); };
  G.recipesFor(el,orb).forEach(R=>{ const cite=String(R.cite||"").split(";")[0];
    try{
      if(R.kind==="c1s"){ R.lines.concat([R.sp2Line,R.piLine]).forEach(l=>{ if(l.ref) add(l.name,(l.key==="sp2"?R.sp2ref:R.ref)-R.beTol,(l.key==="sp2"?R.sp2ref:R.ref)+R.beTol,cite); else if(l.rng) add(l.name,R.ref+l.rng[0],R.ref+l.rng[1],cite); }); }
      else if(R.kind==="o1s"){ R.lines.forEach(l=>{ if(l.rng&&!l.lat) add(l.name,l.rng[0],l.rng[1],cite); else if(l.rng) add(l.name,l.rng[0],l.rng[1],cite); });
        Object.keys(G.O1S_LAT).forEach(m=>{ if(P.has(m)) add(`lattice O (${m})`,G.O1S_LAT[m].be-0.5,G.O1S_LAT[m].be+0.5,G.O1S_LAT[m].src); }); }
      else if(R.kind==="env"){ R.species.forEach(s=>(s.sets||[s]).forEach(t=>{ const b=t.peaks[0][0]; add(s.label,b-R.beTol,b+R.beTol,cite); })); }
      else if(R.species){ R.species.forEach(s=>(s.lines||[]).forEach(l=>{ if(l.be!=null) add(`${s.label} ${l.name}`,l.be-(R.beTol||0.3),l.be+(R.beTol||0.3),cite); })); }
      else if(R.lines){ R.lines.forEach(l=>{ if(l.be!=null) add(`${l.name||l.key} (${l.state||""})`,l.be-(R.beTol||0.6),l.be+(R.beTol||0.6),cite); }); }
    }catch(e){} });
  try{ G.overlapFindings(el,orb,P,{hv,lo:be-1,hi:be+1}).forEach(f=>out.push({name:f.label,lo:null,hi:null,src:"overlap table (SPEC §17.8)"})); }catch(e){}
  return out; }

/* ---------- 기록 내보내기 (Markdown) ---------- */
function logText(L){ const lines=[`# XPS Auto build log — ${L.region||""}`,"",`Started ${L.when||""}. Rules: SPEC §18 (T = ${AB_TH.T} for ΔAICc_eff and ΔBIC_eff, simplify when ΔAICc_eff > −${AB_TH.simp}, min. separation ${AB_TH.minSep} FWHM, N_eff = N(1 − ρ₁)/(1 + ρ₁)).`,"",
  "| # | decision | step | Δχ²ν | ΔAICc_eff | ΔBIC_eff | Δp | ρ₁ | residual run | bounds |","|---|---|---|---|---|---|---|---|---|---|"];
  (L.steps||[]).forEach((s,i)=>{ const e=s.ev||{}; lines.push(`| ${i+1} | ${s.decision||"proposed"} | ${s.label} | ${sg(e.dChi2nu,3)} | ${sg(e.dAICc,1)} | ${sg(e.dBIC,1)} | ${e.dp!=null?e.dp:"—"} | ${fx(e.rho,2)} | ${e.run?fx(e.run.be,2)+" eV (z "+fx(e.run.z,1)+")":"—"} | ${(e.bounds||[]).join("; ")||"—"} |`); });
  lines.push("",`Stopped: ${L.stop||"—"}`); const F=L.final;
  if(F){ lines.push("",`## Result: ${F.established?"established":"MODEL NOT ESTABLISHED"}`); if(!F.established) F.reasons.forEach(r=>lines.push("- "+r));
    lines.push("",`χ²ν ${fx(F.chi2nu,3)} · ρ₁ ${fx(F.rho1,2)} (σ × ${fx(F.infl,2)}) · multi-start same minimum ${F.multi?F.multi.near+"/"+F.multi.N:"—"}`,"","| line | BE (eV) | ± | area % | FWHM | note |","|---|---|---|---|---|---|");
    F.lines.forEach(l=>lines.push(`| ${l.name} | ${fx(l.be,2)} | ${l.separate?fx(l.beSig,2):"—"} | ${l.separate?fx(l.areaPct,1):"not separately determined"} | ${fx(l.fwhm,2)} | ${[l.weak?"weak support":"",l.bounds.length?"at bound: "+l.bounds.join(", "):""].filter(Boolean).join("; ")||""} |`));
    F.groups.forEach(g=>lines.push("",`Identifiable combination: ${g.names.join(" + ")} — area ${fx(g.area,0)} ± ${fx(g.sig,0)} (${fx(100*g.sig/g.area,1)} %), centroid ${fx(g.be,2)} eV`));
    (F.species||[]).forEach(s=>{ if(s.cands.length) lines.push("",`Candidates for the ${fx(s.be,1)} eV line (not assigned): `+s.cands.map(c=>`${c.name}${c.lo!=null?` [${c.lo}–${c.hi}]`:""} (${c.src})`).join("; ")); }); }
  return lines.join("\n"); }

/* Worker 작업: {mode:"step"|"run"|"final", prob, ctx, excluded, first} — 결과는 직렬화 가능 */
function abJob(a){ const C=CORE(); const ctx=a.ctx||{};
  if(a.mode==="step"){ return autoStep(a.prob,ctx,{excluded:a.excluded,first:a.first}); }
  if(a.mode==="run"){ const R=autoBuild(a.prob,ctx,{excluded:a.excluded,first:a.first,final:false,multi:false}); return {steps:R.steps,stop:R.stop,ms:R.ms}; }
  if(a.mode==="final"){ const S=startState(a.prob,{});   // noAdopt: 사용자가 'Adopt' 를 버렸음 → 지금 모형 그대로 분석
    if(!a.noAdopt){ const ad=multiAdopt(S,{}); if(ad) return {step:stepRec(S,ad,ctx,{runs:null}),final:null,fin:null}; }
    const F=finalAnalysis(S,ctx,{}); const fin=F.fin; delete F.fin; return {step:null,final:F,fin:Object.assign(fin,{cons:S.P.cons})}; }
  throw new Error("unknown Auto build mode"); }

const XPS_AUTO={AB_TH,rhoNoise,icEff,cmpIC,caps,boundHits,resRuns,candidates,passes,stepFrom,localStarts,hiddenBound,autoStep,autoBuild,multiStart,finalAnalysis,candidateSpecies,logText,abJob,startState};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_AUTO; else root.XPS_AUTO=XPS_AUTO;
})(typeof self!=="undefined"?self:this);
