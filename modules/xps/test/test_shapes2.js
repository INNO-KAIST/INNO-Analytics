/* CasaXPS Cookbook 2026 선모양 (SPEC_shapes2.md) 시험
   node test/test_shapes2.js   (modules/xps 에서; QUICK=1 이면 회복 시험을 줄임)
   1) 값·넓이·FWHM ↔ 독립 Python 기준 (test/ref/ref_shapes2.json: quad·FFT 합성곱·누적 Simpson·mpmath)
      허용: 값 ≤ 1e-5 (최댓값 대비), 넓이 ≤ 1e-6, FWHM ≤ 1e-5 (상대)
   2) 표기 왕복 · 오류 · 하위 호환(AGT GL(m)T(k) ≡ GLN tail T)
   3) 합성 잡음 스펙트럼 회복: 매개변수가 3σ 안 ≥ 90 % (모양마다 200 회; QUICK 40 회)
   4) 속도: 10 성분 × 2000 채널 평가·피팅 시간 · 4b 퇴화 매개변수(TP wL)까지 풀었을 때 멈춤 규칙(SPEC §3.3b) */
const path=require("path"), fs=require("fs");
const C=require("../src/xps_core.js");
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref/ref_shapes2.json"),"utf8"));
const QUICK=!!process.env.QUICK;
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const F64=a=>Float64Array.from(a), e1=x=>x==null?"—":(+x).toExponential(1);
const label=r=>`${r.p.mod?r.p.mod+" ":""}${r.shape}${r.p.tail?"+"+r.p.tail:""} ${JSON.stringify(Object.fromEntries(Object.entries(r.p).filter(([k])=>!["mod","tail"].includes(k)).map(([k,v])=>[k,+(+v).toPrecision(4)])))}${r.n!=null?" n="+r.n:""}`;
/* 기준 사례 → INNO 성분 (n 이 있으면 casaN 으로 — n → 폭 규칙도 함께 시험) */
function compOf(r){ const p=r.p, c={id:"t",shape:r.shape,be:0,area:1,doublet:null}; Object.assign(c,p);
  if(r.n!=null){ c.casaN=r.n; } if(p.mod==="ST"){ c.stMu=p.stMu; c.stG=p.stG||0; }
  if(c.wG==null&&["TLA","GLN","SGLN","LA","DS"].includes(r.shape)) c.wG=0; if(r.shape==="TP"&&c.wL==null) c.wL=0; return c; }

console.log("== 1 값 · 넓이 · FWHM vs 독립 기준 (scipy quad / FFT / mpmath)");
let worstV=0, worstA=0, worstF=0; const table=[];
REF.cases.forEach(r=>{ const c=compOf(r), sp=C.shapeParams(c,1);
  if(r.n!=null){ const k=r.shape==="TP"?"wL":"wG"; ok(Math.abs(sp[k]-r.p[k])<1e-12,`${label(r)}: CasaXPS n = ${r.n} → ${k} = n·${r.shape==="DS"?"wL":"F"}/144.87`,`${sp[k]} vs ${r.p[k]}`); }
  const ctx=r.win?{lo:r.win[0],hi:r.win[1]}:null;
  const v=C.shapeUnit(r.shape,F64(r.x),sp,ctx), mx=Math.max(...r.v.map(Math.abs)); let ev=0; r.v.forEach((t,i)=>{ ev=Math.max(ev,Math.abs(v[i]-t)/mx); });
  let ea=null;
  if(r.win){ // 창 넓이 1 (DS·SB)
    const n=40001, h=(r.win[1]-r.win[0])/(n-1), xs=new Float64Array(n); for(let i=0;i<n;i++) xs[i]=r.win[0]+i*h; const u=C.shapeUnit(r.shape,xs,sp,ctx);
    let s=u[0]+u[n-1]; for(let i=1;i<n-1;i++) s+=(i%2?4:2)*u[i]; ea=Math.abs(s*h/3-1);
    if(r.shape==="DS"&&r.p.dsl>0) ea=ea; }
  else if(r.area_in!=null&&r.R_in){ const R=r.R_in, kink=r.shape==="TP"&&!(r.p.wL>0);
    const seg=(a,b,n)=>{ const h=(b-a)/n, xs=new Float64Array(n+1); for(let i=0;i<=n;i++) xs[i]=a+i*h; const u=C.shapeUnit(r.shape,xs,sp,{lo:-R,hi:R}); let s=u[0]+u[n]; for(let i=1;i<n;i++) s+=(i%2?4:2)*u[i]; return s*h/3; };
    const W=r.p.F||r.p.wL||1; const area=kink?seg(-r.p.F/2+1e-12,r.p.F/2-1e-12,20000):seg(-R,-6*W,20000)+seg(-6*W,6*W,40000)+seg(6*W,R,20000);
    ea=Math.abs(area-r.area_in); }
  let ef=null; const tf=C.trueFWHM(r.shape,sp); if(r.fwhm!=null) ef=Math.abs(tf.fwhm-r.fwhm)/r.fwhm;
  worstV=Math.max(worstV,ev); if(ea!=null) worstA=Math.max(worstA,ea); if(ef!=null) worstF=Math.max(worstF,ef);
  table.push([label(r),ev,ea,ef]);
  ok(ev<=1e-5&&(ea==null||ea<=1e-6)&&(ef==null||ef<=1e-5),`${label(r)}${r.win?` · 창 [${r.win}]`:""}`,`값 ${e1(ev)} · 넓이 ${e1(ea)} · FWHM ${e1(ef)}${r.fwhm!=null?` (${tf.fwhm.toFixed(6)} vs ${r.fwhm.toFixed(6)})`:""}`); });
console.log(`     최대: 값 ${e1(worstV)} · 넓이 ${e1(worstA)} · FWHM ${e1(worstF)}`);
{ let e=0; REF.tla_mp.forEach(t=>{ const B=C.s2Base("TLA",{F:1,alpha:t.alpha,mu:t.mu,wG:0}); e=Math.max(e,Math.abs(B.left-t.left_frac)); const Z=2*(Math.PI/2+Math.PI/t.mu)*C.tlaI1(t.alpha); e=Math.max(e,Math.abs(Z-t.A)/t.A); });
  ok(e<1e-12,"TLA 정규화 상수(닫힌 꼴 + J(α) 구적) = mpmath 30 자리",e1(e)); }
{ const r=REF.cases.filter(c=>c.p.mod==="ST"); let e=0; r.forEach(c=>{ e=Math.max(e,Math.abs(c.N_numeric/c.N_identity-1)); });
  ok(e<1e-9,"ST 정규화: 수치 ∫S(μ)F = A²(½ + 1/μ) (기준 쪽 확인 — 엔진은 이 항등식을 씀)",e1(e)); }

console.log("== 2 CasaXPS 표기 왕복 · 오류 · 하위 호환");
{ const strs=["TLA(1.5,2,50)","ST(1.7)LA(1,200)","ST(1)LA(1,143)","ST(1.5,30)GL(30)","DS(0.1,50)","DS(0.3,600,40)","GL(30,50)","SGL(30,50)T(2)","GL(30)PHI(1,2)","GL(30)T(2)","TP(1,15,0)","TP(1,1,0)","SB(0)LA(1,43)","SB(0)GL(30,20)T(1.5)","ST(2)TLA(1.2,3,40)"];
  strs.forEach(s=>{ const a=C.casaToComp(s,{fwhm:1.3,be:100,area:1000,strict:true}); const n1=C.compToCasa(a.comp).notation; const b=C.casaToComp(n1,{fwhm:1.3,be:100,area:1000});
    const keys=["shape","mod","tail","casaN","stMu","stN","alpha","mu","m","tk","tj","ta","tb","dsl","beta"]; const same=keys.every(k=>(a.comp[k]==null&&b.comp[k]==null)||Math.abs((+a.comp[k]||0)-(+b.comp[k]||0))<1e-9||a.comp[k]===b.comp[k]);
    const E=Float64Array.from({length:401},(_,i)=>95+0.025*i), va=C.compEval(a.comp,E,{emin:95,emax:105}), vb=C.compEval(b.comp,E,{emin:95,emax:105}); let d=0,m=0; va.forEach((v,i)=>{ d=Math.max(d,Math.abs(v-vb[i])); m=Math.max(m,Math.abs(v)); });
    ok(same&&d/m<1e-6,`'${s}' → ${a.comp.shape}${a.comp.mod?" + "+a.comp.mod:""} → '${n1}' → 같은 성분`,`곡선 차 ${e1(d/m)}`); });
  const bad=["XYZ(1)","ST(1)","GL(30)T(2)PHI(1,1)","TLA(0.4,2)","LA(1,1)T(2)","ST(0)LA(1,100)","DS(1.2,10)","GL(30","TP(0,0,0)"]; let msgs=[];
  bad.forEach(s=>{ let m=null; try{ C.casaToComp(s,{strict:true}); }catch(e){ m=e.message; } msgs.push(m); });
  ok(msgs.every(m=>m&&m.length>20)&&C.casaToComp("XYZ(1)")===null,"알 수 없는/잘못된 표기 → strict: 분명한 오류 문구 · 기본: null (이전 동작)",msgs[0]);
  // GL(m)T(k) — 이전 AGT(tu 'fwhm', TM 0) 와 같은 모양 (정규화만 다름: AGT = 창 넓이, GLN = ℝ 넓이)
  const E=Float64Array.from({length:801},(_,i)=>280+0.025*i), ctx={emin:280,emax:300};
  const ag=C.compEval(C.newComp("AGT",289,{F:1.2,m:0.3,ET:2,TM:0,CT:0,tu:"fwhm"}),E,ctx), gn=C.compEval(Object.assign(C.casaToComp("GL(30)T(2)",{be:289,fwhm:1.2}).comp,{F:1.2}),E,ctx);   // 같은 F 로 (CasaXPS fwhm 1.2 → F 1.22, SPEC_casa_v3 §3)
  let sA=0,sG=0; for(let i=1;i<E.length;i++){ sA+=0.5*(ag[i]+ag[i-1])*0.025; sG+=0.5*(gn[i]+gn[i-1])*0.025; } let d=0; ag.forEach((v,i)=>{ d=Math.max(d,Math.abs(v/sA-gn[i]/sG)); });
  ok(d<1e-6,"GL(30)T(2): 새 GLN(tail T) 모양 = 이전 AGT(tu 'fwhm') 모양 (창 안 넓이로 맞춘 뒤)",e1(d));
  const old=C.newComp("AGT",289,{F:1.2,m:0.3,ET:2,TM:0,CT:0,tu:"fwhm"}); ok(C.compToCasa(old).notation==="GL(30)T(2)"&&C.paramNames(old).join()==="be,area,F,m,ET,TM,CT","저장된 AGT 성분은 그대로 읽힘 (표기·매개변수 이름 그대로)");
  const st=C.newComp("LA",0,{F:1,alpha:1,beta:1,wG:0.3,mod:"ST",stMu:1}); ok(C.paramNames(st).join()==="be,area,F,alpha,beta,wG,stMu,stG"&&C.buildMap({comps:[st],cons:{},bg:{type:"none"}}).fixed.some(f=>f.k==="stG"),"ST 수식어: 매개변수 stMu·stG 추가, γ 0 이면 stG 기본 고정 (ST(μ) = γ 없음)");
  const tp=C.newComp("TP",0,{F:5}); ok(C.buildMap({comps:[tp],cons:{},bg:{type:"none"}}).fixed.some(f=>f.k==="ta"),"TP: α 는 기본 고정 (넓이와 함께 α/β 비만 정해짐)");
  const ns=C.newComp("GLS",0,{mod:"SB"}), nd=C.newComp("DS",0,{}); nd.mod="ST"; ok(C.compMod(ns)==="SB"&&C.compMod(nd)===null&&C.paramNames(nd).join()==="be,area,alpha,wL,wG","수식어는 넓이가 유한한 바탕(GL·SGL·LA·V·TLA·GLN·SGLN)에만 — DS 에는 무시");
}

console.log("== 3 합성 잡음 스펙트럼에서 회복 (Poisson, 3σ 적중률 ≥ 90 %)");
const R=C.rng(777);
const pois=l=>{ if(l>40) return Math.max(0,Math.round(l+Math.sqrt(l)*R.g())); let L=Math.exp(-l),k=0,p=1; do{k++;p*=R.u();}while(p>L); return k-1; };
const RUNS=QUICK?40:200;
const REC=[
  {t:"TLA(α 1.4, μ 2, wG 0.5)",c:{shape:"TLA",F:1.0,alpha:1.4,mu:2,wG:0.5},free:["be","area","F","mu","wG"],lo:280,hi:296,be:285},
  {t:"GL(30,·)T(2) wG 0.4",c:{shape:"GLN",F:1.1,m:0.3,wG:0.4,tail:"T",tk:2},free:["be","area","F","m","wG","tk"],lo:280,hi:294,be:285},
  {t:"SGL(30,·)PHI(1,2) wG 0.3",c:{shape:"SGLN",F:1.1,m:0.3,wG:0.3,tail:"PHI",tj:1,tk:2},free:["be","area","F","wG","tj"],lo:280,hi:294,be:285},
  {t:"ST(1.5)LA(1,1,·) wG 0.5",c:{shape:"LA",F:0.8,alpha:1,beta:1,wG:0.5,mod:"ST",stMu:1.5,stG:0},free:["be","area","F","wG","stMu"],lo:279,hi:295,be:285},
  {t:"ST(1,γ 0.3)GLS(m 0.3)",c:{shape:"GLS",F:1.0,m:0.3,mod:"ST",stMu:1,stG:0.3},free:["be","area","F","stMu","stG"],lo:279,hi:295,be:285},
  {t:"DS(0.15) wL 0.3 wG 0.5 (창 넓이)",c:{shape:"DS",alpha:0.15,wL:0.3,wG:0.5},free:["be","area","alpha","wL","wG"],lo:280,hi:294,be:285},
  {t:"DS(0.2,·,l 20) wL 0.3 wG 0.4",c:{shape:"DS",alpha:0.2,wL:0.3,wG:0.4,dsl:20},free:["be","area","alpha","wL","wG"],lo:280,hi:294,be:285},
  {t:"LA(1,1,·) + SB(0)LA(1,1,·) (계단)",c:{shape:"LA",F:0.8,alpha:1,beta:1,wG:0.5},c2:{shape:"LA",F:0.8,alpha:1,beta:1,wG:0.5,mod:"SB",area:3000},free:["be","area","F"],free2:["area"],lo:279,hi:295,be:285},
  {t:"GLS + TP(1,3) ⊗ L 0.5 (배경 성분)",c:{shape:"GLS",F:1.0,m:0.2},c2:{shape:"TP",F:14,ta:1,tb:3,wL:0.5,area:4000,be:287},free:["be","area","F"],free2:["area","tb","F"],lo:279,hi:295,be:285},
];
for(const T of REC){
  const E=Float64Array.from({length:Math.round((T.hi-T.lo)/0.05)+1},(_,i)=>T.lo+0.05*i), n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const truth=Object.assign({id:"a",be:T.be,area:20000,doublet:null},T.c), comps=[truth]; if(T.c2) comps.push(Object.assign({id:"b",be:T.c2.be||T.be,doublet:null},T.c2));
  const mu0=new Float64Array(n); comps.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) mu0[i]+=v[i]; }); for(let i=0;i<n;i++) mu0[i]+=50;
  const cons={a:{},b:{}}; C.paramNames(truth).forEach(k=>{ if(!T.free.includes(k)) cons.a[k]={mode:"fixed"}; }); if(T.c2) C.paramNames(comps[1]).forEach(k=>{ if(!T.free2.includes(k)) cons.b[k]={mode:"fixed"}; });
  const hits={}, cnt={}; let fails2=0; const t0=Date.now();
  for(let run=0;run<RUNS;run++){
    const y=Float64Array.from(mu0,pois);
    const st={lo:T.lo,hi:T.hi,bg:{type:"active",c0:50,k:0,s:0,_auto:false,cons:{c0:{mode:"free"},k:{mode:"fixed"},s:{mode:"fixed"}}},comps:C.clone(comps),cons:C.clone(cons)};
    st.comps.forEach(c=>{ (c.id==="a"?T.free:T.free2).forEach(k=>{ if(k==="area") c[k]*=1+0.1*(R.u()-0.5); else if(k==="be") c[k]+=0.05*(R.u()-0.5); else if(c[k]>0) c[k]*=1+0.06*(R.u()-0.5); }); });
    let P,r,cov; try{ P=C.makeProblem({be:E,y,c:1},st); r=C.fitOnce(P,{}); cov=C.covariance(P,r); }catch(e){ fails2++; continue; }
    r.map.free.forEach((f,j)=>{ if(f.ci<0) return; const key=(f.ci?"2:":"")+f.k, tv=C.getP(comps[f.ci],f.k), v=C.getP(P.comps[f.ci],f.k), s=cov.sig[j];
      cnt[key]=(cnt[key]||0)+1; if(s>0&&Math.abs(v-tv)<=3*s) hits[key]=(hits[key]||0)+1; });
  }
  const rate=Object.keys(cnt).map(k=>[k,(hits[k]||0)/RUNS]); const worst=Math.min(...rate.map(x=>x[1]));
  ok(worst>=0.9&&fails2===0,`${T.t}: ${RUNS} 회, 모든 자유 매개변수 |p̂ − p| ≤ 3σ ≥ 90 %`,rate.map(([k,v])=>`${k} ${(100*v).toFixed(0)}%`).join(" · ")+` · ${((Date.now()-t0)/RUNS).toFixed(0)} ms/피팅`);
}

console.log("== 4 속도 (10 성분 × 2000 채널)");
{ const E=Float64Array.from({length:2000},(_,i)=>270+0.01*i), ctx={emin:E[0],emax:E[1999]}; const kinds=[
    ["TLA",{F:1,alpha:1.3,mu:2,wG:0.5}],["GLN",{F:1,m:0.3,wG:0.4,tail:"T",tk:2}],["SGLN",{F:1,m:0.3,casaN:50,tail:"PHI",tj:1,tk:2}],["LA",{F:0.8,alpha:1,beta:1,casaN:143,mod:"ST",stMu:1}],
    ["GLS",{F:1,m:0.3,mod:"ST",stMu:1.5,stG:0.3}],["DS",{alpha:0.2,wL:0.3,wG:0.5}],["DS",{alpha:0.3,wL:0.2,casaN:600,dsl:40}],["LA",{F:0.8,alpha:1,beta:1,wG:0.5,mod:"SB"}],["TP",{F:15,ta:1,tb:3,wL:0.5}],["V",{wG:0.6,wL:0.3,mod:"ST",stMu:2,stG:0}]];
  kinds.forEach(([sh,o])=>{ const c=Object.assign({id:"x",shape:sh,be:280,area:1,doublet:null},o); const t0=process.hrtime.bigint(); let k=0; for(;k<50;k++){ c.be=280+1e-7*k; C.compEval(c,E,ctx); }   /* 위치를 바꿔 캐시를 피함 (실제 계산 시간) */ const ms=Number(process.hrtime.bigint()-t0)/1e6/k; console.log(`     ${sh}${o.mod?"+"+o.mod:""}${o.tail?"+"+o.tail:""}: ${ms.toFixed(2)} ms / 성분 평가`); });
  const comps=kinds.map(([sh,o],i)=>Object.assign({id:"c"+i,shape:sh,be:272+1.6*i,area:5000,doublet:null},o)); const mu=new Float64Array(2000); comps.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<2000;i++) mu[i]+=v[i]; });
  const y=Float64Array.from(mu,v=>Math.max(0,Math.round(v+100+Math.sqrt(v+100)*R.g())));
  const cons={}; comps.forEach(c=>{ cons[c.id]={}; C.paramNames(c).forEach(k=>{ if(!["be","area",c.shape==="DS"||c.shape==="V"?"wL":"F"].includes(k)) cons[c.id][k]={mode:"fixed"}; }); });
  const P=C.makeProblem({be:E,y,c:1},{lo:E[0],hi:E[1999],bg:{type:"active",c0:100,k:0,s:0,_auto:false,cons:{k:{mode:"fixed"},s:{mode:"fixed"}}},comps:C.clone(comps).map(c=>{ c.be+=0.02; return c; }),cons});
  const t0=Date.now(); const r=C.fitOnce(P,{}); const tf=Date.now()-t0; const t1=Date.now(); const mc=C.monteCarlo(P,Object.assign({},r),5,{seed:3}); const tm=(Date.now()-t1)/5;
  ok(tf<5000&&tm<3000&&r.chi2/2000<1.2,`10 성분(새 모양 모두) × 2000 채널 피팅 ${(tf/1000).toFixed(1)} s (≤ 5 s) (${r.iter} 반복) · Monte Carlo ${(tm/1000).toFixed(1)} s/회 (≤ 3 s; 자유 = BE·넓이·주 폭 — CasaXPS 처럼 모양 매개변수 고정)`,`χ²/N ${(r.chi2/2000).toFixed(2)}`);
  /* 4b 퇴화 매개변수까지 푼 경우 (SPEC §3.3b · SPEC_shapes2 §3): TP 의 wL 도 자유 → wL 이 0.0015·F 아래로 가면 모형이 wL 을 보지 못하고
     사다리꼴 가장자리가 채널보다 날카로워 χ² 가 계단처럼 끊긴다. 옛 LM 은 감쇠에 갇혀 기어감(QUICK 자료 450 + 452 반복 · 66 s).
     멈춤 규칙: ≤ 5 s 에 끝나고, 멈췄으면(stalled) wL 을 '잘 정해지지 않음'으로 표시. χ² 참고값 = 옛 LM (XPS_NOSTALL=1 로 잰 값) */
  { const cons2=C.clone(cons); cons2.c8.wL={mode:"free"};
    const P2=C.makeProblem({be:E,y,c:1},{lo:E[0],hi:E[1999],bg:{type:"active",c0:100,k:0,s:0,_auto:false,cons:{k:{mode:"fixed"},s:{mode:"fixed"}}},comps:C.clone(comps).map(c=>{ c.be+=0.02; return c; }),cons:cons2});
    const t2=Date.now(); const r2=C.fitOnce(P2,{}); const tf2=Date.now()-t2; const fl=((r2.stall||{}).flags||[]).map(q=>`${q.id}|${q.k}(${q.why.join("+")})`);
    const REF=QUICK?{chi2:2013.5711,it:"450+452",s:65.8}:{chi2:1940.0021,it:"15+7",s:3.1};
    const rel=(r2.chi2-REF.chi2)/REF.chi2, okStop=r2.stop==="stalled"?fl.some(t=>t.startsWith("c8|wL")):r2.stop==="converged";
    ok(okStop&&tf2<5000&&rel<2.5e-3,`TP wL 까지 자유(퇴화): ${(tf2/1000).toFixed(1)} s (≤ 5 s) · ${r2.iter1}+${r2.iter} 반복 · 정지 ${r2.stop}${fl.length?" · 표시 "+fl.join(", "):""}`,
      `χ² ${r2.chi2.toFixed(2)} vs 옛 LM ${REF.chi2} (${REF.it} 반복, ${REF.s} s): ${rel>=0?"+":""}${rel.toExponential(1)} (허용 < 2.5e-3)`); } }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
