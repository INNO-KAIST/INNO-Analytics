/* CasaXPS/Avantage 대응 기능 시험 (SPEC §2 spline · §3.3b Simplex · §4 at% MC · §12 측정 선모양·LCF · §15 여러 스펙트럼)
   node test/test_xps_casa.js   (modules/xps 에서, 약 3–5 분)
   1) 엔진 ↔ 독립 Python 기준(test/ref/ref_casa.json: PchipInterpolator·quad·nnls·least_squares·CubicSpline·Nelder–Mead·numpy MC)
   2) 참값을 아는 합성 데이터(Poisson 잡음)로 회복·불확도 적중(95 % 구간 ≈ 95 %)                                  */
const path=require("path"), fs=require("fs");
const D=p=>path.join(__dirname,p);
const C=require("../src/xps_core.js");
const REF=JSON.parse(fs.readFileSync(D("ref/ref_casa.json"),"utf8"));
const REF1=JSON.parse(fs.readFileSync(D("ref/ref.json"),"utf8"));
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const F64=a=>Float64Array.from(a);
const QUICK=!!process.env.QUICK;
const R=C.rng(4242);
const pois=l=>{ if(l>40) return Math.max(0,Math.round(l+Math.sqrt(l)*R.g())); let L=Math.exp(-l),k=0,p=1; do{k++;p*=R.u();}while(p>L); return k-1; };
const grid=(a,st,n)=>Float64Array.from({length:n},(_,i)=>a+i*st);
const shapeOn=(E,shape,be,o)=>C.compEval(C.newComp(shape,be,Object.assign({area:1},o)),E,{emin:E[0],emax:E[E.length-1]});
const fixedMeas=(id,tpl,be,o)=>Object.assign(C.newComp("MEAS",be,Object.assign({tpl,id},o||{})));

console.log("== 1a PCHIP (scipy PchipInterpolator) · 보간 적분");
REF.pchip.forEach((t,ci)=>{ const x=F64(t.x), v=F64(t.v), m=C.pchipSlopes(x,v); const sc=Math.max(...t.v.map(Math.abs));
  let e=0; t.q.forEach((q,i)=>{ e=Math.max(e,Math.abs(C.interpAt(x,v,m,q)-t.val[i])/sc); });
  let es=0; t.slopes.forEach((s,i)=>{ es=Math.max(es,Math.abs(m[i]-s)/Math.max(1,Math.abs(s))); });
  const ei=Math.abs(C.interpIntegral(t.x,t.v,"pchip")-t.int)/Math.max(1e-12,Math.abs(t.int)), el=Math.abs(C.interpIntegral(t.x,t.v,"linear")-t.int_lin)/Math.max(1e-12,Math.abs(t.int_lin));
  ok(e<1e-12&&es<1e-12&&ei<1e-12&&el<1e-12,`PCHIP 경우 ${ci} (${t.x.length} 마디): 값·기울기·적분 = scipy`,`값 ${e.toExponential(1)} · 기울기 ${es.toExponential(1)} · ∫ ${ei.toExponential(1)} · ∫선형 ${el.toExponential(1)}`); });

console.log("== 1b template (배경 빼기·넓이 1) · Gaussian 넓힘 vs quad");
REF.tpl.forEach(t=>{ const tp=C.makeTemplate(F64(t.E),F64(t.y),{bg:{type:"linear",n:5},shift:t.shift,eref:t.eref,interp:t.kind});
  const sc=Math.max(...t.v); let ex=0, ev=0; t.x.forEach((x,i)=>{ ex=Math.max(ex,Math.abs(tp.x[i]-x)); ev=Math.max(ev,Math.abs(tp.v[i]-t.v[i])/sc); });
  ok(ex<1e-12&&ev<1e-12&&Math.abs(tp.area0-t.area0)/t.area0<1e-12,`${t.kind}: template x·v·넓이 = 기준`,`x ${ex.toExponential(1)} · v ${ev.toExponential(1)}`);
  let ep=0; const pl=C.tplShape(tp,F64(t.X),0); t.plain.forEach((r,i)=>{ ep=Math.max(ep,Math.abs(pl[i]-r)/sc); });
  ok(ep<1e-12,`${t.kind}: 넓힘 없음 = 보간값(마디 밖 0)`,ep.toExponential(1));
  Object.entries(t.broad).forEach(([w,ref])=>{ const b=C.tplShape(tp,F64(t.X),+w); const s2=Math.max(...ref); let e=0; ref.forEach((r,i)=>{ e=Math.max(e,Math.abs(b[i]-r)/s2); });
    ok(e<2e-4,`${t.kind}: ⊗ G(${w} eV) vs quad (최대오차/최댓값 < 2e-4)`,e.toExponential(2)); });
  // 넓힘 뒤에도 넓이 1
  const X=grid(tp.x[0]-5,0.002,Math.round((tp.x[tp.x.length-1]-tp.x[0]+10)/0.002)); const b=C.tplShape(tp,X,0.8); let a=0; for(let i=1;i<X.length;i++) a+=0.5*(b[i]+b[i-1])*0.002;
  ok(Math.abs(a-1)<2e-4,`${t.kind}: ⊗ G 뒤 넓이 1`,a.toFixed(6)); });

console.log("== 1c 보간 오차(해석 모양을 간격 s 로 표본) — SPEC §12.2 표");
{ let worst=0; REF.interp_err.forEach(r=>{ const xs=[]; for(let x=-8;x<=8+r.step/2;x+=r.step) xs.push(+x.toFixed(10));
    const f=r.shape.startsWith("GLS")?(t=>(0.7*C.gauss(t,1)+0.3*C.lorentz(t,1))):(t=>C.voigt(t,0.8,0.3));
    const X=F64(xs), V=Float64Array.from(xs,f), m=C.pchipSlopes(X,V); let mx=0, e=0; for(let i=0;i<=30000;i++){ const t=-7.5+15*i/30000, tv=f(t); mx=Math.max(mx,tv); e=Math.max(e,Math.abs(C.interpAt(X,V,m,t)-tv)); }
    e/=mx; worst=Math.max(worst,Math.abs(e-r.pchip)/r.pchip); });
  ok(worst<0.02,"PCHIP 보간 오차 표 재현(Python 과 상대 2 % 안)",`최대 상대차 ${worst.toExponential(1)}; 간격 0.05 eV·FWHM 1 eV → ${REF.interp_err.find(r=>r.step===0.05).pchip.toExponential(1)}`); }

console.log("== 1d LCF — 선형(이동·넓힘 없음), 정적 직선 배경, 2회 IRLS vs scipy nnls");
let LCFP;
{ const L=REF.lcf_lin, E=F64(L.E), y=F64(L.y);
  const T=L.tpl.map((t,i)=>({kind:"tpl",id:"r"+i,name:"R"+i,interp:"pchip",eref:L.eref,x:t.x,v:t.v,area0:1,raw:null}));
  const prob={E,y,c:L.c,bg:{type:"linear",n:5},comps:T.map((t,i)=>fixedMeas("c"+i,t,L.eref)),cons:{c0:{be:{mode:"fixed"},wG:{mode:"fixed"}},c1:{be:{mode:"fixed"},wG:{mode:"fixed"}}}};
  const r=C.fitOnce(prob,{}); const cov=C.covariance(prob,r); const sm=C.compSummary(prob,r,cov);
  let ea=0; prob.comps.forEach((c,i)=>{ ea=Math.max(ea,Math.abs(c.area-L.areas[i])/L.areas[i]); });
  ok(ea<1e-7&&Math.abs(r.chi2-L.chi2)/L.chi2<1e-8,"넓이 a_i·χ² = nnls (상대 < 1e-7)",`넓이 ${ea.toExponential(1)}, χ² ${r.chi2.toFixed(4)} vs ${L.chi2.toFixed(4)}`);
  let es=0; prob.comps.forEach((c,i)=>{ es=Math.max(es,Math.abs(sm[i].sig.area-L.sig[i])/L.sig[i]); }); ok(es<0.01,"넓이 σ = (JᵀV⁻¹J)⁻¹φ (상대 < 1 %)",(100*es).toFixed(3)+" %");
  let ef=0; sm.forEach((q,i)=>{ ef=Math.max(ef,Math.abs(q.areaPct-L.frac[i])); }); ok(ef<1e-6,"분율 a_i/Σa (%)",ef.toExponential(1));
  const dg=C.lcfDiag(prob,F64(L.v2)); const ev=Math.max(...dg.vif.map((v,i)=>Math.abs(v-L.vif[i])/L.vif[i]));
  ok(Math.abs(dg.cond-L.cond)/L.cond<1e-8&&ev<1e-8,"공선성: 조건수 = numpy.linalg.cond(열 정규화), VIF = diag(corrcoef⁻¹)",`κ ${dg.cond.toFixed(4)} vs ${L.cond.toFixed(4)}, VIF ${dg.vif.map(v=>v.toFixed(3))} vs ${L.vif.map(v=>v.toFixed(3))}`);
  LCFP=prob; }

console.log("== 1e LCF — 이동 δ·넓힘 w 자유 vs scipy least_squares(TRF, 미세 격자 합성곱)");
{ const L=REF.lcf_shift, B=REF.lcf_lin, E=F64(L.E), y=F64(L.y);
  const T=B.tpl.map((t,i)=>({kind:"tpl",id:"r"+i,name:"R"+i,interp:"pchip",eref:B.eref,x:t.x,v:t.v,area0:1,raw:null}));
  const prob={E,y,c:1,bg:{type:"linear",n:5},comps:[fixedMeas("a",T[0],L.start[1],{area:L.start[0],wG:L.start[2]}),fixedMeas("b",T[1],L.start[4],{area:L.start[3]})],cons:{b:{wG:{mode:"fixed"}}}};
  const r=C.fitOnce(prob,{}); const cov=C.covariance(prob,r);
  const get=[["a","area"],["a","be"],["a","wG"],["b","area"],["b","be"]]; let w=0, ws=0;
  get.forEach(([id,k],j)=>{ const c=prob.comps.find(q=>q.id===id); w=Math.max(w,Math.abs(c[k]-L.x[j])/L.sig[j]); const jj=r.map.free.findIndex(f=>f.id===id&&f.k===k); ws=Math.max(ws,Math.abs(cov.sig[jj]-L.sig[j])/L.sig[j]); });
  ok(w<0.05&&Math.abs(r.chi2-L.chi2)/L.chi2<1e-4,"매개변수 |Δ|/σ < 0.05, χ² 상대차 < 1e-4",`최대 ${w.toFixed(4)}σ, χ² ${r.chi2.toFixed(3)} vs ${L.chi2.toFixed(3)}; δ_a ${(prob.comps[0].be-B.eref).toFixed(3)} eV, w_a ${prob.comps[0].wG.toFixed(3)} eV`);
  ok(ws<0.03,"공분산 σ 일치 (상대 < 3 %)",(100*ws).toFixed(2)+" %"); }

console.log("== 1f Spline 배경 vs scipy CubicSpline(bc_type='natural') (+ 끝 매듭 밖 직선 연장)");
REF.spline.forEach((s,i)=>{ const E=F64(s.E), y=F64(s.y); const B=C.bgStatic(E,y,Object.assign({type:"spline"},s.bg)); const rg=Math.max(...s.y)-Math.min(...s.y);
  let e=0; s.B.forEach((b,j)=>{ e=Math.max(e,Math.abs(B[j]-b)/rg); }); let ek=0; B.knots.forEach((k,j)=>{ ek=Math.max(ek,Math.abs(k.e-s.knots[j][0])+Math.abs(k.y-s.knots[j][1])/rg); });
  ok(e<1e-10&&ek<1e-10&&B.knots.length===s.knots.length,`경우 ${i}: 매듭 ${s.knots.length} 개${s.bg.ends===false?" (끝점 매듭 없음 → 연장 구간 포함)":""}`,`최대오차/세기범위 ${e.toExponential(1)}`); });

console.log("== 1g Simplex (Nelder–Mead) vs scipy Nelder–Mead · least_squares, 같은 분산 v");
{ const Sx=REF.simplex, E=F64(Sx.E), y=F64(Sx.y), v=F64(Sx.v); const s0=Sx.start;
  const mk=()=>({E,y,c:1,bg:{type:"linear",n:5},comps:[C.newComp("GLS",s0[0],{F:s0[1],m:s0[2],area:s0[3],id:"a"}),C.newComp("GLS",s0[4],{F:s0[5],m:s0[6],area:s0[7],id:"b"})],cons:{}});
  const vec=P=>[P.comps[0].be,P.comps[0].F,P.comps[0].m,P.comps[0].area,P.comps[1].be,P.comps[1].F,P.comps[1].m,P.comps[1].area];
  const Pn=mk(), mn=C.buildMap(Pn); const t0=Date.now(); const rn=C.nmFit(Pn,mn,v,{}); const tn=Date.now()-t0;
  const Pl=mk(), ml=C.buildMap(Pl); const rl=C.lmFit(Pl,ml,v,{});
  const dz=(a,b)=>Math.max(...a.map((x,i)=>Math.abs(x-b[i])/Sx.sig[i]));
  ok(dz(vec(Pn),Sx.ls_x)<0.01&&Math.abs(rn.chi2-Sx.ls_chi2)/Sx.ls_chi2<1e-9&&rn.conv,"JS Simplex = 참 최솟값(least_squares): |Δ|/σ < 0.01, χ² 상대 < 1e-9",`${dz(vec(Pn),Sx.ls_x).toExponential(1)}σ, χ² ${rn.chi2.toFixed(6)} vs ${Sx.ls_chi2.toFixed(6)}, ${rn.iter} 평가 · ${tn} ms`);
  ok(dz(vec(Pn),Sx.nm_x)<0.01,"JS Simplex = scipy Nelder–Mead(경계, 재시작)",`${dz(vec(Pn),Sx.nm_x).toExponential(1)}σ (scipy NM χ² ${Sx.nm_chi2.toFixed(6)})`);
  ok(dz(vec(Pn),vec(Pl))<0.01&&Math.abs(rn.chi2-rl.chi2)/rl.chi2<1e-9,"JS Simplex = JS LM",`${dz(vec(Pn),vec(Pl)).toExponential(1)}σ`); }

console.log("== 1h 정량 at% Monte Carlo vs numpy (40만 회)");
{ const Q=REF.quantmc; const rows=Q.rows.map((r,i)=>({label:"e"+i,el:"e"+i,area:r.area,sArea:r.sArea,S:r.S,BE:500,hv:1486.6,tf:{model:"none"}}));
  const mc=C.quantMC(rows,{imfp:{model:"none"}},100000,11); const qd=C.quantify(rows,{imfp:{model:"none"}});
  let es=0, ep=0, em=0, ed=0; mc.forEach((m,i)=>{ es=Math.max(es,Math.abs(m.sd-Q.sd[i])/Q.sd[i]); ep=Math.max(ep,Math.abs(m.lo-Q.lo[i])/Q.sd[i],Math.abs(m.hi-Q.hi[i])/Q.sd[i]); em=Math.max(em,Math.abs(m.mean-Q.mean[i])/Q.sd[i]); ed=Math.max(ed,Math.abs(qd[i].sAt-Q.delta[i])/Q.delta[i]); });
  ok(es<0.015&&ep<0.04&&em<0.02,"MC sd·2.5/97.5 % = numpy (sd 상대 < 1.5 %, 백분위·평균 < 0.04 sd)",`sd ${(100*es).toFixed(2)} %, 백분위 ${ep.toFixed(3)} sd, 평균 ${em.toFixed(3)} sd`);
  ok(ed<1e-10,"델타 σ(at%) = numpy 델타",ed.toExponential(1));
  const rel=mc.map((m,i)=>(m.sd-qd[i].sAt)/qd[i].sAt); console.log("     MC sd / 델타 σ − 1:",rel.map(x=>(100*x).toFixed(1)+" %").join(", "),"(넓이 상대오차가 큰 행일수록 델타가 벗어남 — CasaXPS 식 MC 를 쓰는 까닭)"); }

console.log("== 2a 합성(참값): LCF 분율 — 잡음 있는 기준 스펙트럼, 95 % MC 구간 적중");
{ const n=241, E=grid(282,0.05,n); const K=QUICK?16:60, Mmc=QUICK?60:150; const f1=0.6, Atot=30000;
  // 참 모양은 창 안에서 꼬리가 끝나는 것(Gaussian)으로 — Lorentzian 꼬리는 창 끝 배경 평균에 섞여 기준·대상에서 다르게 잘려 분율을 치우치게 한다(불확도 시험과 무관한 편향)
  const s1=shapeOn(E,"GLS",285.0,{F:1.1,m:0}), s2=shapeOn(E,"GLS",285.9,{F:1.5,m:0});
  let inRef=0, inNo=0, inCov=0, t0=Date.now(); const fr=[];
  for(let k=0;k<K;k++){
    const mkRef=s=>Float64Array.from(s,(v,i)=>pois(300+8000*v));   // 기준: 봉우리 높이 ≈ 2.6k count — 기준 잡음이 무시 못 할 정도
    const T1=C.makeTemplate(E,mkRef(s1),{bg:{type:"linear",n:5},name:"A",id:"A"}), T2=C.makeTemplate(E,mkRef(s2),{bg:{type:"linear",n:5},name:"B",id:"B"});
    const y=Float64Array.from(E,(e,i)=>pois(800+Atot*(f1*s1[i]+(1-f1)*s2[i])));
    const prob={E,y,c:1,bg:{type:"linear",n:5},comps:[fixedMeas("a",T1,T1.eref),fixedMeas("b",T2,T2.eref)],cons:{a:{be:{mode:"fixed"},wG:{mode:"fixed"}},b:{be:{mode:"fixed"},wG:{mode:"fixed"}}}};
    const r=C.fitOnce(prob,{}); const sm=C.compSummary(prob,r,C.covariance(prob,r)); fr.push(sm[0].areaPct);
    if(Math.abs(sm[0].areaPct-100*f1)<1.96*sm[0].sig.areaPct) inCov++;
    const res={mu:r.mu,v:r.v,chi2:r.chi2,map:r.map};
    const a=C.mcSummarize(C.monteCarlo(prob,res,Mmc,{seed:100+k}).samples)[0].areaPct, b=C.mcSummarize(C.monteCarlo(prob,res,Mmc,{seed:100+k,refNoise:false}).samples)[0].areaPct;
    if(a.lo<=100*f1&&100*f1<=a.hi) inRef++; if(b.lo<=100*f1&&100*f1<=b.hi) inNo++;
  }
  const m=fr.reduce((s,x)=>s+x,0)/K;
  ok(inRef/K>=0.85,`MC(기준 잡음 포함) 95 % 구간 적중 ${inRef}/${K}`,`목표 ≈ 95 %; 기준 잡음 뺀 MC ${inNo}/${K}, 공분산 ±1.96σ ${inCov}/${K}; 분율 평균 ${m.toFixed(2)} (참 60) · ${((Date.now()-t0)/1000).toFixed(0)} s`);
  ok(inNo<inRef,"기준 잡음을 빼면 구간이 좁아 적중이 떨어진다(기준 잡음을 넣어야 하는 근거)",`${inNo} < ${inRef}`);
  ok(Math.abs(m-60)<3*Math.sqrt(fr.reduce((s,x)=>s+(x-m)*(x-m),0)/(K-1))/Math.sqrt(K)+0.2,"분율 치우침 없음 (평균이 참값 ± 3 SE)",m.toFixed(3)); }

console.log("== 2b 합성(참값): 측정 선모양 성분 — 비대칭 모양(LA) 두 화학 상태, 위치·넓이 비·추가 넓힘 회복");
{ const n=401, E=grid(280,0.05,n); const K=QUICK?10:30; const sh={F:0.9,alpha:1.1,beta:2.2,wG:0.5};
  const meas=shapeOn(E,"LA",284.0,sh), wB=0.6, shB=Object.assign({},sh,{wG:Math.hypot(sh.wG,wB)}), p1=meas, p2=shapeOn(E,"LA",285.5,shB);   // 참: 두 번째 상태 = 같은 모양 ⊗ G(0.6) (Gaussian 합성곱의 합성)
  let inD=0, inR=0, inW=0, t0=Date.now(); const dW=[];
  for(let k=0;k<K;k++){
    const T=C.makeTemplate(E,Float64Array.from(meas,v=>pois(200+150000*v)),{bg:{type:"linear",n:5},id:"T",name:"metal"});
    const y=new Float64Array(n); const A1=12000, A2=6000;
    for(let i=0;i<n;i++) y[i]=pois(1500+20*(E[i]-280)+A1*p1[i]+A2*p2[i]);
    const prob={E,y,c:1,bg:{type:"linear",n:5},comps:[fixedMeas("a",T,T.eref+0.1),fixedMeas("b",T,T.eref+1.3,{wG:0.3})],cons:{a:{wG:{mode:"fixed"}}}};
    const r=C.fitOnce(prob,{}); const cov=C.covariance(prob,r); const map=r.map; const sg=(id,kk)=>cov.sig[map.free.findIndex(f=>f.id===id&&f.k===kk)];
    const a=prob.comps[0], b=prob.comps[1];
    // 위치 차·넓이 비의 σ: 델타(공분산)
    const ia=map.free.findIndex(f=>f.id==="a"&&f.k==="be"), ib=map.free.findIndex(f=>f.id==="b"&&f.k==="be"), np=cov.np; const sD=Math.sqrt(cov.C[ia*np+ia]+cov.C[ib*np+ib]-2*cov.C[ia*np+ib]);
    const ja=map.free.findIndex(f=>f.id==="a"&&f.k==="area"), jb=map.free.findIndex(f=>f.id==="b"&&f.k==="area"); const rr=b.area/a.area; const g1=-rr/a.area, g2=1/a.area; const sR=Math.sqrt(g1*g1*cov.C[ja*np+ja]+g2*g2*cov.C[jb*np+jb]+2*g1*g2*cov.C[ja*np+jb]);
    if(Math.abs(b.be-a.be-1.5)<1.96*sD) inD++; if(Math.abs(rr-A2/A1)<1.96*sR) inR++; if(Math.abs(b.wG-wB)<1.96*sg("b","wG")) inW++; dW.push(b.wG);
  }
  ok(inD/K>=0.8&&inR/K>=0.8&&inW/K>=0.8,`위치 차 1.5 eV·넓이 비 0.5·추가 넓힘 0.6 eV — 95 % 구간 적중`,`ΔBE ${inD}/${K}, 비 ${inR}/${K}, w ${inW}/${K}; w 평균 ${(dW.reduce((s,x)=>s+x,0)/K).toFixed(3)} eV · ${((Date.now()-t0)/1000).toFixed(0)} s`); }

console.log("== 2c 합성(참값): spline 배경 — 휜 배경 위 봉우리 넓이, 매듭 흔들기(모형 민감도)");
{ const n=281, E=grid(524,0.05,n); const K=QUICK?6:16; let good=0, cover=0; const errs=[], zs=[]; let t0=Date.now();
  // 참 배경: 매듭(봉우리 밖)으로 정해지는 매끄러운 곡선. (봉우리 밑의 계단처럼 매듭 밖에서 보이지 않는 모양은 어떤 보간으로도 정할 수 없다 — SPEC §2 spline 주의)
  const truthB=E.map(e=>2500+500*Math.sin((e-524)/4.5)+15*(e-524)); const A=15000; const pk=shapeOn(E,"GLS",531,{F:1.4,m:0});
  for(let k=0;k<K;k++){
    const y=Float64Array.from(E,(e,i)=>pois(truthB[i]+A*pk[i]));
    const model={lo:E[0],hi:E[n-1],bg:{type:"spline",n:5,knots:[{e:527.0},{e:535.0}]},comps:[C.newComp("GLS",531.1,{F:1.3,m:0.3,area:10000,id:"a"})],cons:{}};
    const spec={be:E,y,c:1};
    const vars=C.sensVariants(model,spec,{mode:"oat",bgs:[],ns:[],dw:0.3,dk:0.5}).filter(v=>!v.f.lo&&!v.f.hi&&(v.f.n===5));
    const res=C.sensJob(spec,model,vars.map((v,i)=>Object.assign({},v,{idx:i}))); const sum=C.sensSummary(vars,res);
    const P=C.makeProblem(spec,model); const r=C.fitOnce(P,{}); const sm=C.compSummary(P,r,C.covariance(P,r)); const st=C.fitStats(P,r);
    const ar=sum.comps[0].area, sS=sm[0].sig.area*st.ar1Inflate; errs.push((sm[0].area-A)/A);
    zs.push((sm[0].area-A)/sS); if(Math.abs(sm[0].area-A)<4*sS) good++; if(ar.min-2*sS<=A&&A<=ar.max+2*sS) cover++;
  }
  const me=errs.reduce((s,x)=>s+x,0)/K;
  // 공분산 σ 에는 매듭 높이(±kh 점 평균)의 잡음이 들어 있지 않다(정적 배경 모두 같음) → |z| 가 크다. 그래서 보고는 '가능 범위'(매듭 흔들기 ± 2σ)로 — SPEC §2 spline 주의
  const sdE=Math.sqrt(errs.reduce((s,x)=>s+(x-me)*(x-me),0)/(K-1));
  ok(Math.abs(me)<3*sdE/Math.sqrt(K)+0.003,`넓이 치우침 없음 (평균 오차 < 3 SE) — 공분산 σ 는 매듭 높이 잡음을 빼고 있음(|z| 참고, ${K} 회)`,`평균 ${(100*me).toFixed(2)} %, 오차 ${errs.map(e=>(100*e).toFixed(1)).join(" ")} %, |z| 최대 ${Math.max(...zs.map(Math.abs)).toFixed(1)}`);
  ok(cover/K>=0.8,`가능 범위(매듭 ±0.5 eV 변형의 범위 ± 2σ통계)가 참 넓이를 포함`,`${cover}/${K} · ${((Date.now()-t0)/1000).toFixed(0)} s`); }

console.log("== 2d Simplex vs LM — 기준 구현 ① 의 피팅 문제들(같은 최솟값)");
for(const f of REF1.fits.filter(f=>QUICK?/GLS3_shirley|DS_linear/.test(f.name):true)){
  const P=f.prob; const mk=o=>Object.assign({E:F64(P.E),y:F64(P.y),c:P.c,bg:C.clone(P.bg),comps:C.clone(P.comps),cons:C.clone(P.cons)},o||{});
  const A=mk(), B=mk({optimizer:"simplex"}); const ra=C.fitOnce(A,{}); const t0=Date.now(); const rb=C.fitOnce(B,{}); const ms=Date.now()-t0; const cov=C.covariance(A,ra);
  let w=0, wk=""; ra.map.free.forEach((fr,j)=>{ const a=fr.ci<0?A.bg[fr.k]:C.getP(A.comps[fr.ci],fr.k), b=fr.ci<0?B.bg[fr.k]:C.getP(B.comps[fr.ci],fr.k); const s=cov.sig[j]; if(!(s>0)) return; const d=Math.abs(a-b)/s; if(d>w){ w=d; wk=fr.id+"."+fr.k; } });
  ok(w<0.02&&Math.abs(ra.chi2-rb.chi2)/ra.chi2<1e-6,`${f.name}: Simplex = LM (|Δ|/σ < 0.02, χ² 상대 < 1e-6)`,`${w.toExponential(1)}σ (${wk}), χ² ${rb.chi2.toFixed(4)} vs ${ra.chi2.toFixed(4)}, ${rb.iter} 평가 · ${ms} ms`);
}

console.log("== 2e 합성(참값): 여러 스펙트럼에 같은 모형(깊이 분포·시료마다 다른 대전 이동) → 일괄 피팅");
{ const n=241, Ec=grid(96,0.05,n); const levels=[0.9,0.75,0.55,0.35,0.2,0.1]; const q=[0,0.6,-0.4,1.2,0.25,-0.8];   // 금속 분율 · 대전 이동(측정 = 참 + q)
  const truth=levels.map((fm,k)=>({fm,q:q[k]})); let okBE=0, okF=0, tot=0, t0=Date.now();
  // 모형: 0 번 스펙트럼에서(보정 BE 기준) — Si 금속 99.4·산화물 103.3 GLS
  const model0={lo:96,hi:107.9,bg:{type:"active",n:5},comps:[C.newComp("GLS",99.4,{F:1.0,m:0.2,id:"m"}),C.newComp("GLS",103.3,{F:1.6,m:0.2,id:"o"})],cons:{o:{m:{mode:"link",ref:"m",k:1,b:0}}}};
  truth.forEach((t,k)=>{ const E=Float64Array.from(Ec,e=>e+t.q); const pm=shapeOn(Ec,"GLS",99.4,{F:1.0,m:0.2}), po=shapeOn(Ec,"GLS",103.3,{F:1.6,m:0.2});
    const P=Float64Array.from(Ec,(e,i)=>20000*(t.fm*pm[i]+(1-t.fm)*po[i])), Q=C.cumtrapz(Ec,P); const y=Float64Array.from(P,(v,i)=>pois(1000+0.02*Q[i]+v));   // 참 배경 = 피팅 안 Shirley
    const md=C.propagateModel(model0,0,-t.q,[E[0],E[n-1]]); const out=C.batchJob({be:E,y,c:1},md,{});
    const m=out.summ[0], o=out.summ[1];
    tot+=1; if(Math.abs(m.be-t.q-99.4)<1.96*m.sig.be) okBE++; if(Math.abs(o.be-t.q-103.3)<1.96*o.sig.be) okBE++;   // 보정 BE = 측정 BE − q
    if(Math.abs(m.areaPct-100*t.fm)<1.96*m.sig.areaPct) okF++; });
  ok(okBE>=0.75*2*truth.length&&okF>=0.75*truth.length,`보정 BE·금속 분율 — 95 % 구간 적중(성분별)`,`BE ${okBE}/${2*truth.length}, 분율 ${okF}/${truth.length} · ${((Date.now()-t0)/1000).toFixed(1)} s`);
  const md=C.propagateModel({lo:96,hi:108,bg:{type:"spline",n:5,knots:[{e:100}]},comps:model0.comps,cons:model0.cons,optimizer:"simplex"},0.3,-0.2,[95,107]);
  ok(Math.abs(md.comps[0].be-99.9)<1e-12&&Math.abs(md.bg.knots[0].e-100.5)<1e-12&&md.hi===107&&md.optimizer==="simplex"&&model0.comps[0].be===99.4,"모형 옮기기: BE·매듭 +(Δ_src − Δ_dst), 창은 스펙트럼 끝에서 자름, 원본 그대로",`${md.comps[0].be} · ${md.bg.knots[0].e} · ${md.lo}–${md.hi}`); }

console.log("== 2g 합성(참값): Avantage 꼬리 성분(GL·T) — ET·F·넓이 회복, 95 % 구간");
{ const n=341, E=grid(280,0.05,n), ctx={emin:E[0],emax:E[n-1]}; const K=QUICK?8:30; const tr={F:1.1,m:0.2,ET:0.06,TM:0.2,CT:0.08,ch:0.05}; let inE=0,inA=0,inF=0;
  const P0=C.compEval(C.newComp("AGT",284,Object.assign({area:20000},tr)),E,ctx);
  for(let k=0;k<K;k++){ const y=Float64Array.from(P0,v=>pois(500+v));
    const pb={E,y,c:1,bg:{type:"active",n:5,k:0,s:0,_auto:false,cons:{k:{mode:"fixed"},s:{mode:"fixed"}}},comps:[C.newComp("AGT",284.05,{F:1.0,m:0.25,ET:0.04,TM:0.2,CT:0.08,ch:0.05,id:"a",area:15000})],cons:{a:{TM:{mode:"fixed"},CT:{mode:"fixed"}}}};
    const r=C.fitOnce(pb,{}); const cov=C.covariance(pb,r); const sg=kk=>cov.sig[r.map.free.findIndex(f=>f.k===kk)]; const c=pb.comps[0];
    if(Math.abs(c.ET-tr.ET)<1.96*sg("ET")) inE++; if(Math.abs(c.area-20000)<1.96*sg("area")) inA++; if(Math.abs(c.F-tr.F)<1.96*sg("F")) inF++; }
  ok(inE/K>=0.8&&inA/K>=0.8&&inF/K>=0.8,"ET·넓이(창 안)·F 95 % 구간 적중 (배경 = 상수 c₀ 자유 — 창 끝 직선은 상수 꼬리를 배경으로 먹는다)",`ET ${inE}/${K}, 넓이 ${inA}/${K}, F ${inF}/${K}`); }

console.log("== 2f 기본 경로 불변(새 기능을 쓰지 않으면 결과가 같다)");
{ const n=201, E=grid(80,0.05,n); const y=Float64Array.from(shapeOn(E,"GLS",85,{F:1.2}),v=>pois(300+4000*v));
  const mk=()=>({E,y,c:1,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",85.1,{F:1.0,id:"a"})],cons:{}});
  const A=mk(), B=mk(); const ra=C.fitOnce(A,{}), rb=C.fitOnce(B,{optimizer:"lm"});
  const ma=C.monteCarlo(A,{mu:ra.mu,v:ra.v,chi2:ra.chi2,map:ra.map},5,{seed:3}), mb=C.monteCarlo(A,{mu:ra.mu,v:ra.v,chi2:ra.chi2,map:ra.map},5,{seed:3,refNoise:true});
  ok(ra.chi2===rb.chi2&&JSON.stringify(ma.samples)===JSON.stringify(mb.samples)&&C.finalize(A,ra).lcf===null,"LM 기본·MC(측정 선모양 없음) 결과 동일, lcf 진단 없음");
  const sv=C.sensVariants({lo:80,hi:90,bg:{type:"smart",n:5}},{be:E},{}); ok(sv.every(v=>v.f.kn==null),"모형 민감도: spline 이 아니면 매듭 변형 없음",sv.length+" 변형"); }

console.log("== 3 CasaXPS 2.3.27 실측 대조 (test/ref/casa_readings.json · test/data/casaval) — CasaXPS 규약(bg.conv = \"casa\")");
(async()=>{ const IO=require("../src/xps_io.js"), CR=JSON.parse(fs.readFileSync(D("ref/casa_readings.json"),"utf8"));
  const rd=async f=>(await IO.readFile(f+".vms",new Uint8Array(fs.readFileSync(D("data/casaval/"+f+".vms"))))).spectra[0];
  const {i0,i1}=CR.region; const rnd=(v,c)=>Math.abs(v-c)<=0.05+5e-6*Math.abs(c);   // 화면 값은 소수 한 자리
  for(const [f,t,w,A,mx,mn] of CR.regions){ if(t.startsWith("tougaard")) continue; const s=await rd(f); const R=C.sliceRegion(s.be,s.y,s.be[i0],s.be[i1]);
    const E=F64(R.E), y=F64(R.y); const B=C.bgStatic(E,y,{type:t,conv:"casa",aw:w,_ends:C.casaEnds(s.be,s.y,R.i0,R.i1,w)}); const st=C.casaRegionStats(E,y,B); const dA=Math.abs(st.area-A)/A;
    const info=`Area ${st.area.toFixed(1)} vs ${A} (${dA.toExponential(1)}), Max ${st.maxH.toFixed(1)} vs ${mx}, Min ${st.minH.toFixed(1)} vs ${mn}`;
    if(f==="Syn1"&&t==="shirley"){ console.log(`     알려진 차이(판정 안 함): ${f} ${t} — 기운 창이라 Shirley 고정점이 없음, CasaXPS 의 멈춤 규칙 미확인: ${info}`); continue; }
    if(f==="C1s"&&t==="shirley"){ ok(dA<1.5e-3&&Math.abs(st.maxH-mx)<2&&Math.abs(st.minH-mn)<2,`${f} ${t} Av.W ${w}: 넓이 0.15 %·높이 2 CPS 안 (실측 자료, 남은 차이 미확인)`,info); continue; }
    ok(rnd(st.maxH,mx)&&rnd(st.minH,mn)&&dA<(t==="shirley"?5e-6:1e-4),`${f} ${t} Av.W ${w}: Max·Min Height = CasaXPS, 넓이 상대 < ${t==="shirley"?"5e-6":"1e-4"}`,info); }
  // 같은 데이터, 기본(INNO) 규약은 그대로 — 차이를 기록
  { const s=await rd("Syn2"); const R=C.sliceRegion(s.be,s.y,s.be[i0],s.be[i1]); const E=F64(R.E), y=F64(R.y); const B=C.bgStatic(E,y,{type:"shirley",n:3}); const Q=C.cumtrapz(E,Float64Array.from(y,(v,i)=>v-B[i]));
    console.log(`     참고: Syn2 Shirley, INNO 규약(사다리꼴·끝 3 점) 넓이 ${Q[Q.length-1].toFixed(1)} vs CasaXPS 15089.8`); }
  for(const ft of CR.fits){ const s=await rd(ft.file); const shape=ft.file.includes("GLP")?"GLP":ft.file.includes("LA")?"LA":"GLS";
    const o=shape==="LA"?{F:0.9,alpha:1.53,beta:2.43,wG:0}:{F:1.1,m:0.3}; const cons={a:shape==="LA"?{alpha:{mode:"fixed"},beta:{mode:"fixed"},wG:{mode:"fixed"}}:{m:{mode:"fixed"}}};
    const P=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[i0],hi:s.be[i1],bg:{type:"linear",conv:"casa",aw:1},comps:[C.newComp(shape,285.05,Object.assign({id:"a",area:12000},o))],cons});
    const r=C.fitOnce(P,{}); const q=C.compSummary(P,r,null)[0]; const dA=Math.abs(q.area-ft.area)/ft.area;
    if(shape==="GLP") ok(dA<1e-4&&Math.abs(q.fwhm-ft.fwhm)/ft.fwhm<3e-4&&Math.abs(q.be-ft.pos)<1e-3,"GLP(m 0.3) = CasaXPS GL(30): 넓이, CasaXPS fwhm = 실제 FWHM (상대 < 3e-4)",`넓이 ${q.area.toFixed(1)} vs ${ft.area}, FWHM ${q.fwhm.toFixed(5)} vs ${ft.fwhm}`);
    else if(shape==="LA") ok(dA<2e-4&&Math.abs(P.comps[0].F-1)<1e-3,"LA(α 높은 BE 1.53, β 2.43) = CasaXPS LA(1.53,2.43,0): 같은 모양·꼬리 쪽, 넓이 (상대 < 2e-4)",`넓이 ${q.area.toFixed(1)} vs ${ft.area}; F ${P.comps[0].F.toFixed(4)} · 실제 FWHM ${q.fwhm.toFixed(5)} vs CasaXPS 폭 숫자 ${ft.width} (환산 미확인)`);
    else console.log(`     참고: GLS m 0.3 자료 — CasaXPS SGL(30) 넓이 ${ft.area}, fwhm ${ft.fwhm}, RMS ${ft.rms}: SGL(30) = GLS m ${C.sglToGLS(30).toFixed(4)} (높이 혼합) — 이 판독은 §5 의 CasaXPS 피팅 규약으로 재현`); }
  console.log("== 4 기준 구현 대조: CasaXPS 규약(§2b) · Avantage Shirley(§2c) · Avantage 꼬리(§1.10) · Avantage 투과 함수(§4) · Tougaard T0(§2)");
  for(const c of REF.casa_conv){ const s=await rd(c.file); const R=C.sliceRegion(s.be,s.y,s.be[i0],s.be[i1]); const E=F64(R.E), y=F64(R.y);
    const B=C.bgStatic(E,y,{type:c.type,conv:"casa",aw:c.aw,_ends:C.casaEnds(s.be,s.y,R.i0,R.i1,c.aw)}); const rg=Math.max(...y)-Math.min(...y); let e=0; c.B.forEach((b,i)=>{ e=Math.max(e,Math.abs(B[i]-b)/rg); });
    const st=C.casaRegionStats(E,y,B); const es=Math.max(Math.abs(st.area-c.stats[0])/Math.abs(c.stats[0]),Math.abs(st.maxH-c.stats[1])/rg,Math.abs(st.minH-c.stats[2])/rg);
    ok(e<1e-9&&es<1e-9,`CasaXPS 규약 ${c.file} ${c.type} Av.W ${c.aw}: 배경·Area/Max/Min = 기준 구현`,`${e.toExponential(1)} · ${es.toExponential(1)}`); }
  REF.avantage_shirley.forEach(c=>{ const E=F64(c.E), y=F64(c.y); const B=C.bgStatic(E,y,{type:"shirley",conv:"avantage",n:5,it:c.it,cl:c.cl}); const rg=Math.max(...c.y)-Math.min(...c.y); let e=0; c.B.forEach((b,i)=>{ e=Math.max(e,Math.abs(B[i]-b)/rg); });
    ok(e<1e-10&&(B.fallback||null)===c.fallback,`Avantage Shirley ${c.file} 반복 ${c.it}${c.cl?` · 넓이 변화 ${c.cl}`:""}${c.fallback?" → Linear 로 바뀜":""}`,e.toExponential(1)); });
  REF.agt.forEach(c=>{ const v=C.shapeUnit("AGT",F64(c.x),c.p,{lo:c.lo,hi:c.hi}); const mx=Math.max(...c.v); let e=0; c.v.forEach((r,i)=>{ e=Math.max(e,Math.abs(v[i]-r)/mx); });
    ok(e<1e-6,`Avantage 꼬리 GL·T ${JSON.stringify(c.p)}: 창 안 넓이 1 정규화 = quad`,e.toExponential(1)); });
  { const T=REF.tf_avantage; let e=0; T.rows.forEach(([KE,PE,v])=>{ e=Math.max(e,Math.abs(C.tfValue({model:"avantage",coef:T.coef},KE,PE)-v)/v); }); ok(e<1e-12,"Avantage 투과 함수 PE·10^Σaₖ(log10(KE/PE))ᵏ",e.toExponential(1)); }
  { const T=REF.tougaard_T0, E=F64(T.E), y=F64(T.y), rg=Math.max(...T.y)-Math.min(...T.y);
    T.cases.forEach(c=>{ const B=C.bgStatic(E,y,{type:"tougaard",n:5,C:c.C,D:c.D,u3:c.u3,T0:c.T0}); let e=0; c.B.forEach((b,i)=>{ e=Math.max(e,Math.abs(B[i]-b)/rg); }); ok(e<1e-9,`Tougaard ${c.u3?"U3":"U2"} C ${c.C} D ${c.D} T0 ${c.T0}`,e.toExponential(1)); }); }
  { const a=C.avantageVoigt(0.8,0.3), b=C.voigtFromAvantage(a.fV,a.M); ok(Math.abs(b.wG-0.8)<1e-12&&Math.abs(b.wL-0.3)<1e-12&&Math.abs(C.voigtFromAvantage(1.1,1).wL-1.1)<1e-15,"Avantage Voigt (fV, M) ↔ (fG, fL) 왕복 · M = 1 → fL = fV",JSON.stringify(a)); }
  console.log("== 5 CasaXPS 2.3.27 V2 실측 — CasaXPS 피팅 규약(weight \"casa\"), SGL = 높이 혼합, LA 폭·LA(1,1,n) (SPEC_casa_v2.md · test/ref/ref_casa_v2.py)");
  { const R2=JSON.parse(fs.readFileSync(D("ref/ref_casa_v2.json"),"utf8"));
    const casaFit=async(file,comp,fixed)=>{ const s=await rd(file); const cons={a:{}}; fixed.forEach(k=>cons.a[k]={mode:"fixed"}); comp.id="a";
      const P=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[i0],hi:s.be[i1],bg:{type:"linear",conv:"casa",aw:1},comps:[comp],cons,weight:"casa"});
      const r=C.fitOnce(P,{}); const st=C.fitStats(P,r); return {c:P.comps[0],st,cz:C.compToCasa(P.comps[0])}; };
    const rel=(a,b)=>Math.abs(a-b)/Math.abs(b), e1=x=>x.toExponential(1);
    for(const t of R2.fits){ const k=t.casa_read; let comp,fixed;
      if(t.kind==="GLT"){ comp=C.newComp("AGT",285.05,{area:12000,F:1.1,m:t.q.m,ET:t.q.k,TM:0,CT:0,tu:"fwhm"}); fixed=["m","ET","TM","CT"]; }
      else { comp=C.casaToComp(t.casa,{be:285.05,area:12000,fwhm:t.kind==="LA"?0.9:1.1}).comp; fixed=t.kind==="LA"?["alpha","beta"]:["m"]; }
      const {c,st,cz}=await casaFit(t.file,comp,fixed); const F=c.F;
      const eng=`넓이 ${e1(rel(c.area,t.area))} · F ${e1(rel(F,t.F))} · 위치 ${e1(Math.abs(c.be-t.pos))} eV · RMS ${e1(rel(st.casa.rms,Math.max(t.rms,1e-3)))} · DoF ${st.casa.dof}/${t.dof}`;
      ok(rel(c.area,t.area)<3e-5&&rel(F,t.F)<3e-5&&Math.abs(c.be-t.pos)<3e-5&&(t.rms<1?st.casa.rms<0.05:rel(st.casa.rms,t.rms)<2e-3)&&st.casa.dof===k.dof,`${t.file} ${t.casa}: CasaXPS 규약 피팅 = 기준 구현(scipy); DoF = CasaXPS 판독 (#{y ≠ 배경} − p, SPEC_casa_v3 §2)`,eng);
      const vs=`넓이 ${c.area.toFixed(1)} vs ${k.area} (${e1(rel(c.area,k.area))}), fwhm ${cz.fwhm!=null?cz.fwhm.toFixed(6):"—"} vs ${k.fwhm}, 위치 ${c.be.toFixed(4)} vs ${k.pos}, RMS ${st.casa.rms.toFixed(3)} vs ${k.rms}, DoF ${st.casa.dof} vs ${k.dof}`;
      if(t.kind==="SGL") ok(rel(c.area,k.area)<1e-3&&rel(cz.fwhm,k.fwhm)<1e-3&&Math.abs(c.be-k.pos)<1.5e-4&&rel(st.casa.rms,k.rms)<2e-3&&st.casa.dof===k.dof&&cz.notation===t.casa,`${t.file}: ${t.casa} = GLS m ${c.m.toFixed(4)} (높이 혼합) — 넓이 ≤ 0.1 %, fwhm ≤ 1e-3, RMS ≤ 0.2 %, DoF 같음`,vs);
      else if(t.kind==="LA"&&t.q.a===t.q.b) ok(rel(c.area,k.area)<1.5e-3&&rel(cz.fwhm,k.fwhm)<6e-4&&Math.abs(c.be-k.pos)<1.5e-4&&rel(st.casa.rms,k.rms)<0.05&&st.casa.dof===k.dof&&cz.exact,`${t.file}: ${t.casa} — CasaXPS fwhm = 실제 FWHM (≤ 6e-4), 넓이 ≤ 1.5e-3`,vs);
      else if(t.kind==="LA"){ ok(rel(c.area,k.area)<3e-4&&Math.abs(c.be-k.pos)<1.5e-4&&rel(cz.fwhm,k.fwhm)<6e-4&&Math.abs(cz.asym-k.asym)<1e-3&&cz.exact&&cz.fwhmExact===true,`${t.file}: ${t.casa} — 넓이·위치 같음, CasaXPS fwhm = 2·min(l, h), asym = 1 − l/h, 최댓값 = 격자 k·F/40 포물선 꼭짓점 (≤ 6e-4 · 1e-3, 정함; SPEC_casa_v4 §3)`,vs+` · asym ${cz.asym.toFixed(4)} vs ${k.asym}`); }
      else { ok(rel(c.area,k.area)<3e-5&&Math.abs(c.be-k.pos)<1e-4&&rel(F,1.2)<1e-4&&cz.notation==="GL(30)T(2)"&&rel(cz.fwhm,k.fwhm)<6e-4&&Math.abs(cz.asym-k.asym)<1e-3&&cz.exact&&cz.fwhmExact===true,`${t.file}: GL(30)T(2) = AGT(tu fwhm, ET 2, TM 0) — 넓이·위치·F 1.2 회복; CasaXPS fwhm·asym 격자 포물선 규칙 (≤ 6e-4 · 1e-3, 정함; SPEC_casa_v4 §3)`,vs+` · asym ${cz.asym.toFixed(4)} vs ${k.asym}`);
        console.log(`     남은 차이(CasaXPS 쪽 수치 모양 오차): RMS ${k.rms} (여기 ${st.casa.rms.toExponential(1)}) — SPEC_casa_v3 §2`); } }
    for(const t of R2.scan){ const k=t.casa_read; const {c,st,cz}=await casaFit("V2_Voigt_wG0.6_wL0.4",C.casaToComp(`LA(1,1,${t.n})`,{be:285.05,area:12000,fwhm:0.8}).comp,["alpha","beta"]);
      ok(rel(c.area,t.area)<3e-5&&rel(c.F,t.F)<1e-4&&rel(st.casa.rms,t.rms)<2e-3,`LA(1,1,${t.n}) on Voigt(wG 0.6, wL 0.4): 기준 구현 (wG = n·F/144.87)`,`넓이 ${e1(rel(c.area,t.area))} · F ${e1(rel(c.F,t.F))} · RMS ${e1(rel(st.casa.rms,t.rms))}`);
      ok(rel(cz.fwhm,k.fwhm)<2.5e-4&&rel(c.area,k.area)<1.5e-3&&rel(st.casa.rms,k.rms)<0.03&&cz.notation===`LA(1,1,${t.n})`,`LA(1,1,${t.n}) vs CasaXPS: fwhm ≤ 2.5e-4, 넓이 ≤ 1.5e-3, RMS ≤ 3 %`,`fwhm ${cz.fwhm.toFixed(6)} vs ${k.fwhm}, 넓이 ${c.area.toFixed(1)} vs ${k.area}, RMS ${st.casa.rms.toFixed(1)} vs ${k.rms}`); }
    { const s=await rd("C1s"); const R=C.sliceRegion(s.be,s.y,s.be[i0],s.be[i1]); const E=F64(R.E), y=F64(R.y);
      const st=(cb)=>C.casaRegionStats(E,y,C.bgStatic(E,y,Object.assign({},cb.bg,{_ends:C.casaEnds(s.be,s.y,R.i0,R.i1,cb.bg.aw)})));
      const u3=C.casaBg("U 3 Tougaard","0, 551, 436, 3",1), u3d=C.casaBg("U 3 Tougaard","",1), q=st(u3);
      ok(u3.bg.type==="linear"&&!u3.exact&&u3d.bg.type==="linear"&&rnd(q.maxH,14451.9)&&rnd(q.minH,-482.0)&&rel(q.area,39564.8)<1e-4,"CasaXPS 'U 3 Tougaard' 기본 Cross Section (B 0, 551, 436, 3) → Linear 로 해석, exact=false (미검증)",`Area ${q.area.toFixed(1)} · Max ${q.maxH.toFixed(1)} · Min ${q.minH.toFixed(1)} vs 39564.8 · 14451.9 · −482.0`);
      const u3b=C.casaBg("U 3 Tougaard","1200, 551, 436, 3",1), u2=C.casaBg("U 2 Tougaard","2866,1643,1,0",1);
      ok(u3b.bg.type==="tougaard"&&u3b.bg.u3&&u3b.bg.C===551&&u3b.bg.D===436&&u3b.bg.T0===3&&!u3b.exact&&u2.bg.type==="tougaard"&&!u2.bg.u3&&u2.bg.C===1643&&!u2.exact&&C.casaBg("Shirley",null,3).bg.aw===3,"CasaXPS 배경 표기: 명시 B > 0 → INNO Tougaard(C, D, T0), '정확하지 않음' 표시",JSON.stringify(u3b.bg.casaTougaard)); }
    { let e=0; [0,5,19.5,30,50,80,100].forEach(p=>{ e=Math.max(e,Math.abs(C.glsToSGL(C.sglToGLS(p))-p)); }); const a=C.casaToComp("LA(50)",{fwhm:1}).comp, b=C.casaToComp("LA(1,1283)",{fwhm:1}).comp, d=C.casaToComp("LA(1.53,2.43)",{fwhm:1}).comp, g=C.casaToComp("GL(30)",{fwhm:1.15343}).comp;
      ok(e<1e-12&&a.casaN===700.5&&b.casaN===1283&&b.alpha===1&&d.beta===2.43&&d.casaN===0&&Math.abs(g.F-1.2)<5e-4&&C.casaToComp("XYZ(1)")===null,"표기 해석: SGL ↔ GLS 왕복 · LA(m) = LA(1,1,1401−14.01m) · LA(α,n) · LA(α,β) · GL(p) fwhm → F",`SGL 왕복 ${e.toExponential(1)} · GL(30) fwhm 1.15343 → F ${g.F.toFixed(5)}`); }
  }
})().then(()=>{
console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
}).catch(e=>{ console.error(e); process.exit(2); });
