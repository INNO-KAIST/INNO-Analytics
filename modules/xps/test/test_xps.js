/* XPS 모듈 회귀·검증 시험  —  node test/test_xps.js  (modules/xps 에서)
   1) 엔진 ↔ 독립 Python 기준 구현(test/ref/ref.json: scipy quad·voigt_profile·wofz·least_squares) 대조
   2) 참값을 아는 합성 데이터로 불확도 보정(coverage) 시험
   3) 파일 읽기 회귀(같은 측정의 txt·xls·vms 동일성, VGD 축·계수)
   4) 라이브러리 무결성 · 서베이 원소 판별(실측은 반례 감시용)            */
const path=require("path"), fs=require("fs");
const D=p=>path.join(__dirname,p);
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js"), LIB=require("../src/xps_lib.js");
const REF=JSON.parse(fs.readFileSync(D("ref/ref.json"),"utf8"));
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond){ pass++; } else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!cond) console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const F64=a=>Float64Array.from(a);

(async()=>{
console.log("== 1a Faddeeva / erfcx vs scipy.special.wofz / erfcx");
{ let mx=0,me=0; REF.wofz.forEach(r=>{ if(r[0]==="erfcx"){ const v=C.erfcx(r[1]); me=Math.max(me,Math.abs(v-r[2])/r[2]); return; } const [x,y,wr,wi]=r; const w=C.faddeeva(x,y); mx=Math.max(mx,Math.hypot(w[0]-wr,w[1]-wi)/Math.hypot(wr,wi)); });
  ok(mx<2e-6,"Faddeeva w(z) 상대오차 < 2e-6 (64 점, Im z ≥ 0)","max "+mx.toExponential(2)); ok(me<2e-6,"erfcx 상대오차 < 2e-6","max "+me.toExponential(2)); }

console.log("== 1b 선모양 (넓이 1 정규화·합성곱) vs quad 기준");
REF.shapes.forEach(s=>{ const v=C.shapeUnit(s.shape,F64(s.x),s.p,{lo:s.lo,hi:s.hi}); const mx=Math.max(...s.v.map(Math.abs)); let e=0; s.v.forEach((r,i)=>e=Math.max(e,Math.abs(v[i]-r)/mx));
  const tol=(s.shape==="GLS"||s.shape==="GLP"||s.shape==="V")?1e-6:1e-3; ok(e<tol,`${s.shape} ${JSON.stringify(s.p)}`,`최대오차/최댓값 ${e.toExponential(2)} (허용 ${tol})`); });
// 넓이 1 (ℝ) — DS 는 영역 정의라 제외
[["GLS",{F:1.2,m:0.3}],["GLP",{F:1.2,m:0.6}],["V",{wG:0.7,wL:0.4}],["LA",{F:0.9,alpha:1.1,beta:2.2,wG:0.4}],["GLT",{F:1,m:0.3,T:0.4,tau:1.2}]].forEach(([s,p])=>{
  const R=4000,N=800001; const x=new Float64Array(N); for(let i=0;i<N;i++) x[i]=-R+2*R*i/(N-1); const v=C.shapeUnit(s,x,p,{lo:-30,hi:30}); let a=0; for(let i=1;i<N;i++) a+=0.5*(v[i]+v[i-1])*(x[i]-x[i-1]);
  const tail=s==="GLS"||s==="V"||(s==="GLP")?2/(Math.PI*R/(p.F||p.wL)) :0; ok(Math.abs(a-1)<3e-3,`${s} ℝ 넓이 = 1`,`±${R} eV 적분 ${a.toFixed(5)}`); });

console.log("== 1c 배경 vs 독립 구현");
REF.bgs.forEach(b=>{ const B=C.bgStatic(F64(b.E),F64(b.y),b.bg); const rg=Math.max(...b.y)-Math.min(...b.y); let e=0; b.B.forEach((r,i)=>e=Math.max(e,Math.abs(B[i]-r)/rg));
  ok(e<1e-6,`배경 ${JSON.stringify(b.bg)}`,`최대오차/세기범위 ${e.toExponential(2)}`); });

console.log("== 1d 피팅 (LM·sin/√ 변환) vs scipy least_squares(TRF·경계)");
for(const f of REF.fits){
  const P=f.prob; const prob={E:F64(P.E),y:F64(P.y),c:P.c,bg:C.clone(P.bg),comps:C.clone(P.comps),cons:C.clone(P.cons)};
  const t0=Date.now(); const r=C.fitOnce(prob,{}); const ms=Date.now()-t0; const cov=C.covariance(prob,r);
  const map=r.map; const R=f.ref; let worst=0, worstK="", sigE=0;
  R.free.forEach(([ci,k],j)=>{ const jj=map.free.findIndex(q=>q.ci===ci&&q.k===k); if(jj<0){ worst=Infinity; worstK=ci+"."+k+" 없음"; return; }
    const xj=ci<0?prob.bg[k]:C.getP(prob.comps[ci],k); const d=Math.abs(xj-R.x[j])/(R.sig[j]||1); if(d>worst){ worst=d; worstK=(ci<0?"bg":prob.comps[ci].id)+"."+k; }
    const s=cov.sig[jj]; if(s!=null) sigE=Math.max(sigE,Math.abs(s-R.sig[j])/R.sig[j]); });
  const dchi=Math.abs(r.chi2-R.chi2)/R.chi2;
  ok(worst<0.05&&dchi<1e-4,`${f.name}: 매개변수 일치 |Δ|/σ < 0.05, χ² 상대차 < 1e-4`,`최대 ${worst.toFixed(4)}σ (${worstK}), χ² ${r.chi2.toFixed(3)} vs ${R.chi2.toFixed(3)}, ${ms} ms, LM ${r.iter}회`);
  ok(sigE<0.02,`${f.name}: 공분산 σ 일치 (상대 < 2 %)`,`최대 ${(100*sigE).toFixed(2)} %`);
}

console.log("== 1e TPP-2M vs 독립 구현");
{ let e=0; REF.tpp2m.forEach(t=>{ e=Math.max(e,Math.abs(C.tpp2m(t.E,t.mat)-t.lam)/t.lam); }); ok(e<1e-12,"TPP-2M 식 일치","max rel "+e.toExponential(2));
  const si=C.tpp2m(1000,{rho:2.33,Nv:4,M:28.086,Eg:1.1}); ok(si>20&&si<32,"Si 1000 eV IMFP 가 문헌 범위(2–3 nm) 안",si.toFixed(2)+" Å"); }

console.log("== 2 불확도 보정: 참값을 아는 합성 데이터 (active Shirley = 참 배경 모형)");
{ const R=C.rng(2026); const pois=l=>{ if(l>40) return Math.max(0,Math.round(l+Math.sqrt(l)*R.g())); let L=Math.exp(-l),k=0,p=1; do{k++;p*=R.u();}while(p>L); return k-1; };
  const n=381, E=Float64Array.from({length:n},(_,i)=>279+i*0.05);
  const truth=[C.newComp("GLS",284.8,{F:1.1,m:0.3,area:9000,id:"a"}),C.newComp("GLS",286.3,{F:1.3,m:0.3,area:2500,id:"b"}),C.newComp("GLS",288.8,{F:1.2,m:0.3,area:1500,id:"c"})];
  const ctx={emin:E[0],emax:E[n-1]}; const Pt=new Float64Array(n); truth.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) Pt[i]+=v[i]; }); const Q=C.cumtrapz(E,Pt); const mu=Float64Array.from(Pt,(p,i)=>2000+0.08*Q[i]+p);
  const K=+process.env.NCOV||120, cf=0.3; const z={be:[],area:[],F:[]}; let inside={be:0,area:0,F:0}, tot=0, t0=Date.now();
  for(let k=0;k<K;k++){ const y=Float64Array.from(mu,m=>pois(m*cf)/cf);
    const st=truth.map(c=>Object.assign(C.clone(c),{be:c.be+0.15*R.g(),F:c.F*(1+0.1*R.g())}));
    const prob={E,y,c:cf,bg:{type:"active",n:5},comps:st,cons:{b:{m:{mode:"link",ref:"a",k:1}},c:{m:{mode:"link",ref:"a",k:1}}}};
    const r=C.fitOnce(prob,{}); const cov=C.covariance(prob,r); const sm=C.compSummary(prob,r,cov);
    prob.comps.forEach((c,i)=>{ ["be","area"].forEach(q=>{ const s=sm[i].sig[q]; const zz=(c[q]-truth[i][q])/s; z[q].push(zz); if(Math.abs(zz)<1.96) inside[q]++; }); const jF=r.map.free.findIndex(f=>f.ci===i&&f.k==="F"); const zF=(c.F-truth[i].F)/cov.sig[jF]; z.F.push(zF); if(Math.abs(zF)<1.96) inside.F++; tot++; });
  }
  const stat=a=>{ const m=a.reduce((s,v)=>s+v,0)/a.length; return {m,sd:Math.sqrt(a.reduce((s,v)=>s+(v-m)*(v-m),0)/(a.length-1))}; };
  ["be","area","F"].forEach(q=>{ const s=stat(z[q]); const cvg=inside[q]/tot; const se=Math.sqrt(0.95*0.05/tot);
    ok(Math.abs(cvg-0.95)<3.5*se+0.01&&Math.abs(s.m)<0.25&&Math.abs(s.sd-1)<0.15,`95 % 구간 적중률 (${q}) — ${K}개 데이터 × 3 성분`,`${(100*cvg).toFixed(1)} %, z 평균 ${s.m.toFixed(3)}, z 표준편차 ${s.sd.toFixed(3)}`); });
  console.log("    ("+((Date.now()-t0)/1000).toFixed(1)+" s)"); }

console.log("== 2b 다중 시작이 국소 최소에서 벗어나는지");
{ const R=C.rng(7); const n=301, E=Float64Array.from({length:n},(_,i)=>525+i*0.05);
  const truth=[C.newComp("GLS",529.6,{F:1.0,m:0.2,area:30000,id:"a"}),C.newComp("GLS",531.4,{F:1.4,m:0.2,area:12000,id:"b"})]; const ctx={emin:E[0],emax:E[n-1]};
  const mu=new Float64Array(n).fill(20000); truth.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) mu[i]+=v[i]; }); const y=Float64Array.from(mu,m=>m+Math.sqrt(m/0.3)*R.g());
  const bad=[C.newComp("GLS",532.9,{F:0.4,m:0.2,id:"a"}),C.newComp("GLS",533.6,{F:0.4,m:0.2,id:"b"})];
  const prob={E,y,c:0.3,bg:{type:"linear",n:5},comps:bad,cons:{}};
  const r1=C.fitOnce(Object.assign({},prob,{comps:C.clone(bad)}),{});
  const b=C.fitMulti(Object.assign({},prob,{comps:C.clone(bad)}),{starts:16,seed:3});
  ok(b.r.chi2p1<=r1.chi2p1+1e-6,"다중 시작 χ² ≤ 단일 시작 χ² (1차 목적함수)",`단일 ${r1.chi2p1.toFixed(1)} → 다중 ${b.r.chi2p1.toFixed(1)} (같은 해 비율 ${(100*b.r.multi.frac).toFixed(0)} %)`);
  const be=b.P.comps.map(c=>c.be).sort((x,y)=>x-y); ok(Math.abs(be[0]-529.6)<0.05&&Math.abs(be[1]-531.4)<0.1,"다중 시작이 참 위치 근처로",be.map(v=>v.toFixed(3)).join(", ")); }

console.log("== 3 파일 읽기 회귀");
{ const d=D("data/bct600/"); const x=(await IO.readFile("a.xls",new Uint8Array(fs.readFileSync(d+"1_D.I._on_BCT_600C_5h.xls")))).spectra;
  ok(x.length===7,"Avantage xls: 6 영역 + survey",x.map(s=>s.region).join(", "));
  for(const r of ["C1s","O1s","Co2p","Co3p","Ta4f","Ba4d"]){
    const t=(await IO.readFile("a.txt",new Uint8Array(fs.readFileSync(d+"_D.I._on_BCT__600_C_5h___"+r+".txt")))).spectra[0];
    const v=(await IO.readFile("a.vms",new Uint8Array(fs.readFileSync(d+"_D.I._on_BCT__600_C_5h___"+r+".txt.vms")))).spectra[0];
    const s=x.find(q=>q.region.replace(" ","")===r); let m=0,e=0; for(let i=0;i<t.n;i++){ m=Math.max(m,Math.abs(t.y[i]-s.y[i]),Math.abs(t.y[i]-v.y[i])); e=Math.max(e,Math.abs(t.be[i]-v.be[i]),Math.abs(t.be[i]-s.be[i])); }
    ok(t.n===v.n&&t.n===s.n&&m===0&&e<1e-9,`${r}: txt = xls = VAMAS (세기 완전 일치, BE < 1e-9 eV)`,`n ${t.n}`); }
  const sv=x.find(s=>s.reg&&s.reg.kind==="survey"); ok(sv&&Math.abs(sv.step-1)<1e-9&&sv.n===1361,"survey 1 eV 간격 1361 점");
  const k=D("data/k2312/"); const ex={"C1s_Scan.VGD":[381,279.08,298.08,0.2],"O1s_Scan.VGD":[401,525.08,545.08,0.15],"Ba4d_Scan.VGD":[341,78.08,95.08,0.25],"Co3p_Scan.VGD":[401,52.08,72.08,0.35],"Ta4d_Scan.VGD":[521,222.08,248.08,0.35]};
  for(const [f,[n,a,b,c]] of Object.entries(ex)){ const s=(await IO.readFile(f,new Uint8Array(fs.readFileSync(k+f)))).spectra[0];
    ok(s.n===n&&Math.abs(s.be[0]-a)<1e-6&&Math.abs(s.be[n-1]-b)<1e-6&&Math.abs(s.c-c)<1e-9&&s.hv===1486.68&&s.meta.passEnergy===50,`VGD ${f}: 축·dwell×scans·hν·PE`,`n ${s.n}, ${s.be[0].toFixed(2)}–${s.be[n-1].toFixed(2)}, c ${s.c}, PE ${s.meta.passEnergy}`);
    let mx=0; for(let i=0;i<n;i++) mx=Math.max(mx,Math.abs(s.y[i]*s.c-s.yRaw[i])); ok(mx<1e-9,`VGD ${f}: CPS × c = VGData`); }
  const s=(await IO.readFile("C1s_Scan.VGD",new Uint8Array(fs.readFileSync(k+"C1s_Scan.VGD")))).spectra[0];
  const T=[200,600,1000,1400].map(KE=>C.tfValue(s.meta.tf,KE,50)); ok(T.every((t,i)=>i===0||t<T[i-1]),"VGD transmission 계수 → T(KE) 단조 감소 (CAE 모드 물리)",T.map(t=>t.toFixed(2)).join(" > "));
}

console.log("== 4 라이브러리");
{ ok(Object.keys(LIB.core).length===92,"핵심 준위 92 원소"); ok(LIB.core.C.lv["1s"]===284.2&&LIB.core.O.lv["1s"]===543.1,"C 1s 284.2 · O 1s 543.1 (X-ray Data Booklet, 원소 상태)");
  ok(LIB.scofield.al["C 1s"]===1&&LIB.scofield.al["O 1s"]===2.93&&LIB.scofield.al["Co 2p3/2"]===12.62,"Scofield Al: C 1s 1, O 1s 2.93, Co 2p3/2 12.62");
  ok(Math.abs(LIB.scofield.al["Ni 2p"]-(LIB.scofield.al["Ni 2p1/2"]+LIB.scofield.al["Ni 2p3/2"]))<1e-9,"doublet 합 = 두 subshell 합 (Ni 2p 재계산)",String(LIB.scofield.al["Ni 2p"]));
  ok((LIB.core.Co.soWarn||[]).includes("3p"),"Booklet 의 Co 3p 순서 역전 표시");
  ok(LIB.auger.length>200&&LIB.auger.some(a=>a.el==="O"&&a.ke===511),"Auger (Avantage KE) 포함");
  for(const [nm,fp] of [["BCT 600 °C (xls)","data/bct600/1_D.I._on_BCT_600C_5h.xls"],["EtOH D.I. (VGD)","data/k2312/Survey.VGD"]]){
    const sv=(await IO.readFile(fp.split("/").pop(),new Uint8Array(fs.readFileSync(D(fp))))).spectra.find(s=>s.reg.kind==="survey");
    const pk=C.findPeaks(F64(sv.be),F64(sv.y),{width:1.5,snr:4,merge:1.5}); const cand=C.identify(pk,LIB,{hv:sv.hv,shift:0,tol:2.0}); const acc=cand.filter(C.acceptDefault).map(c=>c.el);
    const co=cand.find(c=>c.el==="Co");
    ok(["Ba","O","C","Ta"].every(e=>acc.includes(e))&&!!co,`실측 survey ${nm}(반례 감시): Ba·O·C·Ta 채택, Co 는 후보에 있음(2p 가 Ba 3d 와 겹쳐 survey 만으로는 확정 불가 — 영역 스펙트럼으로)`,"채택 "+acc.join(" ")+" · Co 점수 "+(co?co.score.toFixed(2):"없음")+" · 미배정 봉우리 "+cand.unassigned.length); }
}

console.log("== 5 서베이 원소 판별 — 합성 survey 150 개 (참값 원소 3–8 개, 화학 이동 −1.5…+3.5 eV, 검출 가능 = 주선 SNR > 15)");
{ const B=require("./survey_bench.js"); const r=B.run(150,"greedy",{});
  ok(r.recall>=0.95,"검출 가능한 참 원소의 채택 비율 ≥ 95 %",(100*r.recall).toFixed(1)+" %");
  ok(r.fpPerSpectrum<=0.2,"없는 원소 채택 ≤ 0.2 개/스펙트럼",r.fpPerSpectrum.toFixed(2)+" 개, 정밀도 "+(100*r.precisionLike).toFixed(1)+" %, 흔한 오채택 "+r.top_fp.slice(0,4).map(x=>x[0]+"×"+x[1]).join(" ")); }

console.log("== 6 두 번째 독립 기준 구현(test/ref/xps_ref.py — SPEC 만 보고 JS 를 보지 않고 작성, fixtures/*.json) 대조");
{ const FX=n=>JSON.parse(fs.readFileSync(D("ref/fixtures/"+n),"utf8"));
  // 6a 선모양 41 종
  const LS=FX("lineshapes.json").cases; const lsE={a:0,c:0}; let fwE=0;
  LS.forEach(cs=>{ const x=F64(cs.x); const ref=cs.shape==="DS"&&cs.f_region_norm?cs.f_region_norm:cs.f; const ctx=cs.region?{lo:cs.region[0],hi:cs.region[1]}:{lo:x[0],hi:x[x.length-1]};
    const v=C.shapeUnit(cs.shape,x,cs.params,ctx); const mx=Math.max(...ref.map(Math.abs)); let e=0; ref.forEach((r,i)=>e=Math.max(e,Math.abs(v[i]-r)/mx));
    const an=cs.shape==="GLS"||cs.shape==="GLP"||((cs.shape==="LA"||cs.shape==="DS")&&!(cs.params.wG>0)); if(cs.shape==="V") lsE.v=Math.max(lsE.v||0,e); else if(an) lsE.a=Math.max(lsE.a,e); else lsE.c=Math.max(lsE.c,e);
    if(cs.fwhm_true!=null&&cs.shape!=="DS"){ const t=C.trueFWHM(cs.shape,cs.params); if(t.fwhm) fwE=Math.max(fwE,Math.abs(t.fwhm-cs.fwhm_true)/cs.fwhm_true); } });
  ok(lsE.a<1e-10,`선모양 해석형(GLS·GLP·LA wG=0·DS wG=0) ${LS.filter(c=>c.shape==="GLS"||c.shape==="GLP"||((c.shape==="LA"||c.shape==="DS")&&!(c.params.wG>0))).length} 종 < 1e-10`,"최대 "+lsE.a.toExponential(2));
  ok((lsE.v||0)<2e-6,"Voigt (Weideman vs wofz) < 2e-6","최대 "+(lsE.v||0).toExponential(2));
  ok(lsE.c<2e-5,"합성곱형(DS·LA wG>0·GLT) < 2e-5 (SPEC 허용 1e-3)","최대 "+lsE.c.toExponential(2));
  ok(fwE<2e-3,"실제 FWHM (수치 반치폭) 상대차 < 2e-3","최대 "+fwE.toExponential(2));
  // 6b 배경 47 종 (실측 8 영역 × 5 + 합성 7)
  const BG=FX("backgrounds.json").cases; let bmax={}; let bad=[];
  BG.forEach(cs=>{ const bg={type:cs.type,n:cs.n_end}; if(cs.type==="tougaard"){ bg.C=cs.params.C; bg.D=cs.params.D||0; bg.u3=cs.params.kind==="U3"; }
    const B=C.bgStatic(F64(cs.E),F64(cs.y),bg); const dy=Math.abs(cs.endpoints.y_hi-cs.endpoints.y_lo)||1; const ymx=Math.max(...cs.y.map(Math.abs));
    let e=0; cs.B.forEach((r,i)=>e=Math.max(e,Math.abs(B[i]-r)/(cs.type==="shirley"||cs.type==="smart"?dy:ymx)));
    const tol=cs.type==="shirley"||cs.type==="smart"?2e-6:1e-9; bmax[cs.type]=Math.max(bmax[cs.type]||0,e); if(!(e<tol)) bad.push(cs.id+" "+e.toExponential(2)); });
  ok(!bad.length,`배경 ${BG.length} 종 (Shirley·Smart < 2e-6·Δy, Linear·Tougaard < 1e-9)`,Object.entries(bmax).map(([k,v])=>k+" "+v.toExponential(1)).join(", ")+(bad.length?" | "+bad.slice(0,4).join("; "):""));
  // 6c 피팅 10 문제 (합성 7 + 실측 3)
  const FT=FX("fits.json").cases;
  for(const cs of FT){
    const M=cs.model, idOf=n=>n; const comps=[], cons={};
    const put=(cid,k,P,setv)=>{ if(P.mode==="fixed"){ setv(P.value); (cons[cid]=cons[cid]||{})[k]={mode:"fixed"}; } else if(P.mode==="link"){ (cons[cid]=cons[cid]||{})[k]={mode:"link",ref:idOf(P.ref),k:P.k!=null?P.k:1,b:P.b!=null?P.b:0}; } else { setv(P.value); const c={mode:"free"}; if(P.lo!=null) c.lo=P.lo; if(P.hi!=null) c.hi=P.hi; if(k==="be"&&P.lo==null){ c.lo=-Infinity; c.hi=Infinity; } (cons[cid]=cons[cid]||{})[k]=c; } };
    M.components.forEach(m=>{ const c={id:m.name,name:m.name,shape:m.shape,be:0,area:1,doublet:null};
      put(m.name,"be",m.E0,v=>c.be=v); put(m.name,"area",m.A,v=>c.area=v);
      C.SHAPES[m.shape].p.forEach(k=>{ if(m[k]) put(m.name,k,m[k],v=>c[k]=v); });
      if(m.doublet){ c.doublet={dE:m.doublet.delta.value!=null?m.doublet.delta.value:1,r:m.doublet.r.value!=null?m.doublet.r.value:0.5,q:m.doublet.q.value!=null?m.doublet.q.value:1};
        put(m.name,"dE",m.doublet.delta,v=>c.doublet.dE=v); put(m.name,"r",m.doublet.r,v=>c.doublet.r=v); put(m.name,"q",m.doublet.q,v=>c.doublet.q=v); }
      comps.push(c); });
    // 연결된 값 초기화(기준의 값)
    const b=M.background; let bg;
    if(b.type==="active"){ bg={type:"active",n:5,c0:b.c0.value,k:b.k.value,s:b.s?b.s.value||0:0,cons:{}}; ["c0","k","s"].forEach(k=>{ const P=b[k]; if(!P) return; bg.cons[k]=P.mode==="fixed"?{mode:"fixed"}:Object.assign({mode:"free"},P.lo!=null?{lo:P.lo}:{},P.hi!=null?{hi:P.hi}:{}); }); }
    else { bg={type:b.type,n:b.n_end||5}; if(b.type==="tougaard"){ bg.C=b.params.C; bg.D=b.params.D||0; bg.u3=b.params.kind==="U3"; } }
    const prob={E:F64(cs.E),y:F64(cs.y),c:cs.c,bg,comps,cons};
    // (가) 최적화 비교: 1차 목적함수(분산 = 평활 데이터, 두 구현이 같은 정의)의 전역 최솟값. 2차 분산 v=μ₁ 은 1차 해에 따라 달라 2차 χ² 는 직접 비교하지 않는다
    const t0=Date.now(); const best=C.fitMulti(prob,{starts:24,seed:11,irls:false}); const ms=Date.now()-t0;
    const R=cs.result; const conv=/DS|LA|GLT/.test(JSON.stringify(M.components.map(m=>m.shape)));
    const tolX=conv?0.2:0.02, tolC=conv?5e-3:1e-6;
    const dchi=(best.r.chi2-R.chi2_pass1)/R.chi2_pass1;
    ok(dchi<tolC,`${cs.id}: 1차 목적함수 최솟값 JS ≤ 기준 (상대 ${tolC})`,`${best.r.chi2.toFixed(3)} vs ${R.chi2_pass1.toFixed(3)} (${dchi.toExponential(1)})${dchi<-Math.max(tolC,1e-6)?" — JS 가 더 낮음(기준이 1차 최솟값을 놓침)":""}, 다중 시작 24 · ${ms} ms`);
    const vref=F64(R.v);
    // (나) 기준 해에서 같은 v_ref 로 LM 을 시작하면 그 자리에 머무는가 + 공분산 σ 일치
    const P2=C.clone(Object.assign({},prob,{E:undefined,y:undefined})); P2.E=prob.E; P2.y=prob.y;
    P2.comps.forEach(c=>{ const pr=R.params[c.id]||{}; Object.entries(pr).forEach(([k,v])=>{ const kk={E0:"be",A:"area",delta:"dE"}[k]||k; if(v!=null) C.setP(c,kk,v); }); });
    if(R.params.bg&&P2.bg.type==="active") Object.entries(R.params.bg).forEach(([k,v])=>{ if(v!=null) P2.bg[k]=v; });
    const map=C.buildMap(P2); C.applyLinks(P2,map); const r2=C.lmFit(P2,map,vref,{}); r2.v=vref; r2.map=map; const cov2=C.covariance(P2,r2);
    let worst=0, wk="", sE=0; R.free.forEach((fk,j)=>{ const [cn,pn]=fk.split("."); const k={E0:"be",A:"area",delta:"dE"}[pn]||pn;
      let xj, jj; if(cn==="bg"){ xj=P2.bg[k]; jj=map.free.findIndex(f=>f.ci<0&&f.k===k); } else { const ci=P2.comps.findIndex(c=>c.id===cn); xj=C.getP(P2.comps[ci],k); jj=map.free.findIndex(f=>f.ci===ci&&f.k===k); }
      const sg=R.sigma_vec[j]; const d=Math.abs(xj-R.theta[j])/(sg||1); if(d>worst){ worst=d; wk=fk; } if(jj>=0&&cov2.sig[jj]!=null&&sg) sE=Math.max(sE,Math.abs(cov2.sig[jj]-sg)/sg); });
    const dchi2=Math.abs(r2.chi2-R.chi2)/R.chi2;
    ok(worst<tolX&&dchi2<tolC,`${cs.id}: 기준 해에서 시작한 LM 이 그 자리 (|Δθ| < ${tolX}σ, χ² 상대 < ${tolC})`,`최대 ${worst.toFixed(4)}σ (${wk}), χ² 차 ${dchi2.toExponential(1)}`);
    ok(sE<0.02,`${cs.id}: 공분산 σ 상대차 < 2 % (같은 해·같은 가중)`,`최대 ${(100*sE).toFixed(3)} %`);
  }
  // 6d 정량 12 조합 (T 4종 × λ 3종)
  const Q=FX("quant.json"); let qe=0;
  Q.toy_cases.forEach(tc=>{ const tf=tc.T.type==="none"?{model:"none"}:tc.T.type==="power"?{model:"power",x:tc.T.x,pe:tc.PE}:tc.T.type==="thermo"?{model:"thermo",coef:tc.T.a,pe:tc.PE}:{model:"table",table:tc.T.KE.map((k,i)=>[k,tc.T.T[i]])};
    const im=tc.lam.type==="none"?{model:"none"}:tc.lam.type==="power"?{model:"power",n:tc.lam.n}:{model:"tpp2m",mat:{rho:tc.lam.rho,Nv:tc.lam.Nv,M:tc.lam.M,Eg:tc.lam.Eg}};
    const rows=Q.toy_peaks.map(pk=>({label:pk.name,area:pk.A,sArea:pk.sigma_A,S:pk.S,BE:pk.BE,hv:tc.hnu,PE:tc.PE,tf}));
    const out=C.quantify(rows,{imfp:im}); out.forEach((o,i)=>{ const R=tc.result[i]; [["T","T"],["L","lam"],["at","at"],["sAt","sigma_at"]].forEach(([a,b])=>{ qe=Math.max(qe,Math.abs(o[a]-R[b])/Math.abs(R[b])); }); }); });
  ok(qe<1e-10,"정량 12 조합 (T·λ·at %·σ_at %) 상대차 < 1e-10","최대 "+qe.toExponential(2));
  let te=0; Q.tpp2m_materials.forEach(m=>m.E.forEach((e,i)=>{ te=Math.max(te,Math.abs(C.tpp2m(e,{rho:m.rho,Nv:m.Nv,M:m.M,Eg:m.Eg})-m.lambda_A[i])/m.lambda_A[i]); }));
  ok(te<1e-12,`TPP-2M ${Q.tpp2m_materials.length} 물질 (Au·Si·SiO₂·Ta₂O₅·Co₃O₄·BaTiO₃·PMMA)`,"최대 "+te.toExponential(2)); }

ok(IO.pathSample("C:\\DATA\\K231213C\\Experiment2\\X-Ray018 400um - FG  ON\\W1\\Group\\C1s Scan.VGD")==="K231213C / Experiment2 / X-Ray018 400um - FG ON / W1"
   &&IO.pathSample("C:\\DATA\\K231213C\\Experiment2\\X-Ray018 400um - FG  ON\\W1\\Survey.VGD")==="K231213C / Experiment2 / X-Ray018 400um - FG ON / W1"
   &&IO.pathSample("C:\\DATA\\K240101A\\Experiment1\\X-Ray018 400um - FG  ON\\W1\\Group\\C1s Scan.VGD")!==IO.pathSample("C:\\DATA\\K231213C\\Experiment2\\X-Ray018 400um - FG  ON\\W1\\Group\\C1s Scan.VGD"),
   "VGD 시료 키 = 자료 / 실험 / 측정 조건 / 지점 (Group 폴더·Survey 위치 무관, 다른 실험의 같은 지점 이름은 다른 시료)");
console.log("== 4b 실측 자동 제안 반례 감시 — K2312 Ba 4d (VGD): 작은 봉우리가 연결 기준이 되어 창 밖·넓이 0 으로 남던 결함(10-16)");
{ const sp=(await IO.readFile("Ba4d_Scan.VGD",new Uint8Array(fs.readFileSync(D("data/k2312/Ba4d_Scan.VGD"))))).spectra[0];
  const E=F64(sp.be), y=F64(sp.y), bg={type:"smart",n:5};
  const s=C.suggestModel(E,y,{bg,el:"Ba",orb:"4d",lib:LIB});
  const rf=C.refineByResidual({E,y,c:sp.c,bg,comps:s.comps,cons:s.cons},{}); const P=rf.prob; const tot=P.comps.reduce((a,c)=>a+c.area,0);
  const inWin=P.comps.every(c=>c.be>=E[0]&&c.be<=E[E.length-1]), alive=P.comps.every(c=>c.area/tot>0.05);
  const anyBad=Object.values(P.cons).some(o=>Object.values(o).some(q=>q.mode==="link"&&!P.comps.some(c=>c.id===q.ref)));
  ok(inWin&&alive&&!anyBad&&P.comps.every(c=>c.doublet),`자동 제안 + 잔차 보강: 모든 성분이 창 안·넓이 몫 > 5 %·연결 대상 존재 (doublet ${P.comps.length} 개, χ²ν ${rf.stats.chi2nu.toFixed(2)}, ${rf.log.filter(t=>/채택|뺌/.test(t)).join("; ")})`,
    P.comps.map(c=>`${c.be.toFixed(2)} A${(100*c.area/tot).toFixed(1)}%`).join(" | ")); }

{ // removeComp: 연결 기준을 빼면 남은 것 중 가장 큰 성분이 기준과 그 제약을 물려받는다
  const a=C.newComp("GLS",80,{id:"a",area:1,doublet:{dE:2,r:0.5}}), b=C.newComp("GLS",85,{id:"b",area:50,doublet:{dE:2,r:0.5}}), c=C.newComp("GLS",87,{id:"c",area:30,doublet:{dE:2,r:0.5}});
  const P={comps:[a,b,c],cons:{a:{dE:{mode:"free",lo:1.5,hi:2.5}},b:{dE:{mode:"link",ref:"a",k:1,b:0},m:{mode:"link",ref:"a",k:1,b:0}},c:{dE:{mode:"link",ref:"a",k:1,b:0},m:{mode:"link",ref:"a",k:1,b:0}}}};
  C.removeComp(P,"a");
  ok(P.comps.length===2&&P.cons.b.dE.mode==="free"&&P.cons.b.dE.lo===1.5&&P.cons.c.dE.ref==="b"&&!P.cons.b.m&&P.cons.c.m.ref==="b","removeComp: 기준 성분을 빼면 가장 큰 성분이 Δ 의 자유 범위를 물려받고 나머지는 거기에 연결",JSON.stringify(P.cons)); }

console.log("== 7 원소 독립성 — Z 3–92 의 주요 궤도 전부(원소마다 센 것 2 개): 이름 해석 → 라이브러리 → 자동 제안 → 피팅 → 참값 회복");
{ const A=require("./all_elements.js"); const r=A.run({seed:31});
  ok(r.parse.length===0,`파일 이름 해석 ${r.n} 궤도 × 5 형식 ("Fe2p", "Fe 2p Scan", "fe2p", "[S-1] Fe2p.txt", "Fe_2p")`,r.parse.slice(0,3).join("; "));
  ok(r.rsf.length===0,`감도계수(Scofield Al · ALTHERMO1) ${r.n} 궤도 모두 있음`,r.rsf.slice(0,5).join(", ")+(r.mgGap.length?` · Mg Kα 표 빈칸 ${r.mgGap.length}: ${r.mgGap.join(", ")} (Mg 광원일 때 직접 입력)`:""));
  ok(r.dbl.length<=3&&r.dbl.every(x=>/^(Ho|Er|Tm) 4d$/.test(x)),"p·d·f 궤도의 doublet Δ (라이브러리) — 빠진 것은 다중항이 지배하는 Ho·Er·Tm 4d 뿐",r.dbl.join(", ")||"모두 있음");
  ok(r.sug.length===0,`단일 화학 상태 ${r.nSug} 경우(active·smart 배경): 자동 제안 → 피팅 → BE ±0.06 eV, 넓이 z<4.5(active)/±5 %(smart, 창 ≤20 eV)·±10 %(20–40 eV), Δ ±0.15 eV`,r.sug.slice(0,4).join(" | "));
  ok(r.two.length===0,`화학 상태 2 개 ${r.nTwo - r.nTwoFixed} 경우(분해되는 doublet·s 궤도): 흔든 출발 → 피팅 → BE·넓이 z < 4`,r.two.slice(0,4).join(" | "));
  ok(r.added.length===0,`잔차 보강(refineByResidual)이 참 1 성분 스펙트럼 ${r.nSug} 경우에 가짜 성분을 더하지 않음 (판단 배경: active Shirley + 기울기, z≥4 · ΔBIC≥10)`,r.added.slice(0,4).join(" | "));
  ok(r.auto2.length===0,`화학 상태 2 개 ${r.nAuto2} 경우를 자동 제안 + 잔차 보강 + 죽은 성분 덜기만으로: 성분 수 2 · BE ±0.15 eV · 넓이 ±10 %`,r.auto2.slice(0,4).join(" | "));
  ok(r.relaxed.length===0,`doublet 넓이 비 r·폭 배율 q 를 참이 통계 값인 ${r.nSug} 경우에 잘못 풀지 않음 (ΔBIC ≥ 10)`,r.relaxed.slice(0,4).join(" | "));
  ok(r.rpos.length<=Math.floor(0.05*r.nRpos),`양성 대조: 참 r = 1.2 × 통계 값인 두 상태 ${r.nRpos} 경우(분해되는 doublet) 중 95 % 이상에서 r 을 풀고 ±5 % 로 회복, 성분 수 2`,`허용 밖 ${r.rpos.length}: `+r.rpos.slice(0,4).join(" | "));
  console.log(`    정보: 넓은 창(>40 eV)의 smart 배경 허용 밖 ${r.sugWideSmart.length} 건(넓이 오차 최대 ≈14 %, 화면 경고) · 분해 안 되는 doublet(Δ<1.2 eV, Δ 는 범위 안에서 자유) 두 상태 ${r.nTwoFixed} 중 ${r.twoFixed.length} 건이 절대 허용 밖: ${r.twoFixed.join(" | ")} · ${r.times.toFixed(1)} s`); }

console.log("== 8 화학식 해석 · 화학량론 비교 · Auger parameter (SPEC §10)");
{ const P=s=>C.parseFormula(s,IO.SYMS); const eq=(a,b)=>JSON.stringify(Object.entries(a).sort())===JSON.stringify(Object.entries(b).sort());
  const okF=eq(P("NaLuF4"),{Na:1,Lu:1,F:4})&&eq(P("Na(Lu0.92Tm0.08)F4"),{Na:1,Lu:0.92,Tm:0.08,F:4})&&eq(P("Na3LuCl6"),{Na:3,Lu:1,Cl:6})
    &&eq(P("[Ba0.9Sr0.1]TiO3"),{Ba:0.9,Sr:0.1,Ti:1,O:3})&&eq(P("Ca10(PO4)6(OH)2"),{Ca:10,P:6,O:26,H:2})&&eq(P("CO"),{C:1,O:1})&&eq(P("Co3O4"),{Co:3,O:4});
  const bad=["Xx2","Na(LuF4","NaLu)F4","Na0F4","","Na-F","nalu"].filter(s=>{ try{ P(s); return true; }catch(e){ return false; } });
  ok(okF&&!bad.length,"화학식 해석: 괄호·대괄호·소수·중첩 괄호, 대소문자로 C+O vs Co 구분, 잘못된 기호·괄호·0 개·빈 값은 오류",bad.join(" | "));
  // 비교: NaLuF4 에서 Na 0.8 배, F 0.9 배인 측정 + 화학식 밖 C, 화학식에 있지만 측정 안 된 Tm
  const f=P("Na(Lu0.95Tm0.05)F4"); const rows=[{el:"Na",nrm:0.8,sNrm:0.02},{el:"Lu",nrm:0.95,sNrm:0.03},{el:"F",nrm:3.6,sNrm:0.05},{el:"C",nrm:2,sNrm:0.1}];
  const R=C.stoichCompare(rows,f,"Lu"); const it=e=>R.items.find(x=>x.el===e);
  const Wm=0.8+0.95+3.6, Nf=1+0.95+4;
  const okR=Math.abs(it("Na").R-(0.8/Wm)/(1/Nf))<1e-12&&Math.abs(it("F").ratioRel-0.9)<1e-12&&Math.abs(it("Na").ratioRel-0.8)<1e-12&&it("Lu").ratioRel===1
    &&R.outside.join()==="C"&&R.missing.join()==="Tm"&&R.ref==="Lu";
  // σ: 델타 방법 vs Monte Carlo(넓이 독립 정규 잡음 20000 회)
  const g=C.rng(7); const xs={Na:[],Lu:[],F:[]}, rr={Na:[],F:[]};
  for(let k=0;k<20000;k++){ const w={Na:0.8+0.02*g.g(),Lu:0.95+0.03*g.g(),F:3.6+0.05*g.g()}; const W=w.Na+w.Lu+w.F; ["Na","Lu","F"].forEach(e=>xs[e].push(w[e]/W)); rr.Na.push(w.Na/w.Lu/(1/0.95)); rr.F.push(w.F/w.Lu/(4/0.95)); }
  const sd=a=>{ const m=a.reduce((s,x)=>s+x,0)/a.length; return Math.sqrt(a.reduce((s,x)=>s+(x-m)*(x-m),0)/(a.length-1)); };
  const relErr=Math.max(...["Na","Lu","F"].map(e=>Math.abs(it(e).sx/sd(xs[e])-1)),Math.abs(it("Na").sRatioRel/sd(rr.Na)-1),Math.abs(it("F").sRatioRel/sd(rr.F)-1));
  ok(okR&&relErr<0.03,"화학량론 비교: 화학식 원소끼리 재정규화·측정/기대·기준 원소 비 정확, 화학식 밖(C)·측정 없음(Tm) 분리, σ(델타) vs Monte Carlo 20000 회 차 < 3 %",`최대 상대 차 ${(100*relErr).toFixed(2)} %`);
  const R2=C.stoichCompare(rows.concat([{el:"Na",nrm:0.7,sNrm:0.02}]),f,"Lu"); ok(R2.dup.join()==="Na","같은 원소 두 행이면 dup 으로 알림",JSON.stringify(R2.dup));
  // Auger parameter: 대전 이동 δ 와 무관해야 한다 (합성: Na 1s BE 1071.50, Na KLL KE 989.80 → α′ 2061.30)
  const hv=1486.68, mk=(c0,F)=>{ const E=Float64Array.from({length:401},(_,i)=>c0-10+i*0.05); return E; };
  let worst=0; const st={};
  for(const d of [0,1.73,4.37,-0.6]){
    // 중심을 격자 점에서 비킨다(+0.017 / −0.023 eV) — 포물선 다듬기의 부호가 틀리면 잡히게
    const Ep=mk(1071.5+d), Ea=mk(hv-989.8+d); const cP=1071.517+d, cA=hv-989.823+d;
    const pk=(E,c,F)=>{ const comp=C.newComp("GLS",c,{F,area:5000}); return C.compEval(comp,E,{emin:E[0],emax:E[E.length-1]}); };
    const yp=Float64Array.from(pk(Ep,cP,1.6),v=>v+800), ya=Float64Array.from(pk(Ea,cA,2.4),v=>v+1200);
    const a0=C.augerParam(C.peakMaxPos(Ep,yp,Ep[0],Ep[400]).be,C.peakMaxPos(Ea,ya,Ea[0],Ea[400]).be,hv).alpha; worst=Math.max(worst,Math.abs(a0-2061.34));
    const yn=y=>Float64Array.from(y,v=>{ const l=v*0.5; return Math.max(0,Math.round(l+Math.sqrt(l)*g.g()))/0.5; });
    const al=[]; for(let k=0;k<100;k++) al.push(C.augerParam(C.peakMaxPos(Ep,yn(yp),Ep[0],Ep[400]).be,C.peakMaxPos(Ea,yn(ya),Ea[0],Ea[400]).be,hv).alpha);
    const m=al.reduce((s,x)=>s+x,0)/al.length; st[d]={m,se:sd(al)/Math.sqrt(al.length),sd:sd(al)}; }
  const zs=Object.values(st).map(o=>Math.abs(o.m-2061.34)/o.se);
  ok(worst<0.005&&Math.max(...zs)<4,"Auger parameter α′ = BE_p + hν − BE_A(겉보기): 대전 이동 δ = 0·1.73·4.37·−0.6 eV 에서 같은 값 (잡음 없음 < 5 meV; Poisson 잡음 100 회 평균이 참값과 4 SE 안 — 치우침 없음)",
    `잡음 없음 ${(1000*worst).toFixed(2)} meV · 잡음 한 번의 sd ≈ ${(1000*st[0].sd).toFixed(0)} meV (최댓값 위치 방식) · |z| 최대 ${Math.max(...zs).toFixed(2)}`);
  const pm=C.peakMaxPos(mk(100),Float64Array.from({length:401},(_,i)=>i),95,105);
  const Ed=Float64Array.from({length:401},(_,i)=>110-i*0.05), yd=C.compEval(C.newComp("GLS",100.023,{F:1.5,area:5000}),Float64Array.from(Ed).reverse(),{emin:90,emax:110}).reverse();
  const pd=C.peakMaxPos(Ed,yd,90.5,100.0), pd2=C.peakMaxPos(Ed,yd,95,105);
  ok(pm===null&&pd===null&&pd2&&Math.abs(pd2.be-100.023)<0.005,"최댓값이 구간 끝이면 null — 오름·내림차순 에너지 모두, 내림차순에서도 위치 정확",`내림차순 ${pd2&&pd2.be.toFixed(4)}`);
  // '보이지 않는 주선' 규칙의 반례(검토): Au 4f 84 eV + 91 eV 작은 단일 봉우리(짝 94.7 은 창 밖) — 87 eV 에 가짜 doublet 을 만들면 안 됨
  { const Lb=LIB; const inf=C.doubletInfo(Lb,"Au","4f"); const E=Float64Array.from({length:301},(_,i)=>78+i*0.05);
    const cs=[C.newComp("GLS",84,{F:1.0,m:0.3,area:20000}),C.newComp("GLS",91,{F:1.0,m:0.3,area:1500})]; cs[0].doublet={dE:inf.dE,r:inf.r,q:1};
    const P=new Float64Array(E.length); cs.forEach(c=>{ const q=C.compEval(c,E,{emin:E[0],emax:E[E.length-1]}); for(let i=0;i<E.length;i++) P[i]+=q[i]; });
    const s=C.suggestModel(E,Float64Array.from(P,v=>v+500),{bg:{type:"smart",n:5},el:"Au",orb:"4f",lib:Lb});
    ok(s.comps.length===2&&s.comps.some(c=>Math.abs(c.be-84)<0.1&&c.doublet)&&s.comps.some(c=>Math.abs(c.be-91)<0.1),"자동 제안: 앞 doublet 의 짝·꼬리를 '보이지 않는 주선'으로 착각하지 않음 (Au 4f + 91 eV 단일 봉우리)",s.comps.map(c=>c.be.toFixed(2)+(c.doublet?"D":"")).join(" | ")); } }

console.log("== 9 모형 민감도 (SPEC §11) — 분석 선택(배경·끝점·창 끝)을 바꿔 다시 맞춘 흔들림");
{ const be=Float64Array.from({length:401},(_,i)=>80+i*0.05);
  const v1=C.sensVariants({lo:82,hi:96,bg:{type:"smart",n:5}},{be},{mode:"oat"}), v2=C.sensVariants({lo:82,hi:96,bg:{type:"smart",n:5}},{be},{mode:"full"});
  const v3=C.sensVariants({lo:80,hi:100,bg:{type:"active",n:5}},{be},{mode:"oat"});
  ok(v1.length===9&&v2.length===45&&v3.length===5&&v1[0].f.bg==="smart"&&v1[0].lo===82&&v1[0].hi===96&&v3.every(v=>v.lo>=80&&v.hi<=100),
    "변형 목록: 빠르게 9 · 꼼꼼히 45 (Smart, 창 82–96), 창이 스펙트럼 끝이면 바깥 변형을 빼고 중복 없이 5, 0 번 = 지금 설정",`${v1.length}/${v2.length}/${v3.length}`);
  // 주효과: 배경 종류만 값을 바꾸는 가짜 결과 → 끝점·창 효과 0 (배경이 바뀌면 끝점도 바뀌는 혼동이 섞이지 않음)
  const fake=v1.map((v,i)=>({i,ok:true,chi2nu:1,comps:[{id:"a",be:84,area:1,areaPct:v.f.bg==="active"?40:v.f.bg==="shirley"?35:30,fwhm:1}]}));
  const sf=C.sensSummary(v1,fake).comps[0].areaPct.eff;
  ok(Math.abs(sf.bg-10)<1e-9&&sf.n===0&&sf.lo===0&&sf.hi===0,"주효과: 배경 효과 10 %p, 끝점·창 끝 효과 0 (요인 혼동 없음)",JSON.stringify(sf));
  // 끝점 점 수만 값을 바꾸는 가짜 결과 → 배경 효과 0 (지금 배경의 n=3·10 변형이 배경 효과에 섞이지 않음)
  const fake2=v1.map((v,i)=>({i,ok:true,chi2nu:1,comps:[{id:"a",be:84,area:1,areaPct:v.f.n===3?40:30,fwhm:1}]}));
  const sf2=C.sensSummary(v1,fake2).comps[0].areaPct.eff;
  ok(sf2.bg===0&&Math.abs(sf2.n-10)<1e-9,"주효과: 끝점 점 수만 바꾸면 배경 효과 0 · 끝점 효과 10 %p",JSON.stringify(sf2));
  // 성분 판정: 지금 창 밖에 있던 성분(꼬리)은 그대로, 변형 창이 잘라 낸 성분만 건너뜀
  { const E=Float64Array.from({length:401},(_,i)=>80+i*0.05); const y=Float64Array.from(E,(e,i)=>1000+5000*Math.exp(-0.5*((e-88)/0.6)**2)+3000*Math.exp(-0.5*((e-90.2)/0.6)**2));
    const cA=C.newComp("GLS",88,{F:1.4,area:7000,id:"a"}), cB=C.newComp("GLS",90.2,{F:1.4,area:4000,id:"b"}), cT=C.newComp("GLS",81.0,{F:1.4,area:10,id:"t"});
    const m1={lo:82,hi:92,bg:{type:"smart",n:5},comps:[cA,cB,cT],cons:{}}; const r1=C.sensJob({be:E,y,c:null},m1,C.sensVariants(m1,{be:E},{mode:"oat"}));
    const m2={lo:84,hi:90.4,bg:{type:"smart",n:5},comps:[cA,cB],cons:{}}; const vs2=C.sensVariants(m2,{be:E},{mode:"oat"}); const r2=C.sensJob({be:E,y,c:null},m2,vs2);
    const cut=vs2.map((v,i)=>v.hi<90.2?i:-1).filter(i=>i>=0);
    ok(r1[0].ok&&cut.length>0&&cut.every(i=>!r2[i].ok)&&r2[0].ok,"성분 판정: 지금 창 밖의 꼬리 성분이 있어도 기준이 실패하지 않고, 변형 창이 잘라 낸 성분만 건너뜀",`기준 ${r1[0].ok} · 잘린 ${cut.length} 개`); }
  // 피팅 안 배경의 기울기를 '고정'한 사용자 모형: 0 번 변형이 그 제약 그대로(검토 반례 — 예전엔 기울기를 풀어 1.8σ 어긋남)
  { const R=C.rng(5); const lo=522,hi=540,st=0.05,n=Math.round((hi-lo)/st)+1; const E=Float64Array.from({length:n},(_,i)=>lo+i*st);
    const tr=[C.newComp("GLS",530.0,{F:1.3,m:0.3,area:20000,id:"a"}),C.newComp("GLS",531.4,{F:1.4,m:0.3,area:10000,id:"b"})];
    const ctx={emin:E[0],emax:E[n-1]}; const P=new Float64Array(n); tr.forEach(c=>{ const q=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=q[i]; });
    const Q=C.cumtrapz(E,P); const y=Float64Array.from(P,(p,i)=>{ const m=3000+25*(E[i]-lo)+0.08*Q[i]+p; return Math.max(0,Math.round(m*0.3+Math.sqrt(m*0.3)*R.g()))/0.3; });
    const spec={be:E,y,c:0.3}; const model={lo:524,hi:538,bg:{type:"active",n:5,cons:{s:{mode:"fixed"}},s:0},comps:tr.map(t=>C.clone(t)),cons:{}};
    const Pb=C.makeProblem(spec,model); const rb=C.fitOnce(Pb,{}); const sm=C.compSummary(Pb,rb,C.covariance(Pb,rb)); model.comps=Pb.comps;
    const vs=C.sensVariants(model,spec,{mode:"oat"}); const rs=C.sensJob(spec,model,vs); const b0=rs[0].comps.find(c=>c.id==="b").areaPct;
    ok(Math.abs(b0-sm[1].areaPct)<0.01*sm[1].sig.areaPct,"기울기 고정 피팅 안 배경: 0 번 변형 = 사용자 모형 (차이 < 0.01σ)",`차 ${(Math.abs(b0-sm[1].areaPct)/sm[1].sig.areaPct).toExponential(1)} σ`); }
  // 참값을 아는 합성(겹친 두 성분 2:1, 참 배경 = 피팅 안 Shirley + 기울기)에서 보고 구간이 참값을 포함하는 비율
  const covRun=(N,bg,seed)=>{ const R=C.rng(seed); let inS=0,inE=0,base0=0,dmax=0;
    for(let k=0;k<N;k++){ const lo=522,hi=540,st=0.05,n=Math.round((hi-lo)/st)+1; const E=Float64Array.from({length:n},(_,i)=>lo+i*st);
      const tr=[C.newComp("GLS",530.0,{F:1.3,m:0.3,area:20000,id:"a"}),C.newComp("GLS",531.4,{F:1.4,m:0.3,area:10000,id:"b"})];
      const ctx={emin:E[0],emax:E[n-1]}; const P=new Float64Array(n); tr.forEach(c=>{ const q=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=q[i]; });
      const Q=C.cumtrapz(E,P); const y=Float64Array.from(P,(p,i)=>{ const m=3000+25*(E[i]-lo)+0.08*Q[i]+p; return Math.max(0,Math.round(m*0.3+Math.sqrt(m*0.3)*R.g()))/0.3; });
      const spec={be:E,y,c:0.3}; const comps=tr.map(t=>{ const c=C.clone(t); c.be+=0.1*R.g(); c.area*=0.8; return c; });
      const model={lo:524,hi:538,bg:C.clone(bg),comps,cons:{a:{F:{mode:"free",lo:0.3,hi:4}},b:{F:{mode:"free",lo:0.3,hi:4},m:{mode:"link",ref:"a",k:1,b:0}}}};
      const Pb=C.makeProblem(spec,model); const rb=C.fitOnce(Pb,{}); const sm=C.compSummary(Pb,rb,C.covariance(Pb,rb)); const bA=sm[1].areaPct, sS=sm[1].sig.areaPct; model.comps=Pb.comps;
      const vs=C.sensVariants(model,spec,{mode:"oat"}); const res=C.sensJob(spec,model,vs); const s=C.sensSummary(vs,res).comps.find(c=>c.id==="b").areaPct;
      dmax=Math.max(dmax,Math.abs(s.base-bA)/sS); if(Math.abs(s.base-bA)<0.01*sS) base0++; const t=100/3; if(Math.abs(bA-t)<2*sS) inS++; if(t>=s.min-2*sS&&t<=s.max+2*sS) inE++; }
    return {inS,inE,base0,N,dmax}; };
  const w=covRun(40,{type:"smart",n:5},11), c=covRun(20,{type:"active",n:5,cons:{s:{mode:"free"}}},12);
  ok(w.base0===40&&c.base0===20,"0 번 변형(지금 설정)이 직접 피팅과 같은 최솟값 — 차이 < 0.01σ (60/60; 다시 맞춘 수치 오차만)",`${w.base0+c.base0}/60, 최대 ${Math.max(w.dmax,c.dmax).toExponential(1)} σ`);
  ok(w.inS/w.N<=0.3&&w.inE/w.N>=0.85,`기준 배경이 틀렸을 때(Smart, 참 = 피팅 안 Shirley + 기울기) 참 넓이 % 포함: 통계 ±2σ 만 ${w.inS}/${w.N} → 가능 범위 ${w.inE}/${w.N}`);
  ok(c.inS/c.N>=0.85&&c.inE/c.N>=0.85,`기준 배경이 맞을 때 참 넓이 % 포함: 통계 ±2σ ${c.inS}/${c.N} · 가능 범위 ${c.inE}/${c.N} (보수적이되 과하지 않음)`); }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
})().catch(e=>{ console.error(e); process.exit(2); });
