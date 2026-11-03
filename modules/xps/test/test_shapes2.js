/* CasaXPS Cookbook 2026 선모양 (SPEC_shapes2.md) 시험
   node test/test_shapes2.js   (modules/xps 에서; QUICK=1 이면 회복 시험을 줄임)
   1) 값·넓이·FWHM ↔ 독립 Python 기준 (test/ref/ref_shapes2.json: quad·FFT 합성곱·누적 Simpson·mpmath)
      허용: 값 ≤ 1e-5 (최댓값 대비), 넓이 ≤ 1e-6, FWHM ≤ 1e-5 (상대)
   2) 표기 왕복 · 오류 · 하위 호환(AGT GL(m)T(k) ≡ GLN tail T)
   3) 합성 잡음 스펙트럼 회복: 매개변수가 3σ 안 ≥ 90 % (모양마다 200 회; QUICK 40 회)
   4) 속도: 10 성분 × 2000 채널 평가·피팅 시간 · 4b TP wL 까지 풀었을 때 같은 프로세스 옛 LM 과 비교 + 구조적 갇힘 앙상블에서 멈춤 규칙(SPEC §3.3b)이 구해 냄
   4c) 멈춤 규칙이 최종 χ² 를 나쁘게 만들지 않음 — 무작위 합성 스펙트럼 고정 표본(규칙 끔·켬 비교) */
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

/* 1b (v1.9.12) 좁은 구조 vs 적응 구적 — 엔진과 따로 쓴 기준(Cookbook 식 그대로, 정규화·Φ·합성곱 모두 적응 Gauss–Kronrod): 허용 3e-7 (최댓값 대비)
   큰 α TLA (s0 = F 였던 v1.9.11: 1e-5–1e-4), 큰 γ 의 꺾인 ST (꺾임 보정 조각 4σ: 3.8e-4), α 20 LA 의 넓은 창 (층 격자 — SPEC_shapes2 §3) */
console.log("== 1b 좁은 구조 (큰 α TLA · 큰 γ 꺾인 ST · 층 격자) vs 적응 구적");
{ const PI=Math.PI, FW=1/(2*Math.sqrt(2*Math.LN2));
  const XG=[0.991455371120812639,0.949107912342758525,0.864864423359769073,0.741531185599394440,0.586087235467691130,0.405845151377397167,0.207784955007898468,0];
  const WK=[0.022935322010529225,0.063092092629978553,0.104790010322250184,0.140653259715525919,0.169004726639267903,0.190350578064785410,0.204432940075298892,0.209482141084727828], WG=[0.129484966168869693,0.279705391489276668,0.381830050505118945,0.417959183673469388];
  const gk=(f,a,b)=>{ const c=(a+b)/2,h=(b-a)/2,fc=f(c); let k=fc*WK[7],g=fc*WG[3]; for(let j=0;j<7;j++){ const d=h*XG[j],t=f(c-d)+f(c+d); k+=WK[j]*t; if(j&1) g+=WG[(j-1)>>1]*t; } return [k*h,Math.abs((k-g)*h)]; };
  const quad=(f,a,b,tol,brk)=>{ const P=[a].concat((brk||[]).filter(t=>t>a&&t<b).sort((u,v)=>u-v),[b]), L=[]; for(let i=0;i+1<P.length;i++){ const r=gk(f,P[i],P[i+1]); L.push([P[i],P[i+1],r[0],r[1]]); }
    for(let it=0;it<5000;it++){ let e=0,im=0,S=0; L.forEach((q,i)=>{ e+=q[3]; S+=q[2]; if(q[3]>L[im][3]) im=i; }); if(e<=tol*Math.max(1e-300,Math.abs(S))) break;
      const q=L[im],m=(q[0]+q[1])/2,r1=gk(f,q[0],m),r2=gk(f,m,q[1]); L.splice(im,1,[q[0],m,r1[0],r1[1]],[m,q[1],r2[0],r2[1]]); } return L.reduce((s,q)=>s+q[2],0); };
  const qinf=(f,a,sg)=>quad(t=>{ const u=1-t; return f(a+sg*t/u)/(u*u); },0,1,1e-14)*sg;   // ∫_a^{a+sg·∞}
  const erfc=z=>{ if(z<0) return 2-erfc(-z); if(z<3){ let s=0,t=z,k=0; const z2=z*z; for(;;){ const q=t/(2*k+1); s+=q; k++; t*=-z2/k; if(Math.abs(q)<1e-18*Math.abs(s)||k>300) break; } return 1-2/Math.sqrt(PI)*s; }
    let f=0; for(let k=60;k>=1;k--) f=k/2/(z+f); return Math.exp(-z*z)/Math.sqrt(PI)/(z+f); };
  const gau=(t,s)=>Math.exp(-0.5*(t/s)*(t/s))/(s*Math.sqrt(2*PI));
  const l6=(T,i0,h,z)=>{ const t=(z-i0)/h; let b=Math.floor(t)-2; b=Math.max(0,Math.min(T.length-6,b)); let s=0; for(let i=0;i<6;i++){ let w=1; for(let j=0;j<6;j++) if(j!==i) w*=(t-b-j)/(i-j); s+=w*T[b+i]; } return s; };
  const K15=(f,a,b)=>{ const c=(a+b)/2,h=(b-a)/2; let s=f(c)*WK[7]; for(let j=0;j<7;j++) s+=WK[j]*(f(c-h*XG[j])+f(c+h*XG[j])); return s*h; };
  // Cookbook 정의 (u = x/F = −x_c): TLA = (1+4u²)^−α·(π/2 + atan 2u + π/μ), LA = (1+4u²)^−(u ≥ 0 ? α : β)
  const CASES=[
    {t:"TLA α 20 μ 2 F 1 wG 0.05",sh:"TLA",p:{F:1,alpha:20,mu:2,wG:0.05}},
    {t:"TLA α 20 μ 0.3 F 0.67 wG 0.2",sh:"TLA",p:{F:0.67,alpha:20,mu:0.3,wG:0.2}},
    {t:"TLA α 5 μ 300 F 0.67 wG 0.01",sh:"TLA",p:{F:0.67,alpha:5,mu:300,wG:0.01}},
    {t:"SB TLA α 12 μ 2 F 1 wG 0.1",sh:"TLA",p:{F:1,alpha:12,mu:2,wG:0.1,mod:"SB"}},
    {t:"ST(1.5) TLA α 15 μ 1 F 0.8 wG 0.05",sh:"TLA",p:{F:0.8,alpha:15,mu:1,wG:0.05,mod:"ST",stMu:1.5,stG:0}},
    {t:"ST(3.09, γ 1.1) LA(3.5, 12.4) F 0.7 (꺾임 + 큰 γ)",sh:"LA",p:{F:0.7,alpha:3.5,beta:12.4,wG:0,mod:"ST",stMu:3.09,stG:1.1}},
    {t:"SB LA(20, 2.5) F 0.14 wG 0.09 · 창 ±13 eV (층 격자)",sh:"LA",p:{F:0.14,alpha:20,beta:2.5,wG:0.09,mod:"SB"},W:13},
    {t:"LA(15, 15) F 0.2 wG 0.05 · 창 ±10 eV (층 격자)",sh:"LA",p:{F:0.2,alpha:15,beta:15,wG:0.05},W:10}];
  let worst=0;
  for(const T of CASES){ const p=T.p, F=p.F, a=p.alpha, b=p.beta, mu=p.mu, s1=(p.wG||0)*FW, s2=p.mod==="ST"?(p.stG||0)*FW:0, muI=p.mod==="ST"?1/p.stMu:0;
    const raw=T.sh==="TLA"?x=>{ const u=x/F; return Math.pow(1+4*u*u,-a)*(PI/2+Math.atan(2*u)+PI/mu); }:x=>{ const u=x/F; return Math.pow(1+4*u*u,-(u>=0?a:b)); };
    const sn=F*Math.min(1,Math.sqrt(Math.pow(2,1/Math.max(a,b||0))-1)), br=[0]; [0.25,0.5,1,2,4,8,16,32].forEach(k=>br.push(-k*sn,k*sn));
    const Z=quad(raw,-40*sn,40*sn,1e-15,br)+qinf(raw,40*sn,sn)+qinf(t=>raw(-t),40*sn,sn), f=x=>raw(x)/Z;
    const W=T.W||Math.max(10*F,12*(s1+s2)), lo=-W, hi=W;
    let Phi=null; if(p.mod){ const R=W+12*(s1+s2)+sn, hc=sn/256, n=Math.ceil(R/hc)+4, P=new Float64Array(2*n+1); P[0]=qinf(t=>f(-t),n*hc,sn); for(let i=1;i<=2*n;i++) P[i]=P[i-1]+K15(f,(i-1-n)*hc,(i-n)*hc); Phi=z=>l6(P,-n*hc,hc,z); }
    const c1=(g,y)=>s1>0?quad(t=>gau(y-t,s1)*g(t),y-9.5*s1,y+9.5*s1,1e-15,br.concat([y])):g(y);
    const u=p.mod==="SB"?y=>c1(Phi,y):p.mod==="ST"?y=>c1(f,y)*(c1(Phi,y)+muI)/(0.5+muI):y=>c1(f,y);
    const val=s2>0?y=>quad(t=>gau(y-t,s2)*u(t),y-9.5*s2,y+9.5*s2,1e-14,br.concat([y])):u;
    const xs=[]; const c=6*sn+9*(s1+s2); for(let i=0;i<=100;i++) xs.push(-c+2*c*i/100); for(let i=0;i<=40;i++) xs.push(lo+(hi-lo)*i/40); xs.sort((q,r)=>q-r);
    let ref=xs.map(val);
    if(p.mod==="SB"){ const Pn=z=>0.5*erfc(-z/(s1*Math.SQRT2)); const N=quad(t=>Phi(t)*(Pn(hi-t)-Pn(lo-t)),lo-9.5*s1,hi+9.5*s1,1e-14,br.concat([lo,hi,lo+9.5*s1,hi-9.5*s1])); ref=ref.map(v=>v/N); }
    const cp=Object.assign({id:"t",shape:T.sh,be:0,area:1,doublet:null},p), v=C.shapeUnit(T.sh,F64(xs),C.shapeParams(cp,1),{lo,hi}); let mx=0,e=0; ref.forEach((r,i)=>{ mx=Math.max(mx,Math.abs(r)); }); ref.forEach((r,i)=>{ e=Math.max(e,Math.abs(v[i]-r)/mx); });
    worst=Math.max(worst,e); ok(e<=3e-7,`${T.t}`,`최댓값 대비 ${e1(e)} (허용 3e-7)`); }
  console.log(`     최대 ${e1(worst)}`); }

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
  /* 4b 퇴화 매개변수까지 푼 경우 (SPEC §3.3b · SPEC_shapes2 §3): TP 의 wL 도 자유. 판정은 저장한 기준값 없이 **같은 프로세스에서 옛 LM(규칙 끔 = XPS_NOSTALL)을 함께 돌려** 비교한다
     (v1.9.11: 예전엔 옛 LM χ² 를 숫자로 저장 — 이 문제의 LM 경로는 혼돈적이라 모형이 1e-8 만 바뀌어도 다시 재야 했다):
     ① 최종 χ² 새 ≤ 옛·(1 + 1e-7) (예전과 같은 부등식, 기준은 언제나 지금 모형의 옛 LM) ② 갈래를 하나도 풀지 않았으면 χ²·매개변수 비트 동일(예전: seed 4 만 저장값과 같음)
     ③ 새 ≤ 5 s (주 자료) / 10 s (seed). 옛 LM 을 시험 안에서 돌릴 수 있는 것은 v1.9.11 에서 TP 의 wL 문턱이 없어져 이 문제가 더는 퇴화가 아니기 때문
     (옛 LM 0.6–1.6 s; v1.9.10 모형에서는 wL < 0.0015·F 에서 기울기 0 → 옛 LM 이 2차 500 번 기어감, 27–52 s). 갈래가 실제로 구해 내는지는 아래 '구조적 갇힘 앙상블'에서 본다 */
  const oldNew=mk=>{ const Pa=mk(); C._hook.stall=false; const ta=Date.now(); let a; try{ a=C.fitOnce(Pa,{}); } finally{ C._hook.stall=true; } const tA=Date.now()-ta;
    const Pb=mk(), tb=Date.now(), b=C.fitOnce(Pb,{}), tB=Date.now()-tb; const pa=C.getVec(Pa,a.map), pb=C.getVec(Pb,b.map);
    return {a,b,tA,tB,rel:b.chi2/a.chi2-1,ident:a.chi2===b.chi2&&JSON.stringify(pa)===JSON.stringify(pb)}; };
  const verdict=(o,lim)=>o.rel<=1e-7&&o.tB<lim&&(o.b.trig||o.ident);
  const desc=o=>`χ² ${o.b.chi2.toFixed(2)} vs 옛 LM ${o.a.chi2.toFixed(2)} (같은 프로세스, ${o.a.iter1}+${o.a.iter} 반복, ${(o.tA/1000).toFixed(1)} s): ${o.rel>=0?"+":""}${o.rel.toExponential(1)} (허용 ≤ +1e-7)${o.b.trig?"":o.ident?" · 갈래 없음 → 비트 같음":" · 갈래 없는데 비트 다름"}`;
  { const cons2=C.clone(cons); cons2.c8.wL={mode:"free"};
    const mk=yy=>()=>C.makeProblem({be:E,y:yy,c:1},{lo:E[0],hi:E[1999],bg:{type:"active",c0:100,k:0,s:0,_auto:false,cons:{k:{mode:"fixed"},s:{mode:"fixed"}}},comps:C.clone(comps).map(c=>{ c.be+=0.02; return c; }),cons:cons2});
    const o=oldNew(mk(y)), r2=o.b, fl=((r2.stall||{}).flags||[]).map(q=>`${q.id}|${q.k}(${q.why.join("+")})`);
    const br=r2.branch?` · 갈래 ${r2.branch.kind}${r2.branch.reset.length?" ("+r2.branch.reset.map(q=>q.id+"|"+q.k).join(", ")+" 되돌림)":""}`:"";
    ok(verdict(o,5000),`TP wL 까지 자유: ${(o.tB/1000).toFixed(1)} s (≤ 5 s) · ${r2.iter1}+${r2.iter} 반복 · 정지 ${r2.stop}${br}${fl.length?" · 표시 "+fl.join(", "):""}`,desc(o));
    /* 같은 문제·다른 잡음(전체 모드만): 잡음 = C.rng(seed) 의 Poisson 근사 */
    if(!QUICK) [3,4,5,6].forEach(sd=>{ const Rs=C.rng(sd), ys=Float64Array.from(mu,v=>Math.max(0,Math.round(v+100+Math.sqrt(v+100)*Rs.g())));
      const q=oldNew(mk(ys)), rs=q.b;
      ok(verdict(q,10000),`  같은 문제·잡음 seed ${sd}: ${(q.tB/1000).toFixed(1)} s (≤ 10 s) · ${rs.iter1}+${rs.iter} 반복 · 갈래 ${rs.trig}${rs.branch?" ("+rs.branch.kind+")":""}`,desc(q)); }); }
  /* 4b 구조적 갇힘 앙상블 (전체 모드만): 폭이 한계에 붙어 LM 이 감쇠에 갇히는 무작위 문제(test/stall_gen.js seed 279 · 423 · 150 — rescue 로 χ² 21–93 % 낮아진 seed)를
     자료에 1e-9 상대 섭동(C.rng(9000+k), k = 0…7)을 줘 24 개로. 옛 LM 이 갇히는지 자체가 혼돈적이다(실측: k = 0…6 에서 옛 LM 이 갇혀 갈래가 구한 것
     v1.9.11 모형 5/21, v1.9.10 모형 4/21 — 나머지는 옛 LM 이 덫을 비켜 가 갈래가 필요 없음). 그래서 'seed X 는 구해져야 한다' 대신
     ① 모든 피팅: 새 ≤ 옛·(1 + 1e-7), 갈래 없으면 비트 동일, 새 ≤ 10 s ② 앙상블에서 갈래(rescue)가 채택되어 옛보다 1e-3 넘게 낮은 피팅 ≥ 1
     (모형이 바뀌어 경로가 다시 섞여도, 실측 비율 p ≈ 0.21 이면 24 개 모두 실패할 확률 0.79²⁴ ≈ 0.4 %; 지금 모형에서는 결정적) */
  if(!QUICK){ const G=require("./stall_gen.js"); let n=0, rescued=[], bad=[], tmax=0;
    for(const sd of [279,423,150]) for(let k=0;k<8;k++){ const g0=G.gen(C,sd), Rk=C.rng(9000+k), yk=k?Float64Array.from(g0.spec.y,v=>v*(1+1e-9*Rk.g())):g0.spec.y;
      const q=oldNew(()=>C.makeProblem(Object.assign({},g0.spec,{y:yk}),g0.st)); n++; tmax=Math.max(tmax,q.tB);
      if(!verdict(q,10000)) bad.push(`${sd}/${k} ${q.rel.toExponential(1)} ${q.tB} ms${q.b.trig||q.ident?"":" 비트 다름"}`);
      if(q.b.branch&&q.b.branch.kind==="rescue"&&q.rel<-1e-3) rescued.push(`${sd}/${k} ${q.rel.toExponential(1)} (${q.b.branch.reset.map(t=>t.id+"|"+t.k).join(",")})`); }
    ok(!bad.length&&rescued.length>=1,`  구조적 갇힘 앙상블 ${n} 피팅: 모두 새 ≤ 옛·(1+1e-7) · 갈래 없으면 비트 같음 · ≤ 10 s (최대 ${(tmax/1000).toFixed(1)} s) · rescue 로 구해진 것 ${rescued.length} (≥ 1)`,
      (rescued.join(" · ")||"없음")+(bad.length?" · 실패 "+bad.join(", "):"")); } }

/* 4c 멈춤 규칙이 최종 χ² 를 나쁘게 만들지 않음 (SPEC §3.3b): 무작위 합성 스펙트럼(test/stall_gen.js) 의 고정 표본을 같은 프로세스에서
   규칙 끔(옛 LM)·켬으로 피팅 → 최종 χ² 새 ≤ 옛·(1 + 1e-7), 갈래를 풀지 않은 피팅은 비트 단위로 같음. 표본 = 갈래가 생기는 seed(바닥에 못 닿아 현직이 이김 —
   rescue·kick, 1차·2차) + 갇혔다가 grace 안에 끝나는 seed + 갇히지 않는 seed. 3989·3132 = 현직이 1차를 끝까지 간 뒤 바닥의 갈래 채택 → 2차를 두 끝점에서(옛 2차 = 옛 최종). 1차 현직을 도중에 버리는 경우는 v1.9.10 모형의 §4b 잡음 seed 5 (v1.9.11 모형에서는 그 문제가 퇴화가 아니어서 갈래 없음 — 구조적 갇힘은 §4b 앙상블).
   전체 1100 개는 node test/stall_neverworse.js (2 코어 약 40 분) */
console.log("== 4c 멈춤 규칙 — 최종 χ² 가 옛 LM 보다 나빠지지 않음 (고정 표본)");
{ const G=require("./stall_gen.js");
  const seeds=QUICK?[3,4,5,8,38,722,1018,418,340,3989]:[3,4,5,8,38,722,1018,418,340,150,420,910,521,60,98,3989,3132];
  let worse=[], notId=[], n=0, nb=0, maxR=-Infinity, ta=0, tb=0;
  for(const s of seeds){ const g=G.gen(C,s); C._hook.stall=false; const a=G.run(C,g); C._hook.stall=true; const b=G.run(C,g); if(!a||!b||a.rec.err||b.rec.err) continue;
    n++; ta+=a.rec.ms; tb+=b.rec.ms; const rel=b.rec.chi2/a.rec.chi2-1; maxR=Math.max(maxR,rel); if(rel>1e-7) worse.push(`${s} ${rel.toExponential(2)}`);
    if(b.rec.trig) nb++; else if(b.rec.chi2!==a.rec.chi2||JSON.stringify(b.rec.p)!==JSON.stringify(a.rec.p)) notId.push(s); }
  ok(!worse.length&&!notId.length&&n===seeds.length,`${n} 피팅(갈래 ${nb}): 최종 χ² 새 ≤ 옛·(1+1e-7) 모두 · 갈래 없는 ${n-nb} 개 비트 같음`,
    `최대 새/옛 − 1 = ${maxR.toExponential(2)} · 시간 옛 ${(ta/1000).toFixed(1)} s → 새 ${(tb/1000).toFixed(1)} s${worse.length?" · 나빠짐 "+worse.join(", "):""}${notId.length?" · 비트 다름 "+notId.join(","):""}`); }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.map(s=>" - "+s).join("\n")); process.exit(1); }
