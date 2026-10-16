/* ============================================================
   INNO Analytics — Energy Diagram 화면
   ED_CORE(물리·데이터) 위에 SVG 도식·편집·툴팁을 얹는다.
   화면에 남기는 글은 출처와 함정뿐 — 조작법은 눌러 보면 알 수 있게.
   ============================================================ */
(function(){
"use strict";
const C = window.ED_CORE;
const KEY = "innoanalytics_ed_v1";
const esc = C.escHTML;
const $q = (r,s)=>r.querySelector(s);
const uid = () => Math.random().toString(36).slice(2,9);
const fmt = (v,d)=> v==null||!isFinite(v) ? "—" : Number(v).toLocaleString("en-US",{maximumFractionDigits:d==null?0:d,minimumFractionDigits:d==null?0:d});
const note = (msg,kind)=>{ try{ if(typeof toast==="function") toast(msg,kind); }catch(e){} };

/* ---------------- 상태 ---------------- */
const DEF = () => ({v:1,formula:"",cols:[],host:{key:null,name:"",w:360,src:"",est:true,T:295,user:false},
  lasers:[],lasersUser:false,procs:[],show:{stark:true,allLabels:false,emis:"main",gsa:true,nums:false,unit:"cm"},
  Emax:"auto",mode:"view",sel:null,warn:[],layers:[]});
let S = load();
function load(){
  try{ const s=JSON.parse(localStorage.getItem(KEY)||"null"); if(s&&s.v===1){ const o=Object.assign(DEF(),s,{mode:"view",sel:null}); /* 예전 판에서 루프를 여러 번 눌러 쌓인 똑같은 과정은 하나만 남긴다 */ const u=[]; (o.procs||[]).forEach(p=>{ if(!u.some(q=>sameProc(q,p))) u.push(p); }); o.procs=u; return o; } }catch(e){}
  return DEF();
}
function save(){ try{ const o=Object.assign({},S); delete o.sel; localStorage.setItem(KEY,JSON.stringify(o)); }catch(e){} }

/* 편집 단계(메모리) */
let PICK = null;          /* {kind:'et'|'ph'|'em', steps:[{c,id}]} */
let ROOT = null, TIP = null, MODEL = null;

/* ---------------- 프리셋 (문헌 표준 메커니즘만) ---------------- */
const PRESETS = [
  {name:"β-NaYF₄:20%Yb,2%Er — 980 nm upconversion",f:"β-NaYF4:20%Yb,2%Er",l:[980],
   p:[["PH",0,"2F7/2","2F5/2"],["ET",[0,"2F5/2","2F7/2"],[1,"4I15/2","4I11/2"]],["ET",[0,"2F5/2","2F7/2"],[1,"4I11/2","4F7/2"]],
      ["MPR",1,"4F7/2","4S3/2"],["MPR",1,"4I11/2","4I13/2"],["ET",[0,"2F5/2","2F7/2"],[1,"4I13/2","4F9/2"]]],
   ref:"Auzel, Chem. Rev. 104, 139 (2004); Haase & Schäfer, Angew. Chem. Int. Ed. 50, 5808 (2011)"},
  {name:"β-NaYF₄:25%Yb,0.3%Tm — 980 nm upconversion",f:"β-NaYF4:25%Yb,0.3%Tm",l:[980],
   p:[["PH",0,"2F7/2","2F5/2"],["ET",[0,"2F5/2","2F7/2"],[1,"3H6","3H5"]],["MPR",1,"3H5","3F4"],
      ["ET",[0,"2F5/2","2F7/2"],[1,"3F4","3F2"]],["MPR",1,"3F2","3H4"],["ET",[0,"2F5/2","2F7/2"],[1,"3H4","1G4"]],
      ["ET",[1,"3F2","3H6"],[1,"3H4","1D2"]],["ET",[0,"2F5/2","2F7/2"],[1,"1D2","3P2"]],["MPR",1,"3P2","1I6"]],
   ref:"Auzel (2004); Wang & Liu, Chem. Soc. Rev. 38, 976 (2009). ¹D₂ 는 Tm–Tm 과정 ³F₂+³H₄→³H₆+¹D₂ 로 채워진다는 해석이 일반적 — 문헌은 흔히 ‘cross-relaxation’이라 부르지만 acceptor 가 여기 상태(³H₄)에서 올라가므로 이 도구의 엄밀한 분류로는 ETU"},
  {name:"β-NaYF₄:8%Tm — 1064 nm photon avalanche(PA)",f:"β-NaYF4:8%Tm",l:[1064],
   p:[["PH",0,"3H6","3H5"],["MPR",0,"3H5","3F4"],["PH",0,"3F4","3F2"],["MPR",0,"3F2","3H4"],["ET",[0,"3H4","3F4"],[0,"3H6","3F4"]]],
   ref:"Lee et al., Nature 589, 230 (2021) — 3준위 모델(³H₆, ³F₄, ³H₄), 약한 non-resonant GSA + resonant ESA + cross-relaxation 루프"},
  {name:"β-NaYF₄:20%Yb,2%Ho — 980 nm",f:"β-NaYF4:20%Yb,2%Ho",l:[980],
   p:[["PH",0,"2F7/2","2F5/2"],["ET",[0,"2F5/2","2F7/2"],[1,"5I8","5I6"]],["MPR",1,"5I6","5I7"],
      ["ET",[0,"2F5/2","2F7/2"],[1,"5I6","5F4"]],["ET",[0,"2F5/2","2F7/2"],[1,"5I7","5F5"]]],
   ref:"Auzel (2004); Yb→Ho ⁵I₆ 전달은 포논 보조(불일치 ≈1,700 cm⁻¹)"},
  {name:"β-NaYF₄:Nd,Yb — 808 nm (composition only)",f:"β-NaYF4:5%Nd,10%Yb",l:[808],p:[["PH",0,"4I9/2","4F5/2"],["MPR",0,"4F5/2","4F3/2"]],ref:""},
  {name:"β-NaGdF₄:Yb,Tm@NaGdF₄:Tb — Gd migration (composition only)",f:"β-NaGdF4:49%Yb,1%Tm@NaGdF4:15%Tb",l:[980],p:[],ref:"Wang et al., Nat. Mater. 10, 968 (2011)"},
  {name:"Y₂O₃:5%Eu³⁺ (composition only)",f:"Y2O3:5%Eu",l:[394],p:[],ref:""},
  {name:"CaF₂:Eu²⁺ (composition only)",f:"CaF2:Eu2+",l:[365],p:[],ref:""},
  {name:"K₂SiF₆:Mn⁴⁺ (composition only)",f:"K2SiF6:Mn4+",l:[455],p:[],ref:""},
  {name:"Al₂O₃:Cr³⁺ ruby (composition only)",f:"Al2O3:Cr3+",l:[532],p:[],ref:""}
];

/* ---------------- 화학식 적용 ---------------- */
function applyFormula(str,opt){
  opt=opt||{};
  const r=C.parseFormula(str);
  S.formula=r.input; S.warn=r.warn.slice(); S.layers=r.layers.map(l=>({tag:l.tag,hostKey:l.hostKey,phase:l.phase,host:l.hostInfo?l.hostInfo.name:null}));
  const old=S.cols; const cols=[];
  r.ions.forEach(ion=>{
    if(!ion.known){
      if(ion.inactive===false) S.warn.push(`${C.ionLabel(ion.key)}: 준위 데이터 없음 (지원하지 않는 이온/산화수)`);
      return;
    }
    const prev=old.find(c=>c.ion===ion.key);
    cols.push({uid:prev?prev.uid:uid(),ion:ion.key,conc:ion.layers.map(x=>x.conc),layers:ion.layers.map(x=>x.role),
      fromHost:ion.layers.some(x=>x.fromHost),notes:ion.notes,p:prev?prev.p:null,band:prev?prev.band:null,qSrc:ion.layers[0].qSrc});
  });
  /* 비활성 이온 안내는 생략(La, Lu, Y 등은 호스트) */
  S.cols=cols; if(!opt.keepRef) S.presetRef="";
  const core=r.layers[0];
  if(core && !S.host.user){
    if(core.hostInfo && core.hostInfo.w){ S.host=Object.assign(S.host,{key:core.hostKey,name:core.hostInfo.name+(core.phase==="α"?" (α, 입방)":""),w:core.hostInfo.w,src:core.hostInfo.src,est:!!core.hostInfo.est}); }
    else { const a=core.anion&&C.ANION_PH[core.anion]; S.host=Object.assign(S.host,{key:core.hostKey,name:core.hostKey||"",w:a?a[0]:400,src:a?a[1]:"음이온 정보 없음 — 값을 입력하세요 (추정)",est:true}); }
  }
  if(core && core.phase==="α") S.warn.push("α(입방)상 — 준위값은 LaF₃ 기준이고, α-상은 β-상보다 자리 대칭이 높아 방출 세기비가 다를 수 있음");
  if(!S.lasersUser || opt.resetLasers) S.lasers=C.defaultLasers(cols.map(c=>c.ion));
  /* 사라진 열을 가리키는 과정은 버린다 */
  const ids=new Set(cols.map(c=>c.uid));
  S.procs=S.procs.filter(p=>ids.has(p.a?p.a.c:p.c) && (!p.d||ids.has(p.d.c)));
  S.sel=null; PICK=null;
  save();
}
function addIonColumn(key){
  if(!C.IONS[key]) return;
  S.cols.push({uid:uid(),ion:key,conc:[null],layers:["-"],fromHost:false,notes:[],p:null,band:null,qSrc:"주기율표"});
  if(!S.lasersUser && !S.lasers.length) S.lasers=C.defaultLasers(S.cols.map(c=>c.ion));
  save(); rerender();
}
function loadPreset(i){
  const P=PRESETS[i]; if(!P) return;
  S.lasersUser=false; S.host.user=false; S.procs=[];
  applyFormula(P.f,{resetLasers:true});
  S.lasers=P.l.slice();
  P.p.forEach(x=>{
    if(x[0]==="PH"){ const c=S.cols[x[1]]; if(c) S.procs.push({uid:uid(),kind:"PH",c:c.uid,f:x[2],t:x[3],nm:null}); }
    else if(x[0]==="MPR"){ const c=S.cols[x[1]]; if(c) S.procs.push({uid:uid(),kind:"MPR",c:c.uid,f:x[2],t:x[3]}); }
    else if(x[0]==="ET"){ const cd=S.cols[x[1][0]], ca=S.cols[x[2][0]]; if(cd&&ca) S.procs.push({uid:uid(),kind:"ET",d:{c:cd.uid,f:x[1][1],t:x[1][2]},a:{c:ca.uid,f:x[2][1],t:x[2][2]},coop:false,type:null}); }
  });
  S.presetRef=P.ref||"";
  save(); rerender();
}

/* ---------------- 모델 계산 ---------------- */
function ctx(){ return {hw:+S.host.w||null,T:+S.host.T||295}; }
function buildModel(){
  const cx=ctx();
  const cols=S.cols.map(c=>{
    const lv=C.annotate(C.levelsFor(c),cx.hw,(C.IONS[c.ion]||{}).kind);
    const byId={}; lv.forEach(l=>byId[l.id]=l);
    return {c,def:C.IONS[c.ion],lv,byId};
  });
  const byUid={}; cols.forEach((m,i)=>{ m.i=i; byUid[m.c.uid]=m; });
  /* 과정 해석 */
  const procs=[];
  S.procs.forEach((p,k)=>{
    const r=resolveProc(p,byUid,cx); if(r) { r.n=k+1; procs.push(r); }
  });
  /* 자동 광흡수 후보 */
  const photons=[];
  const userPH=new Set(procs.filter(p=>p.kind==="PH"&&!p.broken).map(p=>p.m.c.uid+"|"+p.F.id+"|"+p.T.id));
  const extraOf=m=>{ const a=[]; procs.forEach(p=>{ if(p.broken) return; if(p.kind==="ET"){ if(p.ma===m) a.push(p.At.id); if(p.md===m) a.push(p.Dt.id); } else if(p.m===m) a.push(p.T.id); }); return a; };
  if(S.show.gsa) cols.forEach(m=>C.photonCandidates(m.c,m.lv,S.lasers,cx,extraOf(m)).forEach(x=>{ if(!userPH.has(m.c.uid+"|"+x.from.id+"|"+x.to.id)) photons.push(Object.assign({m},x)); }));
  /* population이 쌓이는 준위: 과정·자동 광흡수의 도착 준위와 그 multiphonon 연쇄 끝 */
  const pop=new Set(); const addP=(m,l)=>{ if(!l) return; (l.relax||[l.id]).forEach(id=>pop.add(m.c.uid+"|"+id)); };
  procs.forEach(p=>{ if(p.broken) return; if(p.kind==="ET"){ addP(p.ma,p.At); addP(p.md,p.Dt); } else addP(p.m,p.T); });
  photons.forEach(x=>{ if(x.c.grade<=2) addP(x.m,x.to); });
  /* 대표 방출 */
  const emis=[];
  cols.forEach(m=>{
    if(S.show.emis==="none") return;
    if(S.show.emis==="all"){
      m.lv.forEach(u=>{ if(u.ground||!(u.regime==="rad")) return;
        m.lv.forEach(l=>{ if(l.E>=u.E-300) return; const nm=C.cm2nm(u.E-l.E); if(nm<200||nm>2600) return; emis.push({m,up:u,lo:l,nm,src:"계산: 복사 우세 준위(p≥"+C.P_RAD+")에서 아래 모든 준위"}); }); });
      return;
    }
    let list=[];
    if(m.def.kind==="d") list=C.dMain(m.c,m.lv).filter(x=>x.kind==="em").map(x=>({up:[x.up],lo:x.lo,src:x.note}));
    else list=C.curatedFor(m.c.ion).filter(x=>x.kind==="em"||x.kind==="both");
    list.forEach(x=>{
      const u=m.byId[x.up[0]]; let l=m.byId[x.lo];
      const colPop=[...pop].some(k=>k.startsWith(m.c.uid+"|"));
      if(S.show.emis==="main" && colPop && !x.up.some(id=>pop.has(m.c.uid+"|"+id))) return;
      if(x.lo==="5d") return;
      if(!u||!l) return;
      const E=u.kind==="band"? (u.E-(u.stokes||0))-l.E : u.E-l.E;
      if(!(E>0) || C.cm2nm(E)>3000) return;     /* mid-IR(예: Ce³⁺ ²F₇/₂→²F₅/₂)는 도식에서 제외 */
      emis.push({m,up:u,lo:l,nm:C.cm2nm(E),typ:x.nm,src:x.src,pair:x.up.length>1?x.up:null});
    });
  });
  /* PA 루프 */
  const pa=[]; cols.forEach(m=>{ if(m.def.kind!=="ln3") return; C.paLoops(m.i,m.lv,S.lasers,cx).slice(0,1).forEach(r=>pa.push(Object.assign({m},r))); });
  /* 에너지 축 상한 */
  let top=0;
  const bump=E=>{ if(E>top) top=E; };
  procs.forEach(p=>p.levels.forEach(l=>bump(l.hi)));
  emis.forEach(e=>bump(e.up.hi)); photons.forEach(x=>bump(x.to.hi));
  if(!top) cols.forEach(m=>{ const g=m.lv.filter(l=>l.regime==="rad"&&!l.ground); if(g.length) bump(g[Math.min(g.length-1,2)].E); });
  let Emax = S.Emax==="auto" ? Math.min(52000,Math.max(12000,Math.ceil((top*1.06+1500)/1000)*1000)) : +S.Emax;
  MODEL={cols,byUid,procs,photons,emis,pa,Emax,cx};
  return MODEL;
}
function resolveProc(p,byUid,cx){
  if(p.kind==="ET"){
    const md=byUid[p.d.c], ma=byUid[p.a.c]; if(!md||!ma) return null;
    const Df=md.byId[p.d.f], Dt=md.byId[p.d.t], Af=ma.byId[p.a.f], At=ma.byId[p.a.t];
    if(!Df||!Dt||!Af||!At) return {p,broken:true,levels:[]};
    const c=C.classifyET({d:{ion:md.c.ion,from:Df,to:Dt},a:{ion:ma.c.ion,from:Af,to:At},coop:p.coop},cx);
    return {p,kind:"ET",type:p.type||c.type,auto:c.type,c,md,ma,Df,Dt,Af,At,levels:[Df,Dt,Af,At]};
  }
  const m=byUid[p.c]; if(!m) return null;
  const F=m.byId[p.f], T=m.byId[p.t]; if(!F||!T) return {p,broken:true,levels:[]};
  if(p.kind==="PH"){
    /* 여기 파장: 지정 없으면 활성 광원 중 가장 잘 맞는 것 */
    let best=null;
    const cands=p.nm?[p.nm]:S.lasers;
    cands.forEach(nm=>{ const c=C.classifyPhoton(F,T,C.nm2cm(nm),cx); if(!best||c.delta<best.c.delta) best={nm,c}; });
    const need=C.cm2nm(T.E-F.E);
    return {p,kind:"PH",type:F.ground?"GSA":"ESA",m,F,T,nm:best?best.nm:null,c:best?best.c:null,need,levels:[F,T]};
  }
  if(p.kind==="RAD"||p.kind==="MPR"){
    return {p,kind:p.kind,type:p.kind,m,F,T,nm:C.cm2nm(F.E-T.E),levels:[F,T]};
  }
  return null;
}

/* ---------------- 색 ---------------- */
function isDark(){ return document.documentElement.getAttribute("data-theme")==="dark"; }
function palette(forceLight){
  const cs=getComputedStyle(document.documentElement);
  const d=forceLight?false:isDark();
  const v=(n,f)=>forceLight?f:((cs.getPropertyValue(n)||"").trim()||f);
  return {dark:d,bg:v("--surface",d?"#232f2c":"#fff"),ink:v("--ink",d?"#e3ebe2":"#243433"),ink2:v("--ink-2",d?"#c3cfc5":"#4a5a56"),
    muted:v("--muted",d?"#a3b0a6":"#62706b"),faint:v("--faint",d?"#7f8c84":"#8b978f"),grid:v("--grid",d?"#33403b":"#e4e8e2"),
    line:v("--border",d?"#3a4841":"#dde3dc"),accent:v("--accent",d?"#d9a2c9":"#7a3f6e"),
    T:d?{ET:"#5cc9bf",ETU:"#c69be6",CR:"#f0a35e",EM:"#a9b6ad",BET:"#ef8c86",CET:"#8fb3e6",MPR:"#9aa79f",RAD:"#c3cfc5"}
       :{ET:"#0f7c74",ETU:"#7b3fa0",CR:"#b4580f",EM:"#6a7872",BET:"#b23a3a",CET:"#2f5f9e",MPR:"#7b8983",RAD:"#4a5a56"}};
}
const TYPE_SHORT={ET:"ET",ETU:"ETU",CR:"CR",EM:"EM",BET:"BET",CET:"CET",GSA:"GSA",ESA:"ESA",RAD:"RAD",MPR:"MPR"};
/* ed_core.js 가 한국어로 돌려주는 화면 라벨(등급 이름·과정 종류 이름·색 이름)을 화면용 영어로 바꾼다 — core 는 건드리지 않는다 */
const TYPE_NAME={ET:"Energy transfer (ET)",ETU:"Energy transfer upconversion (ETU)",CR:"Cross-relaxation (CR)",EM:"Energy migration (EM)",BET:"Back energy transfer (BET)",CET:"Cooperative energy transfer (CET)",GSA:"Ground-state absorption (GSA)",ESA:"Excited-state absorption (ESA)",RAD:"Radiative transition (emission)",MPR:"Multiphonon relaxation (MPR)"};
const typeName=t=>TYPE_NAME[t]||t;
const gradeName=c=>c.edge?"resonant (Stark edge)":(C.GRADE_EN[c.grade]||c.gradeName);
const BAND_EN={"보라":"Violet","파랑":"Blue","초록":"Green","노랑":"Yellow","주황":"Orange","빨강":"Red","근적외(NIR-I/II)":"NIR (NIR-I/II)","단파 적외(SWIR)":"SWIR"};
const bandEn=nm=>{ const b=C.bandName(nm); return BAND_EN[b]||b; };

/* ---------------- 렌더 진입 ---------------- */
function render(root){
  ROOT=root;
  ensureTip();
  root.innerHTML="";
  const wrap=document.createElement("div"); wrap.className="ed";
  wrap.innerHTML=`
  <div class="card ed-top">
    <h3>Composition <span class="hint">화학식(코어@셸) 또는 주기율표에서 이온 선택</span></h3>
    <div class="ed-formrow">
      <input type="text" id="edF" spellcheck="false" autocomplete="off" aria-label="Formula" placeholder="e.g. β-NaYF4:20%Yb,2%Er@NaYF4 · NaY0.78Yb0.2Er0.02F4 · K2SiF6:Mn4+" value="${esc(S.formula)}">
      <button class="btn primary" id="edApply">Apply</button>
      <select class="ed-sel" id="edPreset" aria-label="Preset"><option value="">Preset…</option>${PRESETS.map((p,i)=>`<option value="${i}">${esc(p.name)}</option>`).join("")}</select>
      <button class="btn" id="edPTbtn" aria-expanded="false" aria-controls="edPTbox">Periodic table</button>
    </div>
    <div id="edPTbox" hidden></div>
    <div class="ed-chips" id="edChips"></div>
    <div class="ed-host" id="edHost"></div>
    <div class="ed-warn" id="edWarn"></div>
  </div>
  <div class="ed-main">
    <div class="card ed-diag">
      <div class="ed-tools">
        <div class="seg" id="edMode" role="tablist" aria-label="Edit mode">
          <button data-m="view">View</button><button data-m="et">Energy transfer</button><button data-m="ph">Absorption</button><button data-m="em">Emission &amp; relaxation</button>
        </div>
        <span class="sp"></span>
        <label class="ed-tog"><input type="checkbox" id="edStark"> Stark width</label>
        <label class="ed-tog"><input type="checkbox" id="edAll"> All labels</label>
        <label class="ed-tog"><input type="checkbox" id="edNums"> Values</label>
        <select class="ed-sel" id="edEmis" aria-label="Emission display" title="대표 방출: 문헌 대표 전이 중 지금 과정으로 population이 쌓이는 준위에서 나오는 것만"><option value="main">Emission: populated levels</option><option value="lit">Emission: all literature-typical</option><option value="all">Emission: all possible</option><option value="none">Hide emission</option></select>
        <select class="ed-sel" id="edEmax" aria-label="Energy axis maximum"><option value="auto">Axis: auto</option>${[15000,20000,25000,30000,40000,52000].map(v=>`<option value="${v}">≤ ${v/1000}k cm⁻¹</option>`).join("")}</select>
        <select class="ed-sel" id="edUnit" aria-label="Axis unit"><option value="cm">cm⁻¹</option><option value="eV">eV</option></select>
        <select class="ed-sel" id="edExport" aria-label="Export"><option value="">Export…</option><option value="svg">SVG (diagram)</option><option value="png">PNG (diagram, 3×)</option><option value="csv">CSV (process table)</option><option value="json">JSON (save diagram)</option><option value="import">Import JSON…</option></select>
        <input type="file" id="edImport" accept=".json,application/json" hidden>
      </div>
      <div class="ed-prompt" id="edPrompt" aria-live="polite"></div>
      <div class="ed-svgwrap" id="edWrap"></div>
      <div class="ed-legend" id="edLegend"></div>
    </div>
    <div class="ed-side">
      <div class="card" id="edLasCard"></div>
      <div class="card" id="edProcCard"></div>
      <div class="card ed-det" id="edDetCard"></div>
      <div class="card" id="edIntCard" hidden></div>
      <div class="card" id="edSugCard"></div>
      <div class="card" id="edSrcCard"></div>
    </div>
  </div>`;
  root.appendChild(wrap);
  bindTop(wrap);
  rerender();
}
function rerender(){
  if(!ROOT||!ROOT.isConnected) return;
  buildModel();
  renderChips(); renderHost(); renderWarn(); renderTools(); drawSVG(); renderLegend();
  renderLasers(); renderProcs(); renderDetail(); renderInter(); renderSug(); renderSrc(); renderPrompt();
}

/* ---------------- 상단 ---------------- */
function bindTop(w){
  const F=$q(w,"#edF");
  $q(w,"#edApply").onclick=()=>{ S.host.user=false; applyFormula(F.value); rerender(); };
  F.addEventListener("keydown",e=>{ if(e.key==="Enter"){ e.preventDefault(); S.host.user=false; applyFormula(F.value); rerender(); } });
  $q(w,"#edPreset").onchange=e=>{ const v=e.target.value; e.target.value=""; if(v!=="") { loadPreset(+v); const f=$q(ROOT,"#edF"); if(f) f.value=S.formula; } };
  const pb=$q(w,"#edPTbtn"), box=$q(w,"#edPTbox");
  pb.onclick=()=>{ const open=box.hidden; box.hidden=!open; pb.setAttribute("aria-expanded",String(open)); if(open) renderPT(box); };
  $q(w,"#edMode").addEventListener("click",e=>{ const b=e.target.closest("button[data-m]"); if(!b) return; setMode(b.dataset.m); });
  $q(w,"#edStark").onchange=e=>{ S.show.stark=e.target.checked; save(); rerender(); };
  $q(w,"#edAll").onchange=e=>{ S.show.allLabels=e.target.checked; save(); rerender(); };
  $q(w,"#edNums").onchange=e=>{ S.show.nums=e.target.checked; save(); rerender(); };
  $q(w,"#edEmis").onchange=e=>{ S.show.emis=e.target.value; save(); rerender(); };
  $q(w,"#edEmax").onchange=e=>{ S.Emax=e.target.value==="auto"?"auto":+e.target.value; save(); rerender(); };
  $q(w,"#edUnit").onchange=e=>{ S.show.unit=e.target.value; save(); rerender(); };
  $q(w,"#edExport").onchange=e=>{ const v=e.target.value; e.target.value=""; if(v==="svg") exportSVG(); else if(v==="png") exportPNG(); else if(v==="csv") exportCSV(); else if(v==="json") exportJSON(); else if(v==="import") $q(ROOT,"#edImport").click(); };
  $q(w,"#edImport").onchange=e=>{ const f=e.target.files&&e.target.files[0]; if(f) importJSON(f); e.target.value=""; };
  w.addEventListener("keydown",e=>{
    if(e.key==="Escape" && PICK && PICK.steps.length){ e.stopPropagation(); PICK.steps=[]; rerender(); }
  },true);
}
function setMode(m){ S.mode=m; PICK= m==="view"?null:{kind:m,steps:[]}; save(); rerender(); }

const PT=[["H",1,1],["He",1,18],["Li",2,1],["Be",2,2],["B",2,13],["C",2,14],["N",2,15],["O",2,16],["F",2,17],["Ne",2,18],
 ["Na",3,1],["Mg",3,2],["Al",3,13],["Si",3,14],["P",3,15],["S",3,16],["Cl",3,17],["Ar",3,18],
 ["K",4,1],["Ca",4,2],["Sc",4,3],["Ti",4,4],["V",4,5],["Cr",4,6],["Mn",4,7],["Fe",4,8],["Co",4,9],["Ni",4,10],["Cu",4,11],["Zn",4,12],["Ga",4,13],["Ge",4,14],["As",4,15],["Se",4,16],["Br",4,17],["Kr",4,18],
 ["Rb",5,1],["Sr",5,2],["Y",5,3],["Zr",5,4],["Nb",5,5],["Mo",5,6],["Tc",5,7],["Ru",5,8],["Rh",5,9],["Pd",5,10],["Ag",5,11],["Cd",5,12],["In",5,13],["Sn",5,14],["Sb",5,15],["Te",5,16],["I",5,17],["Xe",5,18],
 ["Cs",6,1],["Ba",6,2],["Lu",6,3],["Hf",6,4],["Ta",6,5],["W",6,6],["Re",6,7],["Os",6,8],["Ir",6,9],["Pt",6,10],["Au",6,11],["Hg",6,12],["Tl",6,13],["Pb",6,14],["Bi",6,15],["Po",6,16],["At",6,17],["Rn",6,18]];
const LNROW=["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"];
const HOSTEL=new Set(["Li","Na","K","Rb","Cs","Mg","Ca","Sr","Ba","Sc","Y","La","Lu","Al","Ga","In","Si","Ge","Ti","Zr","Zn","Pb","Bi","B","P","V","Nb","Ta","Mo","W","F","Cl","Br","I","O","S","N","Se"]);
function renderPT(box){
  const Z={}; let z=0; const ALL=ELORDER(); ALL.forEach(e=>Z[e]=++z);
  const cell=el=>{ const act=!!C.OX[el]; const cls=act?"act":(HOSTEL.has(el)?"host":"");
    return `<button type="button" data-el="${el}" class="${cls}" title="${act?"Optically active — pick oxidation state":(HOSTEL.has(el)?"Host element":"No level data")}"><span>${el}</span><small>${Z[el]||""}</small></button>`; };
  let h='<div class="ed-pt" role="group" aria-label="Periodic table">';
  for(let r=1;r<=6;r++){ for(let c=1;c<=18;c++){ const e=PT.find(x=>x[1]===r&&x[2]===c); h+= e?cell(e[0]):'<span class="gap"></span>'; } }
  h+='<span class="gap"></span><span class="gap"></span><div class="lnrow" style="display:grid;grid-template-columns:repeat(15,minmax(0,1fr));gap:inherit;margin-top:6px">'+LNROW.map(cell).join("")+'<span></span></div>';
  h+='</div><div class="ed-ptpop" id="edPTpop"><span class="muted">색이 칠해진 원소 = 준위 데이터 보유 (Ln³⁺ 전체, Sm²⁺·Eu²⁺·Tm²⁺·Yb²⁺, Cr³⁺·Mn²⁺·Mn⁴⁺·Ni²⁺)</span></div>';
  box.innerHTML=h;
  box.querySelectorAll("button[data-el]").forEach(b=>b.onclick=()=>{
    box.querySelectorAll("button.on").forEach(x=>x.classList.remove("on")); b.classList.add("on");
    const el=b.dataset.el, pop=$q(box,"#edPTpop"), ox=C.OX[el];
    if(!ox){ pop.innerHTML=`<span class="muted"><b>${el}</b>: ${HOSTEL.has(el)?"호스트 구성 원소 — 광학 준위 없음(화학식에 넣으면 호스트로 인식)":"준위 데이터가 없습니다"}</span>`; return; }
    pop.innerHTML=`<b>${el}</b>`+ox.map(q=>{ const k=C.ionKey(el,q), d=C.IONS[k]; return `<button class="btn sm" data-k="${k}">${C.ionLabel(k)} <span class="muted">${d?d.conf:""}</span></button>`; }).join("")+
      `<span class="muted">Add as column</span>`;
    pop.querySelectorAll("button[data-k]").forEach(x=>x.onclick=()=>addIonColumn(x.dataset.k));
  });
}
function ELORDER(){ return ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split(" "); }

function renderChips(){
  const el=$q(ROOT,"#edChips"); if(!el) return;
  if(!S.cols.length){ el.innerHTML='<span class="muted">No ions</span>'; return; }
  el.innerHTML=S.cols.map((c,i)=>{
    const d=C.IONS[c.ion]; const conc=(c.conc||[]).map((x,k)=>x==null?null:`${fmt(x,x<1?2:(x%1?1:0))}%${c.layers&&c.layers.length>1?"·"+(c.layers[k]==="core"?"core":"shell"):""}`).filter(Boolean).join(", ");
    const est=d.kind==="d"||(d.band&&true)&&d.kind!=="ln3";
    const dup=S.cols.filter(x=>x.ion===c.ion).length>1 ? ` #${S.cols.filter((x,k)=>x.ion===c.ion&&k<=i).length}`:"";
    return `<span class="ed-ion ${d.kind!=="ln3"?"est":""}" data-i="${i}"><b>${C.ionLabel(c.ion)}${dup}</b><span class="c">${esc(conc||d.conf)}</span>
      <button data-a="l" aria-label="Move left" title="Move left" ${i?"":"disabled"}>◂</button><button data-a="r" aria-label="Move right" title="Move right" ${i<S.cols.length-1?"":"disabled"}>▸</button><button data-a="d" aria-label="Duplicate ion column" title="같은 이온 열 하나 더 (동종 이온 간 CR·ETU 를 두 열로 그리기)">⧉</button><button class="x" data-a="x" aria-label="Remove column" title="Remove column">×</button></span>`;
  }).join("");
  el.querySelectorAll(".ed-ion").forEach(ch=>ch.addEventListener("click",e=>{
    const b=e.target.closest("button[data-a]"); if(!b||b.disabled) return; const i=+ch.dataset.i;
    if(b.dataset.a==="l"&&i>0){ [S.cols[i-1],S.cols[i]]=[S.cols[i],S.cols[i-1]]; }
    else if(b.dataset.a==="r"&&i<S.cols.length-1){ [S.cols[i+1],S.cols[i]]=[S.cols[i],S.cols[i+1]]; }
    else if(b.dataset.a==="d"){ const c=S.cols[i]; S.cols.splice(i+1,0,Object.assign({},c,{uid:uid(),conc:c.conc,notes:[]})); }
    else if(b.dataset.a==="x"){ const u=S.cols[i].uid; S.cols.splice(i,1); S.procs=S.procs.filter(p=>p.c!==u&&!(p.d&&(p.d.c===u||p.a.c===u))); }
    save(); rerender();
  }));
}
function renderHost(){
  const el=$q(ROOT,"#edHost"); if(!el) return;
  const H=S.host;
  const layers=(S.layers||[]).map(l=>`${esc(String(l.tag||"").replace(/^코어/,"Core").replace(/^셸/,"Shell"))} ${esc(l.host||l.hostKey||"")}`).join(" · ");
  el.innerHTML=`<span>${layers?`<b>${layers}</b>`:""}</span>
    <label title="최대 포논 에너지 — multiphonon 이완과 포논 보조 판정의 기준">ħω<sub>max</sub> <input type="number" id="edW" min="50" max="1500" step="10" value="${H.w||""}"> cm⁻¹</label>
    ${H.est?'<span class="badge warn" title="'+esc(H.src)+'">Est.</span>':(H.src?'<span class="badge gray" title="'+esc(H.src)+'">Source</span>':"")}
    <label>T <input type="number" id="edT" min="4" max="900" step="1" value="${H.T||295}"> K</label>
    <span class="muted">kT = ${fmt(C.KT_PER_K*(H.T||295))} cm⁻¹</span>`;
  $q(el,"#edW").onchange=e=>{ const v=+e.target.value; if(v>0){ S.host.w=v; S.host.user=true; S.host.est=true; S.host.src="User input"; save(); rerender(); } };
  $q(el,"#edT").onchange=e=>{ const v=+e.target.value; if(v>0){ S.host.T=v; save(); rerender(); } };
}
function renderWarn(){
  const el=$q(ROOT,"#edWarn"); if(!el) return;
  const notes=[]; S.warn.forEach(w=>notes.push(`<div class="note warn">${esc(w)}</div>`));
  S.cols.forEach(c=>(c.notes||[]).forEach(n=>notes.push(`<div class="note">${C.ionLabel(c.ion)} — ${esc(n)}</div>`)));
  if(S.presetRef) notes.push(`<div class="note">Preset mechanism reference: ${esc(S.presetRef)}</div>`);
  el.innerHTML=notes.join("");
}
function renderTools(){
  const w=ROOT;
  w.querySelectorAll("#edMode button").forEach(b=>{ const on=b.dataset.m===S.mode; b.classList.toggle("active",on); b.setAttribute("aria-selected",on?"true":"false"); });
  $q(w,"#edStark").checked=!!S.show.stark; $q(w,"#edAll").checked=!!S.show.allLabels; $q(w,"#edNums").checked=!!S.show.nums;
  $q(w,"#edEmis").value=S.show.emis; $q(w,"#edEmax").value=String(S.Emax); $q(w,"#edUnit").value=S.show.unit||"cm";
}
function renderPrompt(){
  const el=$q(ROOT,"#edPrompt"); if(!el) return;
  if(!PICK){ el.innerHTML= S.cols.length?'<span class="muted">준위·화살표에 마우스를 올리면 수치가 보입니다.</span>':""; return; }
  const n=PICK.steps.length;
  const T={et:["① Donor initial level (excited state)","② Donor final level (same ion, lower)","③ Acceptor initial level (other or same column)","④ Acceptor final level (upper) — colored outline = resonance grade"],
           ph:["① Absorption initial level","② Absorption final level (same ion, upper)"],
           em:["① Initial level","② Final level (same ion, lower)"]}[PICK.kind];
  el.innerHTML=`<span class="step">${T[Math.min(n,T.length-1)]}</span>${n?'<button class="btn sm" id="edUndo">Undo step</button>':""}<span class="muted">Esc = restart</span>`;
  const u=$q(el,"#edUndo"); if(u) u.onclick=()=>{ PICK.steps.pop(); rerender(); };
}

/* ---------------- SVG 도식 ---------------- */
const GEO={top:58,bot:34,axisL:66,strip:12,colMin:150,colMax:400,lab:58};
function drawSVG(){
  const wrap=$q(ROOT,"#edWrap"); if(!wrap) return;
  if(!MODEL.cols.length){ wrap.innerHTML=`<div class="ed-empty"><b>No composition yet</b><br>화학식을 입력하거나 프리셋·주기율표에서 이온을 고르세요.<br><span class="muted">준위값: LaF₃:Ln³⁺ (Carnall et al., J. Chem. Phys. 90, 3443 (1989))</span></div>`; return; }
  wrap.innerHTML=svgMarkup(palette(),Math.max(320,wrap.clientWidth||900),true);
  bindSVG(wrap.firstElementChild);
}
/* 도식 마크업 — 화면(현재 테마, 편집 표시 포함)과 내보내기(밝은 바탕, 편집 표시 없음)가 같은 함수를 쓴다 */
function svgMarkup(P,avail,interactive){
  const M=MODEL;
  const savedPick=PICK; if(!interactive) PICK=null;
  const n=M.cols.length;
  let cw=Math.max(GEO.colMin,Math.min(GEO.colMax,(avail-GEO.axisL-20)/n));
  const W=Math.round(GEO.axisL+n*cw+16);
  const H=Math.max(460,Math.min(760,Math.round((window.innerHeight||900)*0.72)));
  const plotH=H-GEO.top-GEO.bot, Emax=M.Emax;
  const y=E=>GEO.top+(1-E/Emax)*plotH;
  const colX=i=>GEO.axisL+i*cw;
  const out=[];
  out.push(`<svg class="ed-svg ${PICK?"m-"+PICK.kind:""}" id="edSvg" xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="Energy level diagram">`);
  out.push(`<defs>${arrowDefs(P)}</defs>`);
  out.push(`<rect x="0" y="0" width="${W}" height="${H}" fill="${P.bg}"/>`);
  /* 축 */
  const unit=S.show.unit==="eV";
  const step=unit? pickStep(Emax/C.CM_PER_EV,[0.25,0.5,1,2]) : pickStep(Emax,[2000,5000,10000]);
  const maxv=unit?Emax/C.CM_PER_EV:Emax;
  for(let v=0; v<=maxv+1e-9; v+=step){
    const E=unit?v*C.CM_PER_EV:v, yy=y(E);
    out.push(`<line x1="${GEO.axisL-6}" x2="${W-8}" y1="${yy}" y2="${yy}" stroke="${P.grid}" stroke-width="1"/>`);
    out.push(`<text x="${GEO.axisL-24}" y="${yy+4}" font-size="11" text-anchor="end" fill="${P.muted}">${unit?fmt(v,v%1?2:0):fmt(v/1000)}</text>`);
  }
  out.push(`<text x="14" y="${GEO.top+plotH/2}" font-size="11.5" fill="${P.ink2}" transform="rotate(-90 14 ${GEO.top+plotH/2})" text-anchor="middle">${unit?"Energy (eV)":"Energy (10³ cm⁻¹)"}</text>`);
  /* 광자 색 띠: 바닥에서 그 높이까지의 광자 파장 */
  const sx=GEO.axisL-18;
  for(let E=250;E<Emax;E+=250){ const nm=C.cm2nm(E+125); out.push(`<rect x="${sx}" y="${y(E+250)}" width="${GEO.strip-4}" height="${Math.max(0.6,y(E)-y(E+250))+0.4}" fill="${C.wlColor(nm,P.dark)}"/>`); }
  out.push(`<rect x="${sx}" y="${y(Emax)}" width="${GEO.strip-4}" height="${y(250)-y(Emax)}" fill="none" stroke="${P.line}" class="tipable" data-tip="strip"/>`);
  [2000,1000,800,700,600,500,400,300,250,200].forEach(nm=>{ const E=C.nm2cm(nm); if(E>Emax||E<400) return; out.push(`<line x1="${sx+GEO.strip-4}" x2="${sx+GEO.strip}" y1="${y(E)}" y2="${y(E)}" stroke="${P.ink2}"/>`); });
  out.push(`<text x="${sx+4}" y="${GEO.top-12}" font-size="10" fill="${P.muted}" text-anchor="middle">λ</text>`);
  /* 열 */
  const picked=pickState(M);
  M.cols.forEach((m,i)=>out.push(drawColumn(m,i,colX(i),cw,y,P,picked,M)));
  /* 광흡수(자동) · 방출 · 과정 */
  const lanes=allocLanes(M,cw);
  out.push(drawAutoPhotons(M,colX,cw,y,P,lanes));
  out.push(drawEmissions(M,colX,cw,y,P,lanes));
  out.push(drawProcs(M,colX,cw,y,P,lanes));
  out.push(`</svg>`);
  PICK=savedPick;
  return out.join("");
}
function pickStep(max,steps){ for(const s of steps){ if(max/s<=12) return s; } return steps[steps.length-1]; }
function arrowDefs(P){
  const mk=(id,c)=>`<marker id="${id}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="${c}"/></marker>`;
  return Object.entries(P.T).map(([k,c])=>mk("edA_"+k,c)).join("")+mk("edA_ink",P.ink2);
}
function colorMarker(c){ return "edA_c"+c.replace(/[^0-9a-f]/gi,""); }

/* 선택 단계에 따른 준위 강조 */
function pickState(M){
  if(!PICK) return null;
  const st=PICK.steps, n=st.length; const res={cands:new Set(),grade:{},chosen:new Set(st.map(s=>s.c+"|"+s.id))};
  const lv=(c,id)=>M.byUid[c]&&M.byUid[c].byId[id];
  M.cols.forEach(m=>m.lv.forEach(l=>{
    const k=m.c.uid+"|"+l.id; let ok=false;
    if(PICK.kind==="et"){
      if(n===0) ok=!l.ground;
      else if(n===1) ok=m.c.uid===st[0].c && l.E<lv(st[0].c,st[0].id).E-100;
      else if(n===2) ok=true;
      else if(n===3){ const a0=lv(st[2].c,st[2].id); ok=m.c.uid===st[2].c && l.E>a0.E+100;
        if(ok){ const D=lv(st[0].c,st[0].id), Df=lv(st[1].c,st[1].id); const c=C.classifyET({d:{ion:M.byUid[st[0].c].c.ion,from:D,to:Df},a:{ion:m.c.ion,from:a0,to:l}},M.cx); res.grade[k]=c; } }
    } else if(PICK.kind==="ph"){
      if(n===0) ok=true; else ok=m.c.uid===st[0].c && l.E>lv(st[0].c,st[0].id).E+100;
      if(n===1&&ok&&S.lasers.length){ const F=lv(st[0].c,st[0].id); let b=null; S.lasers.forEach(nm=>{ const c=C.classifyPhoton(F,l,C.nm2cm(nm),M.cx); if(!b||c.delta<b.delta) b=c; }); res.grade[k]=b; }
    } else if(PICK.kind==="em"){
      if(n===0) ok=!l.ground; else ok=m.c.uid===st[0].c && l.E<lv(st[0].c,st[0].id).E-100;
    }
    if(ok) res.cands.add(k);
  }));
  return res;
}
function gradeColor(g,P){ return g==null?P.muted:(g<=1?(P.dark?"#6cc68a":"#2a7141"):g===2?(P.dark?"#8fb3e6":"#2f5f9e"):g===3?(P.dark?"#e3b04f":"#8a5f0a"):(P.dark?"#e8847a":"#ab3f35")); }

/* 라벨 표시 여부 */
function keyLevels(m,M){
  const set=new Set([m.lv[0].id]);
  C.curatedFor(m.c.ion).forEach(t=>{ t.up.forEach(u=>set.add(u)); set.add(t.lo); });
  if(m.def.kind==="d") m.lv.forEach(l=>{ if(/^[46]|^3[AT]|^2E$|^1E$/.test(l.id)) set.add(l.id); });
  if(m.def.kind!=="ln3") m.lv.forEach(l=>set.add(l.id));
  M.procs.forEach(p=>p.levels.forEach(l=>{ if((p.md===m||p.ma===m||p.m===m)) set.add(l.id); }));
  M.photons.forEach(x=>{ if(x.m===m){ set.add(x.from.id); set.add(x.to.id); } });
  /* 희소한 이온은 전부 */
  const vis=m.lv.filter(l=>l.E<=M.Emax);
  if(vis.length<=16) vis.forEach(l=>set.add(l.id));
  return set;
}
function drawColumn(m,i,x0,cw,y,P,picked,M){
  const o=[]; const d=m.def;
  const bx0=x0+12, bx1=x0+cw-GEO.lab;
  const conc=(m.c.conc||[]).filter(v=>v!=null).map(v=>`${fmt(v,v<1?2:(v%1?1:0))}%`).join(", ");
  const same=S.cols.filter(c=>c.ion===m.c.ion); const idx=same.length>1?` #${same.indexOf(m.c)+1}`:"";
  o.push(`<text x="${(bx0+bx1)/2+10}" y="22" font-size="15" font-weight="700" text-anchor="middle" fill="${P.ink}">${esc(C.ionLabel(m.c.ion))}${idx}</text>`);
  o.push(`<text x="${(bx0+bx1)/2+10}" y="38" font-size="10.5" text-anchor="middle" fill="${P.muted}">${esc([d.conf,conc,d.kind!=="ln3"?"host-dependent":""].filter(Boolean).join(" · "))}</text>`);
  const keys=keyLevels(m,M);
  const vis=m.lv.filter(l=>l.E<=M.Emax+1);
  /* 라벨 배치(겹침 회피) */
  const labs=[];
  vis.forEach(l=>{ if(S.show.allLabels||keys.has(l.id)) labs.push({l,y:y(l.E)}); });
  labs.sort((a,b)=>b.y-a.y);
  const minGap=11.5; let last=1e9;
  labs.forEach(a=>{ a.ty=Math.min(a.y+3.5,last-minGap); last=a.ty; });
  /* 위로 넘친 경우 되밀기 */
  let lim=GEO.top-8; for(let k=labs.length-1;k>=0;k--){ if(labs[k].ty<lim){ labs[k].ty=lim; } lim=labs[k].ty+minGap; }
  const labOf=new Map(labs.map(a=>[a.l.id,a]));
  vis.forEach(l=>{
    const yy=y(l.E), k=m.c.uid+"|"+l.id;
    const cand=picked&&picked.cands.has(k), chosen=picked&&picked.chosen.has(k), dim=picked&&!cand&&!chosen;
    const g=picked&&picked.grade[k];
    let stroke=P.ink, sw=1.6, op=1;
    if(l.kind==="band"){ /* 밴드는 사각 그라디언트 */
      const y1=y(Math.min(l.hi,M.Emax)), y0=y(l.lo);
      o.push(`<g class="lv${dim?" dim":""}"><rect x="${bx0}" y="${y1}" width="${bx1-bx0}" height="${Math.max(2,y0-y1)}" fill="${P.accent}" opacity="${P.dark?.20:.13}"/>`+
        `<line x1="${bx0}" x2="${bx1}" y1="${yy}" y2="${yy}" stroke="${P.accent}" stroke-width="2" stroke-dasharray="${l.note&&/입력/.test(l.note)?"0":"5 3"}"/>`);
    } else {
      if(l.regime==="mpr"){ stroke=P.faint; sw=1.1; }
      else if(l.regime==="mixed"){ stroke=P.ink2; sw=1.4; }
      else if(l.ground||l.regime==="rad"){ stroke=P.ink; sw=2; }
      o.push(`<g class="lv${dim?" dim":""}">`);
      if(S.show.stark && l.hi-l.lo>6){ const a=y(l.hi), b=y(l.lo); o.push(`<rect x="${bx0}" y="${a}" width="${bx1-bx0}" height="${Math.max(0.8,b-a)}" fill="${stroke}" opacity="${P.dark?.16:.09}"/>`); }
      o.push(`<line x1="${bx0}" x2="${bx1}" y1="${yy}" y2="${yy}" stroke="${stroke}" stroke-width="${sw}" opacity="${op}"/>`);
    }
    if(chosen) o.push(`<rect x="${bx0-3}" y="${yy-4}" width="${bx1-bx0+6}" height="8" rx="3" fill="none" stroke="${P.accent}" stroke-width="2.2"/>`);
    else if(cand){ const gc=g?gradeColor(g.grade,P):P.accent;
      o.push(`<circle cx="${bx0-6}" cy="${yy}" r="${g?4.2:3}" fill="${gc}"/>`);
      if(g) o.push(`<line x1="${bx0}" x2="${bx1}" y1="${yy}" y2="${yy}" stroke="${gc}" stroke-width="2.4" opacity=".75"/>`); }
    const L=labOf.get(l.id);
    if(L){
      const tx=bx1+5, ty=L.ty;
      if(Math.abs(ty-3.5-yy)>2) o.push(`<path d="M${bx1} ${yy} L${tx-1} ${ty-3.5}" stroke="${P.faint}" stroke-width=".8" fill="none"/>`);
      o.push(`<text x="${tx}" y="${ty}" font-size="11" fill="${l.regime==="mpr"?P.muted:P.ink}">${svgTerm(l.t)}</text>`);
    }
    /* 선택·툴팁 영역 */
    const hh=Math.max(8,Math.min(16,(l.kind==="band"?(y(l.lo)-y(Math.min(l.hi,M.Emax))):0)+8));
    const focusable = picked ? cand : (L?true:false);
    o.push(`<rect class="lvhit tipable" data-tip="lv" data-c="${m.c.uid}" data-id="${esc(l.id)}" x="${bx0-2}" y="${yy-hh/2}" width="${bx1-bx0+4+(L?44:0)}" height="${hh}" fill="transparent" ${focusable?'tabindex="0" role="button"':''} aria-label="${esc(C.ionLabel(m.c.ion)+" "+C.termText(l.t)+" "+Math.round(l.E)+" cm-1")}"/>`);
    o.push(`</g>`);
  });
  return o.join("");
}
function svgTerm(t){
  const m=/^(\d+)([A-Z])(\d+(?:\/\d)?|\d?)(.*)$/.exec(t);
  if(!m) return esc(t);
  return `<tspan font-size="8" dy="-4">${m[1]}</tspan><tspan dy="4">${m[2]}</tspan>${m[3]?`<tspan font-size="8" dy="2.5">${m[3]}</tspan><tspan dy="-2.5">${esc(m[4])}</tspan>`:esc(m[4])}`;
}

/* 열마다 화살표 차선: 광흡수(왼쪽) → 과정 → 방출(오른쪽) */
function allocLanes(M,cw){
  const lanes={};
  M.cols.forEach(m=>{ lanes[m.c.uid]={items:[]}; });
  const add=(uidc,key)=>{ const L=lanes[uidc]; if(!L) return; if(!L.items.includes(key)) L.items.push(key); };
  M.photons.forEach((x,k)=>add(x.m.c.uid,"ph"+k));
  M.procs.forEach(p=>{
    if(p.broken) return;
    if(p.kind==="ET"){ add(p.md.c.uid,"pd"+p.md.c.uid+"|"+p.Df.id+"|"+p.Dt.id); add(p.ma.c.uid,"pa"+p.p.uid); }
    else add(p.m.c.uid,"p"+p.p.uid);
  });
  M.emis.forEach((e,k)=>add(e.m.c.uid,"em"+k));
  Object.values(lanes).forEach(L=>{ L.n=L.items.length; });
  return {pos:(uidc,key,x0)=>{ const L=lanes[uidc]; const i=L.items.indexOf(key); const bx0=x0+18, bx1=x0+cw-GEO.lab-6; const n=Math.max(1,L.n); return bx0+(i+0.5)*(bx1-bx0)/n; }};
}
function arrow(x,y1,y2,color,opt){
  opt=opt||{};
  const id=opt.marker||"";
  const dash=opt.dash?` stroke-dasharray="${opt.dash}"`:"";
  return `<line x1="${x}" x2="${x}" y1="${y1}" y2="${y2}" stroke="${color}" stroke-width="${opt.w||2}"${dash} marker-end="url(#${id})" stroke-linecap="round"/>`;
}
function coloredMarker(color){ const id="edC"+Math.abs(hash(color)); return {id,def:`<marker id="${id}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="${color}"/></marker>`}; }
function hash(s){ let h=0; for(let i=0;i<s.length;i++) h=(h*31+s.charCodeAt(i))|0; return h; }
function drawAutoPhotons(M,colX,cw,y,P,lanes){
  const o=[]; const defs=new Map();
  M.photons.forEach((x,k)=>{
    const col=C.wlColor(x.nm,P.dark); const mk=coloredMarker(col); defs.set(mk.id,mk.def);
    const xx=lanes.pos(x.m.c.uid,"ph"+k,colX(x.m.i));
    const g=x.c.grade;
    const dash=g===0?null:(g<=2?"5 3":"2 3");
    const y1=y(x.from.E), y2=y(x.to.E)+3;
    o.push(`<g class="proc tipable" data-tip="auto" data-k="${k}" tabindex="-1">`+
      `<line x1="${xx}" x2="${xx}" y1="${y1}" y2="${y2}" stroke="transparent" stroke-width="10"/>`+
      arrow(xx,y1,y2,col,{marker:mk.id,w:g>=3?1.6:2.6,dash})+`</g>`);
    if(S.show.nums) o.push(`<text x="${xx+3}" y="${(y1+y2)/2}" font-size="9.5" fill="${P.muted}">${x.nm}</text>`);
  });
  return `<defs>${[...defs.values()].join("")}</defs>`+o.join("");
}
function drawEmissions(M,colX,cw,y,P,lanes){
  const o=[]; const defs=new Map();
  M.emis.forEach((e,k)=>{
    const col=C.wlColor(e.nm,P.dark); const mk=coloredMarker(col); defs.set(mk.id,mk.def);
    const xx=lanes.pos(e.m.c.uid,"em"+k,colX(e.m.i));
    const Eu=e.up.kind==="band"?e.up.E-(e.up.stokes||0):e.up.E;
    const y1=y(Eu), y2=y(e.lo.E)-3;
    const faint=S.show.emis==="all";
    o.push(`<g class="proc tipable" data-tip="em" data-k="${k}"><line x1="${xx}" x2="${xx}" y1="${y1}" y2="${y2}" stroke="transparent" stroke-width="10"/>`+
      arrow(xx,y1,y2,col,{marker:mk.id,w:faint?1.3:2.2})+`</g>`);
    if(S.show.nums) o.push(`<text x="${xx+3}" y="${y2-6}" font-size="9.5" fill="${P.muted}">${fmt(e.nm)}</text>`);
  });
  return `<defs>${[...defs.values()].join("")}</defs>`+o.join("");
}
function squiggle(x,y1,y2,color){
  /* multiphonon 이완: 물결 화살표 */
  const n=Math.max(3,Math.round((y2-y1)/7)); let d=`M${x} ${y1}`; const h=(y2-y1)/n;
  for(let i=0;i<n;i++){ const s=i%2?-3.2:3.2; d+=` Q${x+s} ${y1+h*(i+.5)} ${x} ${y1+h*(i+1)}`; }
  return `<path d="${d}" fill="none" stroke="${color}" stroke-width="1.6" marker-end="url(#edA_MPR)"/>`;
}
function drawProcs(M,colX,cw,y,P,lanes){
  const o=[]; const defs=new Map();
  M.procs.forEach(p=>{
    if(p.broken) return;
    const sel=S.sel===p.p.uid;
    if(p.kind==="ET"){
      const col=P.T[p.type]||P.T.ET;
      const dkey=p.md.c.uid+"|"+p.Df.id+"|"+p.Dt.id;
      const xd=lanes.pos(p.md.c.uid,"pd"+dkey,colX(p.md.i)), xa=lanes.pos(p.ma.c.uid,"pa"+p.p.uid,colX(p.ma.i));
      const dShared=M.procs.filter(q=>q.kind==="ET"&&!q.broken&&q.md.c.uid+"|"+q.Df.id+"|"+q.Dt.id===dkey);
      const dFirst=dShared[0]===p, dMixed=new Set(dShared.map(q=>q.type)).size>1;
      const yd1=y(p.Df.E), yd2=y(p.Dt.E)-3, ya1=y(p.Af.E), ya2=y(p.At.E)+3;
      const w=sel?3:2.1;
      const mdY=(yd1+yd2)/2, maY=(ya1+ya2)/2;
      const cxm=(xd+xa)/2, cym=Math.min(mdY,maY)-Math.max(18,Math.abs(xa-xd)*0.12);
      const sameX=Math.abs(xa-xd)<2;
      o.push(`<g class="proc tipable${sel?" sel":""}" data-tip="proc" data-u="${p.p.uid}" tabindex="0" role="button" aria-label="${esc(typeName(p.type))} ${p.n}">`);
      o.push(`<line x1="${xd}" x2="${xd}" y1="${yd1}" y2="${yd2}" stroke="transparent" stroke-width="12"/><line x1="${xa}" x2="${xa}" y1="${ya1}" y2="${ya2}" stroke="transparent" stroke-width="12"/>`);
      if(dFirst||sel){ const dc=dMixed&&!sel?P.ink2:col; o.push(arrow(xd,yd1,yd2,dc,{marker:dMixed&&!sel?"edA_ink":"edA_"+(P.T[p.type]?p.type:"ET"),w,dash:"6 3"})); }
      o.push(arrow(xa,ya1,ya2,col,{marker:"edA_"+(P.T[p.type]?p.type:"ET"),w,dash:"6 3"}));
      if(p.p.coop) o.push(arrow(xd-7,yd1,yd2,col,{marker:"edA_"+p.type,w:w*0.8,dash:"6 3"}));
      if(!sameX) o.push(`<path d="M${xd} ${mdY} Q${cxm} ${cym} ${xa} ${maY}" fill="none" stroke="${col}" stroke-width="${sel?1.8:1.2}" stroke-dasharray="2 3"/>`);
      const bx=sameX?xd+9:cxm, by=sameX?mdY:(mdY+maY)/4+cym/2;
      o.push(`<circle cx="${bx}" cy="${by}" r="8.5" fill="${P.bg}" stroke="${col}" stroke-width="${sel?2.2:1.4}"/><text x="${bx}" y="${by+3.6}" font-size="10" font-weight="700" text-anchor="middle" fill="${col}">${p.n}</text>`);
      if(S.show.nums) o.push(`<text x="${bx+11}" y="${by+3.5}" font-size="9.5" fill="${P.muted}">${p.type} ${p.c.dE>0?"+":""}${fmt(p.c.dE)}</text>`);
      o.push(`</g>`);
    } else if(p.kind==="PH"){
      const nm=p.nm||p.need; const col=C.wlColor(nm,P.dark); const mk=coloredMarker(col); defs.set(mk.id,mk.def);
      const xx=lanes.pos(p.m.c.uid,"p"+p.p.uid,colX(p.m.i));
      const g=p.c?p.c.grade:4; const dash=!p.c?"1 3":(g===0?null:(g<=2?"5 3":"2 3"));
      const y1=y(p.F.E), y2=y(p.T.E)+3;
      o.push(`<g class="proc tipable${sel?" sel":""}" data-tip="proc" data-u="${p.p.uid}" tabindex="0" role="button" aria-label="${p.type} ${p.n}"><line x1="${xx}" x2="${xx}" y1="${y1}" y2="${y2}" stroke="transparent" stroke-width="12"/>`+
        arrow(xx,y1,y2,col,{marker:mk.id,w:sel?3.4:2.6,dash})+
        `<text x="${xx-4}" y="${(y1+y2)/2}" font-size="9.5" font-weight="700" text-anchor="end" fill="${P.muted}">${p.n}</text></g>`);
    } else {
      const xx=lanes.pos(p.m.c.uid,"p"+p.p.uid,colX(p.m.i));
      const y1=y(p.F.E), y2=y(p.T.E)-3;
      let body;
      if(p.kind==="MPR") body=squiggle(xx,y1,y2,P.T.MPR);
      else { const col=C.wlColor(p.nm,P.dark); const mk=coloredMarker(col); defs.set(mk.id,mk.def); body=arrow(xx,y1,y2,col,{marker:mk.id,w:sel?3:2.3}); }
      o.push(`<g class="proc tipable${sel?" sel":""}" data-tip="proc" data-u="${p.p.uid}" tabindex="0" role="button" aria-label="${p.kind} ${p.n}"><line x1="${xx}" x2="${xx}" y1="${y1}" y2="${y2}" stroke="transparent" stroke-width="12"/>${body}</g>`);
    }
  });
  return `<defs>${[...defs.values()].join("")}</defs>`+o.join("");
}

/* ---------------- 상호작용 ---------------- */
function bindSVG(svg){
  if(!svg) return;
  // 마우스: 올리면 풍선. 손가락(touch): hover 가 없으므로 탭 = 풍선 (손가락이 떨어질 때의 pointerleave 로 바로 닫히지 않게 touch 는 무시)
  svg.addEventListener("pointermove",e=>{ if(e.pointerType==="touch") return; const t=e.target.closest(".tipable"); if(t) showTip(t,e.clientX,e.clientY); else hideTip(); });
  svg.addEventListener("pointerleave",e=>{ if(e.pointerType==="touch") return; hideTip(); });
  let TD=null;
  svg.addEventListener("pointerdown",e=>{ if(e.pointerType!=="touch") return; TD={x:e.clientX,y:e.clientY,t:e.target.closest(".tipable"),time:Date.now()}; });
  svg.addEventListener("pointerup",e=>{ if(e.pointerType!=="touch"||!TD) return; const d=Math.hypot(e.clientX-TD.x,e.clientY-TD.y); const tap=d<12&&Date.now()-TD.time<700; const t=TD.t; TD=null;
    if(!tap) return; const x=e.clientX, y=e.clientY; // 탭 뒤 click 이 도식을 다시 그리며 풍선을 닫으므로 그다음에 연다
    if(t) setTimeout(()=>{ showTip(t,x,y); TIP_TOUCH=true; },30); else hideTip(); });
  svg.addEventListener("pointercancel",()=>{ TD=null; });
  svg.addEventListener("focusin",e=>{ const t=e.target.closest(".tipable"); if(t){ const r=t.getBoundingClientRect(); showTip(t,r.right,r.top+r.height/2); } });
  svg.addEventListener("focusout",hideTip);
  const act=e=>{
    const lv=e.target.closest(".lvhit"); const pr=e.target.closest(".proc[data-u]");
    if(PICK && lv){ onPick(lv.dataset.c,lv.dataset.id); return; }
    if(pr){ S.sel=pr.dataset.u; rerender(); return; }
    if(lv){ S.sel="lv:"+lv.dataset.c+"|"+lv.dataset.id; rerender(); return; }
  };
  svg.addEventListener("click",act);
  svg.addEventListener("keydown",e=>{
    if(e.key==="Enter"||e.key===" "){ e.preventDefault(); act(e); }
    if((e.key==="Delete"||e.key==="Backspace") && S.sel && !String(S.sel).startsWith("lv:")){ e.preventDefault(); delProc(S.sel); }
  });
}
function onPick(c,id){
  const M=MODEL; const st=PICK.steps; const k=c+"|"+id;
  const ps=pickState(M); if(!ps.cands.has(k)){ note("이 단계에서 고를 수 없는 준위입니다"); return; }
  st.push({c,id});
  if(PICK.kind==="et" && st.length===4){
    const p={uid:uid(),kind:"ET",d:{c:st[0].c,f:st[0].id,t:st[1].id},a:{c:st[2].c,f:st[2].id,t:st[3].id},coop:false,type:null};
    S.procs.push(p); S.sel=p.uid; PICK.steps=[]; save();
  } else if(PICK.kind==="ph" && st.length===2){
    const p={uid:uid(),kind:"PH",c:st[0].c,f:st[0].id,t:st[1].id,nm:null}; S.procs.push(p); S.sel=p.uid; PICK.steps=[]; save();
  } else if(PICK.kind==="em" && st.length===2){
    const m=MODEL.byUid[st[0].c], F=m.byId[st[0].id], T=m.byId[st[1].id];
    /* 간격이 multiphonon 영역이면 MPR 로, 아니면 복사 전이로 시작 (나중에 바꿀 수 있음) */
    const hw=MODEL.cx.hw; const gap=Math.max(0,F.lo-T.hi);
    const kind=(hw && gap/hw<C.P_MPR)?"MPR":"RAD";
    const p={uid:uid(),kind,c:st[0].c,f:st[0].id,t:st[1].id}; S.procs.push(p); S.sel=p.uid; PICK.steps=[]; save();
  }
  rerender();
  /* 다음 후보로 초점 */
  const f=ROOT.querySelector('.lvhit[tabindex="0"]'); if(f && document.activeElement && document.activeElement.classList && document.activeElement.classList.contains("lvhit")) f.focus();
}
function sameProc(a,b){ if(a.kind!==b.kind) return false;
  if(a.kind==="ET") return a.d.c===b.d.c&&a.d.f===b.d.f&&a.d.t===b.d.t&&a.a.c===b.a.c&&a.a.f===b.a.f&&a.a.t===b.a.t;
  return a.c===b.c&&a.f===b.f&&a.t===b.t&&(a.kind!=="PH"||(+a.nm||0)===(+b.nm||0)); }
function delProc(u){ S.procs=S.procs.filter(p=>p.uid!==u); if(S.sel===u) S.sel=null; save(); rerender(); }

/* ---------------- 툴팁 ---------------- */
function ensureTip(){ TIP=document.getElementById("edTip"); if(!TIP){ TIP=document.createElement("div"); TIP.id="edTip"; TIP.className="ed-tip"; TIP.setAttribute("role","tooltip"); document.body.appendChild(TIP); } }
function hideTip(){ if(TIP) TIP.classList.remove("on"); TIP_TOUCH=false; }
/* 손가락으로 연 풍선: 도식 밖을 누르거나 쪽을 스크롤하면 닫음 */
let TIP_TOUCH=false;
if(typeof document!=="undefined"){ document.addEventListener("pointerdown",e=>{ if(!TIP_TOUCH) return; if(e.target.closest&&(e.target.closest(".tipable")||e.target.closest("#edTip"))) return; hideTip(); },true);
  window.addEventListener("scroll",()=>{ if(TIP_TOUCH) hideTip(); },{passive:true}); }
function showTip(t,x,yv){
  const html=tipHTML(t); if(!html){ hideTip(); return; }
  TIP.innerHTML=html; TIP.classList.add("on");
  const r=TIP.getBoundingClientRect(); let L=x+14, T=yv+14;
  if(L+r.width>window.innerWidth-8) L=x-r.width-14; if(T+r.height>window.innerHeight-8) T=Math.max(8,yv-r.height-10);
  TIP.style.left=Math.max(8,L)+"px"; TIP.style.top=T+"px";
}
const E3 = E => `${fmt(E)} cm⁻¹ · ${fmt(C.cm2eV(E),3)} eV`;
function tipHTML(t){
  const M=MODEL, k=t.dataset.tip;
  if(k==="strip") return `<h4>Photon color strip</h4>바닥(0)에서 그 높이까지 전이할 때의 광자 파장 색. 눈금: 2000·1000·800·700·600·500·400·300·250 nm`;
  if(k==="lv"){
    const m=M.byUid[t.dataset.c]; const l=m&&m.byId[t.dataset.id]; if(!l) return "";
    return levelTip(m,l);
  }
  if(k==="auto"){ const x=M.photons[+t.dataset.k]; if(!x) return "";
    return photonTip(x.m,x.from,x.to,x.nm,x.c,"Auto candidate (based on active sources) — 활성 광원 기준"); }
  if(k==="em"){ const e=M.emis[+t.dataset.k]; if(!e) return "";
    return `<h4>${esc(C.ionLabel(e.m.c.ion))} ${C.termHTML(e.pair?e.pair.join("+"):e.up.t)} → ${C.termHTML(e.lo.t)} · Emission</h4>
      <table><tr><td class="k">λ (barycenter)</td><td><b>${fmt(e.nm,1)} nm</b> · ${bandEn(e.nm)}</td></tr>
      ${e.typ?`<tr><td class="k">Literature value</td><td>${e.typ} nm</td></tr>`:""}
      <tr><td class="k">Energy</td><td>${E3(C.nm2cm(e.nm))}</td></tr>
      ${e.up.kind==="f"&&e.lo.kind==="f"?`<tr><td class="k">Stark window</td><td>${fmt(C.cm2nm(e.up.hi-e.lo.lo),1)}–${fmt(C.cm2nm(Math.max(1,e.up.lo-e.lo.hi)),1)} nm</td></tr>`:""}</table>
      <div class="s">${esc(e.src||"")}${e.up.kind==="band"?" · 밴드 방출은 흡수 에너지에서 Stokes 이동을 뺀 값":""}</div>`; }
  if(k==="proc"){ const p=M.procs.find(x=>x.p.uid===t.dataset.u); if(!p) return ""; return procTip(p); }
  return "";
}
function levelTip(m,l){
  const g=l.ground?"Ground level":(l.regime==="rad"?"Radiative-dominant (metastable)":l.regime==="mixed"?"Radiative / multiphonon competition":l.regime==="mpr"?"Multiphonon-dominant":"—");
  const srcN={e:"Measured (all Stark levels observed, Carnall 1989)",c:"Calculated (mean of Carnall 1989 calculated Stark levels)",q:"Recalculated (Carnall parameters, qlanth Hamiltonian)",f:"Free-ion calculation (strong J mixing)",n:"NIST free ion",ts:"Tanabe–Sugano (entered Dq·B·C)",band:"Host-dependent band (parameter)"}[l.src]||"";
  return `<h4>${esc(C.ionLabel(m.c.ion))} ${C.termHTML(l.t)}</h4><table>
   <tr><td class="k">Energy</td><td><b>${E3(l.E)}</b></td></tr>
   ${!l.ground?`<tr><td class="k">λ from ground</td><td>${fmt(C.cm2nm(l.E-m.lv[0].E),1)} nm</td></tr>`:""}
   ${l.hi-l.lo>1?`<tr><td class="k">${l.kind==="f"?"Stark range":"Width (approx.)"}</td><td>${fmt(l.lo)}–${fmt(l.hi)} cm⁻¹</td></tr>`:""}
   ${!l.ground&&l.gap!=null?`<tr><td class="k">Gap below</td><td>${fmt(l.gap)} cm⁻¹ (edge) · p = ${l.p!=null?fmt(l.p,1):"—"} phonons</td></tr>`:""}
   <tr><td class="k">Relaxation</td><td>${g}${l.cascade&&l.cascade!==l.id?` → cascades to ${C.termHTML(m.byId[l.cascade]?m.byId[l.cascade].t:l.cascade)}`:""}</td></tr>
   ${l.comp?`<tr><td class="k">Composition</td><td>${esc(l.comp)}</td></tr>`:""}
  </table><div class="s">${esc(srcN)}${l.note?" · "+esc(l.note):""}</div>`;
}
function photonTip(m,F,T,nm,c,sub){
  if(!c) return `<h4>${esc(C.ionLabel(m.c.ion))} ${C.termHTML(F.t)} → ${C.termHTML(T.t)}</h4>필요 파장 ${fmt(C.cm2nm(T.E-F.E),1)} nm — 활성 광원 없음 (no active source)`;
  return `<h4>${F.ground?"GSA":"ESA"} · ${esc(C.ionLabel(m.c.ion))} ${C.termHTML(F.t)} → ${C.termHTML(T.t)}</h4><table>
   <tr><td class="k">Source</td><td><b>${nm} nm</b> (${fmt(C.nm2cm(nm))} cm⁻¹)</td></tr>
   <tr><td class="k">Transition window</td><td>${fmt(c.w.lo)}–${fmt(c.w.hi)} cm⁻¹ (center ${fmt(c.w.E)}, ${fmt(C.cm2nm(c.w.E),1)} nm)</td></tr>
   <tr><td class="k">Grade</td><td><span class="gr g${c.grade}">${gradeName(c)}</span> δ = ${fmt(c.delta)} cm⁻¹${c.nPh!=null&&c.delta>0?` (${fmt(c.nPh,1)} ħω)`:""}</td></tr>
   ${c.delta>0?`<tr><td class="k">Phonon</td><td>${c.endo?`absorption required (anti-Stokes) · exp(−δ/kT) = ${c.boltz.toExponential(1)}`:"emission (Stokes, excess energy)"}</td></tr>`:""}
  </table><div class="s">${esc(sub||"")}</div>`;
}
function procTip(p){
  if(p.broken) return "<h4>Unresolved</h4>이 과정이 가리키는 준위가 현재 열에 없습니다";
  if(p.kind==="ET"){
    const c=p.c;
    return `<h4>${p.n}. ${esc(typeName(p.type))}${p.type!==p.auto?` <span class="k">(auto: ${p.auto})</span>`:""}</h4><table>
     <tr><td class="k">donor</td><td>${esc(C.ionLabel(p.md.c.ion))} ${C.termHTML(p.Df.t)}→${C.termHTML(p.Dt.t)} · ${fmt(c.dw.E)} cm⁻¹${p.p.coop?" ×2":""} (${fmt(C.cm2nm(c.dw.E),0)} nm)</td></tr>
     <tr><td class="k">acceptor</td><td>${esc(C.ionLabel(p.ma.c.ion))} ${C.termHTML(p.Af.t)}→${C.termHTML(p.At.t)} · ${fmt(c.aw.E)} cm⁻¹ (${fmt(C.cm2nm(c.aw.E),0)} nm)</td></tr>
     <tr><td class="k">ΔE (center)</td><td><b>${c.dE>0?"+":""}${fmt(c.dE)} cm⁻¹</b> · ${fmt(c.dE/C.CM_PER_EV*1000,0)} meV</td></tr>
     <tr><td class="k">Grade</td><td><span class="gr g${c.grade}">${gradeName(c)}</span>${c.delta>0?` δ = ${fmt(c.delta)} cm⁻¹ (${c.nPh!=null?fmt(c.nPh,1)+" ħω":""})`:""}</td></tr>
     ${c.delta>0?`<tr><td class="k">Phonon</td><td>${c.endo?`endothermic — phonon absorption required, exp(−δ/kT) = ${c.boltz.toExponential(1)}`:"exothermic — excess energy emitted as phonons"}</td></tr>`:""}
    </table><div class="s">Stark windows: donor ${fmt(c.dw.lo)}–${fmt(c.dw.hi)}, acceptor ${fmt(c.aw.lo)}–${fmt(c.aw.hi)} cm⁻¹. ${esc((C.TYPE_INFO[p.type]||{}).desc||"")}</div>`;
  }
  if(p.kind==="PH") return photonTip(p.m,p.F,p.T,p.nm,p.c,p.p.nm?"Specified source":"Best-matching active source");
  return `<h4>${p.n}. ${p.kind==="MPR"?"Multiphonon relaxation":"Radiative transition"} · ${esc(C.ionLabel(p.m.c.ion))} ${C.termHTML(p.F.t)} → ${C.termHTML(p.T.t)}</h4>
    <table><tr><td class="k">Gap</td><td>${fmt(p.F.E-p.T.E)} cm⁻¹ · ${fmt(p.nm,1)} nm</td></tr>
    <tr><td class="k">Edge gap</td><td>${fmt(Math.max(0,p.F.lo-p.T.hi))} cm⁻¹${MODEL.cx.hw?` = ${fmt(Math.max(0,p.F.lo-p.T.hi)/MODEL.cx.hw,1)} ħω`:""}</td></tr></table>`;
}

/* ---------------- 옆 패널 ---------------- */
function renderLasers(){
  const el=$q(ROOT,"#edLasCard"); if(!el) return; const P=palette();
  const pa=MODEL.pa;
  el.innerHTML=`<h3>Excitation sources <span class="hint">GSA·ESA 자동 후보의 기준</span></h3>
   <div class="ed-lasers">${S.lasers.map((nm,i)=>`<span class="ed-laser"><i style="background:${C.wlColor(nm,P.dark)}"></i>${nm} nm<button data-i="${i}" aria-label="Remove ${nm} nm">×</button></span>`).join("")}
   <input type="number" id="edLasIn" min="150" max="5000" step="1" placeholder="Add nm" aria-label="Add excitation wavelength (nm)">
   <label class="ed-tog"><input type="checkbox" id="edGsa" ${S.show.gsa?"checked":""}> Show auto candidates</label></div>
   ${pa.map(r=>`<div class="ed-pa"><b>Photon avalanche loop candidate · ${esc(C.ionLabel(r.m.c.ion))} · ${r.nm} nm</b><br>
     GSA ${C.termHTML(r.gsa.to.t)} <span class="gr g${r.gsa.c.grade}">${gradeName(r.gsa.c)}</span> ·
     ESA ${C.termHTML(r.M.t)}→${C.termHTML(r.U.t)} <span class="gr g${r.esa.grade}">${gradeName(r.esa)}</span>${r.U.id!==r.H.id?` ⇝ ${C.termHTML(r.H.t)}`:""} ·
     CR ${C.termHTML(r.H.t)}+${C.termHTML(MODEL.cols[r.m.i].lv[0].t)}→${C.termHTML(r.X.t)}+${C.termHTML(r.Y.t)} <span class="gr g${r.cr.grade}">${gradeName(r.cr)}</span>
     <div class="muted" style="margin-top:4px">준위 구조상 필요 조건만 점검 — 실제 avalanche 여부는 σ<sub>ESA</sub>/σ<sub>GSA</sub>, CR 속도(농도), metastable 수명에 달림 (Joubert, Opt. Mater. 11, 181 (1999); Lee et al., Nature 589, 230 (2021)).</div>
     <div class="row" style="margin-top:6px">${S.procs.some(p=>p.grp==="pa:"+r.m.c.uid+"|"+r.nm)?`<button class="btn sm" data-parm="${r.m.c.uid}|${r.nm}">Remove this loop from diagram</button><span class="muted">In diagram</span>`:`<button class="btn sm" data-pa="${r.m.c.uid}|${r.nm}">Add loop to diagram</button>`}</div></div>`).join("")}`;
  el.querySelectorAll(".ed-laser button").forEach(b=>b.onclick=()=>{ S.lasers.splice(+b.dataset.i,1); S.lasersUser=true; save(); rerender(); });
  const inp=$q(el,"#edLasIn");
  inp.addEventListener("keydown",e=>{ if(e.key==="Enter"){ const v=Math.round(+inp.value); if(v>=150&&v<=5000&&!S.lasers.includes(v)){ S.lasers.push(v); S.lasers.sort((a,b)=>a-b); S.lasersUser=true; save(); rerender(); } } });
  $q(el,"#edGsa").onchange=e=>{ S.show.gsa=e.target.checked; save(); rerender(); };
  el.querySelectorAll("button[data-pa]").forEach(b=>b.onclick=()=>{
    const [cu,nm]=b.dataset.pa.split("|"); const r=MODEL.pa.find(x=>x.m.c.uid===cu&&String(x.nm)===nm); if(!r) return;
    const g=MODEL.byUid[cu].lv[0]; const grp="pa:"+cu+"|"+nm; let added=0, dup=0;
    const add=o=>{ if(S.procs.some(q=>sameProc(q,o))){ dup++; return; } S.procs.push(Object.assign({uid:uid(),grp},o)); added++; };
    add({kind:"PH",c:cu,f:g.id,t:r.gsa.to.id,nm:+nm});
    if(r.gsa.to.cascade&&r.gsa.to.cascade!==r.gsa.to.id) add({kind:"MPR",c:cu,f:r.gsa.to.id,t:r.gsa.to.cascade});
    else if(r.gsa.to.id!==r.M.id) add({kind:"MPR",c:cu,f:r.gsa.to.id,t:r.M.id});
    add({kind:"PH",c:cu,f:r.M.id,t:r.U.id,nm:+nm});
    if(r.U.id!==r.H.id) add({kind:"MPR",c:cu,f:r.U.id,t:r.H.id});
    add({kind:"ET",d:{c:cu,f:r.H.id,t:r.X.id},a:{c:cu,f:g.id,t:r.Y.id},coop:false,type:null});
    save(); rerender(); note(added?(`루프 과정 ${added}개를 추가했습니다`+(dup?` (이미 있던 ${dup}개는 건너뜀)`:"")):"이 루프의 과정은 이미 모두 도식에 있습니다");
  });
  el.querySelectorAll("button[data-parm]").forEach(b=>b.onclick=()=>{ const grp="pa:"+b.dataset.parm; const n=S.procs.filter(p=>p.grp===grp).length; S.procs=S.procs.filter(p=>p.grp!==grp); if(S.sel&&!S.procs.some(p=>p.uid===S.sel)) S.sel=null; save(); rerender(); note(`루프 과정 ${n}개를 지웠습니다`); });
}
function procLabel(p){
  const L=m=>esc(C.ionLabel(m.c.ion));
  if(p.broken) return "Unresolved";
  if(p.kind==="ET") return `${L(p.md)} ${C.termHTML(p.Df.t)}→${C.termHTML(p.Dt.t)} ⇒ ${p.md===p.ma?"":L(p.ma)+" "}${C.termHTML(p.Af.t)}→${C.termHTML(p.At.t)}`;
  return `${L(p.m)} ${C.termHTML(p.F.t)}→${C.termHTML(p.T.t)}${p.kind==="PH"&&p.nm?` · ${p.nm} nm`:""}`;
}
function typeBadge(t){ const P=palette(); const c=P.T[t]||(t==="GSA"||t==="ESA"?P.accent:P.ink2); return `<span class="tb" style="background:${c}">${TYPE_SHORT[t]||t}</span>`; }
function renderProcs(){
  const el=$q(ROOT,"#edProcCard"); if(!el) return;
  const ps=MODEL.procs;
  el.innerHTML=`<h3>Processes <span class="hint">${ps.length}</span>${ps.length?'<button class="btn xs" id="edProcClr" style="margin-left:auto" title="도식에 넣은 과정을 모두 지운다 (준위·이온은 그대로)">Clear all</button>':""}</h3>`+(ps.length?`<table class="ed-tbl"><thead><tr><th>#</th><th>Type</th><th>Process</th><th class="num">ΔE</th><th>Grade</th><th></th></tr></thead><tbody>`+
    ps.map(p=>{ const c=p.c; const dE=p.kind==="ET"?c.dE:(p.kind==="PH"&&c?c.dE:null);
      const gr=(p.kind==="ET"||p.kind==="PH")&&c?`<span class="gr g${c.grade}">${gradeName(c)}</span>`:(p.kind==="PH"?'<span class="gr g4">No source</span>':"");
      return `<tr data-u="${p.p.uid}" class="${S.sel===p.p.uid?"sel":""}"><td>${p.n}</td><td>${typeBadge(p.type)}</td><td>${procLabel(p)}</td><td class="num">${dE==null?"":(dE>0?"+":"")+fmt(dE)}</td><td>${gr}</td><td><button class="btn xs" data-pdel="${p.p.uid}" aria-label="Delete process" title="Delete">×</button></td></tr>`; }).join("")+
    `</tbody></table><div class="muted" style="margin-top:6px">ΔE = donor(또는 광자) − acceptor 전이 에너지 [cm⁻¹]. 양수면 남는 에너지를 포논으로 방출.</div>`:
    `<div class="muted">‘Energy transfer’·‘Absorption’·‘Emission &amp; relaxation’ 모드에서 준위를 차례로 누르거나, 아래 ‘Path suggestions’에서 추가하세요.</div>`);
  el.querySelectorAll("tbody tr").forEach(tr=>tr.onclick=e=>{ if(e.target.closest("button[data-pdel]")) return; S.sel=tr.dataset.u; rerender(); });
  el.querySelectorAll("button[data-pdel]").forEach(b=>b.onclick=()=>delProc(b.dataset.pdel));
  const clr=$q(el,"#edProcClr"); if(clr) clr.onclick=()=>{ if(!confirm(`도식의 과정 ${S.procs.length}개를 모두 지울까요? (준위·이온은 그대로)`)) return; S.procs=[]; S.sel=null; save(); rerender(); note("과정을 모두 지웠습니다"); };
}
function renderDetail(){
  const el=$q(ROOT,"#edDetCard"); if(!el) return;
  const sel=S.sel;
  if(!sel){ el.hidden=true; return; } el.hidden=false;
  if(String(sel).startsWith("lv:")){
    const [cu,id]=sel.slice(3).split("|"); const m=MODEL.byUid[cu]; const l=m&&m.byId[id];
    if(!l){ el.hidden=true; return; }
    el.innerHTML=`<h3>Level</h3>`+levelTip(m,l).replace("<h4>","<div style='font-weight:700;margin-bottom:4px'>").replace("</h4>","</div>")+dParamsHTML(m);
    bindDParams(el,m); return;
  }
  const p=MODEL.procs.find(x=>x.p.uid===sel); if(!p){ el.hidden=true; return; }
  let body="";
  if(p.kind==="ET"){
    body=procTip(p).replace(/<h4>.*?<\/h4>/,"")+
    `<div class="row"><label class="muted">Type <select class="ed-sel" id="edType">${["auto","ET","ETU","CR","EM","BET","CET"].map(t=>`<option value="${t}" ${(t==="auto"&&!p.p.type)||t===p.p.type?"selected":""}>${t==="auto"?"Auto ("+p.auto+")":t}</option>`).join("")}</select></label>
     <label class="ed-tog"><input type="checkbox" id="edCoop" ${p.p.coop?"checked":""}> Cooperative (2 donors)</label></div>`;
  } else if(p.kind==="PH"){
    body=procTip(p).replace(/<h4>.*?<\/h4>/,"")+`<div class="row"><label class="muted">Source <select class="ed-sel" id="edPhNm"><option value="">Auto (best match)</option>${S.lasers.map(nm=>`<option value="${nm}" ${p.p.nm===nm?"selected":""}>${nm} nm</option>`).join("")}</select></label></div>`;
  } else {
    body=procTip(p).replace(/<h4>.*?<\/h4>/,"")+`<div class="row"><div class="seg" id="edRk"><button data-k="RAD" class="${p.kind==="RAD"?"active":""}">Radiative</button><button data-k="MPR" class="${p.kind==="MPR"?"active":""}">Multiphonon</button></div></div>`;
  }
  el.innerHTML=`<h3>${p.n}. ${typeBadge(p.type)} <span style="font-weight:500">${procLabel(p)}</span></h3>${body}
   <div class="row"><textarea id="edNote" placeholder="Note (reference etc.)">${esc(p.p.note||"")}</textarea></div>
   <div class="row right"><button class="btn sm" id="edDel">Delete</button></div>`;
  const ty=$q(el,"#edType"); if(ty) ty.onchange=e=>{ p.p.type=e.target.value==="auto"?null:e.target.value; save(); rerender(); };
  const co=$q(el,"#edCoop"); if(co) co.onchange=e=>{ p.p.coop=e.target.checked; save(); rerender(); };
  const pn=$q(el,"#edPhNm"); if(pn) pn.onchange=e=>{ p.p.nm=e.target.value?+e.target.value:null; save(); rerender(); };
  const rk=$q(el,"#edRk"); if(rk) rk.onclick=e=>{ const b=e.target.closest("button[data-k]"); if(!b) return; p.p.kind=b.dataset.k; save(); rerender(); };
  $q(el,"#edNote").onchange=e=>{ p.p.note=e.target.value; save(); };
  $q(el,"#edDel").onclick=()=>delProc(p.p.uid);
}
function dParamsHTML(m){
  const d=m.def;
  if(d.kind==="d"){ const p=Object.assign({},d.p,m.c.p||{});
    return `<div class="ed-dparam"><b>${C.ionLabel(m.c.ion)} crystal field</b>
      <label>Dq <input type="number" data-p="Dq" value="${p.Dq}" step="10"></label><label>B <input type="number" data-p="B" value="${p.B}" step="5"></label><label>C <input type="number" data-p="C" value="${p.C}" step="10"></label>
      <span class="muted">Dq/B ${fmt(p.Dq/p.B,2)} · C/B ${fmt(p.C/p.B,2)} · β = B/B₀ ${fmt(p.B/d.B0,2)}</span></div>
      <div class="muted" style="margin-top:4px">기본값: ${esc(d.preset)}${d.est?" (추정)":""}. 팔면체·스핀-궤도 무시. 격자 범위 Dq/B 0.05–4.0, C/B 3.4–5.6.</div>`; }
  if(d.band){ const b=Object.assign({},d.band,m.c.band||{});
    return `<div class="ed-dparam"><b>${esc(b.label)} band</b><label>First absorption <input type="number" data-b="E" value="${b.E==null?"":b.E}" step="100"> cm⁻¹</label>
      <label>Stokes <input type="number" data-b="stokes" value="${b.stokes||0}" step="100"></label>
      <span class="muted">${b.E?fmt(C.cm2nm(b.E),0)+" nm":"Value needed"}</span></div><div class="muted" style="margin-top:4px">${esc(b.note||"")}</div>`; }
  return "";
}
function bindDParams(el,m){
  el.querySelectorAll("input[data-p]").forEach(i=>i.onchange=()=>{ const v=+i.value; if(!(v>0)) return; m.c.p=Object.assign({},m.def.p,m.c.p||{},{[i.dataset.p]:v}); save(); rerender(); });
  el.querySelectorAll("input[data-b]").forEach(i=>i.onchange=()=>{ const v=i.value===""?null:+i.value; m.c.band=Object.assign({},m.c.band||{},{[i.dataset.b]:v}); save(); rerender(); });
}
let SUG=null, SUGF="all", SUGSIG="";
const sugSig=()=>JSON.stringify([S.cols.map(c=>c.uid+c.ion+JSON.stringify(c.p||0)+JSON.stringify(c.band||0)),S.host.w,S.host.T,MODEL&&MODEL.Emax]);
function renderSug(){
  const el=$q(ROOT,"#edSugCard"); if(!el) return;
  if(SUG && SUGSIG!==sugSig()) SUG=null;       /* 조성·호스트가 바뀌면 제안은 무효 */
  el.innerHTML=`<h3>Path suggestions <span class="hint">Stark 창 겹침 ≤ 2ħω, population 가능한 준위만</span></h3>
    <div class="row"><button class="btn sm" id="edSugRun">${SUG?"Find again":"Find"}</button>
    <div class="seg" id="edSugF">${["all","ET","ETU","CR","BET","EM"].map(t=>`<button data-f="${t}" class="${SUGF===t?"active":""}">${t==="all"?"All":t}</button>`).join("")}</div></div>
    <div class="ed-sug" id="edSugList"></div>`;
  $q(el,"#edSugRun").onclick=()=>{ runSug(); renderSug(); };
  $q(el,"#edSugF").onclick=e=>{ const b=e.target.closest("button[data-f]"); if(!b) return; SUGF=b.dataset.f; renderSug(); };
  const list=$q(el,"#edSugList");
  if(!SUG){ list.innerHTML=`<div class="muted">donor = 복사 우세·경쟁 준위, acceptor 출발 = 바닥·metastable 준위. 같은 이온 쌍도 포함.</div>`; return; }
  const items=SUG.filter(s=>SUGF==="all"||s.c.type===SUGF);
  list.innerHTML=items.length?items.map((s,i)=>{ const md=MODEL.cols[s.di], ma=MODEL.cols[s.ai];
    return `<div class="it">${typeBadge(s.c.type)}<div class="tx">${esc(C.ionLabel(md.c.ion))} ${C.termHTML(s.Di.t)}→${C.termHTML(s.Df.t)} ⇒ ${esc(C.ionLabel(ma.c.ion))} ${C.termHTML(s.Ai.t)}→${C.termHTML(s.Af.t)}
      <div class="mt"><span class="gr g${s.c.grade}">${gradeName(s.c)}</span> ΔE ${s.c.dE>0?"+":""}${fmt(s.c.dE)} cm⁻¹${s.c.endo?" · endothermic":""}</div></div>
      <button class="btn sm" data-s="${SUG.indexOf(s)}" aria-label="Add to diagram">＋</button></div>`; }).join(""):`<div class="muted">None of this type</div>`;
  list.querySelectorAll("button[data-s]").forEach(b=>b.onclick=()=>{
    const s=SUG[+b.dataset.s]; const md=MODEL.cols[s.di], ma=MODEL.cols[s.ai];
    let acol=ma.c.uid;
    /* 같은 이온 쌍이면 두 번째 열이 있으면 그쪽에 그린다 */
    if(s.di===s.ai){ const twin=S.cols.find(c=>c.ion===md.c.ion&&c.uid!==md.c.uid); if(twin) acol=twin.uid; }
    const p={uid:uid(),kind:"ET",d:{c:md.c.uid,f:s.Di.id,t:s.Df.id},a:{c:acol,f:s.Ai.id,t:s.Af.id},coop:false,type:null};
    S.procs.push(p); S.sel=p.uid; save(); rerender();
  });
}
function runSug(){
  const M=MODEL; SUGSIG=sugSig(); SUG=C.suggestET(M.cols.map(m=>m.c),M.cols.map(m=>m.lv),M.cx,{Emax:M.Emax,limit:80,maxGrade:2});
}
/* ---------- 이온 간 상호작용 (서로 다른 이온이 2종 이상일 때) ---------- */
let INT_SEL=null, INT_MAXN=6;
function renderInter(){
  const el=$q(ROOT,"#edIntCard"); if(!el) return; const M=MODEL;
  const ions=[...new Set(M.cols.map(m=>m.c.ion))]; if(ions.length<2){ el.hidden=true; el.innerHTML=""; return; } el.hidden=false;
  const all=C.interIon(M.cols.map(m=>m.c),M.cols.map(m=>m.lv),M.cx,{Emax:40000,maxGrade:4});   /* 도식 축 상한과 무관하게 40,000 cm⁻¹ 까지 본다 */
  const hw=M.cx.hw, T=M.cx.T||295; const keep=all.filter(x=>x.N==null||x.N<=INT_MAXN);
  /* 방향별 행렬: 행 = donor 출발 준위, 열 = acceptor 출발 준위, 칸 = δ 가 가장 작은 조합 */
  const dirs={}; keep.forEach(x=>{ const k=x.di+">"+x.ai; (dirs[k]=dirs[k]||[]).push(x); });
  const cell=x=>`<span class="gr g${x.c.grade}">${x.N!=null?(x.c.delta?`N ${fmt(x.N,1)}`:"Resonant"):esc(gradeName(x.c))}</span>${x.c.endo?`<span class="muted"> ↑${x.boltz.toExponential(0)}</span>`:""}`;
  const mats=Object.entries(dirs).map(([k,L])=>{ const [di,ai]=k.split(">").map(Number); const md=M.cols[di], ma=M.cols[ai];
    const rows=[...new Map(L.map(x=>[x.Di.id,x.Di])).values()].sort((a,b)=>b.E-a.E), cols=[...new Map(L.map(x=>[x.Ai.id,x.Ai])).values()].sort((a,b)=>a.E-b.E);
    const best=(r,c)=>L.filter(x=>x.Di.id===r.id&&x.Ai.id===c.id)[0];
    return `<div class="ed-int-dir"><div class="h"><b>${esc(C.ionLabel(md.c.ion))}</b> → <b>${esc(C.ionLabel(ma.c.ion))}</b> <span class="muted">donor 여기 준위(행)가 acceptor 준위(열)를 올린다${L.some(x=>x.c.type==="BET")?" · back transfer 포함":""}</span></div>
      <div class="tblwrap"><table class="ed-tbl ed-int"><thead><tr><th>${esc(C.ionLabel(md.c.ion))} donor ∖ ${esc(C.ionLabel(ma.c.ion))} acceptor</th>${cols.map(c=>`<th>${C.termHTML(c.t)}${c.ground?' <span class="muted">(ground)</span>':""}</th>`).join("")}</tr></thead>
      <tbody>${rows.map(r=>`<tr><th>${C.termHTML(r.t)} <span class="muted">${fmt(r.E)}</span></th>${cols.map(c=>{ const x=best(r,c); return x?`<td class="${INT_SEL===k+"|"+r.id+"|"+c.id?"sel":""}"><button class="ed-int-c" data-int="${k}|${r.id}|${c.id}" title="${esc(`${C.ionLabel(md.c.ion)} ${x.Di.t}→${x.Df.t} ⇒ ${C.ionLabel(ma.c.ion)} ${x.Ai.t}→${x.Af.t} · ${x.c.type} · δ ${Math.round(x.c.delta)} cm⁻¹`)}">${cell(x)}</button></td>`:`<td class="muted">—</td>`; }).join("")}</tr>`).join("")}</tbody></table></div></div>`; }).join("");
  /* 선택한 칸(없으면 전체 상위)의 조합 목록 */
  let list=keep, cap="Top 12 by smallest δ";
  if(INT_SEL){ const [k,r,c]=INT_SEL.split("|"); const [di,ai]=k.split(">").map(Number); list=keep.filter(x=>x.di===di&&x.ai===ai&&x.Di.id===r&&x.Ai.id===c); cap="All combinations in the selected cell"; }
  const items=list.slice(0,INT_SEL?40:12).map(x=>{ const md=M.cols[x.di], ma=M.cols[x.ai]; return `<div class="it">${typeBadge(x.c.type)}<div class="tx">${esc(C.ionLabel(md.c.ion))} ${C.termHTML(x.Di.t)}→${C.termHTML(x.Df.t)} ⇒ ${esc(C.ionLabel(ma.c.ion))} ${C.termHTML(x.Ai.t)}→${C.termHTML(x.Af.t)}
      <div class="mt"><span class="gr g${x.c.grade}">${esc(gradeName(x.c))}</span> δ ${fmt(x.c.delta)} cm⁻¹${x.N!=null?` · ${fmt(x.N,1)} phonons`:""} · ΔE ${x.c.dE>0?"+":""}${fmt(x.c.dE)}${x.c.endo?` · endothermic, detailed balance ×${x.boltz.toExponential(1)} (${T} K)`:""}</div></div>
      <button class="btn sm" data-intadd="${all.indexOf(x)}" aria-label="Add to diagram">＋</button></div>`; }).join("");
  const accC=(()=>{ const m=M.cols.find(q=>q.c.ion!==M.cols[0].c.ion); const v=m&&(m.c.conc||[]).find(q=>q!=null); return v!=null?v:""; })();
  const nre=/Lu/.test(S.host&&S.host.name||"")?1.42e22:/Y/.test(S.host&&S.host.name||"")?1.37e22:"";
  el.innerHTML=`<h3>Inter-ion Interactions <span class="hint">between different ions · required phonons ≤ <select class="ed-sel" id="edIntN">${[2,4,6,8].map(n=>`<option ${n===INT_MAXN?"selected":""}>${n}</option>`).join("")}</select></span></h3>
    ${mats||`<div class="muted">포논 ${INT_MAXN}개 안쪽의 조합이 없습니다.</div>`}
    <div class="muted" style="margin:6px 0 4px">칸 = 그 두 준위 사이에서 에너지가 가장 잘 맞는 조합. N = Stark 창 간격 δ / ħω${hw?` (ħω ${fmt(hw)} cm⁻¹)`:""}, ↑ = endothermic(포논 흡수)의 detailed balance 인자. 칸을 누르면 아래에 모든 조합이 나온다.</div>
    <div class="ed-sug"><div class="muted" style="margin-bottom:4px">${cap}${INT_SEL?` <button class="btn xs" id="edIntClr">Show all</button>`:""}</div>${items||'<div class="muted">None</div>'}</div>
    <details class="ed-src" style="margin-top:8px"><summary>Absolute transfer probability (Förster)</summary>
      <div class="muted" style="margin:6px 0">위 표는 <b>에너지 조건만</b> 본 상대 순위다. 실제 전달 속도는 (1) 두 전이의 선세기와 선모양의 겹침, (2) 이온 간 거리 분포에 각각 수 자릿수씩 좌우되는데, 이 둘은 준위 에너지만으로 정해지지 않는다.
        예: Yb→Tm ³H₆→³H₅ 는 δ 가 1,000 cm⁻¹ 를 넘는 phonon-assisted 전달인데도 실측에서 효율적이다. 그래서 여기서 '확률'을 내면 신뢰도가 낮다.
        측정된(또는 문헌의 같은 호스트) microparameter C<sub>DA</sub> 가 있으면 아래에서 정적 Förster 평균으로 전달 효율을 계산한다.</div>
      <div class="row" style="flex-wrap:wrap;gap:8px"><label class="ed-f">C<sub>DA</sub> [cm⁶ s⁻¹] <input type="number" id="edCda" step="any" placeholder="e.g. 1e-40"></label>
        <label class="ed-f">Donor intrinsic lifetime τ<sub>D</sub> [ms] <input type="number" id="edTauD" step="any" placeholder="acceptor-free sample"></label>
        <label class="ed-f">Acceptor concentration [%] <input type="number" id="edNApct" step="any" value="${accC}"></label>
        <label class="ed-f">Cation site density [cm⁻³] <input type="number" id="edNsite" step="any" value="${nre}" placeholder="e.g. NaYF₄ 1.37e22"></label></div>
      <div id="edFr" class="muted" style="margin-top:6px"></div>
      <div class="muted" style="margin-top:4px">식: R₀⁶ = C<sub>DA</sub>·τ<sub>D</sub>, c₀ = 3/(2π<sup>3/2</sup>R₀³), γ = (√π/2)·N<sub>A</sub>/c₀, η = √π γ e<sup>γ²</sup> erfc γ — dipole–dipole, acceptor 무작위 분포, donor 간 이동 없음 (Förster 1949; Inokuti & Hirayama 1965). Yb³⁺ 처럼 donor 간 energy migration 이 빠르면 이 값은 <b>하한</b>이다.</div>
    </details>`;
  $q(el,"#edIntN").onchange=e=>{ INT_MAXN=+e.target.value; INT_SEL=null; renderInter(); };
  el.querySelectorAll("button[data-int]").forEach(b=>b.onclick=()=>{ INT_SEL=INT_SEL===b.dataset.int?null:b.dataset.int; renderInter(); });
  const clr=$q(el,"#edIntClr"); if(clr) clr.onclick=()=>{ INT_SEL=null; renderInter(); };
  el.querySelectorAll("button[data-intadd]").forEach(b=>b.onclick=()=>{ const x=all[+b.dataset.intadd]; const md=M.cols[x.di], ma=M.cols[x.ai];
    const p={uid:uid(),kind:"ET",d:{c:md.c.uid,f:x.Di.id,t:x.Df.id},a:{c:ma.c.uid,f:x.Ai.id,t:x.Af.id},coop:false,type:null};
    if(S.procs.some(q=>q.kind==="ET"&&sameProc(q,p))) return; S.procs.push(p); S.sel=p.uid; save(); rerender(); });
  const calc=()=>{ const Cda=+$q(el,"#edCda").value, tau=+$q(el,"#edTauD").value/1e3, pct=+$q(el,"#edNApct").value, ns=+$q(el,"#edNsite").value; const o=$q(el,"#edFr");
    if(!(Cda>0&&tau>0&&pct>0&&ns>0)){ o.textContent="네 값을 모두 넣으면 계산합니다."; return; } const r=C.forsterEff(Cda,tau,pct/100*ns);
    o.innerHTML=`R₀ = <b>${fmt(r.R0*1e7,2)} nm</b> · c₀ = ${r.c0.toExponential(2)} cm⁻³ · γ = ${fmt(r.gamma,3)} · transfer efficiency η = <b>${fmt(100*r.eta,1)} %</b> · mean transfer rate ≈ η/((1−η)τ<sub>D</sub>) = ${r.eta<1?(r.eta/((1-r.eta)*tau)).toExponential(2):"∞"} s⁻¹`; };
  ["#edCda","#edTauD","#edNApct","#edNsite"].forEach(id=>$q(el,id).oninput=calc); calc();
}
function renderSrc(){
  const el=$q(ROOT,"#edSrcCard"); if(!el) return; const R=window.ED_DATA.refs;
  el.innerHTML=`<details class="ed-src"><summary>Data sources &amp; assumptions</summary><ul>
   <li><b>Ln³⁺ levels</b>: ${esc(R.ln3)} 준위선 = multiplet barycenter, 옅은 띠 = Stark 범위. 다른 호스트에서는 수십–수백 cm⁻¹ 이동.</li>
   <li><b>Multiphonon criterion</b>: p = (아래 준위와의 가장자리 간격)/ħω<sub>max</sub>. p&lt;${C.P_MPR} multiphonon 우세(옅은 선), ${C.P_MPR}–${C.P_RAD} 경쟁, ≥${C.P_RAD} 복사 우세(굵은 선) — 경험칙이며 호스트·온도에 따라 경계가 움직임.</li>
   <li><b>Resonance criterion</b>: donor·acceptor 전이의 Stark 창(가능한 Stark–Stark 조합 전체)이 겹치면 resonant. 중심 불일치가 ħω 보다 크면 ‘Stark 끝’. 겹치지 않으면 간격 δ 를 ħω 로 나눠 1-·2-·3–4 포논, 그 이상 non-resonant. endothermic(포논 흡수 필요)은 exp(−δ/kT) 를 따로 표시. 상온에서 위쪽 Stark 준위의 population은 Boltzmann 분포로 줄어들므로 창 가장자리 resonance는 실제보다 유리하게 보일 수 있음.</li>
   <li><b>Automatic type assignment</b>: acceptor 가 바닥에서 출발 → donor 도 바닥으로 가면 ET, 중간 준위로 가면 CR, sensitizer 이온(Yb³⁺·Nd³⁺·Ce³⁺)으로 되돌아가면 BET. acceptor 가 여기 상태에서 출발 → ETU. 같은 전이의 맞교환 → EM.</li>
   <li><b>dⁿ ions</b>: ${esc(R.ts)} 선폭은 근사(스핀 허용 ±1,200 cm⁻¹).</li>
   <li><b>Divalent · 5d</b>: ${esc(R.nist)} ${esc(R.dorenbos)} 5d 밴드는 호스트 의존 매개변수.</li>
   <li><b>Phonon energy</b>: ${esc(S.host.src||"—")}</li></ul></details>`;
}
function renderLegend(){
  const el=$q(ROOT,"#edLegend"); if(!el) return; const P=palette();
  const ln=(c,d,w)=>`<svg width="26" height="10"><line x1="1" x2="25" y1="5" y2="5" stroke="${c}" stroke-width="${w||2}" ${d?`stroke-dasharray="${d}"`:""}/></svg>`;
  el.innerHTML=`<span>${ln(P.ink,"",2)}Radiative-dominant</span><span>${ln(P.ink2,"",1.4)}Competing</span><span>${ln(P.faint,"",1.1)}Multiphonon-dominant</span>
   <span>${ln(C.wlColor(980,P.dark),"",2.6)}Absorption, resonant</span><span>${ln(C.wlColor(980,P.dark),"5 3",2.6)}Phonon-assisted</span><span>${ln(C.wlColor(980,P.dark),"2 3",1.6)}≥3 phonons</span>
   ${["ET","ETU","CR","BET","EM"].map(t=>`<span>${ln(P.T[t],"6 3",2)}${t}</span>`).join("")}<span><svg width="26" height="10"><path d="M1 5 Q4 1 7 5 T13 5 T19 5 T25 5" fill="none" stroke="${P.T.MPR}" stroke-width="1.5"/></svg>MPR</span>`;
}

/* ---------------- 내보내기 ---------------- */
function download(name,blob){ const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },500); }
function fileBase(){ return "EnergyDiagram_"+(S.formula||S.cols.map(c=>c.ion).join("-")||"ions").replace(/β/g,"b").replace(/α/g,"a").replace(/%/g,"pct").replace(/[^\w.+-]+/g,"_").replace(/_+/g,"_").slice(0,60); }
function svgString(){
  const cur=$q(ROOT,"#edSvg"); if(!cur||!MODEL.cols.length) return null;
  /* 내보내기는 테마와 무관하게 흰 바탕 · 편집 표시 없이 다시 그린다 */
  const tmp=document.createElement("div"); tmp.innerHTML=svgMarkup(palette(true),+cur.getAttribute("width"),false);
  const cl=tmp.firstElementChild;
  cl.querySelectorAll(".lvhit").forEach(n=>n.remove());
  cl.querySelectorAll("[tabindex]").forEach(n=>n.removeAttribute("tabindex"));
  cl.setAttribute("font-family","DM Sans, Noto Sans KR, Arial, sans-serif");
  return '<?xml version="1.0" encoding="UTF-8"?>\n'+new XMLSerializer().serializeToString(cl);
}
function exportSVG(){ const s=svgString(); if(s) download(fileBase()+".svg",new Blob([s],{type:"image/svg+xml"})); }
function exportPNG(){
  const s=svgString(); if(!s) return; const svg=$q(ROOT,"#edSvg"); const W=+svg.getAttribute("width"), H=+svg.getAttribute("height");
  const img=new Image(); const url=URL.createObjectURL(new Blob([s],{type:"image/svg+xml"}));
  img.onload=()=>{ const k=3, cv=document.createElement("canvas"); cv.width=W*k; cv.height=H*k; const g=cv.getContext("2d"); g.scale(k,k); g.drawImage(img,0,0); URL.revokeObjectURL(url);
    cv.toBlob(b=>download(fileBase()+".png",b),"image/png"); };
  img.onerror=()=>{ URL.revokeObjectURL(url); note("PNG 변환 실패 — SVG 로 내보내세요","warn"); };
  img.src=url;
}
function exportCSV(){
  const rows=[["#","type","auto_type","donor_ion","donor_from","donor_to","acceptor_ion","acceptor_from","acceptor_to","donor_E_cm1","acceptor_E_cm1","dE_cm1","delta_cm1","n_phonon","grade","endothermic","laser_nm","transition_nm","note"]];
  MODEL.procs.forEach(p=>{
    if(p.broken) return;
    if(p.kind==="ET") rows.push([p.n,p.type,p.auto,p.md.c.ion,p.Df.t,p.Dt.t,p.ma.c.ion,p.Af.t,p.At.t,Math.round(p.c.dw.E),Math.round(p.c.aw.E),Math.round(p.c.dE),Math.round(p.c.delta),p.c.nPh==null?"":p.c.nPh.toFixed(2),C.GRADE_EN[p.c.grade]+(p.c.edge?" (Stark edge)":""),p.c.endo?"yes":"no","","",p.p.note||""]);
    else if(p.kind==="PH") rows.push([p.n,p.type,p.type,p.m.c.ion,p.F.t,p.T.t,"","","",p.c?Math.round(C.nm2cm(p.nm)):"",p.c?Math.round(p.c.w.E):"",p.c?Math.round(p.c.dE):"",p.c?Math.round(p.c.delta):"",p.c&&p.c.nPh!=null?p.c.nPh.toFixed(2):"",p.c?C.GRADE_EN[p.c.grade]:"no laser",p.c&&p.c.endo?"yes":"no",p.nm||"",p.need?p.need.toFixed(1):"",p.p.note||""]);
    else rows.push([p.n,p.kind,p.kind,p.m.c.ion,p.F.t,p.T.t,"","","",Math.round(p.F.E-p.T.E),"","","","","","","",p.nm.toFixed(1),p.p.note||""]);
  });
  const csv="﻿"+`# ${S.formula||""} | host ${S.host.name||""} hw=${S.host.w} cm-1 T=${S.host.T} K | levels: Carnall et al. JCP 90, 3443 (1989) LaF3 barycenters\n`+rows.map(r=>r.map(v=>{ const s=String(v==null?"":v); return /[",\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s; }).join(",")).join("\n");
  download(fileBase()+"_processes.csv",new Blob([csv],{type:"text/csv"}));
}
function exportJSON(){ const o=Object.assign({},S); delete o.sel; download(fileBase()+".json",new Blob([JSON.stringify(o,null,1)],{type:"application/json"})); }
function importJSON(f){
  const r=new FileReader(); r.onload=()=>{ try{ const o=JSON.parse(r.result); if(!o||o.v!==1||!Array.isArray(o.cols)) throw new Error("형식이 다릅니다");
    o.cols=o.cols.filter(c=>C.IONS[c.ion]); S=Object.assign(DEF(),o,{mode:"view",sel:null}); PICK=null; save(); const f=$q(ROOT,"#edF"); if(f) f.value=S.formula; rerender(); note("도식을 불러왔습니다"); }
    catch(e){ note("불러오기 실패: "+e.message,"warn"); } }; r.readAsText(f);
}

/* 창 크기 변화 → 다시 그림 */
let _rz=null; window.addEventListener("resize",()=>{ clearTimeout(_rz); _rz=setTimeout(()=>{ if(ROOT&&ROOT.isConnected&&MODEL) drawSVG(); },150); });

window.ED_UI={render,rerender,_state:()=>S,_model:()=>MODEL,_applyFormula:(f)=>{applyFormula(f);rerender();},_loadPreset:loadPreset,PRESETS};
})();
