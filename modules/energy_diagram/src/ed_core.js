/* ============================================================
   INNO Analytics — Energy Diagram 핵심 로직 (UI 없음, node 로 시험 가능)

   원칙
   - 준위값은 전부 출처가 있는 표(ED_DATA)에서만 온다. 호스트 의존 값(5d 밴드, dⁿ 이온의
     Dq·B·C, 포논 에너지)은 '추정' 표지와 함께 사용자가 바꿀 수 있는 매개변수로 둔다.
   - 에너지 전달의 'resonant'은 두 전이의 Stark 창(가능한 Stark–Stark 조합 전체 폭)이 겹치는지로
     판정하고, 겹치지 않으면 남는 간격 δ 를 최대 포논 에너지 ħω 로 나눈 포논 수로 등급을 매긴다.
   - endothermic(δ 를 포논 흡수로 메워야 하는) 과정은 Boltzmann 인자 exp(−δ/kT)를 따로 보인다.
   ============================================================ */
(function(root){
"use strict";
const DATA = () => root.ED_DATA;
const KT_PER_K = 0.6950348;                 /* cm⁻¹ / K */
const HC_EV = 1239.84198;                    /* eV·nm */
const CM_PER_EV = 8065.544;

/* ---------- 단위 ---------- */
const nm2cm = nm => 1e7/nm;
const cm2nm = E => 1e7/E;
const cm2eV = E => E/CM_PER_EV;

/* ---------- 이온 목록 ---------- */
const LN3 = ["Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"];
const SUPN = {"1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹","0":"⁰"};
const sup = s => String(s).replace(/\d/g,c=>SUPN[c]);
/* 광학 활성 중심. kind: ln3 | ln2 | d */
const IONS = {};
LN3.forEach(el=>{ IONS[el+"3+"]={el,q:3,kind:"ln3",conf:"4f"+sup(DATA_N(el))}; });
function DATA_N(el){ const d=DATA(); return d&&d.ln3[el]?d.ln3[el].N:0; }
/* 2가 란탄족: 4f 준위는 NIST 자유 이온(호스트 이동 수십 cm⁻¹ 수준), 5d 는 호스트 의존 매개변수 */
IONS["Sm2+"]={el:"Sm",q:2,kind:"ln2",conf:"4f⁶",
  f:[["7F0",0],["7F1",293],["7F2",814],["7F3",1493],["7F4",2278],["7F5",3132],["7F6",4027],["5D0",14533]],
  fsrc:"NIST ASD, Sm III 자유 이온", band:{E:14800,stokes:0,need:false,label:"4f⁵5d",
    note:"SrF₂:Sm²⁺ 에서 4f⁵5d 가 ⁵D₀ 바로 위(675.6 nm 선) — 호스트에 따라 ⁵D₀ 아래로도 내려감(CaF₂: d–f 방출 710 nm)."},
  emit:[["5D0","7F1","SrF₂:Sm²⁺ 5D0→7F1 696.7 nm (MD)"],["5D0","7F0","5D0→7F0 683.8 nm (Oh 자리에서는 금지)"]]};
IONS["Eu2+"]={el:"Eu",q:2,kind:"ln2",conf:"4f⁷",f:[["8S7/2",0]],fsrc:"NIST ASD, Eu III",
  band:{E:24680,stokes:1130,need:false,label:"4f⁶5d",note:"기본값 CaF₂:Eu²⁺ 첫 f–d 흡수 3.06 eV · 방출 2.92 eV (arXiv:1709.00311 Table I). 다른 호스트는 값을 바꿀 것."},
  emit:[["5d","8S7/2","4f⁶5d → 4f⁷ (광대역)"]]};
IONS["Tm2+"]={el:"Tm",q:2,kind:"ln2",conf:"4f¹³",f:[["2F7/2",0],["2F5/2",8774]],fsrc:"NIST ASD, Tm III (호스트에서 수백 cm⁻¹ 이동 가능)",
  band:null, emit:[["2F5/2","2F7/2","²F₅/₂→²F₇/₂ (MF₂ 호스트 ≈1.12–1.18 µm)"]]};
IONS["Yb2+"]={el:"Yb",q:2,kind:"ln2",conf:"4f¹⁴",f:[["1S0",0]],fsrc:"4f¹⁴ 닫힌 껍질",
  band:{E:null,stokes:0,need:true,label:"4f¹³5d",note:"호스트 값 입력 필요 (자유 이온 최저 4f¹³5d = 33,386 cm⁻¹)."},
  emit:[["5d","1S0","4f¹³5d → 4f¹⁴"]]};
IONS["Ce3+"].band={E:40000,stokes:3000,need:false,label:"5d",
  note:"기본값: LaF₃ 첫 4f→5d 흡수 ≈250 nm (1차 출처 미확인, 근사). α-NaYF₄ 는 300 nm 흡수·332 nm 방출(Chin. Opt. Lett. 14, 061601 (2016))."};
IONS["Ce3+"].emit=[["5d","2F5/2","5d→²F₅/₂"],["5d","2F7/2","5d→²F₇/₂"]];
/* dⁿ 이온 — 팔면체 Tanabe–Sugano. 기본값은 표지된 참조 호스트 */
IONS["Cr3+"]={el:"Cr",q:3,kind:"d",n:3,conf:"3d³",B0:918,
  p:{Dq:1800,B:619,C:3248},preset:"Al₂O₃:Cr³⁺ (루비) — 표준 흡수대 18,000/24,400 cm⁻¹와 R₁ 694.3 nm 에서 역산",
  est:false};
IONS["Mn4+"]={el:"Mn",q:4,kind:"d",n:3,conf:"3d³",B0:1065,
  p:{Dq:2198,B:768,C:3456},preset:"K₂SiF₆:Mn⁴⁺ — ⁴T₂ 455 nm·²E ZPL 620.5 nm(Sijbom 2016)에서 Dq·C, B 는 C/B=4.5 가정",est:true};
IONS["Mn2+"]={el:"Mn",q:2,kind:"d",n:5,conf:"3d⁵",B0:960,
  p:{Dq:800,B:800,C:3200},preset:"일반 예시값 — 호스트의 Dq·B·C 를 입력할 것",est:true};
IONS["Ni2+"]={el:"Ni",q:2,kind:"d",n:8,conf:"3d⁸",B0:1041,
  p:{Dq:800,B:900,C:4050},preset:"일반 예시값 — 호스트의 Dq·B·C 를 입력할 것",est:true};
/* 광학적으로 불활성(가시–근적외 4f/3d 준위 없음) */
const INACTIVE = {La3:"4f⁰",Lu3:"4f¹⁴",Y3:"4d⁰",Sc3:"3d⁰",Ce4:"4f⁰ (전하이동 밴드만)",Gd3:null};

/* 원소별 선택 가능한 산화수 (주기율표) */
const OX = {Ce:[3],Pr:[3],Nd:[3],Pm:[3],Sm:[2,3],Eu:[2,3],Gd:[3],Tb:[3],Dy:[3],Ho:[3],Er:[3],Tm:[2,3],Yb:[2,3],
            Cr:[3],Mn:[2,4],Ni:[2]};
const ionKey = (el,q) => `${el}${q}+`;
const ionLabel = k => { const m=/^([A-Z][a-z]?)(\d)\+$/.exec(k); return m? m[1]+sup(m[2])+"⁺" : k; };

/* ---------- 호스트 ---------- */
/* 포논 에너지: host_phonons.json 의 출처값. key 는 정규화한 화학식 */
const HOSTS = [
  {key:"NaYF4",name:"β-NaYF₄",w:360,src:"PCCP 25, 11986 (2023) ← Wu et al., JPCC 122, 9611 (2018)",sub:3},
  {key:"NaGdF4",name:"β-NaGdF₄",w:360,est:true,src:"출처 미확보 — 동형 β-NaYF₄ 값(360) 차용",sub:3},
  {key:"NaLuF4",name:"β-NaLuF₄",w:360,est:true,src:"출처 미확보 — 동형 β-NaYF₄ 값(360) 차용",sub:3},
  {key:"NaYbF4",name:"β-NaYbF₄",w:360,est:true,src:"출처 미확보 — 동형 β-NaYF₄ 값 차용",sub:3},
  {key:"NaErF4",name:"β-NaErF₄",w:360,est:true,src:"출처 미확보 — 동형 β-NaYF₄ 값 차용",sub:3},
  {key:"LiYF4",name:"LiYF₄",w:460,src:"PCCP 25, 11986 (2023) ← Miller et al., JCP 52, 4172 (1970)",sub:3},
  {key:"LiLuF4",name:"LiLuF₄",w:460,est:true,src:"출처 미확보 — 동형 LiYF₄ 값 차용",sub:3},
  {key:"LaF3",name:"LaF₃",w:350,src:"PCCP 25, 11986 (2023) ← Weber, Phys. Rev. 157, 262 (1967)",sub:3},
  {key:"CaF2",name:"CaF₂",w:320,src:"PCCP 25, 11986 (2023) ← Mead & Wilkinson (1977), Raman T2g. IR LO 포논은 ≈460",sub:2},
  {key:"SrF2",name:"SrF₂",w:284,src:"PCCP 25, 11986 (2023) ← Mead & Wilkinson (1977)",sub:2},
  {key:"BaF2",name:"BaF₂",w:240,src:"PCCP 25, 11986 (2023) ← Mead & Wilkinson (1977)",sub:2},
  {key:"PbF2",name:"PbF₂",w:257,src:"PCCP 25, 11986 (2023) ← Zhang et al., Opt. Lett. 39, 3942 (2014)",sub:2},
  {key:"KPb2Cl5",name:"KPb₂Cl₅",w:203,src:"arXiv:0707.1434 (저포논 Er³⁺ 재료 리뷰)",sub:2},
  {key:"KPb2Br5",name:"KPb₂Br₅",w:138,src:"arXiv:0707.1434",sub:2},
  {key:"Y2O3",name:"Y₂O₃",w:600,src:"Light Sci. Appl. (2020) 10.1038/s41377-020-00346-z",sub:3},
  {key:"Y3Al5O12",name:"YAG",w:860,src:"Light Sci. Appl. (2020) 10.1038/s41377-020-00346-z",sub:3},
  {key:"YVO4",name:"YVO₄",w:950,src:"Ermeneux et al., PRB 61, 3915 (2000)",sub:3},
  {key:"YBO3",name:"YBO₃",w:1050,src:"Light Sci. Appl. (2020)",sub:3},
  {key:"Al2O3",name:"Al₂O₃",w:null,src:"",sub:3},
  {key:"K2SiF6",name:"K₂SiF₆",w:null,src:"",sub:4}
];
/* 음이온 계열 일반 범위 — 출처 호스트가 없을 때의 '추정' 기본값 */
const ANION_PH = {F:[400,"불화물 일반 범위 250–500 cm⁻¹ (추정)"],Cl:[250,"염화물 일반 범위 (추정)"],Br:[180,"브롬화물 (추정)"],
                  I:[150,"요오드화물 (추정)"],O:[700,"산화물 일반 범위 500–1100 cm⁻¹ (추정)"],S:[400,"황화물 (추정)"]};

/* ---------- 원소 정보 ---------- */
const VAL = {H:1,Li:1,Na:1,K:1,Rb:1,Cs:1,Ag:1,Tl:1,Be:2,Mg:2,Ca:2,Sr:2,Ba:2,Zn:2,Cd:2,Pb:2,Mn:2,
  B:3,Al:3,Ga:3,In:3,Sc:3,Y:3,La:3,Ce:3,Pr:3,Nd:3,Pm:3,Sm:3,Eu:3,Gd:3,Tb:3,Dy:3,Ho:3,Er:3,Tm:3,Yb:3,Lu:3,Bi:3,Cr:3,Fe:3,
  Si:4,Ge:4,Ti:4,Zr:4,Hf:4,Sn:4,Th:4,P:5,V:5,Nb:5,Ta:5,Mo:6,W:6,Ni:2,Co:2,Cu:2};
const ANIONS = {F:-1,Cl:-1,Br:-1,I:-1,O:-2,S:-2,Se:-2,N:-3};
const ACTIVE_EL = new Set(["Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Cr","Mn","Ni"]);
/* 치환 자리의 주 원소 — 활성 원소가 호스트 식에 들어 있을 때 호스트 키를 만드는 데 쓴다 */
const SITE_HOST = {3:["Y","Lu","Gd","La","Sc","Al","Ga","In","Bi"],2:["Ca","Sr","Ba","Mg","Zn","Pb","Cd"],4:["Si","Ge","Ti","Zr","Sn","Hf"]};

/* ---------- 화학식 파서 ---------- */
const SUB_MAP = {"₀":"0","₁":"1","₂":"2","₃":"3","₄":"4","₅":"5","₆":"6","₇":"7","₈":"8","₉":"9",
                 "⁰":"0","¹":"1","²":"2","³":"3","⁴":"4","⁵":"5","⁶":"6","⁷":"7","⁸":"8","⁹":"9","⁺":"+","⁻":"-","−":"-","–":"-","：":":","，":","};
function normalize(s){
  return String(s||"").replace(/[₀-₉⁰-⁹⁺⁻−–：，¹²³]/g,c=>SUB_MAP[c]||c).replace(/\s+/g," ").trim();
}
/* 원소 기호 목록(파서 검증용) */
const ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "+
 "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "+
 "Fr Ra Ac Th Pa U").split(" ");
const ELSET = new Set(ELEMENTS);

function parseComposition(str){
  /* "Na(Y0.78Yb0.2)F4" → {Na:1, Y:0.78, Yb:0.2, F:4} */
  const out={}; let s=str;
  const re=/\(([^()]*)\)(\d*\.?\d*)/;
  let guard=0;
  while(re.test(s) && guard++<8){
    s=s.replace(re,(m,inner,mul)=>{
      const k=mul?parseFloat(mul):1; const sub=parseComposition(inner); let t="";
      Object.entries(sub).forEach(([e,c])=>{ t+=e+(+(c*k).toFixed(6)); });
      return t;
    });
  }
  const tok=/([A-Z][a-z]?)(\d*\.?\d*)/g; let m, bad=[];
  const rest=s.replace(tok,"");
  if(rest.replace(/[\s\-·.]/g,"").length) bad.push(rest.trim());
  while((m=tok.exec(s))){
    if(!ELSET.has(m[1])){ bad.push(m[1]); continue; }
    out[m[1]]=(out[m[1]]||0)+(m[2]?parseFloat(m[2]):1);
  }
  if(bad.length) Object.defineProperty(out,"_bad",{value:bad,enumerable:false});
  return out;
}

function parseDopants(str){
  /* "20%Yb,2%Er" · "Yb20%/Er2%" · "Yb(20%)" · "Yb3+,Er3+" · "Eu2+" · "0.2Yb" */
  const out=[], warn=[];
  const parts=str.split(/[,;\/&]+|\s+and\s+/).map(x=>x.trim()).filter(Boolean);
  for(let p of parts){
    p=p.replace(/mol\s*%/gi,"%").replace(/at\.?\s*%/gi,"%");
    /* 원소 기호 뒤 전하와 농도를 분리: "Er3+2%" 같은 조합까지 허용 */
    const m=/^(\d+(?:\.\d+)?\s*%?)?\s*([A-Z][a-z]?)\s*(?:(\d?)([+-])|([+-])(\d))?\s*(?:\(?\s*(\d+(?:\.\d+)?)\s*(%?)\s*\)?)?\s*$/.exec(p);
    if(!m || !ELSET.has(m[2])){ warn.push(`해석 못함: "${p}"`); continue; }
    const el=m[2];
    let q=null;
    if(m[4]) q=(m[3]?+m[3]:1)*(m[4]==="+"?1:-1);
    else if(m[5]) q=(+m[6])*(m[5]==="+"?1:-1);
    let conc=null, pre=m[1]?m[1].trim():"", post=m[7]||"", raw=null, pct=false;
    if(pre){ raw=parseFloat(pre); pct=/%/.test(pre); }
    else if(post){ raw=parseFloat(post); pct=!!m[8]; }
    if(raw!=null) conc=pct?raw:(raw<=1?raw*100:raw);
    out.push({el,q,conc,_raw:raw,_pct:pct});
  }
  /* 한 식 안에서 % 없는 숫자는 한 가지로 읽는다: 모두 ≤1 이면 분율(×100), 하나라도 >1 이거나 % 가 붙은 값이 있으면 모두 %.
     (예: "Yb1.5,Er0.02" 를 1.5 % 와 2 % 로 섞어 읽지 않는다 → 1.5 %, 0.02 %) */
  const bare=out.filter(d=>d._raw!=null&&!d._pct);
  if(bare.length&&(bare.some(d=>d._raw>1)||out.some(d=>d._pct))){
    const chg=bare.filter(d=>d._raw<=1&&d._raw>0);
    bare.forEach(d=>d.conc=d._raw);
    if(chg.length) warn.push(`분율·% 혼용 — % 없는 숫자를 모두 % 로 읽음 (${chg.map(d=>d.el+" "+d._raw+" %").join(", ")}). 분율이면 모두 ≤1 로, %면 % 를 붙이세요`);
  }
  out.forEach(d=>{ delete d._raw; delete d._pct; });
  return {dopants:out,warn};
}

function findHost(key){ return HOSTS.find(h=>h.key===key)||null; }
function canonicalFormula(comp){
  /* 음이온을 뒤로, 나머지는 입력 순서 유지 */
  const cat=Object.keys(comp).filter(e=>!(e in ANIONS)), an=Object.keys(comp).filter(e=>e in ANIONS);
  const f=c=>{ const r=Math.round(c*1000)/1000; return Math.abs(r-1)<1e-9?"":String(r); };
  return cat.concat(an).map(e=>e+f(comp[e])).join("");
}

/* 화학식 전체: "β-NaYF4:20%Yb,2%Er@NaYF4" — 층은 @ 로 구분(코어@셸@…) */
function parseFormulaFull(input){
  const raw=normalize(input);
  const res={input:raw,layers:[],ions:[],warn:[]};
  if(!raw) return res;
  const parts=raw.split("@").map(x=>x.trim()).filter(Boolean);
  parts.forEach((L,li)=>{
    const one=parseLayer(L,li,res.warn);
    res.layers.push(one);
  });
  /* 이온 목록 — 층별로 합치되 같은 이온은 하나로(층 정보는 배열로) */
  const map=new Map();
  res.layers.forEach(ly=>ly.dopants.forEach(dp=>{
    const k=ionKey(dp.el,dp.q);
    const known=!!IONS[k];
    const inactive=!known && (INACTIVE[dp.el+dp.q]!==undefined);
    if(!map.has(k)) map.set(k,{key:k,el:dp.el,q:dp.q,known,inactive,layers:[],notes:[]});
    const it=map.get(k);
    it.layers.push({role:ly.role,conc:dp.conc,fromHost:dp.fromHost,qSrc:dp.qSrc});
    if(dp.note) it.notes.push(dp.note);
  }));
  /* Gd 호스트(NaGdF4 등): Gd³⁺ 는 UV 준위(⁶P)를 가진 활성 이온 — 에너지 이동 매개로 중요 */
  res.ions=[...map.values()];
  res.host=res.layers.length?res.layers[0]:null;
  return res;
}
function parseLayer(L,li,warn){
  let phase=null, s=L;
  const ph=/^(α|β|alpha|beta|hex|cubic)\s*-?\s*(?=[A-Z])/i.exec(s);
  if(ph){ const t=ph[1].toLowerCase(); phase=(t==="α"||t==="alpha"||t==="cubic")?"α":"β"; s=s.slice(ph[0].length); }
  const idx=s.indexOf(":");
  const hostStr=idx>=0?s.slice(0,idx):s, dopStr=idx>=0?s.slice(idx+1):"";
  const comp=parseComposition(hostStr);
  const role=li===0?"core":(li===1?"shell":`shell${li}`);
  const tag=li===0?"코어":`셸${li>1?li:""}`;
  if(comp._bad&&comp._bad.length) warn.push(`${tag} 식에서 모르는 부분: ${comp._bad.join(", ")}`);
  const siteOf=el=> (el in VAL)? VAL[el] : null;
  const cations=Object.keys(comp).filter(e=>!(e in ANIONS));
  const actInHost=cations.filter(e=>ACTIVE_EL.has(e));
  const siteSum={};
  cations.forEach(e=>{ const v=siteOf(e); if(v) siteSum[v]=(siteSum[v]||0)+comp[e]; });
  const dopants=[];
  actInHost.forEach(e=>{
    const v=siteOf(e)||3;
    const frac=comp[e]/(siteSum[v]||comp[e]);
    dopants.push({el:e,q:null,conc:+(frac*100).toFixed(3),fromHost:true,site:v});
  });
  const hostComp=Object.assign({},comp);
  actInHost.forEach(e=>{
    const v=siteOf(e)||3; const hostEl=cations.find(c=>!ACTIVE_EL.has(c)&&siteOf(c)===v);
    if(hostEl){ hostComp[hostEl]+=hostComp[e]; delete hostComp[e]; }
  });
  Object.keys(hostComp).forEach(e=>{ const r=Math.round(hostComp[e]); if(Math.abs(hostComp[e]-r)<0.021&&r>0) hostComp[e]=r; });
  let hostKey=canonicalFormula(hostComp);
  let hostInfo=findHost(hostKey);
  if(!hostInfo && actInHost.length){ const k2=canonicalFormula(comp); const h2=findHost(k2); if(h2){ hostInfo=h2; hostKey=k2; } }
  let qsum=0, known=Object.keys(comp).length>0;
  Object.entries(comp).forEach(([e,c])=>{ const v=(e in ANIONS)?ANIONS[e]:siteOf(e); if(v==null) known=false; else qsum+=v*c; });
  if(known && Math.abs(qsum)>0.05) warn.push(`${tag} ${hostKey}: 통상 원자가로 전하 합 ${qsum>0?"+":""}${+qsum.toFixed(2)} — 식을 확인하세요.`);
  const sites=Object.keys(siteSum).map(Number);
  if(dopStr){ const d=parseDopants(dopStr); d.warn.forEach(w=>warn.push(w)); d.dopants.forEach(x=>dopants.push(Object.assign({fromHost:false},x))); }
  const anion=Object.keys(comp).find(e=>e in ANIONS)||null;
  /* 농도 범위 점검 */
  dopants.forEach(dp=>{ if(dp.fromHost||dp.conc==null) return;
    if(!(dp.conc>0)) warn.push(`${tag} ${dp.el}: 농도 ${dp.conc} % — 0 이하`);
    else if(dp.conc>100) warn.push(`${tag} ${dp.el}: 농도 ${+dp.conc.toFixed(3)} % — 100 % 초과`); });
  dopants.forEach(dp=>{
    if(dp.q!=null){ dp.qSrc="입력"; return; }
    const ox=OX[dp.el]||[VAL[dp.el]||3];
    if(dp.fromHost&&dp.site){ dp.q=ox.includes(dp.site)?dp.site:ox[ox.length-1]; dp.qSrc="자리 원자가"; return; }
    if(dp.el==="Mn"){ dp.q=sites.includes(4)?4:2; dp.qSrc="치환 자리";
      if(!sites.includes(4)&&!sites.includes(2)) dp.note="Mn 이 들어갈 자리의 원자가가 분명하지 않음 — 산화수를 확인하세요"; return; }
    if(dp.el==="Ni"){ dp.q=2; dp.qSrc="기본"; return; }
    if(dp.el==="Cr"){ dp.q=3; dp.qSrc="기본"; return; }
    const div=sites.length && !sites.includes(3) && sites.includes(2);
    if(dp.el==="Eu" && div){ dp.q=2; dp.qSrc="2가 자리"; dp.note="2가 양이온 자리 — Eu²⁺ 로 두었음. 산화 분위기 합성이면 Eu³⁺(전하 보상)일 수 있음"; return; }
    if(["Sm","Yb","Tm"].includes(dp.el) && div) dp.note=`2가 양이온 자리 — 환원 조건이면 ${dp.el}²⁺ 가능`;
    if(div && dp.el!=="Eu" && ox.includes(3)) dp.note=(dp.note?dp.note+". ":"")+"3가 이온이 2가 자리에 들어가 전하 보상 결함이 생김 → 자리 대칭·준위 분열이 LaF₃ 와 다를 수 있음";
    dp.q=ox.includes(3)?3:ox[0]; dp.qSrc="기본";
  });
  /* 같은 자리(같은 원자가)에 들어가는 도펀트 합 */
  const sumQ={}; dopants.forEach(dp=>{ if(dp.fromHost||!(dp.conc>0)||dp.conc>100) return; const k=dp.q==null?3:dp.q; (sumQ[k]=sumQ[k]||[]).push(dp); });
  Object.values(sumQ).forEach(L=>{ const t=L.reduce((a,d)=>a+d.conc,0); if(L.length>1&&t>100+1e-9) warn.push(`${tag} 같은 자리 도펀트 합 ${+t.toFixed(3)} % (${L.map(d=>d.el).join("+")}) — 100 % 초과`); });
  return {index:li,role,tag,phase,hostRaw:hostStr,comp,hostKey,hostInfo,anion,sites,dopants};
}

/* ---------- Tanabe–Sugano 보간 ---------- */
const TS_NAMES = {
  d3:{"4A2#0":"4A2","4T2#0":"4T2","2E#0":"2E","2T1#0":"2T1","4T1#0":"4T1(F)","2T2#0":"2T2","4T1#1":"4T1(P)"},
  d5:{"6A1#0":"6A1","4T1#0":"4T1(G)","4T2#0":"4T2(G)","4A1#0":"4A1(G)","4E#0":"4E(G)","4T2#1":"4T2(D)","4E#1":"4E(D)","4T1#1":"4T1(P)","2T2#0":"2T2(I)","2A2#0":"2A2","2T1#0":"2T1(I)"},
  d8:{"3A2#0":"3A2","3T2#0":"3T2","1E#0":"1E","3T1#0":"3T1(F)","1A1#0":"1A1","1T2#0":"1T2","3T1#1":"3T1(P)"}};
function tsEval(n,B,C,Dq){
  const T=DATA().ts, g=T.g["d"+n]; if(!g) return [];
  const x=Math.max(T.DQ[0],Math.min(T.DQ[T.DQ.length-1],Dq/B));
  const y=Math.max(T.CB[0],Math.min(T.CB[T.CB.length-1],C/B));
  const fi=(arr,v)=>{ let i=0; while(i<arr.length-2 && arr[i+1]<v) i++; return [i,(v-arr[i])/(arr[i+1]-arr[i])]; };
  const [qi,qt]=fi(T.DQ,x), [ci,ct]=fi(T.CB,y);
  const out=[];
  g.labels.forEach((lab,k)=>{
    const E=g.E[k];
    const v=(1-ct)*((1-qt)*E[ci][qi]+qt*E[ci][qi+1])+ct*((1-qt)*E[ci+1][qi]+qt*E[ci+1][qi+1]);
    /* 기울기 dE/dDq (선폭 근사용) */
    const v2=(1-ct)*E[ci][qi+1]+ct*E[ci+1][qi+1], v1=(1-ct)*E[ci][qi]+ct*E[ci+1][qi];
    const slope=(v2-v1)/100/(T.DQ[qi+1]-T.DQ[qi]);   /* (E/B)/(Dq/B) = dE/dDq */
    out.push({lab,E:v/100*B,slope});
  });
  const clamped=(Dq/B!==x)||(C/B!==y);
  const m=Math.min(...out.map(o=>o.E)); out.forEach(o=>o.E-=m);
  return {levels:out.sort((a,b)=>a.E-b.E),clamped,x,y};
}

/* ---------- 준위 만들기 ---------- */
function termHTML(t){
  /* "4I15/2" → ⁴I₁₅/₂ 형태의 HTML, "4T1(F)" → ⁴T₁(F), "5d" 등은 그대로 */
  const m=/^(\d+)([A-Z])(\d+(?:\/\d)?|\d?)(.*)$/.exec(t);
  if(!m) return escHTML(t);
  const [_,mult,L,J,rest]=m;
  return `<sup>${mult}</sup>${L}${J?`<sub>${J}</sub>`:""}${escHTML(rest)}`;
}
function termText(t){
  const SUB={"0":"₀","1":"₁","2":"₂","3":"₃","4":"₄","5":"₅","6":"₆","7":"₇","8":"₈","9":"₉","/":"/"};
  const m=/^(\d+)([A-Z])(\d+(?:\/\d)?|\d?)(.*)$/.exec(t);
  if(!m) return t;
  return sup(m[1])+m[2]+m[3].replace(/[\d\/]/g,c=>SUB[c])+m[4];
}
function escHTML(s){ return String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c])); }

/* col: {ion:"Tm3+", p:{Dq,B,C}, band:{E,stokes,width}} */
function levelsFor(col){
  const def=IONS[col.ion]; if(!def) return [];
  const out=[];
  if(def.kind==="ln3"){
    const d=DATA().ln3[def.el];
    d.lv.forEach((r,i)=>out.push({id:r[0],t:r[1],E:r[2],lo:Math.min(r[3],r[2]),hi:Math.max(r[4],r[2]),src:r[5],dom:r[6],comp:r[7],note:r[8],kind:"f",ground:i===0}));
  } else if(def.kind==="ln2"){
    def.f.forEach((r,i)=>out.push({id:r[0],t:r[0],E:r[1],lo:r[1],hi:r[1],src:"n",kind:"f",ground:i===0,note:def.fsrc}));
  } else if(def.kind==="d"){
    const p=Object.assign({},def.p,col.p||{});
    const r=tsEval(def.n,p.B,p.C,p.Dq);
    const names=TS_NAMES["d"+def.n];
    r.levels.forEach((o,i)=>{
      const nm=names[o.lab]||o.lab;
      /* 선폭 근사: 바닥과 기울기가 같은(스핀 뒤집기) 준위는 좁고, 다르면 넓다 */
      const g0=r.levels[0].slope, broad=Math.abs(o.slope-g0)>1.5;
      const hw=i===0?0:(broad?1200:120);
      out.push({id:nm,t:nm,E:o.E,lo:o.E-hw,hi:o.E+hw,src:"ts",kind:"d",ground:i===0,broad,slope:o.slope-g0,
        note:(broad?"스핀 허용·전자 배치 변화 → 넓은 띠(±1,200 cm⁻¹ 근사)":"스핀 뒤집기(같은 배치) → 좁은 선")+(r.clamped?" · 격자 범위 밖 값은 가장자리로 고정됨":"")});
    });
  }
  /* 5d / 4fⁿ⁻¹5d 밴드 */
  const b=def.band?Object.assign({},def.band,col.band||{}):null;
  if(b && b.E!=null && b.E>0){
    const w=b.width||4000;
    out.push({id:"5d",t:b.label||"5d",E:b.E,lo:b.E-600,hi:b.E+w,src:"band",kind:"band",ground:false,stokes:b.stokes||0,
      note:(col.band&&col.band.E!=null?"사용자 입력값. ":"")+(b.note||"")});
  }
  out.sort((a,b)=>a.E-b.E);
  return out;
}

/* ---------- 간격·multiphonon ---------- */
const P_MPR=5, P_RAD=8;     /* p<5: multiphonon 우세, 5≤p<8: 경쟁, p≥8: 복사 우세 (경험칙) */
function gapInfo(levels,i,hw){
  if(i<=0) return {gap:null,p:null,regime:"ground"};
  const L=levels[i];
  /* 바로 아래 준위(띠 준위 무시하지 않음) */
  let j=i-1; while(j>0 && levels[j].E>=L.E) j--;
  const lower=levels[j];
  const gapEdge=Math.max(0,L.lo-lower.hi), gapBary=L.E-lower.E;
  const p=hw?gapEdge/hw:null;
  const regime=p==null?"unknown":(p<P_MPR?"mpr":(p<P_RAD?"mixed":"rad"));
  return {gap:gapEdge,gapBary,p,regime,lowerId:lower.id};
}
function annotate(levels,hw,kind){
  levels.forEach((L,i)=>Object.assign(L,gapInfo(levels,i,hw)));
  /* dⁿ 이온의 넓은 띠 준위: 에너지 간격 법칙이 아니라 계간 교차·Stokes 이완으로 가장 낮은 들뜬 준위에
     빠르게 모인다(Kasha 형). 가장 낮은 들뜬 준위만 방출 준위로 둔다. */
  if(kind==="d"){ levels.forEach((L,i)=>{ if(i===0) return; L.regime= i===1?"rad":"mpr"; L.p=null; L.lowerId= i===1?levels[0].id:levels[i-1].id; }); }
  /* multiphonon으로 빠르게 흘러내리는 끝(연쇄) */
  levels.forEach((L,i)=>{
    let k=i, guard=0;
    while(k>0 && levels[k].regime==="mpr" && guard++<60){ const lid=levels[k].lowerId; k=levels.findIndex(x=>x.id===lid); if(k<0){k=0;break;} }
    L.cascade=levels[k].id;
    /* population 흐름: multiphonon 우세와 경쟁 준위를 모두 따라 내려간 끝 (population 가능성 판단용) */
    let r=i; guard=0;
    while(r>0 && (levels[r].regime==="mpr"||levels[r].regime==="mixed") && guard++<60){ const lid=levels[r].lowerId; const nr=levels.findIndex(x=>x.id===lid); if(nr<0) break; r=nr; }
    L.relax=[]; let q=i; guard=0;
    while(guard++<60){ L.relax.push(levels[q].id); if(q===r||q<=0) break; const lid=levels[q].lowerId; const nq=levels.findIndex(x=>x.id===lid); if(nq<0||nq===q) break; q=nq; }
  });
  return levels;
}

/* ---------- 전이 창 ---------- */
function win(up,down,mult){
  const k=mult||1;
  return {E:k*(up.E-down.E),lo:k*(up.lo-down.hi),hi:k*(up.hi-down.lo)};
}
const GRADES=["resonant","1-포논","2-포논","3–4 포논","non-resonant"];
const GRADE_EN=["resonant","1-phonon","2-phonon","3-4 phonon","off-resonant"];
function gradeOf(delta,hw){
  if(delta<=0) return 0;
  if(!hw) return delta<=200?1:(delta<=600?3:4);
  if(delta<=hw) return 1; if(delta<=2*hw) return 2; if(delta<=4*hw) return 3; return 4;
}

/* 에너지 전달: donor Di→Df (아래로), acceptor Ai→Af (위로) */
const SENSITIZERS=new Set(["Yb3+","Nd3+","Ce3+"]);
const MIGRATORS=new Set(["Yb3+","Gd3+"]);
function classifyET(et,ctx){
  /* et: {d:{col,from,to}, a:{col,from,to}, coop} — from/to 는 level 객체 */
  const {d,a}=et, hw=ctx.hw, kT=KT_PER_K*(ctx.T||295);
  const mult=et.coop?2:1;
  const dw=win(d.from,d.to,mult), aw=win(a.to,a.from);
  const dE=dw.E-aw.E;       /* >0: donor 가 남는 에너지 → 포논 방출(발열) */
  let delta;
  if(dw.hi<aw.lo) delta=aw.lo-dw.hi;        /* donor 가 모자람 → 포논 흡수 필요(endothermic) */
  else if(aw.hi<dw.lo) delta=dw.lo-aw.hi;   /* donor 가 남음 → 포논 방출 */
  else delta=0;
  const endo = delta>0 && dw.hi<aw.lo;
  const grade=gradeOf(delta,hw);
  const nPh = hw? delta/hw : null;
  const boltz = endo? Math.exp(-delta/kT) : 1;
  /* 종류 */
  const dIon=d.ion, aIon=a.ion;
  const aGround=a.from.ground, dGround=d.to.ground;
  let type;
  if(et.coop) type="CET";
  else if(dIon===aIon && d.from.id===a.to.id && d.to.id===a.from.id) type="EM";
  else if(aGround && SENSITIZERS.has(aIon) && !SENSITIZERS.has(dIon)) type="BET";
  else if(aGround && dGround) type="ET";
  else if(aGround && !dGround) type="CR";
  else type="ETU";
  /* 창은 겹치지만 중심 불일치가 ħω 보다 크면 Stark 가장자리끼리만 겹친다 — 따로 표시 */
  const edge = grade===0 && hw && Math.abs(dE)>hw;
  return {type,dw,aw,dE,delta,endo,grade,edge,gradeName:edge?"resonant(Stark 끝)":GRADES[grade],nPh,boltz,kT};
}
const TYPE_INFO = {
  ET:{name:"에너지 전달 (ET)",desc:"여기된 donor 가 바닥상태로 돌아가며 바닥상태 acceptor 를 여기"},
  ETU:{name:"에너지 전달 upconversion (ETU)",desc:"이미 여기된 acceptor 를 한 단계 더 올림"},
  CR:{name:"교차 이완 (CR)",desc:"여기 donor 가 중간 준위로 내려가며 바닥상태 acceptor 를 중간 준위로 — 에너지를 나눠 가짐"},
  EM:{name:"에너지 이동 (EM)",desc:"같은 이온 사이의 동일 전이 — resonance가 자동 성립, 에너지가 격자를 따라 이동"},
  BET:{name:"역에너지 전달 (BET)",desc:"활성 이온에서 sensitizer 이온으로 되돌아가는 전달"},
  CET:{name:"협동 에너지 전달 (CET)",desc:"여기된 donor 두 개가 함께 acceptor 하나를 여기 (에너지 2배)"},
  GSA:{name:"바닥상태 흡수 (GSA)",desc:""}, ESA:{name:"여기상태 흡수 (ESA)",desc:""},
  RAD:{name:"복사 전이 (방출)",desc:""}, MPR:{name:"multiphonon 이완 (MPR)",desc:""}
};

/* 광흡수: laser 광자(cm⁻¹) 와 전이 창의 간격 */
function classifyPhoton(from,to,Elas,ctx){
  const w=win(to,from), hw=ctx.hw, kT=KT_PER_K*(ctx.T||295);
  let delta=0, endo=false;
  if(Elas<w.lo){ delta=w.lo-Elas; endo=true; }        /* 광자가 모자람 → anti-Stokes(포논 흡수) */
  else if(Elas>w.hi){ delta=Elas-w.hi; }              /* 광자가 남음 → 포논 방출 */
  const grade=gradeOf(delta,hw);
  const edge = grade===0 && hw && Math.abs(Elas-w.E)>hw;
  return {type:from.ground?"GSA":"ESA",w,delta,endo,grade,edge,gradeName:edge?"resonant(Stark 끝)":GRADES[grade],nPh:hw?delta/hw:null,
          boltz:endo?Math.exp(-delta/kT):1,dE:Elas-w.E};
}

/* 주어진 여기 파장들에 대한 GSA/ESA 후보 */
/* 광흡수 후보 — GSA 에서 시작해 실제로 population이 생길 수 있는 준위에서만 ESA 를 찾는다.
   population: GSA·ESA 도착 준위가 (multiphonon 우세·경쟁 준위를 따라) 흘러내린 준위들 + 가장 낮은 metastable 준위
   + 사용자가 그린 과정의 도착 준위(extra). ESA 는 두 단계까지. */
function photonCandidates(col,levels,lasers,ctx,extra){
  const out=[]; const byId={}; levels.forEach(l=>byId[l.id]=l);
  const best=(s,El,limit,nm)=>{
    const cand=[];
    levels.forEach(f=>{ if(f.E<=s.E+100) return; const c=classifyPhoton(s,f,El,ctx); cand.push({nm,from:s,to:f,c}); });
    cand.sort((x,y)=>x.c.delta-y.c.delta||Math.abs(x.c.dE)-Math.abs(y.c.dE));
    const keep=cand.filter(x=>x.c.grade<=1);
    if(!keep.length && cand.length && cand[0].c.grade<=limit) keep.push(cand[0]);
    return keep;
  };
  const lowMeta=levels.find(L=>!L.ground && L.regime==="rad" && L.kind!=="band");
  lasers.forEach(nm=>{
    const El=nm2cm(nm); const seen=new Set();
    const g=best(levels[0],El,3,nm); g.forEach(k=>out.push(k));
    let front=new Set();
    const reach=l=>{ (l.relax||[l.id]).forEach(id=>{ const L=byId[id]; if(L && !L.ground) front.add(id); }); };
    g.forEach(k=>reach(k.to));
    if(lowMeta) front.add(lowMeta.id);
    (extra||[]).forEach(id=>{ if(byId[id]) reach(byId[id]); });
    for(let gen=0; gen<2; gen++){
      const next=new Set();
      front.forEach(id=>{
        if(seen.has(id)) return; seen.add(id);
        const s=byId[id]; if(!s||s.ground) return;
        if(s.regime!=="rad" && s.kind!=="band") return;       /* 수명이 짧은(경쟁·multiphonon) 준위는 ESA 출발로 보지 않음 */
        best(s,El,1,nm).forEach(k=>{ if(k.c.grade>1) return; out.push(k); if(k.c.grade===0) (k.to.relax||[k.to.id]).forEach(x=>next.add(x)); });
      });
      front=next;
    }
  });
  /* 중복 제거 */
  const u=new Set(); return out.filter(x=>{ const k=x.nm+"|"+x.from.id+"|"+x.to.id; if(u.has(k)) return false; u.add(k); return true; });
}

/* photon avalanche(PA) 루프 후보 — 문헌의 필요 조건(Joubert 1999; Lee 2021)을 준위 구조만으로 점검
   ① GSA 가 약하거나 non-resonant  ② metastable 준위 M 에서의 ESA 가 resonant  ③ ESA 도착 준위가 흘러내린 준위 H 에서
   H + 바닥 → M' + M'' cross-relaxation이 있고 두 결과가 모두 M 으로 cascade(한 광자로 M 이 둘) */
function paLoops(colIdx,levels,lasers,ctx){
  const res=[];
  const byId=id=>levels.find(l=>l.id===id);
  const ground=levels[0];
  lasers.forEach(nm=>{
    const El=nm2cm(nm);
    /* GSA 최선 */
    let gBest=null;
    levels.forEach(f=>{ if(f.ground) return; const c=classifyPhoton(ground,f,El,ctx); if(!gBest||c.delta<gBest.c.delta) gBest={to:f,c}; });
    levels.forEach(M=>{
      if(M.ground || M.regime!=="rad" || M.kind==="band") return;
      levels.forEach(U=>{
        if(U.E<=M.E) return;
        const e=classifyPhoton(M,U,El,ctx); if(e.grade>1) return;
        const H=byId(U.cascade)||U;
        if(H.id===M.id||H.ground) return;
        /* cross-relaxation H→X, 바닥→Y (X,Y 는 흘러내려 M) */
        levels.forEach(X=>{
          if(X.E>=H.E||X.ground) return;
          if((byId(X.cascade)||X).id!==M.id) return;
          levels.forEach(Y=>{
            if(Y.ground||Y.E>=H.E) return;
            if((byId(Y.cascade)||Y).id!==M.id) return;
            const cr=classifyET({d:{ion:"x",from:H,to:X},a:{ion:"x",from:ground,to:Y}},ctx);
            if(cr.grade>2) return;
            const ratio=gBest? gBest.c.grade : 4;
            if(ratio<2) return;              /* GSA 가 resonant이면 avalanche의 전제(약한 GSA)가 아니다 */
            res.push({nm,M,U,H,X,Y,esa:e,cr,gsa:gBest,
              gsaWeak: ratio>=2,
              score:(ratio>=2?0:10)+e.grade*2+cr.grade});
          });
        });
      });
    });
  });
  res.sort((a,b)=>a.score-b.score);
  /* 같은 (nm,M,H) 중복 제거 */
  const seen=new Set();
  return res.filter(r=>{ const k=r.nm+"|"+r.M.id+"|"+r.U.id+"|"+r.X.id+"|"+r.Y.id; if(seen.has(k)) return false; seen.add(k); return true; });
}

/* 에너지 전달 경로 제안 */
function suggestET(cols,levelsByCol,ctx,opt){
  opt=opt||{}; const maxGrade=opt.maxGrade==null?2:opt.maxGrade;
  const out=[];
  const Emax=opt.Emax||1e9;
  cols.forEach((cd,i)=>cols.forEach((ca,j)=>{
    const LD=levelsByCol[i], LA=levelsByCol[j];
    const donors=LD.filter(L=>!L.ground && (L.regime==="rad"||L.regime==="mixed"||L.kind==="band") && L.E<=Emax);
    const accStart=LA.filter(L=>L.ground || L.regime==="rad");
    donors.forEach(Di=>{
      LD.forEach(Df=>{
        if(Df.E>=Di.E-200) return;
        accStart.forEach(Ai=>{
          LA.forEach(Af=>{
            if(Af.E<=Ai.E+200 || Af.E>Emax) return;
            const same=cd.ion===ca.ion;
            const mirror=same && Di.id===Af.id && Df.id===Ai.id;
            if(mirror && !(MIGRATORS.has(cd.ion) && Df.ground)) return;
            if(mirror && i!==j) return;   /* 이동은 같은 이온 열 하나로 대표 */
            const c=classifyET({d:{ion:cd.ion,from:Di,to:Df},a:{ion:ca.ion,from:Ai,to:Af}},ctx);
            if(c.grade>maxGrade) return;
            /* 같은 이온·같은 열의 CR/ETU 는 i===j 에서만(중복 방지) */
            if(same && i!==j) return;
            out.push({di:i,ai:j,Di,Df,Ai,Af,c,
              /* 점수: 등급 → endothermic → population(높은 준위일수록 드묾) → 불일치 순으로 불리 */
              score:c.grade*1000+(c.edge?400:0)+(c.endo?300:0)+Di.E*0.02+Ai.E*0.05+Math.abs(c.dE)*0.1});
          });
        });
      });
    });
  }));
  out.sort((a,b)=>a.score-b.score);
  return out.slice(0,opt.limit||60);
}

/* ---------- 이온 간 상호작용 (서로 다른 이온) ----------
   서로 다른 두 이온 A, B 에 대해 A 의 여기 준위(donor)가 B 의 바닥·metastable 준위(acceptor)를 올리는 모든 조합을
   에너지 불일치로 분류한다. 양방향(A→B, B→A: back transfer 포함)을 모두 본다.
   정량하는 것: Stark 창 간격 δ, 필요한 포논 수 N = δ/ħω, endothermic 이면 detailed balance exp(−δ/kT).
   정량하지 않는 것: 두 전이의 세기(스펙트럼 겹침의 크기 = JO 선세기·선모양)와 이온 간 거리 분포.
   이 둘이 절대 속도를 몇 자릿수씩 바꾸므로(예: Yb→Tm 은 δ≈1,600 cm⁻¹ 로 비공명인데도 효율적) 여기서 내는 값은 확률이 아니다.
   절대 확률은 측정된 microparameter C_DA 가 있을 때만 forsterEff 로 계산한다. */
function interIon(cols,levelsByCol,ctx,opt){
  opt=opt||{}; const Emax=opt.Emax||1e9, maxGrade=opt.maxGrade==null?4:opt.maxGrade; const hw=ctx.hw||null, kT=KT_PER_K*(ctx.T||295);
  const out=[];
  cols.forEach((cd,i)=>cols.forEach((ca,j)=>{
    if(i===j||cd.ion===ca.ion) return;
    const LD=levelsByCol[i], LA=levelsByCol[j];
    const donors=LD.filter(L=>!L.ground&&(L.regime==="rad"||L.regime==="mixed"||L.kind==="band")&&L.E<=Emax);
    const accStart=LA.filter(L=>L.ground||L.regime==="rad");
    donors.forEach(Di=>LD.forEach(Df=>{ if(Df.E>=Di.E-200) return;
      accStart.forEach(Ai=>LA.forEach(Af=>{ if(Af.E<=Ai.E+200||Af.E>Emax) return;
        const c=classifyET({d:{ion:cd.ion,from:Di,to:Df},a:{ion:ca.ion,from:Ai,to:Af}},ctx); if(c.grade>maxGrade) return;
        out.push({di:i,ai:j,Di,Df,Ai,Af,c,N:hw?c.delta/hw:null,boltz:c.endo?Math.exp(-c.delta/kT):1}); })); }));
  }));
  /* 같은 (donor 출발, acceptor 출발) 칸 안에서는 δ 가 작은 순 */
  out.sort((a,b)=>a.c.delta-b.c.delta||Math.abs(a.c.dE)-Math.abs(b.c.dE));
  return out;
}
/* Förster 정적 평균 (dipole–dipole, acceptor 무작위 분포, donor 간 이동 없음):
   R0⁶ = C_DA·τ_D, c0 = 3/(2π^{3/2} R0³), γ = (√π/2)·N_A/c0, 전달 효율 η = √π γ e^{γ²} erfc(γ).
   Förster, Z. Naturforsch. 4a, 321 (1949); Inokuti & Hirayama, J. Chem. Phys. 43, 1978 (1965). */
function erfcx(x){ if(x<0) return 2*Math.exp(x*x)-erfcx(-x);
  if(x<2){ const t=1/(1+0.5*x); const y=t*Math.exp(-x*x-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277))))))))); return y*Math.exp(x*x); }
  let f=x; for(let k=60;k>=1;k--) f=x+(k/2)/f; return 1/(Math.sqrt(Math.PI)*f); }
function forsterEff(CDA,tauD,NA){ if(!(CDA>0&&tauD>0&&NA>0)) return null; const R0=Math.pow(CDA*tauD,1/6); const c0=3/(2*Math.pow(Math.PI,1.5)*R0*R0*R0);
  const g=Math.sqrt(Math.PI)/2*NA/c0; return {R0,c0,gamma:g,eta:Math.sqrt(Math.PI)*g*erfcx(g)}; }

/* ---------- 파장 → 색 (가시광은 스펙트럼 색, UV·NIR·IR 은 구분 색) ---------- */
function wlColor(nm,dark){
  let r=0,g=0,b=0;
  if(nm<380){ return dark?"#b9a3ff":"#5b3fc4"; }                /* UV */
  if(nm>1000) return dark?"#c9a58a":"#7a5a44";                   /* IR */
  if(nm>750) return dark?"#e0828a":"#8f2d38";                    /* NIR */
  if(nm<440){ r=-(nm-440)/(440-380); b=1; }
  else if(nm<490){ g=(nm-440)/(490-440); b=1; }
  else if(nm<510){ g=1; b=-(nm-510)/(510-490); }
  else if(nm<580){ r=(nm-510)/(580-510); g=1; }
  else if(nm<645){ r=1; g=-(nm-645)/(645-580); }
  else { r=1; }
  let f=1; if(nm<420) f=0.3+0.7*(nm-380)/40; else if(nm>700) f=0.3+0.7*(750-nm)/50;
  /* 흰 바탕 대비 확보: 라이트는 조금 어둡게, 다크는 밝게 */
  const k=dark?1.0:0.78;
  const c=v=>Math.round(255*Math.pow(Math.max(0,Math.min(1,v*f)),0.8)*k + (dark?40:0));
  return `rgb(${Math.min(255,c(r))},${Math.min(255,c(g))},${Math.min(255,c(b))})`;
}
function bandName(nm){
  if(nm<280) return "UV-C"; if(nm<315) return "UV-B"; if(nm<400) return "UV-A";
  if(nm<450) return "보라"; if(nm<495) return "파랑"; if(nm<570) return "초록"; if(nm<590) return "노랑";
  if(nm<620) return "주황"; if(nm<750) return "빨강"; if(nm<1400) return "근적외(NIR-I/II)"; if(nm<3000) return "단파 적외(SWIR)";
  return "mid-IR";
}

/* 문헌 대표 전이(방출·흡수) */
function curatedFor(ionKey){
  const def=IONS[ionKey]; if(!def) return [];
  const extra=def.emit?def.emit.map(e=>({up:[e[0]],lo:e[1],kind:"em",nm:null,src:e[2]})):[];
  if(def.kind==="ln3") return DATA().trans.filter(t=>t[0]===def.el).map(t=>({up:t[1],lo:t[2],kind:t[3],nm:t[4],src:t[5]||"Carnall 1989 바리센터로 계산한 λ; 대표값은 문헌 관례"})).concat(extra);
  if(def.emit) return def.emit.map(e=>({up:[e[0]],lo:e[1],kind:"em",nm:null,src:e[2]}));
  return [];
}
/* dⁿ 이온의 대표 전이(장 세기에 따라 달라짐) */
function dMain(col,levels){
  const def=IONS[col.ion]; const f=id=>levels.find(l=>l.id===id);
  const out=[]; const g=levels[0];
  if(def.n===3){
    const E2=f("2E"), T2=f("4T2");
    if(E2&&T2){
      if(E2.E<T2.E-1500) out.push({up:"2E",lo:g.id,kind:"em",note:"강한 장(²E < ⁴T₂): 좁은 R선 스핀 금지 방출"});
      else out.push({up:"4T2",lo:g.id,kind:"em",note:"약한 장(⁴T₂ ≲ ²E): 넓은 띠 방출 — 실제 봉우리는 Stokes 이동만큼 장파장"});
    }
    ["4T2","4T1(F)","4T1(P)"].forEach(t=>{ if(f(t)) out.push({up:t,lo:g.id,kind:"abs",note:"스핀 허용 흡수대"}); });
  } else if(def.n===5){
    if(g.id==="6A1"){ if(f("4T1(G)")) out.push({up:"4T1(G)",lo:"6A1",kind:"em",note:"⁴T₁(G)→⁶A₁ 스핀 금지 방출 — Stokes 이동 큼"});
      ["4T1(G)","4T2(G)","4A1(G)","4T2(D)"].forEach(t=>{ if(f(t)) out.push({up:t,lo:"6A1",kind:"abs",note:"스핀 금지 흡수(약함)"}); }); }
  } else if(def.n===8){
    if(f("3T2")) out.push({up:"3T2",lo:g.id,kind:"em",note:"³T₂→³A₂ 근적외 넓은 띠 — 실제 봉우리는 Stokes 이동만큼 장파장"});
    ["3T2","3T1(F)","3T1(P)"].forEach(t=>{ if(f(t)) out.push({up:t,lo:g.id,kind:"abs",note:"스핀 허용 흡수대"}); });
  }
  return out;
}

/* 이온별 대표 여기 파장 (조성에 sensitizer 이온이 있으면 그쪽) */
function defaultLasers(ionKeys){
  const s=new Set(ionKeys);
  if(s.has("Yb3+")) return [980];
  if(s.has("Nd3+")) return [808];
  if(s.has("Tm3+")) return [1064];
  if(s.has("Er3+")) return [980];
  if(s.has("Ho3+")) return [1150];
  if(s.has("Pr3+")) return [445];
  if(s.has("Eu3+")) return [394];
  if(s.has("Tb3+")) return [378];
  if(s.has("Sm3+")) return [405];
  if(s.has("Dy3+")) return [453];
  if(s.has("Gd3+")) return [273];
  if(s.has("Ce3+")) return [254];
  if(s.has("Eu2+")) return [365];
  if(s.has("Mn4+")||s.has("Cr3+")) return [455];
  if(s.has("Mn2+")) return [405];
  return [];
}

const API={nm2cm,cm2nm,cm2eV,KT_PER_K,CM_PER_EV,IONS,OX,INACTIVE,HOSTS,ANION_PH,ionKey,ionLabel,sup,
  normalize,parseComposition,parseDopants,parseFormula:parseFormulaFull,findHost,canonicalFormula,
  tsEval,levelsFor,annotate,gapInfo,win,gradeOf,GRADES,GRADE_EN,classifyET,classifyPhoton,interIon,forsterEff,erfcx,photonCandidates,paLoops,
  suggestET,wlColor,bandName,curatedFor,dMain,defaultLasers,termHTML,termText,escHTML,TYPE_INFO,SENSITIZERS,P_MPR,P_RAD};
root.ED_CORE=API;
if(typeof module!=="undefined") module.exports=API;
})(typeof window!=="undefined"?window:globalThis);
