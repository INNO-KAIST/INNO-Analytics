/* ============================================================
   INNO Analytics — XPS 탭 (읽기 · 배경 · 피팅 · 불확도 · 정량 · 서베이)
   계산은 XPS_CORE (Web Worker 에서 실행, 취소 가능). 파일 읽기는 XPS_IO. 원소 자료는 XPS_LIB.
   수식 정의: modules/xps/SPEC.md
   ============================================================ */
(function(){
"use strict";
const C=window.XPS_CORE, IO=window.XPS_IO;
const LIB=()=>window.XPS_LIB||{core:{},scofield:{al:{},mg:{}},althermo1:{},auger:[],meta:{sources:{}}};
const esc=s=>String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const $q=(r,s)=>r.querySelector(s), $a=(r,s)=>[...r.querySelectorAll(s)];
/* 알림: XPS 화면이 있으면 그래프 아래 알림 줄(최근 것 + 이전 6 개), 없으면 앱의 toast. act = {label, fn} 이면 버튼 */
let MSGS=[];
const note=(m,k,act)=>{ try{ const el=ROOT&&ROOT.querySelector&&ROOT.querySelector("#xpsMsg");
  if(el){ MSGS.unshift({m:String(m),k:k||"",act:act||null,t:Date.now()}); MSGS=MSGS.slice(0,7); drawMsg(); return; }
  if(typeof toast==="function") toast(m,k); else console.log(m); }catch(e){} };
function drawMsg(){ const el=ROOT&&$q(ROOT,"#xpsMsg"); if(!el) return; if(!MSGS.length){ el.innerHTML=""; el.hidden=true; return; } el.hidden=false;
  const m=MSGS[0], old=MSGS.slice(1);
  el.innerHTML=`<div class="cur ${esc(m.k)}"><span class="ic">${m.k==="warn"?"!":m.k==="ok"?"✓":"i"}</span><span class="tx">${esc(m.m)}</span>${m.act?`<button class="btn xs" id="xpsMsgA">${esc(m.act.label)}</button>`:""}<button class="x" id="xpsMsgX" title="Close">×</button></div>
    ${old.length?`<details class="old"><summary>Earlier messages (${old.length})</summary>${old.map(o=>`<div class="${esc(o.k)}">${esc(o.m)}</div>`).join("")}</details>`:""}`;
  $q(el,"#xpsMsgX").onclick=()=>{ MSGS.shift(); drawMsg(); };
  const b=$q(el,"#xpsMsgA"); if(b) b.onclick=()=>{ const f=m.act&&m.act.fn; m.act=null; drawMsg(); if(f) f(); };
}
const f=(v,d)=>v==null||!isFinite(v)?"—":(+v).toFixed(d==null?2:d);
const g=(v,d)=>{ if(v==null||!isFinite(v)) return "—"; const a=Math.abs(v); if(a!==0&&(a<1e-3||a>=1e6)) return v.toExponential((d||3)-1); return (+v.toPrecision(d||4)).toLocaleString("en-US",{maximumFractionDigits:6}); };
const pm=(v,s,d)=>v==null||!isFinite(v)?"—":`${f(v,d)}${s!=null&&isFinite(s)?` <span class="muted xs">± ${f(s,d)}</span>`:""}`;
/* 결과 표: σ 가 null = 고정(또는 고정된 것에 연결)된 값 → '± —' (0.000 이 아니라) */
const pmF=(v,s,d,fixed)=>v==null||!isFinite(v)?"—":fixed?`${f(v,d)} <span class="muted xs" title="Fixed (or linked to a fixed value) — no statistical uncertainty">± —</span>`:pm(v,s,d);
const SLOTS=8;
const ACCEPT=".vgd,.vms,.npl,.xls,.xlsx,.txt,.csv,.dat,.xy,.avg,.json";
const ORB_J={s:null,p:["3/2","1/2"],d:["5/2","3/2"],f:["7/2","5/2"]};

/* ---------------- 상태 ---------------- */
function defState(){ return {v:1,specs:[],sel:null,proc:{},samples:{},quant:{lib:"scofield",imfp:{model:"power",n:0.6,mat:{rho:5,Nv:24,M:150,Eg:3}},tf:{model:"auto",x:1},custom:{},rows:{}},ui:{tab:"fit",showRef:true,resid:"norm",overlay:[],ovNorm:"max",ovOffset:0,selComp:null,addMode:false,survey:{},ws:"spectra"}}; }
let S=defState(), ROOT=null, PV={}, LOADED=false, WORK=null, RUN=null;
const DBN="innoanalytics_xps", DBS="proj";
function idb(){ return new Promise((res,rej)=>{ try{ const q=indexedDB.open(DBN,1); q.onupgradeneeded=()=>q.result.createObjectStore(DBS); q.onsuccess=()=>res(q.result); q.onerror=()=>rej(q.error); }catch(e){ rej(e); } }); }
async function dbGet(){ try{ const db=await idb(); return await new Promise(res=>{ const t=db.transaction(DBS,"readonly").objectStore(DBS).get("cur"); t.onsuccess=()=>res(t.result||null); t.onerror=()=>res(null); }); }catch(e){ return null; } }
async function dbPut(v){ try{ const db=await idb(); await new Promise(res=>{ const t=db.transaction(DBS,"readwrite"); t.objectStore(DBS).put(v,"cur"); t.oncomplete=res; t.onerror=res; }); }catch(e){} }
let _svT=null;
const _snap=()=>JSON.parse(JSON.stringify(Object.assign({},S),(k,v)=>ArrayBuffer.isView(v)?{__f64:Array.from(v)}:v));
function save(){ clearTimeout(_svT); _svT=setTimeout(()=>{ _svT=null; dbPut(_snap()); },600); }
/* 탭을 닫거나 숨길 때 미뤄 둔 저장을 바로 */
function flushSave(){ if(_svT){ clearTimeout(_svT); _svT=null; dbPut(_snap()); } }
try{ window.addEventListener("pagehide",flushSave); document.addEventListener("visibilitychange",()=>{ if(document.visibilityState==="hidden") flushSave(); }); }catch(e){}
function revive(o){ return JSON.parse(JSON.stringify(o),(k,v)=>v&&typeof v==="object"&&Array.isArray(v.__f64)?Float64Array.from(v.__f64):v); }

/* ---------------- 접근자 ---------------- */
const cur=()=>S.specs.find(s=>s.uid===S.sel)||null;
const byUid=u=>S.specs.find(s=>s.uid===u)||null;
function pr(sp){ sp=sp||cur(); if(!sp) return null; let p=S.proc[sp.uid]; if(!p){ p=S.proc[sp.uid]=defProc(sp); } return p; }
function defProc(sp){ const sv=sp.reg&&sp.reg.kind==="survey", bs=!!(sp.meta&&sp.meta.bgSub); return {lo:sp.be[0],hi:sp.be[sp.n-1],weight:bs?"const":"poisson",bg:{type:sv||bs?"none":"smart",n:5,C:1643,D:0,u3:false,c0:null,k:null,s:0,cons:{s:{mode:"fixed"}}},comps:[],cons:{},result:null,hist:[],models:[],ver:0}; }
const NOSMP="(시료 이름 없음)";   /* 시료 키(상태) — 바꾸지 않음. 화면 표시는 NOSMP_L */
const NOSMP_L="(no sample name)";
const FMTN=t=>t==="텍스트"?"Text":t;   /* 형식 이름의 화면 표시 (저장값은 그대로) */
const SK=sp=>sp&&sp.sample||NOSMP;
/* 시료 키 표시: "자료 / 실험 / 조건 / W1" → 이름 W1, 맥락 "자료 · 실험 · 조건" */
function smpDisp(k){ if(!k||k===NOSMP) return {name:NOSMP_L,ctx:""}; const a=String(k||NOSMP).split(" / "); return {name:a[a.length-1],ctx:a.slice(0,-1).join(" · ")}; }
function smpOf(sp){ const k=SK(sp); if(!S.samples[k]) S.samples[k]={shift:0,target:284.8,refNote:""}; return S.samples[k]; }
const sampleSpecs=k=>S.specs.filter(s=>SK(s)===k);
const isC1s=s=>!!(s&&s.reg&&s.reg.el==="C"&&s.reg.orb==="1s");
const shOf=sp=>sp?smpOf(sp).shift||0:0;
function slotFor(p){ const used=new Set(p.comps.map(c=>c.slot)); for(let i=0;i<SLOTS;i++) if(!used.has(i)) return i; return SLOTS; }
const colorOf=c=>c.slot!=null&&c.slot<SLOTS?`var(--s${c.slot+1})`:"var(--faint)";
function touch(p){ p.ver=(p.ver||0)+1; if(p.result) p.result.stale=true; save(); }

/* ============================================================
   파일 읽기
   ============================================================ */
async function loadFiles(list,lopt){
  const arr=[...list]; if(!arr.length) return;
  let nNew=0, nDup=0, nImg=0; const errs=[], warns=[], axq=[];
  const imgMode=!!(lopt&&lopt.img)||S.ui.ws==="img"; if(imgMode) await imgLoadDB();
  // 이미 열린 같은 스펙트럼(같은 시료·영역·파일·형식·점 수·첫 BE·첫 세기 · VAMAS 블록 번호 · 여러 열 텍스트의 열 번호)은 다시 넣지 않는다
  const sig=s=>[SK(s),s.region,s.file,s.format,s.n,(+s.be[0]).toFixed(4),(+s.y[0]).toPrecision(8),s.meta&&s.meta.vamas?s.meta.vamas.blockIndex:"",s.meta&&s.meta.mcol?s.meta.mcol.col:""].join("|");
  const have=new Set(S.specs.map(sig)); const isDup=s=>{ const k=sig(s); if(have.has(k)) return true; have.add(k); return false; };
  for(const file of arr){
    try{
      const u8=new Uint8Array(await file.arrayBuffer());
      if(/\.json$/i.test(file.name)){ importProject(new TextDecoder().decode(u8)); continue; }
      const mc=multiText(file.name,u8,imgMode);
      if(mc&&mc.kind==="img"){ stackFromText(mc.text,file.name); nImg++; continue; }
      if(mc&&mc.kind==="series"){ mc.spectra.forEach(sp=>{ if(isDup(sp)){ nDup++; return; } S.specs.push(sp); pr(sp); smpOf(sp); nNew++; }); continue; }
      const r=await IO.readFile(file.name,u8,{});
      if(r.vamas&&/^(MAP|SEM)/i.test(r.vamas.mode||"")){ const pos=r.spectra.filter(s=>s.meta&&s.meta.vamas&&s.meta.vamas.x!=null&&s.meta.vamas.y!=null);
        if(pos.length>=4){ await imgLoadDB(); nImg+=stacksFromVamas(pos,file.name).length; r.spectra=r.spectra.filter(s=>!pos.includes(s)); if(!r.spectra.length) continue; } }
      r.spectra.forEach(sp=>{ if(isDup(sp)){ nDup++; return; } S.specs.push(sp); pr(sp); smpOf(sp); nNew++;
        const ac=axisCheck(sp); if(ac){ sp.warn.push(ac.msg); if(ac.fits&&isText(sp)) axq.push(sp); }
        if(sp.warn.length) warns.push(sp.name+": "+sp.warn.join(" / ")); });
      if(r.tables&&r.tables.length) r.tables.forEach(t=>{ if(/peak/i.test(t.name)&&t.rows.length>1) warns.push(file.name+": Avantage '"+t.name+"' 표가 있습니다(Info 탭)"); });
    }catch(e){ errs.push(file.name+": "+(e.message||e)); }
  }
  if(nNew&&(!S.sel||!cur())) S.sel=S.specs[S.specs.length-nNew].uid;
  if(nImg){ S.ui.ws="img"; RES.img=null; RES.roi=null; note(`영상 stack ${nImg}개를 열었습니다 — Imaging`,"ok"); }
  else if(nNew>=2&&S.ui.ws==="spectra"&&S.specs.slice(-nNew).some(s=>s.meta&&(s.meta.mcol||(s.meta.vamas&&s.meta.vamas.expVars&&s.meta.vamas.expVars.length)))) note("여러 단계(깊이·각도) 자료입니다 — 위의 Series · ARXPS 에서 묶어 볼 수 있습니다","",{label:"Open Series",fn:()=>{ S.ui.ws=(S.specs.slice(-1)[0].meta.vamas&&/angle/i.test(((S.specs.slice(-1)[0].meta.vamas.expVars||[])[0]||{}).label||""))?"arxps":"series"; save(); rerender(); }});
  save(); rerender();
  if(nNew) note(`스펙트럼 ${nNew}개를 읽었습니다`+(nDup?` (이미 열린 ${nDup}개는 건너뜀)`:"")+(warns.length?` — 확인할 점 ${warns.length}건(Info 탭)`:""),"ok");
  else if(nDup) note(`이미 열려 있는 스펙트럼입니다 (${nDup}개) — 새로 넣지 않았습니다`);
  errs.forEach(e=>note(e,"warn"));
  axq.forEach(sp=>{ const ac=axisCheck(sp); if(ac) note(`${sp.name}: ${ac.msg}`,"warn",{label:`Read as ${ac.other}`,fn:()=>applyAxis(sp,ac.other,sp.hv)}); });
}
const isText=sp=>sp&&sp.format==="텍스트";
/* 여러 열 텍스트(에너지 + 스펙트럼 여러 열): 열 이름마다 숫자(식각 시간·각도 …)가 있으면 series 스펙트럼, Imaging 화면이면 영상 stack(열 = 화소).
   열 이름에 숫자가 없는 여러 열(Avantage 내보내기: Counts · Backgnd. · Envelope)은 지금처럼 2열 스펙트럼으로 읽는다 */
function multiText(name,u8,imgMode){
  if(/\.(vms|npl|avg|vgd|xls|xlsx|json)$/i.test(name)||(u8[0]===0xD0&&u8[1]===0xCF)||(u8[0]===0x50&&u8[1]===0x4B)) return null;
  const text=new TextDecoder().decode(u8); if(/^\s*VAMAS/i.test(text.slice(0,200))) return null;
  let r; try{ r=IO.textMultiParse(text,name); }catch(e){ return null; }
  if(imgMode&&r.Y.length>=4) return {kind:"img",text};
  if(r.Y.length<3||r.coord.some(v=>v==null)||new Set(r.coord).size!==r.coord.length) return null;
  const base=name.replace(/^.*[\\/]/,"").replace(/\.[^.]+$/,"");
  const reg=[base].concat(r.header).map(h=>IO.parseRegionName(h)).find(q=>q&&q.el&&q.orb);
  const m=/^(.*?)[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*(.*)$/.exec(r.labels[0])||[]; const lab=(m[1]||"").trim()||"Coordinate", unit=(m[2]||"").trim();
  const spectra=r.Y.map((y,j)=>IO.mkSpec({name:`${reg?reg.label:base} · ${r.labels[j]}`,region:reg?reg.label:base,sample:base,file:name,format:"텍스트",hv:r.hv,hvSrc:"기본값(Al Kα 1486.6) — 텍스트에는 없음",be:r.E,y,
    meta:{mcol:{label:lab,unit,coord:r.coord[j],col:j+1,of:r.Y.length},axis:r.axis,header:r.header,commaDecimal:r.commaDecimal},warn:r.warn.slice()}));
  return {kind:"series",spectra}; }
/* 축 점검: 이름에서 읽은 원소·궤도의 라이브러리 BE 가 자료 범위(±8 eV: 대전 이동 여유)에 없으면 경고. 반대 축(KE↔BE)으로 읽으면 들어오는지(fits)도 */
function axisCheck(sp){
  const r=sp&&sp.reg||{}; if(r.kind!=="region"||!r.el||!r.orb) return null;
  const lv=((LIB().core||{})[r.el]||{}).lv||{}; const bes=Object.entries(lv).filter(([k])=>k.startsWith(r.orb)).map(([,v])=>v); if(!bes.length) return null;
  const M=8, lo=sp.be[0], hi=sp.be[sp.n-1], inR=(a,b)=>bes.some(v=>v>=a-M&&v<=b+M);
  if(inR(lo,hi)) return null;
  const ax=(sp.meta&&sp.meta.axis)||"BE", other=ax==="KE"?"BE":"KE", fits=inR(sp.hv-hi,sp.hv-lo);
  const lib=`${r.el} ${r.orb} ${bes.map(v=>f(v,1)).join("/")} eV`;
  return {other,fits,msg:fits?`BE ${f(lo,0)}–${f(hi,0)} eV 에 ${lib}(라이브러리)가 없습니다 — ${other} 축 자료로 보입니다`:`BE ${f(lo,0)}–${f(hi,0)} eV 에 ${lib}(라이브러리)가 없습니다 — 원소·궤도·축을 확인하세요`};
}
/* 텍스트 자료의 축(BE/KE)·hν 를 바꿔 다시 읽기 — BE 가 통째로 바뀌므로 이 스펙트럼의 피팅·보기·survey 판별은 처음부터 */
function applyAxis(sp,axis,hv){
  if(!sp||!byUid(sp.uid)) return;
  let ch=false; try{ ch=IO.reaxis(sp,axis,hv); }catch(e){ note("축을 바꾸지 못했습니다: "+(e.message||e),"warn"); return; }
  if(!ch){ note("바뀐 것이 없습니다"); return; }
  sp.warn=(sp.warn||[]).filter(w=>!/라이브러리\)가 없습니다/.test(w)); const ac=axisCheck(sp); if(ac) sp.warn.push(ac.msg);
  S.proc[sp.uid]=defProc(sp); delete PV[sp.uid]; delete S.ui.survey[sp.uid]; save(); rerender();
  note(`${sp.name}: ${axis} 축, hν ${f(sp.hv,2)} eV 로 다시 읽었습니다 (BE ${f(sp.be[0],1)}–${f(sp.be[sp.n-1],1)} eV)${ac?" — "+ac.msg:""}`,ac?"warn":"ok");
}

/* ============================================================
   보이는 모형 (피팅 전에도 지금 성분으로 계산)
   ============================================================ */
function region(sp,p){ return C.sliceRegion(sp.be,sp.y,p.lo,p.hi); }
/* 배경 설정 + CasaXPS 규약이면 끝점 평균(창 밖 채널 포함, 전체 스펙트럼에서) */
function bgOf(sp,p,R){ if(p.bg.conv!=="casa") return p.bg; R=R||region(sp,p); if(!R) return p.bg; return Object.assign({},p.bg,{_ends:C.casaEnds(sp.be,sp.y,R.i0,R.i1,p.bg.aw!=null?p.bg.aw:1)}); }
function liveBg(p,E,y,sp,R){
  if(p.bg.type!=="active") return C.bgStatic(E,y,sp?bgOf(sp,p,R):p.bg);
  return null;
}
function curModel(sp,p){
  sp=sp||cur(); p=p||pr(sp); if(!sp) return null;
  const R=region(sp,p); if(!R||R.E.length<3) return null;
  const E=Float64Array.from(R.E), y=Float64Array.from(R.y), n=E.length, ctx={emin:E[0],emax:E[n-1]};
  const parts=p.comps.map(c=>{ try{ return C.compEval(c,E,ctx); }catch(e){ return new Float64Array(n); } });
  const P=new Float64Array(n); parts.forEach(v=>{ for(let i=0;i<n;i++) P[i]+=v[i]; });
  let B=liveBg(p,E,y,sp,R);
  if(!B){ // active
    const Q=C.cumtrapz(E,P); const ep=C.endpoints(E,y,p.bg.n||5);
    let c0=p.bg.c0, k=p.bg.k; if(c0==null||k==null){ c0=ep.ylo; k=Q[n-1]>0?Math.max(0,(ep.yhi-ep.ylo)/Q[n-1]):0; if(!p.comps.length){ B=C.bgShirley(E,y,p.bg.n||5,false); } }
    if(!B){ B=new Float64Array(n); for(let i=0;i<n;i++) B[i]=c0+k*Q[i]+(p.bg.s||0)*(E[i]-E[0]); }
  }
  const env=Float64Array.from(B,(b,i)=>b+P[i]);
  let v=null; if(p.result&&!p.result.stale&&p.result.v&&p.result.v.length===n) v=p.result.v;
  if(!v){ if(p.weight==="const") v=Array.from(C.constVar(y)); else { const ys=C.smoothQ(E,y,C.smoothHW(E)); const c=sp.c||1; v=Array.from(ys,t=>Math.max(t,1e-9)/c); } }
  const phi=p.result&&!p.result.stale&&p.result.stats?p.result.stats.chi2nu:null;
  const resid=Float64Array.from(y,(t,i)=>(t-env[i])/Math.sqrt(v[i]*(phi||1)));
  return {E,y,B,parts,env,resid,i0:R.i0,n};
}

/* ============================================================
   그래프
   ============================================================ */
const PW=920, PH=400, RH=96, ML=66, MR=14, MT=12, MB=34;
function plotView(sp){
  let v=PV[sp.uid]; const sh=shOf(sp);
  if(!v){ v=PV[sp.uid]={x0:sp.be[0]+sh,x1:sp.be[sp.n-1]+sh,auto:true}; }
  return v;
}
function niceTicks(a,b,n){ const span=b-a; if(!(span>0)) return [a]; const st0=span/(n||6); const p=Math.pow(10,Math.floor(Math.log10(st0))); const m=st0/p; const st=(m<1.5?1:m<3?2:m<7?5:10)*p; const out=[]; for(let t=Math.ceil(a/st)*st;t<=b+1e-9*st;t+=st) out.push(+t.toFixed(10)); return out; }
function fmtTick(v){ const a=Math.abs(v); if(a>=1e5) return (v/1000).toFixed(0)+"k"; if(a>=1e4) return (v/1000).toFixed(a>=1e5?0:1)+"k"; return String(+v.toPrecision(6)); }
let GEO=null;
/* 좁은 화면(휴대폰): viewBox(920 폭)가 줄어 글자가 4–5 px 가 되므로 글자 배율을 키우고(화면에서 ≈ 8–9 px) 눈금 수·아래 여백을 맞춘다 */
function fontScale(host){ const w=(host&&host.clientWidth)||PW; return w>=700?1:Math.max(1,Math.min(2.2,+(0.75*PW/w).toFixed(2))); }
function drawPlot(){
  const host=ROOT&&$q(ROOT,"#xpsPlot"); if(!host) return;
  const sp=cur();
  if(!sp){ host.innerHTML=`<div class="xps-empty">Drop files here or click <b>Open files</b>.<br><span class="sm">VGD (Avantage raw) · VAMAS (.vms) · Avantage export (.xls/.xlsx) · .avg · 2-column text (.txt/.csv)</span></div>`; return; }
  if(S.ui.overlay.length>1){ drawOverlay(host); return; }
  const p=pr(sp), sh=shOf(sp), v=plotView(sp), M=curModel(sp,p);
  const fs=fontScale(host), mb=Math.round(MB*fs);
  const W=PW, H=PH+RH, pw=W-ML-MR, ph=PH-MT-mb;
  // 보이는 x 범위의 데이터로 y 범위
  let ylo=Infinity, yhi=-Infinity;
  for(let i=0;i<sp.n;i++){ const x=sp.be[i]+sh; if(x<v.x0||x>v.x1) continue; ylo=Math.min(ylo,sp.y[i]); yhi=Math.max(yhi,sp.y[i]); }
  if(M) for(let i=0;i<M.n;i++){ const x=M.E[i]+sh; if(x<v.x0||x>v.x1) continue; yhi=Math.max(yhi,M.env[i]); ylo=Math.min(ylo,M.B[i]); }
  { const pq=procYRange(sp); if(pq){ ylo=Math.min(ylo,pq[0]); yhi=Math.max(yhi,pq[1]); } }
  if(!isFinite(ylo)){ ylo=0; yhi=1; }
  const pad=(yhi-ylo)*0.06||1; let y0=ylo-pad, y1=yhi+pad*1.6;
  if(v.y0!=null&&!v.auto){ y0=v.y0; y1=v.y1; }
  const X=be=>ML+(v.x1-be)/(v.x1-v.x0)*pw, Y=t=>MT+(y1-t)/(y1-y0)*ph;
  GEO={X,Y,v,y0,y1,pw,ph,sh,sp,p,M,inv:px=>v.x1-(px-ML)/pw*(v.x1-v.x0),invY:py=>y1-(py-MT)/ph*(y1-y0)};
  const o=[`<svg class="xps-svg" id="xpsSvg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(sp.region)} spectrum" style="--xfs:${fs}">`];
  o.push(`<defs><clipPath id="xpsClip"><rect x="${ML}" y="${MT}" width="${pw}" height="${ph}"/></clipPath><clipPath id="xpsClipR"><rect x="${ML}" y="${PH}" width="${pw}" height="${RH-22}"/></clipPath></defs>`);
  // 격자·축
  const xt=niceTicks(v.x0,v.x1,fs>1.3?5:8), yt=niceTicks(y0,y1,fs>1.3?4:5);
  xt.forEach(t=>o.push(`<line class="gd" x1="${X(t)}" x2="${X(t)}" y1="${MT}" y2="${MT+ph}"/>`));
  yt.forEach(t=>o.push(`<line class="gd" x1="${ML}" x2="${ML+pw}" y1="${Y(t)}" y2="${Y(t)}"/><text x="${ML-6}" y="${Y(t)+4}" text-anchor="end">${fmtTick(t)}</text>`));
  o.push(`<rect class="ax" x="${ML}" y="${MT}" width="${pw}" height="${ph}" fill="none"/>`);
  if(fs<=1.3) o.push(`<text class="al" x="${14}" y="${MT+ph/2}" transform="rotate(-90 14 ${MT+ph/2})" text-anchor="middle">${sp.c?"CPS":"Intensity (CPS)"}</text>`);   // 좁은 화면: 큰 눈금 글자와 겹쳐 뺌
  o.push(`<g clip-path="url(#xpsClip)">`);
  // 영역 밖 가림
  const rl=p.lo+sh, rh=p.hi+sh;
  if(rh<v.x1) o.push(`<rect class="veil" x="${ML}" y="${MT}" width="${Math.max(0,X(rh)-ML)}" height="${ph}"/>`);
  if(rl>v.x0) o.push(`<rect class="veil" x="${X(rl)}" y="${MT}" width="${Math.max(0,ML+pw-X(rl))}" height="${ph}"/>`);
  // 참고선(라이브러리)
  if(S.ui.showRef) refLines(sp).forEach(r=>{ const x=X(r.be+0); if(x<ML||x>ML+pw) return; o.push(`<line class="ref" x1="${x}" x2="${x}" y1="${MT}" y2="${MT+10}"/><text class="ref" x="${x+2}" y="${MT+20}">${esc(r.lab)}</text>`); });
  // 성분
  if(M){
    p.comps.forEach((c,ci)=>{ const a=M.parts[ci]; if(!a) return; let d="", d2=""; for(let i=0;i<M.n;i++){ const x=X(M.E[i]+sh), yy=Y(M.B[i]+a[i]); d+=(i?"L":"M")+x.toFixed(1)+","+yy.toFixed(1); } for(let i=M.n-1;i>=0;i--){ d2+="L"+X(M.E[i]+sh).toFixed(1)+","+Y(M.B[i]).toFixed(1); }
      const col=colorOf(c); o.push(`<path d="${d}${d2}Z" fill="${col}" fill-opacity="${S.ui.selComp===c.id?0.32:0.18}" stroke="none"/><path d="${d}" fill="none" stroke="${col}" stroke-width="${S.ui.selComp===c.id?2.2:1.5}"/>`); });
    // 배경
    let db=""; for(let i=0;i<M.n;i++) db+=(i?"L":"M")+X(M.E[i]+sh).toFixed(1)+","+Y(M.B[i]).toFixed(1); o.push(`<path class="bgl" d="${db}"/>`);
  }
  // 데이터
  let nIn=0; for(let i=0;i<sp.n;i++){ const x=sp.be[i]+sh; if(x>=v.x0&&x<=v.x1) nIn++; }
  if(nIn<=900){ let s=""; for(let i=0;i<sp.n;i++){ const x=sp.be[i]+sh; if(x<v.x0-1||x>v.x1+1) continue; s+=`<circle class="datp" cx="${X(x).toFixed(1)}" cy="${Y(sp.y[i]).toFixed(1)}" r="1.7"/>`; } o.push(s); }
  else { let d=""; let st=false; for(let i=0;i<sp.n;i++){ const x=sp.be[i]+sh; if(x<v.x0-1||x>v.x1+1) continue; d+=(st?"L":"M")+X(x).toFixed(1)+","+Y(sp.y[i]).toFixed(1); st=true; } o.push(`<path class="dat" d="${d}"/>`); }
  // 포락선
  if(M&&p.comps.length){ let d=""; for(let i=0;i<M.n;i++) d+=(i?"L":"M")+X(M.E[i]+sh).toFixed(1)+","+Y(M.env[i]).toFixed(1); o.push(`<path class="env" d="${d}"/>`); }
  // Process 탭: 디컨볼루션 · 다시 표본화 · 평활/미분 · Fermi edge 겹쳐 그리기
  o.push(procOverlay(sp,X,Y,y0,y1));
  // 성분 표지(꼭대기) — 끌어서 위치·세기
  if(M) p.comps.forEach((c,ci)=>{ const a=M.parts[ci]; let im=0; for(let i=1;i<M.n;i++) if(a[i]>a[im]) im=i; const be=c.be; const x=X(be+sh); const bgAt=interp(M.E,M.B,be); const hgt=compHeight(c); const yy=Y(bgAt+hgt); c._mk={im};
    o.push(`<circle class="mk ${S.ui.selComp===c.id?"sel":""}" data-mk="${esc(c.id)}" cx="${x.toFixed(1)}" cy="${yy.toFixed(1)}" r="6" fill="${colorOf(c)}"/>`);
    if(p.comps.length<=6) o.push(`<text class="cl" x="${x.toFixed(1)}" y="${(yy-10).toFixed(1)}" text-anchor="middle">${esc(c.name||("#"+(ci+1)))}</text>`); });
  // spline 매듭 (끌기 · 두 번 클릭으로 지우기) — 데이터 점 위에 그려야 잡힌다
  if(M&&p.bg.type==="spline") o.push(knotMarks(sp,p,M,X,Y,sh));
  // 서베이 판별 표지
  const sv=S.ui.survey[sp.uid]; if(sv&&sv.marks) sv.marks.forEach(m=>{ const x=X(m.be+sh); if(x<ML||x>ML+pw) return; o.push(`<line class="ref" x1="${x}" x2="${x}" y1="${Y(m.y)-6}" y2="${Y(m.y)-20}"/><text class="pid" x="${x}" y="${Y(m.y)-24}" text-anchor="middle">${esc(m.lab)}</text>`); });
  o.push(`</g>`);
  // 영역 손잡이
  [["lo",rl],["hi",rh]].forEach(([k,xv])=>{ if(xv<v.x0||xv>v.x1) return; const x=X(xv); o.push(`<line class="hdl" x1="${x}" x2="${x}" y1="${MT}" y2="${MT+ph}"/><line class="hdlhit" data-hdl="${k}" x1="${x}" x2="${x}" y1="${MT}" y2="${MT+ph}"/>`); });
  // x 축 글자
  xt.forEach(t=>o.push(`<text x="${X(t)}" y="${MT+ph+4+Math.round(11*fs)}" text-anchor="middle">${fmtTick(t)}</text>`));
  o.push(`<text class="al" x="${ML+pw/2}" y="${MT+ph+mb-4}" text-anchor="middle">Binding energy (eV)${sh?` — charge-corrected (${sh>0?"+":""}${f(sh,2)} eV)`:""}</text>`);
  // 잔차
  if(M&&p.comps.length){
    const r=M.resid; let rmax=0; for(let i=0;i<M.n;i++){ const x=M.E[i]+sh; if(x>=v.x0&&x<=v.x1) rmax=Math.max(rmax,Math.abs(r[i])); } rmax=Math.max(3,Math.min(rmax*1.1,50));
    const r0=PH+4, rh2=RH-26, RY=t=>r0+rh2/2-t/rmax*rh2/2;
    o.push(`<rect class="ax" x="${ML}" y="${r0}" width="${pw}" height="${rh2}" fill="none"/><line class="rz" x1="${ML}" x2="${ML+pw}" y1="${RY(2)}" y2="${RY(2)}"/><line class="rz" x1="${ML}" x2="${ML+pw}" y1="${RY(-2)}" y2="${RY(-2)}"/><line class="gd" x1="${ML}" x2="${ML+pw}" y1="${RY(0)}" y2="${RY(0)}"/>`);
    o.push(`<text x="${ML-6}" y="${RY(2)+4}" text-anchor="end">+2</text><text x="${ML-6}" y="${RY(-2)+4}" text-anchor="end">−2</text>${fs<=1.3?`<text class="al" x="${14}" y="${r0+rh2/2}" transform="rotate(-90 14 ${r0+rh2/2})" text-anchor="middle" style="font-size:10.5px">Residual/σ</text>`:""}`);
    let d=""; for(let i=0;i<M.n;i++){ const x=X(M.E[i]+sh); if(x<ML||x>ML+pw) continue; d+=(d?"L":"M")+x.toFixed(1)+","+RY(Math.max(-rmax,Math.min(rmax,r[i]))).toFixed(1); }
    o.push(`<path d="${d}" fill="none" stroke="var(--dat)" stroke-width="1" clip-path="url(#xpsClipR)"/>`);
  }
  o.push(`<g id="xpsOver"></g></svg>`);
  host.innerHTML=`<div class="xps-plotwrap">${o.join("")}<div class="xps-tip" id="xpsTip" hidden></div></div>`+legendHTML(sp,p)+procLegend(sp);
  bindPlot(host);
}
function legendHTML(sp,p){
  const items=[`<span><i style="background:var(--dat);border-radius:50%"></i>Data</span>`];
  if(p.comps.length) items.push(`<span><span class="ln"></span>Envelope</span>`);
  items.push(`<span><span class="ln bg"></span>Background · ${esc(C.BG_TYPES[p.bg.type]||p.bg.type)}</span>`);
  p.comps.forEach((c,i)=>items.push(`<span><i style="background:${colorOf(c)}"></i>${esc(c.name||("#"+(i+1)))}</span>`));
  return `<div class="xps-legend">${items.join("")}</div>`;
}
function interp(E,Y,x){ const n=E.length; if(x<=E[0]) return Y[0]; if(x>=E[n-1]) return Y[n-1]; let lo=0,hi=n-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(E[m]<=x) lo=m; else hi=m; } const t=(x-E[lo])/(E[hi]-E[lo]); return Y[lo]+t*(Y[hi]-Y[lo]); }
function compHeight(c){ try{ const t=C.trueFWHM(c.shape,C.shapeParams(c,1)); const A=c.doublet?c.area/(1+c.doublet.r):c.area; return A*t.ymax; }catch(e){ return 0; } }
function refLines(sp){
  const L=LIB(), r=sp.reg||{}; if(!r.el||!L.core[r.el]) return [];
  const lv=L.core[r.el].lv||{}; const out=[];
  Object.entries(lv).forEach(([k,be])=>{ if(k.startsWith(r.orb)) out.push({be,lab:r.el+" "+k}); });
  return out;
}
function drawOverlay(host){
  const list=S.ui.overlay.map(byUid).filter(Boolean); if(!list.length){ host.innerHTML=""; return; }
  let x0=Infinity,x1=-Infinity; list.forEach(sp=>{ const sh=shOf(sp); x0=Math.min(x0,sp.be[0]+sh); x1=Math.max(x1,sp.be[sp.n-1]+sh); });
  const fs=fontScale(host), mb=Math.round(MB*fs);
  const W=PW, H=PH, pw=W-ML-MR, ph=PH-MT-mb; const norm=S.ui.ovNorm, off=+S.ui.ovOffset||0;
  const series=list.map((sp,k)=>{ const sh=shOf(sp); let a=Array.from(sp.y); let mn=Math.min(...a), mx=Math.max(...a);
    if(norm==="max") a=a.map(v=>(v-mn)/(mx-mn||1)); else if(norm==="area"){ let s=0; for(let i=1;i<sp.n;i++) s+=0.5*(a[i]+a[i-1]-2*mn)*(sp.be[i]-sp.be[i-1]); a=a.map(v=>(v-mn)/(s||1)); }
    return {sp,x:Array.from(sp.be,b=>b+sh),y:a.map(v=>v+k*off*(norm==="none"?mx:1))}; });
  let y0=Infinity,y1=-Infinity; series.forEach(s=>s.y.forEach(v=>{ y0=Math.min(y0,v); y1=Math.max(y1,v); })); const pad=(y1-y0)*0.05||1; y0-=pad; y1+=pad;
  const X=be=>ML+(x1-be)/(x1-x0)*pw, Y=t=>MT+(y1-t)/(y1-y0)*ph;
  const o=[`<svg class="xps-svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="Overlay" style="--xfs:${fs}">`];
  niceTicks(x0,x1,fs>1.3?5:8).forEach(t=>o.push(`<line class="gd" x1="${X(t)}" x2="${X(t)}" y1="${MT}" y2="${MT+ph}"/><text x="${X(t)}" y="${MT+ph+4+Math.round(11*fs)}" text-anchor="middle">${fmtTick(t)}</text>`));
  o.push(`<rect class="ax" x="${ML}" y="${MT}" width="${pw}" height="${ph}" fill="none"/>`);
  series.forEach((s,k)=>{ let d=""; s.x.forEach((x,i)=>d+=(i?"L":"M")+X(x).toFixed(1)+","+Y(s.y[i]).toFixed(1)); o.push(`<path d="${d}" fill="none" stroke="var(--s${(k%SLOTS)+1})" stroke-width="1.6"/>`); });
  o.push(`<text class="al" x="${ML+pw/2}" y="${MT+ph+mb-4}" text-anchor="middle">Binding energy (eV, per-sample charge correction applied)</text></svg>`);
  host.innerHTML=o.join("")+`<div class="xps-legend">${series.map((s,k)=>`<span><i style="background:var(--s${(k%SLOTS)+1})"></i><span title="${esc(SK(s.sp))}">${esc(smpDisp(SK(s.sp)).name)}</span> · ${esc(s.sp.region)}</span>`).join("")}</div>`;
}

/* ---------------- 그래프 조작 ---------------- */
let DRAG=null, KNOT2=null;
function svgPt(svg,e){ const pt=svg.createSVGPoint(); pt.x=e.clientX; pt.y=e.clientY; const m=svg.getScreenCTM(); return m?pt.matrixTransform(m.inverse()):{x:0,y:0}; }
function bindPlot(host){
  const svg=$q(host,"#xpsSvg"); if(!svg) return; const tip=$q(host,"#xpsTip");
  svg.addEventListener("wheel",e=>{ if(!GEO) return; e.preventDefault(); const q=svgPt(svg,e); const v=GEO.v; const xc=GEO.inv(q.x); const k=e.deltaY>0?1.2:1/1.2; v.x0=xc-(xc-v.x0)*k; v.x1=xc+(v.x1-xc)*k; clampView(GEO.sp,v); drawPlot(); },{passive:false});
  svg.addEventListener("dblclick",e=>{ const sp=cur(); if(!sp) return; const t=e.target; if(t&&t.dataset&&t.dataset.knot!=null){ delKnot(+t.dataset.knot); return; } delete PV[sp.uid]; drawPlot(); });
  svg.addEventListener("pointerdown",e=>{
    if(!GEO) return; const q=svgPt(svg,e); const t=e.target;
    if(t.dataset&&t.dataset.hdl){ DRAG={kind:"hdl",k:t.dataset.hdl,lo:GEO.p.lo,hi:GEO.p.hi}; svg.setPointerCapture(e.pointerId); return; }
    if(t.dataset&&t.dataset.knot!=null){ const i=+t.dataset.knot, k=(GEO.p.bg.knots||[])[i]; if(!k) return;
      // 두 번 누르기 = 지우기 (pointer capture 때문에 dblclick 의 대상이 svg 가 되므로 직접 잰다)
      const now=Date.now(); if(KNOT2&&KNOT2.i===i&&now-KNOT2.t<450){ KNOT2=null; delKnot(i); return; } KNOT2={i,t:now};
      pushHist(GEO.p,"매듭 끌기"); DRAG={kind:"knot",k,e0:k.e,y0:k.y}; svg.setPointerCapture(e.pointerId); return; }
    if(t.dataset&&t.dataset.mk){ const c=GEO.p.comps.find(c=>c.id===t.dataset.mk); if(c){ S.ui.selComp=c.id; DRAG={kind:"mk",c,x:q.x,y:q.y,be:c.be,area:c.area,h:compHeight(c)}; pushHist(GEO.p,"끌기"); svg.setPointerCapture(e.pointerId); renderPanel(); } return; }
    if(q.x<ML||q.x>ML+GEO.pw||q.y<MT||q.y>MT+GEO.ph) return;
    if(S.ui.knotMode){ addKnot(GEO.inv(q.x)-GEO.sh); return; }
    if(S.ui.addMode){ addCompAt(GEO.inv(q.x)-GEO.sh,GEO.invY(q.y)); return; }
    DRAG={kind:"zoom",x:q.x}; svg.setPointerCapture(e.pointerId);
  });
  svg.addEventListener("pointermove",e=>{
    if(!GEO) return; const q=svgPt(svg,e);
    if(!DRAG){ hover(svg,tip,q,e); return; }
    if(DRAG.kind==="hdl"){ const be=GEO.inv(q.x)-GEO.sh; const p=GEO.p; if(DRAG.k==="lo") p.lo=Math.min(be,p.hi-0.5); else p.hi=Math.max(be,p.lo+0.5); clampRegion(GEO.sp,p); touch(p); drawPlot(); syncRegionInputs(); return; }
    if(DRAG.kind==="knot"){ DRAG.k.e=GEO.inv(q.x)-GEO.sh; DRAG.k.y=GEO.invY(q.y); touch(GEO.p); drawPlot(); return; }
    if(DRAG.kind==="mk"){ const c=DRAG.c; const dx=GEO.inv(q.x)-GEO.inv(DRAG.x); c.be=DRAG.be+dx; const M=GEO.M; const bgAt=M?interp(M.E,M.B,c.be):0; const hNew=Math.max(1e-9,GEO.invY(q.y)-bgAt); if(DRAG.h>0) c.area=Math.max(1e-9,DRAG.area*hNew/DRAG.h); touch(GEO.p); drawPlot(); return; }
    if(DRAG.kind==="zoom"){ let ov=$q(svg,"#xpsOver"); const a=Math.min(DRAG.x,q.x), b=Math.max(DRAG.x,q.x); ov.innerHTML=`<rect class="zr" x="${a}" y="${MT}" width="${b-a}" height="${GEO.ph}"/>`; }
  });
  const up=e=>{
    if(!DRAG) return; const q=svgPt(svg,e); const d=DRAG; DRAG=null;
    if(d.kind==="zoom"){ if(Math.abs(q.x-d.x)>6){ const a=GEO.inv(Math.max(d.x,q.x)), b=GEO.inv(Math.min(d.x,q.x)); const v=GEO.v; v.x0=a; v.x1=b; clampView(GEO.sp,v); } else { selectNearComp(q); } drawPlot(); }
    if(d.kind==="knot"&&d.k.e===d.e0&&d.k.y===d.y0){ const h=GEO.p.hist; if(h&&h.length) h.pop(); return; }   // 누르기만 = 바뀐 것 없음
    if(d.kind==="mk"||d.kind==="hdl"||d.kind==="knot"){ renderPanel(); drawList(); }
  };
  svg.addEventListener("pointerup",up); svg.addEventListener("pointercancel",up);
  svg.addEventListener("pointerleave",()=>{ if(!DRAG&&tip) tip.hidden=true; const ov=$q(svg,"#xpsOver"); if(ov&&!DRAG) ov.innerHTML=""; });
  bindTouch(host);
}
/* ---------------- 터치 (window.INNO_TOUCH) ----------------
   두 손가락 = BE 축 확대(핀치 중심 고정)·이동, 한 손가락 끌기 = BE 축 이동, 두 번 탭 = 전체 보기(데스크톱의 두 번 클릭).
   탭 = 그 자리 값 풍선(마우스를 올렸을 때와 같음)·가까운 성분 선택, 'Add component' 를 누른 뒤의 탭 = 그 자리에 성분 놓기.
   길게 누르기 = 풍선이 손가락을 따라감. 성분 꼭대기 표지·영역 손잡이는 한 손가락으로 그대로 끈다(둘째 손가락이 오면 끌기 취소). */
function bindTouch(host){
  if(!window.INNO_TOUCH||host._touch) return;
  const svgOf=()=>$q(host,"#xpsSvg");
  const tipAt=p=>{ const svg=svgOf(); if(!svg||!GEO) return; hover(svg,$q(host,"#xpsTip"),svgPt(svg,p)); };
  host._touch=INNO_TOUCH.attach(host,{
    filter:e=>!!(e.target&&e.target.closest&&e.target.closest("#xpsSvg")),
    touchAction:false,                                   // touch-action:none 은 CSS 가 그래프(svg)에만 — 범례 위에서는 쪽이 스크롤됨
    oneFinger:e=>{ const t=e.target; return t&&t.dataset&&(t.dataset.hdl||t.dataset.mk||t.dataset.knot!=null)?"pass":"pan"; },
    getView:()=>GEO?{uid:GEO.sp.uid,x0:GEO.v.x0,x1:GEO.v.x1}:null,
    setView:(v0,g)=>{ const svg=svgOf(); if(!v0||!svg||!GEO||GEO.sp.uid!==v0.uid) return;
      const q0=svgPt(svg,{clientX:g.clientX0,clientY:g.clientY0}), q=svgPt(svg,{clientX:g.clientX,clientY:g.clientY});
      const span0=v0.x1-v0.x0, be0=v0.x1-(q0.x-ML)/GEO.pw*span0, span=span0/(g.kind==="pinch"?g.scaleX:1);
      const v=GEO.v; v.x1=be0+(q.x-ML)/GEO.pw*span; v.x0=v.x1-span; clampView(GEO.sp,v); drawPlot(); },
    onPassCancel:cancelDrag,
    onTap:p=>{ const svg=svgOf(); if(!svg||!GEO) return; const q=svgPt(svg,p);
      if(q.x<ML||q.x>ML+GEO.pw||q.y<MT||q.y>MT+GEO.ph){ tipAt(p); return; }
      if(S.ui.knotMode){ addKnot(GEO.inv(q.x)-GEO.sh); return; }
      if(S.ui.addMode){ addCompAt(GEO.inv(q.x)-GEO.sh,GEO.invY(q.y)); return; }
      selectNearComp(q); drawPlot(); tipAt(p); },
    onDoubleTap:()=>{ const sp=cur(); if(!sp) return; delete PV[sp.uid]; drawPlot(); },
    doubleTap:()=>!S.ui.addMode&&!S.ui.knotMode,                         // 성분 놓기는 탭 즉시
    onLongPress:tipAt, onLongPressMove:tipAt
  });
}
/* 손가락 끌기를 둘째 손가락이 가져갈 때: 끌기 전 값으로 */
function cancelDrag(){ const d=DRAG; DRAG=null; if(!d||!GEO) return;
  if(d.kind==="mk"){ d.c.be=d.be; d.c.area=d.area; const h=GEO.p.hist; if(h&&h.length) h.pop(); touch(GEO.p); }
  else if(d.kind==="hdl"&&d.lo!=null){ GEO.p.lo=d.lo; GEO.p.hi=d.hi; touch(GEO.p); syncRegionInputs(); }
  else if(d.kind==="knot"){ d.k.e=d.e0; if(d.y0==null) delete d.k.y; else d.k.y=d.y0; const h=GEO.p.hist; if(h&&h.length) h.pop(); touch(GEO.p); }
  drawPlot(); }
function clampView(sp,v){ const sh=shOf(sp); const a=sp.be[0]+sh, b=sp.be[sp.n-1]+sh; if(v.x1-v.x0<0.2){ const c=(v.x0+v.x1)/2; v.x0=c-0.1; v.x1=c+0.1; } if(v.x1-v.x0>b-a){ v.x0=a; v.x1=b; } if(v.x0<a){ v.x1+=a-v.x0; v.x0=a; } if(v.x1>b){ v.x0-=v.x1-b; v.x1=b; } }
function clampRegion(sp,p){ p.lo=Math.max(sp.be[0],Math.min(p.lo,sp.be[sp.n-1])); p.hi=Math.max(sp.be[0],Math.min(p.hi,sp.be[sp.n-1])); }
function selectNearComp(q){ const M=GEO.M; if(!M) return; const be=GEO.inv(q.x)-GEO.sh; const yv=GEO.invY(q.y); let best=null; GEO.p.comps.forEach((c,ci)=>{ const a=M.parts[ci]; const B=interp(M.E,M.B,be), val=interp(M.E,a,be); if(yv>=B&&yv<=B+val&&(!best||val>best.v)) best={c,v:val}; }); S.ui.selComp=best?best.c.id:null; renderPanel(); }
function hover(svg,tip,q){
  const ov=$q(svg,"#xpsOver"); if(q.x<ML||q.x>ML+GEO.pw||q.y<MT||q.y>PH+RH){ tip.hidden=true; ov.innerHTML=""; return; }
  const be=GEO.inv(q.x), sp=GEO.sp, sh=GEO.sh;
  // 가장 가까운 데이터 점
  let k=0, d=Infinity; for(let i=0;i<sp.n;i++){ const dd=Math.abs(sp.be[i]+sh-be); if(dd<d){ d=dd; k=i; } }
  const x=GEO.X(sp.be[k]+sh); ov.innerHTML=`<line class="xh" x1="${x}" x2="${x}" y1="${MT}" y2="${MT+GEO.ph}"/>`;
  const rows=[`<div><span class="k">BE</span> <b>${f(sp.be[k]+sh,2)}</b> eV${sh?` <span class="k">(measured ${f(sp.be[k],2)})</span>`:""}</div>`,`<div><span class="k">Data</span> ${g(sp.y[k],5)}</div>`];
  const M=GEO.M; if(M){ const j=k-M.i0; if(j>=0&&j<M.n){ rows.push(`<div><span class="k">Background</span> ${g(M.B[j],5)}</div>`); if(GEO.p.comps.length){ rows.push(`<div><span class="k">Envelope</span> ${g(M.env[j],5)} <span class="k">Residual/σ</span> ${f(M.resid[j],2)}</div>`); GEO.p.comps.forEach((c,ci)=>rows.push(`<div><i style="background:${colorOf(c)}"></i>${esc(c.name||"#"+(ci+1))} ${g(M.parts[ci][j],4)}</div>`)); } } }
  tip.innerHTML=rows.join(""); tip.hidden=false;
  const wrap=svg.getBoundingClientRect(); const rx=(x/PW)*wrap.width; tip.style.left=(rx+14+180>wrap.width?rx-14-tip.offsetWidth:rx+14)+"px"; tip.style.top="8px";
}

/* ============================================================
   성분 추가 · 자동 제안 · 라이브러리
   ============================================================ */
function pushHist(p,why){ p.hist=p.hist||[]; p.hist.push({why,comps:C.clone(p.comps),cons:C.clone(p.cons),bg:C.clone(p.bg),lo:p.lo,hi:p.hi}); if(p.hist.length>25) p.hist.shift(); }
function undo(){ const p=pr(); if(!p||!p.hist||!p.hist.length){ note("되돌릴 것이 없습니다","warn"); return; } const h=p.hist.pop(); p.comps=h.comps; p.cons=h.cons; p.bg=h.bg; p.lo=h.lo; p.hi=h.hi; touch(p); rerender(); note("되돌렸습니다: "+h.why); }
function defName(sp,p,doublet){ const r=sp.reg||{}; const mu=r.el&&C.multipletInfo(r.el,r.orb); const base=r.el?(r.el+" "+(r.orb||"")+(doublet||mu||!ORB_J[(r.orb||"s").slice(-1)]?"":ORB_J[(r.orb||"s").slice(-1)][0])):"Component"; const used=new Set(p.comps.map(c=>c.name||"")); if(!used.has(base)) return base; for(let i=1;i<26;i++){ const nm=base+" "+String.fromCharCode(65+i); if(!used.has(nm)) return nm; } return base+" "+(p.comps.length+1); }
function addCompAt(be,yv){
  const sp=cur(), p=pr(sp), M=curModel(sp,p); pushHist(p,"성분 추가");
  const bgAt=M?interp(M.E,M.B,be):0; const h=Math.max(1,(yv!=null?yv:(M?interp(M.E,M.y,be):1))-bgAt);
  const r=sp.reg||{}; const dbl=S.ui.addDoublet&&soInfo(r.el,r.orb);
  const c=C.newComp(S.ui.addShape||"GLS",be,{F:1.2,el:r.el||"",line:dbl?r.orb:(r.orb||""),ch:chOf(sp)});
  if(dbl&&dbl.multiplet) note("다중항 준위에 doublet 을 넣었습니다 — 껍질이 닫힌 상태("+(dbl.multiplet.closed||"해당 없음")+")가 아니면 doublet 모형은 맞지 않습니다","warn");
  if(dbl){ c.doublet={dE:dbl.dE||1,r:dbl.r,q:1}; p.cons[c.id]={dE:C.dEConstraint(dbl),r:{mode:"fixed"},q:{mode:"fixed"}}; if(!dbl.dE) note("Δ(스핀–궤도 간격)를 라이브러리에서 쓸 수 없어 1 eV 로 두고 풀어 둡니다 — 문헌값을 넣으세요","warn"); }
  c.name=defName(sp,p,!!dbl); c.slot=slotFor(p);
  const hh=compHeight(Object.assign({},c,{area:1})); c.area=h/(hh||1)*(dbl?(1+dbl.r):1);
  p.comps.push(c); S.ui.selComp=c.id; S.ui.addMode=false; touch(p); rerender();
}
/* 다중항 준위 안내 (피팅·결과 탭) */
function multBanner(sp){ const r=sp.reg||{}; const mu=r.el&&C.multipletInfo(r.el,r.orb); if(!mu) return "";
  return `<div class="note" style="margin:4px 0 6px"><b>다중항 준위</b> — ${esc(mu.why)}. Auto suggest 는 넓은 성분들로 모양(envelope)만 맞춥니다: <b>성분 하나하나를 화학 상태로 해석하지 말고</b>, 넓이 합을 정량에 쓰세요. 화학 상태는 같은 조건의 표준 시료 스펙트럼 모양과 비교해 판단합니다.${mu.closed?` 예외: ${esc(mu.closed)} 처럼 껍질이 닫힌 상태는 일반 doublet 으로 맞출 수 있습니다.`:""}</div>`; }
function soInfo(el,orb){ const d=C.doubletInfo(LIB(),el,orb); if(!d) return null; return d; }
function soNote(so){ if(!so) return ""; const a=so.alt||{}; const both=a.avantage!=null&&a.booklet!=null;
  return `Δ ${so.dE!=null?so.dE+" eV":"—"} (${so.src}${both?`; Booklet ${a.booklet} eV`:""}) · 넓이 비 ${f(so.r,3)} — ${so.dE==null?"직접 넣으세요":so.resolvable?"분해되는 doublet 이라 Δ 를 데이터로 정함(±범위 안)":"분해되지 않는 doublet 이라 Δ 고정"+(both&&so.spread>0.1?` — 두 출처가 ${f(so.spread,2)} eV 다름: 문헌값을 확인하세요`:"")}`; }
/* spArg: 대상 스펙트럼(없으면 지금 보는 것) — 계산 중 다른 스펙트럼으로 옮겨도 그 스펙트럼에 적용 */
async function autoSuggest(spArg){
  const sp=spArg&&spArg.uid?spArg:cur(); if(!sp) return; const p=pr(sp); const R=region(sp,p); if(!R) return;
  if(RUN){ note("이미 계산 중입니다","warn"); return; }
  const r=sp.reg||{};
  if(!r.el) note("이 스펙트럼의 원소·궤도를 모릅니다(파일 이름에서 못 읽음) — Info 탭에서 지정하면 doublet·감도계수를 씁니다","warn");
  const s=C.suggestModel(Float64Array.from(R.E),Float64Array.from(R.y),{bg:p.bg,el:r.el,orb:r.orb,lib:LIB(),shape:S.ui.addShape||"GLS"});
  if(!s.comps.length){ note("뚜렷한 봉우리를 찾지 못했습니다 — 그래프에서 직접 추가하세요","warn"); return; }
  pushHist(p,"자동 제안");
  p.comps=[]; p.cons=s.cons;
  s.comps.forEach(c=>{ c.name=defName(sp,p,!!c.doublet); c.slot=slotFor(p); p.comps.push(c); });
  touch(p); rerender();
  const soTxt=s.mult?" — 다중항 준위라 doublet 대신 넓은 성분으로 모양만 맞췄습니다(성분 하나하나는 화학 상태가 아님, 정량은 넓이 합)":s.info&&s.info.dE&&p.comps.some(c=>c.doublet)?" — "+soNote(s.info):"";
  // 2 단계: 잔차 보강(덜 설명된 봉우리 추가 · 넓이 0 으로 간 성분 덜기) — 판단은 피팅 안 Shirley + 기울기 배경, z≥4 · ΔBIC≥10
  let prob; try{ prob=probOf(sp,p); }catch(e){ note(`성분 ${p.comps.length}개를 제안했습니다${soTxt}. 출발점일 뿐이니 확인하세요`); return; }
  setProg("자동 제안: 잔차 보강 중…",0.3); const ver=p.ver;
  try{
    const [o]=await pool([{cmd:"auto",prob,opt:{}}]); hideProg(); if(!o) throw new Error("결과 없음");
    if(p.ver!==ver||S.proc[sp.uid]!==p){ note(`자동 제안: 계산 중에 ${sp.region} 성분이 바뀌어 잔차 보강 결과를 버렸습니다`,"warn"); return; }
    const ch=o.added||o.removed||(o.relaxed&&o.relaxed.length);
    if(ch){
      const keep=new Map(p.comps.map(c=>[c.id,c]));
      const olds=o.comps.filter(c=>keep.has(c.id)).map(c=>{ const k=keep.get(c.id); c.name=k.name; c.slot=k.slot; return c; });
      p.comps=olds; o.comps.filter(c=>!keep.has(c.id)).forEach(c=>{ c.name=defName(sp,p,!!c.doublet); c.slot=slotFor(p); p.comps.push(c); });
      p.cons=o.cons; }
    else { const m=new Map(o.comps.map(c=>[c.id,c])); p.comps.forEach(c=>{ const q=m.get(c.id); if(q) Object.keys(q).forEach(k=>{ if(k!=="name"&&k!=="slot"&&k!=="color") c[k]=q[k]; }); }); }
    touch(p); rerender();
    note(`성분 ${p.comps.length}개를 제안했습니다${soTxt}${ch?` · 잔차 보강: ${o.log.filter(t=>/채택|뺌/.test(t)).map(t=>t.replace(/^→ 채택: /,"")).join("; ")}`:""}. 출발점일 뿐이니 확인한 뒤 'Fit' 을 누르세요`);
  }catch(e){ hideProg(); RUN=null; note(`성분 ${p.comps.length}개를 제안했습니다${soTxt} (잔차 보강 생략: ${e.message||e}). 출발점일 뿐이니 확인하세요`,"warn"); }
}
/* C 1s adventitious 템플릿 (출발값; 간격은 문헌의 흔한 값 ±0.4 eV 범위로만 풀어 둠) */
function templateC1s(){
  const sp=cur(), p=pr(sp); pushHist(p,"C 1s 템플릿");
  const R=region(sp,p); const E=Float64Array.from(R.E), y=Float64Array.from(R.y); const B=C.bgShirley(E,y,5,true);
  let im=0; for(let i=1;i<E.length;i++) if(y[i]-B[i]>y[im]-B[im]) im=i; const be0=E[im], h0=y[im]-B[im];
  const defs=[["C–C/C–H",0,1],["C–O",1.5,0.2],["C=O",3.0,0.08],["O–C=O",4.0,0.1]];
  p.comps=[]; p.cons={};
  defs.forEach(([nm,d,rel],i)=>{ const c=C.newComp("GLS",be0+d,{F:1.2,el:"C",line:"1s",name:nm}); c.slot=slotFor(p); const hh=compHeight(Object.assign({},c,{area:1})); c.area=h0*rel/(hh||1); p.comps.push(c);
    if(i){ const f0=p.comps[0]; p.cons[c.id]={be:{mode:"free",lo:be0+d-0.4,hi:be0+d+0.4},F:{mode:"link",ref:f0.id,k:1,b:0},m:{mode:"link",ref:f0.id,k:1,b:0}}; } });
  touch(p); rerender(); note("C 1s 템플릿: 간격 +1.5 / +3.0 / +4.0 eV 는 흔히 쓰는 출발값입니다(문헌마다 ±0.3 eV). 시료에 맞게 확인하세요");
}
function delComp(id){ const p=pr(); pushHist(p,"성분 삭제"); p.comps=p.comps.filter(c=>c.id!==id); delete p.cons[id]; Object.values(p.cons).forEach(o=>Object.values(o).forEach(cs=>{ if(cs.mode==="link"&&cs.ref===id) cs.mode="free"; })); touch(p); rerender(); }

/* ============================================================
   피팅 실행 (Web Worker)
   ============================================================ */
function workerURL(){
  if(WORK&&WORK.url) return WORK.url;
  const el=document.getElementById("xps-core-src"); const src=el?el.textContent:null;
  if(!src) return null;
  /* 엔진 순서: core → multi → arxps → image → proc (각 파일이 self 에 전역을 단다). xj = Series·ARXPS·Imaging·Processing 작업 */
  const more=["xps-multi-src","xps-arxps-src","xps-image-src","xps-proc-src"].map(id=>{ const q=document.getElementById(id); return q?q.textContent:""; }).join("\n;\n");
  const body=src+"\n;\n"+more+"\nconst XJ=("+xjobs.toString()+")(self);\n"+`
self.onmessage=function(e){ const d=e.data; const C=self.XPS_CORE; try{
  let out;
  if(d.cmd==="fit") out=C.fitJob(d.prob,d.starts,d.opt);
  else if(d.cmd==="final") out=C.finalJob(d.prob,d.opt);
  else if(d.cmd==="mc") out=C.mcJob(d.prob,d.res,d.M,d.seed,d.opt);
  else if(d.cmd==="auto") out=C.autoJob(d.prob,d.opt);
  else if(d.cmd==="sens") out=C.sensJob(d.spec,d.model,d.variants);
  else if(d.cmd==="batch") out=C.batchJob(d.spec,d.model,d.opt);
  else if(d.cmd==="xj") out=XJ[d.job](d.args,(fr,t)=>self.postMessage({id:d.id,prog:fr,txt:t||""}));
  self.postMessage({id:d.id,ok:true,out});
}catch(err){ self.postMessage({id:d.id,ok:false,err:String(err&&err.message||err)}); } };`;
  const url=URL.createObjectURL(new Blob([body],{type:"text/javascript"})); WORK={url}; return url;
}
function nWorkers(){ const hc=(navigator&&navigator.hardwareConcurrency)||2; return Math.max(1,Math.min(16,hc-1)); }
function pool(jobs,onDone){
  const url=workerURL(); const K=Math.min(nWorkers(),jobs.length);
  return new Promise((resolve,reject)=>{
    if(!url){ // 주 스레드 (예외: 주입 안 된 시험 환경)
      try{ const outs=jobs.map(j=>{ const o=j.cmd==="fit"?C.fitJob(j.prob,j.starts,j.opt):j.cmd==="final"?C.finalJob(j.prob,j.opt):j.cmd==="auto"?C.autoJob(j.prob,j.opt):j.cmd==="sens"?C.sensJob(j.spec,j.model,j.variants):j.cmd==="batch"?C.batchJob(j.spec,j.model,j.opt):C.mcJob(j.prob,j.res,j.M,j.seed,j.opt); onDone&&onDone(); return o; }); resolve(outs); }catch(e){ reject(e); } return; }
    const ws=[...Array(K)].map(()=>new Worker(url)); RUN={ws,cancel:()=>{ ws.forEach(w=>w.terminate()); reject(new Error("취소")); }};
    const outs=new Array(jobs.length); let next=0, done=0, failed=null;
    const feed=w=>{ if(next>=jobs.length) return; const id=next++; w.postMessage(Object.assign({id},jobs[id])); };
    ws.forEach(w=>{ w.onmessage=e=>{ const d=e.data; if(!d.ok){ failed=d.err; } else outs[d.id]=d.out; done++; onDone&&onDone(); if(done===jobs.length){ ws.forEach(x=>x.terminate()); RUN=null; if(failed&&outs.every(o=>!o)) reject(new Error(failed)); else resolve(outs); } else feed(w); };
      w.onerror=e=>{ e.preventDefault&&e.preventDefault(); failed=e.message||"Worker 오류"; done++; if(done===jobs.length){ ws.forEach(x=>x.terminate()); RUN=null; reject(new Error(failed)); } else feed(w); };
      feed(w); });
  });
}
function probOf(sp,p){ const P=C.makeProblem(sp,{lo:p.lo,hi:p.hi,bg:p.bg,comps:p.comps,cons:p.cons,weight:p.weight,optimizer:p.optimizer}); if(p.bg.type==="active"){ /* c0,k: 현재 값이 없으면 엔진이 초기화 */ } return P; }
function setProg(txt,frac){ const el=ROOT&&$q(ROOT,"#xpsProg"); if(!el) return; el.hidden=false; $q(el,".t").textContent=txt; const ind=frac==null||!isFinite(frac); $q(el,".bar").classList.toggle("ind",ind); $q(el,".bar i").style.width=ind?"30%":Math.round(100*Math.min(1,frac))+"%"; }
function hideProg(){ const el=ROOT&&$q(ROOT,"#xpsProg"); if(el) el.hidden=true; }
async function runFit(kind,spArg){
  const sp=spArg&&spArg.uid?spArg:cur(); if(!sp) return; const p=pr(sp);
  if(!p.comps.length){ note("성분이 없습니다 — 'Auto suggest' 또는 그래프에서 '+ Add component'","warn"); return; }
  if(RUN){ note("이미 계산 중입니다","warn"); return; }
  let prob; try{ prob=probOf(sp,p); }catch(e){ note(e.message,"warn"); return; }
  pushHist(p,kind==="multi"?"다중 시작 피팅":"피팅"); const hEnt=p.hist[p.hist.length-1];
  const t0=performance.now(); const ver=p.ver;   // 계산 중 사용자가 성분·영역을 바꾸면 결과를 덮어쓰지 않는다
  try{
    const N=kind==="multi"?Math.max(2,+S.ui.starts||24):1;
    const starts=C.startPoints(prob,N,(+S.ui.seed||1));
    const K=Math.min(nWorkers(),N); const chunks=[...Array(K)].map(()=>[]); starts.forEach((s,i)=>chunks[i%K].push(s));
    let dn=0; setProg(N>1?`다중 시작 0/${K} 묶음`:"피팅 중…",0.02);
    const outs=await pool(chunks.map(st=>({cmd:"fit",prob,starts:st,opt:{}})),()=>{ dn++; setProg(`시작점 묶음 ${dn}/${K}`,dn/(K+1)); });
    const ok=outs.filter(Boolean); if(!ok.length) throw new Error("피팅 실패");
    const best=ok.reduce((a,b)=>b.chi2<a.chi2?b:a); const chis=[].concat(...ok.map(o=>o.chis));
    prob.comps=best.comps; prob.bg=best.bg;
    setProg("불확도 계산…",K/(K+1));
    const [fin]=await pool([{cmd:"final",prob,opt:{}}]);
    const near=chis.filter(x=>x-best.chi2<1).length;
    fin.multi=N>1?{N,chis,near,frac:near/N}:null;
    if(p.ver!==ver||S.proc[sp.uid]!==p){ hideProg(); const hi=(p.hist||[]).indexOf(hEnt); if(hi>=0) p.hist.splice(hi,1); note(`피팅 중에 ${sp.region} 의 모형이 바뀌어 결과를 버렸습니다 — 다시 Fit 하세요`,"warn"); return; }
    fin.iter0=best.iter||0;   // 시작점 묶음 피팅의 반복 수(최종 다듬기 fin.iter 와 따로)
    p.comps=fin.comps; p.bg=Object.assign(p.bg,{c0:fin.bg.c0,k:fin.bg.k,s:fin.bg.s});
    p.result=Object.assign(fin,{when:new Date().toISOString(),ms:Math.round(performance.now()-t0),stale:false,lo:p.lo,hi:p.hi,bgType:p.bg.type,mc:null});
    p.ver++; save(); hideProg(); rerender();
    const st=fin.stats; note(`피팅 완료 — χ²ν ${f(st.chi2nu,3)}, ${fin.conv?"수렴":"최대 반복 도달"}${fin.multi?`, ${N}개 시작 중 ${near}개가 같은 최솟값`:""}`,fin.conv?"ok":"warn");
  }catch(e){ hideProg(); RUN=null; note("피팅 중단: "+(e.message||e),"warn"); }
}
async function runMC(){
  const sp=cur(), p=pr(sp); if(!p||!p.result||p.result.stale){ note("먼저 피팅하세요 (수정 후에는 다시 피팅)","warn"); return; }
  if(RUN){ note("이미 계산 중입니다","warn"); return; }
  const M=Math.max(20,+S.ui.mcN||200); let prob; try{ prob=probOf(sp,p); }catch(e){ note(e.message,"warn"); return; }
  const K=Math.min(nWorkers(),Math.ceil(M/10)); const per=Math.ceil(M/K); let dn=0;
  setProg(`Monte Carlo 0/${M}`,0.01);
  try{
    const res={mu:p.result.mu,v:p.result.v,chi2:p.result.chi2};
    const outs=await pool([...Array(K)].map((_,k)=>({cmd:"mc",prob,res,M:Math.min(per,M-k*per),seed:1000+k*7919,opt:{}})),()=>{ dn++; setProg(`Monte Carlo 묶음 ${dn}/${K}`,dn/K); });
    const samples=[].concat(...outs.filter(Boolean).map(o=>o.samples)); const fails=outs.reduce((s,o)=>s+(o?o.fails:0),0);
    p.result.mc={summ:C.mcSummarize(samples),n:samples.length,fails,rho:outs[0]?outs[0].rho:null,phi:outs[0]?outs[0].phi:null,when:new Date().toISOString(),
      areaSum:samples.length?areaSumStats(samples):null};
    save(); hideProg(); S.ui.tab="result"; rerender(); note(`Monte Carlo ${samples.length}회 완료${fails?` (실패 ${fails})`:""}`,"ok");
  }catch(e){ hideProg(); RUN=null; note("Monte Carlo 중단: "+(e.message||e),"warn"); }
}
function areaSumStats(samples){ const a=samples.map(s=>s.reduce((t,q)=>t+q.area,0)); const m=a.reduce((x,y)=>x+y,0)/a.length; return {mean:m,sd:Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/Math.max(1,a.length-1))}; }
function stopRun(){ if(RUN){ try{ RUN.cancel(); }catch(e){} RUN=null; hideProg(); note("계산을 중단했습니다"); } }


/* ============================================================
   기준 스펙트럼 · 측정 선모양 · LCF (SPEC §12)
   ============================================================ */
const REFS=()=>S.refs||(S.refs=[]);
function refCompat(sp,t){ const r=sp.reg||{}; return !r.el||!t.el||(t.el===r.el&&(t.orb||"")===(r.orb||"")); }
/* 지금 스펙트럼의 피팅 창 − 지금 배경 → 넓이 1 template (보정 BE). active 배경은 피팅 결과의 배경(없으면 Smart), None 은 0 */
function storeRef(sp){
  const p=pr(sp); const R=region(sp,p); if(!R||R.E.length<8){ note("영역 안 데이터가 너무 적습니다","warn"); return null; }
  const E=Float64Array.from(R.E), y=Float64Array.from(R.y), t=p.bg.type;
  const o={shift:shOf(sp),interp:S.ui.tplInterp==="linear"?"linear":"pchip",name:`${smpDisp(SK(sp)).name} · ${sp.region}`,noise:p.weight==="const"?"const":"poisson",c:sp.c||null};
  let bgT=t;
  if(t==="active"){ const r=p.result; if(r&&!r.stale&&r.B&&r.B.length===E.length) o.B=r.B; else { o.bg={type:"smart",n:p.bg.n||5}; bgT="smart"; } }
  else if(t==="none") o.B=new Float64Array(E.length);
  else o.bg=C.clone(bgOf(sp,p,R));
  let tp; try{ tp=C.makeTemplate(E,y,o); }catch(e){ note(e.message,"warn"); return null; }
  Object.assign(tp,{src:sp.uid,sample:SK(sp),region:sp.region,el:(sp.reg||{}).el||"",orb:(sp.reg||{}).orb||"",bgType:C.BG_TYPES[bgT]||bgT,when:new Date().toISOString()});
  REFS().push(tp); save(); renderPanel(); note(`기준 저장: ${tp.name} — 배경(${tp.bgType})을 빼고 넓이 1, 기준 위치 ${f(tp.eref,2)} eV`,"ok"); return tp;
}
function addRefComp(sp,t,lcf){
  const p=pr(sp), sh=shOf(sp), be=t.eref-sh; const M=curModel(sp,p);
  const c=C.newComp("MEAS",be,{tpl:C.clone(t),el:t.el||(sp.reg||{}).el||"",line:t.orb||(sp.reg||{}).orb||""});
  let vm=0; t.v.forEach(v=>{ if(v>vm) vm=v; }); const h=M?Math.max(1,interp(M.E,M.y,be)-interp(M.E,M.B,be)):1;
  c.area=Math.max(1e-6,h/(vm||1)/(lcf?2:1)); c.name=(t.name||"Ref").slice(0,40); c.slot=slotFor(p);
  p.cons[c.id]={be:lcf?{mode:"fixed"}:{mode:"free",lo:be-1.5,hi:be+1.5},wG:{mode:"fixed"}};
  p.comps.push(c); S.ui.selComp=c.id; return c;
}
function lcfAll(sp){
  const p=pr(sp); const have=new Set(p.comps.filter(c=>c.shape==="MEAS"&&c.tpl).map(c=>c.tpl.id));
  const L=REFS().filter(t=>refCompat(sp,t)&&!have.has(t.id)); if(!L.length){ note("더할 기준이 없습니다 — 같은 원소·준위의 스펙트럼에서 'Store as reference'","warn"); return; }
  pushHist(p,"LCF"); L.forEach(t=>addRefComp(sp,t,true)); touch(p); rerender();
  note(`LCF: 기준 ${L.length}개를 고정 위치로 넣었습니다 — Fit 후 Area % 가 분율. 이동·넓힘은 BE·Broadening 을 Free 로`,"ok");
}
function refHTML(sp,p){
  const L=REFS(), ok=L.filter(t=>refCompat(sp,t)), other=L.length-ok.length;
  const rows=ok.map(t=>`<tr><td title="${esc(`${t.sample} · ${t.region} · background ${t.bgType} · ${t.interp==="linear"?"linear":"PCHIP"} · ${t.x.length} pts`)}">${esc(t.name)}</td><td class="num">${f(t.eref,2)}</td>
    <td style="white-space:nowrap"><button class="btn xs" data-refl="${t.id}" title="Add as a fixed-position reference (LCF: area only)">+ LCF</button><button class="btn xs" data-refm="${t.id}" title="Add as a measured line-shape component with free position (CasaXPS 'line shape from spectrum')">+ Shape</button><button class="btn xs" data-refx="${t.id}" title="Delete this reference (press twice)">×</button></td></tr>`).join("");
  return `<details class="xps-d" id="xpsRefD" ${S.ui.refOpen?"open":""}><summary title="Reference spectra: linear-combination fitting (Avantage reference levels) and measured line shapes (CasaXPS)">Reference spectra · LCF <span class="muted xs">(${ok.length}${other?` + ${other} other`:""})</span></summary><div class="b">
   <div class="xps-row"><button class="btn xs" id="xpsRefS" title="Store this spectrum's fit range, minus its current background and normalized to area 1, as a reference. Use it in any sample with the same element and level">Store as reference</button>
    <select class="xps-sel" id="xpsRefI" title="How a stored reference is read between its points — PCHIP: shape-preserving cubic (scipy PchipInterpolator); Linear. Interpolation error ≈ 6e-4 of the peak at 0.05 eV steps / 1 eV FWHM"><option value="pchip" ${S.ui.tplInterp!=="linear"?"selected":""}>PCHIP</option><option value="linear" ${S.ui.tplInterp==="linear"?"selected":""}>Linear</option></select>
    <button class="btn xs primary" id="xpsLcf" ${ok.length?"":"disabled"} title="Add every compatible reference at a fixed position (fractions = Area %). Free BE / Broadening for shift / extra width">LCF: add all</button></div>
   ${ok.length?`<div class="tw"><table class="xps-tbl"><thead><tr><th>Reference</th><th class="num" title="Charge-corrected BE of the reference maximum">BE max</th><th></th></tr></thead><tbody>${rows}</tbody></table></div>`:`<div class="muted xs">같은 원소·준위의 기준 스펙트럼이 없습니다.</div>`}
  </div></details>`;
}
function bindRef(el,sp,p){
  const d=$q(el,"#xpsRefD"); if(!d) return; d.ontoggle=()=>{ S.ui.refOpen=d.open; };
  $q(el,"#xpsRefS").onclick=()=>storeRef(sp);
  $q(el,"#xpsRefI").onchange=e=>{ S.ui.tplInterp=e.target.value; save(); };
  $q(el,"#xpsLcf").onclick=()=>lcfAll(sp);
  $a(el,"[data-refl],[data-refm]").forEach(b=>b.onclick=()=>{ const id=b.dataset.refl||b.dataset.refm; const t=REFS().find(x=>x.id===id); if(!t) return; pushHist(p,"기준 추가"); addRefComp(sp,t,!!b.dataset.refl); touch(p); rerender(); });
  $a(el,"[data-refx]").forEach(b=>b.onclick=()=>{ if(!confirmInline("refx"+b.dataset.refx,"다시 누르면 기준을 지웁니다(이미 넣은 성분은 남음)")) return; S.refs=REFS().filter(x=>x.id!==b.dataset.refx); save(); renderPanel(); });
}

/* ============================================================
   Spline 매듭 (그래프)
   ============================================================ */
function knotMarks(sp,p,M,X,Y,sh){
  const E=M.E, y=M.y; let o="";
  C.splineKnots(E,y,Object.assign({},p.bg,{knots:[]})).forEach(k=>{ o+=`<circle class="knot0" cx="${X(k.e+sh).toFixed(1)}" cy="${Y(k.y).toFixed(1)}" r="3.5"><title>Window-end knot (endpoint average)</title></circle>`; });
  (p.bg.knots||[]).forEach((k,i)=>{ const q=C.splineKnots(E,y,{ends:false,knots:[k],kh:p.bg.kh})[0]; if(!q) return; const x=X(q.e+sh), yy=Y(q.y);
    o+=`<rect class="knot${k.y!=null?" man":""}" data-knot="${i}" x="${(x-5).toFixed(1)}" y="${(yy-5).toFixed(1)}" width="10" height="10"><title>Knot ${i+1}: ${f(q.e+sh,2)} eV${k.y!=null?" (dragged height)":" (auto height)"} — drag to move, double-click to remove</title></rect>`; });
  return o;
}
function addKnot(be){ const sp=cur(), p=pr(sp); if(p.bg.type!=="spline") return; pushHist(p,"매듭 추가"); (p.bg.knots=p.bg.knots||[]).push({e:+(+be).toFixed(4)}); S.ui.knotMode=false; touch(p); rerender(); }
function delKnot(i){ const p=pr(); if(!p||!p.bg.knots||!p.bg.knots[i]) return; pushHist(p,"매듭 지움"); p.bg.knots.splice(i,1); touch(p); rerender(); }

/* ============================================================
   여러 스펙트럼에 같은 모형 (SPEC §15)
   ============================================================ */
function batchCands(sp){ return S.specs.filter(s=>s.uid!==sp.uid&&s.reg&&s.reg.kind!=="survey"); }
function batchSel(sp){ const c=batchCands(sp); if(!S.ui.batSel||S.ui.batSrc!==sp.uid){ S.ui.batSrc=sp.uid; S.ui.batSel={}; const r=sp.reg||{}; c.forEach(s=>{ if(r.el&&s.reg.el===r.el&&s.reg.orb===r.orb) S.ui.batSel[s.uid]=true; }); } return c.filter(s=>S.ui.batSel[s.uid]); }
function batchPanelHTML(sp,p){
  const c=batchCands(sp); if(!c.length) return ""; const sel=new Set(batchSel(sp).map(s=>s.uid));
  const groups={}; c.forEach(s=>{ (groups[SK(s)]=groups[SK(s)]||[]).push(s); });
  return `<details class="xps-d" id="xpsBatD" ${S.ui.batOpen?"open":""}><summary title="Copy this peak model — components, constraints, background (incl. knots), fit range, weighting, optimizer — to other spectra in charge-corrected BE, then fit them all (e.g. a depth profile)">Apply model to other spectra <span class="muted xs">(${sel.size} selected)</span></summary><div class="b">
   <div class="xps-batl">${Object.entries(groups).map(([k,a])=>`<div class="muted xs" style="margin-top:3px">${esc(smpDisp(k).name)}</div>${a.map(s=>`<label class="xps-f" style="margin-right:8px"><input type="checkbox" data-bat="${s.uid}" ${sel.has(s.uid)?"checked":""}> ${esc(s.region)}</label>`).join("")}`).join("")}</div>
   <div class="xps-row" style="margin-top:6px"><button class="btn xs" id="xpsBatA" ${sel.size?"":"disabled"} title="Replace the model of the selected spectra (undo per spectrum)">Apply</button><button class="btn xs primary" id="xpsBatF" ${sel.size&&p.comps.length?"":"disabled"} title="Apply, then fit every selected spectrum in parallel (areas restarted by NNLS). Results → Results tab, Batch results">Apply &amp; fit all</button></div>
  </div></details>`;
}
function bindBatchPanel(el,sp,p){
  const d=$q(el,"#xpsBatD"); if(!d) return; d.ontoggle=()=>{ S.ui.batOpen=d.open; };
  $a(el,"[data-bat]").forEach(i=>i.onchange=()=>{ S.ui.batSel[i.dataset.bat]=i.checked; renderPanel(); });
  $q(el,"#xpsBatA").onclick=()=>runBatch(sp,false); $q(el,"#xpsBatF").onclick=()=>runBatch(sp,true);
}
async function runBatch(sp,fitAll){
  const p=pr(sp); if(!p.comps.length){ note("성분이 없습니다","warn"); return; } if(RUN){ note("이미 계산 중입니다","warn"); return; }
  const tg=batchSel(sp); if(!tg.length){ note("적용할 스펙트럼을 고르세요","warn"); return; }
  const model={lo:p.lo,hi:p.hi,bg:p.bg,comps:p.comps,cons:p.cons,weight:p.weight,optimizer:p.optimizer}; const shS=shOf(sp);
  tg.forEach(t=>{ const q=pr(t); pushHist(q,"모형 적용 ("+sp.region+")"); const md=C.propagateModel(model,shS,shOf(t),[t.be[0],t.be[t.n-1]]);
    q.lo=md.lo; q.hi=md.hi; q.bg=md.bg; q.comps=md.comps; q.cons=md.cons; q.weight=md.weight; q.optimizer=md.optimizer; touch(q); });
  if(!fitAll){ save(); rerender(); note(`모형을 ${tg.length}개 스펙트럼에 적용했습니다(대전 보정 BE 기준) — 각 스펙트럼에서 Fit`,"ok"); return; }
  const jobs=[], idx=[]; tg.forEach(t=>{ const q=pr(t); if(!(q.hi-q.lo>0.5)) return; idx.push(t); jobs.push({cmd:"batch",spec:{be:t.be,y:t.y,c:t.c||null},model:{lo:q.lo,hi:q.hi,bg:q.bg,comps:q.comps,cons:q.cons,weight:q.weight,optimizer:q.optimizer},opt:{}}); });
  const vers=idx.map(t=>pr(t).ver); let dn=0; const t0=performance.now(); setProg(`Batch fit 0/${jobs.length}`,0.01);
  try{
    const outs=await pool(jobs,()=>{ dn++; setProg(`Batch fit ${dn}/${jobs.length}`,dn/jobs.length); });
    let nOk=0; const bad=[];
    idx.forEach((t,i)=>{ const q=pr(t), o=outs[i]; if(!o){ bad.push(t.region); return; } if(q.ver!==vers[i]) return;
      q.comps=o.comps; q.bg=Object.assign(q.bg,{c0:o.bg.c0,k:o.bg.k,s:o.bg.s});
      q.result=Object.assign(o,{when:new Date().toISOString(),ms:Math.round(performance.now()-t0),stale:false,lo:q.lo,hi:q.hi,bgType:q.bg.type,mc:null,batch:true}); q.ver++; nOk++; });
    S.ui.batch={src:sp.uid,uids:idx.map(t=>t.uid),when:new Date().toISOString()};
    save(); hideProg(); S.ui.tab="result"; rerender(); note(`Batch fit ${nOk}/${tg.length}${bad.length?` — 실패: ${bad.join(", ")}`:""} · ${((performance.now()-t0)/1000).toFixed(1)} s`,bad.length?"warn":"ok");
  }catch(e){ hideProg(); RUN=null; note("Batch fit 중단: "+(e.message||e),"warn"); }
}
/* 일괄 결과 표: 스펙트럼마다 χ²ν · 성분(이름으로 맞춤) BE·FWHM·넓이·% · 그 원소의 at %(시료 정량) */
function batchRows(){
  const b=S.ui.batch; if(!b) return null; const list=[b.src].concat(b.uids).map(byUid).filter(Boolean); if(!list.length) return null;
  const names=[]; list.forEach(s=>{ const r=fitted(s); if(r) r.comps.forEach(c=>{ const nm=c.name||c.id; if(!names.includes(nm)) names.push(nm); }); });
  const atCache={};
  const atOf=s=>{ const k=SK(s); if(!(k in atCache)){ try{ const rows=quantRows(k), used=rows.filter(r=>r.inc&&r.S>0); const Q=used.length?C.quantify(used.map(r=>({label:r.el,el:r.el,area:r.area,sArea:r.sArea,S:r.S,BE:r.BE,hv:r.hv,PE:r.PE,tf:r.tf})),{imfp:S.quant.imfp}):[]; atCache[k]={used,Q}; }catch(e){ atCache[k]={used:[],Q:[]}; } }
    const a=atCache[k]; const i=a.used.findIndex(r=>r.sp.uid===s.uid); return i>=0?a.Q[i]:null; };
  return {names,rows:list.map(s=>{ const r=fitted(s), sh=shOf(s); return {s,sh,r,at:r?atOf(s):null,comp:nm=>{ if(!r) return null; const i=r.comps.findIndex(c=>(c.name||c.id)===nm); return i<0?null:{c:r.comps[i],q:r.summ[i]}; }}; })};
}
function batchHTML(sp){
  const B=batchRows(); if(!B) return ""; const b=S.ui.batch;
  const head=`<tr><th>Sample · region</th><th class="num">χ²ν</th>${B.names.map(n=>`<th class="num">${esc(n)}<div class="muted xs">BE · FWHM · %</div></th>`).join("")}<th class="num" title="at % of this spectrum's element in its sample (Quantification tab settings)">at %</th></tr>`;
  const body=B.rows.map(x=>`<tr${x.s.uid===sp.uid?` class="tot"`:""}><td><span title="${esc(SK(x.s))}">${esc(smpDisp(SK(x.s)).name)}</span> · ${esc(x.s.region)}</td><td class="num">${x.r?f(x.r.stats.chi2nu,2):`<span class="badge wr" title="No current fit (modified or failed)">—</span>`}</td>
    ${B.names.map(n=>{ const v=x.comp(n); return `<td class="num">${v?`${f(v.q.be+x.sh,2)} · ${f(v.q.fwhm,2)}<div class="muted xs">${f(v.q.areaPct,1)}${v.q.sig&&v.q.sig.areaPct!=null?" ± "+f(v.q.sig.areaPct,1):""} %</div>`:"—"}</td>`; }).join("")}
    <td class="num">${x.at&&x.at.at!=null?`${f(x.at.at,2)}${x.at.sAt?`<div class="muted xs">± ${f(x.at.sAt,2)}</div>`:""}`:"—"}</td></tr>`).join("");
  return `<details class="xps-d" open><summary>Batch results <span class="muted xs">(${B.rows.length} spectra · ${esc(String(b.when).slice(0,16).replace("T"," "))})</span></summary><div class="b"><div class="tw"><table class="xps-tbl xps-bat">${head}${body}</table></div>
   <div class="xps-row"><button class="btn xs" id="xpsBatCsv">Batch CSV</button><span class="muted xs">BE 는 시료마다의 대전 보정 값. ± = 공분산 1σ</span></div></div></details>`;
}
function bindBatch(el,sp){ const b=$q(el,"#xpsBatCsv"); if(b) b.onclick=exportBatch; }
function exportBatch(){ const B=batchRows(); if(!B) return;
  const cols=[]; B.names.forEach(n=>["BE_corr_eV","sBE","FWHM_eV","area","sArea","area_pct","s_area_pct"].forEach(k=>cols.push(n+"_"+k)));
  const lines=[["sample","region","element","chi2nu","charge_shift_eV",...cols,"at_pct","s_at_pct"].map(csvE).join(",")];
  B.rows.forEach(x=>{ const v=[SK(x.s),x.s.region,(x.s.reg||{}).el||"",x.r?x.r.stats.chi2nu.toFixed(4):"",x.sh];
    B.names.forEach(n=>{ const c=x.comp(n); if(!c){ v.push("","","","","","",""); return; } const q=c.q, sg=q.sig||{}; v.push((q.be+x.sh).toFixed(3),sg.be!=null?sg.be.toFixed(4):"",q.fwhm!=null?q.fwhm.toFixed(3):"",q.area.toFixed(2),sg.area!=null?sg.area.toFixed(2):"",q.areaPct!=null?q.areaPct.toFixed(3):"",sg.areaPct!=null?sg.areaPct.toFixed(3):""); });
    v.push(x.at&&x.at.at!=null?x.at.at.toFixed(3):"",x.at&&x.at.sAt?x.at.sAt.toFixed(3):""); lines.push(v.map(csvE).join(",")); });
  download("XPS_batch_"+new Date().toISOString().slice(0,10)+".csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"})); }


/* CasaXPS 규약의 영역 숫자(Area · Max · Min Height) — CasaXPS 영역 표와 바로 비교 */
function casaRegionTxt(sp,p){ if(p.bg.conv!=="casa"||!(p.bg.type==="linear"||p.bg.type==="shirley")) return ""; const R=region(sp,p); if(!R||R.E.length<5) return "";
  const E=Float64Array.from(R.E), y=Float64Array.from(R.y), B=C.bgStatic(E,y,bgOf(sp,p,R)); const st=C.casaRegionStats(E,y,B);
  return `<div class="xps-casa xs" title="Same definitions as the CasaXPS Regions table (trapezoid of data − background over channels Start+1 … End−1; Max/Min of data − background there)">CasaXPS region: Area ${g(st.area,6)} · Max ${f(st.maxH,1)} · Min ${f(st.minH,1)}</div>`; }
/* 같은 모양의 CasaXPS 표기 (SPEC §13): GLP = GL(100m), CasaXPS 의 fwhm = 실제 FWHM. LA 는 모양·꼬리 쪽이 같고 폭 숫자의 환산은 미확인 */
function casaEqTxt(c,sm){ return casaEq1(c,sm)+avEqTxt(c,sm); }
/* Avantage 표기 (SPEC §13b, Thermo Avantage manual: Product = Sherwood GL 곱, β ≈ ½FWHM → F = 2β; Convolution 은 (fV, M) 표시) */
function avEqTxt(c,sm){
  if(c.shape==="GLP") return `<div class="muted xs" title="Avantage G-L Product: exp(−(1−M)ln2·x²/β²)/(1+M·x²/β²) with β = F/2, M = L/G Mix/100 — same function. Avantage's peak-table FWHM is the measured (true) FWHM">Avantage: Product, L/G Mix ${+(100*c.m).toFixed(1)} %, β ${f(c.F/2,4)} eV</div>`;
  if(c.shape==="V"){ const a=C.avantageVoigt(c.wG,c.wL); return `<div class="muted xs" title="Avantage G-L Convolution shows fV (Olivero–Longbothum FWHM) and a mix M; the fL/fG(M) formula is an image in the manual — M = fL/(fL+fG) assumed, pending Avantage check">Avantage: Convolution fV ${f(a.fV,4)} eV, L/G Mix ≈ ${f(100*a.M,1)} % (assumed)</div>`; }
  if(c.shape==="AGT") return `<div class="muted xs" title="Avantage tail function on a G-L Product: H[GL + (1−GL)T], T = TM·CT + (1−TM)·exp(−ET·Dx), Dx in channels (here ${f(c.ch,3)} eV per channel)">Avantage: Product + tail, ET ${g(c.ET,3)} · TM ${g(c.TM,3)} · CT ${g(c.CT,3)}</div>`;
  return ""; }
function chOf(sp){ return sp&&sp.n>1?Math.abs((sp.be[sp.n-1]-sp.be[0])/(sp.n-1)):0.05; }
function casaEq1(c,sm){ let q; try{ q=C.compToCasa(c); }catch(e){ return ""; }
  if(!q||!q.notation) return "";
  const fw=c.shape==="GLP"&&sm&&sm.fwhm?sm.fwhm:q.fwhm;
  return `<div class="muted xs" style="margin-top:4px" title="${esc(q.note+(q.exact?"":" — approximate"))}">CasaXPS: ${esc(q.notation)}${fw!=null?`, fwhm ${f(fw,4)}`:""}${q.exact?"":" ≈"}</div>`; }
/* ============================================================
   화면 뼈대
   ============================================================ */
function render(root){
  ROOT=root;
  if(!LOADED){ LOADED=true; root.innerHTML=`<div class="card"><div class="xps-empty">XPS 작업을 불러오는 중…</div></div>`; dbGet().then(o=>{ if(o&&o.v===1){ const r=revive(o); S=Object.assign(defState(),r); S.ui=Object.assign(defState().ui,r.ui||{}); S.quant=Object.assign(defState().quant,r.quant||{}); } if(ROOT===root) build(root); }); return; }
  build(root);
}
function build(root){
  ROOT=root;
  root.innerHTML=`<div class="xps">
    <div class="card xps-top">
      <div class="grp"><label class="btn primary sm" style="cursor:pointer" title="VGD · VAMAS (.vms) · Avantage xls/xlsx/avg · Text. Multiple files or folders can also be dropped onto the page">Open files<input type="file" id="xpsFin" multiple accept="${ACCEPT}" style="display:none"></label><span class="muted sm">또는 화면에 끌어다 놓기</span></div>
      <div class="seg xps-ws-sw" id="xpsWs" role="tablist" aria-label="XPS workspace"></div>
      <div class="grp" style="margin-left:auto"><button class="btn sm" id="xpsSaveP" title="Save spectra, fits and charge correction as one JSON file">Save work</button><label class="btn sm" style="cursor:pointer" title="Open a saved JSON work file">Open work<input type="file" id="xpsLoadP" accept=".json" style="display:none"></label><button class="btn sm" id="xpsClear" title="Remove the list and all work (press twice to confirm)">Clear all</button></div>
      <div class="xps-chg" id="xpsChg"></div>
    </div>
    <div class="xps-main">
      <div class="card xps-list" id="xpsList"></div>
      <div class="card xps-center"><div class="xps-plotbar" id="xpsBar"></div><div id="xpsPlot"></div>
        <div class="xps-prog" id="xpsProg" hidden><span class="t"></span><div class="bar"><i></i></div><button class="btn xs" id="xpsStop">Stop</button></div>
        <div class="xps-msg" id="xpsMsg" hidden></div></div>
      <div class="card xps-side"><div class="seg xps-tabs" id="xpsTabs"></div><div id="xpsPanel"></div></div>
    </div><div class="xps-ws" id="xpsWsBox" hidden></div></div>`;
  const x=$q(root,".xps");
  $q(root,"#xpsFin").onchange=e=>{ loadFiles(e.target.files); e.target.value=""; };
  $q(root,"#xpsLoadP").onchange=e=>{ const fl=e.target.files[0]; if(fl) fl.text().then(importProject); e.target.value=""; };
  $q(root,"#xpsSaveP").onclick=exportProject;
  $q(root,"#xpsClear").onclick=()=>{ if(!S.specs.length&&!IMS.stacks.length) return; if(!confirmInline("xpsClear","정말 지울까요? 다시 누르면 지웁니다")) return; const ui=S.ui; S=defState(); S.ui.tab=ui.tab; S.ui.ws=ui.ws; PV={}; Object.keys(RES).forEach(k=>RES[k]=null); Object.keys(PRES).forEach(k=>delete PRES[k]); IMS.stacks=[]; imgSave(); save(); rerender(); };
  $q(root,"#xpsStop").onclick=stopRun;
  ["dragenter","dragover"].forEach(ev=>x.addEventListener(ev,e=>{ if(e.dataTransfer&&[...(e.dataTransfer.types||[])].includes("Files")){ e.preventDefault(); e.stopPropagation(); x.classList.add("over"); } }));
  x.addEventListener("dragleave",e=>{ if(e.target===x) x.classList.remove("over"); });
  x.addEventListener("drop",e=>{ if(e.dataTransfer&&e.dataTransfer.files&&e.dataTransfer.files.length){ e.preventDefault(); e.stopPropagation(); x.classList.remove("over"); loadFiles(e.dataTransfer.files); } });
  rerender();
}
const _cf={};
function confirmInline(k,msg){ const now=Date.now(); if(_cf[k]&&now-_cf[k]<4000){ delete _cf[k]; return true; } _cf[k]=now; note(msg,"warn"); return false; }
function rerender(){ if(!ROOT||!$q(ROOT,".xps")) return; drawWsSw(); drawList(); drawBar(); drawChg(); drawTabs(); drawPlot(); renderPanel(); drawWS(); drawMsg(); }

/* ---------------- 왼쪽 목록 ---------------- */
function drawList(){
  const el=$q(ROOT,"#xpsList"); if(!el) return;
  if(!S.specs.length){ el.innerHTML=`<div class="xps-drop" id="xpsDrop">Drop files here or<br><b>click to open</b><br><span class="xs">한 측정 지점의 파일을 <b>C 1s 와 함께</b> 한꺼번에</span></div>
    <ol class="xps-steps"><li><b>Open files</b> — Survey·C 1s·관심 원소를 같이. 같은 시료끼리 자동으로 묶입니다</li><li><b>Charge correction</b> — 위의 <b>Calibrate to C 1s</b> 한 번(C 1s 를 자동으로 피팅해 맞춤)</li><li><b>Peak fit</b> — 스펙트럼을 고르고 <b>Auto suggest</b> → <b>Fit</b></li><li><b>Quantification · export</b> — Quantification·Results 탭</li></ol>`; $q(el,"#xpsDrop").onclick=()=>$q(ROOT,"#xpsFin").click(); return; }
  const groups={}; S.specs.forEach(s=>{ const k=SK(s); (groups[k]=groups[k]||[]).push(s); });
  const ov=new Set(S.ui.overlay);
  el.innerHTML=Object.entries(groups).map(([k,arr])=>{ const d=smpDisp(k); const cs=chgState(arr[0]); const ren=S.ui.renSmp===k;
    const badge=cs.st==="none"?(arr.some(isC1s)?`<span class="badge wr" title="Not charge-corrected — use ‘Calibrate to C 1s’ above">Uncorrected</span>`:`<span class="badge" title="Not charge-corrected (no C 1s in this sample)">Uncorrected</span>`)
      :cs.st==="ok"||cs.st==="manual"?`<span class="badge ok" title="Charge correction ${esc(smpOf(arr[0]).refNote||"")}">Δ ${cs.sh>0?"+":""}${f(cs.sh,2)}</span>`:`<span class="badge wr" title="${esc(cs.why||"")}">Δ ${cs.sh>0?"+":""}${f(cs.sh,2)} ⟳</span>`;
    return `<div class="smp"><div class="h" title="${esc(k)}">${ren?`<input type="text" class="xps-in" id="xpsRenI" value="${esc(k===NOSMP?"":k)}" placeholder="Sample name"><button class="btn xs" id="xpsRenOk">OK</button>`
      :`<span class="snm">${esc(d.name)}</span><button class="ren" data-ren="${esc(k)}" title="Rename sample — renaming to an existing name merges them">✎</button>${badge}`}</div>${d.ctx&&!ren?`<div class="ctx" title="${esc(k)}">${esc(d.ctx)}</div>`:""}${arr.map(s=>{ const p=S.proc[s.uid]; const r=p&&p.result; const st=r?(r.stale?`<span class="st stale">Modified</span>`:`<span class="st ok">χ²ν ${f(r.stats.chi2nu,2)}</span>`):(p&&p.comps.length?`<span class="st">${p.comps.length} comp.</span>`:"");
      return `<div class="it ${s.uid===S.sel?"sel":""}" data-u="${s.uid}"><input type="checkbox" data-ov="${s.uid}" ${ov.has(s.uid)?"checked":""} title="Add to overlay"><span class="nm">${esc(s.region)}</span><span class="fm">${esc(FMTN(s.format).replace("Avantage ",""))}</span>${st}<button class="x" data-del="${s.uid}" title="Remove from list">×</button></div>`; }).join("")}</div>`; }).join("")
    +`<div class="row" style="margin-top:6px;gap:6px;flex-wrap:wrap">${S.ui.overlay.length>1?`<span class="badge inf">Overlay ×${S.ui.overlay.length}</span> <select class="xps-sel" id="xpsOvN" title="Normalization"><option value="max">Max = 1</option><option value="area">Area = 1</option><option value="none">Raw intensity</option></select> <span class="xps-f">Offset <input type="number" class="w56" id="xpsOvO" step="0.1" value="${S.ui.ovOffset}"></span> <button class="btn xs" id="xpsOvX">Off</button>`:`<span class="muted xs">체크하면 겹쳐 봅니다 (2개 이상)</span>`}</div>`;
  $a(el,"[data-ren]").forEach(b=>b.onclick=()=>{ S.ui.renSmp=b.dataset.ren; drawList(); const i=$q(el,"#xpsRenI"); if(i){ i.focus(); i.select(); } });
  const ri=$q(el,"#xpsRenI"); if(ri){ const go=()=>{ const o=S.ui.renSmp; S.ui.renSmp=null; renameSample(o,ri.value); }; $q(el,"#xpsRenOk").onclick=go; ri.onkeydown=e=>{ if(e.key==="Enter") go(); if(e.key==="Escape"){ S.ui.renSmp=null; drawList(); } }; }
  $a(el,".it").forEach(d=>d.onclick=e=>{ if(e.target.closest("button")||e.target.tagName==="INPUT") return; S.sel=d.dataset.u; S.ui.selComp=null; save(); rerender(); });
  $a(el,"[data-del]").forEach(b=>b.onclick=()=>{ const u=b.dataset.del; if(!confirmInline("del"+u,"다시 누르면 목록에서 뺍니다")) return; S.specs=S.specs.filter(s=>s.uid!==u); delete S.proc[u]; S.ui.overlay=S.ui.overlay.filter(x=>x!==u); if(S.sel===u) S.sel=S.specs[0]?S.specs[0].uid:null; save(); rerender(); });
  $a(el,"[data-ov]").forEach(b=>b.onchange=()=>{ const u=b.dataset.ov; S.ui.overlay=b.checked?[...new Set(S.ui.overlay.concat(u))]:S.ui.overlay.filter(x=>x!==u); save(); rerender(); });
  const on=$q(el,"#xpsOvN"); if(on){ on.value=S.ui.ovNorm; on.onchange=()=>{ S.ui.ovNorm=on.value; save(); drawPlot(); }; }
  const oo=$q(el,"#xpsOvO"); if(oo) oo.onchange=()=>{ S.ui.ovOffset=+oo.value||0; save(); drawPlot(); };
  const ox=$q(el,"#xpsOvX"); if(ox) ox.onclick=()=>{ S.ui.overlay=[]; save(); rerender(); };
}
/* ---------------- 그래프 위 막대 ---------------- */
function drawBar(){
  const el=$q(ROOT,"#xpsBar"); if(!el) return; const sp=cur();
  if(!sp||S.ui.overlay.length>1){ el.innerHTML=S.ui.overlay.length>1?`<span class="ttl">Overlay</span>`:""; return; }
  const p=pr(sp);
  el.innerHTML=`<span class="ttl">${esc(sp.region)}</span><span class="muted" title="${esc(SK(sp))}">${esc(smpDisp(SK(sp)).name)}</span>
    <span class="sp"></span>
    <button class="btn xs ${S.ui.addMode?"primary":""}" id="xpsAdd" title="Click, then click on the graph where the component goes (vertical position = height). Click again to cancel">${S.ui.addMode?"Click the graph to place":"+ Add component"}</button>
    <select class="xps-sel" id="xpsAddS" title="Line shape of the new component (GLS = Avantage L/G Mix Sum)">${Object.entries(C.SHAPES).filter(([k])=>k!=="MEAS").map(([k,v])=>`<option value="${k}" ${S.ui.addShape===k?"selected":""}>${esc(k)}</option>`).join("")}</select>
    <label class="xps-f" title="Add as a spin–orbit doublet"><input type="checkbox" id="xpsAddD" ${S.ui.addDoublet?"checked":""} ${soInfo(sp.reg.el,sp.reg.orb)?"":"disabled"}> Doublet</label>
    <label class="xps-f" title="Library positions of this element's levels (charge-corrected) as vertical lines"><input type="checkbox" id="xpsRef" ${S.ui.showRef?"checked":""}> Reference lines</label>
    <button class="btn xs" id="xpsFitV" title="Zoom to the fit range">Zoom to fit range</button><button class="btn xs" id="xpsAll" title="Show the full spectrum (double-clicking the graph also works)">Reset view</button>`;
  $q(el,"#xpsAdd").onclick=()=>{ S.ui.addMode=!S.ui.addMode; drawBar(); };
  $q(el,"#xpsAddS").onchange=e=>{ S.ui.addShape=e.target.value; save(); };
  $q(el,"#xpsAddD").onchange=e=>{ S.ui.addDoublet=e.target.checked; save(); };
  $q(el,"#xpsRef").onchange=e=>{ S.ui.showRef=e.target.checked; save(); drawPlot(); };
  $q(el,"#xpsAll").onclick=()=>{ delete PV[sp.uid]; drawPlot(); };
  $q(el,"#xpsFitV").onclick=()=>{ const v=plotView(sp); const sh=shOf(sp); const w=(p.hi-p.lo)*0.04; v.x0=p.lo+sh-w; v.x1=p.hi+sh+w; clampView(sp,v); drawPlot(); };
}
/* ---------------- 대전 보정 ----------------
   시료마다 Δ 하나. 기준(sm.ref) = {uid, compId, be, name, file, target, how}: 기준 스펙트럼의 피팅된 성분.
   피팅은 측정 BE 로 하고 화면·표·내보내기에만 Δ 를 더한다 → 피팅 순서와 무관. 기준 성분의 BE 가 다시 피팅으로 바뀌면 '다시 맞추기' 표시. */
function chgCands(k){ return sampleSpecs(k).filter(isC1s); }
/* 기준 성분 자동 선택: 넓이가 C 성분 합의 10 % 이상인 것 중 가장 낮은 BE (보통 C–C/C–H; 작은 탄화물·잡음 성분은 건너뜀) */
function autoRefComp(comps){ const cc=comps.filter(c=>!c.el||c.el==="C"); const tot=cc.reduce((s,c)=>s+Math.max(0,c.area),0)||1;
  return cc.filter(c=>c.area>=0.1*tot).sort((a,b)=>a.be-b.be)[0]||cc.sort((a,b)=>a.be-b.be)[0]||null; }
function fitted(s){ const p=pr(s); return p.result&&!p.result.stale?p.result:null; }
function chgState(sp){
  const sm=smpOf(sp), sh=sm.shift||0;
  if(sm.ref){ const s=byUid(sm.ref.uid); if(!s) return {st:"stale",sh,why:"기준 스펙트럼이 목록에서 빠졌습니다"};
    const r=fitted(s); const c=r&&r.comps.find(c=>c.id===sm.ref.compId);
    if(!r) return {st:"stale",sh,why:"기준 C 1s 의 피팅이 바뀌었거나 지워졌습니다 — 다시 맞추세요"};
    if(!c) return {st:"stale",sh,why:"기준 성분이 없어졌습니다 — 다시 맞추세요"};
    if(Math.abs(c.be-sm.ref.be)>0.005) return {st:"stale",sh,why:`기준 성분 BE 가 ${f(sm.ref.be,2)} → ${f(c.be,2)} eV 로 바뀌었습니다 — 다시 맞추세요`};
    return {st:"ok",sh}; }
  if(sh) return {st:"manual",sh};
  return {st:"none",sh:0};
}
function setShift(sp,sh,ref,noteTxt){ const sm=smpOf(sp); sm.shift=+(+sh).toFixed(3); sm.ref=ref||null; sm.refNote=noteTxt||""; Object.keys(PV).forEach(k=>delete PV[k]); save(); rerender(); }
/* 결과·정량 탭 위에 붙이는 한 줄 */
function chgBanner(sp){ const cs=chgState(sp); if(cs.st==="ok"||cs.st==="manual") return "";
  if(cs.st==="none") return `<div class="note ${chgCands(SK(sp)).length?"warn":""}" style="margin:4px 0 6px">대전 보정 전 BE 입니다${chgCands(SK(sp)).length?" — 위의 <b>Calibrate to C 1s</b> 를 누르세요":" (이 시료에 C 1s 없음 — 위 <b>Reference</b> 에서 다른 기준이나 Δ 를 넣을 수 있습니다)"}.</div>`;
  return `<div class="note warn" style="margin:4px 0 6px">${esc(cs.why)}.</div>`; }
function drawChg(){
  const el=$q(ROOT,"#xpsChg"); if(!el) return; const sp=cur(); if(!sp){ el.innerHTML=""; el.hidden=true; return; } el.hidden=false;
  const k=SK(sp), sm=smpOf(sp), d=smpDisp(k), cs=chgState(sp), cands=chgCands(k);
  const refS=cands.find(s=>s.uid===sm.refUid)||cands[0]||null;
  const tg=sm.target!=null?sm.target:284.8;
  const sgn=v=>(v>0?"+":"")+f(v,2);
  let chip, btn;
  if(cs.st==="ok"){ chip=`<span class="chip ok" title="${esc(sm.refNote)}"><b>Δ ${sgn(cs.sh)} eV</b><span class="sub">${esc(sm.refNote)}</span></span>`; btn=cands.length?`<button class="btn xs" id="xpsChgGo" title="Re-read the reference C 1s and recompute Δ">Recalibrate</button>`:""; }
  else if(cs.st==="manual"){ chip=`<span class="chip ok"><b>Δ ${sgn(cs.sh)} eV</b><span class="sub">${esc(sm.refNote||"직접 입력")}</span></span>`; btn=cands.length?`<button class="btn sm" id="xpsChgGo">Calibrate to C 1s</button>`:""; }
  else if(cs.st==="stale"){ chip=`<span class="chip wr"><b>Δ ${sgn(cs.sh)} eV</b><span class="sub">${esc(cs.why)}</span></span>`; btn=cands.length?`<button class="btn primary sm" id="xpsChgGo">Recalibrate</button>`:""; }
  else { chip=`<span class="chip ${cands.length?"wr":""}"><b>Not corrected</b><span class="sub">${cands.length?"BE 가 측정값 그대로입니다":"이 시료에 C 1s 가 없습니다 — Reference 에서 다른 기준·Δ"}</span></span>`;
    btn=cands.length?`<button class="btn primary sm" id="xpsChgGo" title="If C 1s is not fitted yet: auto suggest → fit first, then set the C–C component to ${tg} eV">Calibrate to C 1s</button>`:""; }
  const open=!!S.ui.chgOpen;
  let set="";
  if(open){
    const rr=refS?fitted(refS):null; const comps=rr?rr.comps:(refS?pr(refS).comps:[]); const auto=autoRefComp(comps);
    const selC=pr(sp).comps.find(c=>c.id===S.ui.selComp);
    set=`<div class="set">
      <div class="r"><span class="k">Reference spectrum</span>${cands.length?`<select class="xps-sel" id="xcSpec">${cands.map(s=>`<option value="${s.uid}" ${refS&&s.uid===refS.uid?"selected":""}>${esc(s.region)} · ${esc(FMTN(s.format).replace("Avantage ",""))} · ${esc(s.file||s.name)}${fitted(s)?"":" (not fitted)"}</option>`).join("")}</select>${cands.length>1?`<span class="muted xs">이 시료에 C 1s 가 ${cands.length} 개 — 기준으로 쓸 것을 고르세요</span>`:""}`:`<span class="muted sm">이 시료에 C 1s 가 없습니다. 같은 측정의 C 1s 파일을 함께 열면 여기에 나타납니다.</span>`}</div>
      ${cands.length?`<div class="r"><span class="k">Reference component</span><select class="xps-sel" id="xcComp"><option value="">Auto — ${auto?esc((auto.name||"Component")+" "+f(auto.be,2)+" eV (measured)"):"set after fitting"} (lowest BE among components with ≥ 10 % area)</option>${comps.map(c=>`<option value="${c.id}" ${sm.refComp===c.id?"selected":""}>${esc(c.name||"Component")} · ${f(c.be,2)} eV (measured)</option>`).join("")}</select>${rr?"":`<span class="muted xs">C 1s 가 아직 피팅되지 않았습니다 — Calibrate to C 1s 를 누르면 자동으로 피팅합니다</span>`}</div>`:""}
      <div class="r"><span class="k">Reference BE</span><input type="number" step="0.01" class="w64" id="xcTg" value="${tg}"> eV <button class="btn xs" data-tg="284.8" title="Thermo/Avantage convention">284.8</button><button class="btn xs" data-tg="285.0" title="Other common convention">285.0</button><button class="btn xs" data-tg="84.0" title="Au 4f7/2 — when gold is used as the reference">Au 84.0</button><span class="muted xs">논문에는 쓴 값을 밝히세요</span></div>
      <div class="r"><span class="k">Other reference</span><button class="btn xs" id="xcSel" ${selC?"":"disabled"} title="Use the component selected in the current spectrum (e.g. Au 4f7/2 → 84.0, metal substrate) as the reference BE">Use selected component${selC?` (${esc(selC.name||"Component")} ${f(selC.be,2)})`:""}</button>${selC?"":`<span class="muted xs">성분 목록이나 그래프에서 성분을 먼저 고르세요</span>`}</div>
      <div class="r"><span class="k">Manual Δ</span><input type="number" step="0.01" class="w64" id="xcSh" value="${+(sm.shift||0).toFixed(3)}"> eV <button class="btn xs" id="xcShOk">Apply</button><button class="btn xs" id="xcOff" ${sm.shift||sm.ref?"":"disabled"}>Clear correction</button></div>
      <div class="muted xs" style="line-height:1.55">Δ 는 시료(<b>${esc(d.name)}</b>)의 모든 스펙트럼에 같이 적용됩니다. 피팅은 측정 BE 로 하고 그래프·표·정량·내보내기에만 Δ 를 더하므로 <b>피팅 순서와 무관</b>합니다 — 넓이·FWHM·at% 는 Δ 와 무관. adventitious C 1s 기준은 시료 일함수에 따라 수백 meV 까지 달라질 수 있습니다(Greczynski &amp; Hultman, Prog. Mater. Sci. 2020).</div>
    </div>`; }
  el.innerHTML=`<div class="bar"><span class="lbl">Charge Correction</span><span class="smpn" title="${esc(k)}">${esc(d.name)}</span>${chip}${btn||""}<button class="btn xs ${open?"on":""}" id="xpsChgSet" title="Reference spectrum, component, reference BE, manual Δ">Reference ${open?"▴":"▾"}</button></div>${set}`;
  const go=$q(el,"#xpsChgGo"); if(go) go.onclick=()=>calibrateC1s(sp);
  $q(el,"#xpsChgSet").onclick=()=>{ S.ui.chgOpen=!S.ui.chgOpen; drawChg(); };
  if(!open) return;
  // 설정을 바꾸면: 기준 C 1s 가 이미 피팅돼 있고 지금 C 1s 기준(또는 아직 보정 전)이면 곧바로 다시 맞춘다. 다른 기준(선택 성분)·직접 입력 상태는 건드리지 않음
  const reapply=()=>{ const rs=chgCands(k).find(s=>s.uid===sm.refUid)||chgCands(k)[0]; const viaC1s=!sm.ref||isC1s(byUid(sm.ref.uid));
    if(rs&&fitted(rs)&&viaC1s&&(sm.ref||cs.st==="none")) calibrateC1s(sp); else drawChg(); };
  const xs=$q(el,"#xcSpec"); if(xs) xs.onchange=()=>{ sm.refUid=xs.value; sm.refComp=null; save(); reapply(); };
  const xc=$q(el,"#xcComp"); if(xc) xc.onchange=()=>{ sm.refComp=xc.value||null; save(); reapply(); };
  const xt=$q(el,"#xcTg"); xt.onchange=()=>{ const v=+xt.value; if(isFinite(v)&&v>0){ sm.target=v; save(); if(sm.ref&&!isC1s(byUid(sm.ref.uid))) drawChg(); else reapply(); } };
  $a(el,"[data-tg]").forEach(b=>b.onclick=()=>{ sm.target=+b.dataset.tg; save(); if(sm.ref&&!isC1s(byUid(sm.ref.uid))) drawChg(); else reapply(); });
  const xsel=$q(el,"#xcSel"); if(xsel) xsel.onclick=()=>shiftFromSel(sp);
  $q(el,"#xcShOk").onclick=()=>{ const v=+$q(el,"#xcSh").value; if(!isFinite(v)) return; setShift(sp,v,null,"직접 입력"); note(`대전 보정 Δ = ${sgn(v)} eV (직접 입력) — 시료 ${d.name} 전체`,"ok"); };
  $q(el,"#xcOff").onclick=()=>{ setShift(sp,0,null,""); note("대전 보정을 해제했습니다 — BE 가 측정값 그대로입니다"); };
}
function shiftFromSel(sp){
  const p=pr(sp); const c=p&&p.comps.find(x=>x.id===S.ui.selComp);
  if(!c){ note("먼저 기준으로 쓸 성분을 그래프나 성분 목록에서 선택하세요","warn"); return; }
  const sm=smpOf(sp); const tg=sm.target!=null?sm.target:284.8; const r=fitted(sp);
  const ref=r&&r.comps.find(x=>x.id===c.id)?{uid:sp.uid,compId:c.id,be:c.be,name:c.name||"",file:sp.file||"",target:tg}:null;
  // C 1s 의 성분을 명시적으로 고른 것이면 기준 성분으로 기억 — 'Recalibrate' 가 같은 성분을 쓰게
  if(isC1s(sp)&&ref){ sm.refUid=sp.uid; sm.refComp=c.id; }
  setShift(sp,tg-c.be,ref,`${sp.region} '${c.name||"Component"}' ${f(c.be,2)} → ${tg} eV${r?"":" (피팅 전 값)"}`);
  const w=isC1s(sp)&&r?refWarn(r.comps,c):"";
  note(`대전 보정 Δ = ${sm.shift>0?"+":""}${f(sm.shift,2)} eV — ${sm.refNote}${w}`,w?"warn":"ok");
}
/* 기준 성분이 '가장 낮은 BE 의 주 C 1s 성분'(autoRefComp)이 아니면 한 줄 경고 */
function refWarn(comps,c){ const a=autoRefComp(comps); if(!a||a.id===c.id) return "";
  return ` · 주의: 기준이 가장 낮은 BE 의 주 C 성분(‘${a.name||"Component"}’ ${f(a.be,2)} eV, 측정)이 아닙니다`; }
/* C 1s 로 맞추기: 기준 C 1s 가 피팅돼 있지 않으면 그 스펙트럼을 (화면을 옮기지 않고) 자동 제안 → 피팅한 뒤, 기준 성분을 기준 BE 로.
   기준 성분 = Reference 에서 고른 것(sm.refComp) 또는 자동(autoRefComp). 그래프에서 클릭해 둔 성분은 쓰지 않는다 —
   선택 성분을 쓰려면 Reference ▸ 'Use selected component' (명시적). */
async function calibrateC1s(sp){
  if(RUN){ note("이미 계산 중입니다","warn"); return; }
  const k=SK(sp), sm=smpOf(sp), cands=chgCands(k);
  const s=cands.find(x=>x.uid===sm.refUid)||cands[0];
  if(!s){ note("이 시료에 C 1s 스펙트럼이 없습니다 — 같은 측정의 C 1s 파일을 함께 열거나, Reference 에서 다른 기준·Δ 를 쓰세요","warn"); S.ui.chgOpen=true; drawChg(); return; }
  let didFit=false;
  if(!fitted(s)){
    setProg(`C 1s 자동 피팅 중 (${s.file||s.name})…`,0.05);
    if(!pr(s).comps.length) await autoSuggest(s);
    if(!pr(s).comps.length){ hideProg(); note("C 1s 에서 봉우리를 찾지 못해 맞추지 못했습니다 — Δ 를 직접 넣으세요","warn"); return; }
    await runFit("one",s); didFit=true;
    if(!fitted(s)){ note("C 1s 피팅이 끝나지 않아 대전 보정을 하지 못했습니다","warn"); return; }
  }
  const r=fitted(s);
  const chosen=sm.refComp?r.comps.find(x=>x.id===sm.refComp):null;
  const lostChoice=!!(sm.refComp&&!chosen);
  const c=chosen||autoRefComp(r.comps);
  if(!c){ note("C 1s 피팅에 C 성분이 없습니다 — Reference 에서 성분을 고르세요","warn"); S.ui.chgOpen=true; drawChg(); return; }
  if(lostChoice) sm.refComp=null;
  const tg=sm.target!=null?sm.target:284.8;
  const how=chosen?"고른":"가장 낮은 BE";
  const nm=c.name&&c.name!=="C 1s"?`‘${c.name}’ `:"";
  setShift(s,tg-c.be,{uid:s.uid,compId:c.id,be:c.be,name:c.name||"",file:s.file||"",target:tg},`C 1s ${nm}${f(c.be,2)} → ${f(tg,2)} eV · ${s.file||s.name} · ${how} 성분`);
  const w=chosen?refWarn(r.comps,c):"";
  const msg=`대전 보정 Δ = ${sm.shift>0?"+":""}${f(sm.shift,2)} eV — 기준 ${c.name||"C 1s 성분"} ${f(c.be,2)} eV(측정) → ${f(tg,2)} eV, 시료 ${smpDisp(k).name} 전체`
    +(didFit?` · C 1s 를 자동 제안 → 피팅 (성분 ${r.comps.length}, χ²ν ${f(r.stats.chi2nu,2)})`:"")+(lostChoice?" · 고른 기준 성분이 없어져 자동으로 골랐습니다":"")+w;
  note(msg,w||lostChoice?"warn":"ok",S.sel!==s.uid?{label:"Show C 1s fit",fn:()=>{ S.sel=s.uid; S.ui.selComp=null; save(); rerender(); }}:null);
}
/* 시료 이름 바꾸기(같은 이름이 있으면 합침) · 스펙트럼 하나를 다른 시료로 옮기기 */
function renameSample(oldK,newName){
  const nk=String(newName||"").trim()||NOSMP; if(nk===oldK){ drawList(); return; }
  const merge=!!S.specs.some(s=>SK(s)===nk);
  S.specs.forEach(s=>{ if(SK(s)===oldK) s.sample=nk===NOSMP?"":nk; });
  if(!merge&&S.samples[oldK]) S.samples[nk]=S.samples[oldK];
  delete S.samples[oldK]; AUG().entries.forEach(e=>{ if(e.sample===oldK) e.sample=nk; }); save(); rerender();
  note(merge?`시료 ‘${smpDisp(oldK).name}’ 를 ‘${smpDisp(nk).name}’ 에 합쳤습니다 — 대전 보정은 ‘${smpDisp(nk).name}’ 것을 씁니다`:`시료 이름을 ‘${smpDisp(nk).name}’ 로 바꿨습니다`,"ok");
}
function moveSpec(sp,newName){
  const nk=String(newName||"").trim()||NOSMP; if(nk===SK(sp)) return; const oldK=SK(sp);
  sp.sample=nk===NOSMP?"":nk; smpOf(sp); if(!sampleSpecs(oldK).length) delete S.samples[oldK];
  // Auger 항목: 두 스펙트럼이 모두 새 시료에 있으면 항목도 옮긴다(아니면 남겨 두고 카드에 '다른 시료' 경고)
  AUG().entries.forEach(e=>{ if(e.sample!==oldK) return; const a=byUid(e.p.uid), b=byUid(e.a.uid); if(a&&b&&SK(a)===nk&&SK(b)===nk) e.sample=nk; });
  save(); rerender(); note(`${sp.region} 을 시료 ‘${smpDisp(nk).name}’ 로 옮겼습니다`,"ok");
}

/* ---------------- 탭 ---------------- */
const TABS=[["fit","Peak Fit"],["result","Results"],["quant","Quantification"],["auger","Auger α′"],["survey","Survey"],["proc","Process"],["info","Info"]];
function drawTabs(){ const el=$q(ROOT,"#xpsTabs"); if(!el) return; el.innerHTML=TABS.map(([k,t])=>`<button data-t="${k}" class="${S.ui.tab===k?"active":""}">${t}</button>`).join(""); $a(el,"button").forEach(b=>b.onclick=()=>{ const was=S.ui.tab; S.ui.tab=b.dataset.t; save(); drawTabs(); renderPanel(); if(was==="proc"||S.ui.tab==="proc") drawPlot(); }); }
function renderPanel(){
  const el=ROOT&&$q(ROOT,"#xpsPanel"); if(!el) return; const sp=cur();
  if(!sp){ el.innerHTML=`<div class="muted sm" style="padding:8px">Select a spectrum.</div>`; return; }
  try{
    if(S.ui.tab==="fit") panelFit(el,sp);
    else if(S.ui.tab==="result") panelResult(el,sp);
    else if(S.ui.tab==="quant") panelQuant(el,sp);
    else if(S.ui.tab==="survey") panelSurvey(el,sp);
    else if(S.ui.tab==="auger") panelAuger(el,sp);
    else if(S.ui.tab==="proc") panelProc(el,sp);
    else panelInfo(el,sp);
  }catch(e){ console.error(e); el.innerHTML=`<div class="note warn">화면 오류: ${esc(e.message||e)}</div>`; }
}
function syncRegionInputs(){ const sp=cur(); if(!sp||!ROOT) return; const p=pr(sp), sh=shOf(sp); const a=$q(ROOT,"#xpsLo"), b=$q(ROOT,"#xpsHi"); if(a) a.value=f(p.lo+sh,2); if(b) b.value=f(p.hi+sh,2); }

/* ============================================================
   피팅 패널
   ============================================================ */
const CMLAB={free:"Free",fixed:"Fixed",link:"Linked"};
function panelFit(el,sp){
  const p=pr(sp), sh=shOf(sp), bg=p.bg;
  const r=p.result;
  const bgOpts=Object.entries(C.BG_TYPES).map(([k,v])=>`<option value="${k}" ${bg.type===k?"selected":""}>${esc(v)}</option>`).join("");
  el.innerHTML=`
   <div class="xps-h">Fit Range &amp; Background</div>
   <div class="xps-row"><span class="xps-f">BE <input type="number" step="0.05" id="xpsLo" value="${f(p.lo+sh,2)}"> – <input type="number" step="0.05" id="xpsHi" value="${f(p.hi+sh,2)}"> eV</span><button class="btn xs" id="xpsRgA" title="Use the whole spectrum as the fit range">Full range</button><button class="btn xs" id="xpsRgV" title="Use the range currently shown on the graph as the fit range">Visible range</button></div>
   <div class="xps-row"><select class="xps-sel" id="xpsBg">${bgOpts}</select>${bg.type!=="active"&&bg.type!=="none"&&!(bg.conv==="casa"&&(bg.type==="linear"||bg.type==="shirley"))?`<span class="xps-f" title="Number of points averaged at each endpoint">Endpoints <input type="number" class="w56" min="1" max="50" id="xpsBgN" value="${bg.n||5}"> pts</span>`:""}
     ${bg.type==="tougaard"?`<select class="xps-sel" id="xpsTgU"><option value="u2" ${!bg.u3?"selected":""}>U2 (C=1643 eV²)</option><option value="u3" ${bg.u3?"selected":""}>U3 (C, D)</option></select>${bg.u3?`<span class="xps-f">C <input type="number" class="w64" id="xpsTgC" value="${bg.C}"> D <input type="number" class="w64" id="xpsTgD" value="${bg.D}"></span>`:""}<span class="xps-f" title="Loss-onset threshold T0 (eV): cross-section = 0 for losses ≤ T0 (CasaXPS U(x; B, C, D, T0), Fairley et al. 2023). 0 = off">T0 <input type="number" class="w56" step="0.5" min="0" id="xpsTgT" value="${bg.T0||0}"></span>`:""}
     ${bg.type==="active"?`<label class="xps-f"><input type="checkbox" id="xpsBgS" ${bg.cons&&bg.cons.s&&bg.cons.s.mode!=="fixed"?"checked":""}> Fit slope too</label>`:""}
     ${bg.type==="linear"||bg.type==="shirley"?`<select class="xps-sel" id="xpsConv" title="Region convention — INNO: endpoint averages of n points at the window ends, trapezoid Shirley (SPEC §2). CasaXPS: endpoints at the 2nd / 2nd-last channel averaged over ± Av. width channels, rectangle Shirley sum — reproduces CasaXPS 2.3.27 region Area / Max / Min (SPEC §2b)"><option value="">INNO</option><option value="casa" ${bg.conv==="casa"?"selected":""}>CasaXPS</option>${bg.type==="shirley"?`<option value="avantage" ${bg.conv==="avantage"?"selected":""}>Avantage</option>`:""}</select>${bg.conv==="avantage"&&bg.type==="shirley"?`<span class="xps-f" title="Avantage 'Iterations' (1 = simple single pass; manual recommends 10). Avantage switches to Linear if the high-BE end is not higher than the low-BE end">Iter. <input type="number" class="w56" min="1" max="200" id="xpsAvIt" value="${bg.it!=null?bg.it:10}"></span><span class="xps-f" title="Avantage 'Convergence limit': stop when the relative change of the peak area between iterations is at most this (0 = run all iterations)">Conv. <input type="number" class="w64" step="any" min="0" id="xpsAvCl" value="${bg.cl!=null?bg.cl:0}"></span>`:""}${bg.conv==="casa"?`<span class="xps-f" title="CasaXPS 'Av. Width' (channels averaged on each side of the endpoint channel)">Av. width <input type="number" class="w56" min="0" max="25" id="xpsAw" value="${bg.aw!=null?bg.aw:1}"></span>`:""}`:""}
     ${bg.type==="spline"?`<span class="xps-f" title="Knots = the two window ends (endpoint averages) + your knots. A knot's height is the mean of ±2 data points unless you drag it. Natural cubic spline through the knots">Knots <b id="xpsKnN">${(bg.knots||[]).length}</b>+2</span><button class="btn xs ${S.ui.knotMode?"primary":""}" id="xpsKnA" title="Then click the graph where the background is (outside peaks). Drag a knot to move it; double-click a knot to remove it">${S.ui.knotMode?"Click the graph…":"+ Knot"}</button><button class="btn xs" id="xpsKnR" title="Return dragged knot heights to the local data mean">Auto heights</button><button class="btn xs" id="xpsKnC" title="Remove all knots (straight line between the window ends)">Clear</button>`:""}</div>
   <div class="muted xs">${bgHelp(bg.type)}</div>${casaRegionTxt(sp,p)}${bgWarn(sp,p)}
   <hr class="xps-hr">
   <div class="xps-h">Components <span class="muted sm">(${p.comps.length})</span></div>${multBanner(sp)}
   <div class="xps-row"><button class="btn sm primary" id="xpsAuto" title="Find peaks to build a starting model, then refine from the residual (add components, tune the doublet area ratio; judged by BIC). The result is only a starting point — check it, then press Fit">Auto suggest</button>${sp.reg.el==="C"&&sp.reg.orb==="1s"?`<button class="btn sm" id="xpsTplC">C 1s template</button>`:""}<button class="btn sm" id="xpsTplS" title="Save the current component set as a template">Save template</button><select class="xps-sel" id="xpsTplL"><option value="">Load template…</option>${tplList().map((t,i)=>`<option value="${i}">${esc(t.name)}</option>`).join("")}</select><button class="btn sm" id="xpsUndo" title="Undo the last change">Undo</button></div>
   <div id="xpsComps">${p.comps.map((c,ci)=>compCard(sp,p,c,ci)).join("")||`<div class="muted sm" style="padding:6px 2px">성분이 없습니다. <b>Auto suggest</b> 또는 그래프 위 <b>+ Add component</b> 로 직접 놓으세요.</div>`}</div>
   ${refHTML(sp,p)}
   <hr class="xps-hr">
   <div class="xps-h">Fit</div>
   <div class="xps-row"><button class="btn primary sm" id="xpsFit" title="Least-squares fit of the current components (LM + 2-pass IRLS weighting)">Fit</button><button class="btn sm" id="xpsMulti" title="Start from many points spread evenly within the bounds of the nonlinear parameters → keep the lowest χ²">Multi-start</button><span class="xps-f"><input type="number" class="w56" id="xpsStarts" min="2" max="400" value="${S.ui.starts||24}"> runs</span>
     <select class="xps-sel" id="xpsOpt" title="Optimizer — Marquardt: Levenberg–Marquardt (default, fast). Simplex: Nelder–Mead, derivative-free (slower; for rough χ² surfaces). Both reach the same minimum on well-posed fits"><option value="lm" ${p.optimizer!=="simplex"?"selected":""}>Marquardt</option><option value="simplex" ${p.optimizer==="simplex"?"selected":""}>Simplex</option></select>
     <select class="xps-sel" id="xpsW" title="Weighting — Poisson: variance ∝ intensity (raw counts). Constant: equal weight for every point (background-subtracted or processed data)"><option value="poisson" ${p.weight!=="const"?"selected":""}>Poisson weight</option><option value="const" ${p.weight==="const"?"selected":""}>Constant weight</option></select>
     <button class="btn sm" id="xpsMC" title="Re-add noise to the best model and refit → 95 % confidence intervals (uncertainty for reporting)">Monte Carlo uncertainty</button><span class="xps-f"><input type="number" class="w56" id="xpsMCN" min="20" max="5000" value="${S.ui.mcN||200}"> runs</span></div>
   ${r?fitStrip(r,sp):`<div class="muted xs">Weighting: ${p.weight==="const"?"constant (척도는 잔차에서 추정)":sp.c?`Poisson (${esc(sp.cSrc)})`:"quasi-Poisson (척도는 잔차에서 추정)"} · LM${p.weight==="const"?"":" + 2-pass IRLS"}</div>`}
   <div class="xps-row" style="margin-top:6px"><button class="btn xs" id="xpsSnap" ${r&&!r.stale?"":"disabled"} title="Save this fit to the model-comparison list">Save to model comparison</button><span class="muted xs">저장한 모형은 Results 탭에서 AICc·BIC 로 비교</span></div>
   ${batchPanelHTML(sp,p)}`;
  // 영역
  const lo=$q(el,"#xpsLo"), hi=$q(el,"#xpsHi");
  const setRg=()=>{ pushHist(p,"영역"); p.lo=Math.min(+lo.value,+hi.value)-sh; p.hi=Math.max(+lo.value,+hi.value)-sh; clampRegion(sp,p); touch(p); drawPlot(); drawList(); };
  lo.onchange=setRg; hi.onchange=setRg;
  $q(el,"#xpsRgA").onclick=()=>{ pushHist(p,"영역"); p.lo=sp.be[0]; p.hi=sp.be[sp.n-1]; touch(p); rerender(); };
  $q(el,"#xpsRgV").onclick=()=>{ const v=plotView(sp); pushHist(p,"영역"); p.lo=v.x0-sh; p.hi=v.x1-sh; clampRegion(sp,p); touch(p); rerender(); };
  $q(el,"#xpsBg").onchange=e=>{ pushHist(p,"배경"); bg.type=e.target.value; if(bg.type==="active"){ bg.c0=null; bg.k=null; bg.s=0; bg.cons=bg.cons||{s:{mode:"fixed"}}; } if(bg.type==="spline") bg.knots=bg.knots||[]; S.ui.knotMode=false; touch(p); rerender(); };
  const kA=$q(el,"#xpsKnA"); if(kA){ kA.onclick=()=>{ S.ui.knotMode=!S.ui.knotMode; S.ui.addMode=false; renderPanel(); drawBar(); };
    $q(el,"#xpsKnR").onclick=()=>{ pushHist(p,"매듭 높이"); (bg.knots||[]).forEach(k=>{ delete k.y; }); touch(p); drawPlot(); };
    $q(el,"#xpsKnC").onclick=()=>{ pushHist(p,"매듭 지움"); bg.knots=[]; touch(p); rerender(); }; }
  $q(el,"#xpsOpt").onchange=e=>{ pushHist(p,"최적화"); p.optimizer=e.target.value==="simplex"?"simplex":undefined; touch(p); renderPanel(); };
  bindRef(el,sp,p); bindBatchPanel(el,sp,p);
  const cv=$q(el,"#xpsConv"); if(cv) cv.onchange=()=>{ pushHist(p,"영역 규약"); if(cv.value==="casa"){ bg.conv="casa"; if(bg.aw==null) bg.aw=1; } else if(cv.value==="avantage"){ bg.conv="avantage"; if(bg.it==null) bg.it=10; } else delete bg.conv; touch(p); rerender(); };
  const ait=$q(el,"#xpsAvIt"); if(ait) ait.onchange=()=>{ bg.it=Math.max(1,Math.min(200,Math.round(+ait.value||10))); touch(p); rerender(); };
  const acl=$q(el,"#xpsAvCl"); if(acl) acl.onchange=()=>{ bg.cl=Math.max(0,+acl.value||0); touch(p); rerender(); };
  const aw=$q(el,"#xpsAw"); if(aw) aw.onchange=()=>{ bg.aw=Math.max(0,Math.min(25,Math.round(+aw.value||0))); touch(p); rerender(); };
  const bn=$q(el,"#xpsBgN"); if(bn) bn.onchange=()=>{ bg.n=Math.max(1,Math.min(50,+bn.value||5)); touch(p); drawPlot(); };
  const tu=$q(el,"#xpsTgU"); if(tu) tu.onchange=()=>{ bg.u3=tu.value==="u3"; if(bg.u3&&!bg.D){ bg.C=551; bg.D=436; } if(!bg.u3){ bg.C=1643; bg.D=0; } touch(p); rerender(); };
  const tT=$q(el,"#xpsTgT"); if(tT) tT.onchange=()=>{ const v=+tT.value; if(v>0) bg.T0=v; else delete bg.T0; touch(p); drawPlot(); };
  const tc=$q(el,"#xpsTgC"), td=$q(el,"#xpsTgD"); if(tc) tc.onchange=()=>{ bg.C=+tc.value; touch(p); drawPlot(); }; if(td) td.onchange=()=>{ bg.D=+td.value; touch(p); drawPlot(); };
  const bs=$q(el,"#xpsBgS"); if(bs) bs.onchange=()=>{ bg.cons=bg.cons||{}; bg.cons.s={mode:bs.checked?"free":"fixed"}; if(!bs.checked) bg.s=0; touch(p); };
  $q(el,"#xpsAuto").onclick=()=>autoSuggest(sp);
  $q(el,"#xpsW").onchange=e=>{ pushHist(p,"가중"); p.weight=e.target.value; touch(p); renderPanel(); drawPlot(); };
  const tpc=$q(el,"#xpsTplC"); if(tpc) tpc.onclick=templateC1s;
  $q(el,"#xpsTplS").onclick=saveTpl; $q(el,"#xpsTplL").onchange=e=>{ if(e.target.value!=="") loadTpl(+e.target.value); };
  $q(el,"#xpsUndo").onclick=undo;
  $q(el,"#xpsFit").onclick=()=>runFit("one",sp);
  $q(el,"#xpsMulti").onclick=()=>runFit("multi",sp);
  $q(el,"#xpsMC").onclick=runMC;
  $q(el,"#xpsStarts").onchange=e=>{ S.ui.starts=Math.max(2,Math.min(400,+e.target.value||24)); save(); };
  $q(el,"#xpsMCN").onchange=e=>{ S.ui.mcN=Math.max(20,Math.min(5000,+e.target.value||200)); save(); };
  $q(el,"#xpsSnap").onclick=()=>snapshotModel(sp,p);
  bindComps(el,sp,p);
}
/* 배경 경고: Tougaard 음의 세기(0 으로 둠), 정적 배경의 넓이 편향 */
function bgWarn(sp,p){
  const out=[];
  if(p.bg.type==="shirley"&&p.bg.conv==="avantage"){ const R=region(sp,p); if(R&&R.E.length>8){ const B=C.bgStatic(Float64Array.from(R.E),Float64Array.from(R.y),p.bg); if(B&&B.fallback) out.push("Avantage 규약: 높은 BE 끝이 낮은 BE 끝보다 낮아 Avantage 처럼 Linear 로 바꿨습니다."); } }
  if(p.bg.type==="shirley"||p.bg.type==="smart"){ const R=region(sp,p); if(R&&R.E.length>8){ const B=C.bgStatic(Float64Array.from(R.E),Float64Array.from(R.y),bgOf(sp,p,R)); if(B&&B.converged===false) out.push("Shirley 반복이 수렴하지 않았습니다(200 회) — 배경이 기울어 Shirley 고정점이 없는 경우입니다. ‘Shirley (in-fit)’ + ‘Fit slope too’, 또는 Linear 를 쓰세요."); } }
  if((p.bg.type==="shirley"||p.bg.type==="smart"||p.bg.type==="linear")&&p.hi-p.lo>20) out.push("창이 넓습니다("+f(p.hi-p.lo,0)+" eV): 끝점으로 정하는 배경은 창이 넓을수록 모양 오차가 쌓여 넓이가 틀립니다(합성 143 경우, Smart: 20–40 eV 최대 ≈8 %, 40 eV 넘으면 최대 ≈14 %; 20 eV 이하는 ≈3 %). ‘Shirley (in-fit)’ + ‘Fit slope too’을 권합니다.");
  if(p.bg.type==="tougaard"){ const R=region(sp,p); if(R&&R.E.length>8){ const B=C.bgStatic(Float64Array.from(R.E),Float64Array.from(R.y),p.bg); if(B&&B.negClamped) out.push("높은 BE 끝이 낮은 BE 끝보다 낮아 Tougaard 세기가 음수가 됩니다 — 비물리적이라 0 으로 두었습니다(평평한 배경). 영역을 봉우리 뒤로 30 eV 이상 넓히거나 다른 배경을 쓰세요."); } }
  if(p.bg.type!=="active"&&p.bg.type!=="none"&&p.comps.some(c=>(c.shape==="GLS"||c.shape==="GLP"||c.shape==="V")&&((c.m||0)>0.3||(c.wL||0)>0.4))) out.push("Lorentzian 꼬리가 큰 성분 + 끝점으로 정하는 배경: 창 밖·끝점 평균 안으로 들어간 꼬리 넓이가 빠져 넓이가 몇 % 작게 나옵니다(합성 시험 −2…−5 %). 창을 넓히거나 ‘Shirley (in-fit)’를 쓰세요.");
  return out.map(t=>`<div class="note warn" style="margin-top:4px">${esc(t)}</div>`).join("");
}
function bgHelp(t){
  return {linear:"끝점 평균을 잇는 직선.",
    shirley:"Proctor–Sherwood 반복 Shirley: 각 BE 의 배경이 그보다 낮은 BE 쪽 봉우리 넓이에 비례.",
    smart:"Avantage 의 Smart 와 같은 목적: Shirley 반복 중 배경이 (평활한) 데이터보다 높아지지 않게 자른다. 끝이 수평일 필요 없음.",
    tougaard:"Tougaard universal cross-section. 세기 s 는 높은 BE 끝점에서 배경 = 데이터가 되게 정함. 영역을 넓게(봉우리 뒤 30 eV 이상 권장).",
    spline:"매듭(창 양끝 평균 + 직접 놓은 점)을 지나는 자연 3차 spline. 매듭은 봉우리 밖 배경 위에 두세요.",
    active:"배경 = c₀ + k·∫(성분 합) (+ 기울기). 배경과 봉우리를 함께 피팅 — 겹친 영역에서 가장 일관된 방법.",
    none:"배경 없음 (서베이·이미 뺀 데이터)."}[t]||"";
}
function fitStrip(r,sp){
  const st=r.stats; const at=(r.cov&&r.cov.atB||[]).map((b,i)=>b?r.free[i]:null).filter(Boolean);
  // 이름은 그 피팅 당시의 성분 목록(r.comps)에서 — 피팅 뒤 성분을 바꿨으면 지금 목록에 없을 수 있다
  const nm=id=>{ if(id==="bg") return "배경"; const c=(r.comps||[]).find(c=>c.id===id)||pr(sp).comps.find(c=>c.id===id); return c?(c.name||"Component"):"(지운 성분)"; };
  const atTxt=a=>{ const c=(r.comps||[]).find(c=>c.id===a.id); let val=null; try{ val=a.id==="bg"?(r.bg||{})[a.k]:(c?C.getP(c,a.k):null); }catch(e){}
    const add=a.k==="be"?shOf(sp):0; const side=val!=null&&a.lo!=null&&a.hi!=null?(Math.abs(val-a.lo)<=Math.abs(val-a.hi)?"하한":"상한"):"";
    const bnd=side==="하한"?a.lo:side==="상한"?a.hi:null;
    return `${nm(a.id)} ${(C.PDEF[a.k]||{lab:a.k}).lab}${val!=null?` = ${g(val+add,4)}`:""}${side?` (${side} ${g(bnd+add,4)})`:""}`; };
  const parts=[
    `<div class="${st.chi2nu>2?"warn":""}">χ²ν<b>${f(st.chi2nu,3)}</b></div>`,
    `<div>R²<b>${st.R2!=null?f(st.R2,5):"—"}</b></div>`,
    `<div class="${Math.abs(st.rho1)>0.3?"warn":""}" title="Neighbour correlation of the normalized residuals">Residual ρ₁<b>${f(st.rho1,2)}</b></div>`,
    `<div class="${Math.abs(st.runsZ)>3?"warn":""}" title="Runs-test z of the residual signs (|z|>3 = structured residual)">runs z<b>${f(st.runsZ,1)}</b></div>`,
    `<div>AICc<b>${f(st.AICc,1)}</b></div>`,
    `<div class="${r.conv?"ok":"warn"}" title="${(r.iter0||0)+r.iter?"LM iterations (start-point fit + final refinement)":"The starting values were already at the minimum (e.g. after Auto suggest) — nothing to change"}">${r.conv?"Converged":"Not converged"}<b>${(r.iter0||0)+r.iter?((r.iter0||0)+r.iter)+" iter":"at start"}</b></div>`];
  if(r.multi) parts.push(`<div class="${r.multi.frac<0.2?"warn":""}" title="Fraction of starts reaching the minimum (Δχ²<1)">Same minimum<b>${r.multi.near}/${r.multi.N}</b></div>`);
  if(r.lcf) parts.push(`<div class="${r.lcf.warn?"warn":""}" title="Collinearity of the reference shapes on this data: scaled condition number κ (Belsley) and VIF per reference (${r.lcf.vif.map(v=>f(v,1)).join(", ")}). κ > 30 or VIF > 10 = references too alike — fractions poorly determined">Collinearity κ<b>${f(r.lcf.cond,1)}</b></div>`);
  if(r.optimizer==="simplex") parts.push(`<div title="Nelder–Mead simplex; count = χ² evaluations">Simplex<b>${r.iter} eval.</b></div>`);
  let w="";
  // 피팅 뒤 모형을 바꿨으면(stale) 위 숫자와 경고는 이전 모형의 것 — 경고는 숨기고 다시 피팅하라고만 한다
  if(r.stale){ w+=`<div class="note warn" style="margin-top:4px">피팅 뒤 성분·영역을 바꿨습니다 — 위 숫자는 <b>이전 모형</b>의 값입니다. <b>Fit</b> 을 다시 누르세요.</div>`; return `<div class="strip stale">${parts.join("")}</div>${w}`; }
  if(at.length) w+=`<div class="note warn" style="margin-top:4px"><b>한계에 닿은 매개변수</b>: ${at.map(a=>esc(atTxt(a))).join(", ")}<br>이 값은 데이터가 아니라 한계가 정한 값이라 σ 를 계산하지 않습니다. 넓이가 하한 0 이면 그 성분을 빼고 다시 피팅하세요. Δ·FWHM·BE 등이 물리적 한계에 붙었다면 모형(성분 수·배경·모양)이 데이터와 맞지 않는다는 신호입니다.</div>`;
  if(r.lcf&&r.lcf.warn) w+=`<div class="note warn" style="margin-top:4px">기준 스펙트럼끼리 너무 닮았습니다(κ ${f(r.lcf.cond,0)}) — 분율의 불확도가 큽니다. 기준을 줄이거나 Monte Carlo 로 구간을 확인하세요.</div>`;
  if(st.chi2nu>3) w+=`<div class="note warn" style="margin-top:4px">χ²ν 가 큽니다 — 성분이 부족하거나 배경·모양이 맞지 않을 수 있습니다(잔차 그림 확인).</div>`;
  if(r.cov&&r.cov.singular) w+=`<div class="note warn" style="margin-top:4px">공분산 행렬이 특이에 가깝습니다 — 서로 구별되지 않는 매개변수가 있습니다(Results 탭의 Parameter Correlation Matrix).</div>`;
  return `<div class="strip">${parts.join("")}</div>${w}`;
}
function compCard(sp,p,c,ci){
  const sel=S.ui.selComp===c.id, sh=shOf(sp); const r=p.result&&!p.result.stale?p.result:null;
  const sm=r&&r.summ&&r.summ[ci]?r.summ[ci]:null;
  const atb=new Set(); if(r&&r.cov) r.cov.atB.forEach((b,i)=>{ if(b&&r.free[i]&&r.free[i].id===c.id) atb.add(r.free[i].k); });
  const names=C.paramNames(c);
  const field=k=>{ const cs=(p.cons[c.id]||{})[k]||{mode:"free"}; let v=C.getP(c,k); if(k==="be") v+=sh; const d=k==="area"?0:k==="be"||k==="dE"?3:3;
    const sig=sm&&sm.sig?(k==="be"?sm.sig.be:k==="area"?sm.sig.area:null):null;
    const lab=c.shape==="MEAS"&&k==="wG"?"Broadening":c.shape==="MEAS"&&k==="be"?"BE (max)":c.shape==="LA"&&k==="alpha"?"α (hi BE)":c.shape==="LA"&&k==="beta"?"β (lo BE)":c.shape==="DS"&&k==="alpha"?"α (asym.)":(C.PDEF[k]||{lab:k}).lab;
    const tip=c.shape==="MEAS"&&k==="wG"?"Extra Gaussian broadening (FWHM, eV) applied to the reference shape. Fixed = off":c.shape==="MEAS"&&k==="be"?"Position of the reference maximum. Fixed = LCF (no shift); Free = shift δ fitted":c.shape==="LA"&&(k==="alpha"||k==="beta")?"Lorentzian power-law exponent — smaller = longer tail on that side (α: high-BE side, β: low-BE side; metals: α < β)":((C.PDEF[k]||{}).lab||k)+" "+((C.PDEF[k]||{}).u||"");
    return `<div class="pf ${atb.has(k)?"atb":""}"><span class="pl" title="${esc(tip)}">${esc(lab)}</span><input type="number" step="${k==="area"?"any":k==="be"?"0.01":"0.01"}" data-p="${k}" value="${k==="area"?Math.round(v):(+v).toFixed(d)}" ${cs.mode==="link"?"disabled":""}><button class="cm ${cs.mode}" data-cm="${k}" title="Free → Fixed → Linked">${CMLAB[cs.mode]}</button></div>`; };
  const open=sel||S.ui.openAll;
  const head=`<div class="ch" data-sel="${c.id}"><span class="sw" style="background:${colorOf(c)}"></span><input class="nmi" data-nm="${c.id}" value="${esc(c.name||"")}" placeholder="Name"><span class="sum">${c.shape==="MEAS"?"Ref":esc(c.shape)}${c.doublet?" ·D":""} · ${f(c.be+sh,2)} eV${sm&&sm.fwhm?` · ${f(sm.fwhm,2)} eV`:""}</span><button class="btn xs" data-dl="${c.id}" title="Delete">×</button></div>`;
  if(!open) return `<div class="comp" data-c="${c.id}">${head}</div>`;
  const ced=S.ui.ced&&S.ui.ced.id===c.id?consEditor(p,c,S.ui.ced.k):"";
  const isM=c.shape==="MEAS", tp=c.tpl||{};
  const shapeSel=isM?`<span class="badge inf" title="${esc(`Measured line shape from ${tp.name||"a reference spectrum"} — background ${tp.bgType||""}, area-normalized, ${tp.interp==="linear"?"linear":"PCHIP"} interpolation. Reference position ${f(tp.eref,2)} eV (charge-corrected)`)}">Measured · ${esc(tp.name||"reference")}</span>`
    :`<select class="xps-sel" data-shape="${c.id}">${Object.entries(C.SHAPES).filter(([k])=>k!=="MEAS").map(([k,v])=>`<option value="${k}" ${c.shape===k?"selected":""}>${esc(v.name)}</option>`).join("")}</select>`;
  const so=soInfo(c.el||sp.reg.el,c.line||sp.reg.orb);
  return `<div class="comp ${sel?"sel":""}" data-c="${c.id}">${head}<div class="cb">
    <div class="xps-row">${shapeSel}${isM?"":`<label class="xps-f"><input type="checkbox" data-dbl="${c.id}" ${c.doublet?"checked":""}> Doublet</label>`}
    <span class="xps-f" title="Element this component belongs to (for quantification)">Element <input type="text" class="xps-in w56" data-el="${c.id}" value="${esc(c.el||"")}"></span><label class="xps-f" title="Include in the quantification sum"><input type="checkbox" data-q="${c.id}" ${c.noQuant?"":"checked"}> Quantify</label></div>
    <div class="pg">${names.map(field).join("")}${ced}</div>
    ${c.doublet&&so?`<div class="muted xs" style="margin-top:4px">${esc(soNote(so))}. 화학 상태에 따라 Δ 는 조금 다를 수 있음</div>`:""}
    ${casaEqTxt(c,sm)}
    ${isM?`<div class="muted xs" style="margin-top:4px">Shift δ = ${f(c.be+sh-tp.eref,3)} eV</div>`:""}
    ${sm?`<div class="muted xs" style="margin-top:4px">True FWHM ${f(sm.fwhm,3)}${sm.sig.fwhm!=null?" ± "+f(sm.sig.fwhm,3):""} eV · height ${g(sm.height,4)} · area ${f(sm.areaPct,1)} %${sm.sig.areaPct!=null?" ± "+f(sm.sig.areaPct,1):""}</div>`:""}
  </div></div>`;
}
function consEditor(p,c,k){
  const cs=(p.cons[c.id]||{})[k]||{mode:"free"}; const sh=shOf(cur());
  if(cs.mode==="link"){ const others=p.comps.filter(o=>o.id!==c.id&&(k!=="dE"&&k!=="r"&&k!=="q"||o.doublet)&&C.paramNames(o).includes(k));
    const isOff=k==="be";
    return `<div class="ced" data-ced="${k}"><b>${esc((C.PDEF[k]||{lab:k}).lab)}</b> = <select class="xps-sel" data-lr>${others.map(o=>`<option value="${o.id}" ${cs.ref===o.id?"selected":""}>${esc(o.name||o.id)}</option>`).join("")}</select> ${isOff?"+":"×"} <input type="number" class="w64" step="any" data-lk value="${isOff?(cs.b!=null?cs.b:0):(cs.k!=null?cs.k:1)}"> ${isOff?"eV":""}<button class="btn xs" data-cok>Apply</button></div>`; }
  if(cs.mode==="free"){ const [dl,dh]=C.defBounds(c,k); const lo=cs.lo!=null?cs.lo:dl, hi=cs.hi!=null?cs.hi:dh; const add=k==="be"?sh:0;
    return `<div class="ced" data-ced="${k}"><b>${esc((C.PDEF[k]||{lab:k}).lab)}</b> bounds <input type="number" class="w64" step="any" data-blo value="${isFinite(lo)?+(lo+add).toFixed(4):""}" placeholder="−∞"> – <input type="number" class="w64" step="any" data-bhi value="${isFinite(hi)?+(hi+add).toFixed(4):""}" placeholder="∞"> <button class="btn xs" data-cok>Apply</button><span class="muted xs">${k==="be"?"비우면 지금 값 ±3 eV":""}</span></div>`; }
  return "";
}
function bindComps(el,sp,p){
  const sh=shOf(sp);
  $a(el,"[data-sel]").forEach(h=>h.onclick=e=>{ if(e.target.closest("button,input,select")) return; const id=h.dataset.sel; S.ui.selComp=S.ui.selComp===id?null:id; S.ui.ced=null; renderPanel(); drawPlot(); });
  $a(el,"[data-nm]").forEach(i=>i.onchange=()=>{ const c=p.comps.find(c=>c.id===i.dataset.nm); c.name=i.value; touch(p); drawPlot(); drawList(); });
  $a(el,"[data-dl]").forEach(b=>b.onclick=()=>delComp(b.dataset.dl));
  $a(el,".comp[data-c]").forEach(card=>{
    const c=p.comps.find(x=>x.id===card.dataset.c); if(!c) return;
    $a(card,"input[data-p]").forEach(inp=>inp.onchange=()=>{ const k=inp.dataset.p; let v=+inp.value; if(!isFinite(v)) return; if(k==="be") v-=sh; pushHist(p,"값 수정"); C.setP(c,k,v); propagate(p); touch(p); drawPlot(); });
    $a(card,"[data-cm]").forEach(b=>b.onclick=()=>{ const k=b.dataset.cm; const cs=(p.cons[c.id]=p.cons[c.id]||{}); const curm=(cs[k]||{mode:"free"}).mode;
      const order=["free","fixed","link"]; let nx=order[(order.indexOf(curm)+1)%3];
      if(nx==="link"&&!p.comps.some(o=>o.id!==c.id&&C.paramNames(o).includes(k))) nx="free";
      const prev=cs[k]||{}; cs[k]=Object.assign({},prev,{mode:nx});
      if(nx==="link"){ const ref=p.comps.find(o=>o.id!==c.id&&C.paramNames(o).includes(k)); cs[k].ref=prev.ref||ref.id; if(k==="be"){ cs[k].k=1; cs[k].b=+(c.be-(p.comps.find(o=>o.id===cs[k].ref)||ref).be).toFixed(3); } else { const rv=C.getP(p.comps.find(o=>o.id===cs[k].ref)||ref,k); cs[k].k=rv?+(C.getP(c,k)/rv).toPrecision(4):1; cs[k].b=0; } }
      S.ui.ced=nx==="fixed"?null:{id:c.id,k}; propagate(p); touch(p); renderPanel(); drawPlot(); });
    const shp=$q(card,"[data-shape]"); if(shp) shp.onchange=()=>{ pushHist(p,"모양"); const i=p.comps.indexOf(c); const n=C.convertShape(c,shp.value); if(n.shape==="AGT") n.ch=chOf(sp); n.slot=c.slot; n.noQuant=c.noQuant; p.comps[i]=n; const cs=p.cons[c.id]||{}; Object.keys(cs).forEach(k=>{ if(!C.paramNames(n).includes(k)) delete cs[k]; }); touch(p); renderPanel(); drawPlot(); };
    const db=$q(card,"[data-dbl]"); if(db) db.onchange=()=>{ pushHist(p,"doublet"); if(db.checked){ const so=soInfo(c.el||sp.reg.el,c.line||sp.reg.orb)||{dE:null,r:0.5}; if(so.multiplet) note("다중항 준위에 doublet 을 켰습니다 — 껍질이 닫힌 상태("+(so.multiplet.closed||"해당 없음")+")가 아니면 doublet 모형은 맞지 않습니다","warn"); c.doublet={dE:so.dE||1,r:so.r||0.5,q:1}; c.area*=1+c.doublet.r; const cs=(p.cons[c.id]=p.cons[c.id]||{}); cs.dE=C.dEConstraint(so.dE?so:null); cs.r={mode:"fixed"}; cs.q={mode:"fixed"}; if(!so.dE) note("이 궤도의 Δ 를 라이브러리에서 찾지 못했습니다 — 직접 넣으세요","warn"); } else { c.area/=1+c.doublet.r; c.doublet=null; const cs=p.cons[c.id]||{}; delete cs.dE; delete cs.r; delete cs.q; } touch(p); renderPanel(); drawPlot(); };
    const ei=$q(card,"[data-el]"); if(ei) ei.onchange=()=>{ c.el=ei.value.trim(); touch(p); };
    const qi=$q(card,"[data-q]"); if(qi) qi.onchange=()=>{ c.noQuant=!qi.checked; touch(p); };
    const ced=$q(card,".ced"); if(ced){ const k=ced.dataset.ced; const ok=$q(ced,"[data-cok]"); ok.onclick=()=>{ const cs=(p.cons[c.id]=p.cons[c.id]||{}); const o=cs[k]||{mode:"free"};
      if(o.mode==="link"){ o.ref=$q(ced,"[data-lr]").value; const val=+$q(ced,"[data-lk]").value; if(k==="be"){ o.k=1; o.b=isFinite(val)?val:0; } else { o.k=isFinite(val)?val:1; o.b=0; } }
      else { const a=$q(ced,"[data-blo]").value, b=$q(ced,"[data-bhi]").value; const add=k==="be"?sh:0; o.lo=a===""?null:+a-add; o.hi=b===""?null:+b-add; }
      cs[k]=o; propagate(p); touch(p); S.ui.ced=null; renderPanel(); drawPlot(); }; }
  });
}
function propagate(p){ try{ const prob={comps:p.comps,cons:p.cons,bg:p.bg}; const map=C.buildMap(prob); C.applyLinks(prob,map); }catch(e){ note(e.message,"warn"); } }
function snapshotModel(sp,p){ const r=p.result; if(!r||r.stale) return; p.models=p.models||[]; const nm=`Model ${p.models.length+1} (${p.comps.length} comp., ${C.BG_TYPES[p.bg.type]})`;
  p.models.push({name:nm,when:new Date().toISOString(),lo:p.lo,hi:p.hi,bg:C.clone(p.bg),comps:C.clone(p.comps),cons:C.clone(p.cons),stats:{N:r.stats.N,p:r.stats.p,chi2:r.stats.chi2,chi2nu:r.stats.chi2nu,AICc:r.stats.AICc,BIC:r.stats.BIC,rho1:r.stats.rho1},c:sp.c||null});
  save(); note("모형 비교 목록에 저장: "+nm,"ok"); }
/* 템플릿 (브라우저에 저장) */
const TKEY="innoanalytics_xps_tpl_v1";
function tplList(){ try{ return JSON.parse(localStorage.getItem(TKEY)||"[]"); }catch(e){ return []; } }
/* 템플릿은 대전 보정된 BE 로 저장(bev:2): 성분 BE 와 'free' BE 한계를 함께 +Δ, 불러올 때 함께 −Δ (C.shiftModelBE) → 저장·불러오기는 어떤 Δ 에서도 항등.
   bev 가 없는 옛 템플릿은 한계만 측정 BE 로 저장돼 있어(결함) 한계를 불러온 성분 중심으로 다시 둔다(폭 유지). */
function saveTpl(){ const sp=cur(), p=pr(sp); if(!p.comps.length){ note("저장할 성분이 없습니다","warn"); return; } const L=tplList(); const nm=`${sp.region} · ${p.comps.length} comp. · ${new Date().toISOString().slice(0,10)}`; const sh=shOf(sp);
  const m=C.shiftModelBE(p.comps,p.cons,+sh);
  const bgc=C.clone(p.bg); if(bgc.knots) bgc.knots.forEach(k=>{ k.e+=sh; });
  L.push({name:nm,region:sp.region,bev:2,comps:m.comps,cons:m.cons,bg:bgc,weight:p.weight||"poisson",optimizer:p.optimizer}); try{ localStorage.setItem(TKEY,JSON.stringify(L.slice(-40))); }catch(e){} renderPanel(); note("템플릿 저장: "+nm+" (BE 는 대전 보정 값으로 저장)","ok"); }
function tplApply(t,sh){
  if(t.bev===2) return C.shiftModelBE(t.comps,t.cons,-sh);
  const m=C.shiftModelBE(t.comps,{},-sh); m.cons=C.clone(t.cons||{}); let moved=0;
  Object.entries(m.cons).forEach(([id,o])=>{ const b=o&&o.be, c=m.comps.find(x=>x.id===id); if(!b||!c||b.mode==="link"||b.lo==null||b.hi==null||!isFinite(b.lo)||!isFinite(b.hi)) return;
    const w0=c.be-(b.lo+b.hi)/2; if(Math.abs(w0)<1e-9) return; b.lo+=w0; b.hi+=w0; moved++; });
  m.legacy=moved; return m; }
function loadTpl(i){ const t=tplList()[i]; if(!t) return; const sp=cur(), p=pr(sp), sh=shOf(sp); pushHist(p,"템플릿");
  const m=tplApply(t,sh); p.comps=m.comps; p.cons=m.cons; p.bg=Object.assign(p.bg,{type:t.bg.type,n:t.bg.n,C:t.bg.C,D:t.bg.D,u3:t.bg.u3}); if(t.weight) p.weight=t.weight;
  if(t.bg.knots) p.bg.knots=t.bg.knots.map(k=>Object.assign({},k,{e:k.e-(t.bev===2?sh:0)})); if(t.optimizer) p.optimizer=t.optimizer;
  touch(p); rerender(); note("템플릿을 불러왔습니다 — 넓이는 피팅 때 다시 맞춰집니다"+(m.legacy?` (옛 형식: BE 범위 ${m.legacy}개를 성분 중심으로 다시 맞춤)`:"")); }

/* ============================================================
   결과 패널
   ============================================================ */
function panelResult(el,sp){
  const p=pr(sp), r=p.result, sh=shOf(sp);
  if(!r){ el.innerHTML=`<div class="muted sm" style="padding:6px">No fit yet. <span class="xs">Peak Fit 탭에서 Fit 을 누르세요.</span></div>${modelsHTML(p)}`; bindModels(el,sp,p); return; }
  const mc=r.mc&&r.mc.summ;
  const rows=[];
  r.comps.forEach((c,i)=>{ const s=r.summ[i]||{sig:{}}; const m=mc?mc[i]:null; const col=colorOf(p.comps[i]||c);
    const lg=c.shape==="MEAS"?`${esc((c.tpl&&c.tpl.name)||"")} · δ ${f(c.be+sh-(c.tpl?c.tpl.eref:0),2)} eV${c.wG>1e-6?` · w ${f(c.wG,2)}`:""}`:c.shape==="GLS"||c.shape==="GLP"||c.shape==="GLT"?`${Math.round(100*c.m)} %`:c.shape==="V"?`G ${f(c.wG,2)}/L ${f(c.wL,2)}`:c.shape==="DS"?`α ${f(c.alpha,3)}`:c.shape==="LA"?`α ${f(c.alpha,2)} β ${f(c.beta,2)}`:"";
    const nm1=c.name||("#"+(i+1));
    const be=s.be+sh, sb=m?m.be.sd:s.sig.be; const fx=k=>s.sig&&s.sig[k]===null;
    rows.push(`<tr><td><i class="sw" style="background:${col}"></i>${esc(nm1)}${c.doublet?` <span class="muted xs">(pair sum)</span>`:""}</td><td class="num">${pmF(be,sb,2,fx("be"))}</td><td class="num">${pmF(s.fwhm,m?m.fwhm.sd:s.sig.fwhm,2,fx("fwhm"))}</td><td class="num">${g(s.height,4)}</td><td class="num">${pmF(s.area,m?m.area.sd:s.sig.area,0,fx("area"))}</td><td class="num">${pm(s.areaPct,m?m.areaPct.sd:s.sig.areaPct,1)}</td><td>${esc(c.shape)} ${lg}</td></tr>`);
    if(c.doublet){ const j=(c.line||sp.reg.orb||"").slice(-1); const J=ORB_J[j]||["",""]; const base=nm1.replace(/(1\/2|3\/2|5\/2|7\/2)$/,"");
      rows.push(`<tr class="sub"><td style="padding-left:18px">${esc(base.includes(J[0])?base:base+J[0])}</td><td class="num">${f(be,2)}</td><td class="num">${f(s.fwhm,2)}</td><td class="num"></td><td class="num">${f(s.area1,0)}</td><td></td><td class="muted xs">major</td></tr>`);
      rows.push(`<tr class="sub"><td style="padding-left:18px">${esc((base.includes(J[0])?base.replace(J[0],J[1]):base+J[1]))}</td><td class="num">${f(s.be2+sh,2)}</td><td class="num">${f(s.fwhm*c.doublet.q,2)}</td><td class="num"></td><td class="num">${f(s.area2,0)}</td><td></td><td class="muted xs">Δ ${f(c.doublet.dE,2)}, ratio ${f(c.doublet.r,3)}</td></tr>`); }
  });
  const tot=r.summ.reduce((s,q)=>s+q.area,0);
  rows.push(`<tr class="tot"><td>Total</td><td></td><td></td><td></td><td class="num">${pm(tot,r.mc&&r.mc.areaSum?r.mc.areaSum.sd:r.sumSig&&r.sumSig.all,0)}</td><td class="num">100</td><td></td></tr>`);
  const st=r.stats;
  el.innerHTML=`${chgBanner(sp)}${multBanner(sp)}
   <div class="xps-h">Peak Table ${r.stale?`<span class="badge wr">Stale fit</span>`:""} <span class="badge ${mc?"ok":"inf"}" ${r.comps.some(c=>c.shape==="MEAS")?`title="${mc?"Monte Carlo includes the noise of the reference spectra (templates rebuilt from noisy copies)":"Covariance σ does not include reference-spectrum noise — run Monte Carlo"}"`:""}>${mc?`± = Monte Carlo SD (${r.mc.n} runs${r.comps.some(c=>c.shape==="MEAS")?", incl. ref. noise":""})`:"± = covariance 1σ"}</span></div>
   <div class="tw"><table class="xps-tbl"><thead><tr><th>Name</th><th class="num">BE (eV)</th><th class="num">FWHM</th><th class="num">Height (CPS)</th><th class="num">Area (CPS·eV)</th><th class="num"${r.comps.some(c=>c.shape==="MEAS")?` title="For reference (LCF) components: fraction of the reference sum"`:""}>Area %${r.comps.some(c=>c.shape==="MEAS")?" (fraction)":""}</th><th>Shape</th></tr></thead><tbody>${rows.join("")}</tbody></table></div>
   <div class="muted xs" style="margin-top:4px">BE 는 대전 보정 ${sh?`(Δ ${sh>0?"+":""}${f(sh,2)} eV)`:"없음"} 값. FWHM 은 모양 함수의 실제 반치폭. 넓이 = 성분 함수의 ℝ 적분(DS 는 영역 안 적분).${(r.stats.ar1Inflate>1.15)?` 잔차 이웃상관 ρ₁=${f(st.rho1,2)} → 공분산 σ 는 약 ×${f(st.ar1Inflate,2)} 과소평가일 수 있음(Monte Carlo 는 이를 반영).`:""}</div>
   ${mc?mcHTML(r,sp):""}
   <div class="xps-h">Goodness of Fit</div>
   <div class="strip">
     <div>N · p<b>${st.N} · ${st.p}</b></div><div>χ²<b>${g(st.chi2,5)}</b></div><div>χ²ν (φ)<b>${f(st.chi2nu,3)}</b></div><div>R²<b>${st.R2!=null?f(st.R2,5):"—"}</b></div>
     <div>Residual RMS/σ<b>${f(st.rmsNorm,3)}</b></div><div>Durbin–Watson<b>${f(st.DW,2)}</b></div><div>ρ₁<b>${f(st.rho1,3)}</b></div><div>runs z<b>${f(st.runsZ,2)}</b></div>
     <div>AICc<b>${f(st.AICc,2)}</b></div><div>BIC<b>${f(st.BIC,2)}</b></div><div>Iterations · time<b>${(r.iter0||0)+r.iter} · ${r.ms} ms</b></div>${r.multi?`<div>Multi-start<b>${r.multi.near}/${r.multi.N}</b></div>`:""}
   </div>
   <div class="muted xs">Weighting: ${p.weight==="const"?"constant (모든 점 같은 무게, 척도 φ 는 잔차에서 추정)":sp.c?`Poisson, 1 CPS = ${g(sp.c,4)} count (${esc(sp.cSrc)})`:"quasi-Poisson (척도 φ 는 잔차에서 추정)"}. ${sp.c&&p.weight!=="const"?`φ≈1 이면 잡음이 Poisson 과 일치. 지금 φ=${f(st.chi2nu,2)}.`:""}</div>
   ${sensHTML(sp,p)}
   ${corrHTML(r,p)}
   ${modelsHTML(p)}
   ${batchHTML(sp)}
   <div class="xps-h">Export</div>
   <div class="xps-row"><button class="btn sm" id="xpsCsv">Curves CSV</button><button class="btn sm" id="xpsPkCsv">Peak table CSV</button><button class="btn sm" id="xpsSvgE">Figure SVG</button><button class="btn sm" id="xpsPngE">Figure PNG</button></div>`;
  $q(el,"#xpsCsv").onclick=()=>exportCurves(sp); $q(el,"#xpsPkCsv").onclick=()=>exportPeaks(sp); $q(el,"#xpsSvgE").onclick=()=>exportFig("svg"); $q(el,"#xpsPngE").onclick=()=>exportFig("png");
  bindSens(el,sp,p); bindModels(el,sp,p); bindBatch(el,sp);
}
function mcHTML(r,sp){ const sh=shOf(sp); const mc=r.mc.summ; return `<details class="xps-d" open><summary>Monte Carlo 95 % interval (${r.mc.n} runs, AR(1) ρ=${f(r.mc.rho,2)}, φ=${f(r.mc.phi,2)})</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th>Name</th><th class="num">BE 2.5–97.5 %</th><th class="num">FWHM</th><th class="num">Area %</th></tr></thead><tbody>${r.comps.map((c,i)=>{ const m=mc[i]; if(!m) return ""; return `<tr><td>${esc(c.name||"#"+(i+1))}</td><td class="num">${f(m.be.lo+sh,2)} – ${f(m.be.hi+sh,2)}</td><td class="num">${m.fwhm?f(m.fwhm.lo,2)+" – "+f(m.fwhm.hi,2):"—"}</td><td class="num">${f(m.areaPct.lo,1)} – ${f(m.areaPct.hi,1)}</td></tr>`; }).join("")}</tbody></table>${r.mc.fails?`<div class="muted xs">재피팅 실패 ${r.mc.fails}회는 빼고 계산</div>`:""}</div></details>`; }
function corrHTML(r,p){
  const c=r.cov; if(!c||!c.np) return ""; const np=c.np; const lab=i=>{ const fr=r.free[i]; const comp=r.comps.find(q=>q.id===fr.id); return (fr.id==="bg"?"bg":(comp?comp.name||comp.id:fr.id))+"·"+((C.PDEF[fr.k]||{lab:fr.k}).lab); };
  const hi=[]; for(let a=0;a<np;a++) for(let b=0;b<a;b++){ const v=c.corr[a*np+b]; if(isFinite(v)&&Math.abs(v)>0.9) hi.push({a,b,v}); }
  const cw=Math.max(10,Math.min(22,300/np)), lw=150, W=lw+np*cw+4, H=np*cw+4;
  const cell=v=>{ if(!isFinite(v)) return "var(--chip)"; const t=Math.min(1,Math.abs(v)); return v>0?`color-mix(in srgb,var(--s1) ${Math.round(t*100)}%,var(--surface))`:`color-mix(in srgb,var(--s2) ${Math.round(t*100)}%,var(--surface))`; };
  let s=`<svg class="heat" viewBox="0 0 ${W} ${H}" width="${W}">`;
  for(let a=0;a<np;a++){ s+=`<text x="${lw-4}" y="${a*cw+cw*0.7+2}" text-anchor="end" style="font-size:${Math.min(10.5,cw*0.6)}px;fill:var(--ink2)">${esc(lab(a))}</text>`; for(let b=0;b<np;b++){ const v=c.corr[a*np+b]; s+=`<rect x="${lw+b*cw}" y="${a*cw+2}" width="${cw-1}" height="${cw-1}" fill="${cell(v)}"><title>${esc(lab(a))} ↔ ${esc(lab(b))}: ${isFinite(v)?v.toFixed(3):"—"}</title></rect>`; } }
  s+=`</svg>`;
  return `<details class="xps-d" ${hi.length?"open":""}><summary>Parameter Correlation Matrix ${hi.length?`<span class="badge wr">|ρ|>0.9: ${hi.length} pairs</span>`:""}</summary><div class="b"><div class="tw">${s}</div><div class="muted xs" style="margin-top:4px">파랑 = 양의 상관, 주황 = 음의 상관. |ρ|>0.9 인 쌍은 데이터만으로 따로 정해지지 않는다는 뜻 — 연결·고정으로 줄이세요.</div>${hi.slice(0,8).map(h=>`<div class="sm">${esc(lab(h.a))} ↔ ${esc(lab(h.b))}: <b>${f(h.v,3)}</b></div>`).join("")}</div></details>`;
}
/* ---------------- 모형 민감도 (결과 탭) ----------------
   지금 성분 구성을 그대로 두고 분석 선택(배경 종류·끝점 점 수·창 끝)을 바꿔 다시 맞춘다 → 각 값의 흔들림.
   보고 구간 = [변형들의 최소 − 2σ통계, 최대 + 2σ통계] (합성 시험에서 기준 배경이 틀렸을 때도 참값을 거의 늘 포함 — SPEC §11). */
const SENS_FAC={bg:"Background type",n:"Endpoint points",lo:"Low-BE window edge",hi:"High-BE window edge",kn:"Spline knot positions"};
function sensOpts(){ const o=S.ui.sens||(S.ui.sens={mode:"oat",dw:0.5,bgs:{smart:true,shirley:true,active:true,tougaard:false},ns:true}); return o; }
function sensVarList(sp,p){ const o=sensOpts(); const wide=p.hi-p.lo>=30;
  const bgs=Object.keys(o.bgs).filter(k=>o.bgs[k]&&(k!=="tougaard"||wide));
  return C.sensVariants({lo:p.lo,hi:p.hi,bg:p.bg},{be:sp.be},{mode:o.mode,dw:+o.dw||0.5,dk:o.dk!=null?+o.dk:0.5,bgs,ns:o.ns?[3,10]:[]}); }
async function runSens(sp,p){
  if(RUN){ note("이미 계산 중입니다","warn"); return; }
  const r=p.result; if(!r||r.stale){ note("먼저 피팅하세요 (수정 뒤에는 다시 피팅)","warn"); return; }
  const variants=sensVarList(sp,p); const model={lo:p.lo,hi:p.hi,bg:C.clone(p.bg),comps:C.clone(r.comps),cons:C.clone(p.cons),weight:p.weight,optimizer:p.optimizer};
  const spec={be:sp.be,y:sp.y,c:sp.c||null};
  const K=Math.min(nWorkers(),variants.length); const chunks=[...Array(K)].map(()=>[]); variants.forEach((v,i)=>chunks[i%K].push(Object.assign({},v,{idx:i})));
  let dn=0; setProg(`모형 민감도 0/${K} 묶음 (변형 ${variants.length}개)`,0.02); const t0=performance.now(); const ver=p.ver;
  try{
    const outs=await pool(chunks.map(ch=>({cmd:"sens",spec,model,variants:ch})),()=>{ dn++; setProg(`모형 민감도 ${dn}/${K} 묶음`,dn/K); });
    const got=[].concat(...outs.filter(Boolean)); const res=variants.map((v,i)=>got.find(x=>x.i===i)||{i,ok:false,err:"계산 묶음 실패"}); const sum=C.sensSummary(variants,res); hideProg();
    if(p.ver!==ver){ note("계산 중에 모형이 바뀌어 결과를 버렸습니다 — 다시 실행하세요","warn"); return; }
    if(!sum){ note("모든 변형이 실패했습니다 — 창·배경 설정을 확인하세요","warn"); return; }
    p.sens={ver,when:new Date().toISOString(),ms:Math.round(performance.now()-t0),opts:C.clone(sensOpts()),variants:variants.map(v=>({lo:v.lo,hi:v.hi,f:v.f})),res,sum};
    save(); renderPanel(); note(`모형 민감도: 변형 ${sum.nOk}/${sum.n}개를 다시 맞췄습니다`,"ok");
  }catch(e){ hideProg(); RUN=null; note("모형 민감도 중단: "+(e.message||e),"warn"); }
}
/* 통계 σ(넓이 %): 공분산 σ × AR(1) 팽창(잔차 이웃상관), Monte Carlo 가 있으면 그 sd 와 큰 쪽 */
function sigStat(r,ci){ const q=r.summ&&r.summ[ci]; const sc=q&&q.sig&&q.sig.areaPct!=null?q.sig.areaPct*((r.stats&&r.stats.ar1Inflate)||1):0;
  const m=r.mc&&r.mc.summ&&r.mc.summ[ci]&&r.mc.summ[ci].areaPct; return Math.max(sc,m&&m.sd!=null?m.sd:0); }
function sensStrip(vals,base,lo,hi){ const W=150,H=16; if(!(hi>lo)){ lo-=0.5; hi+=0.5; } const X=v=>4+(v-lo)/(hi-lo)*(W-8);
  return `<svg class="sens-strip" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}"><line x1="4" x2="${W-4}" y1="${H/2}" y2="${H/2}" class="ax"/>${vals.map(v=>`<circle cx="${X(v).toFixed(1)}" cy="${H/2}" r="2.6" fill="var(--s1)" opacity=".55"/>`).join("")}${base!=null?`<line x1="${X(base).toFixed(1)}" x2="${X(base).toFixed(1)}" y1="2" y2="${H-2}" stroke="var(--ink)" stroke-width="2"/>`:""}</svg>`; }
function sensHTML(sp,p){
  const o=sensOpts(), r=p.result, sh=shOf(sp); const wide=p.hi-p.lo>=30;
  let nVar=0; try{ nVar=sensVarList(sp,p).length; }catch(e){}
  const set=`<div class="xps-row"><span class="xps-f">Mode <select class="xps-sel" id="xsnM"><option value="oat" ${o.mode==="oat"?"selected":""}>Fast — one at a time</option><option value="full" ${o.mode==="full"?"selected":""}>Thorough — all combinations</option></select></span>
     <span class="xps-f" title="Move the low- and high-BE window edges in and out by this much">Window edge ± <input type="number" class="w56" step="0.1" min="0.1" id="xsnW" value="${o.dw}"> eV</span>
     <span class="xps-f">Background ${["smart","shirley","active","tougaard"].map(k=>`<label class="xps-f" ${k===p.bg.type?`title="Current background — always included (baseline)"`:k==="tougaard"&&!wide?`title="Tougaard only when the window is ≥ 30 eV"`:""}><input type="checkbox" data-snb="${k}" ${o.bgs[k]||k===p.bg.type?"checked":""} ${k===p.bg.type||(k==="tougaard"&&!wide)?"disabled":""}> ${esc(C.SENS_BG_LABEL[k])}</label>`).join(" ")}</span>
     <label class="xps-f" title="Also vary the endpoint-average points (3 and 10) for static backgrounds"><input type="checkbox" id="xsnN" ${o.ns?"checked":""}> Endpoint points</label>
     ${p.bg.type==="spline"&&(p.bg.knots||[]).length?`<span class="xps-f" title="Move each of your spline knots by ± this much (one at a time; auto heights are re-taken at the new place)">Knot ± <input type="number" class="w56" step="0.1" min="0.05" id="xsnK" value="${o.dk!=null?o.dk:0.5}"> eV</span>`:""}
     <button class="btn sm primary" id="xsnRun" ${r&&!r.stale?"":"disabled"}>Run sensitivity (${nVar} fits)</button></div>`;
  const S_=p.sens; let body="";
  if(S_&&S_.ver!==p.ver) body=`<div class="note warn">모형(성분·영역·배경)을 바꾼 뒤의 결과가 아닙니다 — 다시 실행하세요.</div>`;
  else if(S_){ const sum=S_.sum; const rc=r.comps; const base=S_.res.find(x=>x.i===0&&x.ok);
    const rows=sum.comps.map(cs=>{ const c=rc.find(x=>x.id===cs.id); const ci=rc.indexOf(c); const q=r.summ&&r.summ[ci]; const a=cs.areaPct; if(!a) return "";
      const sS=sigStat(r,ci); const lo=Math.max(0,a.min-2*sS), hi=Math.min(100,a.max+2*sS);
      const half=(a.max-a.min)/2; const dom=half>2*sS; const top=Object.entries(a.eff||{}).sort((x,y)=>y[1]-x[1])[0];
      const vals=S_.res.filter(x=>x.ok).map(x=>{ const cc=x.comps.find(z=>z.id===cs.id); return cc?cc.areaPct:null; }).filter(v=>v!=null);
      const be=cs.be, fw=cs.fwhm;
      return `<tr><td><i class="sw" style="background:${c?colorOf(c):"var(--faint)"}"></i><b>${esc(c?c.name||"Component":cs.id)}</b></td>
        <td class="num">${f(q&&q.areaPct!=null?q.areaPct:a.base,1)}<div class="muted xs" title="Covariance σ × residual-correlation correction; the larger of that and the Monte Carlo SD if present">± ${f(sS,1)} stat.</div></td>
        <td class="num">${f(a.min,1)} – ${f(a.max,1)}<div class="muted xs">sd ${f(a.sd,1)}</div></td>
        <td class="num"><b>${f(lo,1)} – ${f(hi,1)}</b><div class="xs"><span class="badge ${dom?"wr":"inf"}" title="${dom?"Half the spread from analysis choices exceeds the statistical 2σ — always state the background and window choices in the report":"Statistical 2σ exceeds the spread from analysis choices"}">${dom?"Analysis-dominated":"Noise-dominated"}</span></div></td>
        <td>${sensStrip(vals,a.base,Math.min(a.min,a.base!=null?a.base:a.min)-2*sS,Math.max(a.max,a.base!=null?a.base:a.max)+2*sS)}<div class="muted xs">${top&&top[1]>0?`Largest effect: ${esc(SENS_FAC[top[0]])} (${f(top[1],1)} %p)`:""}</div></td>
        <td class="num">${be?`${f(be.min+sh,2)} – ${f(be.max+sh,2)}`:"—"}<div class="muted xs">FWHM ${fw?`${f(fw.min,2)}–${f(fw.max,2)}`:"—"}</div></td></tr>`; }).join("");
    const fails=S_.res.filter(x=>!x.ok); const why=[...new Set(fails.map(x=>x.err))];
    const ex=sum.comps.map(cs=>{ const c=rc.find(x=>x.id===cs.id), a=cs.areaPct; if(!a||!c) return ""; const ci=rc.indexOf(c); const sS=sigStat(r,ci);
      const bv=r.summ[ci]&&r.summ[ci].areaPct!=null?r.summ[ci].areaPct:a.base;
      return `${c.name||"Component"} ${f(bv,1)} % (stat. σ ${f(sS,1)}; analysis range ${f(a.min,1)}–${f(a.max,1)} %; plausible range ${f(Math.max(0,a.min-2*sS),1)}–${f(Math.min(100,a.max+2*sS),1)} %)`; }).filter(Boolean).join(" · ");
    body=(sum.base?"":`<div class="note warn">기준(0 번, 지금 설정) 변형이 실패했습니다 — 범위를 믿기 어렵습니다. 변형 목록에서 이유를 확인하세요.</div>`)+`<div class="tw"><table class="xps-tbl sens"><thead><tr><th>Component</th><th class="num">Area % baseline</th><th class="num" title="Minimum–maximum over the variants">Analysis range</th><th class="num" title="[min − 2σ stat., max + 2σ stat.] — the interval to report">Plausible range</th><th>Values per variant (│ = baseline)</th><th class="num">BE · FWHM range</th></tr></thead><tbody>${rows}</tbody></table></div>
      <div class="muted xs" style="margin-top:4px">변형 ${sum.nOk}/${sum.n}개 성공${fails.length?` (건너뜀 ${fails.length}: ${esc(why.join(", "))})`:""} · χ²ν ${f(sum.chi2nu.min,2)}–${f(sum.chi2nu.max,2)} · ${(S_.ms/1000).toFixed(1)} s. 변형들은 같은 무게(창·배경이 달라 χ²·BIC 로 서로 비교할 수 없음).</div>
      <div class="note" style="margin-top:6px;font-size:12px"><b>Example for reporting</b>: ${esc(ex)}</div>
      <details class="xps-d"><summary>${sum.n} variants in detail</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th>#</th><th>Background</th><th>Window (corrected BE)</th><th class="num">χ²ν</th>${rc.map(c=>`<th class="num">${esc(c.name||"Component")} %</th>`).join("")}</tr></thead><tbody>
      ${S_.variants.map((v,i)=>{ const x=S_.res.find(z=>z.i===i); return `<tr${i===0?` class="tot"`:""}><td>${i===0?"Baseline":i}</td><td>${esc(C.SENS_BG_LABEL[v.f.bg]||(v.f.bg==="spline"?"Spline":v.f.bg))}${v.f.n!=null?` · ${v.f.n} pts`:""}${v.f.kn?` · knot ${parseInt(v.f.kn)+1} ${v.f.kn.slice(-1)==="-"?"−":"+"}`:""}</td><td>${f(v.lo+sh,2)} – ${f(v.hi+sh,2)}</td><td class="num">${x&&x.ok?f(x.chi2nu,2):`<span class="badge wr" title="${esc(x?x.err:"")}">Skipped</span>`}</td>${rc.map(c=>{ const cc=x&&x.ok&&x.comps.find(z=>z.id===c.id); return `<td class="num">${cc?f(cc.areaPct,1):"—"}</td>`; }).join("")}</tr>`; }).join("")}
      </tbody></table><div class="xps-row"><button class="btn xs" id="xsnCsv">Variants CSV</button></div></div></details>`;
  }
  return `<div class="xps-h">Model Sensitivity <span class="muted xs">분석 선택에 따라 값이 얼마나 흔들리나</span></div>
   <div class="muted xs" style="line-height:1.55">성분 구성은 그대로 두고 배경 종류·끝점·창 끝을 바꿔 다시 맞춥니다. 위 표의 ± 는 잡음(통계)만이고, 여기 <b>Plausible range</b> = [Analysis range − 2σ, + 2σ] 가 분석자의 선택까지 넣은 보고용 구간입니다. 합성 시험(40 회, 참 배경이 후보 중에 있을 때): 지금 배경이 틀렸으면 통계 ±2σ 는 참값을 12 % 만, Plausible range 는 93 % 포함 — 참 배경이 후보 밖이면 보장되지 않음.</div>
   ${set}${body}`;
}
function bindSens(el,sp,p){
  const o=sensOpts(); const m=$q(el,"#xsnM"); if(!m) return;
  m.onchange=()=>{ o.mode=m.value; save(); renderPanel(); };
  $q(el,"#xsnW").onchange=e=>{ const v=+e.target.value; o.dw=isFinite(v)&&e.target.value!==""?Math.max(0.1,Math.min(v,5)):0.5; save(); renderPanel(); };
  $a(el,"[data-snb]").forEach(c=>c.onchange=()=>{ o.bgs[c.dataset.snb]=c.checked; save(); renderPanel(); });
  $q(el,"#xsnN").onchange=e=>{ o.ns=e.target.checked; save(); renderPanel(); };
  const kk=$q(el,"#xsnK"); if(kk) kk.onchange=()=>{ const v=+kk.value; o.dk=isFinite(v)&&kk.value!==""?Math.max(0.05,Math.min(v,5)):0.5; save(); renderPanel(); };
  $q(el,"#xsnRun").onclick=()=>runSens(sp,p);
  const cv=$q(el,"#xsnCsv"); if(cv) cv.onclick=()=>exportSens(sp,p);
}
function exportSens(sp,p){ const S_=p.sens; if(!S_) return; const sh=shOf(sp); const rc=p.result.comps;
  const lines=[`# INNO Analytics XPS model sensitivity · ${SK(sp)} · ${sp.region} · charge shift ${sh} eV · ${S_.when}`,["variant","background","endpoint_n","lo_corr_eV","hi_corr_eV","ok","chi2nu",...rc.map(c=>(c.name||c.id)+"_areaPct"),...rc.map(c=>(c.name||c.id)+"_BEcorr"),...rc.map(c=>(c.name||c.id)+"_FWHM")].map(csvE).join(",")];
  S_.variants.forEach((v,i)=>{ const x=S_.res.find(z=>z.i===i); const get=(c,k)=>{ const cc=x&&x.ok&&x.comps.find(z=>z.id===c.id); return cc&&cc[k]!=null&&isFinite(cc[k])?(k==="be"?(cc.be+sh).toFixed(3):(+cc[k]).toFixed(4)):""; };
    lines.push([i,v.f.bg+(v.f.kn?" knot"+v.f.kn:""),v.f.n==null?"":v.f.n,(v.lo+sh).toFixed(3),(v.hi+sh).toFixed(3),x&&x.ok?1:0,x&&x.ok?x.chi2nu.toFixed(4):"",...rc.map(c=>get(c,"areaPct")),...rc.map(c=>get(c,"be")),...rc.map(c=>get(c,"fwhm"))].map(csvE).join(",")); });
  download(safe(SK(sp)+"_"+sp.region)+"_sensitivity.csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"})); }

function modelsHTML(p){ const L=p.models||[]; if(!L.length) return ""; const same=L.filter(m=>Math.abs(m.lo-p.lo)<1e-6&&Math.abs(m.hi-p.hi)<1e-6); const amin=Math.min(...L.map(m=>m.stats.AICc)), bmin=Math.min(...L.map(m=>m.stats.BIC));
  return `<details class="xps-d" open><summary>Model Comparison (${L.length})</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th>Model</th><th class="num">p</th><th class="num">χ²ν</th><th class="num">ΔAICc</th><th class="num">ΔBIC</th><th></th></tr></thead><tbody>${L.map((m,i)=>`<tr><td>${esc(m.name)}${Math.abs(m.lo-p.lo)>1e-6||Math.abs(m.hi-p.hi)>1e-6?` <span class="badge wr">Different range</span>`:""}</td><td class="num">${m.stats.p}</td><td class="num">${f(m.stats.chi2nu,3)}</td><td class="num">${f(m.stats.AICc-amin,1)}</td><td class="num">${f(m.stats.BIC-bmin,1)}</td><td><button class="btn xs" data-mload="${i}">Load</button> <button class="btn xs" data-mdel="${i}">×</button></td></tr>`).join("")}</tbody></table><div class="muted xs" style="margin-top:4px">같은 영역·같은 데이터에서만 비교 의미가 있습니다. ΔAICc &lt; 2 는 사실상 구별 불가, &gt; 10 은 강한 근거. 성분을 늘려 χ² 가 조금 줄어도 AICc 가 커지면 과적합입니다.${same.length<L.length?" (영역이 다른 모형 포함)":""}</div></div></details>`; }
function bindModels(el,sp,p){ $a(el,"[data-mload]").forEach(b=>b.onclick=()=>{ const m=p.models[+b.dataset.mload]; pushHist(p,"모형 불러오기"); p.comps=C.clone(m.comps); p.cons=C.clone(m.cons); p.bg=C.clone(m.bg); p.lo=m.lo; p.hi=m.hi; touch(p); S.ui.tab="fit"; rerender(); note("불러왔습니다 — 다시 피팅하면 결과가 갱신됩니다"); }); $a(el,"[data-mdel]").forEach(b=>b.onclick=()=>{ p.models.splice(+b.dataset.mdel,1); save(); renderPanel(); }); }

/* ============================================================
   정량 패널
   ============================================================ */
/* 광원 → 감도계수 표: Al Kα(1486.6) · Mg Kα(1253.6) 만 Scofield 표가 있다. 그 밖(Ag Lα, Cr Kα, 방사광)은 직접 입력 */
function srcOf(hv){ if(Math.abs(hv-1486.6)<3) return "al"; if(Math.abs(hv-1253.6)<3) return "mg"; return null; }
function rsfLookup(el,line,src,hv){
  const L=LIB(); const key=el+" "+line; const s=srcOf(hv||1486.6);
  if(src==="custom") return S.quant.custom[key]!=null?{v:+S.quant.custom[key],src:"Custom"}:null;
  if(src==="althermo1"){ if(s!=="al") return null; return L.althermo1&&L.althermo1[key]!=null?{v:L.althermo1[key],src:"ALTHERMO1"}:null; }
  if(!s) return null; const t=(L.scofield&&L.scofield[s])||{}; return t[key]!=null?{v:t[key],src:"Scofield "+(s==="al"?"Al":"Mg")}:null;
}
/* 이 원소에서 감도계수가 있는 선 전부(궤도 합 먼저, Scofield 큰 순) — 어떤 원소든 */
function lineChoices(el,orb){
  const t=(LIB().scofield||{}).al||{}; const ks=Object.keys(t).filter(k=>k.startsWith(el+" ")).map(k=>k.slice(el.length+1));
  const tot=ks.filter(k=>k.length===2).sort((a,b)=>t[el+" "+b]-t[el+" "+a]), sub=ks.filter(k=>k.length>2).sort((a,b)=>t[el+" "+b]-t[el+" "+a]);
  const out=[]; if(orb) out.push(orb); tot.concat(sub).forEach(k=>{ if(!out.includes(k)) out.push(k); }); return out.length?out:[orb||""];
}
/* 한 원소의 '주선'(Scofield Al 이 가장 큰 궤도 합 — 예: Ba 3d, Co 2p) 과 그 subshell 라이브러리 BE */
function mainLine(el){ const L=LIB(), t=(L.scofield||{}).al||{}, lv=((L.core||{})[el]||{}).lv||{};
  const ks=Object.keys(t).filter(k=>k.startsWith(el+" ")&&k.length-el.length-1===2).map(k=>k.slice(el.length+1)).filter(o=>Object.keys(lv).some(q=>q.startsWith(o)));
  if(!ks.length) return null; ks.sort((a,b)=>t[el+" "+b]-t[el+" "+a]); const o=ks[0];
  return {orb:o,key:el+" "+o,bes:Object.entries(lv).filter(([q])=>q.startsWith(o)).map(([q,v])=>({sub:q,be:v}))}; }
/* 같은 원소의 여러 선(예: Co 2p·Co 3p) 중 '센 선' 순위 값: 행의 선 이름으로 Scofield Al */
function lineStrength(el,line){ const t=((LIB().scofield||{}).al)||{}; const v=t[el+" "+line]; if(v!=null) return v;
  const o=(line||"").replace(/(1\/2|3\/2|5\/2|7\/2)$/,""); return t[el+" "+o]!=null?t[el+" "+o]:0; }
function quantRows(sample){
  const rows=[]; const skipped=[];
  sampleSpecs(sample).filter(s=>s.reg&&s.reg.kind==="region").forEach(sp=>{
    const p=pr(sp); const r=p.result&&!p.result.stale?p.result:null; const sh=shOf(sp);
    // 원소별 묶음
    const els={};
    if(r){ r.comps.forEach((c,i)=>{ if(c.noQuant) return; const e=c.el||sp.reg.el||"?"; (els[e]=els[e]||[]).push(i); }); }
    const keys=Object.keys(els).length?Object.keys(els):(sp.reg.el?[sp.reg.el]:[]);
    keys.forEach(e=>{
      const id=sp.uid+"|"+e; const ov=S.quant.rows[id]||{};
      let area=null, sArea=null, src="", BE=null;
      if(r&&els[e]){ const idx=els[e]; area=idx.reduce((s,i)=>s+r.summ[i].area,0); BE=idx.reduce((s,i)=>s+r.summ[i].be*r.summ[i].area,0)/(area||1);
        if(r.mc&&r.mc.summ&&idx.length===r.comps.length&&r.mc.areaSum) sArea=r.mc.areaSum.sd; else if(idx.length===r.comps.length&&r.sumSig) sArea=r.sumSig.all; else sArea=Math.sqrt(idx.reduce((s,i)=>s+Math.pow(r.summ[i].sig.area||0,2),0));
        src=`Fit, ${idx.length} comp.`; }
      else { const R=region(sp,p); if(!R||R.E.length<=8){ skipped.push({sp,el:e,n:R?R.E.length:0}); return; } if(R&&R.E.length>8){ const E=Float64Array.from(R.E), y=Float64Array.from(R.y); const B=C.bgStatic(E,y,p.bg.type==="active"?{type:"shirley",n:p.bg.n}:bgOf(sp,p,R))||C.bgShirley(E,y,5,false); const Q=C.cumtrapz(E,Float64Array.from(y,(v,i)=>v-B[i])); area=Q[Q.length-1]; let w=0,ws=0; for(let i=0;i<E.length;i++){ const t=Math.max(0,y[i]-B[i]); w+=t*E[i]; ws+=t; } BE=ws>0?w/ws:(E[0]+E[E.length-1])/2; src="Region integral (above background)"; } }
      if(area==null) return;
      const cl=r&&els[e]?(r.comps[els[e][0]].line||""):"";
      const line=ov.line||(e===sp.reg.el?sp.reg.orb:(cl?cl.replace(/(1\/2|3\/2|5\/2|7\/2)$/,""):""));
      const lk=rsfLookup(e,line||"",S.quant.lib,sp.hv);
      const S_=ov.rsf!=null?+ov.rsf:(lk?lk.v:null);
      const tf=tfFor(sp);
      rows.push({id,sp,el:e,line,area,sArea,BE,src,S:S_,Ssrc:ov.rsf!=null?"Manual":(lk?lk.src:"None"),hv:sp.hv,PE:sp.meta&&sp.meta.passEnergy,tf,incSet:ov.inc,sh});
    });
  });
  // 기본 포함: 원소마다 가장 센 선(Scofield) 하나만 — 같은 원소의 다른 선(Co 2p 와 Co 3p)이 둘 다 Σ 에 들어가지 않게. 사용자가 체크한 것은 그대로
  const best={}; rows.forEach(r=>{ const v=lineStrength(r.el,r.line); if(!best[r.el]||v>best[r.el].v) best[r.el]={v,id:r.id}; });
  rows.forEach(r=>{ r.weaker=best[r.el].id!==r.id; r.inc=r.incSet!=null?r.incSet!==false:!r.weaker; });
  rows.skipped=skipped;
  return rows;
}
/* 창 안에 다른 원소의 주선: 포함된 원소 + 이 시료 survey 에서 채택한 원소의 주선 BE(라이브러리, 원소 상태)가 행의 피팅 창(보정 BE) 안이면 경고.
   예: Co 2p 창(775–810 eV) 안의 Ba 3d5/2 780 eV — 넓이에 Ba 가 섞인다 */
function interferences(sample,rows){
  const els=new Set(rows.filter(r=>r.inc).map(r=>r.el));
  sampleSpecs(sample).forEach(s=>{ const st=S.ui.survey[s.uid]; if(st&&st.acc) Object.keys(st.acc).forEach(e=>{ if(st.acc[e]) els.add(e); }); });
  const out=[], seen=new Set();
  rows.filter(r=>r.inc).forEach(r=>{ const p=pr(r.sp), lo=p.lo+r.sh, hi=p.hi+r.sh;
    els.forEach(e=>{ if(e===r.el) return; const m=mainLine(e); if(!m) return; const hit=m.bes.filter(b=>b.be>=lo&&b.be<=hi); if(!hit.length) return;
      const k=r.id+"|"+e; if(seen.has(k)) return; seen.add(k);
      out.push(`${r.el} ${r.line||""} 창(${f(lo,0)}–${f(hi,0)} eV) 안에 ${e} ${hit.map(b=>b.sub).join("/")} 주선(${hit.map(b=>f(b.be,1)).join(", ")} eV) — 넓이에 ${e} 가 섞였을 수 있습니다`); }); });
  return out;
}
function tfFor(sp){
  const q=S.quant.tf; const own=sp.meta&&sp.meta.tf;
  const av=t=>S.quant.tfAv&&t&&t.model==="thermo"?Object.assign({},t,{model:"avantage"}):t;   // Avantage 설명서 형태(log10 · ×PE, SPEC §4) — 선택
  if(q.model==="auto") return av(own?own:(S.quant.tfProfile?S.quant.tfProfile:{model:"none"}));
  if(q.model==="power") return {model:"power",x:+q.x||1,pe:sp.meta&&sp.meta.passEnergy};
  if(q.model==="profile") return av(S.quant.tfProfile||{model:"none"});
  return {model:"none"};
}
function panelQuant(el,sp){
  const smpK=SK(sp);
  const rows=quantRows(smpK);
  const used=rows.filter(r=>r.inc&&r.S>0);
  const qo={imfp:S.quant.imfp};
  const Q=used.length?C.quantify(used.map(r=>({label:r.el,el:r.el,area:r.area,sArea:r.sArea,S:r.S,BE:r.BE,hv:r.hv,PE:r.PE,tf:r.tf})),qo):[];
  const qi=new Map(); used.forEach((r,i)=>qi.set(r.id,Q[i]));
  // at% Monte Carlo (CasaXPS 식): 넓이를 N(A, σ_A) 에서 2000 번 → 95 % 구간. 넓이 상대오차가 크면 델타 ±σ 와 다름
  const QM=used.length?C.quantMC(used.map(r=>({label:r.el,el:r.el,area:r.area,sArea:r.sArea,S:r.S,BE:r.BE,hv:r.hv,PE:r.PE,tf:r.tf})),qo,2000,7):[]; const qmi=new Map(); used.forEach((r,i)=>qmi.set(r.id,QM[i]));
  const anyVgdTF=S.specs.some(s=>s.meta&&s.meta.tf);
  const noTF=used.some(r=>!r.tf||r.tf.model==="none");
  const peMix=new Set(used.map(r=>r.PE).filter(v=>v!=null)).size>1;
  const dup=Object.entries(used.reduce((m,r)=>{ m[r.el]=(m[r.el]||0)+1; return m; },{})).filter(([k,v])=>v>1).map(([k])=>k);
  const nds=used.filter(r=>{ const q=qi.get(r.id); return q&&q.nd; }).map(r=>r.el+" "+(r.line||""));
  const intf=interferences(smpK,rows);
  const skipped=rows.skipped||[];
  const noSrc=[...new Set(rows.filter(r=>!srcOf(r.hv)).map(r=>f(r.hv,1)))];
  const dsIn=used.filter(r=>{ const res=pr(r.sp).result; return res&&!res.stale&&res.comps.some(c=>c.shape==="DS"&&!c.noQuant&&(c.el||r.sp.reg.el)===r.el); }).map(r=>r.el+" "+(r.line||"")).join(", ");
  const im=S.quant.imfp, tq=S.quant.tf;
  el.innerHTML=`
   <div class="xps-h">Sample <b>${esc(smpDisp(smpK).name)}</b> <span class="muted xs">${esc(smpDisp(smpK).ctx)}</span></div>${chgBanner(sp)}
   <div class="xps-row"><span class="xps-f">RSF library <select class="xps-sel" id="xqLib"><option value="scofield">Scofield 1976 (Al Kα)</option><option value="althermo1">Avantage ALTHERMO1</option><option value="custom">Custom values</option></select></span></div>
   <div class="xps-row"><span class="xps-f">Transmission <select class="xps-sel" id="xqTf"><option value="auto">From file (VGD/avg), else none</option><option value="power">(PE/KE)^x</option><option value="none">None</option></select></span>${tq.model==="power"?`<span class="xps-f">x <input type="number" class="w56" step="0.05" id="xqTfx" value="${tq.x}"></span>`:""}${tq.model==="auto"&&S.specs.some(s=>s.meta&&s.meta.tf)?`<label class="xps-f" title="Read the VGD/avg transmission coefficients as in the Avantage manual: log10[A/(PE·XSF)] = a + b·x + c·x² + d·x³, x = log10(KE/PE) → T = PE·10^poly. Off = earlier natural-log reading (SPEC §4)"><input type="checkbox" id="xqTfAv" ${S.quant.tfAv?"checked":""}> Avantage form (log10)</label>`:""}</div>
   <div class="xps-row"><span class="xps-f">IMFP <select class="xps-sel" id="xqIm"><option value="power">KE^n</option><option value="tpp2m">TPP-2M (enter material)</option><option value="none">None</option></select></span>${im.model==="power"?`<span class="xps-f">n <input type="number" class="w56" step="0.05" id="xqImn" value="${im.n}"></span>`:""}
     ${im.model==="tpp2m"?`<span class="xps-f" title="Density (g/cm³)">ρ <input type="number" class="w56" step="0.1" id="xqR" value="${im.mat.rho}"></span><span class="xps-f" title="Valence electrons per molecule (or atom)">Nv <input type="number" class="w56" id="xqNv" value="${im.mat.Nv}"></span><span class="xps-f" title="Molar mass (g/mol)">M <input type="number" class="w64" id="xqM" value="${im.mat.M}"></span><span class="xps-f" title="Band gap (eV)">Eg <input type="number" class="w56" step="0.1" id="xqEg" value="${im.mat.Eg}"></span>`:""}</div>
   <div class="tw" style="margin-top:6px"><table class="xps-tbl xq-tbl"><thead><tr><th></th><th>Element · line</th><th class="num">at %</th><th class="num">BE</th><th class="num">Area</th><th class="num">RSF (S)</th><th class="num">T</th><th class="num">λ</th></tr></thead><tbody>
   ${rows.map(r=>{ const q=qi.get(r.id); const atC=`<td class="num xq-at">${q&&q.nd?`<span class="badge wr" title="Area ≤ 0 or below 2σ — not detected; excluded from the sum">n.d.</span>`:`<b>${q?f(q.at,2):"—"}</b>${q&&q.sAt?`<div class="muted xs">± ${f(q.sAt,2)}</div>`:""}${q&&q.sAt&&qmi.get(r.id)?`<div class="muted xs xq-mc" title="Monte Carlo 95 % interval: 2000 draws of every area within its σ (CasaXPS-style). Differs from ±σ when area errors are large">${f(qmi.get(r.id).lo,2)}–${f(qmi.get(r.id).hi,2)}</div>`:""}`}</td>`;
     return `<tr><td><input type="checkbox" data-qinc="${esc(r.id)}" ${r.inc?"checked":""} ${r.weaker&&r.incSet==null?`title="Another line of ${esc(r.el)} is stronger (Scofield) — excluded by default"`:""}></td><td><b>${esc(r.el)}</b> <select class="xps-sel" data-qline="${esc(r.id)}" style="padding:1px 3px;font-size:11.5px">${lineChoices(r.el,r.sp.reg.el===r.el?r.sp.reg.orb:r.line).filter(Boolean).map(l=>`<option ${l===r.line?"selected":""}>${esc(l)}</option>`).join("")}${r.sp.reg.el!==r.el?`<option value="">Manual</option>`:""}</select><div class="muted xs">${esc(r.sp.region)} · ${esc(r.src)}${r.weaker&&!r.inc?" · weaker line":""}</div></td>${atC}
     <td class="num">${f(r.BE+r.sh,2)}</td><td class="num">${g(r.area,5)}${r.sArea?`<div class="muted xs">± ${g(r.sArea,3)}</div>`:""}</td>
     <td class="num"><input type="number" class="w56" step="any" data-qrsf="${esc(r.id)}" value="${r.S!=null?+(+r.S).toPrecision(5):""}" placeholder="?"><div class="muted xs">${esc(r.Ssrc)}</div></td>
     <td class="num">${q?g(q.T,4):"—"}</td><td class="num">${q?g(q.L,4):"—"}</td></tr>`; }).join("")}
   ${skipped.map(k=>`<tr><td></td><td colspan="7"><b>${esc(k.el||"?")}</b> <span class="muted xs">${esc(k.sp.region)}</span> <span class="badge wr" title="The fit range has ${k.n} data points (≥ 9 needed). Widen it in Peak Fit">Fit range too narrow</span></td></tr>`).join("")}
   ${!rows.length&&!skipped.length?`<tr><td colspan="8" class="muted">No region spectra in this sample.</td></tr>`:""}
   </tbody></table></div>
   <div class="muted xs" style="margin-top:6px;line-height:1.55">at% = (A/(S·T·λ)) / Σ. A: 배경 위 넓이(피팅이 있으면 성분 합, 원소별 — 성분 카드의 'Element'·'Quantify' 표시로 나눔). σ 는 넓이 불확도만 전파(감도계수·T·λ 의 계통 오차 제외).</div>
   ${noTF?`<div class="note warn" style="margin-top:6px">Transmission 보정이 없는 행이 있습니다 (txt/xls/VAMAS 에는 장비 함수가 없음). 같은 장비의 VGD 를 하나 열어 'Use VGD transmission as default' 를 누르거나 (PE/KE)^x 를 고르세요. Avantage 결과와 비교할 땐 같은 보정을 써야 합니다.</div>`:""}
   ${anyVgdTF?`<div class="xps-row"><button class="btn xs" id="xqProf" title="Also apply the transmission coefficients read from VGD to the other files (txt/xls)">Use VGD transmission as default</button>${S.quant.tfProfile?`<span class="badge inf">Transmission saved (PE ${esc(S.quant.tfProfile.pe)})</span>`:""}</div><div class="muted xs">VGD 의 transmission 계수 해석 ln T = Σaₖ(ln(KE/PE))ᵏ 은 계수 값의 거동(높은 KE/PE 에서 T∝(KE/PE)^≈−0.9)으로 추정한 것입니다 — Avantage 의 같은 데이터 at% 와 한 번 대조해 확인하세요.</div>`:""}
   ${peMix?`<div class="note warn" style="margin-top:6px">pass energy 가 다른 영역이 섞여 있습니다 — transmission 보정 없이는 비교할 수 없습니다.</div>`:""}
   ${dup.length?`<div class="note warn" style="margin-top:6px">같은 원소가 여러 행: ${esc(dup.join(", "))} — 한 원소당 한 선만 체크하세요(예: Co 2p 와 Co 3p 중 하나).</div>`:""}
   ${nds.length?`<div class="note warn" style="margin-top:6px">검출 안 됨(넓이 ≤ 0 또는 2σ 미만): ${esc(nds.join(", "))} — Σ 에서 뺐습니다.</div>`:""}
   ${intf.map(t=>`<div class="note warn" style="margin-top:6px">${esc(t)}</div>`).join("")}
   ${noSrc.length?`<div class="note warn" style="margin-top:6px">광원 hν = ${esc(noSrc.join(", "))} eV 에 대한 감도계수 표가 없습니다(Scofield 표는 Al Kα·Mg Kα 만). 감도계수를 직접 넣으세요.</div>`:""}
   ${dsIn?`<div class="note warn" style="margin-top:6px">DS(Doniach–Šunjić) 성분이 있습니다(${esc(dsIn)}): DS 넓이는 피팅 창 안의 적분이라 창 폭에 따라 달라집니다. 같은 창 규칙으로 잰 시료끼리만 at% 를 비교하세요.</div>`:""}
   <div class="muted xs" style="margin-top:4px">ALTHERMO1 은 Scofield 대비 선마다 다르지만(5–95 % 범위 −26…+49 %) KE 와 상관이 없다(log 기울기 0.01) — T·λ 가 들어 있지 않은 값이므로 T·λ 보정을 함께 적용한다.</div>
   ${stoichHTML(sp,used,Q)}
   <div class="xps-row" style="margin-top:6px"><button class="btn sm" id="xqCsv">Quantification CSV</button></div>`;
  $q(el,"#xqLib").value=S.quant.lib; $q(el,"#xqTf").value=tq.model==="profile"?"auto":tq.model; $q(el,"#xqIm").value=im.model;
  $q(el,"#xqLib").onchange=e=>{ S.quant.lib=e.target.value; save(); renderPanel(); };
  $q(el,"#xqTf").onchange=e=>{ S.quant.tf.model=e.target.value; save(); renderPanel(); };
  const tav=$q(el,"#xqTfAv"); if(tav) tav.onchange=()=>{ S.quant.tfAv=tav.checked; save(); renderPanel(); };
  const tx=$q(el,"#xqTfx"); if(tx) tx.onchange=()=>{ S.quant.tf.x=+tx.value; save(); renderPanel(); };
  $q(el,"#xqIm").onchange=e=>{ S.quant.imfp.model=e.target.value; save(); renderPanel(); };
  const imn=$q(el,"#xqImn"); if(imn) imn.onchange=()=>{ S.quant.imfp.n=+imn.value; save(); renderPanel(); };
  [["xqR","rho"],["xqNv","Nv"],["xqM","M"],["xqEg","Eg"]].forEach(([id,k])=>{ const i=$q(el,"#"+id); if(i) i.onchange=()=>{ S.quant.imfp.mat[k]=+i.value; save(); renderPanel(); }; });
  $a(el,"[data-qinc]").forEach(i=>i.onchange=()=>{ const o=(S.quant.rows[i.dataset.qinc]=S.quant.rows[i.dataset.qinc]||{}); o.inc=i.checked; save(); renderPanel(); });
  $a(el,"[data-qline]").forEach(i=>i.onchange=()=>{ const o=(S.quant.rows[i.dataset.qline]=S.quant.rows[i.dataset.qline]||{}); o.line=i.value; delete o.rsf; save(); renderPanel(); });
  $a(el,"[data-qrsf]").forEach(i=>i.onchange=()=>{ const o=(S.quant.rows[i.dataset.qrsf]=S.quant.rows[i.dataset.qrsf]||{}); if(i.value==="") delete o.rsf; else o.rsf=+i.value; save(); renderPanel(); });
  const pf=$q(el,"#xqProf"); if(pf) pf.onclick=()=>{ const s=S.specs.find(x=>x.meta&&x.meta.tf); S.quant.tfProfile=C.clone(s.meta.tf); save(); renderPanel(); note("VGD 의 transmission 계수를 기본 장비 함수로 저장했습니다 (PE "+s.meta.tf.pe+" eV)","ok"); };
  $q(el,"#xqCsv").onclick=()=>exportQuant(smpK,rows,qi,qmi);
  bindStoich(el,sp);
}

/* ---------------- 화학량론 비교 (정량 탭 아래) ---------------- */
const ZOF=e=>{ const i=(IO.SYMS||[]).indexOf(e); return i<0?0:i+1; };
function stoichHTML(sp,used,Q){
  const sm=smpOf(sp); const fs=sm.formula||"";
  let body="", fm=null, err="";
  if(fs){ try{ fm=C.parseFormula(fs,IO.SYMS); }catch(e){ err=e.message; } }
  if(fm&&used.length){
    // 핵심 계산과 같은 걸러내기(넓이 정규화 값이 유한하고 양수인 행만, 같은 원소는 첫 행) — 화면과 계산이 다른 행을 쓰지 않게
    const rows=used.map((r,i)=>({el:r.el,nrm:Q[i].nrm,sNrm:Q[i].sNrm,KE:Q[i].KE,L:Q[i].L,line:r.line,nd:Q[i].nd})).filter(r=>!r.nd&&isFinite(r.nrm)&&r.nrm>0);
    const rowOf=e=>rows.find(r=>r.el===e);
    const inF=Object.keys(fm).filter(e=>rowOf(e));
    const refDef=inF.slice().sort((a,b)=>ZOF(b)-ZOF(a))[0]; const R=C.stoichCompare(rows,fm,sm.stRef&&inF.includes(sm.stRef)?sm.stRef:refDef); const ref=R.ref;
    const im=S.quant.imfp; const Lref=rowOf(ref)&&rowOf(ref).L;
    const depth=e=>{ const r=rowOf(e); if(!r) return "—"; if(im.model==="tpp2m") return `${f(3*r.L/10,1)} nm`; if(im.model==="power"&&Lref) return `× ${f(r.L/Lref,2)}`; return "—"; };
    const fmtR=(v,s)=>`<b class="${Math.abs(v-1)>Math.max(0.1,s!=null?2*s:0)?"wr":""}">${f(v,2)}</b>${s!=null?`<span class="muted xs"> ± ${f(s,2)}</span>`:""}`;
    body=`<div class="tw"><table class="xps-tbl"><thead><tr><th>Element · line</th><th class="num">KE</th><th class="num" title="% renormalized over the elements that are both in the formula and measured">Measured %</th><th class="num">Formula %</th><th class="num">Meas./formula</th><th class="num" title="(X/ref) measured ÷ (X/ref) formula">÷ ${esc(ref||"")}</th><th class="num" title="${im.model==="tpp2m"?"Information depth 3λ (TPP-2M)":"λ ratio to the reference element — are the same depths probed?"}">${im.model==="tpp2m"?"Depth 3λ":"λ/λ(ref)"}</th></tr></thead><tbody>
      ${R.items.map(it=>{ const r=rowOf(it.el); return `<tr><td><b>${esc(it.el)}</b> <span class="muted xs">${esc(r.line||"")}</span></td><td class="num">${f(r.KE,0)}</td><td class="num">${f(100*it.x,2)}${it.sx!=null?`<div class="muted xs">± ${f(100*it.sx,2)}</div>`:""}</td><td class="num">${f(100*it.ex,2)}</td><td class="num">${fmtR(it.R,it.sR)}</td><td class="num">${it.el===ref?`<span class="muted">ref</span>`:fmtR(it.ratioRel,it.sRatioRel)}</td><td class="num">${depth(it.el)}</td></tr>`; }).join("")}
    </tbody></table></div>`;
    if(!R.items.length) return `<hr class="xps-hr"><div class="xps-h">Stoichiometry Comparison</div><div class="row"><span class="xps-f">Formula <input type="text" class="xps-in" id="xstF" style="width:220px" value="${esc(fs)}"></span><button class="btn xs" id="xstOk">Compare</button></div><div class="note warn">화학식의 원소 중 정량 표에서 체크된 것이 없습니다.</div>`;
    const kes=R.items.map(it=>rowOf(it.el).KE); const kr=Math.max(...kes)/Math.min(...kes);
    const notes=[];
    if(R.dup.length) notes.push(`같은 원소가 여러 행: ${R.dup.join(", ")} — 정량 표에서 한 원소당 한 선만 체크하세요 (첫 행만 썼습니다).`);
    if(R.missing.length) notes.push(`화학식에 있지만 측정(체크)되지 않은 원소: ${R.missing.join(", ")} — 비교에서 빼고 나머지끼리 정규화했습니다.`);
    if(R.outside.length) notes.push(`화학식 밖 원소 ${R.outside.join(", ")}(오염·리간드 등)는 비교에서 뺐습니다 — 정량 표의 at% 에는 들어 있습니다.`);
    if(kr>1.8) notes.push(`선들의 운동에너지가 ${f(Math.min(...kes),0)}–${f(Math.max(...kes),0)} eV 로 크게 달라 보는 깊이가 다릅니다(λ ∝ KE^≈0.75). 표면 조성이 깊이에 따라 다르면 측정/화학식 비가 1 에서 벗어나는 것이 '조성'이 아니라 '깊이 분포' 때문일 수 있습니다 — 같은 원소의 KE 가 다른 두 선(예: Na 1s 와 Na 2s)을 비교하면 구분됩니다.`);
    if(im.model==="none") notes.push("IMFP 보정이 꺼져 있어 깊이 비를 표시하지 않습니다.");
    body+=notes.map(t=>`<div class="note warn" style="margin-top:4px">${esc(t)}</div>`).join("")
      +`<div class="xps-row"><span class="xps-f">Reference element <select class="xps-sel" id="xstRef">${inF.map(e=>`<option ${e===ref?"selected":""}>${esc(e)}</option>`).join("")}</select></span><span class="muted xs">굵은 주황 = 1 에서 0.1 이상, σ 의 2 배 이상 벗어남. σ 는 넓이 불확도만(감도계수·T·λ 의 계통 오차 제외 — 감도계수 오차만으로도 흔히 ±10–20 %)</span></div>`;
  } else if(fm) body=`<div class="muted sm">정량 표에서 체크된 원소가 없습니다.</div>`;
  return `<hr class="xps-hr"><div class="xps-h">Stoichiometry Comparison <span class="muted xs">측정 비율 ↔ 화학식</span></div>
   <div class="xps-row"><span class="xps-f">Formula <input type="text" class="xps-in" id="xstF" style="width:220px" value="${esc(fs)}" placeholder="e.g. NaLuF4, Na(Lu0.92Tm0.08)F4, Na3LuCl6"></span><button class="btn xs" id="xstOk">Compare</button>${fs?`<button class="btn xs" id="xstX">Clear</button>`:""}</div>
   ${err?`<div class="note warn">화학식을 읽지 못했습니다: ${esc(err)}</div>`:""}${body}`;
}
function bindStoich(el,sp){ const sm=smpOf(sp); const inp=$q(el,"#xstF"); if(!inp) return;
  const go=()=>{ sm.formula=inp.value.trim(); save(); renderPanel(); };
  $q(el,"#xstOk").onclick=go; inp.onkeydown=e=>{ if(e.key==="Enter") go(); };
  const x=$q(el,"#xstX"); if(x) x.onclick=()=>{ sm.formula=""; save(); renderPanel(); };
  const rf=$q(el,"#xstRef"); if(rf) rf.onchange=()=>{ sm.stRef=rf.value; save(); renderPanel(); }; }

/* ============================================================
   Auger parameter 패널 (SPEC §10) — α′ = BE(광전자) + KE(Auger) = BE_p + hν − BE_A(겉보기). 대전 보정과 무관.
   항목 = {id, sample, el, p:{uid,comp,lo,hi}, a:{uid,comp,lo,hi}} — comp 가 성분 id 면 피팅된 그 성분의 BE(σ 있음), 비면 구간 [lo,hi](측정 BE) 의 최댓값.
   기준점(문헌·표준 시료)은 사용자가 넣는다 — 앱에 문헌 α′ 표를 싣지 않음(검증된 자료가 없어서).
   ============================================================ */
function AUG(){ if(!S.auger||typeof S.auger!=="object") S.auger={entries:[],refs:[]}; S.auger.entries=S.auger.entries||[]; S.auger.refs=S.auger.refs||[]; return S.auger; }
const augId=()=>"a"+Math.random().toString(36).slice(2,8);
function augSide(s,side){
  if(!side.uid) return {err:"스펙트럼을 고르세요"};
  const sp=byUid(side.uid); if(!sp) return {err:"스펙트럼이 목록에 없습니다"};
  if(side.comp){ const r=fitted(sp); if(!r) return {err:`${sp.region}: 피팅 결과가 없거나 수정 뒤 다시 피팅하지 않았습니다`};
    const i=r.comps.findIndex(c=>c.id===side.comp); if(i<0) return {err:`${sp.region}: 고른 성분이 없어졌습니다 — 다시 고르세요`};
    const q=r.summ&&r.summ[i]; return {sp,be:r.comps[i].be,s:q&&q.sig&&q.sig.be!=null?q.sig.be:null,how:`성분 ‘${r.comps[i].name||"Component"}’`}; }
  const lo=side.lo!=null?side.lo:sp.be[0], hi=side.hi!=null?side.hi:sp.be[sp.n-1];
  const m=C.peakMaxPos(sp.be,sp.y,Math.min(lo,hi),Math.max(lo,hi)); if(!m) return {err:`${sp.region}: 구간 ${f(lo,1)}–${f(hi,1)} eV 의 끝이 가장 높습니다 — 봉우리를 감싸도록 구간을 고치세요`};
  return {sp,be:m.be,s:null,how:`최댓값(${f(Math.min(lo,hi),1)}–${f(Math.max(lo,hi),1)})`};
}
/* survey 에서 Auger 선을 찾을 기본 구간: 라이브러리(Avantage KE)의 가장 센 Auger 선 ±6 eV (측정 BE; 대전 이동 Δ 만큼 옮김) */
function augWin(elx,s){ const al=(LIB().auger||[]).filter(a=>a.el===elx).sort((a,b)=>b.rsf-a.rsf)[0]; if(!al||!s.hv) return null;
  const be=s.hv-al.ke-((S.samples[SK(s)]||{}).shift||0); return [+(be-6).toFixed(1),+(be+6).toFixed(1)]; }
function augCalc(e){
  const P=augSide(e,e.p), A=augSide(e,e.a); if(P.err||A.err) return {err:(P.err?"Photoelectron line: "+P.err:"Auger line: "+A.err)};
  // 쓸 수 없는 짝: Auger 쪽이 광전자 영역 스펙트럼 / 광전자 쪽이 Auger 영역 / survey 가 아닌 같은 스펙트럼 / 두 위치가 같음(α′ = hν 가 됨)
  const kA=A.sp.reg&&A.sp.reg.kind, kP=P.sp.reg&&P.sp.reg.kind;
  if(kA==="region"&&A.sp.reg.el) return {err:`Auger line 쪽이 광전자 영역 스펙트럼(${A.sp.region})입니다 — Auger 영역(예: ${e.el} KLL) 또는 survey 를 고르세요`};
  if(kP==="auger") return {err:`Photoelectron line 쪽이 Auger 영역 스펙트럼(${P.sp.region})입니다`};
  if(P.sp.uid===A.sp.uid&&kA!=="survey") return {err:"두 선이 같은 스펙트럼입니다 — 광전자 영역과 Auger 영역(또는 survey 의 다른 구간)을 고르세요"};
  if(Math.abs(P.be-A.be)<1) return {err:`두 위치가 같습니다(${f(P.be,2)} / ${f(A.be,2)} eV) — 같은 봉우리를 두 번 고르면 α′ = hν 가 됩니다. 구간을 고치세요`};
  const hv=A.sp.hv; const res=C.augerParam(P.be,A.be,hv,P.s,A.s);
  const sh=shOf(P.sp);     // 대전 보정은 광전자선 스펙트럼이 속한 시료의 것
  const warn=[]; if(SK(P.sp)!==SK(A.sp)) warn.push(`두 스펙트럼이 다른 시료입니다(${smpDisp(SK(P.sp)).name} / ${smpDisp(SK(A.sp)).name}) — α′ 는 같은 측정의 두 선이어야 대전이 지워집니다`); if(/기본값/.test(A.sp.hvSrc||"")) warn.push(`Auger 스펙트럼의 hν 가 파일에 없어 ${A.sp.hv} eV 로 둠 — 실제 광원 에너지와 다르면 α′ 가 그만큼 틀림`);
  if(P.sp.hv&&A.sp.hv&&Math.abs(P.sp.hv-A.sp.hv)>0.05) warn.push(`두 스펙트럼의 hν 가 다름(${P.sp.hv} / ${A.sp.hv}) — 같은 측정인지 확인`);
  if(A.sp.reg&&A.sp.reg.kind==="region"&&!A.sp.reg.el) warn.push(`Auger 스펙트럼(${A.sp.region})의 종류를 모릅니다 — Auger 영역이 맞는지 확인`);
  return {alpha:res.alpha,s:res.s,beP:P.be,keA:res.keA,beC:P.be+sh,keC:res.keA-sh,howP:P.how,howA:A.how,P,A,warn};
}
/* 자동 짝짓기: 이 시료의 Auger 영역 스펙트럼(파일 이름 'Na KLL' 등)마다 같은 원소의 광전자 영역 스펙트럼을 찾는다 */
function augAutoPair(k){ const A=AUG(); const specs=sampleSpecs(k); let n=0;
  specs.filter(s=>s.reg&&s.reg.kind==="auger"&&s.reg.el).forEach(as=>{
    if(A.entries.some(e=>e.sample===k&&e.a.uid===as.uid)) return;
    const ps=specs.filter(s=>s.reg&&s.reg.kind==="region"&&s.reg.el===as.reg.el); if(!ps.length) return;
    const sc=(LIB().scofield||{}).al||{}; ps.sort((x,y)=>(sc[y.reg.el+" "+y.reg.orb]||0)-(sc[x.reg.el+" "+x.reg.orb]||0)); const p=ps[0];
    const big=s=>{ const r=fitted(s); if(!r||!r.comps.length) return ""; return r.comps.slice().sort((a,b)=>b.area-a.area)[0].id; };
    A.entries.push({id:augId(),sample:k,el:as.reg.el,p:{uid:p.uid,comp:big(p),lo:null,hi:null},a:{uid:as.uid,comp:big(as),lo:null,hi:null}}); n++; });
  return n; }
function panelAuger(el,sp){
  const A=AUG(), k=SK(sp), d=smpDisp(k), specs=sampleSpecs(k);
  const ents=A.entries.filter(e=>e.sample===k);
  const nAug=specs.filter(s=>s.reg&&s.reg.kind==="auger").length;
  // 고른 것이 없어졌으면 '— 없어짐 —' 을 선택된 자리표시로 둔다(아니면 첫 항목이 선택된 것처럼 보이고, 하나뿐이면 change 가 일어나지 않아 바꿀 수 없음)
  const GONE="__gone";
  const specOpt=sel=>(specs.some(s=>s.uid===sel)?"":`<option value="${GONE}" selected>${sel?"— 없어짐: 다시 고르세요 —":"— Select —"}</option>`)+specs.map(s=>`<option value="${s.uid}" ${s.uid===sel?"selected":""}>${esc(s.region)} · ${esc(FMTN(s.format).replace("Avantage ",""))}</option>`).join("");
  const compOpt=(uid,sel)=>{ const s=byUid(uid); const r=s&&fitted(s); const has=!sel||(r&&r.comps.some(c=>c.id===sel));
    return (has?"":`<option value="${GONE}" selected>— 고른 성분 없어짐 —</option>`)+`<option value="" ${!sel?"selected":""}>Maximum (range)</option>`+(r?r.comps.map(c=>`<option value="${c.id}" ${c.id===sel?"selected":""}>Component ${esc(c.name||"Component")} · ${f(c.be,2)}</option>`).join(""):""); };
  const sideHTML=(e,side,key,lab)=>{ const s=byUid(side.uid); return `<div class="r"><span class="k">${lab}</span><select class="xps-sel" data-ax="${e.id}|${key}|uid">${specOpt(side.uid)}</select><select class="xps-sel" data-ax="${e.id}|${key}|comp">${compOpt(side.uid,side.comp)}</select>${side.comp?"":`<span class="xps-f">Range <input type="number" step="0.1" class="w64" data-ax="${e.id}|${key}|lo" value="${side.lo!=null?side.lo:(s?f(s.be[0],1):"")}"> – <input type="number" step="0.1" class="w64" data-ax="${e.id}|${key}|hi" value="${side.hi!=null?side.hi:(s?f(s.be[s.n-1],1):"")}"> eV <span class="muted xs">(measured BE)</span></span>`}</div>`; };
  const cards=ents.map(e=>{ const r=augCalc(e);
    return `<div class="aug-card"><div class="h"><b>${esc(e.el)}</b> ${r.err?`<span class="badge wr">Not computed</span>`:`<span class="al">α′ = <b>${f(r.alpha,2)}</b>${r.s!=null?` <span class="muted xs">± ${f(r.s,2)}</span>`:""} eV</span>`}<button class="x" data-axdel="${e.id}" title="Delete this entry">×</button></div>
      ${sideHTML(e,e.p,"p","Photoelectron line")}${sideHTML(e,e.a,"a","Auger line")}
      ${r.err?`<div class="note warn">${esc(r.err)}</div>`:`<div class="muted xs">BE_p ${f(r.beP,2)} (${esc(r.howP)}) · KE_A ${f(r.keA,2)} (${esc(r.howA)}) — 측정값. 대전 보정 적용 좌표: BE ${f(r.beC,2)}, KE ${f(r.keC,2)}${r.s==null?" · σ 는 두 쪽 모두 피팅된 성분일 때만 (최댓값 방식의 위치 불확도는 잡음에 따라 0.1 eV 정도)":""}</div>${r.warn.map(w=>`<div class="note warn">${esc(w)}</div>`).join("")}`}
    </div>`; }).join("");
  // Wagner plot: 원소 하나 — 모든 시료의 항목 + 기준점
  const els=[...new Set(A.entries.map(e=>e.el).concat(A.refs.map(r=>r.el)))];
  const pe=A.plotEl&&els.includes(A.plotEl)?A.plotEl:(ents[0]?ents[0].el:els[0]);
  el.innerHTML=`<div class="xps-h">Auger Parameter α′ <span class="muted xs">Sample ${esc(d.name)}</span></div>
   <div class="muted xs" style="line-height:1.55">α′ = BE(광전자선) + KE(Auger 선). 두 선이 대전으로 같은 만큼 움직이므로 <b>대전 보정과 무관</b> — 대전이 큰 절연체(플루오라이드·클로라이드 등)의 화학 환경 비교에 씁니다. 같은 측정(같은 광원·같은 대전 조건)의 두 스펙트럼을 고르세요. 위치는 피팅된 성분(σ 있음) 또는 구간 최댓값.</div>
   <div class="xps-row"><button class="btn sm primary" id="xaAuto" ${nAug?"":"disabled"} title="Pair each Auger region of this sample (KLL, LMM … in the file name) with the photoelectron region of the same element">Auto-pair</button><button class="btn sm" id="xaAdd" ${specs.length?"":"disabled"}>+ Add manually</button>${nAug?"":`<span class="muted xs">이 시료에 Auger 영역 스펙트럼이 없습니다 — Na KLL·F KLL 같은 영역을 측정해 함께 열거나, + Add manually 에서 survey 를 고르고 구간을 넣으세요</span>`}</div>
   ${cards||`<div class="muted sm" style="padding:4px 2px">No entries.</div>`}
   <hr class="xps-hr"><div class="xps-h">Wagner Plot ${els.length?`<select class="xps-sel" id="xaEl">${els.map(e=>`<option ${e===pe?"selected":""}>${esc(e)}</option>`).join("")}</select>`:""}</div>
   ${pe?wagnerSVG(pe):`<div class="muted sm">항목이나 기준점이 생기면 그려집니다.</div>`}
   <div class="xps-h">Reference Points <span class="muted xs">문헌·표준 시료 (대전 보정된 BE, Fermi 기준 KE)</span></div>
   <div class="tw"><table class="xps-tbl"><thead><tr><th>Element</th><th>Name</th><th class="num">BE</th><th class="num">KE</th><th class="num">α′</th><th></th></tr></thead><tbody>
   ${A.refs.map(r=>`<tr><td>${esc(r.el)}</td><td>${esc(r.label)}</td><td class="num">${f(r.BE,2)}</td><td class="num">${f(r.KE,2)}</td><td class="num">${f(r.BE+r.KE,2)}</td><td><button class="x" data-ardel="${r.id}">×</button></td></tr>`).join("")}
   <tr><td><input type="text" class="xps-in w56" id="xrEl" placeholder="Na" value="${esc(pe||"")}"></td><td><input type="text" class="xps-in" id="xrLb" placeholder="e.g. NaF standard" style="width:110px"></td><td class="num"><input type="number" step="0.01" class="w64" id="xrBE"></td><td class="num"><input type="number" step="0.01" class="w64" id="xrKE"></td><td></td><td><button class="btn xs" id="xrAdd">Add</button></td></tr>
   </tbody></table></div>
   <div class="muted xs">기준값은 출처를 이름에 적어 두세요(예: “NaF, NIST SRD 20”). 문헌 BE 는 대개 C 1s 284.8 또는 Fermi 준위 기준 — 이 시료의 대전 보정 기준과 같은지 확인.</div>`;
  const ax=$a(el,"[data-ax]");
  ax.forEach(i=>i.onchange=()=>{ const [id,key,fld]=i.dataset.ax.split("|"); const e=A.entries.find(x=>x.id===id); if(!e) return; const side=e[key];
    if(i.value===GONE) return;
    if(fld==="uid"){ side.uid=i.value; side.comp=""; side.lo=null; side.hi=null; const s=byUid(i.value);
      if(key==="p"&&s&&s.reg&&s.reg.el) e.el=s.reg.el;
      if(key==="a"&&s&&s.reg&&s.reg.kind==="survey"){ const w=augWin(e.el,s); if(w){ side.lo=w[0]; side.hi=w[1]; } } } else if(fld==="comp") side.comp=i.value; else { const v=+i.value; side[fld]=i.value===""||!isFinite(v)?null:v; }
    save(); renderPanel(); });
  $a(el,"[data-axdel]").forEach(b=>b.onclick=()=>{ A.entries=A.entries.filter(x=>x.id!==b.dataset.axdel); save(); renderPanel(); });
  $a(el,"[data-ardel]").forEach(b=>b.onclick=()=>{ A.refs=A.refs.filter(x=>x.id!==b.dataset.ardel); save(); renderPanel(); });
  const au=$q(el,"#xaAuto"); if(au) au.onclick=()=>{ const n=augAutoPair(k); save(); renderPanel(); note(n?`Auger 항목 ${n}개를 만들었습니다`:"새로 짝지을 것이 없습니다 (이미 있거나, 같은 원소의 광전자 영역이 없음)",n?"ok":""); };
  // 직접 추가: 광전자 쪽 = 지금 보는 광전자 영역(아니면 첫 광전자 영역), Auger 쪽 = 같은 원소 Auger 영역 → 다른 Auger 영역 → survey. 없으면 비워 두고 고르게 한다(광전자 영역으로 채우지 않음)
  const ad=$q(el,"#xaAdd"); if(ad) ad.onclick=()=>{ const isP=s=>s.reg&&s.reg.kind==="region"; const s0=isP(sp)?sp:specs.find(isP)||null;
    const elx=(s0&&s0.reg&&s0.reg.el)||"";
    const sa=specs.find(s=>s.reg&&s.reg.kind==="auger"&&s.reg.el===elx)||specs.find(s=>s.reg&&s.reg.kind==="auger")||specs.find(s=>s.reg&&s.reg.kind==="survey")||null;
    let lo=null,hi=null;
    if(sa&&sa.reg.kind==="survey"&&elx){ const w=augWin(elx,sa); if(w){ lo=w[0]; hi=w[1]; } }
    A.entries.push({id:augId(),sample:k,el:elx||"?",p:{uid:s0?s0.uid:"",comp:"",lo:null,hi:null},a:{uid:sa?sa.uid:"",comp:"",lo,hi}}); save(); renderPanel();
    if(!sa) note("이 시료에 Auger 영역·survey 스펙트럼이 없습니다 — Auger line 에서 스펙트럼을 고르세요","warn"); };
  const xe=$q(el,"#xaEl"); if(xe) xe.onchange=()=>{ A.plotEl=xe.value; save(); renderPanel(); };
  $q(el,"#xrAdd").onclick=()=>{ const raw=$q(el,"#xrEl").value.trim(); const e=IO.SYMS.find(s=>s.toLowerCase()===raw.toLowerCase()); const be=+$q(el,"#xrBE").value, ke=+$q(el,"#xrKE").value;
    if(!e){ note("원소 기호를 확인하세요","warn"); return; } if(!(be>0&&ke>0)){ note("BE 와 KE 를 넣으세요","warn"); return; }
    A.refs.push({id:augId(),el:e,label:$q(el,"#xrLb").value.trim()||"Reference",BE:be,KE:ke}); A.plotEl=e; save(); renderPanel(); };
}
function wagnerSVG(elx){
  const A=AUG(); const pts=[];
  A.entries.filter(e=>e.el===elx).forEach(e=>{ const r=augCalc(e); if(!r.err) pts.push({x:r.beC,y:r.keC,lab:smpDisp(e.sample).name,kind:"m",cur:e.sample===SK(cur()),s:r.s}); });
  A.refs.filter(r=>r.el===elx).forEach(r=>pts.push({x:r.BE,y:r.KE,lab:r.label,kind:"r"}));
  if(!pts.length) return `<div class="muted sm">${esc(elx)} 의 점이 없습니다.</div>`;
  const W=400,H=300,ml=52,mr=14,mt=12,mb=40; const xs=pts.map(p=>p.x), ys=pts.map(p=>p.y);
  // 같은 eV/px 눈금(α′ 일정선이 45°): 그림 칸이 가로로 더 길어 가로 범위를 그 비만큼 넓힌다
  const pw=W-ml-mr, ph=H-mt-mb; let x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  const cx=(x0+x1)/2, cy=(y0+y1)/2; const hx0=(x1-x0)/2, hy0=(y1-y0)/2; const padE=Math.max(1,0.25*Math.max(x1-x0,y1-y0,2));
  const k=Math.max((hx0+padE)/pw,(hy0+padE)/ph);   // eV/px
  x0=cx-k*pw; x1=cx+k*pw; y0=cy-k*ph; y1=cy+k*ph;
  const X=x=>ml+(x1-x)/(x1-x0)*(W-ml-mr), Y=y=>mt+(y1-y)/(y1-y0)*(H-mt-mb);   // BE 는 오른쪽으로 감소(관례)
  const step=[0.5,1,2,5,10,20,50].find(s=>(x1-x0)/s<=8)||100;
  let g=`<svg class="xps-svg wag" viewBox="0 0 ${W} ${H}" width="100%"><rect x="${ml}" y="${mt}" width="${W-ml-mr}" height="${H-mt-mb}" fill="none" class="ax"/>`;
  for(let v=Math.ceil(x0/step)*step;v<=x1;v+=step) g+=`<line class="gd" x1="${X(v)}" x2="${X(v)}" y1="${mt}" y2="${H-mb}"/><text x="${X(v)}" y="${H-mb+14}" text-anchor="middle">${+v.toFixed(1)}</text>`;
  for(let v=Math.ceil(y0/step)*step;v<=y1;v+=step) g+=`<line class="gd" x1="${ml}" x2="${W-mr}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${ml-4}" y="${Y(v)+4}" text-anchor="end">${+v.toFixed(1)}</text>`;
  // α′ 일정선: KE = α′ − BE
  const aMin=x0+y0, aMax=x1+y1; for(let a=Math.ceil(aMin/step)*step;a<=aMax;a+=step){ const xa=Math.max(x0,a-y1), xb=Math.min(x1,a-y0); if(xb<=xa) continue;
    g+=`<line x1="${X(xa)}" y1="${Y(a-xa)}" x2="${X(xb)}" y2="${Y(a-xb)}" stroke="var(--info)" stroke-width="1" stroke-dasharray="4 4" opacity=".55"/><text x="${X(xa)-3}" y="${Y(a-xa)+11}" text-anchor="end" class="ref" font-size="10">α′ ${+a.toFixed(1)}</text>`; }
  pts.forEach(p=>{ g+=p.kind==="r"?`<rect x="${X(p.x)-4}" y="${Y(p.y)-4}" width="8" height="8" fill="none" stroke="var(--ink2)" stroke-width="1.5"/>`:`<circle cx="${X(p.x)}" cy="${Y(p.y)}" r="${p.cur?5:4}" fill="${p.cur?"var(--accent)":"var(--s1)"}" stroke="var(--surface)" stroke-width="1.5"/>`;
    g+=`<text x="${X(p.x)+7}" y="${Y(p.y)+4}" class="cl">${esc(p.lab)}</text>`; });
  g+=`<text x="${(ml+W-mr)/2}" y="${H-6}" text-anchor="middle" class="al">BE (photoelectron line, charge-corrected) eV</text><text transform="translate(12 ${(mt+H-mb)/2}) rotate(-90)" text-anchor="middle" class="al">KE (Auger line) eV</text></svg>`;
  return `<div class="wagwrap">${g}</div><div class="muted xs">● Samples in this project (highlighted = current) · □ Reference point · dashed = constant α′. Equal x/y scale, so α′ lines run at 45°.</div>`;
}

/* ============================================================
   서베이 패널
   ============================================================ */
function panelSurvey(el,sp){
  const isSv=sp.reg&&sp.reg.kind==="survey"; const st=S.ui.survey[sp.uid]||(S.ui.survey[sp.uid]={snr:4,tol:2,cands:null,acc:{}}); if(!(st.snr>0)) st.snr=4; if(!(st.tol>0)) st.tol=2;
  const L=LIB(); const libOK=L.core&&Object.keys(L.core).length>0;
  el.innerHTML=`
   <div class="xps-h">Element Identification ${isSv?"":`<span class="badge wr">Not a survey spectrum</span>`}</div>
   ${libOK?"":`<div class="note warn">원소 라이브러리(XPS_LIB)가 비어 있습니다.</div>`}
   <div class="xps-row"><span class="xps-f">SNR ≥ <input type="number" class="w56" id="xsSnr" value="${st.snr}"></span><span class="xps-f">Tolerance ± <input type="number" class="w56" step="0.1" id="xsTol" value="${st.tol}"> eV</span><button class="btn sm primary" id="xsRun">Find peaks &amp; identify</button></div>
   <div class="muted xs">대전 보정(Δ)을 적용한 BE 로 라이브러리(Avantage 원소표 위치 — 없으면 X-ray Data Booklet, Scofield 단면적 가중, Auger KE)와 대조합니다. 가장 큰 미설명 봉우리부터 차례로, 그 봉우리를 ‘센 선’으로 설명하면서 다른 선들도 같은 화학 이동 δ(−4…+6 eV)로 맞는 원소를 채택합니다. 합성 survey 150 개 시험: 검출 가능한 원소 97.2 % 채택, 없는 원소 오채택 0.05 개/스펙트럼 — 후보일 뿐이니 영역 스펙트럼으로 확인하세요. 방사성 원소·H·He 는 기본 제외.</div>
   <div id="xsOut">${st.cands?surveyHTML(sp,st):""}</div>`;
  // 비우거나 잘못 넣으면 기본값(SNR 4 · ±2 eV)으로 — 화면도 그 값으로 되돌림. 상태와 함께 저장
  $q(el,"#xsSnr").onchange=e=>{ const v=+e.target.value; st.snr=e.target.value!==""&&v>0?v:4; e.target.value=st.snr; save(); };
  $q(el,"#xsTol").onchange=e=>{ const v=+e.target.value; st.tol=e.target.value!==""&&v>0?v:2; e.target.value=st.tol; save(); };
  $q(el,"#xsRun").onclick=()=>{ const pk=C.findPeaks(sp.be,sp.y,{width:Math.max(1.2,sp.step*2),snr:st.snr,merge:Math.max(1.2,sp.step*1.5)}); st.peaks=pk; st.cands=C.identify(pk,LIB(),{hv:sp.hv,shift:shOf(sp),tol:st.tol,E:sp.be,y:sp.y,width:Math.max(1.2,sp.step*2)});
    st.acc={}; st.cands.slice(0,12).forEach(c=>{ if(C.acceptDefault(c)) st.acc[c.el]=true; }); updMarks(sp,st); save(); renderPanel(); drawPlot(); };
  bindSurvey(el,sp,st);
}
function surveyHTML(sp,st){
  const top=st.cands.slice(0,25);
  const q=surveyQuant(sp,st);
  return `<div class="muted xs" style="margin:4px 0">Peaks ${st.peaks.length} · accepted ${st.cands.filter(c=>c.accepted).length} · candidates ${st.cands.length}${st.cands.unassigned&&st.cands.unassigned.length?` · unassigned peaks ${st.cands.unassigned.length} (${st.cands.unassigned.slice(0,8).map(v=>f(v,0)).join(", ")}${st.cands.unassigned.length>8?" …":""} eV — 위성·손실·Auger 일 수 있음)`:""}</div>
   <div class="tw"><table class="xps-tbl pid-tbl"><thead><tr><th></th><th>Element</th><th class="num">Score</th><th class="num" title="Chemical shift relative to the library position">δ (eV)</th><th>Matched lines (observed BE)</th><th>Missing</th></tr></thead><tbody>
   ${top.map(c=>`<tr class="${c.overlap?"muted":""}"><td><input type="checkbox" data-acc="${esc(c.el)}" ${st.acc[c.el]?"checked":""}></td><td><b>${esc(c.el)}</b>${c.accepted?` <span class="badge ok" title="Accepted by the assignment that explains the largest unexplained peaks first">Accepted</span>`:c.possible?` <span class="badge inf" title="High score, but the main line coincides with a peak of another element — check with a region spectrum">Possible · overlap</span>`:c.overlap?` <span class="badge wr" title="No peak is explained by this element alone">Overlap</span>`:""}</td><td class="num">${f(c.score,2)}</td><td class="num">${c.shift!=null?(c.shift>0?"+":"")+f(c.shift,1):"—"}</td><td>${c.found.filter(x=>!x.auger).slice(0,4).map(x=>`${esc(x.line)} <span class="muted">${f(x.obs,1)}</span>`).join(", ")}${c.found.some(x=>x.auger)?` <span class="muted xs">+ Auger ${c.found.filter(x=>x.auger).length}</span>`:""}</td><td class="muted">${esc(c.missing.slice(0,3).join(", "))}</td></tr>`).join("")}
   </tbody></table></div>
   ${q?`<div class="xps-h">Survey Semi-quantification (main line)</div><div class="tw"><table class="xps-tbl"><thead><tr><th>Element · line</th><th class="num">Window (eV)</th><th class="num">Area</th><th class="num">S</th><th class="num">at %</th></tr></thead><tbody>${q.map(r=>`<tr><td>${esc(r.key)}</td><td class="num">${f(r.lo,1)}–${f(r.hi,1)}</td><td class="num">${g(r.area,4)}</td><td class="num">${g(r.S,4)}</td><td class="num"><b>${f(r.at,1)}</b></td></tr>`).join("")}</tbody></table></div><div class="muted xs">주선 둘레를 Shirley 로 빼고 적분, Quantification 탭의 T·λ 설정을 씀. 겹치는 선·위성은 고려하지 않으므로 대략값 — 정밀 정량은 영역 스펙트럼으로.</div>`:""}`;
}
function bindSurvey(el,sp,st){ $a(el,"[data-acc]").forEach(i=>i.onchange=()=>{ st.acc[i.dataset.acc]=i.checked; updMarks(sp,st); save(); renderPanel(); drawPlot(); }); }
function updMarks(sp,st){ const sh=shOf(sp); st.marks=[]; (st.cands||[]).filter(c=>st.acc[c.el]).forEach(c=>c.found.filter(x=>!x.auger).slice(0,3).forEach(x=>{ const i=sp.be.findIndex(b=>b>=x.obs-sh); st.marks.push({be:x.obs-sh,y:i>=0?sp.y[i]:0,lab:x.line}); })); }
function surveyQuant(sp,st){
  const L=LIB(); const acc=(st.cands||[]).filter(c=>st.acc[c.el]); if(!acc.length) return null; const sh=shOf(sp); const sc=L.scofield[srcOf(sp.hv)||"al"]||{};
  const rows=[];
  acc.forEach(c=>{ // 주선: Scofield 가장 큰 선 중 관측된 것 → 궤도 전체(두 subshell) 창
    const lines=c.found.filter(x=>!x.auger).map(x=>({...x,s:sc[x.line]||0})).sort((a,b)=>b.s-a.s); if(!lines.length) return;
    const m=lines[0]; const sub=m.line.split(" ")[1]; const orb=sub.replace(/(1\/2|3\/2|5\/2|7\/2)$/,""); const key=c.el+" "+orb;
    const lv=(L.core[c.el]||{}).lv||{}; const sk=Object.keys(lv).filter(k=>k.startsWith(orb)); const subs=sk.map(k=>lv[k]); if(!subs.length) return;
    // 관측 BE 와 같은 기준의 Booklet 값(subshell 이면 그 값, 궤도 전체면 Scofield 가중 평균) → 창 위치
    const ref=lv[sub]!=null?lv[sub]:sk.reduce((t,k)=>t+lv[k]*(sc[c.el+" "+k]||1),0)/sk.reduce((t,k)=>t+(sc[c.el+" "+k]||1),0); const off=m.obs-ref;
    const lo=Math.min(...subs)+off-4, hi=Math.max(...subs)+off+6;
    const R=C.sliceRegion(sp.be,sp.y,lo-sh,hi-sh); if(!R||R.E.length<5) return;
    const E=Float64Array.from(R.E), y=Float64Array.from(R.y); const B=C.bgShirley(E,y,2,true); const Q=C.cumtrapz(E,Float64Array.from(y,(v,i)=>v-B[i]));
    const S_=sc[key]!=null?sc[key]:subs.reduce((s,_,i)=>s,0)||sc[m.line];
    rows.push({key,lo,hi,area:Q[Q.length-1],S:S_,BE:m.obs-sh});
  });
  const Q=C.quantify(rows.filter(r=>r.S>0&&r.area>0).map(r=>({area:r.area,S:r.S,BE:r.BE,hv:sp.hv,PE:sp.meta&&sp.meta.passEnergy,tf:tfFor(sp)})),{imfp:S.quant.imfp});
  let k=0; return rows.filter(r=>r.S>0&&r.area>0).map(r=>Object.assign(r,{at:Q[k++].at}));
}

/* ============================================================
   정보 패널
   ============================================================ */
function panelInfo(el,sp){
  const m=sp.meta||{}; const nz=C.noiseEstimate(sp.be,sp.y);
  const kv=(k,v)=>v==null||v===""?"":`<tr><td class="muted">${esc(k)}</td><td>${v}</td></tr>`;
  const cExp=sp.c?1/sp.c:null;
  const noiseTxt=nz?`κ = ${g(nz.kappa,3)} CPS${sp.c?` (Poisson: 1/c = ${g(cExp,3)} → ratio ${f(nz.kappa/cExp,2)})`:` → effective counting time ≈ ${g(1/nz.kappa,3)} s/pt`}`:"—";
  const L=LIB();
  const orbs=["1s","2s","2p","3s","3p","3d","4s","4p","4d","4f","5s","5p","5d","6s","6p"];
  const regEd=`<div class="xps-h">Element · Orbital ${sp.reg&&sp.reg.el?`<span class="badge ok">${esc(sp.reg.el+" "+(sp.reg.orb||""))}</span>`:`<span class="badge wr">Unknown</span>`}</div>
    <div class="xps-row"><span class="xps-f">Element <input type="text" class="xps-in w56" id="xiEl" value="${esc(sp.reg&&sp.reg.el||"")}" placeholder="e.g. Fe"></span><span class="xps-f">Orbital <select class="xps-sel" id="xiOrb">${orbs.map(o=>`<option ${sp.reg&&sp.reg.orb===o?"selected":""}>${o}</option>`).join("")}</select></span><button class="btn xs" id="xiSet">Apply</button></div>
    <div class="muted xs">파일 이름에서 읽지 못했거나 틀렸을 때 지정합니다. doublet Δ·넓이 비, Reference lines, RSF 가 이 값을 씁니다.</div>`;
  const smpEd=`<div class="xps-h">Sample Grouping</div>
    <div class="xps-row"><span class="xps-f">Sample of this spectrum <input type="text" class="xps-in w120" id="xiSmp" list="xiSmpL" value="${esc(sp.sample||"")}" placeholder="${NOSMP_L}" style="width:220px"></span><datalist id="xiSmpL">${[...new Set(S.specs.map(SK))].map(k=>`<option value="${esc(k)}">`).join("")}</datalist><button class="btn xs" id="xiSmpOk">Move</button></div>
    <div class="muted xs">대전 보정 Δ 는 시료마다 하나입니다. 자동 묶음이 틀렸으면(예: C 1s 가 다른 묶음에 들어감) 같은 시료 이름으로 옮기세요. 목록의 ✎ 로 묶음 전체 이름을 바꿀 수도 있습니다.</div>`;
  const ax=(m.axis)||"BE", ac=axisCheck(sp);
  const axEd=isText(sp)?`<div class="xps-h">Energy Axis</div>
    <div class="xps-row"><span class="xps-f" title="How the first column of the text file is read">Column 1 <select class="xps-sel" id="xiAx"><option value="BE" ${ax==="BE"?"selected":""}>Binding energy</option><option value="KE" ${ax==="KE"?"selected":""}>Kinetic energy</option></select></span><span class="xps-f" title="X-ray energy — used for KE ↔ BE and for quantification">hν <input type="number" step="0.01" class="w64" id="xiHv" value="${+(+sp.hv).toFixed(2)}"> eV</span><button class="btn xs" id="xiAxOk" title="Re-read the file with this axis and hν (fit of this spectrum is reset)">Re-read</button></div>
    ${ac?`<div class="note warn" style="margin-top:4px">${esc(ac.msg)}${ac.fits?` <button class="btn xs" id="xiAxS">Read as ${ac.other}</button>`:""}</div>`:""}`:(ac?`<div class="note warn" style="margin-top:6px">${esc(ac.msg)}</div>`:"");
  el.innerHTML=regEd+smpEd+axEd+`<div class="xps-h">File · Source</div><div class="tw"><table class="xps-tbl">
    ${kv("File",esc(sp.file))}${kv("Format",esc(FMTN(sp.format)))}${kv("Region",esc(sp.region))}${kv("Sample",`<span title="${esc(SK(sp))}">${esc(smpDisp(SK(sp)).name)}</span> <span class="muted xs">${esc(smpDisp(SK(sp)).ctx)}</span>`)}
    ${kv("Points · step",`${sp.n} · ${f(sp.step,4)} eV`)}${kv("BE range",`${f(sp.be[0],2)} – ${f(sp.be[sp.n-1],2)} eV (measured)`)}
    ${kv("X-ray energy (hν)",`${f(sp.hv,2)} eV <span class="muted xs">${esc(sp.hvSrc)}</span>`)}
    ${kv("Count conversion c",sp.c?`${g(sp.c,4)} s <span class="muted xs">${esc(sp.cSrc)}</span>`:`<span class="muted">모름 — 잡음 척도는 잔차에서 추정</span>`)}
    ${kv("Noise (variance/intensity)",noiseTxt)}
    ${kv("pass energy",m.passEnergy!=null?m.passEnergy+" eV":"")}${kv("dwell · scans",m.dwell!=null?`${m.dwell} s × ${m.scans}`:"")}${kv("lens mode",esc(m.lensMode||""))}${kv("spot",m.spot&&m.spot[0]?`${m.spot[0]} × ${m.spot[1]} µm`:"")}
    ${kv("flood gun",esc(m.floodGun||""))}${kv("Analyzer work function",m.workFunction!=null?m.workFunction+" eV":"")}${kv("Instrument",esc(m.instrument||""))}${kv("Avantage",esc(m.software||""))}
    ${kv("Acquired",esc(m.start?m.start.replace("T"," ").slice(0,19)+" UTC":(m.when||"")))}${kv("Original path",`<span class="mono">${esc(m.origPath||"")}</span>`)}
    ${kv("transmission",m.tf?`${m.tf.coef.map(v=>g(v,6)).join(", ")} <div class="muted xs">${esc(m.tf.src)}</div>`:"")}
    </table></div>
    ${sp.warn&&sp.warn.length?`<div class="note warn" style="margin-top:6px">${sp.warn.map(esc).join("<br>")}</div>`:""}
    ${sp.extra&&sp.extra.avantageCols?`<div class="note" style="margin-top:6px">Avantage 내보내기의 추가 열: ${esc(Object.keys(sp.extra.avantageCols).join(", "))}</div>`:""}
    <div class="xps-h">Library Sources</div><div class="muted xs" style="line-height:1.6">${Object.entries((L.meta&&L.meta.sources)||{}).map(([k,v])=>`<b>${esc(k)}</b>: ${esc(v)}`).join("<br>")||"—"}</div>
    <div class="xps-h">Shortcuts</div><div class="muted xs" style="line-height:1.6">그래프: 휠 = 확대/축소 · 끌기 = 구간 확대 · 두 번 누르기 = 전체 · 점선 손잡이 = 피팅 영역 · 색 동그라미 끌기 = 성분 위치·높이 · 성분 면적 누르기 = 선택</div>`;
  const xm=$q(el,"#xiSmpOk"); if(xm) xm.onclick=()=>moveSpec(sp,$q(el,"#xiSmp").value);
  const xa=$q(el,"#xiAxOk"); if(xa) xa.onclick=()=>{ const hv=+$q(el,"#xiHv").value; if(!(hv>100)){ note("hν 를 확인하세요 (eV)","warn"); return; } applyAxis(sp,$q(el,"#xiAx").value,hv); };
  const xas=$q(el,"#xiAxS"); if(xas&&ac) xas.onclick=()=>applyAxis(sp,ac.other,sp.hv);
  const xs=$q(el,"#xiSet"); if(xs) xs.onclick=()=>{ const raw=$q(el,"#xiEl").value.trim(); const e=IO.SYMS.find(s=>s.toLowerCase()===raw.toLowerCase()); if(!e){ note("원소 기호를 확인하세요 (예: Fe, Co, Ta)","warn"); return; } const o=$q(el,"#xiOrb").value; if(!(LIB().core[e]||{lv:{}}).lv||!Object.keys((LIB().core[e]||{}).lv||{}).some(k=>k.startsWith(o))) note(e+" 에는 라이브러리에 "+o+" 준위가 없습니다 — 그래도 지정합니다","warn");
    sp.reg={kind:"region",el:e,orb:o,j:null,label:e+" "+o,manual:true}; sp.region=sp.reg.label; sp.warn=(sp.warn||[]).filter(w=>!/라이브러리\)가 없습니다/.test(w)); save(); rerender(); note("원소·궤도를 "+sp.reg.label+" 로 지정했습니다","ok"); };
}

/* ============================================================
   내보내기 · 프로젝트
   ============================================================ */
function download(name,blob){ const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },500); }
const safe=s=>String(s||"xps").replace(/[\\/:*?"<>|\s]+/g,"_").slice(0,80);
const csvE=v=>{ const s=String(v==null?"":v); return /[",\n]/.test(s)?`"${s.replace(/"/g,'""')}"`:s; };
function exportCurves(sp){ const p=pr(sp), M=curModel(sp,p), sh=shOf(sp); if(!M) return;
  const head=["BE_measured_eV","BE_corrected_eV","data_CPS","background","envelope","residual_over_sigma",...p.comps.map((c,i)=>(c.name||"comp"+(i+1)))];
  const lines=[`# INNO Analytics XPS · ${sp.sample} · ${sp.region} · file ${sp.file} · bg ${p.bg.type} · charge shift ${sh} eV (${smpOf(sp).refNote||(sh?"":"not charge-corrected")})`,head.map(csvE).join(",")];
  for(let i=0;i<M.n;i++) lines.push([M.E[i].toFixed(4),(M.E[i]+sh).toFixed(4),M.y[i],M.B[i].toFixed(4),M.env[i].toFixed(4),M.resid[i].toFixed(4),...M.parts.map(a=>a[i].toFixed(4))].join(","));
  download(safe(sp.sample+"_"+sp.region)+"_curves.csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"})); }
function exportPeaks(sp){ const k=SK(sp); const lines=["sample,region,component,shape,BE_corr_eV,sBE,FWHM_eV,sFWHM,height_CPS,area_CPSeV,sArea,area_pct,sArea_pct,doublet_dE,doublet_ratio,chi2nu,uncertainty"];
  sampleSpecs(k).forEach(s=>{ const p=pr(s), r=p.result; if(!r) return; const sh=shOf(s); const mc=r.mc&&r.mc.summ;
    r.comps.forEach((c,i)=>{ const q=r.summ[i], m=mc?mc[i]:null; lines.push([s.sample,s.region,c.name,c.shape,(q.be+sh).toFixed(3),q.sig.be===null?"":f(m?m.be.sd:q.sig.be,3),f(q.fwhm,3),q.sig.fwhm===null?"":f(m?m.fwhm.sd:q.sig.fwhm,3),g(q.height,6),q.area.toFixed(2),q.sig.area===null?"":f(m?m.area.sd:q.sig.area,2),f(q.areaPct,2),f(m?m.areaPct.sd:q.sig.areaPct,2),c.doublet?c.doublet.dE:"",c.doublet?c.doublet.r:"",f(r.stats.chi2nu,4),m?"MonteCarlo":"covariance"].map(csvE).join(",")); }); });
  { const sm=S.samples[k]||{}; lines.push(`# charge shift ${sm.shift||0} eV · reference: ${sm.refNote||(sm.shift?"manual":"not charge-corrected (BE as measured)")}`); }
  download(safe(k)+"_peaks.csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"})); }
function exportQuant(k,rows,qi,qmi){ const lines=["sample,element,line,region,BE_corr_eV,area,sArea,S,S_source,T,lambda,at_pct,s_at_pct,at_pct_mc95_lo,at_pct_mc95_hi,included"];
  rows.forEach(r=>{ const q=qi.get(r.id), m=qmi&&qmi.get(r.id); lines.push([k,r.el,r.line,r.sp.region,(r.BE+r.sh).toFixed(2),r.area,r.sArea||"",r.S,r.Ssrc,q?q.T:"",q?q.L:"",q&&q.at!=null?q.at.toFixed(3):(q&&q.nd?"n.d.":""),q&&q.sAt?q.sAt.toFixed(3):"",m&&q&&q.sAt?m.lo.toFixed(3):"",m&&q&&q.sAt?m.hi.toFixed(3):"",r.inc].map(csvE).join(",")); });
  lines.push(`# S: ${S.quant.lib}; T: ${S.quant.tf.model}; IMFP: ${JSON.stringify(S.quant.imfp)}`);
  { const sm=S.samples[k]||{}; lines.push(`# charge shift ${sm.shift||0} eV · reference: ${sm.refNote||(sm.shift?"manual":"not charge-corrected (BE as measured)")}`); }
  download(safe(k)+"_quant.csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"})); }
function figSVG(){ const svg=ROOT&&$q(ROOT,"#xpsSvg"); if(!svg) return null; const cl=svg.cloneNode(true); const cs=getComputedStyle(ROOT.querySelector(".xps"));
  // 변수·클래스를 값으로 굳힌다
  const vars=["--s1","--s2","--s3","--s4","--s5","--s6","--s7","--s8","--dat","--env","--bgl","--grid","--line2","--muted","--ink2","--ink","--surface","--page","--accent","--info","--faint","--chip"];
  const map={}; vars.forEach(v=>map[v]=cs.getPropertyValue(v).trim());
  const fix=s=>s.replace(/var\((--[a-z0-9-]+)\)/g,(m,v)=>map[v]||"#888");
  const st={".ax":"stroke:var(--line2);stroke-width:1;fill:none",".gd":"stroke:var(--grid);stroke-width:1","text":"fill:var(--muted);font-size:11px;font-family:DM Sans,Noto Sans KR,Arial,sans-serif","text.al":"fill:var(--ink2);font-size:12px",".dat":"fill:none;stroke:var(--dat);stroke-width:1.2",".datp":"fill:var(--dat);opacity:.75",".env":"fill:none;stroke:var(--env);stroke-width:2",".bgl":"fill:none;stroke:var(--bgl);stroke-width:1.5;stroke-dasharray:5 4",".veil":"fill:var(--page);opacity:.62",".hdl":"display:none",".hdlhit":"display:none",".mk":"display:none",".ref":"stroke:var(--info);stroke-width:1;opacity:.7","text.ref":"fill:var(--info);font-size:10.5px","text.pid":"fill:var(--ink);font-size:11px;font-weight:600","text.cl":"fill:var(--ink2);font-size:11px;font-weight:600",".rz":"stroke:var(--line2);stroke-dasharray:3 3"};
  const css=Object.entries(st).map(([k,v])=>`${k}{${fix(v)}}`).join("\n");
  cl.querySelectorAll("[fill],[stroke]").forEach(n=>{ ["fill","stroke"].forEach(a=>{ const v=n.getAttribute(a); if(v&&v.includes("var(")) n.setAttribute(a,fix(v)); }); });
  const o=cl.querySelector("#xpsOver"); if(o) o.remove();
  cl.insertAdjacentHTML("afterbegin",`<style>${css}</style><rect width="100%" height="100%" fill="${map["--surface"]||"#fff"}"/>`);
  cl.setAttribute("xmlns","http://www.w3.org/2000/svg"); cl.setAttribute("width",PW); cl.setAttribute("height",PH+RH);
  return new XMLSerializer().serializeToString(cl); }
function exportFig(kind){ const sp=cur(); const s=figSVG(); if(!s){ note("그림이 없습니다","warn"); return; } const nm=safe(sp.sample+"_"+sp.region);
  if(kind==="svg"){ download(nm+".svg",new Blob([s],{type:"image/svg+xml"})); return; }
  const img=new Image(); img.onload=()=>{ const k=3; const cv=document.createElement("canvas"); cv.width=PW*k; cv.height=(PH+RH)*k; const ctx=cv.getContext("2d"); ctx.scale(k,k); ctx.drawImage(img,0,0); cv.toBlob(b=>download(nm+".png",b),"image/png"); }; img.src="data:image/svg+xml;charset=utf-8,"+encodeURIComponent(s); }
function exportProject(){ if(!S.specs.length){ note("저장할 것이 없습니다","warn"); return; } const o=JSON.stringify({app:"INNO Analytics XPS",v:1,saved:new Date().toISOString(),state:S},(k,v)=>ArrayBuffer.isView(v)?{__f64:Array.from(v)}:v);
  download("XPS_project_"+new Date().toISOString().slice(0,10)+".json",new Blob([o],{type:"application/json"})); }
function importProject(txt){ try{ const o=revive(JSON.parse(txt)); if(!o||o.app!=="INNO Analytics XPS"||!o.state) throw new Error("INNO Analytics XPS 프로젝트 파일이 아닙니다"); S=Object.assign(defState(),o.state); S.ui=Object.assign(defState().ui,o.state.ui||{}); S.quant=Object.assign(defState().quant,o.state.quant||{}); PV={}; save(); rerender(); note("프로젝트를 열었습니다 ("+S.specs.length+" 스펙트럼)","ok"); }catch(e){ note("열기 실패: "+(e.message||e),"warn"); } }


/* ============================================================
   v1.8 — Series · ARXPS · Imaging · Processing (SPEC_multi.md · SPEC_proc.md)
   계산: XPS_MULTI · XPS_ARXPS · XPS_IMAGE · XPS_PROC (DOM 없음). 긴 계산(PCA·MCR·MaxEnt·영상 NNLS·디컨볼루션)은 Worker 의 xj 작업.
   상태: 설정만 S.ser · S.arx · S.imu · S.prc (IndexedDB "cur"), 영상 stack 은 IndexedDB "img"(구조 복제), 결과는 RES(메모리).
   ============================================================ */
const MU=()=>window.XPS_MULTI, ARX_=()=>window.XPS_ARXPS, IMG=()=>window.XPS_IMAGE, PRC=()=>window.XPS_PROC;
const defs=(o,d)=>{ Object.keys(d).forEach(k=>{ if(o[k]===undefined) o[k]=d[k]; }); return o; };
const RES={ser:null,mv:null,arx:null,img:null,pca:null,roi:null};
const PRES={};
const numv=(v,d)=>{ const x=parseFloat(v); return isFinite(x)?x:d; };
const trapzA=(E,y)=>{ let s=0; for(let i=1;i<E.length;i++) s+=0.5*(y[i]+y[i-1])*(E[i]-E[i-1]); return s; };
const unitArea=(E,y)=>{ const a=trapzA(E,y); return Float64Array.from(y,v=>a>0?v/a:v); };
const sliceWin=(E,lo,hi)=>{ let i0=0,i1=E.length-1; if(lo!=null) while(i0<E.length&&E[i0]<lo) i0++; if(hi!=null) while(i1>=0&&E[i1]>hi) i1--; return [i0,i1]; };

/* ---------- Worker 작업 (이 함수의 글이 Worker 안에도 들어간다 — 바깥 변수를 쓰지 말 것) ---------- */
function xjobs(G){
  const M=G.XPS_MULTI, AR=G.XPS_ARXPS, IM=G.XPS_IMAGE, P=G.XPS_PROC;
  function winStack(st,lo,hi){ const r=IM.idxRange(st.E,lo,hi), i0=r[0], i1=r[1], N=i1-i0+1, np=st.nx*st.ny, out=new Float32Array(np*N);
    for(let p=0;p<np;p++) out.set(st.data.subarray(p*st.nE+i0,p*st.nE+i1+1),p*N); return Object.assign({},st,{nE:N,E:st.E.slice(i0,i1+1),data:out}); }
  return {
    serPCA(a,prog){ prog(0.05,"PCA");
      const P1=M.pca(a.raw,{scale:a.scale||"poisson",center:"none"});
      const nMax=Math.max(1,Math.min(a.nMax||6,Math.min(a.raw.length,a.raw[0].length)-1));
      prog(0.25,"noise test (MLPCA)"); const NF=M.noiseFactors(a.raw,P1,a.cfac||1,{method:"mlpca",nMax,full:true});
      prog(0.9,"target-test PCA"); const P2=M.pca(a.proc,{scale:"none",center:"none"});
      return {eig:P1.eig,varPct:P1.varPct,mal:P1.malinowski,warn:P1.warn,NF,P2}; },
    serMCR(a){ return M.mcrAls(a.cols,a.n,{init:a.init,closure:a.closure,maxIter:a.maxIter||3000}); },
    serRef(a){ return M.fitReferences(a.cols,a.refs,a.opt); },
    arxMaxent(a){ const r=AR.maxentProfile(a.th,a.L,a.I,a.sI,a.opt); delete r.setup; return r; },
    imgArea(a){ return IM.areaMap(a.st,a.opt); },
    imgAreas(a,prog){ return a.items.map((it,i)=>{ prog(i/a.items.length,it.name||""); return IM.areaMap(it.st,it.opt); }); },
    imgPCA(a,prog){ const st=a.win?winStack(a.st,a.win[0],a.win[1]):a.st; prog(0.05,"PCA"); const P=IM.pcaStack(st,{scale:a.scale,center:"none",nKeep:a.nKeep});
      prog(0.5,"noise test (MLPCA subset)"); let NF=null; try{ NF=IM.noiseFactorsStack(st,P,{method:"mlpca",nMax:Math.min(a.nMax||8,st.nE-1)}); }catch(e){ NF={error:String(e.message||e)}; }
      return {P,NF,E:st.E}; },
    imgRecon(a){ const st=a.win?winStack(a.st,a.win[0],a.win[1]):a.st; return IM.pcaReconstruct(st,a.P,a.n); },
    imgNNLS(a){ const sb=a.bgsub?IM.subtractBackgroundStack(a.st,{lo:a.lo,hi:a.hi,bg:a.bgsub,n:3}):winStack(a.st,a.lo,a.hi);
      return IM.nnlsMap(sb,a.refs,{bg:a.bg||0,weight:a.weight}); },
    deconv(a){ const r=P.deconvolve(a.E,a.y,a.opt); delete r.chi2Hist; delete r.resid; return r; }
  };
}
let _XJM=null; const XJM=()=>_XJM||(_XJM=xjobs(window));
/* 하나의 Worker 로 xj 작업: 진행(작업이 알려 주면 비율, 아니면 흐르는 막대 + 경과 시간) · Stop = Worker 종료 */
function xrun(job,args,o){
  o=o||{}; if(RUN) return Promise.reject(new Error("이미 계산 중입니다"));
  const url=workerURL(); const label=o.label||job;
  if(!url){ try{ return Promise.resolve(XJM()[job](args,()=>{})); }catch(e){ return Promise.reject(e); } }
  return new Promise((resolve,reject)=>{
    const w=new Worker(url), t0=performance.now(); let fr=null, txt="";
    const show=()=>setProg(`${label}${txt?" — "+txt:""} · ${((performance.now()-t0)/1000).toFixed(1)} s`,fr);
    const tick=setInterval(show,250); show();
    const end=()=>{ clearInterval(tick); try{ w.terminate(); }catch(e){} RUN=null; hideProg(); };
    RUN={ws:[w],cancel:()=>{ end(); reject(new Error("취소")); }};
    w.onmessage=e=>{ const d=e.data; if(d.prog!==undefined){ fr=d.prog; txt=d.txt||""; return; } end(); if(d.ok) resolve(d.out); else reject(new Error(d.err)); };
    w.onerror=e=>{ e.preventDefault&&e.preventDefault(); end(); reject(new Error(e.message||"Worker 오류")); };
    w.postMessage({id:1,cmd:"xj",job,args},o.transfer||[]);
  });
}
/* 여러 Worker 에 나눈 같은 작업(영상 NNLS 화소 묶음): 진행 = 끝난 묶음 수 */
function xpool(job,argsList,o){
  o=o||{}; if(RUN) return Promise.reject(new Error("이미 계산 중입니다"));
  const url=workerURL(); const label=o.label||job, n=argsList.length;
  if(!url){ try{ return Promise.resolve(argsList.map(a=>XJM()[job](a,()=>{}))); }catch(e){ return Promise.reject(e); } }
  const K=Math.min(nWorkers(),n);
  return new Promise((resolve,reject)=>{
    const ws=[...Array(K)].map(()=>new Worker(url)), t0=performance.now(), outs=new Array(n); let next=0, done=0, dead=false;
    const show=()=>setProg(`${label} ${done}/${n} · ${((performance.now()-t0)/1000).toFixed(1)} s`,done/n);
    const tick=setInterval(show,250); show();
    const end=()=>{ dead=true; clearInterval(tick); ws.forEach(w=>{ try{ w.terminate(); }catch(e){} }); RUN=null; hideProg(); };
    RUN={ws,cancel:()=>{ end(); reject(new Error("취소")); }};
    const feed=w=>{ if(next>=n) return; const id=next++; w.postMessage({id,cmd:"xj",job,args:argsList[id]}); };
    ws.forEach(w=>{ w.onmessage=e=>{ const d=e.data; if(d.prog!==undefined||dead) return; if(!d.ok){ end(); reject(new Error(d.err)); return; } outs[d.id]=d.out; done++; if(done===n){ end(); resolve(outs); } else feed(w); };
      w.onerror=e=>{ e.preventDefault&&e.preventDefault(); end(); reject(new Error(e.message||"Worker 오류")); }; feed(w); });
  });
}
const failNote=(what,e)=>{ const m=String(e&&e.message||e); if(m==="취소") return; note(what+": "+m,"warn"); };

/* ---------- 작은 HTML 조각 ---------- */
const nIn=(k,v,o)=>{ o=o||{}; return `<input type="number"${k?` data-st="${esc(k)}"`:""} value="${v==null||v===""||!isFinite(v)?"":esc(v)}" class="${o.cls||"w64"}"${o.step?` step="${o.step}"`:` step="any"`}${o.min!=null?` min="${o.min}"`:""}${o.max!=null?` max="${o.max}"`:""}${o.ph?` placeholder="${esc(o.ph)}"`:""}${o.title?` title="${esc(o.title)}"`:""}${o.id?` id="${o.id}"`:""}>`; };
const sIn=(k,cur,opts,o)=>{ o=o||{}; return `<select class="xps-sel" data-st="${esc(k)}"${o.title?` title="${esc(o.title)}"`:""}${o.id?` id="${o.id}"`:""}${o.dis?" disabled":""}>${opts.map(([v,t])=>`<option value="${esc(v)}" ${String(cur)===String(v)?"selected":""}>${esc(t)}</option>`).join("")}</select>`; };
const cIn=(k,on,lab,o)=>{ o=o||{}; return `<label class="xps-f"${o.title?` title="${esc(o.title)}"`:""}><input type="checkbox" data-st="${esc(k)}" ${on?"checked":""}${o.id?` id="${o.id}"`:""}> ${lab}</label>`; };
/* data-st="obj|key|key" → S 의 그 자리. 숫자 칸은 숫자(빈칸 = null), 체크는 참/거짓 */
function bindSt(el,after){ $a(el,"[data-st]").forEach(n=>{ n.onchange=()=>{ const path=n.dataset.st.split("|"); let o=S; for(let i=0;i<path.length-1;i++){ if(o[path[i]]==null||typeof o[path[i]]!=="object") o[path[i]]={}; o=o[path[i]]; }
  const k=path[path.length-1]; o[k]=n.type==="checkbox"?n.checked:n.type==="number"?(n.value===""?null:+n.value):n.value; save(); if(after) after(k,n); }; }); }

/* ---------- 그래프(일반) ---------- */
const CW=640, CH=300;
function chartSVG(o){
  const W=o.w||CW, H=o.h||CH, fs=o.fs||1; const ml=Math.round((o.ml||60)*Math.min(fs,1.5)), mr=14, mt=12, mb=Math.round(36*fs);
  const pw=W-ml-mr, ph=H-mt-mb, ylog=!!(o.y&&o.y.log);
  const SER=(o.series||[]).filter(s=>s&&s.x&&s.x.length);
  const ty=v=>ylog?Math.log10(Math.max(v,1e-300)):v;
  let x0=o.x&&o.x.min, x1=o.x&&o.x.max, y0=o.y&&o.y.min, y1=o.y&&o.y.max;
  if(x0==null||x1==null){ let a=Infinity,b=-Infinity; SER.forEach(s=>{ for(let i=0;i<s.x.length;i++){ const v=s.x[i]; if(isFinite(v)){ a=Math.min(a,v); b=Math.max(b,v); } } }); if(x0==null) x0=a; if(x1==null) x1=b; }
  if(y0==null||y1==null){ let a=Infinity,b=-Infinity; SER.forEach(s=>{ for(let i=0;i<s.y.length;i++){ let v=s.y[i]; if(!isFinite(v)) continue; const e=s.err?s.err[i]:0; const lo=s.band?s.band.lo[i]:v-(isFinite(e)?e:0), hi=s.band?s.band.hi[i]:v+(isFinite(e)?e:0);
      if(ylog){ if(v>0){ a=Math.min(a,v); b=Math.max(b,v); } } else { a=Math.min(a,lo,v); b=Math.max(b,hi,v); } } }); if(!isFinite(a)){ a=0; b=1; }
    if(ylog){ y0=y0==null?Math.pow(10,Math.floor(Math.log10(a))):y0; y1=y1==null?Math.pow(10,Math.ceil(Math.log10(b))):y1; }
    else { const pad=(b-a)*0.06||Math.abs(b)*0.1||1; if(y0==null) y0=o.y&&o.y.zero&&a>=0?0:a-pad; if(y1==null) y1=b+pad; } }
  if(!(x1>x0)){ x0-=0.5; x1+=0.5; } if(!(ty(y1)>ty(y0))){ y1=ylog?y0*10:y0+1; }
  const rev=!!(o.x&&o.x.rev);
  const X=v=>ml+(rev?(x1-v):(v-x0))/(x1-x0)*pw, Y=v=>mt+(ty(y1)-ty(v))/(ty(y1)-ty(y0))*ph;
  const out=[`<svg class="xc-svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.title||"chart")}" style="--xfs:${fs}" xmlns="http://www.w3.org/2000/svg">`];
  const xt=o.xTicks||niceTicks(x0,x1,fs>1.3?4:7).map(v=>({v,t:fmtTick(v)}));
  let yt; if(ylog){ yt=[]; for(let e=Math.ceil(Math.log10(y0)-1e-9);e<=Math.floor(Math.log10(y1)+1e-9);e++) yt.push({v:Math.pow(10,e),t:"1e"+e}); } else yt=niceTicks(y0,y1,fs>1.3?4:5).map(v=>({v,t:fmtTick(v)}));
  xt.forEach(t=>out.push(`<line class="gd" x1="${X(t.v).toFixed(1)}" x2="${X(t.v).toFixed(1)}" y1="${mt}" y2="${mt+ph}"/><text x="${X(t.v).toFixed(1)}" y="${mt+ph+4+Math.round(11*fs)}" text-anchor="middle">${esc(t.t)}</text>`));
  yt.forEach(t=>out.push(`<line class="gd" x1="${ml}" x2="${ml+pw}" y1="${Y(t.v).toFixed(1)}" y2="${Y(t.v).toFixed(1)}"/><text x="${ml-5}" y="${(Y(t.v)+4).toFixed(1)}" text-anchor="end">${esc(t.t)}</text>`));
  out.push(`<rect class="ax" x="${ml}" y="${mt}" width="${pw}" height="${ph}" fill="none"/>`);
  out.push(`<clipPath id="${o.id||"xc"}Clip"><rect x="${ml}" y="${mt}" width="${pw}" height="${ph}"/></clipPath><g clip-path="url(#${o.id||"xc"}Clip)">`);
  (o.vlines||[]).forEach(v=>{ if(!(v.x>=Math.min(x0,x1)&&v.x<=Math.max(x0,x1))) return; out.push(`<line x1="${X(v.x).toFixed(1)}" x2="${X(v.x).toFixed(1)}" y1="${mt}" y2="${mt+ph}" stroke="${v.c||"var(--accent)"}" stroke-width="1.2" stroke-dasharray="4 3"/>${v.t?`<text x="${(X(v.x)+3).toFixed(1)}" y="${mt+12}" style="fill:${v.c||"var(--accent)"}">${esc(v.t)}</text>`:""}`); });
  (o.shade||[]).forEach(v=>{ const a=X(v.a), b=X(v.b); out.push(`<rect x="${Math.min(a,b).toFixed(1)}" y="${mt}" width="${Math.abs(b-a).toFixed(1)}" height="${ph}" fill="var(--faint)" opacity=".12"/>`); });
  SER.forEach(s=>{ const col=s.c||"var(--s1)";
    if(s.band){ let d="", d2=""; for(let i=0;i<s.x.length;i++){ const h=s.band.hi[i]; if(!isFinite(h)) continue; d+=(d?"L":"M")+X(s.x[i]).toFixed(1)+","+Y(h).toFixed(1); } for(let i=s.x.length-1;i>=0;i--){ const l=s.band.lo[i]; if(!isFinite(l)) continue; d2+="L"+X(s.x[i]).toFixed(1)+","+Y(ylog?Math.max(l,y0):l).toFixed(1); } if(d) out.push(`<path d="${d}${d2}Z" fill="${col}" fill-opacity=".18" stroke="none"/>`); }
    if(s.line!==false){ let d="", pen=false; for(let i=0;i<s.x.length;i++){ const v=s.y[i]; if(!isFinite(v)||(ylog&&!(v>0))){ pen=false; continue; } d+=(pen?"L":"M")+X(s.x[i]).toFixed(1)+","+Y(v).toFixed(1); pen=true; } if(d) out.push(`<path d="${d}" fill="none" stroke="${col}" stroke-width="${s.w||1.6}"${s.dash?` stroke-dasharray="${s.dash}"`:""} stroke-linejoin="round"/>`); }
    if(s.err) for(let i=0;i<s.x.length;i++){ const v=s.y[i], e=s.err[i]; if(!isFinite(v)||!isFinite(e)||!(e>0)) continue; const x=X(s.x[i]).toFixed(1); out.push(`<path d="M${x},${Y(v-e).toFixed(1)}V${Y(v+e).toFixed(1)}M${(+x-3).toFixed(1)},${Y(v-e).toFixed(1)}h6M${(+x-3).toFixed(1)},${Y(v+e).toFixed(1)}h6" stroke="${col}" stroke-width="1" fill="none"/>`); }
    if(s.pts) for(let i=0;i<s.x.length;i++){ const v=s.y[i]; if(!isFinite(v)||(ylog&&!(v>0))) continue; out.push(`<circle cx="${X(s.x[i]).toFixed(1)}" cy="${Y(v).toFixed(1)}" r="${s.r||3}" fill="${s.hollow?"var(--surface)":col}" stroke="${col}" stroke-width="1.2"/>`); }
    if(s.labels) s.labels.forEach((t,i)=>{ const v=s.y[i]; if(!isFinite(v)||!t) return; out.push(`<text x="${(X(s.x[i])+5).toFixed(1)}" y="${(Y(v)-5).toFixed(1)}" style="fill:${col}">${esc(t)}</text>`); });
  });
  out.push(`</g>`);
  if(o.x&&o.x.label) out.push(`<text class="al" x="${ml+pw/2}" y="${H-5}" text-anchor="middle">${esc(o.x.label)}</text>`);
  if(o.y&&o.y.label&&fs<=1.3) out.push(`<text class="al" x="13" y="${mt+ph/2}" transform="rotate(-90 13 ${mt+ph/2})" text-anchor="middle">${esc(o.y.label)}</text>`);
  out.push(`</svg>`);
  return out.join("");
}
function chartLegend(o){ const L=(o.series||[]).filter(s=>s&&s.name); if(!L.length) return ""; return `<div class="xps-legend">${L.map(s=>`<span><i style="background:${s.c||"var(--s1)"}${s.dash?";opacity:.6":""}"></i>${esc(s.name)}</span>`).join("")}</div>`; }
/* host 폭에 맞춘 글자 배율로 그린다(휴대폰에서 글자가 너무 작지 않게) */
function putChart(host,o){ if(!host) return; const w=host.clientWidth||CW; o.fs=w>=560?1:Math.max(1,Math.min(2.2,+(0.8*(o.w||CW)/w).toFixed(2))); host.innerHTML=chartSVG(o)+chartLegend(o); host._chart=o; }
/* 그림 내보내기: var(--…) 를 값으로 굳힌 SVG, 또는 3× PNG */
function svgFreeze(svg){ const cl=svg.cloneNode(true); const cs=getComputedStyle(ROOT.querySelector(".xps")); const memo={};
  const fix=s=>s.replace(/var\((--[a-z0-9-]+)\)/g,(m,v)=>memo[v]||(memo[v]=cs.getPropertyValue(v).trim()||"#888"));
  cl.querySelectorAll("*").forEach(n=>{ ["fill","stroke","style"].forEach(a=>{ const v=n.getAttribute(a); if(v&&v.includes("var(")) n.setAttribute(a,fix(v)); }); });
  const css=fix(`.gd{stroke:var(--grid);stroke-width:1}.ax{stroke:var(--line2);stroke-width:1;fill:none}text{fill:var(--muted);font-size:11px;font-family:DM Sans,Noto Sans KR,Arial,sans-serif}text.al{fill:var(--ink2);font-size:12px}`);
  cl.insertAdjacentHTML("afterbegin",`<style>${css}</style><rect width="100%" height="100%" fill="${fix("var(--surface)")}"/>`);
  const vb=(cl.getAttribute("viewBox")||"0 0 640 300").split(/\s+/).map(Number); cl.setAttribute("width",vb[2]); cl.setAttribute("height",vb[3]); cl.style.removeProperty("--xfs");
  return {s:new XMLSerializer().serializeToString(cl),w:vb[2],h:vb[3]}; }
function exportChart(host,name,kind){ const svg=host&&host.querySelector("svg"); if(!svg){ note("그림이 없습니다","warn"); return; } const F=svgFreeze(svg);
  if(kind==="svg"){ download(safe(name)+".svg",new Blob([F.s],{type:"image/svg+xml"})); return; }
  const img=new Image(); img.onload=()=>{ const k=3, cv=document.createElement("canvas"); cv.width=F.w*k; cv.height=F.h*k; const c=cv.getContext("2d"); c.scale(k,k); c.drawImage(img,0,0); cv.toBlob(b=>download(safe(name)+".png",b),"image/png"); }; img.src="data:image/svg+xml;charset=utf-8,"+encodeURIComponent(F.s); }
const csvDL=(name,lines)=>download(safe(name)+".csv",new Blob(["﻿"+lines.join("\n")],{type:"text/csv"}));
/* 그림 카드: 머리(제목 + SVG·PNG·CSV) + 그림 자리 */
const figCard=(id,title,o)=>{ o=o||{}; return `<div class="xws-fig" id="${id}"><div class="xws-fh"><span class="t"${o.tip?` title="${esc(o.tip)}"`:""}>${title}</span>${o.extra||""}<span class="sp"></span><button class="btn xs" data-fx="svg" title="Download the figure as SVG">SVG</button><button class="btn xs" data-fx="png" title="Download the figure as PNG (3×)">PNG</button>${o.csv?`<button class="btn xs" data-fx="csv" title="Download the numbers as CSV">CSV</button>`:""}</div><div class="xws-plot"></div>${o.below||""}</div>`; };
function bindFig(root,id,name,csvFn){ const c=$q(root,"#"+id); if(!c) return null; const host=$q(c,".xws-plot"); $a(c,"[data-fx]").forEach(b=>b.onclick=()=>{ const k=b.dataset.fx; if(k==="csv"){ if(csvFn) csvFn(); } else if(c._map) exportMapPNG(c._map,name); else exportChart(host,name,k); }); return host; }
const levelColor=(j,n)=>`hsl(${Math.round(215-200*(n>1?j/(n-1):0))} 70% 46%)`;

/* ---------- 색 지도 ---------- */
const CMAPS={viridis:[[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[109,205,89],[180,222,44],[253,231,37]],
  magma:[[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],[229,80,100],[251,135,97],[254,194,135],[252,253,191],[252,253,191]],
  gray:[[0,0,0],[255,255,255]],
  diverge:[[33,102,172],[103,169,207],[209,229,240],[247,247,247],[253,219,199],[239,138,98],[178,24,43]]};
function cmapRGB(name,t){ const L=CMAPS[name]||CMAPS.viridis; t=Math.max(0,Math.min(1,t)); const x=t*(L.length-1), i=Math.min(L.length-2,Math.floor(x)), f=x-i; return L[i].map((v,k)=>Math.round(v+(L[i+1][k]-v)*f)); }
const cmapCSS=name=>`linear-gradient(90deg,${(CMAPS[name]||CMAPS.viridis).map((c,i,a)=>`rgb(${c.join(",")}) ${Math.round(100*i/(a.length-1))}%`).join(",")})`;
function mapRange(v,robust){ const a=[]; for(let i=0;i<v.length;i++) if(isFinite(v[i])) a.push(v[i]); if(!a.length) return [0,1]; a.sort((p,q)=>p-q);
  const q=t=>a[Math.min(a.length-1,Math.max(0,Math.round(t*(a.length-1))))]; let lo=robust?q(0.01):a[0], hi=robust?q(0.99):a[a.length-1]; if(!(hi>lo)){ hi=lo+(Math.abs(lo)*0.01||1); } return [lo,hi]; }
/* m = {nx, ny, v, s?, label, unit, cmap, aspect} → 캔버스(화소 그대로 확대) + 색 막대 */
function mapHTML(m){ const ar=m.aspect||m.nx/m.ny; return `<div class="xim-box"><div class="xim-wrap" style="aspect-ratio:${ar.toFixed(4)}"><canvas class="xim" width="${m.nx}" height="${m.ny}"></canvas><div class="xim-roi" hidden></div></div>
  <div class="xim-cb"><span class="lo"></span><i style="background:${cmapCSS(m.cmap)}"></i><span class="hi"></span><span class="u">${esc(m.unit||"")}</span></div><div class="xim-tip muted xs">&nbsp;</div></div>`; }
function paintMap(box,m){ const cv=$q(box,"canvas.xim"); if(!cv) return; const [lo,hi]=m.range||mapRange(m.v,m.robust!==false); const ctx=cv.getContext("2d"); const im=ctx.createImageData(m.nx,m.ny);
  for(let p=0;p<m.nx*m.ny;p++){ const v=m.v[p], o=p*4; if(!isFinite(v)){ im.data[o]=im.data[o+1]=im.data[o+2]=128; im.data[o+3]=60; continue; } const c=cmapRGB(m.cmap,(v-lo)/(hi-lo)); im.data[o]=c[0]; im.data[o+1]=c[1]; im.data[o+2]=c[2]; im.data[o+3]=255; }
  ctx.putImageData(im,0,0); $q(box,".xim-cb .lo").textContent=g(lo,3); $q(box,".xim-cb .hi").textContent=g(hi,3); m.range=[lo,hi]; }
function exportMapPNG(m,name){ const k=Math.max(1,Math.min(32,Math.floor(900/Math.max(m.nx,m.ny)))); const W=m.nx*k, H=m.ny*k, cb=34; const cv=document.createElement("canvas"); cv.width=W; cv.height=H+cb; const c=cv.getContext("2d");
  c.fillStyle="#fff"; c.fillRect(0,0,W,H+cb); const src=document.createElement("canvas"); src.width=m.nx; src.height=m.ny; const box={querySelector:()=>null}; const [lo,hi]=m.range||mapRange(m.v,true);
  const sc=src.getContext("2d"), im=sc.createImageData(m.nx,m.ny); for(let p=0;p<m.nx*m.ny;p++){ const v=m.v[p], o=p*4; if(!isFinite(v)){ im.data[o+3]=0; continue; } const q=cmapRGB(m.cmap,(v-lo)/(hi-lo)); im.data[o]=q[0]; im.data[o+1]=q[1]; im.data[o+2]=q[2]; im.data[o+3]=255; } sc.putImageData(im,0,0);
  c.imageSmoothingEnabled=false; c.drawImage(src,0,0,W,H); const gr=c.createLinearGradient(60,0,W-60,0); (CMAPS[m.cmap]||CMAPS.viridis).forEach((q,i,a)=>gr.addColorStop(i/(a.length-1),`rgb(${q.join(",")})`)); c.fillStyle=gr; c.fillRect(60,H+8,Math.max(10,W-120),12);
  c.fillStyle="#222"; c.font="12px sans-serif"; c.textAlign="right"; c.fillText(g(lo,3),56,H+19); c.textAlign="left"; c.fillText(g(hi,3)+" "+(m.unit||""),W-56,H+19); void box;
  cv.toBlob(b=>download(safe(name)+".png",b),"image/png"); }

/* ---------- 기준 스펙트럼 후보 (TFA · 기준 NNLS · ARXPS · 영상 상태 지도 공용) ----------
   ① 저장한 기준(Store as reference, 넓이 1 template)  ② 선택한 스펙트럼(목록에서 고른 것)의 피팅 성분  ③ (영상) ROI 기준
   ev(Ecorr) → 대전 보정 BE 격자 위의 세기. 쓰기 전에 그 격자에서 넓이 1 로 맞춘다 */
function refCands(el,orb,extra){
  const out=[]; const ok=t=>!el||!t.el||(t.el===el&&(!orb||!t.orb||t.orb===orb));
  REFS().filter(ok).forEach(t=>out.push({id:"t:"+t.id,name:t.name||"Reference",kind:"Stored",ev:E=>C.tplShape(t,Float64Array.from(E,e=>e-t.eref),0)}));
  const s=cur(), r=s&&fitted(s);
  if(r&&(!el||!s.reg.el||s.reg.el===el)){ const sh=shOf(s); r.comps.forEach(c=>out.push({id:"c:"+s.uid+":"+c.id,name:`${c.name||c.id} · ${smpDisp(SK(s)).name} ${s.region}`,kind:"Fit",ev:E=>{ const Em=Float64Array.from(E,e=>e-sh); return C.compEval(c,Em,{emin:Em[0],emax:Em[Em.length-1]}); }})); }
  (extra||[]).forEach(x=>out.push(x));
  return out;
}
function refPickHTML(path,cands,sel,o){ o=o||{}; if(!cands.length) return `<div class="muted xs" title="Fit a spectrum of this element in the Spectra view and keep it selected, or use ‘Store as reference’ (Peak Fit tab)">${o.empty||"기준 없음 — Spectra 에서 같은 원소 스펙트럼을 피팅하거나 Store as reference"}</div>`;
  return `<div class="xws-refs">${cands.map(c=>`<label class="xps-f" title="${esc(c.kind)}: ${esc(c.name)}"><input type="checkbox" data-st="${esc(path+"|"+c.id)}" ${sel&&sel[c.id]?"checked":""}> <span class="badge ${c.kind==="Fit"?"inf":c.kind==="ROI"?"ok":""}">${esc(c.kind)}</span> ${esc(c.name)}</label>`).join("")}</div>`; }
const chosenRefs=(cands,sel)=>cands.filter(c=>sel&&sel[c.id]);

/* ============================================================
   데이터 묶음 · 좌표 (Series · ARXPS 공용)
   ============================================================ */
function dataSets(){ const out=[]; const G={};
  S.specs.forEach(s=>{ if(s.reg&&s.reg.kind==="survey") return; const k=SK(s); const g=G[k]||(G[k]={}); g[s.region]=(g[s.region]||0)+1; });
  Object.entries(G).forEach(([k,g])=>{ const n=Math.max(...Object.values(g)); if(n>=2) out.push({k,label:`${smpDisp(k).name} · ${Object.keys(g).length} region${Object.keys(g).length>1?"s":""} × ${n}`}); });
  if(S.ui.overlay.length>=2) out.push({k:"__checked",label:`Checked spectra (${S.ui.overlay.length})`});
  return out; }
function coordModes(specs){ const o=[["auto","Auto"]]; const vm=specs.map(s=>s.meta&&s.meta.vamas).find(Boolean);
  if(vm&&vm.expVars) vm.expVars.forEach((v,i)=>o.push(["var"+i,`${v.label}${v.unit?` (${v.unit})`:""}`]));
  if(specs.some(s=>s.meta&&s.meta.vamas&&s.meta.vamas.takeOffPolar!=null)) o.push(["tilt","Take-off angle (deg)"]);
  if(specs.some(s=>s.meta&&s.meta.mcol)) o.push(["mcol","Column header"]);
  o.push(["index","Order loaded"],["manual","Manual list"]); return o; }
function resolveMode(specs,mode){ if(mode&&mode!=="auto") return mode; const vm=specs.map(s=>s.meta&&s.meta.vamas).find(Boolean);
  if(vm&&vm.expVars&&vm.expVars.length) return "var0"; if(specs.some(s=>s.meta&&s.meta.mcol)) return "mcol";
  if(vm&&new Set(specs.map(s=>s.meta&&s.meta.vamas&&s.meta.vamas.takeOffPolar)).size>1) return "tilt"; return "index"; }
function coordVal(s,i,mode,man){ const vm=s.meta&&s.meta.vamas, mc=s.meta&&s.meta.mcol;
  if(mode.startsWith("var")&&vm){ const v=parseFloat((vm.expVals||[])[+mode.slice(3)]); return isFinite(v)?v:i; }
  if(mode==="tilt"&&vm&&vm.takeOffPolar!=null) return +vm.takeOffPolar;
  if(mode==="mcol"&&mc&&mc.coord!=null) return mc.coord;
  if(mode==="manual"&&man&&isFinite(man[i])) return man[i];
  return i; }
function coordLabel(specs,mode){ const vm=specs.map(s=>s.meta&&s.meta.vamas).find(Boolean), mc=specs.map(s=>s.meta&&s.meta.mcol).find(Boolean);
  if(mode.startsWith("var")&&vm&&vm.expVars[+mode.slice(3)]) return {label:vm.expVars[+mode.slice(3)].label,unit:vm.expVars[+mode.slice(3)].unit||""};
  if(mode==="tilt") return {label:"Take-off angle",unit:"deg"}; if(mode==="mcol"&&mc) return {label:mc.label||"Coordinate",unit:mc.unit||""};
  if(mode==="manual") return {label:"Coordinate",unit:""}; return {label:"Level",unit:"#"}; }
/* st = {set, coord, manual} → {specs, regions:[{region, el, orb, specs (좌표 순), coord}], label, unit} */
function buildSet(st){ const sets=dataSets(); if(!sets.find(x=>x.k===st.set)) st.set=sets.length?sets[0].k:null; if(!st.set) return null;
  const specs=(st.set==="__checked"?S.ui.overlay.map(byUid).filter(Boolean):sampleSpecs(st.set)).filter(s=>!(s.reg&&s.reg.kind==="survey"));
  const mode=resolveMode(specs,st.coord); const man=String(st.manual||"").split(/[,;\s]+/).filter(Boolean).map(Number);
  const byR={}; specs.forEach(s=>{ (byR[s.region]=byR[s.region]||[]).push(s); });
  const regions=Object.entries(byR).filter(([,a])=>a.length>=2).map(([region,a])=>{ const pts=a.map((s,i)=>({s,x:coordVal(s,i,mode,man),i})); pts.sort((p,q)=>p.x-q.x||p.i-q.i);
    return {region,el:(a[0].reg||{}).el||"",orb:(a[0].reg||{}).orb||"",specs:pts.map(p=>p.s),coord:pts.map(p=>p.x)}; });
  const lab=coordLabel(specs,mode); return {key:st.set,specs,regions,mode,label:lab.label,unit:lab.unit,modes:coordModes(specs)}; }
/* 한 region series 를 대전 보정 BE 의 공통 격자로 */
function serGrid(rs){ const ps=rs.specs.map(s=>{ const sh=shOf(s); return {be:Float64Array.from(s.be,e=>e+sh),y:s.y,step:s.step}; }); const g=MU().commonGrid(ps); g.c=rs.specs[0].c||null; return g; }
const xLab=(L,u)=>L+(u&&u!=="#"?` (${u})`:"");
const setPickHTML=(pfx,st,D)=>{ const sets=dataSets(); return `<div class="xps-row"><span class="xws-k">Data</span>${sIn(pfx+"|set",st.set,sets.map(x=>[x.k,x.label]),{title:"A sample with ≥ 2 spectra of a region (depth profile, angle series), or the spectra checked in the Spectra list"})}</div>
  <div class="xps-row"><span class="xws-k">Coordinate</span>${sIn(pfx+"|coord",st.coord,D?D.modes:[["auto","Auto"]],{title:"Order of the levels: VAMAS experimental variable (etch time / level / angle), take-off angle, multi-column header, load order, or a manual list"})}${st.coord==="manual"?`<input type="text" class="xps-in w120" data-st="${pfx}|manual" value="${esc(st.manual||"")}" placeholder="0, 30, 60 …" title="One value per spectrum in load order, separated by commas">`:""}</div>`; };

/* ============================================================
   Series (깊이 분포 · 각도 · 여러 스펙트럼) — SPEC_multi §1–2
   ============================================================ */
function SER(){ S.ser=S.ser||{}; defs(S.ser,{set:null,coord:"auto",manual:"",rate:null,reg:{},show:"at",sig:true,idc:false,hw:2});
  S.ser.mv=S.ser.mv||{}; defs(S.ser.mv,{region:null,lo:null,hi:null,bgsub:"shirley",scale:"poisson",nMax:6,n:null,method:"mcr",init:"simplisma",closure:false,refs:{},rbg:0,shift:0,frac:true}); return S.ser; }
function serDefS(rs){ const s0=rs.specs[0]; const lk=rsfLookup(rs.el,rs.orb||"",S.quant.lib,s0.hv); return lk?lk.v:null; }
function serDefLam(rs,R){ const s0=rs.specs[0]; const KE=s0.hv-(R.lo+R.hi)/2; const mat=(S.quant.imfp&&S.quant.imfp.mat)||{rho:5,Nv:24,M:150,Eg:3}; let L=C.tpp2m(KE,mat)/10; const th=s0.meta&&s0.meta.vamas&&s0.meta.vamas.takeOffPolar; if(th!=null&&isFinite(th)) L*=Math.cos(th*Math.PI/180); return L; }
function serReg(rs){ const st=SER(); let r=st.reg[rs.region]; if(!r){ const s0=rs.specs[0], p=pr(s0), sh=shOf(s0); r=st.reg[rs.region]={inc:true,mode:"int",lo:+(p.lo+sh).toFixed(2),hi:+(p.hi+sh).toFixed(2),bg:"shirley",S:null,lam:null}; } return r; }
function serActive(D){ const mv=SER().mv; if(!D||!D.regions.length) return null; let rs=D.regions.find(r=>r.region===mv.region); if(!rs){ rs=D.regions[0]; mv.region=rs.region; } return rs; }

function serCompute(quiet){
  const D=buildSet(SER()); if(!D||!D.regions.length){ if(!quiet) note("Series 로 쓸 스펙트럼이 없습니다 (한 region 이 2개 이상)","warn"); return false; }
  const st=SER(); const rows=[], areas=[], meta=[]; let coord=null;
  for(const rs of D.regions){ const R=serReg(rs); if(!R.inc) continue; const s0=rs.specs[0];
    if(coord&&rs.coord.length!==coord.length){ note(`${rs.region}: 단계 수(${rs.coord.length})가 다른 region(${coord.length})과 다릅니다 — 빼고 계산`,"warn"); continue; }
    const S_=R.S!=null&&R.S!==""?+R.S:serDefS(rs); if(!(S_>0)){ note(`${rs.region}: 감도계수 S 가 없습니다 — 표의 S 칸에 넣으세요`,"warn"); return false; }
    const base={el:rs.el||rs.region,S:S_,hv:s0.hv,PE:s0.meta&&s0.meta.passEnergy,tf:tfFor(s0)};
    if(R.mode==="fit"){ const fits=rs.specs.map(fitted); const miss=fits.map((r,j)=>r?null:j).filter(v=>v!=null);
      if(miss.length){ note(`${rs.region}: 피팅이 없는 단계 ${miss.length}개 — ‘Apply model & fit all levels’ 를 먼저`,"warn"); return false; }
      const names=[]; fits.forEach(r=>r.comps.forEach(c=>{ const nm=c.name||c.id; if(!c.noQuant&&!names.includes(nm)) names.push(nm); }));
      names.forEach(nm=>{ const a=fits.map((r,j)=>{ const i=r.comps.findIndex(c=>(c.name||c.id)===nm); if(i<0) return {area:0,sArea:null}; const q=r.summ[i]; return {area:q.area,sArea:q.sig&&q.sig.area!=null?q.sig.area:null,be:q.be+shOf(rs.specs[j])}; });
        const be=a.reduce((t,q)=>t+(q.be||0),0)/a.length; rows.push(Object.assign({},base,{label:`${rs.region} · ${nm}`,BE:be})); areas.push(a); meta.push({rs,R}); }); }
    else { const g=serGrid(rs); let a; try{ a=MU().integrateSeries(g.E,g.Y,{lo:R.lo,hi:R.hi,bg:{type:R.bg,n:5},cfac:g.c}); }catch(e){ note(`${rs.region}: ${e.message}`,"warn"); return false; }
      rows.push(Object.assign({},base,{label:rs.region,BE:(R.lo+R.hi)/2})); areas.push(a); meta.push({rs,R}); }
    coord=coord||rs.coord; }
  if(!rows.length){ note("포함된 region 이 없습니다","warn"); return false; }
  const comp=MU().compositionSeries(rows,areas,{imfp:S.quant.imfp});
  const rate=+st.rate>0?+st.rate:null; const x=coord.map(v=>rate?v*rate:v);
  const out={key:D.key,coord,x,depth:!!rate,label:rate?"Depth":D.label,unit:rate?"nm":D.unit,rows:rows.map(r=>r.label),S:rows.map(r=>r.S),areas,comp,idc:null,when:Date.now()};
  if(st.idc){ if(!rate) note("정보 깊이 보정은 깊이 눈금이 필요합니다 — Sputter rate 를 넣으세요","warn");
    else { const X=rows.map((r,k)=>Float64Array.from(comp.map(lv=>lv[k].nrm))); const lam=meta.map(m=>m.R.lam>0?+m.R.lam:serDefLam(m.rs,m.R));
      try{ const r=MU().infoDepthCorrect(Float64Array.from(x),X,lam,{hw:Math.max(1,Math.round(+st.hw||2))}); out.idc={at:r.at,lam,warn:r.warn}; }catch(e){ note("정보 깊이 보정: "+e.message,"warn"); } } }
  RES.ser=out; return true;
}
async function serApplyFit(rs){
  const src=(cur()&&rs.specs.includes(cur())&&pr(cur()).comps.length)?cur():rs.specs.find(s=>pr(s).comps.length);
  if(!src){ note(`${rs.region}: 성분이 있는 단계가 없습니다 — Spectra 에서 한 단계를 피팅하세요`,"warn"); return; }
  if(!fitted(src)){ await runFit("one",src); if(!fitted(src)) return; }
  S.ui.batSel={}; rs.specs.forEach(s=>{ if(s.uid!==src.uid) S.ui.batSel[s.uid]=true; }); S.ui.batSrc=src.uid;
  const ws=S.ui.ws; await runBatch(src,true); S.ui.ws=ws; S.ui.tab="fit"; serReg(rs).mode="fit"; save(); rerender();
}

/* 다변량 자료: 활성 region 의 창(보정 BE). raw = 창 안 원자료(잡음 시험·요인 수), proc = 배경을 뺀 자료(TFA·MCR·기준 NNLS) */
function mvData(D){ const rs=serActive(D); if(!rs) return null; const mv=SER().mv, R=serReg(rs); const lo=mv.lo!=null?+mv.lo:R.lo, hi=mv.hi!=null?+mv.hi:R.hi; const g=serGrid(rs);
  const [i0,i1]=sliceWin(g.E,lo,hi); if(i1-i0<7) throw new Error("창 안 점이 8개 미만입니다");
  const E=g.E.slice(i0,i1+1), raw=g.Y.map(y=>y.slice(i0,i1+1));
  const proc=mv.bgsub!=="none"?MU().subtractBackground(g.E,g.Y,{lo,hi,bg:{type:mv.bgsub,n:5}}).cols:raw.map(v=>Float64Array.from(v));
  return {rs,E,raw,proc,c:g.c,lo,hi,resampled:g.resampled}; }
async function serPCA(){
  const D=buildSet(SER()); let d; try{ d=mvData(D); }catch(e){ note(e.message,"warn"); return; } if(!d){ note("Series 가 없습니다","warn"); return; }
  if(d.raw.length<3){ note("PCA 에는 단계가 3개 이상 필요합니다","warn"); return; }
  const mv=SER().mv;
  try{ const o=await xrun("serPCA",{raw:d.raw,proc:d.proc,cfac:d.c||1,scale:mv.scale,nMax:Math.max(1,Math.min(12,+mv.nMax||6))},{label:"PCA"});
    RES.mv={d,o,res:null,key:D.key+"|"+d.rs.region};
    const nn=o.NF.nNoise; const nDef=nn==null?Math.min(3,d.raw.length-1):Math.max(1,nn); mv.n=nDef; save(); drawWS();
    note(`PCA: noise-limited factors ${nn==null?"> "+o.NF.table.length:nn} · IND ${o.mal.nIND} · F-test ${o.mal.nF}${d.c?"":" — cfac 를 몰라 count 로 가정"}`,"ok");
  }catch(e){ failNote("PCA",e); }
}
function ttRows(){ const M=RES.mv; if(!M) return []; const mv=SER().mv, n=Math.max(1,Math.min(+mv.n||1,M.o.P2.U.length-1)); const d=M.d;
  const cands=refCands(d.rs.el,d.rs.orb); return cands.map(c=>{ let x=c.ev(d.E); x=unitArea(d.E,x); try{ const t=MU().targetTest(M.o.P2,n,x); return {name:c.name,kind:c.kind,t}; }catch(e){ return {name:c.name,kind:c.kind,err:e.message}; } }); }
async function serResolve(){
  const M=RES.mv; if(!M){ note("먼저 PCA","warn"); return; } const d=M.d, mv=SER().mv; const n=Math.max(1,Math.min(+mv.n||1,d.proc.length));
  try{
    if(mv.method==="mcr"){ const o=await xrun("serMCR",{cols:d.proc,n,init:mv.init,closure:!!mv.closure},{label:"MCR-ALS"}); M.res={kind:"mcr",n,o}; note(`MCR-ALS: ${n} comp., LOF ${f(o.lof,3)} %, ${o.iter} iter.`,"ok"); }
    else { const cands=chosenRefs(refCands(d.rs.el,d.rs.orb),mv.refs); if(cands.length<1){ note("기준을 고르세요","warn"); return; }
      const refs=cands.map(c=>unitArea(d.E,c.ev(d.E))); const sub=mv.bgsub!=="none";
      const opt={E:d.E,weight:sub?"none":"poisson",cfac:d.c||1,bg:sub?Math.min(+mv.rbg||0,2):Math.max(1,Math.min(+mv.rbg||1,2)),maxShift:+mv.shift>0?+mv.shift:0,shiftStep:Math.abs(d.E[1]-d.E[0])/2};
      const o=await xrun("serRef",{cols:d.proc,refs,opt},{label:"Reference fit"}); M.res={kind:"ref",names:cands.map(c=>c.name),refs,o,weight:opt.weight};
      note(`Reference fit: ${cands.length} refs × ${o.length} levels, median χ²ν ${f(o.map(q=>q.phi).sort((a,b)=>a-b)[o.length>>1],2)}`,"ok"); }
    drawWS();
  }catch(e){ failNote("Resolve",e); }
}
/* 성분 분포(단계마다): MCR → 양(Σ 정규화한 S 의 C) 또는 분율, 기준 → 분율 ± σ (σ 는 √φ 로 키움) */
function resolveProfiles(){ const M=RES.mv; if(!M||!M.res) return null; const r=M.res, nL=M.d.proc.length;
  if(r.kind==="mcr"){ const K=r.o.C.length; const amt=r.o.C.map(c=>Float64Array.from(c)); const tot=Float64Array.from({length:nL},(_,j)=>amt.reduce((s,a)=>s+a[j],0));
    return {names:amt.map((_,k)=>"MCR "+(k+1)),amt,frac:amt.map(a=>Float64Array.from(a,(v,j)=>tot[j]>0?100*v/tot[j]:NaN)),sFrac:null,sAmt:null,spectra:r.o.S}; }
  const fits=r.o, K=r.names.length; const sc=q=>Math.sqrt(r.weight==="none"?q.phi:Math.max(q.phi,1));
  return {names:r.names,amt:[...Array(K)].map((_,k)=>Float64Array.from(fits,q=>q.coef[k])),sAmt:[...Array(K)].map((_,k)=>Float64Array.from(fits,q=>q.sig[k]*sc(q))),
    frac:[...Array(K)].map((_,k)=>Float64Array.from(fits,q=>q.fraction[k])),sFrac:[...Array(K)].map((_,k)=>Float64Array.from(fits,q=>q.sFraction[k]*sc(q))),spectra:r.refs,shift:fits.map(q=>q.shift),phi:fits.map(q=>q.phi)}; }

function wsSeries(el){
  const st=SER(), D=buildSet(st); const mv=st.mv;
  if(!D){ el.innerHTML=`<div class="card xws-empty"><b>Series</b><div class="muted sm" style="margin-top:6px">한 시료의 같은 region 스펙트럼이 2개 이상이면(깊이 분포 VAMAS·각도 series·여러 열 텍스트) 여기서 묶입니다. 목록에서 체크한 스펙트럼으로도 만들 수 있습니다.</div></div>`; return; }
  const rsA=serActive(D); const R0=rsA?serReg(rsA):null;
  const regRows=D.regions.map(rs=>{ const R=serReg(rs), k="ser|reg|"+rs.region, dS=serDefS(rs);
    return `<tr><td>${cIn(k+"|inc",R.inc,`<b>${esc(rs.region)}</b>`,{title:"Include in the composition"})}<div class="muted xs">${rs.specs.length} levels</div></td>
     <td>${sIn(k+"|mode",R.mode,[["int","Integrate"],["fit","Fit model"]],{title:"Integrate: background-subtracted window area per level. Fit model: component areas of each level's fit (Apply model & fit all levels)"})}${R.mode==="fit"?`<button class="btn xs" data-sfit="${esc(rs.region)}" title="Fit the selected level (or the first level with components), copy the model to every level and fit all (Peak Fit → Apply model to other spectra)">Apply &amp; fit all</button>`:""}</td>
     <td class="nw">${nIn(k+"|lo",R.lo,{cls:"w64",step:0.1,title:"Window low BE (charge-corrected)"})}–${nIn(k+"|hi",R.hi,{cls:"w64",step:0.1,title:"Window high BE"})}</td>
     <td>${sIn(k+"|bg",R.bg,[["linear","Linear"],["shirley","Shirley"],["smart","Smart"]],{title:"Background under the window (SPEC §2)"})}</td>
     <td>${nIn(k+"|S",R.S,{cls:"w56",ph:dS!=null?g(dS,3):"?",title:`Sensitivity factor S (empty = library ${S.quant.lib}${dS!=null?": "+g(dS,4):": none"})`})}</td>
     ${st.idc?`<td>${nIn(k+"|lam",R.lam,{cls:"w56",ph:f(serDefLam(rs,R),2),title:"Information depth λe = L·cosθ (nm) for the escape-depth correction (empty = TPP-2M with the Quantification material × cos take-off)"})}</td>`:""}</tr>`; }).join("");
  const nT=RES.mv&&RES.mv.o?RES.mv.o.P2.U.length:0;
  el.innerHTML=`<div class="card xws-side">
    <div class="xps-h">Series</div>${setPickHTML("ser",st,D)}
    <div class="xps-row"><span class="xws-k" title="Converts the coordinate to depth: z = coordinate × rate (e.g. nm/s for etch time). Empty = plot against the coordinate">Sputter rate</span>${nIn("ser|rate",st.rate,{cls:"w64",ph:"nm/unit",title:"Depth per coordinate unit (nm/s, nm/level …). Needed for the information-depth correction"})}${cIn("ser|idc",st.idc,"Info-depth corr.",{title:"Removes the exponential information-depth response from the profile: c = I − λe·dI/dz (SPEC_multi §1.5, from Hofmann's MRI model). Assumes constant λe, no sputter mixing/roughness; amplifies noise"})}${st.idc?nIn("ser|hw",st.hw,{cls:"w56",min:1,max:10,step:1,title:"Half-width (points) of the local quadratic derivative"}):""}</div>
    <div class="tw"><table class="xps-tbl xws-reg"><thead><tr><th>Region</th><th>Area from</th><th>Window (eV)</th><th>BG</th><th title="Sensitivity factor">S</th>${st.idc?`<th title="Information depth λe (nm)">λe</th>`:""}</tr></thead><tbody>${regRows}</tbody></table></div>
    <div class="xps-row"><button class="btn sm primary" id="xsrRun" title="Areas → composition per level (Quantification settings: S library, transmission, IMFP), σ by the delta method">Compute composition</button></div>
    <hr class="xps-hr"><div class="xps-h" title="PCA, factor count, target testing and pure components on one region (SPEC_multi §2)">Multivariate</div>
    <div class="xps-row"><span class="xws-k">Region</span>${sIn("ser|mv|region",mv.region,D.regions.map(r=>[r.region,r.region]))}<span class="xps-f">${nIn("ser|mv|lo",mv.lo!=null?mv.lo:(R0?R0.lo:null),{cls:"w64",step:0.1,title:"Window low BE (empty = region window)"})}–${nIn("ser|mv|hi",mv.hi!=null?mv.hi:(R0?R0.hi:null),{cls:"w64",step:0.1,title:"Window high BE"})}</span></div>
    <div class="xps-row"><span class="xws-k" title="Subtracted from every level before target testing / MCR / reference fit. The noise test always uses the raw window (Poisson statistics)">Background</span>${sIn("ser|mv|bgsub",mv.bgsub,[["shirley","Subtract Shirley"],["linear","Subtract linear"],["none","Keep"]])}${sIn("ser|mv|scale",mv.scale,[["poisson","Poisson scaling"],["none","No scaling"],["area","Area scaling"]],{title:"PCA preprocessing for the eigenvalues / IND / F-test (Keenan & Kotula Poisson scaling recommended for counts)"})}</div>
    <div class="xps-row"><button class="btn sm primary" id="xsrPca" title="Eigenvalues, Malinowski IND / F-test and the noise-limited factor count (MLPCA weighted low-rank fit, Poisson)">PCA &amp; factor count</button><span class="xps-f" title="Largest factor count tested">max ${nIn("ser|mv|nMax",mv.nMax,{cls:"w56",min:1,max:12,step:1})}</span></div>
    <div class="xps-row"><span class="xws-k">Factors n</span>${nIn("ser|mv|n",mv.n,{cls:"w56",min:1,max:12,step:1,title:"Number of components for target testing and resolution (default: noise-limited count of the raw window; a level-dependent background counts as one factor — lower n by 1 if it was subtracted)"})}${nT?"":`<span class="muted xs">PCA 후</span>`}</div>
    <div class="xps-row"><span class="xws-k">Resolve</span>${sIn("ser|mv|method",mv.method,[["mcr","MCR-ALS"],["ref","Reference fit (NNLS)"]],{title:"MCR-ALS: non-negative pure spectra + profiles (rotational ambiguity!). Reference fit: non-negative fit of each level to chosen reference spectra"})}
      ${mv.method==="mcr"?sIn("ser|mv|init",mv.init,[["simplisma","SIMPLISMA start"],["opa","OPA start"],["ittfa","ITTFA start"]],{title:"Initial estimate (SIMPLISMA is the validated default)"})+cIn("ser|mv|closure",mv.closure,"Closure",{title:"Force the profiles of each level to sum to 1"}):
       sIn("ser|mv|rbg",mv.rbg,[["0","No bg term"],["1","+ constant"],["2","+ linear"]],{title:"Sign-free background terms in the reference fit (with the raw window at least a constant is used)"})+`<span class="xps-f" title="Common energy shift searched per level (± eV, 0 = off)">shift ± ${nIn("ser|mv|shift",mv.shift,{cls:"w56",min:0,max:3,step:0.05})}</span>`}</div>
    ${mv.method==="ref"?refPickHTML("ser|mv|refs",rsA?refCands(rsA.el,rsA.orb):[],mv.refs):""}
    <div class="xps-row"><button class="btn sm primary" id="xsrRes" ${RES.mv?"":"disabled"} title="Run in a Worker (Stop to cancel)">Resolve components</button></div>
  </div>
  <div class="card xws-main"><div class="xws-pg" id="xwsProgSlot"></div>
    ${figCard("xsrFigC","Composition vs "+esc(RES.ser?xLab(RES.ser.label,RES.ser.unit):xLab(D.label,D.unit)),{csv:true,tip:"at % (or normalized intensity n = A/(S·T·λ)) per level, ± 1σ (Poisson, delta method)",extra:` ${sIn("ser|show",st.show,[["at","at %"],["nrm","Normalized n"]])} ${cIn("ser|sig",st.sig,"± σ")}`})}
    <div id="xsrTbl"></div>
    ${rsA?figCard("xsrFigS",`${esc(rsA.region)} — all levels`,{tip:"Charge-corrected spectra of every level (blue → red = first → last); dashed = multivariate window"}):""}
    ${RES.mv&&RES.mv.key===D.key+"|"+(rsA&&rsA.region)?figCard("xsrFigP","Eigenvalues · factor count",{csv:true,tip:"Scree plot (log). Table: Malinowski RE/IND/F-test (homoscedastic assumption) and the MLPCA noise test χ²/dof (Poisson, recommended)"})+`<div id="xsrPcaT"></div><div id="xsrTT"></div>`:""}
    ${RES.mv&&RES.mv.res?figCard("xsrFigR","Component profiles",{csv:true,extra:` ${cIn("ser|mv|frac",mv.frac,"Fraction %",{title:"Fraction of the summed component amounts at each level (otherwise amounts)"})}`})+figCard("xsrFigQ","Pure-component spectra",{tip:RES.mv.res.kind==="mcr"?"MCR-ALS spectra (unit sum). Rotational ambiguity: other non-negative solutions may fit equally well":"Reference spectra used (unit area)"}):""}
  </div>`;
  bindSt(el,(k)=>{ if(k==="set"||k==="coord"){ RES.ser=null; RES.mv=null; } if(k==="region"){ RES.mv=null; mv.lo=null; mv.hi=null; } if(["show","sig","frac"].includes(k)){ drawWS(); return; } if(RES.ser&&["inc","mode","lo","hi","bg","S","lam","rate","idc","hw"].includes(k)) serCompute(true); drawWS(); });
  $a(el,"[data-sfit]").forEach(b=>b.onclick=()=>{ const rs=D.regions.find(r=>r.region===b.dataset.sfit); if(rs) serApplyFit(rs); });
  $q(el,"#xsrRun").onclick=()=>{ if(serCompute()) drawWS(); };
  $q(el,"#xsrPca").onclick=serPCA; $q(el,"#xsrRes").onclick=serResolve;
  // 그림
  const R=RES.ser&&RES.ser.key===D.key?RES.ser:null;
  const hC=bindFig(el,"xsrFigC","series_composition",()=>serCSV(R));
  if(hC){ if(!R) hC.innerHTML=`<div class="xps-empty sm">Compute composition</div>`;
    else { const k=st.show; const ser=R.rows.map((nm,r)=>({name:nm,c:`var(--s${r%8+1})`,x:R.x,y:R.comp.map(lv=>k==="at"?lv[r].at:lv[r].nrm),err:st.sig?R.comp.map(lv=>k==="at"?lv[r].sAt:lv[r].sNrm):null,pts:true}));
      if(R.idc&&k==="at") R.rows.forEach((nm,r)=>ser.push({name:nm+" (info-depth corr.)",c:`var(--s${r%8+1})`,x:R.x,y:Array.from(R.idc.at[r]),dash:"5 4",w:1.4}));
      putChart(hC,{id:"xsrC",x:{label:xLab(R.label,R.unit)},y:{label:k==="at"?"at %":"n = A/(S·T·λ)",zero:true},series:ser}); }
    const tb=$q(el,"#xsrTbl"); if(tb&&R) tb.innerHTML=`<details class="xps-d" id="xsrTblD"><summary>Table <span class="muted xs">(${R.x.length} levels)</span></summary><div class="b tw"><table class="xps-tbl"><thead><tr><th class="num">${esc(xLab(R.label,R.unit))}</th>${R.rows.map(n=>`<th class="num">${esc(n)}</th>`).join("")}</tr></thead><tbody>${R.x.map((x,j)=>`<tr><td class="num">${g(x,4)}</td>${R.rows.map((_,r)=>{ const q=R.comp[j][r]; return `<td class="num">${q.at!=null?f(q.at,2):(q.nd?"n.d.":"—")}${q.sAt?` <span class="muted xs">± ${f(q.sAt,2)}</span>`:""}</td>`; }).join("")}</tr>`).join("")}</tbody></table>${R.idc?`<div class="muted xs" title="${esc(R.idc.warn.join(" · "))}">Info-depth correction λe = ${R.idc.lam.map(v=>f(v,2)).join(", ")} nm</div>`:""}</div></details>`; }
  const hS=bindFig(el,"xsrFigS","series_spectra"); if(hS&&rsA){ const gq=serGrid(rsA); const n=gq.Y.length; const lo=mv.lo!=null?+mv.lo:R0.lo, hi=mv.hi!=null?+mv.hi:R0.hi;
    putChart(hS,{id:"xsrS",h:260,x:{label:"Binding energy (eV)",rev:true},y:{label:"Intensity"},vlines:[{x:lo},{x:hi}],series:gq.Y.map((y,j)=>({x:gq.E,y,c:levelColor(j,n),w:1.1}))}); }
  if(RES.mv&&RES.mv.key===D.key+"|"+(rsA&&rsA.region)) drawMV(el);
}
function drawMV(el){ const M=RES.mv, o=M.o, mv=SER().mv;
  const hP=bindFig(el,"xsrFigP","series_pca",()=>{ const L=["n,eigenvalue,var_pct,RE,IND,F_SL,noise_chi2,noise_dof,noise_limit,within_noise"]; o.eig.forEach((e,i)=>{ const m=o.mal.table[i]||{}, nf=o.NF.table.find(t=>t.n===i+1)||{}; L.push([i+1,e,o.varPct[i],m.RE,m.IND,m.SL,nf.chi2,nf.dof,nf.limit,nf.within].join(",")); }); csvDL("series_pca",L); });
  if(hP){ const K=Math.min(o.eig.length,15); const nn=o.NF.nNoise; putChart(hP,{id:"xsrP",h:230,x:{label:"Factor",min:0.5,max:K+0.5},xTicks:[...Array(K)].map((_,i)=>({v:i+1,t:String(i+1)})),y:{label:"Eigenvalue",log:true},vlines:nn!=null?[{x:nn+0.5,t:"noise",c:"var(--warn)"}]:[],series:[{x:[...Array(K)].map((_,i)=>i+1),y:Array.from(o.eig).slice(0,K),pts:true,c:"var(--s1)"}]}); }
  const t=$q(el,"#xsrPcaT"); if(t){ const nf=o.NF; const rows=nf.table.map(q=>{ const m=q.n>0?o.mal.table[q.n-1]:null; return `<tr${q.n===nf.nNoise?` class="tot"`:""}><td class="num">${q.n}</td><td class="num">${m?f(o.varPct[q.n-1],2):"—"}</td><td class="num">${m?g(m.IND,3):"—"}</td><td class="num">${m?g(m.SL,2):"—"}</td><td class="num">${f(q.ratio,3)}</td><td>${q.within?`<span class="badge ok">within noise</span>`:`<span class="badge">above</span>`}</td></tr>`; }).join("");
    t.innerHTML=`<div class="strip"><div class="${nf.nNoise==null?"warn":"ok"}" title="Smallest n whose MLPCA reconstruction reaches χ² ≤ dof + 2√(2dof) (Poisson; SPEC_multi §2.3)"><b>${nf.nNoise==null?"> "+(nf.table.length-1):nf.nNoise}</b>noise-limited n</div><div title="Malinowski indicator function minimum"><b>${o.mal.nIND}</b>IND</div><div title="Malinowski F-test, 5 % (under-counts for Poisson data with a large intensity range)"><b>${o.mal.nF}</b>F-test 5 %</div></div>
     <details class="xps-d"><summary>Factor table</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th class="num">n</th><th class="num">var %</th><th class="num">IND</th><th class="num" title="F-test significance level">SL</th><th class="num" title="MLPCA χ²/dof">χ²/dof</th><th></th></tr></thead><tbody>${rows}</tbody></table></div></details>`; }
  const tt=$q(el,"#xsrTT"); if(tt){ const L=ttRows(); tt.innerHTML=`<details class="xps-d" open><summary title="Malinowski target testing on the ${mv.bgsub!=="none"?"background-subtracted":"raw"} window with n = ${mv.n} factors: SPOIL < 3 acceptable, 3–6 moderate, > 6 rejected">Target testing <span class="muted xs">(n = ${mv.n})</span></summary><div class="b">${L.length?`<div class="tw"><table class="xps-tbl"><thead><tr><th>Target</th><th class="num">SPOIL</th><th class="num" title="Apparent error of the target">AET</th><th></th></tr></thead><tbody>${L.map(r=>`<tr><td><span class="badge ${r.kind==="Fit"?"inf":""}">${esc(r.kind)}</span> ${esc(r.name)}</td><td class="num">${r.t?f(r.t.SPOIL,2):"—"}</td><td class="num">${r.t?g(r.t.AET,3):"—"}</td><td>${r.t?`<span class="badge ${r.t.SPOIL<3?"ok":r.t.SPOIL<=6?"wr":"rk"}">${esc(r.t.verdict)}</span>`:`<span class="muted xs">${esc(r.err||"")}</span>`}</td></tr>`).join("")}</tbody></table></div>`:`<div class="muted xs">시험할 기준 없음 — Spectra 에서 피팅한 스펙트럼을 선택하거나 Store as reference</div>`}</div></details>`; }
  const P=resolveProfiles(); if(!P) return;
  const hR=bindFig(el,"xsrFigR","series_profiles",()=>{ const R=RES.ser; const x=M.d.rs.coord; const L=[["level_coord",...P.names.flatMap(n=>[n+"_amount",n+"_s_amount",n+"_fraction_pct",n+"_s_fraction_pct"]),...(P.shift?["shift_eV","chi2nu"]:[])].map(csvE).join(",")];
    x.forEach((v,j)=>L.push([v,...P.names.flatMap((_,k)=>[P.amt[k][j],P.sAmt?P.sAmt[k][j]:"",P.frac[k][j],P.sFrac?P.sFrac[k][j]:""]),...(P.shift?[P.shift[j],P.phi[j]]:[])].join(","))); void R; csvDL("series_profiles",L); });
  if(hR){ const fr=mv.frac; putChart(hR,{id:"xsrR",x:{label:xLab(buildSet(SER()).label,buildSet(SER()).unit)},y:{label:fr?"Fraction %":"Amount",zero:true},series:P.names.map((nm,k)=>({name:nm,c:`var(--s${k%8+1})`,x:M.d.rs.coord,y:Array.from(fr?P.frac[k]:P.amt[k]),err:(fr?P.sFrac:P.sAmt)?Array.from((fr?P.sFrac:P.sAmt)[k]):null,pts:true}))}); }
  const hQ=bindFig(el,"xsrFigQ","series_pure_spectra"); if(hQ) putChart(hQ,{id:"xsrQ",h:240,x:{label:"Binding energy (eV)",rev:true},y:{label:"Intensity (unit area)"},series:P.spectra.map((s,k)=>({name:P.names[k],c:`var(--s${k%8+1})`,x:M.d.E,y:Array.from(s)}))});
  if(M.res.kind==="mcr"&&hQ){ hQ.insertAdjacentHTML("beforeend",`<div class="muted xs" title="${esc((M.res.o.warn||[]).join(" "))}">LOF ${f(M.res.o.lof,3)} % · ${M.res.o.iter} iterations · rotational ambiguity ⓘ</div>`); }
}
function serCSV(R){ if(!R){ note("먼저 Compute composition","warn"); return; } const L=[[xLab(R.label,R.unit),...(R.depth?["coordinate"]:[]),...R.rows.flatMap(n=>[n+"_area",n+"_s_area",n+"_n",n+"_at_pct",n+"_s_at_pct",...(R.idc?[n+"_at_pct_infodepth"]:[])])].map(csvE).join(",")];
  R.x.forEach((x,j)=>L.push([x,...(R.depth?[R.coord[j]]:[]),...R.rows.flatMap((_,k)=>{ const q=R.comp[j][k], a=R.areas[k][j]; return [a.area,a.sArea!=null?a.sArea:"",q.nrm,q.at!=null?q.at:(q.nd?"n.d.":""),q.sAt!=null?q.sAt:"",...(R.idc?[R.idc.at[k][j]]:[])]; })].join(",")));
  L.push(`# S: ${R.rows.map((n,k)=>n+"="+R.S[k]).join("; ")} · T: ${S.quant.tf.model} · IMFP: ${JSON.stringify(S.quant.imfp)}`); csvDL("series_composition",L); }

/* ============================================================
   ARXPS (SPEC_multi §4) — 각도 series 의 세기 → 덮개층 두께 · 다층 · 상대 깊이 · MaxEnt
   ============================================================ */
function ARX(){ S.arx=S.arx||{}; const a=S.arx; defs(a,{set:null,coord:"auto",manual:"",from:"normal",src:"ref",region:null,lo:null,hi:null,bgsub:"shirley",refs:{},wins:[],
  o:0,s:1,Lo:null,Ls:null,Rinf:1,fitR:false,ml:null,nSet:1,J:20,zmax:null,sfl:1,al:null});
  a.ml=a.ml||{N:1,layers:[{name:"Overlayer",d0:1,fix:false},{name:"Layer 2",d0:1,fix:false},{name:"Layer 3",d0:1,fix:false}],sp:[]};
  a.al=a.al||{kind:"inorganic",M:60,rho:2.2,n:3,omega:0,how:"sd"}; return a; }
function arxAngles(D){ const a=ARX(); return Float64Array.from(D.regions[0].coord,v=>a.from==="surface"?90-v:v); }
/* 세기 표: I[k][a] (넓이), σ[k][a] — 기준 NNLS / 창 넓이 / 피팅 성분 */
function arxIntensities(){
  const a=ARX(), D=buildSet(a); if(!D||!D.regions.length) throw new Error("각도 series 가 없습니다");
  const th=arxAngles(D), nA=th.length; let names=[], I=[], sI=[], info="";
  if(a.src==="ref"){ const rs=D.regions.find(r=>r.region===a.region)||D.regions[0]; a.region=rs.region; const g=serGrid(rs);
    const lo=a.lo!=null?+a.lo:g.E[0], hi=a.hi!=null?+a.hi:g.E[g.E.length-1]; const sub=a.bgsub!=="none";
    let E, cols; if(sub){ const sb=MU().subtractBackground(g.E,g.Y,{lo,hi,bg:{type:a.bgsub,n:5}}); E=sb.E; cols=sb.cols; } else { const [i0,i1]=sliceWin(g.E,lo,hi); E=g.E.slice(i0,i1+1); cols=g.Y.map(y=>y.slice(i0,i1+1)); }
    const cands=chosenRefs(refCands(rs.el,rs.orb),a.refs); if(cands.length<2) throw new Error("기준을 2개 이상 고르세요 (덮개층·기판)");
    const refs=cands.map(c=>unitArea(E,c.ev(E))); const fits=MU().fitReferences(cols,refs,{E,weight:sub?"none":"poisson",cfac:g.c||1,bg:sub?0:1});
    names=cands.map(c=>c.name); I=cands.map((_,k)=>Array.from(fits,q=>q.coef[k])); sI=cands.map((_,k)=>Array.from(fits,q=>q.sig[k]*Math.sqrt(sub?q.phi:Math.max(q.phi,1))));
    info=`Reference fit on ${rs.region} ${f(E[0],1)}–${f(E[E.length-1],1)} eV (${sub?"background subtracted, unweighted":"Poisson-weighted + constant"}) · median χ²ν ${f(fits.map(q=>q.phi).sort((x,y)=>x-y)[fits.length>>1],2)}`; }
  else if(a.src==="win"){ if(!a.wins.length) throw new Error("창을 하나 이상 더하세요");
    a.wins.forEach(w=>{ const rs=D.regions.find(r=>r.region===w.region)||D.regions[0]; if(rs.coord.length!==nA) throw new Error(rs.region+": 각도 수가 다릅니다"); const gq=serGrid(rs);
      const ar=MU().integrateSeries(gq.E,gq.Y,{lo:+w.lo,hi:+w.hi,bg:{type:w.bg||"shirley",n:5},cfac:gq.c}); names.push(w.name||`${rs.region} ${w.lo}–${w.hi}`);
      I.push(ar.map(q=>q.area)); sI.push(ar.map(q=>q.sArea!=null?q.sArea:0.01*Math.abs(q.area))); });
    info="Window areas"+(D.regions[0].specs[0].c?"":" — σ: cfac 를 몰라 1 % 로 둠"); }
  else { D.regions.forEach(rs=>{ const fits=rs.specs.map(fitted); if(fits.some(r=>!r)) return; const nm=[]; fits.forEach(r=>r.comps.forEach(c=>{ const k=c.name||c.id; if(!nm.includes(k)) nm.push(k); }));
      nm.forEach(k=>{ names.push(`${rs.region} · ${k}`); I.push(fits.map(r=>{ const i=r.comps.findIndex(c=>(c.name||c.id)===k); return i<0?0:r.summ[i].area; })); sI.push(fits.map(r=>{ const i=r.comps.findIndex(c=>(c.name||c.id)===k); const s=i<0?null:r.summ[i].sig&&r.summ[i].sig.area; return s>0?s:0.01*Math.max(1e-9,i<0?1:r.summ[i].area); })); }); });
    if(!names.length) throw new Error("모든 각도가 피팅된 region 이 없습니다 — Series 의 ‘Apply & fit all’ 또는 Peak Fit 의 ‘Apply model to other spectra’");
    info="Fitted component areas"; }
  const out={th:Array.from(th),names,I,sI,info,key:D.key,label:D.label};
  // 종마다 운동에너지(감쇠 길이 계산용): 기준·성분은 region 의 BE 중심
  RES.arx={tab:out,ov:null,ml:null,rdp:null,me:null}; return out;
}
function arxKE(k){ const a=ARX(), D=buildSet(a); if(!D) return null; const rs=a.src==="win"?(D.regions.find(r=>r.region===(a.wins[k]||{}).region)||D.regions[0]):(D.regions.find(r=>r.region===a.region)||D.regions[0]);
  const s0=rs.specs[0]; const mid=a.src==="win"&&a.wins[k]?(+a.wins[k].lo+ +a.wins[k].hi)/2:(s0.be[0]+s0.be[s0.n-1])/2+shOf(s0); return s0.hv-mid; }
/* 감쇠 길이 계산기: Seah–Dench (원소·무기·유기) 또는 TPP-2M IMFP × (1 − 0.735ω) */
function arxAL(KE){ const al=ARX().al; if(!(KE>0)) return null;
  try{ if(al.how==="imfp"){ const mat=(S.quant.imfp&&S.quant.imfp.mat)||{rho:5,Nv:24,M:150,Eg:3}; return ARX_().ealFromImfp(C.tpp2m(KE,mat)/10,+al.omega||0); }
    if(al.kind==="organic") return ARX_().seahDench(KE,{kind:"organic",rhoKgM3:(+al.rho||1)*1000});
    return ARX_().seahDench(KE,{kind:al.kind,a:ARX_().atomSize(+al.M,+al.rho,+al.n||1)}); }catch(e){ return null; } }
function arxFitOv(){ const T=RES.arx&&RES.arx.tab; if(!T) return; const a=ARX(); const o=+a.o, s=+a.s; if(o===s||!T.I[o]||!T.I[s]){ note("덮개층·기판 종을 다르게 고르세요","warn"); return; }
  if(!(+a.Lo>0&&+a.Ls>0)){ note("감쇠 길이 L_o · L_s 를 넣으세요 (Estimate 단추)","warn"); return; }
  const R=T.th.map((_,i)=>T.I[o][i]/T.I[s][i]); const sR=T.th.map((_,i)=>R[i]*Math.hypot(T.sI[o][i]/T.I[o][i],T.sI[s][i]/T.I[s][i]));
  const opt={Lo:+a.Lo,Ls:+a.Ls,Rinf:+a.Rinf||1,fitRinf:!!a.fitR}; const r=ARX_().fitOverlayer(T.th,R,sR,opt); const h=ARX_().hillPlot(T.th,R,sR,{Lo:opt.Lo,Rinf:r.Rinf});
  RES.arx.ov={R,sR,r,h,opt}; note(`Overlayer d = ${f(r.d,3)} ± ${f(r.sdScaled,3)} nm (χ²ν ${f(r.phi,2)})`,"ok"); drawWS(); }
function arxFitMl(){ const T=RES.arx&&RES.arx.tab; if(!T) return; const a=ARX(), ml=a.ml, N=Math.max(1,Math.min(3,+ml.N||1));
  const sp=T.names.map((nm,k)=>{ const q=ml.sp[k]||{}; return {name:nm,layer:q.layer!=null?+q.layer:(k===0?0:N),Iinf:+q.Iinf>0?+q.Iinf:1,L:+q.L>0?+q.L:arxAL(arxKE(k))}; });
  if(sp.some(q=>!(q.L>0))){ note("종마다 L 을 넣으세요","warn"); return; }
  const model={layers:ml.layers.slice(0,N).map(l=>({name:l.name,d:+l.d0>0?+l.d0:1})),species:sp.map(q=>({name:q.name,layer:q.layer,Iinf:q.Iinf,L:[q.L,q.L,q.L]}))};
  try{ const r=ARX_().fitMultilayer(model,T.th,T.I,T.sI,{fix:ml.layers.slice(0,N).map(l=>!!l.fix),d0:model.layers.map(l=>l.d)}); RES.arx.ml={r,model,N}; note(`Multilayer: ${r.d.slice?Array.from(r.d).map(v=>f(v,2)).join(" / "):""} nm (χ²ν ${f(r.phi,2)})`,"ok"); drawWS(); }catch(e){ note("Multilayer: "+e.message,"warn"); } }
function arxRDP(){ const T=RES.arx&&RES.arx.tab; if(!T) return; try{ RES.arx.rdp=ARX_().relativeDepth(T.th,T.I,T.sI,{nSet:Math.max(1,Math.min(Math.floor(T.th.length/2),+ARX().nSet||1))}); drawWS(); }catch(e){ note("RDP: "+e.message,"warn"); } }
async function arxMaxent(){ const T=RES.arx&&RES.arx.tab; if(!T) return; const a=ARX(), ml=a.ml; if(T.names.length<2){ note("MaxEnt 에는 종이 2개 이상 필요합니다","warn"); return; }
  const L=T.names.map((_,k)=>{ const q=ml.sp[k]||{}; return +q.L>0?+q.L:arxAL(arxKE(k)); }); if(L.some(v=>!(v>0))){ note("종마다 L 을 넣으세요 (Multilayer 표)","warn"); return; }
  /* MaxEnt 모형은 종마다 같은 bulk 세기를 가정(c 의 합 = 1) → 세기를 I∞(다층 표, 상대값)로 나눠 넣는다 */
  const Iinf=T.names.map((_,k)=>{ const q=ml.sp[k]||{}; return +q.Iinf>0?+q.Iinf:1; }); const fl=Math.max(0,+a.sfl||0)/100;   // σ 바닥(모형 오차): 계수 σ 와 제곱합
  const I=T.I.map((r,k)=>r.map(v=>v/Iinf[k])), sI=T.sI.map((r,k)=>r.map((v,i)=>Math.hypot(v/Iinf[k],fl*I[k][i])));
  try{ const r=await xrun("arxMaxent",{th:T.th,L,I,sI,opt:{J:Math.max(5,Math.min(60,+a.J||20)),zmax:+a.zmax>0?+a.zmax:undefined}},{label:"MaxEnt"}); RES.arx.me={r,L,Iinf,fl}; note(`MaxEnt: α ${g(r.alpha,3)}, χ² ${f(r.chi2,1)} / N ${r.N}${r.notReached?" — 잡음 수준에 닿지 못함":""}`,r.notReached?"warn":"ok"); drawWS(); }catch(e){ failNote("MaxEnt",e); } }

function wsArxps(el){
  const a=ARX(), D=buildSet(a);
  if(!D){ el.innerHTML=`<div class="card xws-empty"><b>ARXPS</b><div class="muted sm" style="margin-top:6px">각도 series(VAMAS Angle 변수 · 블록마다 take-off 각)를 열면 여기서 덮개층 두께·다층·상대 깊이·MaxEnt 깊이 분포를 구합니다.</div></div>`; return; }
  const rs=D.regions.find(r=>r.region===a.region)||D.regions[0]; a.region=rs.region; const th=arxAngles(D);
  const T=RES.arx&&RES.arx.tab&&RES.arx.tab.key===D.key?RES.arx.tab:null; const names=T?T.names:[];
  const spOpts=names.map((n,k)=>[k,n]);
  const srcBlock=a.src==="ref"?`<div class="xps-row"><span class="xws-k">Region</span>${sIn("arx|region",a.region,D.regions.map(r=>[r.region,r.region]))}<span class="xps-f">${nIn("arx|lo",a.lo,{cls:"w64",step:0.1,ph:"lo",title:"Window low BE (empty = full)"})}–${nIn("arx|hi",a.hi,{cls:"w64",step:0.1,ph:"hi",title:"Window high BE"})}</span>${sIn("arx|bgsub",a.bgsub,[["shirley","− Shirley"],["linear","− linear"],["none","Keep bg"]],{title:"Background subtracted before the reference fit (then unweighted, σ scaled by √χ²ν)"})}</div>${refPickHTML("arx|refs",refCands(rs.el,rs.orb),a.refs)}`
    :a.src==="win"?`<div class="tw"><table class="xps-tbl"><thead><tr><th>Name</th><th>Region</th><th>Window (eV)</th><th>BG</th><th></th></tr></thead><tbody>${a.wins.map((w,i)=>`<tr><td><input type="text" class="xps-in w90" data-st="arx|wins|${i}|name" value="${esc(w.name||"")}"></td><td>${sIn(`arx|wins|${i}|region`,w.region,D.regions.map(r=>[r.region,r.region]))}</td><td class="nw">${nIn(`arx|wins|${i}|lo`,w.lo,{cls:"w56",step:0.1})}–${nIn(`arx|wins|${i}|hi`,w.hi,{cls:"w56",step:0.1})}</td><td>${sIn(`arx|wins|${i}|bg`,w.bg||"shirley",[["linear","Lin."],["shirley","Shirley"],["smart","Smart"]])}</td><td><button class="x" data-wdel="${i}" title="Remove">×</button></td></tr>`).join("")}</tbody></table></div><button class="btn xs" id="xarWAdd">+ Window</button>`
    :`<div class="muted xs">각 각도의 피팅 결과(성분 넓이 ± σ). Series 의 ‘Apply &amp; fit all’ 로 만듭니다.</div>`;
  const ml=a.ml, N=Math.max(1,Math.min(3,+ml.N||1));
  el.innerHTML=`<div class="card xws-side">
    <div class="xps-h">ARXPS</div>${setPickHTML("arx",a,D)}
    <div class="xps-row"><span class="xws-k">Angles</span>${sIn("arx|from",a.from,[["normal","from normal"],["surface","from surface"]],{title:"Engine convention: emission angle θ from the surface normal. Choose ‘from surface’ if the coordinate is a take-off angle measured from the surface (θ = 90° − value)"})}<span class="muted xs" title="${esc(Array.from(th).map(v=>f(v,1)).join(", "))}">${th.length} angles, ${f(Math.min(...th),0)}–${f(Math.max(...th),0)}°</span></div>
    <div class="xps-row"><span class="xws-k">Intensities</span>${sIn("arx|src",a.src,[["ref","Reference fit"],["win","Windows"],["fit","Fitted components"]],{title:"How species intensities per angle are obtained: non-negative fit to reference shapes (overlapping lines, e.g. Si⁰/Si⁴⁺), window areas, or the per-angle peak fits"})}</div>
    ${srcBlock}
    <div class="xps-row"><button class="btn sm primary" id="xarGet">Get intensities</button></div>
    ${T?`<div class="muted xs" title="${esc(T.info)}">${esc(T.names.length+" species × "+T.th.length+" angles")} ⓘ</div>`:""}
    <hr class="xps-hr"><div class="xps-h" title="Attenuation length L (nm): Seah–Dench 1979 (elements / inorganic / organic) or TPP-2M IMFP × (1 − 0.735ω) (NIST SRD 82). The Cumpson–Seah CS2 formula is not implemented — enter L directly if you have it">Attenuation length</div>
    <div class="xps-row">${sIn("arx|al|how",a.al.how,[["sd","Seah–Dench"],["imfp","TPP-2M IMFP"]])}${a.al.how==="sd"?sIn("arx|al|kind",a.al.kind,[["element","Element"],["inorganic","Inorganic"],["organic","Organic"]])+(a.al.kind==="organic"?`<span class="xps-f" title="Density g/cm³">ρ ${nIn("arx|al|rho",a.al.rho,{cls:"w56"})}</span>`:`<span class="xps-f" title="Molar mass (g/mol), density (g/cm³), atoms per formula unit → atom size a">M ${nIn("arx|al|M",a.al.M,{cls:"w56"})} ρ ${nIn("arx|al|rho",a.al.rho,{cls:"w56"})} n ${nIn("arx|al|n",a.al.n,{cls:"w56",step:1})}</span>`):`<span class="xps-f" title="Single-scattering albedo ω (material; TPP-2M uses the Quantification tab material)">ω ${nIn("arx|al|omega",a.al.omega,{cls:"w56",step:0.05,min:0,max:1})}</span>`}</div>
  </div>
  <div class="card xws-main"><div class="xws-pg" id="xwsProgSlot"></div>
    ${T?`<details class="xps-d" id="xarID" open><summary>Intensities</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th class="num">θ (°)</th>${T.names.map(n=>`<th class="num">${esc(n.length>28?n.slice(0,27)+"…":n)}</th>`).join("")}</tr></thead><tbody>${T.th.map((t,i)=>`<tr><td class="num">${f(t,1)}</td>${T.names.map((_,k)=>`<td class="num">${g(T.I[k][i],4)} <span class="muted xs">± ${g(T.sI[k][i],2)}</span></td>`).join("")}</tr>`).join("")}</tbody></table><div class="xps-row"><button class="btn xs" id="xarICsv">CSV</button></div></div></details>`:`<div class="xps-empty sm">Get intensities</div>`}
    ${T?`<details class="xps-d" id="xarOvD" open><summary title="Single overlayer (Hill et al. 1976): R = I_o/I_s = R∞(1 − e^(−d/(L_o cosθ)))·e^(d/(L_s cosθ)); weighted least squares over all angles">Overlayer thickness</summary><div class="b">
      <div class="xps-row"><span class="xps-f">Overlayer ${sIn("arx|o",a.o,spOpts)}</span><span class="xps-f">Substrate ${sIn("arx|s",a.s,spOpts)}</span></div>
      <div class="xps-row"><span class="xps-f" title="Attenuation length of overlayer electrons in the overlayer (nm)">L<sub>o</sub> ${nIn("arx|Lo",a.Lo,{cls:"w56",step:0.05})}</span><span class="xps-f" title="Attenuation length of substrate electrons in the overlayer (nm)">L<sub>s</sub> ${nIn("arx|Ls",a.Ls,{cls:"w56",step:0.05})}</span><button class="btn xs" id="xarEst" title="Fill L_o and L_s from the calculator (left) at each species' kinetic energy">Estimate</button>
       <span class="xps-f" title="R∞ = I_o∞ / I_s∞ (bulk intensities: densities × λ × sensitivity ratio)">R∞ ${nIn("arx|Rinf",a.Rinf,{cls:"w56",step:0.05})}</span>${cIn("arx|fitR",a.fitR,"fit R∞")}<button class="btn sm primary" id="xarOv">Fit</button></div>
      <div id="xarOvR"></div></div></details>
      ${figCard("xarFigR","R(θ) and model",{csv:true})}${figCard("xarFigH","Hill plot",{tip:"ln(1 + R/R∞) vs 1/cosθ; slope = d/L — exact only when L_o = L_s (use the fit otherwise)"})}
    <details class="xps-d" id="xarMlD"><summary title="Up to 3 layers on a substrate; per-angle instrument factor profiled out (SPEC_multi §4.3)">Multilayer</summary><div class="b">
      <div class="xps-row"><span class="xps-f">Layers ${sIn("arx|ml|N",N,[[1,"1"],[2,"2"],[3,"3"]])}</span>${ml.layers.slice(0,N).map((l,i)=>`<span class="xps-f"><input type="text" class="xps-in w90" data-st="arx|ml|layers|${i}|name" value="${esc(l.name)}"> d₀ ${nIn(`arx|ml|layers|${i}|d0`,l.d0,{cls:"w56",step:0.1,title:"Start (or fixed) thickness, nm"})}${cIn(`arx|ml|layers|${i}|fix`,l.fix,"fix")}</span>`).join("")}</div>
      <div class="tw"><table class="xps-tbl"><thead><tr><th>Species</th><th>In layer</th><th class="num" title="Attenuation length (nm) of this species' electrons, same in every layer (empty = calculator)">L (nm)</th><th class="num" title="Relative bulk intensity of the pure species (only ratios matter)">I∞</th></tr></thead><tbody>${names.map((n,k)=>{ const q=ml.sp[k]||(ml.sp[k]={}); return `<tr><td>${esc(n)}</td><td>${sIn(`arx|ml|sp|${k}|layer`,q.layer!=null?q.layer:(k===0?0:N),[...ml.layers.slice(0,N).map((l,i)=>[i,l.name]),[N,"Substrate"]])}</td><td class="num">${nIn(`arx|ml|sp|${k}|L`,q.L,{cls:"w56",ph:f(arxAL(arxKE(k)),2)})}</td><td class="num">${nIn(`arx|ml|sp|${k}|Iinf`,q.Iinf,{cls:"w56",ph:"1"})}</td></tr>`; }).join("")}</tbody></table></div>
      <div class="xps-row"><button class="btn sm primary" id="xarMl">Fit layers</button></div><div id="xarMlR"></div></div></details>
      ${figCard("xarFigD","Relative depth plot",{csv:true,tip:"RDP = ln(Σ grazing / Σ normal); larger = nearer the surface. Only order and differences are meaningful (common angle factor)",extra:` <span class="xps-f" title="Angles in each of the normal / grazing sets">set ${nIn("arx|nSet",a.nSet,{cls:"w56",min:1,step:1})}</span>`})}
      ${figCard("xarFigM","MaxEnt depth profile",{csv:true,tip:"Livesey & Smith maximum-entropy profile, α by the discrepancy criterion χ² = N. Uses L and I∞ (bulk intensity) of each species from the Multilayer table. ARXPS has only a few independent depth parameters (Cumpson 1995): this is the smoothest profile consistent with the data, interfaces are smeared (≈ 0.2 nm pulled toward the surface in tests), nothing is learned below ≈ 3·L·cosθ_max, and angles > 65° are biased (no elastic scattering)",extra:` <span class="xps-f">J ${nIn("arx|J",a.J,{cls:"w56",min:5,max:60,step:1,title:"Number of layers"})}</span><span class="xps-f" title="Deepest layer top (nm); empty = 3·max L·max cosθ">z<sub>max</sub> ${nIn("arx|zmax",a.zmax,{cls:"w56",step:0.5})}</span><span class="xps-f" title="Relative uncertainty added in quadrature to the counting σ (model error: straight-line attenuation, L values, discretised layers). With counting σ alone (≈ 0.1 %) the discrepancy level χ² = N is usually not reachable even for exact synthetic data (a sharp interface between layer tops leaves ≈ 0.3 % misfit)">σ floor ${nIn("arx|sfl",a.sfl,{cls:"w56",min:0,step:0.5})} %</span><button class="btn xs primary" id="xarMe">Run</button>`})}`:""}
  </div>`;
  bindSt(el,k=>{ if(["set","coord","manual","src","region","from"].includes(k)) RES.arx=null; drawWS(); });
  $q(el,"#xarGet").onclick=()=>{ try{ arxIntensities(); arxRDP(); drawWS(); }catch(e){ note(e.message,"warn"); } };
  const wa=$q(el,"#xarWAdd"); if(wa) wa.onclick=()=>{ a.wins.push({name:"",region:rs.region,lo:+(rs.specs[0].be[0]+shOf(rs.specs[0])).toFixed(1),hi:+(rs.specs[0].be[rs.specs[0].n-1]+shOf(rs.specs[0])).toFixed(1),bg:"shirley"}); save(); drawWS(); };
  $a(el,"[data-wdel]").forEach(b=>b.onclick=()=>{ a.wins.splice(+b.dataset.wdel,1); save(); drawWS(); });
  if(!T) return;
  $q(el,"#xarICsv").onclick=()=>{ const L=[["theta_deg",...T.names.flatMap(n=>[n,"s_"+n])].map(csvE).join(",")]; T.th.forEach((t,i)=>L.push([t,...T.names.flatMap((_,k)=>[T.I[k][i],T.sI[k][i]])].join(","))); L.push("# "+T.info); csvDL("arxps_intensities",L); };
  $q(el,"#xarEst").onclick=()=>{ const Lo=arxAL(arxKE(+a.o)), Ls=arxAL(arxKE(+a.s)); if(!(Lo>0&&Ls>0)){ note("계산기 값을 확인하세요","warn"); return; } a.Lo=+Lo.toFixed(3); a.Ls=+Ls.toFixed(3); save(); drawWS(); };
  $q(el,"#xarOv").onclick=arxFitOv; $q(el,"#xarMl").onclick=arxFitMl; $q(el,"#xarMe").onclick=arxMaxent;
  const OV=RES.arx.ov, ovR=$q(el,"#xarOvR");
  if(OV&&ovR){ const r=OV.r; ovR.innerHTML=`<div class="strip"><div class="ok" title="± σ scaled by √max(χ²ν,1)"><b>${f(r.d,3)} ± ${f(r.sdScaled,3)}</b>d (nm)</div><div><b>${f(r.phi,2)}</b>χ²ν (${r.dof} dof)</div><div title="${a.fitR?"fitted":"fixed"}"><b>${f(r.Rinf,3)}${r.sRinf!=null?" ± "+f(r.sRinf,3):""}</b>R∞</div><div title="Hill-plot slope × L_o (biased when L_o ≠ L_s)"><b>${f(OV.h.d,3)}</b>Hill d (nm)</div></div>${r.warn&&r.warn.length?`<div class="note warn xs">${esc(r.warn.join(" · "))}</div>`:""}
     <details class="xps-d"><summary>Per angle</summary><div class="b tw"><table class="xps-tbl"><thead><tr><th class="num">θ</th><th class="num">R</th><th class="num">σR</th><th class="num" title="Exact single-angle solution">d(θ) nm</th></tr></thead><tbody>${T.th.map((t,i)=>`<tr><td class="num">${f(t,1)}</td><td class="num">${g(OV.R[i],4)}</td><td class="num">${g(OV.sR[i],2)}</td><td class="num">${f(r.perAngle[i],3)}</td></tr>`).join("")}</tbody></table></div></details>`; }
  const hR=bindFig(el,"xarFigR","arxps_ratio",()=>{ if(!OV) return; const L=["theta_deg,R,sR,R_model,d_single_angle_nm"]; T.th.forEach((t,i)=>L.push([t,OV.R[i],OV.sR[i],ARX_().overlayerR(OV.r.d,t,Object.assign({},OV.opt,{Rinf:OV.r.Rinf})),OV.r.perAngle[i]].join(","))); L.push(`# d=${OV.r.d} sd=${OV.r.sdScaled} Lo=${OV.opt.Lo} Ls=${OV.opt.Ls} Rinf=${OV.r.Rinf}`); csvDL("arxps_overlayer",L); });
  if(hR){ if(!OV) hR.innerHTML=`<div class="xps-empty sm">Fit</div>`; else { const tt=[]; for(let t=0;t<=Math.max(...T.th)+2;t+=0.5) tt.push(t); const q=Object.assign({},OV.opt,{Rinf:OV.r.Rinf});
    putChart(hR,{id:"xarR",h:250,x:{label:"θ from normal (deg)"},y:{label:"R = I_o / I_s",zero:true},series:[{name:"Measured",x:T.th,y:OV.R,err:OV.sR,pts:true,line:false,c:"var(--s1)"},{name:`Model d = ${f(OV.r.d,2)} nm`,x:tt,y:tt.map(t=>ARX_().overlayerR(OV.r.d,t,q)),c:"var(--s2)"}]}); } }
  const hH=bindFig(el,"xarFigH","arxps_hill"); if(hH){ if(!OV) hH.innerHTML=""; else { const h=OV.h; const xs=[0,Math.max(...h.x)*1.05]; putChart(hH,{id:"xarH",h:220,x:{label:"1 / cos θ",min:0.9},y:{label:"ln(1 + R/R∞)",zero:true},series:[{name:"Data",x:h.x,y:h.y,pts:true,line:false,c:"var(--s1)"},{name:`slope ${f(h.slope,3)} → d = ${f(h.d,2)} nm`,x:xs,y:xs.map(v=>v*h.slope),c:"var(--s3)",dash:"5 4"}]}); } }
  const MLr=RES.arx.ml, mR=$q(el,"#xarMlR");
  if(MLr&&mR){ const r=MLr.r; mR.innerHTML=`<div class="strip">${MLr.model.layers.map((l,i)=>`<div class="ok"><b>${f(r.d[i],3)}${isFinite(r.sdScaled[i])?" ± "+f(r.sdScaled[i],3):" (fixed)"}</b>${esc(l.name)} (nm)</div>`).join("")}<div><b>${f(r.phi,2)}</b>χ²ν (${r.dof} dof)</div></div>`;
    const id="xarMlF"; mR.insertAdjacentHTML("beforeend",figCard(id,"Measured vs model (÷ per-angle factor)")); const hM=bindFig(mR,id,"arxps_multilayer");
    putChart(hM,{id:"xarMl",h:230,x:{label:"θ (deg)"},y:{label:"I / G(θ)",zero:true},series:T.names.flatMap((n,k)=>[{name:n,x:T.th,y:T.th.map((_,i)=>T.I[k][i]/r.G[i]),pts:true,line:false,c:`var(--s${k%8+1})`},{x:T.th,y:T.th.map((_,i)=>r.model[i][k]),c:`var(--s${k%8+1})`,dash:"4 3"}])}); }
  const RD=RES.arx.rdp; const hD=bindFig(el,"xarFigD","arxps_rdp",()=>{ if(!RD) return; const L=["species,RDP,s_RDP,relative"]; RD.rows.forEach(q=>L.push([csvE(T.names[q.k]),q.rdp,q.s,q.rel].join(","))); csvDL("arxps_rdp",L); });
  if(hD&&RD){ const ord=RD.order; putChart(hD,{id:"xarD",h:220,x:{label:"",min:0.5,max:ord.length+0.5},xTicks:ord.map((k,i)=>({v:i+1,t:T.names[k].length>14?T.names[k].slice(0,13)+"…":T.names[k]})),y:{label:"RDP (surface ↑)"},series:[{x:ord.map((_,i)=>i+1),y:ord.map(k=>RD.rows[k].rdp),err:ord.map(k=>RD.rows[k].s),pts:true,line:false,r:4,c:"var(--s1)"}]}); }
  const ME=RES.arx.me; const hM2=bindFig(el,"xarFigM","arxps_maxent",()=>{ if(!ME) return; const r=ME.r; const L=[["z_top_nm",...T.names].map(csvE).join(",")]; Array.from(r.z).forEach((z,j)=>L.push([z,...r.c.map(c=>c[j])].join(","))); L.push(`# alpha=${r.alpha} chi2=${r.chi2} N=${r.N} notReached=${r.notReached} L=${ME.L.join(";")} Iinf=${ME.Iinf.join(";")} sigma_floor=${100*ME.fl}%`); csvDL("arxps_maxent",L); });
  if(hM2){ if(!ME) hM2.innerHTML=`<div class="xps-empty sm">Run (Worker)</div>`; else { const r=ME.r; const zz=[], st=[]; Array.from(r.z).forEach((z,j)=>{ const z2=j<r.z.length-1?r.z[j+1]:z+r.dz; zz.push(z,z2); st.push(j,j); });
    putChart(hM2,{id:"xarM",h:250,x:{label:"Depth (nm)"},y:{label:"Fraction",min:0,max:1},series:r.c.map((c,k)=>({name:T.names[k],x:zz,y:st.map(j=>c[j]),c:`var(--s${k%8+1})`}))});
    hM2.insertAdjacentHTML("beforeend",`<div class="strip"><div class="${r.notReached?"warn":"ok"}" title="Discrepancy criterion χ²(α) = N (data minus profiled angle factors)"><b>${g(r.alpha,3)}</b>α</div><div class="${r.notReached?"warn":""}"><b>${f(r.chi2,1)} / ${r.N}</b>χ² / N${r.notReached?" — not reached":""}</div><div title="${esc((r.warn||[]).join(" · "))}"><b>ⓘ</b>limits</div></div>`); } }
}

/* ============================================================
   Imaging (SPEC_multi §5) — 분광 영상 stack: 창·넓이·at% 지도, ROI 스펙트럼, binning, PCA, 화소별 NNLS 화학 상태 지도
   stack 은 IndexedDB "img" 에 구조 복제로(JSON 아님 — 256×256×100 도 그대로)
   ============================================================ */
let IMS={stacks:[],loaded:false};
async function imgLoadDB(){ if(IMS.loaded) return; IMS.loaded=true; try{ const db=await idb(); const v=await new Promise(res=>{ const t=db.transaction(DBS,"readonly").objectStore(DBS).get("img"); t.onsuccess=()=>res(t.result||null); t.onerror=()=>res(null); }); if(v&&Array.isArray(v.stacks)) IMS.stacks=v.stacks; }catch(e){} }
let _imT=null; function imgSave(){ clearTimeout(_imT); _imT=setTimeout(async()=>{ try{ const db=await idb(); await new Promise(res=>{ const t=db.transaction(DBS,"readwrite"); t.objectStore(DBS).put({stacks:IMS.stacks},"img"); t.oncomplete=res; t.onerror=res; }); }catch(e){} },400); }
function IMU(){ S.imu=S.imu||{}; defs(S.imu,{sel:null,map:"area",lo:null,hi:null,bg:"shirley",wmode:"sum",cmap:"viridis",sig:false,robust:true,roi:null,refs:{},bgsub:"shirley",nbg:0,weight:"none",show:"f0",pcaScale:"poisson",pcaK:8,pcaN:3,score:0,at:{},atShow:0,bin:2}); S.imgRefs=S.imgRefs||[]; return S.imu; }
const curStack=()=>{ const u=IMU(); let st=IMS.stacks.find(s=>s.id===u.sel); if(!st&&IMS.stacks.length){ st=IMS.stacks[0]; u.sel=st.id; } return st||null; };
function addStack(st,meta){ const o=Object.assign(st,{id:"st"+Math.random().toString(36).slice(2,8)},meta||{}); const r=IO.parseRegionName(o.region||o.name||""); o.el=o.el||r.el||""; o.orb=o.orb||r.orb||""; IMS.stacks.push(o); IMU().sel=o.id; imgSave(); return o; }
/* VAMAS MAP 블록(meta.vamas.x/y) → region 마다 stack */
function stacksFromVamas(spectra,file){ const byR={}; spectra.forEach(s=>{ (byR[s.region]=byR[s.region]||[]).push(s); }); const out=[];
  Object.entries(byR).forEach(([reg,a])=>{ const st=IMG().fromSpectra(a); const vm=a[0].meta.vamas; out.push(addStack(st,{name:`${reg} · ${(a[0].sample||file||"").slice(0,40)}`,region:reg,file,fov:vm.fov||null,src:`VAMAS ${vm.mode} ${st.nx}×${st.ny}`,missing:st.missing})); }); return out; }
/* 여러 열 텍스트: 에너지 + 화소마다 한 열(행 우선). nx 는 머리줄 없이 알 수 없어 정사각이면 √, 아니면 한 줄(나중에 nx 를 바꿈) */
function stackFromText(text,file){ const r=IO.textMultiParse(text,file); const np=r.Y.length; let nx=Math.round(Math.sqrt(np)); if(nx*nx!==np) nx=np; const ny=np/nx; const nE=r.E.length;
  const data=new Float32Array(np*nE); r.Y.forEach((y,p)=>data.set(y,p*nE)); const st=IMG().makeStack({nx,ny,E:r.E,data,cfac:null,name:file});
  return addStack(st,{name:file.replace(/^.*[\\/]/,""),region:(r.header.map(h=>IO.parseRegionName(h)).find(q=>q.el)||{label:""}).label,file,src:`Text ${np} columns`}); }
function reshapeStack(st,nx){ nx=Math.round(nx); const np=st.nx*st.ny; if(!(nx>0)||np%nx){ note(`화소 수 ${np} 가 ${nx} 로 나누어떨어지지 않습니다`,"warn"); return; } st.nx=nx; st.ny=np/nx; st.x=null; st.y=null; imgSave(); RES.img=null; RES.pca=null; drawWS(); }
const sameGrid=(a,b)=>a&&b&&a.nx===b.nx&&a.ny===b.ny;
const stackInfo=st=>`${st.nx}×${st.ny} px · ${st.nE} ch · ${f(st.E[0],1)}–${f(st.E[st.nE-1],1)} eV${st.cfac?"":" · cfac ?"}`;
function imgWin(st){ const u=IMU(); const lo=u.lo!=null?+u.lo:st.E[0], hi=u.hi!=null?+u.hi:st.E[st.nE-1]; return [Math.min(lo,hi),Math.max(lo,hi)]; }
function imgRefCands(st){ const u=IMU(); const roi=S.imgRefs.filter(r=>!st.el||!r.el||r.el===st.el).map(r=>({id:"r:"+r.id,name:r.name,kind:"ROI",ev:E=>{ const o=new Float64Array(E.length); for(let i=0;i<E.length;i++) o[i]=interp(r.E,r.y,E[i]); return o; }})); void u; return refCands(st.el,st.orb,roi); }
async function imgRun(){
  const st=curStack(); if(!st){ note("stack 을 여세요","warn"); return; } const u=IMU(); const [lo,hi]=imgWin(st); const np=st.nx*st.ny;
  try{
    if(u.map==="window"){ const v=IMG().windowImage(st,lo,hi,u.wmode); RES.img={sid:st.id,kind:"window",maps:[{name:`Σ ${f(lo,1)}–${f(hi,1)} eV`,v,unit:u.wmode==="mean"?"mean":"sum"}]}; }
    else if(u.map==="area"){ const r=await xrun("imgArea",{st,opt:{lo,hi,bg:u.bg,n:3}},{label:"Area map"}); RES.img={sid:st.id,kind:"area",maps:[{name:`Area ${f(lo,1)}–${f(hi,1)} eV (${u.bg})`,v:r.area,s:r.sArea,unit:"area"}],conv:r.converged}; }
    else if(u.map==="at"){ const items=IMS.stacks.filter(q=>sameGrid(q,st)&&(u.at[q.id]||{}).inc!==false); if(items.length<2){ note("at % 지도에는 같은 화소 격자의 stack 이 2개 이상 필요합니다","warn"); return; }
      const cfg=items.map(q=>{ const c=u.at[q.id]||{}; const lk=rsfLookup(q.el,q.orb||"",S.quant.lib,1486.6); return {q,lo:c.lo!=null?+c.lo:q.E[0],hi:c.hi!=null?+c.hi:q.E[q.nE-1],bg:c.bg||"shirley",S:c.S>0?+c.S:(lk?lk.v:null)}; });
      const bad=cfg.filter(c=>!(c.S>0)); if(bad.length){ note(`S 가 없습니다: ${bad.map(c=>c.q.name).join(", ")}`,"warn"); return; }
      const ar=await xrun("imgAreas",{items:cfg.map(c=>({st:c.q,opt:{lo:c.lo,hi:c.hi,bg:c.bg,n:3},name:c.q.name}))},{label:"at % maps"});
      const A=IMG().atomicMaps(cfg.map((c,k)=>({label:c.q.el||c.q.name,area:ar[k].area,sArea:ar[k].sArea,S:c.S})),{});
      RES.img={sid:st.id,kind:"at",maps:A.labels.map((l,k)=>({name:`${l} at %`,v:A.at[k],s:A.sAt[k],unit:"at %"}))}; }
    else if(u.map==="nnls"){ const cands=chosenRefs(imgRefCands(st),u.refs); if(!cands.length){ note("기준을 고르세요","warn"); return; }
      const [i0,i1]=IMG().idxRange(st.E,lo,hi); const Ew=st.E.slice(i0,i1+1); if(Ew.length<6){ note("창이 너무 좁습니다","warn"); return; }
      const refs=cands.map(c=>unitArea(Ew,c.ev(Ew))); const sub=u.bgsub!=="none"; const weight=sub?"none":u.weight; const K=Math.min(nWorkers(),Math.max(1,Math.ceil(np/2048)));
      const chunk=Math.ceil(np/K); const args=[]; for(let k=0;k<K;k++){ const a=k*chunk, b=Math.min(np,a+chunk); if(b<=a) continue; args.push({st:{nx:b-a,ny:1,nE:st.nE,E:st.E,data:st.data.slice(a*st.nE,b*st.nE),cfac:st.cfac},lo,hi,bgsub:sub?u.bgsub:null,refs,bg:+u.nbg||0,weight}); }
      const outs=await xpool("imgNNLS",args,{label:"State maps"}); const cat=key=>refs.map((_,k)=>{ const v=new Float64Array(np); let o=0; outs.forEach(r=>{ v.set(r[key][k],o); o+=r[key][k].length; }); return v; });
      const coef=cat("coef"), sig=cat("sig"), frac=cat("fraction"); const red=new Float64Array(np); { let o=0; outs.forEach(r=>{ red.set(r.redchi,o); o+=r.redchi.length; }); }
      const sc=Float64Array.from(red,v=>weight==="none"?Math.sqrt(v):Math.sqrt(Math.max(v,1))); const sig2=sig.map(s=>Float64Array.from(s,(v,p)=>v*sc[p]));
      const maps=[]; cands.forEach((c,k)=>{ maps.push({name:`${c.name} · fraction %`,v:frac[k],unit:"%",key:"f"+k}); maps.push({name:`${c.name} · amount`,v:coef[k],s:sig2[k],unit:"area",key:"c"+k}); }); maps.push({name:"Reduced χ²",v:red,unit:"χ²ν",key:"x"});
      RES.img={sid:st.id,kind:"nnls",maps,note:outs[0].note}; if(!maps.find(m=>m.key===u.show)) u.show="f0"; }
    else if(u.map==="pca"){ const r=await xrun("imgPCA",{st,win:[lo,hi],scale:u.pcaScale,nKeep:Math.max(2,Math.min(20,+u.pcaK||8)),nMax:Math.max(2,Math.min(12,+u.pcaK||8))},{label:"PCA"});
      RES.pca={sid:st.id,P:r.P,NF:r.NF,E:r.E,win:[lo,hi]}; RES.img={sid:st.id,kind:"pca",maps:r.P.scores.map((v,k)=>({name:`Score ${k+1} (${f(r.P.varPct[k],1)} %)`,v,unit:"score",key:"p"+k}))};
      if(r.NF&&r.NF.nNoise!=null){ u.pcaN=r.NF.nNoise; } note(`PCA: noise-limited factors ${r.NF&&r.NF.nNoise!=null?r.NF.nNoise:"?"} (MLPCA on ${r.NF&&r.NF.nSub||"?"} px) · IND ${r.P.malinowski.nIND}`,"ok"); }
    save(); drawWS();
  }catch(e){ failNote("Map",e); }
}
async function imgRecon(){ const st=curStack(), P=RES.pca; if(!st||!P||P.sid!==st.id){ note("먼저 PCA","warn"); return; } const n=Math.max(1,Math.min(+IMU().pcaN||3,P.P.U.length));
  try{ const r=await xrun("imgRecon",{st,win:P.win,P:P.P,n},{label:"Reconstruct"}); addStack(r,{name:`${st.name} · PCA${n}`,region:st.region,el:st.el,orb:st.orb,src:`PCA reconstruction (${n} factors) of ${st.name}`}); RES.img=null; note(`PCA ${n} 요인으로 다시 만든 stack 을 더했습니다 (잡음 감소)`,"ok"); drawWS(); }catch(e){ failNote("Reconstruct",e); } }
function imgBin(){ const st=curStack(); if(!st) return; const b=Math.max(2,+IMU().bin||2); if(st.nx<b||st.ny<b){ note("영상이 너무 작습니다","warn"); return; }
  const r=IMG().binStack(st,b,"sum"); addStack(r,{name:`${st.name} · bin${b}`,region:st.region,el:st.el,orb:st.orb,src:`${b}×${b} sum of ${st.name}`}); RES.img=null; drawWS(); }
function roiCompute(){ const st=curStack(), u=IMU(); if(!st||!u.roi){ RES.roi=null; return; } const q=u.roi; const r={x0:Math.max(0,Math.min(q.x0,q.x1)),x1:Math.min(st.nx-1,Math.max(q.x0,q.x1)),y0:Math.max(0,Math.min(q.y0,q.y1)),y1:Math.min(st.ny-1,Math.max(q.y0,q.y1))};
  RES.roi=Object.assign(IMG().roiSpectrum(st,r),{rect:r,sid:st.id}); }
function roiToSpectrum(){ const R=RES.roi, st=curStack(); if(!R||!st) return; const r=R.rect; const nm=`${st.region||"ROI"}`;
  const sp=IO.mkSpec({name:nm,region:st.region||"ROI",sample:`${st.name} ROI x${r.x0}-${r.x1} y${r.y0}-${r.y1}`,file:st.file||st.name,format:"ROI",be:R.E,y:R.y,c:st.cfac?st.cfac*R.n:null,cSrc:st.cfac?`stack cfac × ${R.n} px`:null,meta:{roi:{stack:st.name,rect:r,n:R.n}}});
  S.specs.push(sp); pr(sp); smpOf(sp); S.sel=sp.uid; save(); note(`ROI 평균 스펙트럼을 목록에 더했습니다 (${R.n} px) — Spectra 에서 피팅·기준 저장`,"ok",{label:"Open in Spectra",fn:()=>{ S.ui.ws="spectra"; S.ui.tab="fit"; save(); rerender(); }}); }
function roiAsRef(){ const R=RES.roi, st=curStack(); if(!R||!st) return; const [lo,hi]=imgWin(st); const [i0,i1]=IMG().idxRange(R.E,lo,hi); const E=R.E.slice(i0,i1+1), y=R.y.slice(i0,i1+1);
  const B=C.bgStatic(E,y,{type:"shirley",n:3})||C.bgShirley(E,y,3,false); const z=unitArea(E,Float64Array.from(y,(v,i)=>Math.max(0,v-B[i])));
  const r=R.rect; S.imgRefs.push({id:"r"+Math.random().toString(36).slice(2,7),name:`ROI ${S.imgRefs.length+1} (x${r.x0}–${r.x1}, y${r.y0}–${r.y1})`,el:st.el,E:Array.from(E),y:Array.from(z)}); save(); note("ROI 를 기준으로 저장했습니다 (Shirley 를 빼고 넓이 1)","ok"); drawWS(); }

function wsImaging(el){
  const u=IMU(); if(!IMS.loaded){ el.innerHTML=`<div class="card xws-empty">…</div>`; imgLoadDB().then(()=>drawWS()); return; }
  const st=curStack();
  if(!st){ el.innerHTML=`<div class="card xws-empty"><b>Imaging</b><div class="muted sm" style="margin:6px 0 10px">VAMAS MAP 블록(화소마다 스펙트럼) 또는 여러 열 텍스트(에너지 + 화소마다 한 열)를 여세요.</div><label class="btn primary sm" style="cursor:pointer">Open stack<input type="file" id="ximFin" multiple accept=".vms,.npl,.txt,.csv,.dat" style="display:none"></label></div>`;
    $q(el,"#ximFin").onchange=e=>{ loadFiles(e.target.files); e.target.value=""; }; return; }
  const [lo,hi]=imgWin(st); const same=IMS.stacks.filter(q=>sameGrid(q,st));
  const mapOpts=[["window","Window image"],["area","Peak area"],["at","at % (several regions)"],["nnls","Chemical state (NNLS)"],["pca","PCA scores"]];
  let extra="";
  if(u.map==="window") extra=sIn("imu|wmode",u.wmode,[["sum","Sum"],["mean","Mean"]]);
  if(u.map==="area") extra=sIn("imu|bg",u.bg,[["none","No bg"],["linear","Linear"],["shirley","Shirley"]],{title:"Per-pixel background; σ map from Poisson statistics (needs cfac)"});
  if(u.map==="at") extra=`<div class="tw"><table class="xps-tbl"><thead><tr><th>Stack</th><th>Window (eV)</th><th>BG</th><th>S</th></tr></thead><tbody>${same.map(q=>{ const c=u.at[q.id]||(u.at[q.id]={}); const lk=rsfLookup(q.el,q.orb||"",S.quant.lib,1486.6); return `<tr><td>${cIn(`imu|at|${q.id}|inc`,c.inc!==false,esc(q.region||q.name))}</td><td class="nw">${nIn(`imu|at|${q.id}|lo`,c.lo,{cls:"w56",ph:f(q.E[0],1)})}–${nIn(`imu|at|${q.id}|hi`,c.hi,{cls:"w56",ph:f(q.E[q.nE-1],1)})}</td><td>${sIn(`imu|at|${q.id}|bg`,c.bg||"shirley",[["linear","Lin."],["shirley","Shirley"]])}</td><td>${nIn(`imu|at|${q.id}|S`,c.S,{cls:"w56",ph:lk?g(lk.v,3):"?"})}</td></tr>`; }).join("")}</tbody></table></div>`;
  if(u.map==="nnls") extra=`${sIn("imu|bgsub",u.bgsub,[["shirley","− Shirley"],["linear","− linear"],["none","Keep bg"]],{title:"Per-pixel background subtraction before the fit (then unweighted — subtracted data are not Poisson)"})}${u.bgsub==="none"?sIn("imu|weight",u.weight,[["none","Unweighted"],["poisson","Poisson IRLS"]])+sIn("imu|nbg",u.nbg,[["0","No bg term"],["1","+ const"],["2","+ linear"]]):""}${refPickHTML("imu|refs",imgRefCands(st),u.refs,{empty:"기준 없음 — ROI 를 그리고 ‘Use as reference’, 또는 Spectra 의 피팅 성분"})}`;
  if(u.map==="pca") extra=`${sIn("imu|pcaScale",u.pcaScale,[["poisson","Poisson scaling"],["none","No scaling"]])}<span class="xps-f" title="Factors kept / tested">k ${nIn("imu|pcaK",u.pcaK,{cls:"w56",min:2,max:20,step:1})}</span>`;
  const R=RES.img&&RES.img.sid===st.id?RES.img:null;
  const maps=R?R.maps:[]; let mi=R?Math.max(0,maps.findIndex(m=>m.key===u.show)):0; if(R&&R.kind==="pca") mi=Math.max(0,Math.min(maps.length-1,+u.score||0)); if(R&&(R.kind==="at")) mi=Math.max(0,Math.min(maps.length-1,+u.atShow||0));
  const M=maps[mi]; const showSig=u.sig&&M&&M.s;
  el.innerHTML=`<div class="card xws-side">
    <div class="xps-h">Imaging</div>
    <div class="xps-row">${sIn("imu|sel",st.id,IMS.stacks.map(q=>[q.id,q.name]),{title:"Spectral image stack"})}<label class="btn xs" style="cursor:pointer" title="VAMAS MAP (.vms) or multi-column text (energy + one column per pixel)">Open<input type="file" id="ximFin" multiple accept=".vms,.npl,.txt,.csv,.dat" style="display:none"></label><button class="btn xs" id="ximDel" title="Remove this stack (press twice)">×</button></div>
    <div class="muted xs" title="${esc(st.src||"")}">${esc(stackInfo(st))}${st.missing?` · ${st.missing} px missing`:""}</div>
    <div class="xps-row"><span class="xps-f" title="Pixels per row (re-shape a text stack)">nx ${nIn("",st.nx,{cls:"w56",step:1,min:1,id:"ximNx"})}</span><span class="xps-f" title="Sum b×b pixels (Poisson statistics preserved) → new stack">Bin ${sIn("imu|bin",u.bin,[["2","2×2"],["3","3×3"],["4","4×4"]])}</span><button class="btn xs" id="ximBin">Bin</button>${st.cfac?"":`<span class="xps-f" title="Counts per data unit (dwell × scans). Needed for σ maps and Poisson weights">cfac ${nIn("",null,{cls:"w56",ph:"1",id:"ximCf"})}</span>`}</div>
    <hr class="xps-hr">
    <div class="xps-row"><span class="xws-k">Map</span>${sIn("imu|map",u.map,mapOpts)}</div>
    <div class="xps-row"><span class="xws-k">Window</span>${nIn("imu|lo",u.lo,{cls:"w64",step:0.1,ph:f(st.E[0],1),title:"Low BE"})}–${nIn("imu|hi",u.hi,{cls:"w64",step:0.1,ph:f(st.E[st.nE-1],1),title:"High BE"})}</div>
    <div class="xps-row">${extra}</div>
    <div class="xps-row"><button class="btn sm primary" id="ximRun" title="Long maps run in Workers (Stop to cancel)">Make map</button>${u.map==="pca"&&RES.pca&&RES.pca.sid===st.id?`<span class="xps-f" title="Rebuild the stack from n factors (noise reduction) → new stack">n ${nIn("imu|pcaN",u.pcaN,{cls:"w56",min:1,step:1})}</span><button class="btn xs" id="ximRec">Reconstruct</button>`:""}</div>
    <hr class="xps-hr">
    <div class="xps-row"><span class="xws-k">Display</span>${sIn("imu|cmap",u.cmap,[["viridis","Viridis"],["magma","Magma"],["gray","Gray"]])}${cIn("imu|robust",u.robust,"1–99 %",{title:"Colour range from the 1st–99th percentile (otherwise min–max)"})}${cIn("imu|sig",u.sig,"σ map",{title:"Show the 1σ uncertainty map instead of the value (area, at %, NNLS amount)"})}</div>
    ${R&&maps.length>1?`<div class="xps-row"><span class="xws-k">Show</span>${R.kind==="pca"?sIn("imu|score",mi,maps.map((m,i)=>[i,m.name])):R.kind==="at"?sIn("imu|atShow",mi,maps.map((m,i)=>[i,m.name])):sIn("imu|show",M?M.key:"",maps.map(m=>[m.key,m.name]))}</div>`:""}
  </div>
  <div class="card xws-main"><div class="xws-pg" id="xwsProgSlot"></div>
    ${M?figCard("ximFigM",esc(showSig?M.name+" · σ":M.name),{csv:true,tip:"Drag on the map to select a rectangular ROI; its mean spectrum appears below",below:mapHTML({nx:st.nx,ny:st.ny,cmap:u.cmap,unit:M.unit,aspect:st.fov&&st.fov[0]>0&&st.fov[1]>0?st.fov[0]/st.fov[1]:st.nx/st.ny})}):`<div class="xps-empty sm">Make map</div>`}
    ${R&&R.kind==="nnls"&&R.note?`<div class="muted xs">${esc(R.note)}</div>`:""}
    ${RES.pca&&RES.pca.sid===st.id&&R&&R.kind==="pca"?figCard("ximFigP","Eigenvalues",{csv:true,tip:"Scree plot (log) and the MLPCA noise-limited factor count (Poisson, strided pixel subset)"})+`<div id="ximPcaT"></div>`:""}
    ${figCard("ximFigR","ROI spectrum",{csv:true,tip:"Mean spectrum of the ROI pixels ± 1σ of the mean (Poisson)",extra:` <button class="btn xs" id="ximRoiS" ${RES.roi?"":"disabled"} title="Add the ROI mean spectrum to the Spectra list (fit it, store it as a reference …)">Add to spectra</button><button class="btn xs" id="ximRoiR" ${RES.roi?"":"disabled"} title="Use the ROI spectrum (Shirley removed over the map window, unit area) as a reference for state maps">Use as reference</button>`})}
  </div>`;
  bindSt(el,k=>{ if(["sel"].includes(k)){ RES.img=null; RES.roi=null; u.roi=null; } if(["cmap","robust","sig","show","score","atShow"].includes(k)){ drawWS(); return; } drawWS(); });
  $q(el,"#ximFin").onchange=e=>{ loadFiles(e.target.files,{img:true}); e.target.value=""; };
  $q(el,"#ximDel").onclick=()=>{ if(!confirmInline("ximDel","다시 누르면 stack 을 지웁니다")) return; IMS.stacks=IMS.stacks.filter(q=>q.id!==st.id); u.sel=null; RES.img=null; RES.roi=null; RES.pca=null; imgSave(); save(); drawWS(); };
  $q(el,"#ximNx").onchange=e=>reshapeStack(st,+e.target.value);
  const cf=$q(el,"#ximCf"); if(cf) cf.onchange=()=>{ const v=+cf.value; if(v>0){ st.cfac=v; imgSave(); drawWS(); } };
  $q(el,"#ximBin").onclick=imgBin; $q(el,"#ximRun").onclick=imgRun; const rc=$q(el,"#ximRec"); if(rc) rc.onclick=imgRecon;
  $q(el,"#ximRoiS").onclick=roiToSpectrum; $q(el,"#ximRoiR").onclick=roiAsRef;
  if(M){ const card=$q(el,"#ximFigM"); const box=$q(card,".xim-box"); const m={nx:st.nx,ny:st.ny,v:showSig?M.s:M.v,cmap:u.cmap,unit:M.unit,robust:u.robust};
    paintMap(box,m); card._map=m; bindFig(el,"ximFigM",`map_${M.name}`,()=>{ const L=["x_index,y_index,x,y,value,sigma"]; for(let y=0;y<st.ny;y++) for(let x=0;x<st.nx;x++){ const p=y*st.nx+x; L.push([x,y,st.x?st.x[x]:"",st.y?st.y[y]:"",M.v[p],M.s?M.s[p]:""].join(",")); } L.unshift("# "+M.name+" · "+st.name); csvDL("map_"+M.name,L); });
    bindROI(box,st,M); }
  const hP=bindFig(el,"ximFigP","image_pca",()=>{ const P=RES.pca.P; const L=["n,eigenvalue,var_pct"]; Array.from(P.eig).forEach((e,i)=>L.push([i+1,e,P.varPct[i]].join(","))); csvDL("image_pca",L); });
  if(hP&&RES.pca){ const P=RES.pca.P, NF=RES.pca.NF; const K=Math.min(P.eig.length,15); const nn=NF&&NF.nNoise;
    putChart(hP,{id:"ximP",h:220,x:{label:"Factor",min:0.5,max:K+0.5},xTicks:[...Array(K)].map((_,i)=>({v:i+1,t:String(i+1)})),y:{label:"Eigenvalue",log:true},vlines:nn!=null?[{x:nn+0.5,t:"noise",c:"var(--warn)"}]:[],series:[{x:[...Array(K)].map((_,i)=>i+1),y:Array.from(P.eig).slice(0,K),pts:true,c:"var(--s1)"}]});
    const t=$q(el,"#ximPcaT"); if(t) t.innerHTML=`<div class="strip"><div class="${nn==null?"warn":"ok"}" title="MLPCA noise test on ${NF&&NF.nSub||"?"} strided pixels (Poisson; SPEC_multi §5)"><b>${nn==null?(NF&&NF.error?"—":"?"):nn}</b>noise-limited n</div><div><b>${P.malinowski.nIND}</b>IND</div><div><b>${f(P.varPct[0],1)} %</b>factor 1</div></div>${NF&&NF.error?`<div class="note warn xs">${esc(NF.error)}</div>`:""}`; }
  const hR=bindFig(el,"ximFigR","roi_spectrum",()=>{ const r=RES.roi; if(!r) return; const L=["BE_eV,mean,sigma_mean,sum"]; for(let i=0;i<r.E.length;i++) L.push([r.E[i],r.y[i],r.sigma?r.sigma[i]:"",r.sum[i]].join(",")); L.unshift(`# ROI ${JSON.stringify(r.rect)} n=${r.n} · ${st.name}`); csvDL("roi_spectrum",L); });
  if(hR){ const r=RES.roi&&RES.roi.sid===st.id?RES.roi:null; if(!r) hR.innerHTML=`<div class="xps-empty sm">지도 위를 끌어 ROI 를 고르세요</div>`;
    else putChart(hR,{id:"ximR",h:230,x:{label:"Binding energy (eV)",rev:true},y:{label:"Mean intensity"},vlines:[{x:lo},{x:hi}],series:[{name:`ROI ${r.n} px`,x:r.E,y:r.y,band:r.sigma?{lo:Array.from(r.y,(v,i)=>v-r.sigma[i]),hi:Array.from(r.y,(v,i)=>v+r.sigma[i])}:null,c:"var(--s1)"}]}); }
}
/* ROI: 지도 위에서 끌기(마우스·터치) → 사각형. 끌지 않고 움직이면 그 화소의 값 */
function bindROI(box,st,M){ const cv=$q(box,"canvas.xim"), ov=$q(box,".xim-roi"), tip=$q(box,".xim-tip"), u=IMU(); let d0=null;
  const pix=e=>{ const r=cv.getBoundingClientRect(); return {x:Math.max(0,Math.min(st.nx-1,Math.floor((e.clientX-r.left)/r.width*st.nx))),y:Math.max(0,Math.min(st.ny-1,Math.floor((e.clientY-r.top)/r.height*st.ny)))}; };
  const showRect=q=>{ if(!q){ ov.hidden=true; return; } const x0=Math.min(q.x0,q.x1), x1=Math.max(q.x0,q.x1), y0=Math.min(q.y0,q.y1), y1=Math.max(q.y0,q.y1); ov.hidden=false; Object.assign(ov.style,{left:(100*x0/st.nx)+"%",top:(100*y0/st.ny)+"%",width:(100*(x1-x0+1)/st.nx)+"%",height:(100*(y1-y0+1)/st.ny)+"%"}); };
  showRect(u.roi&&RES.roi&&RES.roi.sid===st.id?u.roi:null);
  cv.onpointerdown=e=>{ e.preventDefault(); const p=pix(e); d0=p; try{ cv.setPointerCapture(e.pointerId); }catch(_){} showRect({x0:p.x,y0:p.y,x1:p.x,y1:p.y}); };
  cv.onpointermove=e=>{ const p=pix(e); if(d0){ showRect({x0:d0.x,y0:d0.y,x1:p.x,y1:p.y}); return; } const i=p.y*st.nx+p.x; tip.textContent=`x ${p.x}, y ${p.y}: ${g(M.v[i],4)}${M.s?" ± "+g(M.s[i],2):""} ${M.unit||""}`; };
  cv.onpointerup=e=>{ if(!d0) return; const p=pix(e); u.roi={x0:d0.x,y0:d0.y,x1:p.x,y1:p.y}; d0=null; roiCompute(); save(); drawWS(); };
  cv.onpointercancel=()=>{ d0=null; };
}

/* ============================================================
   Processing (SPEC_proc) — 한 스펙트럼: 디컨볼루션 · 다시 표본화 · 평활/미분 · Fermi edge · 기준 봉우리 PSF
   결과는 그래프 위에 겹쳐 그림(Process 탭일 때). 디컨볼루션 결과는 보기·시작값용 — 정량은 원자료 피팅으로
   ============================================================ */
function PRS(){ S.prc=S.prc||{}; return defs(S.prc,{method:"jansson",psf:"gauss",fwhm:0.8,wG:0.6,wL:0.3,mpsf:null,stored:null,iter:null,mc:6,rstep:null,rmeth:"auto",sm:"sg",win:7,ord:2,smF:0.3,der:0,T:300,wLref:0.3}); }
/* 지금 스펙트럼의 피팅 창(측정 BE) — 고른 격자가 아니면 같은 평균 간격으로 다시 표본화(자동) */
function procData(sp){ const p=pr(sp); const R=region(sp,p); if(!R||R.E.length<16) throw new Error("피팅 창 안 점이 16개 미만입니다");
  let E=Float64Array.from(R.E), y=Float64Array.from(R.y), note_=""; const gi=PRC().gridInfo(E);
  if(!gi.uniform){ const r=PRC().resample(E,y,{step:gi.step,method:"auto"}); const ok=[]; for(let i=0;i<r.E.length;i++) if(isFinite(r.y[i])) ok.push(i); E=Float64Array.from(ok.map(i=>r.E[i])); y=Float64Array.from(ok.map(i=>r.y[i])); note_=`고르지 않은 격자 → ${f(gi.step,3)} eV 로 다시 표본화(${r.method})`; }
  const v=sp.c?Float64Array.from(y,t=>Math.max(t,1/sp.c)/sp.c):null; return {E,y,v,note:note_,step:Math.abs(E[1]-E[0])}; }
function psfSpec(){ const P=PRS();
  if(P.psf==="gauss") return {type:"gauss",fwhm:+P.fwhm};
  if(P.psf==="voigt") return {type:"voigt",wG:+P.wG,wL:+P.wL};
  if(P.psf==="stored"){ if(!P.stored) throw new Error("저장된 PSF 가 없습니다 — Fermi edge 또는 기준 봉우리에서 ‘Use as PSF’"); return P.stored; }
  if(P.psf==="measured"){ const s=byUid(P.mpsf); if(!s) throw new Error("PSF 로 쓸 스펙트럼을 고르세요"); const R=region(s,pr(s)); if(!R||R.E.length<5) throw new Error("PSF 스펙트럼의 창이 너무 좁습니다"); return {type:"measured",E:Array.from(R.E),y:Array.from(R.y)}; }
  throw new Error("PSF?"); }
async function procDeconv(sp){ let d; try{ d=procData(sp); }catch(e){ note(e.message,"warn"); return; } const P=PRS(); let psf; try{ psf=psfSpec(); }catch(e){ note(e.message,"warn"); return; }
  const opt={method:P.method,psf,noise:d.v||"auto",mc:Math.max(0,Math.min(40,Math.round(+P.mc||0))),seed:12345}; if(+P.iter>0) opt.iter=Math.round(+P.iter);
  try{ const r=await xrun("deconv",{E:d.E,y:d.y,opt},{label:"Deconvolution"}); const X=PRES[sp.uid]=PRES[sp.uid]||{}; X.dec={E:r.E,o:r.o,sig:r.noise?r.noise.sigma:null,k:r.k,stop:r.stop,res:r.resolution,ng:r.noise?r.noise.gain:null,rough:r.noise?r.noise.roughness:null,area:r.areaRatio,edge:r.edgePoints,warn:r.warnings.concat(d.note?[d.note]:[]),psf:r.kernel,method:r.method,ver:pr(sp).ver};
    drawPlot(); renderPanel(); note(`Deconvolution: k = ${r.k} (${r.stop}), resolution ×${f(r.resolution.gain,2)}${r.noise?`, noise ×${f(r.noise.gain,2)}`:""}`,r.warnings.length?"warn":"ok"); }catch(e){ failNote("Deconvolution",e); } }
function procResample(sp,add){ let d; try{ d=procData(sp); }catch(e){ note(e.message,"warn"); return; } const P=PRS(); const step=+P.rstep>0?+P.rstep:d.step;
  let r; try{ r=PRC().resample(d.E,d.y,{step,method:P.rmeth}); }catch(e){ note(e.message,"warn"); return; } const ok=[]; for(let i=0;i<r.E.length;i++) if(isFinite(r.y[i])) ok.push(i);
  const E=Float64Array.from(ok.map(i=>r.E[i])), y=Float64Array.from(ok.map(i=>r.y[i]));
  if(add){ const ns=IO.mkSpec({name:sp.name,region:sp.region,sample:sp.sample,file:sp.file,format:sp.format,hv:sp.hv,hvSrc:sp.hvSrc,be:E,y,c:null,meta:Object.assign({},sp.meta,{resampled:{step,method:r.method,from:sp.uid}}),warn:[`Resampled ${f(d.step,3)} → ${f(step,3)} eV (${r.method}); output points are correlated`]});
    ns.region=sp.region+` (Δ${+step.toFixed(3)})`; S.specs.push(ns); pr(ns); S.sel=ns.uid; save(); rerender(); note(`다시 표본화한 스펙트럼을 더했습니다 (${r.method}, ${E.length} 점) — 점끼리 상관이 있어 χ² 가 작게 나옵니다`,"ok"); return; }
  (PRES[sp.uid]=PRES[sp.uid]||{}).rs={E,y,method:r.method,warn:r.warnings,step}; drawPlot(); renderPanel(); }
function procSmooth(sp){ let d; try{ d=procData(sp); }catch(e){ note(e.message,"warn"); return; } const P=PRS(); let y;
  try{ const der=+P.der||0; if(der>0) y=PRC().derivative(d.E,d.y,{window:Math.max(5,(Math.round(+P.win)|1)),order:Math.max(der,Math.min(+P.ord||2,6)),deriv:der});
    else y=P.sm==="gauss"?PRC().smooth(d.E,d.y,{method:"gauss",fwhmEV:+P.smF}):PRC().smooth(d.E,d.y,{method:"sg",window:Math.max(3,(Math.round(+P.win)|1)),order:Math.max(0,Math.min(+P.ord||2,6))}); }catch(e){ note(e.message,"warn"); return; }
  (PRES[sp.uid]=PRES[sp.uid]||{}).sm={E:d.E,y,der:+P.der||0}; drawPlot(); renderPanel(); }
function procFermi(sp){ let d; try{ d=procData(sp); }catch(e){ note(e.message,"warn"); return; } const P=PRS();
  try{ const r=PRC().fermiFit(d.E,d.y,{T:+P.T||300,axis:"BE",noise:d.v||"auto"}); (PRES[sp.uid]=PRES[sp.uid]||{}).fermi={E:d.E,model:r.model,EF:r.EF,res:r.resolution,red:r.redchi2,warn:r.warnings}; drawPlot(); renderPanel();
    note(`Fermi edge: E_F ${f(r.EF.value+shOf(sp),3)} ± ${f(r.EF.sigma,3)} eV, resolution ${f(r.resolution.fwhm,3)} ± ${f(r.resolution.sigma,3)} eV`,"ok"); }catch(e){ note("Fermi fit: "+e.message,"warn"); } }
function procPeak(sp){ let d; try{ d=procData(sp); }catch(e){ note(e.message,"warn"); return; } const P=PRS();
  try{ const r=PRC().peakPSF(d.E,d.y,{wL:+P.wLref||0,noise:d.v||"auto"}); (PRES[sp.uid]=PRES[sp.uid]||{}).peak={E:d.E,model:r.model,p:r.p,s:r.sigma,red:r.redchi2,warn:r.warnings}; drawPlot(); renderPanel(); }catch(e){ note("Reference peak: "+e.message,"warn"); } }
function usePSF(o,src){ const P=PRS(); P.stored=Object.assign({},o,{src}); P.psf="stored"; save(); renderPanel(); note(`PSF: ${src}`,"ok"); }
function panelProc(el,sp){
  const P=PRS(), X=PRES[sp.uid]||{}, sh=shOf(sp), p=pr(sp); const specs=S.specs.filter(s=>s.uid!==sp.uid&&!(s.reg&&s.reg.kind==="survey"));
  const D=X.dec, stale=D&&D.ver!==p.ver;
  $a(el,"details[id]").forEach(d=>DOPEN[d.id]=d.open);
  el.innerHTML=`<div class="xps-h">Processing <span class="muted xs">· fit range ${f(p.lo+sh,2)}–${f(p.hi+sh,2)} eV</span></div>
  <details class="xps-d" id="xprDecD" open><summary title="Iterative deconvolution with a noise-based stop (discrepancy principle, SPEC_proc §5). Deconvolved spectra are for display, peak finding and start values only — take positions and areas from a fit of the raw data">Deconvolution <span class="badge wr" title="Display / start values only. Partial restoration pushes blended maxima apart (≈ +0.08 eV for a 0.7 eV doublet) and splits areas with bias; total area is conserved">display only</span></summary><div class="b">
    <div class="xps-row">${sIn("prc|method",P.method,[["jansson","Jansson"],["rl","Richardson–Lucy"],["landweber","Landweber"]],{title:"Jansson (reblurred, bounded), Richardson–Lucy (Poisson ML, positive, flux-conserving), projected Landweber"})}
      ${sIn("prc|psf",P.psf,[["gauss","Gaussian PSF"],["voigt","Voigt PSF"],["measured","Measured (spectrum)"],["stored","From Fermi edge / ref. peak"]],{title:"Instrument point-spread function"})}</div>
    <div class="xps-row">${P.psf==="gauss"?`<span class="xps-f">FWHM ${nIn("prc|fwhm",P.fwhm,{cls:"w56",step:0.05,min:0})} eV</span>`:P.psf==="voigt"?`<span class="xps-f">w<sub>G</sub> ${nIn("prc|wG",P.wG,{cls:"w56",step:0.05})} w<sub>L</sub> ${nIn("prc|wL",P.wL,{cls:"w56",step:0.05})} eV</span>`:P.psf==="measured"?sIn("prc|mpsf",P.mpsf,[["","(choose spectrum)"]].concat(specs.map(s=>[s.uid,`${smpDisp(SK(s)).name} · ${s.region}`])),{title:"Its fit range, minus a linear baseline, centred → PSF"}):`<span class="muted xs">${P.stored?esc(P.stored.src||"")+` (${P.stored.type} ${f(P.stored.fwhm||P.stored.wG,3)} eV)`:"Fermi edge 또는 기준 봉우리에서 ‘Use as PSF’"}</span>`}
      <span class="xps-f" title="Fixed iteration count (empty = stop at the noise level)">k ${nIn("prc|iter",P.iter,{cls:"w56",min:1,step:1,ph:"auto"})}</span><span class="xps-f" title="Monte-Carlo re-runs for the noise band and noise gain">MC ${nIn("prc|mc",P.mc,{cls:"w56",min:0,max:40,step:1})}</span><button class="btn sm primary" id="xprDec">Deconvolve</button>${D?`<button class="btn xs" id="xprDecX">Clear</button><button class="btn xs" id="xprDecC">CSV</button>`:""}</div>
    ${D?`<div class="strip${stale?" stale":""}"><div class="${D.res.gain>1.05?"ok":""}" title="FWHM of the PSF / FWHM of the restored point source at the same k (SPEC_proc §5.5)"><b>×${f(D.res.gain,2)}</b>resolution (${f(D.res.fwhmPSF,2)} → ${f(D.res.fwhmEff,2)} eV)</div><div class="${D.ng>3?"warn":""}" title="√(Σσ_o² / Σvar): noise amplification from ${PRS().mc} Monte-Carlo re-runs"><b>${D.ng!=null?"×"+f(D.ng,2):"—"}</b>noise</div><div title="Iterations and stop rule"><b>${D.k}</b>${esc(D.stop)}</div><div title="Area out / area in"><b>${f(D.area,4)}</b>area ratio</div></div>${D.warn.length?`<div class="note warn xs" title="${esc(D.warn.join(" · "))}">${esc(D.warn[0])}${D.warn.length>1?` (+${D.warn.length-1})`:""}</div>`:""}`:""}
  </div></details>
  <details class="xps-d" id="xprFerD"><summary title="Fermi-edge fit: linear DOS × Fermi–Dirac (T) ⊗ Gaussian → E_F and instrumental resolution (SPEC_proc §6); fit range = current fit range">Fermi edge · resolution</summary><div class="b">
    <div class="xps-row"><span class="xps-f">T ${nIn("prc|T",P.T,{cls:"w56",step:1})} K</span><button class="btn sm" id="xprFer">Fit edge</button>${X.fermi?`<button class="btn xs" id="xprFerU" title="Use the fitted Gaussian width as the deconvolution PSF">Use as PSF</button>`:""}</div>
    ${X.fermi?`<div class="strip"><div class="ok"><b>${f(X.fermi.EF.value+sh,3)} ± ${f(X.fermi.EF.sigma,3)}</b>E<sub>F</sub> (eV)</div><div class="ok"><b>${f(X.fermi.res.fwhm,3)} ± ${f(X.fermi.res.sigma,3)}</b>resolution FWHM</div><div><b>${f(X.fermi.red,2)}</b>χ²ν</div></div>${X.fermi.warn.length?`<div class="note warn xs">${esc(X.fermi.warn.join(" · "))}</div>`:""}`:""}
    <div class="xps-row" style="margin-top:6px"><span class="xps-f" title="Lorentzian (lifetime) width of the reference line, fixed (literature)">w<sub>L</sub> ${nIn("prc|wLref",P.wLref,{cls:"w56",step:0.01})} eV</span><button class="btn sm" id="xprPk" title="Voigt fit of a narrow reference line (e.g. Ag 3d5/2): Gaussian part = instrument PSF (SPEC_proc §7)">Fit reference peak</button>${X.peak?`<button class="btn xs" id="xprPkU">Use as PSF</button>`:""}</div>
    ${X.peak?`<div class="strip"><div class="ok"><b>${f(X.peak.p.wG,3)} ± ${f(X.peak.s.wG,3)}</b>w<sub>G</sub> (eV)</div><div><b>${f(X.peak.p.E0+sh,3)}</b>E₀</div><div><b>${f(X.peak.red,2)}</b>χ²ν</div></div><div class="muted xs" title="${esc(X.peak.warn.join(" · "))}">ⓘ asymmetric lines: w<sub>G</sub> is an upper estimate</div>`:""}
  </div></details>
  <details class="xps-d" id="xprRsD"><summary title="Flux-conserving rebin (down-sampling), PCHIP / Lanczos / linear (up-sampling) with exact variance propagation (SPEC_proc §2)">Resample</summary><div class="b">
    <div class="xps-row"><span class="xps-f">Step ${nIn("prc|rstep",P.rstep,{cls:"w64",step:0.005,ph:f(chOf(sp),3)})} eV</span>${sIn("prc|rmeth",P.rmeth,[["auto","Auto"],["rebin","Rebin"],["pchip","PCHIP"],["lanczos","Lanczos"],["linear","Linear"]])}<button class="btn sm" id="xprRs">Preview</button><button class="btn xs" id="xprRsA" title="Add the resampled spectrum to the list (points are correlated: χ² of a fit is then too small)">Add as spectrum</button></div>
    ${X.rs&&X.rs.warn.length?`<div class="note warn xs">${esc(X.rs.warn.join(" · "))}</div>`:""}
  </div></details>
  <details class="xps-d" id="xprSmD"><summary title="Savitzky–Golay (= scipy savgol_filter) or Gaussian smoothing; SG derivatives. For display and peak finding only">Smooth · derivative</summary><div class="b">
    <div class="xps-row">${sIn("prc|der",P.der,[["0","Smooth"],["1","1st derivative"],["2","2nd derivative"]])}${+P.der?"":sIn("prc|sm",P.sm,[["sg","Savitzky–Golay"],["gauss","Gaussian"]])}${P.sm==="gauss"&&!+P.der?`<span class="xps-f">FWHM ${nIn("prc|smF",P.smF,{cls:"w56",step:0.05})} eV</span>`:`<span class="xps-f">window ${nIn("prc|win",P.win,{cls:"w56",min:3,step:2})} order ${nIn("prc|ord",P.ord,{cls:"w56",min:0,max:6,step:1})}</span>`}<button class="btn sm" id="xprSm">Show</button>${X.sm||X.rs?`<button class="btn xs" id="xprSmX">Clear</button>`:""}</div>
  </div></details>`;
  $a(el,"details[id]").forEach(d=>{ if(d.id in DOPEN) d.open=DOPEN[d.id]; d.addEventListener("toggle",()=>{ DOPEN[d.id]=d.open; }); });
  bindSt(el,()=>renderPanel());
  $q(el,"#xprDec").onclick=()=>procDeconv(sp);
  const dx=$q(el,"#xprDecX"); if(dx) dx.onclick=()=>{ delete X.dec; drawPlot(); renderPanel(); };
  const dc=$q(el,"#xprDecC"); if(dc) dc.onclick=()=>{ const L=["BE_measured_eV,BE_corrected_eV,deconvolved,sigma"]; for(let i=0;i<D.E.length;i++) L.push([D.E[i],D.E[i]+sh,D.o[i],D.sig?D.sig[i]:""].join(",")); L.unshift(`# ${D.method} k=${D.k} (${D.stop}) PSF ${D.psf.type} FWHM ${f(D.psf.fwhm,4)} eV · resolution ×${f(D.res.gain,3)} · noise ×${f(D.ng,3)} · display/start values only`); csvDL(`${sp.sample}_${sp.region}_deconvolved`,L); };
  $q(el,"#xprFer").onclick=()=>procFermi(sp); const fu=$q(el,"#xprFerU"); if(fu) fu.onclick=()=>usePSF({type:"gauss",fwhm:X.fermi.res.fwhm},`Fermi edge · ${smpDisp(SK(sp)).name} ${sp.region}`);
  $q(el,"#xprPk").onclick=()=>procPeak(sp); const pu=$q(el,"#xprPkU"); if(pu) pu.onclick=()=>usePSF({type:"gauss",fwhm:X.peak.p.wG},`Ref. peak wG · ${smpDisp(SK(sp)).name} ${sp.region}`);
  $q(el,"#xprRs").onclick=()=>procResample(sp,false); $q(el,"#xprRsA").onclick=()=>procResample(sp,true);
  $q(el,"#xprSm").onclick=()=>procSmooth(sp); const sx=$q(el,"#xprSmX"); if(sx) sx.onclick=()=>{ delete X.sm; delete X.rs; drawPlot(); renderPanel(); };
}
/* 주 그래프 위 겹쳐 그리기 (Process 탭에서만) */
function procYRange(sp){ const X=S.ui.tab==="proc"&&PRES[sp.uid]; if(!X) return null; let lo=Infinity, hi=-Infinity; const sh=shOf(sp), v=plotView(sp);
  const acc=(E,y,s)=>{ for(let i=0;i<E.length;i++){ const x=E[i]+sh; if(x<v.x0||x>v.x1||!isFinite(y[i])) continue; const e=s?s[i]:0; lo=Math.min(lo,y[i]-e); hi=Math.max(hi,y[i]+e); } };
  if(X.dec) acc(X.dec.E,X.dec.o,X.dec.sig); if(X.rs) acc(X.rs.E,X.rs.y); if(X.sm&&!X.sm.der) acc(X.sm.E,X.sm.y); if(X.fermi) acc(X.fermi.E,X.fermi.model); if(X.peak) acc(X.peak.E,X.peak.model);
  return isFinite(lo)?[lo,hi]:null; }
function procOverlay(sp,X_,Y_,y0,y1){ const X=S.ui.tab==="proc"&&PRES[sp.uid]; if(!X) return ""; const sh=shOf(sp), o=[]; const P=(E,y)=>{ let d=""; for(let i=0;i<E.length;i++){ if(!isFinite(y[i])) continue; d+=(d?"L":"M")+X_(E[i]+sh).toFixed(1)+","+Y_(y[i]).toFixed(1); } return d; };
  if(X.dec){ const D=X.dec; if(D.sig){ let a="", b=""; for(let i=0;i<D.E.length;i++) a+=(i?"L":"M")+X_(D.E[i]+sh).toFixed(1)+","+Y_(D.o[i]+D.sig[i]).toFixed(1); for(let i=D.E.length-1;i>=0;i--) b+="L"+X_(D.E[i]+sh).toFixed(1)+","+Y_(D.o[i]-D.sig[i]).toFixed(1); o.push(`<path d="${a}${b}Z" fill="var(--s7)" fill-opacity=".2" stroke="none"/>`); }
    const R=Math.min(D.edge||0,D.E.length>>1); if(R>0){ [[0,R],[D.E.length-1-R,D.E.length-1]].forEach(([i,j])=>{ const a=X_(D.E[i]+sh), b=X_(D.E[j]+sh); o.push(`<rect x="${Math.min(a,b).toFixed(1)}" y="${Y_(y1).toFixed(1)}" width="${Math.abs(b-a).toFixed(1)}" height="${(Y_(y0)-Y_(y1)).toFixed(1)}" fill="var(--faint)" opacity=".13"><title>Edge points (depend on the reflection assumption)</title></rect>`); }); }
    o.push(`<path d="${P(D.E,D.o)}" fill="none" stroke="var(--s7)" stroke-width="1.8"/>`); }
  if(X.rs) X.rs.E.forEach((e,i)=>o.push(`<circle cx="${X_(e+sh).toFixed(1)}" cy="${Y_(X.rs.y[i]).toFixed(1)}" r="2.2" fill="var(--s3)"/>`));
  if(X.sm){ if(!X.sm.der) o.push(`<path d="${P(X.sm.E,X.sm.y)}" fill="none" stroke="var(--s2)" stroke-width="1.6"/>`);
    else { let mx=0; X.sm.y.forEach(v=>mx=Math.max(mx,Math.abs(v))); const mid=(y0+y1)/2, k=mx>0?0.45*(y1-y0)/mx:1; o.push(`<line x1="${X_(X.sm.E[0]+sh)}" x2="${X_(X.sm.E[X.sm.E.length-1]+sh)}" y1="${Y_(mid)}" y2="${Y_(mid)}" stroke="var(--s5)" stroke-width=".8" stroke-dasharray="2 3"/><path d="${P(X.sm.E,Array.from(X.sm.y,v=>mid+k*v))}" fill="none" stroke="var(--s5)" stroke-width="1.5" stroke-dasharray="5 3"/>`); } }
  if(X.fermi) o.push(`<path d="${P(X.fermi.E,X.fermi.model)}" fill="none" stroke="var(--s4)" stroke-width="2"/><line x1="${X_(X.fermi.EF.value+sh)}" x2="${X_(X.fermi.EF.value+sh)}" y1="${Y_(y1)}" y2="${Y_(y0)}" stroke="var(--s4)" stroke-dasharray="4 3"/>`);
  if(X.peak) o.push(`<path d="${P(X.peak.E,X.peak.model)}" fill="none" stroke="var(--s4)" stroke-width="2"/>`);
  return o.join(""); }
function procLegend(sp){ const X=S.ui.tab==="proc"&&PRES[sp.uid]; if(!X) return ""; const it=[];
  if(X.dec) it.push(`<span title="Display / start values only"><i style="background:var(--s7)"></i>Deconvolved ± σ (${esc(X.dec.method)}, k ${X.dec.k})</span>`); if(X.rs) it.push(`<span><i style="background:var(--s3);border-radius:50%"></i>Resampled (${esc(X.rs.method)})</span>`);
  if(X.sm) it.push(`<span><i style="background:var(--s5)"></i>${X.sm.der?X.sm.der+(X.sm.der===1?"st":"nd")+" derivative (scaled)":"Smoothed"}</span>`); if(X.fermi) it.push(`<span><i style="background:var(--s4)"></i>Fermi-edge model</span>`); if(X.peak) it.push(`<span><i style="background:var(--s4)"></i>Voigt reference fit</span>`);
  return it.length?`<div class="xps-legend">${it.join("")}</div>`:""; }

/* ============================================================
   작업 공간 전환: Spectra (기존) · Series · ARXPS · Imaging
   ============================================================ */
const DOPEN={};   /* 작업 공간의 접는 칸(details[id]) 열림 상태 — 다시 그려도 유지 */
const WSL=[["spectra","Spectra"],["series","Series"],["arxps","ARXPS"],["img","Imaging"]];
function drawWsSw(){ const el=ROOT&&$q(ROOT,"#xpsWs"); if(!el) return; const w=S.ui.ws||"spectra";
  el.innerHTML=WSL.map(([k,t])=>`<button data-ws="${k}" class="${w===k?"active":""}" title="${esc({spectra:"Single spectra: fit, quantify, process",series:"Depth profiles and other multi-spectrum sets: composition vs level, PCA, target testing, MCR-ALS, reference fit",arxps:"Angle-resolved XPS: overlayer thickness, multilayer, relative depth, MaxEnt profile",img:"Spectral images: window / area / at % maps, ROI spectra, PCA, chemical-state maps"}[k])}">${t}</button>`).join("");
  $a(el,"button").forEach(b=>b.onclick=()=>{ S.ui.ws=b.dataset.ws; save(); rerender(); }); }
function placeShared(){ const ws=S.ui.ws||"spectra"; const pg=$q(ROOT,"#xpsProg"), mg=$q(ROOT,"#xpsMsg"); if(!pg||!mg) return;
  const slot=ws==="spectra"?$q(ROOT,".xps-center"):$q(ROOT,"#xwsProgSlot"); if(!slot) return; if(pg.parentNode!==slot) slot.appendChild(pg); if(mg.parentNode!==slot) slot.appendChild(mg); }
function drawWS(){ if(!ROOT) return; const ws=S.ui.ws||"spectra"; const main=$q(ROOT,".xps-main"), box=$q(ROOT,"#xpsWsBox"); if(!main||!box) return;
  main.hidden=ws!=="spectra"; box.hidden=ws==="spectra";
  // 진행 줄·알림 줄을 지금 화면으로 옮긴 뒤 다시 그린다 (다시 그리면 자리가 없어지므로 먼저 빼 둔다)
  const pg=$q(ROOT,"#xpsProg"), mg=$q(ROOT,"#xpsMsg"); const center=$q(ROOT,".xps-center"); if(pg&&center&&pg.parentNode!==center){ center.appendChild(pg); center.appendChild(mg); }
  if(ws!=="spectra"){ const y=window.scrollY; $a(box,"details[id]").forEach(d=>DOPEN[d.id]=d.open); try{ if(ws==="series") wsSeries(box); else if(ws==="arxps") wsArxps(box); else wsImaging(box); }catch(e){ console.error(e); box.innerHTML=`<div class="card"><div class="note warn">화면 오류: ${esc(e.message||e)}</div></div>`; }
    $a(box,"details[id]").forEach(d=>{ if(d.id in DOPEN) d.open=DOPEN[d.id]; d.addEventListener("toggle",()=>{ DOPEN[d.id]=d.open; }); }); window.scrollTo(0,y); }
  else box.innerHTML="";
  placeShared(); drawMsg(); }

window.XPS_UI={render,rerender,_state:()=>S,_setState:o=>{ S=Object.assign(defState(),o); rerender(); },_load:loadFiles,_cur:cur,_pr:pr,_fit:runFit,_mc:runMC,_auto:autoSuggest,_shiftC1s:()=>calibrateC1s(cur()),_calib:sp=>calibrateC1s(sp||cur()),_saveTpl:saveTpl,_loadTpl:loadTpl,_tplList:tplList,_quantRows:k=>quantRows(k),_interf:k=>interferences(k,quantRows(k)),_augCalc:e=>augCalc(e),_applyAxis:(u,a,hv)=>applyAxis(byUid(u),a,hv),_axisCheck:u=>axisCheck(byUid(u)),_rename:renameSample,_move:(u,k)=>moveSpec(byUid(u),k),_geo:()=>GEO,_storeRef:sp=>storeRef(sp||cur()),_refs:()=>REFS(),_addRef:(id,lcf)=>{ const t=REFS().find(x=>x.id===id); if(t){ const sp=cur(), p=pr(sp); pushHist(p,"기준 추가"); addRefComp(sp,t,!!lcf); touch(p); rerender(); } },_lcfAll:()=>lcfAll(cur()),_batch:(uids,fit)=>{ S.ui.batSel={}; uids.forEach(u=>S.ui.batSel[u]=true); return runBatch(cur(),fit!==false); },_addKnot:be=>addKnot(be),
  /* v1.8 Series · ARXPS · Imaging · Processing (시험용 손잡이) */
  _res:()=>RES,_pres:()=>PRES,_ims:()=>IMS,_drawWS:()=>drawWS(),_ser:{compute:q=>serCompute(q),pca:()=>serPCA(),resolve:()=>serResolve(),build:()=>buildSet(SER()),tt:()=>ttRows(),profiles:()=>resolveProfiles(),applyFit:r=>{ const D=buildSet(SER()); const rs=D&&D.regions.find(x=>x.region===r); return rs?serApplyFit(rs):null; }},
  _arx:{get:()=>arxIntensities(),ov:()=>arxFitOv(),ml:()=>arxFitMl(),rdp:()=>arxRDP(),me:()=>arxMaxent(),al:k=>arxAL(arxKE(k))},
  _img:{run:()=>imgRun(),recon:()=>imgRecon(),bin:()=>imgBin(),roi:r=>{ IMU().roi=r; roiCompute(); drawWS(); return RES.roi; },roiRef:()=>roiAsRef(),roiSpec:()=>roiToSpectrum(),cands:()=>{ const st=curStack(); return st?imgRefCands(st).map(c=>({id:c.id,name:c.name,kind:c.kind})):[]; }},
  _proc:{dec:sp=>procDeconv(sp||cur()),fermi:sp=>procFermi(sp||cur()),peak:sp=>procPeak(sp||cur()),resample:(sp,add)=>procResample(sp||cur(),add),smooth:sp=>procSmooth(sp||cur())},
  _refCands:(el,orb)=>refCands(el,orb).map(c=>({id:c.id,name:c.name,kind:c.kind}))};
})();
