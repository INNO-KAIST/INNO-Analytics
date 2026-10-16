/* ============================================================
   INNO Analytics — Image 탭 UI (ImageJ 를 본뜬 영상 분석)
   계산 = IMG_CORE (Web Worker), 파일 = IMG_IO. 이 파일들은 inject_all.py 가 하나의 IIFE 로 묶는다.
   ============================================================ */
const IC=window.IMG_CORE, IO=window.IMG_IO;
const esc=s=>String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const $q=(r,s)=>r&&r.querySelector(s), $a=(r,s)=>r?[...r.querySelectorAll(s)]:[];
const clamp=(v,a,b)=>v<a?a:v>b?b:v;
const fmt=(v,d)=>{ if(v==null||v!==v) return "—"; if(!isFinite(v)) return v>0?"∞":"−∞"; const a=Math.abs(v); if(a!==0&&(a<1e-3||a>=1e6)) return v.toExponential((d==null?3:d)); return (+v.toFixed(d==null?3:d)).toString(); };
const g4=v=>{ if(v==null||v!==v) return "—"; if(!isFinite(v)) return String(v); const a=Math.abs(v); if(a!==0&&(a<1e-3||a>=1e6)) return v.toExponential(3); return String(+v.toPrecision(5)); };
const uid=()=>Math.random().toString(36).slice(2,9)+Date.now().toString(36).slice(-3);
const sleep=ms=>new Promise(r=>setTimeout(r,ms));

/* ---------------- 상태 ----------------
   S 는 JSON 으로 저장되는 부분(영상 메타데이터·ROI·결과·매크로·화면). 화소는 PX(id → Float32Array[]) 에 따로 둔다 */
function defState(){ return {v:1,images:[],sel:null,rois:{},roiSel:{},results:[],particles:{},summary:[],plots:[],
  measure:{cols:["Area","Mean","StdDev","Min","Max","X","Y","Perim.","Feret","MinFeret","Circ.","AR","Round","Solidity","ECD"]},
  macro:{rec:true,steps:[]},
  ui:{panel:"adjust",tool:"rect",bottom:"results",bottomOpen:true,proc:"gaussian",views:{},thr:{method:"Default",dark:true,manual:false,lower:0,upper:255,show:true},
      cthr:{space:"HSB",min:[0,0,0],max:[255,255,255],pass:[true,true,true],method:"Default",dark:false},
      pa:{minSize:0,maxSize:"",unitPx:false,minCirc:0,maxCirc:1,exclude:false,holes:false,redirect:"",show:"outlines",addRoi:false},
      dist:{col:"ECD",bins:0,log:false},scalebar:{on:true,pos:"br",len:0},lw:1}}; }
let S=defState(); const PX=new Map(); let ROOT=null, LOADED=false;
const UNDO=new Map(); // id → [{planes, meta}]
const cur=()=>S.images.find(i=>i.id===S.sel)||null;
const byId=id=>S.images.find(i=>i.id===id)||null;
const planesOf=im=>PX.get(im.id)||[];
const pIndex=(im,c,z,t)=>c+im.nc*(z+im.nz*t);
function planeAt(im,c,z,t){ const P=planesOf(im); const k=pIndex(im,c==null?im.c:c,z==null?im.z:z,t==null?im.t:t); return {w:im.w,h:im.h,bit:im.rgb?8:im.bit,data:P[k]}; }
const curPlane=()=>{ const im=cur(); return im?planeAt(im):null; };
const isStack=im=>im&&(im.nz*im.nt>1);
const calOf=im=>im&&im.cal&&im.cal.pw>0?im.cal:{pw:1,ph:1,unit:"px"};
const unitSq=im=>{ const u=calOf(im).unit; return u==="px"?"px²":u+"²"; };

/* ---------------- 저장 (IndexedDB) ---------------- */
const DBN="innoanalytics_image", ST="kv";
let _db=null; function idb(){ if(_db) return _db; _db=new Promise((res,rej)=>{ try{ const q=indexedDB.open(DBN,1); q.onupgradeneeded=()=>q.result.createObjectStore(ST); q.onsuccess=()=>{ const db=q.result; db.onclose=()=>{ _db=null; }; db.onversionchange=()=>{ db.close(); _db=null; }; res(db); }; q.onerror=()=>{ _db=null; rej(q.error); }; }catch(e){ _db=null; rej(e); } }); return _db; } // 연결 하나를 다시 씀
async function dbKeys(){ try{ const db=await idb(); return await new Promise(r=>{ const t=db.transaction(ST,"readonly").objectStore(ST).getAllKeys(); t.onsuccess=()=>r(t.result||[]); t.onerror=()=>r([]); }); }catch(e){ return []; } }
async function dbGet(k){ try{ const db=await idb(); return await new Promise(r=>{ const t=db.transaction(ST,"readonly").objectStore(ST).get(k); t.onsuccess=()=>r(t.result); t.onerror=()=>r(null); }); }catch(e){ return null; } }
async function dbPut(k,v){ try{ const db=await idb(); await new Promise((r,j)=>{ const t=db.transaction(ST,"readwrite"); t.objectStore(ST).put(v,k); t.oncomplete=r; t.onerror=()=>j(t.error); t.onabort=()=>j(t.error); }); return true; }catch(e){ return false; } }
async function dbDel(k){ try{ const db=await idb(); await new Promise(r=>{ const t=db.transaction(ST,"readwrite"); t.objectStore(ST).delete(k); t.oncomplete=r; t.onerror=r; }); }catch(e){} }
let _svT=null; const DIRTY=new Set(); let _quotaWarned=false;
function save(){ clearTimeout(_svT); _svT=setTimeout(flushSave,700); }
async function flushSave(){ clearTimeout(_svT); _svT=null; const snap=JSON.parse(JSON.stringify(S));
  await dbPut("state",snap);
  for(const id of [...DIRTY]){ DIRTY.delete(id); const P=PX.get(id); if(!P) { await dbDel("px:"+id); continue; } const im=byId(id); const ok=await dbPut("px:"+id,{v:im?im.ver:0,P}); // 판 번호와 함께 → 불러올 때 맞는지 확인
    if(!ok) DIRTY.add(id); // 다음 저장 때 다시
    if(!ok&&!_quotaWarned){ _quotaWarned=true; note("브라우저 저장 공간이 부족해 일부 영상 화소를 저장하지 못했습니다 — 새로 고치면 그 영상은 다시 열어야 합니다","warn"); } } }
function markDirty(id){ DIRTY.add(id); save(); }
try{ window.addEventListener("pagehide",()=>{ if(_svT||DIRTY.size) flushSave(); }); document.addEventListener("visibilitychange",()=>{ if(document.visibilityState==="hidden"&&(_svT||DIRTY.size)) flushSave(); }); }catch(e){}
async function loadSaved(){
  const st=await dbGet("state"); if(!st||st.v!==1) return;
  const d=defState(); S=Object.assign(d,st); S.ui=Object.assign(d.ui,st.ui||{}); S.ui.thr=Object.assign(d.ui.thr,(st.ui||{}).thr||{}); S.ui.pa=Object.assign(d.ui.pa,(st.ui||{}).pa||{}); S.ui.dist=Object.assign(d.ui.dist,(st.ui||{}).dist||{}); S.ui.scalebar=Object.assign(d.ui.scalebar,(st.ui||{}).scalebar||{});
  S.measure=Object.assign(d.measure,st.measure||{}); S.macro=Object.assign(d.macro,st.macro||{});
  const keep=[]; let stale=0; for(const im of S.images){ const R=await dbGet("px:"+im.id); const P=R&&!Array.isArray(R)&&R.P?R.P:R; if(R&&R.P&&R.v!==im.ver){ stale++; continue; } // 화소가 기록보다 오래됨(저장 실패) → 버림
    if(P&&P.length===im.nc*im.nz*im.nt){ PX.set(im.id,P); keep.push(im); } }
  if(stale) setTimeout(()=>note(`영상 ${stale}개는 마지막 작업의 화소가 저장되지 않아 열지 않았습니다 — 다시 여세요`,"warn"),60);
  { const ids=new Set(keep.map(x=>"px:"+x.id)); for(const k of await dbKeys()) if(typeof k==="string"&&k.startsWith("px:")&&!ids.has(k)) dbDel(k); } // 닫힌 영상의 화소 지움
  const lost=S.images.length-keep.length; S.images=keep; if(S.sel&&!byId(S.sel)) S.sel=keep[0]?keep[0].id:null;
  if(lost) setTimeout(()=>note(`저장되지 않은 영상 ${lost}개는 다시 열어야 합니다`,"warn"),50);
}

/* ---------------- 알림 줄 ---------------- */
let MSGS=[];
function note(m,k,act){ MSGS.unshift({m:String(m),k:k||"",act:act||null}); MSGS=MSGS.slice(0,8); drawMsg(); if(!ROOT||!$q(ROOT,"#imMsg")){ if(typeof toast==="function") toast(m,k==="warn"?"warn":k); } }
function drawMsg(){ const el=ROOT&&$q(ROOT,"#imMsg"); if(!el) return; if(!MSGS.length){ el.hidden=true; el.innerHTML=""; return; } el.hidden=false; const m=MSGS[0], old=MSGS.slice(1);
  el.innerHTML=`<div class="cur ${esc(m.k)}"><span class="ic">${m.k==="warn"?"!":m.k==="ok"?"✓":"i"}</span><span class="tx">${esc(m.m)}</span>${m.act?`<button class="btn xs" id="imMsgA">${esc(m.act.label)}</button>`:""}<button class="x" id="imMsgX" title="닫기">×</button></div>${old.length?`<details class="old"><summary>이전 알림 ${old.length}</summary>${old.map(o=>`<div class="${esc(o.k)}">${esc(o.m)}</div>`).join("")}</details>`:""}`;
  $q(el,"#imMsgX").onclick=()=>{ MSGS.shift(); drawMsg(); }; const b=$q(el,"#imMsgA"); if(b) b.onclick=()=>{ const f=m.act.fn; m.act=null; drawMsg(); f(); }; }
const _cf={}; function confirm2(k,msg){ const n=Date.now(); if(_cf[k]&&n-_cf[k]<4000){ delete _cf[k]; return true; } _cf[k]=n; note(msg,"warn"); return false; }

/* ---------------- 진행 표시 ---------------- */
let BUSY=null;
function setProg(t,f){ const el=ROOT&&$q(ROOT,"#imProg"); if(!el) return; el.hidden=false; $q(el,".t").textContent=t; $q(el,".bar i").style.width=Math.round(100*clamp(f||0,0,1))+"%"; }
function hideProg(){ const el=ROOT&&$q(ROOT,"#imProg"); if(el) el.hidden=true; }

/* ---------------- Worker RPC ----------------
   IMG_CORE 함수 이름과 인자(plane 들)를 보내면 결과를 받는다. 여러 장(스택)은 일꾼 여러 개가 나눠 맡는다 */
let WURL=null; const WK=[]; let JOBSEQ=0;
function workerURL(){ if(WURL) return WURL; const el=document.getElementById("img-core-src"); if(!el) return null;
  const body=el.textContent+`
self.onmessage=function(e){ const d=e.data; const IC=self.IMG_CORE; IC._prog=f=>self.postMessage({id:d.id,prog:f});
  try{ const out=IC[d.fn].apply(null,d.args); const tr=[]; const walk=(o,dep)=>{ if(!o||dep>4) return; if(ArrayBuffer.isView(o)){ if(o.buffer&&!tr.includes(o.buffer)) tr.push(o.buffer); return; } if(typeof o==="object") for(const k in o) walk(o[k],dep+1); }; walk(out,0);
    self.postMessage({id:d.id,ok:true,out},tr); }catch(err){ self.postMessage({id:d.id,ok:false,err:String(err&&err.message||err)}); } };`;
  WURL=URL.createObjectURL(new Blob([body],{type:"text/javascript"})); return WURL; }
const nW=()=>Math.max(1,Math.min(6,((navigator&&navigator.hardwareConcurrency)||2)-1));
function getWorker(){ const url=workerURL(); if(!url) return null; let w=WK.find(x=>!x.busy); if(!w&&WK.length<nW()){ w={w:new Worker(url),busy:false}; WK.push(w); } return w||null; }
function killWorkers(){ WK.forEach(x=>{ try{ x.w.terminate(); }catch(e){} }); WK.length=0; }
/* call(fn, args, onProg) → Promise. 일꾼이 없으면(시험 환경) 주 스레드 */
function call(fn,args,onProg){
  const gen=CANCEL_GEN;
  return new Promise((res,rej)=>{
    const tryRun=()=>{ if(gen!==CANCEL_GEN){ rej(new Error("취소")); return; } const w=getWorker();
      if(!workerURL()){ try{ IC._prog=onProg||null; const o=IC[fn].apply(null,args); IC._prog=null; res(o); }catch(e){ rej(e); } return; }
      if(!w){ setTimeout(tryRun,15); return; }
      w.busy=true; const id=++JOBSEQ; const cancel=()=>{ try{ w.w.terminate(); }catch(e){} const i=WK.indexOf(w); if(i>=0) WK.splice(i,1); rej(new Error("취소")); };
      CANCELS.add(cancel);
      w.w.onmessage=e=>{ const d=e.data; if(d.id!==id) return; if(d.prog!=null){ onProg&&onProg(d.prog); return; } CANCELS.delete(cancel); w.w.onmessage=w.w.onerror=null; w.busy=false; if(d.ok) res(d.out); else rej(new Error(d.err)); }; // 끝난 일의 처리기를 떼어 큰 배열을 놓아 줌
      w.w.onerror=e=>{ e.preventDefault&&e.preventDefault(); CANCELS.delete(cancel); w.w.onmessage=w.w.onerror=null; w.busy=false; rej(new Error(e.message||"Worker 오류")); };
      w.w.postMessage({id,fn,args}); };
    tryRun();
  });
}
const CANCELS=new Set();
let CANCEL_GEN=0; // 중단할 때마다 늘림: 기다리던 일·batch 고리도 이 값을 보고 멈춤
function cancelAll(){ CANCEL_GEN++; for(const c of [...CANCELS]) c(); CANCELS.clear(); hideProg(); } // BUSY 는 그 작업이 끝나며 직접 풂

/* ---------------- 영상 추가·삭제·되돌리기 ---------------- */
function addImage(d,o){
  o=o||{}; const id=uid(); const nc=d.nc||1;
  const im={id,name:d.name,w:d.w,h:d.h,bit:d.bit,rgb:!!d.rgb,nc:d.rgb?3:nc,nz:d.nz||1,nt:d.nt||1,cal:d.cal||null,calCandidates:d.calCandidates||[],info:d.info||[],warn:d.warn||[],
    desc:d.desc||"",metaText:d.metaText||"",databar:d.databar||null,invertLut:!!d.invertLut,c:0,z:0,t:0,composite:(nc>1&&!d.rgb),parent:o.parent||null,ops:o.ops||[],created:new Date().toISOString(),
    file:o.file||d.file||d.name,src:o.src||"file",valueUnit:d.valueUnit||"",siScale:d.siScale>0?d.siScale:null,siUnit:d.siUnit||"",ver:1,kind:o.kind||"",fft:o.fft||null,
    chNames:Array.isArray(d.chNames)&&d.chNames.length===nc&&!d.rgb?d.chNames.slice():null,zStep:d.zStep>0?d.zStep:0};
  im.disp=[...Array(im.rgb?1:im.nc)].map((_,c)=>({min:0,max:255,lut:(!im.rgb&&d.chLut&&d.chLut[c]&&LUTS.includes(d.chLut[c]))?d.chLut[c]:defaultLut(im,c)})); // 파일의 채널 색(nd2·czi) 우선
  PX.set(id,d.planes);
  if(d.display&&d.display[0]&&!im.rgb) { im.disp[0].min=d.display[0].min; im.disp[0].max=d.display[0].max; } else resetDisplay(im); // ImageJ 와 같이 열 때는 min–max (Auto 는 사용자가)
  if(o.after){ const i=S.images.findIndex(x=>x.id===o.after); S.images.splice(i+1,0,im); } else S.images.push(im);
  S.sel=id; markDirty(id); return im;
}
function defaultLut(im,c){ if(im.rgb) return "rgb"; if(im.nc>1){ return ["Green","Magenta","Cyan","Red","Yellow","Blue"][c%6]; } return im.invertLut?"Inverted":"Grays"; }
function closeImage(id,force){ const i=S.images.findIndex(x=>x.id===id); if(i<0) return; if(BUSY&&!force){ note("계산 중에는 영상을 닫을 수 없습니다 — 끝나거나 중단한 뒤에","warn"); return; } S.images.splice(i,1); PX.delete(id); UNDO.delete(id); ctDrop(id); delete S.rois[id]; delete S.roiSel[id]; delete S.particles[id]; delete S.ui.views[id];
  DISP.key=null; DISP.cv=null; DISPZ.key=null; DISPZ.cv=null; OVL.key=null; OVL.thr=null; THUMB.delete(id); ADJ.key=""; ADJ.v=null; // 표시 버퍼도 놓아 줌
  DIRTY.add(id); if(S.sel===id) S.sel=S.images[Math.min(i,S.images.length-1)]?S.images[Math.min(i,S.images.length-1)].id:null; save(); }
const UNDO_MAX=8, UNDO_BYTES=1.2e9; // 영상마다 8 단계, 모든 영상 합쳐 1.2 GB
function pushUndo(im,label){ const L=UNDO.get(im.id)||[]; const e={label,planes:planesOf(im),macroAt:null,calChanged:false,
    meta:{w:im.w,h:im.h,bit:im.bit,rgb:im.rgb,nc:im.nc,nz:im.nz,nt:im.nt,cal:im.cal,chNames:im.chNames||null,zStep:im.zStep||0,composite:im.composite,databar:im.databar||null,valueUnit:im.valueUnit||"",siScale:im.siScale||null,siUnit:im.siUnit||"",disp:JSON.parse(JSON.stringify(im.disp)),c:im.c,z:im.z,t:im.t,ops:im.ops.slice(),kind:im.kind}};
  L.push(e); while(L.length>UNDO_MAX) L.shift(); UNDO.set(im.id,L);
  // 전체 메모리 상한: 가장 오래된 항목부터 버림 (지금 넣은 것은 남김)
  // 전체 크기: 같은 배열(한 장만 바뀐 스택의 나머지 장 등)과 지금 보이는 화소는 한 번만·세지 않음
  const live=new Set(); for(const P of PX.values()) for(const p of P) live.add(p.buffer); const extra=x=>{ let b=0; for(const p of x.planes){ if(!live.has(p.buffer)&&!SEEN.has(p.buffer)){ SEEN.add(p.buffer); b+=p.byteLength; } } return b; };
  let SEEN=new Set(); let tot=0; for(const [,LL] of UNDO) for(const x of LL) tot+=extra(x);
  if(tot>UNDO_BYTES){ for(const [id,LL] of UNDO){ while(LL.length&&tot>UNDO_BYTES&&LL[0]!==e){ const o=LL.shift(); SEEN=new Set(); tot=0; for(const [,L2] of UNDO) for(const x of L2) tot+=extra(x); void o; } void id; } }
  return e; }
function undo(){ const im=cur(); if(!im) return; if(BUSY){ note("계산 중에는 되돌릴 수 없습니다 — 끝나거나 중단한 뒤에","warn"); return; } const L=UNDO.get(im.id); if(!L||!L.length){ note("되돌릴 작업이 없습니다"); return; } const u=L.pop();
  const keepCal=im.cal; PX.set(im.id,u.planes); Object.assign(im,u.meta); if(!u.calChanged) im.cal=keepCal; im.ver++; markDirty(im.id);
  // 이 되돌리기에 해당하는 매크로 단계만 지운다
  if(u.macroAt!=null&&S.macro.steps[u.macroAt]&&S.macro.steps[u.macroAt]._u===u.uid) S.macro.steps.splice(u.macroAt,1);
  invalidate(im); save(); rerender(); note(`되돌림: ${u.label}`); }
/* 화소 교체 (제자리 연산) */
function setPlanes(im,planes,meta,label,step){
  if(!byId(im.id)) throw new Error("영상이 닫혔습니다"); // 계산 중 닫힌 영상에 결과를 쓰지 않음
  const u=pushUndo(im,label); u.uid=uid(); const cal0=im.cal, nc0=im.nc; PX.set(im.id,planes); if(meta) Object.assign(im,meta); if(im.nc!==nc0||im.rgb) im.chNames=null; /* 채널 수가 바뀌면 이름은 더 맞지 않음 */ im.ver++; u.calChanged=im.cal!==cal0;
  if(step){ im.ops.push(step); const n=S.macro.steps.length; record(Object.assign({},step,{_u:u.uid}),im); if(S.macro.steps.length>n) u.macroAt=n; }
  if(meta&&(meta.bit!=null||meta.nc!=null)){ im.disp=[...Array(im.rgb?1:im.nc)].map((_,c)=>({min:0,max:255,lut:(im.disp[c]&&im.disp[c].lut)||defaultLut(im,c)})); }
  if(!meta||!meta.keepDisplay) resetDisplay(im);
  markDirty(im.id); invalidate(im);
}

/* ---------------- 표시 범위 ---------------- */
function autoDisplay(im,quiet){ // ImageJ 'Auto' 처럼 포화 0.35 % (8-bit 은 0..255 아닌 실제 분포)
  if(im.rgb){ im.disp[0].min=0; im.disp[0].max=255; return; }
  for(let c=0;c<im.nc;c++){ const p=planeAt(im,c); let r; try{ r=IC.saturatedRange(p,0.35); }catch(e){ r=IC.minMax(p.data); }
    if(!(r[1]>r[0])){ const m=IC.minMax(p.data); r=[m[0],m[1]>m[0]?m[1]:m[0]+1]; }
    if(IC.isBinary(p)) r=[0,255];
    im.disp[c].min=r[0]; im.disp[c].max=r[1]; }
  void quiet;
}
function resetDisplay(im){ for(let c=0;c<(im.rgb?1:im.nc);c++){ if(im.rgb||im.bit===8){ im.disp[c].min=0; im.disp[c].max=255; continue; } const m=IC.minMax(planeAt(im,c).data); im.disp[c].min=m[0]; im.disp[c].max=m[1]>m[0]?m[1]:m[0]+1; } }
