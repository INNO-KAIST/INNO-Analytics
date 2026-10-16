/* ============================================================
   INNO Analytics — Diffraction Phase DB (DOM 없음)
   층(layer) 셋:  pack  = 앱에 들어 있는 기본 상 (DIFF_PACK, 읽기 전용)
                 lab   = 연구실 공유 Firestore (projectId inno-analytics-base) — 내려받아 IndexedDB 에 캐시
                 local = 이 브라우저에만 (IndexedDB "inno-diffraction")
   · Firebase SDK 는 로그인·동기화할 때만 import() 로 불러온다 (없어도 pack·local·캐시로 동작).
   · Firestore 규칙(게시본)을 그대로 따른다: 새 항목 version 1·hidden false·createdBy = uid,
     고치기 = 이전 내용을 history/{옛 version} 에 만들고 version+1, 지우기 = hidden true, verified 는 관리자만.
   · 판정 재현: refOf(id) → {id, layer, version, updatedAt, hash} 를 결과에 저장, getVersion(ref) 로 그 내용을 다시 찾는다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const DC=G.DIFF_CORE||(typeof require!=="undefined"?require("./diff_core.js"):null);
const DB={};

/* ---------------- 설정 ---------------- */
const FB_VER="10.12.2", FB_BASE="https://www.gstatic.com/firebasejs/"+FB_VER+"/";
const EXPECT_PROJECT="inno-analytics-base";
const DEFAULT_CONFIG=Object.freeze({apiKey:"AIzaSyBIpiaPfrXwG6OQe7qQ0qJYIfcXBUZRCzg",authDomain:"inno-analytics-base.firebaseapp.com",projectId:"inno-analytics-base",storageBucket:"inno-analytics-base.firebasestorage.app",messagingSenderId:"394709749281",appId:"1:394709749281:web:3d0e3ac1cbced468e49d60"});
const IDB_NAME="inno-diffraction", IDB_VER=1;
const CONTENT=["name","formula","cell","sg","symops","atoms","source","note"];
const VER_META=["verified","verifiedBy","verifiedByEmail","verifiedAt","verifiedVersion","hidden","promotedFrom"];
const LOCAL_HIST_MAX=30;
const POPUP_FALLBACK=new Set(["auth/popup-blocked","auth/operation-not-supported-in-this-environment","auth/internal-error","auth/web-storage-unsupported"]);
const FILE_MSG="Google 로그인은 웹 주소(https://inno-kaist.github.io/INNO-Analytics/)로 열었을 때만 됩니다. 파일(file://)로 연 앱에서는 Built-in·Local 과 마지막으로 내려받은 Lab 캐시만 씁니다.";

/* ---------------- 상태 ---------------- */
const MEM={local:new Map(),lab:new Map(),cal:new Map(),meta:{}};
const PACK={src:null,version:null,recs:new Map()};
const PH=new Map();                                 // id → {sig, phase}
const FBS={sdk:null,app:null,auth:null,fs:null,user:null,isMember:null,isAdmin:false,admins:null,error:null,connecting:null,unsub:null};
let IDB_OVERRIDE=null, SDK_OVERRIDE=null, ENV={}, dbp=null, INIT=null, STORAGE="memory", LAB_IGNORED=null;
const listeners=new Set();

/* ---------------- 작은 도구 ---------------- */
const clone=o=>o==null?o:JSON.parse(JSON.stringify(o));
const now=()=>Date.now();
function rid(){ let s=""; try{ const a=new Uint8Array(8); (G.crypto||require("crypto").webcrypto).getRandomValues(a); s=[...a].map(x=>(x%36).toString(36)).join(""); }catch(e){ s=Math.random().toString(36).slice(2,10); } return now().toString(36)+s; }
function err(code,message,extra){ const e=new Error(message); e.code=code; e._own=true; if(extra) Object.assign(e,extra); return e; }
const tsMs=v=>{ if(v==null) return null; if(typeof v==="number") return v; if(typeof v.toMillis==="function") return v.toMillis(); if(v instanceof Date) return v.getTime(); if(typeof v.seconds==="number") return v.seconds*1000+(v.nanoseconds||0)/1e6; return null; };
function emit(type,extra){ const ev=Object.assign({type},extra||{}); for(const cb of [...listeners]){ try{ cb(ev); }catch(e){ if(G.console) console.error(e); } } }
DB.onChange=function(cb){ listeners.add(cb); return ()=>listeners.delete(cb); };
/* FNV-1a 32 bit — 내용(CONTENT 필드) 지문 */
function contentHash(r){ const s=JSON.stringify(canon(contentOf(r))); let h=0x811c9dc5; for(let i=0;i<s.length;i++){ h^=s.charCodeAt(i); h=Math.imul(h,0x01000193)>>>0; } return ("0000000"+h.toString(16)).slice(-8); }
function canon(o){ if(Array.isArray(o)) return o.map(canon); if(o&&typeof o==="object"){ const r={}; for(const k of Object.keys(o).sort()) r[k]=canon(o[k]); return r; } return o; }
const num=v=>{ if(v===""||v==null) return null; const x=+v; return isFinite(x)?x:null; };

/* PhaseRec 의 내용만 (정규화, undefined·NaN 없음 — Firestore 가 undefined 를 받지 않는다) */
function contentOf(r){ r=r||{}; const o={};
  o.name=String(r.name||"").trim(); o.formula=String(r.formula||"").trim();
  const c=r.cell||{}; const a=num(c.a); o.cell={a, b:num(c.b)??a, c:num(c.c)??a, al:num(c.al)??90, be:num(c.be)??90, ga:num(c.ga)??90};
  const sg=r.sg||{}; o.sg={number:num(sg.number)!=null?Math.round(+sg.number):null, symbol:String(sg.symbol||"").trim()}; if(sg.hall) o.sg.hall=String(sg.hall);
  if(Array.isArray(r.symops)&&r.symops.length) o.symops=r.symops.map(String);
  if(Array.isArray(r.atoms)&&r.atoms.length) o.atoms=r.atoms.map(at=>{ const q={label:String(at.label||at.el||""),el:String(at.el||""),x:num(at.x),y:num(at.y),z:num(at.z)}; q.occ=num(at.occ)??1; q.B=num(at.B)??0; return q; });
  const s=r.source||{}; o.source={kind:["COD","DOI","paper","measured","manual"].includes(s.kind)?s.kind:(s.codId?"COD":s.doi?"DOI":"manual")};
  for(const k of ["codId","doi","citation","note"]) if(s[k]!=null&&String(s[k]).trim()) o.source[k]=String(s[k]).trim();
  if(r.note!=null&&String(r.note).trim()) o.note=String(r.note).trim();
  return clone(o); }
DB.contentOf=contentOf;
DB.contentHash=contentHash;

/* ---------------- IndexedDB (없으면 메모리만) ---------------- */
function idbFactory(){ return IDB_OVERRIDE||(typeof indexedDB!=="undefined"?indexedDB:null); }
function openDB(){ if(dbp) return dbp; const f=idbFactory(); if(!f) return (dbp=Promise.resolve(null));
  dbp=new Promise(res=>{ try{ const q=f.open(IDB_NAME,IDB_VER);
    q.onupgradeneeded=()=>{ const d=q.result; for(const s of ["local","lab","cal"]) if(!d.objectStoreNames.contains(s)) d.createObjectStore(s,{keyPath:"id"}); if(!d.objectStoreNames.contains("meta")) d.createObjectStore("meta"); };
    q.onsuccess=()=>res(q.result); q.onerror=()=>res(null); }catch(e){ res(null); } });
  return dbp; }
async function tx(store,mode,fn){ const d=await openDB(); if(!d) return null; return new Promise((res,rej)=>{ try{ const t=d.transaction(store,mode); const os=t.objectStore(store); const out=fn(os); t.oncomplete=()=>res(out&&out.result!==undefined?out.result:out); t.onerror=()=>rej(t.error); t.onabort=()=>rej(t.error||err("idb-abort","IndexedDB 저장이 취소되었습니다")); }catch(e){ rej(e); } }); }
const idbPut=(s,v,k)=>tx(s,"readwrite",os=>{ k===undefined?os.put(v):os.put(v,k); });
const idbDel=(s,k)=>tx(s,"readwrite",os=>{ os.delete(k); });
const idbReplace=(s,list,keep)=>tx(s,"readwrite",os=>{ if(keep){ const q=os.openCursor(); q.onsuccess=()=>{ const c=q.result; if(!c) { for(const v of list) os.put(v); return; } if(!keep(c.value)) c.delete(); c.continue(); }; } else { os.clear(); for(const v of list) os.put(v); } });
async function idbAll(s){ const d=await openDB(); if(!d) return []; return new Promise(res=>{ try{ const q=d.transaction(s,"readonly").objectStore(s).getAll(); q.onsuccess=()=>res(q.result||[]); q.onerror=()=>res([]); }catch(e){ res([]); } }); }
async function idbMeta(){ const d=await openDB(); if(!d) return {}; return new Promise(res=>{ try{ const os=d.transaction("meta","readonly").objectStore("meta"); const out={}; const q=os.openCursor(); q.onsuccess=()=>{ const c=q.result; if(!c) return res(out); out[c.key]=c.value; c.continue(); }; q.onerror=()=>res(out); }catch(e){ res({}); } }); }
async function setMeta(k,v){ MEM.meta[k]=v; try{ await idbPut("meta",v,k); }catch(e){} }

/* ---------------- 기본 상 (DIFF_PACK) ---------------- */
function loadPack(){ let P=G.DIFF_PACK;
  if(!P&&typeof module!=="undefined"&&module.exports&&typeof require!=="undefined"&&!PACK.tried){ PACK.tried=true; try{ const m=require("./diff_pack.js"); P=G.DIFF_PACK||(m&&m.phases?m:null); }catch(e){} }
  if(!P||!Array.isArray(P.phases)||PACK.src===P) return; PACK.src=P; PACK.version=P.version||null; PACK.recs=new Map();
  for(const r of P.phases){ if(!r||!r.id) continue; PACK.recs.set(r.id,Object.assign(clone(r),{hidden:false})); }
  emit("pack",{version:PACK.version,count:PACK.recs.size}); }
DB.packInfo=()=>{ loadPack(); return {version:PACK.version,count:PACK.recs.size}; };

/* ---------------- 시작 ---------------- */
DB.init=function(){ if(INIT) return INIT; INIT=(async()=>{ loadPack(); const d=await openDB(); STORAGE=d?"indexeddb":"memory";
  if(d){ const [loc,lab,cal,meta]=await Promise.all([idbAll("local"),idbAll("lab"),idbAll("cal"),idbMeta()]);
    MEM.meta=meta||{}; for(const r of loc) MEM.local.set(r.id,r);
    const pid=currentConfig().projectId;
    if(MEM.meta.cacheProject&&MEM.meta.cacheProject!==pid){ LAB_IGNORED={cacheProject:MEM.meta.cacheProject,count:lab.length}; }   // 다른 프로젝트 캐시는 섞지 않는다
    else for(const r of lab) MEM.lab.set(r.id,r);
    for(const c of cal) if(c.layer!=="lab"||!LAB_IGNORED) MEM.cal.set(c.id,c); }
  emit("init",{storage:STORAGE}); return DB; })(); return INIT; };
DB.ready=()=>DB.init();
DB.storage=()=>STORAGE;

/* ---------------- 읽기 ---------------- */
const layerMap=L=>L==="pack"?(loadPack(),PACK.recs):L==="lab"?MEM.lab:L==="local"?MEM.local:null;
const tag=(r,L)=>Object.assign(clone(r),{layer:L});
function findRaw(id){ loadPack(); for(const L of ["pack","lab","local"]){ const r=layerMap(L).get(id); if(r) return [r,L]; } return [null,null]; }
DB.get=function(id){ const [r,L]=findRaw(id); return r?tag(r,L):null; };
DB.layerOf=id=>findRaw(id)[1];
/* 원소: 화학식에서, 없으면 원자 자리에서 */
function elementsOf(r){ const f=DC.parseFormula(r.formula||""); let el=Object.keys(f).filter(e=>DC.elementOf(e)===e); if(!el.length&&r.atoms) el=[...new Set(r.atoms.map(a=>DC.elementOf(a.el||a.label)).filter(Boolean))]; return el; }
DB.elementsOf=elementsOf;
function parseEls(e){ if(!e) return null; const arr=Array.isArray(e)?e:String(e).split(/[\s,;]+/); const out=[]; for(let s of arr){ s=String(s).trim(); if(!s) continue; s=s[0].toUpperCase()+s.slice(1).toLowerCase(); if(DC.elementOf(s)===s&&!out.includes(s)) out.push(s); } return out; }
DB.parseElements=parseEls;
function matchEls(r,o){ const want=parseEls(o.elements); if(!want||!want.length) return true; const have=elementsOf(r); if(!have.length) return false;
  if(o.elementMode==="all"||o.elementMode==="contains") return want.every(e=>have.includes(e));
  return have.every(e=>want.includes(e)); }                   // 기본: 상의 원소 ⊆ 고른 원소 (REDI 의 EDS 원소)
/* 화학식 비교: 원소별 비(가장 작은 수로 나눔)가 같으면 같은 식 — 'TiO2' = 'O2 Ti', 'NaYF4' = 'Na1.5 Y1.5 F6' */
function ratioKey(f){ const k=Object.keys(f).filter(e=>DC.elementOf(e)===e).sort(); if(!k.length) return ""; const m=Math.min(...k.map(e=>f[e])); return k.map(e=>e+(+(f[e]/m).toFixed(3))).join(""); }
function matchQ(r,q){ if(!q) return true; const raw=String(q).trim(); q=raw.toLowerCase(); if(!q) return true; const qs=q.replace(/\s+/g,"");
  const s=r.source||{}; const hay=[r.name,r.id,s.codId,s.doi,r.sg&&r.sg.symbol].map(x=>String(x||"").toLowerCase());
  if(hay.some(h=>h.includes(q))||String(r.formula||"").replace(/\s+/g,"").toLowerCase().includes(qs)||(r.sg&&String(r.sg.number)===q)) return true;
  if(/^(\s*[A-Z][a-z]?\d*\.?\d*)+\s*$/.test(raw)){ const a=ratioKey(DC.parseFormula(raw)); return !!a&&a===ratioKey(DC.parseFormula(r.formula||"")); }
  return false; }
DB.list=function(o){ o=o||{}; const layers=o.layers||["pack","lab","local"]; const out=[];
  for(const L of layers){ const m=layerMap(L); if(!m) continue; for(const r of m.values()){ if(r.hidden&&!o.includeHidden) continue; if(!matchEls(r,o)||!matchQ(r,o.query)) continue; out.push(tag(r,L)); } }
  const ord={pack:0,lab:1,local:2}; out.sort((a,b)=>String(a.name||a.formula).localeCompare(String(b.name||b.formula),"en",{sensitivity:"base"})||ord[a.layer]-ord[b.layer]);
  return out; };
DB.counts=()=>({pack:(loadPack(),PACK.recs.size),lab:MEM.lab.size,labHidden:[...MEM.lab.values()].filter(r=>r.hidden).length,local:MEM.local.size,cal:MEM.cal.size});
/* 계산용 상 (DC.phaseFromRecord, 버전이 바뀔 때만 다시) */
DB.phase=function(id){ const [r,L]=findRaw(id); if(!r) return null; const sig=L+"|"+(r.version||0)+"|"+(r.updatedAt||0)+"|"+(L==="pack"?PACK.version:""); const c=PH.get(id); if(c&&c.sig===sig) return c.phase;
  const p=DC.phaseFromRecord(r); PH.set(id,{sig,phase:p}); return p; };
/* 판정에 쓴 상의 기록: 결과에 함께 저장한다 */
DB.refOf=function(id){ const [r,L]=findRaw(id); if(!r) return null; const ref={id,layer:L,version:r.version!=null?r.version:null,updatedAt:r.updatedAt!=null?r.updatedAt:null,hash:contentHash(r),name:r.name||r.formula||id};
  if(L==="pack") ref.packVersion=PACK.version; if(L==="lab") ref.projectId=MEM.meta.cacheProject||currentConfig().projectId; return ref; };
/* refOf 로 저장한 판본을 다시 찾기: {rec, exact} — exact=false 면 지금 내용과 다를 수 있음 */
DB.getVersion=async function(ref){ if(!ref||!ref.id) return null; await DB.init(); const [r,L]=findRaw(ref.id);
  const ok=x=>x&&(!ref.hash||contentHash(x)===ref.hash);
  if(r&&(ref.version==null||r.version===ref.version)&&ok(r)) return {rec:tag(r,L),exact:true};
  if(L==="local"&&r&&Array.isArray(r.history)){ const h=r.history.find(x=>x.version===ref.version); if(h) return {rec:Object.assign(clone(h.content),{id:ref.id,version:h.version,updatedAt:h.updatedAt,layer:"local"}),exact:ok(h.content)}; }
  if((L==="lab"||ref.layer==="lab")&&ref.version!=null&&FBS.isMember){ try{ const s=await getDocS(docRef("phases",ref.id,"history",String(ref.version))); if(s.exists()){ const d=fromFS(ref.id,s.data()); return {rec:tag(d,"lab"),exact:ok(d)}; } }catch(e){} }
  if(r) return {rec:tag(r,L),exact:false,note:"저장된 판본을 찾지 못해 현재 내용을 돌려줍니다"};
  return null; };

/* 중복: 모든 층(숨긴 것 제외)에서 DC.sameRecord */
DB.findDuplicate=function(rec,o){ o=o||{}; const ex=new Set(o.exclude||[]); if(rec&&rec.id) ex.add(rec.id); const c=contentOf(rec); if(!(c.cell.a>0)) return null;
  for(const r of DB.list({layers:o.layers})){ if(ex.has(r.id)) continue; try{ if(DC.sameRecord(c,r)) return r; }catch(e){} } return null; };
DB.validate=rec=>DC.validateRecord(contentOf(rec));

/* ---------------- CIF ---------------- */
DB.importCIF=function(text,o){ o=o||{}; const out=[]; let blocks;
  try{ blocks=DC.cifBlocks(text).filter(b=>b.items["_cell_length_a"]!=null); }catch(e){ return [{error:"CIF 를 읽지 못했습니다: "+(e.message||e),file:o.fileName||""}]; }
  if(!blocks.length) return [{error:"CIF 에 결정 구조 블록(_cell_length_a)이 없습니다",file:o.fileName||""}];
  for(const b of blocks){ const it={block:b.name,file:o.fileName||""};
    try{ const ph=DC.phaseFromCIF(b); const rec=DC.recordFromPhase(ph); const I=b.items;
      const src={codId:ph.source.codId||"",doi:ph.source.doi||""};
      const j=[I["_journal_name_full"]||I["_journal_name_abbrev"],I["_journal_year"],I["_journal_volume"],I["_journal_page_first"]].map(x=>String(x||"").replace(/['"]/g,"").trim()).filter(x=>x&&x!=="?"&&x!==".");
      const au=(b.loops.find(l=>l.keys.includes("_publ_author_name"))||{rows:[]}).rows.map(r=>r[0]);
      if(j.length) src.citation=(au.length?au[0].replace(/['"]/g,"")+(au.length>1?" et al., ":", "):"")+j.join(", ");
      src.kind=src.codId?"COD":src.doi?"DOI":src.citation?"paper":"manual";
      src.note="CIF"+(o.fileName?" "+o.fileName:"")+" · data_"+b.name;
      rec.source=src; if(!rec.sg.number&&rec.sg.symbol){ const e=DC.sgBySymbol(rec.sg.symbol); if(e) rec.sg.number=e.number; }
      it.rec=contentOf(rec); it.validation=DC.validateRecord(it.rec);
      it.duplicate=DB.findDuplicate(it.rec)||null;
      if(!it.duplicate){ const prev=out.find(p=>p.rec&&DC.sameRecord(p.rec,it.rec)); if(prev) it.duplicate=Object.assign({layer:"import",id:"",name:prev.rec.name},prev.rec); }
    }catch(e){ it.error=String(e.message||e); }
    out.push(it); }
  return out; };

/* ---------------- Local ---------------- */
/* 저장: {ok, rec, validation} 또는 {ok:false, reason:"invalid"|"duplicate", validation, duplicate} */
DB.saveLocal=async function(rec,o){ o=o||{}; await DB.init(); const content=contentOf(rec); const v=DC.validateRecord(content); if(!v.ok) return {ok:false,reason:"invalid",validation:v};
  const old=rec&&rec.id?MEM.local.get(rec.id):null; if(rec&&rec.id&&!old&&!/^local:/.test(rec.id)) throw err("bad-id","Local 항목의 id 는 'local:' 로 시작해야 합니다");
  if(!o.force&&(!old||o.checkDuplicate)){ const dup=DB.findDuplicate(content,{exclude:[rec&&rec.id].concat(o.exclude||[])}); if(dup) return {ok:false,reason:"duplicate",duplicate:dup,validation:v}; }
  const t=now(); const r=Object.assign(content,{id:old?old.id:(rec.id||"local:"+rid()),version:old?(old.version|0)+1:1,createdAt:old?old.createdAt:t,updatedAt:t});
  for(const k of ["copiedFrom","promotedTo","importedFrom"]) if(rec[k]||(old&&old[k])) r[k]=rec[k]||old[k];
  if(old){ if(contentHash(old)===contentHash(r)&&!o.forceVersion) return {ok:true,rec:tag(old,"local"),validation:v,unchanged:true};
    r.history=[{version:old.version|0,updatedAt:old.updatedAt,content:contentOf(old)}].concat(old.history||[]).slice(0,LOCAL_HIST_MAX); }
  MEM.local.set(r.id,r); PH.delete(r.id); await idbPut("local",r); emit("local",{ids:[r.id]}); return {ok:true,rec:tag(r,"local"),validation:v}; };
DB.removeLocal=async function(id){ await DB.init(); if(!MEM.local.has(id)) return false; MEM.local.delete(id); PH.delete(id); await idbDel("local",id); emit("local",{ids:[id],removed:true}); return true; };
/* Built-in·Lab → Local 사본 (고쳐 쓰려고) */
DB.copyToLocal=async function(id){ const r=DB.get(id); if(!r) throw err("not-found","상을 찾지 못했습니다: "+id); return DB.saveLocal(Object.assign(contentOf(r),{copiedFrom:id}),{force:true}); };

/* ---------------- 내보내기·가져오기 ---------------- */
const pub=r=>{ const o=clone(r); if(o.layer!=="local") delete o.history; return o; };   // Local 은 지난 판본까지 (재현용)
DB.exportJSON=function(ids,o){ o=o||{}; const recs=ids&&ids.length?ids.map(DB.get).filter(Boolean):DB.list({includeHidden:true,layers:o.layers||["lab","local"]});
  const doc={format:"inno-diffraction-phasedb",v:1,exportedAt:new Date().toISOString(),projectId:currentConfig().projectId,packVersion:PACK.version,phases:recs.map(pub)};
  if(!(ids&&ids.length)&&o.calibrations!==false) doc.calibrations=DB.listCalibrations();
  return JSON.stringify(doc,null,1); };
const csvCell=v=>{ const s=v==null?"":String(v); return /[",\n\r]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s; };
DB.exportCSV=function(ids,o){ o=o||{}; const recs=ids&&ids.length?ids.map(DB.get).filter(Boolean):DB.list({includeHidden:true,layers:o.layers||["pack","lab","local"]});
  const H=["id","layer","name","formula","sg_number","sg_symbol","system","a_nm","b_nm","c_nm","alpha_deg","beta_deg","gamma_deg","n_atoms","source_kind","cod_id","doi","citation","verified","hidden","version","updated"];
  const rows=recs.map(r=>{ const c=r.cell||{}, s=r.source||{}, n=r.sg&&r.sg.number; return [r.id,r.layer,r.name,r.formula,n||"",r.sg&&r.sg.symbol||"",n?DC.crystalSystem(+n):"",c.a,c.b,c.c,c.al,c.be,c.ga,(r.atoms||[]).length,s.kind||"",s.codId||"",s.doi||"",s.citation||"",r.verified?1:0,r.hidden?1:0,r.version==null?"":r.version,r.updatedAt?new Date(r.updatedAt).toISOString():""]; });
  return [H].concat(rows).map(r=>r.map(csvCell).join(",")).join("\r\n")+"\r\n"; };
/* 가져오기 → Local 층. local: id 는 그대로(덮어쓰기), 다른 층 항목은 이미 있으면 건너뛰고 없으면 Local 사본 */
DB.importJSON=async function(text){ await DB.init(); let d; try{ d=typeof text==="string"?JSON.parse(text):text; }catch(e){ throw err("bad-json","JSON 을 읽지 못했습니다: "+e.message); }
  const phases=Array.isArray(d)?d:(d&&Array.isArray(d.phases)?d.phases:null); if(!phases) throw err("bad-format","Phase DB 백업 형식이 아닙니다 (phases 배열이 없음)");
  const rep={added:0,updated:0,unchanged:0,skipped:[],errors:[],calibrations:0};
  for(const r of phases){ try{ const id=String(r&&r.id||""); const isLocal=/^local:/.test(id)||r.layer==="local";
      if(isLocal&&id){ const had=MEM.local.has(id); const res=await DB.saveLocal(Object.assign({},r,{id}),{force:true}); if(!res.ok){ rep.errors.push({id,name:r.name,errors:res.validation.errors}); continue; }
        if(had&&res.unchanged) rep.unchanged++; else if(had) rep.updated++; else rep.added++;
        if(!had){ const raw=MEM.local.get(id); if(r.version) raw.version=r.version; if(r.createdAt) raw.createdAt=r.createdAt; if(r.updatedAt) raw.updatedAt=r.updatedAt; if(Array.isArray(r.history)) raw.history=r.history.slice(0,LOCAL_HIST_MAX); await idbPut("local",raw); }
        continue; }
      if(id&&findRaw(id)[0]){ rep.skipped.push({id,name:r.name,reason:"이미 있음"}); continue; }
      const dup=DB.findDuplicate(r); if(dup){ rep.skipped.push({id,name:r.name,reason:"중복: "+(dup.name||dup.id)+" ("+dup.layer+")"}); continue; }
      const res=await DB.saveLocal(Object.assign(contentOf(r),{importedFrom:id||undefined}),{force:true}); if(res.ok) rep.added++; else rep.errors.push({id,name:r.name,errors:res.validation.errors});
    }catch(e){ rep.errors.push({id:r&&r.id,name:r&&r.name,errors:[String(e.message||e)]}); } }
  if(d&&Array.isArray(d.calibrations)) for(const c of d.calibrations){ if(!c||!c.id) continue; if(MEM.cal.has(c.id)) continue; const lc=Object.assign(clone(c),{layer:"local"}); if(c.layer==="lab"){ lc.importedFrom=c.id; lc.id="local:"+rid(); } MEM.cal.set(lc.id,lc); await idbPut("cal",lc); rep.calibrations++; }
  emit("local",{import:rep}); return rep; };

/* ============================================================
   Lab (Firestore)
   ============================================================ */
function currentConfig(){ const o=MEM.meta.fbConfig; return o&&o.projectId?o:DEFAULT_CONFIG; }
DB.config=()=>Object.assign({},currentConfig());
/* 설정 바꾸기 (null = 기본). 다른 projectId 면 상태에 위험 경고, 쓰기는 막는다 */
DB.setConfig=async function(cfg){ await DB.init(); await setMeta("fbConfig",cfg&&cfg.projectId?Object.assign({},cfg):null); FBS.app=FBS.auth=FBS.fs=null; FBS.user=null; FBS.isMember=null; FBS.isAdmin=false; FBS.connecting=null; if(FBS.unsub){ try{ FBS.unsub(); }catch(e){} FBS.unsub=null; }
  const pid=currentConfig().projectId; if(MEM.meta.cacheProject&&MEM.meta.cacheProject!==pid){ LAB_IGNORED={cacheProject:MEM.meta.cacheProject,count:MEM.lab.size}; MEM.lab.clear(); } emit("status"); };
const isFile=()=>ENV.protocol!=null?ENV.protocol==="file:":(typeof location!=="undefined"&&location.protocol==="file:");
const online=()=>ENV.online!=null?ENV.online:(typeof navigator!=="undefined"&&navigator.onLine===false?false:true);
const F=()=>FBS.sdk.firestore;
const docRef=(...p)=>F().doc(FBS.fs,...p);
const colRef=(...p)=>F().collection(FBS.fs,...p);
const getDocS=r=>(F().getDocFromServer||F().getDoc)(r);     // 서버에서만 (오프라인 캐시가 빈 목록을 돌려주지 않게)
const getDocsS=r=>(F().getDocsFromServer||F().getDocs)(r);
const userOf=u=>u?{uid:u.uid,email:String(u.email||"").toLowerCase(),emailVerified:!!u.emailVerified,name:u.displayName||""}:null;
function fbErr(e,ctx){ const c=e&&e.code||""; const pre=ctx?ctx+": ":"";
  if(c==="permission-denied"||c==="firestore/permission-denied") return err("permission-denied",pre+"권한이 없습니다 (Firestore 규칙에 막힘)",{cause:e});
  if(c==="unavailable"||c==="firestore/unavailable"||c==="auth/network-request-failed") return err("unavailable",pre+"서버에 연결할 수 없습니다 (오프라인?) — 캐시를 그대로 씁니다",{cause:e});
  if(c==="not-found") return err("not-found",pre+"문서가 없습니다",{cause:e});
  if(c==="auth/popup-closed-by-user"||c==="auth/cancelled-popup-request") return err(c,"로그인 창이 닫혔습니다",{cause:e});
  if(c==="auth/unauthorized-domain") return err(c,"이 주소는 Firebase 에 허용된 도메인이 아닙니다 — https://inno-kaist.github.io/INNO-Analytics/ 로 여세요",{cause:e});
  return err(c||"firebase",pre+(e&&e.message||String(e)),{cause:e}); }
function fromFS(id,d){ const r=Object.assign({},d); for(const k of ["createdAt","updatedAt","verifiedAt","archivedAt","addedAt"]) if(k in r) r[k]=tsMs(r[k]); r.id=id; return clone(r); }

async function loadSDK(){ if(FBS.sdk) return FBS.sdk; if(SDK_OVERRIDE) return (FBS.sdk=SDK_OVERRIDE);
  const [app,auth,firestore]=await Promise.all(["firebase-app.js","firebase-auth.js","firebase-firestore.js"].map(f=>import(FB_BASE+f)));
  return (FBS.sdk={app,auth,firestore}); }
/* SDK 불러오기 + 앱 시작 + 이전 로그인 되살리기. 실패해도 예외 대신 상태(error) 로 */
DB.connect=async function(){ await DB.init(); if(FBS.app&&FBS.auth) return DB.firebaseStatus(); if(FBS.connecting) return FBS.connecting;
  FBS.connecting=(async()=>{ try{
      if(isFile()){ FBS.error=err("file-protocol",FILE_MSG); emit("status"); return DB.firebaseStatus(); }
      let s; try{ s=await loadSDK(); }catch(e){ FBS.error=err("sdk-load","Firebase 를 불러오지 못했습니다 (오프라인이거나 www.gstatic.com 에 못 닿음) — Built-in·Local·캐시로 계속합니다"); emit("status"); return DB.firebaseStatus(); }
      const cfg=currentConfig(); const name="inno-diffraction-"+cfg.projectId; const apps=s.app.getApps?s.app.getApps():[];
      FBS.app=apps.find(a=>a.name===name)||s.app.initializeApp(cfg,name); FBS.auth=s.auth.getAuth(FBS.app); FBS.fs=s.firestore.getFirestore(FBS.app);
      try{ await s.auth.getRedirectResult(FBS.auth); }catch(e){ FBS.error=fbErr(e,"로그인(redirect)"); }
      await new Promise(res=>{ let first=true; FBS.unsub=s.auth.onAuthStateChanged(FBS.auth,u=>{ const nu=userOf(u); const changed=(nu&&nu.uid)!==(FBS.user&&FBS.user.uid); FBS.user=nu;
          if(first){ first=false; res(); return; } if(changed) refreshMembership().then(()=>emit("status")); }); });
      if(FBS.user) await refreshMembership(); emit("status"); return DB.firebaseStatus();
    }finally{ FBS.connecting=null; } })();
  return FBS.connecting; };
async function refreshMembership(){ const u=FBS.user; FBS.isMember=null; FBS.isAdmin=false; FBS.admins=null; if(!u) return;
  if(!u.emailVerified){ FBS.isMember=false; FBS.error=err("email-unverified","이 Google 계정의 이메일이 확인되지 않았습니다 ("+u.email+")"); return; }
  try{ const snap=await getDocS(docRef("config","admins")); const em=(snap.exists()&&Array.isArray(snap.data().emails))?snap.data().emails.map(x=>String(x).toLowerCase()):[];
    FBS.isMember=true; FBS.admins=em; FBS.isAdmin=em.includes(u.email); FBS.error=null; }
  catch(e){ const c=e&&e.code||"";
    if(c==="permission-denied"){ FBS.isMember=false; FBS.error=err("not-member","구성원으로 등록되지 않았습니다 — 관리자에게 요청하세요 ("+u.email+")"); }
    else { FBS.error=fbErr(e,"구성원 확인"); const lm=MEM.meta.lastMember; if(lm&&lm.email===u.email){ FBS.isMember=lm.isMember; FBS.isAdmin=lm.isAdmin; FBS.admins=lm.admins||null; } return; } }
  await setMeta("lastMember",{email:u.email,isMember:FBS.isMember,isAdmin:FBS.isAdmin,admins:FBS.admins,at:now()}); }
DB.firebaseStatus=function(){ const cfg=currentConfig(); const u=FBS.user;
  return {projectId:cfg.projectId,expectedProjectId:EXPECT_PROJECT,configOk:cfg.projectId===EXPECT_PROJECT,configOverridden:!!(MEM.meta.fbConfig&&MEM.meta.fbConfig.projectId),
    sdkLoaded:!!FBS.sdk,connected:!!(FBS.app&&FBS.auth),fileProtocol:isFile(),online:online(),
    user:u?Object.assign({},u):null,isMember:FBS.isMember,isAdmin:!!FBS.isAdmin,admins:FBS.admins?FBS.admins.slice():null,
    lastSync:MEM.meta.lastSync||null,cacheProject:MEM.meta.cacheProject||null,labCount:MEM.lab.size,ignoredCache:LAB_IGNORED?Object.assign({},LAB_IGNORED):null,
    wasSignedIn:!!MEM.meta.wasSignedIn,storage:STORAGE,
    error:FBS.error?{code:FBS.error.code,message:FBS.error.message}:null}; };
DB.signIn=async function(){ await DB.init(); if(isFile()) throw err("file-protocol",FILE_MSG); await DB.connect(); if(!FBS.auth) throw FBS.error||err("not-connected","Firebase 에 연결하지 못했습니다");
  const A=FBS.sdk.auth; const prov=new A.GoogleAuthProvider(); if(prov.setCustomParameters) prov.setCustomParameters({prompt:"select_account"});
  try{ const r=await A.signInWithPopup(FBS.auth,prov); FBS.user=userOf(r&&r.user||FBS.auth.currentUser); }
  catch(e){ if(POPUP_FALLBACK.has(e&&e.code)){ await setMeta("wasSignedIn",true); await A.signInWithRedirect(FBS.auth,prov); return Object.assign(DB.firebaseStatus(),{redirect:true}); } throw fbErr(e,"로그인"); }
  await refreshMembership(); await setMeta("wasSignedIn",true); emit("status");
  if(FBS.isMember){ try{ await DB.sync(); }catch(e){} }
  return DB.firebaseStatus(); };
DB.signOut=async function(o){ o=o||{}; if(FBS.auth){ try{ await FBS.sdk.auth.signOut(FBS.auth); }catch(e){} } FBS.user=null; FBS.isMember=null; FBS.isAdmin=false; FBS.admins=null; FBS.error=null; await setMeta("wasSignedIn",false);
  if(o.clearCache){ MEM.lab.clear(); for(const [k,c] of MEM.cal) if(c.layer==="lab") MEM.cal.delete(k); await idbReplace("lab",[]); await idbReplace("cal",[],c=>c.layer!=="lab"); await setMeta("lastSync",null); emit("lab"); }
  emit("status"); };
function needMember(){ if(!FBS.app||!FBS.auth) throw FBS.error&&FBS.error.code==="file-protocol"?err("file-protocol",FILE_MSG):err("not-connected","Lab DB 에 연결되어 있지 않습니다 — Sign in 을 누르세요");
  if(!FBS.user) throw err("not-signed-in","로그인하지 않았습니다 — Lab DB 는 Google 로그인 후 씁니다");
  if(FBS.isMember===false) throw err("not-member","구성원으로 등록되지 않았습니다 — 관리자에게 요청하세요 ("+FBS.user.email+")"); }
function needWrite(){ needMember(); const pid=currentConfig().projectId; if(pid!==EXPECT_PROJECT) throw err("wrong-project","연결된 Firebase 프로젝트가 '"+pid+"' 입니다 (기대값 '"+EXPECT_PROJECT+"') — 쓰기를 막았습니다");
  if(FBS.isMember!==true) throw err("not-member","구성원 확인이 아직 안 됐습니다 — Sync 를 다시 누르세요"); }
function needAdmin(){ needWrite(); if(!FBS.isAdmin) throw err("not-admin","관리자만 할 수 있습니다"); }
function wrapF(ctx){ return e=>{ if(e&&e._own) throw e; throw fbErr(e,ctx); }; }

/* 내려받기: phases + calibrations → IndexedDB 캐시. 실패하면 캐시를 그대로 둔다 */
DB.sync=async function(){ await DB.connect(); needMember();
  let ph,cal; try{ [ph,cal]=await Promise.all([getDocsS(colRef("phases")),getDocsS(colRef("calibrations"))]); }
  catch(e){ const x=fbErr(e,"동기화"); if(x.code==="permission-denied"){ await refreshMembership(); emit("status"); throw FBS.isMember===false?err("not-member",FBS.error.message):x; } FBS.error=x; emit("status"); throw x; }
  if((ph.metadata&&ph.metadata.fromCache)||(cal.metadata&&cal.metadata.fromCache)){ const x=err("unavailable","서버 대신 Firestore 캐시에서 왔습니다 — 동기화하지 않았습니다"); FBS.error=x; emit("status"); throw x; }
  const recs=ph.docs.map(d=>fromFS(d.id,d.data())); const cals=cal.docs.map(d=>Object.assign(fromFS(d.id,d.data()),{layer:"lab"}));
  MEM.lab.clear(); for(const r of recs) MEM.lab.set(r.id,r); for(const [k,c] of MEM.cal) if(c.layer==="lab") MEM.cal.delete(k); for(const c of cals) MEM.cal.set(c.id,c); PH.clear();
  await idbReplace("lab",recs); await idbReplace("cal",cals,c=>c.layer!=="lab");
  LAB_IGNORED=null; await setMeta("cacheProject",currentConfig().projectId); await setMeta("lastSync",now()); FBS.error=null;
  emit("lab",{sync:true}); emit("status"); return {phases:recs.length,hidden:recs.filter(r=>r.hidden).length,calibrations:cals.length,lastSync:MEM.meta.lastSync}; };

function cacheLab(id,data){ const t=now(); const r=fromFS(id,data); for(const k of ["createdAt","updatedAt","verifiedAt"]) if(k in data&&r[k]==null) r[k]=t; MEM.lab.set(id,r); PH.delete(id); return idbPut("lab",r).then(()=>r); }
const stripU=o=>{ for(const k of Object.keys(o)) if(o[k]===undefined) delete o[k]; return o; };

/* 새 Lab 항목: 검사 → 중복 확인(모든 층) → version 1 */
DB.addLab=async function(rec,o){ o=o||{}; await DB.init(); needWrite(); const content=contentOf(rec); const v=DC.validateRecord(content); if(!v.ok) return {ok:false,reason:"invalid",validation:v};
  if(!o.force){ const dup=DB.findDuplicate(content,{exclude:o.exclude}); if(dup) return {ok:false,reason:"duplicate",duplicate:dup,validation:v}; }
  const S=F(); const u=FBS.user; const id="lab:"+rid();
  const data=Object.assign(content,{id,createdBy:u.uid,createdByEmail:u.email,createdAt:S.serverTimestamp(),updatedAt:S.serverTimestamp(),updatedBy:u.uid,updatedByEmail:u.email,version:1,hidden:false,verified:false});
  if(rec.promotedFrom) data.promotedFrom=String(rec.promotedFrom);
  if(FBS.isAdmin&&rec.verified===true) Object.assign(data,{verified:true,verifiedBy:u.uid,verifiedByEmail:u.email,verifiedAt:S.serverTimestamp(),verifiedVersion:1});
  await S.setDoc(docRef("phases",id),data).catch(wrapF("Lab 에 추가"));
  const r=await cacheLab(id,data); emit("lab",{ids:[id]}); return {ok:true,rec:tag(r,"lab"),validation:v}; };
/* 트랜잭션: 옛 내용 → history/{옛 version} (create), 새 내용 version+1.
   ch.content = 새 내용(없으면 옛 내용 그대로), ch.meta(old) = 바꿀 부가 필드 (hidden·verified…) */
async function labMutate(id,ch,why){ needWrite(); const S=F(); const u=FBS.user; const ref=docRef("phases",id); let next=null;
  await S.runTransaction(FBS.fs,async t=>{ const snap=await t.get(ref); if(!snap.exists()) throw err("not-found","Lab 에 그 항목이 없습니다: "+id); const old=snap.data();
    if(!FBS.isAdmin&&old.createdBy!==u.uid) throw err("not-owner","이 항목은 만든 사람("+(old.createdByEmail||"?")+")이나 관리자만 고칠 수 있습니다");
    next=ch.content?contentOf(ch.content):contentOf(old); for(const k of VER_META) if(old[k]!==undefined) next[k]=old[k];
    if(ch.meta) Object.assign(next,ch.meta(old)||{});
    if(!FBS.isAdmin) for(const k of ["verified","verifiedBy","verifiedByEmail","verifiedAt","verifiedVersion"]){ if(old[k]===undefined) delete next[k]; else next[k]=old[k]; }
    Object.assign(next,{id,createdBy:old.createdBy,createdByEmail:old.createdByEmail,createdAt:old.createdAt,version:(old.version|0)+1,updatedAt:S.serverTimestamp(),updatedBy:u.uid,updatedByEmail:u.email});
    if(next.hidden===undefined) next.hidden=false; stripU(next);
    const hist=stripU(Object.assign({},old,{archivedAt:S.serverTimestamp(),archivedBy:u.uid,archivedByEmail:u.email,archivedReason:why||"edit"}));
    t.set(docRef("phases",id,"history",String(old.version|0)),hist);
    t.set(ref,next); }).catch(wrapF("Lab 고치기"));
  const r=await cacheLab(id,next); emit("lab",{ids:[id]}); return tag(r,"lab"); }
DB.updateLab=async function(id,rec,o){ o=o||{}; const content=contentOf(rec); const v=DC.validateRecord(content); if(!v.ok) return {ok:false,reason:"invalid",validation:v};
  if(o.checkDuplicate){ const dup=DB.findDuplicate(content,{exclude:[id]}); if(dup) return {ok:false,reason:"duplicate",duplicate:dup,validation:v}; }
  const r=await labMutate(id,{content},"edit"); return {ok:true,rec:r,validation:v}; };
DB.hideLab=(id,hidden)=>labMutate(id,{meta:()=>({hidden:hidden!==false})},hidden===false?"unhide":"hide");
DB.setVerified=function(id,on){ try{ needAdmin(); }catch(e){ return Promise.reject(e); } const u=FBS.user; const S=F();
  return labMutate(id,{meta:old=>on?{verified:true,verifiedBy:u.uid,verifiedByEmail:u.email,verifiedAt:S.serverTimestamp(),verifiedVersion:(old.version|0)+1}
    :{verified:false,verifiedBy:u.uid,verifiedByEmail:u.email,verifiedAt:S.serverTimestamp(),verifiedVersion:null}},on?"verify":"unverify"); };
/* Local → Lab */
DB.promoteLocal=async function(id,o){ o=o||{}; await DB.init(); const r=MEM.local.get(id); if(!r) throw err("not-found","Local 항목이 없습니다: "+id);
  const res=await DB.addLab(Object.assign(contentOf(r),{promotedFrom:id}),{force:o.force,exclude:[id]}); if(!res.ok) return res;
  r.promotedTo=res.rec.id; await idbPut("local",r); emit("local",{ids:[id]}); return res; };
/* Lab 항목의 지난 판본 (서버) */
DB.labHistory=async function(id){ await DB.connect(); needMember(); const s=await getDocsS(colRef("phases",id,"history")).catch(wrapF("판본 기록"));
  return s.docs.map(d=>{ const x=fromFS(id,d.data()); return {version:x.version!=null?x.version:+d.id,updatedAt:x.updatedAt,updatedByEmail:x.updatedByEmail||x.createdByEmail||"",archivedAt:x.archivedAt,archivedByEmail:x.archivedByEmail||"",archivedReason:x.archivedReason||"",rec:x}; }).sort((a,b)=>b.version-a.version); };
/* 판정 기록이 여전히 지금 내용과 같은지 */
DB.verifyState=function(id){ const r=DB.get(id); if(!r||r.layer!=="lab") return null; if(!r.verified) return {state:"unverified"}; if(r.verifiedVersion!=null&&r.version>r.verifiedVersion) return {state:"edited",verifiedVersion:r.verifiedVersion,version:r.version}; return {state:"verified",verifiedVersion:r.verifiedVersion}; };

/* ---------------- 보정(calibrations) ---------------- */
DB.listCalibrations=function(o){ o=o||{}; const layers=o.layers||["lab","local"]; return [...MEM.cal.values()].filter(c=>layers.includes(c.layer||"local")).map(clone).sort((a,b)=>(b.updatedAt||0)-(a.updatedAt||0)); };
const CAL_FIELDS=["name","microscope","kV","cameraLength","mode","scale","ellip","phi","sigmaRel","standard","standardName","sourceFile","note"];
function calContent(c){ const o={}; for(const k of CAL_FIELDS){ const v=c[k]; if(v===undefined||v===null||v==="") continue; o[k]=["kV","cameraLength","scale","ellip","phi","sigmaRel"].includes(k)?num(v):String(v); } if(c.standardRef) o.standardRef=clone(c.standardRef); return stripU(o); }
DB.saveCalibration=async function(cal,o){ o=o||{}; await DB.init(); const c=calContent(cal||{}); if(!c.name) throw err("invalid","보정 이름이 필요합니다"); if(!(c.scale>0)) throw err("invalid","scale (1/nm/px) 은 양수여야 합니다");
  if(c.standard&&!c.standardRef){ const ref=DB.refOf(c.standard); if(ref) c.standardRef=ref; }
  if(!o.lab){ const old=cal.id?MEM.cal.get(cal.id):null; if(old&&old.layer==="lab") throw err("bad-layer","Lab 보정은 {lab:true} 로 고칩니다"); const t=now();
    const r=Object.assign(c,{id:old?old.id:"local:"+rid(),layer:"local",createdAt:old?old.createdAt:t,updatedAt:t}); MEM.cal.set(r.id,r); await idbPut("cal",r); emit("cal",{ids:[r.id]}); return clone(r); }
  needWrite(); const S=F(); const u=FBS.user; const old=cal.id&&/^lab:/.test(cal.id)?MEM.cal.get(cal.id):null;
  let id, data;
  if(old){ if(!FBS.isAdmin&&old.createdBy!==u.uid) throw err("not-owner","이 보정은 만든 사람이나 관리자만 고칠 수 있습니다"); id=old.id;
    data=Object.assign(c,{id,createdBy:old.createdBy,createdByEmail:old.createdByEmail||"",updatedAt:S.serverTimestamp(),updatedBy:u.uid,updatedByEmail:u.email}); if(old.createdAt) data.createdAt=new Date(old.createdAt); }
  else { id="lab:"+rid(); data=Object.assign(c,{id,createdBy:u.uid,createdByEmail:u.email,createdAt:S.serverTimestamp(),updatedAt:S.serverTimestamp()}); }
  await S.setDoc(docRef("calibrations",id),data).catch(wrapF("보정 저장"));
  const t=now(); const r=Object.assign(fromFS(id,data),{layer:"lab"}); if(r.createdAt==null) r.createdAt=t; if(r.updatedAt==null) r.updatedAt=t; MEM.cal.set(id,r); await idbPut("cal",r); emit("cal",{ids:[id]}); return clone(r); };

/* ---------------- 구성원·관리자 (관리자) ---------------- */
const EMAIL_RE=/^[^\s@]+@[^\s@]+\.[^\s@]+$/;
DB.listMembers=async function(){ await DB.connect(); needMember(); const s=await getDocsS(colRef("members")).catch(wrapF("구성원 목록"));
  return s.docs.map(d=>Object.assign({email:d.id},fromFS(d.id,d.data()))).map(m=>{ delete m.id; return m; }).sort((a,b)=>a.email.localeCompare(b.email)); };
DB.addMember=async function(email,name){ needAdmin(); const e=String(email||"").trim().toLowerCase(); if(!EMAIL_RE.test(e)) throw err("invalid","이메일 형식이 아닙니다: "+email); const S=F(); const ref=docRef("members",e);
  const snap=await getDocS(ref).catch(wrapF("구성원 확인"));
  if(snap.exists()) await S.updateDoc(ref,{name:String(name||snap.data().name||""),active:true,changedBy:FBS.user.email,changedAt:S.serverTimestamp()}).catch(wrapF("구성원 추가"));
  else await S.setDoc(ref,{name:String(name||""),active:true,addedBy:FBS.user.email,addedAt:S.serverTimestamp()}).catch(wrapF("구성원 추가"));
  emit("members"); return {email:e,name:String(name||""),active:true}; };
DB.setMemberActive=async function(email,on){ needAdmin(); const S=F(); const e=String(email).toLowerCase(); await S.updateDoc(docRef("members",e),{active:!!on,changedBy:FBS.user.email,changedAt:S.serverTimestamp()}).catch(wrapF("구성원 상태")); emit("members"); return true; };
DB.listAdmins=async function(){ await DB.connect(); needMember(); const s=await getDocS(docRef("config","admins")).catch(wrapF("관리자 목록")); const em=s.exists()&&Array.isArray(s.data().emails)?s.data().emails.slice():[]; FBS.admins=em.map(x=>String(x).toLowerCase()); return em; };
async function editAdmins(fn){ needAdmin(); const S=F(); const ref=docRef("config","admins"); let out;
  await S.runTransaction(FBS.fs,async t=>{ const s=await t.get(ref); const cur=s.exists()&&Array.isArray(s.data().emails)?s.data().emails.slice():[]; out=fn(cur); if(!out.length) throw err("last-admin","마지막 관리자 이메일은 지울 수 없습니다"); t.set(ref,Object.assign({},s.exists()?s.data():{},{emails:out})); }).catch(wrapF("관리자 목록"));
  FBS.admins=out.map(x=>String(x).toLowerCase()); FBS.isAdmin=FBS.admins.includes(FBS.user.email); emit("status"); return out; }
DB.addAdmin=email=>{ const e=String(email||"").trim().toLowerCase(); if(!EMAIL_RE.test(e)) return Promise.reject(err("invalid","이메일 형식이 아닙니다")); return editAdmins(cur=>cur.map(x=>x.toLowerCase()).includes(e)?cur:cur.concat([e])); };
DB.removeAdmin=email=>{ const e=String(email||"").trim().toLowerCase(); return editAdmins(cur=>{ const n=cur.filter(x=>x.toLowerCase()!==e); if(!n.length) throw err("last-admin","마지막 관리자 이메일은 지울 수 없습니다"); return n; }); };

/* ---------------- 시험용 ---------------- */
DB._setIDB=f=>{ IDB_OVERRIDE=f; };
DB._setSDK=s=>{ SDK_OVERRIDE=s; FBS.sdk=null; };
DB._setEnv=e=>{ ENV=Object.assign({},e||{}); };
DB._reset=async function(){ try{ const d=await dbp; if(d) d.close(); }catch(e){} dbp=null; INIT=null; STORAGE="memory"; LAB_IGNORED=null; MEM.local.clear(); MEM.lab.clear(); MEM.cal.clear(); MEM.meta={}; PH.clear();
  if(FBS.unsub){ try{ FBS.unsub(); }catch(e){} } Object.assign(FBS,{sdk:null,app:null,auth:null,fs:null,user:null,isMember:null,isAdmin:false,admins:null,error:null,connecting:null,unsub:null}); PACK.src=null; PACK.tried=false; PACK.recs=new Map(); PACK.version=null; };
DB.EXPECT_PROJECT=EXPECT_PROJECT; DB.DEFAULT_CONFIG=DEFAULT_CONFIG; DB.FB_BASE=FB_BASE;

if(typeof module!=="undefined"&&module.exports) module.exports=DB;
G.DIFF_DB=DB;
})();
