/* ============================================================
   INNO Analytics — Diffraction · Phase DB 화면 (DIFF_DB_UI.render(container, opts))
   자료는 DIFF_DB (src/diff_db.js), 계산은 DIFF_CORE. 화면 뼈대는 영어, 설명·알림은 한국어 (docs/UI-English-style.md).
   opts: {onUse(ids, refs)?  — "Use selected" 단추 (Diffraction 판정 후보로 넘기기), autoConnect?: 이전에 로그인했으면 SDK 를 불러 되살림 (기본 true)}
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const DB=()=>G.DIFF_DB, DC=()=>G.DIFF_CORE;
const esc=s=>String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const $q=(r,s)=>r&&r.querySelector(s), $a=(r,s)=>r?[...r.querySelectorAll(s)]:[];
const LAYER={pack:"Built-in",lab:"Lab",local:"Local"};
const KINDS=["COD","DOI","paper","measured","manual"];
const fx=(v,d)=>v==null||!isFinite(+v)?"—":(+v).toFixed(d==null?4:d);
const when=t=>{ if(!t) return "—"; const d=new Date(t); if(isNaN(d)) return "—"; const p=n=>String(n).padStart(2,"0"); return `${d.getFullYear()}-${p(d.getMonth()+1)}-${p(d.getDate())} ${p(d.getHours())}:${p(d.getMinutes())}`; };
const ago=t=>{ if(!t) return "never"; const s=(Date.now()-t)/1000; if(s<60) return "just now"; if(s<3600) return Math.round(s/60)+" min ago"; if(s<86400) return Math.round(s/3600)+" h ago"; return when(t); };
const PREF="inno-diff-db-ui";

/* ---------------- 알림 ---------------- */
function say(m,k){ try{ if(typeof G.toast==="function"){ G.toast(m,k==="err"?"err":k); return; } }catch(e){}
  if(typeof document==="undefined") return; let box=document.querySelector(".ddb-toasts"); if(!box){ box=document.createElement("div"); box.className="ddb-toasts"; document.body.appendChild(box); }
  const t=document.createElement("div"); t.className="t "+(k||""); t.textContent=m; box.appendChild(t); setTimeout(()=>t.remove(),k==="err"?6000:3200); }
const errMsg=e=>String(e&&e.message||e);
function download(name,text,mime){ const b=new Blob([text],{type:mime||"application/json"}); const u=URL.createObjectURL(b); const a=document.createElement("a"); a.href=u; a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(u); a.remove(); },500); }
const stamp=()=>new Date().toISOString().slice(0,10).replace(/-/g,"");

/* ---------------- 화면 상태 ---------------- */
let ROOT=null, OPTS={}, OFF=null, MEMBERS=null, HIST={}, BUSY=false;
const U={layers:["pack","lab","local"],q:"",els:[],elMode:"subset",showHidden:false,sel:null,checked:[],limit:300};
try{ const p=JSON.parse(localStorage.getItem(PREF)||"null"); if(p&&typeof p==="object"){ if(Array.isArray(p.layers)) U.layers=p.layers.filter(x=>LAYER[x]); if(Array.isArray(p.els)) U.els=p.els; if(p.elMode) U.elMode=p.elMode; } }catch(e){}
function savePref(){ try{ localStorage.setItem(PREF,JSON.stringify({layers:U.layers,els:U.els,elMode:U.elMode})); }catch(e){} }

function render(container,opts){ ROOT=container; OPTS=opts||{}; if(OFF){ OFF(); OFF=null; }
  container.innerHTML=`<div class="ddb">
    <div class="card ddb-status" id="ddbStatus"></div>
    <div class="card ddb-bar" id="ddbBar"></div>
    <div class="ddb-main">
      <div class="card ddb-listc" id="ddbList"></div>
      <div class="card ddb-detail" id="ddbDetail"></div>
    </div>
    <div class="card" id="ddbCal" style="margin-top:12px"></div>
    <div class="card" id="ddbMem" hidden></div>
  </div>`;
  const D=DB(); if(!D){ $q(container,"#ddbStatus").innerHTML=`<div class="note danger">DIFF_DB 가 없습니다 (diff_db.js 를 먼저 넣어야 합니다)</div>`; return; }
  OFF=D.onChange(ev=>{ if(!ROOT||!ROOT.isConnected) return; if(ev.type==="status"){ drawStatus(); drawDetail(); drawMembers(); } else if(ev.type==="cal") drawCal(); else { drawStatus(); drawBar(); drawList(); drawDetail(); drawCal(); } });
  drawAll();
  D.init().then(()=>{ if(ROOT!==container) return; drawAll(); const st=D.firebaseStatus();
    if(OPTS.autoConnect!==false&&st.wasSignedIn&&!st.fileProtocol) D.connect().then(s=>{ if(s&&s.isMember&&(!s.lastSync||Date.now()-s.lastSync>10*60*1000)) D.sync().catch(()=>{}); }).catch(()=>{}); });
}
function drawAll(){ drawStatus(); drawBar(); drawList(); drawDetail(); drawCal(); drawMembers(); }

/* ---------------- 상태 (Lab DB) ---------------- */
function drawStatus(){ const el=$q(ROOT,"#ddbStatus"); if(!el) return; const D=DB(); const s=D.firebaseStatus(); const c=D.counts();
  const role=!s.user?"":s.isMember===null?"Checking…":s.isAdmin?"Admin":s.isMember?"Member":"Not a member";
  const dot=s.user&&s.isMember?"on":s.user?"warn":"";
  const notes=[];
  if(!s.configOk) notes.push(`<div class="note danger"><b>Firebase 프로젝트가 다릅니다:</b> 지금 설정은 <span class="mono">${esc(s.projectId)}</span>, 이 앱의 공유 DB 는 <span class="mono">${esc(s.expectedProjectId)}</span> 입니다. 다른 프로젝트에 몰래 쓰지 않도록 쓰기를 막았습니다. <button class="btn sm" data-a="resetcfg">Reset config</button></div>`);
  if(s.ignoredCache) notes.push(`<div class="note warn">다른 프로젝트(<span class="mono">${esc(s.ignoredCache.cacheProject)}</span>)에서 받은 Lab 캐시 ${s.ignoredCache.count} 개는 섞지 않고 숨겼습니다.</div>`);
  if(s.storage==="memory") notes.push(`<div class="note warn">이 브라우저에서 IndexedDB 를 쓸 수 없어 Local 항목과 Lab 캐시가 창을 닫으면 사라집니다 (사생활 보호 창?). Export JSON 으로 따로 보관하세요.</div>`);
  if(s.fileProtocol) notes.push(`<div class="note warn">Google 로그인은 웹 주소(<span class="mono">https://inno-kaist.github.io/INNO-Analytics/</span>)로 열었을 때만 됩니다. 파일로 연 지금은 Built-in·Local 과 마지막으로 받은 Lab 캐시만 씁니다.</div>`);
  else if(s.error&&s.error.code!=="file-protocol") notes.push(`<div class="note ${s.error.code==="not-member"||s.error.code==="email-unverified"?"danger":"warn"}">${esc(s.error.message)}</div>`);
  let line;
  if(!s.user){ line=`<span><span class="dot ${s.connected?"warn":""}"></span><span class="k">Lab DB</span><span class="v" id="ddbStLine">로그인하지 않음 — Built-in·Local${c.lab?` 과 Lab 캐시(${c.lab} 개)`:""} 만 씁니다</span></span>`; }
  else line=`<span><span class="dot ${dot}"></span><span class="k">Account</span><span class="v brk" id="ddbStLine">${esc(s.user.email)}</span></span><span><span class="k">Role</span><span class="v">${esc(role)}</span></span>`;
  el.innerHTML=`<div class="ddb-h"><h3>Phase DB</h3><span class="muted">Built-in ${c.pack} · Lab ${c.lab}${c.labHidden?` (hidden ${c.labHidden})`:""} · Local ${c.local}</span></div>
    <div class="ddb-st">${line}
      <span><span class="k">Project</span><span class="v mono ${s.configOk?"":"bad"}" id="ddbProj">${esc(s.projectId)}</span></span>
      <span><span class="k">Last sync</span><span class="v" title="${esc(when(s.lastSync))}">${esc(ago(s.lastSync))}</span></span>
      <span class="acts">${s.user?`<button class="btn sm" data-a="sync" ${s.isMember===false?"disabled":""}>Sync</button><button class="btn sm" data-a="signout">Sign out</button>`
        :`<button class="btn sm primary" data-a="signin" ${s.fileProtocol?"disabled title=\"웹 주소로 열어야 로그인됩니다\"":""}>Sign in with Google</button>`}</span>
    </div><div class="ddb-notes">${notes.join("")}</div>`;
  /* 단추에 손이 가면 SDK 를 미리 불러 둔다 — 누를 때 바로 팝업을 열어야 브라우저가 막지 않는다 */
  const si=$q(el,"button[data-a=signin]"); if(si&&!si.disabled){ const pre=()=>{ D.connect().catch(()=>{}); }; si.addEventListener("pointerenter",pre,{once:true}); si.addEventListener("touchstart",pre,{once:true,passive:true}); si.addEventListener("focus",pre,{once:true}); }
  el.onclick=async e=>{ const b=e.target.closest("button[data-a]"); if(!b||BUSY) return; const a=b.dataset.a; BUSY=true; b.disabled=true;
    try{ if(a==="signin"){ const r=await D.signIn(); if(r.redirect) say("팝업이 막혀 로그인 페이지로 이동합니다…"); else if(r.isMember) say("로그인했습니다 — Lab DB 를 받았습니다","ok"); else if(r.isMember===false) say("구성원으로 등록되지 않았습니다 — 관리자에게 요청하세요 ("+r.user.email+")","warn"); }
      else if(a==="signout"){ await D.signOut(); say("로그아웃했습니다 — Lab 캐시는 이 기기에 남습니다"); }
      else if(a==="sync"){ const r=await D.sync(); say(`동기화: 상 ${r.phases} 개 (숨김 ${r.hidden}) · 보정 ${r.calibrations} 개`,"ok"); }
      else if(a==="resetcfg"){ await D.setConfig(null); say("기본 Firebase 설정으로 돌렸습니다","ok"); }
    }catch(x){ say(errMsg(x),"err"); } finally{ BUSY=false; drawStatus(); } }; }

/* ---------------- 도구 막대 ---------------- */
function drawBar(){ const el=$q(ROOT,"#ddbBar"); if(!el) return; const D=DB(); const c=D.counts(); const s=D.firebaseStatus();
  el.innerHTML=`<div class="grp"><span class="lbl">Layers</span><div class="seg ddb-layers">${["pack","lab","local"].map(L=>`<button data-l="${L}" class="${U.layers.includes(L)?"active":""}" aria-pressed="${U.layers.includes(L)}">${LAYER[L]} <span class="n">${L==="lab"?c.lab-(U.showHidden?0:c.labHidden):c[L]}</span></button>`).join("")}</div></div>
    <div class="grp q"><input type="search" id="ddbQ" placeholder="Search name, formula, COD id" value="${esc(U.q)}" aria-label="Search"></div>
    <div class="grp"><span class="lbl">Elements</span><span class="chips" id="ddbChips">${U.els.map(e=>`<span class="chip">${esc(e)}<button data-x="${esc(e)}" title="Remove ${esc(e)}" aria-label="Remove ${esc(e)}">×</button></span>`).join("")}</span>
      <input type="text" class="ddb-elin" id="ddbEl" placeholder="e.g. Na Y F" aria-label="Add elements" title="원소 기호를 띄어 쓰고 Enter (EDS 에서 본 원소)">
      <select id="ddbElMode" aria-label="Element match" title="Only these = 상의 원소가 모두 고른 원소 안에 (REDI 방식). Contains all = 고른 원소를 모두 포함">
        <option value="subset" ${U.elMode==="subset"?"selected":""}>Only these</option><option value="all" ${U.elMode==="all"?"selected":""}>Contains all</option></select>
      ${U.els.length?`<button class="btn sm" id="ddbElClr">Clear</button>`:""}</div>
    <label class="chk"><input type="checkbox" id="ddbHid" ${U.showHidden?"checked":""}>Show hidden</label>
    <div class="grp" style="margin-left:auto">
      <button class="btn sm primary" id="ddbAdd">Add phase</button>
      <label class="btn sm" style="cursor:pointer" title="Phase DB 백업(JSON)을 Local 로 가져오기">Import JSON<input type="file" id="ddbImp" accept=".json,application/json" hidden></label>
      <button class="btn sm" id="ddbExpJ" title="고른 항목(없으면 Lab·Local 전부)을 JSON 으로">Export JSON</button>
      <button class="btn sm" id="ddbExpC" title="목록 요약 CSV">Export CSV</button></div>`;
  $a(el,".ddb-layers button").forEach(b=>b.onclick=()=>{ const L=b.dataset.l; U.layers=U.layers.includes(L)?U.layers.filter(x=>x!==L):U.layers.concat([L]); savePref(); drawBar(); drawList(); });
  const q=$q(el,"#ddbQ"); let tq=null; q.oninput=()=>{ clearTimeout(tq); tq=setTimeout(()=>{ U.q=q.value; drawList(); },150); };
  const ei=$q(el,"#ddbEl"); const addEls=()=>{ const v=D.parseElements(ei.value); const bad=String(ei.value).split(/[\s,;]+/).filter(Boolean).filter(t=>!v.some(x=>x.toLowerCase()===t.toLowerCase()));
      if(bad.length) say("모르는 원소 기호: "+bad.join(" "),"warn"); if(v.length){ U.els=[...new Set(U.els.concat(v))]; savePref(); drawBar(); drawList(); const n=$q(ROOT,"#ddbEl"); if(n) n.focus(); } else ei.value=""; };
  ei.onkeydown=e=>{ if(e.key==="Enter"||e.key===","){ e.preventDefault(); addEls(); } else if(e.key==="Backspace"&&!ei.value&&U.els.length){ U.els.pop(); savePref(); drawBar(); drawList(); $q(ROOT,"#ddbEl").focus(); } };
  ei.onblur=()=>{ if(ei.value.trim()) addEls(); };
  $a(el,"#ddbChips button").forEach(b=>b.onclick=()=>{ U.els=U.els.filter(x=>x!==b.dataset.x); savePref(); drawBar(); drawList(); });
  $q(el,"#ddbElMode").onchange=e=>{ U.elMode=e.target.value; savePref(); drawList(); };
  const cl=$q(el,"#ddbElClr"); if(cl) cl.onclick=()=>{ U.els=[]; savePref(); drawBar(); drawList(); };
  $q(el,"#ddbHid").onchange=e=>{ U.showHidden=e.target.checked; drawBar(); drawList(); };
  $q(el,"#ddbAdd").onclick=()=>openEditor(null);
  $q(el,"#ddbImp").onchange=async e=>{ const f=e.target.files&&e.target.files[0]; e.target.value=""; if(!f) return;
    try{ const r=await D.importJSON(await f.text()); say(`가져오기: 추가 ${r.added} · 갱신 ${r.updated} · 같음 ${r.unchanged} · 건너뜀 ${r.skipped.length}${r.errors.length?` · 오류 ${r.errors.length}`:""}${r.calibrations?` · 보정 ${r.calibrations}`:""}`,r.errors.length?"warn":"ok"); }catch(x){ say(errMsg(x),"err"); } };
  $q(el,"#ddbExpJ").onclick=()=>{ const ids=U.checked.filter(id=>D.get(id)); download(`phasedb_${ids.length?ids.length+"sel":"all"}_${stamp()}.json`,D.exportJSON(ids.length?ids:null)); };
  $q(el,"#ddbExpC").onclick=()=>{ const ids=U.checked.filter(id=>D.get(id)); download(`phasedb_${stamp()}.csv`,"﻿"+D.exportCSV(ids.length?ids:currentList().map(r=>r.id)),"text/csv"); };
  void s; }
function currentList(){ return DB().list({layers:U.layers,includeHidden:U.showHidden,query:U.q,elements:U.els,elementMode:U.elMode}); }

/* ---------------- 목록 ---------------- */
function statusBadge(r){ if(r.hidden) return `<span class="badge gray">Hidden</span>`; if(r.layer!=="lab") return r.layer==="pack"&&r.verified?`<span class="badge ok">Verified</span>`:""; const v=DB().verifyState(r.id);
  return v&&v.state==="verified"?`<span class="badge ok">Verified</span>`:v&&v.state==="edited"?`<span class="badge warn" title="v${esc(v.verifiedVersion)} 에서 검증, 그 뒤 v${esc(v.version)} 로 편집">Edited</span>`:`<span class="badge gray">Unverified</span>`; }
const srcShort=s=>{ s=s||{}; return s.codId?"COD "+s.codId:s.doi?"DOI":s.kind||"—"; };
function drawList(){ const el=$q(ROOT,"#ddbList"); if(!el) return; const D=DB(); const L=currentList(); const shown=L.slice(0,U.limit); const chk=new Set(U.checked);
  if(U.sel&&!D.get(U.sel)) U.sel=null;
  if(!L.length){ const any=D.counts(); const none=!any.pack&&!any.lab&&!any.local;
    el.innerHTML=`<div class="empty"><b>${none?"No phases yet":"No matching phases"}</b>${none?"Add phase 로 CIF 를 넣거나 매개변수로 만드세요. Lab DB 는 로그인 후 Sync 로 받습니다.":"검색어·원소·층(Layers) 조건을 바꿔 보세요."}</div>`; return; }
  el.innerHTML=`<div class="tblwrap"><table class="tbl ddb-tbl"><thead><tr><th class="cb"><label title="Select all shown"><input type="checkbox" id="ddbAll" aria-label="Select all"></label></th><th class="l">Name</th><th class="l">Formula</th><th>SG</th><th>System</th><th>a (nm)</th><th>b (nm)</th><th>c (nm)</th><th class="l">Source</th><th>Status</th><th>Layer</th></tr></thead><tbody>
    ${shown.map(r=>{ const n=r.sg&&+r.sg.number; const c=r.cell||{}; return `<tr data-id="${esc(r.id)}" class="${r.id===U.sel?"sel":""} ${r.hidden?"hid":""}">
      <td class="cb"><label><input type="checkbox" data-c="${esc(r.id)}" ${chk.has(r.id)?"checked":""} aria-label="Select ${esc(r.name)}"></label></td>
      <td class="l nm" title="${esc(r.name)}">${esc(r.name||"—")}</td><td class="l fm" title="${esc(r.formula)}">${esc(r.formula||"—")}</td>
      <td title="${esc(r.sg&&r.sg.symbol)}">${n?`${n} <span class="muted">${esc(r.sg.symbol||"")}</span>`:esc(r.sg&&r.sg.symbol||"—")}</td><td>${n?esc(DC().crystalSystem(n)):"—"}</td>
      <td class="num">${fx(c.a)}</td><td class="num">${fx(c.b)}</td><td class="num">${fx(c.c)}</td>
      <td class="l">${esc(srcShort(r.source))}</td><td>${statusBadge(r)}</td><td><span class="lay ${r.layer}">${LAYER[r.layer]}</span></td></tr>`; }).join("")}
    </tbody></table></div>
    <div class="ddb-foot"><span>${L.length} phase${L.length===1?"":"s"}${L.length>shown.length?` · showing ${shown.length}`:""}${U.checked.length?` · ${U.checked.length} selected`:""}</span><span class="sp"></span>
      ${L.length>shown.length?`<button class="btn sm" id="ddbMore">Show more</button>`:""}
      ${OPTS.onUse?`<button class="btn sm primary" id="ddbUse" ${U.checked.length?"":"disabled"}>Use selected (${U.checked.length})</button>`:""}
      ${U.checked.length?`<button class="btn sm" id="ddbUnchk">Clear selection</button>`:""}</div>`;
  const tb=$q(el,"tbody");
  tb.onclick=e=>{ if(e.target.closest("label")||e.target.matches("input")) return; const tr=e.target.closest("tr[data-id]"); if(!tr) return; U.sel=tr.dataset.id; $a(tb,"tr.sel").forEach(x=>x.classList.remove("sel")); tr.classList.add("sel"); drawDetail();
    if(G.matchMedia&&G.matchMedia("(max-width:1000px)").matches){ const d=$q(ROOT,"#ddbDetail"); if(d&&d.scrollIntoView) d.scrollIntoView({behavior:"smooth",block:"start"}); } };
  tb.onchange=e=>{ const id=e.target.dataset.c; if(!id) return; U.checked=e.target.checked?[...new Set(U.checked.concat([id]))]:U.checked.filter(x=>x!==id); drawList(); };
  $q(el,"#ddbAll").onchange=e=>{ const ids=shown.map(r=>r.id); U.checked=e.target.checked?[...new Set(U.checked.concat(ids))]:U.checked.filter(x=>!ids.includes(x)); drawList(); };
  const m=$q(el,"#ddbMore"); if(m) m.onclick=()=>{ U.limit+=300; drawList(); };
  const us=$q(el,"#ddbUse"); if(us) us.onclick=()=>{ const ids=U.checked.filter(id=>D.get(id)); OPTS.onUse(ids,ids.map(D.refOf)); };
  const uc=$q(el,"#ddbUnchk"); if(uc) uc.onclick=()=>{ U.checked=[]; drawList(); }; }

/* ---------------- 상세 ---------------- */
const hklStr=h=>h.map(x=>x<0?"−"+(-x):String(x)).join(Math.max(...h.map(Math.abs))>9?" ":"");
function repHkl(hkls){ return hkls.slice().sort((a,b)=>{ const na=a.filter(x=>x<0).length, nb=b.filter(x=>x<0).length; if(na!==nb) return na-nb; for(let i=0;i<3;i++) if(a[i]!==b[i]) return b[i]-a[i]; return 0; })[0]; }
const RING_CACHE=new Map();
function rings(id){ const D=DB(); const r=D.get(id); const key=id+"|"+(r&&r.version)+"|"+(r&&r.updatedAt); if(RING_CACHE.has(key)) return RING_CACHE.get(key);
  let out; try{ const ph=D.phase(id); const L=DC().ringList(ph,10); const Im=Math.max(0,...L.map(x=>x.I||0)); out={rings:L.map(x=>({d:x.d,g:x.g,m:x.m,hkl:repHkl(x.hkls),I:x.I!=null&&Im>0?100*x.I/Im:null})),hasI:!!(ph.atoms&&ph.atoms.length),V:ph.cell.V}; }catch(e){ out={error:errMsg(e)}; }
  if(RING_CACHE.size>60) RING_CACHE.clear(); RING_CACHE.set(key,out); return out; }
function canEdit(r,s){ if(r.layer==="local") return true; if(r.layer!=="lab"||!s.user||!s.isMember||!s.configOk) return false; return s.isAdmin||r.createdBy===s.user.uid; }
function drawDetail(){ const el=$q(ROOT,"#ddbDetail"); if(!el) return; const D=DB(); const r=U.sel?D.get(U.sel):null;
  if(!r){ el.innerHTML=`<div class="empty"><b>Phase Detail</b>목록에서 상을 고르면 격자·공간군·원자·출처와 링(d, hkl) 표가 여기에 나옵니다.</div>`; return; }
  const s=D.firebaseStatus(); const n=r.sg&&+r.sg.number; const c=r.cell||{}; const src=r.source||{}; const R=rings(r.id); const ed=canEdit(r,s); const member=!!(s.user&&s.isMember&&s.configOk);
  const vs=r.layer==="lab"?D.verifyState(r.id):null; const els=D.elementsOf(r);
  const acts=[];
  if(r.layer==="local"){ acts.push(`<button class="btn sm" data-a="edit">Edit</button>`); acts.push(`<button class="btn sm" data-a="promote" ${member?"":`disabled title="Lab 구성원으로 로그인해야 합니다"`}>Copy to Lab</button>`); acts.push(`<button class="btn sm" data-a="del">Delete</button>`); }
  if(r.layer==="lab"){ acts.push(`<button class="btn sm" data-a="edit" ${ed?"":`disabled title="만든 사람이나 관리자만 고칠 수 있습니다"`}>Edit</button>`); acts.push(`<button class="btn sm" data-a="hide" ${ed?"":"disabled"}>${r.hidden?"Unhide":"Hide"}</button>`);
    if(s.isAdmin&&s.configOk) acts.push(`<button class="btn sm" data-a="verify">${r.verified?"Unverify":"Verify"}</button>`); }
  if(r.layer!=="local") acts.push(`<button class="btn sm" data-a="copy" title="고쳐 쓰려고 Local 에 사본을 만듭니다">Copy to Local</button>`);
  acts.push(`<button class="btn sm" data-a="json">Export JSON</button>`);
  const link=(u,t)=>`<a href="${esc(u)}" target="_blank" rel="noopener">${esc(t)}</a>`;
  const ringRows=R.rings?R.rings.map((x,i)=>`<tr><td>${i+1}</td><td class="num">${x.d.toFixed(5)}</td><td class="num">${x.g.toFixed(4)}</td><td class="mono">{${esc(hklStr(x.hkl))}}</td><td class="num">${x.m}</td><td class="num">${x.I==null?"—":x.I<0.05?"&lt;0.1":x.I.toFixed(1)}</td></tr>`).join(""):"";
  const hist=r.layer==="local"?(r.history||[]).map(h=>({version:h.version,updatedAt:h.updatedAt,by:"",reason:"",rec:h.content})):(HIST[r.id]||null);
  el.innerHTML=`<div class="ttl"><h3>${esc(r.name||r.formula||r.id)}</h3><span class="lay ${r.layer}">${LAYER[r.layer]}</span>${statusBadge(r)}</div>
    <div class="acts">${acts.join("")}</div>
    ${vs&&vs.state==="edited"?`<div class="note warn" style="margin-bottom:8px">v${esc(vs.verifiedVersion)} 에서 검증된 뒤 v${esc(vs.version)} 로 편집되었습니다 — 관리자가 다시 확인해야 합니다.</div>`:""}
    ${r.promotedTo?`<div class="note" style="margin-bottom:8px">Lab 에 올렸습니다: <span class="mono">${esc(r.promotedTo)}</span></div>`:""}
    <dl class="kv">
      <dt>Formula</dt><dd>${esc(r.formula||"—")}${els.length?` <span class="muted">(${esc(els.join(" "))})</span>`:""}</dd>
      <dt>Space group</dt><dd>${n?`#${n} ${esc(r.sg.symbol||"")}`:esc(r.sg&&r.sg.symbol||"—")}${r.symops?` <span class="muted">· ${r.symops.length} ops from CIF</span>`:""}</dd>
      <dt>System</dt><dd>${n?esc(DC().crystalSystem(n)):"—"}</dd>
      <dt>a, b, c (nm)</dt><dd class="mono">${fx(c.a,5)}, ${fx(c.b,5)}, ${fx(c.c,5)}</dd>
      <dt>α, β, γ (°)</dt><dd class="mono">${fx(c.al,3)}, ${fx(c.be,3)}, ${fx(c.ga,3)}</dd>
      <dt>Volume (nm³)</dt><dd class="mono">${R.V?R.V.toFixed(5):"—"}</dd>
      <dt>Source</dt><dd>${esc(src.kind||"—")}${src.codId?` · ${link("https://www.crystallography.net/cod/"+encodeURIComponent(src.codId)+".html","COD "+src.codId)}`:""}${src.doi?` · ${link("https://doi.org/"+src.doi,"doi:"+src.doi)}`:""}${src.citation?`<br><span class="muted">${esc(src.citation)}</span>`:""}${src.note?`<br><span class="muted">${esc(src.note)}</span>`:""}</dd>
      ${r.note?`<dt>Note</dt><dd>${esc(r.note)}</dd>`:""}
      <dt>ID</dt><dd class="mono brk">${esc(r.id)}</dd>
      <dt>Version</dt><dd>${r.version!=null?"v"+esc(r.version):"—"}${r.updatedAt?` <span class="muted">· ${esc(when(r.updatedAt))}</span>`:""}${r.layer==="pack"&&D.packInfo().version?` <span class="muted">· pack ${esc(D.packInfo().version)}</span>`:""}</dd>
      ${r.layer==="lab"?`<dt>Created by</dt><dd class="brk">${esc(r.createdByEmail||"—")} <span class="muted">· ${esc(when(r.createdAt))}</span></dd>${r.updatedByEmail&&r.updatedByEmail!==r.createdByEmail?`<dt>Updated by</dt><dd class="brk">${esc(r.updatedByEmail)}</dd>`:""}${r.verified?`<dt>Verified by</dt><dd class="brk">${esc(r.verifiedByEmail||"—")} <span class="muted">· v${esc(r.verifiedVersion)}</span></dd>`:""}`:""}
      ${r.copiedFrom?`<dt>Copied from</dt><dd class="mono brk">${esc(r.copiedFrom)}</dd>`:""}
    </dl>
    <div class="sec"><h4>Rings <span class="muted">g ≤ 10 1/nm · d, hkl, multiplicity, intensity</span></h4>
      ${R.error?`<div class="note danger">링을 계산하지 못했습니다: ${esc(R.error)}</div>`:`${R.hasI?"":`<div class="note" style="margin-bottom:6px">원자 좌표가 없어 세기(I)는 계산하지 않습니다 — 소멸칙은 공간군으로만 적용합니다.</div>`}
      <div class="tblwrap"><table class="tbl" id="ddbRings"><thead><tr><th>#</th><th>d (nm)</th><th>g (1/nm)</th><th>hkl</th><th>m</th><th>I (%)</th></tr></thead><tbody>${ringRows}</tbody></table></div>
      <div class="muted xs" style="margin-top:4px">I = 운동학 |F|² × 다중도 합 (전자 산란 인자, Debye–Waller 포함), 가장 센 링 = 100.</div>`}</div>
    ${r.atoms&&r.atoms.length?`<details class="sec"><summary>Atoms (${r.atoms.length})</summary><div class="tblwrap"><table class="tbl"><thead><tr><th class="l">Label</th><th>El</th><th>x</th><th>y</th><th>z</th><th>Occ</th><th>B (Å²)</th></tr></thead><tbody>
      ${r.atoms.map(a=>`<tr><td class="l">${esc(a.label)}</td><td>${esc(a.el)}</td><td class="num">${fx(a.x,5)}</td><td class="num">${fx(a.y,5)}</td><td class="num">${fx(a.z,5)}</td><td class="num">${fx(a.occ,3)}</td><td class="num">${fx(a.B,3)}</td></tr>`).join("")}</tbody></table></div></details>`:""}
    ${r.layer!=="pack"?`<details class="sec" id="ddbHistSec" ${r.layer==="lab"||(hist&&hist.length)?"open":""}><summary>History</summary>
      ${r.layer==="lab"&&!hist?`<div class="row" style="margin:4px 0"><button class="btn sm" data-a="hist" ${member?"":`disabled title="Lab 구성원으로 로그인해야 합니다"`}>Load history</button><span class="muted">지난 판본은 서버의 history 에 있습니다</span></div>`:""}
      ${hist&&hist.length?`<div class="tblwrap"><table class="tbl" id="ddbHist"><thead><tr><th>Version</th><th>Updated</th><th class="l">By</th><th class="l">Change</th><th>a (nm)</th><th>c (nm)</th></tr></thead><tbody>
        ${hist.map(h=>`<tr><td>v${esc(h.version)}</td><td>${esc(when(h.updatedAt))}</td><td class="l brk">${esc(h.by||h.updatedByEmail||"")}</td><td class="l">${esc(h.reason||h.archivedReason||"")}</td><td class="num">${fx(h.rec&&h.rec.cell&&h.rec.cell.a,5)}</td><td class="num">${fx(h.rec&&h.rec.cell&&h.rec.cell.c,5)}</td></tr>`).join("")}</tbody></table></div>
        <div class="muted xs" style="margin-top:4px">표의 각 줄은 그 판본이 바뀌기 전의 내용입니다. 판정 결과에 저장한 refOf(version) 로 다시 찾을 수 있습니다.</div>`:hist?`<div class="muted">지난 판본이 없습니다 (v1).</div>`:""}</details>`:""}`;
  el.onclick=async e=>{ const b=e.target.closest("button[data-a]"); if(!b||BUSY) return; const a=b.dataset.a; BUSY=true;
    try{ if(a==="edit") openEditor(r);
      else if(a==="json") download(`phase_${String(r.name||r.id).replace(/[^\w.-]+/g,"_")}.json`,D.exportJSON([r.id]));
      else if(a==="copy"){ const x=await D.copyToLocal(r.id); if(x.ok){ U.sel=x.rec.id; say("Local 에 사본을 만들었습니다","ok"); } }
      else if(a==="promote"){ let x=await D.promoteLocal(r.id); if(!x.ok&&x.reason==="duplicate"){ if(confirm(`Lab 에 같은 상이 이미 있습니다: ${x.duplicate.name} (${LAYER[x.duplicate.layer]||x.duplicate.layer}, ${x.duplicate.id}).\n그래도 올릴까요?`)) x=await D.promoteLocal(r.id,{force:true}); else { BUSY=false; return; } }
        if(x.ok){ U.sel=x.rec.id; say("Lab 에 올렸습니다 (v1, 검증 전)","ok"); } else say(x.reason==="invalid"?"검사 오류: "+x.validation.errors.join(" / "):"올리지 못했습니다","err"); }
      else if(a==="del"){ if(confirm(`Local 항목 '${r.name}' 을 지울까요?\n이 상을 쓴 판정 결과는 refOf 로 다시 찾을 수 없게 됩니다.`)){ await D.removeLocal(r.id); U.sel=null; say("지웠습니다"); } }
      else if(a==="hide"){ const x=await D.hideLab(r.id,!r.hidden); say(x.hidden?"숨겼습니다 (Lab 에는 판본으로 남습니다)":"다시 보이게 했습니다","ok"); }
      else if(a==="verify"){ const x=await D.setVerified(r.id,!r.verified); say(x.verified?`검증했습니다 (v${x.version})`:"검증을 풀었습니다","ok"); }
      else if(a==="hist"){ HIST[r.id]=(await D.labHistory(r.id)).map(h=>Object.assign(h,{by:h.updatedByEmail,reason:""})); }
    }catch(x){ say(errMsg(x),"err"); } finally{ BUSY=false; drawDetail(); } }; }

/* ---------------- 추가·편집 대화상자 ---------------- */
const SYS_HINT={Cubic:"a = b = c, α = β = γ = 90°",Tetragonal:"a = b, α = β = γ = 90°",Orthorhombic:"α = β = γ = 90°",Hexagonal:"a = b, α = β = 90°, γ = 120°",Trigonal:"hexagonal axes: a = b, γ = 120° · rhombohedral (:R): a = b = c, α = β = γ",Monoclinic:"unique axis b: α = γ = 90°",Triclinic:"no constraints"};
function openEditor(rec){ const D=DB(); const s=D.firebaseStatus(); const member=!!(s.user&&s.isMember&&s.configOk); const editing=!!(rec&&rec.id);
  const layer0=editing?rec.layer:(member&&U.layers.includes("lab")&&!U.layers.includes("local")?"lab":"local");
  const M={mode:editing?"params":"cif",target:layer0,blocks:[],rec:rec?D.contentOf(rec):{name:"",formula:"",cell:{a:null,b:null,c:null,al:90,be:90,ga:90},sg:{number:null,symbol:""},source:{kind:"manual"}},keepOps:!!(rec&&rec.symops),dupOk:false,id:editing?rec.id:null,layer:editing?rec.layer:null};
  if(!M.rec.atoms) M.rec.atoms=[];
  const ov=document.createElement("div"); ov.className="overlay center"; ov.innerHTML=`<div class="modal ddb-modal" role="dialog" aria-modal="true" aria-labelledby="ddbEdT"><button class="closex" aria-label="Close" title="Close">✕</button><h2 id="ddbEdT">${editing?"Edit Phase":"Add Phase"}</h2><div id="ddbEdBody"></div></div>`;
  document.body.appendChild(ov); const close=()=>{ ov.remove(); document.removeEventListener("keydown",onKey); }; const onKey=e=>{ if(e.key==="Escape") close(); }; document.addEventListener("keydown",onKey);
  $q(ov,".closex").onclick=close; ov.addEventListener("mousedown",e=>{ if(e.target===ov) close(); });
  const body=$q(ov,"#ddbEdBody");
  const targetSeg=()=>editing?`<span class="muted">Layer: <span class="lay ${M.layer}">${LAYER[M.layer]}</span></span>`:`<span class="lbl muted">Save to</span><div class="seg" id="ddbTgt"><button data-t="local" class="${M.target==="local"?"active":""}">Local</button><button data-t="lab" class="${M.target==="lab"?"active":""}" ${member?"":`disabled title="Lab 구성원으로 로그인해야 합니다"`}>Lab</button></div>`;
  function draw(){ body.innerHTML=`${editing?"":`<div class="seg" id="ddbMode"><button data-m="cif" class="${M.mode==="cif"?"active":""}">From CIF</button><button data-m="params" class="${M.mode==="params"?"active":""}">From parameters</button></div>`}
      <div id="ddbPane"></div>
      <div class="foot">${targetSeg()}<span class="sp"></span><button class="btn" id="ddbCancel">Cancel</button><button class="btn primary" id="ddbSave">Save</button></div>`;
    $a(body,"#ddbMode button").forEach(b=>b.onclick=()=>{ M.mode=b.dataset.m; draw(); });
    $a(body,"#ddbTgt button").forEach(b=>b.onclick=()=>{ if(b.disabled) return; M.target=b.dataset.t; $a(body,"#ddbTgt button").forEach(x=>x.classList.toggle("active",x===b)); refresh(); });
    $q(body,"#ddbCancel").onclick=close; $q(body,"#ddbSave").onclick=save;
    if(M.mode==="cif") drawCif(); else drawParams(); }
  /* --- CIF --- */
  function drawCif(){ const p=$q(body,"#ddbPane");
    p.innerHTML=`<div class="row" style="margin-bottom:8px"><label class="btn sm" style="cursor:pointer">Choose CIF files<input type="file" id="ddbCifF" accept=".cif,.txt,chemical/x-cif" multiple hidden></label><span class="muted">또는 아래에 CIF 글을 붙여 넣고 Parse (여러 data_ 블록 가능)</span></div>
      <label class="f">CIF text<textarea id="ddbCifT" spellcheck="false" placeholder="data_..."></textarea></label>
      <div class="row" style="margin:8px 0"><button class="btn sm" id="ddbParse">Parse</button><span class="muted" id="ddbCifN"></span></div>
      <div id="ddbBlocks"></div>`;
    $q(p,"#ddbCifF").onchange=async e=>{ const fs=[...(e.target.files||[])]; e.target.value=""; for(const f of fs){ try{ addBlocks(D.importCIF(await f.text(),{fileName:f.name})); }catch(x){ say(f.name+": "+errMsg(x),"err"); } } drawBlocks(); };
    $q(p,"#ddbParse").onclick=()=>{ const t=$q(p,"#ddbCifT").value; if(!t.trim()){ say("CIF 글이 비어 있습니다","warn"); return; } addBlocks(D.importCIF(t,{fileName:"pasted"})); drawBlocks(); };
    drawBlocks(); }
  function addBlocks(list){ for(const it of list){ it.use=!!(it.rec&&it.validation&&it.validation.ok&&!it.duplicate); M.blocks.push(it); } }
  function drawBlocks(){ const el=$q(body,"#ddbBlocks"); if(!el) return; const n=M.blocks.filter(b=>b.use).length; const cn=$q(body,"#ddbCifN"); if(cn) cn.textContent=M.blocks.length?`${M.blocks.length} block${M.blocks.length>1?"s":""} · ${n} selected`:"";
    el.innerHTML=M.blocks.map((b,i)=>{ if(!b.rec) return `<div class="blk"><span></span><div><div class="nm">${esc(b.file||"CIF")}${b.block?" · data_"+esc(b.block):""}</div></div><span></span><div class="msg"><div class="note danger">${esc(b.error)}</div></div></div>`;
      const r=b.rec, c=r.cell, v=b.validation;
      return `<div class="blk" data-i="${i}"><input type="checkbox" data-u="${i}" ${b.use?"checked":""} ${v.ok?"":"disabled"} aria-label="Use ${esc(r.name)}">
        <div><div class="nm">${esc(r.name)} <span class="muted">${esc(r.formula)}</span></div><div class="sub">#${esc(r.sg.number||"?")} ${esc(r.sg.symbol)} · a ${fx(c.a)} b ${fx(c.b)} c ${fx(c.c)} nm · ${fx(c.al,2)}° ${fx(c.be,2)}° ${fx(c.ga,2)}° · ${(r.atoms||[]).length} atoms${r.source.codId?" · COD "+esc(r.source.codId):""}</div><div class="sub">${esc(b.file)} · data_${esc(b.block)}</div></div>
        <button class="btn sm" data-e="${i}" title="이 블록을 매개변수 편집기로">Edit…</button>
        <div class="msg">${v.errors.length?`<div class="note danger">${v.errors.map(esc).join("<br>")}</div>`:""}${v.warnings.length?`<div class="note warn">${v.warnings.map(esc).join("<br>")}</div>`:""}${b.duplicate?`<div class="note warn">중복 후보: <b>${esc(b.duplicate.name||b.duplicate.formula)}</b> (${esc(LAYER[b.duplicate.layer]||"이 가져오기 안")}${b.duplicate.id?" · "+esc(b.duplicate.id):""}) — 그래도 넣으려면 체크하세요.</div>`:""}</div></div>`; }).join("");
    $a(el,"input[data-u]").forEach(x=>x.onchange=()=>{ M.blocks[+x.dataset.u].use=x.checked; drawBlocks(); });
    $a(el,"button[data-e]").forEach(x=>x.onclick=()=>{ const b=M.blocks[+x.dataset.e]; M.rec=D.contentOf(b.rec); if(!M.rec.atoms) M.rec.atoms=[]; M.keepOps=!!M.rec.symops; M.mode="params"; M.fromBlock=+x.dataset.e; draw(); });
    refresh(); }
  /* --- 매개변수 --- */
  function drawParams(){ const p=$q(body,"#ddbPane"); const r=M.rec, c=r.cell, src=r.source||{};
    const inp=(k,lab,v,o)=>`<label class="f ${o&&o.cls||""}">${lab}<input type="${o&&o.type||"text"}" ${o&&o.type==="number"?`step="any" inputmode="decimal"`:""} data-k="${k}" value="${esc(v==null?"":v)}" ${o&&o.ph?`placeholder="${esc(o.ph)}"`:""}></label>`;
    p.innerHTML=`<div class="fgrid">
        ${inp("name","Name",r.name,{cls:"w2",ph:"e.g. β-NaYF4"})}${inp("formula","Formula",r.formula,{cls:"w2",ph:"e.g. Na1.5 Y1.5 F6"})}
        ${inp("sgn","SG number",r.sg.number,{type:"number",ph:"1–230"})}${inp("sgs","SG symbol",r.sg.symbol,{ph:"e.g. Fm-3m"})}
        <label class="f">System<input type="text" id="ddbSys" disabled value=""></label>
        ${inp("a","a (nm)",c.a,{type:"number"})}${inp("b","b (nm)",c.b,{type:"number"})}${inp("c","c (nm)",c.c,{type:"number"})}
        ${inp("al","α (°)",c.al,{type:"number"})}${inp("be","β (°)",c.be,{type:"number"})}${inp("ga","γ (°)",c.ga,{type:"number"})}
      </div><div class="hint" id="ddbSysHint"></div>
      ${r.symops?`<label class="chk" style="margin-bottom:8px"><input type="checkbox" id="ddbKeep" ${M.keepOps?"checked":""}>Keep CIF symmetry operators (${r.symops.length})</label>`:""}
      <details class="sec" ${r.atoms.length?"open":""}><summary>Atoms (${r.atoms.length}) <span class="muted">— 없으면 세기 없이 소멸칙만</span></summary>
        <div class="tblwrap"><table class="tbl atoms"><thead><tr><th>Label</th><th>El</th><th>x</th><th>y</th><th>z</th><th>Occ</th><th>B (Å²)</th><th></th></tr></thead><tbody>
        ${r.atoms.map((a,i)=>`<tr>${["label","el","x","y","z","occ","B"].map(k=>`<td><input data-at="${i}" data-f="${k}" value="${esc(a[k]==null?"":a[k])}" ${k==="label"||k==="el"?"":`inputmode="decimal"`} aria-label="${k}"></td>`).join("")}<td><button class="btn sm" data-rm="${i}" aria-label="Remove atom" title="Remove atom">×</button></td></tr>`).join("")}
        </tbody></table></div><div class="row" style="margin-top:6px"><button class="btn sm" id="ddbAddAt">Add atom</button></div></details>
      <div class="fgrid" style="margin-top:10px">
        <label class="f">Source<select data-k="kind">${KINDS.map(k=>`<option ${src.kind===k?"selected":""}>${k}</option>`).join("")}</select></label>
        ${inp("codId","COD id",src.codId,{ph:"e.g. 9008565"})}${inp("doi","DOI",src.doi,{cls:"w2",ph:"10.xxxx/…"})}
        ${inp("citation","Citation",src.citation,{cls:"wall"})}${inp("snote","Source note",src.note,{cls:"wall"})}${inp("note","Note",r.note,{cls:"wall"})}
      </div>`;
    p.oninput=e=>{ const t=e.target; if(t.dataset.at!=null){ const a=r.atoms[+t.dataset.at]; a[t.dataset.f]=t.value; refresh(); return; } const k=t.dataset.k; if(!k) return; readField(k,t.value); if(k==="sgn"||k==="sgs"||["a","al"].includes(k)) constrain(k); refresh(); };
    p.onchange=e=>{ const t=e.target; if(t.id==="ddbKeep"){ M.keepOps=t.checked; refresh(); } else if(t.dataset.k==="kind"){ r.source.kind=t.value; refresh(); } };
    p.onclick=e=>{ const rm=e.target.closest("button[data-rm]"); if(rm){ r.atoms.splice(+rm.dataset.rm,1); drawParams(); return; } if(e.target.id==="ddbAddAt"){ r.atoms.push({label:"",el:"",x:0,y:0,z:0,occ:1,B:0}); drawParams(); const ins=$a(p,"input[data-f=label]"); if(ins.length) ins[ins.length-1].focus(); } };
    constrain(null); refresh(); }
  function readField(k,v){ const r=M.rec; const n=v===""?null:+v;
    if(k==="name"||k==="formula") r[k]=v; else if(k==="note") r.note=v;
    else if(k==="sgn"){ r.sg.number=n; delete r.sg.hall; if(n>=1&&n<=230){ try{ r.sg.symbol=DC().sgByNumber(n).short; const si=$q(body,"input[data-k=sgs]"); if(si) si.value=r.sg.symbol; }catch(e){} } if(r.symops&&M.keepOps){ M.keepOps=false; const kp=$q(body,"#ddbKeep"); if(kp) kp.checked=false; } }
    else if(k==="sgs"){ r.sg.symbol=v; delete r.sg.hall; let e=null; try{ e=DC().sgBySymbol(v); }catch(x){} if(e){ r.sg.number=e.number; const ni=$q(body,"input[data-k=sgn]"); if(ni) ni.value=e.number; } }
    else if(["a","b","c","al","be","ga"].includes(k)) r.cell[k]=n;
    else if(["codId","doi","citation"].includes(k)) r.source[k]=v; else if(k==="snote") r.source.note=v; }
  /* 결정계 제약: 묶인 칸은 잠그고 값을 맞춘다 */
  function constrain(){ const r=M.rec, c=r.cell; const n=+r.sg.number; const sys=n>=1&&n<=230?DC().crystalSystem(n):""; const sy=$q(body,"#ddbSys"); if(sy) sy.value=sys||"—";
    const rh=sys==="Trigonal"&&/:\s*r$/i.test(String(r.sg.symbol||"")); const lock={b:0,c:0,al:0,be:0,ga:0};
    const set=(k,v)=>{ c[k]=v; lock[k]=1; };
    if(sys==="Cubic"){ set("b",c.a); set("c",c.a); set("al",90); set("be",90); set("ga",90); }
    else if(sys==="Tetragonal"){ set("b",c.a); set("al",90); set("be",90); set("ga",90); }
    else if(sys==="Orthorhombic"){ set("al",90); set("be",90); set("ga",90); }
    else if(sys==="Hexagonal"||(sys==="Trigonal"&&!rh)){ set("b",c.a); set("al",90); set("be",90); set("ga",120); }
    else if(rh){ set("b",c.a); set("c",c.a); set("be",c.al); set("ga",c.al); }
    else if(sys==="Monoclinic"){ set("al",90); set("ga",90); }
    for(const k of ["b","c","al","be","ga"]){ const i=$q(body,`input[data-k=${k}]`); if(!i) continue; i.disabled=!!lock[k]; if(lock[k]) i.value=c[k]==null?"":c[k]; }
    const h=$q(body,"#ddbSysHint"); if(h) h.textContent=sys?`${sys}: ${SYS_HINT[rh?"Trigonal":sys]}`:"공간군 번호나 기호를 넣으면 결정계에 맞게 격자 칸이 묶입니다."; }
  function buildRec(){ const r=D.contentOf(Object.assign({},M.rec,{symops:M.keepOps?M.rec.symops:undefined,atoms:(M.rec.atoms||[]).filter(a=>String(a.el||a.label||"").trim()).map(a=>Object.assign({},a,{el:String(a.el||"").trim()||(DC().elementOf(a.label)||"")}))})); if(editing) r.id=M.id; return r; }
  let tDup=null;
  function refresh(){ const sv=$q(body,"#ddbSave"); if(!sv) return; let v=$q(body,"#ddbVal");
    if(M.mode==="cif"){ if(v) v.remove(); const n=M.blocks.filter(b=>b.use).length; sv.disabled=!n; sv.textContent=n?`Save ${n} phase${n>1?"s":""}`:"Save"; return; }
    if(!v){ v=document.createElement("div"); v.id="ddbVal"; v.className="val"; $q(body,"#ddbPane").appendChild(v); }
    const rec=buildRec(); const val=DC().validateRecord(rec);
    clearTimeout(tDup); M.dup=null; tDup=setTimeout(()=>{ M.dup=val.ok?D.findDuplicate(rec,{exclude:[M.id]}):null; paint(); },120);
    const paint=()=>{ const d=M.dup; v.innerHTML=`${val.errors.length?`<div class="note danger" id="ddbErr"><b>${val.errors.length} error${val.errors.length>1?"s":""}</b> — 저장할 수 없습니다<ul>${val.errors.map(x=>`<li>${esc(x)}</li>`).join("")}</ul></div>`:""}
        ${val.warnings.length?`<div class="note warn" id="ddbWarn"><b>Check</b><ul>${val.warnings.map(x=>`<li>${esc(x)}</li>`).join("")}</ul></div>`:""}
        ${d?`<div class="note warn" id="ddbDup">같은 상이 이미 있습니다: <b>${esc(d.name||d.formula)}</b> (${esc(LAYER[d.layer]||d.layer)} · <span class="mono">${esc(d.id)}</span>, #${esc(d.sg&&d.sg.number)} a = ${fx(d.cell&&d.cell.a)} nm)
          ${editing?`<br><span class="muted">이미 있는 항목을 고치는 중이라 그대로 저장할 수 있습니다.</span>`:`<label class="chk"><input type="checkbox" id="ddbDupOk" ${M.dupOk?"checked":""}>Save anyway</label>`}</div>`:""}
        ${!val.errors.length&&!val.warnings.length&&!d?`<div class="note okay">검사 통과</div>`:""}`;
      const dk=$q(v,"#ddbDupOk"); if(dk) dk.onchange=()=>{ M.dupOk=dk.checked; paint(); };
      sv.disabled=!val.ok||(!!d&&!M.dupOk&&!editing); sv.textContent=d&&M.dupOk&&!editing?"Save anyway":"Save"; };   // 편집은 중복이어도 막지 않음 (안내만)
    paint(); }
  async function save(){ const sv=$q(body,"#ddbSave"); if(sv.disabled) return; sv.disabled=true;
    try{
      if(M.mode==="cif"){ const use=M.blocks.filter(b=>b.use&&b.rec); let ok=0, dup=0; const errs=[]; let last=null;
        for(const b of use){ const force=!!b.duplicate; const res=M.target==="lab"?await D.addLab(b.rec,{force}):await D.saveLocal(b.rec,{force}); if(res.ok){ ok++; last=res.rec.id; b.use=false; b.saved=true; } else if(res.reason==="duplicate"){ dup++; b.duplicate=res.duplicate; b.use=false; } else errs.push(b.rec.name+": "+res.validation.errors.join(" / ")); }
        M.blocks=M.blocks.filter(b=>!b.saved); if(last) U.sel=last;
        say(`${LAYER[M.target]} 에 ${ok} 개 저장${dup?` · 중복 ${dup} 개는 건너뜀 (체크하면 그래도 저장)`:""}${errs.length?` · 오류 ${errs.length}`:""}`,errs.length||dup?"warn":"ok");
        if(!M.blocks.length&&!errs.length) close(); else drawBlocks(); return; }
      const rec=buildRec(); const force=!!M.dupOk; let res;
      if(editing&&M.layer==="lab") res=await D.updateLab(M.id,rec);
      else if(editing) res=await D.saveLocal(rec,{force:true});
      else res=M.target==="lab"?await D.addLab(rec,{force}):await D.saveLocal(rec,{force});
      if(res.ok){ U.sel=res.rec.id; if(M.fromBlock!=null&&M.blocks[M.fromBlock]) M.blocks.splice(M.fromBlock,1); say(editing?(res.unchanged?"바뀐 것이 없습니다":`저장했습니다 (v${res.rec.version})`):`${LAYER[res.rec.layer]} 에 추가했습니다`,"ok"); close(); drawList(); drawDetail(); return; }
      if(res.reason==="duplicate"){ M.dup=res.duplicate; M.dupOk=false; refresh(); say("같은 상이 이미 있습니다 — 확인 후 Save anyway","warn"); }
      else { refresh(); say("검사 오류: "+res.validation.errors.join(" / "),"err"); }
    }catch(x){ say(errMsg(x),"err"); } finally{ if(sv.isConnected) refresh(); } }
  draw(); return ov; }

/* ---------------- 보정 ---------------- */
function drawCal(){ const el=$q(ROOT,"#ddbCal"); if(!el) return; const D=DB(); const L=D.listCalibrations();
  el.innerHTML=`<details class="sec" style="margin:0" ${L.length&&L.length<=8?"open":""}><summary>Calibrations (${L.length}) <span class="muted">— 카메라 길이 보정 (Diffraction 화면에서 만듦), 읽기 전용</span></summary>
    ${L.length?`<div class="tblwrap"><table class="tbl" id="ddbCalT"><thead><tr><th class="l">Name</th><th class="l">Microscope</th><th>kV</th><th>Camera length (mm)</th><th>Mode</th><th>Scale (1/nm/px)</th><th>Ellip.</th><th>φ (°)</th><th>σ rel</th><th class="l">Standard</th><th>Layer</th><th class="l">By</th><th>Updated</th></tr></thead><tbody>
      ${L.map(c=>{ const st=c.standard?D.get(c.standard):null; return `<tr><td class="l">${esc(c.name)}</td><td class="l">${esc(c.microscope||"—")}</td><td class="num">${esc(c.kV??"—")}</td><td class="num">${esc(c.cameraLength??"—")}</td><td>${esc(c.mode||"—")}</td><td class="num">${c.scale!=null?(+c.scale).toPrecision(6):"—"}</td><td class="num">${c.ellip!=null?(+c.ellip).toFixed(4):"—"}</td><td class="num">${c.phi!=null?(+c.phi).toFixed(1):"—"}</td><td class="num">${c.sigmaRel!=null?(+c.sigmaRel).toExponential(1):"—"}</td><td class="l">${esc(st?st.name:c.standardRef?c.standardRef.name:c.standard||"—")}</td><td><span class="lay ${c.layer==="lab"?"lab":"local"}">${c.layer==="lab"?"Lab":"Local"}</span></td><td class="l brk">${esc(c.createdByEmail||"")}</td><td>${esc(when(c.updatedAt))}</td></tr>`; }).join("")}
      </tbody></table></div>`:`<div class="muted" style="padding:6px 0">저장된 보정이 없습니다. Lab 보정은 로그인 후 Sync 로 받습니다.</div>`}</details>`; }

/* ---------------- 구성원 (관리자) ---------------- */
function drawMembers(){ const el=$q(ROOT,"#ddbMem"); if(!el) return; const D=DB(); const s=D.firebaseStatus(); if(!(s.isAdmin&&s.user)){ el.hidden=true; el.innerHTML=""; return; } el.hidden=false;
  const admins=s.admins||[];
  el.innerHTML=`<details class="sec" style="margin:0" id="ddbMemD" ${MEMBERS?"open":""}><summary>Members <span class="muted">— 관리자만 보입니다</span></summary>
    <div class="row" style="margin:6px 0"><button class="btn sm" data-m="load">${MEMBERS?"Reload":"Load members"}</button><span class="muted">구성원(members) 문서가 active 일 때만 Lab DB 를 읽고 쓸 수 있습니다.</span></div>
    ${MEMBERS?`<div class="tblwrap"><table class="tbl" id="ddbMemT"><thead><tr><th class="l">Email</th><th class="l">Name</th><th>Status</th><th class="l">Added by</th><th></th></tr></thead><tbody>
      ${MEMBERS.length?MEMBERS.map(m=>`<tr><td class="l brk">${esc(m.email)}</td><td class="l">${esc(m.name||"")}</td><td>${m.active?`<span class="badge ok">Active</span>`:`<span class="badge gray">Inactive</span>`}</td><td class="l brk">${esc(m.addedBy||"")}</td><td><button class="btn sm" data-act="${esc(m.email)}" data-on="${m.active?0:1}">${m.active?"Deactivate":"Activate"}</button></td></tr>`).join(""):`<tr><td colspan="5" class="muted">구성원이 없습니다</td></tr>`}
      </tbody></table></div>`:""}
    <div class="fgrid" style="display:grid;grid-template-columns:repeat(auto-fill,minmax(160px,1fr));gap:8px;margin-top:8px;align-items:end">
      <label class="f">Email<input type="text" id="ddbMemE" inputmode="email" autocomplete="off" placeholder="ALLOWED_EMAIL_4@example.com"></label><label class="f">Name<input type="text" id="ddbMemN"></label><div><button class="btn sm primary" data-m="add">Add member</button></div></div>
    <div class="sec"><h4>Admins <span class="muted">config/admins</span></h4>
      <div class="row">${admins.map(a=>`<span class="chip">${esc(a)}<button data-rmadm="${esc(a)}" ${admins.length<=1?`disabled title="마지막 관리자는 지울 수 없습니다"`:`title="Remove admin"`} aria-label="Remove admin ${esc(a)}">×</button></span>`).join("")}</div>
      <div class="row" style="margin-top:6px"><input type="text" id="ddbAdmE" inputmode="email" placeholder="admin email" aria-label="Admin email" style="max-width:240px"><button class="btn sm" data-m="addadm">Add admin</button></div></div></details>`;
  el.onclick=async e=>{ const b=e.target.closest("button"); if(!b||b.disabled||BUSY) return; BUSY=true;
    try{ if(b.dataset.m==="load"){ MEMBERS=await D.listMembers(); }
      else if(b.dataset.m==="add"){ const em=$q(el,"#ddbMemE").value, nm=$q(el,"#ddbMemN").value; await D.addMember(em,nm); say("구성원을 추가했습니다: "+em.trim().toLowerCase(),"ok"); MEMBERS=await D.listMembers(); }
      else if(b.dataset.act){ await D.setMemberActive(b.dataset.act,b.dataset.on==="1"); MEMBERS=await D.listMembers(); }
      else if(b.dataset.m==="addadm"){ await D.addAdmin($q(el,"#ddbAdmE").value); say("관리자를 추가했습니다","ok"); }
      else if(b.dataset.rmadm){ if(confirm(`관리자에서 빼기: ${b.dataset.rmadm}?`)){ await D.removeAdmin(b.dataset.rmadm); say("관리자에서 뺐습니다","ok"); } }
    }catch(x){ say(errMsg(x),"err"); } finally{ BUSY=false; drawMembers(); } }; }

G.DIFF_DB_UI={render,openEditor:r=>openEditor(r||null),_state:()=>U,rerender:()=>{ if(ROOT) drawAll(); }};
})();
