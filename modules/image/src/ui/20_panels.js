/* ============================================================
   오른쪽 패널: Adjust · Process · Threshold · Binary · Analyze · Measure · Macro · Info
   ============================================================ */
const PANELS=[["adjust","Adjust"],["process","Process"],["threshold","Threshold"],["binary","Binary"],["analyze","Analyze"],["measure","Measure"],["macro","Macro"],["info","Info"]];

/* ---------- 매개변수 입력 칸 (등록부에서 만든다) ---------- */
const FORM={}; // op id → 마지막 값
function formVals(o){ const v=FORM[o.id]||(FORM[o.id]={}); o.params.forEach(p=>{ if(v[p.k]===undefined) v[p.k]=p.type==="sel"&&Array.isArray(p.def)?p.def[0]:p.def; }); return v; }
function formHTML(o,im){ const v=formVals(o);
  return o.params.filter(p=>!(p.stackOnly&&!isStack(im))&&(!p.showIf||p.showIf(v))).map(p=>{ const id=`imf_${o.id}_${p.k}`; const u=p.unit?`<span class="u">${esc(p.unit)}</span>`:"";
    if(p.type==="bool") return `<label class="im-chk"><input type="checkbox" id="${id}" data-k="${p.k}" ${v[p.k]?"checked":""}> ${esc(p.label)}</label>`;
    if(p.type==="sel") return `<label class="im-f"><span>${esc(p.label)}</span><select id="${id}" data-k="${p.k}" class="im-sel">${p.opts.map(([a,b])=>`<option value="${esc(a)}" ${String(v[p.k])===String(a)?"selected":""}>${esc(b)}</option>`).join("")}</select></label>`;
    if(p.type==="img"){ const L=S.images.filter(x=>x.id!==im.id); if(L.length&&!L.some(x=>x.id===v[p.k])) v[p.k]=L[0].id; /* 보이는 첫 항목 = 값 (고르지 않고 실행해도 화면과 같게) */
      return `<label class="im-f"><span>${esc(p.label)}</span><select id="${id}" data-k="${p.k}" class="im-sel">${L.map(x=>`<option value="${x.id}" ${v[p.k]===x.id?"selected":""}>${esc(x.name)} (${x.w}×${x.h})</option>`).join("")}</select></label>`; }
    if(p.type==="text") return `<label class="im-f col"><span>${esc(p.label)}</span><textarea id="${id}" data-k="${p.k}" rows="3" class="mono">${esc(v[p.k])}</textarea></label>`;
    return `<label class="im-f"><span>${esc(p.label)}</span><input type="number" id="${id}" data-k="${p.k}" value="${esc(v[p.k])}" ${p.min!=null?`min="${p.min}"`:""} ${p.max!=null?`max="${p.max}"`:""} step="${p.step||(p.type==="int"?1:"any")}">${u}</label>`; }).join("");
}
function bindForm(el,o,onChange){ $a(el,"[data-k]").forEach(inp=>{ const p=o.params.find(q=>q.k===inp.dataset.k); const ev=inp.tagName==="SELECT"||inp.type==="checkbox"?"change":"input";
  inp.addEventListener(ev,()=>{ const v=formVals(o); v[p.k]=p.type==="bool"?inp.checked:(p.type==="num"||p.type==="int")?(inp.value===""?"":+inp.value):inp.value; if(p.type==="sel"||p.type==="bool") onChange(true); else onChange(false); }); }); }
async function runOp(id,params,im){
  const o=OPS[id]; im=im||cur(); if(!im){ note("이미지를 먼저 여세요","warn"); return; } if(BUSY){ note("계산 중입니다 — 끝나면 다시 누르세요","warn"); return; }
  const P=Object.assign({},params||formVals(o),{__op:id}); _pvSeq++; clearTimeout(_pvT); clearPreview(); const pvWas=S.ui.preview; BUSY=id;
  try{ if(pvWas&&o.kind==="plane") S.ui.preview=false; const r=await runWithUnit(o,im,P); BUSY=null; if(r&&r.id) S.sel=r.id; save(); rerender(); if(pvWas&&o.kind==="plane") note("적용했습니다 — Preview 는 꺼 두었습니다 (적용한 결과에 한 번 더 미리 보이지 않도록)"); return r; }
  catch(e){ BUSY=null; hideProg(); if(String(e.message)!=="취소") note(o.label+": "+(e.message||e),"warn"); rerender(); }
}

/* ---------- 패널 그리기 ---------- */
function drawTabs(){ const el=$q(ROOT,"#imTabs"); if(!el) return; el.innerHTML=PANELS.map(([k,l])=>`<button data-p="${k}" class="${S.ui.panel===k?"active":""}">${l}</button>`).join("");
  $a(el,"button").forEach(b=>b.onclick=()=>{ S.ui.panel=b.dataset.p; clearPreview(); save(); drawTabs(); renderPanel(); draw(); }); }
function renderPanel(){ const el=$q(ROOT,"#imPanel"); if(!el) return; if(S.ui.panel==="drop") S.ui.panel="analyze"; /* 예전 Drop 탭 → Analyze 안 */ const im=cur(); if(S.ui.panel!=="process"&&PREVIEW.plane){ PREVIEW.plane=false; PREVIEW.planes=null; PREVIEW.key=""; DISP.key=null; }
  if(!im&&S.ui.panel!=="macro"){ el.innerHTML=`<div class="im-empty sm">이미지를 열면 여기에서 조정·처리·분석합니다.</div>`; return; }
  ({adjust:panelAdjust,process:panelProcess,threshold:panelThreshold,binary:panelBinary,analyze:panelAnalyze,measure:panelMeasure,macro:panelMacro,info:panelInfo}[S.ui.panel]||panelAdjust)(el,im); }

/* ---------- 히스토그램 SVG ---------- */
function histSVG(h,o){ o=o||{}; const W=o.w||380, H=o.h||110, n=h.h.length; let mx=0; for(let i=0;i<n;i++) mx=Math.max(mx,o.log?Math.log1p(h.h[i]):h.h[i]); if(!mx) mx=1;
  const bw=W/n; let path=""; for(let i=0;i<n;i++){ const v=(o.log?Math.log1p(h.h[i]):h.h[i])/mx*(H-4); path+=`M${(i*bw).toFixed(1)},${H}V${(H-v).toFixed(1)}`; }
  const xOf=val=>{ const span=h.max-h.min||1; return clamp((val-h.min)/span*W,0,W); };
  let marks=""; if(o.lo!=null){ const a=xOf(o.lo), b=xOf(o.hi); marks+=`<rect x="${Math.min(a,b)}" y="0" width="${Math.abs(b-a)}" height="${H}" class="sel"/><line x1="${a}" x2="${a}" y1="0" y2="${H}" class="mk"/><line x1="${b}" x2="${b}" y1="0" y2="${H}" class="mk"/>`; }
  if(o.dmin!=null){ const a=xOf(o.dmin), b=xOf(o.dmax); marks+=`<polyline points="${a},${H} ${b},0" class="ramp"/><line x1="${a}" x2="${a}" y1="0" y2="${H}" class="dm"/><line x1="${b}" x2="${b}" y1="0" y2="${H}" class="dm"/>`; }
  return `<svg class="im-hist" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none" width="100%" height="${H}" data-min="${h.min}" data-max="${h.max}"><path d="${path}" class="bars" stroke-width="${Math.max(1,bw*0.9)}"/>${marks}</svg>`; }
function planeHist(im,c){ const p=planeAt(im,c); if(p.bit===8&&!im.rgb) return IC.histogram(p,null); return IC.histogram(p,null,256); }

/* ===== Adjust ===== */
// 히스토그램·통계는 (이미지, 판, 채널, z, t) 마다 한 번만 — 패널을 다시 그릴 때마다 큰 이미지 전체를 다시 세지 않게
const ADJ={key:"",v:null};
function adjStats(im,c){ const k=[im.id,im.ver,c,im.z,im.t].join("|"); if(ADJ.key===k) return ADJ.v; const p=im.rgb?curPlaneOf(im):planeAt(im,c);
  const v={h:im.rgb?IC.histogram(p,null):planeHist(im,c),st:IC.stats(p,null,null,null,false,{noMedian:true})}; ADJ.key=k; ADJ.v=v; return v; }
function panelAdjust(el,im){
  const c=im.rgb?0:im.c, ds=im.disp[c]; const {h,st}=adjStats(im,c);
  const lutOpts=LUTS.map(n=>`<button class="lut ${ds.lut===n?"on":""}" data-l="${n}" title="${n}" style="background:${lutCss(n)}"></button>`).join("");
  const cal=im.cal; const cs=calOf(im);
  el.innerHTML=`
  <div class="im-h">Brightness / Contrast ${im.nc>1&&!im.rgb?`<span class="sp"></span><select class="im-sel" id="imCh">${[...Array(im.nc).keys()].map(k=>`<option value="${k}" ${k===im.c?"selected":""}>${im.chNames&&im.chNames[k]?esc(im.chNames[k]):"Channel "+(k+1)}</option>`).join("")}</select><label class="im-chk"><input type="checkbox" id="imComp" ${im.composite?"checked":""}> Composite</label>`:""}</div>
  ${histSVG(h,{dmin:im.rgb?0:ds.min,dmax:im.rgb?255:ds.max,log:S.ui.histLog})}
  <div class="im-row xs muted keep"><span>min ${g4(st.min)}</span><span>max ${g4(st.max)}</span><span>mean ${g4(st.mean)}</span><span>SD ${g4(st.stdDev)}</span><span class="sp"></span><label class="im-chk"><input type="checkbox" id="imHLog" ${S.ui.histLog?"checked":""}> log</label></div>
  ${im.rgb?`<div class="note sm">RGB 이미지는 표시 범위를 바꾸지 않습니다. 세기 분석은 Type ▸ RGB → 8-bit 또는 Split channels 후에.</div>`:`
  <div class="im-slide"><span>Minimum</span><input type="range" id="imDmin" min="${h.min}" max="${h.max}" step="${(h.max-h.min)/1000||1}" value="${ds.min}"><input type="number" id="imDminN" value="${+ds.min.toPrecision(6)}" step="any"></div>
  <div class="im-slide"><span>Maximum</span><input type="range" id="imDmax" min="${h.min}" max="${h.max}" step="${(h.max-h.min)/1000||1}" value="${ds.max}"><input type="number" id="imDmaxN" value="${+ds.max.toPrecision(6)}" step="any"></div>
  <div class="im-row"><button class="btn sm" id="imAuto" title="Auto: 0.35 % saturated (as in ImageJ)">Auto</button><button class="btn sm" id="imReset">Reset</button>${im.bit!==32?`<button class="btn sm" id="imApplyLut" title="Apply the display range to the pixel values (undoable)">Apply LUT</button>`:""}<span class="muted xs">표시만 바꿉니다 — 측정값은 그대로</span></div>
  <div class="im-h">LUT</div><div class="im-luts">${lutOpts}</div>`}
  <div class="im-h">Type <span class="muted xs">${im.rgb?"RGB":im.bit+"-bit"}${im.nc>1&&!im.rgb?` · ${im.nc} channels`:""}${isStack(im)?` · ${im.nz}×${im.nt} slices`:""}</span></div>
  <div class="im-row">${im.rgb?`<button class="btn sm" data-type="rgb2gray">RGB → 8-bit</button><button class="btn sm" data-type="split">Split channels</button>`:`${["8","16","32"].filter(b=>+b!==im.bit).map(b=>`<button class="btn sm" data-type="${b}">${b}-bit</button>`).join("")}${im.nc>1?`<button class="btn sm" data-type="split">Split channels</button>`:""}`}</div>
  <div class="im-h">Set Scale <span class="muted xs">${cal?`${g4(cs.pw)} ${esc(cs.unit)}/px`:"uncalibrated (px)"}</span></div>
  ${cal?`<div class="xs muted keep">Source: ${esc(cal.src||"user-defined")}</div>`:`<div class="note warn sm">화소 크기가 없습니다. 파일 메타데이터에 없으면 <b>Set Scale 도구</b>(도구 막대의 자 모양)로 축척 막대 위에 선을 그어 지정하세요.</div>`}
  ${(im.calCandidates||[]).map((k,i)=>`<div class="note sm">Candidate: ${g4(k.pw)} ${esc(k.unit)}/px — ${esc(k.src)} <button class="btn xs" data-cand="${i}">Apply</button><br><span class="xs muted">TIFF/PNG 해상도 태그는 기본값(72 dpi 등)인 경우가 많아 자동으로 쓰지 않았습니다.</span></div>`).join("")}
  <div class="im-row"><label class="im-f"><span>Distance in pixels</span><input type="number" id="imScPx" step="any" value="${esc(SCALE.px||"")}"></label><label class="im-f"><span>Known distance</span><input type="number" id="imScD" step="any" value="${esc(SCALE.d||"")}"></label>
    <label class="im-f"><span>Unit</span><select class="im-sel" id="imScU">${["nm","µm","mm","Å","1/nm","px"].map(u=>`<option ${(SCALE.u||cs.unit)===u?"selected":""}>${u}</option>`).join("")}</select></label></div>
  <div class="im-row"><button class="btn sm primary" id="imScSet">Set</button><label class="im-chk"><input type="checkbox" id="imScAll"> All open images of the same size</label>${cal?`<button class="btn sm" id="imScClr">Remove scale</button>`:""}</div>
  <div class="im-row"><label class="im-chk"><input type="checkbox" id="imSb" ${S.ui.scalebar.on?"checked":""}> Scale bar</label><select class="im-sel" id="imSbPos">${[["br","Lower right"],["bl","Lower left"],["tr","Upper right"],["tl","Upper left"]].map(([a,b])=>`<option value="${a}" ${S.ui.scalebar.pos===a?"selected":""}>${b}</option>`).join("")}</select><label class="im-f"><span>Length</span><input type="number" id="imSbLen" step="any" value="${S.ui.scalebar.len||""}" placeholder="auto"></label></div>
  ${im.databar?`<div class="note warn sm">이미지 아래 ${im.databar.rows} 줄은 장비의 데이터 막대입니다 (${esc(im.databar.src)}). 분석 전에 잘라내세요. <button class="btn xs" id="imBar">Crop</button></div>`:""}`;
  const q=s=>$q(el,s);
  if(!im.rgb){ const setD=(k,v,now)=>{ ds[k]=+v; if(ds.min>ds.max){ if(k==="min") ds.max=ds.min; else ds.min=ds.max; } invalidate(); if(now) draw(); else drawSoon(); save(); }; // min=max 는 ImageJ 처럼 이진 표시
    q("#imDmin").oninput=e=>{ setD("min",e.target.value); q("#imDminN").value=+(+e.target.value).toPrecision(6); };
    q("#imDmax").oninput=e=>{ setD("max",e.target.value); q("#imDmaxN").value=+(+e.target.value).toPrecision(6); };
    q("#imDminN").onchange=e=>{ setD("min",e.target.value,true); renderPanel(); drawSide(); }; q("#imDmaxN").onchange=e=>{ setD("max",e.target.value,true); renderPanel(); drawSide(); };
    q("#imDmin").onchange=q("#imDmax").onchange=()=>{ renderPanel(); drawSide(); };
    q("#imAuto").onclick=()=>{ autoDisplay(im); invalidate(); draw(); renderPanel(); save(); }; q("#imReset").onclick=()=>{ resetDisplay(im); invalidate(); draw(); renderPanel(); save(); };
    if(q("#imApplyLut")) q("#imApplyLut").onclick=()=>runOp("applyLut",{stack:isStack(im)});
    $a(el,".lut").forEach(b=>b.onclick=()=>{ ds.lut=b.dataset.l; invalidate(); draw(); renderPanel(); save(); }); }
  if(q("#imCh")) q("#imCh").onchange=e=>{ im.c=+e.target.value; invalidate(); rerender(); };
  if(q("#imComp")) q("#imComp").onchange=e=>{ im.composite=e.target.checked; invalidate(); draw(); save(); };
  q("#imHLog").onchange=e=>{ S.ui.histLog=e.target.checked; renderPanel(); };
  $a(el,"[data-type]").forEach(b=>b.onclick=()=>runOp("type",Object.assign({},formVals(OPS.type),{to:b.dataset.type})));
  $a(el,"[data-cand]").forEach(b=>b.onclick=()=>{ const k=im.calCandidates[+b.dataset.cand]; setCal(im,{pw:k.pw,ph:k.pw,unit:k.unit,src:k.src+" (사용자가 적용)"}); });
  q("#imScSet").onclick=()=>{ const px=+q("#imScPx").value, d=+q("#imScD").value, u=q("#imScU").value; if(!(px>0&&d>0)){ note("Distance in pixels 와 Known distance 를 넣으세요","warn"); return; }
    const cal={pw:d/px,ph:d/px,unit:u,src:`Set Scale: ${g4(px)} px = ${g4(d)} ${u}`}; if(u==="px"){ removeCal(im); return; }
    const L=q("#imScAll").checked?S.images.filter(x=>x.w===im.w&&x.h===im.h):[im]; L.forEach(x=>setCal(x,cal,true)); record({op:"setScale",params:{px,d,u},label:`Set Scale ${g4(px)} px = ${g4(d)} ${u}`},im); note(`보정: ${g4(cal.pw)} ${u}/px (${L.length}개 이미지)`,"ok"); rerender(); };
  if(q("#imScClr")) q("#imScClr").onclick=()=>removeCal(im);
  q("#imSb").onchange=e=>{ S.ui.scalebar.on=e.target.checked; save(); draw(); }; q("#imSbPos").onchange=e=>{ S.ui.scalebar.pos=e.target.value; save(); draw(); }; q("#imSbLen").onchange=e=>{ S.ui.scalebar.len=+e.target.value||0; save(); draw(); };
  if(q("#imBar")) q("#imBar").onclick=()=>{ runOp("cropbar",{rows:im.databar.rows}).then(()=>{ im.databar=null; rerender(); }); };
}
const SCALE={};
function openSetScale(im,L){ SCALE.px=+L.toFixed(3); S.ui.panel="adjust"; drawTabs(); renderPanel(); const e=$q(ROOT,"#imScD"); if(e){ e.focus(); } note(`선 길이 ${L.toFixed(2)} px — Known distance 와 단위를 넣고 Set`); }
function setCal(im,cal,quiet){ im.cal=cal; save(); if(!quiet){ note(`보정: ${g4(cal.pw)} ${cal.unit}/px`,"ok"); rerender(); } }
function removeCal(im){ im.cal=null; save(); rerender(); note("보정을 지웠습니다 (px)"); }
OPS.setScale={id:"setScale",label:"Set Scale",hidden:true,params:[],run:(im,P)=>{ setCal(im,{pw:P.d/P.px,ph:P.d/P.px,unit:P.u,src:`Set Scale (macro): ${P.px} px = ${P.d} ${P.u}`},true); },ij:P=>`run("Set Scale...", "distance=${P.px} known=${P.d} unit=${P.u==="µm"?"micron":P.u}");`};

/* ===== Process ===== */
function panelProcess(el,im){
  const list=Object.values(OPS).filter(o=>!o.hidden&&o.group&&GROUPS.find(g=>g[0]===o.group));
  if(!OPS[S.ui.proc]||OPS[S.ui.proc].hidden) S.ui.proc="gaussian"; const o=OPS[S.ui.proc];
  const prevAble=o.kind==="plane";
  el.innerHTML=`<div class="im-proc"><div class="im-oplist">${GROUPS.map(([g,gl])=>{ const L=list.filter(x=>x.group===g); if(!L.length) return ""; return `<div class="g">${gl}</div>${L.map(x=>`<button class="op ${x.id===o.id?"on":""}" data-op="${x.id}">${esc(x.label)}</button>`).join("")}`; }).join("")}</div>
    <div class="im-opform"><div class="im-h">${esc(o.label)}</div>${o.help?`<div class="xs muted im-help">${esc(o.help)}</div>`:""}<div id="imForm">${formHTML(o,im)}</div>
    <div class="im-row">${prevAble?`<label class="im-chk"><input type="checkbox" id="imPrev" ${S.ui.preview?"checked":""}> Preview</label>`:""}<span class="sp"></span><button class="btn primary sm" id="imApply">${o.kind==="new"?"Create":"Apply"}</button></div>
    <div class="xs muted">ImageJ: <code>${esc(o.ij?o.ij(Object.assign({},formVals(o))):"")}</code></div></div></div>`;
  $a(el,"[data-op]").forEach(b=>b.onclick=()=>{ S.ui.proc=b.dataset.op; clearPreview(); save(); renderPanel(); });
  bindForm(el,o,(re)=>{ if(re) renderPanel(); else { const c=$q(el,"code"); if(c&&o.ij) c.textContent=o.ij(formVals(o)); } schedulePreview(); });
  const pv=$q(el,"#imPrev"); if(pv) pv.onchange=e=>{ S.ui.preview=e.target.checked; save(); if(S.ui.preview) schedulePreview(); else clearPreview(); };
  $q(el,"#imApply").onclick=()=>runOp(o.id);
  if(S.ui.preview&&prevAble) schedulePreview();
}
let _pvT=null, _pvSeq=0;
function schedulePreview(){ clearTimeout(_pvT); const o=OPS[S.ui.proc]; if(!S.ui.preview||!o||o.kind!=="plane"||S.ui.panel!=="process") return; _pvT=setTimeout(runPreview,220); }
async function runPreview(){ const im=cur(); const o=OPS[S.ui.proc]; if(!im||!o) return; const seq=++_pvSeq; const P=Object.assign({},formVals(o),{__op:o.id,stack:false,__pv:true});
  try{ const r=await o.run(im,P); if(seq!==_pvSeq||!r||!r.planes||cur()!==im) return;
    const res={}; (im.rgb?[0,1,2]:[im.c]).forEach(c=>{ res[c]=r.planes[pIndex(im,c,im.z,im.t)]; });
    PREVIEW.imId=im.id; PREVIEW.planes=res; PREVIEW.plane=true; PREVIEW.c=im.c; PREVIEW.key="pv"+seq; DISP.key=null; draw(); }
  catch(e){ if(String(e.message)!=="취소") note("Preview: "+e.message,"warn"); }
}
const PVMODE=false; // (이전 전역 미리 보기 표시 — 이제 매개변수 __pv 로 전달)

/* ===== Threshold ===== */
const THR_LOCAL={key:"",mask:null};
async function updateAutoThr(im){ const T=S.ui.thr; if(T.manual||im.rgb) return; const t=await call("autoThreshold",[planeAt(im),T.method,!!T.dark]); T.lower=t.lower; T.upper=t.upper; T.level=t.level; OVL.key=null; if(t.none) note("이미지가 일정해 문턱을 구할 수 없습니다 (ImageJ 와 같음)","warn"); }
/* ===== Color Threshold (RGB) — ImageJ Image ▸ Adjust ▸ Color Threshold ===== */
/* 지금 장의 색 채널 캐시. 변환은 한 번에 하나만 Worker 로 — 도는 동안 바뀐 요청은 마지막 것만 이어서 (빠르게 바꿔도 일이 쌓이지 않음) */
const CT={key:"",ch:null,hist:null,busy:"",want:null};
const ctKey=im=>[im.id,im.ver,im.z,im.t,S.ui.cthr.space].join("|");
const ctReady=im=>!!CT.ch&&CT.key===ctKey(im);
function ctEnsure(im){ const k=ctKey(im); if(CT.key===k) return true; if(CT.busy){ CT.want=im.id; return false; } CT.busy=k; CT.want=null; const P=planesOf(im);
  const u=IC.rgbU8(P[pIndex(im,0,im.z,im.t)],P[pIndex(im,1,im.z,im.t)],P[pIndex(im,2,im.z,im.t)]);
  const next=()=>{ CT.busy=""; const w=CT.want&&byId(CT.want); CT.want=null; const c=cur(); if(c&&c.rgb&&S.ui.panel==="threshold"){ if(ctEnsure(c)){ OVL.key=null; renderPanel(); draw(); } else if(!CT.busy) renderPanel(); } void w; };
  call("colorHist",[u[0],u[1],u[2],S.ui.cthr.space]).then(r=>{ CT.key=k; CT.ch=r.ch; CT.hist=r.hist; next(); })
    .catch(e=>{ if(String(e.message)!=="취소") note("Color Threshold: "+e.message,"warn"); next(); }); return false; }
function ctDrop(id){ if(CT.key.startsWith(id+"|")){ CT.key=""; CT.ch=null; CT.hist=null; } }
/* 선택 화소 수 — 큰 이미지(> 4 M 화소)는 끄는 동안 듬성듬성 세어 '≈', 놓으면 정확히 */
function ctCount(rough){ const N=CT.ch[0].length; const st=rough&&N>4e6?Math.ceil(N/1e6):1; let m=0; for(let i=0;i<N;i+=st) if(ctInside(i)) m++; return st>1?Math.round(m*N/Math.ceil(N/st)):m; }
function ctCountText(n,rough){ const N=CT.ch[0].length; const ap=rough&&N>4e6?"≈ ":""; return `<b>${ap}${(100*n/N).toFixed(2)} %</b> selected (${ap}${n.toLocaleString()} px)`; }
function ctInside(i){ const C=S.ui.cthr; for(let k=0;k<3;k++){ const v=CT.ch[k][i]; const inn=v>=C.min[k]&&v<=C.max[k]; if(C.pass[k]?!inn:inn) return false; } return true; }
function panelColorThreshold(el,im){
  const C=S.ui.cthr; const ready=ctEnsure(im); const nm=IC.COLOR_CHANNEL_NAMES[C.space];
  const n=ready?ctCount(false):0;
  const hue=C.space==="HSB"?`<div class="im-huebar"></div>`:"";
  el.innerHTML=`<div class="im-h">Color Threshold <span class="muted xs">ImageJ Image ▸ Adjust ▸ Color Threshold</span></div>
  <div class="im-row"><label class="im-f"><span>Color space</span><select class="im-sel" id="ctSp">${IC.COLOR_SPACES.map(s=>`<option ${s===C.space?"selected":""}>${s}</option>`).join("")}</select></label>
    <label class="im-chk"><input type="checkbox" id="ctShow" ${S.ui.thr.show?"checked":""}> Show in red</label></div>
  ${!ready?`<div class="muted sm">색 변환 중…</div>`:[0,1,2].map(k=>`<div class="ct-ch" data-k="${k}"><div class="im-row xs"><b>${esc(nm[k])}</b><span class="sp"></span><label class="im-chk"><input type="checkbox" class="ctPass" data-k="${k}" ${C.pass[k]?"checked":""}> Pass</label><span class="muted">${C.pass[k]?"inside range":"outside range (stop)"}</span></div>
    ${histSVG({h:Array.from(CT.hist[k]),min:0,max:255},{lo:C.min[k],hi:C.max[k],log:S.ui.histLog,h:64})}${k===0?hue:""}
    <div class="im-slide"><span>min</span><input type="range" class="ctMin" data-k="${k}" min="0" max="255" step="1" value="${C.min[k]}"><input type="number" class="ctMinN" data-k="${k}" min="0" max="255" step="1" value="${C.min[k]}"></div>
    <div class="im-slide"><span>max</span><input type="range" class="ctMax" data-k="${k}" min="0" max="255" step="1" value="${C.max[k]}"><input type="number" class="ctMaxN" data-k="${k}" min="0" max="255" step="1" value="${C.max[k]}"></div></div>`).join("")}
  <div class="im-row"><label class="im-f"><span>Auto</span><select class="im-sel" id="ctM">${IC.THRESHOLD_METHODS.map(m=>`<option ${m===C.method?"selected":""}>${m}</option>`).join("")}</select></label>
    <label class="im-chk"><input type="checkbox" id="ctD" ${C.dark?"checked":""}> Dark background</label><button class="btn sm" id="ctAuto" ${ready?"":"disabled"} title="Auto threshold on ${C.space==="HSB"?"Brightness":C.space==="RGB"?"all three channels":"the first channel"} (as in ImageJ)">Auto</button>
    <button class="btn sm" id="ctSample" ${ready?"":"disabled"} title="Set the ranges from the colors inside the ROI (ImageJ Sample)">Sample (ROI)</button><button class="btn sm" id="ctReset">Reset</button></div>
  <div class="im-row xs"><span id="ctCount">${ready?ctCountText(n,false):""}</span></div>
  <div class="im-row">${isStack(im)?`<label class="im-chk"><input type="checkbox" id="ctSt" ${C.stack?"checked":""}> All slices</label>`:""}<label class="im-chk"><input type="checkbox" id="ctNew" ${C.newImage!==false?"checked":""}> New image</label><span class="sp"></span><button class="btn primary sm" id="ctApply" ${ready?"":"disabled"}>Create mask</button></div>
  <div class="xs muted">Hue 는 0–255 (빨강 0·255, 초록 85, 파랑 170). 빨강처럼 양 끝을 감싸는 색은 Hue 의 Pass 를 끄고(stop) 가운데를 제외하세요. 결과 mask 로 Analyze Particles (Redirect 로 원본 측정).</div>`;
  const q=s=>$q(el,s), qa=s=>$a(el,s);
  const upd=(full)=>{ OVL.key=null; save(); if(full){ renderPanel(); draw(); return; } drawSoon(); if(ctReady(im)){ const c=q("#ctCount"); if(c) c.innerHTML=ctCountText(ctCount(true),true); } };
  q("#ctSp").onchange=e=>{ C.space=e.target.value; C.min=[0,0,0]; C.max=[255,255,255]; C.pass=[true,true,true]; upd(true); };
  q("#ctShow").onchange=e=>{ S.ui.thr.show=e.target.checked; save(); draw(); };
  const setR=(cls,k,v,full)=>{ v=clamp(Math.round(+v),0,255); if(cls==="min"){ C.min[k]=v; if(C.max[k]<v){ C.max[k]=v; const o=q(`.ctMax[data-k="${k}"]`), p=q(`.ctMaxN[data-k="${k}"]`); if(o) o.value=v; if(p) p.value=v; } } else { C.max[k]=v; if(C.min[k]>v){ C.min[k]=v; const o=q(`.ctMin[data-k="${k}"]`), p=q(`.ctMinN[data-k="${k}"]`); if(o) o.value=v; if(p) p.value=v; } } /* ImageJ 처럼 min ≤ max */ const a=q(`.ct${cls==="min"?"Min":"Max"}[data-k="${k}"]`), b=q(`.ct${cls==="min"?"Min":"Max"}N[data-k="${k}"]`); if(a) a.value=v; if(b) b.value=v; upd(full); };
  qa(".ctMin").forEach(x=>{ x.oninput=()=>setR("min",+x.dataset.k,x.value,false); x.onchange=()=>renderPanel(); });
  qa(".ctMax").forEach(x=>{ x.oninput=()=>setR("max",+x.dataset.k,x.value,false); x.onchange=()=>renderPanel(); });
  qa(".ctMinN").forEach(x=>x.onchange=()=>setR("min",+x.dataset.k,x.value,true)); qa(".ctMaxN").forEach(x=>x.onchange=()=>setR("max",+x.dataset.k,x.value,true));
  qa(".ctPass").forEach(x=>x.onchange=()=>{ C.pass[+x.dataset.k]=x.checked; upd(true); });
  q("#ctM").onchange=e=>{ C.method=e.target.value; save(); }; q("#ctD").onchange=e=>{ C.dark=e.target.checked; save(); };
  q("#ctAuto").onclick=()=>{ if(!ctReady(im)){ note("색 변환이 끝난 뒤에 다시 누르세요","warn"); return; } const a=IC.colorAuto(CT.hist,C.space,C.method,C.dark); for(const k in a){ C.min[k]=a[k][0]; C.max[k]=a[k][1]; } upd(true); };
  q("#ctSample").onclick=()=>{ if(!ctReady(im)){ note("색 변환이 끝난 뒤에 다시 누르세요","warn"); return; } const m=roisMask(im,false); if(!m){ note("Sample: 먼저 영역 ROI 를 그리세요 (사각형·타원·다각형)","warn"); return; } const idx=[]; for(let i=0;i<m.length;i++) if(m[i]) idx.push(i);
    try{ const s=IC.colorSample(CT.ch,idx,C.space); C.min=s.min; C.max=s.max; C.pass[0]=s.pass0; upd(true); note(`Sample: ${idx.length.toLocaleString()} px${C.space==="HSB"&&!s.pass0?" — Hue 가 0/255 를 감싸 stop 으로":""}`,"ok"); }catch(e){ note(e.message,"warn"); } };
  q("#ctReset").onclick=()=>{ C.min=[0,0,0]; C.max=[255,255,255]; C.pass=[true,true,true]; upd(true); };
  if(q("#ctSt")) q("#ctSt").onchange=e=>{ C.stack=e.target.checked; save(); }; q("#ctNew").onchange=e=>{ C.newImage=e.target.checked; save(); };
  q("#ctApply").onclick=()=>runOp("colorThreshold",{space:C.space,min:C.min.slice(),max:C.max.slice(),pass:C.pass.slice(),stack:!!C.stack,newImage:C.newImage!==false});
}
function panelThreshold(el,im){
  const T=S.ui.thr; if(im.rgb){ panelColorThreshold(el,im); return; }
  const p=planeAt(im); const h=IC.histogram(p,null,256); const [mn,mx]=[h.min,h.max];
  if(T.forImg!==im.id+"|"+im.ver+"|"+im.c+"|"+im.z+"|"+im.t+"|"+T.method+"|"+T.dark&&!T.manual){ T.forImg=im.id+"|"+im.ver+"|"+im.c+"|"+im.z+"|"+im.t+"|"+T.method+"|"+T.dark; updateAutoThr(im).then(()=>{ renderPanel(); draw(); }); }
  const lo=clamp(+T.lower,mn,mx), hi=clamp(+T.upper,mn,mx); let n=0; for(let i=0;i<p.data.length;i++){ const v=p.data[i]; if(v>=T.lower&&v<=T.upper) n++; }
  const L=T.local||{method:"Phansalkar",radius:15,p1:"",p2:""};
  el.innerHTML=`
  <div class="im-h">Global threshold</div>
  ${histSVG(h,{lo,hi,log:S.ui.histLog})}
  <div class="im-row"><label class="im-f"><span>Method</span><select class="im-sel" id="imTm" ${T.manual?"disabled":""}>${IC.THRESHOLD_METHODS.map(m=>`<option ${m===T.method?"selected":""}>${m}</option>`).join("")}</select></label>
    <label class="im-chk"><input type="checkbox" id="imTd" ${T.dark?"checked":""}> Dark background</label><label class="im-chk"><input type="checkbox" id="imTman" ${T.manual?"checked":""}> Manual</label></div>
  <div class="im-slide"><span>Lower</span><input type="range" id="imTl" min="${mn}" max="${mx}" step="${(mx-mn)/1000||1}" value="${lo}" ${T.manual?"":"disabled"}><input type="number" id="imTlN" value="${+(+T.lower).toPrecision(7)}" step="any" ${T.manual?"":"disabled"}></div>
  <div class="im-slide"><span>Upper</span><input type="range" id="imTu" min="${mn}" max="${mx}" step="${(mx-mn)/1000||1}" value="${hi}" ${T.manual?"":"disabled"}><input type="number" id="imTuN" value="${+(+T.upper).toPrecision(7)}" step="any" ${T.manual?"":"disabled"}></div>
  <div class="im-row xs"><span><b>${(100*n/p.data.length).toFixed(2)} %</b> selected (${n.toLocaleString()} px)</span>${n/p.data.length>0.5&&!T.manual?`<span class="note warn xs" style="flex-basis:100%">절반이 넘게 선택됨 — 입자가 아니라 바탕을 고른 것 같습니다. ${T.dark?"밝은 바탕에 어두운 입자(TEM 명시야·SEM 반사상 아님)면 <b>Dark background</b> 를 끄세요":"어두운 바탕에 밝은 입자(형광·STEM HAADF)면 <b>Dark background</b> 를 켜세요"}.</span>`:""}${T.level!=null&&!T.manual?`<span class="muted">level ${T.level}/255${im.bit!==8?" (histogram scaled to 8-bit)":""}</span>`:""}<span class="sp"></span><label class="im-chk"><input type="checkbox" id="imTshow" ${T.show?"checked":""}> Show in red</label></div>
  <div class="im-row"><button class="btn sm" id="imTall" title="Show all 16 methods side by side (ImageJ Try all)">Try all methods</button><span class="sp"></span><label class="im-chk"><input type="checkbox" id="imTnew" ${T.newImage!==false?"checked":""}> New image (mask)</label>${isStack(im)?`<label class="im-chk"><input type="checkbox" id="imTst" ${T.stack?"checked":""}> All slices</label>`:""}<button class="btn primary sm" id="imTapply">Convert to Mask</button></div>
  <div id="imTry"></div>
  <hr class="im-hr"><div class="im-h">Auto Local Threshold <span class="muted xs">8-bit · 원형 창</span></div>
  <div class="im-row"><label class="im-f"><span>Method</span><select class="im-sel" id="imLm">${IC.LOCAL_METHODS.map(m=>`<option ${m===L.method?"selected":""}>${m}</option>`).join("")}</select></label><label class="im-f"><span>Radius</span><input type="number" id="imLr" value="${L.radius}" min="1"><span class="u">px</span></label></div>
  <div class="im-row"><label class="im-f"><span>Parameter 1</span><input type="number" id="imL1" value="${esc(L.p1)}" step="any" placeholder="default"></label><label class="im-f"><span>Parameter 2</span><input type="number" id="imL2" value="${esc(L.p2)}" step="any" placeholder="default"></label></div>
  <div class="xs muted">${esc(LOCAL_HELP[L.method]||"")}</div>
  <div class="im-row"><button class="btn sm" id="imLprev">Preview</button>${T.localPreview?`<button class="btn sm" id="imLoff">Preview off</button>`:""}<span class="sp"></span><button class="btn primary sm" id="imLapply" ${im.bit!==8?"disabled title='8-bit images only'":""}>Apply local</button></div>`;
  const q=s=>$q(el,s);
  q("#imTm").onchange=e=>{ T.method=e.target.value; save(); renderPanel(); }; q("#imTd").onchange=e=>{ T.dark=e.target.checked; save(); renderPanel(); };
  q("#imTman").onchange=e=>{ T.manual=e.target.checked; save(); renderPanel(); };
  const setT=(k,v,fromN)=>{ T[k]=+v; OVL.key=null; if(fromN){ draw(); renderPanel(); } else drawSoon(); save(); };
  q("#imTl").oninput=e=>{ setT("lower",e.target.value); q("#imTlN").value=+(+e.target.value).toPrecision(7); }; q("#imTu").oninput=e=>{ setT("upper",e.target.value); q("#imTuN").value=+(+e.target.value).toPrecision(7); };
  q("#imTl").onchange=()=>renderPanel(); q("#imTu").onchange=()=>renderPanel();
  q("#imTlN").onchange=e=>setT("lower",e.target.value,true); q("#imTuN").onchange=e=>setT("upper",e.target.value,true);
  q("#imTshow").onchange=e=>{ T.show=e.target.checked; save(); draw(); };
  q("#imTnew").onchange=e=>{ T.newImage=e.target.checked; save(); }; if(q("#imTst")) q("#imTst").onchange=e=>{ T.stack=e.target.checked; save(); };
  q("#imTapply").onclick=()=>runOp("threshold",{method:T.method,dark:!!T.dark,manual:!!T.manual,lower:+T.lower,upper:+T.upper,stack:!!T.stack,newImage:T.newImage!==false});
  q("#imTall").onclick=()=>tryAll(im,q("#imTry"));
  const LL=()=>({method:q("#imLm").value,radius:+q("#imLr").value||15,p1:q("#imL1").value,p2:q("#imL2").value});
  q("#imLm").onchange=()=>{ T.local=LL(); save(); renderPanel(); }; ["#imLr","#imL1","#imL2"].forEach(s=>q(s).onchange=()=>{ T.local=LL(); save(); });
  q("#imLprev").onclick=async()=>{ if(im.bit!==8){ note("Local threshold 는 8-bit 에서만 — Adjust ▸ Type ▸ 8-bit","warn"); return; } T.local=LL(); setProg("Local threshold…",0.3);
    try{ const m=await call("localThreshold",[planeAt(im),T.local.method,{radius:T.local.radius,p1:T.local.p1===""?null:+T.local.p1,p2:T.local.p2===""?null:+T.local.p2,dark:T.dark}]); hideProg(); THR_LOCAL.mask=m.data; THR_LOCAL.key=JSON.stringify(T.local)+im.ver; T.localPreview=THR_LOCAL.key; OVL.key=null; draw(); renderPanel(); }catch(e){ hideProg(); note(e.message,"warn"); } };
  if(q("#imLoff")) q("#imLoff").onclick=()=>{ T.localPreview=""; OVL.key=null; draw(); renderPanel(); };
  q("#imLapply").onclick=()=>{ T.local=LL(); runOp("threshold",{dark:!!T.dark,stack:!!T.stack,newImage:T.newImage!==false,local:T.local}).then(()=>{ T.localPreview=""; const m=cur(); if(!m) return; const pl=curPlaneOf(m); if(!IC.isBinary(pl)) return; let n=0; for(let i=0;i<pl.data.length;i++) if(pl.data[i]>0) n++; const f=n/pl.data.length; if(f>0.5) note(`Local 문턱 결과 화소의 ${(100*f).toFixed(0)} % 가 선택됨 — 바탕까지 고른 것 같습니다. 먼저 Process ▸ Gaussian Blur(σ 1–2) 를 하거나 Dark background·방법·반지름을 바꿔 보세요`,"warn"); }); };
  // 히스토그램 끌기 (manual)
  const hs=q(".im-hist"); hs.onpointerdown=e=>{ if(!T.manual){ T.manual=true; } const r=hs.getBoundingClientRect(); const f=x=>mn+(x-r.left)/r.width*(mx-mn); const a=f(e.clientX); const k=Math.abs(a-T.lower)<Math.abs(a-T.upper)?"lower":"upper"; T[k]=a; OVL.key=null; draw();
    hs.setPointerCapture(e.pointerId); hs.onpointermove=ev=>{ T[k]=clamp(f(ev.clientX),mn,mx); OVL.key=null; draw(); }; hs.onpointerup=()=>{ hs.onpointermove=null; save(); renderPanel(); }; };
}
const LOCAL_HELP={Niblack:"t = mean + k·SD − c   (P1 = k, 기본 0.2; P2 = c, 기본 0)",Sauvola:"t = mean·(1 + k·(SD/r − 1))   (P1 = k 0.5, P2 = r 128)",Phansalkar:"저대비 이미지용 Sauvola 변형 (P1 = k 0.25, P2 = r 0.5)",Bernsen:"대비 < P1(15) 이면 중간값 기준, 아니면 (max+min)/2",Mean:"화소 > 지역 평균 − P1",Median:"화소 > 지역 중앙값 − P1",MidGrey:"화소 > (max+min)/2 − P1",Contrast:"지역 max·min 중 가까운 쪽"};
async function tryAll(im,host){ host.innerHTML=`<div class="xs muted">계산 중…</div>`; const p=planeAt(im); const T=S.ui.thr;
  const all=await call("tryAllThresholds",[p,!!T.dark]); const sc=Math.max(1,Math.ceil(Math.max(im.w,im.h)/150)); const w=Math.ceil(im.w/sc), h=Math.ceil(im.h/sc);
  host.innerHTML=`<div class="im-try">${IC.THRESHOLD_METHODS.map(m=>`<button class="t" data-m="${m}" title="${m}: [${g4(all[m].lower)}, ${g4(all[m].upper)}]"><canvas width="${w}" height="${h}"></canvas><span>${m}</span></button>`).join("")}</div><div class="xs muted">클릭하면 그 방법을 고릅니다.</div>`;
  $a(host,".t").forEach(b=>{ const r=all[b.dataset.m]; const cv=$q(b,"canvas"), g=cv.getContext("2d"), id=g.createImageData(w,h); for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const v=p.data[(y*sc)*im.w+x*sc]; const on=v>=r.lower&&v<=r.upper; const j=(y*w+x)*4; id.data[j]=on?255:20; id.data[j+1]=on?255:20; id.data[j+2]=on?255:20; id.data[j+3]=255; } g.putImageData(id,0,0);
    b.onclick=()=>{ T.method=b.dataset.m; T.manual=false; save(); renderPanel(); }; }); }

/* ===== Binary ===== */
function panelBinary(el,im){
  const bin=!im.rgb&&IC.isBinary(planeAt(im)); const v=formVals(OPS.binary);
  el.innerHTML=`${bin?"":`<div class="note warn sm">이진(0/255) 8-bit 이미지가 필요합니다. <b>Threshold ▸ Convert to Mask</b> 먼저. <button class="btn xs" id="imGoThr">Go to Threshold</button></div>`}
  <div class="im-h">Morphology <span class="muted xs">물체 = 255 (Black background)</span></div>
  <div class="im-row"><label class="im-f"><span>Iterations</span><input type="number" id="imBi" value="${v.iterations}" min="1" step="1"></label><label class="im-f"><span>Count</span><input type="number" id="imBc" value="${v.count}" min="1" max="8" step="1"></label>${isStack(im)?`<label class="im-chk"><input type="checkbox" id="imBst" ${v.stack?"checked":""}> All slices</label>`:""}</div>
  <div class="im-grid">${[["erode","Erode","Erode object edges by 1 px (needs ≥ count background neighbours)"],["dilate","Dilate","Grow objects by 1 px"],["open","Open","Erode → Dilate: removes small protrusions and specks"],["close","Close","Dilate → Erode: fills small gaps"]].map(([k,l,t])=>`<button class="btn sm" data-b="${k}" title="${t}" ${bin?"":"disabled"}>${l}</button>`).join("")}</div>
  <div class="im-h">Shape</div>
  <div class="im-grid">${[["fill","Fill Holes","Fill holes inside objects"],["outline","Outline","1 px inner outline of objects"],["skeleton","Skeletonize","1 px skeleton"],["watershed","Watershed","Split touching particles (EDM maxima, tolerance 0.5)"]].map(([k,l,t])=>`<button class="btn sm" data-b="${k}" title="${t}" ${bin?"":"disabled"}>${l}</button>`).join("")}</div>
  <div class="im-h">Distance</div>
  <div class="im-grid"><button class="btn sm" data-b="edm" ${bin?"":"disabled"} title="Distance from each object pixel to the nearest background (8-bit, clipped at 255)">Distance Map</button><button class="btn sm" data-b="edm32" ${bin?"":"disabled"} title="New 32-bit image (not clipped)">Distance Map (32-bit)</button><button class="btn sm" data-b="uep" ${bin?"":"disabled"} title="EDM maxima = one point per object, value = inscribed-circle radius">Ultimate Points</button></div>
  <div class="im-row"><button class="btn sm" id="imBinv">Invert</button><span class="xs muted">ImageJ 와 같은 규칙(이미지 밖 = 배경, 8-연결). Watershed 는 볼록한 입자가 붙은 경우에 맞습니다.</span></div>`;
  const q=s=>$q(el,s); if(q("#imGoThr")) q("#imGoThr").onclick=()=>{ S.ui.panel="threshold"; rerender(); };
  q("#imBi").onchange=e=>{ v.iterations=+e.target.value||1; }; q("#imBc").onchange=e=>{ v.count=clamp(+e.target.value||1,1,8); }; if(q("#imBst")) q("#imBst").onchange=e=>{ v.stack=e.target.checked; };
  $a(el,"[data-b]").forEach(b=>b.onclick=()=>{ const k=b.dataset.b; runOp("binary",Object.assign({},v,{fn:k==="edm32"?"edm":k,edm32:k==="edm32"})); });
  q("#imBinv").onclick=()=>runOp("invert",{stack:!!v.stack});
}

/* ===== Info ===== */
function panelInfo(el,im){ const c=calOf(im); const tot=planesOf(im).reduce((a,p)=>a+p.byteLength,0);
  el.innerHTML=`<div class="im-h">${esc(im.name)}</div>
  <table class="im-kv"><tr><td>Size</td><td>${im.w} × ${im.h} px${im.cal?` (${g4(im.w*c.pw)} × ${g4(im.h*c.ph)} ${esc(c.unit)})`:""}</td></tr><tr><td>Type</td><td>${im.rgb?"RGB":im.bit+"-bit"}${im.nc>1&&!im.rgb?`, ${im.nc} channels`:""}${isStack(im)?`, ${im.nz} slices × ${im.nt} frames`:""}</td></tr>
  <tr><td>Pixel size</td><td>${im.cal?`${g4(c.pw)} × ${g4(c.ph)} ${esc(c.unit)}`:"none"}${im.cal?`<div class="xs muted">${esc(im.cal.src||"")}</div>`:""}</td></tr>${im.chNames?`<tr><td>Channels</td><td>${im.chNames.map((n,i)=>`<span class="lutdot" style="background:${lutCss(im.disp[i]&&im.disp[i].lut||"Grays")}"></span>${esc(n)}`).join(" · ")}</td></tr>`:""}${im.zStep>0&&im.nz>1?`<tr><td>Z step</td><td>${g4(im.zStep)} µm</td></tr>`:""}<tr><td>Memory</td><td>${(tot/1048576).toFixed(1)} MB</td></tr><tr><td>Source</td><td>${esc(im.file||"")}${im.parent&&byId(im.parent)?` ← ${esc(byId(im.parent).name)}`:""}</td></tr>
  ${(im.info||[]).map(([k,v])=>`<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join("")}</table>
  ${(im.warn||[]).map(w=>`<div class="note warn sm">${esc(w)}</div>`).join("")}
  <div class="im-h">Processing history <span class="muted xs">${im.ops.length}</span></div>${im.ops.length?`<ol class="im-ops">${im.ops.map(o=>`<li>${esc(o.label||o.op)}</li>`).join("")}</ol>`:`<div class="xs muted">unprocessed</div>`}
  <div class="im-row"><button class="btn sm" id="imRen">Rename</button><button class="btn sm" id="imToMacro" ${im.ops.length?"":"disabled"} title="Convert the operations applied to this image into a macro">History → Macro</button></div>
  ${im.metaText||im.desc?`<details><summary class="sm">Original file metadata</summary><pre class="im-pre">${esc((im.metaText||im.desc).slice(0,20000))}</pre></details>`:""}`;
  $q(el,"#imRen").onclick=()=>{ const n=prompt("New name",im.name); if(n){ im.name=n; save(); rerender(); } };
  $q(el,"#imToMacro").onclick=()=>{ if(S.macro.steps.length&&!confirm2("tomac",`지금 매크로 ${S.macro.steps.length} 단계를 이 이미지의 작업으로 바꿀까요? 한 번 더`)) return; S.macro.steps=im.ops.map(o=>Object.assign({},o)); S.ui.panel="macro"; save(); rerender(); note(`${im.ops.length} 단계를 매크로로 가져왔습니다`,"ok"); };
}
