/* ============================================================
   화면 뼈대 · 파일 열기 · 영상 목록 · 도구 막대 · 스택 막대 · 단축키 · 외부 API
   ============================================================ */
const svgI=d=>`<svg viewBox="0 0 24 24" width="17" height="17" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="${d}"/></svg>`;
function render(root){
  ROOT=root;
  if(!LOADED){ LOADED=true; root.innerHTML=`<div class="card"><div class="im-empty">Image 작업을 불러오는 중…</div></div>`; loadSaved().then(()=>{ if(ROOT===root) build(root); }); return; }
  build(root);
}
function build(root){
  ROOT=root;
  root.innerHTML=`<div class="imx">
    <div class="card im-top">
      <div class="grp"><label class="btn primary sm" style="cursor:pointer" title="TIFF·ND2·CZI·PNG·JPEG·DM3/DM4·MRC·EMD·GWY·Nanoscope·ARDF·IBW·SXM·NPY·텍스트 … 여러 개·폴더째 끌어다 놓아도 됩니다">파일 열기<input type="file" id="imFin" multiple accept="${IO.ACCEPT}" hidden></label>
        <button class="btn sm" id="imFmts" title="읽는 형식">형식</button><span class="muted sm">또는 화면에 끌어다 놓기</span></div>
      <div class="grp" style="margin-left:auto"><button class="btn sm" id="imUndo" title="Ctrl+Z">되돌리기</button><button class="btn sm" id="imSaveT" title="ImageJ 호환 TIFF (보정·스택 포함)">TIFF 저장</button><button class="btn sm" id="imSaveP" title="보이는 그대로(LUT·겹쳐 그린 것·축척 막대) 원래 해상도 PNG">PNG 내보내기</button><button class="btn sm" id="imClose" title="현재 영상 닫기">닫기</button><button class="btn sm" id="imClear" title="모두 지우기 (두 번 눌러 확인)">모두 지우기</button></div>
    </div>
    <div class="im-main">
      <div class="card im-list" id="imList"></div>
      <div class="card im-center">
        <div class="im-toolbar" id="imTools"></div>
        <div class="im-canvaswrap" id="imCanvasWrap"><canvas id="imCanvas"></canvas><canvas id="imOverlay"></canvas><div class="im-drop" id="imDropHint" hidden>여기에 놓기</div></div>
        <div class="im-stackbar" id="imStack"></div>
        <div class="im-status" id="imStatus"></div>
        <div class="im-prog" id="imProg" hidden><span class="t"></span><div class="bar"><i></i></div><button class="btn xs" id="imStop">중단</button></div>
        <div class="im-msg" id="imMsg" hidden></div>
        <div class="im-profile" id="imProfile" hidden></div>
      </div>
      <div class="card im-side"><div class="seg im-tabs" id="imTabs"></div><div id="imPanel"></div></div>
    </div>
    <div class="card im-bottom" id="imBottom"></div>
    <div class="im-modal" id="imModal" hidden></div>
  </div>`;
  const x=$q(root,".imx"), q=s=>$q(root,s);
  q("#imFin").onchange=e=>{ loadFiles(e.target.files); e.target.value=""; };
  q("#imFmts").onclick=showFormats;
  q("#imUndo").onclick=undo; q("#imSaveT").onclick=saveTIFF; q("#imSaveP").onclick=exportPNG;
  q("#imClose").onclick=()=>{ const im=cur(); if(im){ closeImage(im.id); rerender(); } };
  q("#imClear").onclick=()=>{ if(!S.images.length) return; if(BUSY){ note("계산 중에는 지울 수 없습니다 — 끝나거나 중단한 뒤에","warn"); return; } if(!confirm2("imclr","열린 영상·ROI·결과를 모두 지울까요? 다시 누르면 지웁니다")) return; S.images.map(i=>i.id).forEach(id=>{ PX.delete(id); DIRTY.add(id); }); const ui=S.ui; S=defState(); S.ui.panel=ui.panel; save(); rerender(); };
  q("#imStop").onclick=()=>{ cancelAll(); note("계산을 중단했습니다"); };
  ["dragenter","dragover"].forEach(ev=>x.addEventListener(ev,e=>{ if(e.dataTransfer&&[...(e.dataTransfer.types||[])].includes("Files")){ e.preventDefault(); e.stopPropagation(); x.classList.add("over"); } }));
  x.addEventListener("dragleave",e=>{ if(e.target===x) x.classList.remove("over"); });
  x.addEventListener("drop",e=>{ if(e.dataTransfer&&e.dataTransfer.files&&e.dataTransfer.files.length){ e.preventDefault(); e.stopPropagation(); x.classList.remove("over"); loadFiles(e.dataTransfer.files); } });
  const ov=q("#imOverlay"); ov.addEventListener("pointerdown",onDown); ov.addEventListener("pointermove",onMove); ov.addEventListener("pointerup",onUp); ov.addEventListener("dblclick",onDbl);
  ov.addEventListener("pointerleave",()=>{ HOVER=null; drawStatus(); }); ov.addEventListener("wheel",onWheel,{passive:false}); ov.addEventListener("contextmenu",e=>e.preventDefault());
  try{ new ResizeObserver(()=>draw()).observe(q("#imCanvasWrap")); }catch(e){}
  rerender();
}
function rerender(){ if(!ROOT||!$q(ROOT,".imx")) return; drawSide(); drawTools(); drawStack(); drawTabs(); renderPanel(); draw(); drawBottom(); drawMsg(); }

/* ---------- 왼쪽: 영상 목록 ---------- */
const THUMB=new Map();
function thumbOf(im){ const k=im.id+"|"+im.ver+"|"+JSON.stringify(im.disp)+"|"+im.composite; const t=THUMB.get(im.id); if(t&&t.k===k) return t.url;
  try{ const sc=Math.max(1,Math.floor(Math.max(im.w,im.h)/112)); const src=renderRegion(im,{x:0,y:0,w:im.w,h:im.h,sc}); const s=56/Math.max(im.w,im.h); const c=document.createElement("canvas"); c.width=Math.max(1,Math.round(im.w*s)); c.height=Math.max(1,Math.round(im.h*s)); // 작은 표본에서 바로 (전체 표시를 다시 만들지 않음)
    const g=c.getContext("2d"); g.imageSmoothingQuality="high"; g.drawImage(src,0,0,c.width,c.height); const url=c.toDataURL(); THUMB.set(im.id,{k,url}); return url; }catch(e){ return ""; } }
function drawSide(){ const el=ROOT&&$q(ROOT,"#imList"); if(!el) return;
  if(!S.images.length){ el.innerHTML=`<div class="im-dropbox" id="imDrop">여기에 영상 파일을 놓거나<br><b>눌러서 열기</b><br><span class="xs">TIFF · ND2 · CZI · PNG · DM3 · MRC · EMD · GWY · Nanoscope · ARDF · IBW · SXM …</span></div>
    <ol class="im-steps0"><li><b>열기</b> — SEM·TEM·AFM·광학 영상. 화소 크기는 메타데이터에서 자동으로</li><li><b>보정 확인</b> — 없으면 Set Scale 도구로 축척 막대에 선</li><li><b>Threshold</b> → <b>Analyze Particles</b> — 크기 분포·면적 분율</li><li><b>Macro</b> — 같은 처리를 다른 영상·파일에 batch</li></ol>`; $q(el,"#imDrop").onclick=()=>$q(ROOT,"#imFin").click(); return; }
  const depth=im=>{ let d=0, p=im.parent; while(p&&byId(p)&&d<4){ d++; p=byId(p).parent; } return d; };
  el.innerHTML=`<div class="im-lh"><b>Images</b><span class="muted xs">${S.images.length}</span></div>`+S.images.map(im=>`<div class="it ${im.id===S.sel?"sel":""}" data-id="${im.id}" style="padding-left:${6+depth(im)*12}px" title="${esc(im.name)}"><img src="${thumbOf(im)}" alt=""><div class="tx"><div class="nm">${esc(im.name)}</div><div class="xs muted">${im.w}×${im.h} · ${im.rgb?"RGB":im.bit+"-bit"}${isStack(im)?` · ${im.nz*im.nt}`:""}${im.nc>1&&!im.rgb?` · ${im.nc}ch`:""}${im.cal?` · ${g4(im.cal.pw)} ${esc(im.cal.unit)}`:` · <span class="warnc">no scale</span>`}</div></div><button class="x" data-x="${im.id}" title="닫기">×</button></div>`).join("");
  $a(el,".it").forEach(d=>d.onclick=e=>{ if(e.target.dataset.x){ closeImage(e.target.dataset.x); rerender(); return; } if(S.sel!==d.dataset.id){ S.sel=d.dataset.id; clearPreview(); invalidate(); save(); rerender(); } }); }

/* ---------- 도구 막대 ---------- */
function drawTools(){ const el=ROOT&&$q(ROOT,"#imTools"); if(!el) return; const im=cur();
  el.innerHTML=`<div class="tools">${TOOLS.map(([k,l,d])=>`<button class="tb ${S.ui.tool===k?"on":""}" data-t="${k}" title="${esc(l)}">${svgI(d)}</button>`).join("")}</div>
    <div class="sp"></div>
    ${im?`<span class="im-name">${esc(im.name)}</span><button class="btn xs" id="imZo" title="−">−</button><button class="btn xs" id="imZ1" title="100 %">1:1</button><button class="btn xs" id="imZf" title="화면에 맞추기">Fit</button><button class="btn xs" id="imZi" title="+">+</button><button class="btn xs" id="imSn" title="선택 해제 (Esc)">Select none</button>`:""}`;
  $a(el,"[data-t]").forEach(b=>b.onclick=()=>{ S.ui.tool=b.dataset.t; save(); drawTools(); if(S.ui.tool==="wand"&&S.ui.panel!=="threshold"&&cur()&&!IC.isBinary(curPlaneOf(cur()))) note("Wand: Threshold 탭이 열려 있으면 그 범위로, 아니면 클릭한 값(± Measure 탭의 tolerance)으로 따라갑니다"); });
  if(!im) return; const q=s=>$q(el,s);
  q("#imZo").onclick=()=>zoomStep(im,-1); q("#imZi").onclick=()=>zoomStep(im,1); q("#imZf").onclick=()=>{ fitView(im); draw(); };
  q("#imZ1").onclick=()=>{ const v=viewOf(im); const {W,H}=canvasSize(); const cx=v.x+W/v.zoom/2, cy=v.y+H/v.zoom/2; v.zoom=1; v.x=cx-W/2; v.y=cy-H/2; draw(); };
  q("#imSn").onclick=()=>{ setActive(im,null); };
}
/* ---------- 스택 막대 ---------- */
let PLAY=null;
function drawStack(){ const el=ROOT&&$q(ROOT,"#imStack"); if(!el) return; const im=cur(); if(!im||(!isStack(im)&&(im.nc<2||im.rgb))){ el.hidden=true; el.innerHTML=""; return; } el.hidden=false;
  const cn=c=>im.chNames&&im.chNames[c]?` ${esc(im.chNames[c])}`:"";
  const row=(k,lab,n,v)=>n>1?`<div class="sr"><b>${lab}</b><input type="range" min="1" max="${n}" value="${v+1}" data-d="${k}"><span class="v">${v+1}/${n}${k==="c"?cn(v):""}</span></div>`:"";
  el.innerHTML=`${!im.rgb?row("c","c",im.nc,im.c):""}${row("z","z",im.nz,im.z)}${row("t","t",im.nt,im.t)}${isStack(im)?`<button class="btn xs" id="imPlay">${PLAY?"■ Stop":"▶ Play"}</button><label class="im-f"><span>fps</span><input type="number" id="imFps" value="${S.ui.fps||7}" min="1" max="60" style="width:48px"></label>`:""}`;
  $a(el,"input[type=range]").forEach(r=>r.oninput=()=>{ im[r.dataset.d]=+r.value-1; $q(r.parentNode,".v").textContent=`${r.value}/${r.max}${r.dataset.d==="c"&&im.chNames&&im.chNames[+r.value-1]?" "+im.chNames[+r.value-1]:""}`; clearPreview(); invalidate(); draw(); });
  $a(el,"input[type=range]").forEach(r=>r.onchange=()=>{ save(); renderPanel(); });
  const pb=$q(el,"#imPlay"); if(pb) pb.onclick=()=>{ if(PLAY){ clearInterval(PLAY); PLAY=null; drawStack(); return; } const k=im.nz>1?"z":"t", n=im.nz>1?im.nz:im.nt; PLAY=setInterval(()=>{ if(cur()!==im){ clearInterval(PLAY); PLAY=null; return; } im[k]=(im[k]+1)%n; invalidate(); draw(); const r=$q(el,`input[data-d="${k}"]`); if(r){ r.value=im[k]+1; $q(r.parentNode,".v").textContent=`${im[k]+1}/${n}`; } },1000/(S.ui.fps||7)); drawStack(); };
  const fp=$q(el,"#imFps"); if(fp) fp.onchange=e=>{ S.ui.fps=clamp(+e.target.value||7,1,60); save(); if(PLAY){ clearInterval(PLAY); PLAY=null; pb.click(); } }; // 재생 중 fps 바꾸면 바로 반영
}

/* ---------- 파일 열기 ---------- */
async function loadFiles(list){
  const arr=[...list]; if(!arr.length) return;
  const base=n=>n.replace(/\.[^.]+$/,"").toLowerCase();
  const imgNames=new Set(arr.filter(f=>!/\.txt$/i.test(f.name)).map(f=>base(f.name)));
  const side=arr.filter(f=>/\.txt$/i.test(f.name)&&imgNames.has(base(f.name)));
  const sidecars={}; for(const f of side) sidecars[f.name.toLowerCase()]=IO._txt(new Uint8Array(await f.arrayBuffer())); // UTF-8 아니면 latin1 (µ)
  const files=arr.filter(f=>!side.includes(f));
  let n=0; const errs=[];
  for(let k=0;k<files.length;k++){ const f=files[k]; setProg(`읽는 중 ${k+1}/${files.length}: ${f.name}`,k/files.length);
    try{ const u8=new Uint8Array(await f.arrayBuffer()); let r;
      try{ r=await IO.readFile(f.name,u8,{sidecars,choose:askImport,onProgress:p=>setProg(`읽는 중: ${f.name}`,(k+p)/files.length)}); if(r.cancel){ note(`${f.name}: 열기 취소`); continue; } }
      catch(e){ if(e.needRaw){ hideProg(); const o=await askRaw(f.name,u8.length); if(!o) continue; r=IO.readRaw(f.name,u8,o); } else throw e; }
      for(const d of r.images){ const im=addImage(d,{file:f.name}); n++; if(im.cal) im.info.unshift(["Scale",`${g4(im.cal.pw)} ${im.cal.unit}/px`]); (d.warn||[]).forEach(w=>errs.push(`${d.name}: ${w}`)); }
      (r.warn||[]).forEach(w=>errs.push(`${f.name}: ${w}`));
    }catch(e){ errs.push(`${f.name}: ${e.message||e}`); }
    await sleep(0); }
  hideProg(); save(); rerender();
  if(n){ const im=cur(); note(`영상 ${n}개를 열었습니다${side.length?` (옆 파일 ${side.length}개 메타데이터 사용)`:""}${im&&!im.cal?" — 화소 크기가 없습니다: Adjust ▸ Set Scale":""}`,im&&!im.cal?"warn":"ok"); }
  errs.slice(0,6).forEach(e=>note(e,"warn"));
}
function showFormats(){ modal(`<div class="im-h">읽는 형식</div><table class="im-kv sm">${IO.FORMATS.map(([a,b,c])=>`<tr><td><b>${esc(a)}</b><div class="xs muted">.${esc(b.split(",").join(" ."))}</div></td><td>${esc(c)}</td></tr>`).join("")}</table>
  <div class="xs muted">화소 크기 출처: ImageJ·OME 메타데이터, SEM 장비 태그(Zeiss 34118, FEI/Thermo 34680·34682, TESCAN 50431), Hitachi·JEOL 옆 .txt(같은 이름으로 함께 열기), DM·MRC·EMD 보정, Nikon ND2 dCalibration, Zeiss CZI Scaling, Igor wave scaling, Nanonis SCAN_RANGE. TIFF·PNG 해상도 태그는 기본값인 경우가 많아 '후보'로만 보입니다.</div>
  <div class="im-row"><span class="sp"></span><button class="btn sm" data-close>닫기</button></div>`); }
function modal(html){ const m=$q(ROOT,"#imModal"); m.hidden=false; m.innerHTML=`<div class="box">${html}</div>`; m.onclick=e=>{ if(e.target===m||e.target.dataset.close!=null){ m.hidden=true; m.innerHTML=""; } }; return m; }
/* 큰 현미경 파일(nd2·czi): 위치/scene·Z·T 고르기 — 작은 파일은 묻지 않고 전부 */
let IMPORT_Q=Promise.resolve(); // 고르기 창은 한 번에 하나 (열기 중 또 끌어 놓아도 앞 파일이 사라지지 않게)
function askImport(o){ const r=IMPORT_Q.then(()=>askImport1(o)); IMPORT_Q=r.catch(()=>{}); return r; }
function askImport1(o){ const posB=o.posBytes||[...Array(o.positions)].map(()=>o.bytesPerPlane*o.planesPer*o.nz*o.nt); const total=posB.reduce((a,b)=>a+b,0); const LIM=1.6e9;
  if(o.positions<=6&&total<=6e8&&o.nz*o.nt<=400) return Promise.resolve({}); // 위치 6 개 이하·작은 파일은 묻지 않고 전부
  hideProg();
  return new Promise(res=>{ const gb=b=>b>=1e9?(b/1e9).toFixed(2)+" GB":b>=1e6?(b/1e6).toFixed(0)+" MB":Math.max(1,Math.round(b/1e3))+" KB";
    const rng=(id,n,lab)=>n>1?`<div class="im-row"><b style="width:18px">${lab}</b><label class="im-f"><span>from</span><input type="number" id="${id}a" value="1" min="1" max="${n}" style="width:64px"></label><label class="im-f"><span>to</span><input type="number" id="${id}b" value="${n}" min="1" max="${n}" style="width:64px"></label><label class="im-f"><span>step</span><input type="number" id="${id}s" value="1" min="1" max="${n}" style="width:56px"></label><span class="xs muted">/ ${n}</span></div>`:"";
    const pos=o.positions>1?`<div class="im-h sm">${o.kind==="czi"?"Scene":"XY position"} <span class="muted xs">${o.positions}</span><span class="sp"></span><button class="btn xs" id="ipAll">모두</button><button class="btn xs" id="ipNone">없음</button></div>
      <div class="im-list ip-list">${[...Array(o.positions).keys()].map(i=>`<label class="im-chk"><input type="checkbox" class="ipP" value="${i}" ${i<64?"checked":""}> ${esc(o.posNames[i]||String(i+1))}</label>`).join("")}</div>`:"";
    const m=modal(`<div class="im-h">Import — ${esc(o.name)}</div>
      <div class="xs muted">${o.posBytes&&o.positions>1?"scene 마다 크기 다름":`${o.w}×${o.h} px`} · C ${o.nc} · Z ${o.nz} · T ${o.nt}${o.positions>1?` · ${o.kind==="czi"?"scene":"위치"} ${o.positions}`:""} — 전부 열면 ${gb(total)} (32-bit 메모리)</div>
      ${pos}${rng("ipZ",o.nz,"Z")}${rng("ipT",o.nt,"T")}
      <div class="im-row"><span class="xs" id="ipInfo"></span><span class="sp"></span><button class="btn sm" data-close>취소</button><button class="btn primary sm" id="ipOk">Open</button></div>`);
    const g=q=>$q(m,q); const seq=(id,n)=>{ if(n<=1) return [0]; const a=clamp(Math.round(+g(`#${id}a`).value||1),1,n), b=clamp(Math.round(+g(`#${id}b`).value||n),a,n), st=Math.max(1,Math.round(+g(`#${id}s`).value||1)); const r=[]; for(let i=a;i<=b;i+=st) r.push(i-1); return r; };
    const pick=()=>{ const p=o.positions>1?$a(m,".ipP").filter(x=>x.checked).map(x=>+x.value):[0]; return {p,z:seq("ipZ",o.nz),t:seq("ipT",o.nt)}; };
    const cnt=(list,n)=>n<=1?1:list.filter(v=>v<n).length; /* 위치마다 Z·T 수가 다르면(posNz·posNt) 그 위치 안의 것만 셈 */
    const bytesOf=r=>o.posNz||o.posNt?r.p.reduce((a,i)=>{ const nz=o.posNz?o.posNz[i]:o.nz, nt=o.posNt?o.posNt[i]:o.nt; return a+posB[i]*cnt(r.z,nz)/nz*cnt(r.t,nt)/nt; },0):r.p.reduce((a,i)=>a+posB[i],0)*r.z.length*r.t.length/(o.nz*o.nt);
    const upd=()=>{ const r=pick(); const b=bytesOf(r); const el=g("#ipInfo"); el.textContent=`${r.p.length} × Z ${r.z.length} × T ${r.t.length} → ${gb(b)}${b>LIM?" — 너무 큼(한도 "+gb(LIM)+")":""}`; el.className="xs "+(b>LIM||!r.p.length?"warn":"ok"); g("#ipOk").disabled=b>LIM||!r.p.length; };
    $a(m,"input").forEach(x=>x.oninput=upd); if(g("#ipAll")) g("#ipAll").onclick=()=>{ $a(m,".ipP").forEach(x=>x.checked=true); upd(); }; if(g("#ipNone")) g("#ipNone").onclick=()=>{ $a(m,".ipP").forEach(x=>x.checked=false); upd(); };
    // 처음 값: 한도를 넘으면 Z(없으면 T) step 을 키워 맞춘다
    { const b0=posB.slice(0,Math.min(64,o.positions)).reduce((a,b)=>a+b,0); if(b0>LIM){ const k=o.nz>1?"ipZ":o.nt>1?"ipT":null; if(k){ const n=k==="ipZ"?o.nz:o.nt; g(`#${k}s`).value=Math.min(n,Math.max(1,Math.ceil(b0/LIM))); } } } // 처음 체크된 위치(≤64) 기준으로 한도에 맞는 가장 촘촘한 간격
    upd(); let done=false; const fin=v=>{ if(done) return; done=true; m.hidden=true; m.innerHTML=""; res(v); };
    g("#ipOk").onclick=()=>fin(pick()); m.onclick=e=>{ if(e.target===m||e.target.dataset.close!=null) fin(null); }; }); }
function askRaw(name,size){ return new Promise(res=>{ const m=modal(`<div class="im-h">Import Raw — ${esc(name)} <span class="muted xs">${size.toLocaleString()} bytes</span></div>
  <div class="im-row"><label class="im-f"><span>Type</span><select class="im-sel" id="rwT">${[["u8","8-bit"],["u16","16-bit unsigned"],["i16","16-bit signed"],["u32","32-bit unsigned"],["i32","32-bit signed"],["f32","32-bit float"],["f64","64-bit float"]].map(([a,b])=>`<option value="${a}">${b}</option>`).join("")}</select></label><label class="im-chk"><input type="checkbox" id="rwLE" checked> Little-endian</label></div>
  <div class="im-row"><label class="im-f"><span>Width</span><input type="number" id="rwW" value="512"></label><label class="im-f"><span>Height</span><input type="number" id="rwH" value="512"></label><label class="im-f"><span>Images</span><input type="number" id="rwN" value="1"></label></div>
  <div class="im-row"><label class="im-f"><span>Offset</span><input type="number" id="rwO" value="0"></label><label class="im-f"><span>Gap</span><input type="number" id="rwG" value="0"></label><span class="xs muted" id="rwInfo"></span></div>
  <div class="im-row"><span class="sp"></span><button class="btn sm" data-close>취소</button><button class="btn primary sm" id="rwOk">Open</button></div>`);
  const g=s=>$q(m,s); const info=()=>{ const sz={u8:1,u16:2,i16:2,u32:4,i32:4,f32:4,f64:8}[g("#rwT").value]; const need=+g("#rwO").value+(+g("#rwW").value)*(+g("#rwH").value)*sz*(+g("#rwN").value)+(+g("#rwG").value)*((+g("#rwN").value)-1); g("#rwInfo").textContent=`필요 ${need.toLocaleString()} B${need>size?" — 파일보다 큼":""}`; };
  $a(m,"input,select").forEach(e=>e.oninput=info); info();
  g("#rwOk").onclick=()=>{ const o={type:g("#rwT").value,le:g("#rwLE").checked,w:+g("#rwW").value,h:+g("#rwH").value,n:+g("#rwN").value||1,offset:+g("#rwO").value||0,gap:+g("#rwG").value||0}; m.hidden=true; res(o); };
  m.addEventListener("click",e=>{ if(e.target===m||e.target.dataset.close!=null) res(null); }); }); }

/* ---------- 저장·내보내기 ---------- */
function saveTIFF(){ const im=cur(); if(!im) return; const t=tiffOf(im); downloadBlob(t.name,new Blob([t.data],{type:"image/tiff"})); note(`${t.name} — ImageJ/Fiji 에서 보정·스택 정보와 함께 열립니다`,"ok"); }
function exportPNG(){ const im=cur(); if(!im) return; const c=document.createElement("canvas"); c.width=im.w; c.height=im.h; const g=c.getContext("2d"); g.imageSmoothingEnabled=false;
  if(fullScale(im)===1){ const D=buildDisplay(im); g.drawImage(D.cv,0,0); } else { const hs=Math.max(1,Math.floor(8e6/im.w)); for(let y=0;y<im.h;y+=hs){ const st=renderRegion(im,{x:0,y,w:im.w,h:Math.min(hs,im.h-y),sc:1}); g.drawImage(st,0,y); } } // 큰 영상도 원래 해상도 (띠로 나눠)
  const v={zoom:1,x:0,y:0}; drawThreshold(im,g,v); drawParticles(im,g,v); if(S.ui.showAll!==false) roiList(im).forEach(r=>drawRoi(g,v,r,false,im)); const a=activeRoi(im); if(a) drawRoi(g,v,a,false,im); drawScaleBar(im,g,v,im.w,im.h);
  c.toBlob(b=>downloadBlob(im.name.replace(/\.[^.]+$/,"")+"_view.png",b),"image/png"); }

/* ---------- 단축키 ---------- */
function onKey(e){
  if(!ROOT||!document.body.contains(ROOT)||!$q(ROOT,".imx")) return;
  { const m=$q(ROOT,"#imModal"); if(e.key==="Escape"&&m&&!m.hidden){ e.preventDefault(); const b=$q(m,"[data-close]"); if(b) b.click(); else { m.hidden=true; m.innerHTML=""; } return; } } // Esc = 창 닫기(취소)
  const t=e.target; if(t&&(t.tagName==="INPUT"||t.tagName==="TEXTAREA"||t.tagName==="SELECT"||t.isContentEditable)) return;
  const im=cur(); const k=e.key; const mod=e.ctrlKey||e.metaKey;
  if(k===" "&&!SPACE){ SPACE=true; e.preventDefault(); return; }
  if(mod&&!e.shiftKey&&k.toLowerCase()==="z"){ e.preventDefault(); undo(); return; } // 다시 하기(redo)는 없음
  if(!im) return;
  if(mod&&k.toLowerCase()==="m"){ e.preventDefault(); runOp("measure",{all:false}); return; }
  if(mod&&e.shiftKey&&k.toLowerCase()==="a"||k==="Escape"){ if(DRAG&&(DRAG.mode==="poly"||DRAG.mode==="angle")){ DRAG=null; draw(); return; } e.preventDefault(); setActive(im,null); return; }
  if(mod&&e.shiftKey&&k.toLowerCase()==="d"){ e.preventDefault(); runOp("duplicate",{slice:false}); return; }
  if(!mod&&k.toLowerCase()==="t"){ roiAdd(im); return; }
  if(k==="+"||k==="="){ zoomStep(im,1); return; } if(k==="-"||k==="_"){ zoomStep(im,-1); return; }
  const a=activeRoi(im); if(a&&k.startsWith("Arrow")){ e.preventDefault(); const d=e.shiftKey?10:1; const [dx,dy]={ArrowLeft:[-d,0],ArrowRight:[d,0],ArrowUp:[0,-d],ArrowDown:[0,d]}[k]; a.pts=a.pts.map(p=>[p[0]+dx,p[1]+dy]); save(); draw(); return; }
  if(isStack(im)&&(k===","||k===".")){ const key=im.nz>1?"z":"t", n=im.nz>1?im.nz:im.nt; im[key]=(im[key]+(k==="."?1:n-1))%n; invalidate(); drawStack(); draw(); }
}
try{ document.addEventListener("keydown",onKey); document.addEventListener("keyup",e=>{ if(e.key===" ") SPACE=false; }); window.addEventListener("blur",()=>{ SPACE=false; }); }catch(e){}

window.IMG_UI={render,rerender,_load:loadFiles,_state:()=>S,_px:PX,_cur:cur,_run:runOp,_ops:OPS,_steps:()=>S.macro.steps,_runSteps:runSteps,_addImage:addImage,_setActive:(r)=>setActive(cur(),r),_flush:flushSave,_toIJM:()=>toIJM(S.macro.steps),_busy:()=>BUSY,_undo:undo,_close:closeImage,_cancel:cancelAll,_msgs:()=>MSGS,_wand:wandAt};
