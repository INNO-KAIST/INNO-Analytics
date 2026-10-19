/* ============================================================
   Macro: 자동 기록 → 단계 편집 → 다른 이미지·파일에 batch 적용 → ImageJ .ijm 로 내보내기
   ============================================================ */
let REPLAY=false;
function record(step,im){ if(!S.macro.rec||REPLAY||PVMODE) return; S.macro.steps.push(Object.assign({},step,{img:im?im.id:null,t:Date.now()})); save(); if(S.ui.panel==="macro") renderPanel(); }
async function runSteps(steps,im,onStep){
  let target=im; REPLAY=true; const gen=CANCEL_GEN;
  try{ for(let i=0;i<steps.length;i++){ if(gen!==CANCEL_GEN) throw new Error("취소"); const s=steps[i]; const o=OPS[s.op]; if(!o) throw new Error("알 수 없는 단계: "+s.op); onStep&&onStep(i,s);
      const P=Object.assign({},s.params,{__op:s.op}); if(s.op==="calc"&&P.bName){ const b=S.images.find(x=>x.name===P.bName); if(b) P.b=b.id; }
      if(s.op==="analyze"&&P.redirectName&&P.redirect&&P.redirect!=="__parent"&&P.redirect!=="__root"){ const b=S.images.find(x=>x.name===P.redirectName); if(!b) throw new Error(`Redirect 이미지 '${P.redirectName}' 이 열려 있지 않습니다`); P.redirect=b.id; }
      S.sel=target.id; const r=await runWithUnit(o,target,P); if(r&&r.id) target=r; }
  } finally { REPLAY=false; }
  return target;
}
const stepText=s=>s.label||(OPS[s.op]?OPS[s.op].label:s.op);
function panelMacro(el,im){
  const M=S.macro; const ims=S.images;
  el.innerHTML=`<div class="im-h">Recorder <span class="sp"></span><label class="im-chk"><input type="checkbox" id="mcRec" ${M.rec?"checked":""}> Recording</label></div>
  <div class="xs muted">적용한 작업이 순서대로 쌓입니다(ImageJ Recorder 와 같은 역할). 필요 없는 단계를 지우고 다른 이미지·파일에 그대로 적용하세요.</div>
  <ol class="im-steps">${M.steps.length?M.steps.map((s,i)=>`<li><span class="tx">${esc(stepText(s))}</span><span class="b"><button class="x" data-up="${i}" title="Move up">↑</button><button class="x" data-dn="${i}" title="Move down">↓</button><button class="x" data-del="${i}" title="Delete">×</button></span></li>`).join(""):`<li class="xs muted">아직 기록된 단계가 없습니다.</li>`}</ol>
  <div class="im-row"><button class="btn sm" id="mcClr" ${M.steps.length?"":"disabled"}>Clear</button><button class="btn sm" id="mcSave" ${M.steps.length?"":"disabled"}>Save (.json)</button><label class="btn sm" style="cursor:pointer">Open<input type="file" id="mcOpen" accept=".json" hidden></label><button class="btn sm" id="mcIjm" ${M.steps.length?"":"disabled"} title="Re-run the same processing in ImageJ/Fiji">Export .ijm</button></div>
  <hr class="im-hr"><div class="im-h">Run</div>
  <div class="im-row"><button class="btn primary sm" id="mcRun" ${M.steps.length&&im?"":"disabled"}>Run on current image</button><span class="xs muted">${im?esc(im.name):"no image"}</span></div>
  <div class="im-h">Batch <span class="muted xs">열린 이미지 또는 파일 여러 개</span></div>
  <div class="im-batch">${ims.length?ims.map(x=>`<label class="im-chk"><input type="checkbox" data-bt="${x.id}" ${BATCH.sel.has(x.id)?"checked":""}> ${esc(x.name)} <span class="xs muted">${x.w}×${x.h}</span></label>`).join(""):`<div class="xs muted">No open images</div>`}</div>
  <div class="im-row"><button class="btn sm" id="mcAll">Select all</button><button class="btn primary sm" id="mcBatch" ${M.steps.length?"":"disabled"}>Run on selected images</button></div>
  <div class="im-row"><label class="btn sm" style="cursor:pointer" title="Opens each file, applies the steps and keeps only the result tables">Choose files and run batch<input type="file" id="mcFiles" multiple accept="${IO.ACCEPT}" hidden></label><label class="im-chk"><input type="checkbox" id="mcKeep" ${BATCH.keep?"checked":""}> Keep processed images in the list</label><label class="im-chk"><input type="checkbox" id="mcZip" ${BATCH.zip?"checked":""}> Save result images as TIFF (zip)</label></div>
  <div class="xs muted">batch 결과는 Summary·Results 표에 이미지 이름과 함께 쌓입니다. 문턱은 기록한 <b>방법</b>으로 이미지마다 다시 구합니다(수동 범위였다면 그 값 그대로).</div>`;
  const q=s=>$q(el,s);
  q("#mcRec").onchange=e=>{ M.rec=e.target.checked; save(); };
  $a(el,"[data-del]").forEach(b=>b.onclick=()=>{ M.steps.splice(+b.dataset.del,1); save(); renderPanel(); });
  $a(el,"[data-up]").forEach(b=>b.onclick=()=>{ const i=+b.dataset.up; if(i>0){ const t=M.steps[i]; M.steps[i]=M.steps[i-1]; M.steps[i-1]=t; save(); renderPanel(); } });
  $a(el,"[data-dn]").forEach(b=>b.onclick=()=>{ const i=+b.dataset.dn; if(i<M.steps.length-1){ const t=M.steps[i]; M.steps[i]=M.steps[i+1]; M.steps[i+1]=t; save(); renderPanel(); } });
  q("#mcClr").onclick=()=>{ if(confirm2("mcclr","매크로 단계를 모두 지울까요? 한 번 더")){ M.steps=[]; save(); renderPanel(); } };
  q("#mcSave").onclick=()=>downloadText("INNO_image_macro.json",JSON.stringify({kind:"INNO Analytics image macro",v:1,steps:M.steps.map(s=>({op:s.op,params:s.params,label:s.label}))},null,1),"application/json");
  q("#mcOpen").onchange=e=>{ const f=e.target.files[0]; if(!f) return; f.text().then(t=>{ try{ const o=JSON.parse(t); if(!Array.isArray(o.steps)) throw 0; const bad=o.steps.filter(s=>!OPS[s.op]); if(bad.length) throw new Error("모르는 단계: "+bad.map(s=>s.op).join(", ")); M.steps=o.steps; save(); renderPanel(); note(`매크로 ${o.steps.length} 단계를 읽었습니다`,"ok"); }catch(err){ note("매크로 파일을 읽지 못했습니다"+(err&&err.message?": "+err.message:""),"warn"); } }); e.target.value=""; };
  q("#mcIjm").onclick=()=>downloadText("INNO_image_macro.ijm",toIJM(M.steps),"text/plain");
  q("#mcRun").onclick=async()=>{ if(BUSY) return; BUSY="macro"; try{ const t=await runSteps(M.steps,im,(i,s)=>setProg(`매크로 ${i+1}/${M.steps.length}: ${stepText(s)}`,i/M.steps.length)); hideProg(); S.sel=t.id; note("매크로 실행 완료","ok"); }catch(e){ hideProg(); note("매크로 중단: "+e.message,"warn"); } BUSY=null; save(); rerender(); };
  $a(el,"[data-bt]").forEach(cb=>cb.onchange=()=>{ cb.checked?BATCH.sel.add(cb.dataset.bt):BATCH.sel.delete(cb.dataset.bt); });
  q("#mcAll").onclick=()=>{ ims.forEach(x=>BATCH.sel.add(x.id)); renderPanel(); };
  q("#mcKeep").onchange=e=>{ BATCH.keep=e.target.checked; }; q("#mcZip").onchange=e=>{ BATCH.zip=e.target.checked; };
  q("#mcBatch").onclick=()=>batchImages([...BATCH.sel].map(byId).filter(Boolean));
  q("#mcFiles").onchange=e=>{ const fl=[...e.target.files]; e.target.value=""; batchFiles(fl); };
}
const BATCH={sel:new Set(),keep:false,zip:false};
async function batchImages(list){ if(!list.length){ note("이미지를 고르세요","warn"); return; } if(BUSY) return; BUSY="batch"; const outs=[]; let ok=0; const gen=CANCEL_GEN;
  try{ for(let k=0;k<list.length;k++){ if(gen!==CANCEL_GEN){ note(`Batch 중단 — ${k}/${list.length} 에서`,"warn"); break; } const im=list[k]; setProg(`Batch ${k+1}/${list.length}: ${im.name}`,k/list.length); try{ const t=await runSteps(S.macro.steps,im); outs.push(t); ok++; }catch(e){ note(`${im.name}: ${e.message}`,"warn"); } } }
  finally{ hideProg(); BUSY=null; }
  if(BATCH.zip) zipImages(outs); note(`Batch 완료 — ${ok}/${list.length} 이미지`,"ok"); S.ui.bottom="summary"; S.ui.bottomOpen=true; save(); rerender(); }
async function batchFiles(files){ if(!files.length||BUSY) return; BUSY="batch"; let ok=0, nImg=0; const zipList=[]; const gen=CANCEL_GEN;
  // 옆 파일(.txt: Hitachi·JEOL 보정)은 이미지로 열지 않고 같은 이름 이미지에 넘김 (열기와 같은 규칙)
  const base=n=>n.replace(/\.[^.]+$/,"").toLowerCase(); const imgNames=new Set(files.filter(f=>!/\.txt$/i.test(f.name)).map(f=>base(f.name)));
  const side=files.filter(f=>/\.txt$/i.test(f.name)&&imgNames.has(base(f.name))); const sidecars={}; for(const f of side) sidecars[f.name.toLowerCase()]=IO._txt(new Uint8Array(await f.arrayBuffer()));
  const list=files.filter(f=>!side.includes(f));
  try{ for(let k=0;k<list.length;k++){ if(gen!==CANCEL_GEN){ note(`Batch 중단 — ${k}/${list.length} 파일에서`,"warn"); break; } const f=list[k]; setProg(`Batch ${k+1}/${list.length}: ${f.name}`,k/list.length);
      const before=new Set(S.images.map(x=>x.id));
      try{ const u8=new Uint8Array(await f.arrayBuffer()); const r=await IO.readFile(f.name,u8,{sidecars});
        for(const d of r.images){ if(gen!==CANCEL_GEN) break; nImg++; try{ const im=addImage(d,{file:f.name}); const t=await runSteps(S.macro.steps,im); if(BATCH.zip) zipList.push(tiffOf(t)); ok++; }catch(e){ if(String(e.message)==="취소") break; note(`${d.name}: ${e.message}`,"warn"); } } }
      catch(e){ note(`${f.name}: ${e.message}`,"warn"); }
      finally{ if(!BATCH.keep){ S.images.filter(x=>!before.has(x.id)).map(x=>x.id).forEach(id=>closeImage(id,true)); } } // 실패해도 이 파일에서 생긴 이미지는 닫음
      await sleep(0); } }
  finally{ hideProg(); BUSY=null; }
  if(BATCH.zip&&zipList.length) downloadBlob("batch_results.zip",makeZip(zipList)); note(`Batch 완료 — ${ok}/${nImg} 이미지 (${list.length} 파일${side.length?`, 옆 .txt ${side.length}`:""}). 결과는 Summary 표`,ok===nImg?"ok":"warn"); S.ui.bottom="summary"; S.ui.bottomOpen=true; save(); rerender(); }
function tiffOf(im){ return {name:im.name.replace(/\.[^.]+$/,"")+".tif",data:IO.writeTIFF({w:im.w,h:im.h,bit:im.bit,rgb:im.rgb,nc:im.rgb?1:im.nc,nz:im.nz,nt:im.nt,planes:planesOf(im),cal:im.cal,zStep:im.zStep},im.rgb?{}:{min:im.disp[0].min,max:im.disp[0].max})}; }
function zipImages(list){ if(!list.length) return; downloadBlob("batch_results.zip",makeZip(list.map(tiffOf))); }
/* ZIP (무압축 store) */
const CRC=(()=>{ const t=new Uint32Array(256); for(let n=0;n<256;n++){ let c=n; for(let k=0;k<8;k++) c=c&1?0xEDB88320^(c>>>1):c>>>1; t[n]=c>>>0; } return t; })();
function crc32(u8){ let c=0xFFFFFFFF; for(let i=0;i<u8.length;i++) c=CRC[(c^u8[i])&255]^(c>>>8); return (c^0xFFFFFFFF)>>>0; }
function makeZip(files){ const enc=new TextEncoder(); const parts=[], cent=[]; let off=0; const seen={};
  for(const f of files){ let nm=f.name; if(seen[nm]){ seen[nm]++; nm=nm.replace(/(\.[^.]+)?$/,`_${seen[nm]}$1`); } else seen[nm]=1; const n=enc.encode(nm), d=f.data, crc=crc32(d);
    const h=new DataView(new ArrayBuffer(30)); h.setUint32(0,0x04034b50,true); h.setUint16(4,20,true); h.setUint16(6,0x0800,true); h.setUint32(14,crc,true); h.setUint32(18,d.length,true); h.setUint32(22,d.length,true); h.setUint16(26,n.length,true);
    parts.push(new Uint8Array(h.buffer),n,d);
    const c=new DataView(new ArrayBuffer(46)); c.setUint32(0,0x02014b50,true); c.setUint16(4,20,true); c.setUint16(6,20,true); c.setUint16(8,0x0800,true); c.setUint32(16,crc,true); c.setUint32(20,d.length,true); c.setUint32(24,d.length,true); c.setUint16(28,n.length,true); c.setUint32(42,off,true);
    cent.push(new Uint8Array(c.buffer),n); off+=30+n.length+d.length; }
  const csz=cent.reduce((a,b)=>a+b.length,0); const e=new DataView(new ArrayBuffer(22)); e.setUint32(0,0x06054b50,true); e.setUint16(8,files.length,true); e.setUint16(10,files.length,true); e.setUint32(12,csz,true); e.setUint32(16,off,true);
  return new Blob([...parts,...cent,new Uint8Array(e.buffer)],{type:"application/zip"}); }
/* ImageJ 매크로로 */
function toIJM(steps){
  const body=steps.map(s=>{ const o=OPS[s.op]; let line=""; try{ line=o&&o.ij?o.ij(Object.assign({},s.params)):""; }catch(e){ line=""; } return `// ${stepText(s)}\n${line||"// (ImageJ 에 같은 명령 없음)"}`; }).join("\n");
  return `// INNO Analytics → ImageJ macro (${new Date().toISOString().slice(0,10)})\n// ImageJ 1.54 기준 명령. 한 이미지에: 이미지를 연 뒤 아래 '처리' 부분만 실행.\n// 폴더 전체: 이 파일 그대로 실행 → 폴더 고르기.\nsetOption("BlackBackground", true);\n`+
`dir = getDirectory("이미지 폴더");\nlist = getFileList(dir);\nsetBatchMode(true);\nfor (i=0; i<list.length; i++) {\n  if (File.isDirectory(dir+list[i])) continue;\n  if (!matches(toLowerCase(list[i]), ".*\\\\.(tif|tiff|png|jpg|jpeg|bmp|gif|dm3|dm4|nd2|czi|lif|oib|oir|ome\\\\.tif)$")) continue; // 이미지 파일만\n  open(dir+list[i]);\n  // ---- 처리 ----\n${body.split("\n").map(l=>"  "+l).join("\n")}\n  // -------------\n  run("Close All");\n}\nsetBatchMode(false);\n`;
}
