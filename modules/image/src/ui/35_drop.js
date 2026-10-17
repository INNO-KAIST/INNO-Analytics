/* ============================================================
   Drop 패널 — 접촉각(sessile drop) · 펜던트 드랍(pendant drop) 표면장력
   계산: core/88_drop.js (IMG_CORE.contactAngle / pendantDrop, Worker)
   영역 = 사각형 ROI (없으면 영상 전체) · 기준선 = 자동(기판 윗면 RANSAC) 또는 선 ROI 로 지정
   결과: S.drops (아래 'Drop' 표), 덧그림 S.dropFit[영상 id]
   ============================================================ */
const DROP_M=[["yl","Young–Laplace"],["poly","Tangent (polynomial)"],["ellipse","Ellipse"],["circle","Circle (spherical cap)"]];
const DROP_COL={yl:"#ff3ccf",poly:"#2bdc4a",ellipse:"#ffa000",circle:"#4ea6ff"};
const MM_OF={m:1e3,meter:1e3,cm:10,mm:1,"µm":1e-3,"μm":1e-3,um:1e-3,micron:1e-3,microns:1e-3,nm:1e-6,"Å":1e-7,A:1e-7,pm:1e-9};
const mmPerPxOf=im=>{ const c=calOf(im); const f=MM_OF[String(c.unit).trim()]; return im.cal&&f?c.pw*f:null; };
const calAniso=im=>{ const c=calOf(im); return im.cal&&Math.abs(c.ph/c.pw-1)>0.005; };
const dropUi=()=>S.ui.drop||(S.ui.drop=defState().ui.drop);
/* 기판·용액 밀도 기본값: 물 20 °C 998.2, 공기 1.204 kg/m³ (CRC Handbook) — 이 값은 바꿔 쓰는 자리 */
function dropFrames(im,P){ if(P.frames!=="all"||!isStack(im)) return [{z:im.z,t:im.t,label:""}]; const L=[];
  if(im.nt>1){ for(let t=0;t<im.nt;t++) L.push({z:im.z,t,label:`t${t+1}`}); } else for(let z=0;z<im.nz;z++) L.push({z,t:im.t,label:`${z+1}`}); return L; }
function dropRegion(im,P){ if(P.roi&&P.roi.length===4){ const [x0,y0,x1,y1]=P.roi; return {x0,y0,x1,y1}; } return null; }

op("drop",{label:"Drop Shape",group:"meas",hidden:true,kind:"measure",params:[],
  run:async(im,P)=>{ /* 색 영상은 밝기(0.299 R + 0.587 G + 0.114 B)로 */
    const fr=dropFrames(im,P); const roi=dropRegion(im,P); const mmpp=P.pxMM>0?+P.pxMM:mmPerPxOf(im); const rows=[]; let last=null, lastF=null, fails=0; const runId=uid(); const pre=[];
    if(P.imgW&&(P.imgW!==im.w||P.imgH!==im.h)&&(P.roi||P.baseline)) pre.push(`기록한 영역·기준선은 ${P.imgW}×${P.imgH} 영상의 좌표입니다 — 이 영상(${im.w}×${im.h})에서는 자리가 맞지 않을 수 있습니다`);
    if(calAniso(im)) pre.push("화소가 가로·세로 크기가 달라(비등방) 형상 맞춤이 맞지 않습니다 — 먼저 Scale 로 정사각 화소로");
    for(let k=0;k<fr.length;k++){ const f=fr[k]; if(fr.length>1) setProg(`Drop ${k+1}/${fr.length}`,k/fr.length);
      const pl=im.rgb?lumaPlane(im,f.z,f.t):planeAt(im,im.c,f.z,f.t);
      let r; try{
        if(P.mode==="pendant") r=await call("pendantDrop",[pl,{roi,needleMM:+P.needleMM>0?+P.needleMM:null,pxSize:mmpp,rhoIn:+P.rhoIn,rhoOut:+P.rhoOut,g:+P.g||9.80665}]);
        else r=await call("contactAngle",[pl,{roi,baseline:P.baseline?{p1:P.baseline[0],p2:P.baseline[1]}:null,methods:P.methods,pxSize:mmpp,rhoIn:+P.rhoIn,rhoOut:+P.rhoOut,g:+P.g||9.80665}]);
      }catch(e){ fails++; rows.push({run:runId,Image:im.name,Frame:f.label,Mode:P.mode==="pendant"?"Pendant":"Sessile",Note:"실패: "+(e.message||e)}); if(fr.length===1) throw e; continue; }
      if(pre.length) r.warn=pre.concat(r.warn||[]); const rr=dropRows(im,f.label,P,r,mmpp); rr.forEach(x=>{ x.run=runId; x._k=k+1; }); rows.push(...rr);
      if(f.z===im.z&&f.t===im.t||!last){ last=r; lastF=f; } }
    hideProg(); if(!S.drops) S.drops=[]; S.drops.push(...rows);
    if(last) S.dropFit[im.id]={run:runId,mode:P.mode,ov:last.overlay,contour:last.contour,edge:last.edge,warn:last.warn||[],z:lastF.z,t:lastF.t,ver:im.ver,frame:lastF.label}; /* 결과를 만든 장 — 그 장에서만 덧그림 */
    if(fr.length>1){ const lab=P.mode==="pendant"?"Pendant":DROP_M.find(x=>x[0]===dropMainMethod(P))[1]; const ok=rows.filter(r=>r.Method===lab&&(P.mode==="pendant"?r["γ"]!=null:r["θ mean"]!=null));
      const y=P.mode==="pendant"?ok.map(r=>r["γ"]):ok.map(r=>r["θ mean"]); if(ok.length>1) S.plots.push({id:uid(),kind:"xy",title:`${P.mode==="pendant"?"Surface tension":"Contact angle ("+lab+")"} — ${im.name}`,xl:"frame",yl:P.mode==="pendant"?"γ (mN/m)":"θ (°)",series:[{name:P.mode==="pendant"?"γ":"θ",x:ok.map(r=>r._k),y}]}); }
    record({op:"drop",params:Object.assign(clean(P),{imgW:im.w,imgH:im.h}),label:P.mode==="pendant"?"Pendant drop":"Contact angle"},im);
    S.ui.bottom="drop"; S.ui.bottomOpen=true; save();
    const w=last&&last.warn&&last.warn.length?" · "+last.warn[0]:"";
    if(P.mode==="pendant"){ const r0=last; note(r0&&r0.gamma!=null?`γ = ${r0.gamma.toFixed(2)} mN/m (맞춤 ± ${r0.gammaSd.toFixed(2)}${r0.gammaScaleSd===r0.gammaScaleSd?`, 바늘 눈금 ± ${r0.gammaScaleSd.toFixed(2)}`:""}) · Bo ${r0.Bo.toFixed(4)} · Wo ${r0.Wo.toFixed(3)}${fr.length>1?` (${fr.length-fails}/${fr.length} 장)`:""}${w}`:`Bo = ${r0?r0.Bo.toFixed(4):"—"} — 바늘 지름(mm)을 넣으면 γ 를 계산합니다${w}`,w?"warn":"ok"); }
    else { const m=last&&last.methods[dropMainMethod(P)]; note(m?`θ = ${m.left.toFixed(1)}° / ${m.right.toFixed(1)}° (${DROP_M.find(x=>x[0]===dropMainMethod(P))[1]})${fr.length>1?` (${fr.length-fails}/${fr.length} 장)`:""}${w}`:"접촉각을 계산하지 못했습니다",w||!m?"warn":"ok"); } },
  ij:P=>`// Drop Shape (${P.mode==="pendant"?"pendant drop":"contact angle"}): ImageJ 기본 명령 없음 — DropSnake / LB-ADSA / Pendent_Drop 플러그인에 해당`,
  help:"접촉각(Young–Laplace · 접선 · 타원 · 원)과 펜던트 드랍 표면장력(Young–Laplace, 바늘 지름으로 눈금)."});
function dropMainMethod(P){ const m=P.methods||[]; return ["yl","poly","ellipse","circle"].find(k=>m.includes(k))||"yl"; }
function lumaPlane(im,z,t){ const R=planeAt(im,0,z,t).data, G=planeAt(im,1,z,t).data, B=planeAt(im,2,z,t).data; const o=new Float32Array(R.length); for(let i=0;i<o.length;i++) o[i]=0.299*R[i]+0.587*G[i]+0.114*B[i]; return {w:im.w,h:im.h,bit:32,data:o}; }
function dropRows(im,frame,P,r,mmpp){ const base={Image:im.name,Frame:frame}; const W=(r.warn||[]).join(" / ");
  if(P.mode==="pendant") return [Object.assign(base,{Mode:"Pendant",Method:"Pendant","γ":r.gamma,"γ SD":r.gammaSd,"γ scale SD":r.gammaScaleSd,Bo:r.Bo,Wo:r.Wo,"V (mm³)":r.V,"A (mm²)":r.A,"R0 (mm)":r.R0mm,"Needle (px)":r.needle.diam,"Needle (mm)":r.needleMM,"Tilt (°)":r.fit.omega,"rms (px)":r.rms,Scale:r.scaleFrom==="needle"?"needle":r.scaleFrom==="image"?"image calibration":"none","Δρ":r.drho,g:r.g,Note:W})];
  const out=[]; for(const [k,l] of DROP_M){ const m=r.methods[k]; if(!m) continue; const row=Object.assign({},base,{Mode:"Sessile",Method:l,"θ left":m.left,"θ right":m.right,"θ mean":m.mean,"θ SD":m.thetaSd,Baseline:r.baseline.auto?`auto ${r.baseline.angleDeg.toFixed(2)}°`:`line ${r.baseline.angleDeg.toFixed(2)}°`,Note:W});
    if(k==="yl"){ Object.assign(row,{Bo:m.Bo,"rms (px)":m.rms}); if(mmpp){ Object.assign(row,{"V (mm³)":r.V,"Base Ø (mm)":r.baseDiameterMM,"Height (mm)":r.heightMM,"R0 (mm)":r.R0mm}); if(r.gammaEst!=null){ row["γ est."]=r.gammaEst; row["γ est. SD"]=r.gammaEstSd; } } row["Δρ"]=+P.rhoIn-(+P.rhoOut); row.g=+P.g||9.80665; }
    if(k==="circle"&&mmpp&&r.Vcap!=null) row["V (mm³)"]=r.Vcap; out.push(row); } return out; }

/* ---------- 패널 ---------- */
function panelDrop(el,im){ const D=dropUi(); const a=activeRoi(im); const roiR=a&&a.type==="rect"?roiBounds(a):null; const line=a&&a.type==="line"?a:null; const mmpp=mmPerPxOf(im); const fit=S.dropFit[im.id]; const B=(S.dropBase||{})[im.id];
  const pend=D.mode==="pendant";
  el.innerHTML=`<div class="im-h">Drop Shape <span class="muted xs">ImageJ DropSnake · LB-ADSA · Pendent Drop 대응</span></div>
  <div class="im-row"><div class="seg" id="dpMode"><button data-m="sessile" class="${pend?"":"active"}">Contact angle</button><button data-m="pendant" class="${pend?"active":""}">Pendant drop</button></div></div>
  <div class="xs muted">${pend?"방울이 바늘 아래에 매달린 영상. 영역(사각형 ROI)의 <b>위쪽 테두리가 바늘을 가로지르게</b> 잡으세요 — 바늘 폭으로 눈금을 정합니다.":"기판 위 방울의 옆모습. 영역(사각형 ROI)에 방울과 그 아래 기판 윗면이 들어가게 잡으세요. 거울 같은 기판(반사상)은 기준선을 직접 그으세요."}</div>
  <div class="im-row"><span class="xs">Region: <b>${roiR?`rect ${Math.round(roiR.w)}×${Math.round(roiR.h)} @ (${Math.round(roiR.x)}, ${Math.round(roiR.y)})`:"whole image"}</b></span><span class="xs muted">${roiR?"":"— 사각형 ROI 를 그리면 그 안만"}</span></div>
  ${pend?`<div class="im-row"><label class="im-f"><span>Needle Ø</span><input type="number" id="dpNd" value="${esc(D.needleMM)}" step="any" min="0" placeholder="mm"><span class="u">mm</span></label><span class="xs muted">${mmpp?`없으면 영상 눈금 ${g4(mmpp*1000)} µm/px 사용`:"바늘 바깥 지름 (예: 22 G = 0.718 mm)"}</span></div>`
  :`<div class="im-row"><span class="xs">Baseline: <b>${B?`line (${B[0].map(v=>v.toFixed(0)).join(", ")}) → (${B[1].map(v=>v.toFixed(0)).join(", ")})`:"auto (substrate edge)"}</b></span><span class="sp"></span><button class="btn xs" id="dpBl" ${line?"":"disabled"} title="Draw a straight line along the substrate surface first">Use line ROI</button><button class="btn xs" id="dpBlA" ${B?"":"disabled"}>Auto</button></div>
  <div class="im-row">${DROP_M.map(([k,l])=>`<label class="im-chk"><input type="checkbox" data-dm="${k}" ${D.methods.includes(k)?"checked":""}> <i class="im-dot" style="background:${DROP_COL[k]}"></i>${esc(l)}</label>`).join("")}</div>`}
  <details class="im-det" ${D.physOpen?"open":""} id="dpPhys"><summary class="xs">Densities · gravity ${pend?"":"(volume · Bond 수 용)"}</summary>
   <div class="im-row"><label class="im-f"><span>ρ drop</span><input type="number" id="dpRi" value="${esc(D.rhoIn)}" step="any"><span class="u">kg/m³</span></label><label class="im-f"><span>ρ outside</span><input type="number" id="dpRo" value="${esc(D.rhoOut)}" step="any"><span class="u">kg/m³</span></label></div>
   <div class="im-row"><label class="im-f"><span>g</span><input type="number" id="dpG" value="${esc(D.g)}" step="any"><span class="u">m/s²</span></label><span class="xs muted">물 20 °C 998.2 · 공기 1.204 · 표준 중력 9.80665</span></div></details>
  ${isStack(im)?`<div class="im-row"><label class="im-f"><span>Frames</span><select class="im-sel" id="dpFr"><option value="current" ${D.frames!=="all"?"selected":""}>Current</option><option value="all" ${D.frames==="all"?"selected":""}>All ${im.nt>1?im.nt:im.nz} (time series)</option></select></label></div>`:""}
  ${!mmpp&&!pend?`<div class="xs muted">영상 눈금이 없어 부피·지름은 빼고 각만 계산합니다 (Adjust ▸ Set Scale).</div>`:""}
  <div class="im-row"><span class="sp"></span><button class="btn sm" id="dpClr" ${fit?"":"disabled"}>Clear overlay</button><button class="btn primary sm" id="dpRun">${pend?"Measure γ":"Measure θ"}</button></div>
  ${fit?dropCard(im,fit):""}
  <div class="xs muted">${pend?"Young–Laplace 맞춤(Rotenberg 1983 · Berry 2015). γ = Δρ g R0² / Bo. Worthington 수 Wo = Δρ g V / (π γ D) 가 1 에 가까울수록 정밀.":"Young–Laplace(축대칭, 중력 포함)가 기준. 원은 Bo ≈ 0(작은 방울)에서만 맞고, 접선법은 접촉점 둘레만 봅니다."}</div>`;
  const q=s=>$q(el,s); const up=()=>save();
  $a(el,"#dpMode button").forEach(b=>b.onclick=()=>{ D.mode=b.dataset.m; up(); renderPanel(); draw(); });
  if(q("#dpNd")) q("#dpNd").onchange=e=>{ D.needleMM=e.target.value===""?"":+e.target.value; up(); };
  q("#dpRi").onchange=e=>{ D.rhoIn=+e.target.value; up(); }; q("#dpRo").onchange=e=>{ D.rhoOut=+e.target.value; up(); }; q("#dpG").onchange=e=>{ D.g=+e.target.value||9.80665; up(); };
  q("#dpPhys").ontoggle=e=>{ D.physOpen=e.target.open; up(); };
  if(q("#dpFr")) q("#dpFr").onchange=e=>{ D.frames=e.target.value; up(); };
  $a(el,"[data-dm]").forEach(cb=>cb.onchange=()=>{ const s=new Set(D.methods); cb.checked?s.add(cb.dataset.dm):s.delete(cb.dataset.dm); D.methods=DROP_M.map(x=>x[0]).filter(k=>s.has(k)); up(); draw(); });
  if(q("#dpBl")) q("#dpBl").onclick=()=>{ const ln=activeRoi(im); if(!ln||ln.type!=="line"){ note("먼저 기판 윗면을 따라 직선(Line) ROI 를 그으세요","warn"); return; } const [p,r]=ln.pts; const L=p[0]<=r[0]?[p.slice(),r.slice()]:[r.slice(),p.slice()]; (S.dropBase||(S.dropBase={}))[im.id]=L; up(); renderPanel(); draw(); note("기준선을 정했습니다 — 선 ROI 를 지워도 남습니다","ok"); };
  if(q("#dpBlA")) q("#dpBlA").onclick=()=>{ if(S.dropBase) delete S.dropBase[im.id]; up(); renderPanel(); draw(); };
  q("#dpClr").onclick=()=>{ delete S.dropFit[im.id]; up(); renderPanel(); draw(); };
  q("#dpRun").onclick=()=>{ if(!pend&&!D.methods.length){ note("방법을 하나 이상 고르세요","warn"); return; }
    const a2=activeRoi(im), rr=a2&&a2.type==="rect"?roiBounds(a2):null; // 누른 때의 ROI (그린 뒤 옮겼어도)
    const P={mode:D.mode,methods:D.methods.slice(),needleMM:D.needleMM,rhoIn:D.rhoIn,rhoOut:D.rhoOut,g:D.g,frames:D.frames,roi:rr?[rr.x,rr.y,rr.x+rr.w-1,rr.y+rr.h-1]:null,baseline:!pend&&B?B.map(p=>p.slice()):null};
    runOp("drop",P); };
}
function dropCard(im,F){ const rows=(S.drops||[]).filter(r=>r.run===F.run&&(F.frame==null||r.Frame===F.frame)); const ov=F.ov||{};
  const pend=F.mode==="pendant"; const last=pend?rows.filter(r=>r.Mode==="Pendant").slice(-1)[0]:null; const ses=!pend?rows.filter(r=>r.Mode==="Sessile"):[];
  const stale=F.ver!==im.ver||F.z!==im.z||F.t!==im.t; void stale;
  return `<div class="im-card">${F.ver!==im.ver?`<div class="note warn xs">이 결과 뒤에 화소가 바뀌었습니다(필터 등) — 다시 Measure 하세요.</div>`:F.z!==im.z||F.t!==im.t?`<div class="xs muted">결과는 ${esc(F.frame||"다른")} 장에서 — 그 장으로 가면 덧그림이 보입니다.</div>`:""}
  ${pend&&last?`<div class="im-stats"><div><span>γ</span><b>${last["γ"]!=null?`${fmt(last["γ"],2)} mN/m <span class="xs muted">± ${fmt(last["γ SD"],2)} fit${last["γ scale SD"]===last["γ scale SD"]&&last["γ scale SD"]!=null?` · ± ${fmt(last["γ scale SD"],2)} needle ±0.3 px`:""}</span>`:"— (scale?)"}</b></div><div><span>Bond number</span><b>${fmt(last.Bo,4)}</b></div><div><span>Worthington</span><b>${fmt(last.Wo,3)}</b></div>
    <div><span>Volume</span><b>${fmt(last["V (mm³)"],3)} mm³</b></div><div><span>Surface area</span><b>${fmt(last["A (mm²)"],3)} mm²</b></div><div><span>R0 (apex)</span><b>${fmt(last["R0 (mm)"],4)} mm</b></div><div><span>Tilt</span><b>${fmt(last["Tilt (°)"],2)}°</b></div><div><span>Fit rms</span><b>${fmt(last["rms (px)"],3)} px</b></div></div>`:""}
  ${!pend&&ses.length?`<table class="im-mini"><tr><th></th><th>Left</th><th>Right</th><th>Mean</th></tr>${ses.map(r=>{ const k=(DROP_M.find(x=>x[1]===r.Method)||[""])[0]; return `<tr><td><i class="im-dot" style="background:${DROP_COL[k]||"#999"}"></i>${esc(r.Method)}</td><td>${fmt(r["θ left"],1)}°</td><td>${fmt(r["θ right"],1)}°</td><td><b>${fmt(r["θ mean"],1)}°</b>${r["θ SD"]===r["θ SD"]&&r["θ SD"]!=null?` <span class="xs muted">± ${fmt(r["θ SD"],2)} fit</span>`:""}</td></tr>`; }).join("")}</table>`:""}
  ${(F.warn||[]).map(w=>`<div class="note warn xs">${esc(w)}</div>`).join("")}</div>`; }

/* ---------- 덧그림 ---------- */
function drawDrop(im,ox,v){ const F=S.dropFit&&S.dropFit[im.id]; const D=dropUi(); const B=(S.dropBase||{})[im.id];
  ox.save(); const path=(pts,close)=>{ ox.beginPath(); pts.forEach((p,i)=>{ const [sx,sy]=toScr(v,p[0],p[1]); i?ox.lineTo(sx,sy):ox.moveTo(sx,sy); }); if(close) ox.closePath(); };
  if(S.ui.panel==="drop"&&D.mode!=="pendant"&&B){ ox.setLineDash([6,4]); ox.strokeStyle="#ffe600"; ox.lineWidth=1.5; path(B); ox.stroke(); ox.setLineDash([]); }
  if(!F||!F.ov||F.z!==im.z||F.t!==im.t||F.ver!==im.ver){ ox.restore(); return; } /* 화소가 바뀌면(필터) 옛 덧그림을 그리지 않음 */ const ov=F.ov;
  if(F.contour&&F.contour.length){ ox.strokeStyle="rgba(0,220,255,.85)"; ox.lineWidth=1; ox.beginPath(); let pen=false; F.contour.forEach((p,i)=>{ const [sx,sy]=toScr(v,p[0],p[1]); if(F.edge&&F.edge[i]){ pen=false; return; } pen?ox.lineTo(sx,sy):ox.moveTo(sx,sy); pen=true; }); ox.stroke(); }
  if(ov.baseline){ ox.strokeStyle="#ffe600"; ox.lineWidth=1.5; path(ov.baseline); ox.stroke(); }
  const show=k=>F.mode==="pendant"||D.methods.includes(k);
  const pr={circle:0,ellipse:1,poly:2,yl:3}; (ov.curves||[]).slice().sort((a,b)=>(pr[a.m]||0)-(pr[b.m]||0)).forEach(c=>{ /* 주 방법(YL)을 맨 위에 */ const pts=Array.isArray(c)?c:c.pts, m=Array.isArray(c)?"yl":c.m; if(!show(m)) return; ox.strokeStyle=DROP_COL[m]||"#ff3ccf"; ox.lineWidth=2; path(pts); ox.stroke(); });
  (ov.needle||[]).forEach(l=>{ ox.strokeStyle="#ffe600"; ox.lineWidth=2; path(l); ox.stroke(); });
  if(ov.apex){ const [sx,sy]=toScr(v,ov.apex[0],ov.apex[1]); ox.fillStyle="#ff3ccf"; ox.beginPath(); ox.arc(sx,sy,3.5,0,7); ox.fill(); }
  ox.font="600 12px DM Sans, system-ui, sans-serif"; ox.textBaseline="middle";
  const tg=(ov.tangents||[]).filter(t=>show(t.m)); const main=["yl","poly","ellipse","circle"].find(k=>tg.some(t=>t.m===k));
  tg.forEach(t=>{ const [x0,y0]=toScr(v,t.p[0],t.p[1]); const L=t.len*v.zoom; ox.strokeStyle=DROP_COL[t.m]; ox.lineWidth=t.m===main?2:1.2; ox.beginPath(); ox.moveTo(x0,y0); ox.lineTo(x0+t.d[0]*L,y0+t.d[1]*L); ox.stroke();
    if(t.m===main){ const tx=x0+t.d[0]*L*1.08+(t.side==="L"?-6:6), ty=y0+t.d[1]*L*1.08; const s=t.theta.toFixed(1)+"°"; ox.textAlign=t.side==="L"?"right":"left"; ox.lineWidth=3; ox.strokeStyle="rgba(0,0,0,.7)"; ox.strokeText(s,tx,ty); ox.fillStyle=DROP_COL[t.m]; ox.fillText(s,tx,ty); } });
  ox.restore(); }

/* ---------- 아래 표 ---------- */
function bDrop(B){ const R=S.drops||[]; const keys=["Image","Frame","Mode","Method","θ left","θ right","θ mean","θ SD","γ","γ SD","γ scale SD","Bo","Wo","V (mm³)","A (mm²)","R0 (mm)","Base Ø (mm)","Height (mm)","Needle (mm)","Tilt (°)","rms (px)","γ est.","γ est. SD","Δρ","g","Baseline","Scale","Note"];
  const cols=keys.filter(k=>R.some(r=>r[k]!=null&&r[k]!=="")).map(k=>({k,label:({"γ":"γ (mN/m)","γ SD":"γ fit SD (mN/m)","γ scale SD":"γ needle-scale SD (mN/m)","γ est.":"γ est. (mN/m)","γ est. SD":"γ est. fit SD (mN/m)","θ left":"θ left (°)","θ right":"θ right (°)","θ mean":"θ mean (°)","θ SD":"θ fit SD (°)","Δρ":"Δρ (kg/m³)","g":"g (m/s²)"})[k]||k}));
  B.innerHTML=`<div class="im-row"><button class="btn xs" id="drCsv" ${R.length?"":"disabled"}>CSV</button><button class="btn xs" id="drCopy" ${R.length?"":"disabled"} title="Copy as tab-separated text (paste into Excel or Origin)">Copy</button><button class="btn xs" id="drClr" ${R.length?"":"disabled"}>Clear</button><span class="xs muted">Drop 패널의 측정이 한 줄씩 쌓입니다 (접촉각은 방법마다 한 줄)</span></div>${tableHTML(R,cols,{empty:"Drop 패널에서 Measure 를 누르세요."})}`;
  $q(B,"#drCsv").onclick=()=>downloadText("Drop.csv",csvOf(R,cols)); $q(B,"#drCopy").onclick=()=>copyText(tsvOf(R,cols));
  $q(B,"#drClr").onclick=()=>{ if(confirm2("drclr","Drop 표를 비울까요? 한 번 더 누르면 비웁니다")){ S.drops=[]; save(); drawBottom(); } }; }
