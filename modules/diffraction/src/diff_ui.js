/* ============================================================
   INNO Analytics — Diffraction 탭 (보기 키 "lattice" 그대로). 모드: FFT (HRTEM) · SAED · Phase DB
   FFT = 기존 Lattice FFT 화면(renderLattice) 을 그대로. SAED = 스폿·링 자동 읽기 + 보정 + 상 찾기 + 수동 비교.
   계산은 DIFF_WORK (Web Worker). 상 목록은 DIFF_DB (있으면) → 없으면 DIFF_PACK.
   화면 글자: 영어 뼈대 · 한국어 설명 (docs/UI-English-style.md). 휴대폰: 카드 쌓기, 표는 자체 스크롤, INNO_TOUCH 로 확대·이동·탭.
   ============================================================ */
(function(){
"use strict";
const DC=()=>window.DIFF_CORE, DI=()=>window.DIFF_INDEX, DW=()=>window.DIFF_WORK;
const $q=(s,r)=>(r||document).querySelector(s), $all=(s,r)=>[...(r||document).querySelectorAll(s)];
const E=s=>String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const say=(m,k)=>{ if(typeof window.toast==="function") window.toast(m,k); else console.log(m); };
const f=(x,n)=>x==null||!isFinite(x)?"–":(+x).toFixed(n);
const hklS=h=>h?h.map(v=>v<0?"<span class=\"ov\">"+(-v)+"</span>":String(v)).join(""):"";
const hklT=h=>h?"("+h.join(" ")+")":"";
const uvwS=z=>"["+z.map(v=>v<0?"<span class=\"ov\">"+(-v)+"</span>":String(v)).join("")+"]";
const LS_KEY="inno.diff.v1";
const P=(()=>{ let p={}; try{ p=JSON.parse(localStorage.getItem(LS_KEY)||"{}")||{}; }catch(e){ p={}; } return Object.assign({pattern:"auto",calSrc:"auto",kappaPct:100,tolPct:3,elements:"",elMode:"subset",layers:{pack:true,lab:true,local:true},maxIdx:5,radial:"auto",topN:0,showMask:true,showPred:true,distPct:1,rankCal:"auto",manual:{h1:"1 1 1",h2:"2 0 0"}},p); })();
const save=()=>{ try{ localStorage.setItem(LS_KEY,JSON.stringify(P)); }catch(e){} };
/* 화면 상태 (파일마다 새로) */
let S=null;
const fresh=()=>({key:null,name:"",img:null,w:0,h:0,meta:{},res:null,type:null,busy:false,err:null,center:null,lat:0,spotsEdit:null,cands:null,sel:-1,rcands:null,rsel:-1,
  view:null,tool:"view",disp:null,maskCv:null,gamma:"log",lohi:[0.5,99.95],barVal:null,barLen:null,manualScale:null,savedCal:null,mres:null,progress:""});

/* ---------- 원본 영상 ---------- */
function tagsOf(im){ return (im&&im.md&&im.md.__dm)||{}; }
function readSource(){ const rec=typeof curRec==="function"?curRec():null; if(!rec||!rec.images||!rec.images.length) return null; const c=typeof loadImageData==="function"?loadImageData():null; if(!c) return null;
  const im=c.im||{}; const key=rec.name+"|"+(im.uid||"")+"|"+c.w+"x"+c.h; return {key,name:rec.name+(rec.images.length>1?" · "+(im.title||""):""),c,im}; }
function metaOf(im,c){ const T=tagsOf(im); const g=k=>T[k]; const num=v=>{ const x=Array.isArray(v)?+v[0]:+v; return isFinite(x)?x:null; };
  const unit=String(im.dmUnit||"").trim(); const recip=/^(1\/nm|nm\^-1|nm-1|1\/Å|1\/A|1\/um|1\/µm)$/i.test(unit.replace(/\s/g,""));
  let scale=null; if(recip&&im.dmScale>0){ const u=unit.replace(/\s/g,"").toLowerCase(); scale=im.dmScale*(u==="1/å"||u==="1/a"?10:(u==="1/um"||u==="1/µm"?1e-3:1)); }
  const opMode=g("Microscope Info/Operation Mode")||g("Microscope Info/Imaging Mode")||""; const kV=num(g("Microscope Info/Voltage")); const ind=num(g("Microscope Info/Indicated Magnification"));
  const camPx=num(g("Acquisition/Device/CCD/Pixel Size (um)"))||num(g("Acquisition/Frame/CCD/Pixel Size (um)")); const bin=num(g("Acquisition/Parameters/High Level/Binning"))||num(g("DataBar/Binning"))||1;
  const scope=g("Microscope Info/Microscope")||g("Microscope Info/Name")||"";
  const diffMode=/diff/i.test(String(opMode))||recip;
  return {unit,recip,tagScale:scale,dmScale:im.dmScale,opMode:String(opMode),kV:kV?kV/1000:null,camLen:diffMode&&ind?ind:null,camPx,bin,scope:String(scope),diffMode,eightBit:c&&c.max<=255&&c.min>=0,rawMax:c&&c.max}; }
/* 전자 파장 (nm), 상대론 보정 */
const lambda=kV=>{ const V=kV*1e3; return 1.226426/Math.sqrt(V*(1+0.97848e-6*V)); };

/* ---------- 보정 ---------- */
function calib(){ const m=S.meta; let src=P.calSrc; let scale=null, ellip=1, phi=0, label="", ok=true;
  if(src==="auto") src=m.tagScale?"tag":(S.res&&S.res.scaleBar?"bar":"none");
  if(src==="tag"){ scale=m.tagScale; label="File tag"; if(!scale){ ok=false; } }
  else if(src==="bar"){ const L=S.barLen!=null?S.barLen:(S.res&&S.res.scaleBar&&S.res.scaleBar.lenPx); const v=S.barVal; scale=L>0&&v>0?v/L:null; label="Scale bar"; if(!scale) ok=false; }
  else if(src==="saved"){ const c=S.savedCal; if(c){ scale=c.scale; ellip=c.ellip||1; phi=c.phi||0; label="Saved: "+c.name; } else ok=false; }
  else if(src==="manual"){ scale=S.manualScale; label="Manual"; if(!(scale>0)) ok=false; }
  else { label="None"; }
  if(scale) scale*=P.kappaPct/100;
  // 순위에 보정을 쓸지: 표준 시료로 만든 보정(Saved)·직접 넣은 값만 기본으로 믿는다. 파일 태그·스케일바는 카메라 길이 공칭값이라 수 % 틀리기 쉬움 → 보정 없이 순위, κ 만 표시
  const trusted=P.rankCal==="on"?true:P.rankCal==="off"?false:(src==="saved"||src==="manual");
  return {src,scale:scale||null,ellip,phi,label,ok,err:P.tolPct/100,rank:!!(scale&&trusted)}; }
/* 보정 찌그러짐이 있으면 격자 벡터를 고쳐 '찌그러짐 없는 px' 로 (평균 축척 유지) */
function correctedFit(L,cal){ if(!L||!cal||!(cal.ellip>1.0001)) return L; const A=DI().affineFromParams(1,cal.ellip,cal.phi); const T=A.T; const ap=v=>[T[0][0]*v[0]+T[0][1]*v[1],T[1][0]*v[0]+T[1][1]*v[1]];
  const a=ap(L.R1), b=ap(L.R2), c=[a[0]+b[0],a[1]+b[1]]; const len=v=>Math.hypot(v[0],v[1]); const ang=(u,v)=>Math.acos(Math.max(-1,Math.min(1,(u[0]*v[0]+u[1]*v[1])/(len(u)*len(v)))))*180/Math.PI;
  return Object.assign({},L,{R1:a,R2:b,R3:c,len:[len(a),len(b),len(c)],ang12:ang(a,b),ang13:ang(a,c),ang23:ang(b,c)}); }

/* 링: 방위 적분 반지름은 타원의 단축 기준 → 단축 방향 축척 = 평균 축척 × √(장단축비) */
const ringScale=cal=>cal.scale&&S.res&&S.res.rings?cal.scale*Math.sqrt(S.res.rings.ratio||1):cal.scale;
/* 스폿 벡터(px) → |g| (1/nm): 저장한 보정의 찌그러짐까지 */
function gOf(v,cal){ if(!cal.scale) return null; if(cal.ellip>1.0001){ const T=DI().affineFromParams(1,cal.ellip,cal.phi).T; return Math.hypot(T[0][0]*v[0]+T[0][1]*v[1],T[1][0]*v[0]+T[1][1]*v[1])*cal.scale; } return Math.hypot(v[0],v[1])*cal.scale; }
/* ---------- 상 목록 ---------- */
async function records(){ const els=P.elements.trim()?P.elements.trim().split(/[\s,;]+/).filter(Boolean):null; let list=[];
  if(window.DIFF_DB&&DIFF_DB.list){ try{ if(DIFF_DB.init&&!DIFF_DB._inited){ await DIFF_DB.init(); DIFF_DB._inited=true; } list=await Promise.resolve(DIFF_DB.list({elements:els,elementMode:P.elMode,layers:Object.keys(P.layers).filter(k=>P.layers[k])})); }catch(e){ say("Phase DB 를 읽지 못했습니다: "+(e.message||e),"warn"); list=[]; } }
  else if(window.DIFF_PACK){ list=DIFF_PACK.phases.map(r=>Object.assign({layer:"pack"},r)); if(els){ const set=new Set(els); list=list.filter(r=>{ const k=Object.keys(DC().parseFormula(r.formula)); return P.elMode==="contains"?els.every(e=>k.includes(e)):k.every(e=>set.has(e)); }); } }
  return (list||[]).filter(r=>!r.hidden); }
const refOf=id=>{ try{ return window.DIFF_DB&&DIFF_DB.refOf?DIFF_DB.refOf(id):{id}; }catch(e){ return {id}; } };

/* ---------- 그리기 ---------- */
/* 표시용 캔버스: 긴 변 2048 넘으면 fd×fd 블록의 최댓값으로 줄여(작은 스폿이 사라지지 않게) — 휴대폰 메모리 (4096² RGBA 64 MB → 16 MB) */
const dispF=()=>Math.max(1,Math.ceil(Math.max(S.w,S.h)/2048));
function buildDisplay(){ const {img,w,h}=S; const n=w*h; const st=Math.max(1,Math.floor(n/300000)); const v=new Float32Array(Math.ceil(n/st)); let m=0; for(let i=0;i<n;i+=st){ const x=img[i]; if(x===x) v[m++]=x; } const u=v.subarray(0,m).sort();
  const lo=u[Math.floor(S.lohi[0]/100*(m-1))], hi=u[Math.floor(S.lohi[1]/100*(m-1))]; const fd=dispF(); const W=Math.ceil(w/fd), H=Math.ceil(h/fd);
  const cv=document.createElement("canvas"); cv.width=W; cv.height=H; const ctx=cv.getContext("2d"); const id=ctx.createImageData(W,H); const d=id.data;
  const L=S.gamma==="log"; const span=Math.max(1e-12,hi-lo); const lg=Math.log1p(span);
  for(let Y=0;Y<H;Y++) for(let X=0;X<W;X++){ let mx=-Infinity; for(let j=0;j<fd;j++){ const y=Y*fd+j; if(y>=h) break; const row=y*w; for(let i=0;i<fd;i++){ const x=X*fd+i; if(x>=w) break; const q=img[row+x]; if(q>mx) mx=q; } }
    let t=(mx-lo)/span; if(!(t>0)) t=0; else if(t>1) t=1; if(L) t=Math.log1p(t*span)/lg; const g=(t*255)|0; const k=(Y*W+X)*4; d[k]=d[k+1]=d[k+2]=g; d[k+3]=255; }
  ctx.putImageData(id,0,0); S.disp=cv; S.maskCv=null; }
function buildMask(){ const m=S.res&&S.res.mask; if(!m) return; const {w,h}=S; const fd=dispF(); const W=Math.ceil(w/fd), H=Math.ceil(h/fd); const cv=document.createElement("canvas"); cv.width=W; cv.height=H; const ctx=cv.getContext("2d"); const id=ctx.createImageData(W,H); const d=id.data;
  for(let Y=0;Y<H;Y++) for(let X=0;X<W;X++){ const x=Math.min(w-1,X*fd), y=Math.min(h-1,Y*fd); if(m[y*w+x]){ const k=(Y*W+X)*4; d[k]=60; d[k+1]=130; d[k+2]=255; d[k+3]=70; } } ctx.putImageData(id,0,0); S.maskCv=cv; }
function fitView(cv){ const W=cv.clientWidth||600, H=cv.clientHeight||600; const z=Math.min(W/S.w,H/S.h); S.view={zoom:z,x:-(W/z-S.w)/2,y:-(H/z-S.h)/2}; }
function draw(){ const cv=$q("#dfCv"); if(!cv||!S||!S.disp) return; const dpr=window.devicePixelRatio||1; const W=cv.clientWidth, H=cv.clientHeight; if(cv.width!==Math.round(W*dpr)||cv.height!==Math.round(H*dpr)){ cv.width=Math.round(W*dpr); cv.height=Math.round(H*dpr); }
  if(!S.view) fitView(cv); const v=S.view; const ctx=cv.getContext("2d"); ctx.setTransform(dpr,0,0,dpr,0,0); ctx.fillStyle="#0b0f0e"; ctx.fillRect(0,0,W,H);
  ctx.imageSmoothingEnabled=v.zoom<2; ctx.setTransform(dpr*v.zoom,0,0,dpr*v.zoom,-dpr*v.zoom*v.x,-dpr*v.zoom*v.y); ctx.drawImage(S.disp,0,0,S.w,S.h);
  if(P.showMask&&S.res&&S.res.mask){ if(!S.maskCv) buildMask(); if(S.maskCv) ctx.drawImage(S.maskCv,0,0,S.w,S.h); }
  ctx.setTransform(dpr,0,0,dpr,0,0); const X=x=>(x-v.x)*v.zoom, Y=y=>(y-v.y)*v.zoom; const big=window.INNO_TOUCH&&INNO_TOUCH.touchUI&&INNO_TOUCH.touchUI();
  const r=S.res; if(!r) return; ctx.lineWidth=1.5; ctx.font="12px system-ui, sans-serif";
  if(r.scaleBar){ const b=r.scaleBar; ctx.strokeStyle="#ff4d6d"; ctx.strokeRect(X(b.x0)-2,Y(b.y0)-2,(b.x1-b.x0+1)*v.zoom+4,(b.y1-b.y0+1)*v.zoom+4); }
  if(r.scaleLabel&&r.scaleLabel.box){ const b=r.scaleLabel.box; ctx.strokeStyle=r.scaleLabel.conf>=0.9?"#3ddc84":"#ffb020"; ctx.setLineDash([4,3]); ctx.strokeRect(X(b.x0)-3,Y(b.y0)-3,(b.x1-b.x0+1)*v.zoom+6,(b.y1-b.y0+1)*v.zoom+6); ctx.setLineDash([]); }
  if(S.type==="spots"&&r.spots){ const L=curLat(); const list=S.spotsEdit||r.spots.list; const idx=new Map(); if(L) L.spots.forEach(s=>idx.set(s.x.toFixed(2)+","+s.y.toFixed(2),s));
    for(const s of list){ const q=idx.get(s.x.toFixed(2)+","+s.y.toFixed(2)); const rad=Math.max(big?9:6,(s.radius||0)*v.zoom+3); ctx.strokeStyle=q&&q.indexed?"#3ddc84":"#ffb020"; ctx.lineWidth=s.saturated?2.5:1.5; ctx.beginPath(); ctx.arc(X(s.x),Y(s.y),rad,0,7); ctx.stroke(); }
    const c=L?L.center:r.spots.center; if(c){ ctx.strokeStyle="#ff3b3b"; ctx.lineWidth=2; const s=big?12:8; ctx.beginPath(); ctx.moveTo(X(c.x)-s,Y(c.y)); ctx.lineTo(X(c.x)+s,Y(c.y)); ctx.moveTo(X(c.x),Y(c.y)-s); ctx.lineTo(X(c.x),Y(c.y)+s); ctx.stroke(); }
    if(L){ const arrow=(vec,col,lab)=>{ const x0=X(L.center.x), y0=Y(L.center.y), x1=X(L.center.x+vec[0]), y1=Y(L.center.y+vec[1]); ctx.strokeStyle=col; ctx.fillStyle=col; ctx.lineWidth=2.2; ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x1,y1); ctx.stroke();
        const a=Math.atan2(y1-y0,x1-x0); ctx.beginPath(); ctx.moveTo(x1,y1); ctx.lineTo(x1-10*Math.cos(a-0.4),y1-10*Math.sin(a-0.4)); ctx.lineTo(x1-10*Math.cos(a+0.4),y1-10*Math.sin(a+0.4)); ctx.fill(); ctx.font="bold 13px system-ui"; ctx.fillText(lab,x1+6,y1-6); };
      arrow(L.R1,"#00d1ff","R1"); arrow(L.R2,"#ff7ad9","R2"); arrow(L.R3,"#ffe066","R3"); }
    const C=S.cands&&S.cands[S.sel]; if(C&&P.showPred&&L){ ctx.font=(big?"13px":"11px")+" ui-monospace,monospace";
      for(const row of C.verify.rows){ const x=X(row.x), y=Y(row.y); ctx.fillStyle=row.weak?"#ffb020":"#e8fff0"; ctx.fillText(row.h.join(" "),x+7,y+14); }
      for(const m of C.verify.missing){ const q=predictPx(L,C,m.h); if(!q) continue; const x=X(q.x), y=Y(q.y); ctx.strokeStyle="#ff4d6d"; ctx.setLineDash([3,3]); ctx.strokeRect(x-6,y-6,12,12); ctx.setLineDash([]); } } }
  if(S.type==="rings"&&r.rings){ const e=r.rings.ellipse; const C=S.rcands&&S.rcands[S.rsel]; const peaks=r.rings.peaks.filter(p=>!p.broad&&!p.nearStop);
    const ell=(rm,col,dash)=>{ ctx.strokeStyle=col; ctx.setLineDash(dash||[]); ctx.lineWidth=1.4; ctx.beginPath(); ctx.ellipse(X(e.cx),Y(e.cy),rm*e.a/e.b*v.zoom,rm*v.zoom,e.phi,0,7); ctx.stroke(); ctx.setLineDash([]); };
    for(const p of peaks) ell(p.r,"rgba(61,220,132,.85)");
    if(C&&P.showPred){ ctx.font="12px ui-monospace,monospace"; for(const q of C.rings){ const rm=q.g/C.scale; ell(rm,"rgba(255,122,217,.7)",[5,4]); const ang=-Math.PI/4+e.phi; ctx.fillStyle="#ff7ad9"; ctx.fillText(q.hkl.map(Math.abs).sort((a,b)=>b-a).join(""),X(e.cx+rm*Math.cos(ang)),Y(e.cy+rm*Math.sin(ang))); } }
    ctx.strokeStyle="#ff3b3b"; ctx.lineWidth=2; const s=8; ctx.beginPath(); ctx.moveTo(X(e.cx)-s,Y(e.cy)); ctx.lineTo(X(e.cx)+s,Y(e.cy)); ctx.moveTo(X(e.cx),Y(e.cy)-s); ctx.lineTo(X(e.cx),Y(e.cy)+s); ctx.stroke(); } }
function predictPx(L,C,h){ const z=DI().predictPx({center:L.center,R1:L.R1,R2:L.R2},{h1:C.h1,h2:C.h2},h); return z; }
const curLat=()=>S.res&&S.res.spots&&S.res.spots.lattices[S.lat]||null;

/* ---------- 분석 ---------- */
async function analyze(){ if(!S.img||S.busy) return; const my=S; S.busy=true; S.err=null; S.cands=null; S.rcands=null; S.sel=-1; S.rsel=-1; S.progress="패턴을 읽는 중… (큰 영상은 몇 초)"; paint();
  try{ const copy=new Float32Array(S.img); const r=await DW().run("analyze",{img:copy,w:S.w,h:S.h,opt:{pattern:P.pattern,radial:P.radial,topN:P.topN,maxLattices:3,center:S.center}},[copy.buffer]);
    if(S!==my) return; // 분석 중 다른 파일로 바뀜 — 버림
    S.res=r; S.type=P.pattern==="auto"?r.type:P.pattern; S.spotsEdit=null; S.lat=0; S.maskCv=null; if(S.barLen==null&&r.scaleBar) S.barLen=r.scaleBar.lenPx;
    const lb=r.scaleLabel; if(S.barVal==null&&lb&&lb.recip&&lb.valueInvNm>0){ S.barVal=+lb.valueInvNm.toPrecision(6); S.barFromOCR=true; } // 그림 속 이름표에서 읽은 값 (확인 가능)
    if(DW().mode==="main") say("Worker 를 쓸 수 없어 주 화면에서 계산했습니다","warn"); }
  catch(e){ if(S!==my) return; S.err=String(e.message||e); }
  if(S.res&&!S.err) S.idBusy=true; S.busy=false; S.progress=""; paint(); if(S.res&&!S.err) identify(); }
async function refit(){ if(!S.res||!S.res.spots) return; const my=S; const c=S.center||(S.res.spots.center); if(!c) return; S.busy=true; paint();
  try{ const list=(S.spotsEdit||S.res.spots.list); const sub=P.topN>0?list.slice().sort((a,b)=>b.I-a.I).slice(0,P.topN):list; const r=await DW().run("refit",{spots:sub,center:c,size:Math.max(S.w,S.h),w:S.w,h:S.h,hint:!!S.centerHint,radial:P.radial,max:3}); S.centerHint=false; if(S!==my) return;
    S.res.spots.lattices=r.lattices; S.res.spots.nUnindexed=r.nUnindexed; S.res.spots.center=Object.assign({},r.center||c,{source:S.center?"user":(c.source||"")}); if(S.center&&r.center) S.center={x:r.center.x,y:r.center.y}; S.lat=0; S.cands=null; S.sel=-1; }catch(e){ S.err=String(e.message||e); }
  S.idBusy=true; S.busy=false; paint(); identify(); }
async function identify(){ if(!S.res) return; const my=S; S.idBusy=true; const recs=await records(); if(S!==my) return; if(!recs.length){ S.cands=[]; S.rcands=[]; S.idBusy=false; paint(); return; } const cal=calib(); paint();
  try{ if(S.type==="spots"){ const L=curLat(); if(!L){ S.cands=[]; } else { const Lc=correctedFit(L,cal); const mk=S.res.mask?new Uint8Array(S.res.mask):null; const r=await DW().run("rankSpots",{fit:Lc,records:recs,opt:{scale:cal.rank?cal.scale:null,scaleErr:cal.err,distortion:P.distPct/100,maxIdx:P.maxIdx,gmax:Math.max(15,cal.scale?Lc.len[2]*cal.scale*1.5:15),w:S.w,h:S.h,mask:mk,posFit:{center:L.center,R1:L.R1,R2:L.R2},nVerify:80,dynamic:P.intModel!=="kin",kV:S.meta.kV||200}},mk?[mk.buffer]:[]);
        if(S!==my) return; S.cands=r.cands.slice(0,40); S.badRecs=r.bad; S.sel=S.cands.length?0:-1; } }
    else { const pk=S.res.rings.peaks; const gmax=cal.scale?Math.max(...pk.map(p=>p.r))*ringScale(cal)*1.05:20; const r=await DW().run("matchRings",{peaks:pk,records:recs,opt:{scale:cal.rank?ringScale(cal):null,scaleErr:cal.err},gShow:gmax});
      if(S!==my) return; S.rcands=r.cands.slice(0,30); S.badRecs=r.bad; S.rsel=S.rcands.length?0:-1; } }
  catch(e){ S.err=String(e.message||e); }
  S.idBusy=false; paint(); }

/* ---------- 화면 ---------- */
/* 모드: 사용자가 고른 것(파일마다) → 없으면 파일 종류로: 역공간 교정·회절 모드 태그 → SAED, 실공간 영상 → FFT, 파일 없음 → SAED(열기 안내) */
const MODE={key:null,mode:null};
function curMode(){ const src=readSource(); const key=src?src.key:"(none)"; if(MODE.key===key&&MODE.mode) return MODE.mode; if(!src) return "saed";
  const m=metaOf(src.im,src.c); if(m.recip||m.diffMode) return "saed"; if(m.unit&&/^(nm|µm|um|Å|A|pm)$/i.test(m.unit)&&!/diff|saed|sad\b/i.test(src.name)) return "fft";
  // 교정 정보 없음: 회절도는 대부분 어둡고(중앙값 ≪ 최댓값) 몇 점만 밝다, HRTEM 은 중간 회색 → 중앙값/범위 < 0.15 이면 SAED
  if(/diff|saed|\bsad/i.test(src.name)) return "saed"; const c=src.c; if(c&&c.cum&&c.total){ const half=c.total/2; let k=0; while(k<c.NB-1&&c.cum[k]<half) k++; return k/(c.NB-1)<0.15?"saed":"fft"; } return "saed"; }
function render(root){ if(!S) S=fresh(); const mode=curMode(); P.mode=mode;
  root.innerHTML=`<div class="card df-head"><h3>Diffraction <span class="hint">회절도 읽기 · 상 찾기 · 보정</span></h3>
    <div class="seg" id="dfMode" role="tablist">${[["fft","FFT (HRTEM)"],["saed","SAED"],["db","Phase DB"]].map(([k,t])=>`<button data-m="${k}" class="${mode===k?"active":""}">${t}</button>`).join("")}</div></div><div id="dfBody"></div>`;
  $all("#dfMode button",root).forEach(b=>b.onclick=()=>{ const src=readSource(); MODE.key=src?src.key:"(none)"; MODE.mode=b.dataset.m; render(root); });
  const body=$q("#dfBody",root);
  if(P.mode==="fft"){ const rec=typeof curRec==="function"?curRec():null; if(!rec||!rec.images||!rec.images.length){ body.innerHTML=`<div class="card"><div class="note">FFT 는 HRTEM 영상이 필요합니다. 왼쪽에서 파일을 여세요.</div></div>`; return; } if(typeof renderLattice==="function") renderLattice(body); return; }
  if(P.mode==="db"){ if(window.DIFF_DB_UI&&DIFF_DB_UI.render) DIFF_DB_UI.render(body,{}); else body.innerHTML=`<div class="card"><div class="note">Phase DB 화면을 불러오지 못했습니다.</div></div>`; return; }
  renderSaed(body); }
function renderSaed(body){ const src=readSource();
  if(!src){ body.innerHTML=`<div class="card"><h3>SAED</h3><div class="note">회절 영상(.dm3 · .dm4 · .tif)을 여세요. 왼쪽 파일 목록에 추가되고 여기서 바로 분석합니다.</div>
      <div class="row" style="margin-top:10px"><button class="btn primary" id="dfOpen">Open file…</button></div></div>`; $q("#dfOpen",body).onclick=openFile; return; }
  if(S.key!==src.key){ S=fresh(); S.key=src.key; S.name=src.name; const c=src.c; S.w=c.w; S.h=c.h; S.img=Float32Array.from(c.data); S.meta=metaOf(src.im,c); buildDisplay(); S.view=null; setTimeout(analyze,0); }
  body.innerHTML=`<div class="df-grid">
    <div class="df-left">
      <div class="card df-viewer"><div class="row df-vbar">
        <div class="seg" id="dfTool">${[["view","Pan"],["center","Set center"],["spots","Edit spots"]].map(([k,t])=>`<button data-t="${k}" class="${S.tool===k?"active":""}">${t}</button>`).join("")}</div>
        <button class="btn sm" id="dfFit" title="Fit to window">Fit</button>
        <label class="df-ck"><input type="checkbox" id="dfMaskCk" ${P.showMask?"checked":""}> Mask</label>
        <label class="df-ck"><input type="checkbox" id="dfPredCk" ${P.showPred?"checked":""}> Indices</label>
        <select id="dfGamma" title="Display">${[["log","Log"],["lin","Linear"]].map(([k,t])=>`<option value="${k}" ${S.gamma===k?"selected":""}>${t}</option>`).join("")}</select>
      </div><div class="df-stage"><canvas id="dfCv"></canvas>${S.busy?`<div class="df-busy">${E(S.progress||"계산 중…")}</div>`:""}</div>
      <div class="muted df-hint" id="dfHint">${toolHint()}</div></div>
      <div id="dfProfile"></div>
    </div>
    <div class="df-right">
      ${sourceCard(src)}
      ${calCard()}
      ${patternCard()}
      <div id="dfResult">${S.err?`<div class="card"><div class="note danger">${E(S.err)}</div></div>`:""}${resultCards()}</div>
    </div></div>`;
  wire(body); draw(); if(S.type==="rings") drawProfile(); }
function toolHint(){ if(S.tool==="center") return "영상을 눌러 직접빔(중심)을 지정합니다. 지정한 중심으로 격자를 다시 맞춥니다."; if(S.tool==="spots") return "스폿을 누르면 빼고, 빈 곳을 누르면 그 자리의 밝은 점을 더합니다."; return "끌어서 이동 · 휠/두 손가락으로 확대 · 두 번 눌러 맞춤"; }
function sourceCard(src){ const m=S.meta; const warn=[];
  if(m.unit&&!m.recip&&m.dmScale>0&&(m.diffMode||/diff/i.test(S.name))) warn.push(`파일의 교정 단위가 <b>${E(m.unit)}</b> (실공간) 입니다 — 회절 패턴이면 교정값이 잘못 저장된 것입니다. 스케일바·저장한 보정·표준 시료로 보정하세요.`);
  else if(m.unit&&!m.recip&&m.dmScale>0) warn.push(`교정 단위가 <b>${E(m.unit)}</b> 입니다 — 회절 영상이 아니라면 FFT (HRTEM) 모드를 쓰세요.`);
  if(m.eightBit) warn.push("8비트 표시용 사본입니다 — 포화된 스폿은 덩어리 중심으로 위치를 잡고 세기는 하한값입니다.");
  const L=calib(); const camL=m.kV&&m.camPx&&L.scale?(m.camPx*1e-6*m.bin)/(lambda(m.kV)*1e-9*L.scale*1e9)*1e3:null;
  return `<div class="card"><h3>Source</h3><div class="df-kv">
    <span>File</span><b class="df-trunc" title="${E(S.name)}">${E(S.name)}</b><span>Size</span><b>${S.w} × ${S.h}</b>
    ${m.opMode?`<span>Mode</span><b>${E(m.opMode)}</b>`:""}${m.kV?`<span>Voltage</span><b>${f(m.kV,0)} kV · λ ${f(lambda(m.kV)*1e3,3)} pm</b>`:""}
    ${m.camLen?`<span>Camera length</span><b>${f(m.camLen,0)} mm (indicated)</b>`:""}${camL?`<span>Effective L</span><b>${f(camL,0)} mm <span class="muted">(화소 ${f(m.camPx*m.bin,1)} µm · 지금 보정 기준)</span></b>`:""}
    ${m.scope?`<span>Microscope</span><b>${E(m.scope)}</b>`:""}<span>File calibration</span><b>${m.dmScale?`${E(m.dmScale.toPrecision(6))} ${E(m.unit)}/px`:"없음"}</b></div>
    ${warn.map(w=>`<div class="note warn" style="margin-top:8px">${w}</div>`).join("")}
    <div class="row" style="margin-top:10px"><button class="btn sm" id="dfOpen">Open file…</button></div></div>`; }
function calCard(){ const m=S.meta, c=calib(); const bar=S.res&&S.res.scaleBar; const opts=[["auto","Auto"],["tag","File tag"],["bar","Scale bar"],["saved","Saved"],["manual","Manual"],["none","None"]];
  return `<div class="card"><h3>Calibration <span class="hint">${c.scale?`${c.scale.toPrecision(5)} 1/nm/px · ${E(c.label)}`:"보정 없음 — 길이비·각도로만 비교"}</span></h3>
    <div class="seg df-wrap" id="dfCal">${opts.map(([k,t])=>`<button data-c="${k}" class="${P.calSrc===k?"active":""}">${t}</button>`).join("")}</div>
    <div class="df-calbody">
    ${c.src==="tag"?(m.tagScale?`<div class="muted">파일 태그 ${E(m.tagScale.toPrecision(6))} 1/nm/px. 표준 시료로 확인한 값과 다르면 κ 로 고치세요.</div>`:`<div class="note warn">파일에 역공간(1/nm) 교정값이 없습니다.</div>`):""}
    ${c.src==="bar"?`<div class="row"><label class="f">Bar length (px)<input type="number" id="dfBarLen" value="${S.barLen!=null?S.barLen:""}" step="0.5" min="1" style="width:100px"></label>
       <label class="f">Bar value (1/nm)<input type="number" id="dfBarVal" value="${S.barVal!=null?S.barVal:""}" step="any" min="0" placeholder="e.g. 5" style="width:100px"></label></div>
       <div class="muted">${bar?`그림 속 막대를 찾았습니다 (${bar.lenPx} px${bar.hollow?", 속 빈 테두리 — 선 중심 기준":""}).`:"그림에서 스케일바를 찾지 못했습니다 — 길이를 직접 넣으세요."} ${labelNote()}</div>`:""}
    ${c.src==="saved"?savedCalUI():""}
    ${c.src==="manual"?`<label class="f">Scale (1/nm per px)<input type="number" id="dfManScale" value="${S.manualScale||""}" step="any" min="0" style="width:140px"></label>`:""}
    <div class="row" style="margin-top:8px"><label class="f" title="REDI '측정 오차 %' 와 같은 배율 계수">κ factor (%)<input type="number" id="dfKappa" value="${P.kappaPct}" step="0.1" min="50" max="150" style="width:90px"></label>
      <label class="f">Tolerance (%)<input type="number" id="dfTol" value="${P.tolPct}" step="0.5" min="0.2" max="20" style="width:80px"></label>
      <label class="f" title="렌즈 찌그러짐 허용">Distortion (%)<input type="number" id="dfDist" value="${P.distPct}" step="0.5" min="0" max="15" style="width:80px"></label></div>
    <div class="row"><label class="f" title="후보 순위에 보정값을 쓸지 (Auto = 저장한 표준 보정·직접 넣은 값일 때만)">Rank with calibration<select id="dfRankCal">${[["auto","Auto"],["on","On"],["off","Off"]].map(([k,t])=>`<option value="${k}" ${P.rankCal===k?"selected":""}>${t}</option>`).join("")}</select></label><span class="muted">${c.rank?"순위에 보정 사용":"순위는 보정 없이"}</span></div>
    ${c.ellip>1.0001?`<div class="muted">저장한 보정의 찌그러짐 ${f((c.ellip-1)*100,2)} % (1/nm/px 가 큰 방향 ${f(c.phi*180/Math.PI,1)}°) 을 격자에 적용합니다.</div>`:""}</div></div>`; }
function labelNote(){ const lb=S.res&&S.res.scaleLabel; if(!lb) return S.res&&S.res.scaleBar?"이름표 숫자를 읽지 못했습니다 — 값을 넣으세요.":"";
  const t=`이름표 "${E(lb.text)}" 를 읽었습니다 (확신 ${f(lb.conf*100,0)} %)`; if(!lb.recip) return `<span class="df-kbad">${t} — 실공간 단위입니다. 회절 영상이 아니면 FFT (HRTEM) 모드를 쓰세요.</span>`;
  if(S.barFromOCR&&S.barVal===+lb.valueInvNm.toPrecision(6)) return lb.conf>=0.9?`${t} → 값에 넣었습니다.`:`<span class="df-kbad">${t} → 값에 넣었지만 확신이 낮습니다. 그림과 맞는지 확인하세요.</span>`; return t+"."; }
function savedCalUI(){ const list=S.calList||[]; return list.length?`<label class="f">Saved calibration<select id="dfSaved"><option value="">— choose —</option>${list.map((c,i)=>`<option value="${i}" ${S.savedCal&&S.savedCal.id===c.id?"selected":""}>${E(c.name)} · ${E(c.scale.toPrecision(5))}${c.ellip>1.0001?` · e ${f((c.ellip-1)*100,1)}%`:""}</option>`).join("")}</select></label>`:`<div class="muted">저장한 보정이 없습니다. 표준 시료(Au·Al·Si 등) 패턴을 분석한 뒤 결과에서 'Use as calibration' 을 누르세요.</div>`; }
function patternCard(){ return `<div class="card"><h3>Pattern</h3><div class="row">
    <div class="seg" id="dfPat">${[["auto","Auto"],["spots","Spots"],["rings","Rings"]].map(([k,t])=>`<button data-p="${k}" class="${P.pattern===k?"active":""}">${t}</button>`).join("")}</div>
    <button class="btn primary" id="dfRun" ${S.busy?"disabled":""}>Analyze</button></div>
    ${S.res?`<div class="muted" style="margin-top:6px">${S.type==="rings"?"링(다결정) 패턴":"스폿(단결정) 패턴"}${P.pattern==="auto"?" — 자동 판정":""} · ${f(S.res.ms/1000,1)} s${S.res.spots?` · 스폿 ${S.res.spots.list.length}`:""}</div>`:""}
    ${S.type==="spots"?`<div class="row" style="margin-top:8px"><label class="f" title="센 스폿 N 개만으로 격자 (겹친 상이 있을 때)">Use spots<select id="dfTopN">${[[0,"All"],[8,"Strongest 8"],[12,"Strongest 12"],[16,"Strongest 16"],[24,"Strongest 24"]].map(([k,t])=>`<option value="${k}" ${+P.topN===k?"selected":""}>${t}</option>`).join("")}</select></label>
      <label class="f">Radial distortion<select id="dfRad">${[["auto","Auto"],["on","Fit"],["off","Off"]].map(([k,t])=>`<option value="${k}" ${(P.radial===true?"on":P.radial===false?"off":"auto")===k?"selected":""}>${t}</option>`).join("")}</select></label></div>`:""}</div>`; }
function latCard(){ const L=curLat(); const r=S.res.spots; if(!L) return `<div class="card"><h3>Lattice</h3><div class="note warn">격자를 맞추지 못했습니다 — 중심을 직접 지정하거나(Set center) 'Use spots' 를 바꿔 보세요.${r.center&&r.center.warning?"<br>"+E(r.center.warning):""}</div></div>`;
  const cal=calib(); const Lc=correctedFit(L,cal); const s=cal.scale; const rows=[0,1,2].map(i=>{ const g=s?Lc.len[i]*s:null, sg=s?Math.hypot(L.sig[i]/L.len[i],cal.err)*g:null; return `<tr><td class="l"><b>R${i+1}</b></td><td class="num">${f(L.len[i],2)} ± ${f(L.sig[i],2)}</td><td class="num">${g?f(g,4):"–"}</td><td class="num">${g?f(1/g,4)+" ± "+f(sg/g/g,4):"–"}</td></tr>`; }).join("");
  const lats=r.lattices; const c=r.center;
  return `<div class="card"><h3>Lattice <span class="hint">${L.nIndexed} / ${(S.spotsEdit||r.list).length} spots indexed · σ ${f(L.sigmaPx,2)} px</span></h3>
    ${lats.length>1?`<div class="seg" id="dfLat" style="margin-bottom:8px">${lats.map((q,i)=>`<button data-l="${i}" class="${S.lat===i?"active":""}">Lattice ${i+1}</button>`).join("")}</div>`:""}
    <div class="tblwrap"><table class="tbl"><thead><tr><th class="l">Vector</th><th>|R| (px)</th><th>g (1/nm)</th><th>d (nm)</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div class="df-kv" style="margin-top:8px"><span>∠R1·R2</span><b>${f(Lc.ang12,2)}° ± ${f(L.sigAng,2)}°</b><span>∠R1·R3</span><b>${f(Lc.ang13,2)}°</b><span>∠R2·R3</span><b>${f(Lc.ang23,2)}°</b>
      <span>R2/R1</span><b>${f(Lc.len[1]/Lc.len[0],4)}</b><span>Center</span><b>${f(L.center.x,1)}, ${f(L.center.y,1)} <span class="muted">(${E(c&&c.source==="user"?"직접 지정":c&&c.origin==="stop"?"빔스토퍼 뒤 · Friedel 짝":c&&c.origin==="spot"?"직접빔":"Friedel 짝")})</span></b>
      ${L.radialUsed?`<span>Radial distortion</span><b>${f(L.radialPct,3)} % at R</b>`:""}</div>
    ${c&&c.warning?`<div class="note warn" style="margin-top:8px">${E(c.warning)}</div>`:""}
    ${r.nUnindexed>3?`<div class="muted" style="margin-top:6px">지수 안 붙은 스폿 ${r.nUnindexed} 개 — 이중 회절·제2상·석출물일 수 있습니다${lats.length>1?" (다른 격자 탭 참고)":""}.</div>`:""}
    <div class="muted" style="margin-top:6px">길이비와 각도는 보정과 무관합니다.</div></div>`; }
function ringCard(){ const R=S.res.rings; const cal=calib(); const s=ringScale(cal); const pk=R.peaks;
  const rows=pk.map((p,i)=>{ const g=s?p.r*s:null; const flag=p.nearStop?"beam-stop edge":p.broad?"broad":""; return `<tr class="${flag?"df-dim":""}"><td>${i+1}</td><td class="num">${f(p.r,2)}</td><td class="num">${g?f(g,4):"–"}</td><td class="num">${g?f(1/g,4):"–"}</td><td class="num">${f(p.fwhm,1)}</td><td class="num">${f(p.snr,0)}</td><td class="l muted">${flag}</td></tr>`; }).join("");
  return `<div class="card"><h3>Rings <span class="hint">center ${f(R.center.x,1)}, ${f(R.center.y,1)} · a/b ${f(R.ratio,4)} · ${R.ringFits.length} rings fitted</span></h3>
    <div class="tblwrap" style="max-height:40vh"><table class="tbl"><thead><tr><th>#</th><th>r (px)</th><th>g (1/nm)</th><th>d (nm)</th><th>FWHM (px)</th><th>S/N</th><th class="l">Note</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div class="muted" style="margin-top:6px">반지름은 타원을 원으로 편 단축 기준입니다. 빔스토퍼 가장자리·넓은 봉우리는 상 찾기에서 뺍니다.</div></div>`; }
function idCard(){ const cal=calib(); const els=P.elements;
  const head=`<div class="card"><h3>Phase ID <span class="hint">${cal.scale?"보정값 + 길이비·각도·세기":"보정 없이 — 길이비·각도·세기"}</span></h3>
    <div class="row"><label class="f" style="flex:1;min-width:160px" title="EDS 로 확인한 원소 (예: Na Y F)">Elements<input type="text" id="dfEls" value="${E(els)}" placeholder="e.g. Na Y F  (empty = all)"></label>
      <label class="f">Match<select id="dfElMode"><option value="subset" ${P.elMode==="subset"?"selected":""}>Only these</option><option value="contains" ${P.elMode==="contains"?"selected":""}>Contains all</option></select></label>
      ${S.type==="spots"?`<label class="f" title="동역학 = Bloch-wave (이중 회절·Gjønnes–Moodie 소멸 반영, 상위 12 후보)">Intensities<select id="dfIntM"><option value="dyn" ${P.intModel!=="kin"?"selected":""}>Dynamical</option><option value="kin" ${P.intModel==="kin"?"selected":""}>Kinematic</option></select></label><label class="f">Max zone index<select id="dfMaxIdx">${[3,4,5,6,7].map(k=>`<option ${+P.maxIdx===k?"selected":""}>${k}</option>`).join("")}</select></label>`:""}</div>
    <div class="row df-layers">${[["pack","Built-in"],["lab","Lab"],["local","Local"]].map(([k,t])=>`<label class="df-ck"><input type="checkbox" data-layer="${k}" ${P.layers[k]?"checked":""}> ${t}</label>`).join("")}<button class="btn sm" id="dfId" ${S.idBusy?"disabled":""}>Search</button></div>`;
  if(S.idBusy) return head+`<div class="muted" style="margin-top:8px">후보를 계산하는 중…</div></div>`;
  if(S.type==="spots"){ const C=S.cands; if(!C) return head+"</div>"; if(!C.length) return head+`<div class="note" style="margin-top:8px">후보가 없습니다 — 원소·Phase DB 범위를 넓혀 보세요.</div></div>`;
    const tb=C.map((c,i)=>{ const dErr=cal.scale?c.d.map((d,j)=>{ const L=correctedFit(curLat(),cal); const dm=1/(L.len[j]*cal.scale); return (dm-d)/d*100; }):null;
      return `<tr data-i="${i}" class="${i===S.sel?"df-sel":""}"><td>${i+1}</td><td class="l"><b>${E(c.phaseName)}</b><div class="muted">${E(c.formula||"")} · ${E(c.sg&&c.sg.symbol||"")} ${c.layer?`<span class="df-tag">${E(c.layer)}</span>`:""}</div></td><td>${uvwS(c.zone)}</td>
        <td class="mono">${hklS(c.h1)}<br>${hklS(c.h2)}</td><td class="num">${f(c.d[0],4)}<br>${f(c.d[1],4)}</td>${dErr?`<td class="num">${f(dErr[0],1)}<br>${f(dErr[1],1)}</td>`:""}<td class="num">${f(c.ang12,1)}°</td>
        <td class="num">${f(c.fom,2)}</td><td class="num">${c.hasF?f(c.rho,2):"–"}</td><td class="num">${c.hasF?f(c.forb*100,0)+"%":"–"}</td><td class="num">${c.hasF?f(c.missFrac*100,0)+"%":"–"}</td>${kCell(cal.scale?c.scaleImplied/cal.scale:null,c.scaleImplied,cal)}<td class="num"><b>${f(c.score,2)}</b></td></tr>`; }).join("");
    return head+`<div class="tblwrap" style="margin-top:10px;max-height:52vh"><table class="tbl df-cand"><thead><tr><th>#</th><th class="l">Phase</th><th>Zone</th><th>R1 / R2 (hkl)</th><th>d calc (nm)</th>${cal.scale?"<th>Δd (%)</th>":""}<th>∠</th>
      <th title="길이비·각도(+보정) χ²/자유도">FOM</th><th title="관측 세기 ↔ 예측 세기(동역학 또는 운동학) 순위 상관 — 포화 스폿 제외">ρ</th><th title="운동학적으로 0 인 자리의 관측 세기 비율 (이중 회절 가능)">Forb.</th><th title="예측 강한 반사 중 안 보이는 비율">Miss</th><th title="${cal.scale?"후보가 맞다면 필요한 보정 배율 = 예측 g ÷ 측정 g":"후보가 맞다면 필요한 1/nm/px"}">${cal.scale?"κ":"Scale"}</th><th>Score</th></tr></thead><tbody>${tb}</tbody></table></div>
      <div class="muted" style="margin-top:6px">Score = FOM + 5·Forb + (1 − ρ) + 3·Miss (작을수록 좋음). ${rankNote(cal)}</div>${tieNote(C,"low")}
      ${candDetail()}</div>`; }
  const C=S.rcands; if(!C) return head+"</div>"; if(!C.length) return head+`<div class="note" style="margin-top:8px">후보가 없습니다.</div></div>`;
  const tb=C.map((c,i)=>`<tr data-i="${i}" class="${i===S.rsel?"df-sel":""}"><td>${i+1}</td><td class="l"><b>${E(c.phaseName)}</b><div class="muted">${E(c.formula||"")} · ${E(c.sg&&c.sg.symbol||"")} ${c.layer?`<span class="df-tag">${E(c.layer)}</span>`:""}</div></td>
    <td class="num">${c.pairs.length}</td><td class="num">${f(c.scale,6)}</td>${kCell(cal.scale?c.scale/ringScale(cal):null,null,cal)}<td class="num">${f(c.explained*100,0)}%</td><td class="num">${f(c.missing*100,0)}%</td><td class="num">${f(c.rho,2)}</td><td class="num">${f(c.dev*100,2)}%</td><td class="num"><b>${f(c.score,3)}</b></td></tr>`).join("");
  return head+`<div class="tblwrap" style="margin-top:10px;max-height:52vh"><table class="tbl df-cand"><thead><tr><th>#</th><th class="l">Phase</th><th>Rings</th><th>Scale (1/nm/px)</th><th>κ</th><th title="짝지은 관측 세기 비율">Explained</th><th title="관측 범위 안 강한 예측 링이 빠진 비율">Missing</th><th>ρ</th><th>Mean dev.</th><th>Score</th></tr></thead><tbody>${tb}</tbody></table></div>
    <div class="muted" style="margin-top:6px">Score = 설명한 세기 × (1 − 빠진 링) × exp(−어긋남/허용) × (1+ρ)/2 (클수록 좋음). ${rankNote(cal)}</div>${tieNote(C,"high")}${ringDetail()}</div>`; }
/* 1·2 위가 거의 같으면(점수 차 < 5 %) 알림 — 같은 구조형(예: fcc 금속끼리)은 보정이 1 % 안으로 맞거나 원소가 없으면 가를 수 없다 */
function tieNote(list,better){ if(!list||list.length<2) return ""; const a=list[0].score, b=list[1].score; const close=better==="low"?(b-a)<=0.05*Math.max(Math.abs(a),0.05):(a-b)<=0.05*Math.abs(a);
  return close?`<div class="note warn" style="margin-top:8px">1 위(${E(list[0].phaseName)})와 2 위(${E(list[1].phaseName)})의 점수가 거의 같습니다 — 이 패턴만으로는 가를 수 없습니다. Elements 에 EDS 원소를 넣거나 표준 시료로 만든 보정을 쓰세요.</div>`:""; }
function kCell(k,sc,cal){ if(k==null) return `<td class="num">${sc!=null?f(sc,5):"–"}</td>`; const ok=Math.abs(k-1)<=cal.err; return `<td class="num ${ok?"df-kok":"df-kbad"}" title="${ok?"보정 허용 안":"보정 허용 밖"}">${f(k,3)}</td>`; }
function rankNote(cal){ if(!cal.scale) return "보정 없이 — 모양(길이비·각도)과 세기가 같은 상은 가를 수 없습니다. 원소나 보정으로 좁히세요.";
  return cal.rank?`보정(${E(cal.label)}, 허용 ±${f(cal.err*100,1)} %)을 순위에 씁니다.`:`순위는 보정 없이 매기고 κ 만 보여 줍니다 — ${E(cal.label)} 는 공칭값이라 수 % 틀릴 수 있습니다. 초록 κ = 허용 ±${f(cal.err*100,1)} % 안. 보정을 믿으면 Calibration 의 'Rank' 를 켜세요.`; }
function candDetail(){ const c=S.cands&&S.cands[S.sel]; if(!c) return ""; const L=curLat(); const cal=calib(); const v=c.verify;
  const dy=!!c.dyn; const note=r=>dy?(r.gm?"Gjønnes–Moodie 소멸 — 정대축이면 0, 조금 기울면 생김":r.dz?"ZOLZ 로는 생길 수 없음 (HOLZ·다른 상?)":r.dd?"운동학 금지 — 이중 회절로 생김":""):(r.absent?"공간군 소멸 — 이중 회절?":r.weak?"운동학적으로 약함":"");
  const rows=v.rows.slice().sort((a,b)=>Math.hypot(...a.n)-Math.hypot(...b.n)).slice(0,30).map(r=>`<tr><td class="mono">${hklS(r.h)}</td><td class="num">${f(r.d,4)}</td><td class="num">${r.I==null?"–":f(r.I*100,1)}</td>${dy?`<td class="num">${r.Id==null?"–":f(r.Id*100,1)}</td>`:""}<td class="l muted">${note(r)}</td></tr>`).join("");
  const missL=dy?c.dyn.missing:v.missing;
  return `<div class="df-detail"><h4>${E(c.phaseName)} ${uvwS(c.zone)}</h4>
    <div class="df-kv"><span>R1</span><b class="mono">${hklS(c.h1)} · d ${f(c.d[0],4)} nm</b><span>R2</span><b class="mono">${hklS(c.h2)} · d ${f(c.d[1],4)} nm</b><span>R3</span><b class="mono">${hklS(c.h3)} · d ${f(c.d[2],4)} nm</b>
      <span>∠R1·R2</span><b>${f(c.ang12,2)}° (측정 ${f(correctedFit(L,cal).ang12,2)}°)</b><span>Spots on lattice</span><b>${v.nObs} · 강한 예측 ${v.nStrongPred} 중 빠짐 ${v.missing.length}</b></div>
    ${missL.length?`<div class="muted">빠진 ${dy?"(동역학) ":""}강한 반사: ${missL.slice(0,8).map(m=>`<span class="mono">${hklS(m.h)}</span> (${f(m.I*100,0)}%)`).join(", ")}</div>`:""}
    ${dy?`<div class="muted">세기 모형: 동역학 Bloch-wave (ZOLZ ${c.dyn.nBeams} 빔, ${c.dyn.kV} kV, 두께 ${c.dyn.tRange[0]}–${c.dyn.tRange[1]} nm 평균). 운동학이었다면 ρ ${c.kin&&c.kin.rho!=null?f(c.kin.rho,2):"–"} · Forb ${c.kin&&c.kin.forb!=null?f(c.kin.forb*100,0)+"%":"–"} · Score ${c.kin?f(c.kin.score,2):"–"}.</div>`:""}
    <div class="tblwrap" style="max-height:30vh;margin-top:6px"><table class="tbl"><thead><tr><th>hkl</th><th>d (nm)</th><th>|F|² rel. (%)</th>${dy?"<th>Dyn. I (%)</th>":""}<th class="l">Note</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div class="row" style="margin-top:8px"><button class="btn sm" id="dfCalib">Use as calibration</button><button class="btn sm" id="dfCsv">Export CSV</button><button class="btn sm" id="dfPng">Export PNG</button><button class="btn sm" id="dfJson">Export JSON</button></div></div>`; }
function ringDetail(){ const c=S.rcands&&S.rcands[S.rsel]; if(!c) return ""; const rows=c.pairs.map(p=>`<tr><td class="num">${f(p.r,2)}</td><td class="num">${f(p.d,4)}</td><td class="mono">${hklS(p.hkl.map(Math.abs).sort((a,b)=>b-a))}</td><td class="num">${p.m}</td><td class="num">${f(1/p.gPred,4)}</td><td class="num">${f(p.e*100,2)}</td></tr>`).join("");
  return `<div class="df-detail"><h4>${E(c.phaseName)} <span class="muted">scale ${f(c.scale,6)} ± ${c.scaleErr.toExponential(1)} 1/nm/px</span></h4>
    <div class="tblwrap" style="max-height:30vh"><table class="tbl"><thead><tr><th>r (px)</th><th>d meas (nm)</th><th>hkl</th><th>m</th><th>d calc (nm)</th><th>Δ (%)</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div class="row" style="margin-top:8px"><button class="btn sm" id="dfCalib">Use as calibration</button><button class="btn sm" id="dfCsv">Export CSV</button><button class="btn sm" id="dfPng">Export PNG</button><button class="btn sm" id="dfJson">Export JSON</button></div></div>`; }
function manualCard(){ if(S.type!=="spots") return ""; const L=curLat(); if(!L) return ""; const m=S.mres;
  return `<div class="card"><h3>Manual check <span class="hint">REDI 수동 해석처럼 — 상과 R1·R2 의 hkl 을 넣어 비교</span></h3>
    <div class="row"><label class="f" style="flex:1;min-width:150px">Phase<select id="dfMPh">${(S.recList||[]).map(r=>`<option value="${E(r.id)}" ${S.mPh===r.id?"selected":""}>${E(r.name)} · ${E(r.formula||"")} · ${E(r.sg&&r.sg.symbol||"")}</option>`).join("")}</select></label>
      <label class="f">R1 (h k l)<input type="text" id="dfMh1" value="${E(P.manual.h1)}" style="width:90px"></label><label class="f">R2 (h k l)<input type="text" id="dfMh2" value="${E(P.manual.h2)}" style="width:90px"></label>
      <button class="btn sm" id="dfMGo">Compare</button><button class="btn sm" id="dfMSug" title="측정 d 에 맞는 hkl 쌍 찾기">Suggest hkl</button></div>${m?manualResult(m):""}</div>`; }
function manualResult(m){ if(m.err) return `<div class="note warn" style="margin-top:8px">${E(m.err)}</div>`;
  if(m.sug) return `<div class="tblwrap" style="margin-top:8px;max-height:30vh"><table class="tbl"><thead><tr><th>R1</th><th>R2</th><th>Zone</th><th>d1 / d2 (nm)</th><th>∠ calc</th><th>Δ∠</th></tr></thead><tbody>${m.sug.map(s=>`<tr><td class="mono">${hklS(s.h1)}</td><td class="mono">${hklS(s.h2)}</td><td>${uvwS(s.z)}</td><td class="num">${f(s.d1,4)} / ${f(s.d2,4)}</td><td class="num">${f(s.ang,2)}°</td><td class="num">${f(s.dAng,2)}</td></tr>`).join("")||`<tr><td colspan="6">맞는 쌍이 없습니다</td></tr>`}</tbody></table></div>`;
  return `<div class="tblwrap" style="margin-top:8px"><table class="tbl"><thead><tr><th class="l"></th><th>d calc (nm)</th><th>d meas (nm)</th><th>Δ (%)</th></tr></thead><tbody>${m.rows.map(r=>`<tr><td class="l"><b>${r.k}</b> <span class="mono">${hklS(r.h)}</span></td><td class="num">${f(r.dc,4)}</td><td class="num">${r.dm?f(r.dm,4):"–"}</td><td class="num">${r.dm?f((r.dm-r.dc)/r.dc*100,2):"–"}</td></tr>`).join("")}</tbody></table></div>
    <div class="df-kv" style="margin-top:6px"><span>∠R1·R2</span><b>${f(m.ang,2)}° (측정 ${f(m.angM,2)}°)</b><span>Zone</span><b>${uvwS(m.zone)}</b>${m.kappa?`<span>κ</span><b>${f(m.kappa,4)}</b>`:""}</div>${m.note?`<div class="note warn" style="margin-top:6px">${E(m.note)}</div>`:""}`; }
function resultCards(){ if(!S.res) return ""; if(S.type==="spots") return (S.res.spots?latCard():"")+idCard()+manualCard(); if(S.type==="rings") return ringCard()+idCard(); return ""; }

/* ---------- 링 프로파일 ---------- */
function drawProfile(){ const host=$q("#dfProfile"); if(!host||!S.res||!S.res.rings) return; host.innerHTML=`<div class="card"><h3>Radial Profile <span class="hint">타원 보정 방위 평균</span></h3><canvas id="dfProf" style="width:100%;height:220px;display:block"></canvas></div>`;
  const cv=$q("#dfProf"); const dpr=window.devicePixelRatio||1; const W=cv.clientWidth, H=cv.clientHeight; cv.width=W*dpr; cv.height=H*dpr; const ctx=cv.getContext("2d"); ctx.setTransform(dpr,0,0,dpr,0,0);
  const R=S.res.rings, pr=R.profile; const cal=calib(); const s=ringScale(cal); const xs=pr.r.map(r=>s?r*s:r); const ys=pr.I.map((v,i)=>v-R.background[i]); const r0=Math.max(R.rStop||0,0)*1.1; let ymax=0; for(let i=0;i<ys.length;i++) if(pr.r[i]>r0&&ys[i]>ymax) ymax=ys[i];
  const x0=36,x1=W-8,y0=8,y1=H-26; const xmax=xs[xs.length-1]||1; const X=x=>x0+(x1-x0)*x/xmax, Y=y=>y1-(y1-y0)*Math.max(0,y)/(ymax||1);
  const css=getComputedStyle(document.body); ctx.strokeStyle=css.getPropertyValue("--line")||"#ccc"; ctx.lineWidth=1; ctx.beginPath(); ctx.moveTo(x0,y1); ctx.lineTo(x1,y1); ctx.stroke();
  ctx.fillStyle=css.getPropertyValue("--muted")||"#666"; ctx.font="11px system-ui"; for(let k=0;k<=5;k++){ const xv=xmax*k/5; ctx.fillText(s?xv.toFixed(1):xv.toFixed(0),X(xv)-8,H-10); } ctx.fillText(s?"g (1/nm)":"r (px)",x1-50,H-10);
  ctx.strokeStyle=css.getPropertyValue("--primary")||"#3a7"; ctx.lineWidth=1.3; ctx.beginPath(); let first=true; for(let i=0;i<xs.length;i++){ if(pr.r[i]<r0) continue; const X1=X(xs[i]), Y1=Y(ys[i]); if(first){ ctx.moveTo(X1,Y1); first=false; } else ctx.lineTo(X1,Y1); } ctx.stroke();
  const C=S.rcands&&S.rcands[S.rsel]; if(C){ ctx.strokeStyle="rgba(220,60,160,.8)"; ctx.fillStyle="rgba(220,60,160,1)"; for(const q of C.rings){ const xv=s?q.g*s/C.scale:q.g/C.scale; if(xv>xmax) continue; ctx.beginPath(); ctx.moveTo(X(xv),y1); ctx.lineTo(X(xv),y1-14); ctx.stroke(); ctx.save(); ctx.translate(X(xv)+3,y1-16); ctx.rotate(-Math.PI/2); ctx.fillText(q.hkl.map(Math.abs).sort((a,b)=>b-a).join(""),0,0); ctx.restore(); } } }

/* ---------- 손 조작 ---------- */
function wire(body){ const cv=$q("#dfCv",body);
  $q("#dfOpen",body).onclick=openFile;
  $all("#dfTool button",body).forEach(b=>b.onclick=()=>{ S.tool=b.dataset.t; $all("#dfTool button",body).forEach(x=>x.classList.toggle("active",x===b)); $q("#dfHint",body).textContent=toolHint(); });
  $q("#dfFit",body).onclick=()=>{ fitView(cv); draw(); };
  $q("#dfMaskCk",body).onchange=e=>{ P.showMask=e.target.checked; save(); draw(); }; $q("#dfPredCk",body).onchange=e=>{ P.showPred=e.target.checked; save(); draw(); };
  $q("#dfGamma",body).onchange=e=>{ S.gamma=e.target.value; buildDisplay(); draw(); };
  $all("#dfCal button",body).forEach(b=>b.onclick=async()=>{ P.calSrc=b.dataset.c; save(); if(P.calSrc==="saved") await loadCals(); rerender(); identify(); });
  const num=(id,fn)=>{ const el=$q(id,body); if(el) el.onchange=()=>{ fn(parseFloat(el.value)); save(); rerender(); identify(); }; };
  num("#dfBarLen",v=>S.barLen=v>0?v:null); num("#dfBarVal",v=>{ S.barVal=v>0?v:null; S.barFromOCR=false; }); num("#dfManScale",v=>S.manualScale=v>0?v:null);
  num("#dfKappa",v=>{ if(v>0) P.kappaPct=v; }); num("#dfTol",v=>{ if(v>0) P.tolPct=v; }); num("#dfDist",v=>{ if(v>=0) P.distPct=v; });
  const rk=$q("#dfRankCal",body); if(rk) rk.onchange=()=>{ P.rankCal=rk.value; save(); rerender(); identify(); };
  const sv=$q("#dfSaved",body); if(sv) sv.onchange=()=>{ S.savedCal=(S.calList||[])[+sv.value]||null; rerender(); identify(); };
  $all("#dfPat button",body).forEach(b=>b.onclick=()=>{ P.pattern=b.dataset.p; save(); if(S.res&&P.pattern!=="auto"&&((P.pattern==="spots"&&S.res.spots)||(P.pattern==="rings"&&S.res.rings))){ S.type=P.pattern; rerender(); identify(); } else analyze(); });
  $q("#dfRun",body).onclick=()=>{ S.center=null; analyze(); };
  const tn=$q("#dfTopN",body); if(tn) tn.onchange=()=>{ P.topN=+tn.value; save(); refit(); };
  const rd=$q("#dfRad",body); if(rd) rd.onchange=()=>{ P.radial=rd.value==="on"?true:rd.value==="off"?false:"auto"; save(); refit(); };
  $all("#dfLat button",body).forEach(b=>b.onclick=()=>{ S.lat=+b.dataset.l; S.cands=null; rerender(); identify(); });
  const els=$q("#dfEls",body); if(els) els.onchange=()=>{ P.elements=els.value; save(); identify(); };
  const em=$q("#dfElMode",body); if(em) em.onchange=()=>{ P.elMode=em.value; save(); identify(); };
  const mi=$q("#dfMaxIdx",body); if(mi) mi.onchange=()=>{ P.maxIdx=+mi.value; save(); identify(); };
  const im=$q("#dfIntM",body); if(im) im.onchange=()=>{ P.intModel=im.value; save(); identify(); };
  $all("[data-layer]",body).forEach(c=>c.onchange=()=>{ P.layers[c.dataset.layer]=c.checked; save(); identify(); });
  const idb=$q("#dfId",body); if(idb) idb.onclick=identify;
  $all(".df-cand tbody tr",body).forEach(tr=>tr.onclick=()=>{ if(S.type==="spots") S.sel=+tr.dataset.i; else S.rsel=+tr.dataset.i; rerender(); });
  const cb=$q("#dfCalib",body); if(cb) cb.onclick=useAsCalibration;
  const cs=$q("#dfCsv",body); if(cs) cs.onclick=exportCSV; const pn=$q("#dfPng",body); if(pn) pn.onclick=exportPNG; const js=$q("#dfJson",body); if(js) js.onclick=exportJSON;
  const mg=$q("#dfMGo",body); if(mg){ ensureRecList(body); mg.onclick=()=>manualCompare(body); $q("#dfMSug",body).onclick=()=>manualSuggest(body); }
  // 보기: 마우스
  let drag=null; cv.addEventListener("wheel",e=>{ e.preventDefault(); const v=S.view; const rc=cv.getBoundingClientRect(); const mx=e.clientX-rc.left, my=e.clientY-rc.top; const ix=v.x+mx/v.zoom, iy=v.y+my/v.zoom; const k=Math.exp(-e.deltaY*0.0015); v.zoom=Math.max(0.05,Math.min(40,v.zoom*k)); v.x=ix-mx/v.zoom; v.y=iy-my/v.zoom; draw(); },{passive:false});
  cv.addEventListener("pointerdown",e=>{ if(e.pointerType==="touch") return; drag={x:e.clientX,y:e.clientY,vx:S.view.x,vy:S.view.y,moved:false}; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener("pointermove",e=>{ if(!drag) return; const dx=e.clientX-drag.x, dy=e.clientY-drag.y; if(Math.hypot(dx,dy)>3) drag.moved=true; if(drag.moved&&S.tool==="view"){ S.view.x=drag.vx-dx/S.view.zoom; S.view.y=drag.vy-dy/S.view.zoom; draw(); } });
  cv.addEventListener("pointerup",e=>{ if(!drag) return; const d=drag; drag=null; if(!d.moved){ const rc=cv.getBoundingClientRect(); tapAt(e.clientX-rc.left,e.clientY-rc.top); } });
  cv.addEventListener("dblclick",()=>{ fitView(cv); draw(); });
  if(window.INNO_TOUCH){ if(S._detach) S._detach(); S._detach=INNO_TOUCH.attach(cv,{getView:()=>Object.assign({},S.view),setView:(v0,g)=>{ S.view=INNO_TOUCH.affine(v0,g,{min:0.05,max:40}); draw(); },oneFinger:"pan",
    onTap:p=>tapAt(p.x,p.y),onDoubleTap:()=>{ fitView(cv); draw(); },doubleTap:()=>S.tool==="view"}); }
  if(!S._ro&&window.ResizeObserver){ S._ro=new ResizeObserver(()=>draw()); } if(S._ro){ S._ro.disconnect(); S._ro.observe(cv); } }
/* 다시 그리기는 한 박자 뒤에 한 번만 — change 사건(blur 중)에 입력 칸을 지우면 브라우저 오류가 난다 */
let RR=0; function rerender(){ if(RR) return; RR=setTimeout(()=>{ RR=0; const b=$q("#dfBody"); if(b&&$q("#dfCv",b)!==undefined&&P.mode==="saed") renderSaed(b); },0); }
function tapAt(sx,sy){ if(!S||!S.view||!S.res) return; const ix=S.view.x+sx/S.view.zoom, iy=S.view.y+sy/S.view.zoom;
  if(S.tool==="center"){ if(S.type==="rings"){ S.ringCenter={x:ix,y:iy}; } S.center={x:ix,y:iy}; S.centerHint=true; S.tool="view"; if(S.type==="spots") refit(); else analyzeRingsAt(ix,iy); return; }
  if(S.tool==="spots"&&S.type==="spots"){ const list=(S.spotsEdit||S.res.spots.list).slice(); const tol=Math.max(6,10/S.view.zoom); let bi=-1,bd=1e9; list.forEach((s,i)=>{ const d=Math.hypot(s.x-ix,s.y-iy); if(d<bd){ bd=d; bi=i; } });
    if(bi>=0&&bd<=tol+(list[bi].radius||0)){ list.splice(bi,1); } else { const q=localPeak(ix,iy,Math.max(4,8/S.view.zoom)); if(q) list.push(q); }
    S.spotsEdit=list; refit(); return; }
  if(S.type==="spots"&&S.res.spots){ const L=curLat(); if(!L) return; let b=null,bd=1e9; for(const s of L.spots){ const d=Math.hypot(s.x-ix,s.y-iy); if(d<bd){ bd=d; b=s; } } if(b&&bd<12/S.view.zoom+6){ const C=S.cands&&S.cands[S.sel]; const h=C&&b.indexed?[0,1,2].map(i=>b.n[0]*C.h1[i]+b.n[1]*C.h2[i]):null; const g=Math.hypot(b.x-L.center.x,b.y-L.center.y); const cal=calib();
      const msg=`스폿 (${b.x.toFixed(1)}, ${b.y.toFixed(1)}) · |g| ${g.toFixed(1)} px${cal.scale?` · d ${(1/gOf([b.x-L.center.x,b.y-L.center.y],cal)).toFixed(4)} nm`:""}${b.indexed?` · n = (${b.n.join(", ")})`:" · 지수 없음"}${h?` · hkl ${hklT(h)}`:""}`; say(msg,"ok"); } } }
async function analyzeRingsAt(x,y){ // 지정 중심에서 링 다시 (analyze 의 ringCenter)
  S.busy=true; paint(); try{ const copy=new Float32Array(S.img); const r=await DW().run("analyze",{img:copy,w:S.w,h:S.h,opt:{pattern:"rings",ringCenter:{x,y}}},[copy.buffer]); S.res.rings=r.rings; S.res.mask=r.mask; S.maskCv=null; }catch(e){ S.err=String(e.message||e); } S.busy=false; paint(); identify(); }
function localPeak(x,y,r){ const {img,w,h}=S; let best=null; const R=Math.ceil(r); for(let dy=-R;dy<=R;dy++) for(let dx=-R;dx<=R;dx++){ const xx=Math.round(x)+dx, yy=Math.round(y)+dy; if(xx<1||yy<1||xx>=w-1||yy>=h-1) continue; const v=img[yy*w+xx]; if(!best||v>best.v) best={x:xx,y:yy,v}; }
  if(!best) return null; let sw=0,sx=0,sy=0,mn=Infinity; for(let dy=-3;dy<=3;dy++) for(let dx=-3;dx<=3;dx++){ const v=img[(best.y+dy)*w+best.x+dx]; if(v<mn) mn=v; } for(let dy=-2;dy<=2;dy++) for(let dx=-2;dx<=2;dx++){ const v=img[(best.y+dy)*w+best.x+dx]-mn; if(v>0){ sw+=v; sx+=v*(best.x+dx); sy+=v*(best.y+dy); } }
  return sw>0?{x:sx/sw,y:sy/sw,I:sw,snr:99,manual:true}:null; }
function paint(){ rerender(); }
async function openFile(){ const inp=document.createElement("input"); inp.type="file"; inp.accept=".dm3,.dm4,.tif,.tiff,.emd,.mrc"; inp.onchange=async()=>{ if(inp.files&&inp.files.length&&typeof loadFiles==="function"){ await loadFiles(inp.files); if(typeof state!=="undefined"){ state.view="lattice"; } if(typeof renderView==="function") renderView(); } }; inp.click(); }

/* ---------- 수동 비교 (REDI 수동 해석) ---------- */
async function ensureRecList(body){ if(S.recList) return; S.recList=await records(); if(!S.mPh&&S.recList.length) S.mPh=S.recList[0].id; const sel=$q("#dfMPh",body); if(sel&&!sel.options.length){ rerender(); } }
const parseH=s=>{ const m=String(s).trim().match(/^\(?\s*(-?\d+)[\s,]+(-?\d+)[\s,]+(-?\d+)\s*\)?$/)||String(s).trim().match(/^(-?\d)(-?\d)(-?\d)$/); return m?[+m[1],+m[2],+m[3]]:null; };
function mPhase(body){ const id=$q("#dfMPh",body).value; S.mPh=id; const r=(S.recList||[]).find(x=>x.id===id); return r?DC().phaseFromRecord(r):null; }
function manualCompare(body){ const h1=parseH($q("#dfMh1",body).value), h2=parseH($q("#dfMh2",body).value); P.manual={h1:$q("#dfMh1",body).value,h2:$q("#dfMh2",body).value}; save();
  const ph=mPhase(body); if(!ph||!h1||!h2){ S.mres={err:"상과 hkl 두 개(예: 1 1 1)를 넣으세요"}; rerender(); return; } const dc=DC(); const h3=[h1[0]+h2[0],h1[1]+h2[1],h1[2]+h2[2]]; const z=dc.zoneAxis(h1,h2);
  if(!z.some(Boolean)){ S.mres={err:"두 hkl 이 평행합니다"}; rerender(); return; }
  const L=correctedFit(curLat(),calib()); const s=calib().scale; const dm=i=>s?1/(L.len[i]*s):null; const ab=dc.absent(ph.ops,h1)||dc.absent(ph.ops,h2);
  const rows=[["R1",h1],["R2",h2],["R3",h3]].map(([k,h],i)=>({k,h,dc:dc.dSpacing(ph.cell,h),dm:dm(i)})); const ang=dc.angleG(ph.cell,h1,h2);
  S.mres={rows,ang,angM:L.ang12,zone:z,kappa:s?rows[0].dm/rows[0].dc:null,note:ab?"넣은 hkl 중 공간군 소멸 반사가 있습니다 (이중 회절이 아니면 맞지 않음)":""}; rerender(); }
function manualSuggest(body){ const ph=mPhase(body); const cal=calib(); const L=correctedFit(curLat(),cal); if(!ph){ return; } const dc=DC(); const tol=Math.max(cal.err,0.01);
  const refl=dc.reflections(ph,Math.max(20,(cal.scale?L.len[1]*cal.scale:0)*1.2),{allowed:true}); const out=[];
  if(cal.scale){ const g1=L.len[0]*cal.scale, g2=L.len[1]*cal.scale; const A=refl.filter(r=>Math.abs(r.g-g1)/g1<=tol), B=refl.filter(r=>Math.abs(r.g-g2)/g2<=tol);
    for(const a of A) for(const b of B){ const z=dc.zoneAxis(a.h,b.h); if(!z.some(Boolean)) continue; const ang=dc.angleG(ph.cell,a.h,b.h); const d=Math.abs(ang-L.ang12); if(d>3) continue; out.push({h1:a.h,h2:b.h,z:dc.reduceInt(z),d1:a.d,d2:b.d,ang,dAng:ang-L.ang12}); } }
  else { const ratio=L.len[1]/L.len[0]; for(const a of refl.slice(0,60)) for(const b of refl){ if(Math.abs(b.g/a.g-ratio)/ratio>tol) continue; const z=dc.zoneAxis(a.h,b.h); if(!z.some(Boolean)) continue; const ang=dc.angleG(ph.cell,a.h,b.h); if(Math.abs(ang-L.ang12)>2) continue; out.push({h1:a.h,h2:b.h,z:dc.reduceInt(z),d1:a.d,d2:b.d,ang,dAng:ang-L.ang12}); if(out.length>400) break; } }
  out.sort((p,q)=>Math.abs(p.dAng)-Math.abs(q.dAng)); const seen=new Set(), res=[]; for(const s of out){ const k=s.z.map(Math.abs).sort().join(); if(seen.has(k)) continue; seen.add(k); res.push(s); if(res.length>=15) break; }
  S.mres={sug:res}; rerender(); }

/* ---------- 보정 저장 ---------- */
async function loadCals(){ let list=[]; try{ if(window.DIFF_DB&&DIFF_DB.listCalibrations) list=await Promise.resolve(DIFF_DB.listCalibrations()); }catch(e){ list=[]; } if(!list||!list.length){ try{ list=JSON.parse(localStorage.getItem("inno.diff.cal")||"[]"); }catch(e){ list=[]; } } S.calList=list||[]; }
async function useAsCalibration(){ const m=S.meta; let cal=null, std="", stdId="";
  if(S.type==="spots"){ const c=S.cands&&S.cands[S.sel]; const L=curLat(); if(!c||!L) return; const A=DI().calibrateAffine(L.R1,L.R2,c.g[0],c.g[1],c.ang12); cal={scale:A.scale,ellip:A.ellip,phi:A.phi}; std=c.phaseName+" "+c.zone.join(""); stdId=c.phaseId; }
  else { const c=S.rcands&&S.rcands[S.rsel]; if(!c) return; const e=S.res.rings.ellipse; cal={scale:c.scale/Math.sqrt(e.a/e.b),ellip:e.a/e.b,phi:(e.phi+Math.PI/2)%Math.PI}; std=c.phaseName+" rings"; stdId=c.phaseId; }
  // 링: 반지름은 단축 기준 → 단축 방향 축척이 큼(smax = c.scale), 기하 평균 = c.scale/√e, 큰 축척 방향 = 단축 방향 (= 장축 φ + 90°)
  const name=prompt("보정 이름 (예: Libra200 375mm)",`${m.scope||"TEM"} ${m.camLen?m.camLen+" mm":""} ${m.kV?m.kV+" kV":""}`.trim()); if(name==null) return;
  const rec={id:"cal:"+Date.now().toString(36),name,microscope:m.scope||"",kV:m.kV||null,cameraLength:m.camLen||null,mode:m.opMode||"",scale:cal.scale,ellip:cal.ellip,phi:cal.phi,sigmaRel:null,standard:stdId,standardName:std,sourceFile:S.name,note:"",createdAt:new Date().toISOString()};
  let where="이 기기"; try{ if(window.DIFF_DB&&DIFF_DB.saveCalibration){ const st=DIFF_DB.firebaseStatus?await Promise.resolve(DIFF_DB.firebaseStatus()):null; const lab=!!(st&&st.isMember)&&confirm("연구실 공유 DB 에도 저장할까요?"); await DIFF_DB.saveCalibration(rec,{lab}); where=lab?"이 기기 + 연구실 DB":"이 기기"; } else { const L=JSON.parse(localStorage.getItem("inno.diff.cal")||"[]"); L.push(rec); localStorage.setItem("inno.diff.cal",JSON.stringify(L)); } }
  catch(e){ say("보정 저장 실패: "+(e.message||e),"err"); return; }
  say(`보정 저장 (${where}): ${cal.scale.toPrecision(5)} 1/nm/px · 찌그러짐 ${((cal.ellip-1)*100).toFixed(2)} %`,"ok"); await loadCals(); }

/* ---------- 내보내기 ---------- */
function dl(name,blob){ const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },1000); }
const base=()=>(S.name||"pattern").replace(/\.[^.]+$/,"").replace(/[^\w.-]+/g,"_");
function exportCSV(){ const cal=calib(); let s=""; if(S.type==="spots"){ const L=curLat(); const C=S.cands&&S.cands[S.sel]; s+="# INNO Analytics Diffraction — spots\n# file,"+S.name+"\n# calibration,"+(cal.scale||"")+",1/nm/px,"+cal.label+"\n";
    if(L) s+="# R1_px,"+L.len[0]+"\n# R2_px,"+L.len[1]+"\n# R3_px,"+L.len[2]+"\n# angle12_deg,"+L.ang12+"\n# center,"+L.center.x+","+L.center.y+"\n";
    if(C) s+="# phase,"+C.phaseName+",zone,"+C.zone.join(" ")+",ref,"+JSON.stringify(refOf(C.phaseId)).replace(/,/g,";")+"\n";
    s+="x_px,y_px,intensity,snr,n1,n2,indexed,h,k,l,d_nm\n"; for(const p of (L?L.spots:S.res.spots.list)){ const h=C&&p.indexed?[0,1,2].map(i=>p.n[0]*C.h1[i]+p.n[1]*C.h2[i]):null; const g=L?Math.hypot(p.x-L.center.x,p.y-L.center.y):null; s+=[p.x.toFixed(3),p.y.toFixed(3),p.I.toFixed(1),isFinite(p.snr)?p.snr.toFixed(1):"sat",p.n?p.n[0]:"",p.n?p.n[1]:"",p.indexed?1:0,h?h[0]:"",h?h[1]:"",h?h[2]:"",L&&cal.scale?(1/gOf([p.x-L.center.x,p.y-L.center.y],cal)).toFixed(5):""].join(",")+"\n"; } }
  else { const C=S.rcands&&S.rcands[S.rsel]; s+="# INNO Analytics Diffraction — rings\n# file,"+S.name+"\n# center,"+S.res.rings.center.x+","+S.res.rings.center.y+",a/b,"+S.res.rings.ratio+"\n"; if(C) s+="# phase,"+C.phaseName+",scale_1/nm/px,"+C.scale+"\n";
    s+="r_px,g_1/nm,d_nm,fwhm_px,snr,flag,hkl,d_calc_nm\n"; const sc=C?C.scale:ringScale(cal); for(const p of S.res.rings.peaks){ const q=C&&C.pairs.find(x=>Math.abs(x.r-p.r)<1e-6); s+=[p.r.toFixed(3),sc?(p.r*sc).toFixed(5):"",sc?(1/(p.r*sc)).toFixed(5):"",p.fwhm.toFixed(2),p.snr.toFixed(1),p.nearStop?"beamstop":p.broad?"broad":"",q?q.hkl.map(Math.abs).sort((a,b)=>b-a).join(""):"",q?(1/q.gPred).toFixed(5):""].join(",")+"\n"; } }
  dl(base()+"_diffraction.csv",new Blob([s],{type:"text/csv"})); }
function exportPNG(){ const cv=$q("#dfCv"); if(!cv) return; cv.toBlob(b=>{ if(b) dl(base()+"_diffraction.png",b); }); }
function exportJSON(){ const cal=calib(); const out={app:"INNO Analytics Diffraction",version:1,file:S.name,size:[S.w,S.h],meta:S.meta,calibration:cal,type:S.type};
  if(S.type==="spots"){ out.lattice=curLat(); const C=S.cands&&S.cands[S.sel]; if(C) out.selected=Object.assign({},C,{ref:refOf(C.phaseId)}); out.candidates=(S.cands||[]).slice(0,10).map(c=>({phase:c.phaseName,ref:refOf(c.phaseId),zone:c.zone,h1:c.h1,h2:c.h2,fom:c.fom,rho:c.rho,forb:c.forb,miss:c.missFrac,score:c.score,kappa:c.kappa,scale:c.scaleImplied})); }
  else { out.rings={center:S.res.rings.center,ellipse:S.res.rings.ellipse,peaks:S.res.rings.peaks}; const C=S.rcands&&S.rcands[S.rsel]; if(C) out.selected={phase:C.phaseName,ref:refOf(C.phaseId),scale:C.scale,pairs:C.pairs}; }
  dl(base()+"_diffraction.json",new Blob([JSON.stringify(out,null,1)],{type:"application/json"})); }

window.DIFF_UI={render,_S:()=>S,_P:P,_identify:identify,_analyze:analyze};
})();
