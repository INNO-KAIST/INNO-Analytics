/* ============================================================
   보기: LUT · 표시 버퍼 · 확대/이동 · 겹쳐 그리기(threshold·입자·ROI·축척 막대) · 도구
   ============================================================ */
/* ---------- LUT (ij/plugin/LutLoader 의 fire·ice·spectrum, 원색, HiLo, 지각 균일 색표) ---------- */
const LUTS=["Grays","Inverted","Fire","Ice","Spectrum","HiLo","viridis","magma","inferno","Red","Green","Blue","Cyan","Magenta","Yellow","3-3-2 RGB"];
const _lutCache={};
function interp32(r,g,b){ const n=r.length, out=new Uint8Array(768), sc=n/256; for(let i=0;i<256;i++){ const i1=Math.trunc(i*sc); let i2=i1+1; if(i2===n) i2=n-1; const f=i*sc-i1; out[i]=Math.trunc((1-f)*r[i1]+f*r[i2]); out[256+i]=Math.trunc((1-f)*g[i1]+f*g[i2]); out[512+i]=Math.trunc((1-f)*b[i1]+f*b[i2]); } return out; }
function anchors(A){ const out=new Uint8Array(768), n=A.length-1; for(let i=0;i<256;i++){ const t=i/255*n, k=Math.min(n-1,Math.floor(t)), f=t-k; for(let c=0;c<3;c++) out[c*256+i]=Math.round(A[k][c]*(1-f)+A[k+1][c]*f); } return out; }
function hsb(h){ const i=Math.floor(h*6), f=h*6-i, q=1-f; const r=[1,q,0,0,f,1][i%6], g=[f,1,1,q,0,0][i%6], b=[0,0,f,1,1,q][i%6]; return [Math.round(r*255),Math.round(g*255),Math.round(b*255)]; }
function lutTable(name){
  if(_lutCache[name]) return _lutCache[name]; let t=new Uint8Array(768);
  const prim={Red:[1,0,0],Green:[0,1,0],Blue:[0,0,1],Cyan:[0,1,1],Magenta:[1,0,1],Yellow:[1,1,0]};
  if(name==="Fire") t=interp32([0,0,1,25,49,73,98,122,146,162,173,184,195,207,217,229,240,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255],[0,0,0,0,0,0,0,0,0,0,0,0,0,14,35,57,79,101,117,133,147,161,175,190,205,219,234,248,255,255,255,255],[0,61,96,130,165,192,220,227,210,181,151,122,93,64,35,5,0,0,0,0,0,0,0,0,0,0,0,35,98,160,223,255]);
  else if(name==="Ice") t=interp32([0,0,0,0,0,0,19,29,50,48,79,112,134,158,186,201,217,229,242,250,250,250,250,251,250,250,250,250,251,251,243,230],[156,165,176,184,190,196,193,184,171,162,146,125,107,93,81,87,92,97,95,93,93,90,85,69,64,54,47,35,19,0,4,0],[140,147,158,166,170,176,209,220,234,225,236,246,250,251,250,250,245,230,230,222,202,180,163,142,123,114,106,94,84,64,26,27]);
  else if(name==="Spectrum"){ for(let i=0;i<256;i++){ const c=hsb(i/255); t[i]=c[0]; t[256+i]=c[1]; t[512+i]=c[2]; } }
  else if(name==="viridis") t=anchors([[68,1,84],[71,44,122],[59,81,139],[44,113,142],[33,144,141],[39,173,129],[92,200,99],[170,220,50],[253,231,37]]);
  else if(name==="magma") t=anchors([[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],[229,80,100],[251,135,97],[254,194,135],[252,253,191]]);
  else if(name==="inferno") t=anchors([[0,0,4],[31,12,72],[85,15,109],[136,34,106],[186,54,85],[227,89,51],[249,140,10],[249,201,50],[252,255,164]]);
  else if(name==="3-3-2 RGB"){ for(let i=0;i<256;i++){ t[i]=i&0xe0; t[256+i]=(i<<3)&0xe0; t[512+i]=(i<<6)&0xc0; } }
  else if(prim[name]){ const m=prim[name]; for(let i=0;i<256;i++){ t[i]=i*m[0]; t[256+i]=i*m[1]; t[512+i]=i*m[2]; } }
  else { const inv=name==="Inverted"; for(let i=0;i<256;i++){ const v=inv?255-i:i; t[i]=v; t[256+i]=v; t[512+i]=v; } if(name==="HiLo"){ t[0]=0; t[256]=0; t[512]=255; t[255]=255; t[511]=0; t[767]=0; } }
  _lutCache[name]=t; return t;
}
function lutCss(name){ const t=lutTable(name); const st=[]; for(let i=0;i<=8;i++){ const k=Math.round(i*255/8); st.push(`rgb(${t[k]},${t[256+k]},${t[512+k]}) ${i*12.5}%`); } return `linear-gradient(90deg,${st.join(",")})`; }

/* ---------- 표시 버퍼 ---------- */
const DISP={key:null,cv:null,scale:1};
const MAXDISP=24e6;
function invalidate(im){ if(im) im.ver=(im.ver||0)+1; DISP.key=null; DISPZ.key=null; OVL.key=null; }
function dispKey(im){ return [im.id,im.ver,im.c,im.z,im.t,im.composite,JSON.stringify(im.disp),PREVIEW.key||""].join("|"); }
/* 값 → 0..255 (ImageJ create8BitImage 와 같은 규칙)
   8-bit : LUT[i] = i<min ? 0 : i>max ? 255 : (int)(256·(i−min)/(max−min))   (ByteProcessor.setMinAndMax, min·max 반올림)
   16-bit: (int)(max(0,v−min)·256/(max−min+1)+0.5) ≤ 255              (ShortProcessor, min·max 정수)
   32-bit: (int)(max(0,v−min)·255/(max−min)+0.5) ≤ 255, NaN → 0          (FloatProcessor; max=min 이면 v>min → 255) */
function qMap(bit,mn,mx){
  if(bit===8||bit===16){ const N=bit===8?256:65536, t=new Uint8Array(N);
    if(bit===8){ const a=Math.round(mn), b=Math.round(mx); for(let i=0;i<256;i++){ t[i]=i<a?0:i>b?255:(b>a?Math.min(255,Math.max(0,Math.floor(256*(i-a)/(b-a)))):0); } }
    else { const a=Math.trunc(mn), b=Math.trunc(mx), sc=256/(b-a+1); for(let i=0;i<N;i++){ let v=i-a; if(v<0) v=0; v=Math.floor(v*sc+0.5); t[i]=v>255?255:v; } }
    return v=>v!==v?0:t[v<0?0:v>=N?N-1:v|0]; }
  const sc=255/(mx-mn); return v=>{ if(v!==v) return 0; let x=v-mn; if(x<0) x=0; const q=Math.floor(x*sc+0.5); return q>255||q!==q?(x>0?255:0):q; };
}
/* 이미지 영역 R={x,y,w,h,sc} 를 표시 색으로 (sc>1 이면 sc 칸마다 한 화소) — 화면·썸네일·PNG 내보내기 공용 */
function renderRegion(im,R,cvIn){
  const W=Math.max(1,Math.ceil(R.w/R.sc)), H=Math.max(1,Math.ceil(R.h/R.sc)); let cv=cvIn; if(!cv||cv.width!==W||cv.height!==H){ cv=document.createElement("canvas"); cv.width=W; cv.height=H; }
  const cx=cv.getContext("2d"), id=cx.createImageData(W,H), o=id.data;
  const chans=im.rgb?[0,1,2]:(im.composite&&im.nc>1?[...Array(im.nc).keys()]:[im.c]);
  const add=chans.length>1&&!im.rgb;
  for(let i=0;i<o.length;i+=4) o[i+3]=255;
  const bitOf=im.rgb?8:im.bit;
  chans.forEach(c=>{
    const prev=PREVIEW.plane&&PREVIEW.imId===im.id&&PREVIEW.planes?PREVIEW.planes[c]||null:null;
    const d=prev||planesOf(im)[pIndex(im,c,im.z,im.t)]; const ds=im.rgb?im.disp[0]:im.disp[c];
    const q8=im.rgb?(v=>v<0?0:v>255?255:v|0):qMap(bitOf,ds.min,ds.max);
    const L=im.rgb?null:lutTable(add?ds.lut:(ds.lut||"Grays")); const hilo=!im.rgb&&ds.lut==="HiLo"; const mn=ds.min, mx=ds.max;
    for(let y=0;y<H;y++){ const sy=R.y+y*R.sc; if(sy>=im.h) break; const ro=sy*im.w;
      for(let x=0;x<W;x++){ const sx=R.x+x*R.sc; if(sx>=im.w) break; const v=d[ro+sx]; const q=q8(v); const j=(y*W+x)*4;
        if(im.rgb){ o[j+c]=q; continue; }
        if(hilo&&v===v&&(v<=mn||v>=mx)){ const r=v>=mx?255:0, b=v<=mn?255:0; if(add){ o[j]=Math.min(255,o[j]+r); o[j+2]=Math.min(255,o[j+2]+b); } else { o[j]=r; o[j+1]=0; o[j+2]=b; } continue; }
        if(add){ o[j]=Math.min(255,o[j]+L[q]); o[j+1]=Math.min(255,o[j+1]+L[256+q]); o[j+2]=Math.min(255,o[j+2]+L[512+q]); }
        else { o[j]=L[q]; o[j+1]=L[256+q]; o[j+2]=L[512+q]; } } }
  });
  cx.putImageData(id,0,0); return cv;
}
function fullScale(im){ return im.w*im.h>MAXDISP?Math.ceil(Math.sqrt(im.w*im.h/MAXDISP)):1; }
function buildDisplay(im){
  const key=dispKey(im); if(DISP.key===key&&DISP.cv) return DISP;
  const sc=fullScale(im); DISP.cv=renderRegion(im,{x:0,y:0,w:im.w,h:im.h,sc},DISP.cv); DISP.key=key; DISP.scale=sc; return DISP;
}
/* 아주 큰 이미지(> 2.4×10⁷ 화소)를 확대했을 때: 보이는 영역만 원래 해상도로 (여유 25 %, 영역 안에서 움직이면 다시 그리지 않음) */
const DISPZ={key:null,cv:null,R:null};
function zoomRegion(im,v,W,H){
  const sc=fullScale(im); if(sc===1||v.zoom*sc<=1.0001) return null;
  const [ax,ay]=toImg(v,0,0), [bx,by]=toImg(v,W,H); const vx0=Math.max(0,Math.floor(ax)), vy0=Math.max(0,Math.floor(ay)), vx1=Math.min(im.w,Math.ceil(bx)), vy1=Math.min(im.h,Math.ceil(by)); if(vx1<=vx0||vy1<=vy0) return null;
  const key=dispKey(im); const R0=DISPZ.R;
  if(DISPZ.key===key&&R0&&vx0>=R0.x&&vy0>=R0.y&&vx1<=R0.x+R0.w&&vy1<=R0.y+R0.h) return DISPZ;
  const mx=Math.ceil((vx1-vx0)*0.25), my=Math.ceil((vy1-vy0)*0.25); let x=Math.max(0,vx0-mx), y=Math.max(0,vy0-my), w=Math.min(im.w,vx1+mx)-x, h=Math.min(im.h,vy1+my)-y;
  if(w*h>MAXDISP) return null; const R={x,y,w,h,sc:1}; DISPZ.cv=renderRegion(im,R,DISPZ.cv); DISPZ.R=R; DISPZ.key=key; return DISPZ;
}
/* 미리 보기(필터 매개변수를 바꾸는 동안 결과를 화면에만) */
const PREVIEW={key:"",imId:null,planes:null,plane:false,c:0};
function clearPreview(){ if(PREVIEW.plane){ PREVIEW.plane=false; PREVIEW.planes=null; PREVIEW.key=""; DISP.key=null; draw(); } }

/* ---------- 보기 변환 ---------- */
function viewOf(im){ let v=S.ui.views[im.id]; if(!v){ v=S.ui.views[im.id]={zoom:0,x:0,y:0}; } return v; }
function canvasSize(){ const el=ROOT&&$q(ROOT,"#imCanvasWrap"); if(!el) return {W:800,H:600}; const r=el.getBoundingClientRect(); return {W:Math.max(200,Math.floor(r.width)),H:Math.max(200,Math.floor(r.height))}; }
function fitView(im){ const {W,H}=canvasSize(); const z=Math.min(W/im.w,H/im.h); const v=viewOf(im); v.zoom=z; v.x=(im.w-W/z)/2; v.y=(im.h-H/z)/2; }
function resetView(im){ const v=viewOf(im); v.zoom=0; }
const ZSTEPS=[1/72,1/48,1/32,1/24,1/16,1/12,1/8,1/6,1/4,1/3,1/2,0.75,1,1.5,2,3,4,6,8,12,16,24,32,48,64];
function zoomAt(im,f,sx,sy){ const v=viewOf(im); const ix=v.x+sx/v.zoom, iy=v.y+sy/v.zoom; let z=v.zoom*f; z=clamp(z,1/72,64); v.zoom=z; v.x=ix-sx/z; v.y=iy-sy/z; draw(); }
function zoomStep(im,dir,sx,sy){ const v=viewOf(im); const {W,H}=canvasSize(); if(sx==null){ sx=W/2; sy=H/2; } let z=v.zoom;
  const nz=dir>0?(ZSTEPS.find(s=>s>z*1.001)||64):([...ZSTEPS].reverse().find(s=>s<z/1.001)||1/72); zoomAt(im,nz/z,sx,sy); }
const toImg=(v,sx,sy)=>[v.x+sx/v.zoom,v.y+sy/v.zoom];
const toScr=(v,ix,iy)=>[(ix-v.x)*v.zoom,(iy-v.y)*v.zoom];

/* ---------- 그리기 ---------- */
let DRAW_RAF=0; function drawSoon(){ if(DRAW_RAF) return; DRAW_RAF=requestAnimationFrame(()=>{ DRAW_RAF=0; draw(); }); } // 막대를 끄는 동안 한 프레임에 한 번만
const OVL={key:null,thr:null};
function draw(){
  const im=cur(); const cv=ROOT&&$q(ROOT,"#imCanvas"), ov=ROOT&&$q(ROOT,"#imOverlay"); if(!cv||!ov) return;
  const {W,H}=canvasSize(); const dpr=window.devicePixelRatio||1;
  [cv,ov].forEach(c=>{ if(c.width!==Math.round(W*dpr)||c.height!==Math.round(H*dpr)){ c.width=Math.round(W*dpr); c.height=Math.round(H*dpr); c.style.width=W+"px"; c.style.height=H+"px"; } });
  const cx=cv.getContext("2d"); cx.setTransform(dpr,0,0,dpr,0,0); cx.fillStyle=getComputedStyle(cv).getPropertyValue("--imbg")||"#111"; cx.fillRect(0,0,W,H);
  const ox=ov.getContext("2d"); ox.setTransform(dpr,0,0,dpr,0,0); ox.clearRect(0,0,W,H);
  if(!im){ return; }
  const v=viewOf(im); if(!v.zoom) fitView(im);
  const Z=zoomRegion(im,v,W,H);
  if(Z){ cx.imageSmoothingEnabled=false; const [x0,y0]=toScr(v,Z.R.x,Z.R.y); cx.drawImage(Z.cv,x0,y0,Z.R.w*v.zoom,Z.R.h*v.zoom); }
  else { const D=buildDisplay(im); cx.imageSmoothingEnabled=v.zoom<1; cx.imageSmoothingQuality="high"; const [x0,y0]=toScr(v,0,0); cx.drawImage(D.cv,x0,y0,im.w*v.zoom,im.h*v.zoom); }
  // 화소 격자 (아주 크게 확대했을 때)
  if(v.zoom>=24){ ox.strokeStyle="rgba(255,255,255,.08)"; ox.lineWidth=1; ox.beginPath(); const [ax,ay]=toImg(v,0,0), [bx,by]=toImg(v,W,H);
    for(let x=Math.max(0,Math.floor(ax));x<=Math.min(im.w,Math.ceil(bx));x++){ const s=toScr(v,x,0)[0]; ox.moveTo(s,0); ox.lineTo(s,H); }
    for(let y=Math.max(0,Math.floor(ay));y<=Math.min(im.h,Math.ceil(by));y++){ const s=toScr(v,0,y)[1]; ox.moveTo(0,s); ox.lineTo(W,s); } ox.stroke(); }
  drawThreshold(im,ox,v); drawParticles(im,ox,v); drawDrop(im,ox,v); drawRois(im,ox,v); drawScaleBar(im,ox,v,W,H);
  drawStatus(); drawProfileMini(); drawPolyBar();
}
function thrRange(im){ const T=S.ui.thr; return {lo:+T.lower,hi:+T.upper}; }
function drawThreshold(im,ox,v){
  if(S.ui.panel!=="threshold"||!S.ui.thr.show) return; if(im.rgb&&!ctReady(im)){ ctEnsure(im); return; } const {lo,hi}=thrRange(im);
  const key=[im.id,im.ver,im.c,im.z,im.t,lo,hi,S.ui.thr.localPreview||"",im.rgb?JSON.stringify(S.ui.cthr)+CT.key:""].join("|");
  if(OVL.key!==key){ const sc=im.w*im.h>MAXDISP?Math.ceil(Math.sqrt(im.w*im.h/MAXDISP)):1; const W=Math.ceil(im.w/sc),H=Math.ceil(im.h/sc);
    const c=document.createElement("canvas"); c.width=W; c.height=H; const g=c.getContext("2d"), idd=g.createImageData(W,H), o=idd.data;
    const lp=S.ui.thr.localPreview&&THR_LOCAL.mask&&THR_LOCAL.key===S.ui.thr.localPreview?THR_LOCAL.mask:null; const d=planesOf(im)[pIndex(im,im.c,im.z,im.t)];
    for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const i=(y*sc)*im.w+x*sc; const inn=im.rgb?ctInside(i):lp?lp[i]>0:(d[i]>=lo&&d[i]<=hi); if(inn){ const j=(y*W+x)*4; o[j]=255; o[j+1]=im.rgb?0:38; o[j+2]=im.rgb?0:38; o[j+3]=im.rgb?235:150; } } /* RGB 는 ImageJ Color Threshold 처럼 진한 빨강 (색 이미지 위에서도 보이게) */
    g.putImageData(idd,0,0); OVL.thr=c; OVL.key=key; }
  ox.save(); ox.imageSmoothingEnabled=false; const [x0,y0]=toScr(v,0,0); ox.drawImage(OVL.thr,x0,y0,im.w*v.zoom,im.h*v.zoom); ox.restore();
}
function drawParticles(im,ox,v){
  const pa=S.particles[im.id]; if(!pa||!pa.polys||S.ui.pa.show==="none"&&!S.ui.paPick) return; const sel=new Set(pa.sel||[]), pk=new Set(pa.pick||[]); const st=paStatus(pa);
  ox.save(); ox.lineWidth=1.2; ox.font="600 11px DM Sans, system-ui, sans-serif"; ox.textAlign="center"; ox.textBaseline="middle";
  pa.polys.forEach((pl,i)=>{ const [xs,ys]=pl; ox.beginPath(); for(let k=0;k<xs.length;k++){ const [sx,sy]=toScr(v,xs[k],ys[k]); if(k) ox.lineTo(sx,sy); else ox.moveTo(sx,sy); } ox.closePath();
    // 색: 사용 = 하늘색, 겹침 = 주황 점선, 크기 이상치 = 자홍 점선, 고른 입자 = 노랑 칠, 표에서 고른 행 = 노랑
    const s=sel.has(i), p=pk.has(i), x=st[i]; if(S.ui.pa.show==="masks"||s||p){ ox.fillStyle=p?"rgba(255,220,0,.5)":s?"rgba(255,200,0,.45)":x?"rgba(255,140,60,.22)":"rgba(0,200,255,.28)"; ox.fill(); }
    ox.setLineDash(x&&!p?[4,3]:[]); ox.lineWidth=p?2:1.2; ox.strokeStyle=p||s?"#ffd400":x==="ovl"?"#ff8c3c":x==="out"?"#ff4fd8":"#00d0ff"; ox.stroke(); });
  ox.setLineDash([]);
  if(v.zoom*Math.sqrt(pa.meanArea||25)>12&&S.ui.pa.show!=="masks"){ pa.cent.forEach((c,i)=>{ const [sx,sy]=toScr(v,c[0],c[1]); ox.lineWidth=3; ox.strokeStyle="rgba(0,0,0,.65)"; ox.strokeText(String(i+1),sx,sy); ox.fillStyle="#fff"; ox.fillText(String(i+1),sx,sy); }); }
  ox.restore();
}
/* 축척 막대 (Analyze ▸ Tools ▸ Scale Bar): 화면 폭의 약 1/5 에 맞는 1·2·5 × 10ⁿ 길이 */
function niceLen(x){ const p=Math.pow(10,Math.floor(Math.log10(x))); const m=x/p; return (m<1.5?1:m<3.5?2:m<7.5?5:10)*p; }
function scaleBarSpec(im,viewW){ const c=calOf(im); if(!im.cal) return null; const fov=viewW*c.pw; const auto=niceLen(fov/5); let L=S.ui.scalebar.len>0?+S.ui.scalebar.len:auto;
  if(L/c.pw>0.9*viewW) L=auto; // 지정 길이가 보이는 이미지보다 길면 자동 길이로
  return {L,px:L/c.pw,unit:c.unit,user:S.ui.scalebar.len>0&&L===+S.ui.scalebar.len}; }
/* 눈금 막대 글자: 값이 1–1000 이 되는 단위로 (0.02 µm → 20 nm, 2000 nm → 2 µm, Å → nm) */
function sbLabel(L,u){ const M={pm:1e-12,"Å":1e-10,nm:1e-9,"µm":1e-6,um:1e-6,mm:1e-3,cm:1e-2,m:1}; if(!(u in M)) return `${g4(L)} ${u}`; const m=L*M[u];
  for(const [nu,k] of [["m",1],["mm",1e-3],["µm",1e-6],["nm",1e-9],["Å",1e-10]]) if(m>=k*0.9999||nu==="Å") return `${g4(m/k)} ${nu}`; return `${g4(L)} ${u}`; }
/* 화면·내보내기 공용: 보이는 이미지 영역 안에 두고, 글자·막대 굵기는 그 영역 크기에 비례 */
function drawScaleBar(im,ox,v,W,H){
  if(!S.ui.scalebar.on||!im.cal) return;
  const [sx0,sy0]=toScr(v,0,0), [sx1,sy1]=toScr(v,im.w,im.h); const ax=Math.max(0,sx0), ay=Math.max(0,sy0), bx=Math.min(W,sx1), by=Math.min(H,sy1); const vw=bx-ax, vh=by-ay; if(vw<40||vh<30) return;
  const sb=scaleBarSpec(im,vw/v.zoom); if(!sb) return; const len=sb.px*v.zoom; if(len<8) return;
  const fs=Math.round(Math.max(12,Math.min(96,Math.min(vw,vh)/40))), th=Math.max(4,Math.round(fs/2.6)), pad=Math.round(fs*1.1);
  const pos=S.ui.scalebar.pos, x=pos.includes("l")?ax+pad:bx-pad-len, y=pos.includes("t")?ay+pad+fs+th+4:by-pad;
  ox.save(); ox.fillStyle="rgba(0,0,0,.45)"; ox.fillRect(x-fs/2,y-th-fs-8,len+fs,th+fs+8+fs/2); ox.fillStyle="#fff"; ox.fillRect(x,y-th,len,th);
  ox.font=`600 ${fs}px DM Sans, system-ui, sans-serif`; ox.textAlign="center"; ox.fillText(sbLabel(sb.L,sb.unit),x+len/2,y-th-6); ox.restore();
}

/* ---------- ROI ---------- */
const AREA_T=new Set(["rect","oval","poly","free","traced"]);
const roiList=im=>S.rois[im.id]||(S.rois[im.id]=[]);
const activeRoi=im=>S.roiSel[im.id]||null;
function setActive(im,r){ if(r) S.roiSel[im.id]=r; else delete S.roiSel[im.id]; save(); draw(); drawSide(); if(S.ui.panel==="analyze"&&typeof dropRefresh==="function") dropRefresh(); } // Analyze 안 Drop Shape 칸은 영역·기준선 ROI 를 보여 줌
function clearRois(im){ delete S.roiSel[im.id]; S.rois[im.id]=[]; delete S.particles[im.id]; if(S.dropFit) delete S.dropFit[im.id]; if(S.dropBase) delete S.dropBase[im.id]; } // 좌표가 바뀌는 연산 → 방울 덧그림·기준선도
function roiBounds(r){ const xs=r.pts.map(p=>p[0]), ys=r.pts.map(p=>p[1]); const x=Math.min(...xs), y=Math.min(...ys); return {x,y,w:Math.max(...xs)-x,h:Math.max(...ys)-y}; }
function roiRect(im){ const r=activeRoi(im); if(!r||r.type!=="rect") return null; return roiBounds(r); }
function roiGeom(im,r){ // {mask, rect, poly, kind} — IC.measureRoi 형식
  if(r.type==="rect"){ const b=roiBounds(r); const X0=Math.round(b.x), Y0=Math.round(b.y), X1=Math.round(b.x+b.w), Y1=Math.round(b.y+b.h);
    const x=clamp(X0,0,im.w), y=clamp(Y0,0,im.h); const w=clamp(X1,0,im.w)-x, h=clamp(Y1,0,im.h)-y;
    return {mask:null,rect:{x,y,w,h},bounds:{x:X0,y:Y0,w:X1-X0,h:Y1-Y0},poly:{xs:[X0,X1,X1,X0],ys:[Y0,Y0,Y1,Y1]},kind:"rect"}; }
  if(r.type==="oval"){ const b=roiBounds(r); const X0=Math.round(b.x), Y0=Math.round(b.y), BW=Math.round(b.w), BH=Math.round(b.h); const m=IC.ovalMask(X0,Y0,BW,BH,im.w,im.h); const n=64, xs=[], ys=[]; for(let k=0;k<n;k++){ const a=2*Math.PI*k/n; xs.push(X0+BW/2+BW/2*Math.cos(a)); ys.push(Y0+BH/2+BH/2*Math.sin(a)); } return {mask:m.mask,rect:m.rect,bounds:{x:X0,y:Y0,w:BW,h:BH},poly:{xs,ys},kind:"oval"}; }
  if(r.type==="poly"||r.type==="free"||r.type==="traced"){ const xs=r.pts.map(p=>p[0]), ys=r.pts.map(p=>p[1]); const m=IC.polygonMask(xs,ys,im.w,im.h); return {mask:m.mask,rect:m.rect,poly:{xs,ys},kind:r.type==="free"?"freehand":r.type==="traced"?"traced":"polygon"}; }
  return null;
}
function roisMask(im,includeManager){ const L=(includeManager?roiList(im):[]).concat(activeRoi(im)?[activeRoi(im)]:[]).filter(r=>AREA_T.has(r.type)); if(!L.length) return null; const m=new Uint8Array(im.w*im.h);
  for(const r of L){ const g=roiGeom(im,r); if(!g) continue; const R=g.rect; for(let y=0;y<R.h;y++) for(let x=0;x<R.w;x++) if(!g.mask||g.mask[y*R.w+x]) m[(y+R.y)*im.w+x+R.x]=1; } return m; }
function lineLen(im,r){ const c=calOf(im); let L=0; for(let k=1;k<r.pts.length;k++) L+=Math.hypot((r.pts[k][0]-r.pts[k-1][0])*c.pw,(r.pts[k][1]-r.pts[k-1][1])*c.ph); return L; }
function roiAngle(r){ if(r.type==="line"){ const [a,b]=r.pts; return Math.atan2(a[1]-b[1],b[0]-a[0])*180/Math.PI; } if(r.type==="angle"&&r.pts.length===3){ const [a,b,c]=r.pts; const a1=Math.atan2(a[1]-b[1],a[0]-b[0]), a2=Math.atan2(c[1]-b[1],c[0]-b[0]); let d=Math.abs(a1-a2)*180/Math.PI; if(d>180) d=360-d; return d; } return NaN; }
function pathOf(ox,v,r,closed){ ox.beginPath(); r.pts.forEach((p,k)=>{ const [sx,sy]=toScr(v,p[0],p[1]); if(k) ox.lineTo(sx,sy); else ox.moveTo(sx,sy); }); if(closed) ox.closePath(); }
function drawRoi(ox,v,r,act,im){
  ox.save(); const col=r.color||(act?"#ffd400":"#3fd2ff"); ox.strokeStyle=col; ox.lineWidth=act?1.6:1.2; ox.fillStyle=col;
  if(r.type==="rect"||r.type==="oval"){ const b=roiBounds(r); const [x,y]=toScr(v,b.x,b.y); const w=b.w*v.zoom, h=b.h*v.zoom; ox.beginPath(); if(r.type==="rect") ox.rect(x,y,w,h); else ox.ellipse(x+w/2,y+h/2,Math.abs(w/2),Math.abs(h/2),0,0,2*Math.PI); ox.stroke(); }
  else if(r.type==="poly"||r.type==="free"||r.type==="traced"){ pathOf(ox,v,r,!r.open); ox.stroke(); }
  else if(r.type==="line"||r.type==="angle"){ pathOf(ox,v,r,false); if(r.type==="line"&&(S.ui.lw||1)>1){ ox.save(); ox.globalAlpha=0.25; ox.lineWidth=(S.ui.lw||1)*v.zoom; ox.stroke(); ox.restore(); } ox.stroke(); }
  else if(r.type==="point"){ r.pts.forEach((p,k)=>{ const [sx,sy]=toScr(v,p[0],p[1]); ox.beginPath(); ox.arc(sx,sy,act?4.5:3.5,0,2*Math.PI); ox.fillStyle=r.cls!=null?CLS_COL[r.cls[k]%CLS_COL.length]:col; ox.fill(); ox.lineWidth=1; ox.strokeStyle="#000"; ox.stroke(); }); }
  if(act&&r.type!=="point"){ const hs=TUI()?5:3; handlesOf(r).forEach(([hx,hy])=>{ const [sx,sy]=toScr(v,hx,hy); ox.fillStyle="#fff"; ox.fillRect(sx-hs,sy-hs,2*hs,2*hs); ox.strokeStyle="#000"; ox.lineWidth=1; ox.strokeRect(sx-hs,sy-hs,2*hs,2*hs); }); } // 손가락이면 손잡이를 크게
  if(r.name&&!act){ const b=roiBounds(r); const [sx,sy]=toScr(v,b.x,b.y); ox.font="600 10.5px DM Sans,system-ui"; ox.fillStyle=col; ox.fillText(r.name,sx+2,sy-3); }
  void im; ox.restore();
}
const CLS_COL=["#ff4d4d","#3fd2ff","#ffd400","#7CFF6B","#ff7ad9","#ff9a3c","#9a8cff","#ffffff"];
function handlesOf(r){ if(r.type==="rect"||r.type==="oval"){ const b=roiBounds(r); return [[b.x,b.y],[b.x+b.w,b.y],[b.x+b.w,b.y+b.h],[b.x,b.y+b.h],[b.x+b.w/2,b.y],[b.x+b.w,b.y+b.h/2],[b.x+b.w/2,b.y+b.h],[b.x,b.y+b.h/2]]; } if(r.type==="free"||r.type==="traced") return []; return r.pts; }
function drawRois(im,ox,v){ if(S.ui.showAll!==false) roiList(im).forEach(r=>{ if(r!==activeRoi(im)) drawRoi(ox,v,r,false,im); }); const a=activeRoi(im); if(a) drawRoi(ox,v,a,true,im); if(DRAG&&DRAG.temp) drawRoi(ox,v,DRAG.temp,true,im); }

/* ---------- 도구·포인터 ---------- */
const TOOLS=[["rect","Rectangle","M3 5h18v14H3z"],["oval","Oval","M12 5c5 0 9 3 9 7s-4 7-9 7-9-3-9-7 4-7 9-7z"],["poly","Polygon","M4 18 7 5l9 3 4 9-9 3z"],["free","Freehand","M4 15c3-9 6 4 9-3s5-3 7 1"],
  ["line","Straight line (profile)","M4 19 20 5"],["angle","Angle","M4 19h16M4 19 15 6"],["point","Multi-point (count)","M7 8h.01M15 7h.01M10 15h.01M17 16h.01"],["wand","Wand (tracing)","M4 20 14 10M15 4v3M18 6l-2 2M20 9h-3"],
  ["scale","Set Scale (draw a line on the scale bar)","M3 17h18M3 13v8M21 13v8M8 17v-2M13 17v-2"],["hand","Scrolling tool (Space)","M8 13V6a1.5 1.5 0 0 1 3 0v6V4a1.5 1.5 0 0 1 3 0v8V6a1.5 1.5 0 0 1 3 0v8a6 6 0 0 1-6 6h-1a6 6 0 0 1-5-3l-2-4a1.5 1.5 0 0 1 2.5-1.5L8 13"],["zoom","Zoom (click = in, Alt+click = out)","M10 17a7 7 0 1 0 0-14 7 7 0 0 0 0 14zM21 21l-6-6M7 10h6M10 7v6"]];
let DRAG=null, SPACE=false, HOVER=null, LASTPT="mouse";
const TUI=()=>LASTPT==="touch"||!!(window.INNO_TOUCH&&INNO_TOUCH.touchUI()); // 손가락 화면: 손잡이·허용 거리를 크게
function evPos(e){ const cv=$q(ROOT,"#imOverlay"); const r=cv.getBoundingClientRect(); return [e.clientX-r.left,e.clientY-r.top]; }
function hitHandle(im,r,sx,sy,tol){ tol=tol||5; const v=viewOf(im); const H=handlesOf(r); for(let k=0;k<H.length;k++){ const [hx,hy]=toScr(v,H[k][0],H[k][1]); if(Math.abs(hx-sx)<=tol&&Math.abs(hy-sy)<=tol) return k; } return -1; }
function insideRoi(im,r,ix,iy){ if(r.type==="rect"||r.type==="oval"){ const b=roiBounds(r); if(r.type==="rect") return ix>=b.x&&ix<=b.x+b.w&&iy>=b.y&&iy<=b.y+b.h; const dx=(ix-b.x-b.w/2)/(b.w/2), dy=(iy-b.y-b.h/2)/(b.h/2); return dx*dx+dy*dy<=1; }
  if(r.type==="poly"||r.type==="free"||r.type==="traced"){ let c=false; const P=r.pts; for(let i=0,j=P.length-1;i<P.length;j=i++){ if((P[i][1]>iy)!==(P[j][1]>iy)&&ix<(P[j][0]-P[i][0])*(iy-P[i][1])/(P[j][1]-P[i][1])+P[i][0]) c=!c; } return c; }
  if(r.type==="line"){ const v=viewOf(im); const [a,b]=r.pts; const L=Math.hypot(b[0]-a[0],b[1]-a[1])||1; const t=clamp(((ix-a[0])*(b[0]-a[0])+(iy-a[1])*(b[1]-a[1]))/(L*L),0,1); const d=Math.hypot(ix-(a[0]+t*(b[0]-a[0])),iy-(a[1]+t*(b[1]-a[1]))); return d*v.zoom<(LASTPT==="touch"?12:5); }
  return false; }
const snapI=v=>Math.floor(v); // ImageJ: 화면 → 이미지 좌표는 (int) 내림
function onDown(e){
  const im=cur(); if(!im) return; const [sx,sy]=evPos(e); const v=viewOf(im); const [ix,iy]=toImg(v,sx,sy); const tool=S.ui.tool;
  if(e.button===2) return; // 오른쪽 단추는 그리지 않음
  LASTPT=e.pointerType||"mouse"; const touch=LASTPT==="touch"; if(touch&&!e._tap) touchSnap(im); // 손가락: 둘째 손가락이 오면 되돌릴 수 있게
  try{ $q(ROOT,"#imOverlay").setPointerCapture(e.pointerId); }catch(er){}
  if(e.button===1||SPACE||tool==="hand"){ DRAG={mode:"pan",sx,sy,vx:v.x,vy:v.y}; return; }
  if(tool==="zoom"){ zoomStep(im,e.altKey?-1:1,sx,sy); return; }
  if(S.ui.paPick&&S.particles[im.id]){ pickToggle(im,ix,iy); return; } // 입자 고르기 중: 누른 입자를 고르거나 뺌 (ROI 를 그리지 않음)
  const a=activeRoi(im);
  if(a&&tool!=="point"&&tool!=="wand"&&(tool!=="poly"||!(DRAG&&DRAG.mode==="poly"))){ const h=hitHandle(im,a,sx,sy,touch?14:5); if(h>=0){ DRAG={mode:"handle",h,orig:JSON.parse(JSON.stringify(a.pts)),ix,iy}; return; }
    if(insideRoi(im,a,ix,iy)&&!(e.shiftKey)){ DRAG={mode:"move",orig:JSON.parse(JSON.stringify(a.pts)),ix,iy}; return; } }
  if(tool==="poly"){ if(DRAG&&DRAG.mode==="poly"){ const t=DRAG.temp; const [fx,fy]=toScr(v,t.pts[0][0],t.pts[0][1]); if(t.pts.length>=3&&Math.hypot(fx-sx,fy-sy)<(touch?18:7)){ finishPoly(im); return; } const P=[snapI(ix),snapI(iy)]; t.pts[t.pts.length-1]=P; t.pts.push(P.slice()); draw(); return; } // 터치는 hover 가 없어 고무줄 점을 누른 자리로 먼저 옮김
    DRAG={mode:"poly",temp:{id:uid(),type:"poly",pts:[[snapI(ix),snapI(iy)],[snapI(ix),snapI(iy)]]}}; return; }
  if(tool==="angle"){ if(DRAG&&DRAG.mode==="angle"){ DRAG.temp.pts.push([ix,iy]); if(DRAG.temp.pts.length>=4){ DRAG.temp.pts.length=3; setActive(im,DRAG.temp); DRAG=null; } draw(); return; } DRAG={mode:"angle",temp:{id:uid(),type:"angle",pts:[[ix,iy],[ix,iy]]}}; return; }
  if(tool==="point"){ let r=a&&a.type==="point"?a:null; if(!r||!e.shiftKey&&!S.ui.pointAdd&&false){ r={id:uid(),type:"point",pts:[],cls:[]}; } if(!r) r={id:uid(),type:"point",pts:[],cls:[]};
    if(e.altKey){ let bi=-1,bd=8; r.pts.forEach((p,k)=>{ const [px,py]=toScr(v,p[0],p[1]); const d=Math.hypot(px-sx,py-sy); if(d<bd){ bd=d; bi=k; } }); if(bi>=0){ r.pts.splice(bi,1); r.cls&&r.cls.splice(bi,1); } }
    else { r.pts.push([Math.floor(ix)+0.5,Math.floor(iy)+0.5]); (r.cls=r.cls||[]).push(S.ui.ptCls||0); }
    setActive(im,r); return; }
  if(tool==="wand"){ wandAt(im,Math.floor(ix),Math.floor(iy)); return; }
  const type=tool==="scale"?"line":tool; const p0=type==="line"?[ix,iy]:[snapI(ix),snapI(iy)];
  DRAG={mode:"new",type,temp:{id:uid(),type,pts:type==="free"?[p0]:[p0,p0.slice()]},x0:p0[0],y0:p0[1],shift:e.shiftKey,scale:tool==="scale"};
}
function onMove(e){
  const im=cur(); if(!im) return; const [sx,sy]=evPos(e); const v=viewOf(im); const [ix,iy]=toImg(v,sx,sy); HOVER=[ix,iy];
  if(!DRAG){ drawStatus(); updateCursor(im,sx,sy,ix,iy); return; }
  if(DRAG.mode==="pan"){ v.x=DRAG.vx-(sx-DRAG.sx)/v.zoom; v.y=DRAG.vy-(sy-DRAG.sy)/v.zoom; draw(); return; }
  if(DRAG.mode==="poly"||DRAG.mode==="angle"){ const t=DRAG.temp; t.pts[t.pts.length-1]=DRAG.mode==="poly"?[snapI(ix),snapI(iy)]:[ix,iy]; draw(); return; }
  const a=activeRoi(im);
  if(DRAG.mode==="move"&&a){ let dx=ix-DRAG.ix, dy=iy-DRAG.iy; if(a.type!=="line"&&a.type!=="angle"&&a.type!=="point"){ dx=Math.round(dx); dy=Math.round(dy); } a.pts=DRAG.orig.map(p=>[p[0]+dx,p[1]+dy]); draw(); return; }
  if(DRAG.mode==="handle"&&a){ const h=DRAG.h; if(a.type==="rect"||a.type==="oval"){ const b0=(()=>{ const xs=DRAG.orig.map(p=>p[0]), ys=DRAG.orig.map(p=>p[1]); return [Math.min(...xs),Math.min(...ys),Math.max(...xs),Math.max(...ys)]; })(); let [x0,y0,x1,y1]=b0; const X=snapI(ix), Y=snapI(iy);
      if([0,3,7].includes(h)) x0=X; if([1,2,5].includes(h)) x1=X; if([0,1,4].includes(h)) y0=Y; if([2,3,6].includes(h)) y1=Y; a.pts=[[Math.min(x0,x1),Math.min(y0,y1)],[Math.max(x0,x1),Math.max(y0,y1)]]; }
    else { a.pts[h]=a.type==="line"||a.type==="angle"?[ix,iy]:[snapI(ix),snapI(iy)]; } draw(); return; }
  if(DRAG.mode==="new"){ const t=DRAG.temp;
    if(t.type==="free"){ const q=[snapI(ix),snapI(iy)]; const l=t.pts[t.pts.length-1]; if(q[0]!==l[0]||q[1]!==l[1]) t.pts.push(q); }
    else if(t.type==="line"){ let x=ix, y=iy; if(e.shiftKey){ const dx=x-DRAG.x0, dy=y-DRAG.y0; if(Math.abs(dx)>Math.abs(dy)) y=DRAG.y0; else x=DRAG.x0; } t.pts[1]=[x,y]; }
    else { let X=snapI(ix), Y=snapI(iy); if(e.shiftKey){ const d=Math.max(Math.abs(X-DRAG.x0),Math.abs(Y-DRAG.y0)); X=DRAG.x0+Math.sign(X-DRAG.x0||1)*d; Y=DRAG.y0+Math.sign(Y-DRAG.y0||1)*d; } t.pts=[[Math.min(DRAG.x0,X),Math.min(DRAG.y0,Y)],[Math.max(DRAG.x0,X),Math.max(DRAG.y0,Y)]]; }
    draw(); }
}
function onUp(e){
  const im=cur(); if(!im||!DRAG) return; const m=DRAG.mode;
  if(m==="poly"||m==="angle") return;
  if(m==="new"){ const t=DRAG.temp; DRAG=null;
    if(t.type==="free"){ if(t.pts.length>=3) setActive(im,t); else setActive(im,null); }
    else if(t.type==="line"){ const L=Math.hypot(t.pts[1][0]-t.pts[0][0],t.pts[1][1]-t.pts[0][1]); if(L<0.5){ setActive(im,null); } else { setActive(im,t); if(S.ui.tool==="scale"){ openSetScale(im,L); } } }
    else { const b=roiBounds(t); if(b.w<1||b.h<1) setActive(im,null); else setActive(im,t); }
    return; }
  DRAG=null; save(); draw(); drawSide(); void e;
}
function finishPoly(im){ const t=DRAG.temp; DRAG=null; t.pts.pop(); const L=t.pts.filter((p,k)=>k===0||p[0]!==t.pts[k-1][0]||p[1]!==t.pts[k-1][1]); t.pts=L; if(L.length>=3) setActive(im,t); else setActive(im,null); }
function onDbl(){ const im=cur(); if(im&&DRAG&&DRAG.mode==="poly"){ DRAG.temp.pts.pop(); finishPoly(im); } }
function updateCursor(im,sx,sy,ix,iy){ const el=$q(ROOT,"#imOverlay"); const t=S.ui.tool; let c="crosshair"; if(SPACE||t==="hand") c="grab"; else if(t==="zoom") c="zoom-in";
  else { const a=activeRoi(im); if(a){ if(hitHandle(im,a,sx,sy)>=0) c="pointer"; else if(insideRoi(im,a,ix,iy)) c="move"; } } el.style.cursor=c; }
function onWheel(e){ const im=cur(); if(!im) return; e.preventDefault(); const [sx,sy]=evPos(e); if(e.ctrlKey||e.metaKey||!e.shiftKey){ zoomAt(im,Math.exp(-e.deltaY*0.0018),sx,sy); } else { const v=viewOf(im); v.x+=e.deltaY/v.zoom; draw(); } }
/* ---------- 터치 (window.INNO_TOUCH) ----------
   두 손가락 = 확대(핀치 중심 고정)·이동 — 어느 도구든. 한 손가락 = 그리기 도구(rect·oval·poly·free·line·angle·scale)면 그대로 그리기,
   Hand 면 이동, Point·Wand·Zoom 은 탭 한 번이 클릭 한 번(끌면 이동). 탭 = 그 화소 값(상태 줄), 길게 누르기 = 값 읽기(따라 움직임),
   두 번 탭 = 화면에 맞춤. Point 도구에서 길게 누르기 = 가까운 점 지우기(데스크톱의 Alt+클릭). */
const pickOn=()=>{ const im=cur(); return !!(S.ui.paPick&&im&&S.particles[im.id]); };
const TAP_TOOLS=new Set(["point","wand","zoom"]), PASS_DT=new Set(["rect","oval","line","free","scale"]);
let TSNAP=null;
function touchSnap(im){ const a=activeRoi(im); TSNAP={id:im.id,a,pts:a?JSON.parse(JSON.stringify(a.pts)):null,cls:a&&a.cls?a.cls.slice():null,drag:DRAG?JSON.parse(JSON.stringify(DRAG)):null}; }
function touchRestore(){ const t=TSNAP; TSNAP=null; const im=cur(); if(!im||!t||t.id!==im.id){ DRAG=null; draw(); return; }
  if(t.a){ t.a.pts=t.pts; if(t.cls) t.a.cls=t.cls; S.roiSel[im.id]=t.a; } else delete S.roiSel[im.id];
  DRAG=t.drag; save(); draw(); drawSide(); }
function bindTouch(ov){
  if(!window.INNO_TOUCH||ov._touch) return;
  const at=p=>{ const im=cur(); if(!im) return null; const v=viewOf(im); return {im,ix:v.x+p.x/v.zoom,iy:v.y+p.y/v.zoom}; };
  const fake=(p,alt)=>({clientX:p.clientX,clientY:p.clientY,button:0,pointerId:-1,pointerType:"touch",shiftKey:false,altKey:!!alt,ctrlKey:false,metaKey:false,_tap:true});
  const showVal=p=>{ const q=at(p); if(!q) return; HOVER=[q.ix,q.iy]; drawStatus(); };
  ov._touch=INNO_TOUCH.attach(ov,{
    oneFinger:()=>{ const t=S.ui.tool; return !cur()||t==="hand"||TAP_TOOLS.has(t)||pickOn()?"pan":"pass"; },
    getView:()=>{ const im=cur(); return im?Object.assign({id:im.id},viewOf(im)):null; },
    setView:(v0,g)=>{ const im=cur(); if(!im||!v0||v0.id!==im.id) return; const n=INNO_TOUCH.affine(v0,g,{min:1/72,max:64}); const v=viewOf(im); v.zoom=n.zoom; v.x=n.x; v.y=n.y; draw(); },
    onPassCancel:touchRestore,
    onTap:p=>{ LASTPT="touch"; if(TAP_TOOLS.has(S.ui.tool)||pickOn()) onDown(fake(p)); showVal(p); },
    onDoubleTap:(p,info)=>{ const im=cur(); if(!im) return; if(info&&info.pass&&!PASS_DT.has(S.ui.tool)) return; fitView(im); draw(); },
    doubleTap:()=>!TAP_TOOLS.has(S.ui.tool)&&!pickOn(),
    doubleTapInPass:true,
    onLongPress:p=>{ LASTPT="touch"; const im=cur(); const a=im&&activeRoi(im);
      if(S.ui.tool==="point"&&a&&a.type==="point"){ const v=viewOf(im); if(a.pts.some(q=>{ const [x,y]=toScr(v,q[0],q[1]); return Math.hypot(x-p.x,y-p.y)<16; })){ onDown(fake(p,true)); return; } }
      showVal(p); },
    onLongPressMove:showVal
  });
}
function drawPolyBar(){ const el=ROOT&&$q(ROOT,"#imPolyBar"); if(!el) return;
  const on=!!(DRAG&&(DRAG.mode==="poly"||DRAG.mode==="angle")&&LASTPT==="touch"); if(!on){ if(!el.hidden){ el.hidden=true; el.innerHTML=""; } return; }
  const n=DRAG.temp.pts.length-1, key=DRAG.mode+n; if(!el.hidden&&el.dataset.k===key) return; el.dataset.k=key; el.hidden=false;
  el.innerHTML=DRAG.mode==="poly"?`<span class="t">${n} pts</span><button class="btn sm primary" id="imPolyDone" ${n<3?"disabled":""}>Finish</button><button class="btn sm" id="imPolyEsc">Cancel</button>`:`<span class="t">Angle ${n}/3</span><button class="btn sm" id="imPolyEsc">Cancel</button>`;
  const d=$q(el,"#imPolyDone"); if(d) d.onclick=()=>{ const im=cur(); if(im&&DRAG&&DRAG.mode==="poly") finishPoly(im); };
  $q(el,"#imPolyEsc").onclick=()=>{ DRAG=null; draw(); }; }
/* Wand: 문턱 범위(Threshold 탭)·이진 이미지면 그 범위로, 아니면 클릭한 값과 같은 화소(허용 tolerance)로 8-연결 외곽선 */
let WAND_BUSY=false;
async function wandAt(im,x,y){ if(x<0||y<0||x>=im.w||y>=im.h) return; if(WAND_BUSY){ note("Wand: 앞의 외곽선을 아직 찾는 중입니다"); return; }
  const p=curPlaneOf(im); let lo,hi;
  if(S.ui.panel==="threshold"&&!im.rgb){ ({lo,hi}=thrRange(im)); } else { const v=p.data[y*im.w+x]; const tol=+(S.ui.wandTol||0); lo=v-tol; hi=v+tol; }
  const id=im.id, big=im.w*im.h>4e6; WAND_BUSY=true; let ol;
  try{ if(big) setProg("Wand …",0.1); ol=await call("wand",[{w:p.w,h:p.h,bit:p.bit,data:p.data},x,y,lo,hi]); }
  catch(e){ if(!/취소/.test(e.message)) note("Wand: "+e.message,"warn"); return; }
  finally{ WAND_BUSY=false; if(big) hideProg(); }
  if(!ol){ note("클릭한 화소가 범위 밖입니다"); return; }
  const im2=S.images.find(z=>z.id===id); if(!im2) return; // 그 사이 닫힌 이미지
  setActive(im2,{id:uid(),type:"traced",pts:ol.xs.map((xx,k)=>[xx,ol.ys[k]])});
}

/* ---------- 상태 줄 ---------- */
function drawStatus(){
  const el=ROOT&&$q(ROOT,"#imStatus"); if(!el) return; const im=cur(); if(!im){ el.innerHTML=""; return; } const v=viewOf(im); const c=calOf(im);
  let s=`<span>${Math.round(v.zoom*100)}%</span>`;
  if(HOVER){ const [ix,iy]=HOVER; const x=Math.floor(ix), y=Math.floor(iy);
    if(x>=0&&y>=0&&x<im.w&&y<im.h){ const vals=[]; const chans=im.rgb?[0,1,2]:[...Array(im.nc).keys()]; chans.forEach(ch=>{ const d=planesOf(im)[pIndex(im,ch,im.z,im.t)]; vals.push(g4(d[y*im.w+x])); });
      s+=`<span>x ${x}, y ${y}${im.cal?` <b>(${g4((x+0.5)*c.pw)}, ${g4((y+0.5)*c.ph)} ${c.unit})</b>`:""}</span><span>value ${im.rgb?"RGB ":""}${vals.join(", ")}${im.valueUnit?" "+esc(im.valueUnit):""}</span>`; } }
  const a=activeRoi(im);
  if(a){ if(a.type==="line") s+=`<span class="roi">length ${g4(lineLen(im,a))} ${c.unit}, angle ${roiAngle(a).toFixed(1)}°</span>`;
    else if(a.type==="angle"&&a.pts.length===3) s+=`<span class="roi">angle ${roiAngle(a).toFixed(2)}°</span>`;
    else if(a.type==="point") s+=`<span class="roi">${a.pts.length} points</span>`;
    else { const b=roiBounds(a); s+=`<span class="roi">${a.type} ${Math.round(b.w)}×${Math.round(b.h)} px${im.cal?` (${g4(b.w*c.pw)}×${g4(b.h*c.ph)} ${c.unit})`:""}</span>`; } }
  if(PREVIEW.plane) s+=`<span class="pv">Preview — not applied</span>`;
  el.innerHTML=s;
}
/* 선 ROI 의 프로파일 (캔버스 아래 작은 그래프) */
function profileOf(im,r){ const p=curPlaneOf(im); const [a,b]=r.pts; const pr=IC.lineProfile(p,a[0],a[1],b[0],b[1],S.ui.lw||1); const c=calOf(im); return {x:[...pr.values].map((_,i)=>i*pr.step*c.pw),y:[...pr.values],unit:c.unit}; }
function drawProfileMini(){
  const el=ROOT&&$q(ROOT,"#imProfile"); if(!el) return; const im=cur(); const a=im&&activeRoi(im);
  if(!a||a.type!=="line"){ el.hidden=true; return; } el.hidden=false; const pr=profileOf(im,a);
  el.innerHTML=`<div class="ph"><b>Line profile</b><span class="muted xs">width ${S.ui.lw||1} px · ${pr.y.length} points</span><span class="sp"></span><button class="btn xs" id="imProfAdd">Add to Plots</button><button class="btn xs" id="imProfCsv">CSV</button></div>${svgXY([{name:"profile",x:pr.x,y:pr.y}],{w:640,h:150,xl:`distance (${pr.unit})`,yl:"gray value"})}`;
  $q(el,"#imProfAdd").onclick=()=>{ S.plots.push({id:uid(),kind:"xy",title:"Profile — "+im.name,xl:`distance (${pr.unit})`,yl:"gray value",series:[{name:"profile",x:pr.x,y:pr.y}]}); S.ui.bottom="plots"; S.ui.bottomOpen=true; save(); drawBottom(); };
  $q(el,"#imProfCsv").onclick=()=>downloadText(`profile_${im.name}.csv`,"distance_"+pr.unit+",value\n"+pr.x.map((x,i)=>x+","+pr.y[i]).join("\n"));
}
