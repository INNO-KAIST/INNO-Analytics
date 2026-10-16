/* ============================================================
   Process ▸ Binary  —  Erode/Dilate/Open/Close (count·iterations), Outline, Fill Holes, Skeletonize,
   Distance Map (EDM), Ultimate Points, Watershed
   원본: ij/process/ByteProcessor.java (filter ERODE/DILATE), ij/process/BinaryProcessor.java,
         ij/plugin/filter/Binary.java, ij/plugin/filter/EDM.java
   규칙: Black background (물체 = 255, 배경 = 0). 영상 밖 = 배경 (Prefs.padEdges=false)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const FG=255, BG=0;

function isBinary(p){ if(p.bit!==8) return false; const d=p.data; for(let i=0;i<d.length;i++){ const v=d[i]; if(v!==0&&v!==255) return false; } return true; }
IC.isBinary=isBinary;
function needBinary(p){ if(!isBinary(p)) throw new Error("이진(0/255) 8-bit 영상이 필요합니다 — 먼저 Threshold ▸ Apply (Convert to Mask)"); }

function erodeDilate(src,W,H,count,erode){
  const out=new Float32Array(W*H);
  const g=(x,y)=>(x<0||x>=W||y<0||y>=H)?BG:src[y*W+x];
  const target=erode?BG:FG, other=erode?FG:BG;
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){
    const p5=src[y*W+x];
    if(p5===target){ out[y*W+x]=target; continue; }
    let c=0;
    if(g(x-1,y-1)===target) c++; if(g(x,y-1)===target) c++; if(g(x+1,y-1)===target) c++;
    if(g(x-1,y)===target) c++; if(g(x+1,y)===target) c++;
    if(g(x-1,y+1)===target) c++; if(g(x,y+1)===target) c++; if(g(x+1,y+1)===target) c++;
    out[y*W+x]=c>=count?target:other;
  }
  return out;
}
/* op: erode | dilate | open | close ; o.iterations (1), o.count (1..8) */
IC.binaryMorph=function(p,op,o){
  needBinary(p); o=o||{}; const it=Math.max(1,o.iterations|0||1), cnt=Math.min(8,Math.max(1,o.count|0||1));
  let d=Float32Array.from(p.data); const W=p.w,H=p.h;
  const run=(er)=>{ for(let i=0;i<it;i++) d=erodeDilate(d,W,H,cnt,er); };
  if(op==="erode") run(true); else if(op==="dilate") run(false);
  else if(op==="open"){ run(true); run(false); } else if(op==="close"){ run(false); run(true); }
  else throw new Error("binary 연산: "+op);
  return {w:W,h:H,bit:8,data:d};
};

/* Outline (BinaryProcessor.process OUTLINE): 안쪽 화소(8-이웃이 모두 물체)를 지운다. 가장자리 행·열은 손대지 않음 */
IC.outline=function(p){
  needBinary(p); const W=p.w,H=p.h,s=p.data,out=Float32Array.from(s);
  for(let y=1;y<H-1;y++) for(let x=1;x<W-1;x++){ const i=y*W+x; if(s[i]!==FG) continue;
    if(s[i-W-1]!==BG&&s[i-W]!==BG&&s[i-W+1]!==BG&&s[i-1]!==BG&&s[i+1]!==BG&&s[i+W-1]!==BG&&s[i+W]!==BG&&s[i+W+1]!==BG) out[i]=BG; }
  return {w:W,h:H,bit:8,data:out};
};

/* Fill Holes (Binary.fill, Landini): 가장자리에서 4-연결로 닿는 배경만 남기고 나머지는 물체로 */
IC.fillHoles=function(p){
  needBinary(p); const W=p.w,H=p.h,s=p.data, mark=new Uint8Array(W*H), st=[];
  const push=(x,y)=>{ const i=y*W+x; if(s[i]===BG&&!mark[i]){ mark[i]=1; st.push(i); } };
  for(let y=0;y<H;y++){ push(0,y); push(W-1,y); } for(let x=0;x<W;x++){ push(x,0); push(x,H-1); }
  while(st.length){ const i=st.pop(), x=i%W, y=(i-x)/W; if(x>0) push(x-1,y); if(x<W-1) push(x+1,y); if(y>0) push(x,y-1); if(y<H-1) push(x,y+1); }
  const out=new Float32Array(W*H); for(let i=0;i<W*H;i++) out[i]=mark[i]?BG:FG;
  return {w:W,h:H,bit:8,data:out};
};

/* Skeletonize (BinaryProcessor.skeletonize, 표 기반 thinning) */
const T1=[0,0,0,0,0,0,1,3,0,0,3,1,1,0,1,3,0,0,0,0,0,0,0,0,0,0,2,0,3,0,3,3,
0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,3,0,2,2,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,0,0,0,0,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,0,0,0,0,3,0,0,0,3,0,2,0,
0,0,3,1,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,3,1,3,0,0,1,3,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2,3,0,1,0,0,0,1,0,0,0,0,0,0,0,0,3,3,0,1,0,0,0,0,2,2,0,0,2,0,0,0];
const T2=[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
function thin(px,W,H,pass,table){
  const c=Uint8Array.from(px); let removed=0;
  for(let y=1;y<=H-2;y++){ let o=1+y*W; for(let x=1;x<=W-2;x++,o++){
    if(c[o]===0) continue;
    let idx=0;
    if(c[o-W-1]) idx|=1; if(c[o-W]) idx|=2; if(c[o-W+1]) idx|=4; if(c[o+1]) idx|=8;
    if(c[o+W+1]) idx|=16; if(c[o+W]) idx|=32; if(c[o+W-1]) idx|=64; if(c[o-1]) idx|=128;
    const code=table[idx];
    if(pass&1){ if(code===2||code===3){ px[o]=0; removed++; } } else { if(code===1||code===3){ px[o]=0; removed++; } }
  } }
  return removed;
}
IC.skeletonize=function(p){
  needBinary(p); const W0=p.w,H0=p.h; let edge=false;
  for(let x=0;x<W0&&!edge;x++) if(p.data[x]===FG||p.data[(H0-1)*W0+x]===FG) edge=true;
  for(let y=0;y<H0&&!edge;y++) if(p.data[y*W0]===FG||p.data[y*W0+W0-1]===FG) edge=true;
  const W=edge?W0+2:W0, H=edge?H0+2:H0, px=new Uint8Array(W*H), o=edge?1:0;
  for(let y=0;y<H0;y++) for(let x=0;x<W0;x++) px[(y+o)*W+x+o]=p.data[y*W0+x]===FG?1:0;
  // 가장자리 줄을 배경으로 (skeletonize2 의 moveTo/lineTo)
  for(let x=0;x<W;x++){ px[x]=0; px[(H-1)*W+x]=0; } for(let y=0;y<H;y++){ px[y*W]=0; px[y*W+W-1]=0; }
  let pass=0, rem;
  do{ rem=thin(px,W,H,pass++,T1); rem+=thin(px,W,H,pass++,T1); }while(rem>0);
  do{ rem=thin(px,W,H,pass++,T2); rem+=thin(px,W,H,pass++,T2); }while(rem>0);
  const out=new Float32Array(W0*H0); for(let y=0;y<H0;y++) for(let x=0;x<W0;x++) out[y*W0+x]=px[(y+o)*W+x+o]?FG:BG;
  return {w:W0,h:H0,bit:8,data:out};
};

/* ---------- EDM (EDM.makeFloatEDM) — 두 방향 두 번 훑는 ImageJ 알고리즘 ----------
   물체(≠backgroundValue) 화소의 가장 가까운 배경 화소까지 거리. edgesAreBackground=false 가 ImageJ 기본 */
IC.makeFloatEDM=function(p,backgroundValue,edgesAreBackground){
  const W=p.w,H=p.h,b=p.data, f=new Float32Array(W*H); const bgv=backgroundValue||0;
  for(let i=0;i<W*H;i++) if(b[i]!==bgv) f[i]=3.4028234663852886e38;
  const NO=-1, buf0=new Int32Array(W), buf1=new Int32Array(W);
  const minDist2=(points,pPrev,pDiag,x,y,distSqr)=>{
    const p0=points[x]; let nearest=p0;
    if(p0!==NO){ const x0=p0&0xffff, y0=(p0>>16)&0xffff; const d=(x-x0)*(x-x0)+(y-y0)*(y-y0); if(d<distSqr) distSqr=d; }
    if(pDiag!==p0&&pDiag!==NO){ const x1=pDiag&0xffff, y1=(pDiag>>16)&0xffff; const d=(x-x1)*(x-x1)+(y-y1)*(y-y1); if(d<distSqr){ nearest=pDiag; distSqr=d; } }
    if(pPrev!==pDiag&&pPrev!==NO){ const x1=pPrev&0xffff, y1=(pPrev>>16)&0xffff; const d=(x-x1)*(x-x1)+(y-y1)*(y-y1); if(d<distSqr){ nearest=pPrev; distSqr=d; } }
    points[x]=nearest; return distSqr;
  };
  const edmLine=(offset,y,yDist)=>{
    const eab=yDist!==2147483647; let distSqr=2147483647;
    let points=buf0, pPrev=NO, pDiag=NO, pNext;
    for(let x=0;x<W;x++,offset++){ pNext=points[x];
      if(b[offset]===bgv) points[x]=x|(y<<16);
      else { if(eab) distSqr=(x+1<yDist)?(x+1)*(x+1):yDist*yDist; const d2=minDist2(points,pPrev,pDiag,x,y,distSqr); if(f[offset]>d2) f[offset]=d2; }
      pPrev=points[x]; pDiag=pNext; }
    offset--; points=buf1; pPrev=NO; pDiag=NO;
    for(let x=W-1;x>=0;x--,offset--){ pNext=points[x];
      if(b[offset]===bgv) points[x]=x|(y<<16);
      else { if(eab) distSqr=(W-x<yDist)?(W-x)*(W-x):yDist*yDist; const d2=minDist2(points,pPrev,pDiag,x,y,distSqr); if(f[offset]>d2) f[offset]=d2; }
      pPrev=points[x]; pDiag=pNext; }
  };
  buf0.fill(NO); buf1.fill(NO);
  for(let y=0;y<H;y++) edmLine(y*W,y,edgesAreBackground?y+1:2147483647);
  buf0.fill(NO); buf1.fill(NO);
  for(let y=H-1;y>=0;y--) edmLine(y*W,y,edgesAreBackground?H-y:2147483647);
  for(let i=0;i<W*H;i++) f[i]=Math.sqrt(f[i]);
  return {w:W,h:H,bit:32,data:f};
};
/* Distance Map 명령 (EDM output: "8-bit"(Overwrite) 기본 → (int)(v+0.5) 를 255 로 자름) */
IC.distanceMap=function(p,out){
  needBinary(p); const e=IC.makeFloatEDM(p,0,false);
  if(out==="32-bit") return e;
  const d=new Float32Array(e.data.length); for(let i=0;i<d.length;i++){ let v=Math.floor(e.data[i]+0.5); if(v>255) v=255; d[i]=v; }
  return {w:p.w,h:p.h,bit:8,data:d};
};
/* Watershed (EDM.toWatershed): EDM 의 극대(허용 0.5)를 씨앗으로 MaximumFinder SEGMENTED, 원본과 AND */
IC.watershed=function(p){
  needBinary(p); const e=IC.makeFloatEDM(p,0,false);
  const seg=IC.findMaximaImage(e,0.5,{output:"segmented",isEDM:true});
  const out=new Float32Array(p.data.length); for(let i=0;i<out.length;i++) out[i]=(p.data[i]&&seg.data[i])?255:0;
  return {w:p.w,h:p.h,bit:8,data:out};
};
/* Ultimate Points (EDM 'points'): EDM 극대점 하나씩, 값 = EDM 높이(8-bit 로) */
IC.ultimatePoints=function(p){
  needBinary(p); const e=IC.makeFloatEDM(p,0,false);
  const pts=IC.findMaximaImage(e,0.5,{output:"single",isEDM:true});
  const out=new Float32Array(p.data.length); for(let i=0;i<out.length;i++) if(pts.data[i]){ let v=Math.floor(e.data[i]+0.5); if(v>255) v=255; out[i]=v; }
  return {w:p.w,h:p.h,bit:8,data:out};
};
IC.invert=function(p,min,max){ // Edit ▸ Invert: 8-bit 255−v, 16-bit: (min+max)−v (ImageJ 1.54: 표시 범위 아닌 0..65535? → 실제 화소 min/max 사용), 32-bit: (min+max)−v
  const d=p.data, out=new Float32Array(d.length);
  if(p.bit===8){ for(let i=0;i<d.length;i++) out[i]=255-d[i]; return {w:p.w,h:p.h,bit:8,data:out}; }
  let mn=min, mx=max; if(mn==null||mx==null){ [mn,mx]=IC.minMax(d); }
  for(let i=0;i<d.length;i++) out[i]=mn+mx-d[i];
  return {w:p.w,h:p.h,bit:p.bit,data:p.bit===16?IC.roundTo(out,16):out};
};

})();
