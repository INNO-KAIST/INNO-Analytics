/* ============================================================
   Process ▸ Binary  —  Erode/Dilate/Open/Close (count·iterations), Outline, Fill Holes, Skeletonize,
   Distance Map (EDM), Ultimate Points, Watershed
   원본: ij/process/ByteProcessor.java (filter ERODE/DILATE), ij/process/BinaryProcessor.java,
         ij/plugin/filter/Binary.java, ij/plugin/filter/EDM.java
   규칙: Black background (물체 = 255, 배경 = 0). 이미지 밖 = 배경 (Prefs.padEdges=false)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const FG=255, BG=0;

function isBinary(p){ if(p.bit!==8) return false; const d=p.data; for(let i=0;i<d.length;i++){ const v=d[i]; if(v!==0&&v!==255) return false; } return true; }
IC.isBinary=isBinary;
function needBinary(p){ if(!isBinary(p)) throw new Error("이진(0/255) 8-bit 이미지가 필요합니다 — 먼저 Threshold ▸ Apply (Convert to Mask)"); }

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

/* ===== Split by size (크기 기준 나누기) — 2026-10-19 =====
   거의 같은 크기의 둥근 입자(나노입자 TEM)가 겹쳐 한 덩어리가 된 것을 나눈다. ImageJ Watershed(EDM 극대, 허용 0.5)는 목이 넓게 겹친 덩어리를 못 나눔.
   ① 성분마다 넓이·solidity(넓이/볼록껍질) → 홑입자(solidity ≥ 0.93, 9 px 이상)의 넓이 중앙값 A0, r0 = √(A0/π)
   ② 덩어리(넓이 > 1.35·A0 또는 solidity < 0.88): 씨앗 = EDM 의 국소 극대 중 높이 ≥ 0.45·r0, 서로 ≥ 0.75·r0 떨어진 것(높은 것부터; 허용 0 — 낮은 목도 나눔)
      씨앗 수 < round(넓이/A0) 이면 가장 먼 화소를 씨앗으로 더함
   ③ 씨앗에서 k-means(화소 좌표, 10 번) → 가장 가까운 씨앗으로 나눔(같은 크기 원이 겹친 경계 = 두 중심의 수직이등분선) → 경계에 1 px 선(0)
   반환: 이진 이미지(255/0), info {A0, r0, split: 나눈 덩어리 수, parts} */
IC.splitBySize=function(p,o){ o=o||{}; needBinary(p); const W=p.w, H=p.h, N=W*H, src=p.data;
  const lab=new Int32Array(N).fill(-1), st=new Int32Array(N); const comps=[];
  for(let i0=0;i0<N;i0++){ if(!(src[i0]>0)||lab[i0]>=0) continue; const id=comps.length; let sp=0; st[sp++]=i0; lab[i0]=id; const px=[];
    while(sp){ const i=st[--sp]; px.push(i); const x=i%W, y=(i-x)/W; for(let dy=-1;dy<=1;dy++){ const yy=y+dy; if(yy<0||yy>=H) continue; for(let dx=-1;dx<=1;dx++){ const xx=x+dx; if(xx<0||xx>=W) continue; const j=yy*W+xx; if(src[j]>0&&lab[j]<0){ lab[j]=id; st[sp++]=j; } } } }
    comps.push({id,px}); }
  // solidity: 볼록 껍질(화소 모서리 점들의 monotone chain) 넓이 대비
  const hullArea=px=>{ const P=[]; const rows=new Map(); for(const i of px){ const x=i%W, y=(i-x)/W; const r=rows.get(y); if(!r) rows.set(y,[x,x]); else { if(x<r[0]) r[0]=x; if(x>r[1]) r[1]=x; } }
    for(const [y,[a,b]] of rows){ P.push([a,y],[b+1,y],[a,y+1],[b+1,y+1]); } P.sort((u,v)=>u[0]-v[0]||u[1]-v[1]);
    const cr=(O,A,B)=>(A[0]-O[0])*(B[1]-O[1])-(A[1]-O[1])*(B[0]-O[0]); const lo=[],up=[];
    for(const q of P){ while(lo.length>=2&&cr(lo[lo.length-2],lo[lo.length-1],q)<=0) lo.pop(); lo.push(q); }
    for(let k=P.length-1;k>=0;k--){ const q=P[k]; while(up.length>=2&&cr(up[up.length-2],up[up.length-1],q)<=0) up.pop(); up.push(q); }
    const hu=lo.slice(0,-1).concat(up.slice(0,-1)); let a=0; for(let k=0;k<hu.length;k++){ const A=hu[k], B=hu[(k+1)%hu.length]; a+=A[0]*B[1]-B[0]*A[1]; } return Math.abs(a)/2; };
  for(const c of comps){ c.area=c.px.length; c.sol=c.area>=4?c.area/Math.max(c.area,hullArea(c.px)):1; }
  // 잡음 점(수천 개의 작은 성분)이 중앙값을 차지하지 않게: 넓이 가중 중앙값의 minFrac(기본 0.15) 보다 작은 성분은 홑입자 후보에서 뺌
  const allA=comps.filter(c=>c.area>=9).map(c=>c.area).sort((a,b)=>a-b); let wm=0; { let tot=0; for(const a of allA) tot+=a; let cs=0; for(const a of allA){ cs+=a; if(cs>=tot/2){ wm=a; break; } } }
  const aMin=Math.max(9,(o.minFrac>=0?o.minFrac:0.15)*wm);
  let singles=comps.filter(c=>c.area>=aMin&&c.sol>=0.93).map(c=>c.area).sort((a,b)=>a-b);
  if(singles.length<3) singles=comps.filter(c=>c.area>=aMin).map(c=>c.area).sort((a,b)=>a-b);
  const out=new Float32Array(N); for(let i=0;i<N;i++) out[i]=src[i]>0?255:0;
  if(!singles.length) return {w:W,h:H,bit:8,data:out,info:{A0:NaN,r0:NaN,split:0,parts:0}};
  const A0=o.A0>0?o.A0:singles[singles.length>>1], r0=Math.sqrt(A0/Math.PI);
  const e=IC.makeFloatEDM(p,0,false).data; let nSplit=0, nParts=0;
  // 세기 씨앗(선택, o.gray): 입자 크기 척도의 Laplacian-of-Gaussian 극대 = 원판 중심 (겹쳐 붙어 EDM 에 목이 없어도 찾음).
  // 큰 입자는 f 배 줄인 이미지에서 계산(σ' = σ/f) — 속도. 어두운 입자면 +∇²(G*I), 밝은 입자면 −∇²
  let lseed=null; if(o.gray&&o.gray.data&&o.gray.data.length===N){ const f=Math.max(1,Math.floor(r0/8)), w2=Math.ceil(W/f), h2=Math.ceil(H/f); const sm=new Float32Array(w2*h2), cnt=new Float32Array(w2*h2);
      for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const k=((y/f)|0)*w2+((x/f)|0); sm[k]+=o.gray.data[y*W+x]; cnt[k]++; } for(let k=0;k<sm.length;k++) sm[k]/=cnt[k]||1;
      const sg=Math.max(1,r0/Math.SQRT2/f); const gg=IC.gaussian({w:w2,h:h2,bit:32,data:sm},{sigma:sg}); const L=(gg.plane||gg).data; const sgn=o.dark?1:-1;
      const lap=new Float32Array(w2*h2); for(let y=1;y<h2-1;y++) for(let x=1;x<w2-1;x++){ const k=y*w2+x; lap[k]=sgn*(L[k-1]+L[k+1]+L[k-w2]+L[k+w2]-4*L[k]); }
      lseed=[]; for(let y=1;y<h2-1;y++) for(let x=1;x<w2-1;x++){ const k=y*w2+x, v=lap[k]; if(!(v>0)) continue; let mx=true; for(let dy=-1;dy<=1&&mx;dy++) for(let dx=-1;dx<=1;dx++){ if((dx||dy)&&lap[k+dy*w2+dx]>v){ mx=false; break; } } if(mx){ const X=Math.min(W-1,Math.round((x+0.5)*f-0.5)), Y=Math.min(H-1,Math.round((y+0.5)*f-0.5)); lseed.push([X,Y,v]); } } }
  const seedsOf=new Map(); if(lseed){ for(const q of lseed){ const id=lab[q[1]*W+q[0]]; if(id<0) continue; (seedsOf.get(id)||seedsOf.set(id,[]).get(id)).push(q); } }
  for(const c of comps){ if(!(c.area>1.35*A0||(c.sol<0.88&&c.area>0.6*A0))) continue;
    // 씨앗: EDM 국소 극대(3×3) 중 높이 ≥ 0.45 r0 → 높은 것부터, 서로 0.75 r0 이상
    const cand=[]; for(const i of c.px){ const v=e[i]; if(v<0.45*r0) continue; const x=i%W, y=(i-x)/W; let mx=true;
      for(let dy=-1;dy<=1&&mx;dy++) for(let dx=-1;dx<=1;dx++){ if(!dx&&!dy) continue; const xx=x+dx, yy=y+dy; if(xx<0||yy<0||xx>=W||yy>=H) continue; if(e[yy*W+xx]>v){ mx=false; break; } }
      if(mx) cand.push([x,y,v]); }
    cand.sort((a,b)=>b[2]-a[2]); const seeds=[]; const dmin=0.75*r0;
    const ls=seedsOf.get(c.id); if(ls&&ls.length){ ls.sort((a,b)=>b[2]-a[2]); const top=ls[0][2]; const dl=0.9*r0; for(const q of ls){ if(q[2]<0.25*top) break; if(seeds.every(s=>(s[0]-q[0])**2+(s[1]-q[1])**2>=dl*dl)) seeds.push([q[0],q[1]]); } } // 세기 씨앗이 있으면 그것으로
    else for(const q of cand){ if(seeds.every(s=>(s[0]-q[0])**2+(s[1]-q[1])**2>=dmin*dmin)) seeds.push([q[0],q[1]]); }
    // 넓이로 씨앗을 더함: 오목하면 넓이 ≥ (k − 0.4)·A0, 볼록하면 ≥ (k − 0.25)·A0 일 때 k 개 — 볼록한 큰 홑입자(넓이 < 1.75·A0)는 쪼개지 않음
    const kA=Math.max(seeds.length,Math.floor(c.area/A0+(c.sol<0.9?0.4:0.25)));
    while(seeds.length<kA){ let best=-1,bi=-1; for(const i of c.px){ const x=i%W, y=(i-x)/W; let d=Infinity; for(const s of seeds){ const t=(s[0]-x)**2+(s[1]-y)**2; if(t<d) d=t; } if(seeds.length===0) d=e[i]; if(d>best){ best=d; bi=i; } } if(bi<0) break; seeds.push([bi%W,(bi-bi%W)/W]); }
    if(seeds.length<2) continue;
    // k-means (화소 좌표)
    const asg=new Int32Array(c.px.length); let S=seeds.map(s=>[s[0]+0.5,s[1]+0.5]);
    for(let it=0;it<10;it++){ const sx=new Float64Array(S.length), sy=new Float64Array(S.length), sn=new Float64Array(S.length);
      for(let k=0;k<c.px.length;k++){ const i=c.px[k], x=i%W+0.5, y=(i-i%W)/W+0.5; let bd=Infinity, bj=0; for(let j=0;j<S.length;j++){ const d=(S[j][0]-x)**2+(S[j][1]-y)**2; if(d<bd){ bd=d; bj=j; } } asg[k]=bj; sx[bj]+=x; sy[bj]+=y; sn[bj]++; }
      let mv=0; for(let j=0;j<S.length;j++) if(sn[j]){ const nx=sx[j]/sn[j], ny=sy[j]/sn[j]; mv=Math.max(mv,Math.abs(nx-S[j][0])+Math.abs(ny-S[j][1])); S[j]=[nx,ny]; } if(mv<0.05) break; }
    const part=new Map(); for(let k=0;k<c.px.length;k++) part.set(c.px[k],asg[k]);
    // 경계 1 px: 오른쪽·아래 이웃이 다른 조각이면 이 화소를 0 (8-연결로도 붙지 않게 대각 이웃까지)
    for(const [i,a] of part){ const x=i%W, y=(i-x)/W; for(const [dx,dy] of [[1,0],[0,1],[1,1],[-1,1]]){ const xx=x+dx, yy=y+dy; if(xx<0||yy<0||xx>=W||yy>=H) continue; const b=part.get(yy*W+xx); if(b!=null&&b!==a){ out[i]=0; break; } } }
    nSplit++; nParts+=S.length; }
  return {w:W,h:H,bit:8,data:out,info:{A0,r0,split:nSplit,parts:nParts}};
};
})();
