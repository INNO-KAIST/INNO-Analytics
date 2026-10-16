/* ============================================================
   Process ▸ Find Maxima  —  ij/plugin/filter/MaximumFinder.java 이식 (watershed 분할 포함)
   prominence(=noise tolerance): 더 높은 극대로 가는 능선보다 이 값보다 더 솟아야 극대로 인정
   output: "points"(목록) | "count" | "single"(Single Points 영상) | "within"(Maxima Within Tolerance) | "segmented"
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const f32=Math.fround;
const MAXIMUM=1, LISTED=2, PROCESSED=4, MAX_AREA=8, EQUAL=16, MAX_POINT=32, ELIMINATED=64;
const DX=[0,1,1,1,0,-1,-1,-1], DY=[-1,-1,0,1,1,1,0,-1];
const SQRT2=f32(1.4142135624);
const NO_THRESHOLD=-808080;

function MF(W,H){
  this.W=W; this.H=H;
  this.dir=[-W,-W+1,1,W+1,W,W-1,-1,-W-1];
  let shift=0,mult=1; do{ shift++; mult*=2; }while(mult<W);
  this.xMask=mult-1; this.shift=shift;
}
MF.prototype.isWithin=function(x,y,d){ const xm=this.W-1, ym=this.H-1;
  switch(d){ case 0: return y>0; case 1: return x<xm&&y>0; case 2: return x<xm; case 3: return x<xm&&y<ym; case 4: return y<ym; case 5: return x>0&&y<ym; case 6: return x>0; case 7: return x>0&&y>0; } return false; };
MF.prototype.trueEdmHeight=function(x,y,px){
  const W=this.W, xm=W-1, ym=this.H-1, off=x+y*W, v=px[off];
  if(x===0||y===0||x===xm||y===ym||v===0) return v;
  let trueH=f32(v+f32(0.5*SQRT2)); let ridge=false;
  for(let d=0;d<4;d++){ const d2=(d+4)%8; const v1=px[off+this.dir[d]], v2=px[off+this.dir[d2]]; let h;
    if(v>=v1&&v>=v2){ ridge=true; h=f32((v1+v2)/2); } else h=Math.min(v1,v2);
    h=f32(h+((d%2===0)?1:SQRT2)); if(trueH>h) trueH=h; }
  if(!ridge) trueH=v; return trueH;
};
MF.prototype.findMaxima=function(px,isFloat,tolerance,strict,threshold,outputType,excludeOnEdges,isEDM,roiMask){
  const W=this.W,H=this.H, N=W*H; const types=new Uint8Array(N);
  let gMin=Infinity,gMax=-Infinity; for(let i=0;i<N;i++){ const v=px[i]; if(v<gMin) gMin=v; if(v>gMax) gMax=v; }
  gMin=f32(gMin); gMax=f32(gMax);
  let maxPossible=gMax>gMin; if(strict&&gMax-gMin<=tolerance) maxPossible=false;
  if(threshold!==NO_THRESHOLD) threshold-=(gMax-gMin)*1e-6;
  const excludeEdgesNow=excludeOnEdges&&outputType!=="segmented";
  const maxPts=maxPossible?this.sortedMax(px,types,excludeEdgesNow,isEDM,gMin,gMax,threshold):[];
  let maxSortErr=0; if(isFloat) maxSortErr=f32(1.1*(isEDM?SQRT2/2:(gMax-gMin)/2e9));
  const xy=this.analyze(px,types,maxPts,excludeEdgesNow,isEDM,gMin,tolerance,strict,outputType,maxSortErr,excludeOnEdges,roiMask);
  if(outputType==="points"||outputType==="count") return {points:xy};
  let out;
  if(outputType==="segmented"){
    out=this.make8bit(px,types,isEDM,gMin,gMax,threshold);
    this.cleanupMaxima(out,types,maxPts);
    this.watershedSegment(out);
    if(!isEDM) this.cleanupExtraLines(out);
    for(let i=0;i<N;i++) if(out[i]<255) out[i]=0;
    if(excludeOnEdges) this.deleteEdgeParticles(out,types);
  }else{
    const m=outputType==="single"?MAX_POINT:MAX_AREA; out=new Uint8Array(N); for(let i=0;i<N;i++) out[i]=(types[i]&m)?255:0;
  }
  if(roiMask) for(let i=0;i<N;i++) if(!roiMask[i]) out[i]=0;
  return {image:out};
};
MF.prototype.sortedMax=function(px,types,exEdges,isEDM,gMin,gMax,threshold){
  const W=this.W,H=this.H; const chk=threshold!==NO_THRESHOLD; const list=[];
  for(let y=0;y<H;y++) for(let x=0,i=y*W;x<W;x++,i++){
    const v=px[i]; const vT=isEDM?this.trueEdmHeight(x,y,px):v;
    if(v===gMin) continue; if(exEdges&&(x===0||x===W-1||y===0||y===H-1)) continue; if(chk&&v<threshold) continue;
    let isMax=true; const inner=(y!==0&&y!==H-1)&&(x!==0&&x!==W-1);
    for(let d=0;d<8;d++) if(inner||this.isWithin(x,y,d)){ const vn=px[i+this.dir[d]]; const vnT=isEDM?this.trueEdmHeight(x+DX[d],y+DY[d],px):vn; if(vn>v&&vnT>vT){ isMax=false; break; } }
    if(isMax){ types[i]=MAXIMUM; list.push(i); }
  }
  const vF=f32(2e9/(gMax-gMin));
  const arr=list.map(p=>{ const x=p%W, y=(p-x)/W; const fv=isEDM?this.trueEdmHeight(x,y,px):px[p]; return {v:(Math.trunc(f32(f32(fv-gMin)*vF))||0),p}; }); // Java (int)NaN = 0
  arr.sort((a,b)=>a.v-b.v||a.p-b.p);
  return arr.map(a=>a.p);
};
MF.prototype.analyze=function(px,types,maxPts,exEdges,isEDM,gMin,tolerance,strict,outputType,maxSortErr,excludeOnEdges,roiMask){
  const W=this.W,H=this.H,N=W*H; const pList=new Int32Array(N); const xy=[];
  const disp=outputType==="points"||outputType==="count"; const tol=f32(tolerance);
  const val=(x,y,o)=>isEDM?this.trueEdmHeight(x,y,px):px[o];
  for(let iMax=maxPts.length-1;iMax>=0;iMax--){
    let off0=maxPts[iMax]; if(types[off0]&PROCESSED) continue;
    let x0=off0%W, y0=(off0-x0)/W; let v0=val(x0,y0,off0); let sortErr;
    do{
      pList[0]=off0; types[off0]|=(EQUAL|LISTED); let len=1, li=0;
      let isEdgeMax=(x0===0||x0===W-1||y0===0||y0===H-1); sortErr=false; let maxPossible=true;
      let xE=x0, yE=y0, nE=1;
      do{
        const off=pList[li], x=off%W, y=(off-x)/W; const inner=(y!==0&&y!==H-1)&&(x!==0&&x!==W-1);
        for(let d=0;d<8;d++){ const off2=off+this.dir[d];
          if((inner||this.isWithin(x,y,d))&&(types[off2]&LISTED)===0){
            if(isEDM&&px[off2]<=0) continue;
            if(types[off2]&PROCESSED){ maxPossible=false; break; }
            const x2=x+DX[d], y2=y+DY[d]; const v2=val(x2,y2,off2);
            if(v2>f32(v0+maxSortErr)){ maxPossible=false; break; }
            else if(v2>=f32(v0-tol)){
              if(v2>v0){ sortErr=true; off0=off2; v0=v2; x0=x2; y0=y2; }
              pList[len++]=off2; types[off2]|=LISTED;
              if((x2===0||x2===W-1||y2===0||y2===H-1)&&(strict||v2>=v0)){ isEdgeMax=true; if(exEdges){ maxPossible=false; break; } }
              if(v2===v0){ types[off2]|=EQUAL; xE+=x2; yE+=y2; nE++; }
            }
          }
        }
        li++;
      }while(li<len);
      if(sortErr){ for(li=0;li<len;li++) types[pList[li]]=0; }
      else{
        const reset=~(maxPossible?LISTED:(LISTED|EQUAL)); xE/=nE; yE/=nE; let minD=1e20, nearest=0;
        for(li=0;li<len;li++){ const off=pList[li], x=off%W, y=(off-x)/W; types[off]&=reset; types[off]|=PROCESSED;
          if(maxPossible){ types[off]|=MAX_AREA; if(types[off]&EQUAL){ const d2=(xE-x)*(xE-x)+(yE-y)*(yE-y); if(d2<minD){ minD=d2; nearest=li; } } } }
        if(maxPossible){ const off=pList[nearest]; types[off]|=MAX_POINT;
          if(disp&&!(excludeOnEdges&&isEdgeMax)){ const x=off%W, y=(off-x)/W; if(!roiMask||roiMask[off]) xy.push([x,y]); } }
      }
    }while(sortErr);
  }
  if(maxPts.length===0) types.fill(PROCESSED|MAX_AREA);
  return xy;
};
MF.prototype.make8bit=function(px,types,isEDM,gMin,gMax,threshold){
  const N=this.W*this.H; let minValue;
  if(isEDM){ threshold=0.5; minValue=1; } else minValue=(threshold===NO_THRESHOLD)?gMin:threshold;
  const offset=minValue-(gMax-minValue)*(1/253/2-1e-6); let factor=253/(gMax-minValue); if(isEDM&&factor>1) factor=1;
  const out=new Uint8Array(N);
  for(let i=0;i<N;i++){ const raw=px[i];
    if(threshold!==NO_THRESHOLD&&raw<threshold) out[i]=0;
    else if(types[i]&MAX_AREA) out[i]=255;
    else { const v=1+Math.round((raw-offset)*factor); out[i]=v<1?1:v<=254?v:254; } }
  return out;
};
MF.prototype.cleanupMaxima=function(pix,types,maxPts){
  const W=this.W,H=this.H; const pList=new Int32Array(W*H);
  for(let iMax=maxPts.length-1;iMax>=0;iMax--){
    const off0=maxPts[iMax]; if(types[off0]&(MAX_AREA|ELIMINATED)) continue;
    const level=pix[off0]; let lo=level+1; pList[0]=off0; types[off0]|=LISTED; let len=1,last=1,li=0, saddle=false;
    while(!saddle&&lo>0){ lo--; last=len; li=0;
      do{ const off=pList[li], x=off%W, y=(off-x)/W; const inner=(y!==0&&y!==H-1)&&(x!==0&&x!==W-1);
        for(let d=0;d<8;d++){ const off2=off+this.dir[d];
          if((inner||this.isWithin(x,y,d))&&(types[off2]&LISTED)===0){
            if((types[off2]&MAX_AREA)||((types[off2]&ELIMINATED)&&pix[off2]>=lo)){ saddle=true; break; }
            else if(pix[off2]>=lo&&(types[off2]&ELIMINATED)===0){ pList[len++]=off2; types[off2]|=LISTED; }
          } }
        if(saddle) break; li++;
      }while(li<len);
    }
    for(li=0;li<len;li++) types[pList[li]]&=~LISTED;
    for(li=0;li<last;li++){ const off=pList[li]; pix[off]=lo; types[off]|=ELIMINATED; }
  }
};
MF.prototype.fateTable=function(){
  const t=new Int32Array(256), s=new Array(8);
  for(let item=0;item<256;item++){
    for(let i=0,m=1;i<8;i++,m*=2) s[i]=(item&m)===m;
    for(let i=0,m=1;i<8;i++,m*=2) if(s[(i+4)%8]) t[item]|=m;
    for(let i=0;i<8;i+=2) if(s[i]){ s[(i+1)%8]=true; s[(i+7)%8]=true; }
    let tr=0; for(let i=0;i<8;i++) if(s[i]!==s[(i+1)%8]) tr++;
    if(tr>=4) t[item]=0;
  }
  return t;
};
MF.prototype.processLevel=function(pass,pix,table,levelStart,n,coords,setList){
  const W=this.W,H=this.H, xm=W-1, ym=H-1; let nCh=0,nUn=0;
  for(let i=0,p=levelStart;i<n;i++,p++){
    const c=coords[p], x=c&this.xMask, y=c>>>this.shift, off=x+y*W; let idx=0;
    if(y>0&&pix[off-W]===255) idx^=1; if(x<xm&&y>0&&pix[off-W+1]===255) idx^=2; if(x<xm&&pix[off+1]===255) idx^=4;
    if(x<xm&&y<ym&&pix[off+W+1]===255) idx^=8; if(y<ym&&pix[off+W]===255) idx^=16; if(x>0&&y<ym&&pix[off+W-1]===255) idx^=32;
    if(x>0&&pix[off-1]===255) idx^=64; if(x>0&&y>0&&pix[off-W-1]===255) idx^=128;
    const m=1<<pass; if((table[idx]&m)===m) setList[nCh++]=off; else coords[levelStart+(nUn++)]=c;
  }
  for(let i=0;i<nCh;i++) pix[setList[i]]=255;
  return nCh;
};
MF.prototype.watershedSegment=function(pix){
  const W=this.W,H=this.H,N=W*H; const hist=new Int32Array(256); for(let i=0;i<N;i++) hist[pix[i]]++;
  const size=N-hist[0]-hist[255]; const coords=new Int32Array(Math.max(1,size)); let highest=0, maxBin=0, off=0; const levelStart=new Int32Array(256);
  for(let v=1;v<255;v++){ levelStart[v]=off; off+=hist[v]; if(hist[v]>0) highest=v; if(hist[v]>maxBin) maxBin=hist[v]; }
  const lvOff=new Int32Array(highest+1);
  for(let y=0,i=0;y<H;y++) for(let x=0;x<W;x++,i++){ const v=pix[i]; if(v>0&&v<255){ off=levelStart[v]+lvOff[v]; coords[off]=x|(y<<this.shift); lvOff[v]++; } }
  const setList=new Int32Array(Math.max(1,Math.min(maxBin,Math.trunc((N+2)/3))));
  const table=this.fateTable(); const seq=[7,3,1,5,0,4,2,6];
  for(let level=highest;level>=1;level--){
    let remaining=hist[level], idle=0;
    while(remaining>0&&idle<8){ let dI=0;
      do{ const n=this.processLevel(seq[dI%8],pix,table,levelStart[level],remaining,coords,setList); remaining-=n; if(n>0) idle=0; dI++; }while(remaining>0&&idle++<8);
    }
    if(remaining>0&&level>1){
      let next=level; do next--; while(next>1&&hist[next]===0);
      if(next>0){ let newEnd=levelStart[next]+hist[next];
        for(let i=0,p=levelStart[level];i<remaining;i++,p++){ const c=coords[p], x=c&this.xMask, y=c>>>this.shift, po=x+y*W; let add=false;
          if(x===0||y===0||x===W-1||y===H-1) add=true; else for(let d=0;d<8;d++) if(this.isWithin(x,y,d)&&pix[po+this.dir[d]]===0){ add=true; break; }
          if(add) coords[newEnd++]=c; }
        hist[next]=newEnd-levelStart[next]; }
    }
  }
};
MF.prototype.nRadii=function(pix,x,y){
  const W=this.W,H=this.H, off=x+y*W; let tr=0, prev=true, first=true; const inner=(y!==0&&y!==H-1)&&(x!==0&&x!==W-1);
  for(let d=0;d<8;d++){ let set=prev;
    if(inner||this.isWithin(x,y,d)){ const isSet=pix[off+this.dir[d]]!==255; if((d&1)===0) set=isSet; else if(!isSet) set=false; }
    else set=true;
    if(set&&!prev) tr++; prev=set; if(d===0) first=set; }
  if(first&&!prev) tr++; return tr;
};
MF.prototype.removeLineFrom=function(pix,x,y){
  const W=this.W,H=this.H; pix[x+W*y]=255; let cont;
  do{ cont=false; const inner=(y!==0&&y!==H-1)&&(x!==0&&x!==W-1);
    for(let d=0;d<8;d+=2) if(inner||this.isWithin(x,y,d)){ const v=pix[x+W*y+this.dir[d]];
      if(v!==255&&v!==0){ const nr=this.nRadii(pix,x+DX[d],y+DY[d]); if(nr<=1){ x+=DX[d]; y+=DY[d]; pix[x+W*y]=255; cont=nr===1; break; } } }
  }while(cont);
};
MF.prototype.cleanupExtraLines=function(pix){
  const W=this.W,H=this.H;
  for(let y=0,i=0;y<H;y++) for(let x=0;x<W;x++,i++){ const v=pix[i]; if(v!==255&&v!==0){ const nr=this.nRadii(pix,x,y); if(nr===0) pix[i]=255; else if(nr===1) this.removeLineFrom(pix,x,y); } }
};
MF.prototype.deleteEdgeParticles=function(pix,types){
  const W=this.W,H=this.H; const del=(x,y)=>{ const i=x+y*W; if((types[i]&MAX_AREA)&&pix[i]!==0){ const b=IC.componentFilledBox(pix,W,H,x,y,v=>v===255); for(let yy=0;yy<b.h;yy++) for(let xx=0;xx<b.w;xx++) if(b.m[xx+yy*b.w]) pix[(b.y0+yy)*W+b.x0+xx]=0; } };
  for(let x=0;x<W;x++){ del(x,0); del(x,H-1); } for(let y=1;y<H-1;y++){ del(0,y); del(W-1,y); }
};

/* 8-연결 성분 + 그 안의 구멍(= ImageJ 가 traced polygon 을 채운 것과 같은 화소 집합)
   상자 크기 결과 {x0,y0,w,h,m}: 입자마다 영상 전체 크기 배열을 만들지 않음 (입자 많은 큰 영상에서 수십 배 빠름) */
let CF_STAMP=null, CF_GEN=0;
IC.componentFilledBox=function(pix,W,H,x0,y0,inside){
  if(!CF_STAMP||CF_STAMP.length!==W*H){ CF_STAMP=new Int32Array(W*H); CF_GEN=0; } const gen=++CF_GEN; if(gen>2e9){ CF_STAMP.fill(0); CF_GEN=1; }
  const G0=CF_GEN, st=[x0+y0*W], list=[]; CF_STAMP[x0+y0*W]=G0; let minx=x0,maxx=x0,miny=y0,maxy=y0;
  while(st.length){ const i=st.pop(), x=i%W, y=(i-x)/W; list.push(i); if(x<minx)minx=x; if(x>maxx)maxx=x; if(y<miny)miny=y; if(y>maxy)maxy=y;
    for(let dy=-1;dy<=1;dy++) for(let dx=-1;dx<=1;dx++){ const xx=x+dx, yy=y+dy; if(xx<0||yy<0||xx>=W||yy>=H) continue; const j=xx+yy*W; if(CF_STAMP[j]!==G0&&inside(pix[j])){ CF_STAMP[j]=G0; st.push(j); } } }
  // 바깥(경계 상자 테두리 한 칸 밖에서 4-연결로 닿는 비성분)을 지우고 나머지 = 채운 영역
  const bw=maxx-minx+3, bh=maxy-miny+3, comp=new Uint8Array(bw*bh); for(const i of list){ const x=i%W, y=(i-x)/W; comp[(x-minx+1)+(y-miny+1)*bw]=1; }
  const out=new Uint8Array(bw*bh); const q=[0]; out[0]=1;
  while(q.length){ const i=q.pop(), x=i%bw, y=(i-x)/bw;
    for(let d=0;d<4;d++){ const xx=x+(d===0?-1:d===1?1:0), yy=y+(d===2?-1:d===3?1:0); if(xx<0||yy<0||xx>=bw||yy>=bh) continue; const j=xx+yy*bw; if(out[j]||comp[j]) continue; out[j]=1; q.push(j); } }
  const w=maxx-minx+1, h=maxy-miny+1, m=new Uint8Array(w*h); for(let y=0;y<h;y++) for(let x=0;x<w;x++) if(!out[(x+1)+(y+1)*bw]) m[x+y*w]=1;
  return {x0:minx,y0:miny,w,h,m};
};
IC.componentFilled=function(pix,W,H,x0,y0,inside){ const b=IC.componentFilledBox(pix,W,H,x0,y0,inside); const res=new Uint8Array(W*H); for(let y=0;y<b.h;y++) for(let x=0;x<b.w;x++) if(b.m[x+y*b.w]) res[(b.y0+y)*W+b.x0+x]=1; return res; };

/* p: plane. prominence: 허용치. o: {output, strict(false), excludeEdges(false), light(false), isEDM(false), threshold(null), mask(Uint8Array)} */
IC.findMaxima=function(p,prominence,o){
  o=o||{}; const W=p.w,H=p.h; let px;
  if(o.light){ px=new Float32Array(W*H); if(p.bit===8){ for(let i=0;i<px.length;i++) px[i]=255-p.data[i]; } else { const [mn,mx]=IC.minMax(p.data); for(let i=0;i<px.length;i++) px[i]=f32(mn+mx-p.data[i]); } }
  else px=p.data;
  const mf=new MF(W,H);
  const thr=(o.threshold!=null&&!o.light)?o.threshold:NO_THRESHOLD;
  const isFloat=p.bit===32;
  const r=mf.findMaxima(px,isFloat,+prominence||0,!!o.strict,thr,o.output||"points",!!o.excludeEdges,!!o.isEDM,o.mask||null);
  return r;
};
IC.findMaximaImage=function(p,prominence,o){ const r=IC.findMaxima(p,prominence,o); const d=new Float32Array(r.image.length); for(let i=0;i<d.length;i++) d[i]=r.image[i]; return {w:p.w,h:p.h,bit:8,data:d}; };

})();
