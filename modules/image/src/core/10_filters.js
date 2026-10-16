/* ============================================================
   Process ▸ Filters  —  Gaussian · Rank(Mean/Median/Min/Max/Variance/Outliers/Top Hat)
   · 3×3 (Smooth/Sharpen/Find Edges) · Convolve · Unsharp Mask · Subtract Background
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const f32=Math.fround;

/* ---------- Gaussian Blur  (ij/plugin/filter/GaussianBlur.java) ----------
   accuracy: 8-bit·RGB 0.002, 16/32-bit 0.0002 (ImageJ 기본). 가장자리 밖은 가장 가까운 가장자리 화소.
   sigma > 8.5 이면 ImageJ 처럼 축소→합성곱→확대 경로를 쓴다(결과가 ImageJ 와 같도록). */
function makeGaussianKernel(sigma,accuracy,maxRadius){
  let kRadius=Math.ceil(sigma*Math.sqrt(-2*Math.log(accuracy)))+1;
  if(maxRadius<50) maxRadius=50;
  if(kRadius>maxRadius) kRadius=maxRadius;
  const k0=new Float32Array(kRadius), k1=new Float32Array(kRadius);
  for(let i=0;i<kRadius;i++) k0[i]=Math.exp(-0.5*i*i/sigma/sigma);
  if(kRadius<maxRadius&&kRadius>3){
    let sqrtSlope=Number.MAX_VALUE; let r=kRadius;
    while(r>(kRadius>>1)){ r--; const a=Math.sqrt(k0[r])/(kRadius-r); if(a<sqrtSlope) sqrtSlope=a; else break; }
    for(let r1=r+2;r1<kRadius;r1++) k0[r1]=(kRadius-r1)*(kRadius-r1)*sqrtSlope*sqrtSlope;
  }
  let sum;
  if(kRadius<maxRadius){ sum=k0[0]; for(let i=1;i<kRadius;i++) sum+=2*k0[i]; }
  else sum=sigma*Math.sqrt(2*Math.PI);
  let rsum=0.5+0.5*k0[0]/sum;
  for(let i=0;i<kRadius;i++){ const v=k0[i]/sum; k0[i]=v; rsum-=v; k1[i]=rsum; }
  return [k0,k1];
}
IC.makeGaussianKernel=makeGaussianKernel;
function convolveLine(input,pixels,kernel,readFrom,readTo,writeFrom,writeTo,point0,pointInc){
  const length=input.length, first=input[0], last=input[length-1];
  const kern=kernel[0], kern0=kern[0], kernSum=kernel[1], kRadius=kern.length;
  const firstPart=kRadius<length?kRadius:length;
  let p=point0+writeFrom*pointInc, i=writeFrom;
  for(;i<firstPart;i++,p+=pointInc){
    let result=input[i]*kern0; result+=kernSum[i]*first;
    if(i+kRadius>length) result+=kernSum[length-i-1]*last;
    for(let k=1;k<kRadius;k++){ let v=0; if(i-k>=0) v+=input[i-k]; if(i+k<length) v+=input[i+k]; result+=kern[k]*v; }
    pixels[p]=result;
  }
  const iEndInside=length-kRadius<writeTo?length-kRadius:writeTo;
  for(;i<iEndInside;i++,p+=pointInc){ let result=input[i]*kern0; for(let k=1;k<kRadius;k++) result+=kern[k]*(input[i-k]+input[i+k]); pixels[p]=result; }
  for(;i<writeTo;i++,p+=pointInc){
    let result=input[i]*kern0;
    if(i<kRadius) result+=kernSum[i]*first;
    if(i+kRadius>=length) result+=kernSum[length-i-1]*last;
    for(let k=1;k<kRadius;k++){ let v=0; if(i-k>=0) v+=input[i-k]; if(i+k<length) v+=input[i+k]; result+=kern[k]*v; }
    pixels[p]=result;
  }
}
function makeDownscaleKernel(unitLength){
  const mid=(unitLength*3)>>1, kernel=new Float32Array(3*unitLength);
  for(let i=0;i<=(unitLength>>1);i++){ const x=i/unitLength; const v=f32((0.75-x*x)/unitLength); kernel[mid-i]=v; kernel[mid+i]=v; }
  for(let i=(unitLength>>1)+1;i<((unitLength*3+1)>>1);i++){ const x=i/unitLength; const v=f32((0.125+0.5*(x-1)*(x-2))/unitLength); kernel[mid-i]=v; kernel[mid+i]=v; }
  return kernel;
}
function makeUpscaleKernel(unitLength){
  const kernel=new Float32Array(4*unitLength), mid=2*unitLength; kernel[0]=0;
  for(let i=0;i<unitLength;i++){ const x=i/unitLength; const v=f32(2/3-x*x*(1-0.5*x)); kernel[mid+i]=v; kernel[mid-i]=v; }
  for(let i=unitLength;i<2*unitLength;i++){ const x=i/unitLength; const v=f32((2-x)*(2-x)*(2-x)/6); kernel[mid+i]=v; kernel[mid-i]=v; }
  return kernel;
}
function downscaleLine(pixels,cache,kernel,reduceBy,pixel0,unscaled0,length,pointInc,newLength){
  let p=pixel0+pointInc*(unscaled0-((reduceBy*3)>>1));
  const pLast=pixel0+pointInc*(length-1);
  for(let xout=-1;xout<=newLength;xout++){
    let sum0=0,sum1=0,sum2=0;
    for(let x=0;x<reduceBy;x++,p+=pointInc){
      const v=pixels[p<pixel0?pixel0:(p>pLast?pLast:p)];
      sum0+=v*kernel[x+2*reduceBy]; sum1+=v*kernel[x+reduceBy]; sum2+=v*kernel[x];
    }
    if(xout>0) cache[xout-1]+=sum0;
    if(xout>=0&&xout<newLength) cache[xout]+=sum1;
    if(xout+1<newLength) cache[xout+1]=sum2;
  }
}
function upscaleLine(cache,pixels,kernel,reduceBy,pixel0,unscaled0,writeFrom,writeTo,pointInc){
  let p=pixel0+pointInc*writeFrom;
  for(let xout=writeFrom;xout<writeTo;xout++,p+=pointInc){
    const xin=Math.trunc((xout-unscaled0+reduceBy-1)/reduceBy);
    const x=reduceBy-1-((xout-unscaled0+reduceBy-1)%reduceBy);
    pixels[p]=cache[xin-2]*kernel[x]+cache[xin-1]*kernel[x+reduceBy]+cache[xin]*kernel[x+2*reduceBy]+cache[xin+1]*kernel[x+3*reduceBy];
  }
}
function blur1Direction(pixels,width,height,sigma,accuracy,xDirection){
  const UPSCALE_K_RADIUS=2, MIN_DOWNSCALED_SIGMA=4;
  const length=xDirection?width:height, pointInc=xDirection?1:width, lineInc=xDirection?width:1;
  const lineFrom=0, lineTo=xDirection?height:width, writeFrom=0, writeTo=length;
  const doDownscaling=sigma>2*MIN_DOWNSCALED_SIGMA+0.5;
  const reduceBy=doDownscaling?Math.min(Math.floor(sigma/MIN_DOWNSCALED_SIGMA),length):1;
  const sigmaGauss=doDownscaling?Math.sqrt(sigma*sigma/(reduceBy*reduceBy)-1/3-1/4):sigma;
  const maxLength=doDownscaling?Math.trunc((length+reduceBy-1)/reduceBy)+2*(UPSCALE_K_RADIUS+1):length;
  const gaussKernel=makeGaussianKernel(sigmaGauss,accuracy,maxLength);
  const kRadius=gaussKernel[0].length;
  const readFrom=(writeFrom-kRadius<0)?0:writeFrom-kRadius;
  const readTo=(writeTo+kRadius>length)?length:writeTo+kRadius;
  const newLength=doDownscaling?Math.trunc((readTo-readFrom+reduceBy-1)/reduceBy)+2*(UPSCALE_K_RADIUS+1):length;
  const unscaled0=readFrom-(UPSCALE_K_RADIUS+1)*reduceBy;
  const downK=doDownscaling?makeDownscaleKernel(reduceBy):null, upK=doDownscaling?makeUpscaleKernel(reduceBy):null;
  const cache1=new Float32Array(newLength), cache2=doDownscaling?new Float32Array(newLength):null;
  for(let line=lineFrom;line<lineTo;line++){
    const pixel0=line*lineInc;
    if(doDownscaling){
      downscaleLine(pixels,cache1,downK,reduceBy,pixel0,unscaled0,length,pointInc,newLength);
      convolveLine(cache1,cache2,gaussKernel,0,newLength,1,newLength-1,0,1);
      upscaleLine(cache2,pixels,upK,reduceBy,pixel0,unscaled0,writeFrom,writeTo,pointInc);
    }else{
      let p=pixel0+readFrom*pointInc;
      for(let i=readFrom;i<readTo;i++,p+=pointInc) cache1[i]=pixels[p];
      convolveLine(cache1,pixels,gaussKernel,readFrom,readTo,writeFrom,writeTo,pixel0,pointInc);
    }
  }
}
/* sigma 는 화소 단위. sigmaY 생략 시 sigma. accuracy 생략 시 ImageJ 기본값 */
IC.gaussianFloat=function(data,w,h,sx,sy,accuracy){
  if(sx>0) blur1Direction(data,w,h,sx,accuracy,true);
  if(sy>0) blur1Direction(data,w,h,sy,accuracy,false);
  return data;
};
IC.gaussian=function(p,o){
  const sx=+o.sigma, sy=o.sigmaY!=null?+o.sigmaY:sx;
  const acc=o.accuracy||((p.bit===8||o.rgb)?0.002:0.0002);
  const d=Float32Array.from(p.data);
  IC.gaussianFloat(d,p.w,p.h,sx,sy,acc);
  return {w:p.w,h:p.h,bit:p.bit,data:IC.roundTo(d,p.bit)};
};

/* ---------- Rank filters (ij/plugin/filter/RankFilters.java) ----------
   원형 커널: r2=(int)(r*r)+1, kRadius=(int)sqrt(r2), 줄 dy 의 반폭 dx=(int)sqrt(r2-dy*dy).
   (1.5≤r<1.75 → 1.75, 2.5≤r<2.85 → 2.85 — ImageJ 와 같은 크기)
   가장자리 밖 화소 = 가장 가까운 가장자리 화소. 8/16-bit 결과는 (int)(v+0.5) */
IC.lineRadii=function(radius){
  if(radius>=1.5&&radius<1.75) radius=1.75; else if(radius>=2.5&&radius<2.85) radius=2.85;
  const r2=Math.trunc(radius*radius)+1, kR=Math.trunc(Math.sqrt(r2+1e-10));
  const dx=new Int32Array(2*kR+1); let n=0;
  for(let y=-kR;y<=kR;y++){ const d=y===0?kR:Math.trunc(Math.sqrt(r2-y*y+1e-10)); dx[y+kR]=d; n+=2*d+1; }
  return {kR,dx,n};
};
function padded(p,kR){ // (w+2kR)×(h+2kR), 가장자리 복제
  const W=p.w, H=p.h, PW=W+2*kR, PH=H+2*kR, out=new Float32Array(PW*PH), d=p.data;
  for(let y=0;y<PH;y++){ const sy=Math.min(H-1,Math.max(0,y-kR)); const o=y*PW, so=sy*W;
    for(let x=0;x<PW;x++){ const sx=Math.min(W-1,Math.max(0,x-kR)); out[o+x]=d[so+sx]; } }
  return {d:out,PW,PH};
}
/* van Herk / Gil-Werman 1D 이동 최소·최대 (창 반폭 r) — 길이 n 의 a 를 out 에 */
function slideMinMax(a,off,n,r,out,isMax){
  const w=2*r+1; if(r===0){ for(let i=0;i<n;i++) out[i]=a[off+i]; return; }
  const g=new Float32Array(n), hh=new Float32Array(n);
  for(let i=0;i<n;i++){ const v=a[off+i]; if(i%w===0) g[i]=v; else g[i]=isMax?Math.max(g[i-1],v):Math.min(g[i-1],v); }
  for(let i=n-1;i>=0;i--){ const v=a[off+i]; if(i===n-1||(i+1)%w===0) hh[i]=v; else hh[i]=isMax?Math.max(hh[i+1],v):Math.min(hh[i+1],v); }
  for(let i=r;i<n-r;i++){ const a1=hh[i-r], b1=g[i+r]; out[i]=isMax?(a1>b1?a1:b1):(a1<b1?a1:b1); }
}
function rankMinMax(p,radius,isMax){
  const {kR,dx}=IC.lineRadii(radius); const P=padded(p,kR);
  { const d=P.d, fill=isMax?-Infinity:Infinity; for(let i=0;i<d.length;i++) if(d[i]!==d[i]) d[i]=fill; } /* NaN → 비교에서 빠지게 */ const W=p.w,H=p.h,PW=P.PW;
  const out=new Float32Array(W*H); const tmp=new Float32Array(PW); const acc=new Float32Array(W);
  for(let y=0;y<H;y++){
    acc.fill(isMax?-Infinity:Infinity);
    for(let k=-kR;k<=kR;k++){ const r=dx[k+kR]; slideMinMax(P.d,(y+kR+k)*PW,PW,r,tmp,isMax);
      for(let x=0;x<W;x++){ const v=tmp[x+kR]; if(isMax?v>acc[x]:v<acc[x]) acc[x]=v; } }
    out.set(acc,y*W);
    if((y&63)===0) IC.progress(y/H);
  }
  return out;
}
function rankSums(p,radius){ // 평균·분산용 합 (double). 비유한값(NaN·±∞)은 개수만 세어 그 창의 결과를 ImageJ 처럼 NaN/±∞ 로
  const {kR,dx,n}=IC.lineRadii(radius); const P=padded(p,kR); const W=p.w,H=p.h,PW=P.PW,PH=P.PH;
  let anyNF=false; for(let i=0;i<P.d.length;i++){ const v=P.d[i]; if(v!==v||v===Infinity||v===-Infinity){ anyNF=true; break; } }
  let cN=null,cP=null,cM=null; if(anyNF){ const R1=PW+1; cN=new Int32Array(PH*R1); cP=new Int32Array(PH*R1); cM=new Int32Array(PH*R1);
    for(let y=0;y<PH;y++){ let a=0,b=0,c=0; for(let x=0;x<PW;x++){ const v=P.d[y*PW+x]; if(v!==v) a++; else if(v===Infinity) b++; else if(v===-Infinity) c++; cN[y*R1+x+1]=a; cP[y*R1+x+1]=b; cM[y*R1+x+1]=c; if(v!==v||v===Infinity||v===-Infinity) P.d[y*PW+x]=0; } } }
  const c1=new Float64Array(PW*PH+PH), c2=new Float64Array(PW*PH+PH); // 행별 누적합 (길이 PW+1)
  const R=PW+1;
  for(let y=0;y<PH;y++){ let s=0,t=0; c1[y*R]=0; c2[y*R]=0; for(let x=0;x<PW;x++){ const v=P.d[y*PW+x]; s+=v; t+=v*v; c1[y*R+x+1]=s; c2[y*R+x+1]=t; } }
  const S1=new Float64Array(W*H), S2=new Float64Array(W*H);
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){ let s=0,t=0; const cx=x+kR;
    for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*R; s+=c1[row+cx+r+1]-c1[row+cx-r]; t+=c2[row+cx+r+1]-c2[row+cx-r]; }
    if(anyNF){ let a=0,b=0,c=0; for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*R; a+=cN[row+cx+r+1]-cN[row+cx-r]; b+=cP[row+cx+r+1]-cP[row+cx-r]; c+=cM[row+cx+r+1]-cM[row+cx-r]; }
      if(a||(b&&c)){ s=NaN; t=NaN; } else if(b){ s=Infinity; t=NaN; } else if(c){ s=-Infinity; t=NaN; } }
    S1[y*W+x]=s; S2[y*W+x]=t; }
  return {S1,S2,n};
}
function quickSelect(a,k,n){ let l=0,r=n-1; while(l<r){ const piv=a[(l+r)>>1]; let i=l,j=r; while(i<=j){ while(a[i]<piv) i++; while(a[j]>piv) j--; if(i<=j){ const t=a[i]; a[i]=a[j]; a[j]=t; i++; j--; } } if(j<k) l=i; if(k<i) r=j; } return a[k]; }
function rankMedian(p,radius){
  const {kR,dx,n}=IC.lineRadii(radius); const P=padded(p,kR); const W=p.w,H=p.h,PW=P.PW;
  const out=new Float32Array(W*H), half=n>>1;
  if(p.bit===8){ // 이동 히스토그램 (정수 0..255)
    const hist=new Int32Array(256);
    for(let y=0;y<H;y++){
      hist.fill(0);
      for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR; for(let xx=-r;xx<=r;xx++) hist[P.d[row+xx]]++; }
      let m=0, below=0; // below = m 보다 작은 개수
      const fix=()=>{ while(below+hist[m]<=half){ below+=hist[m]; m++; } while(below>half){ m--; below-=hist[m]; } };
      fix(); out[y*W]=m;
      for(let x=1;x<W;x++){
        for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; const vo=P.d[row-r-1], vi=P.d[row+r];
          hist[vo]--; if(vo<m) below--; hist[vi]++; if(vi<m) below++; }
        fix(); out[y*W+x]=m;
      }
      if((y&63)===0) IC.progress(y/H);
    }
    return out;
  }
  const buf=new Float64Array(n), ab=new Float64Array(n), bl=new Float64Array(n);
  for(let y=0;y<H;y++){ for(let x=0;x<W;x++){ let c=0; for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; for(let xx=-r;xx<=r;xx++) buf[c++]=P.d[row+xx]; }
      // ImageJ getNaNAwareMedian: 추정값 g = 가운데 화소(NaN 이면 결과 NaN). NaN 은 뺀 개수 m, half=⌊m/2⌋.
      // g 보다 큰 값이 half 개보다 많으면 그중 (nAbove−half−1) 번째, 작은 값이 half 개보다 많으면 그중 half 번째, 아니면 g
      const g=P.d[(y+kR)*PW+kR+x]; if(g!==g){ out[y*W+x]=NaN; continue; }
      let m=0,na=0,nb=0; for(let i=0;i<c;i++){ const v=buf[i]; if(v!==v) continue; m++; if(v>g) ab[na++]=v; else if(v<g) bl[nb++]=v; }
      const half=m>>1; out[y*W+x]=!m?NaN:na>half?quickSelect(ab,na-half-1,na):nb>half?quickSelect(bl,half,nb):g; }
    if((y&31)===0) IC.progress(y/H); }
  return out;
}
/* type: mean | median | min | max | variance | outliers | despeckle | tophat | open | close
   o.radius, o.threshold, o.which ('bright'|'dark'), o.light (Top Hat: light background) */
IC.rank=function(p,type,o){
  o=o||{}; let r=+o.radius; if(type==="despeckle"){ type="median"; r=1; }
  const W=p.w,H=p.h; let out;
  if(type==="mean"||type==="variance"){
    const {S1,S2,n}=rankSums(p,r); out=new Float32Array(W*H);
    if(type==="mean") for(let i=0;i<W*H;i++) out[i]=S1[i]/n;
    else { const mx=IC.maxOf(p.bit); for(let i=0;i<W*H;i++){ let v=f32((S2[i]-S1[i]*S1[i]/n)/n); if(v>mx) v=mx; if(v<0) v=0; out[i]=v; } } // NaN 은 비교를 지나 그대로 NaN
  }else if(type==="median") out=rankMedian(p,r);
  else if(type==="min") out=rankMinMax(p,r,false);
  else if(type==="max") out=rankMinMax(p,r,true);
  else if(type==="open"||type==="close"||type==="tophat"){
    // ImageJ: Open = Min→Max, Close = Max→Min.
    // Top Hat: 어두운 배경 = 원본 − Open, 밝은 배경 = 원본 − Close + (255 또는 65535; float 는 0)
    const light=!!o.light; const maxFirst=type==="close"||(type==="tophat"&&light);
    const a=rankMinMax(p,r,maxFirst);
    const b=rankMinMax({w:W,h:H,bit:p.bit,data:a},r,!maxFirst);
    if(type==="tophat"){ out=new Float32Array(W*H); const off=(light&&p.bit!==32)?IC.maxOf(p.bit):0;
      for(let i=0;i<W*H;i++) out[i]=p.data[i]-b[i]+off; }
    else out=b;
  }else if(type==="outliers"){
    const med=rankMedian(p,r), thr=+o.threshold||0, bright=(o.which||"bright")==="bright"; out=Float32Array.from(p.data);
    for(let i=0;i<W*H;i++){ const v=p.data[i]; if(bright?(v-thr>med[i]):(v+thr<med[i])) out[i]=med[i]; }
  }else throw new Error("알 수 없는 rank filter: "+type);
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};

/* ---------- 3×3 필터 (ImageProcessor.smooth/sharpen/findEdges) ----------
   가장자리 = 복제. 8-bit 는 ImageJ 정수 연산 그대로 */
IC.filter3=function(p,type){
  const W=p.w,H=p.h,d=p.data,out=new Float32Array(W*H);
  const g=(x,y)=>d[(y<0?0:y>=H?H-1:y)*W+(x<0?0:x>=W?W-1:x)];
  const b8=p.bit===8||p.bit===16; // 8·16-bit: ImageJ 정수 연산(Byte/ShortProcessor.filter3x3)
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){
    const p1=g(x-1,y-1),p2=g(x,y-1),p3=g(x+1,y-1),p4=g(x-1,y),p5=g(x,y),p6=g(x+1,y),p7=g(x-1,y+1),p8=g(x,y+1),p9=g(x+1,y+1);
    let v;
    if(type==="smooth"){ const s=p1+p2+p3+p4+p5+p6+p7+p8+p9; v=b8?Math.trunc((s+4)/9):s/9; }
    else if(type==="sharpen"){ const s=12*p5-(p1+p2+p3+p4+p6+p7+p8+p9); v=b8?Math.trunc((s+2)/4):s/4; }
    else if(type==="edges"){ const s1=p1+2*p2+p3-p7-2*p8-p9, s2=p1+2*p4+p7-p3-2*p6-p9; v=Math.sqrt(s1*s1+s2*s2); if(b8) v=Math.trunc(v); }
    else throw new Error("3×3 필터 종류: "+type);
    out[y*W+x]=v;
  }
  if(b8){ const M=IC.maxOf(p.bit); for(let i=0;i<out.length;i++){ if(out[i]>M) out[i]=M; if(out[i]<0) out[i]=0; } return {w:W,h:H,bit:p.bit,data:out}; }
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};

/* ---------- Convolve (ij/plugin/filter/Convolver.java) ----------
   kernel: {kw, kh, k:Array(kw*kh)} (홀수 크기), normalize: 합이 0 이 아니면 합으로 나눔. 가장자리 = 복제 */
IC.convolve=function(p,o){
  const {kw,kh}=o, k=o.k; if(!(kw%2&&kh%2)) throw new Error("커널 크기는 홀수여야 합니다");
  let scale=1; if(o.normalize!==false){ let s=0; for(const v of k) s+=v; if(s!==0) scale=1/s; }
  const W=p.w,H=p.h,d=p.data,out=new Float32Array(W*H), uc=kw>>1, vc=kh>>1;
  for(let y=0;y<H;y++){ for(let x=0;x<W;x++){ let s=0, i=0;
      for(let v=-vc;v<=vc;v++){ const yy=y+v<0?0:y+v>=H?H-1:y+v; const row=yy*W;
        for(let u=-uc;u<=uc;u++){ const xx=x+u<0?0:x+u>=W?W-1:x+u; s+=d[row+xx]*k[i++]; } }
      out[y*W+x]=s*scale; }
    if((y&63)===0) IC.progress(y/H); }
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};

/* ---------- Unsharp Mask (ij/plugin/filter/UnsharpMask.java) ----------
   blur(accuracy 0.01) 후 (orig − w·blur)/(1 − w) */
IC.unsharp=function(p,o){
  const sigma=+o.sigma, wt=f32(+o.weight);
  const b=Float32Array.from(p.data); IC.gaussianFloat(b,p.w,p.h,sigma,sigma,0.01);
  const out=new Float32Array(b.length); for(let i=0;i<b.length;i++) out[i]=(p.data[i]-wt*b[i])/(1-wt);
  return {w:p.w,h:p.h,bit:p.bit,data:IC.roundTo(out,p.bit)};
};

/* ---------- Subtract Background — Rolling ball (ij/plugin/filter/BackgroundSubtracter.java) ----------
   ImageJ 와 같이: 반지름 > 10 이면 영상을 줄여(최솟값) 공을 굴리고 쌍선형으로 되돌린다.
   공은 위쪽 일부만 쓴다(arcTrim 24/32/40 %). 사전 3×3 평균(smoothing) 기본 켬.
   o: {radius, light, createBackground, smooth(기본 true)} */
function rollingBall(radius){
  let shrink, trim;
  if(radius<=10){ shrink=1; trim=24; } else if(radius<=30){ shrink=2; trim=24; } else if(radius<=100){ shrink=4; trim=32; } else { shrink=8; trim=40; }
  let sr=radius/shrink; if(sr<1) sr=1;
  const rsq=sr*sr, xtrim=Math.trunc(Math.trunc(trim*sr)/100), half=Math.round(sr-xtrim), width=2*half+1;
  const data=new Float32Array(width*width);
  for(let y=0,p=0;y<width;y++) for(let x=0;x<width;x++,p++){ const xv=x-half, yv=y-half; const t=rsq-xv*xv-yv*yv; data[p]=t>0?Math.sqrt(t):0; }
  return {data,width,shrink};
}
function shrinkImage(px,w,h,s){
  const sw=Math.trunc((w+s-1)/s), sh=Math.trunc((h+s-1)/s), out=new Float32Array(sw*sh);
  for(let ys=0;ys<sh;ys++) for(let xs=0;xs<sw;xs++){ let mn=Number.MAX_VALUE;
    for(let j=0,y=s*ys;j<s&&y<h;j++,y++) for(let k=0,x=s*xs;k<s&&x<w;k++,x++){ const v=px[x+y*w]; if(v<mn) mn=v; }
    out[xs+ys*sw]=mn; }
  return {d:out,w:sw,h:sh};
}
function rollBall(ball,pixels,width,height){
  const zBall=ball.data, bw=ball.width, radius=bw>>1, cache=new Float32Array(width*bw);
  const NEG=-3.4028234663852886e38, POS=3.4028234663852886e38;
  for(let y=-radius;y<height+radius;y++){
    const nextW=(y+radius)%bw, nextR=y+radius;
    if(nextR<height){ cache.set(pixels.subarray(nextR*width,nextR*width+width),nextW*width); for(let x=0,p=nextR*width;x<width;x++,p++) pixels[p]=NEG; }
    let y0=y-radius; if(y0<0) y0=0; const yBall0=y0-y+radius; let yend=y+radius; if(yend>=height) yend=height-1;
    for(let x=-radius;x<width+radius;x++){
      let z=POS; let x0=x-radius; if(x0<0) x0=0; const xBall0=x0-x+radius; let xend=x+radius; if(xend>=width) xend=width-1;
      for(let yp=y0,yb=yBall0;yp<=yend;yp++,yb++){ let cp=(yp%bw)*width+x0; for(let xp=x0,bp=xBall0+yb*bw;xp<=xend;xp++,cp++,bp++){ const zr=f32(cache[cp]-zBall[bp]); if(z>zr) z=zr; } }
      for(let yp=y0,yb=yBall0;yp<=yend;yp++,yb++) for(let xp=x0,p=xp+yp*width,bp=xBall0+yb*bw;xp<=xend;xp++,p++,bp++){ const zm=f32(z+zBall[bp]); if(pixels[p]<zm) pixels[p]=zm; }
    }
    if((y&31)===0) IC.progress(0.1+0.8*(y+radius)/(height+bw));
  }
}
function interpArrays(len,smallLen,s){
  const idx=new Int32Array(len), wt=new Float32Array(len);
  for(let i=0;i<len;i++){ let si=Math.trunc((i-(s>>1))/s); if(si>=smallLen-1) si=smallLen-2; idx[i]=si; const dist=f32(f32((i+0.5)/s)-(si+0.5)); wt[i]=f32(1-dist); }
  return {idx,wt};
}
function enlargeImage(sm,px,w,h,s){
  const X=interpArrays(w,sm.w,s), Y=interpArrays(h,sm.h,s), sp=sm.d, sw=sm.w;
  let line0=new Float32Array(w), line1=new Float32Array(w);
  for(let x=0;x<w;x++) line1[x]=f32(f32(sp[X.idx[x]]*X.wt[x])+f32(sp[X.idx[x]+1]*f32(1-X.wt[x])));
  let ys0=-1;
  for(let y=0;y<h;y++){
    if(ys0<Y.idx[y]){ const t=line0; line0=line1; line1=t; ys0++; const o=(Y.idx[y]+1)*sw;
      for(let x=0;x<w;x++) line1[x]=f32(f32(sp[o+X.idx[x]]*X.wt[x])+f32(sp[o+X.idx[x]+1]*f32(1-X.wt[x]))); }
    const wt=Y.wt[y], w1=f32(1-wt); for(let x=0,p=y*w;x<w;x++,p++) px[p]=f32(f32(line0[x]*wt)+f32(line1[x]*w1));
  }
}
function filter3mean(px,w,h){
  const f3=(p0,len,inc)=>{ let v3=px[p0],v2=v3,v1; for(let i=0,p=p0;i<len;i++,p+=inc){ v1=v2; v2=v3; if(i<len-1) v3=px[p+inc]; const m=f32(f32(v1+v2+v3)*0.333333333); if(m===m) px[p]=m; } };
  for(let y=0;y<h;y++) f3(y*w,w,1);
  for(let x=0;x<w;x++) f3(x,h,w);
}
IC.rollingBackground=function(p,o){
  const W=p.w,H=p.h, radius=+o.radius, invert=!!o.light, smooth=o.smooth!==false;
  const ball=rollingBall(radius);
  if(ball.shrink>1&&(Math.ceil(W/ball.shrink)<2||Math.ceil(H/ball.shrink)<2)) throw new Error(`Subtract Background: 영상(${W}×${H})이 반지름 ${radius} 에 비해 너무 얇습니다 — 반지름을 줄이세요 (ImageJ 도 처리하지 않음)`);
  const bg=Float32Array.from(p.data);
  if(invert) for(let i=0;i<bg.length;i++) bg[i]=-bg[i];
  if(smooth) filter3mean(bg,W,H);
  if(ball.shrink>1){ const sm=shrinkImage(bg,W,H,ball.shrink); rollBall(ball,sm.d,sm.w,sm.h); enlargeImage(sm,bg,W,H,ball.shrink); }
  else rollBall(ball,bg,W,H);
  if(invert) for(let i=0;i<bg.length;i++) bg[i]=-bg[i];
  return bg;
};
IC.subtractBackground=function(p,o){
  const bg=IC.rollingBackground(p,o); const W=p.w,H=p.h, invert=!!o.light;
  if(o.createBackground) return {w:W,h:H,bit:p.bit,data:IC.roundTo(bg,p.bit)};
  const out=new Float32Array(W*H);
  if(p.bit===32){ for(let i=0;i<out.length;i++) out[i]=p.data[i]-bg[i]; return {w:W,h:H,bit:32,data:out}; }
  const mx=p.bit===8?255:65535, off=f32(invert?mx+0.5:0.5);
  for(let i=0;i<out.length;i++){ let v=f32(f32(p.data[i]-bg[i])+off); if(v<0) v=0; if(v>mx) v=mx; out[i]=Math.floor(v); }
  return {w:W,h:H,bit:p.bit,data:out};
};

})();
