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
/* 1D 이동 최댓값·최솟값 (van Herk / Gil-Werman): dst[do+x] = max|min(a[off+x .. off+x+2d]), x<n. g·h 는 길이 ≥ n+2d 작업 버퍼 */
function slideMax(a,off,n,d,dst,dOff,g,h){
  if(d===0){ for(let x=0;x<n;x++) dst[dOff+x]=a[off+x]; return; }
  const w=2*d+1, L=n+2*d;
  for(let s=0;s<L;s+=w){ const e=s+w<L?s+w:L;
    let m=a[off+s]; g[s]=m; for(let i=s+1;i<e;i++){ const v=a[off+i]; if(v>m) m=v; g[i]=m; }
    m=a[off+e-1]; h[e-1]=m; for(let i=e-2;i>=s;i--){ const v=a[off+i]; if(v>m) m=v; h[i]=m; } }
  for(let x=0,j=2*d;x<n;x++,j++){ const u=h[x], v=g[j]; dst[dOff+x]=u>v?u:v; }
}
function slideMin(a,off,n,d,dst,dOff,g,h){
  if(d===0){ for(let x=0;x<n;x++) dst[dOff+x]=a[off+x]; return; }
  const w=2*d+1, L=n+2*d;
  for(let s=0;s<L;s+=w){ const e=s+w<L?s+w:L;
    let m=a[off+s]; g[s]=m; for(let i=s+1;i<e;i++){ const v=a[off+i]; if(v<m) m=v; g[i]=m; }
    m=a[off+e-1]; h[e-1]=m; for(let i=e-2;i>=s;i--){ const v=a[off+i]; if(v<m) m=v; h[i]=m; } }
  for(let x=0,j=2*d;x<n;x++,j++){ const u=h[x], v=g[j]; dst[dOff+x]=u<v?u:v; }
}
/* 원형 커널 최댓값·최솟값 (ImageJ RankFilters 처럼 커널 = 줄 조각들).
   패딩한 줄마다, 커널에 나오는 서로 다른 반폭 d 각각의 1D 이동 극값을 한 번만 구해 (2kR+1) 줄 고리 버퍼에 두고
   출력 줄마다 (2kR+1) 개를 합친다 → 화소당 O(커널 줄 수), 옛 구현(출력 줄마다 모든 커널 줄을 다시 계산)의 결과와 같다 */
IC._rankRingLimit=1<<25; // 고리 버퍼 최대 칸 수 (Float32) — 시험에서 0 으로 두면 버퍼 없는 경로
function rankMinMax(p,radius,isMax){
  const {kR,dx}=IC.lineRadii(radius); const P=padded(p,kR);
  { const d=P.d, fill=isMax?-Infinity:Infinity; for(let i=0;i<d.length;i++) if(d[i]!==d[i]) d[i]=fill; } /* NaN → 비교에서 빠지게 */
  const W=p.w,H=p.h,PW=P.PW,KH=2*kR+1, pd=P.d, slide=isMax?slideMax:slideMin;
  const ds=Array.from(new Set(dx)).sort((a,b)=>a-b), nd=ds.length, di=new Int32Array(KH); for(let k=0;k<KH;k++) di[k]=ds.indexOf(dx[k]);
  const g=new Float32Array(PW), h=new Float32Array(PW), out=new Float32Array(W*H);
  if(KH*nd*W<=IC._rankRingLimit){ // 고리 버퍼: 칸 s = 패딩 줄 py (py ≡ s mod KH) 의 반폭별 이동 극값
    const R=new Float32Array(KH*nd*W);
    const fill=py=>{ const slot=py%KH, src=py*PW; for(let j=0;j<nd;j++){ const d=ds[j]; slide(pd,src+kR-d,W,d,R,(slot*nd+j)*W,g,h); } };
    for(let py=0;py<KH-1;py++) fill(py);
    for(let y=0;y<H;y++){ fill(y+KH-1); const o=y*W;
      { const s0=(((y+kR)%KH)*nd+di[kR])*W; for(let x=0;x<W;x++) out[o+x]=R[s0+x]; } // 가운데 줄
      for(let k=0;k<kR;k++){ const s1=(((y+k)%KH)*nd+di[k])*W, s2=(((y+KH-1-k)%KH)*nd+di[KH-1-k])*W; // 위·아래 짝 줄 (같은 반폭)
        if(isMax){ for(let x=0;x<W;x++){ const a=R[s1+x], b=R[s2+x], c=out[o+x], m=a>b?a:b; out[o+x]=m>c?m:c; } }
        else { for(let x=0;x<W;x++){ const a=R[s1+x], b=R[s2+x], c=out[o+x], m=a<b?a:b; out[o+x]=m<c?m:c; } } }
      if((y&63)===0) IC.progress(y/H); }
  } else { // 아주 큰 커널: 버퍼 없이 출력 줄마다 계산
    const tmp=new Float32Array(W);
    for(let y=0;y<H;y++){ const o=y*W;
      for(let k=0;k<KH;k++){ const d=dx[k]; slide(pd,(y+k)*PW+kR-d,W,d,k?tmp:out,k?0:o,g,h);
        if(k){ if(isMax){ for(let x=0;x<W;x++) if(tmp[x]>out[o+x]) out[o+x]=tmp[x]; } else { for(let x=0;x<W;x++) if(tmp[x]<out[o+x]) out[o+x]=tmp[x]; } } }
      if((y&15)===0) IC.progress(y/H); }
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
function quickSelectI(a,k,n){ let l=0,r=n-1; while(l<r){ const piv=a[(l+r)>>1]; let i=l,j=r; while(i<=j){ while(a[i]<piv) i++; while(a[j]>piv) j--; if(i<=j){ const t=a[i]; a[i]=a[j]; a[j]=t; i++; j--; } } if(j<k) l=i; if(k<i) r=j; } return a[k]; }
/* 값 → 조밀한 순위 (같은 값 = 같은 순위, ±∞ 포함, NaN = −1, −0 = +0). 반환 {rk, U, vals(순위→값)}.
   float32 비트를 정렬 가능한 정수 열쇠로 바꿔 기수 정렬(11 비트 × 3) — Float64 정렬 + 이진 탐색보다 빠르다 */
function denseRanks(d){
  const N=d.length, f32=new Float32Array(1), u32=new Uint32Array(f32.buffer);
  let ia=new Uint32Array(N), ib=new Uint32Array(N), ka=new Uint32Array(N), kb=new Uint32Array(N); const cnt=new Int32Array(2048); let nn=0;
  for(let i=0;i<N;i++){ ia[i]=i; let v=d[i]; if(v!==v){ ka[i]=0xFFFFFFFF; continue; } if(v===0) v=0; f32[0]=v; const u=u32[0]; ka[i]=(u&0x80000000)?(~u>>>0):(u|0x80000000)>>>0; nn++; }
  for(let sh=0;sh<33;sh+=11){ cnt.fill(0); for(let i=0;i<N;i++) cnt[(ka[i]>>>sh)&2047]++;
    let s=0; for(let k=0;k<2048;k++){ const c=cnt[k]; cnt[k]=s; s+=c; }
    for(let i=0;i<N;i++){ const kk=ka[i], j=cnt[(kk>>>sh)&2047]++; kb[j]=kk; ib[j]=ia[i]; }
    let t=ia; ia=ib; ib=t; t=ka; ka=kb; kb=t; }
  const rk=new Int32Array(N), vals=new Float32Array(Math.max(1,nn)); let U=0, prev=-1;
  for(let j=0;j<N;j++){ const kk=ka[j], i=ia[j]; if(kk===0xFFFFFFFF&&d[i]!==d[i]){ rk[i]=-1; continue; }
    if(U===0||kk!==prev){ let v=d[i]; if(v===0) v=0; vals[U++]=v; prev=kk; } rk[i]=U-1; }
  return {rk,U,vals:vals.subarray(0,U)};
}
/* 3×3 중앙값 (반지름 1, Despeckle): 줄마다 세로 세 값을 한 번 정렬(lo·mid·hi)해 세 열이 함께 쓰고,
   중앙값 = med3(max lo, med3 mid, min hi) — 9 개의 참 중앙값. NaN 이 든 창만 medianAt (ImageJ 규칙) */
function median3x3(p,P,out){
  const W=p.w,H=p.h,PW=P.PW,d=P.d, lo=new Float32Array(PW), mi=new Float32Array(PW), hi=new Float32Array(PW);
  for(let y=0;y<H;y++){ const a0=y*PW, a1=a0+PW, a2=a1+PW, o=y*W;
    for(let x=0;x<PW;x++){ const a=d[a0+x], b=d[a1+x], c=d[a2+x]; const l1=a<b?a:b, h1=a<b?b:a; lo[x]=l1<c?l1:c; const h2=l1<c?c:l1; mi[x]=h1<h2?h1:h2; hi[x]=h1<h2?h2:h1; }
    for(let x=0;x<W;x++){
      let L=lo[x]>lo[x+1]?lo[x]:lo[x+1]; L=L>lo[x+2]?L:lo[x+2];
      let U=hi[x]<hi[x+1]?hi[x]:hi[x+1]; U=U<hi[x+2]?U:hi[x+2];
      const q=mi[x], r=mi[x+1], s=mi[x+2], qr=q<r?q:r, QR=q<r?r:q, t=QR<s?QR:s, M=qr>t?qr:t;
      const a=L<M?L:M, A=L<M?M:L, u=A<U?A:U; out[o+x]=a>u?a:u; }
    if((y&127)===0) IC.progress(y/H); }
  let anyNaN=false; for(let i=0;i<d.length;i++) if(d[i]!==d[i]){ anyNaN=true; break; }
  if(anyNaN){ const L=IC.lineRadii(1), B=medBufs(L.n), mark=new Uint8Array(W*H);
    for(let py=0;py<P.PH;py++) for(let px=0;px<PW;px++) if(d[py*PW+px]!==d[py*PW+px])
      for(let y=Math.max(0,py-2);y<=Math.min(H-1,py);y++) for(let x=Math.max(0,px-2);x<=Math.min(W-1,px);x++) mark[y*W+x]=1;
    for(let y=0;y<H;y++) for(let x=0;x<W;x++) if(mark[y*W+x]) out[y*W+x]=medianAt(P,1,L.dx,x,y,B); }
  return out;
}
/* 중앙값 (모든 형식): 값을 순위로 바꿔 이동 히스토그램(Huang) + 세 층 점유 비트맵.
   순위: 8/16-bit 정수 영상은 값 그대로, 아니면 denseRanks. 비트맵 L0 비트 q = (C[q]>0), L1 비트 w = (L0[w]≠0), L2 비트 = (L1 낱말≠0).
   더하고 뺄 때 분기 없이 비트를 고치고, 중앙값 순위 m 을 옮길 때 빈 칸을 건너뛰어 바로 다음·이전 값으로 간다.
   서로 다른 값이 65536 개보다 많으면(32-bit) 히스토그램은 거친 순위 q>>sh 로 (캐시에 들어가게) 돌리고,
   중앙값이 든 거친 칸 안에서 창의 값을 훑어 정확한 순위를 고른다.
   NaN 없는 창은 개수 n 이 늘 홀수(원형 커널) → 중앙값이 하나로 정해져 화소마다 모아 고르는 medianAt 과 같은 값.
   창에 NaN 이 있으면(짝수 개일 때 결과가 추정값에 따라 다름) 그 화소만 medianAt → ImageJ 와 같음. */
function rankMedian(p,radius){
  const {kR,dx,n}=IC.lineRadii(radius); const P=padded(p,kR); const W=p.w,H=p.h,PW=P.PW,PH=P.PH,d=P.d,NP=PW*PH;
  const out=new Float32Array(W*H), half=n>>1;
  if(kR===1&&dx[0]===1&&dx[2]===1) return median3x3(p,P,out);
  let rk, vals=null, U=0, intv=p.bit===8||p.bit===16;
  if(intv) for(let i=0;i<NP;i++){ const v=d[i]; if(!(v>=0&&v<=65535&&v===Math.floor(v))){ intv=false; break; } if(v>=U) U=v+1; }
  if(intv){ rk=new Int32Array(NP); for(let i=0;i<NP;i++) rk[i]=d[i]; }
  else if(n<25){ const B=medBufs(n); for(let y=0;y<H;y++){ for(let x=0;x<W;x++) out[y*W+x]=medianAt(P,kR,dx,x,y,B); if((y&31)===0) IC.progress(y/H); } return out; } // 32-bit 작은 창: 화소마다 고르는 편(옛 경로)이 순위 매기기보다 빠르다
  else { const R=denseRanks(d); rk=R.rk; U=R.U; vals=R.vals; }
  if(!U){ out.fill(NaN); return out; }
  let anyNaN=false; if(!intv) for(let i=0;i<NP;i++) if(rk[i]<0){ anyNaN=true; break; }
  if(U<=1024&&!anyNaN){ // 값 종류가 적으면(8-bit 등) 비트맵 없이 칸을 하나씩 걸어도 충분
    const C=new Int32Array(U+1), KH=2*kR+1, left=new Int32Array(KH), right=new Int32Array(KH);
    for(let y=0;y<H;y++){
      for(let k=0;k<KH;k++){ const r=dx[k], row=(y+k)*PW+kR; left[k]=row-r-1; right[k]=row+r; for(let xx=-r;xx<=r;xx++) C[rk[row+xx]]++; }
      let m=0, below=0;
      for(let x=0;x<W;x++){
        if(x>0) for(let k=0;k<KH;k++){ const qo=rk[left[k]+x], qi=rk[right[k]+x]; C[qo]--; C[qi]++; below+=((qi-m)>>>31)-((qo-m)>>>31); }
        while(below+C[m]<=half){ below+=C[m]; m++; }
        while(below>half){ m--; below-=C[m]; }
        out[y*W+x]=vals?vals[m]:m; }
      for(let k=0;k<KH;k++){ const r=dx[k], row=(y+k)*PW+kR+W-1; for(let xx=-r;xx<=r;xx++) C[rk[row+xx]]--; }
      if((y&63)===0) IC.progress(y/H); }
    return out;
  }
  const fine=rk; let sh=0; while((U-1)>>sh>=65536) sh++;
  if(sh){ const rc=new Int32Array(NP); for(let i=0;i<NP;i++){ const q=fine[i]; rc[i]=q<0?-1:q>>sh; } rk=rc; U=((U-1)>>sh)+1; }
  const C=new Int32Array(U+1), n0=(U>>5)+2, n1=(n0>>5)+2, n2=(n1>>5)+2, L0=new Int32Array(n0), L1=new Int32Array(n1), L2=new Int32Array(n2);
  const add=q=>{ C[q]++; const w=q>>5, w1=w>>5; L0[w]|=1<<(q&31); L1[w1]|=1<<(w&31); L2[w1>>5]|=1<<(w1&31); };
  const rem=q=>{ const c=--C[q], w=q>>5, w1=w>>5;
    const v0=L0[w]&~(((c-1)>>>31)<<(q&31)); L0[w]=v0;                       // c=0 이면 비트 지움
    const v1=L1[w1]&~((((v0|-v0)>>>31)^1)<<(w&31)); L1[w1]=v1;              // 낱말이 비면 위 비트 지움
    L2[w1>>5]&=~((((v1|-v1)>>>31)^1)<<(w1&31)); };
  const nxt=q=>{ // q 이상에서 가장 작은 점유 순위 (있다고 보장될 때만 부름)
    let w=q>>5, b=L0[w]&(-1<<(q&31));
    if(b===0){ const q1=w+1; let w1=q1>>5, b1=L1[w1]&(-1<<(q1&31));
      if(b1===0){ const q2=w1+1; let w2=q2>>5, b2=L2[w2]&(-1<<(q2&31)); while(b2===0){ w2++; b2=L2[w2]; } w1=(w2<<5)|(31-Math.clz32(b2&-b2)); b1=L1[w1]; }
      w=(w1<<5)|(31-Math.clz32(b1&-b1)); b=L0[w]; }
    return (w<<5)|(31-Math.clz32(b&-b)); };
  const prv=q=>{ // q 이하에서 가장 큰 점유 순위 (있다고 보장될 때만 부름)
    let w=q>>5, b=L0[w]&(-1>>>(31-(q&31)));
    if(b===0){ const q1=w-1; let w1=q1>>5, b1=q1<0?0:L1[w1]&(-1>>>(31-(q1&31)));
      if(b1===0){ const q2=w1-1; let w2=q2>>5, b2=q2<0?0:L2[w2]&(-1>>>(31-(q2&31))); while(b2===0){ w2--; b2=L2[w2]; } w1=(w2<<5)|(31-Math.clz32(b2)); b1=L1[w1]; }
      w=(w1<<5)|(31-Math.clz32(b1)); b=L0[w]; }
    return (w<<5)|(31-Math.clz32(b)); };
  const MB=medBufs(n), KH=2*kR+1, left=new Int32Array(KH), right=new Int32Array(KH);
  const OFF=new Int32Array(n), cand=new Int32Array(n); { let c=0; for(let k=0;k<KH;k++) for(let xx=-dx[k];xx<=dx[k];xx++) OFF[c++]=k*PW+xx; }
  const refine=(base,mc,t)=>{ // 거친 칸 mc 의 값들 중 t 번째(0 부터)로 작은 고운 순위
    let c=0; for(let i=0;i<n;i++){ const q=fine[base+OFF[i]]; if((q>>sh)===mc&&q>=0) cand[c++]=q; }
    if(c>32) return quickSelectI(cand,t,c);
    for(let i=1;i<c;i++){ const v=cand[i]; let j=i-1; while(j>=0&&cand[j]>v){ cand[j+1]=cand[j]; j--; } cand[j+1]=v; }
    return cand[t]; };
  for(let y=0;y<H;y++){
    let nan=0;
    for(let k=0;k<KH;k++){ const r=dx[k], row=(y+k)*PW+kR; left[k]=row-r-1; right[k]=row+r; for(let xx=-r;xx<=r;xx++){ const q=rk[row+xx]; if(q<0) nan++; else add(q); } }
    let m=0, below=0; // below = 순위 < m 인 값의 개수
    for(let x=0;x<W;x++){
      if(x>0){
        if(anyNaN){ for(let k=0;k<KH;k++){ const qo=rk[left[k]+x], qi=rk[right[k]+x];
            if(qo<0) nan--; else { rem(qo); below-=(qo-m)>>>31; }
            if(qi<0) nan++; else { add(qi); below+=(qi-m)>>>31; } }
          if(nan){ out[y*W+x]=medianAt(P,kR,dx,x,y,MB); continue; } }
        else for(let k=0;k<KH;k++){ const qo=rk[left[k]+x], qi=rk[right[k]+x]; rem(qo); add(qi); below+=((qi-m)>>>31)-((qo-m)>>>31); }
      } else if(nan){ out[y*W]=medianAt(P,kR,dx,0,y,MB); continue; }
      while(below+C[m]<=half){ below+=C[m]; m=nxt(m+1); }
      while(below>half){ m=prv(m-1); below-=C[m]; }
      out[y*W+x]=sh?vals[refine(y*PW+kR+x,m,half-below)]:vals?vals[m]:m;
    }
    // 줄 끝: 마지막 창을 빼서 개수·비트맵을 0 으로
    for(let k=0;k<KH;k++){ const r=dx[k], row=(y+k)*PW+kR+W-1; for(let xx=-r;xx<=r;xx++){ const q=rk[row+xx]; if(q>=0) rem(q); } }
    if((y&31)===0) IC.progress(y/H);
  }
  return out;
}
const medBufs=n=>({buf:new Float64Array(n),ab:new Float64Array(n),bl:new Float64Array(n)});
/* 16/32-bit 한 화소의 중앙값 — ImageJ getNaNAwareMedian 과 같은 규칙.
   추정값 g = 가운데 화소(NaN 이면 결과 NaN). NaN 은 뺀 개수 m, half=⌊m/2⌋.
   g 보다 큰 값이 half 개보다 많으면 그중 (nAbove−half−1) 번째, 작은 값이 half 개보다 많으면 그중 half 번째, 아니면 g */
function medianAt(P,kR,dx,x,y,B){ const PW=P.PW, d=P.d, buf=B.buf, ab=B.ab, bl=B.bl; let c=0;
  for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; for(let xx=-r;xx<=r;xx++) buf[c++]=d[row+xx]; }
  const g=d[(y+kR)*PW+kR+x]; if(g!==g) return NaN;
  let m=0,na=0,nb=0; for(let i=0;i<c;i++){ const v=buf[i]; if(v!==v) continue; m++; if(v>g) ab[na++]=v; else if(v<g) bl[nb++]=v; }
  const half=m>>1; return !m?NaN:na>half?quickSelect(ab,na-half-1,na):nb>half?quickSelect(bl,half,nb):g; }
IC._medianDirect=function(p,radius){ const {kR,dx}=IC.lineRadii(radius); const P=padded(p,kR), B=medBufs(IC.lineRadii(radius).n), out=new Float32Array(p.w*p.h);
  for(let y=0;y<p.h;y++) for(let x=0;x<p.w;x++) out[y*p.w+x]=medianAt(P,kR,dx,x,y,B); return out; }; // 시험용 기준 (화소마다 모아 고름)
/* Remove Outliers — ImageJ RankFilters.filterLine(OUTLIERS) 의 순서 그대로 (float 연산):
   s = 밝은 쪽 −1 · 어두운 쪽 +1 (보통 LUT). 창 안 v·s 의 최댓값 M(NaN 제외, 최소 −Float.MAX) 에 대해 v·s+thr < M 일 때만 중앙값을 구하고
   v·s+thr < med·s 이면 v → med. 중앙값은 getMedian(NaN 을 빼지 않음, 추정값 = 그 줄에서 마지막으로 구한 중앙값, 줄 처음 = 첫 화소(NaN→0)).
   NaN 없는 창은 개수가 홀수라 추정값과 관계없이 참 중앙값 → 빠른 rankMedian 결과를 쓰고, NaN 이 있는 창만 추정값 규칙을 그대로 따른다. */
function removeOutliers(p,radius,threshold,bright){
  const W=p.w,H=p.h,d=p.data, sg=bright?-1:1, thr=Math.fround(threshold), FMAX=3.4028234663852886e38;
  const med=rankMedian(p,radius); const ext=rankMinMax(p,radius,!bright); // 밝은 쪽: 창 최솟값(→ M = −min), 어두운 쪽: 최댓값
  let anyNaN=false; for(let i=0;i<W*H;i++) if(d[i]!==d[i]){ anyNaN=true; break; }
  let hasNaN=null, P=null, kR=0, dx=null, B=null;
  if(anyNaN){ const ind=new Float32Array(W*H); for(let i=0;i<W*H;i++) ind[i]=d[i]!==d[i]?1:0; hasNaN=rankMinMax({w:W,h:H,bit:32,data:ind},radius,true);
    const L=IC.lineRadii(radius); kR=L.kR; dx=L.dx; P=padded(p,kR); B=medBufs(L.n); }
  const out=new Float32Array(W*H);
  for(let y=0;y<H;y++){ let g=d[y*W]; if(g!==g) g=0;
    for(let x=0;x<W;x++){ const i=y*W+x; let v=d[i]; const M=Math.max(-FMAX,sg*ext[i]); const vs=Math.fround(v*sg+thr);
      if(vs<M){ const m=(hasNaN&&hasNaN[i]>0)?medianGuess(P,kR,dx,x,y,B,g):med[i]; g=m; if(vs<m*sg) v=m; }
      out[i]=v; }
    if((y&63)===0) IC.progress(y/H); }
  return out;
}
/* ImageJ getMedian (NaN 을 빼지 않음 — NaN 은 추정값과 같은 것으로 셈) */
function medianGuess(P,kR,dx,x,y,B,g){ const PW=P.PW, d=P.d, ab=B.ab, bl=B.bl; let na=0,nb=0,n=0;
  for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; for(let xx=-r;xx<=r;xx++){ const v=d[row+xx]; n++; if(v>g) ab[na++]=v; else if(v<g) bl[nb++]=v; } }
  const half=n>>1; return na>half?quickSelect(ab,na-half-1,na):nb>half?quickSelect(bl,half,nb):g; }
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
  }else if(type==="outliers") out=removeOutliers(p,r,+o.threshold||0,(o.which||"bright")==="bright");
  else throw new Error("알 수 없는 rank filter: "+type);
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
   ImageJ 와 같이: 반지름 > 10 이면 이미지를 줄여(최솟값) 공을 굴리고 쌍선형으로 되돌린다.
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
  if(ball.shrink>1&&(Math.ceil(W/ball.shrink)<2||Math.ceil(H/ball.shrink)<2)) throw new Error(`Subtract Background: 이미지(${W}×${H})가 반지름 ${radius} 에 비해 너무 얇습니다 — 반지름을 줄이세요 (ImageJ 도 처리하지 않음)`);
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
