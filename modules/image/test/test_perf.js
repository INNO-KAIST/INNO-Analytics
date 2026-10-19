/* 속도 개선 회귀 시험 — node test/test_perf.js   (modules/image 에서)
   core/60_fft.js (FFT·Bandpass·Inverse FFT) 와 core/10_filters.js (Rank filters) 를 빠르게 바꾼 뒤,
   바꾸기 전 구현(아래 makeOld 에 그대로 얼려 둠 — ImageJ 대조를 통과한 판)과 결과가 같은지 확인하고 시간을 잰다.
     · Rank (min/max/median/outliers/tophat/open/close/despeckle): 8/16-bit 는 비트 단위로 같아야 함, 32-bit 도 같은 값 (NaN·±∞ 포함)
     · FFT 계열: 32-bit 결과의 최대차/범위 < 1e-5, 8/16-bit 출력은 float 반올림 경계 때문에 ±1 이 극히 드물게만 허용
   FULL=1  : 2048² 에서 옛 구현도 돌려 같음 확인 + 옛/새 시간 비교 (오래 걸림, 수 분)
   기본    : 작은·중간 크기(홀수·2 의 거듭제곱 아님 포함)에서 같음 확인 + 2048² 새 구현 시간만 */
"use strict";
const IC=require("./load_core.js");
const FULL=!!process.env.FULL;
let pass=0, fail=0; const fails=[];
function ok(c,name,info){ if(c) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!c) console.log((c?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }

/* ======================================================================
   옛 구현 (2026-10-19 판, 그대로 복사). IC 의 공용 도우미(histogram, roundTo …)만 빌려 쓴다.
   ====================================================================== */
function makeOld(IC0){
const IC={progress(){}, histogram:IC0.histogram, minMax:IC0.minMax, to8bitScaled:IC0.to8bitScaled, maxOf:IC0.maxOf, roundTo:IC0.roundTo};
const f32=Math.fround;
/* ---- 60_fft.js ---- */
function fft1(re,im,inv){
  const n=re.length; for(let i=1,j=0;i<n;i++){ let bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j){ let t=re[i]; re[i]=re[j]; re[j]=t; t=im[i]; im[i]=im[j]; im[j]=t; } }
  for(let len=2;len<=n;len<<=1){ const ang=2*Math.PI/len*(inv?1:-1), wr=Math.cos(ang), wi=Math.sin(ang);
    for(let i=0;i<n;i+=len){ let cr=1,ci=0; for(let k=0;k<len/2;k++){ const a=i+k,b=a+len/2; const xr=re[b]*cr-im[b]*ci, xi=re[b]*ci+im[b]*cr; re[b]=re[a]-xr; im[b]=im[a]-xi; re[a]+=xr; im[a]+=xi; const t=cr*wr-ci*wi; ci=cr*wi+ci*wr; cr=t; } } }
  if(inv) for(let i=0;i<n;i++){ re[i]/=n; im[i]/=n; }
}
IC.fft2=function(re,im,W,H,inv){
  const rr=new Float64Array(W), ri=new Float64Array(W);
  for(let y=0;y<H;y++){ const o=y*W; for(let x=0;x<W;x++){ rr[x]=re[o+x]; ri[x]=im[o+x]; } fft1(rr,ri,inv); for(let x=0;x<W;x++){ re[o+x]=rr[x]; im[o+x]=ri[x]; } }
  const cr=new Float64Array(H), ci=new Float64Array(H);
  for(let x=0;x<W;x++){ for(let y=0;y<H;y++){ cr[y]=re[y*W+x]; ci[y]=im[y*W+x]; } fft1(cr,ci,inv); for(let y=0;y<H;y++){ re[y*W+x]=cr[y]; im[y*W+x]=ci[y]; } }
};
IC.nextPow2=function(n){ let i=1; while(i<n) i*=2; return i; };
function tileMirror(p,S,x0,y0){
  const W=p.w,H=p.h,d=p.data,out=new Float64Array(S*S);
  for(let y=0;y<S;y++){ let sy=y-y0; const ky=Math.floor(sy/H); sy=sy-ky*H; if(ky&1) sy=H-1-sy;
    for(let x=0;x<S;x++){ let sx=x-x0; const kx=Math.floor(sx/W); sx=sx-kx*W; if(kx&1) sx=W-1-sx; out[y*S+x]=d[sy*W+sx]; } }
  return out;
}
IC.bandpass=function(p,o){
  o=Object.assign({large:40,small:3,stripes:"none",tolerance:5,autoscale:true,saturate:true},o||{});
  const W=p.w,H=p.h,maxN=Math.max(W,H); let S=2; while(S<1.5*maxN) S*=2;
  const fL=2*o.large/S, fS=2*o.small/S, sL=fL*fL, sS=fS*fS; let sStr=(100-o.tolerance)/100; sStr*=sStr;
  const x0=Math.round((S-W)/2), y0=Math.round((S-H)/2);
  const re=tileMirror(p,S,x0,y0), im=new Float64Array(S*S);
  IC.fft2(re,im,S,S,false);
  const half=S/2, mode=o.stripes==="horizontal"?1:o.stripes==="vertical"?2:0;
  for(let ky=0;ky<S;ky++){ const j=ky<=half?ky:S-ky; for(let kx=0;kx<S;kx++){ if(kx===0&&ky===0) continue; const c=kx<=half?kx:S-kx;
      const fl=Math.fround(Math.exp(-(j*j)*sL))*Math.fround(Math.exp(-(c*c)*sL)), fs=Math.fround(Math.exp(-(j*j)*sS))*Math.fround(Math.exp(-(c*c)*sS));
      let f=(1-fl)*fs; if(mode===1) f*=(1-Math.fround(Math.exp(-(c*c)*sStr))); else if(mode===2) f*=(1-Math.fround(Math.exp(-(j*j)*sStr)));
      const i=ky*S+kx; re[i]*=f; im[i]*=f; } }
  IC.fft2(re,im,S,S,true);
  const out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=re[(y+y0)*S+x+x0];
  let res={w:W,h:H,bit:32,data:out};
  if(p.bit===32) return res;
  let mn,mx;
  if(o.autoscale){ const h=IC.histogram(res,null); const n=h.n; const thr=o.saturate?Math.trunc(n*1.0/200):0;
    let i=-1,c=0; do{ i++; c+=h.h[i]; }while(!(c>thr)&&i<255); const hmin=i; i=256; c=0; do{ i--; c+=h.h[i]; }while(!(c>thr)&&i>0); const hmax=i;
    if(hmax>hmin){ mn=h.min+hmin*h.binSize; mx=h.min+hmax*h.binSize; } else [mn,mx]=IC.minMax(out);
    if(p.bit===8) return {w:W,h:H,bit:8,data:IC.to8bitScaled(res,mn,mx)};
    const sc=mx-mn===0?1:65535/(mx-mn), d=new Float32Array(W*H); for(let k=0;k<d.length;k++){ let v=(out[k]-mn)*sc; if(v<0) v=0; if(v>65535) v=65535; d[k]=Math.trunc(v+0.5); } return {w:W,h:H,bit:16,data:d};
  }
  const d=new Float32Array(W*H), M=IC.maxOf(p.bit); for(let k=0;k<d.length;k++){ let v=out[k]; if(p.bit===8){ v=v+0.5; if(!(v>=0)) v=0; if(v>255) v=255; d[k]=Math.trunc(v); } else { if(v<0) v=0; if(v>M) v=M; d[k]=Math.trunc(v+0.5); } }
  return {w:W,h:H,bit:p.bit,data:d};
};
IC.powerSpectrum=function(p,o){
  o=o||{}; const W=p.w,H=p.h; const S=IC.nextPow2(Math.max(W,H)); const re=new Float64Array(S*S), im=new Float64Array(S*S);
  let mean=0; for(let i=0;i<W*H;i++) mean+=p.data[i]; mean/=W*H;
  const x0=(S-W)>>1, y0=(S-H)>>1;
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ let v=mean; const sx=x-x0, sy=y-y0; if(sx>=0&&sy>=0&&sx<W&&sy<H) v=p.data[sy*W+sx]; let w=1; if(o.window==="hann"){ w=(0.5-0.5*Math.cos(2*Math.PI*(x+0.5)/S))*(0.5-0.5*Math.cos(2*Math.PI*(y+0.5)/S)); } re[y*S+x]=(v-mean)*w; }
  IC.fft2(re,im,S,S,false);
  const amp=new Float32Array(S*S); let mn=Infinity,mx=-Infinity;
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ const sx=(x+S/2)%S, sy=(y+S/2)%S; const a=Math.log(1+Math.hypot(re[sy*S+sx],im[sy*S+sx])); amp[y*S+x]=a; if(a<mn) mn=a; if(a>mx) mx=a; }
  const out=new Float32Array(S*S); for(let i=0;i<out.length;i++) out[i]=mx>mn?Math.round(255*(amp[i]-mn)/(mx-mn)):0;
  return {w:S,h:S,bit:8,data:out,pad:{x0,y0,S}};
};
IC.fftMaskFilter=function(p,maskF,S,keep){
  const W=p.w,H=p.h; const x0=(S-W)>>1, y0=(S-H)>>1; let mean=0; for(let i=0;i<W*H;i++) mean+=p.data[i]; mean/=W*H;
  const re=new Float64Array(S*S), im=new Float64Array(S*S);
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ let v=mean; const sx=x-x0, sy=y-y0; if(sx>=0&&sy>=0&&sx<W&&sy<H) v=p.data[sy*W+sx]; re[y*S+x]=v; }
  IC.fft2(re,im,S,S,false);
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ const cx=(x+S/2)%S, cy=(y+S/2)%S;
    const mx=(x+S/2)%S, my=(y+S/2)%S; const m=maskF[my*S+mx]; const pass=keep?!!m:!m; if(!pass&&!(x===0&&y===0)){ re[y*S+x]=0; im[y*S+x]=0; } void cx; void cy; }
  IC.fft2(re,im,S,S,true);
  const out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=re[(y+y0)*S+x+x0];
  return {w:W,h:H,bit:32,data:out};
};
/* ---- 10_filters.js (Rank) ---- */
IC.lineRadii=function(radius){
  if(radius>=1.5&&radius<1.75) radius=1.75; else if(radius>=2.5&&radius<2.85) radius=2.85;
  const r2=Math.trunc(radius*radius)+1, kR=Math.trunc(Math.sqrt(r2+1e-10));
  const dx=new Int32Array(2*kR+1); let n=0;
  for(let y=-kR;y<=kR;y++){ const d=y===0?kR:Math.trunc(Math.sqrt(r2-y*y+1e-10)); dx[y+kR]=d; n+=2*d+1; }
  return {kR,dx,n};
};
function padded(p,kR){
  const W=p.w, H=p.h, PW=W+2*kR, PH=H+2*kR, out=new Float32Array(PW*PH), d=p.data;
  for(let y=0;y<PH;y++){ const sy=Math.min(H-1,Math.max(0,y-kR)); const o=y*PW, so=sy*W;
    for(let x=0;x<PW;x++){ const sx=Math.min(W-1,Math.max(0,x-kR)); out[o+x]=d[so+sx]; } }
  return {d:out,PW,PH};
}
function slideMinMax(a,off,n,r,out,isMax){
  const w=2*r+1; if(r===0){ for(let i=0;i<n;i++) out[i]=a[off+i]; return; }
  const g=new Float32Array(n), hh=new Float32Array(n);
  for(let i=0;i<n;i++){ const v=a[off+i]; if(i%w===0) g[i]=v; else g[i]=isMax?Math.max(g[i-1],v):Math.min(g[i-1],v); }
  for(let i=n-1;i>=0;i--){ const v=a[off+i]; if(i===n-1||(i+1)%w===0) hh[i]=v; else hh[i]=isMax?Math.max(hh[i+1],v):Math.min(hh[i+1],v); }
  for(let i=r;i<n-r;i++){ const a1=hh[i-r], b1=g[i+r]; out[i]=isMax?(a1>b1?a1:b1):(a1<b1?a1:b1); }
}
function rankMinMax(p,radius,isMax){
  const {kR,dx}=IC.lineRadii(radius); const P=padded(p,kR);
  { const d=P.d, fill=isMax?-Infinity:Infinity; for(let i=0;i<d.length;i++) if(d[i]!==d[i]) d[i]=fill; } const W=p.w,H=p.h,PW=P.PW;
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
function quickSelect(a,k,n){ let l=0,r=n-1; while(l<r){ const piv=a[(l+r)>>1]; let i=l,j=r; while(i<=j){ while(a[i]<piv) i++; while(a[j]>piv) j--; if(i<=j){ const t=a[i]; a[i]=a[j]; a[j]=t; i++; j--; } } if(j<k) l=i; if(k<i) r=j; } return a[k]; }
function rankMedian(p,radius){
  const {kR,dx,n}=IC.lineRadii(radius); const P=padded(p,kR); const W=p.w,H=p.h,PW=P.PW;
  const out=new Float32Array(W*H), half=n>>1;
  if(p.bit===8){
    const hist=new Int32Array(256);
    for(let y=0;y<H;y++){
      hist.fill(0);
      for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR; for(let xx=-r;xx<=r;xx++) hist[P.d[row+xx]]++; }
      let m=0, below=0;
      const fix=()=>{ while(below+hist[m]<=half){ below+=hist[m]; m++; } while(below>half){ m--; below-=hist[m]; } };
      fix(); out[y*W]=m;
      for(let x=1;x<W;x++){
        for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; const vo=P.d[row-r-1], vi=P.d[row+r];
          hist[vo]--; if(vo<m) below--; hist[vi]++; if(vi<m) below++; }
        fix(); out[y*W+x]=m;
      }
    }
    return out;
  }
  if(n>=25) return rankMedianHist(p,P,kR,dx,n);
  const B=medBufs(n); for(let y=0;y<H;y++){ for(let x=0;x<W;x++) out[y*W+x]=medianAt(P,kR,dx,x,y,B); }
  return out;
}
const medBufs=n=>({buf:new Float64Array(n),ab:new Float64Array(n),bl:new Float64Array(n)});
function medianAt(P,kR,dx,x,y,B){ const PW=P.PW, d=P.d, buf=B.buf, ab=B.ab, bl=B.bl; let c=0;
  for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; for(let xx=-r;xx<=r;xx++) buf[c++]=d[row+xx]; }
  const g=d[(y+kR)*PW+kR+x]; if(g!==g) return NaN;
  let m=0,na=0,nb=0; for(let i=0;i<c;i++){ const v=buf[i]; if(v!==v) continue; m++; if(v>g) ab[na++]=v; else if(v<g) bl[nb++]=v; }
  const half=m>>1; return !m?NaN:na>half?quickSelect(ab,na-half-1,na):nb>half?quickSelect(bl,half,nb):g; }
function rankMedianHist(p,P,kR,dx,n){
  const W=p.w,H=p.h,PW=P.PW,PH=P.PH,d=P.d,NP=PW*PH; const out=new Float32Array(W*H);
  const rk=new Int32Array(NP); let vals=null, U; let int16=p.bit===16;
  if(int16) for(let i=0;i<NP;i++){ const v=d[i]; if(!(v>=0&&v<=65535&&v===Math.floor(v))){ int16=false; break; } }
  if(int16){ U=65536; for(let i=0;i<NP;i++) rk[i]=d[i]; }
  else { const srt=new Float64Array(NP); let m=0; for(let i=0;i<NP;i++) if(d[i]===d[i]) srt[m++]=d[i];
    const so=srt.subarray(0,m).sort(); U=0; for(let k=0;k<m;k++) if(!U||so[k]!==so[U-1]) so[U++]=so[k]; vals=so.subarray(0,U);
    for(let i=0;i<NP;i++){ const v=d[i]; if(v!==v){ rk[i]=-1; continue; } let lo=0,hi=U-1; while(lo<hi){ const mid=(lo+hi)>>1; if(vals[mid]<v) lo=mid+1; else hi=mid; } rk[i]=lo; } }
  if(!U){ out.fill(NaN); return out; }
  const L0=new Int32Array(U+4096), L1=new Int32Array((U>>6)+65), L2=new Int32Array((U>>12)+2); const half=n>>1; const MB=medBufs(n);
  const put=(q,s)=>{ L0[q]+=s; L1[q>>6]+=s; L2[q>>12]+=s; };
  for(let y=0;y<H;y++){
    let nan=0;
    for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR; for(let xx=-r;xx<=r;xx++){ const q=rk[row+xx]; if(q<0) nan++; else put(q,1); } }
    let m=0, below=0;
    for(let x=0;x<W;x++){
      if(x>0) for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; const qo=rk[row-r-1], qi=rk[row+r];
        if(qo<0) nan--; else { put(qo,-1); if(qo<m) below--; }
        if(qi<0) nan++; else { put(qi,1); if(qi<m) below++; } }
      if(nan){ out[y*W+x]=medianAt(P,kR,dx,x,y,MB); continue; }
      while(below+L0[m]<=half){ below+=L0[m]; m++;
        if((m&63)===0) while(below+L1[m>>6]<=half){ below+=L1[m>>6]; m+=64; if((m&4095)===0) while(below+L2[m>>12]<=half){ below+=L2[m>>12]; m+=4096; } } }
      while(below>half){
        if((m&4095)===0&&below-L2[(m>>12)-1]>half){ m-=4096; below-=L2[m>>12]; continue; }
        if((m&63)===0&&below-L1[(m>>6)-1]>half){ m-=64; below-=L1[m>>6]; continue; }
        m--; below-=L0[m]; }
      out[y*W+x]=int16?m:vals[m];
    }
    for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+W-1; for(let xx=-r;xx<=r;xx++){ const q=rk[row+xx]; if(q>=0) put(q,-1); } }
  }
  return out;
}
function removeOutliers(p,radius,threshold,bright){
  const W=p.w,H=p.h,d=p.data, sg=bright?-1:1, thr=Math.fround(threshold), FMAX=3.4028234663852886e38;
  const med=rankMedian(p,radius); const ext=rankMinMax(p,radius,!bright);
  let anyNaN=false; for(let i=0;i<W*H;i++) if(d[i]!==d[i]){ anyNaN=true; break; }
  let hasNaN=null, P=null, kR=0, dx=null, B=null;
  if(anyNaN){ const ind=new Float32Array(W*H); for(let i=0;i<W*H;i++) ind[i]=d[i]!==d[i]?1:0; hasNaN=rankMinMax({w:W,h:H,bit:32,data:ind},radius,true);
    const L=IC.lineRadii(radius); kR=L.kR; dx=L.dx; P=padded(p,kR); B=medBufs(L.n); }
  const out=new Float32Array(W*H);
  for(let y=0;y<H;y++){ let g=d[y*W]; if(g!==g) g=0;
    for(let x=0;x<W;x++){ const i=y*W+x; let v=d[i]; const M=Math.max(-FMAX,sg*ext[i]); const vs=Math.fround(v*sg+thr);
      if(vs<M){ const m=(hasNaN&&hasNaN[i]>0)?medianGuess(P,kR,dx,x,y,B,g):med[i]; g=m; if(vs<m*sg) v=m; }
      out[i]=v; } }
  return out;
}
function medianGuess(P,kR,dx,x,y,B,g){ const PW=P.PW, d=P.d, ab=B.ab, bl=B.bl; let na=0,nb=0,n=0;
  for(let k=-kR;k<=kR;k++){ const r=dx[k+kR], row=(y+kR+k)*PW+kR+x; for(let xx=-r;xx<=r;xx++){ const v=d[row+xx]; n++; if(v>g) ab[na++]=v; else if(v<g) bl[nb++]=v; } }
  const half=n>>1; return na>half?quickSelect(ab,na-half-1,na):nb>half?quickSelect(bl,half,nb):g; }
IC.rank=function(p,type,o){
  o=o||{}; let r=+o.radius; if(type==="despeckle"){ type="median"; r=1; }
  const W=p.w,H=p.h; let out;
  if(type==="median") out=rankMedian(p,r);
  else if(type==="min") out=rankMinMax(p,r,false);
  else if(type==="max") out=rankMinMax(p,r,true);
  else if(type==="open"||type==="close"||type==="tophat"){
    const light=!!o.light; const maxFirst=type==="close"||(type==="tophat"&&light);
    const a=rankMinMax(p,r,maxFirst);
    const b=rankMinMax({w:W,h:H,bit:p.bit,data:a},r,!maxFirst);
    if(type==="tophat"){ out=new Float32Array(W*H); const off=(light&&p.bit!==32)?IC.maxOf(p.bit):0;
      for(let i=0;i<W*H;i++) out[i]=p.data[i]-b[i]+off; }
    else out=b;
  }else if(type==="outliers") out=removeOutliers(p,r,+o.threshold||0,(o.which||"bright")==="bright");
  else throw new Error("old rank: "+type);
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};
return IC;
}
const OLD=makeOld(IC);

/* ====================================================================== 시험 영상 */
let seed=12345; const rnd=()=>{ seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff; };
function mkImg(W,H,bit,kind,opt){
  opt=opt||{}; const d=new Float32Array(W*H), M=bit===8?255:bit===16?65535:1000;
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){ let v;
    if(kind==="random") v=rnd()*M;
    else if(kind==="flat") v=M*0.3;
    else { // structured: 기울기 + 작은 반점들 + 줄무늬 + 잡음
      v=0.2*M*(x/W)+0.1*M*(y/H)+0.15*M*Math.sin(x*0.07)*Math.cos(y*0.05)+0.05*M*rnd();
      const cx=(x*7+y*13)%61, cy=(x*3+y*5)%47; if(cx<3&&cy<3) v+=0.4*M; }
    if(bit!==32){ v=Math.max(0,Math.min(M,Math.round(v))); if(kind==="random"&&bit===8&&rnd()<0.02) v=rnd()<0.5?0:255; }
    else { v=v-200; if(opt.nan&&rnd()<opt.nan) v=NaN; if(opt.inf&&rnd()<opt.inf) v=rnd()<0.5?Infinity:-Infinity; }
    d[y*W+x]=v; }
  if(opt.nanBlock&&bit===32) for(let y=Math.floor(H/3);y<Math.floor(H/3)+6;y++) for(let x=Math.floor(W/2);x<Math.floor(W/2)+9;x++) d[y*W+x]=NaN;
  return {w:W,h:H,bit,data:d};
}
const clone=p=>({w:p.w,h:p.h,bit:p.bit,data:Float32Array.from(p.data)});
/* 비트 단위 같음 (NaN=NaN, +0/-0 는 같은 값으로) */
function sameBits(a,b){ if(a.length!==b.length) return {ok:false,info:`길이 ${a.length}≠${b.length}`}; let bad=0,first=-1;
  for(let i=0;i<a.length;i++){ const u=a[i], v=b[i]; if(u===v||(u!==u&&v!==v)) continue; bad++; if(first<0) first=i; }
  return {ok:bad===0,info:bad?`${bad} 화소 다름 @${first}: new ${a[first]} old ${b[first]}`:"동일"}; }
function relDiff(a,b){ let lo=Infinity,hi=-Infinity,mx=0,nanBad=0;
  for(let i=0;i<b.length;i++){ const v=b[i], u=a[i]; if((u!==u)!==(v!==v)){ nanBad++; continue; } if(v!==v) continue; if(v<lo) lo=v; if(v>hi) hi=v; const d=Math.abs(u-v); if(d>mx) mx=d; }
  return {rel:mx/Math.max(1e-30,hi-lo),nanBad}; }
function intDiff(a,b){ let n=0,mx=0; for(let i=0;i<a.length;i++){ const d=Math.abs(a[i]-b[i]); if(d){ n++; if(d>mx) mx=d; } } return {n,mx}; }
const t=(f)=>{ const t0=process.hrtime.bigint(); const r=f(); return [r,Number(process.hrtime.bigint()-t0)/1e6]; };

/* ====================================================================== 1. Rank filters — 같음 */
console.log("== 1. Rank filters: 새 구현 = 옛 구현 (비트 단위)");
const RANK_CASES=[["min",{}],["max",{}],["median",{}],["outliers",{threshold:null,which:"bright"}],["outliers",{threshold:null,which:"dark"}],["tophat",{}],["tophat",{light:true}],["open",{}],["close",{}],["despeckle",{}]];
const thrOf=bit=>bit===8?20:bit===16?3000:30;
function rankEq(p,radii,label){
  for(const [type,o0] of RANK_CASES) for(const r of (type==="despeckle"?[1]:radii)){
    const o=Object.assign({},o0,{radius:r}); if(o.threshold===null) o.threshold=thrOf(p.bit);
    const a=IC.rank(clone(p),type,o).data, b=OLD.rank(clone(p),type,o).data; const s=sameBits(a,b);
    ok(s.ok,`${label} ${type}${o.which?"/"+o.which:""}${o.light?"/light":""} r=${r}`,s.info); }
}
const RADII=[0.5,1,1.5,2,2.5,3,4,5,7.5,10];
for(const [W,H] of [[37,23],[64,64],[131,97]]) for(const bit of [8,16,32]) for(const kind of ["random","structured"])
  rankEq(mkImg(W,H,bit,kind),RADII,`${W}×${H} ${bit}-bit ${kind}`);
rankEq(mkImg(5,3,8,"random"),[1,2,4,10],"5×3 8-bit (커널>영상)");
rankEq(mkImg(1,40,16,"random"),[1,3],"1×40 16-bit");
rankEq(mkImg(40,1,32,"random"),[1,3],"40×1 32-bit");
rankEq(mkImg(70,50,8,"flat"),[2,5],"70×50 8-bit 평탄");
rankEq(mkImg(97,61,32,"structured",{nan:0.03,inf:0.01}),[1,2,3,5,8],"97×61 32-bit NaN·±∞");
rankEq(mkImg(97,61,32,"random",{nanBlock:true}),[1,2,4],"97×61 32-bit NaN 덩어리");
rankEq(mkImg(80,60,32,"random",{nan:0.6}),[1,3],"80×60 32-bit NaN 60 %");
{ const p=mkImg(90,70,16,"random"); for(let i=0;i<p.data.length;i++) p.data[i]=Math.floor(p.data[i]/4096)*4096+ (i%3); rankEq(p,[2,6,12],"90×70 16-bit 값 몇 개뿐"); }
{ const p=mkImg(90,70,32,"random"); for(let i=0;i<p.data.length;i++) p.data[i]=(i%7)*1e-30*(i%2?1:-1); rankEq(p,[2,5],"90×70 32-bit 작은 값·±0"); }
rankEq(mkImg(400,300,16,"structured"),[2,5,15],"400×300 16-bit");
rankEq(mkImg(400,300,8,"structured"),[3,20],"400×300 8-bit");
rankEq(mkImg(300,200,32,"structured",{nan:0.002}),[2,6],"300×200 32-bit");
// 32-bit, 서로 다른 값 > 65536 개 → 거친 순위 + 칸 안 고르기 경로
{ const p=mkImg(420,330,32,"random"); rankEq(p,[1,2,5],"420×330 32-bit 다른 값 많음");
  const q=clone(p); for(let i=0;i<q.data.length;i++){ if(i%2===0) q.data[i]=5; if(i%97===0) q.data[i]=NaN; } rankEq(q,[3,10],"420×330 32-bit 같은 값 반 + NaN"); }
// 아주 큰 커널에서 쓰는 버퍼 없는 최솟값·최댓값 경로
{ const lim=IC._rankRingLimit; IC._rankRingLimit=0;
  for(const bit of [8,32]){ const p=mkImg(77,59,bit,"structured",bit===32?{nan:0.01}:{}); for(const t of ["min","max","tophat"]) for(const r of [1,4,9]){ const o={radius:r}; const s=sameBits(IC.rank(clone(p),t,o).data,OLD.rank(clone(p),t,o).data); ok(s.ok,`버퍼 없는 경로 ${bit}-bit ${t} r=${r}`,s.info); } }
  IC._rankRingLimit=lim; }

/* ====================================================================== 2. FFT 계열 — 같음 */
console.log("== 2. FFT·Bandpass·Inverse FFT: 새 구현 ≈ 옛 구현");
function fftEq(p,label){
  const BPO=[{},{large:25,small:2,stripes:"horizontal"},{large:30,small:4,stripes:"vertical",autoscale:false,saturate:false},{large:100,small:0,tolerance:20}];
  for(const o of BPO){ const a=IC.bandpass(clone(p),o).data, b=OLD.bandpass(clone(p),o).data; const nm=`${label} bandpass ${JSON.stringify(o)}`;
    if(p.bit===32){ const r=relDiff(a,b); ok(r.rel<1e-5&&!r.nanBad,nm,`최대차/범위 ${r.rel.toExponential(2)}`); }
    else { const r=intDiff(a,b); ok(r.mx<=1&&r.n<=Math.max(2,a.length*1e-4),nm,`다른 화소 ${r.n}, 최대차 ${r.mx}`); } }
  // 32-bit 중간 결과(8/16-bit 입력에서도) 직접 비교
  { const q=Object.assign(clone(p),{bit:32}); const r=relDiff(IC.bandpass(q,{}).data,OLD.bandpass(clone(q),{}).data); ok(r.rel<1e-5,`${label} bandpass float 중간값`,`최대차/범위 ${r.rel.toExponential(2)}`); }
  for(const o of [{},{window:"hann"}]){ const A=IC.powerSpectrum(clone(p),o), B=OLD.powerSpectrum(clone(p),o); const r=intDiff(A.data,B.data);
    ok(A.w===B.w&&r.mx<=1&&r.n<=Math.max(2,A.data.length*1e-4),`${label} power spectrum ${o.window||""}`,`다른 화소 ${r.n}, 최대차 ${r.mx}`); }
  const S=IC.nextPow2(Math.max(p.w,p.h)); const mask=new Uint8Array(S*S); // 비대칭 마스크 (원 하나 + 사각형 하나)
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ const dx=x-S/2-S/9, dy=y-S/2+S/13; if(dx*dx+dy*dy<(S/20)**2||(Math.abs(x-S/2)<S/6&&Math.abs(y-S/2)<4)) mask[y*S+x]=1; }
  for(const keep of [true,false]){ const r=relDiff(IC.fftMaskFilter(clone(p),mask,S,keep).data,OLD.fftMaskFilter(clone(p),mask,S,keep).data); ok(r.rel<1e-5,`${label} inverse FFT mask keep=${keep}`,`최대차/범위 ${r.rel.toExponential(2)}`); }
  { const S2=IC.nextPow2(Math.max(p.w,p.h)); for(const [W2,H2] of [[S2,S2],[S2,Math.max(2,S2/4)],[Math.max(2,S2/2),S2]]){ const re=new Float64Array(W2*H2), im=new Float64Array(W2*H2); for(let i=0;i<re.length;i++){ re[i]=rnd()-0.5; im[i]=rnd()-0.5; }
    for(const inv of [false,true]){ const a=[Float64Array.from(re),Float64Array.from(im)], b=[Float64Array.from(re),Float64Array.from(im)]; IC.fft2(a[0],a[1],W2,H2,inv); OLD.fft2(b[0],b[1],W2,H2,inv);
      const r1=relDiff(a[0],b[0]), r2=relDiff(a[1],b[1]); ok(Math.max(r1.rel,r2.rel)<1e-9,`fft2 ${W2}×${H2} inv=${inv}`,`최대차/범위 ${Math.max(r1.rel,r2.rel).toExponential(2)}`); } } }
}
for(const [W,H] of [[37,23],[100,64],[128,128],[257,190]]) for(const bit of [8,16,32]) fftEq(mkImg(W,H,bit,bit===32?"structured":"random"),`${W}×${H} ${bit}-bit`);
fftEq(mkImg(1000,700,16,"structured"),"1000×700 16-bit");
fftEq(mkImg(1000,700,32,"structured"),"1000×700 32-bit");
fftEq(mkImg(1000,700,8,"structured"),"1000×700 8-bit");
{ const p=mkImg(1,1,32,"random"); ok(IC.bandpass(p,{}).data.length===1,"1×1 bandpass"); }
{ const p=mkImg(3,2,8,"random"); const r=intDiff(IC.bandpass(clone(p),{}).data,OLD.bandpass(clone(p),{}).data); ok(r.mx<=1,"3×2 8-bit bandpass",`다른 ${r.n}`); }
{ let calls=0, last=0, mono=true; IC._prog=f=>{ calls++; if(f<last-1e-9) mono=false; last=f; }; IC.bandpass(mkImg(300,200,8,"random"),{}); IC._prog=null; ok(calls>=10&&mono,"bandpass 진행 막대 호출",`${calls} 회`); }

/* ====================================================================== 3. 시간 (2048²) */
console.log("== 3. 시간 — 2048² (ms)"+(FULL?" · 옛 구현 포함 + 같음 확인":" · 새 구현만 (FULL=1 이면 옛 구현도)"));
const BIG={8:mkImg(2048,2048,8,"structured"),16:mkImg(2048,2048,16,"structured"),32:mkImg(2048,2048,32,"structured")};
const S2048=2048, MASK=new Uint8Array(S2048*S2048); for(let y=0;y<S2048;y++) for(let x=0;x<S2048;x++){ const dx=x-1024-200, dy=y-1024; if(dx*dx+dy*dy<900) MASK[y*S2048+x]=1; }
const OPS=[
  ["Bandpass (기본) 8-bit",   L=>L.bandpass(clone(BIG[8]),{}), "int"],
  ["Bandpass (기본) 32-bit",  L=>L.bandpass(clone(BIG[32]),{}), "f"],
  ["FFT (power spectrum)",    L=>L.powerSpectrum(clone(BIG[16]),{}), "int"],
  ["Inverse FFT (mask)",      L=>L.fftMaskFilter(clone(BIG[32]),MASK,S2048,true), "f"],
  ["Top Hat r5 8-bit",        L=>L.rank(clone(BIG[8]),"tophat",{radius:5}), "eq"],
  ["Median r2 16-bit",        L=>L.rank(clone(BIG[16]),"median",{radius:2}), "eq"],
  ["Median r5 16-bit",        L=>L.rank(clone(BIG[16]),"median",{radius:5}), "eq"],
  ["Median r2 32-bit",        L=>L.rank(clone(BIG[32]),"median",{radius:2}), "eq"],
  ["Median r2 8-bit",         L=>L.rank(clone(BIG[8]),"median",{radius:2}), "eq"],
  ["Remove Outliers r2 16-bit",L=>L.rank(clone(BIG[16]),"outliers",{radius:2,threshold:3000,which:"bright"}), "eq"],
  ["Remove Outliers r2 8-bit",L=>L.rank(clone(BIG[8]),"outliers",{radius:2,threshold:20,which:"bright"}), "eq"],
  ["Minimum r2 16-bit",       L=>L.rank(clone(BIG[16]),"min",{radius:2}), "eq"],
  ["Maximum r2 8-bit",        L=>L.rank(clone(BIG[8]),"max",{radius:2}), "eq"],
  ["Maximum r20 32-bit",      L=>L.rank(clone(BIG[32]),"max",{radius:20}), "eq"],
  ["Despeckle 16-bit",        L=>L.rank(clone(BIG[16]),"despeckle",{}), "eq"],
  ["Despeckle 8-bit",         L=>L.rank(clone(BIG[8]),"despeckle",{}), "eq"],
];
const only=process.env.ONLY?new RegExp(process.env.ONLY,"i"):null;
const rows=[];
for(const [nm,f,cmp] of OPS){ if(only&&!only.test(nm)) continue;
  const [a,tn]=t(()=>f(IC)); let to=null, info="";
  if(FULL){ const [b,tt]=t(()=>f(OLD)); to=tt;
    if(cmp==="eq"){ const s=sameBits(a.data,b.data); ok(s.ok,`2048² ${nm} 같음`,s.info); info=s.info; }
    else if(cmp==="f"){ const r=relDiff(a.data,b.data); ok(r.rel<1e-5,`2048² ${nm} 같음`,`최대차/범위 ${r.rel.toExponential(2)}`); info=`rel ${r.rel.toExponential(1)}`; }
    else { const r=intDiff(a.data,b.data); ok(r.mx<=1&&r.n<=a.data.length*1e-4,`2048² ${nm} 같음`,`다른 화소 ${r.n}`); info=`다른 화소 ${r.n}`; } }
  rows.push([nm,to,tn,info]);
  console.log(`  ${nm.padEnd(28)} ${to==null?"":("옛 "+to.toFixed(0).padStart(7)+"  →  ")}새 ${tn.toFixed(0).padStart(6)} ms${to!=null?`  (×${(to/tn).toFixed(1)})  ${info}`:""}`);
}

console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.slice(0,80).join("\n")); process.exit(1); }
