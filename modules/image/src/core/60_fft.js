/* ============================================================
   Process ▸ FFT — 2D FFT, Bandpass Filter (ij/plugin/filter/FFTFilter.java 정의), Power spectrum,
   주파수 영역 마스크 필터(Inverse FFT)
   ImageJ 는 Hartley 변환(FHT)을 쓰지만 필터가 실수·짝대칭이라 복소 FFT 로 같은 결과를 낸다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;

/* 1D radix-2 FFT (제자리), re/im Float64Array, inv=true 면 역변환(1/N 포함) */
function fft1(re,im,inv){
  const n=re.length; for(let i=1,j=0;i<n;i++){ let bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit; if(i<j){ let t=re[i]; re[i]=re[j]; re[j]=t; t=im[i]; im[i]=im[j]; im[j]=t; } }
  for(let len=2;len<=n;len<<=1){ const ang=2*Math.PI/len*(inv?1:-1), wr=Math.cos(ang), wi=Math.sin(ang);
    for(let i=0;i<n;i+=len){ let cr=1,ci=0; for(let k=0;k<len/2;k++){ const a=i+k,b=a+len/2; const xr=re[b]*cr-im[b]*ci, xi=re[b]*ci+im[b]*cr; re[b]=re[a]-xr; im[b]=im[a]-xi; re[a]+=xr; im[a]+=xi; const t=cr*wr-ci*wi; ci=cr*wi+ci*wr; cr=t; } } }
  if(inv) for(let i=0;i<n;i++){ re[i]/=n; im[i]/=n; }
}
IC.fft1=fft1;
/* 2D FFT (N×M 둘 다 2 의 거듭제곱) */
IC.fft2=function(re,im,W,H,inv){
  const rr=new Float64Array(W), ri=new Float64Array(W);
  for(let y=0;y<H;y++){ const o=y*W; for(let x=0;x<W;x++){ rr[x]=re[o+x]; ri[x]=im[o+x]; } fft1(rr,ri,inv); for(let x=0;x<W;x++){ re[o+x]=rr[x]; im[o+x]=ri[x]; } }
  const cr=new Float64Array(H), ci=new Float64Array(H);
  for(let x=0;x<W;x++){ for(let y=0;y<H;y++){ cr[y]=re[y*W+x]; ci[y]=im[y*W+x]; } fft1(cr,ci,inv); for(let y=0;y<H;y++){ re[y*W+x]=cr[y]; im[y*W+x]=ci[y]; } }
};
IC.nextPow2=function(n){ let i=1; while(i<n) i*=2; return i; };

/* FFTFilter.tileMirror: 크기 S×S 에 (x0,y0) 위치로 넣고 사방을 거울상으로 채운다 */
function tileMirror(p,S,x0,y0){
  const W=p.w,H=p.h,d=p.data,out=new Float64Array(S*S);
  for(let y=0;y<S;y++){ let sy=y-y0; const ky=Math.floor(sy/H); sy=sy-ky*H; if(ky&1) sy=H-1-sy;
    for(let x=0;x<S;x++){ let sx=x-x0; const kx=Math.floor(sx/W); sx=sx-kx*W; if(kx&1) sx=W-1-sx; out[y*S+x]=d[sy*W+sx]; } }
  return out;
}
/* Bandpass Filter. o: {large(40), small(3), stripes:"none"|"horizontal"|"vertical", tolerance(5 %), autoscale(true), saturate(true)} */
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
  // 비례 없이 원래 형식으로 (자르기)
  const d=new Float32Array(W*H), M=IC.maxOf(p.bit); for(let k=0;k<d.length;k++){ let v=out[k]; if(p.bit===8){ v=v+0.5; if(!(v>=0)) v=0; if(v>255) v=255; d[k]=Math.trunc(v); } else { if(v<0) v=0; if(v>M) v=M; d[k]=Math.trunc(v+0.5); } }
  return {w:W,h:H,bit:p.bit,data:d};
};

/* Power spectrum (보기용): 가운데 정렬, ln(1+|F|) 를 0..255 로. 창(window)=Hann 선택 가능 */
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
/* 주파수 마스크 필터: maskF (S×S, 가운데 정렬, 1=통과) 를 적용해 역변환. keep=true 면 마스크 안만 통과, false 면 마스크 안 제거 */
IC.fftMaskFilter=function(p,maskF,S,keep){
  const W=p.w,H=p.h; const x0=(S-W)>>1, y0=(S-H)>>1; let mean=0; for(let i=0;i<W*H;i++) mean+=p.data[i]; mean/=W*H;
  const re=new Float64Array(S*S), im=new Float64Array(S*S);
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ let v=mean; const sx=x-x0, sy=y-y0; if(sx>=0&&sy>=0&&sx<W&&sy<H) v=p.data[sy*W+sx]; re[y*S+x]=v; }
  IC.fft2(re,im,S,S,false);
  for(let y=0;y<S;y++) for(let x=0;x<S;x++){ const cx=(x+S/2)%S, cy=(y+S/2)%S; // (x,y) 원점 정렬 → 가운데 정렬 좌표
    const mx=(x+S/2)%S, my=(y+S/2)%S; const m=maskF[my*S+mx]; const pass=keep?!!m:!m; if(!pass&&!(x===0&&y===0)){ re[y*S+x]=0; im[y*S+x]=0; } void cx; void cy; }
  IC.fft2(re,im,S,S,true);
  const out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=re[(y+y0)*S+x+x0];
  return {w:W,h:H,bit:32,data:out};
};

})();
