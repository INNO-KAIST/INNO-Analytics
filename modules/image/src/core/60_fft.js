/* ============================================================
   Process ▸ FFT — 2D FFT, Bandpass Filter (ij/plugin/filter/FFTFilter.java 정의), Power spectrum,
   주파수 영역 마스크 필터(Inverse FFT)
   ImageJ 는 Hartley 변환(FHT)을 쓰지만 필터가 실수·짝대칭이라 복소 FFT 로 같은 결과를 낸다.
   속도: 크기마다 비트 뒤집기·회전 인자 표를 한 번 만들어 두고, 줄 방향은 연속 메모리에서 바로,
   열 방향은 16 열씩 묶어 모아 변환한다. 입력이 실수인 Bandpass·Power spectrum·Inverse FFT 는
   두 줄을 한 복소 줄로 묶어 변환하고 반쪽 스펙트럼((S/2+1)×S, 열 순서로 저장)만 다뤄 계산·메모리를 절반으로.
   (옛 구현과 결과 차 ~1e-12 상대 — test/test_perf.js 가 확인)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;

/* 크기 n (2 의 거듭제곱) 의 표: 비트 뒤집기 짝, 단계별 회전 인자 w[h+k] = e^{iπk/h} (h=1,2,4…, k<h) */
const TABS=new Map();
function tables(n){
  let t=TABS.get(n); if(t) return t;
  let bits=0; while((1<<bits)<n) bits++;
  const pi=[], pj=[]; for(let i=0;i<n;i++){ let r=0; for(let b=0,x=i;b<bits;b++,x>>=1) r=(r<<1)|(x&1); if(i<r){ pi.push(i); pj.push(r); } }
  const wr=new Float64Array(Math.max(2,n)), wi=new Float64Array(Math.max(2,n));
  for(let h=1;h<n;h<<=1) for(let k=0;k<h;k++){ const a=Math.PI*k/h; wr[h+k]=Math.cos(a); wi[h+k]=Math.sin(a); }
  t={n,si:Int32Array.from(pi),sj:Int32Array.from(pj),wr,wi}; TABS.set(n,t); return t;
}
/* re/im[o .. o+n) 제자리 FFT (크기 조정 없음), inv=true 면 e^{+i…}. 비트 뒤집기 후 radix-2 두 단계씩 묶어(radix-2²) 한 번에 */
function fftCore(re,im,o,n,t,inv){
  const si=t.si, sj=t.sj, ns=si.length;
  for(let k=0;k<ns;k++){ const a=o+si[k], b=o+sj[k]; let v=re[a]; re[a]=re[b]; re[b]=v; v=im[a]; im[a]=im[b]; im[b]=v; }
  const e=o+n, WR=t.wr, WI=t.wi, sg=inv?1:-1;
  let h=1;
  for(;4*h<=n;h<<=2){ const len=4*h, h2=2*h;
    for(let i=o;i<e;i+=len){
      for(let k=0;k<h;k++){
        const tr=WR[h+k], ti=sg*WI[h+k], ur=WR[h2+k], ui=sg*WI[h2+k];
        const a0=i+k, a1=a0+h, a2=a1+h, a3=a2+h;
        const r0=re[a0], i0=im[a0], r1=re[a1], i1=im[a1], r2=re[a2], i2=im[a2], r3=re[a3], i3=im[a3];
        let xr=r1*tr-i1*ti, xi=r1*ti+i1*tr; const b0r=r0+xr, b0i=i0+xi, b1r=r0-xr, b1i=i0-xi;
        xr=r3*tr-i3*ti; xi=r3*ti+i3*tr; const b2r=r2+xr, b2i=i2+xi, b3r=r2-xr, b3i=i2-xi;
        xr=b2r*ur-b2i*ui; xi=b2r*ui+b2i*ur; re[a0]=b0r+xr; im[a0]=b0i+xi; re[a2]=b0r-xr; im[a2]=b0i-xi;
        const pr=b3r*ur-b3i*ui, pq=b3r*ui+b3i*ur; xr=-sg*pq; xi=sg*pr; // b3·u·(±i)
        re[a1]=b1r+xr; im[a1]=b1i+xi; re[a3]=b1r-xr; im[a3]=b1i-xi; } } }
  if(h<n){ // 남은 radix-2 한 단계
    for(let i=o;i<e;i+=2*h) for(let k=0;k<h;k++){ const wr=WR[h+k], wi=sg*WI[h+k]; const a=i+k, b=a+h;
      const br=re[b], bi=im[b], xr=br*wr-bi*wi, xi=br*wi+bi*wr; re[b]=re[a]-xr; im[b]=im[a]-xi; re[a]+=xr; im[a]+=xi; } }
}
/* 1D FFT (제자리), re/im Float64Array, inv=true 면 역변환(1/N 포함) — 옛 API 그대로 */
function fft1(re,im,inv){
  const n=re.length; fftCore(re,im,0,n,tables(n),inv);
  if(inv){ const s=1/n; for(let i=0;i<n;i++){ re[i]*=s; im[i]*=s; } }
}
IC.fft1=fft1;
const BLK=16;
/* 2D FFT (W×H 둘 다 2 의 거듭제곱), 제자리. 역변환은 1/(W·H) 포함. prog=[p0,p1] 이면 진행 막대 */
IC.fft2=function(re,im,W,H,inv,prog){
  const tw=tables(W), th=tables(H);
  const P=prog?(f)=>IC.progress(prog[0]+(prog[1]-prog[0])*f):null;
  for(let y=0;y<H;y++){ fftCore(re,im,y*W,W,tw,inv); if(P&&(y&127)===0) P(0.5*y/H); }
  const br=new Float64Array(BLK*H), bi=new Float64Array(BLK*H);
  for(let x0=0;x0<W;x0+=BLK){ const nb=Math.min(BLK,W-x0);
    for(let y=0;y<H;y++){ const o=y*W+x0; for(let b=0;b<nb;b++){ br[b*H+y]=re[o+b]; bi[b*H+y]=im[o+b]; } }
    for(let b=0;b<nb;b++) fftCore(br,bi,b*H,H,th,inv);
    for(let y=0;y<H;y++){ const o=y*W+x0; for(let b=0;b<nb;b++){ re[o+b]=br[b*H+y]; im[o+b]=bi[b*H+y]; } }
    if(P&&(x0&127)===0) P(0.5+0.5*x0/W); }
  if(inv){ const s=1/(W*H); for(let i=0;i<W*H;i++){ re[i]*=s; im[i]*=s; } }
};
IC.nextPow2=function(n){ let i=1; while(i<n) i*=2; return i; };

/* ---------- 실수 입력 2D FFT, 반쪽 스펙트럼 ----------
   결과 T={re,im}: 길이 (S/2+1)·S, T[c·S+ky] = F(kx=c, ky)  (c=0..S/2). 나머지 반은 F(S−c,S−ky)=conj F(c,ky).
   rowFn(y,buf) 이 실수 줄 y (길이 S) 를 buf 에 채운다. */
function rfft2(rowFn,S,prog){
  const t=tables(S), C=(S>>1)+1, Tr=new Float64Array(C*S), Ti=new Float64Array(C*S);
  const zr=new Float64Array(S), zi=new Float64Array(S), Br=new Float64Array(C*(BLK+1)), Bi=new Float64Array(C*(BLK+1));
  for(let y0=0;y0<S;y0+=BLK){ const nb=Math.min(BLK,S-y0);
    for(let q=0;q<nb;q+=2){ rowFn(y0+q,zr); if(q+1<nb) rowFn(y0+q+1,zi); else zi.fill(0);
      fftCore(zr,zi,0,S,t,false);
      const oa=q*C, ob=oa+C;
      for(let c=0;c<C;c++){ const m=c===0?0:S-c; const ar=zr[c], ai=zi[c], mr=zr[m], mi=zi[m];
        Br[oa+c]=0.5*(ar+mr); Bi[oa+c]=0.5*(ai-mi);          // 줄 a 의 F
        Br[ob+c]=0.5*(ai+mi); Bi[ob+c]=-0.5*(ar-mr); } }     // 줄 b 의 F
    for(let c=0;c<C;c++){ const o=c*S+y0; for(let q=0,oq=c;q<nb;q++,oq+=C){ Tr[o+q]=Br[oq]; Ti[o+q]=Bi[oq]; } }
    if(prog&&(y0&127)===0) prog(0.5*y0/S); }
  for(let c=0;c<C;c++){ fftCore(Tr,Ti,c*S,S,t,false); if(prog&&(c&63)===0) prog(0.5+0.5*c/C); }
  return {re:Tr,im:Ti,S,C};
}
/* rfft2 의 역: T(에르미트 대칭 스펙트럼의 반쪽)를 제자리에서 열 역변환한 뒤, 줄 ya..yb−1 만 실수로 되돌려
   outFn(y,buf) 로 넘긴다 (buf 길이 S, 1/S² 포함). T 는 바뀐다. */
function irfft2(T,ya,yb,outFn,prog){
  const S=T.S, C=T.C, Tr=T.re, Ti=T.im, t=tables(S);
  for(let c=0;c<C;c++){ fftCore(Tr,Ti,c*S,S,t,true); if(prog&&(c&63)===0) prog(0.5*c/C); }
  const zr=new Float64Array(S), zi=new Float64Array(S), Br=new Float64Array(C*(BLK+1)), Bi=new Float64Array(C*(BLK+1)), sc=1/(S*S);
  const n=yb-ya;
  for(let y0=ya;y0<yb;y0+=BLK){ const nb=Math.min(BLK,yb-y0);
    if(nb&1){ Br.fill(0,nb*C,(nb+1)*C); Bi.fill(0,nb*C,(nb+1)*C); } // 홀수 줄이면 짝 없는 마지막 줄의 짝 = 0
    for(let c=0;c<C;c++){ const o=c*S+y0; for(let q=0,oq=c;q<nb;q++,oq+=C){ Br[oq]=Tr[o+q]; Bi[oq]=Ti[o+q]; } }
    for(let q=0;q<nb;q+=2){ const two=q+1<nb, oa=q*C, ob=oa+C;
      // Z[c] = A[c] + i·B[c],  Z[S−c] = conj A[c] + i·conj B[c]
      for(let c=0;c<C;c++){ const ar=Br[oa+c], ai=Bi[oa+c], bR=Br[ob+c], bI=Bi[ob+c];
        zr[c]=ar-bI; zi[c]=ai+bR;
        if(c>0&&c<S-c){ zr[S-c]=ar+bI; zi[S-c]=bR-ai; } }
      fftCore(zr,zi,0,S,t,true);
      for(let x=0;x<S;x++){ zr[x]*=sc; zi[x]*=sc; }
      outFn(y0+q,zr); if(two) outFn(y0+q+1,zi); }
    if(prog&&((y0-ya)&127)===0) prog(0.5+0.5*(y0-ya)/Math.max(1,n)); }
}
IC._rfft2=rfft2; IC._irfft2=irfft2;

/* FFTFilter.tileMirror 의 줄 생성기: 크기 S×S 에 (x0,y0) 위치로 넣고 사방을 거울상으로 채운 영상의 줄 y */
function mirrorRows(p,S,x0,y0){
  const W=p.w,H=p.h,d=p.data, xm=new Int32Array(S);
  for(let x=0;x<S;x++){ let sx=x-x0; const kx=Math.floor(sx/W); sx=sx-kx*W; if(kx&1) sx=W-1-sx; xm[x]=sx; }
  return (y,buf)=>{ let sy=y-y0; const ky=Math.floor(sy/H); sy=sy-ky*H; if(ky&1) sy=H-1-sy; const o=sy*W; for(let x=0;x<S;x++) buf[x]=d[o+xm[x]]; };
}
/* 가운데 (x0,y0) 에 넣고 바깥은 평균 m, 값에서 m 을 빼고(sub) 창 hx·hy 를 곱한 줄 */
function meanPadRows(p,S,x0,y0,m,sub,hx,hy){
  const W=p.w,H=p.h,d=p.data, base=sub?m:0;
  return (y,buf)=>{ const sy=y-y0, inY=sy>=0&&sy<H, wy=hy?hy[y]:1;
    for(let x=0;x<S;x++){ let v=m; const sx=x-x0; if(inY&&sx>=0&&sx<W) v=d[sy*W+sx]; buf[x]=hx?(v-base)*(hx[x]*wy):v-base; } };
}
/* Bandpass Filter. o: {large(40), small(3), stripes:"none"|"horizontal"|"vertical", tolerance(5 %), autoscale(true), saturate(true)} */
IC.bandpass=function(p,o){
  o=Object.assign({large:40,small:3,stripes:"none",tolerance:5,autoscale:true,saturate:true},o||{});
  const W=p.w,H=p.h,maxN=Math.max(W,H); let S=2; while(S<1.5*maxN) S*=2;
  const fL=2*o.large/S, fS=2*o.small/S, sL=fL*fL, sS=fS*fS; let sStr=(100-o.tolerance)/100; sStr*=sStr;
  const x0=Math.round((S-W)/2), y0=Math.round((S-H)/2);
  IC.progress(0);
  const T=rfft2(mirrorRows(p,S,x0,y0),S,f=>IC.progress(0.45*f));
  // 필터 = 분리형 인자의 곱 (옛 구현과 같은 float 연산 순서 → 같은 값). 반쪽 스펙트럼 T[c·S+ky], c = kx ≤ S/2
  const half=S/2, mode=o.stripes==="horizontal"?1:o.stripes==="vertical"?2:0, Tr=T.re, Ti=T.im;
  const EL=new Float64Array(half+1), ES=new Float64Array(half+1), EZ=new Float64Array(half+1);
  for(let j=0;j<=half;j++){ EL[j]=Math.fround(Math.exp(-(j*j)*sL)); ES[j]=Math.fround(Math.exp(-(j*j)*sS)); EZ[j]=1-Math.fround(Math.exp(-(j*j)*sStr)); }
  for(let c=0;c<T.C;c++){ const elc=EL[c], esc=ES[c], ezc=EZ[c], o2=c*S;
    for(let ky=0;ky<S;ky++){ if(c===0&&ky===0) continue; const j=ky<=half?ky:S-ky;
      const fl=EL[j]*elc, fs=ES[j]*esc; let f=(1-fl)*fs; if(mode===1) f*=ezc; else if(mode===2) f*=EZ[j];
      Tr[o2+ky]*=f; Ti[o2+ky]*=f; } }
  IC.progress(0.5);
  const out=new Float32Array(W*H);
  irfft2(T,y0,y0+H,(y,buf)=>{ const o3=(y-y0)*W; for(let x=0;x<W;x++) out[o3+x]=buf[x+x0]; },f=>IC.progress(0.5+0.45*f));
  IC.progress(0.97);
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
  o=o||{}; const W=p.w,H=p.h; const S=IC.nextPow2(Math.max(W,H));
  let mean=0; for(let i=0;i<W*H;i++) mean+=p.data[i]; mean/=W*H;
  const x0=(S-W)>>1, y0=(S-H)>>1;
  let hx=null; if(o.window==="hann"){ hx=new Float64Array(S); for(let x=0;x<S;x++) hx[x]=0.5-0.5*Math.cos(2*Math.PI*(x+0.5)/S); }
  IC.progress(0);
  const T=rfft2(meanPadRows(p,S,x0,y0,mean,true,hx,hx),S,f=>IC.progress(0.85*f)), Tr=T.re, Ti=T.im, h2=S>>1;
  // 반쪽 스펙트럼 T[c·S+ky] 를 16 열씩 묶어 읽고 (메모리 순서), 가운데 정렬 위치와 그 점대칭 위치(|F(−k)|=|F(k)|)에 쓴다
  const amp=new Float32Array(S*S); let mn=Infinity,mx=-Infinity;
  for(let c0=0;c0<T.C;c0+=BLK){ const c1=Math.min(T.C,c0+BLK);
    for(let ky=0;ky<S;ky++){ const yo=((ky+h2)%S)*S, ym=((S-ky)%S+h2)%S*S;
      for(let c=c0;c<c1;c++){ const k=c*S+ky, re=Tr[k], im=Ti[k], a=Math.log(1+Math.sqrt(re*re+im*im));
        amp[yo+(c+h2)%S]=a; if(c>0&&c<h2) amp[ym+(S-c+h2)%S]=a; if(a<mn) mn=a; if(a>mx) mx=a; } } }
  const out=new Float32Array(S*S); for(let i=0;i<out.length;i++) out[i]=mx>mn?Math.round(255*(amp[i]-mn)/(mx-mn)):0;
  IC.progress(1);
  return {w:S,h:S,bit:8,data:out,pad:{x0,y0,S}};
};
/* 주파수 마스크 필터: maskF (S×S, 가운데 정렬, 1=통과) 를 적용해 역변환. keep=true 면 마스크 안만 통과, false 면 마스크 안 제거.
   결과 = Re(IFFT(F·m)). 실수 입력이라 Re(IFFT(F·m)) = IFFT(F·(m(k)+m(−k))/2) → 반쪽 스펙트럼에서 대칭화한 마스크를 곱한다 */
IC.fftMaskFilter=function(p,maskF,S,keep){
  const W=p.w,H=p.h; const x0=(S-W)>>1, y0=(S-H)>>1; let mean=0; for(let i=0;i<W*H;i++) mean+=p.data[i]; mean/=W*H;
  IC.progress(0);
  const T=rfft2(meanPadRows(p,S,x0,y0,mean,false,null,null),S,f=>IC.progress(0.45*f)), Tr=T.re, Ti=T.im, h2=S>>1;
  const pass=(x,y)=>{ if(x===0&&y===0) return 1; const m=maskF[((y+h2)%S)*S+(x+h2)%S]; return (keep?!!m:!m)?1:0; };
  for(let c0=0;c0<T.C;c0+=BLK){ const c1=Math.min(T.C,c0+BLK); // 16 열씩: 마스크를 줄 방향으로 읽게
    for(let ky=0;ky<S;ky++){ const kym=(S-ky)%S;
      for(let c=c0;c<c1;c++){ const f=(pass(c,ky)+pass((S-c)%S,kym))*0.5; if(f!==1){ const k=c*S+ky; Tr[k]*=f; Ti[k]*=f; } } } }
  IC.progress(0.5);
  const out=new Float32Array(W*H);
  irfft2(T,y0,y0+H,(y,buf)=>{ const o3=(y-y0)*W; for(let x=0;x<W;x++) out[o3+x]=buf[x+x0]; },f=>IC.progress(0.5+0.5*f));
  return {w:W,h:H,bit:32,data:out};
};

})();
