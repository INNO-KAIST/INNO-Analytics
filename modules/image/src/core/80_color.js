/* ============================================================
   Image ▸ Adjust ▸ Color Threshold — ImageJ 1.54u ij/plugin/frame/ColorThresholder 와 같은 규칙
   · 색 공간: HSB (java.awt.Color.RGBtoHSB 의 float 계산 → (int)(v×255)), RGB, Lab (ColorSpaceConverter, D65:
     L×2.55, a·b ×1.0625+128 반올림, 0..255 로 자름), YUV (Landini: Y=0.299R+0.587G+0.114B 반올림, U=128+⌊0.493(B−Y)+½⌋, V=128+⌊0.877(R−Y)+½⌋, byte 로 감음)
   · 문턱: 채널마다 [min,max] 와 pass(안쪽 선택)·stop(바깥 선택), 세 채널 조건의 AND → 255 (ImageJ apply 의 8 경우와 같음)
   · Auto: HSB 는 Brightness 만, RGB 는 세 채널 모두, Lab·YUV 는 첫 채널만 (autoSetThreshold) — dark 면 [t+1,255], 아니면 [0,t]
   · Sample: ROI 안 화소의 채널별 최소·최대. HSB 의 Hue 는 쓰인 값의 범위와 가장 긴 빈 구간을 비교해 좁은 쪽을 pass 또는 stop 으로
   검증: test/test_color.js — ImageJ 를 실행해 모든 24-bit 색(2²⁴)의 변환 결과·문턱 결과와 대조.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const f=Math.fround;
// java.awt.Color.RGBtoHSB (OpenJDK) — float 연산 순서 그대로
function rgbToHsb(r,g,b,out){
  let cmax=r>g?r:g; if(b>cmax) cmax=b; let cmin=r<g?r:g; if(b<cmin) cmin=b;
  const bri=f(cmax/255); const sat=cmax!==0?f(f(cmax-cmin)/f(cmax)):0; let hue;
  if(sat===0) hue=0;
  else { const d=f(cmax-cmin); const rc=f(f(cmax-r)/d), gc=f(f(cmax-g)/d), bc=f(f(cmax-b)/d);
    if(r===cmax) hue=f(bc-gc); else if(g===cmax) hue=f(f(2+rc)-bc); else hue=f(f(4+gc)-rc);
    hue=f(hue/6); if(hue<0) hue=f(hue+1); }
  out[0]=hue; out[1]=sat; out[2]=bri; return out;
}
// ij.process.ColorSpaceConverter (D65, sRGB) — double
const M=[[0.4124,0.3576,0.1805],[0.2126,0.7152,0.0722],[0.0193,0.1192,0.9505]], WP=[95.0429,100.0,108.8900];
function lin(c){ c=c/255.0; return c<=0.04045?c/12.92:Math.pow((c+0.055)/1.055,2.4); }
const LIN=new Float64Array(256); for(let i=0;i<256;i++) LIN[i]=lin(i)*100.0; /* 같은 double 값을 미리 (정확히 같음) */
function rgbToLab(R,Gg,B,out){
  const r=LIN[R], g=LIN[Gg], b=LIN[B];
  const X=(r*M[0][0])+(g*M[0][1])+(b*M[0][2]), Y=(r*M[1][0])+(g*M[1][1])+(b*M[1][2]), Z=(r*M[2][0])+(g*M[2][1])+(b*M[2][2]);
  let x=X/WP[0], y=Y/WP[1], z=Z/WP[2];
  x=x>0.008856?Math.pow(x,1.0/3.0):(7.787*x)+(16.0/116.0); y=y>0.008856?Math.pow(y,1.0/3.0):(7.787*y)+(16.0/116.0); z=z>0.008856?Math.pow(z,1.0/3.0):(7.787*z)+(16.0/116.0);
  out[0]=(116.0*y)-16.0; out[1]=500.0*(x-y); out[2]=200.0*(y-z); return out;
}
const cl8=v=>v<0?0:v>255?255:v;
IC.COLOR_SPACES=["HSB","RGB","Lab","YUV"];
IC.COLOR_CHANNEL_NAMES={HSB:["Hue","Saturation","Brightness"],RGB:["Red","Green","Blue"],Lab:["L*","a*","b*"],YUV:["Y","U","V"]};
/* RGB 세 장(0..255 값) → 세 채널 Uint8Array */
IC.colorChannels=function(R,Gp,B,space){
  const n=R.length, c0=new Uint8Array(n), c1=new Uint8Array(n), c2=new Uint8Array(n); const t=[0,0,0];
  if(space==="RGB"){ for(let i=0;i<n;i++){ c0[i]=R[i]; c1[i]=Gp[i]; c2[i]=B[i]; } }
  else if(space==="HSB"){ for(let i=0;i<n;i++){ rgbToHsb(R[i]&255,Gp[i]&255,B[i]&255,t); c0[i]=Math.trunc(t[0]*255.0)&255; c1[i]=Math.trunc(t[1]*255.0)&255; c2[i]=Math.trunc(t[2]*255.0)&255; } }
  else if(space==="Lab"){ for(let i=0;i<n;i++){ rgbToLab(R[i]&255,Gp[i]&255,B[i]&255,t);
      const L1=Math.trunc(t[0]*2.55), a1=Math.floor((1.0625*t[1]+128)+0.5), b1=Math.floor((1.0625*t[2]+128)+0.5); c0[i]=cl8(L1); c1[i]=cl8(a1); c2[i]=cl8(b1); } }
  else if(space==="YUV"){ for(let i=0;i<n;i++){ const r=R[i]&255, g=Gp[i]&255, b=B[i]&255; const yf=(0.299*r+0.587*g+0.114*b);
      c0[i]=Math.floor(yf+0.5)&255; c1[i]=(128+Math.floor((0.493*(b-yf))+0.5))&255; c2[i]=(128+Math.floor((0.877*(r-yf))+0.5))&255; } } /* ImageJ getYUV 는 자르지 않고 byte 로 감는다 */
  else throw new Error("색 공간: "+space);
  return [c0,c1,c2];
};
/* 문턱 → 8-bit mask (안 = 255) */
IC.colorMask=function(ch,o){
  const [c0,c1,c2]=ch, n=c0.length, out=new Float32Array(n); const mn=o.min, mx=o.max, ps=o.pass;
  const inR=(v,k)=>v>=mn[k]&&v<=mx[k]; const ok=(v,k)=>ps[k]?inR(v,k):!inR(v,k);
  for(let i=0;i<n;i++) out[i]=ok(c0[i],0)&&ok(c1[i],1)&&ok(c2[i],2)?255:0;
  return out;
};
IC.colorThreshold=function(R,Gp,B,w,h,o){ const ch=IC.colorChannels(R,Gp,B,o.space); return {w,h,bit:8,data:IC.colorMask(ch,o)}; };
/* 채널별 256 칸 히스토그램 (+ 채널 배열 — UI 가 재사용) */
IC.colorHist=function(R,Gp,B,space){ const ch=IC.colorChannels(R,Gp,B,space); const hs=ch.map(c=>{ const h=new Float64Array(256); for(let i=0;i<c.length;i++) h[c[i]]++; return h; }); return {ch,hist:hs}; };
/* Auto (ColorThresholder.autoSetThreshold) — 바꾸는 채널만 돌려줌 */
IC.colorAuto=function(hist,space,method,dark){
  const which=space==="HSB"?[2]:space==="RGB"?[0,1,2]:[0]; const out={};
  for(const k of which){ const t=IC.thresholdLevel(method,Array.from(hist[k])); out[k]=dark?[t+1,255]:[0,t]; }
  return out;
};
/* Sample (ColorThresholder.sample): ROI 안 채널 값들 → 범위, HSB 면 Hue 의 pass/stop 도 */
IC.colorSample=function(ch,maskIdx,space){
  const mn=[256,256,256], mx=[-1,-1,-1]; const bin=new Uint8Array(256);
  for(const i of maskIdx){ for(let k=0;k<3;k++){ const v=ch[k][i]; if(v<mn[k]) mn[k]=v; if(v>mx[k]) mx[k]=v; } bin[ch[0][i]]=1; }
  if(!maskIdx.length) throw new Error("ROI 안에 화소가 없습니다");
  let pass0=true;
  if(space==="HSB"){ let gap=0,maxgap=0,maxgapst=-1,maxgapen=-1,gapst=0; if(bin[0]===0){ gapst=0; gap=1; }
    for(let i=1;i<256;i++){ if(bin[i]===0){ if(bin[i-1]>0){ gap=1; gapst=i; } else gap++; if(gap>maxgap){ maxgap=gap; maxgapst=gapst; maxgapen=i; } } }
    let lo=0,hi=255; for(let i=0;i<256;i++) if(bin[i]>0){ lo=i; break; } for(let i=255;i>=0;i--) if(bin[i]>0){ hi=i; break; }
    if((hi-lo)<maxgap){ pass0=true; mn[0]=lo; mx[0]=hi; } else { pass0=false; mn[0]=maxgapst; mx[0]=maxgapen; } }
  return {min:mn,max:mx,pass0};
};
/* 화면: RGB 세 장(Float32 0..255) → Uint8 (Worker 로 보낼 때 메모리 1/4) */
IC.rgbU8=function(R,Gp,B){ const n=R.length, o=[new Uint8Array(n),new Uint8Array(n),new Uint8Array(n)]; for(let i=0;i<n;i++){ o[0][i]=R[i]; o[1][i]=Gp[i]; o[2][i]=B[i]; } return o; };
IC._rgbToHsb=rgbToHsb; IC._rgbToLab=rgbToLab;
})();
