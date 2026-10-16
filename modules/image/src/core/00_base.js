/* ============================================================
   INNO Analytics — Image 탭 계산 엔진 (IMG_CORE)
   ImageJ 1.54u (공개 도메인, NIH) 의 알고리즘을 옮겨 같은 결과를 내도록 했다.
   원본 위치는 각 함수 머리말에 적었다. 검증: modules/image/test (ImageJ 실행 결과와 대조).
   화소 자료 표현: plane = {w, h, bit:8|16|32, data:Float32Array}
     - 8/16-bit 는 정수값을 Float32 로 담는다(ImageJ 와 같은 반올림 규칙으로 다시 정수화).
     - RGB 는 UI 층에서 8-bit plane 3 장으로 다룬다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE||(G.IMG_CORE={});

IC.VERSION="1.0.0";
IC.REF="ImageJ 1.54u";

/* ---------- plane 도우미 ---------- */
IC.plane=function(w,h,bit,data){ return {w,h,bit:bit||32,data:data||new Float32Array(w*h)}; };
IC.clonePlane=function(p){ return {w:p.w,h:p.h,bit:p.bit,data:Float32Array.from(p.data)}; };
IC.maxOf=function(bit){ return bit===8?255:bit===16?65535:Infinity; };

/* ImageJ 의 float→정수 쓰기 규칙: (int)(v+0.5) 후 범위로 자르기.
   (ByteProcessor.setPixels(channel, fp), ShortProcessor 동일) — 음수는 0 으로 먼저 자른다 */
IC.roundTo=function(data,bit){
  if(bit===32) return data;
  const mx=bit===8?255:65535;
  for(let i=0;i<data.length;i++){ let v=data[i]+0.5; if(!(v>=0)) v=0; if(v>mx) v=mx; data[i]=Math.floor(v); }
  return data;
};

/* 마스크(Uint8Array, 0/255 또는 0/1) 안 화소만 바꾸고 나머지는 원래 값으로 */
IC.applyMask=function(orig,res,mask){ if(!mask) return res; for(let i=0;i<res.length;i++) if(!mask[i]) res[i]=orig[i]; return res; };

IC.minMax=function(d,mask){ let mn=Infinity,mx=-Infinity; for(let i=0;i<d.length;i++){ if(mask&&!mask[i]) continue; const v=d[i]; if(v!==v||v===Infinity||v===-Infinity) continue; if(v<mn) mn=v; if(v>mx) mx=v; } if(mn===Infinity){ mn=0; mx=0; } return [mn,mx]; }; // ImageJ FloatProcessor.findMinAndMax: NaN·±∞ 제외

/* ---------- 히스토그램 (ImageJ ImageStatistics 규칙) ----------
   8-bit : 0..255 의 256 칸.
   16/32 : [min,max] 를 256 칸으로 (binSize=(max-min)/256, index=(int)((v-min)*256/(max-min)), 마지막 칸으로 자름) */
IC.histogram=function(p,mask,nBins,range){
  const d=p.data; nBins=nBins||256;
  if(p.bit===8&&!range){ const h=new Float64Array(256); let n=0; for(let i=0;i<d.length;i++){ if(mask&&!mask[i]) continue; h[d[i]|0]++; n++; } return {h,min:0,max:256,binSize:1,n}; }
  let [mn,mx]=range||IC.minMax(d,mask); const h=new Float64Array(nBins); let n=0;
  if(!(mx>mn)){ for(let i=0;i<d.length;i++){ if(mask&&!mask[i]) continue; const v=d[i]; if(v===v&&v!==Infinity&&v!==-Infinity){ h[0]++; n++; } } return {h,min:mn,max:mx,binSize:0,n}; }
  const sc=nBins/(mx-mn);
  for(let i=0;i<d.length;i++){ if(mask&&!mask[i]) continue; const v=d[i]; if(v!==v||v===Infinity||v===-Infinity) continue; let k=Math.floor((v-mn)*sc); if(k<0) k=0; if(k>=nBins) k=nBins-1; h[k]++; n++; }
  return {h,min:mn,max:mx,binSize:(mx-mn)/nBins,n};
};

/* ImageJ TypeConverter: 16/32-bit → 8-bit (scale=true), 현재 표시 범위 [min,max] 사용.
   16-bit: v8=(int)((v-min)*256/(max-min+1)+0.5), 32-bit: v8=(int)((v-min)*255/(max-min)+0.5), 0..255 */
IC.to8bitScaled=function(p,min,max){
  const d=p.data, out=new Float32Array(d.length);
  if(p.bit===8){ out.set(d); return out; }
  if(p.bit===16){ const mn=Math.trunc(min), mx=Math.trunc(max); const sc=256/(mx-mn+1);
    for(let i=0;i<d.length;i++){ let v=d[i]-mn; if(v<0) v=0; v=Math.floor(v*sc+0.5); if(v>255) v=255; out[i]=v; } return out; }
  const sc=max>min?255/(max-min):0; // min=max 이면 0 (NaN 이 8-bit 에 들어가지 않게)
  for(let i=0;i<d.length;i++){ let v=d[i]-min; if(!(v>=0)) v=0; v=Math.floor(v*sc+0.5); if(v>255||v!==v) v=v>255?255:0; out[i]=v; }
  return out;
};

/* ---------- 측정 통계 (ImageJ ImageStatistics / ByteStatistics / FloatStatistics) ----------
   mask: Uint8Array(w*h) — null 이면 rect 전체. rect={x,y,w,h} (없으면 영상 전체).
   cal = {pw, ph} (화소 크기). 반환: area, mean, stdDev(n-1), min, max, median, mode, skewness, kurtosis,
   xCentroid, yCentroid (화소 중심 = +0.5), xCM, yCM, rawIntDen, intDen, n */
IC.stats=function(p,mask,rect,cal,local,opt){ // local=true: mask 가 rect 크기(rect 좌표). opt.noMedian: 중앙값·최빈값 생략(빠름)
  const W=p.w, d=p.data; rect=rect||{x:0,y:0,w:p.w,h:p.h}; cal=cal||{pw:1,ph:1};
  let n=0,s1=0,s2=0,s3=0,s4=0,xs=0,ys=0,xsw=0,ysw=0,mn=Infinity,mx=-Infinity;
  const keep=!(opt&&opt.noMedian); const vals=keep?new Float32Array(rect.w*rect.h):null; let nv=0; // 형식 배열 정렬 (JS 배열 + 비교 함수보다 수십 배 빠름)
  for(let y=rect.y;y<rect.y+rect.h;y++){
    for(let x=rect.x;x<rect.x+rect.w;x++){
      const i=y*W+x; if(mask&&!(local?mask[(y-rect.y)*rect.w+(x-rect.x)]:mask[i])) continue; const v=d[i]; if(v!==v||v===Infinity||v===-Infinity) continue; // ImageJ FloatStatistics: ±Float.MAX 안의 값만
      n++; s1+=v; const v2=v*v; s2+=v2; s3+=v*v2; s4+=v2*v2; xs+=x; ys+=y; xsw+=x*v; ysw+=y*v;
      if(v<mn) mn=v; if(v>mx) mx=v; if(keep) vals[nv++]=v;
    }
  }
  const r={n,area:n*cal.pw*cal.ph,min:n?mn:NaN,max:n?mx:NaN};
  if(!n){ Object.assign(r,{mean:NaN,stdDev:NaN,median:NaN,mode:NaN,skewness:NaN,kurtosis:NaN,xCentroid:NaN,yCentroid:NaN,xCM:NaN,yCM:NaN,rawIntDen:0,intDen:0}); return r; }
  const mean=s1/n; r.mean=mean;
  let sd=(n*s2-s1*s1)/n; r.stdDev=(sd>0&&n>1)?Math.sqrt(sd/(n-1)):0;
  const mean2=mean*mean, variance=s2/n-mean2, sdv=Math.sqrt(variance);
  r.skewness=variance>0?((s3-3*mean*s2)/n+2*mean*mean2)/(variance*sdv):NaN;
  r.kurtosis=variance>0?(((s4-4*mean*s3+6*mean2*s2)/n-3*mean2*mean2)/(variance*variance)-3):NaN;
  r.xCentroid=(xs/n+0.5)*cal.pw; r.yCentroid=(ys/n+0.5)*cal.ph;
  r.xCM=s1!==0?(xsw/s1+0.5)*cal.pw:NaN; r.yCM=s1!==0?(ysw/s1+0.5)*cal.ph:NaN;
  r.rawIntDen=s1; r.intDen=r.area*mean;
  if(keep){ const sv=vals.subarray(0,nv).sort(); r.median=IC._median(sv,p.bit); r.mode=IC._mode(sv,p.bit,mn,mx); } else { r.median=NaN; r.mode=NaN; }
  return r;
};
/* 중앙값: 8/16-bit 는 ImageJ 처럼 히스토그램에서 (누적 ≥ n/2 인 첫 값), float 는 정렬 후 가운데(짝수면 두 값 평균) */
IC._median=function(sorted,bit){
  const n=sorted.length; if(!n) return NaN;
  if(bit===8||bit===16){ const half=n/2; let i=0; // 누적 개수가 n/2 를 넘는 첫 값
    let cnt=0; for(i=0;i<n;i++){ cnt++; if(cnt>half) break; } return sorted[Math.min(i,n-1)]; }
  return n%2?sorted[(n-1)>>1]:(sorted[n/2-1]+sorted[n/2])/2;
};
IC._mode=function(sorted,bit,mn,mx){
  const n=sorted.length; if(!n) return NaN;
  if(bit===8||bit===16){ let best=sorted[0],bc=0,cur=sorted[0],c=0; for(let i=0;i<n;i++){ if(sorted[i]===cur) c++; else { if(c>bc){bc=c;best=cur;} cur=sorted[i]; c=1; } } if(c>bc){best=cur;} return best; }
  // float: 256 칸 히스토그램의 최빈 칸 중심 (ImageJ FloatStatistics: histMin+(mode+0.5)*binSize)
  if(!(mx>mn)) return mn; const nb=256, h=new Float64Array(nb), sc=nb/(mx-mn);
  for(const v of sorted){ let k=Math.floor((v-mn)*sc); if(k>=nb) k=nb-1; h[k]++; }
  let m=0; for(let k=1;k<nb;k++) if(h[k]>h[m]) m=k; const bs=(mx-mn)/nb; return mn+m*bs+(bs!==1?bs/2:0); // ImageJ: 칸 폭이 정확히 1 이면 반 칸을 더하지 않음
};

/* 선택 비트 수로 변환 (Image ▸ Type)
   8→16/32: 값 그대로. 16/32→8: 표시 범위로 비례(ImageJ 'Scale when converting' 켬). 32→16: 표시 범위를 0..65535 로.
   scale=false 면 자르기만 한다 */
IC.convertBit=function(p,bit,min,max,scale){
  if(p.bit===bit) return IC.clonePlane(p);
  const d=p.data, n=d.length, out=new Float32Array(n);
  if(bit===32){ out.set(d); return {w:p.w,h:p.h,bit:32,data:out}; }
  if(bit===8){
    if(scale===false){ for(let i=0;i<n;i++){ let v=d[i]; if(!(v>=0)) v=0; if(v>255) v=255; out[i]=p.bit===32?Math.floor(v+0.5)>255?255:Math.floor(v+0.5):Math.floor(v); } return {w:p.w,h:p.h,bit:8,data:out}; }
    return {w:p.w,h:p.h,bit:8,data:IC.to8bitScaled(p,min,max)};
  }
  // → 16
  if(p.bit===8){ out.set(d); return {w:p.w,h:p.h,bit:16,data:out}; }
  if(scale===false){ for(let i=0;i<n;i++){ let v=d[i]; if(!(v>=0)) v=0; if(v>65535) v=65535; out[i]=Math.floor(v+0.5)>65535?65535:Math.floor(v+0.5); } return {w:p.w,h:p.h,bit:16,data:out}; }
  // ImageJ TypeConverter.convertFloatToShort: value=(v-min)*65535/(max-min)+0.5, clamp
  const sc=max>min?65535/(max-min):1;
  for(let i=0;i<n;i++){ let v=(d[i]-min)*sc; if(!(v>=0)) v=0; v=Math.floor(v+0.5); if(v>65535) v=65535; out[i]=v; }
  return {w:p.w,h:p.h,bit:16,data:out};
};

/* RGB → 8-bit 회색. ImageJ 기본은 (R+G+B)/3, 'weighted' 면 0.299/0.587/0.114.
   ColorProcessor.getBrightness/convertToByte: (int)(r*rw+g*gw+b*bw+0.5) */
IC.rgbToGray=function(r,g,b,weighted){
  const n=r.data.length, out=new Float32Array(n);
  const [wr,wg,wb]=weighted?[0.299,0.587,0.114]:[1/3,1/3,1/3];
  for(let i=0;i<n;i++){ out[i]=Math.floor(r.data[i]*wr+g.data[i]*wg+b.data[i]*wb+0.5); if(out[i]>255) out[i]=255; }
  return {w:r.w,h:r.h,bit:8,data:out};
};

/* 진행 알림 (Worker 에서는 postMessage 로 바뀐다) */
IC._prog=null;
IC.progress=function(f){ if(IC._prog) try{ IC._prog(f); }catch(e){} };

})();
