/* ============================================================
   Process ▸ Math · Image Calculator · Enhance Contrast · Equalize · Image ▸ Crop/Rotate/Flip/Scale
   Image ▸ Stacks ▸ Z Project · Registration(위상 상관) · Radial Profile · Orientation(구조 텐서)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const f32=Math.fround;

/* ---------- Process ▸ Math (ImageProcessor.process / ShortProcessor.process / FloatProcessor.process) ----------
   ImageJ 정의 그대로. 8-bit: LUT (Add=(int)value 더함, Multiply=round(v·c), Log=round(ln v·255/ln255),
   Exp=round(exp(v·ln255/255)), Gamma=round(exp(ln(v/255)·γ)·255), Sqrt=round(√v)), 0..255 로 자름.
   16-bit: Log/Exp/Gamma 는 화소 min..max(표시 범위) 기준. 32-bit: 실수 연산 그대로 (log 0 = −∞). Divide = Multiply(1/c) */
IC.math=function(p,op,val){
  const d=p.data, n=d.length, out=new Float32Array(n), bit=p.bit;
  if(bit===32){
    for(let i=0;i<n;i++){ const v=d[i]; let r;
      switch(op){ case "add": r=v+val; break; case "sub": r=v-val; break; case "mul": r=v*val; break; case "div": r=v/val; break;
        case "min": r=v<val?val:v; break; case "max": r=v>val?val:v; break; case "set": r=val; break;
        case "gamma": r=v<=0?0:Math.exp(val*Math.log(v)); break; case "log": r=Math.log(v); break; case "exp": r=Math.exp(v); break;
        case "sqr": r=v*v; break; case "sqrt": r=v<=0?0:Math.sqrt(v); break; case "abs": r=Math.abs(v); break; case "reciprocal": r=v===0?NaN:1/v; break; // ImageJ: 1/0 = NaN
        case "nan": r=(v<val[0]||v>val[1])?NaN:v; break;
        default: throw new Error("Math 연산 (32-bit): "+op); }
      out[i]=r; }
    return {w:p.w,h:p.h,bit:32,data:out};
  }
  const jr=x=>Math.floor(x+0.5); // Java Math.round
  if(op==="div"&&val===0) return {w:p.w,h:p.h,bit,data:Float32Array.from(d)}; // ImageJ: 8/16-bit 을 0 으로 나누기는 바꾸지 않음
  if(op==="div"){ op="mul"; val=1/val; }
  if(op==="sub"){ op="add"; val=-val; }
  const M=bit===8?255:65535; let mn=0,mx=M;
  if(bit===16){ [mn,mx]=IC.minMax(d); mn=Math.trunc(mn); mx=Math.trunc(mx); }
  const range=mx-mn, iv=Math.trunc(val);
  const SC=255/Math.log(255);
  for(let i=0;i<n;i++){ const v=d[i]; let r;
    switch(op){
      case "add": r=v+iv; break; case "mul": r=jr(v*val); break; case "set": r=iv; break;
      case "min": r=v<val?iv:v; break; case "max": r=v>val?iv:v; break;
      case "and": r=v&iv; break; case "or": r=v|iv; break; case "xor": r=v^iv; break;
      case "abs": r=v; break; case "reciprocal": throw new Error("Reciprocal 은 32-bit 이미지에서만 (ImageJ 와 같음)");
      case "gamma": r=bit===8?jr(Math.exp(Math.log(v/255)*val)*255):((range<=0||v===mn)?v:jr(Math.exp(val*Math.log((v-mn)/range))*range+mn)); break;
      case "log": r=bit===8?(v===0?0:jr(Math.log(v)*SC)):(v<=0?0:jr(Math.log(v)*(mx/Math.log(mx)))); break;
      case "exp": r=bit===8?jr(Math.exp(v/SC)):jr(Math.exp(v*(Math.log(mx)/mx))); break;
      case "sqr": r=v*v; break; case "sqrt": r=jr(Math.sqrt(v)); break;
      case "nan": throw new Error("NaN Background 는 32-bit 이미지에서만");
      default: throw new Error("Math 연산: "+op); }
    if(!(r>=0)) r=0; if(r>M) r=M; out[i]=r; }
  return {w:p.w,h:p.h,bit,data:out};
};

/* ---------- Image Calculator (ij/plugin/ImageCalculator) ----------
   op: add sub mul div and or xor min max average difference copy transparent-zero
   float=true 면 32-bit 결과 */
IC.calc=function(a,b,op,float){
  if(a.w!==b.w||a.h!==b.h) throw new Error("두 이미지의 크기가 같아야 합니다");
  const n=a.data.length, out=new Float32Array(n), bit=float?32:a.bit;
  for(let i=0;i<n;i++){ const x=a.data[i], y=b.data[i]; let r;
    switch(op){ case "add": r=x+y; break; case "sub": r=x-y; break; case "mul": r=x*y; break; case "div": r=(float||bit===32)?x/y:(y===0?IC.maxOf(bit):Math.trunc(x/y)); break; // 32-bit 은 언제나 실수 나눗셈 (0/0 = NaN, x/0 = ±∞)
      case "and": r=(x&y)>>>0; break; case "or": r=(x|y)>>>0; break; case "xor": r=(x^y)>>>0; break; case "min": r=Math.min(x,y); break; case "max": r=Math.max(x,y); break;
      case "average": r=(x+y)/2; break; case "difference": r=Math.abs(x-y); break; case "copy": r=y; break; case "transparent-zero": r=y===0?x:y; break;
      default: throw new Error("Image Calculator 연산: "+op); }
    out[i]=r; }
  if(bit===32) return {w:a.w,h:a.h,bit:32,data:out};
  // ImageJ 정수 Blitter: 결과를 0..max 로 자르고 나눗셈·평균은 정수 나눗셈
  const M=IC.maxOf(bit); for(let i=0;i<n;i++){ let r=out[i]; if(op==="div"||op==="average") r=Math.trunc(r); else r=Math.round(r); if(!(r>=0)) r=0; if(r>M) r=M; out[i]=r; }
  return {w:a.w,h:a.h,bit,data:out};
};

/* ---------- Enhance Contrast (ContrastEnhancer) ----------
   saturated (%), normalize(8/16: LUT 로 0..max 펼침, 32: 0..1), equalize (sqrt; classic=true 면 표준).
   반환: normalize 면 새 plane, 아니면 {min,max} 표시 범위만 */
IC.saturatedRange=function(p,saturated,mask){
  // 16-bit 은 ImageJ 처럼 전체 65536 칸 히스토그램(ShortStatistics.histogram16)
  const h=p.bit===16?(()=>{ const hh=new Float64Array(65536); let n=0; for(let i=0;i<p.data.length;i++){ if(mask&&!mask[i]) continue; hh[p.data[i]|0]++; n++; } return {h:hh,min:0,binSize:1,n}; })():IC.histogram(p,mask); const n=h.n; const thr=saturated>0?Math.trunc(n*saturated/200):0; const L=h.h.length;
  let i=-1,c=0; do{ i++; c+=h.h[i]; }while(!(c>thr)&&i<L-1); const hmin=i;
  i=L; c=0; do{ i--; c+=h.h[i]; }while(!(c>thr)&&i>0); const hmax=i;
  if(hmax>hmin) return [h.min+hmin*h.binSize,h.min+hmax*h.binSize];
  return null; // ImageJ: hmax ≤ hmin 이면 바꾸지 않음
};
IC.enhanceContrast=function(p,o){
  o=o||{}; const rr=IC.saturatedRange(p,o.saturated!=null?+o.saturated:0.35,o.mask);
  if(!rr) return o.normalize?{plane:{w:p.w,h:p.h,bit:p.bit,data:Float32Array.from(p.data)},unchanged:true}:{range:null,unchanged:true}; // 일정한 이미지·saturated 과다 → 그대로
  const [mn,mx]=rr;
  if(!o.normalize) return {range:[mn,mx]};
  const d=p.data, out=new Float32Array(d.length);
  if(p.bit===32){ const sc=mx>mn?1/(mx-mn):1; for(let i=0;i<d.length;i++){ let v=d[i]-mn; if(v<0) v=0; v*=sc; if(v>1) v=1; out[i]=v; } return {plane:{w:p.w,h:p.h,bit:32,data:out}}; }
  const M=IC.maxOf(p.bit);
  for(let i=0;i<d.length;i++){ const v=d[i]; out[i]=v<=mn?0:v>=mx?M:Math.trunc(((v-mn)/(mx-mn))*M); }
  return {plane:{w:p.w,h:p.h,bit:p.bit,data:out}};
};
IC.equalize=function(p,classic){
  if(p.bit===32) throw new Error("32-bit 이미지는 Equalize 하지 않습니다 (ImageJ 와 같음)");
  const max=p.bit===8?255:65535, range=max; const hist=new Float64Array(max+1); for(const v of p.data) hist[v|0]++;
  const wv=i=>{ const h=hist[i]; return (h<2||classic)?h:Math.sqrt(h); };
  let sum=wv(0); for(let i=1;i<max;i++) sum+=2*wv(i); sum+=wv(max); const scale=range/sum;
  const lut=new Float64Array(max+1); lut[0]=0; sum=wv(0); for(let i=1;i<max;i++){ const dl=wv(i); sum+=dl; lut[i]=Math.round(sum*scale); sum+=dl; } lut[max]=max;
  const out=new Float32Array(p.data.length); for(let i=0;i<out.length;i++) out[i]=lut[p.data[i]|0];
  return {w:p.w,h:p.h,bit:p.bit,data:out};
};

/* ---------- 기하 ---------- */
IC.crop=function(p,r){ const out=new Float32Array(r.w*r.h); for(let y=0;y<r.h;y++) out.set(p.data.subarray((y+r.y)*p.w+r.x,(y+r.y)*p.w+r.x+r.w),y*r.w); return {w:r.w,h:r.h,bit:p.bit,data:out}; };
IC.flip=function(p,dir){ const W=p.w,H=p.h,out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=dir==="h"?p.data[y*W+(W-1-x)]:p.data[(H-1-y)*W+x]; return {w:W,h:H,bit:p.bit,data:out}; };
IC.rotate90=function(p,cw){ const W=p.w,H=p.h,out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const v=p.data[y*W+x]; if(cw) out[x*H+(H-1-y)]=v; else out[(W-1-x)*H+y]=v; } return {w:H,h:W,bit:p.bit,data:out}; };
/* 임의 각도 회전 (ImageProcessor.rotate: 중심 (w−1)/2, (h−1)/2, 쌍선형, 밖 = bg). 양의 각 = 시계방향 (ImageJ) */
IC.rotate=function(p,deg,o){
  o=o||{}; const W0=p.w,H0=p.h; const a=-deg*Math.PI/180, ca=Math.cos(a), sa=Math.sin(a);
  let W=W0,H=H0,ox=0,oy=0; if(o.enlarge){ W=Math.ceil(Math.abs(W0*ca)+Math.abs(H0*sa)); H=Math.ceil(Math.abs(W0*sa)+Math.abs(H0*ca)); ox=(W-W0)/2; oy=(H-H0)/2; }
  const cx=(W0-1)/2, cy=(H0-1)/2, bg=o.bg||0, out=new Float32Array(W*H), d=p.data;
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const dx=x-ox-cx, dy=y-oy-cy; const xs=dx*ca-dy*sa+cx, ys=dx*sa+dy*ca+cy;
    let v=bg; if(xs>=-0.01&&xs<W0&&ys>=-0.01&&ys<H0){ const xc=Math.min(Math.max(xs,0),W0-1.001), yc=Math.min(Math.max(ys,0),H0-1.001); // ImageJ: [−0.01, w) 안이면 w−1.001 로 자름
      const x0=Math.floor(xc), y0=Math.floor(yc), x1=Math.min(W0-1,x0+1), y1=Math.min(H0-1,y0+1), fx=xc-x0, fy=yc-y0;
      v=o.nearest?d[Math.min(H0-1,Math.round(yc))*W0+Math.min(W0-1,Math.round(xc))]:(d[y0*W0+x0]*(1-fx)+d[y0*W0+x1]*fx)*(1-fy)+(d[y1*W0+x0]*(1-fx)+d[y1*W0+x1]*fx)*fy; }
    out[y*W+x]=v; }
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};
/* 크기 바꾸기: 줄일 때 면적 평균(average when downsizing), 키울 때 쌍선형(화소 중심 정렬) */
IC.resize=function(p,W,H,o){
  o=o||{}; const W0=p.w,H0=p.h,d=p.data,out=new Float32Array(W*H); const sx=W0/W, sy=H0/H;
  if(o.method==="nearest"){ for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=d[Math.min(H0-1,Math.floor((y+0.5)*sy))*W0+Math.min(W0-1,Math.floor((x+0.5)*sx))]; }
  else if(sx>=1&&sy>=1&&o.average!==false){ // 면적 평균
    for(let y=0;y<H;y++){ const ya=y*sy, yb=(y+1)*sy; for(let x=0;x<W;x++){ const xa=x*sx, xb=(x+1)*sx; let s=0,wt=0;
      for(let yy=Math.floor(ya);yy<Math.ceil(yb);yy++){ const wy=Math.min(yb,yy+1)-Math.max(ya,yy); if(wy<=0||yy>=H0) continue; for(let xx=Math.floor(xa);xx<Math.ceil(xb);xx++){ const wx=Math.min(xb,xx+1)-Math.max(xa,xx); if(wx<=0||xx>=W0) continue; s+=d[yy*W0+xx]*wx*wy; wt+=wx*wy; } }
      out[y*W+x]=s/wt; } } }
  else { for(let y=0;y<H;y++) for(let x=0;x<W;x++) out[y*W+x]=IC.bilinear(p,(x+0.5)*sx,(y+0.5)*sy); }
  return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)};
};
/* 이미지 이동 (정수·소수): 쌍선형, 밖 = 0 */
IC.translate=function(p,dx,dy){ const W=p.w,H=p.h,out=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const xs=x-dx, ys=y-dy; if(xs<0||ys<0||xs>W-1||ys>H-1){ out[y*W+x]=0; continue; } out[y*W+x]=IC.bilinear(p,xs+0.5,ys+0.5); } return {w:W,h:H,bit:p.bit,data:IC.roundTo(out,p.bit)}; };

/* ---------- Z Project (ij/plugin/ZProjector) ----------
   method: avg | max | min | sum | sd | median. 결과: avg·sum·sd 는 32-bit, max·min·median 은 원래 bit */
IC.zProject=function(planes,method){
  const n=planes.length, N=planes[0].data.length, out=new Float32Array(N);
  if(method==="median"){ const buf=new Float64Array(n); for(let i=0;i<N;i++){ for(let k=0;k<n;k++) buf[k]=planes[k].data[i]; buf.sort(); out[i]=n%2?buf[(n-1)>>1]:(buf[n/2-1]+buf[n/2])/2; } }
  else if(method==="max"||method==="min"){ out.fill(method==="max"?-3.4028234663852886e38:3.4028234663852886e38); for(let k=0;k<n;k++){ const d=planes[k].data; for(let i=0;i<N;i++){ if(method==="max"?d[i]>out[i]:d[i]<out[i]) out[i]=d[i]; } } } // NaN 은 비교에서 빠짐 (ImageJ)
  else { const s=new Float64Array(N), s2=method==="sd"?new Float64Array(N):null, cnt=method==="avg"&&planes[0].bit===32?new Int32Array(N):null;
    for(const p of planes){ const d=p.data; for(let i=0;i<N;i++){ const v=d[i]; if(cnt){ if(v!==v) continue; cnt[i]++; } s[i]+=v; if(s2) s2[i]+=v*v; } }
    for(let i=0;i<N;i++){ if(method==="sum") out[i]=s[i]; else if(method==="avg") out[i]=s[i]/(cnt?cnt[i]:n); else { const v=n>1?(n*s2[i]-s[i]*s[i])/n/(n-1):0; out[i]=v>0?Math.sqrt(v):0; } } } // 32-bit 평균은 NaN 을 뺀 개수로 (ImageJ doAverageFloatProjection)
  // ImageJ 1.54u: avg·max·min 은 원래 형식으로 (int)(v+0.5), sum·sd 는 32-bit, median 은 아래 시험으로 확인
  const keep=method==="avg"||method==="max"||method==="min"||(method==="median"&&planes[0].bit===8); const bit=keep?planes[0].bit:32; // 16-bit median 은 32-bit (ImageJ)
  return {w:planes[0].w,h:planes[0].h,bit,data:keep&&bit!==32?IC.roundTo(out,bit):out};
};

/* ---------- 위상 상관 이동 추정 (translation) ----------
   Hann 창, 교차 전력 스펙트럼 정규화, 정수 극대 + 3 점 포물선 보간(부화소). 반환 {dx,dy,peak} = b 를 a 에 맞추는 이동 */
IC.phaseCorr=function(a,b){
  for(const q of [a,b]){ const d=q.data; for(let i=0;i<d.length;i++){ const v=d[i]; if(v!==v||v===Infinity||v===-Infinity) throw new Error("정렬: 이미지에 NaN·무한대 화소가 있습니다 — Process ▸ Math ▸ NaN Background 등을 되돌리거나 값으로 바꾼 뒤 정렬하세요"); } }
  const W=a.w,H=a.h, S=IC.nextPow2(Math.max(W,H));
  let ma=0,mb=0; for(let i=0;i<W*H;i++){ ma+=a.data[i]; mb+=b.data[i]; } ma/=W*H; mb/=W*H;
  const hann=(t,n)=>t<=-0.5||t>=n-0.5?0:0.5-0.5*Math.cos(2*Math.PI*(t+0.5)/n);
  const fk=k=>k>S/2?k-S:k;
  // 한 번 추정: a 는 고정 창, b 는 창을 (ox,oy) 만큼 옮겨 내용과 함께 움직이게 → 창 때문에 0 쪽으로 끌리는 치우침 제거
  const once=(ox,oy,coarse,c0)=>{
    const ar=new Float64Array(S*S), ai=new Float64Array(S*S), br=new Float64Array(S*S), bi=new Float64Array(S*S);
    for(let y=0;y<H;y++) for(let x=0;x<W;x++){ ar[y*S+x]=(a.data[y*W+x]-ma)*hann(x,W)*hann(y,H); br[y*S+x]=(b.data[y*W+x]-mb)*hann(x+ox,W)*hann(y+oy,H); }
    IC.fft2(ar,ai,S,S,false); IC.fft2(br,bi,S,S,false);
    let mmax=0; const Rr=new Float64Array(S*S), Ri=new Float64Array(S*S);
    for(let i=0;i<S*S;i++){ const r=ar[i]*br[i]+ai[i]*bi[i], im=ai[i]*br[i]-ar[i]*bi[i]; Rr[i]=r; Ri[i]=im; const m=Math.hypot(r,im); if(m>mmax) mmax=m; }
    const eps=1e-6*mmax; for(let i=0;i<S*S;i++){ const m=Math.hypot(Rr[i],Ri[i])+eps; Rr[i]/=m; Ri[i]/=m; }
    let cx,cy,peak;
    if(coarse){ const xr=Float64Array.from(Rr), xi=Float64Array.from(Ri); IC.fft2(xr,xi,S,S,true); let best=-Infinity,bx=0,by=0; for(let i=0;i<S*S;i++) if(xr[i]>best){ best=xr[i]; bx=i%S; by=(i-bx)/S; } cx=bx>S/2?bx-S:bx; cy=by>S/2?by-S:by; peak=best; }
    else { cx=c0[0]; cy=c0[1]; }
    // 국소 DFT 업샘플링 (Guizar-Sicairos 2008, 분리 계산)
    const grid=(gx,gy,step,n)=>{ const off=[]; for(let i=-n;i<=n;i++) off.push(i*step); const m=off.length; const Tr=new Float64Array(S*m), Ti=new Float64Array(S*m);
      for(let ky=0;ky<S;ky++){ const o=ky*S; for(let j=0;j<m;j++){ const x=gx+off[j]; let sr=0,si=0; for(let kx=0;kx<S;kx++){ const ph=2*Math.PI*fk(kx)*x/S, c=Math.cos(ph), sn=Math.sin(ph); sr+=Rr[o+kx]*c-Ri[o+kx]*sn; si+=Rr[o+kx]*sn+Ri[o+kx]*c; } Tr[ky*m+j]=sr; Ti[ky*m+j]=si; } }
      let bv=-Infinity,bxx=gx,byy=gy; for(let i=0;i<m;i++){ const y=gy+off[i]; for(let j=0;j<m;j++){ let sr=0; for(let ky=0;ky<S;ky++){ const ph=2*Math.PI*fk(ky)*y/S; sr+=Tr[ky*m+j]*Math.cos(ph)-Ti[ky*m+j]*Math.sin(ph); } if(sr>bv){ bv=sr; bxx=gx+off[j]; byy=y; } } }
      return [bxx,byy,bv/(S*S)]; };
    let g=grid(cx,cy,1/8,10); g=grid(g[0],g[1],1/128,10);
    return {dx:g[0],dy:g[1],peak:coarse?peak:g[2]};
  };
  let r=once(0,0,true,null);
  for(let it=0;it<3;it++){ const r2=once(r.dx,r.dy,false,[r.dx,r.dy]); const ch=Math.hypot(r2.dx-r.dx,r2.dy-r.dy); r={dx:r2.dx,dy:r2.dy,peak:r.peak}; if(ch<0.004) break; }
  return r;
};
/* 스택 정렬: ref=0(첫 장) 또는 'previous'(이전 장 기준 누적). 반환 각 장의 (dx,dy) */
IC.registerShifts=function(planes,mode){
  const sh=[{dx:0,dy:0,peak:1}]; let acc={dx:0,dy:0};
  for(let k=1;k<planes.length;k++){ const ref=mode==="previous"?planes[k-1]:planes[0]; const s=IC.phaseCorr(ref,planes[k]);
    if(mode==="previous"){ acc={dx:acc.dx+s.dx,dy:acc.dy+s.dy}; sh.push({dx:acc.dx,dy:acc.dy,peak:s.peak}); } else sh.push(s); IC.progress(k/planes.length); }
  return sh;
};

/* ---------- Radial Profile (방위 평균) ----------
   중심 (cx,cy) 화소 좌표(화소 중심 = +0.5), 반지름 rMax, 칸 폭 1 px. 반환 r(칸 중심), mean, n */
IC.radialProfile=function(p,cx,cy,rMax,step){
  step=step||1; const nb=Math.ceil(rMax/step); const s=new Float64Array(nb), c=new Float64Array(nb); const W=p.w,H=p.h;
  const x0=Math.max(0,Math.floor(cx-rMax)), x1=Math.min(W-1,Math.ceil(cx+rMax)), y0=Math.max(0,Math.floor(cy-rMax)), y1=Math.min(H-1,Math.ceil(cy+rMax));
  for(let y=y0;y<=y1;y++) for(let x=x0;x<=x1;x++){ const r=Math.hypot(x+0.5-cx,y+0.5-cy); const k=Math.floor(r/step); if(k<nb){ s[k]+=p.data[y*W+x]; c[k]++; } }
  const r=[],m=[],n=[]; for(let k=0;k<nb;k++){ r.push((k+0.5)*step); m.push(c[k]?s[k]/c[k]:NaN); n.push(c[k]); }
  return {r,mean:m,n};
};

/* ---------- Orientation (구조 텐서; OrientationJ 와 같은 정의) ----------
   기울기: 중앙 차분, 텐서 성분을 Gaussian(sigma) 로 평균. 방향 θ = ½·atan2(2Jxy, Jyy−Jxx) (가로축 기준, 반시계 +, deg),
   coherency = (λ1−λ2)/(λ1+λ2), energy = λ1+λ2. 반환 화소별 θ·coherency 와 가중 히스토그램(−90..90, 1°) */
IC.orientation=function(p,o){
  o=o||{}; const sigma=+o.sigma||2, W=p.w,H=p.h,d=p.data, N=W*H;
  // 기울기: 4차 중앙 차분 (−f₊₂+8f₊₁−8f₋₁+f₋₂)/12 — 2차 중앙 차분은 가는 결에서 방향이 1° 가량 기운다. 가장자리는 복제
  const gx=new Float32Array(N), gy=new Float32Array(N); const at=(x,y)=>d[(y<0?0:y>=H?H-1:y)*W+(x<0?0:x>=W?W-1:x)];
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){ gx[y*W+x]=(-at(x+2,y)+8*at(x+1,y)-8*at(x-1,y)+at(x-2,y))/12; gy[y*W+x]=(-at(x,y+2)+8*at(x,y+1)-8*at(x,y-1)+at(x,y-2))/12; }
  const jxx=new Float32Array(N), jyy=new Float32Array(N), jxy=new Float32Array(N);
  for(let i=0;i<N;i++){ jxx[i]=gx[i]*gx[i]; jyy[i]=gy[i]*gy[i]; jxy[i]=gx[i]*gy[i]; }
  IC.gaussianFloat(jxx,W,H,sigma,sigma,0.0002); IC.gaussianFloat(jyy,W,H,sigma,sigma,0.0002); IC.gaussianFloat(jxy,W,H,sigma,sigma,0.0002);
  const th=new Float32Array(N), coh=new Float32Array(N), en=new Float32Array(N); const hist=new Float64Array(180); let emax=0;
  for(let i=0;i<N;i++){ const a=jxx[i], b=jyy[i], c=jxy[i]; const tr=a+b, df=Math.sqrt((a-b)*(a-b)+4*c*c); const l1=(tr+df)/2, l2=(tr-df)/2;
    // 이미지 y 축이 아래로 향하므로 수학 좌표로 바꿔 반시계 + 로: 결 방향(기울기에 수직)
    // 기울기 방향 φ = ½·atan2(2Jxy, Jxx−Jyy) (이미지 좌표, y 아래). 결 방향 = φ+90°, 화면 반시계 + 로 바꾸면 −(φ+90°)
    let t=-(0.5*Math.atan2(2*c,a-b)*180/Math.PI+90);
    if(t>90) t-=180; if(t<=-90) t+=180;
    th[i]=t; coh[i]=tr>1e-12?(l1-l2)/tr:0; en[i]=tr; if(tr>emax) emax=tr; }
  const minC=o.minCoherency!=null?+o.minCoherency:0;
  const mk=o.mask&&o.mask.length===N?o.mask:null; // mask[i] = 1 → 통계에서 뺌 (축척 막대·글자, ROI 밖)
  if(mk){ emax=0; for(let i=0;i<N;i++) if(!mk[i]&&en[i]>emax) emax=en[i]; }
  const minE2=(o.minEnergy!=null?+o.minEnergy:0.01)*emax;
  for(let i=0;i<N;i++){ if(mk&&mk[i]||en[i]<minE2||coh[i]<minC) continue; let k=Math.floor(th[i]+90); if(k>=180) k=179; if(k<0) k=0; hist[k]+=coh[i]*en[i]; }
  // 주 방향: 선택한 화소의 텐서 합(구간 나눔 없이), order parameter = 이중각 평균 길이
  let A=0,B=0,Cc=0,sc=0,ss=0,sw=0; for(let i=0;i<N;i++){ if(mk&&mk[i]||en[i]<minE2||coh[i]<minC) continue; A+=jxx[i]; B+=jyy[i]; Cc+=jxy[i]; const w=coh[i]*en[i], t=th[i]*Math.PI/180; sc+=w*Math.cos(2*t); ss+=w*Math.sin(2*t); sw+=w; }
  let dom=-(0.5*Math.atan2(2*Cc,A-B)*180/Math.PI+90); if(dom>90) dom-=180; if(dom<=-90) dom+=180; if(!(A+B>0)||!(emax>0)) dom=NaN; /* 일정한 이미지: 방향 없음 */ const R=sw?Math.hypot(sc,ss)/sw:0;
  let ct=0,cw=0; for(let i=0;i<N;i++){ if(mk&&mk[i]||en[i]<minE2) continue; ct+=coh[i]*en[i]; cw+=en[i]; }
  return {theta:th,coherency:coh,energy:en,hist,dominant:dom,orderParam:R,meanCoherency:cw?ct/cw:0,emax};
};

})();
(function(){ const IC=((typeof self!=="undefined")?self:globalThis).IMG_CORE;
/* Image ▸ Adjust ▸ Apply LUT: 표시 범위 [min,max] 를 0..M 로 펼쳐 화소에 씀 (8/16-bit) */
/* ContrastAdjuster.apply: min·max 정수, i≤min → 0, i≥max → range−1, 사이 = (int)((i−min)/(max−min)·range) (range 256 / 65536) */
IC.applyLut=function(p,min,max){ const range=p.bit===8?256:65536, d=p.data, out=new Float32Array(d.length); const a=Math.trunc(min), b=Math.trunc(max);
  for(let i=0;i<d.length;i++){ const v=d[i]; out[i]=v<=a?0:v>=b?range-1:Math.trunc(((v-a)/(b-a))*range); } return {w:p.w,h:p.h,bit:p.bit,data:out}; };
})();
