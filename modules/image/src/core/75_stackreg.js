/* ============================================================
   Plugins ▸ Registration ▸ StackReg / TurboReg — 강도 기반 피라미드 정합 엔진
   pystackreg 0.2.8 의 C++ 소스(src/TurboReg*.cpp, inc/*.h, src/pymain.cpp)와
   pystackreg.py 의 래퍼 규칙을 연산 순서까지 그대로 옮겼다 (double 계산, float 저장 위치 포함).
   검증: modules/image/test/stackreg_parity.py + test_stackreg.js (pystackreg 결과와 대조).

   TurboReg — C++ port by Gregor Lichtner (pystackreg), original code by Philippe Thevenaz, EPFL.
   P. Thevenaz, U.E. Ruttimann, M. Unser, "A Pyramid Approach to Subpixel Registration Based on Intensity",
   IEEE Transactions on Image Processing, vol. 7, no. 1, pp. 27-41, January 1998.
   You are free to use this software for commercial and non-commercial purposes. However, we expect you to include a citation or acknowledgement whenever you present or publish research results that are based on this software. You are free to modify this software or derive works from it, but you are only allowed to distribute it under the same terms as this license specifies. Additionally, you must include a reference to the research paper above in all software and works derived from this software.
   JavaScript port for INNO Analytics (2026).

   좌표 규칙: plane {w,h,data[y*w+x]} 은 numpy 배열 data.reshape(h,w) 와 같다
   (pymain 이 numpy 의 (행,열) 을 C++ 의 (height,width) 로 넘기므로 x=열, y=행).
   행렬(정규형, pystackreg get_matrix): [x_mov, y_mov, 1]^T = M [x_ref, y_ref, 1]^T.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;

const TRANSLATION=2, RIGID_BODY=3, SCALED_ROTATION=4, AFFINE=6, BILINEAR=8;
const PYRAMID_MIN_SIZE=12;
const FIRST_LAMBDA=1.0, LAMBDA_MAGSTEP=4.0, MANY_ITERATIONS=10, PIXEL_HIGH_PRECISION=0.001, ITERATION_PROGRESSION=2;
const GOLDEN_RATIO=0.5*(Math.sqrt(5.0)-1.0);
const f32=Math.fround;          // C++ putRow/putColumn/(float) 캐스트
const trunc=Math.trunc;         // C 의 (int) (0 쪽으로 자름, 범위 밖에서도 값이 감기지 않게)
const MODES={translation:TRANSLATION,rigid:RIGID_BODY,rigid_body:RIGID_BODY,scaled:SCALED_ROTATION,scaled_rotation:SCALED_ROTATION,affine:AFFINE,bilinear:BILINEAR};
function modeCode(mode){ if(typeof mode==="number"&&[2,3,4,6,8].includes(mode)) return mode; const c=MODES[String(mode).toLowerCase()]; if(!c) throw new Error("StackReg: 알 수 없는 변환 "+mode); return c; }

function getPyramidDepth(sw,sh,tw,th){ let d=1; while(2*PYRAMID_MIN_SIZE<=sw&&2*PYRAMID_MIN_SIZE<=sh&&2*PYRAMID_MIN_SIZE<=tw&&2*PYRAMID_MIN_SIZE<=th){ sw=trunc(sw/2); sh=trunc(sh/2); tw=trunc(tw/2); th=trunc(th/2); d++; } return d; }

/* ===================== TurboRegImage (1D/2D B-spline 연산) ===================== */
const Z3=[Math.sqrt(3.0)-2.0];
const Z7=[-0.5352804307964381655424037816816460718339231523426924148812,-0.122554615192326690515272264359357343605486549427295558490763,-0.0091486948096082769285930216516478534156925639545994482648003];
const H3=[2.0/3.0,1.0/6.0], H7=[151.0/315.0,397.0/1680.0,1.0/42.0,1.0/5040.0], H1=[1.0];

function initialCausal(c,n,z){ // getInitialCausalCoefficientMirrorOffBounds (tolerance=0 → horizon=n)
  let z1=z, zn=Math.pow(z,n); let sum=(1.0+z)*(c[0]+zn*c[n-1]); const horizon=n; zn=zn*zn;
  for(let k=1;k<horizon-1;k++){ z1=z1*z; zn=zn/z; sum=sum+(z1+zn)*c[k]; }
  return sum/(1.0-Math.pow(z,2*n));
}
function samplesToInterpolationCoefficient1D(c,n,degree){
  const z=degree===3?Z3:degree===7?Z7:[]; let lambda=1.0;
  if(n===1) return;
  for(let k=0;k<z.length;k++) lambda*=(1.0-z[k])*(1.0-1.0/z[k]);
  for(let i=0;i<n;i++) c[i]=c[i]*lambda;
  for(let k=0;k<z.length;k++){ const zk=z[k];
    c[0]=initialCausal(c,n,zk);
    for(let i=1;i<n;i++) c[i]=c[i]+zk*c[i-1];
    c[n-1]=zk*c[n-1]/(zk-1.0);
    for(let i=n-2;0<=i;i--) c[i]=zk*(c[i+1]-c[i]);
  }
}
function symmetricFir(h,c,s,n){ // symmetricFirMirrorOffBounds1D, c,s 길이 n
  if(h.length===2){ const h0=h[0],h1=h[1];
    if(2<=n){ s[0]=h0*c[0]+h1*(c[0]+c[1]); for(let i=1;i<n-1;i++) s[i]=h0*c[i]+h1*(c[i-1]+c[i+1]); s[n-1]=h0*c[n-1]+h1*(c[n-2]+c[n-1]); }
    else s[0]=(h0+2.0*h1)*c[0];
  } else if(h.length===4){ const h0=h[0],h1=h[1],h2=h[2],h3=h[3];
    if(6<=n){
      s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2])+h3*(c[2]+c[3]);
      s[1]=h0*c[1]+h1*(c[0]+c[2])+h2*(c[0]+c[3])+h3*(c[1]+c[4]);
      s[2]=h0*c[2]+h1*(c[1]+c[3])+h2*(c[0]+c[4])+h3*(c[0]+c[5]);
      for(let i=3;i<n-3;i++) s[i]=h0*c[i]+h1*(c[i-1]+c[i+1])+h2*(c[i-2]+c[i+2])+h3*(c[i-3]+c[i+3]);
      s[n-3]=h0*c[n-3]+h1*(c[n-4]+c[n-2])+h2*(c[n-5]+c[n-1])+h3*(c[n-6]+c[n-1]);
      s[n-2]=h0*c[n-2]+h1*(c[n-3]+c[n-1])+h2*(c[n-4]+c[n-1])+h3*(c[n-5]+c[n-2]);
      s[n-1]=h0*c[n-1]+h1*(c[n-2]+c[n-1])+h2*(c[n-3]+c[n-2])+h3*(c[n-4]+c[n-3]);
    } else switch(n){
      case 5: s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2])+h3*(c[2]+c[3]);
        s[1]=h0*c[1]+h1*(c[0]+c[2])+h2*(c[0]+c[3])+h3*(c[1]+c[4]);
        s[2]=h0*c[2]+h1*(c[1]+c[3])+(h2+h3)*(c[0]+c[4]);
        s[3]=h0*c[3]+h1*(c[2]+c[4])+h2*(c[1]+c[4])+h3*(c[0]+c[3]);
        s[4]=h0*c[4]+h1*(c[3]+c[4])+h2*(c[2]+c[3])+h3*(c[1]+c[2]); break;
      case 4: s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2])+h3*(c[2]+c[3]);
        s[1]=h0*c[1]+h1*(c[0]+c[2])+h2*(c[0]+c[3])+h3*(c[1]+c[3]);
        s[2]=h0*c[2]+h1*(c[1]+c[3])+h2*(c[0]+c[3])+h3*(c[0]+c[2]);
        s[3]=h0*c[3]+h1*(c[2]+c[3])+h2*(c[1]+c[2])+h3*(c[0]+c[1]); break;
      case 3: s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2])+2.0*h3*c[2];
        s[1]=h0*c[1]+(h1+h2)*(c[0]+c[2])+2.0*h3*c[1];
        s[2]=h0*c[2]+h1*(c[1]+c[2])+h2*(c[0]+c[1])+2.0*h3*c[0]; break;
      case 2: s[0]=(h0+h1+h3)*c[0]+(h1+2.0*h2+h3)*c[1];
        s[1]=(h0+h1+h3)*c[1]+(h1+2.0*h2+h3)*c[0]; break;
      case 1: s[0]=(h0+2.0*(h1+h2+h3))*c[0]; break;
    }
  } // h.length===1 (degree 기본값): C++ 처럼 s 를 건드리지 않음
}
function antiSymmetricFir(h1,c,s,n){ // h={0,h1}
  s[0]=h1*(c[1]-c[0]); for(let i=1;i<n-1;i++) s[i]=h1*(c[i+1]-c[i-1]); s[n-1]=h1*(c[n-1]-c[n-2]);
}
function reduceDual1D(c,nc,s,ns){
  const h0=6.0/16.0,h1=4.0/16.0,h2=1.0/16.0;
  if(2<=ns){
    s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2]);
    for(let i=2,j=1;j<ns-1;i+=2,j++) s[j]=h0*c[i]+h1*(c[i-1]+c[i+1])+h2*(c[i-2]+c[i+2]);
    if(nc===2*ns) s[ns-1]=h0*c[nc-2]+h1*(c[nc-3]+c[nc-1])+h2*(c[nc-4]+c[nc-1]);
    else s[ns-1]=h0*c[nc-3]+h1*(c[nc-4]+c[nc-2])+h2*(c[nc-5]+c[nc-1]);
  } else if(nc===3) s[0]=h0*c[0]+h1*(c[0]+c[1])+h2*(c[1]+c[2]);
  else if(nc===2) s[0]=h0*c[0]+h1*(c[0]+c[1])+2.0*h2*c[1];
}
/* 2D: 행 → (float 저장) → 열 → (float 저장) */
function getBasicFromCardinal2D(card,w,h,degree,basic){
  basic=basic||new Float64Array(w*h); const hl=new Float64Array(w), vl=new Float64Array(h);
  for(let y=0;y<h;y++){ const o=y*w; for(let x=0;x<w;x++) hl[x]=card[o+x]; samplesToInterpolationCoefficient1D(hl,w,degree); for(let x=0;x<w;x++) basic[o+x]=f32(hl[x]); }
  for(let x=0;x<w;x++){ for(let y=0;y<h;y++) vl[y]=basic[y*w+x]; samplesToInterpolationCoefficient1D(vl,h,degree); for(let y=0;y<h;y++) basic[y*w+x]=f32(vl[y]); }
  return basic;
}
function basicToCardinal2D(basic,card,w,h,degree){
  const hh=degree===3?H3:degree===7?H7:H1;
  const hl=new Float64Array(w), hd=new Float64Array(w), vl=new Float64Array(h), vd=new Float64Array(h);
  for(let y=0;y<h;y++){ const o=y*w; for(let x=0;x<w;x++) hl[x]=basic[o+x]; symmetricFir(hh,hl,hd,w); for(let x=0;x<w;x++) card[o+x]=f32(hd[x]); }
  for(let x=0;x<w;x++){ for(let y=0;y<h;y++) vl[y]=card[y*w+x]; symmetricFir(hh,vl,vd,h); for(let y=0;y<h;y++) card[y*w+x]=f32(vd[y]); }
  return card;
}
function cardinalToDual2D(card,dual,w,h,degree){ basicToCardinal2D(getBasicFromCardinal2D(card,w,h,degree),dual,w,h,2*degree+1); }
function dualToCardinal2D(dual,card,w,h,degree){ basicToCardinal2D(getBasicFromCardinal2D(dual,w,h,2*degree+1),card,w,h,degree); }
function getHalfDual2D(full,fw,fh){
  const hw=trunc(fw/2), hh=trunc(fh/2);
  const hl=new Float64Array(fw), hd=new Float64Array(hw), vl=new Float64Array(fh), vd=new Float64Array(hh);
  const demi=new Float64Array(hw*fh), half=new Float64Array(hw*hh);
  for(let y=0;y<fh;y++){ const o=y*fw; for(let x=0;x<fw;x++) hl[x]=full[o+x]; reduceDual1D(hl,fw,hd,hw); for(let x=0;x<hw;x++) demi[y*hw+x]=f32(hd[x]); }
  for(let x=0;x<hw;x++){ for(let y=0;y<fh;y++) vl[y]=demi[y*hw+x]; reduceDual1D(vl,fh,vd,hh); for(let y=0;y<hh;y++) half[y*hw+x]=f32(vd[y]); }
  return half;
}
function coefficientToGradient1D(c,n,tmp){ antiSymmetricFir(0.5,c,tmp,n); c.set(tmp.subarray(0,n)); }
function coefficientToSamples1D(c,n,tmp){ symmetricFir(H3,c,tmp,n); c.set(tmp.subarray(0,n)); }
function coefficientToXYGradient2D(basic,w,h){
  const xg=new Float64Array(w*h), yg=new Float64Array(w*h);
  const hl=new Float64Array(w), hd=new Float64Array(w), vl=new Float64Array(h), tw=new Float64Array(w), th=new Float64Array(h);
  for(let y=0;y<h;y++){ const o=y*w; for(let x=0;x<w;x++) hl[x]=basic[o+x]; hd.set(hl);
    coefficientToGradient1D(hl,w,tw); coefficientToSamples1D(hd,w,tw);
    for(let x=0;x<w;x++){ xg[o+x]=f32(hl[x]); yg[o+x]=f32(hd[x]); } }
  for(let x=0;x<w;x++){
    for(let y=0;y<h;y++) vl[y]=xg[y*w+x]; coefficientToSamples1D(vl,h,th); for(let y=0;y<h;y++) xg[y*w+x]=f32(vl[y]);
    for(let y=0;y<h;y++) vl[y]=yg[y*w+x]; coefficientToGradient1D(vl,h,th); for(let y=0;y<h;y++) yg[y*w+x]=f32(vl[y]);
  }
  return [xg,yg];
}
function imageToXYGradient2D(img,w,h){
  const xg=new Float64Array(w*h), yg=new Float64Array(w*h), hl=new Float64Array(w), vl=new Float64Array(h), tw=new Float64Array(w), th=new Float64Array(h);
  for(let y=0;y<h;y++){ const o=y*w; for(let x=0;x<w;x++) hl[x]=img[o+x]; samplesToInterpolationCoefficient1D(hl,w,3); coefficientToGradient1D(hl,w,tw); for(let x=0;x<w;x++) xg[o+x]=f32(hl[x]); }
  for(let x=0;x<w;x++){ for(let y=0;y<h;y++) vl[y]=img[y*w+x]; samplesToInterpolationCoefficient1D(vl,h,3); coefficientToGradient1D(vl,h,th); for(let y=0;y<h;y++) yg[y*w+x]=f32(vl[y]); }
  return [xg,yg];
}
/* TurboRegImage::init. pyramid[0] = 가장 고운 반쪽 단계, 마지막 = 가장 거친 단계 (std::stack top) */
function makeImage(data,w,h,tr,isTarget,depth,coefOnly){
  const im={w,h,image:Float64Array.from(data),pyramid:[]};
  im.coefficient=getBasicFromCardinal2D(im.image,w,h,3);
  if(coefOnly) return im;
  const coefPyr=()=>{ let full=new Float64Array(w*h), hw=w, hh=h; if(1<depth) basicToCardinal2D(im.coefficient,full,w,h,7);
    for(let d=1;d<depth;d++){ const fw=hw, fh=hh; hw=trunc(hw/2); hh=trunc(hh/2); const half=getHalfDual2D(full,fw,fh);
      im.pyramid.push({w:hw,h:hh,img:getBasicFromCardinal2D(half,hw,hh,7)}); full=half; } };
  if(tr!==BILINEAR){
    if(isTarget) coefPyr();
    else { [im.xGradient,im.yGradient]=imageToXYGradient2D(im.image,w,h);
      let full=new Float64Array(w*h), hw=w, hh=h; if(1<depth) cardinalToDual2D(im.image,full,w,h,3);
      for(let d=1;d<depth;d++){ const fw=hw, fh=hh; hw=trunc(hw/2); hh=trunc(hh/2); const half=getHalfDual2D(full,fw,fh);
        const hi=getBasicFromCardinal2D(half,hw,hh,7); const [xg,yg]=coefficientToXYGradient2D(hi,hw,hh); basicToCardinal2D(hi,hi,hw,hh,3);
        im.pyramid.push({w:hw,h:hh,img:hi,xg,yg}); full=half; } }
  } else {
    if(isTarget){ let full=new Float64Array(w*h), hw=w, hh=h; if(1<depth) cardinalToDual2D(im.image,full,w,h,3);
      for(let d=1;d<depth;d++){ const fw=hw, fh=hh; hw=trunc(hw/2); hh=trunc(hh/2); const half=getHalfDual2D(full,fw,fh);
        const hi=new Float64Array(hw*hh); dualToCardinal2D(half,hi,hw,hh,3); im.pyramid.push({w:hw,h:hh,img:hi}); full=half; } }
    else coefPyr();
  }
  return im;
}

/* ===================== TurboRegMask (clearMask → 전부 1) ===================== */
function getHalfMask2D(F,fw,fh){
  const hw=trunc(fw/2), hh=trunc(fh/2), odd=(2*hw!==fw), M=new Float64Array(hw*hh), A=Math.abs;
  let k=0, n=0;
  for(let y=0;y<hh-1;y++){
    for(let x=0;x<hw-1;x++){ M[k]+=A(F[n++]); M[k]+=A(F[n]); M[++k]+=A(F[n++]); }
    M[k]+=A(F[n++]); M[k++]+=A(F[n++]);
    if(odd) n++;
    for(let x=0;x<hw-1;x++){ M[k-hw]+=A(F[n]); M[k]+=A(F[n++]); M[k-hw]+=A(F[n]); M[k-hw+1]+=A(F[n]); M[k]+=A(F[n]); M[++k]+=A(F[n++]); }
    M[k-hw]+=A(F[n]); M[k]+=A(F[n++]); M[k-hw]+=A(F[n]); M[k++]+=A(F[n++]);
    if(odd) n++;
    k-=hw;
  }
  for(let x=0;x<hw-1;x++){ M[k]+=A(F[n++]); M[k]+=A(F[n]); M[++k]+=A(F[n++]); }
  M[k]+=A(F[n++]); M[k++]+=A(F[n++]);
  if(odd) n++;
  k-=hw;
  for(let x=0;x<hw-1;x++){ M[k]+=A(F[n++]); M[k]+=A(F[n]); M[++k]+=A(F[n++]); }
  M[k]+=A(F[n++]); M[k]+=A(F[n]);
  return M;
}
function makeMask(w,h,depth){
  const mk={mask:new Float64Array(w*h).fill(1.0),pyramid:[]}; let full=mk.mask, hw=w, hh=h;
  for(let d=1;d<depth;d++){ const fw=hw, fh=hh; hw=trunc(hw/2); hh=trunc(hh/2); const half=getHalfMask2D(full,fw,fh); mk.pyramid.push(half); full=half; }
  return mk;
}

/* ===================== TurboRegPointHandler ===================== */
function initialPoints(width,height,tr){ // width/height = 잘린 크기 (C++ 에서 +1)
  width++; height++; const F=Math.floor, C=Math.ceil;
  switch(tr){
    case TRANSLATION: return [[F(width/2.0),F(height/2.0)]];
    case RIGID_BODY: return [[F(width/2.0),F(height/2.0)],[F(width/2.0),F(height/4.0)],[F(width/2.0),F((3.0*height)/4.0)]];
    case SCALED_ROTATION: return [[F(width/4.0),F(height/2.0)],[F((3.0*width)/4.0),F(height/2.0)]];
    case AFFINE: return [[F(width/2.0),F(height/4.0)],[F(width/4.0),F((3.0*height)/4.0)],[F((3.0*width)/4.0),F((3.0*height)/4.0)]];
    case BILINEAR: return [[F(0.25*GOLDEN_RATIO*width),F(0.25*GOLDEN_RATIO*height)],[F(0.25*GOLDEN_RATIO*width),height-C(0.25*GOLDEN_RATIO*height)],
      [width-C(0.25*GOLDEN_RATIO*width),F(0.25*GOLDEN_RATIO*height)],[width-C(0.25*GOLDEN_RATIO*width),height-C(0.25*GOLDEN_RATIO*height)]];
  }
}
const copyPts=p=>p.map(r=>r.slice());

/* ===================== TurboRegTransform ===================== */
function invertGauss(m,n){
  const inv=new Float64Array(n*n);
  for(let i=0;i<n;i++){ let max=m[i*n], absMax=Math.abs(max);
    for(let j=0;j<n;j++){ inv[i*n+j]=0.0; if(absMax<Math.abs(m[i*n+j])){ max=m[i*n+j]; absMax=Math.abs(max); } }
    inv[i*n+i]=1.0/max; for(let j=0;j<n;j++) m[i*n+j]/=max; }
  for(let j=0;j<n;j++){ let max=m[j*n+j], absMax=Math.abs(max), k=j;
    for(let i=j+1;i<n;i++) if(absMax<Math.abs(m[i*n+j])){ max=m[i*n+j]; absMax=Math.abs(max); k=i; }
    if(k!==j){ for(let c=j;c<n;c++){ const t=m[j*n+c]; m[j*n+c]=m[k*n+c]; m[k*n+c]=t; } for(let c=0;c<n;c++){ const t=inv[j*n+c]; inv[j*n+c]=inv[k*n+c]; inv[k*n+c]=t; } }
    for(k=0;k<=j;k++) inv[j*n+k]/=max;
    for(k=j+1;k<n;k++){ m[j*n+k]/=max; inv[j*n+k]/=max; }
    for(let i=j+1;i<n;i++){ const mij=m[i*n+j];
      for(k=0;k<=j;k++) inv[i*n+k]-=mij*inv[j*n+k];
      for(k=j+1;k<n;k++){ m[i*n+k]-=mij*m[j*n+k]; inv[i*n+k]-=mij*inv[j*n+k]; } }
  }
  for(let j=n-1;1<=j;j--) for(let i=j-1;0<=i;i--){ const mij=m[i*n+j];
    for(let k=0;k<=j;k++) inv[i*n+k]-=mij*inv[j*n+k];
    for(let k=j+1;k<n;k++){ m[i*n+k]-=mij*m[j*n+k]; inv[i*n+k]-=mij*inv[j*n+k]; } }
  m.set(inv);
}
function matrixMultiply(m,v,n){ const r=new Float64Array(n); for(let i=0;i<n;i++){ let s=0.0; for(let j=0;j<n;j++) s+=m[i*n+j]*v[j]; r[i]=s; } return r; }
/* getTransformationMatrix(from,to): 2×(1|3|4) 짧은 형 */
function getTM(tr,F,T){
  switch(tr){
    case TRANSLATION: return [[T[0][0]-F[0][0]],[T[0][1]-F[0][1]]];
    case RIGID_BODY: { const angle=Math.atan2(F[2][0]-F[1][0],F[2][1]-F[1][1])-Math.atan2(T[2][0]-T[1][0],T[2][1]-T[1][1]);
      const c=Math.cos(angle), s=Math.sin(angle);
      return [[T[0][0]-c*F[0][0]+s*F[0][1],c,-s],[T[0][1]-s*F[0][0]-c*F[0][1],s,c]]; }
    case SCALED_ROTATION: case AFFINE: { const a=new Float64Array(9);
      a[0]=1.0; a[1]=F[0][0]; a[2]=F[0][1]; a[3]=1.0; a[4]=F[1][0]; a[5]=F[1][1]; a[6]=1.0;
      if(tr===SCALED_ROTATION){ a[7]=F[0][1]-F[1][1]+F[1][0]; a[8]=F[1][0]+F[1][1]-F[0][0]; } else { a[7]=F[2][0]; a[8]=F[2][1]; }
      invertGauss(a,3);
      const v0=tr===SCALED_ROTATION?[T[0][0],T[1][0],T[0][1]-T[1][1]+T[1][0]]:[T[0][0],T[1][0],T[2][0]];
      const v1=tr===SCALED_ROTATION?[T[0][1],T[1][1],T[1][0]+T[1][1]-T[0][0]]:[T[0][1],T[1][1],T[2][1]];
      const m=[[0,0,0],[0,0,0]];
      for(let i=0;i<3;i++){ let s=0.0; for(let j=0;j<3;j++) s+=a[i*3+j]*v0[j]; m[0][i]=s; }
      for(let i=0;i<3;i++){ let s=0.0; for(let j=0;j<3;j++) s+=a[i*3+j]*v1[j]; m[1][i]=s; }
      return m; }
    case BILINEAR: { const a=new Float64Array(16);
      for(let r=0;r<4;r++){ a[r*4]=1.0; a[r*4+1]=F[r][0]; a[r*4+2]=F[r][1]; a[r*4+3]=F[r][0]*F[r][1]; }
      invertGauss(a,4); const m=[[0,0,0,0],[0,0,0,0]];
      for(let row=0;row<2;row++) for(let i=0;i<4;i++){ let s=0.0; for(let j=0;j<4;j++) s+=a[i*4+j]*T[j][row]; m[row][i]=s; }
      return m; }
  }
}

/* B-spline 3차 보간 상태 (C++ 멤버 xIndex/yIndex/xWeight/yWeight/dxWeight/dyWeight) */
const XI=new Int32Array(4), YI=new Int32Array(4), XW=new Float64Array(4), YW=new Float64Array(4), DXW=new Float64Array(4), DYW=new Float64Array(4);
// 색인은 (p, nx, ny) 로만 정해지므로 직전 값과 같으면 다시 계산하지 않는다 (결과 동일, 속도용)
let XP=NaN, XN=-1, YP=NaN, YN=-1, YNY=-1;
function xIndexes(x,nx,tnx){ let p=(0.0<=x)?(trunc(x)+2):(trunc(x)+1); if(p===XP&&nx===XN) return; XP=p; XN=nx;
  for(let k=0;k<4;k++,p--){ let q=(p<0)?(-1-p):p; if(tnx<=q) q-=tnx*trunc(q/tnx); XI[k]=(nx<=q)?(tnx-1-q):q; } }
function yIndexes(y,nx,ny,tny){ let p=(0.0<=y)?(trunc(y)+2):(trunc(y)+1); if(p===YP&&nx===YN&&ny===YNY) return; YP=p; YN=nx; YNY=ny;
  for(let k=0;k<4;k++,p--){ let q=(p<0)?(-1-p):p; if(tny<=q) q-=tny*trunc(q/tny); YI[k]=(ny<=q)?((tny-1-q)*nx):(q*nx); } }
function xWeights(x){ let s=1.0-x; XW[3]=s*s*s/6.0; s=x*x; XW[2]=2.0/3.0-0.5*s*(2.0-x); XW[0]=s*x/6.0; XW[1]=1.0-XW[0]-XW[2]-XW[3]; }
function yWeights(y){ let t=1.0-y; YW[3]=t*t*t/6.0; t=y*y; YW[2]=2.0/3.0-0.5*t*(2.0-y); YW[0]=t*y/6.0; YW[1]=1.0-YW[0]-YW[2]-YW[3]; }
function xDxWeights(x){ const s=1.0-x; DXW[0]=0.5*x*x; XW[0]=x*DXW[0]/3.0; DXW[3]=-0.5*s*s; XW[3]=s*DXW[3]/-3.0;
  DXW[1]=1.0-2.0*DXW[0]+DXW[3]; XW[1]=2.0/3.0+(1.0+x)*DXW[3]; DXW[2]=1.5*x*(x-4.0/3.0); XW[2]=2.0/3.0-(2.0-x)*DXW[0]; }
function yDyWeights(y){ const t=1.0-y; DYW[0]=0.5*y*y; YW[0]=y*DYW[0]/3.0; DYW[3]=-0.5*t*t; YW[3]=t*DYW[3]/-3.0;
  DYW[1]=1.0-2.0*DYW[0]+DYW[3]; YW[1]=2.0/3.0+(1.0+y)*DYW[3]; DYW[2]=1.5*y*(y-4.0/3.0); YW[2]=2.0/3.0-(2.0-y)*DYW[0]; }
/* interpolate: s=Σ_i xW[i]·c[p+xI[i]] (i 순서), t=Σ_j yW[j]·s_j (j 순서) — C++ 와 같은 덧셈 순서 */
function interpolate(c){ const x0=XI[0],x1=XI[1],x2=XI[2],x3=XI[3], w0=XW[0],w1=XW[1],w2=XW[2],w3=XW[3]; let p, t;
  p=YI[0]; t=YW[0]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[1]; t+=YW[1]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[2]; t+=YW[2]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[3]; t+=YW[3]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  return t; }
function interpolateDx(c){ const x0=XI[0],x1=XI[1],x2=XI[2],x3=XI[3], w0=DXW[0],w1=DXW[1],w2=DXW[2],w3=DXW[3]; let p, t;
  p=YI[0]; t=YW[0]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[1]; t+=YW[1]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[2]; t+=YW[2]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[3]; t+=YW[3]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  return t; }
function interpolateDy(c){ const x0=XI[0],x1=XI[1],x2=XI[2],x3=XI[3], w0=XW[0],w1=XW[1],w2=XW[2],w3=XW[3]; let p, t;
  p=YI[0]; t=DYW[0]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[1]; t+=DYW[1]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[2]; t+=DYW[2]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  p=YI[3]; t+=DYW[3]*(w0*c[p+x0]+w1*c[p+x1]+w2*c[p+x2]+w3*c[p+x3]);
  return t; }
const rnd=v=>(0.0<=v)?trunc(v+0.5):trunc(v-0.5);   // xMsk/yMsk
const frac=v=>v-((0.0<=v)?trunc(v):(trunc(v)-1));  // x -= (0<=x)?(int)x:((int)x-1)

/* 평균제곱 + 기울기 + 헤시안. H===null 이면 기울기/헤시안 없는 판 (C++ 의 인자 1·2 개 판) */
function symH(H,n){ for(let i=1;i<n;i++) for(let j=0;j<i;j++) H[i*n+j]=H[j*n+i]; }
function translationMS(S,m,H,g){
  const inImg=S.inImg,inMsk=S.inMsk,outImg=S.outImg,outMsk=S.outMsk,xG=S.xGradient,yG=S.yGradient,inNx=S.inNx,inNy=S.inNy,tnx=S.twiceInNx,tny=S.twiceInNy,outNx=S.outNx,outNy=S.outNy;
  let dx=m[0][0], dy=m[1][0]; const dx0=dx; let ms=0.0, area=0, k=0;
  let g0=0,g1=0,h00=0,h01=0,h11=0;
  xWeights(dx-Math.floor(dx)); yWeights(dy-Math.floor(dy));
  for(let v=0;v<outNy;v++){ const y=dy; dy+=1; let yMsk=rnd(y);
    if(0<=yMsk&&yMsk<inNy){ yMsk*=inNx; yIndexes(y,inNx,inNy,tny); dx=dx0;
      for(let u=0;u<outNx;u++,k++){ const x=dx; dx+=1; const xMsk=rnd(x);
        if(0<=xMsk&&xMsk<inNx&&(outMsk[k]*inMsk[yMsk+xMsk])!==0.0){ area++; xIndexes(x,inNx,tnx);
          const d=outImg[k]-interpolate(inImg); ms+=d*d;
          if(H){ const gx=xG[k], gy=yG[k]; g0+=d*gx; g1+=d*gy; h00+=gx*gx; h01+=gx*gy; h11+=gy*gy; } } }
    } else k+=outNx;
  }
  if(H){ g[0]=g0; g[1]=g1; H[0]=h00; H[1]=h01; H[3]=h11; symH(H,2); }
  return ms/area;
}
function rigidMS(S,m,H,g){
  const inImg=S.inImg,inMsk=S.inMsk,outImg=S.outImg,outMsk=S.outMsk,xG=S.xGradient,yG=S.yGradient,inNx=S.inNx,inNy=S.inNy,tnx=S.twiceInNx,tny=S.twiceInNy,outNx=S.outNx,outNy=S.outNy;
  const m01=m[0][1],m11=m[1][1],m02=m[0][2],m12=m[1][2]; let yx=m[0][0], yy=m[1][0]; let ms=0.0, area=0, k=0;
  let g0=0,g1=0,g2=0,h00=0,h01=0,h02=0,h11=0,h12=0,h22=0;
  for(let v=0;v<outNy;v++){ let x0=yx, y0=yy;
    for(let u=0;u<outNx;u++,k++){ const x=x0, y=y0; const xMsk=rnd(x), yMsk=rnd(y);
      if(0<=xMsk&&xMsk<inNx&&0<=yMsk&&yMsk<inNy&&(outMsk[k]*inMsk[yMsk*inNx+xMsk])!==0.0){ area++;
        xIndexes(x,inNx,tnx); yIndexes(y,inNx,inNy,tny); xWeights(frac(x)); yWeights(frac(y));
        const d=outImg[k]-interpolate(inImg); ms+=d*d;
        if(H){ const gx=xG[k], gy=yG[k]; const dT=gy*u-gx*v;
          g0+=d*dT; g1+=d*gx; g2+=d*gy; h00+=dT*dT; h01+=dT*gx; h02+=dT*gy; h11+=gx*gx; h12+=gx*gy; h22+=gy*gy; } }
      x0+=m01; y0+=m11; }
    yx+=m02; yy+=m12; }
  if(H){ g[0]=g0; g[1]=g1; g[2]=g2; H[0]=h00; H[1]=h01; H[2]=h02; H[4]=h11; H[5]=h12; H[8]=h22; symH(H,3); }
  return ms/area;
}
function scaledMS(S,sp,m,H,g){
  const inImg=S.inImg,inMsk=S.inMsk,outImg=S.outImg,outMsk=S.outMsk,xG=S.xGradient,yG=S.yGradient,inNx=S.inNx,inNy=S.inNy,tnx=S.twiceInNx,tny=S.twiceInNy,outNx=S.outNx,outNy=S.outNy;
  const u1=sp[0][0],u2=sp[1][0],v1=sp[0][1],v2=sp[1][1]; const u12=u1-u2, v12=v1-v2, uv2=u12*u12+v12*v12;
  const c=0.5*(u2*v1-u1*v2)/uv2, c1=u12/uv2, c2=v12/uv2, c3=(uv2-u12*v12)/uv2, c4=(uv2+u12*v12)/uv2;
  const c5=c+u1*c1+u2*c2, c6=c*(u12*u12-v12*v12)/uv2, c7=c1*c4, c8=c1-c2-c1*c2*v12, c9=c1+c2-c1*c2*u12, c0=c2*c3;
  const dgxx0=c1*u2+c2*v2, dgyx0=2.0*c, dgxx1=c5+c6, dgyy1=c5-c6;
  const m01=m[0][1],m11=m[1][1],m02=m[0][2],m12=m[1][2]; let yx=m[0][0], yy=m[1][0]; let ms=0.0, area=0, k=0;
  let g0=0,g1=0,g2=0,g3=0,h00=0,h01=0,h02=0,h03=0,h11=0,h12=0,h13=0,h22=0,h23=0,h33=0;
  for(let v=0;v<outNy;v++){ let x0=yx, y0=yy;
    for(let u=0;u<outNx;u++,k++){ const x=x0, y=y0; const xMsk=rnd(x), yMsk=rnd(y);
      if(0<=xMsk&&xMsk<inNx&&0<=yMsk&&yMsk<inNy&&(outMsk[k]*inMsk[yMsk*inNx+xMsk])!==0.0){ area++;
        xIndexes(x,inNx,tnx); yIndexes(y,inNx,inNy,tny); xWeights(frac(x)); yWeights(frac(y));
        const d=outImg[k]-interpolate(inImg); ms+=d*d;
        if(H){ const gx=xG[k], gy=yG[k];
          const gxx0=u*c1+v*c2-dgxx0, gyx0=v*c1-u*c2+dgyx0, gxy0=-gyx0, gyy0=gxx0;
          const gxx1=v*c8-u*c7+dgxx1, gyx1=-c3*gyx0, gxy1=c4*gyx0, gyy1=dgyy1-u*c9-v*c0;
          const dx0=gx*gxx0+gy*gyx0, dy0=gx*gxy0+gy*gyy0, dx1=gx*gxx1+gy*gyx1, dy1=gx*gxy1+gy*gyy1;
          g0+=d*dx0; g1+=d*dy0; g2+=d*dx1; g3+=d*dy1;
          h00+=dx0*dx0; h01+=dx0*dy0; h02+=dx0*dx1; h03+=dx0*dy1; h11+=dy0*dy0; h12+=dy0*dx1; h13+=dy0*dy1; h22+=dx1*dx1; h23+=dx1*dy1; h33+=dy1*dy1; } }
      x0+=m01; y0+=m11; }
    yx+=m02; yy+=m12; }
  if(H){ g[0]=g0; g[1]=g1; g[2]=g2; g[3]=g3; H[0]=h00; H[1]=h01; H[2]=h02; H[3]=h03; H[5]=h11; H[6]=h12; H[7]=h13; H[10]=h22; H[11]=h23; H[15]=h33; symH(H,4); }
  return ms/(area*uv2/S.targetJacobian);
}
function affineMS(S,sp,m,H,g){
  const inImg=S.inImg,inMsk=S.inMsk,outImg=S.outImg,outMsk=S.outMsk,xG=S.xGradient,yG=S.yGradient,inNx=S.inNx,inNy=S.inNy,tnx=S.twiceInNx,tny=S.twiceInNy,outNx=S.outNx,outNy=S.outNy;
  const u1=sp[0][0],u2=sp[1][0],u3=sp[2][0],v1=sp[0][1],v2=sp[1][1],v3=sp[2][1];
  let uv32=u3*v2-u2*v3, uv21=u2*v1-u1*v2, uv13=u1*v3-u3*v1; const det=uv32+uv21+uv13;
  const u12=(u1-u2)/det, u23=(u2-u3)/det, u31=(u3-u1)/det, v12=(v1-v2)/det, v23=(v2-v3)/det, v31=(v3-v1)/det;
  uv32/=det; uv21/=det; uv13/=det;
  const m01=m[0][1],m11=m[1][1],m02=m[0][2],m12=m[1][2]; let yx=m[0][0], yy=m[1][0]; let ms=0.0, area=0, k=0;
  const gr=new Float64Array(6), A=H?new Float64Array(36):null;
  for(let v=0;v<outNy;v++){ let x0=yx, y0=yy;
    for(let u=0;u<outNx;u++,k++){ const x=x0, y=y0; const xMsk=rnd(x), yMsk=rnd(y);
      if(0<=xMsk&&xMsk<inNx&&0<=yMsk&&yMsk<inNy&&(outMsk[k]*inMsk[yMsk*inNx+xMsk])!==0.0){ area++;
        xIndexes(x,inNx,tnx); yIndexes(y,inNx,inNy,tny); xWeights(frac(x)); yWeights(frac(y));
        const d=outImg[k]-interpolate(inImg); ms+=d*d;
        if(H){ const gx=xG[k], gy=yG[k];
          const q0=u23*v-v23*u+uv32, q1=u31*v-v31*u+uv13, q2=u12*v-v12*u+uv21;
          const dx0=gx*q0, dy0=gy*q0, dx1=gx*q1, dy1=gy*q1, dx2=gx*q2, dy2=gy*q2;
          gr[0]+=d*dx0; gr[1]+=d*dy0; gr[2]+=d*dx1; gr[3]+=d*dy1; gr[4]+=d*dx2; gr[5]+=d*dy2;
          A[0]+=dx0*dx0; A[1]+=dx0*dy0; A[2]+=dx0*dx1; A[3]+=dx0*dy1; A[4]+=dx0*dx2; A[5]+=dx0*dy2;
          A[7]+=dy0*dy0; A[8]+=dy0*dx1; A[9]+=dy0*dy1; A[10]+=dy0*dx2; A[11]+=dy0*dy2;
          A[14]+=dx1*dx1; A[15]+=dx1*dy1; A[16]+=dx1*dx2; A[17]+=dx1*dy2;
          A[21]+=dy1*dy1; A[22]+=dy1*dx2; A[23]+=dy1*dy2;
          A[28]+=dx2*dx2; A[29]+=dx2*dy2; A[35]+=dy2*dy2; } }
      x0+=m01; y0+=m11; }
    yx+=m02; yy+=m12; }
  if(H){ g.set(gr); H.set(A); symH(H,6); }
  return ms/(area*Math.abs(det/S.targetJacobian));
}
function bilinearConstants(tp){
  const u1=tp[0][0],u2=tp[1][0],u3=tp[2][0],u4=tp[3][0],v1=tp[0][1],v2=tp[1][1],v3=tp[2][1],v4=tp[3][1];
  const v12=v1-v2,v13=v1-v3,v14=v1-v4,v23=v2-v3,v24=v2-v4,v34=v3-v4;
  const uv12=u1*u2*v12,uv13=u1*u3*v13,uv14=u1*u4*v14,uv23=u2*u3*v23,uv24=u2*u4*v24,uv34=u3*u4*v34;
  const det=uv12*v34-uv13*v24+uv14*v23+uv23*v14-uv24*v13+uv34*v12;
  return [
    (-uv34*v2+uv24*v3-uv23*v4)/det, (u3*v3*v24-u2*v2*v34-u4*v4*v23)/det, (uv23-uv24+uv34)/det, (u4*v23-u3*v24+u2*v34)/det,
    (uv34*v1-uv14*v3+uv13*v4)/det, (-u3*v3*v14+u1*v1*v34+u4*v4*v13)/det, (-uv13+uv14-uv34)/det, (-u4*v13+u3*v14-u1*v34)/det,
    (-uv24*v1+uv14*v2-uv12*v4)/det, (u2*v2*v14-u1*v1*v24-u4*v4*v12)/det, (uv12-uv14+uv24)/det, (u4*v12-u2*v14+u1*v24)/det,
    (uv23*v1-uv13*v2+uv12*v3)/det, (-u2*v2*v13+u1*v1*v23+u3*v3*v12)/det, (-uv12+uv13-uv23)/det, (-u3*v1+u2*v13+u3*v2-u1*v23)/det ];
}
function bilinearMS(S,m,H,g){ // (inMsk 가 비지 않은 판) difference = interpolate − out
  const inImg=S.inImg,inMsk=S.inMsk,outImg=S.outImg,outMsk=S.outMsk,inNx=S.inNx,inNy=S.inNy,tnx=S.twiceInNx,tny=S.twiceInNy,outNx=S.outNx,outNy=S.outNy;
  const C=H?bilinearConstants(S.targetPoint):null;
  const m01=m[0][1],m11=m[1][1],m02=m[0][2],m12=m[1][2],m03=m[0][3],m13=m[1][3];
  let yx=m[0][0], yy=m[1][0], yxy=0.0, yyy=0.0; let ms=0.0, area=0, k=0;
  const gr=new Float64Array(8), A=H?new Float64Array(64):null, D=new Float64Array(8);
  for(let v=0;v<outNy;v++){ let x0=yx, y0=yy;
    for(let u=0;u<outNx;u++,k++){ const x=x0, y=y0; const xMsk=rnd(x), yMsk=rnd(y);
      if(0<=xMsk&&xMsk<inNx&&0<=yMsk&&yMsk<inNy&&(outMsk[k]*inMsk[yMsk*inNx+xMsk])!==0.0){ area++;
        xIndexes(x,inNx,tnx); yIndexes(y,inNx,inNy,tny);
        if(H){ xDxWeights(frac(x)); yDyWeights(frac(y)); } else { xWeights(frac(x)); yWeights(frac(y)); }
        const d=interpolate(inImg)-outImg[k]; ms+=d*d;
        if(H){ const gx=interpolateDx(inImg), gy=interpolateDy(inImg); const uv=u*v;
          for(let i=0;i<4;i++){ const gi=C[4*i+3]*uv+C[4*i+1]*u+C[4*i+2]*v+C[4*i]; D[2*i]=gx*gi; D[2*i+1]=gy*gi; }
          for(let i=0;i<8;i++){ gr[i]+=d*D[i]; const Di=D[i]; for(let j=i;j<8;j++) A[i*8+j]+=Di*D[j]; } } }
      x0+=m01+yxy; y0+=m11+yyy; }
    yx+=m02; yy+=m12; yxy+=m03; yyy+=m13; }
  if(H){ g.set(gr); H.set(A); symH(H,8); }
  return ms/area;
}
/* 주: bilinear 의 기울기 없는 판(getBilinearMeanSquares(m))은 xWeights/yWeights, 기울기 판은 xDxWeights/yDyWeights 로
   보간 가중치를 만든다 (두 식은 수학적으로 같지만 반올림이 다름) — C++ 그대로 나눔. */

function inverseML(S){ // inverseMarquardtLevenbergOptimization (TRANSLATION / SCALED_ROTATION / AFFINE)
  const tr=S.tr, n=tr, np=n/2, sp=S.sourcePoint, tp=S.targetPoint;
  const attempt=[...Array(np)].map(()=>[0,0]);
  const hessian=new Float64Array(n*n), pseudo=new Float64Array(n*n), gradient=new Float64Array(n);
  const MS=(pts,m,H,g)=>tr===TRANSLATION?translationMS(S,m,H,g):tr===SCALED_ROTATION?scaledMS(S,pts,m,H,g):affineMS(S,pts,m,H,g);
  let m=getTM(tr,sp,tp), ms=0.0, lambda=FIRST_LAMBDA, displacement, iteration=0, update;
  let best=MS(sp,m,hessian,gradient); iteration++;
  do {
    for(let k=0;k<n;k++) pseudo[k*n+k]=(1.0+lambda)*hessian[k*n+k]; // C++ 판은 대각만 갱신 (비대각은 이전 역행렬 값 유지)
    invertGauss(pseudo,n); update=matrixMultiply(pseudo,gradient,n);
    displacement=0.0;
    for(let k=0;k<np;k++){ attempt[k][0]=sp[k][0]-update[2*k]; attempt[k][1]=sp[k][1]-update[2*k+1];
      displacement+=Math.sqrt(update[2*k]*update[2*k]+update[2*k+1]*update[2*k+1]); }
    displacement/=0.5*n;
    m=getTM(tr,attempt,tp); ms=MS(attempt,m,hessian,gradient); iteration++;
    if(ms<best){ best=ms; for(let k=0;k<np;k++){ sp[k][0]=attempt[k][0]; sp[k][1]=attempt[k][1]; } lambda/=LAMBDA_MAGSTEP; }
    else lambda*=LAMBDA_MAGSTEP;
  } while(iteration<(S.maxIterations*S.iterationPower-1)&&S.pixelPrecision<=displacement);
  invertGauss(hessian,n); update=matrixMultiply(hessian,gradient,n);
  for(let k=0;k<np;k++){ attempt[k][0]=sp[k][0]-update[2*k]; attempt[k][1]=sp[k][1]-update[2*k+1]; }
  m=getTM(tr,attempt,tp); ms=MS(attempt,m,null,null);
  if(ms<best) for(let k=0;k<np;k++){ sp[k][0]=attempt[k][0]; sp[k][1]=attempt[k][1]; }
}
function inverseMLRigid(S){ // inverseMarquardtLevenbergRigidBodyOptimization
  const n=3, sp=S.sourcePoint, tp=S.targetPoint;
  const attempt=[[0,0,0],[0,0,0]], hessian=new Float64Array(9), pseudo=new Float64Array(9), gradient=new Float64Array(3);
  let m=getTM(RIGID_BODY,tp,sp), ms=0.0, lambda=FIRST_LAMBDA, displacement, iteration=0, update;
  for(let k=0;k<n;k++){ const a=m[0][0]+tp[k][0]*m[0][1]+tp[k][1]*m[0][2], b=m[1][0]+tp[k][0]*m[1][1]+tp[k][1]*m[1][2]; sp[k][0]=a; sp[k][1]=b; }
  m=getTM(RIGID_BODY,sp,tp);
  let best=rigidMS(S,m,hessian,gradient); iteration++;
  const step=()=>{ const angle=Math.atan2(m[0][2],m[0][1])-update[0];
    attempt[0][1]=Math.cos(angle); attempt[0][2]=Math.sin(angle); attempt[1][1]=-attempt[0][2]; attempt[1][2]=attempt[0][1];
    const c=Math.cos(update[0]), s=Math.sin(update[0]);
    attempt[0][0]=(m[0][0]+update[1])*c-(m[1][0]+update[2])*s; attempt[1][0]=(m[0][0]+update[1])*s+(m[1][0]+update[2])*c; };
  do {
    for(let k=0;k<n;k++) pseudo[k*n+k]=(1.0+lambda)*hessian[k*n+k];
    invertGauss(pseudo,n); update=matrixMultiply(pseudo,gradient,n); step();
    displacement=Math.sqrt(update[1]*update[1]+update[2]*update[2])+0.25*Math.sqrt((S.inNx*S.inNx)+(S.inNy*S.inNy))*Math.abs(update[0]);
    ms=rigidMS(S,attempt,hessian,gradient); iteration++;
    if(ms<best){ best=ms; for(let i=0;i<2;i++) for(let j=0;j<3;j++) m[i][j]=attempt[i][j]; lambda/=LAMBDA_MAGSTEP; }
    else lambda*=LAMBDA_MAGSTEP;
  } while(iteration<(S.maxIterations*S.iterationPower-1)&&S.pixelPrecision<=displacement);
  invertGauss(hessian,n); update=matrixMultiply(hessian,gradient,n); step();
  ms=rigidMS(S,attempt,null,null);
  if(ms<best) for(let i=0;i<2;i++) for(let j=0;j<3;j++) m[i][j]=attempt[i][j];
  for(let k=0;k<n;k++){ const a=(tp[k][0]-m[0][0])*m[0][1]+(tp[k][1]-m[1][0])*m[1][1], b=(tp[k][0]-m[0][0])*m[0][2]+(tp[k][1]-m[1][0])*m[1][2]; sp[k][0]=a; sp[k][1]=b; }
}
function marquardtLevenberg(S){ // BILINEAR
  const n=8, np=4, sp=S.sourcePoint, tp=S.targetPoint;
  const attempt=[[0,0],[0,0],[0,0],[0,0]], hessian=new Float64Array(64), pseudo=new Float64Array(64), gradient=new Float64Array(8);
  let m=getTM(BILINEAR,tp,sp), ms=0.0, lambda=FIRST_LAMBDA, displacement, iteration=0, update;
  let best=bilinearMS(S,m,hessian,gradient); iteration++;
  do {
    for(let k=0;k<n;k++) pseudo[k*n+k]=(1.0+lambda)*hessian[k*n+k];
    invertGauss(pseudo,n); update=matrixMultiply(pseudo,gradient,n);
    displacement=0.0;
    for(let k=0;k<np;k++){ attempt[k][0]=sp[k][0]-update[2*k]; attempt[k][1]=sp[k][1]-update[2*k+1];
      displacement+=Math.sqrt(update[2*k]*update[2*k]+update[2*k+1]*update[2*k+1]); }
    displacement/=0.5*n;
    m=getTM(BILINEAR,tp,attempt); ms=bilinearMS(S,m,hessian,gradient); iteration++;
    if(ms<best){ best=ms; for(let k=0;k<np;k++){ sp[k][0]=attempt[k][0]; sp[k][1]=attempt[k][1]; } lambda/=LAMBDA_MAGSTEP; }
    else lambda*=LAMBDA_MAGSTEP;
  } while(iteration<(S.maxIterations*S.iterationPower-1)&&S.pixelPrecision<=displacement);
  invertGauss(hessian,n); update=matrixMultiply(hessian,gradient,n);
  for(let k=0;k<np;k++){ attempt[k][0]=sp[k][0]-update[2*k]; attempt[k][1]=sp[k][1]-update[2*k+1]; }
  m=getTM(BILINEAR,tp,attempt); ms=bilinearMS(S,m,null,null);
  if(ms<best) for(let k=0;k<np;k++){ sp[k][0]=attempt[k][0]; sp[k][1]=attempt[k][1]; }
}
function scaleLandmarks(S,f){ const np=S.tr===RIGID_BODY?3:S.tr/2; for(let i=0;i<np;i++){ S.sourcePoint[i][0]*=f; S.sourcePoint[i][1]*=f; S.targetPoint[i][0]*=f; S.targetPoint[i][1]*=f; } }
function optimize(S){
  switch(S.tr){
    case TRANSLATION: case SCALED_ROTATION: case AFFINE: inverseML(S); break;
    case RIGID_BODY: inverseMLRigid(S); break;
    case BILINEAR: marquardtLevenberg(S); break;
  }
}
/* TurboRegTransform::doRegistration (accelerated=false) */
function doRegistration(src,srcMsk,tgt,tgtMsk,tr,sourcePoint,targetPoint,depth){
  const S={tr,sourcePoint,targetPoint,pixelPrecision:PIXEL_HIGH_PRECISION,maxIterations:MANY_ITERATIONS,targetJacobian:1.0};
  // targetJacobian: C++ 멤버. 피라미드 단계에서만 계산되고 마지막 전체 해상도 단계에서는 직전 값을 그대로 씀.
  // (깊이 1 이면 C++ 에선 초기화되지 않은 값 — 여기서는 1.0. 이 값은 같은 단계 안의 비교에 상수배로만 들어감)
  S.iterationPower=trunc(Math.pow(ITERATION_PROGRESSION,depth));
  for(let d=1;d<depth;d++) scaleLandmarks(S,0.5);
  for(let L=src.pyramid.length-1;L>=0;L--){
    const si=src.pyramid[L], sm=srcMsk.pyramid[L], ti=tgt.pyramid[L], tm=tgtMsk.pyramid[L];
    S.iterationPower=trunc(S.iterationPower/ITERATION_PROGRESSION);
    if(tr===BILINEAR){ S.inNx=si.w; S.inNy=si.h; S.inImg=si.img; S.inMsk=sm; S.outNx=ti.w; S.outNy=ti.h; S.outImg=ti.img; S.outMsk=tm; }
    else { S.inNx=ti.w; S.inNy=ti.h; S.inImg=ti.img; S.inMsk=tm; S.outNx=si.w; S.outNy=si.h; S.outImg=si.img; S.xGradient=si.xg; S.yGradient=si.yg; S.outMsk=sm; }
    S.twiceInNx=2*S.inNx; S.twiceInNy=2*S.inNy;
    const tp=S.targetPoint;
    if(tr===TRANSLATION) S.targetJacobian=1.0;
    else if(tr===SCALED_ROTATION) S.targetJacobian=(tp[0][0]-tp[1][0])*(tp[0][0]-tp[1][0])+(tp[0][1]-tp[1][1])*(tp[0][1]-tp[1][1]);
    else if(tr===AFFINE) S.targetJacobian=(tp[1][0]-tp[2][0])*tp[0][1]+(tp[2][0]-tp[0][0])*tp[1][1]+(tp[0][0]-tp[1][0])*tp[2][1];
    optimize(S);
    scaleLandmarks(S,2.0);
  }
  S.iterationPower=trunc(S.iterationPower/ITERATION_PROGRESSION);
  if(tr===BILINEAR){ S.inNx=src.w; S.inNy=src.h; S.inImg=src.coefficient; S.inMsk=srcMsk.mask; S.outNx=tgt.w; S.outNy=tgt.h; S.outImg=tgt.image; S.outMsk=tgtMsk.mask; }
  else { S.inNx=tgt.w; S.inNy=tgt.h; S.inImg=tgt.coefficient; S.inMsk=tgtMsk.mask; S.outNx=src.w; S.outNy=src.h; S.outImg=src.image; S.xGradient=src.xGradient; S.yGradient=src.yGradient; S.outMsk=srcMsk.mask; }
  S.twiceInNx=2*S.inNx; S.twiceInNy=2*S.inNy;
  optimize(S);
  return S;
}

/* pymain registerImg: ref/mov = Float64 배열 (w×h, 이미 잘린 것) */
function registerRaw(ref,mov,w,h,tr){
  const depth=getPyramidDepth(w,h,w,h);
  const refImg=makeImage(ref,w,h,tr,true,depth), movImg=makeImage(mov,w,h,tr,false,depth);
  const refPts=initialPoints(w,h,tr), movPts=initialPoints(w,h,tr);
  const refMsk=makeMask(w,h,depth), movMsk=makeMask(w,h,depth);
  const S=doRegistration(movImg,movMsk,refImg,refMsk,tr,copyPts(movPts),copyPts(refPts),depth);
  return {short:getTM(tr,S.targetPoint,S.sourcePoint),refPts,movPts:copyPts(S.sourcePoint)};
}
/* pystackreg _matrix_short_to_long / _matrix_long_to_short */
function shortToLong(m,tr){
  if(tr===TRANSLATION) return [[1,0,m[0][0]],[0,1,m[1][0]],[0,0,1]];
  if(tr===BILINEAR) return [[m[0][1],m[0][2],m[0][3],m[0][0]],[m[1][1],m[1][2],m[1][3],m[1][0]],[0,0,1,0],[0,0,0,1]];
  return [[m[0][1],m[0][2],m[0][0]],[m[1][1],m[1][2],m[1][0]],[0,0,1]];
}
function longToShort(M,tr){
  if(tr===TRANSLATION) return [[+M[0][2]],[+M[1][2]]];
  if(tr===BILINEAR) return [[+M[0][3],+M[0][0],+M[0][1],+M[0][2]],[+M[1][3],+M[1][0],+M[1][1],+M[1][2]]];
  return [[+M[0][2],+M[0][0],+M[0][1]],[+M[1][2],+M[1][0],+M[1][1]]];
}
/* TurboRegTransform::doFinalTransform(sourceImg, m) — translationTransform / affineTransform / bilinearTransform (마스크 없는 판) */
function transformRaw(data,w,h,m){
  const cols=m[0].length, tr=cols===1?TRANSLATION:cols===3?AFFINE:cols===4?BILINEAR:0; if(!tr) throw new Error("StackReg: 행렬 모양이 잘못됨");
  const im=makeImage(data,w,h,tr,false,1,true), c=im.coefficient, out=new Float32Array(w*h);
  const inNx=w, inNy=h, tnx=2*w, tny=2*h; let k=0;
  if(tr===TRANSLATION){
    let dx=m[0][0], dy=m[1][0]; const dx0=dx; xWeights(dx-Math.floor(dx)); yWeights(dy-Math.floor(dy));
    for(let v=0;v<h;v++){ const y=dy; dy+=1; const yMsk=rnd(y);
      if(0<=yMsk&&yMsk<inNy){ yIndexes(y,inNx,inNy,tny); dx=dx0;
        for(let u=0;u<w;u++){ const x=dx; dx+=1; const xMsk=rnd(x);
          if(0<=xMsk&&xMsk<inNx){ xIndexes(x,inNx,tnx); out[k++]=interpolate(c); } else out[k++]=0.0; } }
      else for(let u=0;u<w;u++) out[k++]=0.0; }
    return out;
  }
  const m00=m[0][0],m10=m[1][0],m01=m[0][1],m11=m[1][1],m02=m[0][2],m12=m[1][2], m03=tr===BILINEAR?m[0][3]:0, m13=tr===BILINEAR?m[1][3]:0;
  let yx=m00, yy=m10, yxy=0.0, yyy=0.0;
  for(let v=0;v<h;v++){ let x0=yx, y0=yy;
    for(let u=0;u<w;u++){ const x=x0, y=y0; const xMsk=rnd(x), yMsk=rnd(y);
      if(0<=xMsk&&xMsk<inNx&&0<=yMsk&&yMsk<inNy){ xIndexes(x,inNx,tnx); yIndexes(y,inNx,inNy,tny); xWeights(frac(x)); yWeights(frac(y)); out[k++]=interpolate(c); }
      else out[k++]=0.0;
      if(tr===BILINEAR){ x0+=m01+yxy; y0+=m11+yyy; } else { x0+=m01; y0+=m11; } }
    yx+=m02; yy+=m12; if(tr===BILINEAR){ yxy+=m03; yyy+=m13; } }
  return out; // Float32Array = C++ 의 (float)interpolate() 를 정확히 담음
}

/* ===================== 공개 API ===================== */
function cropData(p){ const w=p.w, h=p.h, cw=w-1, ch=h-1, d=p.data, o=new Float64Array(cw*ch);
  for(let y=0;y<ch;y++) for(let x=0;x<cw;x++) o[y*cw+x]=d[y*w+x]; return o; }
function checkPlane(p,name){ if(!p||!(p.w>=2)||!(p.h>=2)||!p.data||p.data.length!==p.w*p.h) throw new Error("StackReg: "+name+" 크기가 잘못됨"); }

/* StackReg.register(ref, mov): pystackreg 처럼 ref[:-1,:-1], mov[:-1,:-1] 로 잘라 정합 */
IC.turboReg=function(refPlane,movPlane,mode){
  const tr=modeCode(mode); checkPlane(refPlane,"ref"); checkPlane(movPlane,"mov");
  if(refPlane.w!==movPlane.w||refPlane.h!==movPlane.h) throw new Error("StackReg: 두 이미지 크기가 다름");
  if(refPlane.w<16||refPlane.h<16) throw new Error("StackReg: 이미지가 너무 작습니다 (가로·세로 16 px 이상)");
  for(const q of [refPlane,movPlane]){ const d=q.data; for(let i=0;i<d.length;i++){ const v=d[i]; if(v!==v||v===Infinity||v===-Infinity) throw new Error("StackReg: 이미지에 NaN·무한대 화소가 있습니다 — 값으로 바꾼 뒤 정렬하세요"); } }
  const r=registerRaw(cropData(refPlane),cropData(movPlane),refPlane.w-1,refPlane.h-1,tr);
  return {mat:shortToLong(r.short,tr),refPts:r.refPts,movPts:r.movPts};
};
/* StackReg.transform(mov, tmat): 같은 크기의 이미지 (Float32Array) */
IC.turboTransform=function(movPlane,mat,mode){
  const tr=modeCode(mode); checkPlane(movPlane,"mov");
  return transformRaw(movPlane.data,movPlane.w,movPlane.h,longToShort(mat,tr));
};
function matmul(A,B){ const n=A.length, C=[]; for(let i=0;i<n;i++){ C.push([]); for(let j=0;j<n;j++){ let s=0; for(let k=0;k<n;k++) s+=A[i][k]*B[k][j]; C[i].push(s); } } return C; }
const ident=n=>[...Array(n)].map((_,i)=>[...Array(n)].map((_,j)=>i===j?1:0));
/* numpy 평균: float32 배열이면 float32 누적(각 덧셈마다 반올림), 정수(8/16-bit)면 float64 누적. 축 0 방향 순차 합 후 n 으로 나눔 */
function meanPlanes(planes,i0,i1){
  const p0=planes[i0], N=p0.w*p0.h, isF32=!(p0.bit===8||p0.bit===16), s=new Float64Array(N), n=i1-i0;
  for(let i=0;i<N;i++) s[i]=planes[i0].data[i];
  for(let k=i0+1;k<i1;k++){ const d=planes[k].data; if(isF32) for(let i=0;i<N;i++) s[i]=f32(s[i]+d[i]); else for(let i=0;i<N;i++) s[i]+=d[i]; }
  if(isF32) for(let i=0;i<N;i++) s[i]=f32(s[i]/n); else for(let i=0;i<N;i++) s[i]/=n;
  return {w:p0.w,h:p0.h,bit:p0.bit,data:s};
}
/* StackReg.register_stack(img, reference, n_frames, axis=0, moving_average=1) */
IC.stackRegMatrices=function(planes,mode,reference,nFrames){
  const tr=modeCode(mode); reference=reference||"previous"; nFrames=nFrames||1;
  if(tr===BILINEAR&&reference==="previous") throw new Error("StackReg: bilinear 는 reference='previous' 와 함께 쓸 수 없음 (pystackreg 와 같은 제한)");
  const N=planes.length, dim=tr===BILINEAR?4:3, mats=[...Array(N)].map(()=>ident(dim));
  let ref=null, i0=1;
  if(reference==="first"){ ref=meanPlanes(planes,0,Math.min(nFrames,N)); if(nFrames>1) i0=0; } // pystackreg 는 n_frames>1 이어도 첫 장을 맞추지 않음(첫 장이 기준과 다른데 그대로) — 여기서는 첫 장도 평균 기준에 맞춘다 (SPEC §3c)
  else if(reference==="mean"){ ref=meanPlanes(planes,0,N); i0=0; }
  else if(reference!=="previous") throw new Error('StackReg: 알 수 없는 reference "'+reference+'"');
  for(let i=i0;i<N;i++){
    if(reference==="previous") ref=planes[i-1];
    mats[i]=IC.turboReg(ref,planes[i],tr).mat;
    if(reference==="previous"&&i>0) mats[i]=matmul(mats[i],mats[i-1]);
    if(IC.progress) IC.progress((i-i0+1)/(N-i0));
  }
  return mats;
};
/* StackReg.transform_stack(img, tmats) */
IC.stackRegApply=function(planes,mats,mode){
  const tr=modeCode(mode); if(mats.length!==planes.length) throw new Error("StackReg: 행렬 수와 장 수가 다름");
  const out=[]; for(let i=0;i<planes.length;i++){ out.push(IC.turboTransform(planes[i],mats[i],tr)); if(IC.progress) IC.progress((i+1)/planes.length); }
  return out;
};
IC._stackreg={TRANSLATION,RIGID_BODY,SCALED_ROTATION,AFFINE,BILINEAR,getPyramidDepth,registerRaw,transformRaw,shortToLong,longToShort,makeImage};
})();
