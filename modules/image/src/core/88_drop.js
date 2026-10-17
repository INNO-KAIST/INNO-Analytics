/* ============================================================
   Drop Shape — 접촉각(sessile drop)·펜던트 드랍(pendant drop) 계산 엔진  (IMG_CORE.drop*)
   - 윤곽: 가우스 흐림 → Otsu 문턱 → 가장 큰(또는 지정한) 어두운 성분 → marching squares 부화소 등고선
   - Young–Laplace(무차원, apex 곡률 반지름 b 로 나눔): dx/ds = cos φ, dz/ds = sin φ,
       세실(apex 위, z 아래로) dφ/ds = 2 + Bo·z − sin φ/x · 펜던트(apex 아래, z 위로) dφ/ds = 2 − Bo·z − sin φ/x,  Bo = Δρ g b²/γ
     (Rotenberg, Boyce, Neumann, J. Colloid Interface Sci. 93 (1983) 169 — ADSA; Berry et al., J. Colloid Interface Sci. 454 (2015) 226)
   - 맞춤: 윤곽 점 ↔ 이론 곡선 최단 거리 제곱합, Levenberg–Marquardt(수치 야코비안)
   - 접촉각 방법: Young–Laplace · 원(구면 모자) · 타원(Fitzgibbon 1999 직접 맞춤) · 접선(접촉점 근처 2 차 다항식)
   정답 자료: ref/gen_drop.py (scipy DOP853, 이 파일과 다른 구현) · 실측 OpenDrop 예제(test/drop_real, 저장소 제외)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE||(G.IMG_CORE={});
const D2R=Math.PI/180;
const amax=(A,f)=>{ let m=-Infinity; for(const q of A){ const v=f(q); if(v>m) m=v; } return m; }, amin=(A,f)=>{ let m=Infinity; for(const q of A){ const v=f(q); if(v<m) m=v; } return m; }; // 큰 배열에서 Math.max(...) 는 RangeError

/* ---------- Young–Laplace 적분 (RK4, 원점은 급수로 시작) ----------
   mode "sessile"|"pendant". 멈춤: z > zMax 또는 φ ≥ phiMax 또는 x ≤ 0 (돌아옴) 또는 s > sMax.
   반환 {s,x,z,p} (Float64Array, 무차원) */
IC.dropYL=function(Bo,mode,o){ o=o||{}; const sg=mode==="pendant"?-1:1; const h=o.ds||0.002; const zMax=o.zMax==null?Infinity:o.zMax, phiMax=o.phiMax==null?Math.PI*0.999:o.phiMax, sMax=o.sMax||40;
  const f=(x,z,p,out)=>{ const sx=x>1e-9?Math.sin(p)/x:1; out[0]=Math.cos(p); out[1]=Math.sin(p); out[2]=2+sg*Bo*z-sx; };
  const N=Math.ceil(sMax/h)+2; const S=new Float64Array(N), X=new Float64Array(N), Z=new Float64Array(N), P=new Float64Array(N);
  // 급수: x = s − s³/6·(…), 작은 s0 에서 x≈s, z≈s²/2, φ≈s
  let s=1e-5, x=s, z=s*s/2, p=s; let n=0; S[n]=0; X[n]=0; Z[n]=0; P[n]=0; n++;
  const k1=[0,0,0],k2=[0,0,0],k3=[0,0,0],k4=[0,0,0];
  S[n]=s; X[n]=x; Z[n]=z; P[n]=p; n++;
  while(n<N){ f(x,z,p,k1); f(x+h/2*k1[0],z+h/2*k1[1],p+h/2*k1[2],k2); f(x+h/2*k2[0],z+h/2*k2[1],p+h/2*k2[2],k3); f(x+h*k3[0],z+h*k3[1],p+h*k3[2],k4);
    const xn=x+h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), zn=z+h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]), pn=p+h/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
    if(!(xn>0)||!isFinite(xn)||!isFinite(zn)) break; s+=h; x=xn; z=zn; p=pn; S[n]=s; X[n]=x; Z[n]=z; P[n]=p; n++;
    if(z>zMax||p>=phiMax||s>=sMax) break; }
  return {s:S.subarray(0,n),x:X.subarray(0,n),z:Z.subarray(0,n),p:P.subarray(0,n),n}; };

/* 점 (r,z) 와 곡선(z 가 s 에 따라 단조 증가하는 부분)의 최단 거리 — z 로 이분 탐색 후 이웃 선분 검사 */
function curveDist(C,r,z,win){ const Z=C.z, X=C.x, n=C.n; let lo=0, hi=n-1; if(n<2) return Infinity;
  while(hi-lo>1){ const m=(lo+hi)>>1; if(Z[m]<z) lo=m; else hi=m; }
  let best=Infinity; const a=Math.max(0,lo-win), b=Math.min(n-2,lo+win);
  for(let i=a;i<=b;i++){ const x0=X[i], z0=Z[i], dx=X[i+1]-x0, dz=Z[i+1]-z0; const L2=dx*dx+dz*dz; let t=L2>0?((r-x0)*dx+(z-z0)*dz)/L2:0; t=t<0?0:t>1?1:t; const ex=x0+t*dx-r, ez=z0+t*dz-z; const d=ex*ex+ez*ez; if(d<best) best=d; }
  if(z>Z[n-1]){ const ex=X[n-1]-r, ez=Z[n-1]-z; best=Math.min(best,ex*ex+ez*ez); }
  return Math.sqrt(best); }
IC._dropCurveDist=curveDist;
/* z 가 단조 증가하는 앞부분만 (φ < π) */
function monoPrefix(C){ let n=C.n; for(let i=1;i<C.n;i++){ if(C.z[i]<C.z[i-1]){ n=i; break; } } return {x:C.x,z:C.z,p:C.p,s:C.s,n}; }

/* ---------- Otsu 문턱 ---------- */
function otsu(vals){ let mn=Infinity,mx=-Infinity; for(const v of vals){ if(v<mn) mn=v; if(v>mx) mx=v; } if(!(mx>mn)) return mn; const B=256, h=new Float64Array(B);
  for(const v of vals) h[Math.min(B-1,Math.floor((v-mn)/(mx-mn)*B))]++; const N=vals.length; let sum=0; for(let i=0;i<B;i++) sum+=i*h[i];
  let w0=0,s0=0,best=-1,bi=0; for(let i=0;i<B;i++){ w0+=h[i]; if(!w0) continue; const w1=N-w0; if(!w1) break; s0+=i*h[i]; const m0=s0/w0, m1=(sum-s0)/w1; const v=w0*w1*(m0-m1)*(m0-m1); if(v>best){ best=v; bi=i; } }
  return mn+(bi+1)/B*(mx-mn); }
function gauss1(sig){ const r=Math.max(1,Math.ceil(3*sig)); const k=new Float64Array(2*r+1); let s=0; for(let i=-r;i<=r;i++){ k[i+r]=Math.exp(-i*i/(2*sig*sig)); s+=k[i+r]; } for(let i=0;i<k.length;i++) k[i]/=s; return {k,r}; }
function blur(src,w,h,sig){ if(!(sig>0)) return Float64Array.from(src); const {k,r}=gauss1(sig); const t=new Float64Array(w*h), o=new Float64Array(w*h);
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ let s=0; for(let i=-r;i<=r;i++){ const xx=x+i<0?0:x+i>=w?w-1:x+i; s+=k[i+r]*src[y*w+xx]; } t[y*w+x]=s; }
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ let s=0; for(let i=-r;i<=r;i++){ const yy=y+i<0?0:y+i>=h?h-1:y+i; s+=k[i+r]*t[yy*w+x]; } o[y*w+x]=s; }
  return o; }

/* ---------- 윤곽 ----------
   plane {w,h,data}, o.roi {x0,y0,x1,y1}(포함, 없으면 전체), o.sigma(1), o.dark(자동: 가운데가 어두운지), o.pick "largest"|"top"|"bottom"
   반환 {pts:[[x,y]…] (영상 좌표, 닫힌 고리 순서), level, dark, comp:{area,x0,x1,y0,y1}, mask(roi 크기 Uint8Array), roi} */
IC.dropContour=function(plane,o){ o=o||{}; const W=plane.w, H=plane.h, d=plane.data; const R=o.roi||{x0:0,y0:0,x1:W-1,y1:H-1};
  const x0=Math.max(0,Math.floor(R.x0)), y0=Math.max(0,Math.floor(R.y0)), x1=Math.min(W-1,Math.ceil(R.x1)), y1=Math.min(H-1,Math.ceil(R.y1)); const w=x1-x0+1, h=y1-y0+1;
  if(w<8||h<8) throw new Error("영역이 너무 작습니다");
  const src=new Float64Array(w*h); let nNaN=0; for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const v=d[(y+y0)*W+x+x0]; if(v===v) src[y*w+x]=v; else { src[y*w+x]=NaN; nNaN++; } }
  if(nNaN){ const fin=[]; for(let i=0;i<src.length;i+=Math.max(1,src.length>>14)) if(src[i]===src[i]) fin.push(src[i]); fin.sort((a,b)=>a-b); const med=fin.length?fin[fin.length>>1]:0; for(let i=0;i<src.length;i++) if(src[i]!==src[i]) src[i]=med; } // NaN = 영역 중앙값(바탕 쪽) — 0 이면 가짜 어두운 점
  const sm=blur(src,w,h,o.sigma==null?1:o.sigma); const step=Math.max(1,Math.floor(Math.sqrt(w*h/200000))); const samp=[]; for(let i=0;i<sm.length;i+=step) samp.push(sm[i]);
  const level=o.level!=null?o.level:otsu(samp);
  // 방울이 어두운가: 영역 테두리(대부분 바탕 — 바늘·기판은 일부)의 중앙값이 문턱보다 밝으면 방울 = 어두움
  let dark=o.dark; if(dark==null){ const c=[]; for(let x=0;x<w;x+=2){ c.push(sm[x]); c.push(sm[(h-1)*w+x]); } for(let y=0;y<h;y+=2){ c.push(sm[y*w]); c.push(sm[y*w+w-1]); } c.sort((a,b)=>a-b); dark=c[c.length>>1]>level; }
  const inside=i=>dark?sm[i]<level:sm[i]>level;
  // 4-연결 성분
  const lab=new Int32Array(w*h).fill(-1), st=new Int32Array(w*h); const comps=[];
  for(let i0=0;i0<w*h;i0++){ if(lab[i0]>=0||!inside(i0)) continue; const id=comps.length; let sp=0; st[sp++]=i0; lab[i0]=id; let n=0,cx0=w,cx1=0,cy0=h,cy1=0;
    while(sp){ const p=st[--sp]; n++; const x=p%w, y=(p-x)/w; if(x<cx0) cx0=x; if(x>cx1) cx1=x; if(y<cy0) cy0=y; if(y>cy1) cy1=y;
      if(x>0&&lab[p-1]<0&&inside(p-1)){ lab[p-1]=id; st[sp++]=p-1; } if(x<w-1&&lab[p+1]<0&&inside(p+1)){ lab[p+1]=id; st[sp++]=p+1; }
      if(y>0&&lab[p-w]<0&&inside(p-w)){ lab[p-w]=id; st[sp++]=p-w; } if(y<h-1&&lab[p+w]<0&&inside(p+w)){ lab[p+w]=id; st[sp++]=p+w; } }
    comps.push({id,area:n,x0:cx0,x1:cx1,y0:cy0,y1:cy1}); }
  if(!comps.length) throw new Error("방울을 찾지 못했습니다 (문턱 "+level.toFixed(1)+")");
  let pick=comps[0]; for(const c of comps){ const sc=c.area*(o.pick==="top"&&c.y0===0?4:1)*(o.pick==="bottom"&&c.y1===h-1?4:1), sb=pick.area*(o.pick==="top"&&pick.y0===0?4:1)*(o.pick==="bottom"&&pick.y1===h-1?4:1); if(sc>sb) pick=c; }
  // 고른 성분의 구멍 메우기(바깥 테두리만) — 밝은 반사점이 방울 안에 있어도 윤곽이 하나가 되게
  const m=new Uint8Array(w*h); for(let i=0;i<w*h;i++) if(lab[i]===pick.id) m[i]=1;
  { const ext=new Uint8Array(w*h); let sp=0; const push=i=>{ if(!m[i]&&!ext[i]){ ext[i]=1; st[sp++]=i; } };
    for(let x=0;x<w;x++){ push(x); push((h-1)*w+x); } for(let y=0;y<h;y++){ push(y*w); push(y*w+w-1); }
    while(sp){ const p=st[--sp]; const x=p%w, y=(p-x)/w; if(x>0) push(p-1); if(x<w-1) push(p+1); if(y>0) push(p-w); if(y<h-1) push(p+w); }
    for(let i=0;i<w*h;i++) if(!ext[i]) m[i]=1; }
  // 성분 밖은 바탕 값, 성분 안은 원래(흐린) 값 — 성분 경계에서 문턱을 지나는 곳 = 윤곽. 구멍(메운 곳)은 안쪽 값으로
  const far=dark?Math.max(level+1,level+Math.abs(level)*0.5+1):Math.min(level-1,level-Math.abs(level)*0.5-1), inv=dark?-Infinity:Infinity;
  const F=new Float64Array(w*h); for(let i=0;i<w*h;i++){ if(!m[i]){ // 성분 둘레 1 px 은 실제 값(부화소 위치 보존), 더 먼 곳은 바탕
      const x=i%w, y=(i-x)/w; let near=false; for(let dy=-1;dy<=1&&!near;dy++) for(let dx=-1;dx<=1;dx++){ const xx=x+dx, yy=y+dy; if(xx>=0&&yy>=0&&xx<w&&yy<h&&m[yy*w+xx]){ near=true; break; } }
      F[i]=near?(dark?Math.max(sm[i],level+1e-9):Math.min(sm[i],level-1e-9)):far; }
    else F[i]=inside(i)?sm[i]:(dark?level-1e-6:level+1e-6); }
  // 영역 테두리에 닿은 성분: 바깥(영역 밖)을 바탕으로 보아 닫힌 고리를 만든다 — 한 줄 덧댐
  const pw=w+2, ph=h+2, Fp=new Float64Array(pw*ph).fill(far); for(let y=0;y<h;y++) for(let x=0;x<w;x++) Fp[(y+1)*pw+x+1]=F[y*w+x];
  const pts=marchingLoop(Fp,pw,ph,level,dark); if(!pts.length) throw new Error("윤곽을 만들지 못했습니다");
  // 덧댄 줄 위의 점(영역 테두리를 따라 그은 가짜 윤곽)은 표시 — 맞춤에서 빼야 한다
  const edge=pts.map(([x,y])=>x<0.75||y<0.75||x>pw-1.75||y>ph-1.75);
  // 경사 최대 자리로 다듬기(Canny 처럼): 바탕 밝기가 고르지 않으면(비네팅) 한 문턱의 등고선은 가장자리에서 밀린다 → 법선 ±3 px 에서 |∂I/∂n| 최대(포물선 보간)
  const loc=pts.map(([x,y])=>[x-1,y-1]); if(o.refine!==false){ const S=(x,y)=>{ x=x<0?0:x>w-1.001?w-1.001:x; y=y<0?0:y>h-1.001?h-1.001:y; const i=Math.floor(x), j=Math.floor(y), fx=x-i, fy=y-j, k=j*w+i; return sm[k]*(1-fx)*(1-fy)+sm[k+1]*fx*(1-fy)+sm[k+w]*(1-fx)*fy+sm[k+w+1]*fx*fy; };
    for(let q=0;q<loc.length;q++){ if(edge[q]) continue; const [x,y]=loc[q]; if(x<4||y<4||x>w-5||y>h-5) continue; const gx=S(x+1,y)-S(x-1,y), gy=S(x,y+1)-S(x,y-1), gn=Math.hypot(gx,gy); if(!(gn>0)) continue; const nx=gx/gn, ny=gy/gn;
      const T=[]; for(let t=-3;t<=3.001;t+=0.5) T.push(Math.abs(S(x+(t+1)*nx,y+(t+1)*ny)-S(x+(t-1)*nx,y+(t-1)*ny)));
      let bi=0; for(let i=1;i<T.length;i++) if(T[i]>T[bi]) bi=i; if(bi===0||bi===T.length-1) continue; const den=T[bi-1]-2*T[bi]+T[bi+1]; const d=den<0?0.5*(T[bi-1]-T[bi+1])/den:0; const t=-3+(bi+d)*0.5;
      if(Math.abs(t)<=2.5){ loc[q]=[x+t*nx,y+t*ny]; } } }
  const out=loc.map(([x,y])=>[x+x0,y+y0]), raw=pts.map(([x,y])=>[x-1+x0,y-1+y0]);
  return {nNaN,pts:out,raw,edge,level,dark,comp:{area:pick.area,x0:pick.x0+x0,x1:pick.x1+x0,y0:pick.y0+y0,y1:pick.y1+y0},mask:m,roi:{x0,y0,x1,y1}}; };

/* marching squares: level 등고선 중 가장 긴 닫힌 고리 (점 = 선형 보간 교차점) */
function marchingLoop(F,w,h,lv,dark){ const ins=v=>dark?v<lv:v>lv; const segs=new Map(); const key=(x,y,e)=>((y*w+x)*4+e);
  // 칸 (x,y)~(x+1,y+1) 의 네 변: 0 위, 1 오른, 2 아래, 3 왼 — 변 id 는 공유되게 (위 = (x,y) 가로, 왼 = (x,y) 세로)
  const eid=(x,y,e)=>e===0?((y*w+x)*2):e===2?(((y+1)*w+x)*2):e===3?((y*w+x)*2+1):((y*w+x+1)*2+1);
  const pt=(x,y,e)=>{ let ax,ay,bx,by; if(e===0){ ax=x;ay=y;bx=x+1;by=y; } else if(e===1){ ax=x+1;ay=y;bx=x+1;by=y+1; } else if(e===2){ ax=x;ay=y+1;bx=x+1;by=y+1; } else { ax=x;ay=y;bx=x;by=y+1; }
    const va=F[ay*w+ax], vb=F[by*w+bx]; const t=(lv-va)/(vb-va); const tt=isFinite(t)?Math.min(1,Math.max(0,t)):0.5; return [ax+tt*(bx-ax),ay+tt*(by-ay)]; };
  const adj=new Map(), P=new Map(); const link=(a,b)=>{ (adj.get(a)||adj.set(a,[]).get(a)).push(b); (adj.get(b)||adj.set(b,[]).get(b)).push(a); };
  for(let y=0;y<h-1;y++) for(let x=0;x<w-1;x++){ const a=ins(F[y*w+x])?1:0, b=ins(F[y*w+x+1])?1:0, c=ins(F[(y+1)*w+x+1])?1:0, d=ins(F[(y+1)*w+x])?1:0; const k=a|(b<<1)|(c<<2)|(d<<3);
    if(k===0||k===15) continue; const E=[]; // 변에 교차가 있는 것
    if(a!==b) E.push(0); if(b!==c) E.push(1); if(c!==d) E.push(2); if(d!==a) E.push(3);
    for(const e of E){ const id=eid(x,y,e); if(!P.has(id)) P.set(id,pt(x,y,e)); }
    if(E.length===2) link(eid(x,y,E[0]),eid(x,y,E[1]));
    else if(E.length===4){ // 안장점: 칸 평균으로 연결 방식 결정
      const mid=ins((F[y*w+x]+F[y*w+x+1]+F[(y+1)*w+x+1]+F[(y+1)*w+x])/4); const e=E.map(q=>eid(x,y,q));
      if(mid===!!a){ link(e[0],e[1]); link(e[2],e[3]); } else { link(e[0],e[3]); link(e[1],e[2]); } } }
  const seen=new Set(); let best=[];
  for(const s0 of adj.keys()){ if(seen.has(s0)) continue; const loop=[]; let prev=-1, cur=s0; while(cur!=null&&!seen.has(cur)){ seen.add(cur); loop.push(P.get(cur)); const nb=adj.get(cur)||[]; const nx=nb.find(q=>q!==prev&&!seen.has(q)); prev=cur; cur=nx; }
    if(loop.length>best.length) best=loop; }
  return best; }

/* ---------- 작은 선형대수 ---------- */
function solveSym(A,b){ const n=b.length; const M=A.map((r,i)=>r.concat([b[i]]));
  for(let c=0;c<n;c++){ let p=c; for(let r=c+1;r<n;r++) if(Math.abs(M[r][c])>Math.abs(M[p][c])) p=r; [M[c],M[p]]=[M[p],M[c]]; const v=M[c][c]; if(Math.abs(v)<1e-300) return null;
    for(let r=0;r<n;r++){ if(r===c) continue; const f=M[r][c]/v; if(!f) continue; for(let k=c;k<=n;k++) M[r][k]-=f*M[c][k]; } }
  return M.map((r,i)=>r[n]/r[i]); }
function invSym(A){ const n=A.length; const out=[]; for(let j=0;j<n;j++){ const e=new Array(n).fill(0); e[j]=1; const x=solveSym(A,e); if(!x) return null; out.push(x); } return out[0].map((_,i)=>out.map(r=>r[i])); }
/* Levenberg–Marquardt: fn(p) → 잔차 배열, 수치 야코비안(상대 단계 h) */
function lm(fn,p0,o){ o=o||{}; let p=p0.slice(); let r=fn(p); if(!r) return null; let cost=r.reduce((a,v)=>a+v*v,0); let lam=o.lam0||1e-3; const n=p.length, it=o.iter||60, hs=o.h||p.map(v=>1e-4*Math.max(1,Math.abs(v)));
  let J=null, JtJ=null, Jtr=null;
  for(let k=0;k<it;k++){ J=[]; for(let j=0;j<n;j++){ const q=p.slice(); q[j]+=hs[j]; const rq=fn(q); if(!rq) return null; J.push(rq.map((v,i)=>(v-r[i])/hs[j])); }
    JtJ=[...Array(n)].map((_,a)=>[...Array(n)].map((_,b)=>{ let s=0; const A=J[a], B=J[b]; for(let i=0;i<A.length;i++) s+=A[i]*B[i]; return s; }));
    Jtr=[...Array(n)].map((_,a)=>{ let s=0; const A=J[a]; for(let i=0;i<A.length;i++) s+=A[i]*r[i]; return s; });
    let improved=false;
    for(let t=0;t<12;t++){ const M=JtJ.map((row,i)=>row.map((v,j)=>i===j?v*(1+lam)+1e-12:v)); const dp=solveSym(M,Jtr.map(v=>-v)); if(!dp){ lam*=10; continue; }
      const q=p.map((v,i)=>v+dp[i]); if(o.clamp) o.clamp(q); const rq=fn(q); if(!rq){ lam*=10; continue; } const c2=rq.reduce((a,v)=>a+v*v,0);
      if(c2<cost){ const rel=(cost-c2)/Math.max(cost,1e-300); p=q; r=rq; cost=c2; lam=Math.max(lam/3,1e-9); improved=true; if(rel<1e-10) k=it; break; } lam*=4; }
    if(!improved) break; }
  const m=r.length; const s2=cost/Math.max(1,m-n); const C=JtJ?invSym(JtJ):null; const sd=C?C.map((row,i)=>Math.sqrt(Math.max(0,row[i]*s2))):p.map(()=>NaN);
  return {p,cost,rms:Math.sqrt(cost/Math.max(1,m)),sd,cov:C?C.map(row=>row.map(v=>v*s2)):null,m}; }
IC._lm=lm;

/* ---------- 바늘 (펜던트) ----------
   윤곽 점에서 위쪽(바늘) 부분: 줄마다 왼·오른 가장자리 → 위에서부터 폭이 일정한 구간 = 바늘, 두 가장자리에 직선 맞춤
   반환 {diam(px), axisDeg(세로에서 기운 각, 시계 +), yTip, left:{a,b}, right:{a,b}} (x = a + b·y) */
IC.dropNeedle=function(pts,edge,o){ o=o||{}; const P=pts.filter((_,i)=>!(edge&&edge[i])); if(P.length<20) throw new Error("윤곽 점이 너무 적습니다");
  let ymin=Infinity, ymax=-Infinity; for(const [,y] of P){ if(y<ymin) ymin=y; if(y>ymax) ymax=y; }
  const rows=new Map(); for(const [x,y] of P){ const k=Math.round(y); const r=rows.get(k)||{l:Infinity,r:-Infinity}; if(x<r.l) r.l=x; if(x>r.r) r.r=x; rows.set(k,r); }
  const ks=[...rows.keys()].sort((a,b)=>a-b); const top=ks.slice(0,Math.max(5,Math.floor(ks.length*0.08))); const wTop=top.map(k=>rows.get(k).r-rows.get(k).l).sort((a,b)=>a-b); const w0=wTop[wTop.length>>1];
  if(!(w0>2)) throw new Error("바늘을 찾지 못했습니다 — 위쪽에 바늘이 보이게 영역을 잡으세요");
  let yTip=ks[0]; for(const k of ks){ const r=rows.get(k); if(Math.abs((r.r-r.l)-w0)>Math.max(1.5,0.03*w0)) break; yTip=k; }
  const span=yTip-ks[0]; if(span<5) throw new Error("바늘 구간이 너무 짧습니다 ("+span+" px) — 영역 위쪽에 바늘을 더 넣으세요");
  const y1=ks[0]+0.1*span, y2=yTip-0.1*span; const L=[],Rr=[]; for(const [x,y] of P){ if(y<y1||y>y2) continue; const r=rows.get(Math.round(y)); if(!r) continue; if(Math.abs(x-r.l)<1.5) L.push([x,y]); else if(Math.abs(x-r.r)<1.5) Rr.push([x,y]); }
  const fit=Q=>{ let sy=0,sx=0,syy=0,sxy=0; const n=Q.length; for(const [x,y] of Q){ sy+=y; sx+=x; syy+=y*y; sxy+=x*y; } const den=n*syy-sy*sy; const b=den?(n*sxy-sx*sy)/den:0; return {a:(sx-b*sy)/n,b,n}; };
  if(L.length<4||Rr.length<4) throw new Error("바늘 가장자리를 맞추지 못했습니다");
  const l=fit(L), r=fit(Rr); const bm=(l.b+r.b)/2; const ym=(y1+y2)/2; const diam=Math.abs((r.a+r.b*ym)-(l.a+l.b*ym))*Math.cos(Math.atan(bm));
  return {diam,axisDeg:Math.atan(bm)/D2R,yTip,yTop:ks[0],left:l,right:r,cx:((l.a+l.b*ym)+(r.a+r.b*ym))/2,ym}; };

/* ---------- 펜던트 드랍 맞춤 ----------
   pts: 방울 윤곽(바늘 아래), o: {rot0(°), maxPts}. 매개변수 [x0, y0(apex), b(px), Bo, ω(rad)] — 영상 y 아래로, apex 가 가장 아래
   반환 {apex:[x,y], b, Bo, omega(°), rms, sd:{…}, curve:{x,z}(무차원), zTop(무차원), n} */
IC.fitPendant=function(pts,o){ o=o||{}; let Q=pts; const maxP=o.maxPts||900; if(Q.length>maxP){ const st=Q.length/maxP; Q=[...Array(maxP)].map((_,i)=>Q[Math.floor(i*st)]); }
  // 초기값: apex = 가장 아래 점 근처, 축 = rot0, b = 아래 10 % 점에 원 맞춤, Bo = 거친 격자
  const om0=(o.rot0||0)*D2R; let yb=-Infinity, xb=0; for(const [x,y] of Q) if(y>yb){ yb=y; xb=x; }
  let ytop=Infinity; for(const [,y] of Q) if(y<ytop) ytop=y; const Hd=yb-ytop;
  const low=Q.filter(([,y])=>y>yb-0.12*Hd); const circ=circleFit(low); const b0=circ&&circ.r>0?circ.r:Hd/2; const x00=circ?circ.cx:xb, y00=circ?circ.cy+circ.r:yb;
  const zMaxOf=p=>{ let zm=0; const c=Math.cos(p[4]), s=Math.sin(p[4]); for(const [x,y] of Q){ const dx=x-p[0], dy=y-p[1]; const v=-dx*s+dy*c; if(-v>zm) zm=-v; } return zm/p[2]*1.03+0.05; };
  const cache=new Map(); const curveOf=(Bo,zm)=>{ const k=Bo.toPrecision(10)+"|"+zm.toFixed(2); let C=cache.get(k); if(!C){ C=monoPrefix(IC.dropYL(Bo,"pendant",{zMax:zm,ds:0.003,phiMax:Math.PI*0.999})); if(cache.size>64) cache.clear(); cache.set(k,C); } return C; };
  const res=p=>{ if(!(p[2]>2)||!(p[3]>-0.05&&p[3]<1.5)) return null; const zm=Math.ceil(zMaxOf(p)*4)/4; const C=curveOf(p[3],zm); if(C.n<10||C.z[C.n-1]<zm*0.97) return null; const c=Math.cos(p[4]), s=Math.sin(p[4]), b=p[2];
    return Q.map(([x,y])=>{ const dx=x-p[0], dy=y-p[1]; const u=dx*c+dy*s, v=-dx*s+dy*c; return b*curveDist(C,Math.abs(u)/b,-v/b,40); }); };
  let best=null; for(const Bo of [0.02,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]){ const p=[x00,y00,b0,Bo,om0]; const r=res(p); if(!r) continue; const c=r.reduce((a,v)=>a+v*v,0); if(!best||c<best.c) best={c,p}; }
  if(!best) throw new Error("초기값을 찾지 못했습니다");
  const F=lm(res,best.p,{iter:80,h:[0.02,0.02,0.01,1e-4,1e-5],clamp:q=>{ if(q[3]<-0.04) q[3]=-0.04; if(q[3]>1.4) q[3]=1.4; if(q[2]<3) q[2]=3; }});
  const atBound=F&&(F.p[3]<=-0.0399||F.p[3]>=1.399);
  if(!F) throw new Error("맞춤 실패");
  const p=F.p; const zTop=zMaxOf(p)/1.03; const C=IC.dropYL(p[3],"pendant",{zMax:zTop+0.5,ds:0.002}); const Cm=monoPrefix(C);
  return {apex:[p[0],p[1]],b:p[2],Bo:p[3],omega:p[4]/D2R,rms:F.rms,atBound,sd:{x0:F.sd[0],y0:F.sd[1],b:F.sd[2],Bo:F.sd[3],omega:F.sd[4]/D2R},cov:F.cov,curve:{x:Array.from(Cm.x.subarray(0,Cm.n)),z:Array.from(Cm.z.subarray(0,Cm.n))},n:Q.length}; };

/* 회전체 부피·겉넓이 (무차원 곡선, z ≤ zTop 까지) */
IC.dropVolArea=function(curve,zTop){ const X=curve.x, Z=curve.z; let V=0,A=0; for(let i=1;i<X.length;i++){ let x0=X[i-1],z0=Z[i-1],x1=X[i],z1=Z[i]; if(z0>=zTop) break; if(z1>zTop){ const t=(zTop-z0)/(z1-z0); x1=x0+t*(x1-x0); z1=zTop; }
    const xm=(x0+x1)/2; V+=Math.PI*xm*xm*(z1-z0); A+=2*Math.PI*xm*Math.hypot(x1-x0,z1-z0); } return {V,A}; };

/* 곡선에서 높이 zt 에서의 반지름 x (z 단조 부분 선형 보간) */
function xAtZ(C,zt){ const Z=C.z, X=C.x; for(let i=1;i<Z.length;i++) if(Z[i]>=zt){ const t=(zt-Z[i-1])/(Z[i]-Z[i-1]); return X[i-1]+t*(X[i]-X[i-1]); } return NaN; }
function phiAtZ(C,zt){ const Z=C.z, P=C.p; for(let i=1;i<C.n;i++) if(Z[i]>=zt){ const t=(zt-Z[i-1])/(Z[i]-Z[i-1]); return P[i-1]+t*(P[i]-P[i-1]); } return NaN; }

/* ---------- 펜던트 전체: 윤곽 → 바늘 → 맞춤 → 물성 ----------
   o: {roi, sigma, needleMM (바늘 지름 mm, 없으면 pxSize 사용), pxSize (mm/px, 선택), rhoIn, rhoOut (kg/m³), g (m/s²)} */
IC.pendantDrop=function(plane,o){ o=o||{}; const ct=IC.dropContour(plane,{roi:o.roi,sigma:o.sigma,pick:"top"});
  const nd=IC.dropNeedle(ct.pts,ct.edge); const tipMargin=Math.max(2,0.02*(ct.comp.y1-nd.yTip));
  const drop=ct.pts.filter((q,i)=>!ct.edge[i]&&q[1]>nd.yTip+tipMargin);
  const fit=IC.fitPendant(drop,{rot0:-nd.axisDeg});
  let mmPerPx=null, scaleFrom=""; if(o.needleMM>0){ mmPerPx=o.needleMM/nd.diam; scaleFrom="needle"; } else if(o.pxSize>0){ mmPerPx=o.pxSize; scaleFrom="image"; }
  const drho=(o.rhoIn==null?1000:+o.rhoIn)-(o.rhoOut==null?0:+o.rhoOut), g=o.g||9.80665;
  // 바늘 끝 높이(무차원): apex 에서 바늘 끝까지
  const c=Math.cos(fit.omega*D2R), s=Math.sin(fit.omega*D2R); const tipPt=[nd.cx,nd.yTip]; const vTip=-(tipPt[0]-fit.apex[0])*s+(tipPt[1]-fit.apex[1])*c; let zTip=-vTip/fit.b;
  // 바늘이 붙는 높이: 적도 위에서 이론 곡선 반지름 = 바늘 반지름인 곳 (폭 문턱으로 찾은 yTip 보다 정확 — 가장자리 흐림·문턱에 덜 민감). 곡선이 바늘보다 늘 넓으면 yTip 그대로
  { const X=fit.curve.x, Z=fit.curve.z, rN=nd.diam/2/fit.b; let ie=0; for(let i=1;i<X.length;i++){ if(X[i]>X[ie]) ie=i; else if(X[i]<X[ie]*0.999) break; }
    for(let i=ie+1;i<X.length;i++){ if(X[i]<=rN){ const t=(X[i-1]-rN)/(X[i-1]-X[i]); const zz=Z[i-1]+t*(Z[i]-Z[i-1]); if(Math.abs(zz-zTip)<0.15*zTip) zTip=zz; break; } } }
  const VA=IC.dropVolArea(fit.curve,zTip); const out={contour:ct.pts,edge:ct.edge,level:ct.level,needle:nd,fit,zTip,R0px:fit.b,Bo:fit.Bo,rms:fit.rms,scaleFrom,mmPerPx};
  out.xTip=xAtZ(fit.curve,zTip); // 바늘 끝에서 이론 곡선 반지름 (바늘 반지름과 비교 — 맞으면 ≈ diam/2/b)
  out.warn=[]; if(!(drho>0)){ out.warn.push(`Δρ = ${isFinite(drho)?drho:"?"} kg/m³ — 방울이 바깥보다 무거워야 펜던트 γ 를 계산합니다 (밀도 칸 확인)`); mmPerPx=null; }
  if(mmPerPx){ const b=fit.b*mmPerPx*1e-3; // m
    out.R0mm=b*1e3; out.gamma=drho*g*b*b/fit.Bo*1e3; // mN/m
    // 맞춤 SD: γ ∝ b²/Bo, b·Bo 상관(보통 −0.9)까지 넣은 gᵀCg. 눈금 SD: 바늘 가장자리 ±0.3 px 가 지름에 주는 몫(γ ∝ 눈금²) — 따로 보인다
    const C=fit.cov; const gb=2/fit.b, gB=-1/fit.Bo; const v=C?(gb*gb*C[2][2]+gB*gB*C[3][3]+2*gb*gB*C[2][3]):NaN; out.gammaSd=out.gamma*Math.sqrt(Math.max(0,v));
    out.gammaScaleSd=scaleFrom==="needle"?out.gamma*2*(0.3*Math.SQRT2/nd.diam):NaN;
    out.V=VA.V*b*b*b*1e9; out.A=VA.A*b*b*1e6; // mm³, mm²
    const Dn=nd.diam*mmPerPx*1e-3; out.needleMM=nd.diam*mmPerPx; out.Wo=drho*g*(out.V*1e-9)/(Math.PI*out.gamma*1e-3*Dn); }
  out.drho=drho; out.g=g;
  // 화면 덧그림(영상 좌표): 맞춘 이론 윤곽(양쪽, 바늘 끝까지) · 바늘 두 가장자리 · apex
  { const c=Math.cos(fit.omega*D2R), s=Math.sin(fit.omega*D2R), b=fit.b, X=fit.curve.x, Z=fit.curve.z; const half=sg=>{ const L=[]; for(let i=0;i<X.length;i++){ if(Z[i]>zTip+1e-9) break; const u=sg*X[i]*b, v=-Z[i]*b; L.push([fit.apex[0]+u*c-v*s,fit.apex[1]+u*s+v*c]); } return L; };
    const l=half(-1).reverse(), r=half(1); const ln=(f)=>[[f.a+f.b*nd.yTop,nd.yTop],[f.a+f.b*nd.yTip,nd.yTip]];
    out.overlay={curves:[l.concat(r.slice(1))],needle:[ln(nd.left),ln(nd.right)],apex:fit.apex.slice()}; }
  // 경고: Worthington 수가 작으면(방울이 작아 중력 변형이 약함) γ 정밀도가 떨어진다 — Berry et al. 2015
  if(fit.atBound) out.warn.push(`Bond 수가 맞춤 한계(${fit.Bo.toFixed(3)})에 닿았습니다 — 바늘·영역·밝기를 확인하세요. γ 는 믿지 마세요`);
  if(out.Wo!=null&&out.Wo<0.5) out.warn.push(`Worthington 수 ${out.Wo.toFixed(2)} — 0.5 보다 작아 γ 불확도가 큽니다(방울을 바늘에서 떨어지기 직전 크기에 가깝게)`);
  if(fit.rms>1) out.warn.push(`맞춤 잔차 rms ${fit.rms.toFixed(2)} px — 윤곽에 바늘·잡티·반사가 섞였는지 보세요`);
  if(Math.abs(fit.omega)>5) out.warn.push(`기울기 ${fit.omega.toFixed(1)}° — 카메라·바늘 수직을 확인하세요`);
  if(!mmPerPx) out.warn.push("바늘 지름(mm)이나 영상 눈금이 없어 γ·부피는 계산하지 않았습니다 (Bo·모양만)");
  return out; };

/* ---------- 원 맞춤 (Kåsa → 기하 LM) ---------- */
function circleFit(P){ const n=P.length; if(n<3) return null; let mx=0,my=0; for(const [x,y] of P){ mx+=x; my+=y; } mx/=n; my/=n;
  let Suu=0,Svv=0,Suv=0,Suuu=0,Svvv=0,Suvv=0,Svuu=0; for(const [x,y] of P){ const u=x-mx, v=y-my; Suu+=u*u; Svv+=v*v; Suv+=u*v; Suuu+=u*u*u; Svvv+=v*v*v; Suvv+=u*v*v; Svuu+=v*u*u; }
  const s=solveSym([[Suu,Suv],[Suv,Svv]],[(Suuu+Suvv)/2,(Svvv+Svuu)/2]); if(!s) return null; let cx=s[0]+mx, cy=s[1]+my, r=Math.sqrt(s[0]*s[0]+s[1]*s[1]+(Suu+Svv)/n);
  const F=lm(p=>P.map(([x,y])=>Math.hypot(x-p[0],y-p[1])-p[2]),[cx,cy,r],{iter:30}); if(F){ [cx,cy,r]=F.p; return {cx,cy,r,rms:F.rms}; } return {cx,cy,r}; }
IC.circleFit=circleFit;

/* ---------- 타원(일반 원뿔) 직접 맞춤 — Fitzgibbon, Pilu, Fisher 1999 (Halir–Flusser 수치 안정판) ---------- */
function ellipseFit(P){ const n=P.length; if(n<6) return null; let mx=0,my=0,sc=0; for(const [x,y] of P){ mx+=x; my+=y; } mx/=n; my/=n; for(const [x,y] of P) sc+=Math.hypot(x-mx,y-my); sc=sc/n||1;
  const D1=[],D2=[]; for(const [x0,y0] of P){ const x=(x0-mx)/sc, y=(y0-my)/sc; D1.push([x*x,x*y,y*y]); D2.push([x,y,1]); }
  const mul=(A,B)=>{ const r=[]; for(let i=0;i<3;i++){ r.push([]); for(let j=0;j<3;j++){ let s=0; for(let k=0;k<n;k++) s+=A[k][i]*B[k][j]; r[i].push(s); } } return r; };
  const S1=mul(D1,D1), S2=mul(D1,D2), S3=mul(D2,D2); const S3i=invSym(S3); if(!S3i) return null;
  const T=[0,1,2].map(i=>[0,1,2].map(j=>-(S3i[i][0]*S2[j][0]+S3i[i][1]*S2[j][1]+S3i[i][2]*S2[j][2]))); // T = −S3⁻¹ S2ᵀ
  const M0=[0,1,2].map(i=>[0,1,2].map(j=>S1[i][j]+S2[i][0]*T[0][j]+S2[i][1]*T[1][j]+S2[i][2]*T[2][j]));
  const M=[[M0[2][0]/2,M0[2][1]/2,M0[2][2]/2],[-M0[1][0],-M0[1][1],-M0[1][2]],[M0[0][0]/2,M0[0][1]/2,M0[0][2]/2]];
  // 3×3 고유벡터: 특성다항식 실근마다 영공간 → 4ac − b² > 0 인 것
  const ev=eig3(M); let a1=null; for(const v of ev){ if(4*v[0]*v[2]-v[1]*v[1]>0){ a1=v; break; } } if(!a1) return null;
  const a2=[0,1,2].map(i=>T[i][0]*a1[0]+T[i][1]*a1[1]+T[i][2]*a1[2]);
  // 원래 좌표로: x = mx + sc·X
  let [A,B,C]=a1, [Dd,E,Fc]=a2; A/=sc*sc; B/=sc*sc; C/=sc*sc; Dd/=sc; E/=sc;
  const D_=Dd-2*A*mx-B*my, E_=E-2*C*my-B*mx, F_=Fc+A*mx*mx+B*mx*my+C*my*my-Dd*mx-E*my;
  return {A,B,C,D:D_,E:E_,F:F_}; }
function eig3(M){ // 실 고유벡터들 (특성다항식 근 → (M−λI) 의 두 행 외적)
  const a=M, tr=a[0][0]+a[1][1]+a[2][2]; const c2=-tr;
  const m2=a[0][0]*a[1][1]-a[0][1]*a[1][0]+a[0][0]*a[2][2]-a[0][2]*a[2][0]+a[1][1]*a[2][2]-a[1][2]*a[2][1];
  const det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
  const roots=cubicRoots(1,c2,m2,-det); const out=[];
  for(const l of roots){ const B=[[a[0][0]-l,a[0][1],a[0][2]],[a[1][0],a[1][1]-l,a[1][2]],[a[2][0],a[2][1],a[2][2]-l]];
    const cr=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]; const cs=[cr(B[0],B[1]),cr(B[0],B[2]),cr(B[1],B[2])]; let v=cs[0]; for(const q of cs) if(Math.hypot(...q)>Math.hypot(...v)) v=q; const nv=Math.hypot(...v); if(nv>0) out.push(v.map(t=>t/nv)); }
  return out; }
function cubicRoots(a,b,c,d){ b/=a; c/=a; d/=a; const p=c-b*b/3, q=2*b*b*b/27-b*c/3+d; const D=q*q/4+p*p*p/27; const out=[];
  if(D>1e-14){ const u=Math.cbrt(-q/2+Math.sqrt(D)), v=Math.cbrt(-q/2-Math.sqrt(D)); out.push(u+v-b/3); }
  else if(Math.abs(D)<=1e-14){ const u=Math.cbrt(-q/2); out.push(2*u-b/3,-u-b/3); }
  else { const r=Math.sqrt(-p*p*p/27), ph=Math.acos(Math.max(-1,Math.min(1,-q/(2*r)))), m=2*Math.cbrt(r); for(let k=0;k<3;k++) out.push(m*Math.cos((ph+2*Math.PI*k)/3)-b/3); }
  return out; }
IC.ellipseFit=ellipseFit;

/* ---------- 기준선 자동 (세실) ----------
   윤곽(방울 + 기판 한 덩어리)에서 열마다 가장 위 윤곽 y → RANSAC 직선(기판 표면) → 안쪽 점 최소제곱
   반환 {p1:[x,y], p2:[x,y], angleDeg, inliers} — 방울 아래에 기판이 보여야 한다(거울 반사 기판은 수동 기준선) */
IC.dropBaseline=function(pts,edge,o){ o=o||{}; const top=new Map(); pts.forEach(([x,y],i)=>{ if(edge&&edge[i]) return; const k=Math.round(x); const v=top.get(k); if(v==null||y<v[1]) top.set(k,[x,y]); });
  const T=[...top.values()].sort((a,b)=>a[0]-b[0]); if(T.length<20) throw new Error("기준선: 점이 너무 적습니다");
  let best=null; const tol=o.tol||1.5; const N=T.length; let seed=12345; const rnd=()=>((seed=(seed*1103515245+12345)>>>0)/4294967296);
  for(let it=0;it<400;it++){ const a=T[Math.floor(rnd()*N)], b=T[Math.floor(rnd()*N)]; if(Math.abs(a[0]-b[0])<N*0.2) continue; const m=(b[1]-a[1])/(b[0]-a[0]); if(Math.abs(m)>0.3) continue;
    let n=0; for(const q of T) if(Math.abs(q[1]-(a[1]+m*(q[0]-a[0])))<=tol) n++; if(!best||n>best.n) best={n,a,m}; }
  if(!best||best.n<0.15*N) throw new Error("기준선(기판 표면)을 찾지 못했습니다 — 기준선을 직접 그으세요(Line ROI)");
  const In=T.filter(q=>Math.abs(q[1]-(best.a[1]+best.m*(q[0]-best.a[0])))<=tol); let sx=0,sy=0,sxx=0,sxy=0; for(const [x,y] of In){ sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; } const n=In.length; const m=(n*sxy-sx*sy)/(n*sxx-sx*sx), c=(sy-m*sx)/n;
  const x1=T[0][0], x2=T[N-1][0]; return {p1:[x1,c+m*x1],p2:[x2,c+m*x2],angleDeg:Math.atan(m)/D2R,inliers:n}; };

/* ---------- 접촉각 ----------
   o: {roi, sigma, baseline:{p1,p2} (없으면 자동), methods:["yl","circle","ellipse","poly"], polyFrac(0.15), pxSize(mm/px 선택), rhoIn/rhoOut/g}
   좌표: 기준선이 x' 축(왼→오른), 높이 y' 는 기판에서 위로(+). 반환 {contour, baseline, profile(방울 점, 기준선 좌표), methods:{…}, width, height} */
IC.contactAngle=function(plane,o){ o=o||{}; const ct=IC.dropContour(plane,{roi:o.roi,sigma:o.sigma,pick:"bottom"});
  let bl=o.baseline; let auto=false; if(!bl){ bl=IC.dropBaseline(ct.raw,ct.edge); auto=true; }
  const [ax,ay]=bl.p1, [bx,by]=bl.p2; const L=Math.hypot(bx-ax,by-ay); const ux=(bx-ax)/L, uy=(by-ay)/L; // 기준선 방향(영상)
  // 위쪽(+y') = 영상에서 기준선 왼쪽 법선 (영상 y 아래이므로 (uy,−ux))
  const nx=uy, ny=-ux; const toB=([x,y])=>[(x-ax)*ux+(y-ay)*uy,(x-ax)*nx+(y-ay)*ny]; const toI=([s,t])=>[ax+s*ux+t*nx,ay+s*uy+t*ny];
  // 기준선 가까이(접촉점 둘레 쐐기)는 두 가장자리가 겹쳐 경사 최대 다듬기가 흔들린다 → 문턱 등고선 그대로
  const res0={}; const rawB=ct.raw.map(toB); const Hr=amax(rawB.filter((_,i)=>!ct.edge[i]),q=>q[1]); const tKeep=Math.max(6,0.04*Hr);
  const all=ct.pts.map((q,i)=>rawB[i][1]<tKeep?rawB[i]:toB(q)); const Hmax=amax(all.filter((_,i)=>!ct.edge[i]),q=>q[1]); if(!(Hmax>4)) throw new Error("기준선 위에 방울이 없습니다");
  // 방울 안에 바늘(주사기)이 꽂혀 영역 위 테두리까지 이어지면: 위쪽 폭이 일정한 구간 = 바늘 → 그 아래 끝(+여유)보다 위 점은 맞춤에서 뺀다 (접선법은 접촉점 둘레만 보므로 영향 없음)
  let cutT=Infinity; if(ct.comp.y0<=ct.roi.y0+1){ try{ const nd=IC.dropNeedle(ct.raw,ct.edge); const tipT=Math.max(...[[nd.left.a+nd.left.b*nd.yTip,nd.yTip],[nd.right.a+nd.right.b*nd.yTip,nd.yTip]].map(q=>toB(q)[1]));
      if(nd.diam<0.8*(amax(rawB,q=>q[0])-amin(rawB,q=>q[0]))){ cutT=tipT-Math.max(4,0.03*Hr); res0.needle={diam:nd.diam,yTip:nd.yTip,cutT}; } }catch(e){} }
  const minH=Math.max(1.5,0.01*Hmax); res0.cut=rawB.some((q,i)=>ct.edge[i]&&q[1]>Math.max(3,0.03*Hmax)&&q[1]<cutT-2);
  const prof=all.filter((q,i)=>!ct.edge[i]&&q[1]>minH&&q[1]<cutT), profRaw=rawB.filter((q,i)=>!ct.edge[i]&&q[1]>minH&&q[1]<cutT); // 접선법은 접촉점 둘레만 보므로 문턱 등고선(profRaw)
  // 방울 점만: 기준선 위로 솟은 덩어리 — 꼭대기 열 둘레로 이어진 x' 범위 (기판 윗면 점은 높이가 minH 아래라 이미 빠짐)
  if(prof.length<20) throw new Error("방울 윤곽 점이 너무 적습니다");
  const res={nNaN:ct.nNaN,contour:ct.pts,edge:ct.edge,level:ct.level,baseline:{p1:bl.p1,p2:bl.p2,auto,angleDeg:Math.atan2(by-ay,bx-ax)/D2R},height:Hmax,methods:{},needle:res0.needle||null,cut:!!res0.cut,why:{}};
  // 접촉점 근사: 낮은 점들 중 왼쪽·오른쪽 끝
  let xl=Infinity,xr=-Infinity; for(const [s,t] of prof){ if(t<minH+0.1*Hmax){ if(s<xl) xl=s; if(s>xr) xr=s; } } res.width=xr-xl;
  const M=new Set(o.methods||["yl","circle","ellipse","poly"]); const inside=(sx,dir)=>dir; // 표시용
  // 접선 각: 접촉점에서 위로 가는 접선 (dx'/dy' = c) → 왼쪽 θ = atan2(1, c), 오른쪽 θ = atan2(1, −c)
  const angL=c=>Math.atan2(1,c)/D2R, angR=c=>Math.atan2(1,-c)/D2R;
  // --- 다항식(접선) ---
  if(M.has("poly")){ const side=sgn=>polyTangent(profRaw,sgn,xl,xr,Hmax,o); const l=side(-1), r=side(1);
    if(!(l&&r)) res.why.poly="접촉점 둘레 윤곽 점이 모자랍니다("+(!l&&!r?"양쪽":!l?"왼쪽":"오른쪽")+")";
    if(l&&r) res.methods.poly={left:l.theta,right:r.theta,mean:(l.theta+r.theta)/2,contactL:toI([l.x,0]),contactR:toI([r.x,0]),n:l.n+r.n,fitRadius:l.R}; }
  // --- 원 ---
  if(M.has("circle")){ const c=circleFit(prof); if(c&&c.r>0&&Math.abs(c.cy)<c.r){ const th=(90+Math.asin(c.cy/c.r)/D2R); const dx=Math.sqrt(c.r*c.r-c.cy*c.cy);
      res.methods.circle={left:th,right:th,mean:th,center:toI([c.cx,c.cy]),r:c.r,rms:c.rms,contactL:toI([c.cx-dx,0]),contactR:toI([c.cx+dx,0]),cB:[c.cx,c.cy]}; } }
  // --- 타원 ---
  if(M.has("ellipse")){ const e=ellipseFit(prof); if(e){ // y'=0 에서 A x² + D x + F = 0
      const {A,B,C,D,E,F}=e; const disc=D*D-4*A*F; if(disc>0&&A!==0){ const r1=(-D-Math.sqrt(disc))/(2*A), r2=(-D+Math.sqrt(disc))/(2*A); const xL=Math.min(r1,r2), xR=Math.max(r1,r2);
        // 접선: ∇ = (2Ax + By + D, Bx + 2Cy + E) 에 수직 → dx/dy = −(Bx+2Cy+E)/(2Ax+By+D)
        const sl=x=>-(B*x+E)/(2*A*x+D); res.methods.ellipse={left:angL(sl(xL)),right:angR(sl(xR)),contactL:toI([xL,0]),contactR:toI([xR,0]),conic:e};
        res.methods.ellipse.mean=(res.methods.ellipse.left+res.methods.ellipse.right)/2; } } }
  // --- Young–Laplace (축대칭) ---
  if(M.has("yl")){ try{ const Y=fitSessile(prof,Hmax,xl,xr); if(!Y) res.why.yl="초기값·수렴 실패 (방울 윗부분이 잘렸거나 윤곽이 방울이 아님)"; if(Y){ res.methods.yl=Y; Y.contactL=toI([Y.x0-Y.rc,0]); Y.contactR=toI([Y.x0+Y.rc,0]); Y.apexImg=toI([Y.x0,Y.yApex]); } }catch(e){ res.methods.ylError=String(e.message||e); } }
  // 물성(선택): 부피·Bond — 화소 크기 있을 때
  if(o.pxSize>0){ const s=o.pxSize; const Y=res.methods.yl; if(Y){ res.V=Y.V*Math.pow(Y.b*s,3); res.R0mm=Y.b*s; res.baseDiameterMM=2*Y.rc*s; res.heightMM=Y.yApex*s;
      const drho=(o.rhoIn==null?1000:+o.rhoIn)-(o.rhoOut==null?0:+o.rhoOut), g=o.g||9.80665; const relBo=Y.BoSd/Math.max(Y.Bo,1e-9); if(drho>0&&Y.Bo>0.1&&relBo<0.1&&Y.mean>=60&&!Y.atBound){ res.gammaEst=drho*g*Math.pow(Y.b*s*1e-3,2)/Y.Bo*1e3; res.gammaEstSd=res.gammaEst*Math.sqrt(4*Math.pow(Y.bSd/Y.b,2)+relBo*relBo); } } // 세실 형상에서 γ: Bo 가 충분히 크고 잘 정해질 때만 (작은 θ·작은 Bo 에서는 형상에 정보가 없음)
    const Cc=res.methods.circle; if(Cc){ const h=Cc.r+Cc.cB[1], a=Math.sqrt(Math.max(0,Cc.r*Cc.r-Cc.cB[1]*Cc.cB[1])); res.Vcap=Math.PI*h*(3*a*a+h*h)/6*Math.pow(s,3); } }
  // 화면 덧그림(영상 좌표): 기준선 · 방법마다 접촉점의 접선(액체 쪽에서 잰 θ) · YL·원 맞춤 곡선
  { const ov={baseline:[bl.p1.slice(),bl.p2.slice()],tangents:[],curves:[]}; const Lt=Math.max(20,0.35*Hmax);
    const tan=(m,k)=>{ if(!m||!m.contactL) return; const d=(th,sg)=>{ const a=th*D2R; const ds=sg*Math.cos(a), dt=Math.sin(a); return [ds*ux+dt*nx,ds*uy+dt*ny]; };
      const dl=d(m.left,1), dr=d(m.right,-1); ov.tangents.push({m:k,p:m.contactL,d:dl,len:Lt,theta:m.left,side:"L"},{m:k,p:m.contactR,d:dr,len:Lt,theta:m.right,side:"R"}); };
    for(const k of ["yl","poly","ellipse","circle"]) tan(res.methods[k],k);
    const Y=res.methods.yl; if(Y){ const half=sg=>{ const L=[]; for(let i=0;i<Y.curve.x.length;i++){ const z=Y.curve.z[i]*Y.b; if(z>Y.yApex+1e-9) break; L.push(toI([Y.x0+sg*Y.curve.x[i]*Y.b,Y.yApex-z])); } return L; }; ov.curves.push({m:"yl",pts:half(-1).reverse().concat(half(1).slice(1))}); }
    const Cc=res.methods.circle; if(Cc){ const L=[]; const a1=Math.asin(Math.max(-1,Math.min(1,-Cc.cB[1]/Cc.r))); for(let i=0;i<=120;i++){ const a=a1+(Math.PI-2*a1)*i/120; L.push(toI([Cc.cB[0]+Cc.r*Math.cos(a),Cc.cB[1]+Cc.r*Math.sin(a)])); } ov.curves.push({m:"circle",pts:L}); }
    if(res.needle){ const c=res.needle.cutT; ov.needle=[[toI([xl-0.1*res.width,c]),toI([xr+0.1*res.width,c])]]; } // 바늘 때문에 이 선 위 윤곽은 뺌
    res.overlay=ov; }
  res.warn=[]; if(res.cut) res.warn.push("방울 윤곽이 영역(ROI) 테두리에 잘렸습니다 — 방울 전체가 들어가게 영역을 넓히세요. 잘린 채 맞춘 값은 틀릴 수 있습니다");
  for(const [k,v] of Object.entries(res.why)) res.warn.push(({poly:"접선",yl:"Young–Laplace"})[k]+" 없음: "+v);
  if(res.methods.yl&&res.methods.yl.atBound) res.warn.push("Young–Laplace Bond 수가 맞춤 한계에 닿았습니다 — 그 θ 는 믿지 마세요");
  if(res.nNaN) res.warn.push(`영역 안 NaN 화소 ${res.nNaN} 개를 바탕 값으로 채웠습니다`);
  if(res.needle) res.warn.push(`방울 안에 바늘(폭 ${res.needle.diam.toFixed(0)} px)이 있어 그 아래 끝보다 위의 윤곽은 맞춤에서 뺐습니다 — 높이·부피는 바늘 부피를 포함하지 않은 값`); const vals=Object.entries(res.methods).filter(([k,v])=>v&&v.mean!=null).map(([k,v])=>v.mean); if(vals.length>1&&Math.max(...vals)-Math.min(...vals)>5) res.warn.push(`방법 사이 차이 ${(Math.max(...vals)-Math.min(...vals)).toFixed(1)}° — 중력으로 납작한(Bo 큼) 대칭 방울이면 원·타원이 틀리고 Young–Laplace 가 맞습니다. 왼·오른 각이 다르면(비대칭·접촉선 고정) 축대칭 가정이 깨져 Young–Laplace 도 틀릴 수 있으니 접선법 왼·오른 값을 함께 보세요`);
  const Pm=res.methods.poly; if(Pm&&Math.abs(Pm.left-Pm.right)>5) res.warn.push(`접선법 왼·오른 차 ${Math.abs(Pm.left-Pm.right).toFixed(1)}° — 방울이 비대칭(고정·이력)이거나 기준선이 기울었습니다. 축대칭인 Young–Laplace·원·타원은 이 경우 평균만 줍니다`);
  if(Pm&&Pm.mean>150) res.warn.push("θ > 150° 에서는 접촉점 둘레 쐐기가 좁아 접선법 오차가 커집니다 — Young–Laplace 권장");
  if(Y_ok(res)&&res.methods.yl.rms>1) res.warn.push(`Young–Laplace 잔차 rms ${res.methods.yl.rms.toFixed(2)} px — 바늘·반사·기판 그림자가 윤곽에 섞였는지 보세요`);
  if(res.methods.ylError) res.warn.push("Young–Laplace 맞춤 실패: "+res.methods.ylError);
  return res; };
const Y_ok=r=>!!(r.methods&&r.methods.yl);

/* 접선(다항식) 접촉각 — 접촉점 둘레 반지름 R 안의 윤곽 점을 접촉점에서 멀어지는 순으로 늘어놓고
   호 길이 u 의 2 차식 x'(u), t(u) 로 맞춘 뒤 t(u)=0 (기준선) 으로 외삽 → 그 자리 접선 (x', t) 미분 방향으로 θ.
   x(t) 꼴 맞춤과 달리 θ > 90° (방울이 접촉점 위로 넘어감) 에서도 한 값. 바닥 가까이 흐림 구간(t < tMin) 은 뺌. */
function polyTangent(prof,sgn,xl,xr,Hmax,o){ const mid=(xl+xr)/2; let xc=sgn<0?xl:xr, tb=Infinity; const tMin0=Math.max(o.polyTmin||3,0.01*Hmax);
  for(const [s,t] of prof) if((sgn<0?s<mid:s>mid)&&t>=tMin0&&t<tb){ tb=t; xc=s; } // 접촉점 근사 = 그 쪽에서 기준선에 가장 가까운 점 (θ > 90° 면 불룩한 곳보다 안쪽)
  const R=Math.max(8,(o.polyFrac||0.12)*Math.min(Hmax,(xr-xl)/2)); const tMin=Math.max(o.polyTmin||3,0.01*Hmax);
  const Q=prof.filter(([s,t])=>(sgn<0?s<mid:s>mid)&&t>=tMin&&Math.hypot(s-xc,t)<=R+tMin); if(Q.length<6) return null;
  Q.sort((a,b)=>Math.hypot(a[0]-xc,a[1])-Math.hypot(b[0]-xc,b[1]));
  const U=[0]; for(let i=1;i<Q.length;i++) U.push(U[i-1]+Math.hypot(Q[i][0]-Q[i-1][0],Q[i][1]-Q[i-1][1]));
  const fit2=(V)=>{ const S=[[0,0,0],[0,0,0],[0,0,0]], b=[0,0,0]; for(let i=0;i<Q.length;i++){ const u=U[i], v=[1,u,u*u]; for(let a=0;a<3;a++){ b[a]+=v[a]*V(i); for(let c=0;c<3;c++) S[a][c]+=v[a]*v[c]; } } return solveSym(S,b); };
  const cx=fit2(i=>Q[i][0]), ct=fit2(i=>Q[i][1]); if(!cx||!ct) return null;
  // t(u) = 0 의 u (가장 0 에 가까운 근, 보통 음수 = 외삽)
  let u0; if(Math.abs(ct[2])<1e-12) u0=-ct[0]/ct[1]; else { const D=ct[1]*ct[1]-4*ct[2]*ct[0]; if(D<0) u0=-ct[1]/(2*ct[2]); else { const r1=(-ct[1]-Math.sqrt(D))/(2*ct[2]), r2=(-ct[1]+Math.sqrt(D))/(2*ct[2]); u0=Math.abs(r1)<Math.abs(r2)?r1:r2; } }
  const dx=cx[1]+2*cx[2]*u0, dt=ct[1]+2*ct[2]*u0; const x0=cx[0]+cx[1]*u0+cx[2]*u0*u0;
  // u 가 커지면 위로(t 증가) — 접선 방향 (dx, dt), 액체 쪽에서 잰 각: 왼쪽 θ = atan2(dt, dx), 오른쪽 θ = atan2(dt, −dx)
  const th=(sgn<0?Math.atan2(dt,dx):Math.atan2(dt,-dx))/D2R; return {theta:th,x:x0,n:Q.length,R}; }
/* 세실 Young–Laplace 맞춤: 매개변수 [x0, yApex(기준선 위 높이), b, Bo] (기준선 좌표, 축은 기준선에 수직) */
function fitSessile(prof,Hmax,xl,xr){ const Q0=prof; let Q=Q0; if(Q.length>900){ const st=Q.length/900; Q=[...Array(900)].map((_,i)=>Q0[Math.floor(i*st)]); }
  const c0=circleFit(Q); const x00=c0?c0.cx:(xl+xr)/2, b0=c0?c0.r:Hmax/2;
  const cache=new Map(); const curveOf=(Bo,zm)=>{ const k=Bo.toPrecision(10)+"|"+zm.toFixed(2); let C=cache.get(k); if(!C){ C=monoPrefix(IC.dropYL(Bo,"sessile",{zMax:zm,ds:0.003})); if(cache.size>64) cache.clear(); cache.set(k,C); } return C; };
  const res=p=>{ const [x0,ya,b,Bo]=p; if(!(b>2)||!(Bo>-0.05&&Bo<5)||!(ya>1)) return null; const zm=Math.ceil((ya/b*1.02+0.05)*4)/4; const C=curveOf(Bo,zm); if(C.n<10) return null;
    if(C.z[C.n-1]<ya/b){ // φ 가 π 에 닿아 기준선까지 못 감 → 벌점
      return null; }
    return Q.map(([s,t])=>b*curveDist(C,Math.abs(s-x0)/b,(ya-t)/b,40)); };
  let best=null; for(const Bo of [0.0005,0.05,0.15,0.3,0.6,1.0,1.6]){ for(const bf of [0.8,1,1.25]){ const p=[x00,Hmax,b0*bf,Bo]; const r=res(p); if(!r) continue; const c=r.reduce((a,v)=>a+v*v,0); if(!best||c<best.c) best={c,p}; } }
  if(!best) return null;
  const F=lm(res,best.p,{iter:80,h:[0.02,0.02,0.01,1e-4],clamp:q=>{ if(q[3]<-0.04) q[3]=-0.04; if(q[2]<3) q[2]=3; }}); if(!F) return null; const atBound=F.p[3]<=-0.0399;
  const [x0,ya,b,Bo]=F.p; const C=monoPrefix(IC.dropYL(Bo,"sessile",{zMax:ya/b+0.2,ds:0.002})); const zt=ya/b; const th=phiAtZ(C,zt)/D2R; const rc=xAtZ(C,zt)*b; const VA=IC.dropVolArea({x:Array.from(C.x.subarray(0,C.n)),z:Array.from(C.z.subarray(0,C.n))},zt);
  // θ 불확도: Bo·b·yApex 공분산을 수치 미분으로 전파
  let thSd=NaN; if(F.cov){ const thOf=p=>{ const Cq=monoPrefix(IC.dropYL(p[3],"sessile",{zMax:p[1]/p[2]+0.2,ds:0.002})); return phiAtZ(Cq,p[1]/p[2])/D2R; }; const g=[0,0,0,0]; const hs=[0,0.05,0.05,1e-3];
    for(const j of [1,2,3]){ const q=F.p.slice(); q[j]+=hs[j]; const t2=thOf(q); g[j]=(t2-th)/hs[j]; } let v=0; for(let i=0;i<4;i++) for(let j=0;j<4;j++) v+=g[i]*g[j]*F.cov[i][j]; thSd=Math.sqrt(Math.max(0,v)); }
  return {left:th,right:th,mean:th,thetaSd:thSd,atBound,BoSd:F.sd[3],bSd:F.sd[2],x0,yApex:ya,b,Bo,rc,rms:F.rms,V:VA.V,A:VA.A,curve:{x:Array.from(C.x.subarray(0,C.n)),z:Array.from(C.z.subarray(0,C.n))},n:Q.length}; }
IC._fitSessile=fitSessile;
})();
