/* ============================================================
   INNO Analytics — Diffraction: 영상에서 패턴 읽기 (SAED 스폿·링). DOM 없음, Worker·Node 공용.
   · mask: 빔스토퍼(어둡고 큰 덩어리), 영상 위 글자·스케일바(가늘고 밝은 선)
   · spots: 대역 통과(가우스 차) → 국소 최대 → 배경 뺀 무게 중심 (부화소)
   · center: 직접빔 → 없으면 Friedel 짝 투표(세기 닮음 가중, 빔스토퍼 안 우선)
   · lattice: 짧은 벡터 쌍으로 2D 역격자 → 가장 성긴(점유율 높은) 격자 → 모든 지수 스폿으로 최소제곱(중심 포함) → R1·R2·R3, 불확도
   · ring: 광선별 링 위치 → 타원 맞춤(Fitzgibbon) → 타원 보정 방위 적분 → 배경 → 봉우리(부화소)
   · scale bar: 그림 위 가로 막대 길이(px) — 값은 사용자가 확인
   좌표: 화소 (x 오른쪽, y 아래). 역공간 길이 = 화소 × scale (1/nm/px).
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const DP={};

/* ---------- 기본 ---------- */
function gauss1d(s){ const r=Math.max(1,Math.ceil(3*s)); const k=new Float32Array(2*r+1); let sum=0; for(let i=-r;i<=r;i++){ const v=Math.exp(-i*i/(2*s*s)); k[i+r]=v; sum+=v; } for(let i=0;i<k.length;i++) k[i]/=sum; return k; }
/* 가우스 흐림 (가장자리 = 복제), mask 가 있으면 가려진 화소는 빼고 정규화 */
/* 큰 σ: 상자 흐림 세 번 (Wells 1986 / Kovesi 2010 너비) — 2048² 에서 σ 18 도 O(n). 작은 σ 는 정확한 가우스 */
function boxWidths(s,n){ const wI=Math.floor(Math.sqrt(12*s*s/n+1)); const wl=wI%2?wI:wI-1; const wu=wl+2; const m=Math.round((12*s*s-n*wl*wl-4*n*wl-3*n)/(-4*wl-4)); return Array.from({length:n},(_,i)=>i<m?wl:wu); }
function box1(inp,out,w,h,r,horiz){ const L=horiz?w:h, M=horiz?h:w, st=horiz?1:w, lst=horiz?w:1; const inv=1/(2*r+1);
  for(let j=0;j<M;j++){ const b=j*lst; const first=inp[b], last=inp[b+(L-1)*st]; let acc=(r+1)*first; for(let i=0;i<r;i++) acc+=inp[b+Math.min(i,L-1)*st];
    for(let i=0;i<L;i++){ const ip=i+r, im=i-r-1; acc+=(ip<L?inp[b+ip*st]:last)-(im>=0?inp[b+im*st]:first); out[b+i*st]=acc*inv; } } }
function boxBlur(src,w,h,s){ const ws=boxWidths(s,3); let a=Float32Array.from(src), t=new Float32Array(src.length); for(const bw of ws){ const r=(bw-1)>>1; box1(a,t,w,h,r,true); box1(t,a,w,h,r,false); } return a; }
DP.blur=function(a,w,h,s,mask){ if(!(s>0)) return Float32Array.from(a);
  if(s>3){ const src=new Float32Array(a.length); let wt=null; if(mask){ wt=new Float32Array(a.length); for(let i=0;i<a.length;i++){ const m=mask[i]?0:1; src[i]=a[i]*m; wt[i]=m; } } else src.set(a);
    const o=boxBlur(src,w,h,s); if(wt){ const o2=boxBlur(wt,w,h,s); for(let i=0;i<a.length;i++) o[i]=o2[i]>1e-3?o[i]/o2[i]:0; } return o; }
  const k=gauss1d(s), r=(k.length-1)>>1;
  const src=new Float32Array(a.length), wt=mask?new Float32Array(a.length):null; for(let i=0;i<a.length;i++){ const m=mask&&mask[i]?0:1; src[i]=a[i]*m; if(wt) wt[i]=m; }
  const pass=(inp,out,horiz)=>{ for(let y=0;y<h;y++) for(let x=0;x<w;x++){ let s2=0; for(let j=-r;j<=r;j++){ let xx=horiz?x+j:x, yy=horiz?y:y+j; if(xx<0) xx=0; else if(xx>=w) xx=w-1; if(yy<0) yy=0; else if(yy>=h) yy=h-1; s2+=k[j+r]*inp[yy*w+xx]; } out[y*w+x]=s2; } };
  const t=new Float32Array(a.length), o=new Float32Array(a.length); pass(src,t,true); pass(t,o,false);
  if(wt){ const t2=new Float32Array(a.length), o2=new Float32Array(a.length); pass(wt,t2,true); pass(t2,o2,false); for(let i=0;i<a.length;i++) o[i]=o2[i]>1e-3?o[i]/o2[i]:0; }
  return o; };
DP.percentile=function(a,p,mask){ const v=[]; const st=Math.max(1,Math.floor(a.length/400000)); for(let i=0;i<a.length;i+=st) if(!(mask&&mask[i])&&a[i]===a[i]) v.push(a[i]); if(!v.length) return NaN; v.sort((x,y)=>x-y); return v[Math.min(v.length-1,Math.max(0,Math.floor(p/100*(v.length-1))))]; };
/* 연결 성분 (4-연결). sel(i) 참인 화소 → [{px:[...], n, x0,x1,y0,y1}] */
DP.components=function(w,h,sel,minN){ const lab=new Int32Array(w*h).fill(-1); const out=[]; const st=new Int32Array(w*h);
  for(let i=0;i<w*h;i++){ if(lab[i]>=0||!sel(i)) continue; let sp=0; st[sp++]=i; lab[i]=out.length; const px=[]; let x0=w,x1=0,y0=h,y1=0;
    while(sp){ const p=st[--sp]; px.push(p); const x=p%w, y=(p-x)/w; if(x<x0) x0=x; if(x>x1) x1=x; if(y<y0) y0=y; if(y>y1) y1=y;
      if(x>0&&lab[p-1]<0&&sel(p-1)){ lab[p-1]=out.length; st[sp++]=p-1; } if(x<w-1&&lab[p+1]<0&&sel(p+1)){ lab[p+1]=out.length; st[sp++]=p+1; }
      if(y>0&&lab[p-w]<0&&sel(p-w)){ lab[p-w]=out.length; st[sp++]=p-w; } if(y<h-1&&lab[p+w]<0&&sel(p+w)){ lab[p+w]=out.length; st[sp++]=p+w; } }
    out.push({px,n:px.length,x0,x1,y0,y1}); }
  return minN?out.filter(c=>c.n>=minN):out; };
/* 거리 변환(체스판 근사 두 번 훑기) — 성분 안쪽 깊이 */
DP.chamfer=function(w,h,inside){ const INF=1e9, d=new Float32Array(w*h); for(let i=0;i<w*h;i++) d[i]=inside[i]?INF:0; const a=1, b=1.4142;
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const i=y*w+x; if(!d[i]) continue; let v=d[i]; if(x>0) v=Math.min(v,d[i-1]+a); if(y>0){ v=Math.min(v,d[i-w]+a); if(x>0) v=Math.min(v,d[i-w-1]+b); if(x<w-1) v=Math.min(v,d[i-w+1]+b); } d[i]=v; }
  for(let y=h-1;y>=0;y--) for(let x=w-1;x>=0;x--){ const i=y*w+x; if(!d[i]) continue; let v=d[i]; if(x<w-1) v=Math.min(v,d[i+1]+a); if(y<h-1){ v=Math.min(v,d[i+w]+a); if(x<w-1) v=Math.min(v,d[i+w+1]+b); if(x>0) v=Math.min(v,d[i+w-1]+b); } d[i]=v; }
  for(let i=0;i<w*h;i++) if(d[i]>=INF) d[i]=Math.max(w,h); return d; };
const dilate=(w,h,m,r)=>{ if(r<=0) return m; const o=new Uint8Array(m.length); const R=Math.ceil(r), r2=r*r;
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ if(!m[y*w+x]) continue; for(let dy=-R;dy<=R;dy++){ const yy=y+dy; if(yy<0||yy>=h) continue; for(let dx=-R;dx<=R;dx++){ if(dx*dx+dy*dy>r2) continue; const xx=x+dx; if(xx>=0&&xx<w) o[yy*w+xx]=1; } } } return o; };
DP.dilate=dilate;

/* ---------- 덧그림(글자·스케일바) 과 빔스토퍼 ---------- */
/* 8-bit 그림의 덧그림: 포화(≥ 최대의 97 %) 화소의 성분 중 가늘거나(깊이 ≤ 2.5) 길쭉한 것. 스케일바 = 가로로 긴 막대 */
DP.findOverlay=function(img,w,h){ let mx=-Infinity; for(let i=0;i<img.length;i++) if(img[i]>mx) mx=img[i]; const thr=mx*0.97;
  const sel=i=>img[i]>=thr; const comps=DP.components(w,h,sel,6); const inside=new Uint8Array(w*h); for(let i=0;i<w*h;i++) inside[i]=img[i]>=thr?1:0; const dt=DP.chamfer(w,h,inside);
  const mask=new Uint8Array(w*h); const bars=[], text=[];
  for(const c of comps){ const bw=c.x1-c.x0+1, bh=c.y1-c.y0+1; let dmax=0; for(const p of c.px) if(dt[p]>dmax) dmax=dt[p]; const fill=c.n/(bw*bh);
    const isBar=bw>=Math.max(30,0.04*w)&&bh<=Math.max(10,bw/8)&&fill>0.5;               // 꽉 찬 가로 막대
    const isHollowBar=bw>=Math.max(30,0.04*w)&&bh<=Math.max(40,bw/4)&&fill<0.5&&dmax<=3;  // 속 빈 사각 테두리 (REDI 예제)
    const thin=dmax<=2.5&&(fill<0.6||c.n>=40);
    if(isBar||isHollowBar){ bars.push({x0:c.x0,x1:c.x1,y0:c.y0,y1:c.y1,len:bw,hollow:isHollowBar,thick:isHollowBar?Math.round(dmax*2):bh}); for(const p of c.px) mask[p]=1; }
    else if(thin){ text.push({x0:c.x0,x1:c.x1,y0:c.y0,y1:c.y1}); for(const p of c.px) mask[p]=1; } }
  return {mask:dilate(w,h,mask,3),bars,text}; };
/* 빔스토퍼: 배경 중앙값의 25 % 보다 어둡고 큰(영상의 0.2 % 이상) 성분 */
DP.findBeamStop=function(img,w,h,opt){ opt=opt||{}; const med=DP.percentile(img,50); const lo=DP.percentile(img,1); const thr=lo+(med-lo)*(opt.frac||0.25);
  const comps=DP.components(w,h,i=>img[i]<=thr,Math.max(50,Math.round(0.002*w*h))); const mask=new Uint8Array(w*h); const parts=[];
  for(const c of comps){ // 영상 가장자리의 어두운 띠(검출기 바깥)도 같은 방법으로 가린다
    let sx=0,sy=0; for(const p of c.px){ mask[p]=1; sx+=p%w; sy+=(p-p%w)/w; } parts.push({n:c.n,cx:sx/c.n,cy:sy/c.n,x0:c.x0,x1:c.x1,y0:c.y0,y1:c.y1}); }
  return {mask:dilate(w,h,mask,opt.pad!=null?opt.pad:4),parts}; };

/* ---------- 스폿 ---------- */
/* 칸별 MAD 잡음 지도 (가려진 화소 뺌). 쓸 화소가 모자란 칸은 전체 값으로 */
DP.noiseMap=function(d,w,h,mask,blk){ const nx=Math.max(1,Math.round(w/blk)), ny=Math.max(1,Math.round(h/blk)); const bw=w/nx, bh=h/ny; const g=new Float64Array(nx*ny).fill(NaN);
  const madOf=v=>{ if(v.length<16) return NaN; v.sort((a,b)=>a-b); const m=v[v.length>>1]; const r=v.map(x=>Math.abs(x-m)).sort((a,b)=>a-b); return r[r.length>>1]*1.4826; };
  const all=[]; for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){ const v=[]; const x0=Math.floor(i*bw), x1=Math.floor((i+1)*bw), y0=Math.floor(j*bh), y1=Math.floor((j+1)*bh); const st=Math.max(1,Math.floor(Math.sqrt((x1-x0)*(y1-y0)/4000)));
    for(let y=y0;y<y1;y+=st) for(let x=x0;x<x1;x+=st){ const k=y*w+x; if(mask&&mask[k]) continue; v.push(d[k]); } const m=madOf(v); g[j*nx+i]=m; if(m===m&&m>0) all.push(m); }
  all.sort((a,b)=>a-b); const glob=all.length?all[all.length>>1]:1e-9; const floor=Math.max(1e-12,(all[0]||glob)*0.5);
  for(let k=0;k<g.length;k++){ if(!(g[k]===g[k])) g[k]=glob; else if(g[k]<floor) g[k]=floor; }
  const at=(x,y)=>{ const fx=Math.min(nx-1,Math.max(0,x/bw-0.5)), fy=Math.min(ny-1,Math.max(0,y/bh-0.5)); const i0=Math.floor(fx), j0=Math.floor(fy), i1=Math.min(nx-1,i0+1), j1=Math.min(ny-1,j0+1); const tx=fx-i0, ty=fy-j0;
    return (1-tx)*(1-ty)*g[j0*nx+i0]+tx*(1-ty)*g[j0*nx+i1]+(1-tx)*ty*g[j1*nx+i0]+tx*ty*g[j1*nx+i1]; };
  return {grid:g,nx,ny,at,global:glob}; };
/* opt: sigma(스폿 크기, px), bg(배경 σ), k(잡음 배수), minDist, mask */
DP.findSpots=function(img,w,h,opt){ opt=opt||{}; const sz=Math.max(w,h); const s1=opt.sigma||Math.max(1,sz/1024*1.5), s2=opt.bg||s1*6; const mask=opt.mask||null;
  const a=DP.blur(img,w,h,s1,mask), b=DP.blur(img,w,h,s2,mask); const d=new Float32Array(w*h); for(let i=0;i<w*h;i++) d[i]=mask&&mask[i]?0:a[i]-b[i];
  // 잡음: 칸마다 대역 통과 값의 MAD → 쌍선형 보간. Poisson 잡음은 배경 세기를 따르므로, 영상 전체 한 값이면 밝은 가운데에서 거짓 스폿이 생긴다
  const nz=DP.noiseMap(d,w,h,mask,opt.block||Math.max(32,Math.round(8*s2))); const mad=nz.global; const kk=opt.k||6;
  const R=Math.max(2,Math.round(opt.minDist||s1*3)); const cand=[];
  for(let y=R;y<h-R;y++) for(let x=R;x<w-R;x++){ const i=y*w+x; const val=d[i]; if(!(val>0)) continue; if(mask&&mask[i]) continue; if(val<kk*nz.at(x,y)) continue; let isMax=true;
    for(let dy=-R;dy<=R&&isMax;dy++) for(let dx=-R;dx<=R;dx++){ if(!dx&&!dy) continue; const j=i+dy*w+dx; if(d[j]>val||(d[j]===val&&j<i)){ isMax=false; break; } }
    if(isMax) cand.push(i); }
  // 부화소: 창(반지름 2σ+1) 안 (원래값 − 큰 흐림 배경) 의 양수 무게 중심, 세기 = 그 합
  const rw=Math.max(2,Math.round(2*s1+1)); const spots=[];
  for(const i of cand){ const x=i%w, y=(i-x)/w; let sw=0,sx=0,sy=0,peak=-Infinity, sat=0; let masked=false;
    for(let dy=-rw;dy<=rw;dy++) for(let dx=-rw;dx<=rw;dx++){ const xx=x+dx, yy=y+dy; if(xx<0||yy<0||xx>=w||yy>=h) continue; if(dx*dx+dy*dy>rw*rw) continue; const j=yy*w+xx; if(mask&&mask[j]){ masked=true; continue; } const q=img[j]-b[j]; if(img[j]>peak) peak=img[j]; if(q>0){ sw+=q; sx+=q*xx; sy+=q*yy; } }
    if(!(sw>0)) continue; spots.push({x:sx/sw,y:sy/sw,I:sw,peak,snr:d[i]/nz.at(x,y),edgeMasked:masked}); void sat; }
  spots.sort((p,q)=>q.I-p.I); return {spots,noise:mad,k:kk,noiseMap:nz}; };

/* ---------- 중심 ---------- */
/* 직접빔이 보이면 가장 크고 밝은 스폿(영상 가운데 쪽), 아니면 Friedel 짝 중점 투표 */
DP.findCenter=function(spots,w,h,opt){ opt=opt||{}; const stop=opt.beamStop; const top=spots.slice(0,Math.min(80,spots.length));
  if(top.length<2) return null; const tol=opt.tol||Math.max(1.5,Math.max(w,h)/1024*2);
  // 후보 중점: 세기가 비슷한 짝 (작은 쪽/큰 쪽 ≥ 0.33)
  const votes=[]; for(let i=0;i<top.length;i++) for(let j=i+1;j<top.length;j++){ const a=top[i], b=top[j]; const r=Math.min(a.I,b.I)/Math.max(a.I,b.I); if(r<0.33) continue; votes.push({x:(a.x+b.x)/2,y:(a.y+b.y)/2,w:r*Math.sqrt(Math.min(a.I,b.I)),i,j}); }
  if(!votes.length) return null;
  // 격자 칸에 모아 가장 무거운 무리
  const cell=tol; const map=new Map(); for(const v of votes){ const k=Math.round(v.x/cell)+","+Math.round(v.y/cell); const e=map.get(k)||{w:0,sx:0,sy:0,n:0,pairs:new Set()}; e.w+=v.w; e.sx+=v.w*v.x; e.sy+=v.w*v.y; e.n++; e.pairs.add(v.i); e.pairs.add(v.j); map.set(k,e); }
  // 이웃 칸 합치기 (칸 경계에 걸린 무리)
  const keys=[...map.keys()]; let best=null;
  for(const k of keys){ const [cx,cy]=k.split(",").map(Number); let W=0,SX=0,SY=0,N=0; const P=new Set(); for(let dx=-1;dx<=1;dx++) for(let dy=-1;dy<=1;dy++){ const e=map.get((cx+dx)+","+(cy+dy)); if(!e) continue; W+=e.w; SX+=e.sx; SY+=e.sy; N+=e.n; e.pairs.forEach(p=>P.add(p)); }
    let score=W*Math.sqrt(P.size); const x=SX/W, y=SY/W;
    if(stop){ const inStop=stop[Math.round(y)*w+Math.round(x)]; if(inStop) score*=2; }   // 직접빔은 대개 빔스토퍼 뒤
    if(!best||score>best.score) best={x,y,score,n:N,members:P.size}; }
  return best; };

/* ---------- 2D 역격자 ---------- */
const det2=(a,b)=>a[0]*b[1]-a[1]*b[0];
function indexAll(vs,a,b,tol){ const D=det2(a,b); const out=[]; for(let i=0;i<vs.length;i++){ const v=vs[i]; const n1=(v[0]*b[1]-v[1]*b[0])/D, n2=(a[0]*v[1]-a[1]*v[0])/D; const m1=Math.round(n1), m2=Math.round(n2); const rx=v[0]-m1*a[0]-m2*b[0], ry=v[1]-m1*a[1]-m2*b[1]; const res=Math.hypot(rx,ry); out.push({n:[m1,m2],res,ok:res<=tol&&(m1||m2)}); } return out; }
/* 가우스 축약: |a| ≤ |b|, a·b ≤ 0, |a·b| ≤ |a|²/2 → R1=a, R2=b, R3=a+b (짧은 대각선, 둘 사이) */
DP.reduce2=function(a,b){ a=a.slice(); b=b.slice(); for(let it=0;it<50;it++){ if(Math.hypot(...a)>Math.hypot(...b)){ const t=a; a=b; b=t; } const mu=Math.round((a[0]*b[0]+a[1]*b[1])/(a[0]*a[0]+a[1]*a[1])); if(!mu) break; b=[b[0]-mu*a[0],b[1]-mu*a[1]]; }
  if(a[0]*b[0]+a[1]*b[1]>0) b=[-b[0],-b[1]]; if(det2(a,b)<0){ /* 손잡이: a→b 반시계(화면 y 아래이므로 시계) 로 통일 */ } return [a,b]; };
/* spots(px), center(px) → 격자. opt.tol(px), opt.maxCand */
DP.fitLattice=function(spots,center,opt){ opt=opt||{}; const N=spots.length; if(N<3) return null; const tol=opt.tol||Math.max(1.5,(opt.size||1024)/1024*2.5);
  const vs=spots.map(s=>[s.x-center.x,s.y-center.y]); const ws=spots.map(s=>Math.sqrt(Math.max(s.I,0)));
  // 격자 고르기에는 뚜렷한 스폿(SNR ≥ minSnrBasis)만 — 문턱 근처의 거짓 스폿이 촘촘한 가짜 격자를 부르지 않게. 지수 붙이기·최소제곱은 모든 스폿
  const snrMin=opt.minSnrBasis!=null?opt.minSnrBasis:9; let inB=spots.map(s=>!(s.snr<snrMin)); if(inB.filter(Boolean).length<6) inB=spots.map(()=>true);
  const totalW=ws.reduce((a,b,i)=>a+(inB[i]?b:0),0);
  // 짧은 벡터 후보: 길이순 + 세기
  const order=vs.map((v,i)=>({i,r:Math.hypot(v[0],v[1])})).filter(o=>o.r>3*tol&&inB[o.i]).sort((p,q)=>p.r-q.r).slice(0,opt.maxCand||16).map(o=>o.i);
  // 차 벡터도 후보 (중심이 불확실할 때 쓸모)
  const cand=order.map(i=>vs[i]); for(let a=0;a<Math.min(8,order.length);a++) for(let b=a+1;b<Math.min(8,order.length);b++){ const d=[vs[order[a]][0]-vs[order[b]][0],vs[order[a]][1]-vs[order[b]][1]]; if(Math.hypot(...d)>3*tol) cand.push(d); }
  const rmax=Math.max(...vs.map(v=>Math.hypot(v[0],v[1])));
  let best=null; const tried=[];
  for(let p=0;p<cand.length;p++) for(let q=p+1;q<cand.length;q++){ const a=cand[p], b=cand[q]; const ar=Math.abs(det2(a,b)); if(ar<(3*tol)**2) continue; if(Math.abs(det2(a,b))/(Math.hypot(...a)*Math.hypot(...b))<0.17) continue; // 10° 미만은 거의 같은 방향
    const ix=indexAll(vs,a,b,tol); let wsum=0,nIdx=0; for(let i=0;i<N;i++) if(ix[i].ok&&inB[i]){ wsum+=ws[i]; nIdx++; }
    const frac=wsum/totalW; const nPred=Math.PI*rmax*rmax/ar; const occ=Math.min(1,(nIdx+1)/Math.max(1,nPred));
    tried.push({a,b,frac,occ,nIdx,area:ar}); }
  if(!tried.length) return null; const fmax=Math.max(...tried.map(t=>t.frac));
  // 설명하는 세기가 최대의 92 % 이상인 것 중 가장 성긴(점유율 큰) 격자 — 더 촘촘한 격자는 모든 점을 설명하지만 빈 자리가 많다
  const good=tried.filter(t=>t.frac>=0.92*fmax).sort((p,q)=>q.occ-p.occ||q.frac-p.frac); best=good[0];
  let [a,b]=DP.reduce2(best.a,best.b);
  const R=DP.refineLattice(vs,ws,a,b,tol,opt); if(!R) return null; a=R.a; b=R.b; const c=R.c;
  const r3=[a[0]+b[0],a[1]+b[1]]; const cen={x:center.x+c[0],y:center.y+c[1]};
  const idx=spots.map((s,i)=>Object.assign({},s,{n:R.ix[i].n,res:R.ix[i].res,indexed:R.ix[i].ok}));
  const la=Math.hypot(...a), lb=Math.hypot(...b), l3=Math.hypot(...r3);
  const ang=(u,v)=>Math.acos(Math.max(-1,Math.min(1,(u[0]*v[0]+u[1]*v[1])/(Math.hypot(...u)*Math.hypot(...v)))))*180/Math.PI;
  // 불확도: 공분산 s²(XᵀWX)⁻¹ (x·y 같은 설계 행렬, 독립) → 길이는 방향 성분, 각도는 수직 성분
  const C=R.inv, s2=R.s2; const sa=Math.sqrt(C[0][0]*s2), sb=Math.sqrt(C[1][1]*s2), s3=Math.sqrt(Math.max(0,(C[0][0]+C[1][1]+2*C[0][1])*s2));
  const sigAng=180/Math.PI*Math.sqrt(C[0][0]*s2/(la*la)+C[1][1]*s2/(lb*lb));
  return {center:cen,R1:a,R2:b,R3:r3,len:[la,lb,l3],sig:[sa,sb,s3],ang12:ang(a,b),ang13:ang(a,r3),ang23:ang(b,r3),sigAng,sigCenter:Math.sqrt(C[2][2]*s2),
    radial:R.k, radialPct:R.k*Math.max(la,lb)**2*100, radialUsed:R.useK, spots:idx,nIndexed:idx.filter(s=>s.indexed).length,fracI:best.frac,occupancy:best.occ,sigmaPx:Math.sqrt(s2),tol}; };
/* 격자 다듬기: 관측 v = D(c + n1 a + n2 b), D(p) = p(1 + k|p|²) (투사 렌즈 방사 왜곡, opt.radial: "auto"|true|false).
   안쪽부터 지수를 붙이고 모형을 다듬으며 바깥으로 넓힌다 (왜곡으로 바깥 스폿이 몇 px 밀려도 놓치지 않게). 가중 √I. */
DP.refineLattice=function(vs,ws,a,b,tol,opt){ opt=opt||{}; const N=vs.length; const mode=opt.radial==null?"auto":opt.radial; let c=[0,0], k=0, useK=false; let ix=null, inv=null, s2=tol*tol, nn=0;
  const len0=Math.max(Math.hypot(...a),Math.hypot(...b)); const rAll=Math.max(...vs.map(v=>Math.hypot(...v)));
  const undist=(v,kk)=>{ const r2=v[0]*v[0]+v[1]*v[1]; const f=1/(1+kk*r2); return [v[0]*f,v[1]*f]; };
  const solve=(kk,rLim)=>{ const P=vs.map(v=>undist([v[0]-c[0],v[1]-c[1]],kk)).map(p=>[p[0],p[1]]); const I=indexAll(P,a,b,tol);
    const M=[[0,0,0],[0,0,0],[0,0,0]], bx=[0,0,0], by=[0,0,0]; let n=0;
    for(let i=0;i<N;i++){ if(!I[i].ok||Math.hypot(...vs[i])>rLim) continue; const r=[I[i].n[0],I[i].n[1],1]; const w2=ws[i]; const q=undist(vs[i],kk); n++;
      for(let u=0;u<3;u++){ for(let v=0;v<3;v++) M[u][v]+=w2*r[u]*r[v]; bx[u]+=w2*r[u]*q[0]; by[u]+=w2*r[u]*q[1]; } }
    if(n<3) return null; const iv=inv3(M); if(!iv) return null; const sx=mul3(iv,bx), sy=mul3(iv,by); return {a:[sx[0],sy[0]],b:[sx[1],sy[1]],c:[sx[2],sy[2]],inv:iv,n}; };
  const stats=(kk)=>{ // 모형 위치와의 잔차 (관측 공간), 방사 성분
    let ss=0, sw=0, n=0; const I=[]; let A=0,B=0;
    for(let i=0;i<N;i++){ const q=undist(vs[i],kk); const D=det2(a,b); const qq=[q[0]-c[0],q[1]-c[1]]; const n1=Math.round((qq[0]*b[1]-qq[1]*b[0])/D), n2=Math.round((a[0]*qq[1]-a[1]*qq[0])/D);
      const p=[c[0]+n1*a[0]+n2*b[0],c[1]+n1*a[1]+n2*b[1]]; const r2=p[0]*p[0]+p[1]*p[1]; const m=[p[0]*(1+kk*r2),p[1]*(1+kk*r2)]; const ex=vs[i][0]-m[0], ey=vs[i][1]-m[1]; const res=Math.hypot(ex,ey);
      const ok=res<=tol*(1+0.5*Math.sqrt(r2)/rAll)&&(n1||n2); I.push({n:[n1,n2],res,ok});
      if(ok){ ss+=ws[i]*(ex*ex+ey*ey); sw+=ws[i]; n++; const r=Math.sqrt(r2); if(r>0){ const dr=(ex*p[0]+ey*p[1])/r; const rk=r*r2; A+=ws[i]*dr*rk; B+=ws[i]*rk*rk; } } }
    return {I,ss,sw,n,dk:B>0?A/B:0}; };
  // 넓혀 가기: 안쪽 2.2 |R| 에서 시작해 1.5 배씩
  let rLim=Math.max(2.2*len0,1); for(let pass=0;pass<12;pass++){ const S=solve(useK||mode===true?k:0,rLim); if(!S) break; a=S.a; b=S.b; c=S.c; inv=S.inv;
    if(mode!==false){ const st=stats(k); k+=st.dk; } const done=rLim>=rAll; rLim=Math.min(rAll*1.01,rLim*1.5); if(mode===true) useK=true; if(done&&pass>3) break; }
  // 왜곡 항을 쓸지 (auto): 잔차가 뚜렷이 줄고(20 %) k 가 0 과 구별될 때
  const s0=(()=>{ const kk0=0; const S=solve(kk0,rAll*1.01); return S; })();
  const fitWith=kk=>{ const S=solve(kk,rAll*1.01); if(!S) return null; const sv={a,b,c}; a=S.a; b=S.b; c=S.c; const st=stats(kk); const out={S,st,s2:st.n>3?st.ss/st.sw*st.n/(2*(st.n-3)):tol*tol}; a=sv.a; b=sv.b; c=sv.c; return out; };
  let F0=fitWith(0), F1=mode!==false?fitWith(k):null; void s0;
  let F=F0; useK=false;
  if(F1&&F0&&(mode===true||(F1.st.n>=8&&F1.s2<0.64*F0.s2&&F1.st.n>=F0.st.n))) { F=F1; useK=true; }
  if(!F) return null; a=F.S.a; b=F.S.b; c=F.S.c;
  // 최종 기저를 REDI 규약(둔각, R3 = 짧은 대각선)으로 축약한 뒤 그 기저로 한 번 더 풀어 공분산을 그 기저에서 얻는다
  for(let it=0;it<4;it++){ const [ra,rb]=DP.reduce2(a,b); const same=Math.hypot(ra[0]-a[0],ra[1]-a[1])+Math.hypot(rb[0]-b[0],rb[1]-b[1])<1e-9*(1+Math.hypot(...a)); if(same&&it>0) break; a=ra; b=rb; const G=fitWith(useK?k:0); if(!G) break; F=G; a=F.S.a; b=F.S.b; c=F.S.c; }
  inv=F.S.inv; s2=F.s2; ix=F.st.I; nn=F.st.n; void nn;
  // 가중 공분산 정규화: Σw 로 나눈 상대 가중 → (XᵀWX)⁻¹ 에 Σw/n 을 곱해 단위 맞춤
  const scale=F.st.sw/Math.max(1,F.st.n); inv=inv.map(r=>r.map(x=>x*scale));
  return {a,b,c,k:useK?k:0,useK,ix,inv,s2}; };
function inv3(m){ const [a,b,c]=m[0],[d,e,f]=m[1],[g,h,i]=m[2]; const A=e*i-f*h, B=-(d*i-f*g), C=d*h-e*g; const det=a*A+b*B+c*C; if(Math.abs(det)<1e-12) return null;
  return [[A/det,-(b*i-c*h)/det,(b*f-c*e)/det],[B/det,(a*i-c*g)/det,-(a*f-c*d)/det],[C/det,-(a*h-b*g)/det,(a*e-b*d)/det]]; }
const mul3=(m,v)=>[m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2],m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2],m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]];

/* ---------- 링 ---------- */
/* 광선(nRay 개)을 따라 r ∈ [rmin,rmax] 에서 가장 강한 국소 봉우리 위치 → 점들 */
DP.ringPoints=function(img,w,h,c,rmin,rmax,opt){ opt=opt||{}; const nRay=opt.nRay||180; const mask=opt.mask; const pts=[];
  for(let k=0;k<nRay;k++){ const th=2*Math.PI*k/nRay, cx=Math.cos(th), sy=Math.sin(th); let best=-Infinity, br=-1; const prof=[];
    for(let r=rmin;r<=rmax;r+=0.5){ const x=c.x+r*cx, y=c.y+r*sy; const xi=Math.round(x), yi=Math.round(y); if(xi<0||yi<0||xi>=w||yi>=h){ prof.push(NaN); continue; } if(mask&&mask[yi*w+xi]){ prof.push(NaN); continue; } prof.push(bilin(img,w,h,x,y)); }
    // 국소 배경(선형) 뺀 최대
    const n=prof.length; if(n<5) continue; const p0=prof[0], p1=prof[n-1]; if(!(p0===p0&&p1===p1)) continue;
    for(let i=1;i<n-1;i++){ const v=prof[i]; if(v!==v) continue; const base=p0+(p1-p0)*i/(n-1); const q=v-base; if(q>best&&v>=prof[i-1]&&v>=prof[i+1]){ best=q; br=i; } }
    if(br>0&&br<n-1){ const a=prof[br-1], b=prof[br], cc=prof[br+1]; const den=a-2*b+cc; const off=den<0?0.5*(a-cc)/den:0; const r=rmin+(br+off)*0.5; pts.push({x:c.x+r*cx,y:c.y+r*sy,th,q:best}); } }
  return pts; };
function bilin(a,w,h,x,y){ const x0=Math.floor(x), y0=Math.floor(y); const fx=x-x0, fy=y-y0; const g=(xx,yy)=>a[Math.min(h-1,Math.max(0,yy))*w+Math.min(w-1,Math.max(0,xx))]; return (1-fx)*(1-fy)*g(x0,y0)+fx*(1-fy)*g(x0+1,y0)+(1-fx)*fy*g(x0,y0+1)+fx*fy*g(x0+1,y0+1); }
DP.bilin=bilin;
/* 타원 맞춤 (Fitzgibbon–Pilu–Fisher 1999 직접 최소제곱, 정규화) → {cx,cy,a,b,phi} (a ≥ b, phi = a 축 방향) */
DP.fitEllipse=function(pts){ if(pts.length<6) return null; let mx=0,my=0; for(const p of pts){ mx+=p.x; my+=p.y; } mx/=pts.length; my/=pts.length; let sc=0; for(const p of pts) sc+=Math.hypot(p.x-mx,p.y-my); sc/=pts.length; if(!(sc>0)) return null;
  const D=pts.map(p=>{ const x=(p.x-mx)/sc, y=(p.y-my)/sc; return [x*x,x*y,y*y,x,y,1]; });
  const S=Array.from({length:6},()=>new Float64Array(6)); for(const r of D) for(let i=0;i<6;i++) for(let j=0;j<6;j++) S[i][j]+=r[i]*r[j];
  // S = [S1 S2; S2' S3], 제약 4ac − b² = 1
  const S1=[[S[0][0],S[0][1],S[0][2]],[S[1][0],S[1][1],S[1][2]],[S[2][0],S[2][1],S[2][2]]], S2=[[S[0][3],S[0][4],S[0][5]],[S[1][3],S[1][4],S[1][5]],[S[2][3],S[2][4],S[2][5]]], S3=[[S[3][3],S[3][4],S[3][5]],[S[4][3],S[4][4],S[4][5]],[S[5][3],S[5][4],S[5][5]]];
  const S3i=inv3(S3); if(!S3i) return null; const T=mm(S3i,tr(S2)).map(r=>r.map(x=>-x)); const Mm=add(S1,mm(S2,T)); // M = C1⁻¹ (S1 + S2 T)
  const M=[[Mm[2][0]/2,Mm[2][1]/2,Mm[2][2]/2],[-Mm[1][0],-Mm[1][1],-Mm[1][2]],[Mm[0][0]/2,Mm[0][1]/2,Mm[0][2]/2]];
  const ev=eig3(M); let a1=null; for(const v of ev){ const cond=4*v[0]*v[2]-v[1]*v[1]; if(cond>0){ a1=v; break; } } if(!a1) return null;
  const a2=mul3(T,a1); let [A,B,C]=a1, [Dd,E,F]=a2;
  // 정규화 풀기: 좌표 (x−mx)/sc
  const cxn=(2*C*Dd-B*E)/(B*B-4*A*C), cyn=(2*A*E-B*Dd)/(B*B-4*A*C);
  const num=2*(A*E*E+C*Dd*Dd-B*Dd*E+(B*B-4*A*C)*F); const s=Math.sqrt((A-C)**2+B*B);
  const ra=-Math.sqrt(num*(A+C+s))/(B*B-4*A*C), rb=-Math.sqrt(num*(A+C-s))/(B*B-4*A*C);
  let phi=Math.atan2(-B,C-A)/2; let ax=ra, bx=rb; if(!(ax>0&&bx>0)) return null; if(bx>ax){ const t=ax; ax=bx; bx=t; phi+=Math.PI/2; }
  return {cx:mx+cxn*sc,cy:my+cyn*sc,a:ax*sc,b:bx*sc,phi:((phi%Math.PI)+Math.PI)%Math.PI}; };
const tr=m=>m[0].map((_,j)=>m.map(r=>r[j])); const mm=(a,b)=>a.map(r=>b[0].map((_,j)=>r.reduce((s,x,k)=>s+x*b[k][j],0))); const add=(a,b)=>a.map((r,i)=>r.map((x,j)=>x+b[i][j]));
/* 3×3 실수 행렬 고유벡터 (거듭 제곱 + 감축이 아니라: 특성 다항식 근 → 영공간) */
function eig3(M){ const a=M; const p=[-(a[0][0]+a[1][1]+a[2][2]), a[0][0]*a[1][1]-a[0][1]*a[1][0]+a[0][0]*a[2][2]-a[0][2]*a[2][0]+a[1][1]*a[2][2]-a[1][2]*a[2][1], -(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))];
  const roots=cubicRoots(1,p[0],p[1],p[2]); const out=[];
  for(const l of roots){ const m=[[a[0][0]-l,a[0][1],a[0][2]],[a[1][0],a[1][1]-l,a[1][2]],[a[2][0],a[2][1],a[2][2]-l]];
    const c1=[m[0][1]*m[1][2]-m[0][2]*m[1][1],m[0][2]*m[1][0]-m[0][0]*m[1][2],m[0][0]*m[1][1]-m[0][1]*m[1][0]], c2=[m[0][1]*m[2][2]-m[0][2]*m[2][1],m[0][2]*m[2][0]-m[0][0]*m[2][2],m[0][0]*m[2][1]-m[0][1]*m[2][0]], c3=[m[1][1]*m[2][2]-m[1][2]*m[2][1],m[1][2]*m[2][0]-m[1][0]*m[2][2],m[1][0]*m[2][1]-m[1][1]*m[2][0]];
    const v=[c1,c2,c3].sort((x,y)=>Math.hypot(...y)-Math.hypot(...x))[0]; const n=Math.hypot(...v); if(n>0) out.push(v.map(x=>x/n)); }
  return out; }
function cubicRoots(a,b,c,d){ b/=a; c/=a; d/=a; const q=(3*c-b*b)/9, r=(9*b*c-27*d-2*b*b*b)/54; const disc=q*q*q+r*r; const out=[];
  if(disc>0){ const s=Math.cbrt(r+Math.sqrt(disc)), t=Math.cbrt(r-Math.sqrt(disc)); out.push(-b/3+s+t); }
  else { const th=Math.acos(Math.max(-1,Math.min(1,r/Math.sqrt(-q*q*q)))); const m=2*Math.sqrt(-q); out.push(m*Math.cos(th/3)-b/3,m*Math.cos((th+2*Math.PI)/3)-b/3,m*Math.cos((th+4*Math.PI)/3)-b/3); }
  return out; }
/* 타원 보정 방위 적분: 화소를 원으로 펴서(장축 방향으로 b/a 배) 반지름 칸(0.5 px)의 평균 */
DP.azimuthal=function(img,w,h,el,opt){ opt=opt||{}; const mask=opt.mask; const step=opt.step||0.5; const cs=Math.cos(el.phi), sn=Math.sin(el.phi); const k=el.b/el.a; // 장축 성분을 k 배 → 반지름 = 단축 기준
  const rmax=opt.rmax||Math.hypot(w,h); const nb=Math.ceil(rmax/step)+1; const sum=new Float64Array(nb), cnt=new Float64Array(nb);
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const i=y*w+x; if(mask&&mask[i]) continue; const dx=x-el.cx, dy=y-el.cy; const u=dx*cs+dy*sn, v=-dx*sn+dy*cs; const r=Math.hypot(u*k,v); const b=Math.round(r/step); if(b<nb){ sum[b]+=img[i]; cnt[b]++; } }
  const r=[], I=[], n=[]; for(let b=0;b<nb;b++){ if(cnt[b]<3) continue; r.push(b*step); I.push(sum[b]/cnt[b]); n.push(cnt[b]); } return {r,I,n,radiusIs:"minor"}; };
/* 1D 배경: 아래 껍질 (구간 최솟값 + 다듬기) */
DP.background1d=function(y,win){ const n=y.length; const m=new Float64Array(n); const W=Math.max(3,win|0); for(let i=0;i<n;i++){ let v=Infinity; for(let j=Math.max(0,i-W);j<=Math.min(n-1,i+W);j++) if(y[j]<v) v=y[j]; m[i]=v; }
  const s=new Float64Array(n); for(let i=0;i<n;i++){ let a=0,c=0; for(let j=Math.max(0,i-W);j<=Math.min(n-1,i+W);j++){ a+=m[j]; c++; } s[i]=a/c; } for(let i=0;i<n;i++) if(s[i]>y[i]) s[i]=y[i]; return s; };
/* 봉우리: 배경 뺀 프로파일의 국소 최대 (돌출 ≥ k·잡음) + 포물선 부화소 */
DP.peaks1d=function(r,y,bg,opt){ opt=opt||{}; const n=y.length; const d=y.map((v,i)=>v-bg[i]); const diffs=[]; for(let i=1;i<n;i++) diffs.push(Math.abs(d[i]-d[i-1])); diffs.sort((a,b)=>a-b); const noise=(diffs[diffs.length>>1]||1e-9)*1.4826/Math.SQRT2;
  const out=[]; const W=opt.win||4;
  for(let i=W;i<n-W;i++){ let isMax=true; for(let j=i-W;j<=i+W;j++) if(j!==i&&d[j]>d[i]){ isMax=false; break; } if(!isMax) continue; if(d[i]<(opt.k||5)*noise) continue;
    const a=d[i-1], b=d[i], c=d[i+1]; const den=a-2*b+c; const off=den<0?0.5*(a-c)/den:0; const rr=r[i]+off*(r[i+1]-r[i]);
    // 반높이 폭
    let lo=i, hi=i; while(lo>0&&d[lo]>b/2) lo--; while(hi<n-1&&d[hi]>b/2) hi++; out.push({r:rr,h:b,fwhm:r[hi]-r[lo],snr:b/noise}); }
  return {peaks:out,noise}; };

/* 링 중심 자동: 줄인 영상(≤ 512)에서 방위 평균 프로파일이 가장 날카로운 점 (Σ (ΔI)², 빈 1 px) — 거칠게 → 곱게.
   대칭인 링 패턴은 참 중심에서 봉우리가 가장 좁고 높다. 빔스토퍼·덧그림은 mask 로 뺀다 */
DP.ringCenter=function(img,w,h,opt){ opt=opt||{}; const f=Math.max(1,Math.ceil(Math.max(w,h)/(opt.work||384))); const W=Math.floor(w/f), H=Math.floor(h/f); const a=new Float32Array(W*H), m=new Uint8Array(W*H);
  for(let y=0;y<H;y++) for(let x=0;x<W;x++){ let s=0, bad=0; for(let j=0;j<f;j++) for(let i=0;i<f;i++){ const k=(y*f+j)*w+x*f+i; s+=img[k]; if(opt.mask&&opt.mask[k]) bad=1; } a[y*W+x]=s/(f*f); m[y*W+x]=bad; }
  const rmax=Math.floor(Math.hypot(W,H)/2); const sum=new Float64Array(rmax+2), cnt=new Float64Array(rmax+2);
  const obj=(cx,cy)=>{ sum.fill(0); cnt.fill(0); for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const k=y*W+x; if(m[k]) continue; const r=Math.hypot(x-cx,y-cy); const b=Math.floor(r); if(b>rmax) continue; const t=r-b; sum[b]+=a[k]*(1-t); cnt[b]+=1-t; sum[b+1]+=a[k]*t; cnt[b+1]+=t; }
    let o=0, prev=null; const r0=opt.rmin?opt.rmin/f:3; for(let b=Math.ceil(r0);b<=rmax;b++){ if(cnt[b]<8) { prev=null; continue; } const v=sum[b]/cnt[b]; if(prev!=null) o+=(v-prev)*(v-prev); prev=v; } return o; };
  let best={x:opt.x0!=null?opt.x0/f:W/2,y:opt.y0!=null?opt.y0/f:H/2}; best.o=obj(best.x,best.y); const span=opt.span!=null?opt.span/f:Math.min(W,H)/4;
  // 거친 격자(span/6 간격) 한 번, 그다음 간격을 반씩 줄이며 ±2 간격 — 줄인 영상 0.5 px 까지 (원래 영상 f/2 px; 나머지는 타원 맞춤이 다듬는다)
  const st0=Math.max(1,Math.round(span/6)); for(let st=st0;;){ const cx=best.x, cy=best.y; const R=st===st0?span:st*2;
    for(let dy=-R;dy<=R+1e-9;dy+=st) for(let dx=-R;dx<=R+1e-9;dx+=st){ if(!dx&&!dy&&st!==st0) continue; const o=obj(cx+dx,cy+dy); if(o>best.o) best={x:cx+dx,y:cy+dy,o}; }
    if(st<=0.5) break; st=st/2; }
  return {x:best.x*f+(f-1)/2,y:best.y*f+(f-1)/2,factor:f}; };
/* 링 패턴 한 번에: 중심 → 가장 센 링들로 타원 → 타원 보정 방위 적분 → 배경 → 봉우리 */
DP.analyzeRings=function(img,w,h,opt){ opt=opt||{}; const mask=opt.mask; const c0=opt.center||DP.ringCenter(img,w,h,{mask,rmin:opt.rmin});
  let el={cx:c0.x,cy:c0.y,a:1,b:1,phi:0}; const rmax=opt.rmax||Math.min(Math.max(c0.x,w-c0.x),Math.max(c0.y,h-c0.y));
  let az=DP.azimuthal(img,w,h,el,{mask,rmax}); let bg=DP.background1d(az.I,Math.max(6,Math.round(rmax/40))); let pk=DP.peaks1d(az.r,az.I,bg,{k:opt.k||6,win:Math.max(3,Math.round(rmax/200))}).peaks.filter(p=>p.r>(opt.rmin||rmax*0.05));
  const fits=[];
  for(let it=0;it<3;it++){ if(!pk.length) break; const strong=pk.slice().sort((p,q)=>q.h*Math.sqrt(q.r)-p.h*Math.sqrt(p.r)).slice(0,opt.nRings||4); const all=[]; fits.length=0;
    for(const p of strong){ const win=Math.max(3*p.fwhm,0.04*p.r); const pts=DP.ringPoints(img,w,h,{x:el.cx,y:el.cy},Math.max(1,p.r*(el.b?1:1)-win),p.r*(el.a/el.b)+win,{mask,nRay:opt.nRay||360});
      // 광선 점 중 튀는 것 빼기 (반지름 중앙값 ± 3 MAD)
      const rs=pts.map(q=>Math.hypot(q.x-el.cx,q.y-el.cy)); const md=rs.slice().sort((x,y)=>x-y)[rs.length>>1]; const mad=rs.map(r=>Math.abs(r-md)).sort((x,y)=>x-y)[rs.length>>1]*1.4826||1; const keep=pts.filter((q,i)=>Math.abs(rs[i]-md)<=3*mad+0.5);
      const e=DP.fitEllipse(keep); if(e) { fits.push(Object.assign({r:p.r,n:keep.length},e)); all.push(...keep.map(q=>Object.assign({ring:fits.length-1},q))); } }
    if(!fits.length) break;
    // 링마다 중심·찌그러짐이 같아야 한다: 가중(점 수) 평균
    let W=0,cx=0,cy=0,e2c=0,e2s=0; for(const f of fits){ const wgt=f.n; W+=wgt; cx+=wgt*f.cx; cy+=wgt*f.cy; const ecc=f.a/f.b-1; e2c+=wgt*ecc*Math.cos(2*f.phi); e2s+=wgt*ecc*Math.sin(2*f.phi); }
    const ecc=Math.hypot(e2c,e2s)/W, phi=((Math.atan2(e2s,e2c)/2)%Math.PI+Math.PI)%Math.PI; el={cx:cx/W,cy:cy/W,a:1+ecc,b:1,phi};
    az=DP.azimuthal(img,w,h,el,{mask,rmax}); bg=DP.background1d(az.I,Math.max(6,Math.round(rmax/40))); pk=DP.peaks1d(az.r,az.I,bg,{k:opt.k||6,win:Math.max(3,Math.round(rmax/200))}).peaks.filter(p=>p.r>(opt.rmin||rmax*0.05)); }
  // 링마다 흩어짐 → 불확도
  const sd=v=>{ const m=v.reduce((a,b)=>a+b,0)/v.length; return Math.sqrt(v.reduce((a,b)=>a+(b-m)**2,0)/Math.max(1,v.length-1)); };
  const unc=fits.length>1?{cx:sd(fits.map(f=>f.cx))/Math.sqrt(fits.length),cy:sd(fits.map(f=>f.cy))/Math.sqrt(fits.length),ratio:sd(fits.map(f=>f.a/f.b))/Math.sqrt(fits.length)}:null;
  return {ellipse:el,ratio:el.a/el.b,center:{x:el.cx,y:el.cy},profile:az,background:bg,peaks:pk,ringFits:fits,unc,coarse:c0}; };

/* ---------- 스케일바 ---------- */
/* findOverlay 의 막대 중 가장 긴 것(아래쪽 우선). 길이 = 바깥 끝 사이(px). 값은 사용자가 확인 */
DP.scaleBar=function(ov,w,h){ if(!ov||!ov.bars.length) return null; const b=ov.bars.slice().sort((p,q)=>(q.y0>h*0.5)-(p.y0>h*0.5)||q.len-p.len)[0];
  return {x0:b.x0,x1:b.x1,y0:b.y0,y1:b.y1,lenPx:b.hollow?b.len-b.thick:b.len,outerPx:b.len,hollow:b.hollow}; };

if(typeof module!=="undefined"&&module.exports) module.exports=DP;
G.DIFF_PATTERN=DP;
})();
