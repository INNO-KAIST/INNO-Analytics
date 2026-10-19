/* ============================================================
   Auto particles (2026-10-19) — 나노입자 TEM 이미지 한 번에 마스크로. 배율(입자 크기)에 맞춰 스스로 맞춘다.
   ① 1 차: Gaussian σ 2 → Default 문턱(어두운 쪽이 > 50 % 면 밝은 입자) → 성분 넓이의 '넓이 가중' 중앙값 = 입자 하나의 대략 넓이
      (개수 중앙값은 잡음 점 수천 개가 차지해 틀림 — 넓이로 가중하면 실제 입자가 차지)
   ② 2 차: σ = clamp(r/6, 2, 12) (r = 그 넓이의 반지름) 로 다시 흐림 → 같은 방법 문턱 → Fill Holes → Watershed → Split by size
      배율이 높아 입자가 크면 σ 도 커져 지지막(carbon grid) 결·잡음이 사라진다
   ③ Size 최솟값 = 홑입자 넓이 A0 의 ¼ (A0 = splitBySize 가 다시 잰 값)
   ④ 겹침 표시(2026-10-19): Watershed·Split 전 마스크에서 한 덩어리였는데 나눈 뒤 (minPx 이상) 조각이 2 개 이상이면 그 덩어리 화소 = 1
      → 겹쳐 붙은 입자는 '찾되' 크기 분포에서 뺄 수 있게 (경계가 잘려 넓이가 작게 나오고, 겹친 부분이 가려져 크기를 믿기 어려움)
   반환 {mask(plane), ovl(Uint8Array, 1 = 겹친 덩어리), info:{sigma1, sigma, dark, Aw, A0, r0, minPx, split, parts, ovlClumps}}
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis; const IC=G.IMG_CORE;
function comps(bin,W,H){ return IC.componentAreas(bin,W,H); }
function wMedian(a){ const v=a.slice().sort((x,y)=>x-y); let tot=0; for(const x of v) tot+=x; let c=0; for(const x of v){ c+=x; if(c>=tot/2) return x; } return v[v.length-1]||0; }
function mask(p,sigma,dark){ const g=IC.gaussian(p,{sigma}); const gp=g.plane||g; const t=IC.autoThreshold(gp,"Default",dark); return IC.applyThreshold(gp,t.lower,t.upper); }
// 8-연결 성분 번호
function label(bin,W,H){ const lab=new Int32Array(W*H).fill(-1), st=new Int32Array(W*H); const size=[];
  for(let i0=0;i0<W*H;i0++){ if(!bin[i0]||lab[i0]>=0) continue; const id=size.length; let sp=0, n=0; st[sp++]=i0; lab[i0]=id;
    while(sp){ const i=st[--sp]; n++; const x=i%W, y=(i-x)/W; for(let dy=-1;dy<=1;dy++){ const yy=y+dy; if(yy<0||yy>=H) continue; for(let dx=-1;dx<=1;dx++){ const xx=x+dx; if(xx<0||xx>=W) continue; const j=yy*W+xx; if(bin[j]&&lab[j]<0){ lab[j]=id; st[sp++]=j; } } } }
    size.push(n); } return {lab,size}; }
// 나누기 전 덩어리(pre) 하나에 나눈 뒤 조각(post, minPx 이상)이 2 개 이상이면 그 덩어리 = 겹침
function overlapMask(pre,post,W,H,minPx){ const A=label(pre,W,H); const pb=new Uint8Array(W*H); for(let i=0;i<pb.length;i++) pb[i]=post[i]>0?1:0; const B=label(pb,W,H);
  const first=new Int32Array(B.size.length).fill(-1); for(let i=0;i<W*H;i++){ const b=B.lab[i]; if(b>=0&&first[b]<0) first[b]=i; }
  const cnt=new Int32Array(A.size.length); for(let b=0;b<B.size.length;b++){ if(B.size[b]<minPx) continue; const a=A.lab[first[b]]; if(a>=0) cnt[a]++; }
  const ovl=new Uint8Array(W*H); let clumps=0; for(let a=0;a<cnt.length;a++) if(cnt[a]>=2) clumps++;
  for(let i=0;i<W*H;i++){ const a=A.lab[i]; if(a>=0&&cnt[a]>=2) ovl[i]=1; } return {ovl,clumps}; }
IC.overlapMask=overlapMask;
/* Find Maxima (auto) — 입자 중심(2026-10-19). 날 이미지에 prominence 10 이면 TEM 잡음 점 수십만 개가 극대가 된다.
   흐린 사본의 세기 극대만으로는 맞닿은 쌍이 하나로 합쳐져(합성 421 개 중 391) Auto particle mask(배율 맞춤 흐림 → 문턱 → Watershed → Split by size)의
   조각마다 무게중심 하나 — Analyze particles (auto) 와 같은 개수. 원자 열처럼 '세기 극대'가 필요한 경우는 Find(직접)를 쓴다.
   반환 {points:[[x,y]], info:{light, r, sigma, minPx}} */
IC.findMaximaAuto=function(p,o){ o=o||{}; const a=IC.autoParticles(p,o); const W=p.w, H=p.h, d=a.mask.data, N=W*H; const lab=new Int32Array(N).fill(-1), st=new Int32Array(N); const pts=[];
  for(let i0=0;i0<N;i0++){ if(!(d[i0]>0)||lab[i0]>=0) continue; let sp=0, n=0, sx=0, sy=0; st[sp++]=i0; lab[i0]=1;
    while(sp){ const i=st[--sp]; n++; const x=i%W, y=(i-x)/W; sx+=x; sy+=y; for(let dy=-1;dy<=1;dy++){ const yy=y+dy; if(yy<0||yy>=H) continue; for(let dx=-1;dx<=1;dx++){ const xx=x+dx; if(xx<0||xx>=W) continue; const j=yy*W+xx; if(d[j]>0&&lab[j]<0){ lab[j]=1; st[sp++]=j; } } } }
    if(n>=a.info.minPx) pts.push([Math.min(W-1,Math.round(sx/n)),Math.min(H-1,Math.round(sy/n))]); }
  return {points:pts,info:{light:!a.info.dark,r:a.info.r0,sigma:a.info.sigma,minPx:a.info.minPx}};
};
IC.autoParticles=function(p,o){ o=o||{}; const W=p.w, H=p.h;
  // ① 극성·대략 크기
  const g1=IC.gaussian(p,{sigma:2}); const gp1=g1.plane||g1; const t=IC.autoThreshold(gp1,"Default",false); let n=0; for(let i=0;i<gp1.data.length;i++){ const v=gp1.data[i]; if(v>=t.lower&&v<=t.upper) n++; }
  const dark=o.dark!=null?!!o.dark:n/gp1.data.length>0.5;
  const m1=IC.applyThreshold(gp1,...(()=>{ const q=IC.autoThreshold(gp1,"Default",dark); return [q.lower,q.upper]; })());
  const b1=new Uint8Array(W*H); for(let i=0;i<b1.length;i++) b1[i]=m1.data[i]>0?1:0;
  const ar=comps(b1,W,H).filter(a=>a>=9&&a<0.25*W*H); const Aw=ar.length?wMedian(ar):100;
  // 홑입자 넓이: 넓이 가중 중앙값의 15 % 이상이고 볼록(solidity ≥ 0.93)한 성분의 중앙값 — 덩어리(가중 중앙값을 키움)·잡음(개수 중앙값을 줄임) 둘 다 피함
  const s1=IC.splitBySize(IC.fillHoles(m1),{minFrac:0.15}); const A1=s1.info.A0>0?s1.info.A0:Aw; const r=Math.sqrt(A1/Math.PI);
  // ② 크기에 맞춘 흐림
  const sigma=o.sigma>0?o.sigma:Math.max(2,Math.min(12,r/6));
  let m=mask(p,sigma,dark); m=IC.fillHoles(m); const pre=new Uint8Array(W*H); for(let i=0;i<pre.length;i++) pre[i]=m.data[i]>0?1:0; m=IC.watershed(m);
  // 작은 조각(¼ 입자보다 작은 것)은 홑입자 크기 추정에서 빼도록 A0 를 2 차 넓이 가중 중앙값 근처 홑입자에서
  const s=IC.splitBySize(m,{minFrac:0.15,gray:p,dark:!dark}); // 문턱 'dark' = 바탕이 어두움 → 입자는 밝음. splitBySize 의 dark = 입자가 어두움
  const A0=s.info.A0>0?s.info.A0:Aw; const minPx=Math.max(9,Math.round(0.25*A0));
  const ov=overlapMask(pre,s.data,W,H,minPx);
  return {mask:{w:W,h:H,bit:8,data:s.data},ovl:ov.ovl,info:{sigma1:2,sigma,dark,Aw,A1,A0,r0:Math.sqrt(A0/Math.PI),minPx,split:s.info.split,parts:s.info.parts,ovlClumps:ov.clumps}};
};
})();
