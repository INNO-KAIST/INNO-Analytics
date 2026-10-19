/* ============================================================
   Auto particles (2026-10-19) — 나노입자 TEM 이미지 한 번에 마스크로. 배율(입자 크기)에 맞춰 스스로 맞춘다.
   ① 1 차: Gaussian σ 2 → Default 문턱(어두운 쪽이 > 50 % 면 밝은 입자) → 성분 넓이의 '넓이 가중' 중앙값 = 입자 하나의 대략 넓이
      (개수 중앙값은 잡음 점 수천 개가 차지해 틀림 — 넓이로 가중하면 실제 입자가 차지)
   ② 2 차: σ = clamp(r/6, 2, 12) (r = 그 넓이의 반지름) 로 다시 흐림 → 같은 방법 문턱 → Fill Holes → Watershed → Split by size
      배율이 높아 입자가 크면 σ 도 커져 지지막(carbon grid) 결·잡음이 사라진다
   ③ Size 최솟값 = 홑입자 넓이 A0 의 ¼ (A0 = splitBySize 가 다시 잰 값)
   반환 {mask(plane), info:{sigma1, sigma, dark, Aw, A0, r0, minPx, split, parts}}
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis; const IC=G.IMG_CORE;
function comps(bin,W,H){ return IC.componentAreas(bin,W,H); }
function wMedian(a){ const v=a.slice().sort((x,y)=>x-y); let tot=0; for(const x of v) tot+=x; let c=0; for(const x of v){ c+=x; if(c>=tot/2) return x; } return v[v.length-1]||0; }
function mask(p,sigma,dark){ const g=IC.gaussian(p,{sigma}); const gp=g.plane||g; const t=IC.autoThreshold(gp,"Default",dark); return IC.applyThreshold(gp,t.lower,t.upper); }
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
  let m=mask(p,sigma,dark); m=IC.fillHoles(m); m=IC.watershed(m);
  // 작은 조각(¼ 입자보다 작은 것)은 홑입자 크기 추정에서 빼도록 A0 를 2 차 넓이 가중 중앙값 근처 홑입자에서
  const s=IC.splitBySize(m,{minFrac:0.15,gray:p,dark:!dark}); // 문턱 'dark' = 바탕이 어두움 → 입자는 밝음. splitBySize 의 dark = 입자가 어두움
  const A0=s.info.A0>0?s.info.A0:Aw; const minPx=Math.max(9,Math.round(0.25*A0));
  return {mask:{w:W,h:H,bit:8,data:s.data},info:{sigma1:2,sigma,dark,Aw,A1,A0,r0:Math.sqrt(A0/Math.PI),minPx,split:s.info.split,parts:s.info.parts}};
};
})();
