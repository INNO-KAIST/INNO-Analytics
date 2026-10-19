// Auto particle mask 시험 (node test/test_autopart.js) — 2026-10-19 사용자 보고 두 가지를 정답 있는 합성으로:
//  ① 겹친 덩어리(중심 거리 1.3–1.9 R, 2·3·4 개 묶음)를 하나로 셈  ② 고배율(입자 R 60 px)에서 지지막 결(상관 잡음)을 입자로 셈
"use strict";
const IC=require("./load_core.js");
let pass=0, fail=0; const ok=(c,m,d)=>{ c?pass++:fail++; console.log((c?"  ✓ ":"  ✗ ")+m+(d!=null?" — "+d:"")); };
function rng(seed){ let s=seed>>>0; return ()=>((s=(s*1664525+1013904223)>>>0)/4294967296); }
function gauss(r){ return ()=>{ const u=r()||1e-9, v=r(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); }; }
function scene(W,H,R,seed,{clusters,singles,noise,corr}){ const r=rng(seed), g=gauss(r); const C=[]; const img=new Float32Array(W*H).fill(175);
  const free=(x,y,rr)=>x>rr+20&&y>rr+20&&x<W-rr-20&&y<H-rr-20&&C.every(c=>(c[0]-x)**2+(c[1]-y)**2>(2*R+8)**2);
  for(const k of clusters){ for(let tries=0;tries<500;tries++){ const x=r()*W, y=r()*H; if(!free(x,y,3*R)) continue; const grp=[[x,y]];
      while(grp.length<k){ const b=grp[Math.floor(r()*grp.length)], a=r()*2*Math.PI, d=(1.3+0.6*r())*R; const q=[b[0]+d*Math.cos(a),b[1]+d*Math.sin(a)]; if(grp.every(o=>(o[0]-q[0])**2+(o[1]-q[1])**2>(1.25*R)**2)) grp.push(q); }
      C.push(...grp); break; } }
  for(let n=0;n<singles;n++) for(let tries=0;tries<500;tries++){ const x=r()*W, y=r()*H; if(free(x,y,R)){ C.push([x,y]); break; } }
  for(const [cx,cy] of C) for(let y=Math.max(0,Math.floor(cy-R-1));y<=Math.min(H-1,cy+R+1);y++) for(let x=Math.max(0,Math.floor(cx-R-1));x<=Math.min(W-1,cx+R+1);x++){ const d=Math.hypot(x+0.5-cx,y+0.5-cy); if(d<=R) img[y*W+x]-=60*Math.sqrt(1-(d/R)**2)+40; } // 두께 ∝ 현의 길이(겹치면 더 어두움)
  // 잡음: 화소 잡음 + (선택) 상관 잡음 = 지지막 결 (σ corr px 로 흐린 잡음)
  const nz=new Float32Array(W*H); for(let i=0;i<nz.length;i++) nz[i]=g(); let cz=null; if(corr){ const q=IC.gaussian({w:W,h:H,bit:32,data:Float32Array.from(nz.map(()=>g()))},{sigma:corr}); cz=(q.plane||q).data; let s=0; for(const v of cz) s+=v*v; const k=Math.sqrt(cz.length/s); for(let i=0;i<cz.length;i++) cz[i]*=k; }
  for(let i=0;i<img.length;i++) img[i]=Math.max(0,Math.min(255,Math.round(img[i]+noise*nz[i]+(cz?noise*1.2*cz[i]:0))));
  return {p:{w:W,h:H,bit:8,data:img},n:C.length,nCl:C.length-singles}; }
function count(p){ const r=IC.autoParticles(p); const bin=new Uint8Array(r.mask.data.length); for(let i=0;i<bin.length;i++) bin[i]=r.mask.data[i]>0?1:0;
  const a=IC.analyzeParticles(bin,p.w,p.h,null,{pw:1,ph:1},{minSize:r.info.minPx,minCirc:0.5,excludeEdges:false}); const e=a.particles.map(q=>q.meas.ECD).sort((x,y)=>x-y);
  const ov=a.particles.map(q=>r.ovl[Math.floor(q.meas.Y)*p.w+Math.floor(q.meas.X)]?1:0); const eUsed=a.particles.filter((q,i)=>!ov[i]).map(q=>q.meas.ECD), eOv=a.particles.filter((q,i)=>ov[i]).map(q=>q.meas.ECD);
  return {n:a.particles.length,info:r.info,e,nOv:ov.reduce((x,y)=>x+y,0),eUsed,eOv}; }
const avg=a=>a.reduce((x,y)=>x+y,0)/a.length;
console.log("== 1. 겹친 덩어리 (R 20 px, 잡음 σ 25)");
{ const s=scene(1400,1400,20,3,{clusters:[2,2,2,2,3,3,3,3,4,4,4,4,2,3,4,5,5],singles:120,noise:25}); const t0=Date.now(); const c=count(s.p);
  ok(Math.abs(c.n-s.n)<=Math.max(3,0.03*s.n),`개수 ${c.n} ↔ 정답 ${s.n} (± 3 %)`,`${Date.now()-t0} ms, 덩어리 ${c.info.split} 개 나눔`);
  const big=c.e.filter(x=>x>1.4*2*20).length; ok(big===0,"지름 1.4 배 넘는 '합쳐진 입자' 없음",big);
  ok(Math.abs(c.nOv-s.nCl)<=3,`겹친 입자 표시 ${c.nOv} ↔ 정답 ${s.nCl} (± 3) — 덩어리 ${c.info.ovlClumps} 개`);
  ok(avg(c.eUsed)>avg(c.eOv),`겹침을 뺀 ECD 평균 ${avg(c.eUsed).toFixed(1)} px > 겹친 조각 평균 ${avg(c.eOv).toFixed(1)} px (조각은 작게 잼)`);
  const sdv=a=>{ const m=avg(a); return Math.sqrt(a.reduce((x,y)=>x+(y-m)**2,0)/(a.length-1)); };
  ok(sdv(c.eUsed)<sdv(c.e),`겹침을 빼면 ECD 퍼짐이 줄어듦 — SD ${sdv(c.eUsed).toFixed(2)} < 전체 ${sdv(c.e).toFixed(2)} px`); }
console.log("== 2. 고배율 (R 60 px) + 지지막 결(상관 잡음 σ 3 px)");
{ const s=scene(1600,1600,60,7,{clusters:[2,3,2],singles:28,noise:22,corr:3}); const t0=Date.now(); const c=count(s.p);
  ok(Math.abs(c.n-s.n)<=2,`개수 ${c.n} ↔ 정답 ${s.n} (± 2)`,`${Date.now()-t0} ms, σ ${c.info.sigma.toFixed(1)} px`);
  ok(c.e[0]>0.6*120,`가장 작은 입자 ECD ${c.e[0].toFixed(0)} px > 0.6 × 120 (잡음 조각 없음)`); ok(c.info.sigma>=6,"흐림이 입자 크기에 맞춰 커짐",c.info.sigma.toFixed(1));
  ok(Math.abs(c.nOv-s.nCl)<=1,`겹친 입자 표시 ${c.nOv} ↔ 정답 ${s.nCl}`); }
console.log("== 3. 저배율 (R 6 px) — 흐림이 너무 커서 입자를 지우지 않음");
{ const s=scene(900,900,6,11,{clusters:[2,3],singles:300,noise:18}); const c=count(s.p); ok(Math.abs(c.n-s.n)<=Math.max(4,0.04*s.n),`개수 ${c.n} ↔ 정답 ${s.n} (± 4 %)`,`σ ${c.info.sigma.toFixed(1)}`); }
console.log("== 4. Find (auto) — 날 이미지 prominence 10 이면 잡음 극대 수만 개 → 자동은 입자 조각마다 중심 하나");
for(const [R,seed,cfg,tol] of [[20,3,{clusters:[2,2,3,3,4],singles:120,noise:25},0.03],[60,7,{clusters:[2,3,2],singles:28,noise:22,corr:3},0.03],[6,11,{clusters:[2,3],singles:300,noise:18},0.04]]){
  const W=R>=60?1600:R>=20?1400:900; const sc=scene(W,W,R,seed,cfg); const raw=IC.findMaxima(sc.p,10,{}).points.length; const t0=Date.now(); const r=IC.findMaximaAuto(sc.p);
  ok(Math.abs(r.points.length-sc.n)<=Math.max(3,tol*sc.n),`R ${R}: 자동 ${r.points.length} ↔ 정답 ${sc.n} (날 이미지 prominence 10: ${raw})`,`${Date.now()-t0} ms, light ${r.info.light}, σ ${r.info.sigma.toFixed(1)}`);
  ok(r.info.light===true,"밝은 바탕·어두운 입자 → Light background 자동"); }
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail) process.exit(1);
