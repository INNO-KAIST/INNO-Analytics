// 지수 붙이기·상 찾기·보정 대조 — 정답을 아는 합성 회절도 (test/synth.js)
const fs=require("fs"), path=require("path");
const {rng,poisson,PH,synthSpots,DC}=require("./synth.js"); const DP=require("../src/diff_pattern.js"); const DI=require("../src/diff_index.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };
const SEED=+(process.env.SEED||20261017);
const phases=fs.readdirSync(path.join(__dirname,"cif")).map(f=>{ const p=DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif",f),"utf8"))[0]; p.name=f.replace(".cif",""); p.id="t:"+p.name; return p; });
const base={w:512,h:512,phase:PH.Si,gmax:14,sig:1.6,peak:900,g0:12,bg:6,direct:true,stop:true,stopR:14,stopAng:2.2,text:false};
function analyze(o){ const {img,truth}=synthSpots(o); const bs=DP.findBeamStop(img,o.w,o.h); const S=DP.findSpots(img,o.w,o.h,{mask:bs.mask}); const c=DP.findCenter(S.spots,o.w,o.h,{beamStop:bs.mask}); const L=DP.fitLattice(S.spots,c,{size:o.w}); return {img,truth,bs,S,c,L}; }
const sameZone=(a,b)=>{ const x=DC.reduceInt(a.map(Math.abs)).sort((p,q)=>p-q).join(), y=DC.reduceInt(b.map(Math.abs)).sort((p,q)=>p-q).join(); return x===y; };

console.log("== 상 찾기: Si [011] (보정 있음 / 없음)");
for(const [k,z] of [[0,[0,1,1]],[1,[1,1,2]],[2,[0,1,3]]]){
  const o=Object.assign({},base,{seed:SEED+k,zone:z,rot:0.4+k,cx:250.3+k,cy:261.7-k,scale:0.04}); const A=analyze(o);
  if(!A.L){ ok(false,`[${z}] 격자 맞춤 실패`); continue; }
  const R=DI.rankSpots(A.L,phases,{scale:o.scale,scaleErr:0.03,distortion:0.01,maxIdx:5,gmax:15,w:o.w,h:o.h,mask:A.bs.mask}); const t=R.cands[0];
  ok(t&&t.phaseName==="Si_Fd-3m"&&sameZone(t.zone,z),`[${z}] 보정 있음: 1 위 = ${t&&t.phaseName} [${t&&t.zone}] score ${t&&t.score.toFixed(2)}`);
  ok(t&&Math.abs(t.kappa-1)<0.004,`[${z}] κ = ${t&&t.kappa.toFixed(5)} (정답 1)`);
  const F=DI.rankSpots(A.L,phases,{distortion:0.01,maxIdx:5,gmax:15,w:o.w,h:o.h,mask:A.bs.mask}); const iSi=F.cands.findIndex(c=>c.phaseName==="Si_Fd-3m"&&sameZone(c.zone,z));
  // 보정 없이는 모양(길이비·각)이 같은 후보를 가를 수 없다 — Si 보다 위의 후보는 모두 같은 모양이어야 한다
  const ahead=F.cands.slice(0,Math.max(0,iSi)); const sameShape=ahead.every(c=>Math.abs(c.ratio21-F.cands[iSi].ratio21)<0.01&&Math.abs(c.ang12-F.cands[iSi].ang12)<0.5);
  ok(iSi>=0&&iSi<4&&sameShape,`[${z}] 보정 없음: Si 순위 ${iSi+1}, 앞선 후보 ${ahead.map(c=>c.phaseName+"["+c.zone+"]").join(", ")||"없음"} (모두 같은 모양 ${sameShape})`);
  ok(F.cands[iSi]&&Math.abs(F.cands[iSi].scaleImplied/o.scale-1)<0.004,`[${z}] 보정 없음: 함께 알려 주는 축척 ${F.cands[iSi]&&F.cands[iSi].scaleImplied.toFixed(6)} (정답 ${o.scale})`);
}

console.log("== 아핀 보정: Si [001] 에 찌그러짐 5 % (φ 30°) — 표준 시료로 되찾기");
{ const e=1.05, phi=30*Math.PI/180; const Tq=DI.affineFromParams(1,e,phi); const T=[[Tq.T[0][0],Tq.T[0][1]],[Tq.T[1][0],Tq.T[1][1]]];
  const o=Object.assign({},base,{seed:SEED+20,zone:[0,0,1],rot:0.2,cx:256.4,cy:249.2,scale:0.035,T}); const A=analyze(o);
  // 기대: 관측 px = T · (g/s) → 보정 px→1/nm 행렬 = s · T⁻¹, 즉 찌그러짐 e, 장축(큰 축척) 방향 = T 의 단축 방향 (φ+90°)
  const zb=DC.zoneBasis(PH.Si,[0,0,1]); const cal=DI.calibrateAffine(A.L.R1,A.L.R2,zb.g[0],zb.g[1],zb.ang12);
  let dphi=Math.abs(((cal.phi-(phi+Math.PI/2))%Math.PI+Math.PI)%Math.PI); dphi=Math.min(dphi,Math.PI-dphi);
  ok(Math.abs(cal.ellip-e)<0.004,`찌그러짐 ${cal.ellip.toFixed(4)} (정답 ${e})`); ok(dphi<1.5*Math.PI/180,`방향 오차 ${(dphi*180/Math.PI).toFixed(2)}°`);
  ok(Math.abs(cal.scale/o.scale-1)<0.003,`평균 축척 ${cal.scale.toFixed(6)} (정답 ${o.scale})`);
  const C=DI.correctLattice(A.L,cal); ok(Math.abs(C.g[0]/zb.g[0]-1)<0.002&&Math.abs(C.g[1]/zb.g[1]-1)<0.002&&Math.abs(C.ang12-zb.ang12)<0.2,`보정 뒤 g ${C.g.map(x=>x.toFixed(4)).join("/")} ∠ ${C.ang12.toFixed(2)}° (정답 ${zb.g.map(x=>x.toFixed(4)).join("/")} ∠ ${zb.ang12.toFixed(2)}°)`);
  // 저장한 보정(축척·찌그러짐·방향)으로 다시 만든 행렬이 같은 결과를 내는지 (회전 없는 대칭 행렬)
  const cal2=DI.affineFromParams(cal.scale,cal.ellip,cal.phi); const C2=DI.correctLattice(A.L,cal2); ok(Math.abs(C2.g[0]-C.g[0])/C.g[0]<1e-3&&Math.abs(C2.ang12-C.ang12)<0.05,`저장값으로 다시 보정: g ${C2.g.map(x=>x.toFixed(4)).join("/")} ∠ ${C2.ang12.toFixed(2)}°`);
  // 찌그러짐을 모르면: 보정 없는 모양 비교가 거의 정사각 [001] 을 놓칠 수 있다 → 허용 찌그러짐을 넓히면 찾는다
  const R=DI.rankSpots(A.L,phases,{distortion:0.05,maxIdx:4,gmax:15}); ok(R.cands.slice(0,6).some(c=>c.phaseName==="Si_Fd-3m"&&sameZone(c.zone,[0,0,1])),"허용 찌그러짐 5 %: Si [001] 이 상위 6 안");
}

console.log("== 방사 왜곡 (가장자리 1 %)");
{ const r0=230, k=0.01/(r0*r0); const o=Object.assign({},base,{seed:SEED+30,zone:[0,0,1],rot:0.5,cx:256.2,cy:255.7,scale:0.06,k,gmax:24,g0:30}); const A=analyze(o);
  const zb=DC.zoneBasis(PH.Si,[0,0,1]); const want=zb.g[0]/o.scale; const L0=DP.fitLattice(A.S.spots,A.c,{size:o.w,radial:false});
  ok(A.L.radialUsed&&Math.abs(A.L.len[0]-want)<0.06&&Math.abs(A.L.len[1]-want)<0.06,`auto: 왜곡 항 사용 ${A.L.radialUsed}, R1·R2 ${A.L.len.slice(0,2).map(x=>x.toFixed(3)).join("/")} px (정답 ${want.toFixed(3)}; 왜곡 무시 ${L0.len.slice(0,2).map(x=>x.toFixed(3)).join("/")})`);
  ok(Math.abs(A.L.radial/k-1)<0.15,`k ${A.L.radial.toExponential(3)} (정답 ${k.toExponential(3)})`);
  ok(A.L.nIndexed>=L0.nIndexed,`지수 붙은 스폿 ${A.L.nIndexed} ≥ 왜곡 무시 ${L0.nIndexed}`);
}

console.log("== 링: Au 자동 중심·타원 + 상 찾기");
{ const R=rng(SEED+9); const w=800, h=800, cx=402.3, cy=391.7, scale=0.02, ell=1.03, phi=0.6; const rings=DC.ringList(PH.Au,14); const img=new Float32Array(w*h); const A0=600/(rings[0].m*rings[0].I/(rings[0].g**2));
  const cs=Math.cos(phi), sn=Math.sin(phi);
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; const u=dx*cs+dy*sn, v=-dx*sn+dy*cs; const r=Math.hypot(u/ell,v);
    let I=8+400*Math.exp(-r/60); for(const g of rings){ const rr=g.g/scale; if(rr>620) continue; I+=A0*g.m*g.I/(g.g*g.g)*Math.exp(-((r-rr)**2)/(2*1.8*1.8)); } img[y*w+x]=poisson(R,I); } // 방위 평균 높이 ∝ m|F|²/g²
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; if(dx*dx+dy*dy<22*22||(dx<0&&Math.abs(dy)<8)) img[y*w+x]=poisson(R,0.3); }
  const bs=DP.findBeamStop(img,w,h); const A=DP.analyzeRings(img,w,h,{mask:bs.mask});
  ok(Math.hypot(A.center.x-cx,A.center.y-cy)<0.3,`자동 중심 오차 ${Math.hypot(A.center.x-cx,A.center.y-cy).toFixed(3)} px (거친 단계 ${Math.hypot(A.coarse.x-cx,A.coarse.y-cy).toFixed(1)} px)`);
  ok(Math.abs(A.ratio-ell)<0.003,`장단축비 ${A.ratio.toFixed(4)} (정답 ${ell})`);
  const M=DI.matchRings(A.peaks,phases,{}); const t=M[0];
  ok(t&&t.phaseName==="Au_Fm-3m",`보정 없음 1 위 ${t&&t.phaseName} (ρ ${t&&t.rho.toFixed(2)}; 2 위 ${M[1]&&M[1].phaseName} ${M[1]&&M[1].score.toFixed(3)} vs ${t&&t.score.toFixed(3)})`);
  ok(t&&Math.abs(t.scale/scale-1)<0.002,`축척 ${t&&t.scale.toFixed(6)} ± ${t&&t.scaleErr.toExponential(1)} (정답 ${scale})`);
  const M2=DI.matchRings(A.peaks,phases,{scale:scale*1.05,scaleErr:0.1}); ok(M2[0]&&M2[0].phaseName==="Au_Fm-3m"&&Math.abs(M2[0].kappa-1/1.05)<0.003,`보정 5 % 틀림: κ ${M2[0]&&M2[0].kappa.toFixed(4)} (정답 ${(1/1.05).toFixed(4)})`);
}
console.log("== 링: 큰 찌그러짐 (같은 장비 Si 에서 본 7.4 % 와 10 %)");
for(const ell of [1.074,1.10]){ const R=rng(SEED+40); const w=800, h=800, cx=402.3, cy=391.7, scale=0.02, phi=0.6; const rings=DC.ringList(PH.Au,14); const img=new Float32Array(w*h); const A0=600/(rings[0].m*rings[0].I/(rings[0].g**2)); const cs=Math.cos(phi), sn=Math.sin(phi);
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; const u=dx*cs+dy*sn, v=-dx*sn+dy*cs; const r=Math.hypot(u/ell,v); let I=8+400*Math.exp(-r/60); for(const g of rings){ const rr=g.g/scale; if(rr>620) continue; I+=A0*g.m*g.I/(g.g*g.g)*Math.exp(-((r-rr)**2)/(2*1.8*1.8)); } img[y*w+x]=poisson(R,I); }
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-cx, dy=y-cy; if(dx*dx+dy*dy<22*22||(dx<0&&Math.abs(dy)<8)) img[y*w+x]=poisson(R,0.3); }
  const bs=DP.findBeamStop(img,w,h); const A=DP.analyzeRings(img,w,h,{mask:bs.mask}); const M=DI.matchRings(A.peaks,phases,{});
  ok(Math.abs(A.ratio-ell)<0.003&&Math.hypot(A.center.x-cx,A.center.y-cy)<0.5,`e ${ell}: 장단축비 ${A.ratio.toFixed(4)}, 중심 오차 ${Math.hypot(A.center.x-cx,A.center.y-cy).toFixed(2)} px`);
  ok(M[0]&&M[0].phaseName==="Au_Fm-3m"&&Math.abs(M[0].scale/scale-1)<0.003,`e ${ell}: 1 위 ${M[0]&&M[0].phaseName}, 단축 기준 축척 ${M[0]&&M[0].scale.toFixed(6)} (정답 ${scale})`); }
console.log("== 동역학 세기(Bloch, 두께 10–80 nm 평균)로 만든 Si [011]: 이중 회절 (200) 자리가 보임");
{ let DY=null; try{ DY=require("../src/diff_dyn.js"); }catch(e){}
  if(!DY){ console.log("  diff_dyn.js 없음 — 건너뜀"); } else {
    const o=Object.assign({},base,{seed:SEED+50,zone:[0,1,1],rot:0.9,cx:252.3,cy:259.6,scale:0.04,dyn:DY,g0:40,peak:1500}); const A=analyze(o);
    const has200=A.truth.some(t=>Math.abs(t.h[0])===2&&t.h[1]===0&&t.h[2]===0||(t.h[0]===0&&Math.abs(t.h[1])===2&&t.h[2]===0)||(t.h[0]===0&&t.h[1]===0&&Math.abs(t.h[2])===2));
    const K=DI.rankSpots(A.L,phases,{scale:o.scale,scaleErr:0.03,distortion:0.01,maxIdx:5,gmax:15,w:o.w,h:o.h,mask:A.bs.mask,dynamic:false}); const Dn=DI.rankSpots(A.L,phases,{scale:o.scale,scaleErr:0.03,distortion:0.01,maxIdx:5,gmax:15,w:o.w,h:o.h,mask:A.bs.mask,dynamic:true,kV:200});
    const ks=K.cands.find(c=>c.phaseName==="Si_Fd-3m"&&sameZone(c.zone,[0,1,1])), ds=Dn.cands.find(c=>c.phaseName==="Si_Fd-3m"&&sameZone(c.zone,[0,1,1]));
    ok(has200&&ks&&ks.forb>0.005,`운동학 채점은 (200) 이중 회절 스폿을 '금지 자리 세기' ${ks&&(ks.forb*100).toFixed(1)} % 로 벌점`);
    ok(ds&&ds.dyn&&ds.dyn.forb<ks.forb*0.2&&ds.rho>ks.rho,`동역학 채점: 금지 벌점 ${ds&&(ds.dyn.forb*100).toFixed(2)} %, ρ ${ds&&ds.rho.toFixed(3)} (운동학 ${ks&&ks.rho.toFixed(3)})`);
    ok(Dn.cands[0].phaseName==="Si_Fd-3m"&&sameZone(Dn.cands[0].zone,[0,1,1]),`동역학 1 위 ${Dn.cands[0].phaseName} [${Dn.cands[0].zone}] (score ${Dn.cands[0].score.toFixed(3)}; 운동학 점수 ${ks.score.toFixed(3)} → ${ds.score.toFixed(3)})`); } }
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
