// 실측 회절 파일로 반증 시험 (ref/gen_real.py 가 만든 float32; 없으면 건너뜀). 정답을 아는 것만 단언한다:
//  Au 나노입자 링(Libra 200, 375 mm) · Si [001] 스폿(같은 장비, 300 mm, 파일 교정 단위 nm = 잘못 저장) · Fe-Cr ODS SAD (8-bit RGB, 스케일바 5 1/nm)
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js"); const DP=require("../src/diff_pattern.js"); const DI=require("../src/diff_index.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };
const R=path.join(__dirname,"..","ref","real"); const M=fs.existsSync(path.join(R,"real.json"))?JSON.parse(fs.readFileSync(path.join(R,"real.json"),"utf8")):null;
if(!M){ console.log("실측 자료 없음 (ref/gen_real.py) — 건너뜀"); process.exit(0); }
const load=f=>{ const b=fs.readFileSync(path.join(R,f+".f32")); return new Float32Array(b.buffer,b.byteOffset,b.length/4); };
const phases=fs.readdirSync(path.join(__dirname,"cif")).map(f=>{ const p=DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif",f),"utf8"))[0]; p.name=f.replace(".cif",""); p.id="t:"+p.name; return p; });
const fe=DC.phaseFromRecord({id:"t:Fe",name:"Fe_bcc",formula:"Fe",cell:{a:0.28665,b:0.28665,c:0.28665},sg:{number:229},atoms:[{label:"Fe",el:"Fe",x:0,y:0,z:0,occ:1,B:0.35}]}); fe.name="Fe_bcc";
let auScale=null;
const AU="GOLD-NP-DIFF.dm3"; if(M[AU]){ console.log("== Au 링 (DM 교정 "+M[AU].scale+" 1/nm/px)"); const m=M[AU]; const img=load(AU); const bs=DP.findBeamStop(img,m.w,m.h); const A=DP.analyzeRings(img,m.w,m.h,{mask:bs.mask});
  ok(Math.abs(A.ratio-1)<0.01,`타원 찌그러짐 ${((A.ratio-1)*100).toFixed(2)} % (< 1 %)`);
  const T=DI.matchRings(A.peaks,phases.concat([fe]),{}); auScale=T[0].scale;
  ok(T[0].phaseName==="Au_Fm-3m",`보정 없이 1 위 ${T[0].phaseName} (2 위 ${T[1].phaseName} ${T[1].score.toFixed(3)} < ${T[0].score.toFixed(3)})`);
  ok(T[0].pairs.length>=10&&T[0].dev<0.005,`링 ${T[0].pairs.length} 개 짝지음, 평균 어긋남 ${(T[0].dev*100).toFixed(2)} %`);
  ok(Math.abs(T[0].scale/m.scale-1.089)<0.01,`실제 축척 ${T[0].scale.toFixed(6)} = DM 교정 × ${(T[0].scale/m.scale).toFixed(4)} (DM 태그가 약 9 % 작음)`); }
const SI="Si001_300mm_alpha0.1_beta9.4_diff.dm3"; if(M[SI]){ console.log("== Si [001] 스폿 (파일 교정 단위 "+M[SI].unit+")"); const m=M[SI]; const img=load(SI); const bs=DP.findBeamStop(img,m.w,m.h); const S=DP.findSpots(img,m.w,m.h,{mask:bs.mask}); const c=DP.findCenter(S.spots,m.w,m.h,{beamStop:bs.mask}); const L=DP.fitLattice(S.spots,c,{size:m.w});
  ok(c.origin==="stop",`중심 = 빔스토퍼 뒤 Friedel 중심 (${c.x.toFixed(1)}, ${c.y.toFixed(1)})`);
  ok(Math.abs(L.len[0]-336.9)<1.5&&Math.abs(L.len[1]-361.4)<1.5&&Math.abs(L.ang12-90.35)<0.5,`R1·R2 ${L.len.slice(0,2).map(x=>x.toFixed(2)).join("/")} px ∠ ${L.ang12.toFixed(2)}° (직교, 길이 7 % 다름 = 투사 렌즈 찌그러짐)`);
  const Rk=DI.rankSpots(L,phases,{distortion:0.08,maxIdx:5,gmax:15,w:m.w,h:m.h,mask:bs.mask}); const i=Rk.cands.findIndex(x=>x.phaseName==="Si_Fd-3m"&&x.zone.map(Math.abs).sort().join()==="0,0,1");
  ok(i>=0&&i<8,`찌그러짐 허용 8 %: Si [001] 순위 ${i+1}`);
  const zb=DC.zoneBasis(phases.find(p=>p.name==="Si_Fd-3m"),[0,0,1]); const cal=DI.calibrateAffine(L.R1,L.R2,zb.g[0],zb.g[1],zb.ang12);
  ok(Math.abs(cal.ellip-1.073)<0.01,`Si 를 표준으로 한 아핀 보정: 찌그러짐 ${((cal.ellip-1)*100).toFixed(2)} % · 평균 축척 ${cal.scale.toFixed(6)}`);
  if(auScale){ const pred=auScale*375/300; ok(Math.abs(cal.scale/pred-1)<0.03,`같은 장비 교차 확인: Au(375 mm) 축척 × 375/300 = ${pred.toFixed(6)} ↔ Si(300 mm) ${cal.scale.toFixed(6)} (${((cal.scale/pred-1)*100).toFixed(1)} %)`); } }
const FE="FeCr_ODS_SAD_zenodo19217212.tif"; if(M[FE]){ console.log("== Fe-Cr ODS SAD (8-bit, 스케일바)"); const m=M[FE]; const img=load(FE); const ov=DP.findOverlay(img,m.w,m.h); const bs=DP.findBeamStop(img,m.w,m.h); const mask=new Uint8Array(m.w*m.h); for(let i=0;i<mask.length;i++) mask[i]=bs.mask[i]|ov.mask[i];
  const bar=DP.scaleBar(ov,m.w,m.h); ok(bar&&Math.abs(bar.lenPx-177)<=1,`스케일바 ${bar&&bar.lenPx} px`); const scale=5/bar.lenPx;
  const S=DP.findSpots(img,m.w,m.h,{mask}); const c=DP.findCenter(S.spots,m.w,m.h,{beamStop:bs.mask}); ok(c.origin==="spot"&&Math.hypot(c.x-466.5,c.y-531.8)<3,`중심 = 포화된 직접빔 (${c.x.toFixed(1)}, ${c.y.toFixed(1)}) — 반 격자점(${c.alternative?c.alternative.x.toFixed(0)+","+c.alternative.y.toFixed(0):"-"}) 이 아니라`);
  const L=DP.fitLattice(S.spots.slice(0,16),c,{size:m.w}); const Rk=DI.rankSpots(L,[fe].concat(phases),{scale,scaleErr:0.05,distortion:0.02,maxIdx:7,gmax:15,w:m.w,h:m.h,mask}); const t=Rk.cands[0];
  ok(t.phaseName==="Fe_bcc"&&t.zone.map(Math.abs).sort().join()==="1,3,5",`센 16 스폿: 1 위 ${t.phaseName} [${t.zone}] (기지 bcc, [135] 계열)`);
  ok(Math.abs(t.kappa-1)<0.03,`κ = ${t.kappa.toFixed(4)} — 그림 스케일바와 bcc Fe 격자상수가 ${((t.kappa-1)*100).toFixed(1)} % 안에서 맞음`); }
// 자동 판정 (Worker 처리기 그대로): 링 = 광선마다 뚜렷한 봉우리, 스폿 = 뚜렷한 격자
{ require("../src/diff_work.js"); const H=globalThis.DIFF_WORK.handler(globalThis); console.log("== 자동 판정 (스폿/링)");
  for(const [f,want] of [[AU,"rings"],[SI,"spots"],[FE,"spots"]]){ if(!M[f]) continue; const m=M[f]; const r=H("analyze",{img:Float32Array.from(load(f)),w:m.w,h:m.h,opt:{pattern:"auto"}}); ok(r.type===want,`${f}: ${r.type} (${JSON.stringify(r.typeWhy)})`); } }
// 큰 영상: 줄여서 찾고 원래 해상도에서 다듬은 결과 = 원래 해상도 전체 계산 (그리고 4096² 로 키운 영상도)
{ const H=globalThis.DIFF_WORK.handler(globalThis); console.log("== 큰 영상: 줄여 찾기 + 원래 해상도 다듬기");
  const up2=(src,w0)=>{ const w=w0*2, o=new Float32Array(w*w); for(let y=0;y<w;y++) for(let x=0;x<w;x++){ const X=Math.max(0,Math.min(w0-1.001,(x+0.5)/2-0.5)), Y=Math.max(0,Math.min(w0-1.001,(y+0.5)/2-0.5)); const x0=X|0,y0=Y|0,fx=X-x0,fy=Y-y0; o[y*w+x]=(1-fx)*(1-fy)*src[y0*w0+x0]+fx*(1-fy)*src[y0*w0+x0+1]+(1-fx)*fy*src[(y0+1)*w0+x0]+fx*fy*src[(y0+1)*w0+x0+1]; } return o; };
  if(M[SI]){ const m=M[SI]; const A=H("analyze",{img:Float32Array.from(load(SI)),w:m.w,h:m.h,opt:{pattern:"spots",full:true}}), B=H("analyze",{img:Float32Array.from(load(SI)),w:m.w,h:m.h,opt:{pattern:"spots"}});
    const t=Date.now(); const C=H("analyze",{img:up2(load(SI),m.w),w:2*m.w,h:2*m.h,opt:{pattern:"auto"}}); const tc=Date.now()-t; const La=A.spots.lattices[0], Lb=B.spots.lattices[0], Lc=C.spots.lattices[0];
    const rel=(x,y)=>Math.abs(x-y)/y; ok(B.work===2&&rel(Lb.len[0],La.len[0])<1e-3&&rel(Lb.len[1],La.len[1])<1e-3&&Math.abs(Lb.ang12-La.ang12)<0.05,`Si 2048²: 줄여 찾기(×${B.work}) R1·R2 ${Lb.len.slice(0,2).map(x=>x.toFixed(2)).join("/")} = 전체 ${La.len.slice(0,2).map(x=>x.toFixed(2)).join("/")} px (${A.ms} → ${B.ms} ms)`);
    ok(C.work===4&&C.type==="spots"&&rel(Lc.len[0]/2,La.len[0])<1.5e-3&&rel(Lc.len[1]/2,La.len[1])<1.5e-3,`Si 4096²(키운 영상): R1·R2 ÷ 2 = ${(Lc.len[0]/2).toFixed(2)}/${(Lc.len[1]/2).toFixed(2)} px, ${tc} ms`); }
  if(M[AU]){ const m=M[AU]; const A=H("analyze",{img:Float32Array.from(load(AU)),w:m.w,h:m.h,opt:{pattern:"rings",full:true}}), B=H("analyze",{img:Float32Array.from(load(AU)),w:m.w,h:m.h,opt:{pattern:"rings"}});
    const ph=phases.find(p=>p.name==="Au_Fm-3m"); const sa=DI.matchRings(A.rings.peaks,[ph],{})[0].scale, sb=DI.matchRings(B.rings.peaks,[ph],{})[0].scale;
    ok(Math.abs(sb/sa-1)<1e-3&&Math.hypot(B.rings.center.x-A.rings.center.x,B.rings.center.y-A.rings.center.y)<0.5&&Math.abs(B.rings.ratio-A.rings.ratio)<0.001,`Au 2048²: 줄여 찾기 축척 ${sb.toFixed(6)} = 전체 ${sa.toFixed(6)} (${((sb/sa-1)*100).toFixed(3)} %), 중심 차 ${Math.hypot(B.rings.center.x-A.rings.center.x,B.rings.center.y-A.rings.center.y).toFixed(2)} px (${A.ms} → ${B.ms} ms)`);
    const t=Date.now(); const C=H("analyze",{img:up2(load(AU),m.w),w:2*m.w,h:2*m.h,opt:{pattern:"auto"}}); const tc=Date.now()-t; const sc=DI.matchRings(C.rings.peaks,[ph],{})[0].scale*2;
    ok(C.type==="rings"&&Math.abs(sc/sa-1)<2e-3,`Au 4096²(키운 영상): 링, 축척×2 ${sc.toFixed(6)} (${((sc/sa-1)*100).toFixed(3)} %), ${tc} ms`); } }
// 내장 묶음(COD) 으로: 원소 필터(REDI 방식 — 상의 원소 ⊆ EDS 원소)
let PACK=null; try{ require("../src/diff_pack.js"); PACK=globalThis.DIFF_PACK; }catch(e){}
if(PACK&&M[FE]&&M[AU]){ console.log("== 내장 묶음 (COD "+PACK.phases.length+" 상) + 원소 필터");
  const only=els=>PACK.phases.filter(r=>Object.keys(DC.parseFormula(r.formula)).every(e=>els.includes(e))).map(r=>{ const p=DC.phaseFromRecord(r); p.id=r.id; p.name=r.name; return p; });
  { const m=M[AU]; const img=load(AU); const bs=DP.findBeamStop(img,m.w,m.h); const A=DP.analyzeRings(img,m.w,m.h,{mask:bs.mask}); const all=DI.matchRings(A.peaks,only(["Au","Al","Cu","Ni","Pt","Ag","Si","C"]),{}); const T=DI.matchRings(A.peaks,only(["Au","C","O"]),{});
    ok(T[0]&&T[0].phase.id==="pack:Au"&&Math.abs(T[0].scale/m.scale-1.089)<0.01,`Au (원소 Au C O): 1 위 ${T[0]&&T[0].phaseName}, 축척 ÷ DM = ${T[0]&&(T[0].scale/m.scale).toFixed(4)}`);
    const fcc=all.slice(0,3).map(c=>c.phase.id); ok(fcc.every(i=>/pack:(Au|Al|Cu|Ni|Pt|Ag)$/.test(i)),`원소 없이 fcc 금속끼리는 못 가름 (상위 3: ${all.slice(0,3).map(c=>c.phaseName).join(", ")})`); }
  { const m=M[FE]; const img=load(FE); const ov=DP.findOverlay(img,m.w,m.h); const bs=DP.findBeamStop(img,m.w,m.h); const mask=new Uint8Array(m.w*m.h); for(let i=0;i<mask.length;i++) mask[i]=bs.mask[i]|ov.mask[i]; const bar=DP.scaleBar(ov,m.w,m.h);
    const S=DP.findSpots(img,m.w,m.h,{mask}); const c=DP.findCenter(S.spots,m.w,m.h,{beamStop:bs.mask}); const L=DP.fitLattice(S.spots.slice(0,16),c,{size:m.w});
    const Rk=DI.rankSpots(L,only(["Fe","Cr","Y","Ti","O"]),{scale:5/bar.lenPx,scaleErr:0.05,distortion:0.02,maxIdx:7,gmax:15,w:m.w,h:m.h,mask}); const t=Rk.cands[0];
    ok(t&&t.phase.id==="pack:aFe"&&t.zone.map(Math.abs).sort().join()==="1,3,5",`Fe-Cr (원소 Fe Cr Y Ti O): 1 위 ${t&&t.phaseName} [${t&&t.zone}] κ ${t&&t.kappa.toFixed(4)}`); } }
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
