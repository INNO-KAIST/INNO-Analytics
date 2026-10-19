// 감사(2026-10-19)에서 찾은 결함의 회귀 시험 (node test/test_fix.js)
//  1) 링·얼룩 패턴을 스폿으로 순위 매겨도 메모리 폭발 없음 (verifyZone gmax 상한·2D 정대축 격자·얼룩 거부·시간 상한)
//  2) 링 판정: 중심이 치우쳐 가장자리에 잘린 링 · 약한 링 · 8-bit 로그 표시 사본 → rings, 스폿 패턴은 spots 그대로
//  3) 이름표 OCR: 단위 확신이 낮으면 보정값으로 바로 쓰지 않음 (DIFF_OCR.labelUse)
//  4) 공칭 보정 ±15 % 창 (scalePrior) · 회절도다움(saedLikeness: 8-bit 로그 SAED ↔ HRTEM)
// 시험 영상은 모두 여기서 합성한다 (실측·내부 표본 없음)
const fs=require("fs"), path=require("path");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js"); require("../src/diff_dyn.js"); const OCR=require("../src/diff_ocr.js");
const DP=require("../src/diff_pattern.js"); const DI=require("../src/diff_index.js"); const DW=require("../src/diff_work.js"); const H=DW.handler(globalThis);
const {rng,poisson,PH,synthSpots}=require("./synth.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };
const phases=fs.readdirSync(path.join(__dirname,"cif")).map(f=>{ const p=DC.parseCIF(fs.readFileSync(path.join(__dirname,"cif",f),"utf8"))[0]; p.name=f.replace(".cif",""); p.id="t:"+p.name; return p; });
const recs=phases.map(p=>Object.assign(DC.recordFromPhase(p,{id:p.id,name:p.name}),{layer:"pack"}));

/* Au 다결정 링 (Poisson): 중심·타원·세기·배경 */
function rings(o){ const R=rng(o.seed); const {w,h}=o; const img=new Float32Array(w*h); const L=DC.ringList(PH.Au,o.gmax||14); let Im=0; for(const q of L) Im=Math.max(Im,(q.I||1)*q.m/(q.g*q.g));
  for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const dx=x-o.cx, dy=y-o.cy; const rr=Math.hypot(dx,dy); let I=o.bg*(1+30*Math.exp(-rr/40));
    for(const q of L){ const rp=q.g/o.scale; const A=o.peak*((q.I||1)*q.m/(q.g*q.g))/Im*Math.exp(-q.g*q.g/144); I+=A*Math.exp(-((rr-rp)**2)/(2*o.sig*o.sig)); }
    img[y*w+x]=poisson(R,I); } return img; }
/* 8-bit 로그 표시 사본 (GMS 'export as displayed' 비슷하게): log1p → 0.1–99.9 백분위 → 0..255 정수 */
function eightBitLog(img){ const L=Float32Array.from(img,v=>Math.log1p(Math.max(0,v))); const s=Float32Array.from(L).sort(); const lo=s[Math.floor(0.001*s.length)], hi=s[Math.floor(0.999*s.length)];
  return Float32Array.from(L,v=>Math.max(0,Math.min(255,Math.round((v-lo)/(hi-lo)*255)))); }
const base={w:768,h:768,cx:380,cy:390,scale:0.016,sig:2,peak:300,bg:10,seed:31};
const heapMB=()=>process.memoryUsage().heapUsed/1048576;

console.log("== 1. 링·얼룩을 스폿으로 순위 매기기 — 메모리·시간 상한");
{ const o=Object.assign({},base,{cx:170,cy:220,seed:32}); const img=rings(o);
  const A=H("analyze",{img:Float32Array.from(img),w:o.w,h:o.h,opt:{pattern:"spots",maxLattices:1}}); const L=A.spots.lattices[0];
  ok(!!L,`치우친 링에서도 (가짜) 격자가 맞춰짐 — ${L?L.len.slice(0,2).map(x=>x.toFixed(1)).join("/")+" px, 지수 "+L.nIndexed:"없음"}`);
  if(L){ const h0=heapMB(), t0=Date.now(); const R=DI.rankSpots(L,phases,{maxIdx:5,gmax:15,w:o.w,h:o.h,nVerify:20,dynamic:false}); const dt=Date.now()-t0;
    ok(R.refused==="speckle"&&R.cands.length===0&&dt<5000,`얼룩 격자는 순위를 매기지 않음 (refused=${R.refused}, ${dt} ms)`);
    ok(DI.isSpeckle(L,o.w),"DI.isSpeckle → true");
    // 거부를 우회해도(검증 직접) 반사 목록이 상한 안: gmax ≤ 25 1/nm, 시간·힙 증가 작음
    const C=DI.matchSpots(DI.measured(L),[phases.find(p=>p.name==="Au_Fm-3m")],{maxIdx:5,gmax:15,top:5}); let gm=0, t1=Date.now();
    for(const c of C){ const V=DI.verifyZone(L,c,{gmax:15}); gm=Math.max(gm,V.gmax); } const dt1=Date.now()-t1, dh=heapMB()-h0;
    ok(C.length>0&&gm<=25&&dt1<5000&&dh<300,`verifyZone 직접: 예측 범위 ≤ ${gm.toFixed(1)} 1/nm (상한 25), ${dt1} ms, 힙 +${dh.toFixed(0)} MB`); }
  // 순수 얼룩 (결정립 반점이 무작위로 흩어진 영상) — 작은 격자라도 거부, 시간 상한
  { const R=rng(77), w=512, h=512; const im=new Float32Array(w*h); for(let i=0;i<im.length;i++) im[i]=poisson(R,8);
    for(let k=0;k<900;k++){ const x=5+R()*(w-10), y=5+R()*(h-10), a=50+400*R(); for(let yy=Math.floor(y)-3;yy<=Math.ceil(y)+3;yy++) for(let xx=Math.floor(x)-3;xx<=Math.ceil(x)+3;xx++) im[yy*w+xx]+=a*Math.exp(-((xx-x)**2+(yy-y)**2)/2.5); }
    const A2=H("analyze",{img:im,w,h,opt:{pattern:"spots",maxLattices:1}}); const L2=A2.spots.lattices[0]; const t0=Date.now();
    const R2=L2?DI.rankSpots(L2,phases,{maxIdx:5,gmax:15,w,h,nVerify:80,budgetMs:8000}):{cands:[]}; const dt=Date.now()-t0;
    ok(dt<10000&&(!L2||R2.refused==="speckle"||R2.cands.length<=80),`무작위 반점: ${L2?"격자 "+L2.len.slice(0,2).map(x=>x.toFixed(1)).join("/")+" px → "+(R2.refused||R2.cands.length+" 후보"):"격자 없음"}, ${dt} ms`); }
  // 시간 상한: budgetMs 를 아주 작게 → truncated, 예외 없음
  { const o2={w:512,h:512,phase:PH.Si,zone:[0,1,1],gmax:12,cx:250,cy:262,rot:0.3,scale:0.02,sig:1.6,peak:400,g0:9,bg:4,seed:7}; const S=synthSpots(o2);
    const A3=H("analyze",{img:S.img,w:512,h:512,opt:{pattern:"spots"}}); const L3=A3.spots.lattices[0]; const R3=DI.rankSpots(L3,phases,{maxIdx:5,gmax:15,w:512,h:512,budgetMs:0,dynamic:false});
    ok(R3.truncated===true&&Array.isArray(R3.cands),`budgetMs 0 → truncated (${R3.cands.length} 후보)`);
    const R4=DI.rankSpots(L3,phases,{maxIdx:5,gmax:15,w:512,h:512,dynamic:false}); ok(!R4.refused&&R4.cands[0]&&R4.cands[0].phaseName==="Si_Fd-3m",`진짜 Si [011] 격자는 그대로 순위 (1 위 ${R4.cands[0]&&R4.cands[0].phaseName}, ${R4.ms} ms)`); } }

console.log("== 2. 링 / 스폿 자동 판정");
{ const cases=[["가운데 링",Object.assign({},base),"rings"],["치우친 링 (가장자리에 잘림)",Object.assign({},base,{cx:170,cy:220,seed:33}),"rings"],["약한 링 (봉우리 20, 배경 10)",Object.assign({},base,{peak:20,seed:34}),"rings"],
    ["8-bit 로그 표시 링",Object.assign({},base,{seed:35,eight:true}),"rings"],["8-bit 로그 + 치우침",Object.assign({},base,{cx:250,cy:300,seed:36,eight:true}),"rings"]];
  for(const [name,o,want] of cases){ let img=rings(o); if(o.eight) img=eightBitLog(img); const r=H("analyze",{img,w:o.w,h:o.h,opt:{pattern:"auto"}}); ok(r.type===want,`${name} → ${r.type} ${JSON.stringify(r.typeWhy)}`); }
  for(const [name,o] of [["Si [011] 스폿",{}],["Si [011] 빔스토퍼",{stop:true,stopR:14,stopAng:0.8,direct:false}],["Si [011] 잡음",{peak:25,bg:6,seed:9}]]){
    const S=synthSpots(Object.assign({w:512,h:512,phase:PH.Si,zone:[0,1,1],gmax:12,cx:250,cy:262,rot:0.3,scale:0.02,sig:1.6,peak:400,g0:9,bg:4,seed:7},o)); const r=H("analyze",{img:S.img,w:512,h:512,opt:{pattern:"auto"}});
    ok(r.type==="spots",`${name} → ${r.type}`); }
  // 큰 영상 (줄여서 찾고 원래 해상도에서 다시 판정)
  { const o=Object.assign({},base,{w:1400,h:1400,cx:500,cy:560,scale:0.009,seed:37}); const r=H("analyze",{img:rings(o),w:o.w,h:o.h,opt:{pattern:"auto"}}); ok(r.type==="rings"&&r.work>1,`1400² 치우친 링 (work ${r.work}) → ${r.type}`); }
  // 광선 수는 보이는 원호 기준: DP.raysInside
  { const el={cx:50,cy:50,a:1,b:1,phi:0}; const n1=DP.raysInside(400,400,el,30), n2=DP.raysInside(400,400,el,200);
    ok(n1===360&&n2>60&&n2<140,`raysInside: 안쪽 원 ${n1}/360, 모서리에 잘린 원 ${n2}/360 (≈ 90°)`); } }

console.log("== 3. 이름표 단위 확신 → 보정에 바로 쓰지 않음");
{ const U=OCR.labelUse;
  const a=U({text:"5 1/A",value:5,unit:"1/Å",valueInvNm:50,recip:true,conf:0.6,unitConf:0.5}); ok(a&&a.auto===false&&a.opts[0]===50&&a.opts[1]===5,`"5 1/A" conf 0.60 → 자동 사용 안 함, 고를 값 ${a&&a.opts.join(" / ")} 1/nm`);
  const b=U({text:"5 1/nm",value:5,unit:"1/nm",valueInvNm:5,recip:true,conf:0.999,unitConf:0.999}); ok(b&&b.auto===true&&b.value===5,"\"5 1/nm\" conf 0.999 → 자동 사용 5 1/nm");
  const c=U({text:"5 1/nm",value:5,unit:"1/nm",valueInvNm:5,recip:true,conf:0.95,unitConf:0.5}); ok(c&&c.auto===false&&c.opts[1]===50,"conf 0.95 이지만 단위 글자 맞춤이 나쁨(unitConf 0.5) → 확인");
  ok(U({text:"5 nm",value:5,unit:"nm",valueInvNm:null,recip:false,conf:0.99})===null&&U({text:"5 1/nm",value:5,unit:"1/nm",valueInvNm:5,recip:true,conf:0.2})===null,"실공간 단위 · conf < 0.3 → null");
  // 실제 OCR: 합성 이름표 묶음이 있으면 conf ≥ 0.9 결과에 unitConf 가 붙는지
  const tf=path.join(__dirname,"ocr","truth.json"); if(fs.existsSync(tf)){ const r=OCR.readScaleLabel(new Float32Array(64*64),64,64,{x0:2,x1:40,y0:50,y1:52}); ok(r===null||typeof r.unitConf==="number","readScaleLabel 결과에 unitConf"); } }

console.log("== 4. 공칭 보정 ±15 % 창 · 회절도다움");
{ const o={w:512,h:512,phase:PH.Si,zone:[0,1,1],gmax:12,cx:250,cy:262,rot:0.3,scale:0.02,sig:1.6,peak:400,g0:9,bg:4,seed:7}; const S=synthSpots(o);
  const A=H("analyze",{img:S.img,w:512,h:512,opt:{pattern:"spots"}}); const L=A.spots.lattices[0];
  const R=DI.rankSpots(L,phases,{maxIdx:5,gmax:15,w:512,h:512,dynamic:false,scalePrior:{scale:0.02,tol:0.15}}); const ks=R.cands.map(c=>c.scaleImplied/0.02);
  ok(R.cands.length&&ks.every(k=>Math.abs(k-1)<=0.15)&&R.cands[0].phaseName==="Si_Fd-3m",`창 안 후보만 (${R.cands.length} 개, κ ${Math.min(...ks).toFixed(3)}–${Math.max(...ks).toFixed(3)}), 1 위 ${R.cands[0]&&R.cands[0].phaseName}`);
  const R2=DI.rankSpots(L,phases,{maxIdx:5,gmax:15,w:512,h:512,dynamic:false,scalePrior:{scale:0.2,tol:0.15}}); ok(R2.priorEmpty===true&&R2.cands.length>0,"창 안에 아무것도 없으면 창 없이 + priorEmpty");
  const img=rings(Object.assign({},base)); const pk=DP.analyzeRings(img,base.w,base.h,{}).peaks; const M=DI.matchRings(pk,phases,{scalePrior:{scale:base.scale,tol:0.15}});
  ok(M.length&&M.every(m=>Math.abs(m.scale/base.scale-1)<=0.15)&&M[0].phaseName==="Au_Fm-3m",`링: 창 안 후보만, 1 위 ${M[0]&&M[0].phaseName}`);
  const M2=DI.matchRings(pk,phases,{scalePrior:{scale:base.scale*10,tol:0.15}}); ok(M2.priorEmpty===true&&M2.length>0,"링: 10 배 틀린 보정 → 창 없이 + priorEmpty");
  // 회절도다움: 8-bit 로그 링·스폿 ↔ HRTEM 비슷한 격자 무늬·잡음
  const lk=(im,w,h)=>DP.saedLikeness(im,w,h); const q1=lk(eightBitLog(img),base.w,base.h), q2=lk(eightBitLog(S.img),512,512);
  const R0=rng(5); const hr=new Float32Array(512*512); for(let y=0;y<512;y++) for(let x=0;x<512;x++) hr[y*512+x]=Math.round(128+40*Math.sin(x*0.6)*Math.sin(y*0.5)+10*(R0()+R0()+R0()-1.5));
  const q3=lk(hr,512,512), q4=lk(new Float32Array(256*256),256,256);
  ok(q1.r2>=0.1&&q1.decay&&q2.r2>=0.1&&q2.decay,`8-bit 로그 SAED: 링 R² ${q1.r2.toFixed(2)}, 스폿 R² ${q2.r2.toFixed(2)} (≥ 0.1, 안쪽이 밝음)`);
  ok(q3.r2<0.04&&q4.r2<0.04,`HRTEM 비슷한 무늬 R² ${q3.r2.toFixed(3)}, 빈 영상 ${q4.r2.toFixed(3)} (< 0.04)`); }

console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ fails.forEach(f=>console.log("  ✗ "+f)); process.exit(1); }
