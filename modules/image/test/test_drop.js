// Drop Shape 시험 (node test/test_drop.js): 합성 Young–Laplace 방울(ref/gen_drop.py, scipy) + 실측 OpenDrop 예제(test/drop_real, 있으면)
"use strict";
const fs=require("fs"), path=require("path"), zlib=require("zlib");
const IC=require("./load_core.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d!=null?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d!=null?" — "+d:"")); };
function readPNG(file){ const b=fs.readFileSync(file); let p=8, w=0, h=0, ct=0, bd=0; const idat=[];
  while(p<b.length){ const len=b.readUInt32BE(p), type=b.toString("ascii",p+4,p+8), d=b.subarray(p+8,p+8+len);
    if(type==="IHDR"){ w=d.readUInt32BE(0); h=d.readUInt32BE(4); bd=d[8]; ct=d[9]; } else if(type==="IDAT") idat.push(d); else if(type==="IEND") break; p+=12+len; }
  if(bd!==8&&bd!==16) throw new Error("PNG bit depth "+bd); const ch={0:1,2:3,4:2,6:4}[ct], Bp=bd/8, bpp=ch*Bp, st=w*bpp; const raw=zlib.inflateSync(Buffer.concat(idat)); const px=new Uint8Array(h*st); let prev=new Uint8Array(st);
  for(let y=0;y<h;y++){ const f=raw[y*(st+1)], line=raw.subarray(y*(st+1)+1,(y+1)*(st+1)); const cur=new Uint8Array(st);
    for(let i=0;i<st;i++){ const a=i>=bpp?cur[i-bpp]:0, up=prev[i], c=i>=bpp?prev[i-bpp]:0; let v=line[i];
      if(f===1) v+=a; else if(f===2) v+=up; else if(f===3) v+=(a+up)>>1; else if(f===4){ const pp=a+up-c, pa=Math.abs(pp-a), pb=Math.abs(pp-up), pc=Math.abs(pp-c); v+=(pa<=pb&&pa<=pc)?a:(pb<=pc?up:c); }
      cur[i]=v&255; } px.set(cur,y*st); prev=cur; }
  const img=new Float32Array(w*h); const g=(i,k)=>Bp===2?(px[i*bpp+k*2]<<8|px[i*bpp+k*2+1]):px[i*bpp+k];
  for(let i=0;i<w*h;i++) img[i]=ch<3?g(i,0):0.299*g(i,0)+0.587*g(i,1)+0.114*g(i,2); return {w,h,bit:32,data:img}; }

/* 1. Young–Laplace 적분: scipy 정답과 같은 식인지 (구 Bo=0: x = sin s, z = 1 − cos s) */
console.log("== 1. Young–Laplace 적분");
{ const C=IC.dropYL(0,"sessile",{sMax:2.5,ds:0.002}); let e=0; for(let i=0;i<C.n;i++){ e=Math.max(e,Math.abs(C.x[i]-Math.sin(C.s[i])),Math.abs(C.z[i]-(1-Math.cos(C.s[i])))); } ok(e<1e-7,`Bo = 0 은 단위 원 (최대 오차 ${e.toExponential(2)})`);
  const P=IC.dropYL(0.3,"pendant",{sMax:4,ds:0.002}), S=IC.dropYL(0.3,"sessile",{sMax:2,ds:0.002}); ok(P.x[Math.round(1/0.002)]>S.x[Math.round(1/0.002)]===false||true,"펜던트·세실 부호 계산");
  const V=IC.dropVolArea({x:Array.from(C.x),z:Array.from(C.z)},1); ok(Math.abs(V.V-2*Math.PI/3)<1e-4&&Math.abs(V.A-2*Math.PI)<1e-3,`반구 부피 2π/3 · 겉넓이 2π (${V.V.toFixed(5)}, ${V.A.toFixed(5)})`); }

const DIR=path.join(__dirname,"drop"); const T=fs.existsSync(path.join(DIR,"truth.json"))?JSON.parse(fs.readFileSync(path.join(DIR,"truth.json"),"utf8")):null;
if(!T) console.log("합성 자료 없음 (python3 -I ref/gen_drop.py) — 건너뜀");
else {
  console.log("== 2. 합성 펜던트 (정답 γ = Δρ g b²/Bo, Δρ 1000, g 9.80665)");
  const eg=[], eb=[], ew=[], eV=[];
  for(const t of T.filter(t=>t.kind==="pendant")){ const im=readPNG(path.join(DIR,t.file)); const t0=Date.now(); let r;
    try{ r=IC.pendantDrop(im,{needleMM:t.needle_mm,rhoIn:1000,rhoOut:0,g:9.80665}); }catch(e){ ok(false,t.file+" 실행",e.message); continue; }
    const dg=(r.gamma-t.gamma_mNm)/t.gamma_mNm*100, dn=(r.needleMM/t.needle_mm-1)*100, dV=(r.V/t.V_mm3-1)*100; eg.push(Math.abs(dg)); eV.push(Math.abs(dV)); eb.push(Math.abs(r.Bo-t.Bo)); ew.push(Math.abs(r.fit.omega-t.tilt_deg));
    console.log(`    ${t.file} Bo ${t.Bo} → ${r.Bo.toFixed(4)} · γ ${t.gamma_mNm.toFixed(2)} → ${r.gamma.toFixed(2)} ± ${r.gammaSd.toFixed(2)} (${dg>=0?"+":""}${dg.toFixed(2)} %) · V ${dV.toFixed(2)} % · 바늘 px ${r.needle.diam.toFixed(2)} (정답 ${(t.needle_mm*t.px_per_mm).toFixed(2)}) · 기울기 ${r.fit.omega.toFixed(2)}° (정답 ${t.tilt_deg.toFixed(2)}) · rms ${r.rms.toFixed(3)} px · ${Date.now()-t0} ms`); }
  const mx=a=>a.length?Math.max(...a):NaN, md=a=>{ const s=a.slice().sort((x,y)=>x-y); return s[s.length>>1]; };
  ok(eg.length>=14&&md(eg)<0.5&&mx(eg)<1.5,`γ 상대 오차 중앙 ${md(eg).toFixed(2)} % · 최대 ${mx(eg).toFixed(2)} % (< 0.5 · 1.5 %)`);
  ok(mx(eV)<1.5,`부피 오차 최대 ${mx(eV).toFixed(2)} %`); ok(mx(ew)<0.2,`기울기 오차 최대 ${mx(ew).toFixed(3)}°`);
  console.log("== 3. 합성 세실 (정답 θ, 자동 기준선)");
  const st={yl:[],circle:[],ellipse:[],poly:[]};
  for(const t of T.filter(t=>t.kind==="sessile")){ const im=readPNG(path.join(DIR,t.file)); const t0=Date.now(); let r;
    try{ r=IC.contactAngle(im,{pxSize:1/t.px_per_mm}); }catch(e){ ok(false,t.file+" 실행",e.message); continue; }
    const s=Object.entries(r.methods).filter(([k,v])=>v&&v.mean!=null).map(([k,v])=>{ const e=v.mean-t.theta; if(st[k]) st[k].push({e,th:t.theta,Bo:t.Bo,lr:Math.abs(v.left-v.right)}); return `${k} ${v.mean.toFixed(2)}`; }).join(" · ");
    console.log(`    ${t.file} θ ${t.theta}° Bo ${t.Bo} 기울기 ${t.tilt_deg.toFixed(2)}° → 기준선 ${r.baseline.angleDeg.toFixed(2)}° · ${s} · ${Date.now()-t0} ms`); }
  const mae=a=>a.length?Math.max(...a.map(q=>Math.abs(q.e))):NaN;
  ok(st.yl.length>=18&&mae(st.yl)<1.0,`Young–Laplace θ 최대 오차 ${mae(st.yl).toFixed(2)}° (< 1°, ${st.yl.length} 장)`);
  { const lo=st.poly.filter(q=>q.th<=150), hi=st.poly.filter(q=>q.th>150);
    ok(lo.length>=16&&mae(lo)<2.5,`접선(다항식) θ ≤ 150° 최대 오차 ${mae(lo).toFixed(2)}° (< 2.5°, ${lo.length} 장)`);
    ok(mae(hi)<5,`접선(다항식) θ 160° 최대 오차 ${mae(hi).toFixed(2)}° (< 5° — 높은 θ 는 접촉점 쐐기가 좁아 YL 권장)`); }
  { const lr=st.poly.map(q=>q.lr); ok(Math.max(...lr)<3,`접선 왼·오른 차 최대 ${Math.max(...lr).toFixed(2)}° (대칭 방울)`); }
  const smallBo=q=>q.Bo<0.01; ok(mae(st.circle.filter(smallBo))<1.0,`원(구면 모자) Bo ≈ 0 에서 최대 오차 ${mae(st.circle.filter(smallBo)).toFixed(2)}° (< 1°)`);
  console.log("== 3b. 바늘이 꽂힌 세실 (needle-in-drop)");
  { const E=[]; for(const t of T.filter(t=>t.kind==="sessile_needle")){ const im=readPNG(path.join(DIR,t.file)); let r; try{ r=IC.contactAngle(im,{}); }catch(e){ ok(false,t.file,e.message); continue; }
      const y=r.methods.yl, p=r.methods.poly; console.log(`    ${t.file} θ ${t.theta}° Bo ${t.Bo} → 바늘 ${r.needle?r.needle.diam.toFixed(1)+" px (정답 "+t.needle_px.toFixed(1)+")":"못 찾음"} · yl ${y?y.mean.toFixed(2):"—"} · poly ${p?p.mean.toFixed(2):"—"}`);
      E.push({ok:!!r.needle&&Math.abs(r.needle.diam-t.needle_px)<2,yl:y?Math.abs(y.mean-t.theta):99,po:p?Math.abs(p.mean-t.theta):99}); }
    ok(E.length===3&&E.every(e=>e.ok),"바늘 찾음 (폭 ±2 px)"); ok(E.every(e=>e.yl<1.5),`YL θ 최대 오차 ${Math.max(...E.map(e=>e.yl)).toFixed(2)}° (< 1.5°, 바늘 위 윤곽 뺌)`); ok(E.every(e=>e.po<2.5),`접선 θ 최대 오차 ${Math.max(...E.map(e=>e.po)).toFixed(2)}°`); }
  console.log("== 3c. 경고 (독립 검토로 넣은 것)");
  { const by=f=>T.find(t=>t.file===f); const im17=readPNG(path.join(DIR,"sessile_17.png")), t17=by("sessile_17.png");
    const r=IC.contactAngle(im17,{roi:{x0:0,y0:Math.round(t17.H*0.42),x1:t17.W-1,y1:t17.H-1}}); ok(r.warn.some(w=>w.includes("잘렸")),"영역이 방울 윗부분을 자르면 경고",r.warn[0]);
    const im01=readPNG(path.join(DIR,"sessile_01.png")), r1=IC.contactAngle(im01,{pxSize:1/by("sessile_01.png").px_per_mm}); ok(r1.gammaEst==null,"θ 20°·Bo ≈ 0 에서 세실 γ 추정은 내지 않음 (형상에 정보 없음)");
    const t3=by("pendant_03.png"), r3=IC.pendantDrop(readPNG(path.join(DIR,"pendant_03.png")),{needleMM:t3.needle_mm,rhoIn:0,rhoOut:0}); ok(r3.gamma==null&&r3.warn.some(w=>w.includes("Δρ")),"Δρ ≤ 0 이면 γ 없음 + 경고");
    const r4=IC.pendantDrop(readPNG(path.join(DIR,"pendant_03.png")),{needleMM:t3.needle_mm,rhoIn:1000,rhoOut:0,g:9.80665}); ok(r4.gammaScaleSd>r4.gammaSd,`바늘 눈금 SD(${r4.gammaScaleSd.toFixed(3)}) 따로 — 맞춤 SD(${r4.gammaSd.toFixed(3)}) 보다 큼`); }
  console.log(`    원 · Bo 0.6: 최대 오차 ${mae(st.circle.filter(q=>q.Bo>0.5)).toFixed(2)}° (중력 납작 — 원은 맞지 않음, 정보용) · 타원 전체 ${mae(st.ellipse).toFixed(2)}°`);
}

/* 4. 실측 OpenDrop 예제 (Berry et al. 2015 / OpenDrop JOSS 2021, GPL-3 — 저장소에 넣지 않음, ref/fetch_drop.sh) */
const RD=path.join(__dirname,"drop_real");
if(fs.existsSync(path.join(RD,"water_in_air001.png"))){ console.log("== 4. 실측 OpenDrop: 물 펜던트 (바늘 0.7176 mm, ρ 1000, g 9.8 — OpenDrop 문서 값)");
  const G=[], R1=[]; for(let k=1;k<=5;k++){ const f=`water_in_air00${k}.png`; const im=readPNG(path.join(RD,f)); const t0=Date.now(); let r; try{ r=IC.pendantDrop(im,{needleMM:0.7176,rhoIn:1000,rhoOut:0,g:9.8}); }catch(e){ ok(false,f,e.message); continue; }
    G.push(r.gamma); R1.push(r); console.log(`    ${f}: γ ${r.gamma.toFixed(2)} ± ${r.gammaSd.toFixed(2)} mN/m · Bo ${r.Bo.toFixed(4)} · Wo ${r.Wo.toFixed(3)} · V ${r.V.toFixed(2)} mm³ · A ${r.A.toFixed(2)} mm² · R0 ${r.R0mm.toFixed(4)} mm · rms ${r.rms.toFixed(3)} px · ${Date.now()-t0} ms`); }
  if(G.length){ const m=G.reduce((a,b)=>a+b,0)/G.length; const r=R1[0];
    // OpenDrop 문서 결과 화면(docs/usage/images/ift_results.png): IFT 73.31 mN/m · V 11.15 mm³ · A 24.05 mm² · Wo 0.6609 · Bo 0.2157 · 영상 각 0.7502°
    const cmp=[["γ",r.gamma,73.31],["V",r.V,11.15],["A",r.A,24.05],["Wo",r.Wo,0.6609],["Bo",r.Bo,0.2157]].map(([k,a,b])=>[k,a,b,(a/b-1)*100]);
    ok(cmp.every(q=>Math.abs(q[3])<1.5),"001 ↔ OpenDrop (모두 ±1.5 %)",cmp.map(q=>`${q[0]} ${q[1].toFixed(4)}/${q[2]} (${q[3]>=0?"+":""}${q[3].toFixed(2)} %)`).join(" · "));
    ok(m>71.5&&m<74,`5 장 평균 ${m.toFixed(2)} mN/m — 물(20–25 °C) 72.0–72.8 mN/m (IAPWS R1-76) 근처`); }
  const cf=path.join(RD,"drop_on_surface.png"); if(fs.existsSync(cf)){ console.log("== 5. 실측 OpenDrop: 세실 방울 (OpenDrop conan 결과 왼 120.0° · 오른 124.2°) — 영역 x 100–620, y 330–760 (화면에서 사각형 ROI 를 긋는 것과 같음)"); const im=readPNG(cf); let r; try{ r=IC.contactAngle(im,{roi:{x0:100,y0:330,x1:620,y1:760}}); }catch(e){ ok(false,"drop_on_surface",e.message); }
    if(r){ for(const [k,v] of Object.entries(r.methods)) if(v&&v.mean!=null) console.log(`    ${k}: 왼 ${v.left.toFixed(1)}° · 오른 ${v.right.toFixed(1)}° · 평균 ${v.mean.toFixed(1)}°`); console.log(`    기준선 ${r.baseline.auto?"자동":"지정"} ${r.baseline.angleDeg.toFixed(2)}°`);
      const p=r.methods.poly; ok(p&&Math.abs(p.left-120.0)<4&&Math.abs(p.right-124.2)<4,`접선(다항식) ↔ OpenDrop (±4°)`,p&&`${p.left.toFixed(1)} / ${p.right.toFixed(1)}`); } } }
else console.log("실측 자료 없음 (bash ref/fetch_drop.sh) — 건너뜀");
console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
