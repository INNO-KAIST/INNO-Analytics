// 스케일바 이름표 OCR 시험 (node test/test_ocr.js)
//  1) 합성 (ref/gen_ocr_test.py → test/ocr): 조건별 정확도, conf 보정 (conf ≥ 0.9 → 틀림 ≤ 1 %)
//  2) 실측 Fe-Cr ODS SAD: "5 1/nm" (위쪽 "SAD 1" 이 아니라)
//  3) 속도: 2048² 에서 < 50 ms     4) 막대 없음·글자 없음·쓰레기 → null, 예외 없음
//  OCR_DIR=다른폴더 로 dev 묶음 평가, OCR_DUMP=파일 로 표본별 결과 저장 (conf 맞춤용)
const fs=require("fs"), path=require("path"), zlib=require("zlib");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); require("../src/diff_core.js"); const DP=require("../src/diff_pattern.js"); const OCR=require("../src/diff_ocr.js");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d?" — "+d:"")); };

/* 작은 PNG 해독기 (8-bit 회색·RGB·RGBA, 비월 없음) → Float32 회색 (Rec.601) */
function readPNG(file){ const b=fs.readFileSync(file); let p=8, w=0, h=0, ct=0, bd=0; const idat=[];
  while(p<b.length){ const len=b.readUInt32BE(p), type=b.toString("ascii",p+4,p+8), d=b.subarray(p+8,p+8+len);
    if(type==="IHDR"){ w=d.readUInt32BE(0); h=d.readUInt32BE(4); bd=d[8]; ct=d[9]; } else if(type==="IDAT") idat.push(d); else if(type==="IEND") break; p+=12+len; }
  if(bd!==8||![0,2,4,6].includes(ct)) throw new Error("PNG 형식 "+ct+"/"+bd);
  const ch={0:1,2:3,4:2,6:4}[ct], bpp=ch, st=w*ch; const raw=zlib.inflateSync(Buffer.concat(idat)); const px=new Uint8Array(w*h*ch); let prev=new Uint8Array(st);
  for(let y=0;y<h;y++){ const f=raw[y*(st+1)], line=raw.subarray(y*(st+1)+1,(y+1)*(st+1)); const cur=new Uint8Array(st);
    for(let i=0;i<st;i++){ const a=i>=bpp?cur[i-bpp]:0, up=prev[i], c=i>=bpp?prev[i-bpp]:0; let v=line[i];
      if(f===1) v+=a; else if(f===2) v+=up; else if(f===3) v+=(a+up)>>1; else if(f===4){ const pp=a+up-c, pa=Math.abs(pp-a), pb=Math.abs(pp-up), pc=Math.abs(pp-c); v+=(pa<=pb&&pa<=pc)?a:(pb<=pc?up:c); }
      cur[i]=v&255; } px.set(cur,y*st); prev=cur; }
  const img=new Float32Array(w*h); for(let i=0;i<w*h;i++){ if(ch<3) img[i]=px[i*ch]; else img[i]=0.299*px[i*ch]+0.587*px[i*ch+1]+0.114*px[i*ch+2]; }
  return {img,w,h}; }
// 자동 = DIFF_PATTERN.scaleBar 를 힌트로 OCR.readAuto (힌트 막대 + 스스로 찾은 밝은/어두운 막대 중 conf 최고)
const auto=(img,w,h,opt)=>OCR.readAuto(img,w,h,DP.scaleBar(DP.findOverlay(img,w,h),w,h),opt);

/* ---------- 1. 합성 ---------- */
const DIR=process.env.OCR_DIR||path.join(__dirname,"ocr");
const truth=fs.existsSync(path.join(DIR,"truth.json"))?JSON.parse(fs.readFileSync(path.join(DIR,"truth.json"),"utf8")):null;
const dump=[];
if(!truth) console.log("합성 시험 그림 없음 (python3 ref/gen_ocr_test.py) — 건너뜀");
else { console.log("== 1. 합성 이름표 "+truth.length+" 장 ("+DIR+")");
  const gmsRight=new Set(); const by={}; let nOK=0, nOKT=0, nHi=0, nHiWrong=0, nBarMiss=0; const bins=[[0,0.5],[0.5,0.7],[0.7,0.9],[0.9,1.01]].map(r=>({r,n:0,ok:0}));
  const t0=Date.now();
  for(const t of truth){ const {img,w,h}=readPNG(path.join(DIR,t.file));
    const isRight=r=>!!r&&Math.abs(r.value-t.value)<1e-9&&r.unit===t.unit&&(t.valueInvNm==null?r.valueInvNm==null:Math.abs(r.valueInvNm-t.valueInvNm)<1e-9);
    const r=auto(img,w,h,{debug:!!process.env.OCR_DUMP}); const bar=r&&r.bar;
    const barOK=!!bar&&Math.abs(bar.x0-t.bar.x0)<=3&&Math.abs(bar.x1-t.bar.x1)<=3&&Math.abs(bar.y0-t.bar.y0)<=3;
    const rT=OCR.readScaleLabel(img,w,h,t.bar); // 정답 막대를 줬을 때 (OCR 만)
    const right=isRight(r), rightT=isRight(rT);
    const c=by[t.cond]||(by[t.cond]={n:0,ok:0,okT:0,nul:0,hiWrong:0,bar:0}); c.n++; if(right) c.ok++; if(rightT) c.okT++; if(!r) c.nul++; if(r&&!barOK) c.bar++;
    if(right&&t.cond==="gms") gmsRight.add(t.file);
    if(right) nOK++; if(rightT) nOKT++; const conf=r?r.conf:0; if(conf>=0.9){ nHi++; if(!right){ nHiWrong++; c.hiWrong++; } }
    for(const bn of bins) if(conf>=bn.r[0]&&conf<bn.r[1]){ bn.n++; if(right) bn.ok++; }
    if(process.env.OCR_DUMP) dump.push({file:t.file,cond:t.cond,truth:t.text,got:r?r.text:null,right,conf,barOK,dbg:r&&r.debug});
    if(!right&&process.env.OCR_VERBOSE) console.log("    ✗ "+t.file+" "+t.text+" ["+t.font+" "+t.size+" "+t.place+(t.sup?" sup":"")+"] → "+(r?r.text+" conf "+r.conf.toFixed(2):"null")+(r&&!barOK?" (다른 막대)":"")+" | 정답 막대: "+(rT?rT.text+" "+rT.conf.toFixed(2):"null")); }
  const ms=(Date.now()-t0)/truth.length;
  console.log("   조건        n   자동 정확   정답막대 정확  null  다른막대  conf≥0.9 틀림");
  for(const k of Object.keys(by)){ const c=by[k]; console.log("   "+k.padEnd(10)+String(c.n).padStart(4)+(100*c.ok/c.n).toFixed(1).padStart(9)+" %"+(100*c.okT/c.n).toFixed(1).padStart(11)+" %"+String(c.nul).padStart(7)+String(c.bar).padStart(8)+String(c.hiWrong).padStart(10)); }
  console.log("   전체      "+String(truth.length).padStart(4)+(100*nOK/truth.length).toFixed(1).padStart(9)+" %"+(100*nOKT/truth.length).toFixed(1).padStart(11)+" %");
  console.log("   conf 구간: "+bins.map(b=>"["+b.r[0]+","+Math.min(1,b.r[1])+") "+b.n+" 장 "+(b.n?(100*b.ok/b.n).toFixed(1):"-")+" %").join(" · "));
  console.log("   그림당 평균 "+ms.toFixed(1)+" ms (PNG 해독·막대 찾기·정답 막대 한 번 더 포함)");
  const clean=by.clean; if(clean) ok(clean.ok/clean.n>=0.98,`깨끗한 조건 정확(자동 막대) ${(100*clean.ok/clean.n).toFixed(1)} % (≥ 98 %)`);
  ok(nOK/truth.length>=0.95,`전체 정확(자동 막대) ${(100*nOK/truth.length).toFixed(1)} % (≥ 95 %)`);
  ok(nHi>0&&nHiWrong/nHi<=0.01,`conf ≥ 0.9: ${nHi} 장 중 틀림 ${nHiWrong} (${(100*nHiWrong/Math.max(1,nHi)).toFixed(2)} % ≤ 1 %)`);
  // GMS 기본 이름표 꼴(검은 명조 + 흰 그림자, REDI 예제 B3·T13 모양): 큰 글자(대문자 ≥ 20 px)는 대부분 읽고, 못 읽으면 conf 를 낮춰 확인을 묻는다(높은 conf 오답 0)
  const g=by.gms; if(g){ const big=truth.filter(t=>t.cond==="gms"&&t.capH>=20); const okBig=big.filter(t=>gmsRight.has(t.file)).length;
    ok(g.ok/g.n>=0.55,`GMS 꼴 전체 ${(100*g.ok/g.n).toFixed(1)} % (≥ 55 %; 대문자 10–19 px 가 어렵다)`);
    ok(big.length>0&&okBig/big.length>=0.75,`GMS 꼴 대문자 ≥ 20 px: ${okBig}/${big.length} (${(100*okBig/Math.max(1,big.length)).toFixed(1)} % ≥ 75 %)`);
    ok(g.hiWrong===0,`GMS 꼴 conf ≥ 0.9 오답 ${g.hiWrong} (= 0)`); }
  if(process.env.OCR_DUMP) fs.writeFileSync(process.env.OCR_DUMP,JSON.stringify(dump)); }

/* ---------- 2. 실측 ---------- */
const RF=path.join(__dirname,"..","ref","real","FeCr_ODS_SAD_zenodo19217212.tif.f32");
if(fs.existsSync(RF)){ console.log("== 2. 실측 Fe-Cr ODS SAD (1024², 8-bit RGB)"); const b=fs.readFileSync(RF); const img=new Float32Array(b.buffer,b.byteOffset,b.length/4); const w=1024,h=1024;
  const bar=DP.scaleBar(DP.findOverlay(img,w,h),w,h); const t0=process.hrtime.bigint(); const r=OCR.readScaleLabel(img,w,h,bar); const dt=Number(process.hrtime.bigint()-t0)/1e6;
  ok(r&&r.value===5&&r.unit==="1/nm"&&r.valueInvNm===5&&r.recip===true,`"${r&&r.text}" → ${r&&r.value} ${r&&r.unit}, valueInvNm ${r&&r.valueInvNm}, conf ${r&&r.conf.toFixed(3)} (${dt.toFixed(1)} ms)`);
  ok(r&&r.box.y1<bar.y0&&r.box.y0>900&&!/SAD/.test(r.text),`막대 바로 위 이름표 상자 (${r&&[r.box.x0,r.box.y0,r.box.x1,r.box.y1].join(",")}), "SAD 1" 아님`);
  // 막대를 영상 위쪽 "SAD 1" 근처에 있다고 속여도 SAD 1 을 수·단위로 읽지 않는다
  const fake={x0:120,x1:300,y0:118,y1:121}; const r2=OCR.readScaleLabel(img,w,h,fake); ok(!r2||r2.conf<0.5,`가짜 막대를 "SAD 1" 아래에 두면 → ${r2?r2.text+" conf "+r2.conf.toFixed(2):"null"}`); }
else console.log("실측 자료 없음 (ref/gen_real.py) — 건너뜀");

/* ---------- 3. 속도 ---------- */
console.log("== 3. 속도 (2048², 이름표는 막대 둘레만)");
{ const W=2048, H=2048; const img=new Float32Array(W*H); let s=12345; const rnd=()=>((s=(s*1103515245+12345)>>>0)/4294967296);
  for(let i=0;i<W*H;i++) img[i]=20+rnd()*8; for(let k=0;k<400;k++){ const cx=rnd()*W, cy=rnd()*H; for(let y=-6;y<=6;y++) for(let x=-6;x<=6;x++){ const X=Math.round(cx+x), Y=Math.round(cy+y); if(X>=0&&Y>=0&&X<W&&Y<H) img[Y*W+X]+=200*Math.exp(-(x*x+y*y)/8); } }
  // 합성 묶음의 큰 글씨 하나를 아래 왼쪽에 붙인다
  let src=null; if(truth){ const t=truth.find(t=>t.cond==="large")||truth[0]; src=Object.assign(readPNG(path.join(DIR,t.file)),{t}); }
  let bar; if(src){ const ox=40, oy=H-src.h-40; for(let y=0;y<src.h;y++) for(let x=0;x<src.w;x++) img[(oy+y)*W+ox+x]=src.img[y*src.w+x]; bar={x0:src.t.bar.x0+ox,x1:src.t.bar.x1+ox,y0:src.t.bar.y0+oy,y1:src.t.bar.y1+oy}; }
  else { bar={x0:100,x1:700,y0:1950,y1:1965}; for(let y=bar.y0;y<=bar.y1;y++) for(let x=bar.x0;x<=bar.x1;x++) img[y*W+x]=255; }
  OCR.readScaleLabel(img,W,H,bar); const n=10; const t0=process.hrtime.bigint(); let r; for(let i=0;i<n;i++) r=OCR.readScaleLabel(img,W,H,bar); const dt=Number(process.hrtime.bigint()-t0)/1e6/n;
  ok(dt<50,`2048² 한 번 ${dt.toFixed(1)} ms (< 50 ms) → ${r?r.text:"null"}${src?" (정답 "+src.t.text+")":""}`); }

/* ---------- 4. 견고성 ---------- */
console.log("== 4. 막대 없음 · 글자 없음 · 쓰레기");
{ const W=512, H=512; const noise=new Float32Array(W*H); let s=7; const rnd=()=>((s=(s*1103515245+12345)>>>0)/4294967296); for(let i=0;i<W*H;i++) noise[i]=rnd()*255;
  const flat=new Float32Array(W*H).fill(10); const barOnly=Float32Array.from(flat); for(let y=480;y<486;y++) for(let x=50;x<250;x++) barOnly[y*W+x]=255;
  const nanImg=Float32Array.from(barOnly); for(let i=0;i<nanImg.length;i+=3) nanImg[i]=NaN;
  const B={x0:50,x1:249,y0:480,y1:485}; const cases=[["막대 null",()=>OCR.readScaleLabel(barOnly,W,H,null)],["글자 없는 막대",()=>OCR.readScaleLabel(barOnly,W,H,B)],["잡음 영상 + 막대",()=>OCR.readScaleLabel(noise,W,H,B)],
    ["NaN 섞인 영상",()=>OCR.readScaleLabel(nanImg,W,H,B)],["막대가 영상 밖",()=>OCR.readScaleLabel(barOnly,W,H,{x0:600,x1:900,y0:10,y1:12})],["막대 NaN",()=>OCR.readScaleLabel(barOnly,W,H,{x0:NaN,x1:5,y0:1,y1:2})],
    ["빈 배열",()=>OCR.readScaleLabel(new Float32Array(0),W,H,B)],["img null",()=>OCR.readScaleLabel(null,W,H,B)],["평평한 영상",()=>OCR.readScaleLabel(flat,W,H,B)],["findBar 잡음",()=>OCR.findBar(noise,W,H)],["findBar 평평",()=>OCR.findBar(flat,W,H)]];
  for(const [name,f] of cases){ let r, err=null; try{ r=f(); }catch(e){ err=e; } const isNull=r==null||(name.startsWith("잡음")&&r.conf<0.5);
    ok(!err&&(isNull||name.startsWith("findBar")&&r==null),`${name} → ${err?"예외 "+err.message:r==null?"null":r.text+" conf "+r.conf.toFixed(2)}`); } }

console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ fails.forEach(f=>console.log("  ✗ "+f)); process.exit(1); }
