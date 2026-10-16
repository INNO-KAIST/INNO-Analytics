// Olympus/Evident OIR 대조: 정답 = ref/gen_oir_synth.py 가 쓴 합성 파일의 정답 (oirfile 과 Bio-Formats 7.3.1 이 같게 읽은 값,
// 다르면 synth.json 의 truth 에 적은 쪽). 실측 시료가 없어 합성 파일만 — 파일이 없으면: python3 -I ref/gen_oir_synth.py <bf jar> <BFDump 폴더>
const fs=require("fs"), path=require("path"), zlib=require("zlib");
const IO=require("../src/img_io.js"); require("../src/img_lif.js"); require("../src/img_oir.js");
IO.inflate=async(u8,raw)=>new Uint8Array(raw?zlib.inflateRawSync(u8):zlib.inflateSync(u8));
const D=path.join(__dirname,"io","oir"), SJ=path.join(D,"synth.json");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push("✗ "+m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d&&!c?" — "+d:"")); };
const f32=k=>{ const b=fs.readFileSync(path.join(D,k+".f32")); return new Float32Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength)); };
const rd=f=>new Uint8Array(fs.readFileSync(path.join(D,f)));
const allWarn=r=>(r.warn||[]).join(" | ")+" | "+r.images.map(im=>(im.warn||[]).join(" | ")).join(" | ");
const rel=(a,b)=>Math.abs(a-b)<=1e-9*Math.max(Math.abs(a),Math.abs(b));
function cmpPlanes(im,tr,sel){ // sel: {z:[...], t:[...], nz, nt} → 정답에서 그 면만
  const N=im.w*im.h; let bad=0, first=-1; const nc=im.nc;
  const zs=sel?sel.z:[...Array(im.nz).keys()], ts=sel?sel.t:[...Array(im.nt).keys()], NZ=sel?sel.nz:im.nz;
  for(let ti=0;ti<ts.length;ti++) for(let zi=0;zi<zs.length;zi++) for(let c=0;c<nc;c++){ const P=im.planes[c+nc*(zi+zs.length*ti)], q=c+nc*(zs[zi]+NZ*ts[ti]);
    for(let j=0;j<N;j++) if(P[j]!==tr[q*N+j]){ bad++; if(first<0) first=q; } }
  return bad?`다른 화소 ${bad} (첫 면 ${first})`:""; }
(async()=>{
  if(!fs.existsSync(SJ)){ console.log("합성 OIR 없음 — 건너뜀 (python3 -I ref/gen_oir_synth.py <bf jar> <BFDump>)"); return; }
  const spec=JSON.parse(fs.readFileSync(SJ,"utf8"));
  ok(IO.sniff("a.oir",rd("gray12.oir"))==="oir"&&IO.sniff("split_00001",rd("split_00001"))==="oir"&&IO.sniff("pair.poir",rd("pair.poir"))==="poir"&&IO.sniff("x.zip",rd("pair.poir"))!=="poir","판별: OLYMPUSRAWFORMAT 머리(확장자 없는 보조 파일 포함)·.poir ZIP");
  for(const [f,sp] of Object.entries(spec)){
    if(sp.error||sp.images) continue;
    let r, err=null, asked=null;
    try{ r=await IO.readFile(f,rd(f),{choose:async o=>{ asked=o; return sp.choose||{}; }}); }catch(e){ err=e.message; }
    if(err){ ok(false,`${f} (${sp.desc})`,"오류: "+err); continue; }
    const im=r.images[0]; let ex="";
    const sel=sp.choose?{z:sp.choose.z,t:sp.choose.t,nz:sp.nz,nt:sp.nt}:null;
    const wantNz=sel?sel.z.length:sp.nz, wantNt=sel?sel.t.length:sp.nt;
    if(sp.rgb){ if(!(im.rgb&&im.bit===8&&im.nc===3)) ex+=` RGB 8-bit 아님 (rgb ${im.rgb} bit ${im.bit} nc ${im.nc})`; }
    if(im.w!==sp.w||im.h!==sp.h||im.nc!==sp.nc||im.nz!==wantNz||im.nt!==wantNt) ex+=` 크기 앱 ${im.w}×${im.h} c${im.nc} z${im.nz} t${im.nt} / 정답 ${sp.w}×${sp.h} c${sp.nc} z${wantNz} t${wantNt}`;
    else ex+=cmpPlanes(im,f32(sp.key),sel);
    if(im.bit!==sp.bit) ex+=` bit ${im.bit} ≠ ${sp.bit}`;
    if(sp.px){ if(!im.cal) ex+=" 보정 없음"; else { const um=im.cal.pw*IO.unitToM(im.cal.unit)/1e-6; if(!rel(um,sp.px)) ex+=` 화소 ${um} µm ≠ ${sp.px}`; } }
    const wz=sel&&sp.zStep?sp.zStep*(sel.z[1]-sel.z[0]):sp.zStep; if((wz||0)!==0||im.zStep){ if(!rel(im.zStep||0,wz||0)) ex+=` z step ${im.zStep} ≠ ${wz}`; }
    if(sp.chNames&&!sp.rgb&&JSON.stringify(im.chNames)!==JSON.stringify(sp.chNames)) ex+=` 채널 이름 ${JSON.stringify(im.chNames)} ≠ ${JSON.stringify(sp.chNames)}`;
    if(sp.warn&&!allWarn(r).includes(sp.warn)) ex+=` 알림에 '${sp.warn}' 없음 (${allWarn(r)})`;
    if(!sp.warn&&!sp.choose&&allWarn(r).replace(/[ |]/g,"")) ex+=` 뜻밖의 알림: ${allWarn(r)}`;
    if(allWarn(r).includes("아직 검증하지 않은 읽기")) ex+=" 일괄 '미검증' 알림이 남아 있음";
    if(sp.fi){ const fi=(im.info.find(x=>x[0]==="Frame interval")||[])[1]||""; if(!fi.startsWith(sp.fi)) ex+=` frame interval '${fi}' ≠ ${sp.fi}`; }
    if(sp.chLut&&JSON.stringify(im.chLut)!==JSON.stringify(sp.chLut)) ex+=` chLut ${JSON.stringify(im.chLut)} ≠ ${JSON.stringify(sp.chLut)}`;
    if(sp.choose&&!(asked&&asked.kind==="oir"&&asked.nz===sp.nz&&asked.nt===sp.nt&&asked.nc===sp.nc&&asked.w===sp.w)) ex+=` 고르기 창 ${JSON.stringify(asked)}`;
    if(sp.nometa){ if(im.cal) ex+=" 메타데이터 없는데 보정이 있음"; }
    else if(sp.refonly){ if(!/reference$/.test(im.name)) ex+=" 참조 영상 이름 아님"; }
    else if(!(typeof im.metaText==="string"&&im.metaText.includes("lsmimage:imageProperties")&&im.metaText.length<=200001)) ex+=" metaText 없음";
    if(sp.ref){ const ri=r.images[1]; if(!ri) ex+=" 참조 영상 없음"; else { if(ri.w!==sp.ref.w||ri.h!==sp.ref.h||ri.nc!==sp.ref.nc) ex+=` 참조 크기 ${ri.w}×${ri.h} c${ri.nc}`; else { const e2=cmpPlanes(ri,f32(sp.key+"#ref")); if(e2) ex+=" 참조 "+e2; } } }
    else if(r.images.length!==1) ex+=` 영상 ${r.images.length} 개`;
    ok(!ex,`${f} ${sp.w}×${sp.h} c${sp.nc} z${wantNz} t${wantNt} — ${sp.desc} [정답: ${sp.truth.split(" (")[0]}]`,ex); }
  // 오류가 나야 하는 파일
  for(const [f,sp] of Object.entries(spec)){ if(!sp.error) continue; let m=null; try{ await IO.readFile(f,rd(f),{choose:async()=>({})}); }catch(e){ m=e.message; }
    ok(m&&m.includes(sp.error),`${f} (${sp.desc}) → 분명한 오류`,m||"오류 없음"); }
  // POIR: 각 OIR = 단일 파일 정답
  for(const [f,sp] of Object.entries(spec)){ if(!sp.images||sp.budget) continue; const r=await IO.readFile(f,rd(f),{choose:async()=>({})}); let ex="";
    if(r.images.length!==sp.images.length) ex+=` 영상 ${r.images.length} ≠ ${sp.images.length}`;
    else sp.images.forEach((k,i)=>{ const im=r.images[i], t=spec[k]; if(im.w!==t.w||im.nc!==t.nc){ ex+=` ${i} 크기`; return; } const e=cmpPlanes(im,f32(k)); if(e) ex+=` ${i} ${e}`; if(im.name!==`${f} · ${sp.inner[i]}`) ex+=` 이름 '${im.name}'`; if(t.chNames&&JSON.stringify(im.chNames)!==JSON.stringify(t.chNames)) ex+=` 채널 이름 ${im.chNames}`; });
    ok(!ex,`${f}: ${sp.desc}`,ex); }
  { const f="q3_three.poir", sp=spec[f]; if(sp){ const one=64*48*4*2*4*3; const names=[], prog=[];
      const r=await IO.readFile(f,rd(f),{maxBytes:one,choose:async o=>{ names.push(o.name); return {}; },onProgress:p=>prog.push(p)});
      ok(r.images.length===1&&/메모리 한도/.test((r.warn||[]).join(" "))&&names[0]===`${f} · ${sp.inner[0]}`&&prog.every((p,i)=>i===0||p>=prog[i-1])&&Math.max(...prog)<=1,
        "q3_three.poir: 메모리 한도는 묶음 전체 (둘째부터 건너뜀), 고르기 창에 안쪽 이름, 진행 막대가 되돌아가지 않음",`${r.images.length} ${r.warn} ${names} ${prog.slice(0,5)}`);
      let c=0; const r2=await IO.readFile(f,rd(f),{choose:async()=>(++c===2?null:{})});
      ok(!r2.cancel&&r2.images.length===1&&!cmpPlanes(r2.images[0],f32("c2z4big.oir"))&&/취소/.test(r2.warn.join(" ")),"q3_three.poir: 둘째 고르기 창 취소 → 첫 영상은 남김",`${r2.cancel} ${r2.images.length} ${r2.warn}`);
      const r3=await IO.readFile(f,rd(f),{choose:async()=>null}); ok(r3.cancel===true&&!r3.images.length,"q3_three.poir: 첫 고르기 창 취소 → cancel"); } }
  // POIR: 암호화·ZIP64 크기 항목 → 분명한 알림
  { const u8=rd("pair.poir").slice(); const dv=new DataView(u8.buffer); let e=-1; for(let p=u8.length-22;p>=0;p--) if(dv.getUint32(p,true)===0x06054b50){ e=p; break; }
    let cd=dv.getUint32(e+16,true); const ent=[]; for(let k=0;k<3;k++){ ent.push(cd); cd+=46+dv.getUint16(cd+28,true)+dv.getUint16(cd+30,true)+dv.getUint16(cd+32,true); }
    dv.setUint16(ent[0]+8,dv.getUint16(ent[0]+8,true)|1,true); dv.setUint32(ent[1]+20,0xFFFFFFFF,true);
    let m=null; try{ await IO.readFile("pair.poir",u8,{}); }catch(err){ m=err.message; } ok(m&&/암호/.test(m)&&/ZIP64/.test(m),"POIR: 암호화 항목·ZIP64 크기 항목 → 분명한 오류",m); }
  // 취소·고르기·메모리 한도
  { const r=await IO.readFile("c2z4big.oir",rd("c2z4big.oir"),{choose:async()=>null}); ok(r.cancel===true&&r.images.length===0,"고르기 창 취소 → cancel"); }
  { let m=null; try{ await IO.readFile("c2z4big.oir",rd("c2z4big.oir"),{maxBytes:64*48*4*2*4*3-1}); }catch(e){ m=e.message; } ok(m&&/너무 큼/.test(m),"메모리 한도(ctx.maxBytes) 넘음 → 분명한 오류",m);
    const r=await IO.readFile("c2z4big.oir",rd("c2z4big.oir"),{maxBytes:64*48*4*2*4*3,choose:async()=>({t:[1]})}); const tr=f32("c2z4big.oir");
    ok(r.images[0].nt===1&&r.images[0].nz===4&&!cmpPlanes(r.images[0],tr,{z:[0,1,2,3],t:[1],nz:4,nt:3})&&/T 3 장 중 1 장/.test(allWarn(r)),"T 하나만 고르기 (z 는 전부)"); }
  // 잘린 파일: 끝에서 조금씩 자르며 — 오류면 분명한 문구, 아니면 알림
  { const u8=rd("c2z3t2.oir"); let okAll=true, bad=""; for(const cut of [1,9,17,100,400,1500,3000,u8.length-60,u8.length-48]){ let m=null, r=null; try{ r=await IO.readFile("x.oir",u8.slice(0,u8.length-cut),{}); }catch(e){ m=e.message; }
      if(m){ if(!/잘|머리|크기|메타데이터/.test(m)){ okAll=false; bad+=` [${cut}: ${m}]`; } } else if(!/잘|index/.test(allWarn(r))){ okAll=false; bad+=` [${cut}: 알림 없음]`; } }
    ok(okAll,"c2z3t2.oir 끝을 여러 길이로 자름 → 읽고 알림 또는 분명한 오류",bad); }
  // 손상: UID 가 PIXEL 과 짝이 아닌 블록 · 이상한 UID
  { const u8=rd("gray12.oir").slice(); const s=Buffer.from(u8).toString("latin1"); const k=s.indexOf("_0_1_"); u8[k+5]="Z".charCodeAt(0); // GUID 를 망가뜨림
    let m=null; try{ await IO.readFile("bad.oir",u8,{}); }catch(e){ m=e.message; } ok(m&&/UID/.test(m),"화소 블록 이름(UID) 손상 → 분명한 오류",m); }
  // ============ 실측 OIR (ref/real/oir, CC BY 4.0) — 정답 = ref/gen_oir_real.py (oirfile = Bio-Formats [= ome.tif], λ 는 oirfile) ============
  { const RD=path.join(__dirname,"..","ref","real","oir"), TJ=path.join(RD,"truth","truth.json"); const crypto=require("crypto");
    if(!fs.existsSync(TJ)) console.log("실측 OIR·정답 없음 — 건너뜀 (python3 -I ref/gen_oir_real.py ref/real/oir <bf jar> <BFDump>)");
    else { const truth=JSON.parse(fs.readFileSync(TJ,"utf8"));
      const sha=P=>{ const u=new Uint16Array(P.length); for(let i=0;i<P.length;i++) u[i]=P[i]; return crypto.createHash("sha1").update(Buffer.from(u.buffer)).digest("hex"); };
      for(const [f,t] of Object.entries(truth)){ const fp=path.join(RD,f); if(!fs.existsSync(fp)){ console.log("  (없음) "+f); continue; }
        const u8=new Uint8Array(fs.readFileSync(fp)); let r=null, err=null; const t0=Date.now();
        try{ r=await IO.readFile(f,u8,{choose:async()=>({})}); }catch(e){ err=e.message; }
        if(t.error){ ok(err&&err.includes(t.error),`실측 ${f}: 분명한 오류 (${t.decisions[0]})`,err||"오류 없음"); continue; }
        if(err){ ok(false,`실측 ${f}`,"오류: "+err); continue; }
        const im=r.images[0]; let ex="";
        if(im.w!==t.w||im.h!==t.h||im.nc!==t.nc||im.nz!==t.nz||im.nt!==t.nt) ex+=` 크기 앱 ${im.w}×${im.h} c${im.nc} z${im.nz} t${im.nt} / 정답 ${t.w}×${t.h} c${t.nc} z${t.nz} t${t.nt}`;
        else { let bad=0, first=-1; im.planes.forEach((P,q)=>{ if(sha(P)!==t.sha[q]){ bad++; if(first<0) first=q; } }); if(bad) ex+=` 다른 면 ${bad}/${t.sha.length} (첫 ${first})`; }
        if(im.bit!==t.bit) ex+=` bit ${im.bit}`;
        if(JSON.stringify(im.chNames)!==JSON.stringify(t.chNames)) ex+=` 채널 이름 ${JSON.stringify(im.chNames)} ≠ ${JSON.stringify(t.chNames)}`;
        { const um=im.cal?im.cal.pw*IO.unitToM(im.cal.unit)/1e-6:null; if(!rel(um||0,t.px)) ex+=` 화소 ${um} ≠ ${t.px}`; }
        if(t.nz>1?!(Math.abs(im.zStep-t.zStep)<=1e-6*t.zStep):im.zStep) ex+=` z step ${im.zStep} ≠ ${t.zStep}`;
        if(t.fi){ const fi=parseFloat((im.info.find(x=>x[0]==="Frame interval")||[])[1]||"NaN"); if(!(Math.abs(fi-t.fi)<=1e-4*t.fi)) ex+=` frame interval ${fi} ≠ ${t.fi}`; }
        if(t.nl>1&&!/λ/.test((im.info.find(x=>x[0]==="λ")||[])[0]||"")) ex+=" λ 정보 없음";
        if(t.ref){ const ri=r.images[1]; if(!ri||ri.w!==t.ref.w||ri.h!==t.ref.h||ri.nc!==t.ref.nc) ex+=" 참조 영상 크기"; else ri.planes.forEach((P,q)=>{ if(sha(P)!==t.ref.sha[q]) ex+=` 참조 면 ${q} 다름`; }); }
        const W=allWarn(r); const wantAx=f.startsWith("coupe"); // 계획 200 장 중 35 장 (실측 파일의 ZSTACK 축 maxSize 200)
        if(/획득 설정보다 면이 적습니다/.test(W)!==wantAx) ex+=` '획득 설정보다 적음' 알림 ${wantAx?"없음":"잘못 뜸"}: ${W}`;
        ok(!ex,`실측 ${f.slice(0,48)} ${t.w}×${t.h} c${t.nc} z${t.nz} t${t.nt}${t.nl>1?` (λ${t.nl}·T${t.nt0})`:""}${t.ref?" + 참조":""} [${t.truth}] ${Date.now()-t0} ms`,ex); }
      // 실측 파일의 채널 색 (채널 GUID 가 붙은 LUT XML 의 이름)
      for(const [f,want] of [["1202-interval_10sec_sequence_frame.oir",["Green","Red"]],["tiling-sample-brain-section_A01_G001_0001.oir",["Cyan","Red"]]]){ const fp=path.join(RD,f); if(!fs.existsSync(fp)) continue;
        const r=await IO.readFile(f,new Uint8Array(fs.readFileSync(fp)),{}); ok(JSON.stringify(r.images[0].chLut)===JSON.stringify(want),`실측 ${f.slice(0,40)}: 채널 색 = LUT XML 이름 ${want}`,JSON.stringify(r.images[0].chLut)); }
      // 잘린 실측 파일: 끝 5 MB 잘림 → block index 없음 → 차례 읽기 (실측 파일은 METADATA 블록이 앞·끝에 둘)
      { const fp=path.join(RD,"venus stack.oir"); if(fs.existsSync(fp)){ const u8=new Uint8Array(fs.readFileSync(fp)); const t=truth["venus stack.oir"]; let m=null, r=null;
          try{ r=await IO.readFile("v.oir",u8.slice(0,u8.length-5e6),{}); }catch(e){ m=e.message; }
          let good=0; if(r){ const im=r.images[0]; im.planes.forEach((P,q)=>{ if(q<t.sha.length&&sha(P)===t.sha[q]) good++; }); }
          ok(r&&/index/.test(allWarn(r))&&/잘렸/.test(allWarn(r))&&good>=60&&r.images[0].chNames[0]==="CH3","실측 venus stack 끝 5 MB 잘림 → 차례 읽기 (앞쪽 METADATA 블록은 남음) · 앞 면들 = 정답 + 알림",m||`${good} 면 같음; ${allWarn(r).slice(0,200)}`); } } } }
  console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
})();
