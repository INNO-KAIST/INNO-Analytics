// AFM·SPM 형식(GWY·GSF·Nanoscope·ARDF) 읽기 시험
// 정답: ref/gen_afm.py — Gwyddion 2.64 의 가져오기 모듈·gwyfile·magic-afm·손 계산 (앱과 독립)
//   합성: test/io/afm (저장소에 있음)   실측: ref/real/afm (sh ref/fetch_afm.sh 로 받고 python3 -I ref/gen_afm.py test/io/afm ref/real/afm)
const fs=require("fs"), path=require("path");
const IO=require("../src/img_io.js"); require("../src/img_afm.js");
let pass=0, fail=0; const fails=[];
const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push("✗ "+m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d&&!c?" — "+d:"")); };
async function rd(name,u8){ return IO.readFile(name,u8,{}); }
async function err(name,u8){ try{ await rd(name,u8); return null; }catch(e){ return e.message; } }

async function folder(D,label){
  const T=path.join(D,"truth"); if(!fs.existsSync(path.join(T,"truth.json"))){ console.log(`${label}: 정답 없음 — 건너뜀 (${D})`); return; }
  const truth=JSON.parse(fs.readFileSync(path.join(T,"truth.json"),"utf8"));
  console.log(`\n== ${label} (${Object.keys(truth).length} 파일)`);
  for(const f of Object.keys(truth).sort()){
    const t=truth[f], fp=path.join(D,f); if(!fs.existsSync(fp)){ console.log("  (없음) "+f); continue; }
    const u8=new Uint8Array(fs.readFileSync(fp));
    if(t.expect_error){ const m=await err(f,u8); ok(m&&m.includes(t.expect_error),`${f}: 읽지 않는 종류를 분명히 알림`,m||"오류 없이 읽음"); continue; }
    let r; const t0=Date.now(); try{ r=await rd(f,u8); }catch(e){ ok(false,f,"예외: "+e.message); continue; } const ms=Date.now()-t0;
    ok(r.images.length===t.images.length,`${f}: 영상 ${t.images.length} 개`,`앱 ${r.images.length} (${r.images.map(i=>i.name).join(" | ")})`);
    for(const w of t.warn||[]) ok(r.warn.some(x=>x.includes(w)),`${f}: 알림에 '${w}'`,JSON.stringify(r.warn));
    t.images.forEach((ti,j)=>{ const im=r.images[j]; if(!im) return;
      const b=fs.readFileSync(path.join(T,ti.file)); const tr=new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
      const nz=ti.nz||1, N=ti.w*ti.h; let ex="";
      if(im.w!==ti.w||im.h!==ti.h||im.nz!==nz) ex+=` 크기 앱 ${im.w}×${im.h}×${im.nz}, 정답 ${ti.w}×${ti.h}×${nz}`;
      let bad=0, worst=0, wi=-1;
      if(!ex){ const f=im.siScale; if(!(f>0)) ex+=" siScale 없음";
        for(let z=0;z<nz&&!ex;z++){ const P=im.planes[z]; for(let i=0;i<N;i++){ const a=P[i], e=tr[z*N+i]*f; if(a!==a&&e!==e) continue; const d=Math.abs(a-e)/Math.max(Math.abs(e),1e-300); if(!(d<=1.5e-7)&&!(a===0&&Math.abs(e)<1e-38)){ bad++; if(!(d<=worst)){ worst=d; wi=z*N+i; } } } }
        if(bad) ex+=` 다른 값 ${bad} (최악 상대 ${worst.toExponential(2)} @${wi}: 앱 ${im.planes[Math.floor(wi/N)][wi%N]} / 정답 ${tr[wi]*f})`; }
      // Gwyddion 은 괄호 등 '거친' 문자가 든 단위(예 log(Pa))를 버리고 단위 없음으로 둔다 — 앱은 파일의 글자를 이름표로 남긴다 (값 배율은 같음)
      const gwyDrops=!ti.uz&&/[!#$&()*,:;=?@\\[\]_`|{}]/.test(im.siUnit||"");
      if(ti.kind!=="mask"&&(im.siUnit||"")!==(ti.uz||"")&&!gwyDrops) ex+=` 단위 앱 '${im.siUnit}' 정답 '${ti.uz}'`;
      if(ti.uxy==="m"&&ti.xreal>0){ if(!im.cal) ex+=" 보정 없음"; else { const k=IO.unitToM(im.cal.unit); const px=im.cal.pw*k, py=im.cal.ph*k, ex0=ti.xreal/ti.w, ey0=ti.yreal/ti.h;
        if(!(Math.abs(px-ex0)<=1e-9*ex0&&Math.abs(py-ey0)<=1e-9*ey0)) ex+=` 화소 앱 ${px}×${py} m, 정답 ${ex0}×${ey0} m`; } }
      else if(!ti.uxy&&im.cal) ex+=` 정답엔 길이 단위 없음, 앱 보정 ${im.cal.pw} ${im.cal.unit}`;
      if(ti.kind==="brick"&&ti.uzz==="m"){ const zs=ti.zreal/nz*1e6; if(!(Math.abs((im.zStep||0)-zs)<=1e-9*zs)) ex+=` z step 앱 ${im.zStep} µm, 정답 ${zs}`; }
      ok(!ex,`${f}#${j} ${ti.kind} '${ti.title}' ${ti.w}×${ti.h}${nz>1?"×"+nz:""} ${ti.uz||""}${j===0?` (${ms} ms)`:""}`,ex); });
  }
}

(async()=>{
  /* ---------- 단위 해석 (Gwyddion 관례) ---------- */
  console.log("== SI 단위 해석");
  const U=s=>{ const o=IO._unitParse(s); return o.p10+" "+IO._unitStr(o); };
  for(const [s,e] of [["m","0 m"],["nm","-9 m"],["~m","-6 m"],["µm","-6 m"],["μm","-6 m"],["um","-6 m"],["nm/V","-9 m/V"],["deg","0 deg"],["°","0 deg"],["Pa","0 Pa"],["kPa","3 Pa"],["mV","-3 V"],["pA/V","-12 A/V"],["N/m","0 N/m"],["","0 "],["V/V","0 "],["Å","-10 m"],["m^2","0 m^2"],["mm","-3 m"],["cd","0 cd"],["log(Pa)","0 log(Pa)"],["\u00ba","0 deg"],["\u00ba/V","0 deg/V"]])
    ok(U(s)===e,`'${s}' → ${e}`,U(s));
  { const a=IO._unitParse("nm/V"), b=IO._unitParse("V"); const o={p10:a.p10+b.p10,u:Object.assign({},a.u)}; for(const k in b.u){ o.u[k]=(o.u[k]||0)+b.u[k]; if(!o.u[k]) delete o.u[k]; } ok(IO._unitStr(o)==="m"&&o.p10===-9,"nm/V × V = 1e-9 m"); }
  { const d=IO._dispUnit({p10:0,u:{Pa:1}},Float64Array.from([2e6,3e6,5e6])); ok(d.unit==="MPa"&&d.f===1e-6,"Pa 값 수 MPa → MPa", JSON.stringify(d));
    const n=IO._dispUnit({p10:0,u:{N:1}},Float64Array.from([2e-9,4e-9])); ok(n.unit==="nN"&&Math.abs(n.f-1e9)<1,"N 값 수 nN → nN",JSON.stringify(n));
    const m=IO._dispUnit({p10:0,u:{m:1}},Float64Array.from([1e-3])); ok(m.unit==="nm"&&m.f===1e9,"길이는 값과 상관없이 nm"); }
  /* ---------- Nanoscope 값 줄 해석 ---------- */
  console.log("== Nanoscope 머리 값 해석 (Gwyddion parse_value 와 같은 규칙)");
  { const v=IO._nsValue("@2:Z scale","V [Sens. Zsens] (0.006713867 V/LSB) 440.0000 V"); ok(v.soft==="Sens. Zsens"&&v.hscale===0.006713867&&v.hscaleU==="V"&&v.hv===440&&v.hvU==="V","@2:Z scale: V [soft] (hard V/LSB) value V");
    const s=IO._nsValue("@2:Image Data",'S [ZSensor] "Height Sensor"'); ok(s.soft==="ZSensor"&&s.hvStr==="Height Sensor","@2:Image Data: S [soft] \"이름\"");
    const o=IO._nsValue("Scan Size","7 7 ~m"); ok(o.hv===7&&o.hvStr==="7 7 ~m","옛 형식 값");
    const e=IO._nsValue("@Sens. Phase","V 1.000000 "); ok(e.hv===1&&!e.hvU,"soft scale 단위 없음 → 단위 null (Gwyddion: Phase 는 deg, 그 밖은 V)"); }
  /* ---------- CRC32 ---------- */
  ok(IO._crc32(new TextEncoder().encode("123456789"),0,9)===0xCBF43926,"CRC32('123456789') = CBF43926 (zlib 표준)");

  await folder(path.join(__dirname,"io","afm"),"합성 (정답: 손 계산·gwyfile·magic-afm, 모두 Gwyddion 과 대조)");
  await folder(path.join(__dirname,"..","ref","real","afm"),"실측 (정답: Gwyddion 2.64 · gwyfile)");

  /* ---------- 손상·변형 파일은 분명한 오류 ---------- */
  console.log("\n== 손상·변형 파일");
  const D=path.join(__dirname,"io","afm"); const L=n=>new Uint8Array(fs.readFileSync(path.join(D,n)));
  { const u=L("ns_v531.spm"); let m=await err("x.spm",u.slice(0,u.length-100)); ok(m&&/잘렸/.test(m),"Nanoscope 자료가 잘림 → 오류",m);
    const s=Buffer.from(u).toString("latin1").replace(/\r\n/g,"\n"); m=await err("x.spm",new Uint8Array(Buffer.from(s,"latin1"))); ok(m&&/줄 끝/.test(m),"CR LF → LF 로 망가진 Nanoscope → 원인을 알림",m);
    const t=new Uint8Array(u); t[0]=0x3F; m=await err("x.spm",t); ok(m&&/텍스트형/.test(m),"텍스트형(?*File list) → 읽지 않음을 알림",m); }
  { const u=new Uint8Array(L("a_fmap.ardf")); u[16+40]^=1; let m=await err("x.ardf",u); ok(m&&/CRC32/.test(m),"ARDF 목차 1 바이트 손상 → CRC32 오류",m);
    const v=L("a_fmap.ardf"); m=await err("x.ardf",v.slice(0,3000)); ok(m&&/(잘림|파일 밖|넘습니다)/.test(m),"ARDF 가 잘림 → 오류",m); }
  { const u=L("g_multi.gwy"); let m=await err("x.gwy",u.slice(0,u.length-500)); ok(m&&/GWY/.test(m),"GWY 가 잘림 → 오류",m); }
  { const u=L("s_a.gsf"); let m=await err("x.gsf",u.slice(0,u.length-4)); ok(m&&/크기/.test(m),"GSF 자료 크기가 맞지 않음 → 오류",m); }
  { const r=await rd("s_b.gsf",L("s_b.gsf")); ok(r.images[0].warn.some(w=>/NaN 1/.test(w)),"GSF NaN 은 그대로 두고 알림",JSON.stringify(r.images[0].warn)); }
  { const r=await rd("g_multi.gwy",L("g_multi.gwy")); const im=r.images[0]; ok(im.name.includes("높이")&&im.info.some(x=>x[0]==="Offset (x, y)"),"GWY: 한글 제목·offset 정보",im.name);
    ok(/Scan rate: 1.5 Hz/.test(im.metaText),"GWY: meta 를 원문으로 보여 줌"); ok(r.images[1].bit===8&&r.images[1].name.includes("mask"),"GWY: mask 를 8-bit 영상으로");
    const mod=r.images.find(i=>/Modulus/.test(i.name)); ok(mod&&mod.valueUnit==="MPa","GWY: Pa 채널 → MPa 로 표시",mod&&mod.valueUnit); }
  { const r=await rd("a_pad.ardf",L("a_pad.ardf")); ok(r.images[0].warn.some(w=>/비어 있는 줄 1/.test(w))&&/SecondNote/.test(r.images[0].metaText),"ARDF: 빈 IDAT 줄 → NaN·알림, 영상 메모 둘 다",JSON.stringify(r.images[0].warn));
    const m=await rd("ns_misc.spm",L("ns_misc.spm")); ok(m.warn.some(w=>/soft scale 'Sens. Nope'/.test(w))&&m.images[2].valueUnit==="deg","Nanoscope: 없는 soft scale 알림, º/V → deg",JSON.stringify(m.warn)+m.images[2].valueUnit);
    const u=Buffer.from(L("ns_v531.spm")).toString("latin1").replace(/\\Data offset: \d+/,"\\Data offset: n/a"); const e=await err("x.spm",new Uint8Array(Buffer.from(u,"latin1"))); ok(e&&/수가 아닙니다/.test(e),"Nanoscope: Data offset 이 수가 아니면 오류 (머리를 영상으로 읽지 않음)",e);
    const f=Buffer.from(L("ns_fv_magic.spm")).toString("latin1").replace("Operating mode: Force Volume","Operating mode: Force"); const e2=await err("x.spm",new Uint8Array(Buffer.from(f,"latin1"))); ok(e2&&/힘 곡선/.test(e2),"Force file list + Operating mode Force → 힘 곡선 오류",e2); }
  { const r=await rd("g_old.gwy",L("g_old.gwy")); ok(/Operator: Synthetic µm/.test(r.images[0].metaText),"GWYO: '/meta/…' 를 원문으로 (latin1 µ)",JSON.stringify(r.images[0].metaText)); }
  { const r=await rd("a_rect.ardf",L("a_rect.ardf")); const im=r.images[0]; ok(im.w===20&&im.h===12&&Math.abs(im.cal.ph/im.cal.pw-2)<1e-12,"ARDF 직사각(20×12, y 간격 2 배): IDEF 간격 사용 (Gwyddion 은 FastScanSize 로 정사각 가정)",`${im.w}×${im.h} ${im.cal.pw}×${im.cal.ph}`); }

  console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
})();
