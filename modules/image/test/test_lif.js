// 실측 Leica LIF 대조: 정답 = Bio-Formats 7.3.1 이 읽은 값 (liffile 과 같음을 ref/gen_lif.py 가 확인)
// 파일이 없으면 건너뜀: sh ref/fetch_lif.sh ref/real/lif && python3 -I ref/gen_lif.py ref/real/lif <bioformats_package.jar> <BFDump 폴더>
const fs=require("fs"), path=require("path");
const IO=require("../src/img_io.js"); require("../src/img_lif.js");
const D=path.join(__dirname,"..","ref","real","lif"), T=path.join(D,"truth");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push("✗ "+m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d&&!c?" — "+d:"")); };
(async()=>{
  // XML 나무 (속성 따옴표 안의 '>', 엔터티, 주석)
  { const t=IO._xmlTree('<?xml version="1.0"?><A x="1 &gt; 0" y=\'b\'><!-- <B/> --><B n="2"/><C><D>text</D></C></A>'); const A=t.c[0];
    ok(A.tag==="A"&&A.a.x==="1 > 0"&&A.a.y==="b"&&A.c.length===2&&A.c[0].tag==="B"&&A.c[1].c[0].tag==="D","XML 나무: 속성 엔터티·주석·중첩"); }
  if(!fs.existsSync(path.join(T,"truth.json"))){ console.log("실측 LIF·정답 없음 — 건너뜀"); }
  else { const truth=JSON.parse(fs.readFileSync(path.join(T,"truth.json"),"utf8")); const byFile={};
    for(const k of Object.keys(truth)){ const [f,s]=k.split("#"); (byFile[f]=byFile[f]||[]).push([+s,truth[k],k]); }
    for(const f of Object.keys(byFile).sort()){ const fp=path.join(D,f); if(!fs.existsSync(fp)){ console.log("  (없음) "+f); continue; }
      const u8=new Uint8Array(fs.readFileSync(fp)); const want=byFile[f].map(x=>x[0]); const t0=Date.now(); let r;
      try{ r=await IO.readFile(f,u8,{choose:async o=>{ if(o.positions!==byFile[f][0][1].series) throw new Error(`위치 수 ${o.positions} ≠ Bio-Formats series ${byFile[f][0][1].series}`); return {p:want}; }}); }catch(e){ ok(false,f,"예외: "+e.message); continue; }
      const ms=Date.now()-t0; ok(r.images.length===want.length,`${f}: series ${byFile[f][0][1].series} 중 ${want.length} 개 (${ms} ms)`,`앱 ${r.images.length}`);
      byFile[f].forEach(([s,t,k],i)=>{ const im=r.images[i]; if(!im) return; const b=fs.readFileSync(path.join(T,k+".truth")); const tr=new Float32Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
        let ex=""; if(im.w!==t.w||im.h!==t.h||im.nz!==t.nz||im.nc!==t.nc||im.nt!==t.nt) ex+=` 크기 앱 ${im.w}×${im.h} c${im.nc} z${im.nz} t${im.nt} / BF ${t.w}×${t.h} c${t.nc} z${t.nz} t${t.nt}`;
        else { const N=im.w*im.h; let bad=0; im.planes.forEach((P,q)=>{ for(let j=0;j<N;j++) if(P[j]!==tr[q*N+j]) bad++; }); if(bad) ex+=` 다른 화소 ${bad}`; }
        const um=v=>v*IO.unitToM(im.cal.unit)/1e-6;
        if(t.px){ if(!im.cal) ex+=" 보정 없음"; else if(Math.abs(um(im.cal.pw)-t.px)>1e-9*t.px||Math.abs(um(im.cal.ph)-t.py)>1e-9*t.py) ex+=` 화소 ${um(im.cal.pw)}×${um(im.cal.ph)} µm / BF ${t.px}×${t.py}`; }
        if(t.pz&&t.nz>1&&Math.abs((im.zStep||0)-t.pz)>1e-9*t.pz) ex+=` z step ${im.zStep} / BF ${t.pz}`;
        ok(!ex,`${k} ${t.w}×${t.h} c${t.nc} z${t.nz} t${t.nt}${t.tfold&&t.tfold.length>1?` (T = ${t.tfold.join("·")})`:""}`,ex); }); } }
  // 합성 경계 사례 (ref/gen_lif_synth.py — 정답은 Bio-Formats 또는 만든 자료)
  { const SD=path.join(__dirname,"io","lif"), SJ=path.join(SD,"synth.json");
    if(!fs.existsSync(SJ)) console.log("합성 LIF 없음 — 건너뜀 (python3 -I ref/gen_lif_synth.py <bf jar> <BFDump>)");
    else { const spec=JSON.parse(fs.readFileSync(SJ,"utf8"));
      for(const [f,sp] of Object.entries(spec)){ const u8=new Uint8Array(fs.readFileSync(path.join(SD,f))); let r, err=null, asked=null;
        try{ r=await IO.readFile(f,u8,{choose:async o=>{ asked=o; return sp.choose||{}; }}); }catch(e){ err=e.message; }
        if(err){ ok(false,`합성 ${f}`,"오류: "+err); continue; }
        let ex=""; const W=(r.warn||[]).join(" | ")+" | "+r.images.map(im=>(im.warn||[]).join(" | ")).join(" | ");
        if(sp.images!=null&&r.images.length!==sp.images) ex+=` 영상 ${r.images.length} ≠ ${sp.images}`;
        if(sp.warn&&!W.includes(sp.warn)) ex+=` 알림에 '${sp.warn}' 없음 (${W})`;
        if(sp.rgb&&!(r.images[0]&&r.images[0].rgb&&r.images[0].bit===8)) ex+=" RGB 8-bit 아님";
        if(sp.chNames&&JSON.stringify(r.images[0].chNames)!==JSON.stringify(sp.chNames)) ex+=` 채널 이름 ${JSON.stringify(r.images[0].chNames)}`;
        if(sp.names&&JSON.stringify(r.images.map(i=>i.name))!==JSON.stringify(sp.names)) ex+=` 이름 ${JSON.stringify(r.images.map(i=>i.name))}`;
        if(sp.pw!=null){ const c=r.images[0].cal; if(!c||Math.abs(c.pw-sp.pw)>1e-9||Math.abs(c.ph-sp.ph)>1e-9) ex+=` 보정 ${c&&c.pw}×${c&&c.ph}`; }
        if(sp.choose){ if(!asked||asked.nt!==20||!asked.posNt||asked.posNt.join()!=="1,20") ex+=` 고르기 창 nt ${asked&&asked.nt} posNt ${asked&&asked.posNt}`; if(sp.nt&&r.images.map(i=>i.nt).join()!==sp.nt.join()) ex+=` T ${r.images.map(i=>i.nt)}`; }
        if(sp.series&&!ex) sp.series.forEach((t,k)=>{ const im=r.images[k]; if(!im){ ex+=` series ${k} 없음`; return; } const b=fs.readFileSync(path.join(SD,t.key+".f32")); const tr=new Float32Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
          const N=im.w*im.h; let bad=0; if(im.w!==t.w||im.h!==t.h){ ex+=` 크기 ${im.w}×${im.h}`; return; } im.planes.forEach((P,q)=>{ for(let j=0;j<N;j++) if(P[j]!==tr[q*N+j]) bad++; }); if(bad) ex+=` series ${k} 다른 화소 ${bad}`; });
        ok(!ex,`합성 ${f}${sp.series?" = "+(sp.bf?"Bio-Formats":"만든 자료"):""}${sp.warn?" · 알림":""}`,ex); } } }
  // 실측: 채널 이름 = liffile
  for(const [f,names] of [["readlif_new_lasx.lif",["HOECHST 33258.","Nile Red triglycerid"]],["readlif_xz-example.lif",["FITC","TRITC"]],["bioio_tiled.lif",["C1","C2","C3","C4"]]]){ const fp=path.join(D,f); if(!fs.existsSync(fp)) continue;
    const r=await IO.readFile(f,new Uint8Array(fs.readFileSync(fp)),{choose:async()=>({p:[0]})}); ok(JSON.stringify(r.images[0].chNames)===JSON.stringify(names),`${f}: 채널 이름 = liffile (${names.join(", ")})`,JSON.stringify(r.images[0].chNames));
    ok(typeof r.images[0].metaText==="string"&&r.images[0].metaText.startsWith("<Element"),`${f}: metaText = 그 영상 요소의 XML`); }
  // 손상 파일
  { const fp=path.join(D,"readlif_xyzt-example.lif"); if(fs.existsSync(fp)){ const u8=new Uint8Array(fs.readFileSync(fp)); let m=null, r=null; try{ r=await IO.readFile("x.lif",u8.slice(0,u8.length-1000),{choose:async()=>({})}); }catch(e){ m=e.message; }
    const W=r?r.images.map(i=>(i.warn||[]).join(" ")).join(" "):""; ok(r&&r.images.length&&/잘려/.test(W),"실측 LIF 끝 1000 바이트 잘림 → 읽고 잘린 부분은 0 + 알림",m||W);
    m=null; try{ await IO.readFile("x.lif",u8.slice(0,200),{}); }catch(e){ m=e.message; } ok(m&&/잘린/.test(m),"XML 머리 안에서 잘림 → 분명한 오류",m); } }
  console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
})();
