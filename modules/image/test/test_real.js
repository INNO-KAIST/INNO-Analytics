// 실측 nd2·czi 대조: 정답 = Python nd2 패키지 / Zeiss libCZI(pylibCZIrw) 가 읽은 값 (ref/gen_real.py)
// 파일이 없으면 건너뜀: sh ref/fetch_real.sh ref/real && python3 -I ref/gen_real.py ref/real
const fs=require("fs"), path=require("path"), zlib=require("zlib");
const IO=require("../src/img_io.js"); require("../src/zstd_vendor.js"); require("../src/img_micro.js");
IO.inflate=async(u8,raw)=>new Uint8Array(raw?zlib.inflateRawSync(u8):zlib.inflateSync(u8));
const D=path.join(__dirname,"..","ref","real"), T=path.join(D,"truth");
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c){ pass++; } else { fail++; fails.push("✗ "+m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(d&&!c?" — "+d:""));};
(async()=>{
  if(!fs.existsSync(path.join(T,"truth.json"))){ console.log("실측 파일·정답 없음 — 건너뜀"); return; }
  const truth=JSON.parse(fs.readFileSync(path.join(T,"truth.json"),"utf8")); const byFile={};
  for(const k of Object.keys(truth)){ const [f,i]=k.split("#"); (byFile[f]=byFile[f]||[]).push([+i||0,truth[k],k]); }
  for(const f of Object.keys(byFile)){ const fp=path.join(D,f); if(!fs.existsSync(fp)){ console.log("  (없음) "+f); continue; }
    const u8=new Uint8Array(fs.readFileSync(fp)); const t0=byFile[f][0][1]; const t1=Date.now();
    if(t0.expect_error){ let e=null; try{ await IO.readFile(f,u8,{}); }catch(err){ e=err; } ok(e&&/v1\.0/.test(e.message),`${f}: 읽지 않는 형식을 분명히 알림`,e?e.message:"오류 없이 읽음"); continue; }
    let r; try{ r=await IO.readFile(f,u8,{choose:t0.choose?async()=>t0.choose:undefined}); }catch(e){ ok(false,f,"예외: "+(e.stack||e).toString().split("\n").slice(0,2).join(" | ")); continue; }
    const ms=Date.now()-t1;
    ok(r.images.length===byFile[f].length,`${f}: 영상 수 ${byFile[f].length}`,`앱 ${r.images.length}`);
    for(const [i,t,k] of byFile[f]){ const im=r.images[i]; if(!im){ ok(false,k,"영상 없음"); continue; }
      const _b=fs.readFileSync(path.join(T,k.replace(/\//g,"_")+".truth")); const tr=new Float32Array(_b.buffer.slice(_b.byteOffset,_b.byteOffset+_b.byteLength)); const N=im.w*im.h;
      let mx=0, bad=0; const flat=im.planes; const sameShape=flat.length*N===tr.length&&t.shape[1]===im.h&&t.shape[2]===im.w;
      if(sameShape) for(let q=0;q<flat.length;q++){ const P=flat[q]; for(let j=0;j<N;j++){ const d=Math.abs(P[j]-tr[q*N+j]); if(d>0){ bad++; if(d>mx) mx=d; } } }
      let ex="";
      if(t.nc&&im.nc!==t.nc) ex+=` nc ${im.nc}≠${t.nc}`; if(im.nz!==t.nz) ex+=` nz ${im.nz}≠${t.nz}`; if(im.nt!==t.nt) ex+=` nt ${im.nt}≠${t.nt}`; if(t.rgb&&!im.rgb) ex+=" RGB 아님";
      if(t.pw){ if(!im.cal) ex+=" 보정 없음"; else { const kk=IO.unitToM(im.cal.unit)/1e-6; if(Math.abs(im.cal.pw*kk-t.pw)>1e-6*t.pw) ex+=` 화소 ${im.cal.pw} ${im.cal.unit}≠${t.pw} µm`; } } else if(im.cal) ex+=` 정답엔 보정 없음, 앱 ${im.cal.pw} ${im.cal.unit}`;
      if(t.zstep&&!(Math.abs((im.zStep||0)-t.zstep)<1e-6*t.zstep)) ex+=` z step ${im.zStep}≠${t.zstep}`;
      if(t.names&&JSON.stringify(im.chNames)!==JSON.stringify(t.names)) ex+=` 이름 ${JSON.stringify(im.chNames)}≠${JSON.stringify(t.names)}`;
      ok(sameShape&&bad===0&&!ex,`${k}: ${im.w}×${im.h} ${im.rgb?"RGB":im.bit+"-bit"} c${im.nc} z${im.nz} t${im.nt}${t.pw?` · ${+t.pw.toPrecision(5)} µm`:""}${t.names?" · 이름":""} (${ms} ms)`,
        sameShape?`다른 화소 ${bad} (최대 ${mx})${ex}`:`모양 다름: 앱 ${flat.length}×${im.h}×${im.w}, 정답 ${t.shape}`+ex); }
  }
  console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
})();
