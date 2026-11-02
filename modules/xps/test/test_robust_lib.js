/* 원소 라이브러리 재현성·출처 대조 (2026-11 감사 #5): node test/test_robust_lib.js
   1) 발행 Scofield 표(scofield.al/.mg) = 논문 전사(test/ref/scofield1976_transcription.json) 칸 하나하나
   2) doublet 합(예: "Ce 3d") = 인쇄된 부껍질 값의 합(하나만 인쇄되면 그 값), 합 칸은 전사에 따로 없음
   3) 문서대로의 빌드(python3 -I build/build_lib.py)가 커밋된 src/xps_lib.js 를 바이트 단위로 다시 만든다 (XPS_LIB_OUT 로 임시 파일에) */
"use strict";
const fs=require("fs"), path=require("path"), os=require("os"), cp=require("child_process");
const ROOT=path.join(__dirname,".."), LIBF=path.join(ROOT,"src","xps_lib.js");
const L=require(LIBF), T=JSON.parse(fs.readFileSync(path.join(__dirname,"ref","scofield1976_transcription.json"),"utf8"));
let pass=0, fail=0; const fails=[];
function ok(c,name,info){ if(c) pass++; else { fail++; fails.push(name); } if(process.env.V||!c) console.log((c?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }

for(const s of ["al","mg"]){
  const lib=L.scofield[s], tr=T[s]; let same=0, diff=[];
  for(const k in tr){ if(lib[k]===tr[k]) same++; else diff.push(`${k}: lib ${lib[k]} vs paper ${tr[k]}`); }
  ok(diff.length===0&&same===Object.keys(tr).length,`Scofield ${s}: every transcribed subshell equals the library`,`${same} values${diff.length?"; "+diff.slice(0,5).join("; "):""}`);
  const grp={}; for(const k in tr){ const [el,sub]=k.split(" "); if(sub.length>2) (grp[el+" "+sub.slice(0,2)]=grp[el+" "+sub.slice(0,2)]||[]).push(tr[k]); }
  let nT=0, bad=[]; for(const g in grp){ nT++; const t=+grp[g].reduce((a,b)=>a+b,0).toFixed(6); if(lib[g]!==t) bad.push(`${g}: ${lib[g]} vs Σ ${t}`); if(g in tr) bad.push(`${g} given explicitly in the transcription`); }
  ok(bad.length===0,`Scofield ${s}: doublet totals = sum of printed subshells`,`${nT} totals${bad.length?"; "+bad.slice(0,5).join("; "):""}`);
  const extra=Object.keys(lib).filter(k=>!(k in tr)&&!(k in grp)); ok(extra.length===0,`Scofield ${s}: no value without a source`,extra.slice(0,5).join(", "));
  ok(lib["C 1s"]===1,`Scofield ${s}: C 1s = 1`);
}
// 두 부껍질이 모두 인쇄된 doublet 의 합 예: Ce 3d (Al) = 3d3/2 + 3d5/2
ok(L.scofield.al["Ce 3d"]===+(T.al["Ce 3d3/2"]+T.al["Ce 3d5/2"]).toFixed(6),"Ce 3d (Al) = 3d3/2 + 3d5/2",`${L.scofield.al["Ce 3d"]}`);
// 하나만 인쇄된 예: Al Tm 3d — 3d3/2 (1515 eV) 는 1487 eV 위
ok(T.al["Tm 3d3/2"]==null&&L.scofield.al["Tm 3d"]===T.al["Tm 3d5/2"],"Al Tm 3d total = 3d5/2 only (3d3/2 above hν)");

// 빌드 재현: 임시 파일에 커밋본을 복사 → 그 자리에 빌드 → 바이트 비교 (built 날짜는 내용이 같으면 유지)
{ const py=["python3","python"].find(c=>{ try{ cp.execFileSync(c,["--version"],{stdio:"ignore"}); return true; }catch(e){ return false; } });
  if(!py) ok(false,"python3 available for the build check");
  else { const tmp=path.join(fs.mkdtempSync(path.join(os.tmpdir(),"xpslib-")),"xps_lib.js"); fs.copyFileSync(LIBF,tmp);
    let out=""; try{ out=cp.execFileSync(py,["-I",path.join(ROOT,"build","build_lib.py")],{env:Object.assign({},process.env,{XPS_LIB_OUT:tmp}),encoding:"utf8"}); }catch(e){ out=String(e.stdout||"")+String(e.stderr||e.message); }
    const a=fs.readFileSync(LIBF), b=fs.existsSync(tmp)?fs.readFileSync(tmp):Buffer.alloc(0);
    ok(a.equals(b),"build/build_lib.py reproduces the committed src/xps_lib.js byte for byte",`${a.length} vs ${b.length} bytes`+(a.equals(b)?"":"; "+out.split("\n").slice(-3).join(" | ")));
    fs.rmSync(path.dirname(tmp),{recursive:true,force:true}); } }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
