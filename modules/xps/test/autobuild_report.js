/* Auto build 벤치마크 결과 표 (SPEC §18.5): node test/autobuild_report.js <out.jsonl> [...] — autobuild_bench.js 의 줄들을 묶어 구간별 표.
   지표: 선 수 맞음(참 · 해상 가능) · 가짜 선(스펙트럼에 참에 없는 선이 하나라도) · 놓친 선(간격·넓이 비별) · 확립된 경우의 위치/넓이 3σ ·
   해상 불가 참 묶음(간격 < 0.5 FWHM)을 정직하게 다룸 · '확립 안 됨' 비율 · 시간. */
"use strict";
const fs=require("fs");
function load(files){ const R=[]; files.forEach(f=>fs.readFileSync(f,"utf8").split("\n").filter(Boolean).forEach(l=>R.push(JSON.parse(l)))); return R; }
const pct=(a,b)=>b?(100*a/b).toFixed(1)+" %":"—";
function summarize(R){ const ok=R.filter(r=>r.ev); const n=ok.length; const s={n,err:R.length-n};
  s.countOK=ok.filter(r=>r.ev.countOK).length; s.countOKres=ok.filter(r=>r.ev.countOKres).length;
  s.fa=ok.filter(r=>r.ev.falseAdd>0).length; s.faLines=ok.reduce((t,r)=>t+r.ev.falseAdd,0); s.fitLines=ok.reduce((t,r)=>t+r.ev.nFit,0);
  s.faStart=ok.filter(r=>r.startFalseAdd>0).length;
  s.est=ok.filter(r=>r.ev.established).length; const E=ok.filter(r=>r.ev.established&&r.ev.posOK!=null);
  // 기준 B (SPEC §18.4): 확립된 경우, 해상 가능한 참 선에 짝지어진 따로 정해진 선 모두 |Δ| ≤ 3u (u = σ + 숨은 선 한계/3); 엄격(σ 만)도 함께
  const resZ=r=>r.ev.zs.filter(z=>z.res!==false); const EB=E.filter(r=>resZ(r).length);
  s.estPos=EB.length; s.posOK=EB.filter(r=>resZ(r).every(z=>z.zbU<=3)).length; s.areaOK=EB.filter(r=>resZ(r).every(z=>z.zaU<=3)).length;
  s.posStrict=EB.filter(r=>resZ(r).every(z=>z.zb<=3)).length;
  s.faEst=ok.filter(r=>r.ev.established&&r.ev.falseAdd>0).length;
  const U=[].concat(...ok.map(r=>r.ev.unres)); s.unres=U.length; s.unresOK=U.filter(u=>u.ok).length;
  const ms=ok.map(r=>r.ms).sort((a,b)=>a-b); s.msMed=ms.length?ms[ms.length>>1]:null; s.msP95=ms.length?ms[Math.floor(0.95*(ms.length-1))]:null; s.msMax=ms.length?ms[ms.length-1]:null;
  const Z=[].concat(...ok.filter(r=>r.ev.established).map(r=>r.ev.zs.filter(z=>z.res!==false))); s.zPos=Z.length; s.zPos3=Z.filter(z=>z.zb<=3).length; s.zPos3U=Z.filter(z=>z.zbU<=3).length;
  s.recipe=ok.filter(r=>r.meta.recipe).length; s.latlo=ok.filter(r=>r.ev.latloFound!=null); s.recExtra=ok.reduce((t,r)=>t+(r.ev.recipeExtra||0),0);
  return s; }
function row(name,s){ return `| ${name} | ${s.n} | ${pct(s.countOK,s.n)} | ${pct(s.countOKres,s.n)} | ${pct(s.fa,s.n)} (${s.faLines}/${s.fitLines} lines) | ${pct(s.faStart,s.n)} | ${pct(s.est,s.n)} | ${pct(s.posOK,s.estPos)} (${s.estPos}) | ${pct(s.posStrict,s.estPos)} | ${pct(s.areaOK,s.estPos)} | ${pct(s.zPos3U,s.zPos)} / ${pct(s.zPos3,s.zPos)} | ${s.unres?pct(s.unresOK,s.unres)+` (${s.unres})`:"—"} | ${s.msMed} / ${s.msP95} / ${s.msMax} |`; }
const HDR="| bin | spectra | count = truth | count = resolvable | false-add spectra | false-add after Auto suggest | established | B: est., all pos. ≤ 3u | est., all pos. ≤ 3σ (strict) | est., all areas ≤ 3u | est. lines: pos. ≤ 3u / ≤ 3σ | unresolvable pairs handled | ms median / p95 / max |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|";
function missTable(R){ const ok=R.filter(r=>r.ev); const M=[].concat(...ok.filter(r=>r.ev.K>1).map(r=>r.ev.miss)); const sb=[[0,0.5],[0.5,1],[1,1.5],[1.5,3.01],[3.01,1e9]], rb=[[0,0.05],[0.05,0.1],[0.1,0.3],[0.3,1.01]];
  const L=["| nearest neighbour (FWHM) \\ area ratio | "+rb.map(b=>`${b[0]}–${b[1]>1?1:b[1]}`).join(" | ")+" |","|---|"+rb.map(()=>"---").join("|")+"|"];
  sb.forEach(s=>{ L.push(`| ${s[0]}–${s[1]>100?"∞":s[1]>3?3:s[1]} | `+rb.map(b=>{ const a=M.filter(m=>m.sep>=s[0]&&m.sep<s[1]&&m.ratio>=b[0]&&m.ratio<b[1]); return a.length?`${pct(a.filter(m=>m.missed).length,a.length)} (${a.length})`:"—"; }).join(" | ")+" |"); }); return L.join("\n"); }
function report(R){ const out=["### overall",HDR,row("all",summarize(R)),row("in model family (symmetric)",summarize(R.filter(r=>!r.meta.asym))),row("asymmetric truth (outside family)",summarize(R.filter(r=>r.meta.asym)))];
  const by=(title,key,vals)=>{ out.push("",`### by ${title}`,HDR); vals.forEach(v=>{ const sub=R.filter(r=>!r.meta.asym&&key(r)===v); if(sub.length) out.push(row(String(v),summarize(sub))); }); };
  by("closest separation (symmetric only)",r=>r.meta.sepBin,["single","0.3-0.5","0.5-1","1-1.5","1.5-3"]);
  by("noise (symmetric only)",r=>r.meta.noise,["poisson","ar1","smooth"]);
  by("peak height, counts (symmetric only)",r=>r.meta.H,[300,3000,30000]);
  by("number of lines (symmetric only)",r=>r.meta.K,[1,2,3,4,5,6]);
  by("background slope (symmetric only)",r=>r.meta.slope?"slope":"no slope",["slope","no slope"]);
  out.push("","### missed-line rate (symmetric, K ≥ 2): by nearest-neighbour separation and area ratio to the largest line",missTable(R.filter(r=>!r.meta.asym)));
  const RC=R.filter(r=>r.meta.recipe); if(RC.length){ out.push("","### recipe starts (C 1s / O 1s): false-add = lines added by Auto build that are not in truth",HDR); ["c1s","o1s"].forEach(k=>{ const a=RC.filter(r=>r.meta.recipe===k); if(a.length) out.push(row("recipe "+k,summarize(a))); });
    const ll=RC.filter(r=>r.ev&&r.ev.latloFound!=null); out.push("",`O 1s with a low-BE lattice line (truth 0.8–1.2 eV below the lattice line, 10–30 %): found ${ll.filter(r=>r.ev.latloFound).length}/${ll.length}; recipe lines left with no truth line: ${RC.reduce((t,r)=>t+(r.ev?r.ev.recipeExtra||0:0),0)}`); }
  return out.join("\n"); }
module.exports={load,summarize,report};
if(require.main===module) console.log(report(load(process.argv.slice(2))));
