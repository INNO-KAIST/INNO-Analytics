/* stall_stress.js 두 결과(옛 = XPS_NOSTALL, 새) 비교: χ² 비(새/옛) 분포, 나빠진 경우, 시간 (SPEC §3.3b 표) */
"use strict";
const fs=require("fs");
const rd=f=>fs.readFileSync(f,"utf8").split("\n").filter(Boolean).map(JSON.parse);
const A={}, B={}; rd(process.argv[2]).forEach(o=>A[o.seed]=o); rd(process.argv[3]).forEach(o=>B[o.seed]=o);
const rows=Object.keys(A).filter(k=>B[k]&&!A[k].err&&!B[k].err).map(k=>({a:A[k],b:B[k]}));
const q=(arr,p)=>{ const s=arr.slice().sort((x,y)=>x-y); if(!s.length) return NaN; const t=(s.length-1)*p, i=Math.floor(t); return s[i]+(s[Math.min(i+1,s.length-1)]-s[i])*(t-i); };
const fmt=x=>(x>=0?"+":"")+x.toExponential(1);
function report(name,rs){
  if(!rs.length){ console.log(name+": 0"); return; }
  const rel=rs.map(r=>(r.b.chi2-r.a.chi2)/r.a.chi2), reld=rs.map(r=>(r.b.chi2d-r.a.chi2d)/r.a.chi2d);
  const same=rs.filter((r,i)=>Math.abs(rel[i])<=1e-6).length, worse=rs.filter((r,i)=>rel[i]>1e-4), better=rs.filter((r,i)=>rel[i]<-1e-4);
  const ta=rs.reduce((s,r)=>s+r.a.ms,0), tb=rs.reduce((s,r)=>s+r.b.ms,0);
  const ia=rs.reduce((s,r)=>s+r.a.iter+r.a.iter1,0), ib=rs.reduce((s,r)=>s+r.b.iter+r.b.iter1,0);
  console.log(`${name}: ${rs.length} 피팅 · χ²(새/옛 − 1) 백분위 1/5/50/95/99 % = ${[0.01,0.05,0.5,0.95,0.99].map(p=>fmt(q(rel,p))).join(" / ")} · |Δ| ≤ 1e-6 ${same} · 나빠짐 >1e-4 ${worse.length} · 좋아짐 <−1e-4 ${better.length} · stalled ${rs.filter(r=>r.b.stop==="stalled").length} · kick 있음 ${rs.filter(r=>r.b.kicks>0).length} · 시간 ${(ta/1000).toFixed(1)} → ${(tb/1000).toFixed(1)} s · 반복 ${ia} → ${ib}`);
  worse.forEach(r=>console.log(`    나빠짐 seed ${r.a.seed} [${r.a.loose?"loose":"tight"} ${r.a.shapes} · ${r.a.bg} · ${r.a.weight}] χ² ${r.a.chi2.toPrecision(7)} → ${r.b.chi2.toPrecision(7)} (${fmt((r.b.chi2-r.a.chi2)/r.a.chi2)}; χ²_d ${fmt((r.b.chi2d-r.a.chi2d)/r.a.chi2d)}) · 반복 ${r.a.iter1}+${r.a.iter} → ${r.b.iter1}+${r.b.iter} · ${r.a.stop} → ${r.b.stop} ${r.b.flags?r.b.flags.join(","):""}`));
  return {rel,reld};
}
report("전체",rows); report("tight (BE·넓이·주 폭)",rows.filter(r=>!r.a.loose)); report("loose (모양 매개변수까지)",rows.filter(r=>r.a.loose));
const slow=rows.filter(r=>r.a.iter1+r.a.iter>=300);
report("옛 반복 ≥ 300 (기어감)",slow);
slow.sort((x,y)=>y.a.ms-x.a.ms).slice(0,12).forEach(r=>console.log(`    seed ${r.a.seed} ${r.a.loose?"loose":"tight"} ${r.a.shapes}: ${(r.a.ms/1000).toFixed(1)} s (${r.a.iter1}+${r.a.iter}) → ${(r.b.ms/1000).toFixed(1)} s (${r.b.iter1}+${r.b.iter}, ${r.b.stop}) · χ² ${fmt((r.b.chi2-r.a.chi2)/r.a.chi2)} · ${r.b.flags?r.b.flags.join(","):""}`));
