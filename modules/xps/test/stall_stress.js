/* LM 멈춤 규칙(SPEC §3.3b) 무작위 비교 — 같은 합성 스펙트럼(stall_gen.js)을 한 엔진으로 피팅해 χ²·반복·시간을 줄마다 적는다.
   node test/stall_stress.js <엔진.js> <출력.jsonl> [첫 seed=1] [개수=320]   (XPS_NOSTALL=1 이면 멈춤 규칙을 끈 옛 동작)
   비교: node test/stall_stress_cmp.js <옛.jsonl> <새.jsonl>   (옛·새를 한 프로세스에서 번갈아 재려면 test/stall_neverworse.js) */
"use strict";
const fs=require("fs"), path=require("path");
const C=require(path.resolve(process.argv[2])); if(process.env.XPS_NOSTALL) C._hook.stall=false;
const out=process.argv[3]; const s0=+(process.argv[4]||1), N=+(process.argv[5]||320);
const {gen,run}=require("./stall_gen.js");
fs.writeFileSync(out,"");
for(let s=s0;s<s0+N;s++){ const x=run(C,gen(C,s)); if(x) fs.appendFileSync(out,JSON.stringify(x.rec)+"\n"); }
