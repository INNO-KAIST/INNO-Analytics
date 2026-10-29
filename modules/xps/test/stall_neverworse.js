/* 멈춤 규칙(SPEC §3.3b) '절대 나빠지지 않음' 회귀 시험 — 같은 문제를 같은 프로세스에서 규칙 끔(옛 LM, XPS_NOSTALL 과 같음)·켬으로 번갈아 피팅.
   node test/stall_neverworse.js [첫 seed=1] [개수=1100] [기록.jsonl]     (SEEDS=3,38,… 이면 그 seed 만)
   판정: (A) 최종 χ²(IRLS 2차 뒤, 코드가 정의한 목적함수) 새 ≤ 옛·(1 + 1e-7) — 모든 피팅
         (B) 갈래를 하나도 풀지 않은 피팅(trig = 0)은 χ²·매개변수가 비트 단위로 같음
   시간: 옛·새 합(같은 프로세스, seed 마다 순서를 바꿔 잼). 하나라도 어긋나면 종료 코드 1.
   전체(1–1100, 약 25 % loose)는 2 코어에서 반씩 나눠 약 25 분: node test/stall_neverworse.js 1 550 a.jsonl & node test/stall_neverworse.js 551 550 b.jsonl */
"use strict";
const fs=require("fs"), path=require("path");
const C=require(path.join(__dirname,"../src/xps_core.js")), {gen,run}=require("./stall_gen.js");
const s0=+(process.argv[2]||1), N=+(process.argv[3]||1100), out=process.argv[4]||null;
const seeds=process.env.SEEDS?process.env.SEEDS.split(",").map(Number):Array.from({length:N},(_,i)=>s0+i);
if(out) fs.writeFileSync(out,"");
const one=(on,g)=>{ C._hook.stall=on; try{ return run(C,g); } finally{ C._hook.stall=true; } };
let n=0, same=0, better=0, worse=[], untouched=0, notIdent=[], maxRel=-Infinity, ta=0, tb=0, ta1=0, tb1=0, acted=0;
for(const s of seeds){
  const g=gen(C,s); let a,b;
  if(s%2){ a=one(false,g); if(!a) continue; b=one(true,g); } else { b=one(true,g); if(!b) continue; a=one(false,g); }
  a=a.rec; b=b.rec; if(a.err||b.err){ if(!!a.err!==!!b.err) worse.push({s,why:"error "+(a.err||"")+" / "+(b.err||"")}); continue; }
  n++; ta+=a.ms; tb+=b.ms; const rel=b.chi2/a.chi2-1; maxRel=Math.max(maxRel,rel);
  if(b.chi2===a.chi2) same++; else if(rel<-1e-7) better++;
  if(rel>1e-7) worse.push({s,rel,a,b});
  if(!b.trig){ untouched++; if(b.chi2!==a.chi2||JSON.stringify(b.p)!==JSON.stringify(a.p)) notIdent.push(s); } else { acted++; ta1+=a.ms; tb1+=b.ms; }
  if(out) fs.appendFileSync(out,JSON.stringify({seed:s,old:a,new:b})+"\n");
}
const e=x=>(x>=0?"+":"")+x.toExponential(2);
console.log(`${n} 피팅 · χ² 비트 같음 ${same} · 좋아짐(< −1e-7) ${better} · 나빠짐(> +1e-7) ${worse.length} · 최대 새/옛 − 1 = ${e(maxRel)}`);
console.log(`갈래 없음(trig 0) ${untouched}: 비트 단위로 다름 ${notIdent.length}${notIdent.length?" — seed "+notIdent.join(","):""} · 갈래 있음 ${acted}`);
console.log(`시간 옛 ${(ta/1000).toFixed(1)} s → 새 ${(tb/1000).toFixed(1)} s (갈래 있는 피팅 ${(ta1/1000).toFixed(1)} → ${(tb1/1000).toFixed(1)} s)`);
worse.forEach(w=>console.log(`  나빠짐 seed ${w.s}: ${w.why||e(w.rel)} ${w.a?`[${w.a.shapes} · ${w.a.bg} · ${w.a.weight}${w.a.loose?" · loose":""}]`:""}`));
if(worse.length||notIdent.length) process.exit(1);
