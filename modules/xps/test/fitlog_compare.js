/* 피팅 기록 비교 (SPEC §3.3b 검증): node test/fitlog_compare.js <기준 폴더> <새 폴더> [허용 χ² 상대 1e-6] [허용 |Δp|/σ 1e-3]
   두 폴더의 같은 이름 *.jsonl 을 줄 번호(seq)로 맞대어: 같은 문제인지(N·자유 매개변수 이름·Σy) 확인 →
   χ² 상대 차, 매개변수 차 |Δp|/σ_기준 (σ 없으면 |Δp|/(|p|+1e-12) 를 따로), 반복 수, 정지 이유를 비교한다. */
"use strict";
const fs=require("fs"), path=require("path");
const [A,B]=process.argv.slice(2,4); const TC=+(process.argv[4]||1e-6), TP=+(process.argv[5]||1e-3);
/* 성분 id 는 무작위일 수 있어(newComp) 처음 나온 순서 번호로 바꿔 비교 */
const norm=r=>{ const m={}; let k=0; r.key=r.names.map(s=>{ const i=s.lastIndexOf("|"), id=s.slice(0,i); if(!(id in m)) m[id]=k++; return m[id]+"|"+s.slice(i+1); }).join(); return r; };
const rd=f=>fs.existsSync(f)?fs.readFileSync(f,"utf8").split("\n").filter(Boolean).map(l=>norm(JSON.parse(l))):[];
const files=fs.readdirSync(A).filter(f=>f.endsWith(".jsonl"));
const tot={n:0,same:0,better:0,worse:0,mism:0,stall:0,earlier:0,iterOld:0,iterNew:0}; const listed=[];
for(const f of files){
  const a=rd(path.join(A,f)), b=rd(path.join(B,f)); const s={n:0,same:0,better:0,worse:0,mism:0,stall:0,earlier:0,iterOld:0,iterNew:0,maxRel:0,maxDp:0,maxDpRel:0};
  const m=Math.min(a.length,b.length); if(a.length!==b.length) console.log(`  ${f}: 줄 수 다름 ${a.length} vs ${b.length}`);
  for(let i=0;i<m;i++){ const x=a[i], y=b[i]; s.n++;
    if(x.n!==y.n||x.key!==y.key||x.ys!==y.ys){ s.mism++; if(s.mism<=3) console.log(`  ${f}#${i}: 다른 문제(정렬 어긋남)`); continue; }
    const rel=(y.chi2-x.chi2)/Math.max(Math.abs(x.chi2),1e-300), rel1=(y.chi2p1-x.chi2p1)/Math.max(Math.abs(x.chi2p1),1e-300);
    let dp=0, dpr=0; x.p.forEach((v,j)=>{ const sg=x.sig&&x.sig[j]; const d=Math.abs(y.p[j]-v); if(sg>0&&isFinite(sg)) dp=Math.max(dp,d/sg); else dpr=Math.max(dpr,d/(Math.abs(v)+1e-12)); });
    const itO=(x.iter1||0)+(x.iter||0), itN=(y.iter1||0)+(y.iter||0); s.iterOld+=itO; s.iterNew+=itN;
    if(y.stop==="stalled") s.stall++;
    if(itN<itO) s.earlier++;
    const same=Math.abs(rel)<=TC&&Math.abs(rel1)<=TC&&dp<=TP&&dpr<=1e-6;
    if(same) s.same++; else if(rel<-TC||(Math.abs(rel)<=TC&&rel1<-TC)) s.better++; else s.worse++;
    s.maxRel=Math.max(s.maxRel,Math.abs(rel)); s.maxDp=Math.max(s.maxDp,dp); s.maxDpRel=Math.max(s.maxDpRel,dpr);
    if(!same||y.stop==="stalled") listed.push({f,i,np:x.np,n:x.n,rel,rel1,dp,dpr,itO,itN,stopO:x.stop,stopN:y.stop,stop1:y.stop1,flags:y.stall?y.stall.flags.map(q=>q.id+"|"+q.k+":"+q.why.join("+")).join(","):""}); }
  console.log(`${f}: ${s.n} 피팅 · 같음 ${s.same} · 더 좋음 ${s.better} · 나쁨 ${s.worse} · 어긋남 ${s.mism} · stalled ${s.stall} · 반복 줄어듦 ${s.earlier} · 반복 합 ${s.iterOld} → ${s.iterNew} · 최대 |Δχ²|/χ² ${s.maxRel.toExponential(1)} · 최대 |Δp|/σ ${s.maxDp.toExponential(1)} · (σ 없는) 최대 |Δp|/|p| ${s.maxDpRel.toExponential(1)}`);
  Object.keys(tot).forEach(k=>tot[k]+=s[k]);
}
console.log(`합계: ${tot.n} 피팅 · 같음 ${tot.same} · 더 좋음 ${tot.better} · 나쁨 ${tot.worse} · 어긋남 ${tot.mism} · stalled ${tot.stall} · 반복 줄어듦 ${tot.earlier} · 반복 합 ${tot.iterOld} → ${tot.iterNew}`);
if(listed.length){ console.log("같지 않거나 stalled 인 피팅:"); listed.slice(0,+process.env.LIST||60).forEach(q=>console.log(`  ${q.f}#${q.i} np ${q.np} N ${q.n}: Δχ²/χ² ${q.rel.toExponential(2)} (1차 ${q.rel1.toExponential(2)}) · |Δp|/σ ${q.dp.toExponential(1)} · |Δp|/|p|(σ 없음) ${q.dpr.toExponential(1)} · 반복 ${q.itO} → ${q.itN} · ${q.stopO||"-"} → ${q.stopN}${q.stop1?"(1차 "+q.stop1+")":""} ${q.flags}`)); }
