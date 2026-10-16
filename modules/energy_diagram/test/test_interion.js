// 이온 간 상호작용·Förster 효율 회귀 시험
global.window=global; require('../src/ed_data.js'); require('../src/ed_core.js'); const C=ED_CORE; let ok=0,bad=0; const t=(c,m)=>{ if(c) ok++; else { bad++; console.log('FAIL',m); } };
const cols=[{ion:'Yb3+'},{ion:'Tm3+'}], ctx={hw:350,T:295}; const L=cols.map(c=>C.annotate(C.levelsFor(c),ctx.hw,(C.IONS[c.ion]||{}).kind));
const r=C.interIon(cols,L,ctx,{maxGrade:4,Emax:40000});
t(r.every(x=>cols[x.di].ion!==cols[x.ai].ion),'같은 이온 쌍 제외');
const f=(d,a)=>r.find(x=>x.Di.t===d[0]&&x.Df.t===d[1]&&x.Ai.t===a[0]&&x.Af.t===a[1]);
const a=f(['2F5/2','2F7/2'],['3H6','3H5']), b=f(['3H5','3H6'],['2F7/2','2F5/2']);
t(a&&!a.c.endo&&a.c.type==='ET','Yb→Tm 3H6→3H5 발열 ET'); t(b&&b.c.endo&&b.boltz<0.01&&b.c.type==='BET','역방향 BET endothermic');
t(Math.abs(a.c.delta-b.c.delta)<1e-9,'정·역방향 δ 동일');
for(let i=1;i<r.length;i++) if(r[i].c.delta<r[i-1].c.delta-1e-9){ t(false,'δ 정렬'); break; }
const erfc=x=>{ let s=0,n=200000,h=12/n; for(let i=0;i<=n;i++){ const u=x+i*h; s+=(i==0||i==n?0.5:1)*Math.exp(-u*u);} return 2/Math.sqrt(Math.PI)*s*h; };
[0.1,1,1.9,2.1,3,6].forEach(x=>t(Math.abs(C.erfcx(x)/(Math.exp(x*x)*erfc(x))-1)<1e-6,'erfcx '+x));
const lo=C.forsterEff(1e-40,1e-3,1e18), hi=C.forsterEff(1e-40,1e-3,1e24);
t(Math.abs(lo.eta/(Math.sqrt(Math.PI)*lo.gamma)-1)<0.01,'η → √π γ (γ→0)'); t(hi.eta>0.99&&hi.eta<1,'η → 1 (γ→∞)');
console.log(`이온 간 상호작용 시험 ${ok} 통과 / ${bad} 실패`);
