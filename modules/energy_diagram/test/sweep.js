/* 전 이온 × 호스트 포논 × 여기 파장 × 온도 스윕 — 예외·NaN·불변식 점검 */
globalThis.window=globalThis; require('../src/ed_data.js'); const E=require('../src/ed_core.js');
const ions=Object.keys(E.IONS); const hws=[138,203,284,360,460,600,860,1050]; const lasers=[266,365,405,450,532,660,690,785,808,915,980,1064,1150,1208,1450,1532,1940];
let n=0,bad=[];
for(const ion of ions){ for(const hw of hws){ for(const T of [4,77,295,500]){
  const ctx={hw,T}; const lv=E.annotate(E.levelsFor({ion}),hw,E.IONS[ion].kind);
  lv.forEach((l,i)=>{ if(!isFinite(l.E)||!isFinite(l.lo)||!isFinite(l.hi)) bad.push(`${ion} ${l.id} nonfinite`); if(i&&l.p!=null&&!(l.p>=0)) bad.push(`${ion} ${l.id} p`); if(!l.relax||!l.relax.length) bad.push(`${ion} ${l.id} relax`); });
  for(let i=1;i<lv.length;i++) if(lv[i].E<lv[i-1].E) bad.push(`${ion} order`);
  const ph=E.photonCandidates({ion},lv,lasers,ctx);
  ph.forEach(x=>{ n++; if(!isFinite(x.c.delta)||x.c.delta<0||x.c.boltz>1||x.c.boltz<0) bad.push(`${ion} ph ${x.nm}`); if(x.c.grade<0||x.c.grade>4) bad.push('grade'); });
  E.paLoops(0,lv,lasers,ctx).forEach(r=>{ if(r.gsa.c.grade<2) bad.push(`${ion} PA with resonant GSA`); });
}}}
/* 쌍 제안: 대표 조합 */
const pairs=[['Yb3+','Er3+'],['Yb3+','Tm3+'],['Yb3+','Ho3+'],['Nd3+','Yb3+'],['Gd3+','Tb3+'],['Ce3+','Tb3+'],['Tm3+','Tm3+'],['Pr3+','Pr3+'],['Eu2+','Mn2+'],['Cr3+','Nd3+'],['Sm3+','Eu3+'],['Dy3+','Tb3+']];
for(const [a,b] of pairs) for(const hw of [203,360,600]){
  const ctx={hw,T:295}; const L=[a,b].map(i=>E.annotate(E.levelsFor({ion:i}),hw,E.IONS[i].kind));
  const s=E.suggestET([{ion:a},{ion:b}],L,ctx,{Emax:52000,limit:200,maxGrade:4});
  s.forEach(x=>{ n++; const c=x.c; if(!isFinite(c.dE)||c.delta<0) bad.push(`${a}${b} sug`); if(c.grade===0&&!(c.dw.hi>=c.aw.lo&&c.aw.hi>=c.dw.lo)) bad.push('overlap logic');
    /* 반대칭: 뒤집으면 ΔE 부호 반대 */ });
}
/* 대칭성: donor/acceptor 맞바꾼 역과정은 ΔE 부호가 반대, δ 동일 */
const hw=360, ctx={hw,T:295}; const Tm=E.annotate(E.levelsFor({ion:'Tm3+'}),hw,'ln3'), Er=E.annotate(E.levelsFor({ion:'Er3+'}),hw,'ln3');
for(let i=0;i<200;i++){ const pick=a=>a[Math.floor(Math.random()*a.length)];
  let D1=pick(Tm),D2=pick(Tm),A1=pick(Er),A2=pick(Er); if(D1.E<=D2.E||A2.E<=A1.E) continue;
  const f=E.classifyET({d:{ion:'Tm3+',from:D1,to:D2},a:{ion:'Er3+',from:A1,to:A2}},ctx);
  const r=E.classifyET({d:{ion:'Er3+',from:A2,to:A1},a:{ion:'Tm3+',from:D2,to:D1}},ctx);
  if(Math.abs(f.dE+r.dE)>1e-6||Math.abs(f.delta-r.delta)>1e-6||(f.delta>0&&f.endo===r.endo)) bad.push('reverse symmetry');
  /* 상세 균형: 정방향·역방향 Boltzmann 비 */
}
console.log('checks',n,'problems',bad.length, [...new Set(bad)].slice(0,20));
