/* 화학식 농도 점검 — node test/test_formula_warn.js */
const C=require('../src/ed_core.js');
let pass=0,fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
const conc=(r,el)=>{ const i=r.ions.find(x=>x.el===el); return i?i.layers[0].conc:undefined; };
const has=(r,re)=>r.warn.some(w=>re.test(w));
let r=C.parseFormula("NaYF4:20%Yb,2%Er"); ok(!r.warn.length&&conc(r,"Yb")===20&&conc(r,"Er")===2,'normal');
r=C.parseFormula("NaYF4:200%Yb"); ok(has(r,/100 % 초과/),'>100 %');
r=C.parseFormula("NaYF4:0%Yb"); ok(has(r,/0 이하/),'0 %');
r=C.parseFormula("NaYF4:80%Yb,30%Er"); ok(has(r,/같은 자리 도펀트 합 110/),'site sum');
r=C.parseFormula("NaYF4:Yb1.5,Er0.02"); ok(conc(r,"Yb")===1.5&&conc(r,"Er")===0.02&&has(r,/혼용/),'mixed bare numbers → all %');
r=C.parseFormula("NaYF4:0.2Yb,0.02Er"); ok(conc(r,"Yb")===20&&conc(r,"Er")===2&&!r.warn.length,'all fractions');
r=C.parseFormula("NaYF4:20%Yb,0.5Er"); ok(conc(r,"Er")===0.5&&has(r,/혼용/),'explicit % + bare ≤1');
r=C.parseFormula("NaYbF4:2%Er"); ok(!r.warn.length,'host-site ion not summed');
r=C.parseFormula("NaYF4:Yb,Er"); ok(!r.warn.length&&conc(r,"Yb")==null,'no concentration');
r=C.parseFormula("β-NaGdF4:49%Yb,1%Tm@NaGdF4:15%Tb"); ok(!r.warn.length,'core-shell');
console.log(`formula warning tests ${pass} pass / ${fail} fail`); if(fail) process.exit(1);
