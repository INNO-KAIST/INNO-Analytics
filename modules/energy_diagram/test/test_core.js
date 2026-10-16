/* Energy Diagram 핵심 회귀 시험 — node test/test_core.js */
globalThis.window=globalThis; require('../src/ed_data.js'); const E=require('../src/ed_core.js');
let pass=0, fail=0; const ok=(c,m)=>{ if(c) pass++; else { fail++; console.log('FAIL',m); } };
const near=(a,b,t,m)=>ok(Math.abs(a-b)<=t,`${m}: ${a} vs ${b} (±${t})`);
const ctx={hw:360,T:295};
const L=(ion,p)=>E.annotate(E.levelsFor(Object.assign({ion},p||{})),ctx.hw,E.IONS[ion].kind);
const f=(Ls,id)=>Ls.find(l=>l.id===id);
const et=(dI,Di,Df,aI,Ai,Af,c)=>{ const LD=L(dI),LA=L(aI); return E.classifyET({d:{ion:dI,from:f(LD,Di),to:f(LD,Df)},a:{ion:aI,from:f(LA,Ai),to:f(LA,Af)}},c||ctx); };

/* 1. 데이터: 모든 Ln³⁺ 준위에 출처 코드, Stark 범위가 바리센터를 포함 */
Object.entries(window.ED_DATA.ln3).forEach(([el,d])=>d.lv.forEach(r=>{ ok('ecqf'.includes(r[5]),`${el} ${r[0]} src`); ok((r[3]<=r[2]+1&&r[4]>=r[2]-1)||(r[5]==="f"&&Math.min(Math.abs(r[3]-r[2]),Math.abs(r[4]-r[2]))<150),`${el} ${r[0]} stark∋E`); }));
/* 2. 문헌 대표 파장과 바리센터 파장: 2% 이내 */
window.ED_DATA.trans.forEach(t=>{ if(!t[4]) return; const Ls=window.ED_DATA.ln3[t[0]].lv; const u=Ls.find(r=>r[0]===t[1][0]), l=Ls.find(r=>r[0]===t[2]); const nm=1e7/(u[2]-l[2]); ok(Math.abs(nm-t[4])/t[4]<0.02,`${t[0]} ${t[1]}→${t[2]} ${nm.toFixed(1)} vs ${t[4]}`); });
/* 3. 화학식 해석 */
const P=s=>E.parseFormula(s);
let r=P("NaY0.78Yb0.2Er0.02F4"); ok(r.layers[0].hostKey==="NaYF4","host key NaYF4"); ok(r.ions.map(i=>i.key+i.layers[0].conc).join()==="Yb3+20,Er3+2","site fractions");
r=P("β-NaYF₄:Yb³⁺(20%),Er³⁺(2%)"); ok(r.ions.map(i=>i.key+i.layers[0].conc).join()==="Yb3+20,Er3+2","unicode + paren conc"); ok(r.layers[0].phase==="β","phase");
r=P("NaGdF4:Yb,Tm@NaGdF4:Tb"); ok(r.layers.length===2 && r.ions.some(i=>i.key==="Gd3+"&&i.layers.length===2),"core-shell Gd both layers");
ok(P("CaF2:Eu").ions[0].key==="Eu2+","Eu in divalent site → 2+ (with note)"); ok(P("CaF2:Eu3+").ions[0].key==="Eu3+","explicit charge wins");
ok(P("K2SiF6:Mn").ions[0].key==="Mn4+","Mn on Si4+ site"); ok(P("ZnF2:Mn").ions[0].key==="Mn2+","Mn on Zn2+ site");
ok(P("NaYF4:Yb0.2,Er0.02").ions.map(i=>i.layers[0].conc).join()==="20,2","fraction ≤1 → %");
ok(P("NaXF4:Er").warn.length>0,"unknown element warns"); ok(P("NaYF4:Qq").warn.length>0,"bad dopant warns");
ok(P("KPb2Cl5:Nd").layers[0].hostInfo.w===203,"KPb2Cl5 phonon");
/* 4. 기준 메커니즘 판정 (β-NaYF4, ħω 360) */
let c;
c=et('Yb3+','2F5/2','2F7/2','Er3+','4I15/2','4I11/2'); ok(c.type==="ET"&&c.grade===0,"Yb→Er ⁴I₁₁/₂ resonant ET");
c=et('Yb3+','2F5/2','2F7/2','Er3+','4I11/2','4F7/2'); ok(c.type==="ETU"&&c.grade===0,"Yb→Er ⁴F₇/₂ ETU resonant");
c=et('Yb3+','2F5/2','2F7/2','Er3+','4I13/2','4F9/2'); ok(c.type==="ETU"&&c.grade===3&&!c.endo,"Yb→Er ⁴F₉/₂ ETU phonon-assisted exo");
c=et('Er3+','4I11/2','4I15/2','Yb3+','2F7/2','2F5/2'); ok(c.type==="BET","Er→Yb BET");
c=et('Er3+','4F7/2','4I11/2','Yb3+','2F7/2','2F5/2'); ok(c.type==="BET","BET even if donor ends excited");
c=et('Tm3+','3H4','3F4','Tm3+','3H6','3F4'); ok(c.type==="CR"&&c.grade===2,"Tm CR 2-phonon"); near(c.dE,1324,5,"Tm CR ΔE");
c=et('Yb3+','2F5/2','2F7/2','Yb3+','2F7/2','2F5/2'); ok(c.type==="EM"&&c.dE===0,"Yb EM");
c=et('Yb3+','2F5/2','2F7/2','Tm3+','3H6','3H5'); ok(c.type==="ET"&&c.grade>=3&&!c.endo,"Yb→Tm ³H₅ large-mismatch exo");
/* 흡열: 방향을 뒤집으면 흡열 + Boltzmann < 1 */
c=et('Tm3+','3F4','3H6','Tm3+','3H4','3H5'); ok(c.endo===false||c.endo===true,"endo flag exists");
c=et('Er3+','4I11/2','4I15/2','Tm3+','3H6','3H5'); ok(c.grade>0,"Er→Tm mismatch>0");
const c2=E.classifyET({d:{ion:'Tm3+',from:f(L('Tm3+'),'3H5'),to:f(L('Tm3+'),'3H6')},a:{ion:'Er3+',from:f(L('Er3+'),'4I15/2'),to:f(L('Er3+'),'4I11/2')}},ctx);
ok(c2.endo && c2.boltz<1,"endothermic → Boltzmann factor");
/* 협동 */
const Tb=L('Tb3+'),Yb=L('Yb3+');
c=E.classifyET({d:{ion:'Yb3+',from:f(Yb,'2F5/2'),to:f(Yb,'2F7/2')},a:{ion:'Tb3+',from:Tb[0],to:f(Tb,'5D4(1)')},coop:true},ctx); ok(c.type==="CET"&&c.grade<=2,`cooperative Yb-Yb→Tb 5D4 (${c.gradeName})`);
/* 5. 광흡수 · PA */
const Tm=L('Tm3+');
let ph=E.photonCandidates({},Tm,[1064],ctx);
ok(ph.some(x=>x.c.type==="GSA"&&x.to.id==="3H5"&&x.c.grade===3),"Tm 1064 GSA weak nonresonant");
ok(ph.some(x=>x.c.type==="ESA"&&x.from.id==="3F4"&&x.to.id==="3F2"&&x.c.grade===0),"Tm 1064 ESA ³F₄→³F₂ resonant");
let pa=E.paLoops(0,Tm,[1064],ctx); ok(pa.length&&pa[0].M.id==="3F4"&&pa[0].H.id==="3H4","Tm 1064 PA loop found");
ok(E.paLoops(0,L('Er3+'),[980],ctx).length===0,"Er 980 (resonant GSA) is not PA");
ok(E.paLoops(0,Tm,[800],ctx).length===0,"Tm 800 (resonant GSA ³H₄) is not PA");
ph=E.photonCandidates({},L('Er3+'),[980],ctx); ok(ph.some(x=>x.from.id==="4I11/2"&&x.to.id==="4F7/2"),"Er 980 ESA ⁴I₁₁/₂→⁴F₇/₂");
/* 6. 다중포논 판정: ħω 가 커지면 경쟁·MPR 이 늘어난다 */
const reg=(hw)=>E.annotate(E.levelsFor({ion:'Er3+'}),hw,'ln3').filter(l=>l.regime==="rad").length;
ok(reg(203)>=reg(360)&&reg(360)>=reg(600)&&reg(600)>=reg(1000),"more radiative levels in low-phonon host");
/* 7. dⁿ: 루비 역산값 재현 · 닫힌 형식과 보간 비교 */
const cr=E.annotate(E.levelsFor({ion:'Cr3+'}),400,'d');
near(1e7/f(cr,'2E').E,694.3,1.5,"ruby 2E R-line nm"); near(f(cr,'4T2').E,18000,15,"ruby 4T2"); near(f(cr,'4T1(F)').E,24400,40,"ruby 4T1(F)");
let maxErr=0;
for(let k=0;k<300;k++){ const B=500+Math.random()*600, CB=3.5+Math.random()*2, Dq=B*(0.3+Math.random()*3.5);
  const Ls=E.tsEval(3,B,CB*B,Dq).levels; const t2=Ls.find(o=>o.lab==="4T2#0").E, t1=Ls.find(o=>o.lab==="4T1#0").E;
  const rr=0.5*Math.sqrt(225*B*B-180*B*Dq+100*Dq*Dq); maxErr=Math.max(maxErr,Math.abs(t2-10*Dq),Math.abs(t1-(7.5*B+15*Dq-rr))); }
ok(maxErr<40,`TS interpolation vs closed form max err ${maxErr.toFixed(1)} cm-1`);
const mn=E.annotate(E.levelsFor({ion:'Mn2+'}),400,'d'); ok(mn[0].id==="6A1"&&Math.abs(f(mn,'4A1(G)').E-f(mn,'4E(G)').E)<30,"d5 4A1/4E(G) degenerate (10B+5C)");
near(f(mn,'4A1(G)').E,10*800+5*3200,40,"d5 4A1(G)=10B+5C");
/* 8. 제안: Yb/Er 에서 ET Yb→Er ⁴I₁₁/₂ 과 ETU ⁴F₇/₂ 가 상위 */
const sug=E.suggestET([{ion:'Yb3+'},{ion:'Er3+'}],[L('Yb3+'),L('Er3+')],ctx,{Emax:30000,limit:30});
const idx=q=>sug.findIndex(s=>s.di===0&&s.ai===1&&s.Ai.id===q[0]&&s.Af.id===q[1]);
ok(idx(['4I15/2','4I11/2'])>=0&&idx(['4I15/2','4I11/2'])<5,"suggest Yb→Er ET top5"); ok(idx(['4I11/2','4F7/2'])>=0&&idx(['4I11/2','4F7/2'])<15,"suggest Yb→Er ETU top15");
console.log(`pass ${pass}, fail ${fail}`); process.exit(fail?1:0);
