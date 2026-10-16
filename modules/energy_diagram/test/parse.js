globalThis.window=globalThis; require('../src/ed_data.js'); const E=require('../src/ed_core.js');
const cases=["NaYF4:20%Yb,2%Er","β-NaYF₄:Yb³⁺(20%),Er³⁺(2%)","NaY0.78Yb0.2Er0.02F4","Na(Y0.78Yb0.20Er0.02)F4","NaLuF4:20%Tm","NaGdF4:Yb,Tm@NaGdF4:Tb","CaF2:Eu","CaF2:Eu3+","K2SiF6:Mn","ZnF2:Mn","Al2O3:Cr","LiYF4:Tm 3%","KPb2Cl5:Nd","NaYbF4:Tm 1%@NaYF4","Y2O3:Eu 5 mol%","NaYF4:Yb0.2,Er0.02","SrF2:Sm2+","NaYF4:Yb/Tm/Er","α-NaYF4:Ce","Y3Al5O12:Ce","NaYF4:20Yb,2Er","NaXF4:Er","NaYF4:Qq"];
for(const c of cases){ const r=E.parseFormula(c);
 console.log(c.padEnd(30), '| layers:', r.layers.map(l=>`${l.tag}:${l.hostKey}${l.hostInfo?'('+l.hostInfo.w+')':''}${l.phase||''}`).join(' / '),
  '| ions:', r.ions.map(i=>`${i.key}${i.known?'':'?'}${i.layers.map(x=>'['+x.role+' '+(x.conc==null?'-':x.conc)+(x.fromHost?'H':'')+']').join('')}`).join(' '), r.warn.length?'| WARN '+r.warn.join(' ; '):'', r.ions.some(i=>i.notes.length)?'| NOTE '+r.ions.map(i=>i.notes.join(';')).join(' '):'' ); }
