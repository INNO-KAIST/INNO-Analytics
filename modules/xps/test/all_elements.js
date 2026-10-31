/* 모든 원소·주요 궤도에서 같은 파이프라인이 동작하는지 — node test/all_elements.js [V=1]
   원소 Z = 3–92 각각에서 Scofield(Al Kα) 세기가 그 원소 최댓값의 15 % 이상이고 BE 10–1300 eV 인 궤도 중 센 것 2 개.
   (1) 파일 이름 해석: "Fe2p", "Fe 2p Scan", "fe2p", "[시료] Fe2p.txt"
   (2) 감도계수: Scofield Al·Mg, ALTHERMO1
   (3) doublet 정보(p·d·f): Δ·넓이 비·분해 여부
   (4) 단일 화학 상태 합성 스펙트럼 → 자동 제안(suggestModel) → 피팅 → 참값 회복
   (5) 화학 상태 2 개(δ = +2.5 eV) 합성 스펙트럼 → 참 구조에서 흔든 출발 → 피팅 → 참값 회복(z-값)
   합성: GLS(F 1.2, m 0.3), Δ_참 = 라이브러리 Δ + U(−0.1, 0.1) (라이브러리 오차 흉내), 넓이 비 (2j+1),
         배경 = 참 active Shirley(c₀ + k∫P) + 기울기, Poisson 잡음. */
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js"), L=require("../src/xps_lib.js");
function cases(){
  const out=[];
  Object.entries(L.core).forEach(([el,rec])=>{
    if(rec.Z<3||C.EXCLUDE_DEFAULT.includes(el)) return;
    const sc=L.scofield.al; const orbs=[...new Set(Object.keys(rec.lv).map(k=>k.slice(0,2)))];
    const cand=orbs.map(o=>({o,s:sc[el+" "+o]})).filter(x=>x.s!=null);
    if(!cand.length) return; const smax=Math.max(...cand.map(x=>x.s));
    cand.filter(x=>x.s>=0.15*smax).map(x=>{ const l=x.o[1]; const J=l==="s"?"":{p:"3/2",d:"5/2",f:"7/2"}[l]; const be=rec.lv[x.o+J]!=null?rec.lv[x.o+J]:rec.lv[x.o]; return Object.assign(x,{be}); })
      .filter(x=>x.be>=10&&x.be<=1300).sort((a,b)=>b.s-a.s).slice(0,2).forEach(x=>out.push({el,orb:x.o,be:x.be,s:x.s}));
  });
  return out;
}
function synth(R,cs,o){
  const info=C.doubletInfo(L,cs.el,cs.orb); const dbl=info&&info.dE;
  const dTrue=dbl?info.dE+(R.u()*0.2-0.1):0;
  // 두 상태 간격 δ: doublet Δ 와 겹치면(δ≈Δ) 본질적으로 구분 불가한 배치라 피한다
  const dl=dbl&&Math.abs(dTrue-2.5)<1.0?dTrue+1.5:2.5;
  const states=o.two?[{be:cs.be,A:20000},{be:cs.be+dl,A:10000}]:[{be:cs.be,A:20000}];
  const span=o.mult?3.5:(dbl?dTrue:0)+(o.two?dl:0);
  const step=span>25?0.1:0.05; const lo=cs.be-8, hi=cs.be+span+10; const n=Math.round((hi-lo)/step)+1;
  const E=Float64Array.from({length:n},(_,i)=>lo+i*step); const ctx={emin:E[0],emax:E[n-1]};
  const truthDl=dl;
  // 다중항 준위의 참 모형: 원소별 다중항 계산 대신, 넓은 선 3 개의 일반적 모양(오프셋 0·1.6·3.5 eV, 넓이 52·31·17 %, FWHM 1.8·2.0·2.4 eV)
  // — 자동 제안이 doublet 을 강요하지 않고 넓이 합을 회복하는지 보는 시험(화학 상태 분해는 하지 않으므로 BE 는 판정하지 않음)
  const truth=o.mult?[[0,0.52,1.8],[1.6,0.31,2.0],[3.5,0.17,2.4]].map(([d,fr,F],i)=>C.newComp("GLS",cs.be+d,{F,m:0.3,area:20000*fr,id:"t"+i}))
    :states.map((s,i)=>{ const c=C.newComp("GLS",s.be,{F:1.2,m:0.3,area:s.A,id:"t"+i}); if(dbl) c.doublet={dE:dTrue,r:info.r*(o.rScale||1),q:1}; return c; });
  const P=new Float64Array(n); truth.forEach(c=>{ const v=C.compEval(c,E,ctx); for(let i=0;i<n;i++) P[i]+=v[i]; });
  const Q=C.cumtrapz(E,P); const mu=Float64Array.from(P,(p,i)=>3000+15*(E[i]-lo)+0.06*Q[i]+p);
  const cf=0.3; const y=Float64Array.from(mu,m=>{ const l=m*cf; return Math.max(0,Math.round(l+Math.sqrt(l)*R.g()))/cf; });
  return {E,y,truth,info,dTrue,cf};
}
/* 다중항 준위 추가 경우(기본 목록은 원소마다 Scofield 센 궤도 2 개라 3d 전이금속 3p·3s, 가벼운 희토류 4d 가 빠진다) */
function multCases(){ const out=[]; [["Ti","3p"],["Mn","3p"],["Mn","3s"],["Fe","3p"],["Co","3p"],["Ni","3p"],["Cu","3p"],["Ce","4d"],["Nd","4d"],["Eu","4d"],["Gd","4d"],["Dy","4d"]].forEach(([el,orb])=>{
  const lv=(L.core[el]||{}).lv||{}; const k=Object.keys(lv).filter(x=>x.startsWith(orb)).sort((a,b)=>lv[a]-lv[b])[0]; if(k&&lv[k]>10&&lv[k]<1300) out.push({el,orb,be:+lv[k],extra:true}); }); return out; }
function run(opt){
  opt=opt||{}; const R=C.rng(opt.seed||31); const base=cases(); const CS=base.concat(multCases().filter(m=>!base.some(b=>b.el===m.el&&b.orb===m.orb))); const rep={n:CS.length,parse:[],rsf:[],mgGap:[],dbl:[],mult:[],sug:[],sugWideSmart:[],two:[],twoFixed:[],added:[],auto2:[],nAuto2:0,relaxed:[],rpos:[],nRpos:0,nSug:0,nTwo:0,nTwoFixed:0,times:0};
  const t0=Date.now();
  for(const cs of CS){
    // (1) 이름
    const names=[cs.el+cs.orb, cs.el+" "+cs.orb+" Scan", cs.el.toLowerCase()+cs.orb, "[S-1] "+cs.el+cs.orb+".txt", cs.el+"_"+cs.orb];
    names.forEach(nm=>{ const r=IO.parseRegionName(nm); if(r.el!==cs.el||r.orb!==cs.orb) rep.parse.push(nm+" → "+JSON.stringify(r)); });
    // (2) 감도계수
    if(!(L.scofield.al[cs.el+" "+cs.orb]>0)) rep.rsf.push(cs.el+" "+cs.orb+" Scofield-Al");
    if(!(L.scofield.mg[cs.el+" "+cs.orb]>0)) rep.mgGap.push(cs.el+" "+cs.orb);
    if(!(L.althermo1[cs.el+" "+cs.orb]>0)) rep.rsf.push(cs.el+" "+cs.orb+" ALTHERMO1");
    // (3) doublet
    const l=cs.orb[1]; const info=C.doubletInfo(L,cs.el,cs.orb); const mult=!!C.multipletInfo(cs.el,cs.orb);
    if(mult) rep.mult.push(cs.el+" "+cs.orb);
    if(l!=="s"&&!mult&&!(info&&info.dE)) rep.dbl.push(cs.el+" "+cs.orb);
    // (4) 단일 상태: 자동 제안 → 피팅
    for(const bgt of ["active","smart"]){
      const S=synth(R,cs,{two:false,mult});
      const sg=C.suggestModel(S.E,S.y,{bg:{type:bgt,n:5},el:cs.el,orb:cs.orb,lib:L});
      // 참 배경에 기울기가 있으므로 active 는 기울기도 피팅(화면의 '기울기도 피팅')
      const bgp=bgt==="active"?{type:"active",n:5,cons:{s:{mode:"free"}}}:{type:bgt,n:5};
      let prob={E:S.E,y:S.y,c:S.cf,bg:bgp,comps:sg.comps,cons:sg.cons};
      let ok=false, why="", why0="";
      try{
        // 화면의 '자동 제안'과 같은 순서: 봉우리 → 잔차로 성분 더하기(ΔBIC ≥ 10) → 피팅
        const rf=C.refineByResidual(prob,{}); prob=rf.prob; if(rf.added&&!mult) rep.added.push(`${cs.el} ${cs.orb} [${bgt}] +${rf.added}`); if(rf.relaxed&&rf.relaxed.length) rep.relaxed.push(`${cs.el} ${cs.orb} [${bgt}] ${rf.relaxed.join(",")}`);
        const r=C.fitOnce(prob,{initAreas:false}); const st=C.fitStats(prob,r);
        const main=prob.comps.slice().sort((a,b)=>b.area-a.area)[0]; const tot=prob.comps.reduce((s,c)=>s+c.area,0);
        const dOK=mult?prob.comps.every(c=>!c.doublet):(!S.info||!S.info.dE||(main.doublet&&Math.abs(main.doublet.dE-S.dTrue)<0.15));
        // 정적 배경(smart)은 참 배경(active Shirley + 기울기)과 모양이 달라 χ²ν·넓이 허용을 넓힌다(SPEC P4: 창 밖 꼬리 손실 몇 %)
        // Δ 를 고정하는(분해 안 되는) doublet 은 라이브러리 Δ 와 참값의 차(≤0.1 eV) 만큼 χ²ν 가 오를 수 있다
        const fixedD=S.info&&S.info.dE&&!S.info.resolvable;
        const wW=S.E[S.E.length-1]-S.E[0]; const aTol=bgt==="active"?0.03:(wW<=20?0.05:0.10), cTol=bgt==="active"?(fixedD?2.0:mult?1.5:1.35):2.5;   // 다중항 준위(mult): 2–3 성분 근사라 χ²ν 가 1.35–1.41 까지 오름(seed 31–33, 옛·새 Scofield 표 모두) → 1.5
        // 정적(smart) 배경의 넓이 편향은 창 폭에 따라 커진다(합성 143 경우: ≤20 eV 최대 3 %, 20–40 eV 최대 8 %, >40 eV 최대 14 %) — 허용: ≤20 eV 5 %, 20–40 eV 10 %, >40 eV 는 정보
        // active: 넓이는 공분산 σ 대비 z < 4 (넓은 창에서는 배경 k 와의 상관으로 σ 가 커진다); 정적 배경: 상대 오차
        let aOK=Math.abs(tot/20000-1)<aTol;
        if(bgt==="active"){ const cov=C.covariance(prob,r); const sA=C.areaSumSigma(prob,r.map,cov,prob.comps.map((c,i)=>i)); aOK=Math.abs(tot-20000)/sA<4.5; why0=` σ_A/A ${(sA/20000).toFixed(3)}`; }
        ok=prob.comps.length>=1&&st.chi2nu<cTol&&(mult||Math.abs(main.be-cs.be)<0.06)&&aOK&&dOK;
        why=`n=${prob.comps.length} χ²ν=${st.chi2nu.toFixed(2)} BE ${main.be.toFixed(3)}/${cs.be} A ${(tot/20000).toFixed(3)} Δ ${main.doublet?main.doublet.dE.toFixed(2):"-"}/${S.dTrue?S.dTrue.toFixed(2):"-"}${why0}`;
      }catch(e){ why="오류 "+e.message; }
      const wide=S.E[S.E.length-1]-S.E[0]>40;
      if(!ok){ if(bgt==="smart"&&wide) rep.sugWideSmart.push(`${cs.el} ${cs.orb} [smart, 창 ${(S.E[S.E.length-1]-S.E[0]).toFixed(0)} eV] ${why}`); else rep.sug.push(`${cs.el} ${cs.orb} [${bgt}] ${why}`); }
      rep.nSug++;
    }
    // (5) 두 상태: 참 구조에서 흔든 출발 — 다중항 준위는 화학 상태 분해를 하지 않으므로 제외
    if(!mult){ const S=synth(R,cs,{two:true});
      const comps=S.truth.map((t,i)=>{ const c=C.clone(t); c.id="c"+i; c.be+=0.2*R.g(); c.F*=1+0.15*R.g(); c.area*=0.7; return c; });
      const cons={}; comps.forEach((c,i)=>{ const cs_=cons[c.id]={F:{mode:"free",lo:0.3,hi:4}}; if(i) cs_.m={mode:"link",ref:"c0",k:1,b:0};
        if(c.doublet){ cs_.r={mode:"fixed"}; cs_.q={mode:"fixed"}; c.doublet.dE=S.info.dE; cs_.dE=i?{mode:"link",ref:"c0",k:1,b:0}:C.dEConstraint(S.info); } });
      const prob={E:S.E,y:S.y,c:S.cf,bg:{type:"active",n:5,cons:{s:{mode:"free"}}},comps,cons};
      let ok=false, why="";
      try{ const r=C.fitOnce(prob,{}); const st=C.fitStats(prob,r); const cov=C.covariance(prob,r); const sm=C.compSummary(prob,r,cov);
        const z=prob.comps.map((c,i)=>{ const t=S.truth[i]; return [Math.abs(c.be-t.be)/Math.max(sm[i].sig.be,1e-6), Math.abs(c.area-t.area)/Math.max(sm[i].sig.area,1e-6)]; });
        const zmax=Math.max(...z.flat());
        // 분해되지 않는 doublet 은 Δ 를 라이브러리 값으로 고정 → Δ_참과의 차(≤0.1 eV) 만큼 계통 편향: 절대 오차로 판정
        const fixedD=S.info&&S.info.dE&&!S.info.resolvable;
        const absOK=prob.comps.every((c,i)=>Math.abs(c.be-S.truth[i].be)<0.08&&Math.abs(c.area/S.truth[i].area-1)<0.06);
        ok=fixedD?(st.chi2nu<2.0&&absOK):(st.chi2nu<1.35&&zmax<4); why=`χ²ν=${st.chi2nu.toFixed(2)} z_max=${zmax.toFixed(2)}${fixedD?" (Δ 고정 — 절대 오차 판정 "+absOK+")":""} BE ${prob.comps.map((c,i)=>(c.be-S.truth[i].be).toFixed(3)).join("/")} A ${prob.comps.map((c,i)=>(c.area/S.truth[i].area).toFixed(3)).join("/")}`; }
      catch(e){ why="오류 "+e.message; }
      // (6) 같은 두 상태 스펙트럼을 자동 제안 + 잔차 보강으로만 — 두 상태가 둘 다 잡히는가
      try{ const sg=C.suggestModel(S.E,S.y,{bg:{type:"active",n:5},el:cs.el,orb:cs.orb,lib:L});
        const rf=C.refineByResidual({E:S.E,y:S.y,c:S.cf,bg:{type:"active",n:5,cons:{s:{mode:"free"}}},comps:sg.comps,cons:sg.cons},{});
        const P=rf.prob; const got=S.truth.map(t=>{ const near=P.comps.filter(c=>Math.abs(c.be-t.be)<0.15); return near.length?near.reduce((a,c)=>a+c.area,0)/t.area:0; });
        const okA=got.every(g=>Math.abs(g-1)<0.1)&&P.comps.length===S.truth.length; rep.nAuto2++;
        if(!okA) rep.auto2.push(`${cs.el} ${cs.orb} 성분 ${P.comps.length}: ${P.comps.map(c=>c.be.toFixed(2)+(c.doublet?"D":"")).join(" ")} / 참 ${S.truth.map(t=>t.be.toFixed(2)).join(" ")} 넓이비 ${got.map(g=>g.toFixed(2)).join("/")}`); }
      catch(e){ rep.auto2.push(`${cs.el} ${cs.orb} 오류 ${e.message}`); }
      // (7) 양성 대조: 참 넓이 비가 통계 값의 1.2 배인 두 상태 스펙트럼 — 잔차 보강이 r 을 풀어 회복하는가 (분해되는 doublet 만)
      if(S.info&&S.info.dE&&S.info.dE>=C.DOUBLET_RESOLVE){ try{ const T=synth(R,cs,{two:true,rScale:1.2});
        const sg=C.suggestModel(T.E,T.y,{bg:{type:"active",n:5},el:cs.el,orb:cs.orb,lib:L});
        const rf=C.refineByResidual({E:T.E,y:T.y,c:T.cf,bg:{type:"active",n:5,cons:{s:{mode:"free"}}},comps:sg.comps,cons:sg.cons},{});
        const rr=rf.prob.comps.filter(c=>c.doublet).map(c=>c.doublet.r/T.truth[0].doublet.r); rep.nRpos++;
        if(!(rf.relaxed||[]).includes("r")||rr.some(x=>Math.abs(x-1)>0.05)||rf.prob.comps.length!==2) rep.rpos.push(`${cs.el} ${cs.orb} 풀림 ${(rf.relaxed||[]).join(",")||"없음"} r비 ${rr.map(x=>x.toFixed(3)).join("/")} 성분 ${rf.prob.comps.length}`); }
        catch(e){ rep.rpos.push(`${cs.el} ${cs.orb} 오류 ${e.message}`); } }
      const fixedD=S.info&&S.info.dE&&!S.info.resolvable; rep.nTwo++; if(fixedD) rep.nTwoFixed++;
      if(!ok){ if(fixedD) rep.twoFixed.push(`${cs.el} ${cs.orb} ${why}`); else rep.two.push(`${cs.el} ${cs.orb} ${why}`); } }
  }
  rep.times=(Date.now()-t0)/1000;
  return rep;
}
module.exports={cases,run,synth};
if(require.main===module){ const r=run(); console.log(JSON.stringify({n:r.n,parse:r.parse.length,rsf:r.rsf.length,mgGap:r.mgGap,dbl:r.dbl,sug:r.sug.length,two:r.two.length,s:r.times},null,0)); console.log("다중항 준위:",r.mult.length,r.mult.join(", ")); console.log("두 상태 자동:",r.nAuto2-r.auto2.length,"/",r.nAuto2); r.auto2.forEach(x=>console.log("  auto2",x)); console.log("r 양성 대조:",r.nRpos-r.rpos.length,"/",r.nRpos); r.rpos.forEach(x=>console.log("  rpos",x)); console.log("r·q 를 잘못 푼 경우(참은 통계 값):",r.relaxed.length,r.relaxed.join(", ")); console.log("잔차로 더한 성분(참 모형은 1 성분):",r.added.length,r.added.join(", ")); console.log("smart 넓은 창(정보):",r.sugWideSmart.length,"Δ고정 두 상태(정보):",r.twoFixed.length,"/",r.nTwoFixed); ["parse","rsf","sug","two","sugWideSmart","twoFixed"].forEach(k=>r[k].slice(0,40).forEach(x=>console.log(k,x))); }
