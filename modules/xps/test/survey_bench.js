/* 합성 survey 로 원소 판별 성능 측정 (참값을 아는 데이터) — node test/survey_bench.js [N]   (trigPow 0·0.5·1·2 비교)
   원소 3–6 개(무작위, 흔한 원소 풀), 원소별 화학 이동 δ∈[−1.5,+3.5], 선 높이 ∝ 농도·Scofield·KE^−0.4,
   FWHM 1.6–2.4 eV, 각 선 뒤 비탄성 계단, Auger(Avantage rsf), Poisson 잡음. 1 eV 간격 0–1350 eV. */
const C=require("../src/xps_core.js"), L=require("../src/xps_lib.js");
const POOL=["C","O","N","F","Na","Mg","Al","Si","P","S","Cl","K","Ca","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Zr","Nb","Mo","Ag","In","Sn","Sr","Y","Ba","La","Ce","Ta","W","Pt","Au","Pb","Bi"];
function synth(seed){
  const R=C.rng(seed); const hv=1486.6; const n=1361; const E=Float64Array.from({length:n},(_,i)=>-10+i);
  const k=3+Math.floor(R.u()*4); const els=[]; while(els.length<k){ const e=POOL[Math.floor(R.u()*POOL.length)]; if(!els.includes(e)) els.push(e); }
  if(!els.includes("C")&&R.u()<0.8) els.push("C"); if(!els.includes("O")&&R.u()<0.8) els.push("O");
  const mu=new Float64Array(n).fill(0); const base=Float64Array.from(E,e=>3000+40*Math.max(0,e)); const truth=[];
  const sc=L.scofield.al, idp=L.idpos;
  els.forEach(el=>{ const conc=Math.pow(10,-1.3+1.3*R.u()); const d=-1.5+5*R.u(); const rec=L.core[el]; let hmain=0, truthPos=0;
    Object.entries(rec.lv).forEach(([sub,bk])=>{ const s=sc[el+" "+sub]; if(s==null) return; const pos0=idp[el+" "+sub]!=null?idp[el+" "+sub]:(rec.vac?null:bk); if(pos0==null||pos0>hv-30) return;
      const pos=pos0+d, KE=hv-pos, h=2.0e5*conc*s*Math.pow(KE/1000,-0.4), w=1.6+0.8*R.u(), sg=w/2.3548;
      if(h>hmain){ hmain=h; truthPos=pos; }
      for(let i=0;i<n;i++){ const x=E[i]-pos; mu[i]+=h*Math.exp(-0.5*x*x/(sg*sg)); if(x>0) mu[i]+=0.35*h*(1-Math.exp(-x/15)); } });
    (L.auger||[]).filter(a=>a.el===el).forEach(a=>{ const pos=hv-a.ke; if(pos<0||pos>hv-30) return; const h=2.0e5*conc*a.rsf*0.5, sg=3/2.3548; for(let i=0;i<n;i++){ const x=E[i]-pos; mu[i]+=h*Math.exp(-0.5*x*x/(sg*sg)); if(x>0) mu[i]+=0.5*h*(1-Math.exp(-x/25)); } });
    truth.push({el,conc,d,hmain,pos:truthPos});
  });
  const cf=0.02; const y=Float64Array.from(mu,(m,i)=>{ const l=(m+base[i])*cf; return Math.max(0,l+Math.sqrt(l)*R.g())/cf; });
  // 검출 가능성: 주선 높이 / 그 자리 Poisson σ (1 eV 점 하나 기준)
  truth.forEach(t=>{ const i=Math.max(0,Math.min(n-1,Math.round(t.pos+10))); t.snr=t.hmain/Math.sqrt((mu[i]+base[i])/cf); });
  return {E,y,truth,hv};
}
function run(N,mode,opt){
  let tp=0,fp=0,fn=0, det=0; const fpEl={}, fnEl={};
  for(let s=1;s<=N;s++){ const S=synth(1000+s); const pk=C.findPeaks(S.E,S.y,{width:1.5,snr:4,merge:1.5});
    const c=C.identify(pk,L,{hv:S.hv,tol:2,E:S.E,y:S.y,width:1.5,...(opt||{})}); const acc=c.filter(C.acceptDefault).map(x=>x.el);
    // 검출 가능한 참 원소: 주선 높이가 잡음의 8 배 이상(대략 conc > 0.06)
    const T=S.truth.filter(t=>t.snr>15).map(t=>t.el); const Tall=S.truth.map(t=>t.el); det+=T.length;
    T.forEach(e=>{ if(acc.includes(e)) tp++; else { fn++; fnEl[e]=(fnEl[e]||0)+1; } });
    acc.forEach(e=>{ if(!Tall.includes(e)){ fp++; fpEl[e]=(fpEl[e]||0)+1; } });
  }
  return {N,mode,recall:tp/det,precisionLike:tp/(tp+fp),fpPerSpectrum:fp/N,fn,fp,top_fp:Object.entries(fpEl).sort((a,b)=>b[1]-a[1]).slice(0,8),top_fn:Object.entries(fnEl).sort((a,b)=>b[1]-a[1]).slice(0,8)};
}
module.exports={synth,run};
if(require.main===module){ const N=+process.argv[2]||40; for(const tp of [0,0.5,1,2]){ const t0=Date.now(); console.log("trigPow",tp,JSON.stringify(run(N,"greedy",{trigPow:tp})),((Date.now()-t0)/1000).toFixed(1)+"s"); } }
