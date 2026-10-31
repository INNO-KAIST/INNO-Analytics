/* ============================================================
   INNO Analytics — XPS fit guidance (v1.9.3): fit check (diagnosis), region recipes, staged fit
   DOM-free: browser main thread, Web Worker (load after xps_core.js) and node.
   Definitions, thresholds, citations and validation: modules/xps/SPEC.md §17 (if they disagree, fix SPEC first).
   General rules only — nothing here is tuned to one data set; the Ce 3d recipe is literature data (SPEC §17.3).
   ============================================================ */
(function(root){
"use strict";
function getCore(){ if(root&&root.XPS_CORE) return root.XPS_CORE; if(typeof require==="function"){ try{ return require("./xps_core.js"); }catch(e){} } return null; }
let _C=null; const CORE=()=>_C||(_C=getCore());

/* ------------------------------------------------------------
   0) Thresholds (SPEC §17.1) — conservative: a finding must be clear in the data, not a matter of taste
   ------------------------------------------------------------ */
const GUIDE_TH={
  boundRel:1e-3,     // parameter within 1e-3·(hi − lo) of a bound (same rule as the stall-rule 'bound' flag, SPEC §3.3b)
  corr:0.95,         // |r| between two component parameters
  corrBg:0.995,      // |r| when one parameter is the in-fit background (c0·k·s) — area↔k is high by construction on wide windows
  corrMax:3,         // report at most this many pairs
  negShare:0.01,     // area < 1 % of the summed component area → negligible
  relSig:1.0,        // σ(area)/area > 100 % → ill-determined
  merge:0.5,         // |ΔBE| < 0.5·mean FWHM between two free main lines → merging
  edgeComp:0.05,     // fitted component sum at a window end > 5 % of its maximum → truncated
  edgeRise:0.10,     // smoothed data rise towards the outside within edgeLen > 10 % of the net peak height → a peak continues outside
  edgeLen:1.0,       // eV
  bgAbove:0.20,      // background above the smoothed data by more than max(1 σ noise, bgAboveAmp · net peak height) on > 20 % of the window
  bgAboveAmp:0.03,
  blobZ:5,           // residual blob |z| ≥ 5 (noise σ and AR(1) ρ estimated from the data's second differences, independent of the model)
  blobHalf:0.6,      // eV, half width of the residual blob window
  blobAmp:0.03,      // and the blob's mean residual ≥ 3 % of the net peak height (a static-background endpoint error of ~1 % is not a missing line)
};
const SEV_ORDER={bad:0,warn:1,info:2};

/* ------------------------------------------------------------
   1) Complex regions (multiplet / satellite-rich) — element·orbital rule, independent of data
   ------------------------------------------------------------ */
const SAT_TM2P=["Cr","Mn","Fe","Co","Ni","Cu"], SAT_RE3D=["Ce","Pr","Nd"];
function complexInfo(el,orb){
  const C=CORE(); if(!el||!orb) return null; const rc=recipesFor(el,orb).map(r=>r.id); const m=C.multipletInfo(el,orb);
  if(m) return {kind:"multiplet",why:`${el} ${orb}: exchange (multiplet) splitting with the open ${m.kind} shell is comparable to or larger than the spin–orbit splitting — the envelope is many lines, not one doublet per chemical state`,closed:m.closed||null,recipes:rc};
  if(orb==="2p"&&SAT_TM2P.includes(el)) return {kind:"satellite",why:`${el} 2p: multiplet splitting and shake-up satellites of the open 3d shell (closed-shell states such as Cu⁺/Cu⁰ excepted) — one doublet per chemical state does not describe the envelope`,closed:el==="Cu"?"Cu⁺·Cu⁰ (3d¹⁰)":null,recipes:rc};
  if(orb==="3d"&&SAT_RE3D.includes(el)) return {kind:"satellite",why:`${el} 3d: 4f final-state (screening) satellites — several peaks per oxidation state, overlapping across the two spin–orbit halves`,closed:null,recipes:rc};
  if(rc.length) return {kind:"recipe",why:`${el} ${orb}: literature starting model available`,closed:null,recipes:rc};
  return null;
}

/* ------------------------------------------------------------
   2) Region recipes (SPEC §17.3) — data, not code. Positions are charge-corrected BE (C 1s C–C/C–H 284.8 eV).
      line = one spin–orbit pair (main = 3d5/2 'v' line, partner = 3d3/2 'u' line at +ΔE, area ratio r, equal FWHM).
   ------------------------------------------------------------ */
const RECIPES=[{
  id:"ce3d", el:"Ce", orb:"3d", name:"Ce 3d (Burroughs v/u)",
  window:[875,925], bg:{type:"active",n:5,c0:null,k:null,s:0,cons:{s:{mode:"free"}}}, shape:"GLS", m:{v:0.3,lo:0,hi:1},
  dE:{v:18.6,lo:18.3,hi:18.9}, dEs:{"Ce3+":{v:18.2,lo:17.9,hi:18.5}}, r:2/3, beTol:0.6, align:3,
  states:{"Ce4+":"Ce⁴⁺","Ce3+":"Ce³⁺"},
  lines:[
    {key:"v0", name:"v0/u0", state:"Ce3+", be:880.5, F:2.4, Flo:1.2, Fhi:3.5},
    {key:"v",  name:"v/u",   state:"Ce4+", be:882.4, F:2.3, Flo:1.2, Fhi:3.5},
    {key:"v1", name:"v′/u′", state:"Ce3+", be:885.0, F:3.1, Flo:1.5, Fhi:4.5},
    {key:"v2", name:"v″/u″", state:"Ce4+", be:888.9, F:5.0, Flo:2.5, Fhi:6.5},
    {key:"v3", name:"v‴/u‴", state:"Ce4+", be:898.2, F:2.3, Flo:1.2, Fhi:3.5}],
  options:[
    {id:"ce4", label:"CeO₂ only (6)", keys:["v","v2","v3"]},
    {id:"ce3", label:"Ce³⁺ only (4)", keys:["v0","v1"]},
    {id:"mix", label:"Mixed (10)", keys:["v0","v","v1","v2","v3"]}],
  derived:[
    {id:"ce3", label:"Ce³⁺", state:"Ce3+", tip:"Ce³⁺ fraction = (v0 + v′ + u0 + u′) / (all ten lines), areas from the fit (each pair = v + u). ± = covariance σ (statistics only). The fraction depends on the background and on the line shapes; beam reduction under X-rays raises it during the measurement (Isaacs et al. 2023; Morgan 2023)."},
    {id:"u3", label:"u‴", key:"v3", part:"partner", tip:"u‴ fraction = u‴ area / (all ten lines). u‴ (≈ 916.7 eV) has no Ce³⁺ counterpart, so it falls as Ce⁴⁺ is reduced; for fully oxidised CeO₂ Shyu et al. (J. Catal. 114, 23, 1988, Fig. 2; nonlinear Shirley) found ≈13–14 % of the Ce 3d region (read from the graph; ±10 % relative), but the value depends on background and fit model, so compare it only between spectra treated the same way."}],
  cite:"Burroughs et al., J. Chem. Soc. Dalton Trans. 1976, 1686 (v/u notation); Romeo et al., Surf. Interface Anal. 20, 508 (1993); Bêche et al., Surf. Interface Anal. 40, 264 (2008); Morgan, Surf. Interface Anal. 55 (2023), doi:10.1002/sia.7254; Isaacs et al., Appl. Surf. Sci. Adv. 18, 100469 (2023); Paparazzo, Mater. Res. Bull. 46, 323 (2011); Shyu et al., J. Catal. 114, 23 (1988) (u‴ fraction)",
  pitfalls:"Five spin–orbit pairs overlap (v‴ ≈ 898 eV lies under u0 ≈ 899 eV, u ≈ 901 eV under u′), so single-line areas are strongly correlated — report the Ce³⁺ fraction, not single peaks. Symmetric lines cannot follow the asymmetric Ce⁴⁺ peaks; the leftover intensity near 885–888 eV is then taken by v′ (Ce³⁺), so 10-line fits tend to overestimate Ce³⁺ — cross-check with the u‴ fraction, which has no Ce³⁺ counterpart. The Ce³⁺ fraction depends on the background (here in-fit Shirley + slope over 875–925 eV; a static Shirley on a sloping baseline biased it by +5 % on synthetic spectra), the u‴ fraction hardly does. X-rays and UHV reduce CeO₂ during acquisition (measure Ce 3d first, low flux). Positions are for C 1s = 284.8 eV and may sit off by charging; u‴ ≈ 916.7 eV is an internal reference. Lines start within ±0.6 eV; ΔE free per state (Ce⁴⁺ 18.6, 18.3–18.9 eV; Ce³⁺ 18.2, 17.9–18.5 eV; Bêche 2008), u/v area = 2/3 (Bêche: I3d5/2/I3d3/2 = 1.5), equal widths within a pair."
}];
/* spin–orbit Δ per oxidation state (Bêche et al. 2008: ≈18.6 eV for Ce⁴⁺/CeO₂, ≈18.2 eV for Ce³⁺/CePO₄) */
function dEof(R,state){ return (R.dEs&&R.dEs[state])||R.dE; }
function recipesFor(el,orb){ return RECIPES.filter(r=>r.el===el&&r.orb===orb); }
function recipeById(id){ return RECIPES.find(r=>r.id===id)||null; }
/* Build a recipe model on raw (uncorrected) BE: shift = charge correction Δ (BE_corr = BE_raw + Δ).
   o = {shift, be, y (full spectrum, raw BE ascending), opt (option id), align (true: search ±recipe.align eV for the best NNLS match, SPEC §17.3)}
   → {lo, hi, bg, comps, cons, recipe:{id,opt,align,...}} */
function recipeModel(id,o){
  const C=CORE(); const R=recipeById(id); if(!R) throw new Error("unknown recipe "+id);
  o=o||{}; const sh=o.shift||0; const op=R.options.find(x=>x.id===o.opt)||R.options[R.options.length-1];
  const lines=R.lines.filter(l=>op.keys.includes(l.key));
  let lo=R.window[0]-sh, hi=R.window[1]-sh; const clipped=[];
  if(o.be&&o.be.length){ const a=o.be[0], b=o.be[o.be.length-1]; if(lo<a){ clipped.push("low"); lo=a; } if(hi>b){ clipped.push("high"); hi=b; } }
  const mk=ds=>lines.map((l,i)=>{ const id="r"+R.id+"_"+l.key+"_"+Math.random().toString(36).slice(2,6);
    const c=C.newComp(R.shape,l.be-sh+ds,{F:l.F,m:R.m.v,el:R.el,line:R.orb,id,name:l.name});
    c.doublet={dE:dEof(R,l.state).v,r:R.r,q:1}; c.rk=l.key; c.rstate=l.state; c.recipe=R.id; c.area=1; return c; });
  let ds=0, score=null, areas=null;
  if(o.be&&o.y&&o.align!==false){
    const reg=C.sliceRegion(o.be,o.y,lo,hi); if(reg&&reg.E.length>=8){
      const E=Float64Array.from(reg.E), y=Float64Array.from(reg.y), n=E.length, ctx={emin:E[0],emax:E[n-1]};
      const B=C.bgStatic(E,y,{type:"smart",n:5}); const z=Float64Array.from(y,(v,i)=>v-B[i]);
      const ys=C.smoothQ(E,y,C.smoothHW(E)); const w=Float64Array.from(ys,v=>1/Math.sqrt(Math.max(v,1e-9)));
      const trial=d=>{ const cs=mk(d); const cols=cs.map(c=>Float64Array.from(C.compEval(c,E,ctx),(v,i)=>v*w[i])); const rhs=Float64Array.from(z,(v,i)=>v*w[i]);
        const a=C.nnls(cols,rhs); let ss=0; for(let i=0;i<n;i++){ let m=0; for(let j=0;j<cols.length;j++) m+=a[j]*cols[j][i]; const r=rhs[i]-m; ss+=r*r; } return {ss,a}; };
      let best=null; const st=0.05, K=Math.round((R.align||0)/st);
      for(let k=-K;k<=K;k++){ const d=k*st, t=trial(d); if(!best||t.ss<best.ss) best={d,ss:t.ss,a:t.a}; }
      ds=best.d; score=best.ss; areas=best.a; } }
  const comps=mk(ds); if(areas) comps.forEach((c,i)=>{ c.area=Math.max(areas[i],1e-6); });
  const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); const cons={};
  const ancS={}; comps.forEach(c=>{ const a=ancS[c.rstate]; if(!a||c.area>a.area) ancS[c.rstate]=c; });
  comps.forEach((c,i)=>{ const l=lines[i]; const q=cons[c.id]={};
    q.be={mode:"free",lo:+(c.be-R.beTol).toFixed(4),hi:+(c.be+R.beTol).toFixed(4)};
    q.F={mode:"free",lo:l.Flo,hi:l.Fhi};
    q.r={mode:"fixed"}; q.q={mode:"fixed"};
    const dS=dEof(R,l.state), aS=ancS[c.rstate]; q.dE=c===aS?{mode:"free",lo:dS.lo,hi:dS.hi}:{mode:"link",ref:aS.id,k:1,b:0};
    q.m=c===anc?{mode:"free",lo:R.m.lo,hi:R.m.hi}:{mode:"link",ref:anc.id,k:1,b:0}; });
  return {lo,hi,bg:C.clone(R.bg),comps,cons,recipe:{id:R.id,opt:op.id,label:op.label,align:ds,score,clipped,name:R.name}};
}
/* Derived recipe output (SPEC §17.3): fractions of summed areas with covariance σ (delta method on the free area parameters).
   fin = finalize() output; returns [{id,label,v,sig,tip}] in % */
function recipeDerived(fin,meta){
  const R=recipeById(meta&&meta.id); if(!R||!fin||!fin.comps) return [];
  const cs=fin.comps.filter(c=>c.recipe===R.id); if(!cs.length) return [];
  const T=cs.reduce((s,c)=>s+Math.max(0,c.area),0); if(!(T>0)) return [];
  const names=(fin.cov&&fin.cov.names)||[], np=(fin.cov&&fin.cov.np)||names.length, Cm=(fin.cov&&fin.cov.C)||[];
  const aIdx={}; names.forEach((q,j)=>{ if(q.k==="area") aIdx[q.id]=j; });
  return R.derived.map(d=>{
    const part=c=>{ if(d.state) return c.rstate===d.state?1:0; if(d.key) return c.rk===d.key?(d.part==="partner"?(c.doublet?c.doublet.r/(1+c.doublet.r):0):(c.doublet?1/(1+c.doublet.r):1)):0; return 0; };
    const num=cs.reduce((s,c)=>s+part(c)*Math.max(0,c.area),0), f=num/T;
    const g=new Float64Array(np); cs.forEach(c=>{ const j=aIdx[c.id]; if(j==null) return; g[j]=(part(c)-f)/T; });
    let s=0; const ok=true; for(let a=0;a<np;a++) if(g[a]) for(let b=0;b<np;b++) if(g[b]){ const v=Cm[a*np+b]; if(isFinite(v)) s+=g[a]*v*g[b]; }   // areas at a bound count as fixed (SPEC §3.4)
    const present=cs.some(c=>part(c)>0);
    return {id:d.id,label:d.label,v:present?100*f:null,sig:present&&ok?100*Math.sqrt(Math.max(0,s)):null,tip:d.tip};
  });
}

/* ------------------------------------------------------------
   3) Noise statistics from the data alone (SPEC §17.1 R8). Second differences d_i = y_i − (y_{i−1} + y_{i+1})/2 remove the smooth
      signal (peaks are many channels wide). Noise model MA(1) — lag-1 correlation ρ, none beyond (measured on VGD flat regions:
      lag 1 ≈ 0.3–0.5, lag ≥ 2 ≈ 0; Poisson counting noise: ρ = 0). Then var(d) = σ²(1.5 − 2ρ), corr(d_i, d_{i+1}) = (1.75ρ − 1)/(1.5 − 2ρ)
      → ρ = (1.5c + 1)/(1.75 + 2c), clipped to [0, 0.5]; the variance of a sum of W noise points ≈ W·σ²(1 + 2ρ).
   ------------------------------------------------------------ */
function d2Corr(rho){ return (1.75*rho-1)/(1.5-2*rho); }
function noiseMA1(y,v){
  const n=y.length; if(n<8) return {rho:0,s2:1,lrv:1};
  const d=[]; for(let i=1;i<n-1;i++){ const sc=v?Math.sqrt(v[i]):1; d.push((y[i]-(y[i-1]+y[i+1])/2)/sc); }
  let a=0,b=0; for(let i=0;i<d.length;i++){ a+=d[i]*d[i]; if(i) b+=d[i]*d[i-1]; } const c=a>0?b/(a*(d.length-1)/d.length):-2/3;
  const rho=Math.max(0,Math.min(0.5,(1.5*c+1)/(1.75+2*c)));
  const q=d.map(x=>x*x).sort((x,y)=>x-y); const med=q[q.length>>1]/0.4549;   // robust variance of d: median of χ²₁ = 0.4549
  const s2=med/(1.5-2*rho); return {rho,s2,lrv:s2*(1+2*rho),c};
}
/* ------------------------------------------------------------
   4) Fit check (SPEC §17.1). ctx = {E, y, fin, el, orb, spec:{be0, be1} (raw spectrum extent), shift, recipe, dEref (literature Δ), anyRecipe}
      → [{id, sev:"bad"|"warn"|"info", label, tip, act|null, ids:[comp ids]}] sorted by severity
   ------------------------------------------------------------ */
function diagnose(ctx,th){
  const C=CORE(); th=Object.assign({},GUIDE_TH,th||{}); const out=[];
  const fin=ctx.fin; if(!fin||!fin.comps) return out;
  const E=ctx.E, y=ctx.y, n=E.length, sh=ctx.shift||0, step=n>1?(E[n-1]-E[0])/(n-1):0.05;
  const comps=fin.comps, byId={}; comps.forEach((c,i)=>byId[c.id]=i);
  const nm=id=>{ if(id==="bg") return "Background"; const c=comps[byId[id]]; return c?(c.name||("#"+(byId[id]+1))):"?"; };
  const lab=k=>(C.PDEF[k]||{lab:k}).lab;
  const fmt=(v,d)=>v==null||!isFinite(v)?"—":(+v).toFixed(d);
  const val=(id,k)=>{ if(id==="bg") return (fin.bg||{})[k]; const c=comps[byId[id]]; if(!c) return null; try{ return C.getP(c,k); }catch(e){ return null; } };
  const free=fin.free||[], cov=fin.cov||{}, np=cov.np||free.length, corr=cov.corr||[], atB=cov.atB||[];
  const summ=fin.summ||[]; const tot=comps.reduce((s,c)=>s+Math.max(0,c.area||0),0)||1;
  const fw=i=>(summ[i]&&summ[i].fwhm)||comps[i].F||1;
  const isFree=(id,k)=>free.some(f=>f.id===id&&f.k===k);
  const push=f=>{ if(!out.some(x=>x.key===f.key)) out.push(f); };
  const dropAct=id=>({type:"drop",id,label:"Remove "+nm(id)});
  const areaDead=new Set();
  // R3a area at 0 (lower bound) → handled as a component finding
  free.forEach((f,j)=>{ if(f.k!=="area"||f.ci<0) return; const p=val(f.id,"area"); if(p!=null&&(atB[j]||p<=Math.max(0,f.lo)+1e-6*tot)) areaDead.add(f.id); });
  // R1 parameters at / near a bound
  free.forEach((f,j)=>{
    if(f.k==="area") return; const p=val(f.id,f.k); if(p==null||!isFinite(p)) return;
    const L=isFinite(f.lo), H=isFinite(f.hi); if(!L&&!H) return;
    const dl=L?p-f.lo:Infinity, dh=H?f.hi-p:Infinity, ref=L&&H?f.hi-f.lo:Math.abs(p)+1;
    if(Math.min(dl,dh)>th.boundRel*ref) return;
    const side=dl<=dh?"lower":"upper", bv=side==="lower"?f.lo:f.hi, add=f.k==="be"?sh:0;
    const shapeP=["m","alpha","beta","T","tau","ET","TM","CT","mu","stMu","stG","tk","tj","ta","tb"].includes(f.k);
    let sev=shapeP||f.id==="bg"?"info":"warn", act=null, hint="";
    if(["F","wG","wL"].includes(f.k)){
      const ci=byId[f.id], others=comps.filter((c,i)=>i!==ci&&c[f.k]!=null&&isFree(c.id,f.k));
      const anc=others.reduce((a,c)=>!a||c.area>a.area?c:a,null);
      const rc=ctx.recipe?ctx.recipe.lines&&ctx.recipe.lines.find(l=>l.key===comps[ci].rk):null;
      if(rc&&f.k==="F"&&side==="upper"&&p>=rc.Fhi-1e-9){ act=null; hint=" The recipe width range is a typical span (not a literature limit) — a wider peak usually means a missing line, a wrong background or a charging/beam-damage broadening, not a wider component."; }
      else if(anc){ act={type:"link",id:f.id,k:f.k,ref:anc.id,label:`Link ${lab(f.k)} to ${nm(anc.id)}`}; hint=` Linking it to ${nm(anc.id)} (same ratio as now) removes one free width; if the line should be broader (e.g. a satellite), fix it at a literature value instead.`; }
      else if(side==="lower"){ act={type:"fix",id:f.id,k:f.k,v:p,label:`Fix ${lab(f.k)}`}; hint=" A width at its lower limit is narrower than the instrument resolution allows — fix it at the width of a well-resolved line in the same spectrum."; }
    } else if(f.k==="be"){
      const edge=Math.abs(bv-E[0])<=1.5*step||Math.abs(bv-E[n-1])<=1.5*step;
      if(edge){ act=dropAct(f.id); hint=" The line drifted to the window edge — it is not supported by the data inside the window. Remove it, or widen the window if the line is really outside."; }
      else { act={type:"widen",id:f.id,k:"be",by:0.3,label:"Widen range ±0.3 eV"}; hint=" Widen the range only if the literature allows a larger shift (charge-reference error is ±0.2–0.3 eV); otherwise the model is missing a line or the line belongs to another state."; }
    } else if(f.k==="dE"){
      const ref=ctx.recipe&&ctx.recipe.dE?dEof(ctx.recipe,(comps[byId[f.id]]||{}).rstate).v:ctx.dEref;
      if(ref!=null){ act={type:"fix",id:f.id,k:"dE",v:ref,label:`Fix Δ at ${fmt(ref,2)} eV`}; hint=` Fix the spin–orbit splitting at the literature value (${fmt(ref,2)} eV) unless the doublet is well resolved and isolated.`; }
    } else if(f.k==="r"||f.k==="q"){
      const ref=f.k==="q"?1:null; if(ref!=null) act={type:"fix",id:f.id,k:f.k,v:ref,label:`Fix ${lab(f.k)} at 1`};
      hint=" A doublet ratio at its limit usually compensates for an overlapping line; fix it at the statistical value (2j+1 ratio) and look at the residual.";
    } else if(f.id==="bg"){ hint=" The in-fit background parameter reached its limit: the window ends are probably not on background (widen the window or use a static Shirley/Smart background)."; }
    else hint=" A shape parameter at a physical limit is often fine (e.g. L frac. 0 = pure Gaussian); fix it there to remove a free parameter.";
    if(shapeP&&!act) act={type:"fix",id:f.id,k:f.k,v:p,label:`Fix ${lab(f.k)} at ${fmt(p,3)}`};
    push({key:"bound|"+f.id+"|"+f.k,id:"bound",sev,ids:f.id==="bg"?[]:[f.id],
      label:`${nm(f.id)} ${lab(f.k)} at ${side} bound (${fmt(bv+add,f.k==="be"?2:3)})`,
      tip:`${nm(f.id)} ${lab(f.k)} = ${fmt(p+add,4)} is within ${th.boundRel*100} % of the range of its ${side} bound ${fmt(bv+add,4)}. The value is set by the bound, not by the data, so it has no σ and the other parameters are biased.`+hint,act});
  });
  // R2 correlated pairs
  if(np&&corr.length===np*np){
    const pairs=[]; for(let a=0;a<np;a++) for(let b=a+1;b<np;b++){ if(atB[a]||atB[b]) continue; const r=corr[a*np+b]; if(!isFinite(r)) continue; const fa=free[a], fb=free[b]; if(!fa||!fb) continue;
      const lim=fa.id==="bg"||fb.id==="bg"?th.corrBg:th.corr; if(Math.abs(r)>lim) pairs.push({a,b,r}); }
    pairs.sort((x,y)=>Math.abs(y.r)-Math.abs(x.r));
    pairs.slice(0,th.corrMax).forEach(({a,b,r})=>{ const fa=free[a], fb=free[b]; let act=null, hint="", sevOv=null;
      if(fa.k===fb.k&&["F","wG","wL"].includes(fa.k)&&fa.id!==fb.id){ const ia=byId[fa.id], ib=byId[fb.id]; const big=comps[ia].area>=comps[ib].area?fa:fb, sm=big===fa?fb:fa;
        act={type:"link",id:sm.id,k:fa.k,ref:big.id,label:`Link ${lab(fa.k)}s`}; hint=" Linking the two widths (fixed ratio) usually resolves it."; }
      else if(fa.k==="be"&&fb.k==="be"){ const ia=byId[fa.id], ib=byId[fb.id]; const big=comps[ia].area>=comps[ib].area?fa:fb, sm=big===fa?fb:fa;
        act={type:"offset",id:sm.id,ref:big.id,label:"Fix separation"}; hint=" Fixing the separation at the current (or literature) value leaves one position free."; }
      else if((()=>{ const unres=id=>{ const c=comps[byId[id]]; return !!(c&&c.doublet&&c.doublet.dE<C.DOUBLET_RESOLVE); }; return (fa.k==="dE"&&unres(fa.id))||(fb.k==="dE"&&unres(fb.id))||(fa.id===fb.id&&unres(fa.id)); })()){
        const q=fa.k==="dE"?fa:fb.k==="dE"?fb:fa; const ref=ctx.recipe&&ctx.recipe.dE?dEof(ctx.recipe,(comps[byId[q.id]]||{}).rstate).v:ctx.dEref;
        if(ref!=null&&isFree(q.id,"dE")) act={type:"fix",id:q.id,k:"dE",v:ref,label:`Fix Δ at ${fmt(ref,2)} eV`}; hint=" The doublet is not resolved (Δ < 1.2 eV), so its two halves trade off against position and width; fix Δ at the literature value and keep the area ratio fixed."; sevOv="info"; }
      else if(fa.k==="area"&&fb.k==="area") hint=" The two areas trade off against each other — their sum is well determined, each alone is not. Report the sum, or constrain widths/positions from literature.";
      else hint=" One of the two parameters should be fixed or linked (literature value) — the data cannot separate them.";
      const inRec=fa.id!=="bg"&&fb.id!=="bg"&&comps[byId[fa.id]].recipe&&comps[byId[fa.id]].recipe===comps[byId[fb.id]].recipe;
      if(inRec){ act=null; hint=" Expected in this recipe (overlapping lines): report the derived fraction — its ± includes this correlation — not single-line areas."; }
      push({key:"corr|"+fa.id+fa.k+"|"+fb.id+fb.k,id:"corr",sev:inRec?"info":(sevOv||"warn"),ids:[fa.id,fb.id].filter(x=>x!=="bg"),
        label:`${nm(fa.id)} ${lab(fa.k)} ↔ ${nm(fb.id)} ${lab(fb.k)} correlated (r ${fmt(r,3)})`,
        tip:`Correlation coefficient ${fmt(r,4)} from the covariance matrix (|r| > ${fa.id==="bg"||fb.id==="bg"?th.corrBg:th.corr}). Changing one can be compensated by the other with almost no change in χ², so both σ are large and the split between them is arbitrary.`+hint,act}); });
  }
  // R3 negligible / ill-determined components
  comps.forEach((c,i)=>{ const share=Math.max(0,c.area||0)/tot, s=summ[i]&&summ[i].sig?summ[i].sig.area:null, rel=s!=null&&c.area>0?s/c.area:null;
    if(areaDead.has(c.id)) push({key:"dead|"+c.id,id:"dead",sev:"warn",ids:[c.id],label:`${nm(c.id)}: area went to zero`,tip:`The fit drove the area of ${nm(c.id)} to its lower bound 0 — the data do not need this component. Remove it and fit again (a zero-area line makes the covariance singular).`,act:dropAct(c.id)});
    else if(rel!=null&&rel>th.relSig) push({key:"ill|"+c.id,id:"ill",sev:"warn",ids:[c.id],label:`${nm(c.id)}: area ill-determined (±${fmt(100*rel,0)} %)`,tip:`σ(area)/area = ${fmt(100*rel,0)} % (> ${th.relSig*100} %): the data cannot say how large this component is. Typical causes: it overlaps a larger line with free width/position, or it is not there. Constrain it (literature position/width, link the width) or remove it.`,act:dropAct(c.id)});
    else if(share<th.negShare&&comps.length>1) push({key:"neg|"+c.id,id:"neg",sev:"info",ids:[c.id],label:`${nm(c.id)}: negligible (${fmt(100*share,1)} % of area)`,tip:`This component carries ${fmt(100*share,2)} % of the component area (< ${th.negShare*100} %). If it is not required physically, remove it — it adds free parameters without information.`+(c.recipe?" In a recipe model a negligible group means that state is absent; the simpler recipe option may be enough.":""),act:dropAct(c.id)});
  });
  // R4 merging main lines (both positions free and not tied to each other)
  const freeBe=comps.map(c=>isFree(c.id,"be"));
  for(let i=0;i<comps.length;i++) for(let j=i+1;j<comps.length;j++){ if(!freeBe[i]||!freeBe[j]) continue; if(areaDead.has(comps[i].id)||areaDead.has(comps[j].id)) continue;
    const d=Math.abs(comps[i].be-comps[j].be), w=(fw(i)+fw(j))/2; if(!(d<th.merge*w)) continue;
    const sm=comps[i].area<=comps[j].area?comps[i]:comps[j], bg_=sm===comps[i]?comps[j]:comps[i];
    push({key:"merge|"+comps[i].id+"|"+comps[j].id,id:"merge",sev:"warn",ids:[comps[i].id,comps[j].id],label:`${nm(comps[i].id)} and ${nm(comps[j].id)} merge (ΔBE ${fmt(d,2)} eV)`,
      tip:`The two lines are ${fmt(d,3)} eV apart, less than ${th.merge} × their mean FWHM (${fmt(w,2)} eV): they describe one peak. Remove the smaller one (${nm(sm.id)}), or fix its position from literature relative to ${nm(bg_.id)}.`,act:dropAct(sm.id)}); }
  // R5 stop reason
  if(fin.stop==="stalled"||fin.stop==="maxIter"||fin.conv===false){ const fl=((fin.stall&&fin.stall.flags)||[]).map(q=>`${nm(q.id)} ${lab(q.k)} (${q.why.join(", ")})`);
    push({key:"stop",id:"stop",sev:"warn",ids:[],label:fin.stop==="stalled"?"Fit stalled":"Fit did not converge",tip:(fin.stop==="stalled"?"LM stopped early: only heavily damped steps were possible (SPEC §3.3b).":"The iteration limit was reached before χ² settled.")+(fl.length?" Poorly determined: "+fl.join("; ")+".":"")+" With many overlapping lines, 'Staged fit' (heights → positions → all) reaches a minimum more reliably; also fix or link the flagged parameters.",act:comps.length>=3?{type:"staged",label:"Staged fit"}:null}); }
  // R6 window edges
  const P=fin.mu&&fin.B?fin.mu.map((m,i)=>m-fin.B[i]):null; const ys=C.smoothQ(E,y,C.smoothHW(E));
  let net=0; if(fin.B) for(let i=0;i<n;i++) net=Math.max(net,ys[i]-fin.B[i]); let pmax=0; if(P) for(const v of P) pmax=Math.max(pmax,v);
  const ne=Math.max(3,Math.round(0.3/Math.max(step,1e-6))), kL=Math.min(n-1,Math.max(1,Math.round(th.edgeLen/Math.max(step,1e-6))));
  const sideTest=(lowSide)=>{ const idx=lowSide?[...Array(ne).keys()]:[...Array(ne).keys()].map(k=>n-1-k); const pe=P?idx.reduce((s,i)=>s+P[i],0)/ne:0;
    const e0=lowSide?0:n-1, inner=lowSide?[...Array(kL+1).keys()]:[...Array(kL+1).keys()].map(k=>n-1-k); let mn=Infinity; inner.forEach(i=>{ if(ys[i]<mn) mn=ys[i]; });
    const rise=ys[e0]-mn; const byComp=pmax>0&&pe>th.edgeComp*pmax, byData=net>0&&rise>th.edgeRise*net; return {byComp,byData,pe:pmax>0?pe/pmax:0,rise:net>0?rise/net:0}; };
  [["low",true],["high",false]].forEach(([nmS,lowS])=>{ const t=sideTest(lowS); if(!t.byComp&&!t.byData) return;
    const wmax=comps.reduce((m,c,i)=>Math.max(m,fw(i)),1), by=Math.max(2,+(2*wmax).toFixed(1));
    const sp=ctx.spec||{}; const room=lowS?(sp.be0!=null?E[0]-sp.be0:0):(sp.be1!=null?sp.be1-E[n-1]:0);
    const act=room>step*2?{type:"window",side:nmS,by:Math.min(by,room),label:`Widen window (${nmS} BE) by ${fmt(Math.min(by,room),1)} eV`}:null;
    push({key:"edge|"+nmS,id:"edge",sev:"warn",ids:[],label:`Peak cut at the ${nmS}-BE window end`,
      tip:(t.byComp?`The fitted components still have ${fmt(100*t.pe,0)} % of their maximum at the ${nmS}-BE end (> ${th.edgeComp*100} %). `:"")+(t.byData?`The smoothed data rise by ${fmt(100*t.rise,0)} % of the net peak height towards the ${nmS}-BE end within ${th.edgeLen} eV (> ${th.edgeRise*100} %) — signal continues outside. `:"")+"The background endpoint then sits on a peak, which biases every area. Widen the window so both ends are on flat background."+(act?"":" The spectrum has no more data on this side — measure a wider range next time."),act}); });
  // R7 background above the data
  if(fin.B&&ctx.bgType!=="none"){ const nz=CORE().noiseEstimate(E,y); const kap=nz&&nz.kappa>0?nz.kappa:1; let above=0;
    for(let i=0;i<n;i++){ const s=Math.max(Math.sqrt(kap*Math.max(ys[i],1e-9)),th.bgAboveAmp*pmax); if(fin.B[i]>ys[i]+s) above++; }
    const fr=above/n; if(fr>th.bgAbove){ const sw=["shirley","linear","tougaard","spline"].includes(ctx.bgType);
      push({key:"bgabove",id:"bgabove",sev:"warn",ids:[],label:`Background above data on ${fmt(100*fr,0)} % of the window`,tip:`The background is more than max(1 σ noise, ${th.bgAboveAmp*100} % of the peak height) above the smoothed data on ${fmt(100*fr,0)} % of the window (> ${th.bgAbove*100} %): areas there are negative. Usual causes: a window end on a peak or on noise spikes, or a Shirley without a fixed point on a sloping baseline. Use Smart (clipped to the data) or the in-fit Shirley, or move the window ends onto flat background.`,act:sw?{type:"bg",bg:"smart",label:"Use Smart background"}:null}); } }
  // R8 structured residual (blob of the normalised residual; noise from second differences of the data — model-independent)
  if(fin.mu&&fin.v){ const v=fin.v, r=new Float64Array(n); for(let i=0;i<n;i++) r[i]=(y[i]-fin.mu[i])/Math.sqrt(v[i]);
    const nz=noiseMA1(y,v); const w=Math.max(2,Math.round(th.blobHalf/Math.max(step,1e-6)));
    let best=null; for(let i=w;i<n-w;i++){ let s=0, a=0; for(let j=i-w;j<=i+w;j++){ s+=r[j]; a+=y[j]-fin.mu[j]; } const z=s/Math.sqrt(nz.lrv*(2*w+1)), amp=a/(2*w+1);
      if(pmax>0&&Math.abs(amp)<th.blobAmp*pmax) continue; if(!best||Math.abs(z)>Math.abs(best.z)) best={i,z,amp}; }
    if(ctx.debug) ctx.debug.blob={z:best&&best.z,be:best&&E[best.i],rho:nz.rho,s2:nz.s2};
    if(best&&Math.abs(best.z)>=th.blobZ){ const be=E[best.i], pos=best.z>0;
      let act=null; if(pos){ const anc=comps.reduce((a,c)=>!a||c.area>a.area?c:a,null); act={type:"add",be,F:Math.min(3,Math.max(0.8,anc?fw(byId[anc.id]):1.2)),h:Math.max(0,y[best.i]-fin.mu[best.i]),label:`Add a line at ${fmt(be+sh,1)} eV`}; }
      push({key:"resid",id:"resid",sev:"warn",ids:[],label:pos?`Missing intensity near ${fmt(be+sh,1)} eV`:`Model above data near ${fmt(be+sh,1)} eV`,
        tip:`The residual summed over ±${th.blobHalf} eV around ${fmt(be+sh,2)} eV is ${fmt(Math.abs(best.z),1)} σ ${pos?"above":"below"} zero and ${fmt(100*Math.abs(best.amp)/pmax,1)} % of the peak height (thresholds ${th.blobZ} σ and ${th.blobAmp*100} %; noise level and its neighbour correlation ρ = ${fmt(nz.rho,2)} (MA(1)) are estimated from the data's second differences, so the test does not depend on the model). `+(pos?"A line or a satellite is missing there (or a width/position is wrong). Add one and compare χ²/BIC — or use the region recipe if this is a satellite-rich region.":"The model is too high there: check the line shapes (asymmetry/tails) and the background.") ,act}); } }
  // R9 complex region → recipe or reference spectra
  const ci_=complexInfo(ctx.el,ctx.orb);
  if(ci_&&!(ctx.recipe&&ctx.recipe.id)){ const rc=ci_.recipes[0];
    push({key:"complex",id:"complex",sev:"info",ids:[],label:rc?`${ctx.el} ${ctx.orb}: use the region recipe`:`${ctx.el} ${ctx.orb}: multiplet/satellite region`,
      tip:ci_.why+". "+(rc?"The recipe sets literature lines with doublet links (ΔE, area ratio, equal widths), the window and the background; derived output (e.g. Ce³⁺ fraction) is shown after the fit.":"No literature recipe is built in for this region: fit only the envelope and use the summed area, or compare with measured reference spectra (Store as reference → LCF), which carry the multiplet structure.")+(ci_.closed?` Closed-shell exception: ${ci_.closed}.`:""),
      act:rc?{type:"recipe",recipe:rc,label:"Apply recipe"}:null}); }
  out.sort((a,b)=>SEV_ORDER[a.sev]-SEV_ORDER[b.sev]);
  return out;
}
/* Diagnosis straight from a problem + fitOnce result (tests and Worker) */
function diagnoseFit(P,r,o){ const C=CORE(); o=o||{}; const fin=o.fin||C.finalize(P,r); return diagnose(Object.assign({E:P.E,y:P.y,fin,bgType:P.bg&&P.bg.type},o),o.th); }

/* ------------------------------------------------------------
   5) One-click actions on a model {lo, hi, bg, comps, cons} (raw BE) — the UI calls this after pushHist (SPEC §17.2)
   ------------------------------------------------------------ */
function applyAction(model,act,o){
  const C=CORE(); o=o||{}; const cs=model.comps, byId={}; cs.forEach(c=>byId[c.id]=c); model.cons=model.cons||{};
  const q=id=>(model.cons[id]=model.cons[id]||{});
  if(act.type==="drop"){ const P={comps:cs,cons:model.cons}; C.removeComp(P,act.id); model.comps=P.comps; model.cons=P.cons; return "Removed "+(act.label||"").replace(/^Remove /,""); }
  if(act.type==="link"){ const c=byId[act.id], a=byId[act.ref]; if(!c||!a) return null; const ka=C.getP(a,act.k), kc=C.getP(c,act.k); q(act.id)[act.k]={mode:"link",ref:act.ref,k:ka>0?+(kc/ka).toFixed(6):1,b:0}; return "Linked"; }
  if(act.type==="fix"){ const c=byId[act.id]; if(!c) return null; C.setP(c,act.k,act.v); q(act.id)[act.k]={mode:"fixed"}; return "Fixed"; }
  if(act.type==="widen"){ const c=byId[act.id]; if(!c) return null; const cs_=q(act.id)[act.k]||{mode:"free"}; const [dlo,dhi]=C.defBounds(c,act.k); const lo=cs_.lo!=null?cs_.lo:dlo, hi=cs_.hi!=null?cs_.hi:dhi; q(act.id)[act.k]={mode:"free",lo:+(lo-act.by).toFixed(4),hi:+(hi+act.by).toFixed(4)}; return "Widened"; }
  if(act.type==="offset"){ const c=byId[act.id], a=byId[act.ref]; if(!c||!a) return null; q(act.id).be={mode:"link",ref:act.ref,k:1,b:+(c.be-a.be).toFixed(4)}; return "Separation fixed"; }
  if(act.type==="window"){ const be0=o.be0!=null?o.be0:-Infinity, be1=o.be1!=null?o.be1:Infinity; if(act.side==="low") model.lo=Math.max(be0,model.lo-act.by); else model.hi=Math.min(be1,model.hi+act.by); return "Window widened"; }
  if(act.type==="bg"){ model.bg=Object.assign({},model.bg,{type:act.bg}); if(act.bg!=="active"){ delete model.bg.conv; } return "Background changed"; }
  if(act.type==="add"){ const anc=cs.reduce((a,c)=>!a||c.area>a.area?c:a,null); const c=C.newComp(anc?anc.shape:"GLS",act.be,{F:act.F,el:anc?anc.el:"",line:anc?anc.line:""});
    if(anc&&anc.m!=null) c.m=anc.m; const hh=C.trueFWHM(c.shape,C.shapeParams(c,1)).ymax; c.area=Math.max(1e-6,(act.h||1)/(hh||1));
    q(c.id).be={mode:"free",lo:+(act.be-1).toFixed(4),hi:+(act.be+1).toFixed(4)}; q(c.id).F={mode:"free",lo:0.3,hi:Math.max(4,act.F*2)}; if(anc&&anc.m!=null) q(c.id).m={mode:"link",ref:anc.id,k:1,b:0};
    cs.push(c); return "Added"; }
  return null;
}

/* ------------------------------------------------------------
   6) Staged fit (SPEC §17.4): (1) heights only → (2) positions within their ranges, widths linked → (3) the user's constraints
      (normal fitOnce, stall rule included). The normal fit from the same start is run too; the lower χ²_d (data-variance χ², the
      common objective of both) is kept, so the result is never worse than the normal fit.
   ------------------------------------------------------------ */
const WKEYS=["F","wG","wL"];
function cpParams(dst,src){ dst.comps=CORE().clone(src.comps); ["c0","k","s","_khi"].forEach(k=>{ if(src.bg&&src.bg[k]!=null) dst.bg[k]=src.bg[k]; }); }
function stageCons(P,user,stage){
  const C=CORE(); const cons=C.clone(user||{}); const get=(id,k)=>(cons[id]||{})[k];
  if(stage===1){ P.comps.forEach(c=>{ const q=cons[c.id]=cons[c.id]||{}; C.paramNames(c).forEach(k=>{ if(k==="area") return; q[k]={mode:"fixed"}; }); }); return cons; }
  // stage 2: positions as the user set them; widths linked to the anchor's (current ratios); all other shape/doublet parameters fixed
  P.comps.forEach(c=>{ const q=cons[c.id]=cons[c.id]||{}; C.paramNames(c).forEach(k=>{ if(k==="area"||k==="be"||WKEYS.includes(k)) return; q[k]={mode:"fixed"}; }); });
  WKEYS.forEach(k=>{ const has=P.comps.filter(c=>C.paramNames(c).includes(k)&&!C.derivedParam(c,k)&&(get(c.id,k)||{}).mode!=="fixed");
    if(has.length<2) return; const anc=has.reduce((a,c)=>!a||c.area>a.area?c:a,null); const ua=get(anc.id,k);
    cons[anc.id][k]=ua&&ua.mode==="free"?ua:{mode:"free"};
    has.forEach(c=>{ if(c===anc) return; const va=C.getP(anc,k), vc=C.getP(c,k); cons[c.id][k]={mode:"link",ref:anc.id,k:va>0?vc/va:1,b:0}; }); });
  return cons;
}
function stagedFit(prob,opt){
  const C=CORE(); opt=opt||{}; const t0=Date.now();
  const base=()=>{ const P=C.clone(C.stripProb(prob)); P.E=prob.E; P.y=prob.y; return P; };
  const P0=base(); { const m=C.buildMap(P0); C.applyLinks(P0,m); }
  const user=C.clone(prob.cons||{}); const stages=[];
  // normal path from the same start
  const PN=base(); const rN=C.fitOnce(PN,opt); const sN=C.fitStats(PN,rN);
  // stage 1: heights (areas) only
  const P1=base(); cpParams(P1,P0); P1.cons=stageCons(P1,user,1); const r1=C.fitOnce(P1,opt); stages.push({stage:1,chi2:r1.chi2,iter:r1.iter+(r1.iter1||0),stop:r1.stop});
  // stage 2: positions within their ranges, widths linked
  const P2=base(); cpParams(P2,P1); P2.cons=stageCons(P2,user,2); const r2=C.fitOnce(P2,Object.assign({},opt,{initAreas:false})); stages.push({stage:2,chi2:r2.chi2,iter:r2.iter+(r2.iter1||0),stop:r2.stop});
  // stage 3: the user's constraints, normal fit path (stall rule applies)
  const P3=base(); cpParams(P3,P2); P3.cons=user; const r3=C.fitOnce(P3,Object.assign({},opt,{initAreas:false})); const s3=C.fitStats(P3,r3); stages.push({stage:3,chi2:r3.chi2,iter:r3.iter+(r3.iter1||0),stop:r3.stop});
  const useStaged=s3.chi2d<sN.chi2d*(1-1e-9);
  const P=useStaged?P3:PN, r=useStaged?r3:rN;
  return {P,r,staged:{used:useStaged?"staged":"normal",chi2dStaged:s3.chi2d,chi2dNormal:sN.chi2d,chi2Staged:r3.chi2,chi2Normal:rN.chi2,stopStaged:r3.stop,stopNormal:rN.stop,stages,ms:Date.now()-t0}};
}
/* Worker job: staged fit → finalize (same summary as finalJob) + staged report */
function stagedJob(prob,opt){ const C=CORE(); const s=stagedFit(prob,opt); const out=C.finalize(s.P,s.r); out.staged=s.staged; return out; }

const XPS_GUIDE={GUIDE_TH,RECIPES,complexInfo,recipesFor,recipeById,recipeModel,recipeDerived,noiseMA1,d2Corr,diagnose,diagnoseFit,applyAction,stageCons,stagedFit,stagedJob};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_GUIDE; else root.XPS_GUIDE=XPS_GUIDE;
})(typeof self!=="undefined"?self:this);
