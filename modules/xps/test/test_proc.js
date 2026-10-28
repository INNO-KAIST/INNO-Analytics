/* XPS pre-processing engine tests  —  node test/test_proc.js   (from modules/xps; ~1–2 min)
   A) engine vs the independent Python reference test/ref/proc_ref.json (scipy where SPEC_proc requires equality)
   B) synthetic truth: deconvolution of a known doublet, stopping behaviour, Fermi-edge coverage, resampling area
   Set V=1 for every line.  QUICK=1 shortens the coverage studies.                                           */
const path=require("path"), fs=require("fs");
const X=require("../src/xps_proc.js");
const REF=JSON.parse(fs.readFileSync(path.join(__dirname,"ref/proc_ref.json"),"utf8"));
let pass=0, fail=0; const fails=[], NUM={};
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!cond) console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const F64=a=>Float64Array.from(a);
const maxRel=(a,b)=>{ let m=0, s=0; for(let i=0;i<b.length;i++) s=Math.max(s,Math.abs(b[i])); for(let i=0;i<b.length;i++) m=Math.max(m,Math.abs(a[i]-b[i])); return m/(s||1); };
const QUICK=!!process.env.QUICK;
const t0=Date.now(); let tS=t0; const sec=t=>{ const n=Date.now(); if(process.env.V) console.log(`   [${((n-tS)/1000).toFixed(1)} s]`); tS=n; console.log(t); };

sec("== A1 special functions");
{ const r=REF.erf; let e1=0,e2=0; r.x.forEach((x,i)=>{ e1=Math.max(e1,Math.abs(X.erf(x)-r.erf[i])); e2=Math.max(e2,Math.abs(X.erfc(x)-r.erfc[i])/Math.max(r.erfc[i],1e-300)); });
  ok(e1<1e-14,"erf vs scipy.special.erf, abs err < 1e-14","max "+e1.toExponential(2)); ok(e2<1e-13,"erfc vs scipy, rel err < 1e-13","max "+e2.toExponential(2)); NUM.erf=e1; NUM.erfc=e2;
  let ev=0; REF.voigt.forEach(c=>{ ev=Math.max(ev,maxRel(c.x.map(x=>X.voigt(x,c.wG,c.wL)),c.v)); });
  ok(ev<1e-7,"Voigt (Weideman N=40) vs scipy.special.voigt_profile, err/peak < 1e-7","max "+ev.toExponential(2)); NUM.voigt=ev; }

sec("== A2 Savitzky–Golay vs scipy.signal.savgol_filter ("+REF.sg.length+" cases)");
{ let mx=0, worst=""; REF.sg.forEach(c=>{ const y=F64(REF.sg_y[String(c.n)]);
    const v=X.savgol(y,c.w,c.order,{deriv:c.deriv,delta:c.delta,mode:c.mode,cval:c.cval}); const e=maxRel(v,c.out);
    if(e>mx){ mx=e; worst=JSON.stringify({w:c.w,o:c.order,d:c.deriv,m:c.mode,n:c.n}); } });
  ok(mx<=1e-10,"all modes (interp/mirror/nearest/constant/wrap), deriv 0–3, window 5–31, order 0–5: max|Δ|/max|ref| ≤ 1e-10","max "+mx.toExponential(2)+" at "+worst); NUM.sg=mx;
  const vv=F64(REF.sg_var.var); let ve=0; REF.sg_var.cases.forEach(c=>{ const r=X.savgol(F64(REF.sg_y["50"]),c.w,c.order,{mode:c.mode,var:vv}); ve=Math.max(ve,maxRel(r.var,c.var_out)); });
  ok(ve<1e-10,"SG variance propagation Σw²σ² vs scipy impulse responses","max "+ve.toExponential(2)); NUM.sgVar=ve;
  // coefficient sanity: classic 5-point quadratic smoothing (-3,12,17,12,-3)/35 (Savitzky & Golay 1964, Table I)
  const c5=X.savgolCoeffs(5,2,0,1), ref5=[-3,12,17,12,-3].map(v=>v/35); ok(maxRel(c5,ref5)<1e-14,"SG 5-pt quadratic = (−3,12,17,12,−3)/35"); }

sec("== A3 Gaussian smoothing vs scipy.ndimage.gaussian_filter1d");
{ let mx=0; REF.gauss.forEach(c=>{ const v=X.gaussSmooth(F64(c.y),c.sigma,{mode:c.mode}); mx=Math.max(mx,maxRel(v,c.out)); });
  ok(mx<1e-12,REF.gauss.length+" cases (5 modes, σ up to 7 pts incl. kernel wider than the signal)","max "+mx.toExponential(2)); NUM.gauss=mx; }

sec("== A4 resampling");
{ const r=REF.pchip; const R=X.resample(r.x,r.y,{Enew:r.xn,method:"pchip",var:r.var});
  const e=maxRel(R.y,r.out); ok(e<1e-12,"PCHIP vs scipy.interpolate.PchipInterpolator (non-uniform x, flat run, sign changes)","max "+e.toExponential(2)); NUM.pchip=e;
  // first-order variance: compare where PCHIP is differentiable (no zero / sign-changing slope within ±2 intervals)
  const PD=X.pchipDeriv(F64(r.x),F64(r.y)); const bad=k=>{ for(let q=k-2;q<=k+2;q++){ if(q<1||q>=PD.m.length) continue; if(PD.m[q]===0||PD.m[q-1]===0||Math.sign(PD.m[q])!==Math.sign(PD.m[q-1])) return true; } return false; };
  let ev=0, nv=0; const vmx=Math.max(...r.var_out); r.xn.forEach((x,j)=>{ let k=0; while(k<r.x.length-2&&r.x[k+1]<x) k++; if(bad(k)) return; nv++; ev=Math.max(ev,Math.abs(R.var[j]-r.var_out[j])/vmx); });
  ok(ev<1e-6,"PCHIP first-order variance vs finite-difference Jacobian of scipy ("+nv+"/"+r.xn.length+" points away from slope kinks)","max "+ev.toExponential(2)); NUM.pchipVar=ev;
  const P=X.pchipRows(F64(r.x),F64(r.y),F64(r.xn)); const lin=X.applyRows(P.rows,F64(r.y)).y; ok(maxRel(lin,P.val)<1e-12,"PCHIP linearised rows reproduce the value (Euler, 1-homogeneity)");
  const b=REF.rebin; const B=X.resample(b.x,b.y,{Enew:b.xn,method:"rebin",var:b.var});
  const eb=maxRel(B.y,b.out); ok(eb<1e-12,"flux-conserving rebin vs cumulative-integral reference","max "+eb.toExponential(2)); NUM.rebin=eb;
  const ebv=maxRel(B.var,b.var_out); ok(ebv<1e-12,"rebin variance Σw²σ² vs dense reference","max "+ebv.toExponential(2));
  const B2=X.resample(b.x,b.y,{Enew:REF.rebin2.xn,method:"rebin"}); ok(maxRel(B2.y,REF.rebin2.out)<1e-12,"rebin with partly covered edge bins (covered-width normalisation)");
  // area: full-cover rebin conserves Σ y_i ΔE_i exactly
  const E=F64(b.x), ei=X.edgesOf(E); let ain=0; for(let i=0;i<E.length;i++) ain+=b.y[i]*(ei[i+1]-ei[i]);
  const Eo=X.uniformGrid(ei[0]+0.15,ei[E.length]+0.3,0.3); const eo=X.edgesOf(Eo);
  const Eall=F64([...Eo]); const Rb=X.resample(b.x,b.y,{Enew:Eall,method:"rebin"}); let aout=0; for(let j=0;j<Eo.length;j++) if(isFinite(Rb.y[j])) aout+=Rb.y[j]*(Math.min(eo[j+1],ei[E.length])-Math.max(eo[j],ei[0]));
  ok(Math.abs(aout-ain)/ain<1e-13,"rebin conserves channel area exactly","rel "+(Math.abs(aout-ain)/ain).toExponential(2));
  const l=REF.lanczos; const L3=X.resample(l.x,l.y,{Enew:l.xn,method:"lanczos",a:3}), L5=X.resample(l.x,l.y,{Enew:l.xn,method:"lanczos",a:5});
  const el=Math.max(maxRel(L3.y,l.a3),maxRel(L5.y,l.a5)); ok(el<1e-12,"Lanczos a=3,5 vs direct reference","max "+el.toExponential(2)); NUM.lanczos=el;
  // Nyquist alias level formula
  ok(Math.abs(X.aliasLevel(1,X.suggestStep(1,1e-3))-1e-3)<1e-15,"suggestStep inverts aliasLevel");
}

sec("== A5 PSF kernels");
{ let eg=0, ev=0; REF.kernels.forEach(k=>{ const K=X.psfKernel(k.type==="gauss"?{type:"gauss",fwhm:k.fwhm}:{type:"voigt",wG:k.wG,wL:k.wL},k.step);
    if(K.h.length!==k.h.length){ ok(false,"kernel length "+k.type,K.h.length+" vs "+k.h.length); return; }
    const e=maxRel(K.h,k.h); if(k.type==="gauss") eg=Math.max(eg,e); else ev=Math.max(ev,e); });
  ok(eg<1e-13,"Gaussian bin-averaged kernel vs scipy ndtr differences","max "+eg.toExponential(2)); ok(ev<1e-7,"Voigt bin-averaged kernel vs quad(voigt_profile)","max "+ev.toExponential(2)); NUM.kG=eg; NUM.kV=ev;
  const a=REF.adjoint; const K=X.psfKernel({type:"gauss",fwhm:a.fwhm},a.step); const r=F64(a.r);
  const Hr=X.convRefl(r,K.h,K.R,new Float64Array(r.length)), HTr=X.convReflT(r,K.h,K.R,new Float64Array(r.length));
  ok(maxRel(Hr,a.Hr)<1e-13&&maxRel(HTr,a.HTr)<1e-13,"reflect-boundary operator H and adjoint Hᵀ vs dense numpy matrices");
  // asymmetric kernel: <Hx, y> = <x, Hᵀy>
  const h=F64([0.05,0.1,0.4,0.3,0.1,0.05]).subarray(0,5); const hs=F64([0.05,0.15,0.4,0.3,0.1]); const g=X.rng(3); const x=new Float64Array(37).map(()=>g.normal()), y=new Float64Array(37).map(()=>g.normal());
  const Hx=X.convRefl(x,hs,2,new Float64Array(37)), Hty=X.convReflT(y,hs,2,new Float64Array(37)); let s1=0,s2=0; for(let i=0;i<37;i++){ s1+=Hx[i]*y[i]; s2+=x[i]*Hty[i]; }
  ok(Math.abs(s1-s2)<1e-13*Math.abs(s1)+1e-13,"adjoint identity <Hx,y> = <x,Hᵀy> for an asymmetric kernel");
}

sec("== A6 noise model");
{ const n=X.noiseModel(F64(REF.noise.y)); ok(Math.abs(n.kappa-REF.noise.kappa)/REF.noise.kappa<1e-12&&Math.abs(n.s2-REF.noise.s2)/REF.noise.s2<1e-12,"κ and s² vs reference","κ="+n.kappa.toFixed(4)+" (Poisson ⇒ 1)"); }

sec("== A7 deconvolution iterations vs dense numpy reference");
{ const D=REF.deconv; D.cases.forEach(c=>{
    const r=X.deconvolve(D.E,D.y,{psf:{type:"gauss",fwhm:D.psf_fwhm},method:c.method,iter:c.iter==null?undefined:c.iter,reblur:c.reblur,stopRule:c.rule,noise:F64(D.var),mc:0,gainTest:false});
    const e=maxRel(r.o,c.o); ok(r.k===c.k&&e<1e-9,`${c.method}${c.reblur===false?" (no reblur)":""} ${c.iter==null?c.rule+" discrepancy stop ("+r.stop+")":"k="+c.iter}`,`k ${r.k} vs ${c.k}, max|Δo|/max ${e.toExponential(2)}, χ² ${r.chi2.toFixed(3)} vs ${c.chi2.toFixed(3)}`);
    NUM["dc_"+c.method+"_"+(c.iter||"stop_"+c.rule)+(c.reblur===false?"_nr":"")]={k:r.k,err:e}; });
  const pb=X.passband(F64(D.y),X.psfKernel({type:"gauss",fwhm:D.psf_fwhm},D.E[1]-D.E[0]),F64(D.var)); const rp=REF.passband;
  ok(pb.nPerp===rp.nPerp&&pb.mc===rp.mc&&Math.abs(pb.chiPerp-rp.chi2Perp)<1e-9*rp.chi2Perp,"noise-only DCT modes (N⊥, m_c, χ²⊥) vs scipy.fft.dct","N⊥="+pb.nPerp+"/"+D.y.length+", m_c="+pb.mc+", χ²⊥="+pb.chiPerp.toFixed(3)); }

sec("== A8 Fermi edge model vs scipy quad, fit vs scipy least_squares");
{ let mx=0; REF.fermi_model.forEach(c=>{ const v=X.fermiModel(c.E,c.p,{T:c.T,axis:c.axis}); mx=Math.max(mx,maxRel(v,c.y)); });
  ok(mx<1e-9,"model (FD ⊗ G, linear DOS + linear bg; BE & KE; wG ≪ kT and ≫ kT) — closed-form step + GL remainder","max "+mx.toExponential(2)); NUM.fermiModel=mx;
  const f=REF.fermi_fit; const r=X.fermiFit(f.E,f.y,{noise:F64(f.var)});
  let worst=0; ["EF","wG","A","S","b0","b1"].forEach(k=>{ worst=Math.max(worst,Math.abs(r.p[k]-f.p[k])/f.sigma[k]); });
  const es=Math.max(...["EF","wG","A","S","b0","b1"].map(k=>Math.abs(r.sigma[k]-f.sigma[k])/f.sigma[k]));
  ok(worst<1e-4,"fit parameters = least_squares optimum (|Δp|/σ)","max "+worst.toExponential(2)+"; EF "+r.p.EF.toFixed(5)+"±"+r.sigma.EF.toFixed(5)+", wG "+r.p.wG.toFixed(4)+"±"+r.sigma.wG.toFixed(4));
  ok(es<1e-3,"fit σ = sqrt(diag((JᵀWJ)⁻¹)) of least_squares","max rel "+es.toExponential(2)); ok(Math.abs(r.chi2-f.chi2)/f.chi2<1e-8,"χ² equal");
  NUM.fermiFit={dp:worst,dsig:es};
  const p=REF.peak_fit; const q=X.peakPSF(p.E,p.y,{wL:p.wL,noise:F64(p.var)});
  let w2=0; ["E0","A","wG","b0","b1"].forEach(k=>{ w2=Math.max(w2,Math.abs(q.p[k]-p.p[k])/p.sigma[k]); });
  const es2=Math.max(...["E0","A","wG","b0","b1"].map(k=>Math.abs(q.sigma[k]-p.sigma[k])/p.sigma[k]));
  ok(w2<1e-4&&es2<1e-3,"Voigt reference-peak fit (wL fixed) = least_squares","max |Δp|/σ "+w2.toExponential(2)+", σ rel "+es2.toExponential(2)+"; wG "+q.p.wG.toFixed(4)+"±"+q.sigma.wG.toFixed(4)); NUM.peakFit={dp:w2,dsig:es2};
}

/* ================================================================ B: synthetic truth */
const FW2SIG=1/(2*Math.sqrt(2*Math.LN2));
function blurTruth(E,truth,fw){ const K=X.psfKernel({type:"gauss",fwhm:fw},E[1]-E[0]); return X.convRefl(truth,K.h,K.R,new Float64Array(truth.length)); }
function peaksOf(E,o,lo,hi){ // local maxima inside [lo,hi], sorted by height
  const pk=[]; for(let i=1;i<o.length-1;i++) if(E[i]>=lo&&E[i]<=hi&&o[i]>=o[i-1]&&o[i]>o[i+1]){ // parabolic vertex
    const d=(o[i-1]-2*o[i]+o[i+1]); const dx=d!==0?0.5*(o[i-1]-o[i+1])/d:0; pk.push({E:E[i]+dx*(E[1]-E[0]),h:o[i],i}); }
  return pk.sort((a,b)=>b.h-a.h); }

sec("== B1 doublet deconvolution (0.25 eV lines 0.70 eV apart, Gaussian PSF 0.80 eV, Poisson ~4e3 counts peak, 2:1)");
const B1={};
{ const st=0.04, E=X.uniformGrid(-4,4,st), N=E.length, bg=300;
  const p1=-0.35, p2=0.35, wT=0.25, psfFW=0.80;
  const truth=E.map(e=>bg+3000*X.voigt(e-p1,wT*0.8,wT*0.25)+1500*X.voigt(e-p2,wT*0.8,wT*0.25));
  const blur=blurTruth(E,truth,psfFW);
  const pkB=peaksOf(E,blur,-1,1); B1.blurPeaks=pkB.length;
  const methods=["jansson","rl","landweber"], seeds=QUICK?4:12;
  methods.forEach(m=>{ const S={sep:[],ratio:[],area:[],dbl:[],fw1:[],k:[],gain:[],ng:[],resolved:0,valley:[],e1:[],e2:[]};
    for(let s=0;s<seeds;s++){ const g=X.rng(1000+s); const y=blur.map(v=>g.poisson(v));
      const r=X.deconvolve(E,y,{psf:{type:"gauss",fwhm:psfFW},method:m,noise:"poisson",mc:s===0?6:0,gainTest:s===0});
      // highest local maximum on each side of the midpoint (|E| ≤ 1.2 eV), each a genuine maximum (not a window edge)
      const pkL=peaksOf(E,r.o,-1.2,-0.02)[0], pkR=peaksOf(E,r.o,0.02,1.2)[0];
      S.k.push(r.k); if(s===0){ S.gain.push(r.resolution.gain); S.ng.push(r.noise?r.noise.gain:NaN); }
      let vmin=Infinity, im=-1; if(pkL&&pkR) for(let i=pkL.i;i<=pkR.i;i++) if(r.o[i]<vmin){ vmin=r.o[i]; im=i; }
      if(pkL&&pkR&&vmin<pkR.h-0.05*(pkR.h-bg)){ const a=pkL, b=pkR; S.resolved++; S.sep.push(b.E-a.E); S.e1.push(a.E-p1); S.e2.push(b.E-p2);
        S.valley.push((vmin-bg)/(b.h-bg));
        // component areas: split at the valley minimum, background-subtracted, within ±2 eV
        let A1=0,A2=0; for(let i=0;i<N;i++){ if(Math.abs(E[i])>2.0) continue; const v=(r.o[i]-bg)*st; if(i<im) A1+=v; else A2+=v; } S.ratio.push(A1/A2);
        // main-line FWHM = 2 × outer (low-BE side) half width at half maximum above bg
        const hm=bg+(a.h-bg)/2; let l=NaN; for(let i=a.i;i>0;i--) if(r.o[i-1]<hm){ l=E[i-1]+(hm-r.o[i-1])/(r.o[i]-r.o[i-1])*st; break; } S.fw1.push(2*(a.E-l)); }
      let ay=0, ao=0, ad=0, at=0; for(let i=0;i<N;i++){ ay+=y[i]; ao+=r.o[i]; if(Math.abs(E[i])<=2){ ad+=r.o[i]-bg; at+=truth[i]-bg; } } S.area.push(ao/ay); S.dbl.push(ad/at); }
    const mean=a=>a.reduce((s,v)=>s+v,0)/a.length, sd=a=>Math.sqrt(a.reduce((s,v)=>s+(v-mean(a))**2,0)/Math.max(1,a.length-1));
    B1[m]={resolved:S.resolved+"/"+seeds,sep:mean(S.sep),sepSD:sd(S.sep),dE1:mean(S.e1),dE2:mean(S.e2),ratio:mean(S.ratio),ratioSD:sd(S.ratio),dblArea:mean(S.dbl),valley:mean(S.valley),fwhm1:mean(S.fw1),areaMaxDev:Math.max(...S.area.map(v=>Math.abs(v-1))),k:[Math.min(...S.k),Math.max(...S.k)],gain:S.gain[0],noiseGain:S.ng[0]};
    const R=B1[m];
    ok(S.resolved===seeds,`${m}: doublet resolved in every noise realisation (blurred data: ${pkB.length} maximum)`,R.resolved+", k "+R.k.join("–"));
    ok(Math.abs(R.dE1)<0.1&&Math.abs(R.dE2)<0.1&&R.sep>0.7,`${m}: peak positions within 0.10 eV (partial restoration pushes the pair apart, SPEC_proc §5.6)`,"ΔE1 "+R.dE1.toFixed(4)+", ΔE2 "+R.dE2.toFixed(4)+", sep "+R.sep.toFixed(3)+"±"+R.sepSD.toFixed(3)+" (true 0.700)");
    ok(Math.abs(R.dblArea-1)<0.01,`${m}: doublet area (bg-subtracted, ±2 eV) = truth within 1 %`,"ratio "+R.dblArea.toFixed(4)+"; valley-split component ratio "+R.ratio.toFixed(2)+"±"+R.ratioSD.toFixed(2)+" (true 2, informational — quantify by fitting, SPEC_proc §5.6)");
    ok(R.areaMaxDev<2e-3,`${m}: total area Σo/Σy conserved < 0.2 %`,"max dev "+R.areaMaxDev.toExponential(2));
    ok(R.fwhm1<0.70,`${m}: main-line FWHM (2×outer HWHM) narrowed from ~0.86 eV`,"FWHM "+R.fwhm1.toFixed(3)+" eV, gain(point-source) "+R.gain.toFixed(2)+", noise gain "+R.noiseGain.toFixed(2));
  });
  B1.blurFWHM=X.fwhm(E,blur,bg).fwhm;
}

sec("== B2 stopping behaviour (semi-convergence) — discrepancy k vs truth error");
const B2={};
{ const st=0.04, E=X.uniformGrid(-4,4,st), bg=300;
  const truth=E.map(e=>bg+3000*X.voigt(e+0.35,0.2,0.0625)+1500*X.voigt(e-0.35,0.2,0.0625)); const blur=blurTruth(E,truth,0.8);
  const g=X.rng(77); const y=blur.map(v=>g.poisson(v));
  const rmse=o=>{ let s=0; for(let i=0;i<o.length;i++) s+=(o[i]-truth[i])**2; return Math.sqrt(s/o.length); };
  ["jansson","rl"].forEach(m=>{
    const r=X.deconvolve(E,y,{psf:{type:"gauss",fwhm:0.8},method:m,noise:"poisson",mc:8,seed:5});
    const g0=r.noise.gain, ks=[0,Math.round(r.k/2),r.k,4*r.k,16*r.k]; const rows=ks.map(k=>{ const q=k===r.k?r:X.deconvolve(E,y,{psf:{type:"gauss",fwhm:0.8},method:m,noise:"poisson",iter:k,mc:k>r.k?4:8,seed:5,gainTest:false,maxNoiseGain:1.2*g0});
      return {k,chi2N:q.chi2/q.N,rmse:rmse(q.o),noiseGain:q.noise?q.noise.gain:1,rough:q.noise?q.noise.roughness:1,warn:q.warnings.some(w=>/noise amplification/.test(w))}; });
    B2[m]=rows;
    let mono=true; for(let i=1;i<r.chi2Hist.length;i++) if(r.chi2Hist[i]>r.chi2Hist[i-1]*(1+1e-12)) mono=false;
    ok(mono,`${m}: χ² decreases monotonically up to the stop`,"k="+r.k+", χ²/N "+(r.chi2/r.N).toFixed(3));
    ok(rows[2].rmse<rows[0].rmse,`${m}: error vs truth at discrepancy stop < raw data`,rows.map(q=>`k=${q.k}: rmse ${q.rmse.toFixed(1)} ng ${q.noiseGain.toFixed(2)}`).join(" | "));
    ok(rows[4].noiseGain>rows[2].noiseGain*1.5,`${m}: over-iteration (16·k) amplifies noise`,"noise gain "+rows[2].noiseGain.toFixed(2)+" → "+rows[4].noiseGain.toFixed(2)+", roughness "+rows[2].rough.toFixed(2)+" → "+rows[4].rough.toFixed(2));
    ok(rows[4].warn&&!rows[1].warn,`${m}: excessive-noise warning fires when over-iterated (threshold 1.2 × gain at stop), not before`);
  });
  // pure noise, no structure: k = 0 (data consistent with PSF-broadened model)
  const g2=X.rng(9); const fl=E.map(()=>g2.poisson(500));
  const rf=X.deconvolve(E,fl,{psf:{type:"gauss",fwhm:0.8},method:"jansson",noise:"poisson"});
  ok(rf.k<=3,"featureless Poisson data: stops at k ≤ 3","k="+rf.k+" χ²₀/N "+(rf.chi2Hist[0]/rf.N).toFixed(3)); B2.flatK=rf.k;
  // noise 'auto' estimate on the doublet data is close to Poisson (κ≈1)
  const nm=X.noiseModel(y); ok(Math.abs(nm.kappa-1)<0.25,"auto noise estimate κ ≈ 1 on Poisson data","κ="+nm.kappa.toFixed(3)); B2.kappa=nm.kappa;
}

sec("== B3 Fermi edge: E_F and resolution recovery, 95 % interval coverage");
const B3={};
{ const E=X.uniformGrid(-1,1,0.02), tp={EF:0.05,wG:0.35,A:2000,S:-300,b0:50,b1:0}; const mu=X.fermiModel(E,tp,{T:300});
  const n=QUICK?60:250; let cE=0,cW=0,cEa=0,cWa=0, bE=0,bW=0, zE=[], zW=[]; const t1=Date.now();
  for(let s=0;s<n;s++){ const g=X.rng(5000+s); const y=mu.map(v=>g.poisson(v));
    const r=X.fermiFit(E,y,{T:300,noise:"poisson"}); const ra=X.fermiFit(E,y,{T:300});         // poisson (unscaled) and auto (scaled by χ²_red)
    const ze=(r.p.EF-tp.EF)/r.sigma.EF, zw=(r.p.wG-tp.wG)/r.sigma.wG; zE.push(ze); zW.push(zw);
    if(Math.abs(ze)<=1.959964) cE++; if(Math.abs(zw)<=1.959964) cW++;
    if(Math.abs(ra.p.EF-tp.EF)<=1.959964*ra.sigma.EF) cEa++; if(Math.abs(ra.p.wG-tp.wG)<=1.959964*ra.sigma.wG) cWa++;
    bE+=r.p.EF-tp.EF; bW+=r.p.wG-tp.wG; }
  const sd=a=>{ const m=a.reduce((s,v)=>s+v,0)/a.length; return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1)); };
  const lo=QUICK?0.85:0.92, hi=QUICK?1:0.98;
  Object.assign(B3,{n,covEF:cE/n,covWG:cW/n,covEFauto:cEa/n,covWGauto:cWa/n,biasEF:bE/n,biasWG:bW/n,sdzEF:sd(zE),sdzWG:sd(zW),sec:(Date.now()-t1)/1000});
  ok(cE/n>=lo&&cE/n<=hi,`E_F 95 % CI coverage (Poisson weights, n=${n})`,(100*cE/n).toFixed(1)+" %, sd(z) "+B3.sdzEF.toFixed(3)+", bias "+(1000*B3.biasEF).toFixed(2)+" meV");
  ok(cW/n>=lo&&cW/n<=hi,"Gaussian FWHM 95 % CI coverage",(100*cW/n).toFixed(1)+" %, sd(z) "+B3.sdzWG.toFixed(3)+", bias "+(1000*B3.biasWG).toFixed(2)+" meV");
  ok(cEa/n>=lo&&cEa/n<=hi&&cWa/n>=lo&&cWa/n<=hi,"coverage with 'auto' (χ²_red-scaled) covariance",(100*cEa/n).toFixed(1)+" % / "+(100*cWa/n).toFixed(1)+" %");
  ok(Math.abs(B3.biasEF)<0.002&&Math.abs(B3.biasWG)<0.005,"mean bias |ΔE_F| < 2 meV, |ΔwG| < 5 meV");
  // KE axis and narrow resolution
  const Ek=X.uniformGrid(1484.5,1487.5,0.02), tk={EF:1486.0,wG:0.45,A:900,S:-40,b0:30,b1:0}; const g=X.rng(4); const yk=X.fermiModel(Ek,tk,{axis:"KE"}).map(v=>g.poisson(v));
  const rk=X.fermiFit(Ek,yk,{axis:"KE",noise:"poisson"}); ok(Math.abs(rk.p.EF-tk.EF)<3*rk.sigma.EF&&Math.abs(rk.p.wG-tk.wG)<3*rk.sigma.wG,"KE axis edge recovered","EF "+rk.p.EF.toFixed(4)+"±"+rk.sigma.EF.toFixed(4)+", wG "+rk.p.wG.toFixed(4)+"±"+rk.sigma.wG.toFixed(4));
  // PSF from the Fermi fit -> deconvolution kernel
  const K=X.psfKernel(rk.psf,0.02); ok(Math.abs(K.fwhmSampled-rk.p.wG)<0.02,"Fermi-fit PSF → kernel FWHM consistent","kernel "+K.fwhmSampled.toFixed(4));
}

const B4={};
sec("== B4 resampling truth tests");
{ // oversampled Gaussian line rebinned ×4: area conserved, variance /4 for Poisson
  const E=X.uniformGrid(280,290,0.01), mu=E.map(e=>100+5e4*X.gauss(e-285,0.8)); const g=X.rng(8); const y=mu.map(v=>g.poisson(v));
  const r=X.resample(E,y,{step:0.04,method:"rebin",var:y.map(v=>Math.max(v,1))});
  ok(Math.abs(r.areaOut-r.areaIn)/r.areaIn<1e-3,"×4 rebin: trapezoid area preserved < 0.1 %","rel "+(Math.abs(r.areaOut-r.areaIn)/r.areaIn).toExponential(2));
  // aligned grids: an output bin holds 3 whole input channels + 2 halves, Σw² = 3/16 + 2/64 = 0.21875
  let rs=0, cnt=0; for(let j=2;j<r.y.length-2;j++){ const m=100+5e4*X.gauss(r.E[j]-285,0.8); rs+=r.var[j]/m; cnt++; } ok(Math.abs(rs/cnt-0.21875)<0.005,"rebin variance of merged Poisson channels = μ·Σw² (0.21875 for ×4 aligned grids)","mean ratio "+(rs/cnt).toFixed(4));
  // undersampled line up-sampled by pchip: line position/height error small vs linear
  const Ec=X.uniformGrid(280,290,0.25), yc=Ec.map(e=>X.gauss(e-285.07,0.9)); const Ef=X.uniformGrid(280,290,0.01);
  const rp=X.resample(Ec,yc,{Enew:Ef,method:"pchip"}), rl=X.resample(Ec,yc,{Enew:Ef,method:"lanczos",a:4}); let ep=0, el=0; for(let i=0;i<Ef.length;i++){ const t=X.gauss(Ef[i]-285.07,0.9); ep=Math.max(ep,Math.abs(rp.y[i]-t)); el=Math.max(el,Math.abs(rl.y[i]-t)); }
  const pk=X.gauss(0,0.9); ok(ep/pk<0.05&&el/pk<0.005,"up-sampling 0.25→0.01 eV of a 0.9 eV line: max error pchip < 5 % (shape-preserving flattening at the maximum), Lanczos-4 < 0.5 %","pchip "+(100*ep/pk).toFixed(3)+" %, lanczos "+(100*el/pk).toFixed(3)+" %");
  // SG derivative of a Gaussian vs analytic derivative
  const Ed=X.uniformGrid(-3,3,0.01), yd=Ed.map(e=>X.gauss(e,1)); const d1=X.derivative(Ed,yd,{window:21,order:4}); let ed=0, md=0; for(let i=0;i<Ed.length;i++){ const s=1*FW2SIG, t=-Ed[i]/(s*s)*X.gauss(Ed[i],1); ed=Math.max(ed,Math.abs(d1[i]-t)); md=Math.max(md,Math.abs(t)); }
  ok(ed/md<2e-4,"SG first derivative (21 pt, order 4) of a sampled Gaussian vs analytic","rel "+(ed/md).toExponential(2));
}

sec("== B5 measured PSF and reference-peak PSF");
{ const st=0.02, Ek=X.uniformGrid(-2,2,0.013), yk=Ek.map(e=>50+2*e+1e4*X.gauss(e-0.031,0.6));   // off-grid, offset centre, sloped baseline
  const K=X.psfKernel({type:"measured",E:Ek,y:yk},st), G=X.psfKernel({type:"gauss",fwhm:0.6},st);
  let c=0; for(let j=0;j<K.h.length;j++) c+=(j-K.R)*st*K.h[j];
  ok(Math.abs(K.fwhmSampled-0.6)<0.01&&Math.abs(c)<2e-3,"measured PSF: baseline removed, centred, FWHM preserved","FWHM "+K.fwhmSampled.toFixed(4)+", centroid "+c.toExponential(2)+" eV, centre found "+K.info.centroid.toFixed(4));
  const E=X.uniformGrid(-4,4,st), tr=E.map(e=>200+800*X.gauss(e+0.3,0.3)+400*X.gauss(e-0.3,0.3)); const K2=X.psfKernel({type:"gauss",fwhm:0.6},st); const y=X.convRefl(tr,K2.h,K2.R,new Float64Array(E.length));
  const a=X.deconvolve(E,y,{kernel:K,method:"rl",iter:50,noise:1,mc:0,gainTest:false}), b=X.deconvolve(E,y,{kernel:G,method:"rl",iter:50,noise:1,mc:0,gainTest:false});
  let d=0,mx=0; for(let i=0;i<E.length;i++){ d=Math.max(d,Math.abs(a.o[i]-b.o[i])); mx=Math.max(mx,b.o[i]); } ok(d/mx<5e-3,"deconvolution with measured vs analytic Gaussian PSF agree (50 RL iterations)","max rel "+(d/mx).toExponential(2));
  // reference line (Ag 3d5/2-like Voigt, wL fixed) with 'auto' noise → PSF FWHM coverage over 100 trials
  const Ep=X.uniformGrid(366,370,0.025), q={E0:368.25,A:3e4,wG:0.55,wL:0.30,b0:400,b1:20}; const mu=X.peakModel(Ep,[q.E0,q.A,q.wG,q.wL,q.b0,q.b1]);
  let cov=0; const n=QUICK?30:100; for(let s=0;s<n;s++){ const g=X.rng(900+s); const yy=mu.map(v=>g.poisson(v)); const r=X.peakPSF(Ep,yy,{wL:0.30}); if(Math.abs(r.p.wG-q.wG)<=1.959964*r.sigma.wG) cov++; }
  ok(cov/n>=(QUICK?0.85:0.9)&&cov/n<=1,"reference-peak wG 95 % CI coverage ('auto' noise, n="+n+")",(100*cov/n).toFixed(0)+" %"); B4.peakCov=cov/n;
  const bad=X.peakPSF(Ep,mu.map(v=>X.rng(1).poisson(v)),{wL:0.20}); ok(bad.p.wG>q.wG,"mis-specified (too small) wL inflates wG — documented systematic","wG "+bad.p.wG.toFixed(4)+" with wL 0.20 vs true 0.55/0.30");
}

console.log(`\n${pass} passed, ${fail} failed  (${((Date.now()-t0)/1000).toFixed(1)} s)`);
if(fail){ console.log("FAILURES:\n  "+fails.join("\n  ")); }
if(process.env.JSON_OUT) fs.writeFileSync(process.env.JSON_OUT,JSON.stringify({NUM,B1,B2,B3,B4},null,1));
process.exit(fail?1:0);
