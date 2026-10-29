/* CasaXPS Round-4 check of INNO's CasaXPS conversions (SPEC_casa_v4.md, casa_v4_patch.md). Read-only use of xps_core.js.
   node test/ref/check_casa_v4.js [path/to/xps_core.js]     (modules/xps; default ../../src/xps_core.js; ~5 s)
   1) INNO CasaXPS-mode fit (weight "casa", Linear conv casa Av.W 1, region 280.03/297.13 or the read window) of each Round-3/4/5/6 file
      → compToCasa(fitted): fwhm, asym, CasaXPS position (be + posShift); DS area via dsCasaAreaFactor when the core has it
   2) casaToComp(CasaXPS string, CasaXPS pos·fwhm) → width and origin of the fitted INNO component
   Prints deviations; the current core (before casa_v4_patch) and a patched copy can be compared by passing its path. */
const path=require("path"), fs=require("fs");
const CORE=process.argv[2]?path.resolve(process.argv[2]):path.join(__dirname,"../../src/xps_core.js");
global.window=global; const C=require(CORE), IO=require(path.join(__dirname,"../../src/xps_io.js"));
const D=p=>path.join(__dirname,"..",p);
const R=[ // file, Casa string, area, fwhm, pos, rms, dof, asym, region end channel (i1)
 ["K_DS0_100_F1","DS(0,100)",15180.3,0.986552,285.0002,693.123,338,0,361],
 ["K_DS02_100_F1","DS(0.2,100)",27119.4,0.795054,285.0927,2522.3,338,0.2550,361],
 ["K_DS01_0_F1","DS(0.1,0)",21288.6,0.868703,285.0689,1649.87,338,0.1603,361],
 ["K_DS01_100_F1","DS(0.1,100)",20811.2,0.907537,285.0523,1570.19,338,0.1178,361],
 ["K_DS01_100_F1","DS(0.1,100)",19946.5,0.876718,285.0489,2202.23,197,0.1178,220],
 ["K_ST17_LA110_F1","ST(1.7)LA(1,1,0)",14734.8,0.891263,285.0737,596.667,338,0.1529,361],
 ["K_ST3_LA11200_F1","ST(3)LA(1,1,200)",14928.0,0.911423,285.1231,279.251,338,0.0440,361],
 ["K_ST17_LA11200_F1","ST(1.7)LA(1,1,200)",14908.6,0.951278,285.0988,322.28,338,0.0071,361],
 ["K_ST17_GL30_F1_2","ST(1.7)GL(30)",14998.0,1.11634,285.1485,37.373,161,-0.0131,361],
 ["K_LA243_153_0_F1","LA(2.43,1.53,0)",14904.2,0.866617,285.0002,320.787,338,-0.2944,361],
 ["K_LA2_3_30_F1","LA(2,3,30)",14996.5,0.930135,285.0000,19.4075,338,0.1296,361],
 ["K_LA2_3_0_F1","LA(2,3,0)",14996.9,0.892084,285.0000,22.9486,338,0.1936,361],
 ["K_TLA12_3_0_F1","TLA(1.2,3,0)",14898.2,0.911717,285.0004,290.636,338,0.1524,361],
 ["K_TLA12_3_100_F1","TLA(1.2,3,100)",14946.4,0.9551,285.0003,180.058,338,0.0817,361],
 ["Shape_LA153_243","LA(1.53,2.43,0)",14997.0,0.580056,285.0000,15.6509,338,0.2275,361],
 ["V2_LA_153_243_F2","LA(1.53,2.43,0)",14979.3,1.15857,284.9999,60.9276,338,0.2275,361],
 ["V2_GL30T2_F1.2","GL(30)T(2)",14999.8,1.18051,285.0000,6.80617,337,0.2138,361],
 // Round 5
 ["K_ST1_GL0_F1_2","ST(1)GL(0)",14998.2,1.16736,285.1241,28.8428,144,-0.0560,361],
 ["K_ST17_GL0_F1_2","ST(1.7)GL(0)",14999.1,1.1715,285.1595,41.0479,148,-0.0483,361],
 ["K_ST3_GL0_F1_2","ST(3)GL(0)",15000.0,1.18428,285.2315,50.1108,153,-0.0256,361],
 ["K_ST10_GL0_F1_2","ST(10)GL(0)",15002.9,1.1855,285.2791,41.3143,161,0.0226,361],
 ["K_ST17_GL30_F0_6","ST(1.7)GL(30)",14999.0,0.59604,285.0793,49.9146,86,-0.0131,361],
 ["K_ST17_GL30_F2_4","ST(1.7)GL(30)",14995.5,2.38363,285.3171,26.4777,337,-0.0131,361],
 ["K_DS03_0_F1","DS(0.3,0)",37178.9,0.814936,284.9485,3185.66,338,0.4834,361],
 // Round 6
 ["K_ST13_GL0_F1_2","ST(1.3)GL(0)",14998.5,1.16827,285.1573,34.8693,146,-0.0543,361],
 ["K_ST2_GL0_F1_2","ST(2)GL(0)",14999.5,1.17442,285.1935,44.0331,149,-0.0429,361],
 ["K_ST25_GL0_F1_2","ST(2.5)GL(0)",14999.8,1.17948,285.1960,48.4436,151,-0.0339,361],
 ["K_ST4_GL0_F1_2","ST(4)GL(0)",15001.7,1.19204,285.2351,50.5211,155,-0.0121,361],
 ["K_ST6_GL0_F1_2","ST(6)GL(0)",15002.6,1.19598,285.2742,48.7939,158,0.0053,361],
];
const e=x=>(x>=0?"+":"")+x.toExponential(1);
(async()=>{
  console.log("core:",CORE,"| dsCasaAreaFactor:",typeof C.dsCasaAreaFactor==="function");
  for(const [f,casa,A,W,P,rms,dof,asy,i1] of R){
    const s=(await IO.readFile(f+".vms",new Uint8Array(fs.readFileSync(D("data/casaval/"+f+".vms"))))).spectra[0];
    const r0=C.casaToComp(casa,{be:285.05,area:12000,fwhm:W*1.06,strict:true}), comp=r0.comp; comp.id="a";
    const cons={a:{}}; C.paramNames(comp).forEach(q=>{ if(["be","area","F"].includes(q)||(comp.shape==="DS"&&q==="wL")||C.derivedParam(comp,q)) return; cons.a[q]={mode:"fixed"}; });
    const Pb=C.makeProblem({be:s.be,y:s.y,c:1},{lo:s.be[19],hi:s.be[i1],bg:{type:"linear",conv:"casa",aw:1},comps:[comp],cons,weight:"casa"});
    const fr=C.fitOnce(Pb,{}), st=C.fitStats(Pb,fr), c=Pb.comps[0], cz=C.compToCasa(c), Wk=c.shape==="DS"?"wL":"F";
    const cpos=cz.posShift==null?null:c.be+cz.posShift;
    let area=c.area, an=""; if(c.shape==="DS"){ const k=C.dsCasaAreaFactor?C.dsCasaAreaFactor(c,Pb.E[0],Pb.E[Pb.E.length-1]):null; if(k){ area=c.area*k; an=" (×"+k.toFixed(5)+")"; } else an=" (INNO window area)"; }
    const back=C.casaToComp(casa,{be:P,area:A,fwhm:W,strict:true}).comp;
    console.log(`${(f+(i1!==361?"@"+i1:"")).padEnd(20)} ${casa.padEnd(19)} area ${area.toFixed(1)}${an} vs ${A} (${e(area/A-1)}) | fwhm ${cz.fwhm!=null?cz.fwhm.toFixed(6):"—"} vs ${W} (${cz.fwhm!=null?e(cz.fwhm/W-1):"—"})`+
      ` | asym ${cz.asym!=null?cz.asym.toFixed(4):"—"} vs ${asy} | pos ${cpos!=null?cpos.toFixed(4):"null"} vs ${P}${cpos!=null?" ("+e(cpos-P)+" eV)":""} | RMS ${st.casa.rms.toFixed(2)} vs ${rms}, DoF ${st.casa.dof}/${dof}`+
      ` | Casa→INNO ${Wk} ${e(back[Wk]/c[Wk]-1)}, be ${e(back.be-c.be)} eV`);
  }
})();
