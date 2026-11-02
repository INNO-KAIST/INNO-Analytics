/* 견고성 회귀 (2026-11 감사): node test/test_robust_core.js
   피팅 코어(xps_core) — 겹친 한계 · invalid 정지 · 숫자 아닌 세기 · 내림차순 축 · CasaXPS 가중 MC · 긴 연결 사슬 · Shirley 미수렴 대체 · doublet 표시
   읽기(xps_io) — 엄격한 숫자 · mkSpec 표시 · .avg 자리 지킴 · 텍스트 건너뜀 경고 · BOM/UTF-16 · OLE2 잘림/0 채움 · 손상 xlsx/VAMAS 오류 · 시료 이름 괄호 */
"use strict";
const path=require("path"), fs=require("fs");
const C=require("../src/xps_core.js"), IO=require("../src/xps_io.js"), LIB=require("../src/xps_lib.js");
let pass=0, fail=0; const fails=[];
function ok(c,name,info){ if(c) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!c) console.log((c?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const throwsMsg=(f,re)=>{ try{ f(); return null; }catch(e){ return re.test(e.message)?e:null; } };
const synth=(n,f)=>{ const be=[],y=[]; for(let i=0;i<n;i++){ const e=280+i*0.05; be.push(e); y.push(f(e,i)); } return {be,y}; };

(async()=>{
console.log("== 1 겹친 한계 lo = hi → 고정 (감사 r1b·r2)");
{ const sp=synth(200,e=>100+2000*C.gauss(e-285,1.2)*1.2+(e>285?50:0));
  const comp=C.newComp("GLS",285.3,{id:"a",F:1.0,m:0.3});
  const fit=cons=>{ const P=C.makeProblem(sp,{lo:280,hi:290,bg:{type:"shirley",n:5},comps:[comp],cons}); const r=C.fitOnce(P,{}); return {P,r,f:C.finalize(P,r)}; };
  const a=fit({a:{F:{mode:"free",lo:1,hi:1}}}), b=fit({a:{F:{mode:"fixed"}}});
  ok(a.r.iter>0&&a.r.stop==="converged"&&a.r.chi2===b.r.chi2&&a.f.cov.sig.every(s=>s>0),"free F with lo = hi fits like fixed F",`iter ${a.r.iter}, χ² ${a.r.chi2.toFixed(3)} vs fixed ${b.r.chi2.toFixed(3)}, σ ${a.f.cov.sig.map(s=>s.toExponential(1)).join(",")}`);
  const m=C.buildMap(a.P); ok(!m.free.some(f=>f.k==="F")&&m.fixed.some(f=>f.k==="F"&&f.pinned),"buildMap moves lo = hi to fixed");
  const P2=C.makeProblem(sp,{lo:280,hi:290,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",285.3,{id:"a",F:1.4,m:0.3})],cons:{a:{F:{mode:"free",lo:1.1,hi:1.1}}}}); C.buildMap(P2);
  ok(P2.comps[0].F===1.1,"value outside a collapsed range is put on the limit",`F ${P2.comps[0].F}`);
  ok(C.toU(1,1,1)===0&&C.dPdU(0,1,1)===0,"toU/dPdU at lo = hi: 0, not NaN"); }
{ const hv=1486.6, sp=synth(240,e=>200+5000*C.gauss(e-285,1.3)*1.3+3000*C.gauss(e-286.6,1.3)*1.3);
  const comp=(name,pos,area)=>({name,shape:"GL(30)",area:{v:area,lo:0,hi:1e9,link:-1,k:1},fwhm:{v:1.3,lo:1.3,hi:1.3,link:-1,k:1},pos:{v:hv-pos,lo:hv-pos-1,hi:hv-pos+1,link:-1,k:1}});
  const cm={shift:0,regions:[{name:"C 1s",bg:"Shirley",keStart:hv-280.5,keEnd:hv-291.5,aw:2}],comps:[comp("A",284.7,6000),comp("B",286.9,3000)]};
  const m=C.casaModelToINNO(cm,hv,{step:0.05}); const P=C.makeProblem(sp,m); const r=C.fitOnce(P,{}); const f=C.finalize(P,r);
  ok(r.iter>0&&r.stop==="converged"&&Math.abs(P.comps[0].be-285)<0.01&&Math.abs(P.comps[1].be-286.6)<0.01&&f.cov.sig.every(s=>s>0),"CasaXPS model with FWHM lo = hi fits (casaModelToINNO)",`iter ${r.iter}, BE ${P.comps.map(c=>c.be.toFixed(3)).join(", ")}`); }

console.log("== 2 invalid 정지 · 숫자 아닌 세기 · 내림차순 축 (감사 r7)");
{ const sp=synth(200,e=>100+1000*C.gauss(e-285,1.2)); sp.y[50]=NaN;
  const e=throwsMsg(()=>C.makeProblem(sp,{lo:280,hi:290,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",285.4,{id:"a"})],cons:{}}),/숫자가 아닌 세기/);
  ok(!!e&&e.code==="nonFinite"&&e.nonFinite===1,"makeProblem rejects NaN intensity in the region",e&&e.message);
  ok(!!C.makeProblem(sp,{lo:283,hi:290,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",285.4,{id:"a"})],cons:{}}),"NaN outside the region is fine");
  const P=C.makeProblem(synth(200,e=>100+1000*C.gauss(e-285,1.2)),{lo:280,hi:290,bg:{type:"linear",n:5},comps:[C.newComp("GLS",285.4,{id:"a"})],cons:{}}); P.y[10]=Infinity;
  ok(!!throwsMsg(()=>C.fitOnce(P,{}),/숫자가 아닌 세기/),"fitOnce rejects non-finite y on a hand-made problem");
  const P2=C.makeProblem(synth(200,e=>100+1000*C.gauss(e-285,1.2)),{lo:280,hi:290,bg:{type:"linear",n:5},comps:[C.newComp("GLS",285.4,{id:"a"})],cons:{}}); P2.comps[0].F=NaN;
  const r2=C.fitOnce(P2,{initAreas:false}); ok(r2.stop==="invalid"&&r2.conv===false&&r2.iter===0,"non-finite χ² at iteration 0 → stop 'invalid', not 'converged'",`stop ${r2.stop}`);
  const d={be:synth(200,()=>0).be.reverse(),y:synth(200,e=>100).y};
  ok(!!throwsMsg(()=>C.checkData(Float64Array.from(d.be),Float64Array.from(d.y)),/오름차순/),"descending energy axis is rejected by checkData"); }

console.log("== 3 CasaXPS 가중 Monte Carlo (감사 r3)");
{ const R=C.rng(3); const sp=synth(240,e=>{ const mu=200+(e>285?80:0)+5000*C.gauss(e-285,1.3)*1.3; return mu+Math.sqrt(mu)*R.g(); });
  for(const bg of [{type:"tougaard",u3:true,C:551,D:436},{type:"shirley",n:5}]){
    const P=C.makeProblem(sp,{lo:280.5,hi:291.5,bg,comps:[C.newComp("GLS",285,{id:"a",F:1.3,m:0.3})],cons:{},weight:"casa"});
    const r=C.fitOnce(P,{}); const A=P.comps[0].area; const mc=C.mcJob(P,{mu:r.mu,v:r.v,chi2:r.chi2},12,5,{}); const s=C.mcSummarize(mc.samples);
    ok(isFinite(s[0].area.mean)&&Math.abs(s[0].area.mean/A-1)<0.05&&s[0].area.sd<0.1*A,`MC with weight 'casa' (${bg.type}): end channels (v = 1e30) not perturbed`,`area ${A.toFixed(1)}, MC ${s[0].area.mean.toFixed(1)} ± ${s[0].area.sd.toFixed(1)}`); } }

console.log("== 4 연결 사슬 > 8 (감사 r6)");
{ const comps=[], cons={}; const K=12; for(let i=K-1;i>=0;i--){ comps.push(C.newComp("GLS",280,{id:"c"+i,F:1})); if(i>0) cons["c"+i]={be:{mode:"link",ref:"c"+(i-1),k:1,b:1}}; }
  const prob={E:new Float64Array(10),y:new Float64Array(10),bg:{type:"none"},comps,cons}; comps.find(c=>c.id==="c0").be=300;
  const map=C.buildMap(prob); C.applyLinks(prob,map); ok(comps.find(c=>c.id==="c11").be===311,"12-long link chain resolved in one pass",`c11.be ${comps.find(c=>c.id==="c11").be}`);
  const m2={free:map.free,links:map.links,fixed:map.fixed}; comps.find(c=>c.id==="c0").be=200; C.applyLinks(prob,m2); ok(comps.find(c=>c.id==="c11").be===211,"map without linkOrder iterates to convergence");
  cons.c0={be:{mode:"link",ref:"c11",k:1,b:1}}; ok(!!throwsMsg(()=>C.buildMap(prob),/순환/),"link cycle still reported"); }

console.log("== 5 Shirley 미수렴 → 직선 (자동 경로) (감사 r14)");
{ let n=0, fb=0, worst=0; for(let s=1;s<=20;s++){ const R=C.rng(s); const sp={be:[],y:[]}; for(let i=0;i<100;i++){ sp.be.push(280+i*0.1); sp.y.push(1000+i*0.3+Math.sqrt(1000)*R.g()); }
    const model={lo:280,hi:289.9,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",285,{id:"a",area:10})],cons:{}};
    const res=C.sensJob(sp,model,C.sensVariants(model,sp,{}).slice(0,1))[0]; n++;
    const P=C.makeProblem(sp,model); const B=C.bgStatic(P.E,P.y,P.bg);
    if(B.converged===false){ ok(res.bgFallback&&res.bgFallback.to==="linear","sensJob falls back to linear when Shirley does not converge",`seed ${s}`); fb++; }
    else ok(!res.bgFallback,"sensJob keeps a converged Shirley",`seed ${s}`);
    if(res.ok) worst=Math.max(worst,Math.abs(res.comps[0].area)); }
  ok(fb>0,"noise-only windows trigger the fallback",`${fb}/${n}`);
  const sp=synth(240,e=>200+(e>285?80:0)+5000*C.gauss(e-285,1.3)*1.3);
  const P=C.makeProblem(sp,{lo:281,hi:291,bg:{type:"shirley",n:5},comps:[C.newComp("GLS",285,{id:"a",F:1.3})],cons:{}}); ok(C.shirleyFallback(P)===null&&!P.bgFallback,"no fallback on a window with a peak"); }

console.log("== 6 doublet 표시 (감사 #7)");
{ const d=(e,o)=>C.doubletInfo(LIB,e,o);
  const ni=d("Ni","3p"); ok(ni.dE===ni.alt.booklet&&ni.unreliable&&ni.unreliable.includes("sources"),"Ni 3p: Booklet (1.8) instead of Avantage 4.0, flagged",`ΔE ${ni.dE}`);
  const rb=d("Rb","3d"), cr=C.dEConstraint(rb); ok(rb.dE===1&&rb.unreliable&&rb.unreliable.includes("trend")&&cr.hi>=1.5&&cr.lo<=0.7,"Rb 3d: value kept, flagged 'trend', window widened toward Kr 1.2 → Sr 1.8",JSON.stringify(cr));
  ["Xe","Cs","Ba","La"].forEach(e=>{ const x=d(e,"4p"); ok(x&&x.erratic&&x.unreliable.includes("erratic"),`${e} 4p flagged erratic`); });
  ok(!d("Hf","4p").erratic&&!d("Ba","4d").unreliable,"Hf 4p / Ba 4d not flagged");
  const yb=C.multipletInfo("Yb","4d"); ok(yb&&yb.kind==="4f"&&yb.valenceDependent&&/Yb²⁺/.test(yb.valenceDependent.doublet),"Yb 4d: multiplet in core + valenceDependent flag");
  ok(d("La","4d").dE===2.8,"La 4d library ΔE unchanged (2.8; discrepancy with recipe 3.0–3.4 reported in lit_audit.md)");
  // Ba 4p 합성: doublet 을 강요하지 않음
  const sp=synth(400,e=>100+4000*C.gauss(e-(280+8),3)*3); const E=Float64Array.from(sp.be), Y=Float64Array.from(sp.y);
  const sm=C.suggestModel(E,Y,{lib:LIB,el:"Ba",orb:"4p",bg:{type:"linear",n:5}}); ok(sm.comps.length>0&&sm.comps.every(c=>!c.doublet)&&sm.mult&&sm.mult.kind==="erratic","suggestModel: Ba 4p → envelope, no doublet",`${sm.comps.length} comps`); }

console.log("== 7 읽기: 엄격한 숫자 · 표시 (감사 t7_nan · fuzz2)");
{ ok(isNaN(IO.numStrict("1.#QNAN"))&&isNaN(IO.numStrict("1.#INF"))&&isNaN(IO.numStrict("12abc"))&&isNaN(IO.numStrict(""))&&IO.numStrict(" -1.5e3 ")===-1500,"numStrict: whole-cell decimals only ('1.#QNAN', '12abc', '' → NaN)");
  const vms=fs.readFileSync(path.join(__dirname,"data","casaval","C1s.vms"),"latin1"); const L=vms.split(/\r?\n/);
  const sp0=IO.vamasParse(vms,"C1s.vms")[0]; const n=sp0.n; let k=-1, hits=0;
  // 세기 값 줄 하나를 "1.#QNAN" 으로: 끝에서 거꾸로 찾아 첫 숫자 줄(세기)
  for(let i=L.length-1;i>0;i--){ if(/^\s*-?\d/.test(L[i])&&hits++===5){ k=i; break; } }
  const L2=L.slice(); L2[k]="1.#QNAN"; const sp=IO.vamasParse(L2.join("\n"),"C1s.vms")[0];
  ok(sp.meta.nonFinite===1&&sp.warn.some(w=>/숫자가 아닌 세기 1/.test(w))&&Array.from(sp.y).filter(v=>!isFinite(v)).length===1,"VAMAS '1.#QNAN' → NaN + flag (not 1)",`n ${n}`);
  const txt="BE CPS\n"+Array.from({length:20},(_,i)=>`${280+i*0.1} ${i===7?"1.#QNAN":1000+i}`).join("\n");
  const t=IO.textParse(txt,"C1s.txt")[0]; ok(t.n===19&&t.warn.some(w=>/건너뜀/.test(w)),"text: row with '1.#QNAN' intensity skipped with a warning",`n ${t.n}`);
  const avg=["DS_EXT_SUPROPID_TITLE : VT_BSTR = 'C1s'","$SPACEAXES=1","0=   280, 0.1, 6, ENERGY, BE, 'BE', 'eV', 'Binding Energy'","$DATA=","LIST@    0=1000, 1001, #empty#, 1.#QNAN, 1004, 1005"].join("\n");
  const a=IO.avgParse(avg,"C1s.avg")[0]; ok(a.n===6&&a.meta.nonFinite===2&&a.y[4]===1004,".avg keeps the place of non-numeric cells (NaN, flagged) — later channels not shifted",`y ${Array.from(a.y).join(",")}`);
  const m=IO.mkSpec({name:"x",be:[1,2,3],y:[1,NaN,3],format:"t"}); ok(m.meta.nonFinite===1&&m.warn.length===1,"mkSpec flags non-finite y"); }

console.log("== 8 읽기: BOM · UTF-16 · 시료 이름 괄호");
{ const txt="Binding Energy\tCPS\n"+Array.from({length:10},(_,i)=>`${280+i*0.1}\t${1000+i}`).join("\n");
  const b=Buffer.from(txt,"latin1");
  for(const [nm,u] of [["UTF-8 BOM",Buffer.concat([Buffer.from([0xEF,0xBB,0xBF]),b])],["UTF-16LE BOM",Buffer.concat([Buffer.from([0xFF,0xFE]),Buffer.from(txt,"utf16le")])],["UTF-16LE no BOM",Buffer.from(txt,"utf16le")],["UTF-16BE BOM",Buffer.concat([Buffer.from([0xFE,0xFF]),Buffer.from(txt,"utf16le").swap16()])]]){
    let r=null, err=""; try{ r=await IO.readFile("C1s.txt",new Uint8Array(u),{}); }catch(e){ err=e.message; }
    ok(r&&r.spectra[0].n===10&&Math.abs(r.spectra[0].be[0]-280)<1e-9&&r.spectra[0].meta.axis==="BE",`text in ${nm} reads like plain text`,err); }
  const vms=fs.readFileSync(path.join(__dirname,"data","casaval","C1s.vms"),"latin1"); const L=vms.split(/\r?\n/);
  const sp0=IO.vamasParse(vms,"C1s.vms")[0]; const bi=L.indexOf(sp0.meta.blockId); const L2=L.slice(); L2[bi]="[D.I. [600'C] 5h] C1s";
  const sp=IO.vamasParse(L2.join("\n"),"C1s.vms")[0]; ok(sp.sample==="D.I. [600'C] 5h","VAMAS sample name with nested brackets",sp.sample);
  ok(IO.sampleFromName("[A [x] B] Ta4f.txt")==="A [x] B"&&IO.sampleFromName("[D.I. on BCT (600'C 5h)] Ta4f.txt")==="D.I. on BCT (600'C 5h)","sampleFromName: balanced brackets"); }

console.log("== 9 읽기: 손상 파일 — 멈추지 않고 읽을 수 있는 오류");
{ const xls=new Uint8Array(fs.readFileSync(path.join(__dirname,"data","bct600","1_D.I._on_BCT_600C_5h.xls")));
  for(const [nm,u] of [["cut 90 %",xls.slice(0,xls.length*0.9|0)],["cut 50 %",xls.slice(0,xls.length*0.5|0)],["cut 600 B",xls.slice(0,600)],["zero-filled",(()=>{ const z=xls.slice(); z.fill(0,512); return z; })()],["zero-filled after FAT",(()=>{ const z=xls.slice(); z.fill(0,xls.length>>1); return z; })()]]){
    const t=Date.now(); let err=""; try{ await IO.readFile("a.xls",u,{}); }catch(e){ err=e.message; }
    ok(err&&/잘렸|손상|순환/.test(err)&&Date.now()-t<3000,`OLE2 ${nm}: error, no hang`,`${Date.now()-t} ms: ${err.slice(0,80)}`); }
  let r=null; try{ r=await IO.readFile("a.xls",xls,{}); }catch(e){} ok(r&&r.spectra.length===7,"intact xls still reads (7 spectra)");
  const zip=new Uint8Array([0x50,0x4B,3,4,...new Array(60).fill(7)]); let ez=""; try{ await IO.readFile("a.xlsx",zip,{}); }catch(e){ ez=e.message; }
  ok(/xlsx 파일을 읽지 못했습니다/.test(ez),"corrupt xlsx: readable error",ez.slice(0,90));
  const vms=fs.readFileSync(path.join(__dirname,"data","casaval","C1s.vms"),"latin1"); let ev=""; try{ await IO.readFile("C1s.vms",new Uint8Array(Buffer.from(vms.slice(0,vms.length>>1),"latin1")),{}); }catch(e){ ev=e.message; }
  ok(/VAMAS/.test(ev),"truncated VAMAS: readable error",ev.slice(0,90)); }

console.log(`\n${pass} 통과 · ${fail} 실패`); if(fail){ console.log(fails.join("\n")); process.exit(1); }
})();
