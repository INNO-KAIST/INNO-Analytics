/* ============================================================
   INNO Analytics — DRE 탭 (rate equation 시뮬레이션·피팅)
   계산은 DRE_CORE. 피팅은 Web Worker 에서(화면이 멈추지 않게, 취소 가능).
   화면 글: 출처와 함정만.
   ============================================================ */
(function(){
"use strict";
const D=window.DRE_CORE;
const KEY="innoanalytics_dre_v1";
const esc=s=>String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const $q=(r,s)=>r.querySelector(s);
const note=(m,k)=>{ try{ if(typeof toast==="function") toast(m,k); }catch(e){} };
const sci=(v,d)=>{ if(v==null||!isFinite(v)) return "—"; if(v===0) return "0"; const a=Math.abs(v); if(a>=1e-2&&a<1e5) return (+v.toPrecision(d||4)).toLocaleString("en-US",{maximumFractionDigits:6}); const e=Math.floor(Math.log10(a)); const m=v/Math.pow(10,e); return `${m.toFixed((d||4)-1)}×10<sup>${e}</sup>`; };
const sciTxt=(v,d)=>v==null||!isFinite(v)?"":Number(v).toPrecision(d||4);
const PU={"W/cm²":1,"kW/cm²":1e3,"MW/cm²":1e6,"W/m²":1e-4};
/* dre_core.js 가 만들어 주는 한국어 라벨(적합성 제목·과정 이름·예측값 이름)을 화면용 영어로 바꾼다 — core 는 건드리지 않는다 */
const ADQ_TITLE={A0:"Adequacy check",A1:"Residual magnitude",A2:"Systematic residual structure",A3:"Windowed shape comparison (model vs. data)",A4:"Up/down consistency",A5:"Intermediate values between branches",A6:"Time trace"};
const adqTitle=q=>q.key==="A6"?"Time trace"+String(q.title||"").replace(/^시간 트레이스/,""):(ADQ_TITLE[q.key]||q.title);
const enL=s=>String(s==null?"":s).replace(/흡수 속도 /g,"Absorption rate ").replace(/전달 fraction /g,"Transfer fraction ").replace(/\(저여기\)/g,"(low excitation)")
  .replace(/threshold I_th \((.*) 최대 기울기 위치\)/,"Threshold I_th (max-slope position, $1)").replace(/최대 기울기 s_max/,"Max slope s_max")
  .replace(/^흡수 /,"Absorption ").replace(/^감쇠 /,"Decay ").replace(/^전달 /,"Transfer ").replace(/· 상승 /,"· Rise ").replace(/· 감쇠 /,"· Decay ")
  .replace(/(\d)준위/g,"$1-level");

/* ---------------- 상태 ---------------- */
function defState(){
  const T=D.templateLee();
  return {v:1,template:"lee",name:T.name,spec:T.spec,nm:T.nm,beam:"peak",params:T.params,samples:T.samples.map(s=>Object.assign({data:[]},s)),
    pu:"kW/cm²",simRange:[0.1,200],showDown:true,showSlope:true,tab:"param",fit:{starts:40,maxit:150,seed:1,profile:true,span:2,npts:13},result:null,modelNote:"Lee et al., Nature 589, 230 (2021) 공개 코드의 3준위 모델과 값"};
}
let S=load();
function load(){ try{ const s=JSON.parse(localStorage.getItem(KEY)||"null"); if(s&&s.v===1&&s.spec) return Object.assign(defState(),s); }catch(e){} return defState(); }
function save(){ try{ localStorage.setItem(KEY,JSON.stringify(S)); }catch(e){ /* 데이터가 크면 결과를 빼고 다시 */ try{ const o=Object.assign({},S,{result:null}); localStorage.setItem(KEY,JSON.stringify(o)); }catch(e2){} } }
let ROOT=null, MODEL=null, WORKER=null, RUN=null, SIM=null;

/* ---------------- 모델 ---------------- */
function compileModel(){
  try{ MODEL=D.compile(S.spec); MODEL.err=null; }catch(e){ MODEL={err:String(e.message||e)}; }
  return MODEL;
}
function envOf(si){
  const env={}; Object.entries(S.params).forEach(([n,p])=>env[n]=+p.v);
  const smp=S.samples[si||0]; if(smp&&smp.vars) Object.assign(env,smp.vars);
  if(S.result&&S.result.bestBySample&&S.useFit){ Object.assign(env,S.result.bestBySample[si||0]||{}); }
  return env;
}
/* 이름 → 쓰임새(표시용) */
function paramUse(){
  const use={}; (S.spec.procs||[]).forEach(p=>{ const ex=p.type==="abs"?p.sigma:p.type==="dec"?p.k:p.W; (ex.match(/[A-Za-z_]\w*/g)||[]).forEach(n=>{ (use[n]=use[n]||[]).push(p); }); });
  return use;
}

/* ---------------- Energy Diagram 에서 가져오기 ---------------- */
function importFromED(){
  const E=window.ED_UI, C=window.ED_CORE;
  if(!E||!E._model){ note("Energy Diagram 탭을 먼저 열어 도식을 만드세요","warn"); return; }
  const st=E._state(); const M=E._model()||null;
  if(!M||!M.procs||!M.procs.length){ note("Energy Diagram 에 과정이 없습니다 — 광흡수·에너지 전달·이완을 먼저 그리세요","warn"); return; }
  const spId=ion=>ion.replace(/(\d)\+$/,"").replace(/[^A-Za-z]/g,"")+ion.replace(/\D/g,"");
  const sp={}; const levelsUsed={};
  const use=(m,l)=>{ const k=m.c.ion; (levelsUsed[k]=levelsUsed[k]||new Map()).set(l.id,l); sp[k]=m; };
  M.procs.forEach(p=>{ if(p.broken) return; if(p.kind==="ET"){ use(p.md,p.Df); use(p.md,p.Dt); use(p.ma,p.Af); use(p.ma,p.At); } else { use(p.m,p.F); use(p.m,p.T); } });
  const species=[], params={}, procs=[], placeholder=[];
  const pname=(pre,ion,a,b)=>`${pre}_${spId(ion)}_${a.replace(/[^A-Za-z0-9]/g,"")}${b?"_"+b.replace(/[^A-Za-z0-9]/g,""):""}`;
  Object.entries(levelsUsed).forEach(([ion,map])=>{
    const m=sp[ion]; const g=m.lv[0];
    const lv=[g.id,...[...map.values()].filter(l=>l.id!==g.id).sort((a,b)=>a.E-b.E).map(l=>l.id)];
    species.push({id:spId(ion),ion,levels:[...new Set(lv)],cvar:"c_"+spId(ion)});
  });
  const addP=(nm,v,unit,src,lo,hi)=>{ if(!params[nm]) params[nm]={v,unit,lo:lo||v/1e4,hi:hi||v*1e4,src,placeholder:true}; placeholder.push(nm); return nm; };
  const decays=new Set();
  M.procs.forEach(p=>{
    if(p.broken) return;
    if(p.kind==="PH"){ const n=pname("s",p.m.c.ion,p.F.id,p.T.id); const g=p.c?p.c.grade:4;
      addP(n,g<=1?1e-20:1e-23,"cm²",`placeholder — ${p.type} ${p.F.id}→${p.T.id} 단면적 (판정 ${p.c?p.c.gradeName:"광원 없음"})`,1e-27,1e-17);
      procs.push({type:"abs",sp:spId(p.m.c.ion),from:p.F.id,to:p.T.id,sigma:n,note:`${p.type} ${p.F.id}→${p.T.id}`}); }
    else if(p.kind==="MPR"||p.kind==="RAD"){ const n=pname(p.kind==="MPR"?"kmpr":"A",p.m.c.ion,p.F.id,p.T.id);
      addP(n,p.kind==="MPR"?1e5:1e3,"s⁻¹",`placeholder — ${p.kind} ${p.F.id}→${p.T.id}`,1e-1,1e9);
      procs.push({type:"dec",sp:spId(p.m.c.ion),from:p.F.id,to:p.T.id,k:n,note:`${p.kind} ${p.F.id}→${p.T.id}`}); decays.add(p.m.c.ion+"|"+p.F.id); }
    else if(p.kind==="ET"){ const n=pname("W"+p.type,p.md.c.ion,p.Df.id,p.At.id);
      addP(n,1e4,"s⁻¹ per mole fraction",`placeholder — ${p.type} (mean-field, 농도 비례)`,1,1e8);
      procs.push({type:"et",d:[spId(p.md.c.ion),p.Df.id,p.Dt.id],a:[spId(p.ma.c.ion),p.Af.id,p.At.id],W:n,scale:"x",note:`${p.type}`}); }
  });
  /* 감쇠가 없는 들뜬 준위는 바닥으로 총 감쇠를 붙인다 (없으면 population이 쌓여 정상상태가 없다) */
  species.forEach(s=>s.levels.slice(1).forEach(l=>{ if(decays.has(s.ion+"|"+l)) return;
    const n=pname("k",s.ion,l,null); addP(n,1e3,"s⁻¹",`placeholder — ${l} 총 감쇠(바닥으로). 수명 τ = 1/k`,1,1e8);
    procs.push({type:"dec",sp:s.id,from:l,to:s.levels[0],k:n,note:`${l} → 바닥 (총 감쇠)`}); }));
  /* 채널: ED 에 그려진 방출 중 모델 준위에서 나오는 것 */
  const channels=[];
  (M.emis||[]).forEach(e=>{ const s=species.find(q=>q.ion===e.m.c.ion); if(!s||!s.levels.includes(e.up.id)) return;
    const id=`${Math.round(e.nm)} nm (${e.up.t})`; if(!channels.some(c=>c.id===id)) channels.push({id,terms:[[s.id,e.up.id,1]]}); });
  species.forEach(s=>s.levels.slice(1).forEach(l=>{ const id=`${s.id} ${l} population`; if(!channels.some(c=>c.terms[0][1]===l&&c.terms[0][0]===s.id)) channels.push({id,terms:[[s.id,l,1]]}); }));
  const vars={}; species.forEach(s=>{ const col=st.cols.find(c=>c.ion===s.ion); const cc=col&&col.conc&&col.conc.find(v=>v!=null); vars[s.cvar]=cc!=null?cc:1; });
  S.spec={species:species.map(s=>({id:s.id,levels:s.levels,cvar:s.cvar})),procs,channels};
  S.params=params; S.samples=[{id:"S1",vars,data:[]}]; S.nm=(st.lasers&&st.lasers[0])||980;
  S.template="ed"; S.name="Model imported from Energy Diagram — "+(st.formula||species.map(s=>s.ion).join(", "));
  S.modelNote=`placeholder 값 ${Object.keys(params).length}개 — 문헌·측정값으로 바꾸기 전에는 형태만 참고하세요.`;
  S.result=null; save(); rerender(); note("Energy Diagram 의 과정 "+procs.length+"개로 모델을 만들었습니다");
}

/* ---------------- 시뮬레이션 ---------------- */
function simulate(){
  const M=compileModel(); if(M.err) return {err:M.err};
  const [a,b]=S.simRange; const f=PU[S.pu];
  const I=D.logGrid(a*f,b*f,40);
  const out={I,up:[],dn:[],hyst:false,err:null};
  try{
    const env=envOf(0);
    const chIdx=M.channels.map((_,k)=>k);
    chIdx.forEach(k=>{
      const u=D.observe(M,env,I,1,S.nm,k,{beam:S.beam}); const d=D.observe(M,env,I,-1,S.nm,k,{beam:S.beam});
      out.up.push(u.y); out.dn.push(d.y);
      if(u.y.some((v,i)=>Math.abs(v/d.y[i]-1)>0.01)) out.hyst=true;
    });
    const sw=D.sweep(M,env,I,1,S.nm); out.pop=sw.n;
  }catch(e){ out.err=String(e.message||e); }
  return out;
}

/* ---------------- 그래프 (SVG) ---------------- */
function palette(){ const cs=getComputedStyle(document.documentElement); const v=(n,f)=>(cs.getPropertyValue(n)||"").trim()||f; const d=document.documentElement.getAttribute("data-theme")==="dark";
  return {dark:d,ink:v("--ink","#243433"),ink2:v("--ink-2","#4a5a56"),muted:v("--muted","#62706b"),grid:v("--grid","#e4e8e2"),line:v("--border","#dde3dc"),bg:v("--surface","#fff"),accent:v("--accent","#7a3f6e"),
    ser:d?["#e39ad0","#7fc6e8","#f0b26a","#8fd19e","#c9a8f0","#f08c86"]:["#7a3f6e","#1f6fa3","#b4580f","#2a7141","#6a4bb0","#b23a3a"]}; }
function chart(o){
  const P=palette(), W=o.w||640, H=o.h||340, m={l:62,r:16,t:14,b:42};
  const pw=W-m.l-m.r, ph=H-m.t-m.b;
  const xs=[],ys=[]; o.series.forEach(s=>s.x.forEach((x,i)=>{ const y=s.y[i]; if(x>0||!o.xlog) if(isFinite(x)) xs.push(x); if(isFinite(y)&&(!o.ylog||y>0)) ys.push(y); }));
  if(!xs.length||!ys.length) return `<div class="dre-empty">No values to plot</div>`;
  let x0=o.xmin!=null?o.xmin:Math.min(...xs), x1=o.xmax!=null?o.xmax:Math.max(...xs), y0=o.ymin!=null?o.ymin:Math.min(...ys), y1=o.ymax!=null?o.ymax:Math.max(...ys);
  if(o.ylog){ y0=Math.pow(10,Math.floor(Math.log10(y0))); y1=Math.pow(10,Math.ceil(Math.log10(y1))); } else { const pad=(y1-y0)*0.06||1; y0-=pad; y1+=pad; }
  if(o.xlog){ x0=Math.pow(10,Math.floor(Math.log10(x0)*4)/4); x1=Math.pow(10,Math.ceil(Math.log10(x1)*4)/4); }
  const X=v=>m.l+(o.xlog?(Math.log10(v)-Math.log10(x0))/(Math.log10(x1)-Math.log10(x0)):(v-x0)/(x1-x0))*pw;
  const Y=v=>m.t+ph-(o.ylog?(Math.log10(v)-Math.log10(y0))/(Math.log10(y1)-Math.log10(y0)):(v-y0)/(y1-y0))*ph;
  const g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="${esc(o.label||"Chart")}" preserveAspectRatio="xMidYMid meet"><rect x="0" y="0" width="${W}" height="${H}" fill="${P.bg}"/>`];
  const tick=(lg,a,b,axis)=>{ const t=[]; if(lg){ const e0=Math.floor(Math.log10(a)), e1=Math.ceil(Math.log10(b)), st=Math.max(1,Math.ceil((e1-e0)/10)); for(let e=e0;e<=e1;e++){ if((e-e0)%st) continue; const v=Math.pow(10,e); if(v>=a*0.999&&v<=b*1.001) t.push(v); } if(t.length<3) for(let e=Math.floor(Math.log10(a));e<=Math.ceil(Math.log10(b));e++) [2,5].forEach(k=>{ const v=k*Math.pow(10,e); if(v>=a&&v<=b) t.push(v); }); }
    else { const span=b-a, st=Math.pow(10,Math.floor(Math.log10(span/5))); const k=[1,2,5,10].find(k=>span/(k*st)<=(o.nt||7))*st; for(let v=Math.ceil(a/k)*k; v<=b+1e-12; v+=k) t.push(+v.toPrecision(10)); } return t.sort((p,q)=>p-q); };
  const lbl=(v,lg)=> lg? (Math.abs(Math.log10(v)%1)<1e-9? `10<tspan dy="-5" font-size="8">${Math.round(Math.log10(v))}</tspan>` : String(+v.toPrecision(2))) : String(+v.toPrecision(4));
  tick(o.xlog,x0,x1).forEach(v=>{ g.push(`<line x1="${X(v)}" x2="${X(v)}" y1="${m.t}" y2="${m.t+ph}" stroke="${P.grid}"/><text x="${X(v)}" y="${m.t+ph+15}" font-size="10.5" text-anchor="middle" fill="${P.muted}">${lbl(v,o.xlog)}</text>`); });
  tick(o.ylog,y0,y1).forEach(v=>{ g.push(`<line x1="${m.l}" x2="${m.l+pw}" y1="${Y(v)}" y2="${Y(v)}" stroke="${P.grid}"/><text x="${m.l-6}" y="${Y(v)+3.5}" font-size="10.5" text-anchor="end" fill="${P.muted}">${lbl(v,o.ylog)}</text>`); });
  g.push(`<rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}" fill="none" stroke="${P.line}"/>`);
  g.push(`<text x="${m.l+pw/2}" y="${H-6}" font-size="11.5" text-anchor="middle" fill="${P.ink2}">${o.xlabel||""}</text><text x="13" y="${m.t+ph/2}" font-size="11.5" text-anchor="middle" fill="${P.ink2}" transform="rotate(-90 13 ${m.t+ph/2})">${o.ylabel||""}</text>`);
  (o.hlines||[]).forEach(h=>{ if(h.y<y0||h.y>y1) return; g.push(`<line x1="${m.l}" x2="${m.l+pw}" y1="${Y(h.y)}" y2="${Y(h.y)}" stroke="${h.color||P.muted}" stroke-dasharray="4 3"/><text x="${m.l+pw-4}" y="${Y(h.y)-4}" font-size="10" text-anchor="end" fill="${h.color||P.muted}">${esc(h.label||"")}</text>`); });
  (o.vlines||[]).forEach(h=>{ if(h.x<x0||h.x>x1) return; g.push(`<line x1="${X(h.x)}" x2="${X(h.x)}" y1="${m.t}" y2="${m.t+ph}" stroke="${h.color||P.muted}" stroke-dasharray="4 3"/>`); });
  g.push(`<clipPath id="dreClip${o.id||""}"><rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}"/></clipPath><g clip-path="url(#dreClip${o.id||""})">`);
  o.series.forEach((s,k)=>{
    const col=s.color||P.ser[k%P.ser.length];
    if(s.line!==false){ let d=""; s.x.forEach((x,i)=>{ const y=s.y[i]; if(!(isFinite(y)&&(!o.ylog||y>0))){ d+=" "; return; } d+=(d&&!d.endsWith(" ")?"L":"M")+X(x).toFixed(1)+" "+Y(y).toFixed(1); });
      g.push(`<path d="${d.replace(/ M/g," M")}" fill="none" stroke="${col}" stroke-width="${s.width||2}" ${s.dash?`stroke-dasharray="${s.dash}"`:""}/>`); }
    if(s.points){ s.x.forEach((x,i)=>{ const y=s.y[i]; if(!(isFinite(y)&&(!o.ylog||y>0))) return; const hol=s.hollow&&s.hollow[i];
      g.push(`<circle cx="${X(x).toFixed(1)}" cy="${Y(y).toFixed(1)}" r="${s.r||3.2}" fill="${hol?P.bg:col}" stroke="${col}" stroke-width="1.3"><title>${esc((s.label||"")+": "+sciTxt(x)+", "+sciTxt(y))}${hol?" (below detection limit — censored)":""}</title></circle>`); }); }
  });
  g.push(`</g>`);
  if(o.legend!==false){ let lx=m.l+8, ly=m.t+14; o.series.filter(s=>s.label).forEach((s,k)=>{ const col=s.color||P.ser[o.series.indexOf(s)%P.ser.length];
    g.push(`<line x1="${lx}" x2="${lx+18}" y1="${ly-4}" y2="${ly-4}" stroke="${col}" stroke-width="2" ${s.dash?`stroke-dasharray="${s.dash}"`:""}/>${s.points?`<circle cx="${lx+9}" cy="${ly-4}" r="3" fill="${col}"/>`:""}<text x="${lx+23}" y="${ly}" font-size="10.5" fill="${P.ink2}">${esc(s.label)}</text>`); ly+=14; }); }
  g.push(`</svg>`); return g.join("");
}

/* ---------------- 렌더 ---------------- */
function render(root){
  ROOT=root; root.innerHTML="";
  compileModel();
  const w=document.createElement("div"); w.className="dre";
  w.innerHTML=`
  <div class="card dre-top">
    <h3>Model <span class="hint" id="dreName"></span></h3>
    <div class="row">
      <select class="dre-sel" id="dreTpl" aria-label="Model"><option value="">Change model…</option><option value="lee">Tm³⁺ 3-level PA — Lee et al. 2021</option><option value="tm5">Tm³⁺ 5-level PA (³H₆–¹G₄) — Chen 2025 · Lee 2021 values</option><option value="tm6">Tm³⁺ 6-level PA (+¹D₂, blue 450 nm)</option><option value="ed">Import from Energy Diagram</option></select>
      <label class="dre-f">Excitation λ <input type="number" id="dreNm" min="200" max="3000" step="1"> nm</label>
      <label class="dre-f">Beam <select class="dre-sel" id="dreBeam"><option value="peak">Peak intensity (single particle)</option><option value="gauss">Gaussian beam integral (ensemble)</option></select></label>
      <label class="dre-f">Power unit <select class="dre-sel" id="drePU">${Object.keys(PU).map(u=>`<option>${u}</option>`).join("")}</select></label>
    </div>
    <div id="dreModelNote"></div>
  </div>
  <div class="dre-main">
    <div class="dre-left">
      <div class="card"><h3>S-curve <span class="hint" id="dreSimHint"></span></h3>
        <div class="row dre-simctl"><label class="dre-f">Range <input type="number" id="dreR0" step="any"> – <input type="number" id="dreR1" step="any"> <span id="drePU2"></span></label>
          <label class="dre-tog"><input type="checkbox" id="dreDown"> down sweep</label><label class="dre-tog"><input type="checkbox" id="dreSlope"> Slope</label>
          <label class="dre-tog" id="dreUseFitL" hidden><input type="checkbox" id="dreUseFit"> Plot with fit result</label>
          <select class="dre-sel" id="dreExport" aria-label="Export"><option value="">Export…</option><option value="csv">Simulation CSV</option><option value="json">Model &amp; results JSON</option><option value="import">Import JSON…</option></select><input type="file" id="dreImp" accept=".json" hidden></div>
        <div id="dreChart"></div><div id="dreSlopeChart"></div><div id="dreSimMsg"></div>
      </div>
      <div class="card" id="dreDiag" hidden></div>
    </div>
    <div class="dre-right">
      <div class="card">
        <div class="seg" id="dreTabs" role="tablist"><button data-t="param">Parameters</button><button data-t="data">Samples &amp; data</button><button data-t="fit">Fit</button><button data-t="model">Rate equations</button></div>
        <div id="drePanel" style="margin-top:10px"></div>
      </div>
    </div>
  </div>`;
  root.appendChild(w);
  $q(w,"#dreTpl").onchange=e=>{ const v=e.target.value; e.target.value=""; if(v==="lee"){ if(!confirm("현재 모델·매개변수·데이터를 Lee 3준위 모델로 바꿉니다. 계속할까요?")) return; S=Object.assign(defState(),{pu:S.pu,simRange:S.simRange}); save(); rerender(); } else if(v==="tm5"||v==="tm6"){ if(!confirm("현재 모델·매개변수·데이터를 "+(v==="tm5"?"Tm 5준위":"Tm 6준위")+" 모델로 바꿉니다. 계속할까요?")) return; const T=D.templateTm(v==="tm5"?5:6);
      Object.values(T.params).forEach(p=>{ if(/^(placeholder|자리표시)/.test(p.src||"")) p.placeholder=true; });
      S=Object.assign(defState(),{template:v,name:T.name,spec:T.spec,nm:T.nm,params:T.params,samples:T.samples.map(q=>Object.assign({data:[]},q)),pu:S.pu,simRange:[1,1e4],modelNote:T.note}); save(); rerender(); }
  else if(v==="ed"){ if(S.samples.some(s=>s.data.length)&&!confirm("데이터가 있는 현재 모델을 바꿉니다. 계속할까요?")) return; importFromED(); } };
  $q(w,"#dreNm").onchange=e=>{ const v=+e.target.value; if(v>0){ S.nm=v; S.result=null; save(); rerender(); } };
  $q(w,"#dreBeam").onchange=e=>{ S.beam=e.target.value; S.result=null; save(); rerender(); };
  $q(w,"#drePU").onchange=e=>{ const f=PU[S.pu]/PU[e.target.value]; S.simRange=S.simRange.map(v=>+(v*f).toPrecision(4)); S.pu=e.target.value; save(); rerender(); };
  $q(w,"#dreR0").onchange=e=>{ const v=+e.target.value; if(v>0&&v<S.simRange[1]){ S.simRange[0]=v; save(); drawSim(); } };
  $q(w,"#dreR1").onchange=e=>{ const v=+e.target.value; if(v>S.simRange[0]){ S.simRange[1]=v; save(); drawSim(); } };
  $q(w,"#dreDown").onchange=e=>{ S.showDown=e.target.checked; save(); drawSim(); };
  $q(w,"#dreSlope").onchange=e=>{ S.showSlope=e.target.checked; save(); drawSim(); };
  $q(w,"#dreUseFit").onchange=e=>{ S.useFit=e.target.checked; save(); drawSim(); };
  $q(w,"#dreTabs").onclick=e=>{ const b=e.target.closest("button[data-t]"); if(!b) return; S.tab=b.dataset.t; save(); renderPanel(); };
  $q(w,"#dreExport").onchange=e=>{ const v=e.target.value; e.target.value=""; if(v==="csv") exportSimCSV(); else if(v==="json") exportJSON(); else if(v==="import") $q(ROOT,"#dreImp").click(); };
  $q(w,"#dreImp").onchange=e=>{ const f=e.target.files&&e.target.files[0]; if(f) importJSON(f); e.target.value=""; };
  rerender();
}
function rerender(){
  if(!ROOT||!ROOT.isConnected) return;
  compileModel();
  $q(ROOT,"#dreName").textContent=enL(S.name||"");
  $q(ROOT,"#dreNm").value=S.nm; $q(ROOT,"#dreBeam").value=S.beam; $q(ROOT,"#drePU").value=S.pu;
  $q(ROOT,"#dreR0").value=S.simRange[0]; $q(ROOT,"#dreR1").value=S.simRange[1]; $q(ROOT,"#drePU2").textContent=S.pu;
  $q(ROOT,"#dreDown").checked=!!S.showDown; $q(ROOT,"#dreSlope").checked=!!S.showSlope;
  const uf=$q(ROOT,"#dreUseFitL"); uf.hidden=!(S.result&&S.result.bestBySample); $q(ROOT,"#dreUseFit").checked=!!S.useFit;
  const ph=Object.values(S.params).filter(p=>p.placeholder).length;
  $q(ROOT,"#dreModelNote").innerHTML=(MODEL.err?`<div class="note danger">모델 오류: ${esc(MODEL.err)}</div>`:"")+
    (ph?`<div class="note warn">placeholder 값 ${ph}개 — 문헌·측정값으로 바꾸기 전의 곡선은 형태 참고용입니다.</div>`:"")+
    (S.modelNote&&!ph?`<div class="note">${esc(S.modelNote)}</div>`:"");
  drawSim(); renderPanel(); renderDiag();
}
function drawSim(){
  const el=$q(ROOT,"#dreChart"); if(!el) return;
  SIM=simulate(); const msg=$q(ROOT,"#dreSimMsg"), hint=$q(ROOT,"#dreSimHint");
  if(SIM.err){ el.innerHTML=`<div class="note danger">${esc(SIM.err)}</div>`; $q(ROOT,"#dreSlopeChart").innerHTML=""; return; }
  const f=PU[S.pu], X=SIM.I.map(v=>v/f); const P=palette();
  const series=[];
  MODEL.channels.forEach((c,k)=>{
    series.push({x:X,y:SIM.up[k],label:c.id+(S.showDown&&SIM.hyst?" (up)":""),color:P.ser[k%6]});
    if(S.showDown&&SIM.hyst) series.push({x:X,y:SIM.dn[k],label:c.id+" (down)",color:P.ser[k%6],dash:"6 4"});
  });
  /* 데이터 겹치기: 채널별, 배율은 피팅 결과가 있으면 그 배율, 없으면 첫 점 맞춤 */
  S.samples.forEach((smp,si)=>smp.data.forEach((d,di)=>{
    if(d.type!=="scurve") return; const k=MODEL.channels.findIndex(c=>c.id===d.ch); if(k<0) return;
    const scale=d._scale|| (S.result&&S.result.scales&&S.result.scales[si+"|"+di]) || autoScale(d,SIM,k);
    const kc=d.k==null?3:d.k;
    series.push({x:d.I.map(v=>v/f),y:d.y.map(v=>v/scale),points:true,line:false,label:`${smp.id} · ${d.label||d.ch}`,color:P.ser[k%6],hollow:d.y.map(v=>!(v>kc*(d.sbg||0)))});
  }));
  el.innerHTML=chart({id:"s",xlog:true,ylog:true,series,xlabel:`excitation intensity (${S.pu}, ${S.beam==="peak"?"peak":"beam center"})`,ylabel:"Signal (normalized population)",label:"S-curve",w:680,h:360});
  const sl=$q(ROOT,"#dreSlopeChart");
  if(S.showSlope){ const ss=MODEL.channels.map((c,k)=>({x:X,y:D.slopes(SIM.I,SIM.up[k]),label:"d ln y / d ln I — "+c.id,color:P.ser[k%6]}));
    sl.innerHTML=chart({id:"sl",xlog:true,ylog:false,series:ss,xlabel:`excitation intensity (${S.pu})`,ylabel:"Slope s",label:"Slope",w:680,h:220,legend:false}); } else sl.innerHTML="";
  const mx=MODEL.channels.map((c,k)=>{ const s=D.slopes(SIM.I,SIM.up[k]); let i=0; s.forEach((v,j)=>{ if(v>s[i]) i=j; }); return `${esc(c.id)}: s<sub>max</sub> ${s[i].toFixed(1)} @ ${sciTxt(X[i],3)} ${S.pu}`; });
  hint.innerHTML=SIM.hyst?'<span class="badge warn">Hysteresis — up/down branches differ</span>':"";
  msg.innerHTML=`<div class="muted" style="margin-top:4px">${mx.join(" · ")}</div>`;
}
function autoScale(d,sim,k){
  /* 데이터와 모델의 log 차 중앙값 (표시용) */
  const lr=[]; d.I.forEach((I,i)=>{ if(!(d.y[i]>0)) return; let q=0; while(q<sim.I.length-1&&sim.I[q]<I) q++; const m=sim.up[k][q]; if(m>0) lr.push(Math.log(d.y[i]/m)); });
  lr.sort((a,b)=>a-b); return lr.length?Math.exp(lr[Math.floor(lr.length/2)]):1;
}

/* ---------------- 옆 패널 ---------------- */
function renderPanel(){
  const el=$q(ROOT,"#drePanel"); if(!el) return;
  ROOT.querySelectorAll("#dreTabs button").forEach(b=>{ const on=b.dataset.t===S.tab; b.classList.toggle("active",on); b.setAttribute("aria-selected",on?"true":"false"); });
  if(S.tab==="param") return panelParams(el);
  if(S.tab==="data") return panelData(el);
  if(S.tab==="fit") return panelFit(el);
  return panelModel(el);
}
function panelParams(el){
  const use=paramUse();
  const names=Object.keys(S.params).sort((a,b)=>(use[b]?1:0)-(use[a]?1:0));
  el.innerHTML=`<div class="dre-tw"><table class="dre-tbl"><thead><tr><th>Name</th><th>Value</th><th>Unit</th><th title="Hold fixed in fit">Fixed</th><th title="시료마다 따로 피팅(농도 계열 전역 피팅에서 시료별 비복사 손실 등)">Per sample</th><th>Search range</th></tr></thead><tbody>`+
    names.map(n=>{ const p=S.params[n]; const u=use[n];
      return `<tr class="${p.placeholder?"ph":""}"><td><b>${esc(n)}</b>${u?"":' <span class="muted">(unused)</span>'}<div class="muted sm">${esc(p.src||"")}</div></td>
        <td><input type="text" data-n="${esc(n)}" data-f="v" value="${sciTxt(p.v,5)}" aria-label="${esc(n)} value"></td><td class="muted">${esc(p.unit||"")}</td>
        <td><input type="checkbox" data-n="${esc(n)}" data-f="fixed" ${p.fixed?"checked":""} aria-label="${esc(n)} fixed"></td>
        <td><input type="checkbox" data-n="${esc(n)}" data-f="per" ${p.per==="sample"?"checked":""} ${p.fixed?"disabled":""} aria-label="${esc(n)} per sample"></td>
        <td class="rng"><input type="text" data-n="${esc(n)}" data-f="lo" value="${sciTxt(p.lo,3)}" aria-label="${esc(n)} lower bound">–<input type="text" data-n="${esc(n)}" data-f="hi" value="${sciTxt(p.hi,3)}" aria-label="${esc(n)} upper bound"></td></tr>`; }).join("")+
    `</tbody></table></div><div class="muted" style="margin-top:6px">값을 바꾸면 곡선이 바로 다시 그려집니다. 탐색 범위는 log 공간 multi-start box입니다.</div>`;
  el.querySelectorAll("input[data-n]").forEach(i=>i.onchange=()=>{
    const p=S.params[i.dataset.n]; const f=i.dataset.f;
    if(f==="fixed"){ p.fixed=i.checked; } else if(f==="per"){ p.per=i.checked?"sample":"global"; }
    else { const v=parseFloat(i.value); if(!(v>0)){ note("양수를 넣으세요","warn"); i.value=sciTxt(p[f]); return; } p[f]=v; if(f==="v") p.placeholder=false; }
    S.result=null; save(); rerender();
  });
}
function panelModel(el){
  const M=MODEL;
  if(M.err){ el.innerHTML=`<div class="note danger">${esc(M.err)}</div>`; return; }
  const lvName=(sp,l)=>`${esc(sp)} ${esc(l)}`;
  el.innerHTML=`<div class="muted">population은 이온 종마다 합이 1 (첫 준위 = 바닥). φ = I/hν, λ = ${S.nm} nm.</div>
   <div class="dre-tw"><table class="dre-tbl"><thead><tr><th>Process</th><th>Rate term</th><th>Note</th></tr></thead><tbody>${S.spec.procs.map(p=>{
     if(p.type==="abs") return `<tr><td>Absorption ${lvName(p.sp,p.from)}→${esc(p.to)}</td><td>(${esc(p.sigma)})·φ·n<sub>${esc(p.from)}</sub></td><td class="muted">${esc(p.note||"")}</td></tr>`;
     if(p.type==="dec") return `<tr><td>Decay ${lvName(p.sp,p.from)}→${esc(p.to)}</td><td>(${esc(p.k)})·n<sub>${esc(p.from)}</sub></td><td class="muted">${esc(p.note||"")}</td></tr>`;
     return `<tr><td>Transfer ${esc(p.d[0])} ${esc(p.d[1])}→${esc(p.d[2])} ⇒ ${esc(p.a[0])} ${esc(p.a[1])}→${esc(p.a[2])}</td><td>(${esc(p.W)})·n<sub>${esc(p.d[1])}</sub>·n<sub>${esc(p.a[1])}</sub>${p.scale==="x"?" × relative-species mole fraction":""}</td><td class="muted">${esc(p.note||"")}</td></tr>`; }).join("")}</tbody></table></div>
   <div class="muted" style="margin-top:6px">Channels: ${M.channels.map(c=>esc(c.id)).join(", ")}</div>
   <details class="dre-src"><summary>Numerical methods</summary><ul>
    <li>정상상태: 앞 파워의 해에서 Newton(선 탐색). branch가 사라지면(threshold) implicit Euler pseudo-time으로 다음 안정 상태까지 흘러간 뒤 Newton 으로 마무리 — up·down 스윕이 측정처럼 안정 branch를 따라가고, bistable이면 히스테리시스가 나온다.</li>
    <li>동역학: L-안정 TR-BDF2 + 단계 배가 오차 제어. 상승·감쇠 시간은 적분 시간이 아니라 목표 도달로 끝낸다(critical slowing down 대응).</li>
    <li>감도: implicit function theorem dn/dθ = −J⁻¹∂f/∂θ (정확). 피팅 기울기는 차분이 아니다.</li>
    <li>검증: Lee 2021 모델 6개 시료를 scipy(LSODA) 독립 구현과 대조 — 정상상태 상대오차 ≤ 2.5×10⁻¹⁰, 동역학 ≤ 5×10⁻⁷; bistable 매개변수에서 up=낮은 안정 branch, down=높은 안정 branch(전근 탐색 + 안정성 판정과 일치).</li></ul></details>`;
}
/* ---------- 데이터 ---------- */
function parseTable(txt){
  const lines=String(txt).replace(/\r/g,"").split("\n").map(l=>l.trim()).filter(l=>l&&!l.startsWith("#"));
  if(!lines.length) return null;
  const sep=lines[0].includes("\t")?"\t":(lines[0].includes(",")?",":(lines[0].includes(";")?";":/\s+/));
  const rows=lines.map(l=>l.split(sep).map(x=>x.trim()));
  const isNum=v=>v!==""&&isFinite(+v);
  let header=null; if(!rows[0].every(isNum)){ header=rows.shift(); }
  const nc=Math.max(...rows.map(r=>r.length));
  const cols=[...Array(nc).keys()].map(j=>rows.map(r=>r[j]==null||r[j]===""||/^nan$/i.test(r[j])?NaN:+r[j]));
  const strs=[...Array(nc).keys()].map(j=>rows.map(r=>r[j]==null?"":String(r[j])));
  return {header:header||cols.map((_,j)=>"Column "+(j+1)),cols,txt:strs,n:rows.length};
}
let PENDING=null;
function panelData(el){
  const chs=MODEL.err?[]:MODEL.channels.map(c=>c.id);
  const cvars=[...new Set(S.spec.species.map(s=>s.cvar).filter(Boolean))];
  el.innerHTML=`
   ${S.samples.map((s,si)=>`<div class="dre-smp"><div class="row"><b>Sample <input type="text" class="dre-sid" data-si="${si}" value="${esc(s.id)}" aria-label="Sample name"></b>
     ${cvars.map(v=>`<label class="dre-f">${esc(v)} <input type="number" step="any" data-si="${si}" data-var="${esc(v)}" value="${s.vars&&s.vars[v]!=null?s.vars[v]:""}"> %</label>`).join("")}
     ${S.samples.length>1?`<button class="btn xs" data-delsmp="${si}">Delete sample</button>`:""}</div>
     ${s.data.length?`<table class="dre-tbl"><tbody>${s.data.map((d,di)=>`<tr><td>${d.type==="scurve"?"S-curve":d.type==="rise"?"Rise time":d.type==="trace"?"Time trace (diagnostic)":"Decay time"}</td><td>${esc(d.ch)}</td><td>${d.type==="trace"?`${d.t.length} pts · ${sciTxt(d.I[0]/PU[S.pu],3)} ${esc(S.pu)}`:`${d.I.length} pts ${d.type==="scurve"?(d.dir<0?"· down":"· up"):""}`}</td><td>${(()=>{ const L=calLabel(d); return `<span class="gr ${L.lvl==="ok"?"g0":"g3"}" title="${esc(L.txt)}">${L.lvl==="ok"?"Calibrated":"Check calibration"}</span> <button class="btn xs" data-cal="${si}|${di}" title="${esc(L.txt)}">Power calibration</button>`; })()}</td><td class="muted">${d.type==="trace"?"피팅에 쓰지 않음 · 적합성 A6":d.type==="scurve"?`Background σ ${sciTxt(d.sbg||0,3)} · multiplicative noise ${Math.round((d.relErr||0.05)*100)}%`:`Relative error ${Math.round((d.relErr||0.05)*100)}%`}</td><td><button class="btn xs" data-del="${si}|${di}" aria-label="Delete data">×</button></td></tr>${issuesRow(d,si,di)}`).join("")}</tbody></table>`:'<div class="muted">No data</div>'}</div>`).join("")}
   <div class="row" style="margin-top:6px"><button class="btn sm" id="dreAddSmp">Add sample</button> <span class="muted">농도 계열은 시료를 나눠 넣고, 비복사 손실만 '시료별'로 두면 identifiability이 크게 오른다.</span></div>
   <hr class="dre-hr">
   <h4 class="dre-h4">Add data</h4>
   <div class="row"><select class="dre-sel" id="dreDType"><option value="scurve">S-curve (power · signal)</option><option value="rise">Rise time (power · time, s)</option><option value="decay">Decay 1/e time (power · time, s)</option><option value="trace">Time trace (time, s · signal) — diagnostic only</option></select>
     <select class="dre-sel" id="dreDSmp">${S.samples.map((s,si)=>`<option value="${si}">${esc(s.id)}</option>`).join("")}</select>
     <label class="dre-f"><input type="file" id="dreFile" accept=".csv,.txt,.tsv,.dat,.h5,.hdf5"></label></div>
   <div class="muted">.h5 (ScopeFoundry power_scan · ni_daq_pulse_resp) 는 형식을 자동으로 알아보고 up·down·트레이스로 나눈다. 표 파일은 CSV·TSV·TXT·DAT.</div>
   <textarea id="drePaste" class="dre-ta" placeholder="또는 붙여넣기: 첫 열 파워, 다음 열 신호 (쉼표·탭·공백 구분, # 주석 무시, 결측 NaN)"></textarea>
   <div class="row"><button class="btn sm" id="dreParse">Parse</button><button class="btn sm" id="dreSynth" title="현재 매개변수로 참값을 아는 합성 데이터를 만든다 — 검증용">Generate synthetic data</button></div>
   <div id="dreMap"></div>`;
  el.querySelectorAll(".dre-sid").forEach(i=>i.onchange=()=>{ const v=i.value.trim(); if(v&&!S.samples.some((s,k)=>s.id===v&&k!==+i.dataset.si)){ S.samples[+i.dataset.si].id=v; save(); } else i.value=S.samples[+i.dataset.si].id; });
  el.querySelectorAll("input[data-var]").forEach(i=>i.onchange=()=>{ const s=S.samples[+i.dataset.si]; s.vars=s.vars||{}; s.vars[i.dataset.var]=+i.value; S.result=null; save(); drawSim(); });
  el.querySelectorAll("button[data-cal]").forEach(b=>b.onclick=()=>{ const [si,di]=b.dataset.cal.split("|").map(Number); calUI(el,si,di); const m=$q(el,"#dreMap"); if(m&&m.scrollIntoView) m.scrollIntoView({block:"nearest"}); });
  el.querySelectorAll("button[data-fix]").forEach(b=>b.onclick=()=>{ const [si,di,act]=b.dataset.fix.split("|"); fixData(el,+si,+di,act); });
  el.querySelectorAll("button[data-del]").forEach(b=>b.onclick=()=>{ const [si,di]=b.dataset.del.split("|").map(Number); S.samples[si].data.splice(di,1); S.result=null; save(); rerender(); });
  el.querySelectorAll("button[data-delsmp]").forEach(b=>b.onclick=()=>{ if(!confirm("시료와 그 데이터를 삭제할까요?")) return; S.samples.splice(+b.dataset.delsmp,1); S.result=null; save(); rerender(); });
  $q(el,"#dreAddSmp").onclick=()=>{ const n=S.samples.length+1; const vars=Object.assign({},S.samples[0]&&S.samples[0].vars); S.samples.push({id:"S"+n,vars,data:[]}); save(); renderPanel(); };
  const readTxt=(t,name)=>{ PENDING=parseTable(t); if(!PENDING||PENDING.cols.length<2){ $q(el,"#dreMap").innerHTML='<div class="note warn">두 열 이상의 숫자 표가 필요합니다</div>'; return; }
    PENDING.name=name||""; if(isPulseTable(PENDING)){ pulseCsvUI(el,chs,PENDING); return; } mapUI(el,chs); };
  $q(el,"#dreFile").onchange=async e=>{ const f=e.target.files&&e.target.files[0]; if(!f) return;
    if(/\.(h5|hdf5)$/i.test(f.name)){ try{ const H=await readH5(f); h5UI(el,chs,H); }catch(err){ $q(el,"#dreMap").innerHTML=`<div class="note warn">${esc(f.name)}: ${esc(err.message||err)}</div>`; } e.target.value=""; return; }
    const r=new FileReader(); r.onload=()=>readTxt(r.result,f.name); r.readAsText(f); e.target.value=""; };
  $q(el,"#dreParse").onclick=()=>readTxt($q(el,"#drePaste").value);
  $q(el,"#dreSynth").onclick=()=>synthUI(el,chs);
}
/* ---------- 파워 보정 (데이터별) ----------
   브라우저에 넣는 데이터는 대개 파워미터(PM) reading이다. 원본은 d.raw 에 그대로 두고, 보정 d.cal 로 d.I(W/cm²)를 다시 만든다.
   I[W/cm²] = P_PM[W] × C × 10⁸ / A_peak[µm²] × 추가배율   (peak intensity 2P/πw² 와 같다: A_peak = πw²/2)  */
const PWU={"mW":1e-3,"W":1,"µW":1e-6,"nW":1e-9};
const CALPRE=[{id:"nir950oil",name:"NIR · 950 DMSP · Oil",C:3.322959765878627,A:1.2393,src:"A_peak 직접 측정(γ=2.2 보정본)"},
  {id:"blue600oil",name:"Blue 450+475 · 600 DMSP · Oil",C:3.322959765878627,A:0.8715,src:"잠정 — 기존 계산 역산값, 재측정 권장"}];
const CALKEY="innoanalytics_dre_cal_v1";
function calProfiles(){ let u=[]; try{ u=JSON.parse(localStorage.getItem(CALKEY)||"[]"); }catch(e){} return CALPRE.map(p=>Object.assign({builtin:true},p)).concat(u); }
function saveCalProfile(p){ try{ const u=JSON.parse(localStorage.getItem(CALKEY)||"[]").filter(q=>q.name!==p.name); u.push(p); localStorage.setItem(CALKEY,JSON.stringify(u)); }catch(e){} }
function calK(cal){ /* W/cm² per W */ if(!cal) return null; const x=+cal.extra||1; if(cal.kPerMW>0) return cal.kPerMW*1e6*x; if(cal.C>0&&cal.A>0) return 1e8*cal.C/cal.A*x; return null; }
function wheelW(f,deg){ if(f.model==="nd") return f.a*Math.pow(10,-f.k*(deg-f.t0))+f.b; const r=(deg-f.t0)*Math.PI/180*(f.model==="hwp"?2:1); return f.a*Math.cos(r)**2+f.b; }
function calcI(d){
  if(!d.raw) return d.I;
  const R=d.raw, c=d.cal||{};
  if(R.kind==="intensity") return R.P.map(v=>v*PU[R.unit]*(+c.extra||1));
  const k=calK(c); if(!k) return d.I;
  const PW=(c.src==="wheel"&&c.wheelFit&&R.wheel)?R.wheel.map(t=>wheelW(c.wheelFit,t)):R.P.map(v=>v*PWU[R.unit]);
  return PW.map(v=>v*k);
}
function applyCal(d,cal,confirmed){ d.cal=JSON.parse(JSON.stringify(cal)); d.calOK=!!confirmed; const I=calcI(d);
  if(d.type==="scurve"||d.type==="trace"||d.type==="rise"||d.type==="decay"){ /* 양수·유한만 유지되도록 원본 인덱스 그대로 */ d.I=I; } }
function calLabel(d){
  if(!d.raw) return {txt:"intensity entered directly",lvl:"ok"};
  const c=d.cal||{}; const k=calK(c);
  if(d.raw.kind==="intensity") return {txt:`intensity (${d.raw.unit})${(+c.extra||1)!==1?` ×${+c.extra}`:""}`,lvl:"ok"};
  const nm=c.name||(c.kPerMW?`${c.kPerMW} kW/cm²/mW`:(c.C&&c.A?`C ${(+c.C).toFixed(3)} · A ${c.A} µm²`:"no conversion"));
  return {txt:`PM ${d.raw.unit} → ${nm}${c.src==="wheel"?" · reconstructed from wheel angle":""}${(+c.extra||1)!==1?` · ×${+c.extra}`:""}${k?` (${(k/1e6).toPrecision(4)} kW/cm² per mW)`:""}`,lvl:d.calOK?"ok":"warn"};
}
/* 휠 모델 3종을 모두 맞추고 상대 rms 가 가장 작은 것을 고른다 (점 ≥ 3)
   pol: a·cos²(θ−θ₀)+b (선형 편광자 회전) · hwp: a·cos²(2(θ−θ₀))+b (반파장판+편광자) · nd: a·10^(−k(θ−θ₀))+b (연속 ND 휠) */
function fitWheel(pts){
  if(pts.length<3) return null; const n=pts.length; const t0r=pts[0][0];
  const lin=g=>{ let sxx=0,sx=0,sy=0,sxy=0; pts.forEach(([t,p])=>{ const x=g(t); sxx+=x*x; sx+=x; sy+=p; sxy+=x*p; }); const det=n*sxx-sx*sx; if(Math.abs(det)<1e-300) return null; const a=(n*sxy-sx*sy)/det, b=(sy-a*sx)/n; if(!(a>0)) return null;
    const ss=pts.reduce((q,[t,p])=>q+Math.pow((p-(a*g(t)+b))/Math.max(Math.abs(p),1e-30),2),0); return {a,b,ss}; };
  const cands=[];
  for(const model of ["pol","hwp"]){ let best=null; for(let t0=-90;t0<90;t0+=0.25){ const m=model==="hwp"?2:1; const r=lin(t=>Math.cos((t-t0)*Math.PI/180*m)**2); if(r&&(!best||r.ss<best.ss)) best=Object.assign(r,{t0,model}); } if(best) cands.push(best); }
  { let best=null; for(let lk=-4;lk<=0;lk+=0.01){ const k=Math.pow(10,lk); for(const sg of [1,-1]){ const r=lin(t=>Math.pow(10,-sg*k*(t-t0r))); if(r&&(!best||r.ss<best.ss)) best=Object.assign(r,{t0:t0r,k:sg*k,model:"nd"}); } } if(best) cands.push(best); }
  if(!cands.length) return null; cands.sort((x,y)=>x.ss-y.ss); const w=cands[0]; w.rms=Math.sqrt(w.ss/n); w.alt=cands.slice(1).map(c=>({model:c.model,rms:Math.sqrt(c.ss/n)})); return w; }
function calUI(el,si,di){
  const m=$q(el,"#dreMap"); const d=S.samples[si].data[di]; const c=Object.assign({},d.cal||{}); const R=d.raw;
  const profs=calProfiles(); const f=PU[S.pu];
  if(!R){ /* 예전 데이터: intensity만 있다 → 배율만 */
    m.innerHTML=`<div class="dre-map"><b>Power calibration — ${esc(d.label||d.ch)}</b><div class="muted">이 데이터는 원본 파워 없이 intensity만 저장돼 있다. intensity에 배율을 곱하는 것만 할 수 있다.</div>
      <div class="row"><label class="dre-f">Multiplier <input type="number" id="cX" step="any" value="1"></label><button class="btn sm primary" id="cOK">Apply</button></div></div>`;
    $q(m,"#cOK").onclick=()=>{ const x=+$q(m,"#cX").value; if(!(x>0)) return; d.raw={kind:"intensity",unit:"W/cm²",P:d.I.slice()}; applyCal(d,{extra:x},true); S.result=null; save(); rerender(); note("배율을 적용했습니다"); }; return; }
  const isP=R.kind==="power";
  m.innerHTML=`<div class="dre-map"><b>Power calibration — ${esc(d.label||d.ch)}</b>
    <div class="muted">원본: ${isP?`파워미터 ${esc(R.unit)} (${sciTxt(Math.min(...R.P.filter(v=>v>0)),3)}–${sciTxt(Math.max(...R.P),3)})`:`intensity ${esc(R.unit)}`}${R.wheel?` · 휠 각도 ${R.wheel.length}점`:""}. 원본은 그대로 두고 intensity만 다시 계산한다.</div>
    ${isP?`<div class="row"><label class="dre-f">Profile <select class="dre-sel" id="cPre"><option value="">(manual)</option>${profs.map((p,i)=>`<option value="${i}" ${c.name===p.name?"selected":""}>${esc(p.name)}${p.builtin?"":" (saved)"}</option>`).join("")}</select></label></div>
    <div class="row"><label class="dre-f" title="파워미터 reading → 대물렌즈 뒤 시료면 파워의 비 (분배비·투과율 포함)">C (sample/PM) <input type="number" id="cC" step="any" value="${c.C||""}"></label>
      <label class="dre-f" title="peak intensity 면적 A_peak = πw²/2 (w = 1/e² 반경). 빔 이미지 가우시안 피팅값">A_peak <input type="number" id="cA" step="any" value="${c.A||""}"> µm²</label>
      <label class="dre-f" title="또는 w(1/e² 반경)로 입력하면 A_peak = πw²/2">w <input type="number" id="cW" step="any" value=""> µm</label>
      <span class="muted">or</span><label class="dre-f">Direct <input type="number" id="cKd" step="any" value="${c.kPerMW||""}"> kW/cm² per mW</label></div>`:""}
    <div class="row"><label class="dre-f" title="예: 파워미터 파장 설정이 1064 가 아니었으면 R(λ_설정)/R(1064) — 광다이오드 응답도 비. 대물렌즈 투과율 등 추가 인자도 여기에">Extra factor <input type="number" id="cX" step="any" value="${c.extra||1}"></label>
      ${isP?`<label class="dre-f">PM wavelength setting <input type="number" id="cLs" step="1" placeholder="1064"> nm → responsivity R(setting)/R(1064) <input type="number" id="cRr" step="any" placeholder="ratio"></label>`:""}</div>
    ${isP&&R.wheel?`<div class="row"><label class="dre-tog"><input type="checkbox" id="cWh" ${c.src==="wheel"?"checked":""}> Reconstruct power from wheel angle instead of PM record (기록이 깨졌을 때)</label></div>
      <div id="cWhBox" ${c.src==="wheel"?"":"hidden"}><div class="muted">P(θ) = a·cos²(θ−θ₀) + b 를 기준점(각도, ${esc(R.unit)})에 맞춘다. 같은 날 다른 파일의 정상 기록이나 손으로 잰 값을 넣는다. 이 데이터에서 믿을 만한 구간만 쓰려면 아래 버튼.</div>
        <textarea id="cAnc" class="dre-ta" placeholder="angle, power (${esc(R.unit)}) — 한 줄에 하나, 3줄 이상">${c.anchors?c.anchors.map(p=>p.join(", ")).join("\n"):""}</textarea>
        <div class="row"><button class="btn xs" id="cAuto">Fill from upper PM range of this data</button><span class="muted" id="cFit"></span></div></div>`:""}
    <div id="cPrev" class="muted"></div>
    <div class="row"><button class="btn sm primary" id="cOne">Apply to this data</button><button class="btn sm" id="cSmp">Whole sample</button><button class="btn sm" id="cAll">All data</button>
      ${isP?`<label class="dre-f">Profile name <input type="text" id="cNm" value="${esc(c.name&&!CALPRE.some(p=>p.name===c.name)?c.name:"")}" placeholder="e.g. 260702 NIR Oil"></label><button class="btn sm" id="cSave">Save as profile</button>`:""}</div></div>`;
  const read=()=>{ const o={extra:+$q(m,"#cX").value||1};
    if(isP){ const C=+$q(m,"#cC").value, A=+$q(m,"#cA").value, w=+$q(m,"#cW").value, kd=+$q(m,"#cKd").value, rr=+(($q(m,"#cRr")||{}).value||0);
      if(kd>0) o.kPerMW=kd; else { o.C=C; o.A=A>0?A:(w>0?Math.PI*w*w/2:A); }
      if(rr>0) o.extra*=rr;
      const pr=$q(m,"#cPre").value; if(pr!==""){ const p=profs[+pr]; if(p&&Math.abs((p.C||0)-(o.C||0))<1e-9&&Math.abs((p.A||0)-(o.A||0))<1e-9) o.name=p.name; }
      const wh=$q(m,"#cWh"); if(wh&&wh.checked){ o.src="wheel"; const pts=($q(m,"#cAnc").value||"").split(/\n/).map(l=>l.trim().split(/[,\s\t]+/).map(Number)).filter(a=>a.length>=2&&isFinite(a[0])&&isFinite(a[1])&&a[1]>0).map(a=>[a[0],a[1]]);
        o.anchors=pts; const fw=fitWheel(pts.map(([t,p])=>[t,p*PWU[R.unit]])); o.wheelFit=fw?{a:fw.a,b:fw.b,t0:fw.t0,k:fw.k,model:fw.model,rms:fw.rms,alt:fw.alt}:null; } }
    return o; };
  const prev=()=>{ const o=read(); const tmp=Object.assign({},d,{cal:o}); const I=calcI(tmp).filter(v=>v>0&&isFinite(v)); const k=calK(o);
    const fw=o.wheelFit; const fe=$q(m,"#cFit"); const MN={pol:"Polarizer cos²(θ−θ₀)",hwp:"Half-wave plate cos²(2(θ−θ₀))",nd:"Continuous ND 10^(−kθ)"}; if(fe) fe.textContent=o.src==="wheel"?(fw?`최적 모델 ${MN[fw.model]}, 상대 rms ${(100*fw.rms).toFixed(1)} %`+(fw.alt&&fw.alt.length?` (다른 모델: ${fw.alt.map(q=>MN[q.model]+" "+(100*q.rms).toFixed(1)+" %").join(", ")})`:"")+(fw.rms>0.05?" — 5 % 넘음, 기준점이나 모델을 확인":""):"기준점 3개 이상 필요"):"";
    $q(m,"#cPrev").innerHTML=k||R.kind==="intensity"?`Intensity range <b>${sciTxt(Math.min(...I)/f,3)} – ${sciTxt(Math.max(...I)/f,3)} ${esc(S.pu)}</b>${k?` · conversion ${(k/1e6).toPrecision(4)} kW/cm² per mW`:""}${Math.max(...I)>1e9?' <span class="gr g4">비정상적으로 큼</span>':Math.min(...I)<1e-2?' <span class="gr g3">매우 작음 — 단위 확인</span>':""}`:'<span class="gr g3">C·A_peak 또는 환산 계수를 넣으세요</span>'; };
  m.querySelectorAll("input,select,textarea").forEach(i=>i.addEventListener("input",prev));
  const pre=$q(m,"#cPre"); if(pre) pre.onchange=()=>{ const p=profs[+pre.value]; if(p){ $q(m,"#cC").value=p.C; $q(m,"#cA").value=p.A; $q(m,"#cW").value=""; $q(m,"#cKd").value=""; } prev(); };
  const wh=$q(m,"#cWh"); if(wh) wh.onchange=()=>{ $q(m,"#cWhBox").hidden=!wh.checked; prev(); };
  const au=$q(m,"#cAuto"); if(au) au.onclick=()=>{ const pts=R.wheel.map((t,i)=>[t,R.P[i]]).filter(p=>p[1]>0); const mx=Math.max(...pts.map(p=>p[1])); $q(m,"#cAnc").value=pts.filter(p=>p[1]>=0.05*mx).map(p=>`${(+p[0]).toFixed(3)}, ${(+p[1]).toPrecision(5)}`).join("\n"); prev(); };
  const go=scope=>{ const o=read(); if(R.kind==="power"&&!calK(o)){ note("환산 계수를 넣으세요","warn"); return; } if(o.src==="wheel"&&!o.wheelFit){ note("휠 기준점이 3개 이상 필요합니다","warn"); return; }
    if(scope!=="one"){ const cnt=S.samples.reduce((a2,s2,a)=>a2+s2.data.filter(q=>q.raw&&q.raw.kind===R.kind&&(scope==="all"||a===si)).length,0); if(!confirm(`${cnt}개 데이터의 보정을 이 설정으로 바꿉니다${o.src==="wheel"?" (휠 재구성은 휠 각도가 있는 데이터에만)":" (다른 데이터에 걸어 둔 휠 재구성도 풀립니다)"}. 계속할까요?`)) return; }
    const targets=[]; S.samples.forEach((s,a)=>s.data.forEach((q,b)=>{ if(scope==="one"?(a===si&&b===di):scope==="smp"?a===si:true) targets.push(q); }));
    let n=0; targets.forEach(q=>{ if(!q.raw) return; if(q.raw.kind!==R.kind) return; const oo=Object.assign({},o); if(oo.src==="wheel"&&!q.raw.wheel) delete oo.src; applyCal(q,oo,true); n++; });
    if(R.kind==="power"){ S.lastCal=Object.assign({},o); delete S.lastCal.anchors; delete S.lastCal.wheelFit; delete S.lastCal.src; }
    S.result=null; save(); rerender(); note(`${n}개 데이터에 보정을 적용했습니다`); };
  $q(m,"#cOne").onclick=()=>go("one"); $q(m,"#cSmp").onclick=()=>go("smp"); $q(m,"#cAll").onclick=()=>go("all");
  const sv=$q(m,"#cSave"); if(sv) sv.onclick=()=>{ const o=read(); const nm=($q(m,"#cNm").value||"").trim(); if(!nm){ note("프로필 이름을 넣으세요","warn"); return; } if(!(o.C>0&&o.A>0)){ note("C 와 A_peak 가 있어야 저장합니다","warn"); return; } saveCalProfile({name:nm,C:o.C,A:o.A,src:"사용자 저장 "+new Date().toLocaleDateString()}); note("프로필을 저장했습니다"); calUI(el,si,di); };
  prev();
}

/* ---------- HDF5 원시 파일 (ScopeFoundry) ---------- */
const H5PRESET=[["NIR · 950 DMSP · A_peak 1.2393 µm² (measured)",268.13],["Blue 450+475 · 600 DMSP · A_peak 0.8715 µm² (provisional · back-calculated)",381.27],["Manual",null]];
async function readH5(file){
  if(typeof h5wasm==="undefined") throw new Error("HDF5 모듈이 없습니다");
  await h5wasm.ready; const buf=new Uint8Array(await file.arrayBuffer());
  try{ h5wasm.FS.mkdir("/work"); }catch(e){}
  const vp="/work/dre_"+Date.now()+".h5"; h5wasm.FS.writeFile(vp,buf);
  const f=new h5wasm.File(vp,"r");
  const get=p=>{ try{ const d=f.get(p); return d&&d.value!=null?Array.from(d.value):null; }catch(e){ return null; } };
  const attrs=p=>{ try{ const g=f.get(p); if(!g||!g.attrs) return {}; const o={}; Object.keys(g.attrs).forEach(k=>{ try{ o[k]=g.attrs[k].value; }catch(e){} }); return o; }catch(e){ return {}; } };
  const out={name:file.name};
  try{
    const meas=(()=>{ try{ return f.get("measurement").keys(); }catch(e){ return []; } })();
    out.stage=attrs("hardware/PI_xyz_stage/settings");
    if(meas.includes("power_scan")){
      const b="measurement/power_scan/"; const P=get(b+"apd_counter_powers"), R=get(b+"apd_counter_apd_count_rates"), W=get(b+"wheel_positions");
      const d=f.get(b+"apd_counter_apd_count_rates_raw"); let raw=null; if(d){ const v=Array.from(d.value), sh=d.shape; raw=[]; for(let i=0;i<sh[0];i++) raw.push(v.slice(i*sh[1],(i+1)*sh[1])); }
      if(!P||!R) throw new Error("power_scan 에 apd_counter_powers / apd_counter_apd_count_rates 가 없습니다");
      const st=attrs(b+"settings"); let turn=P.length-1; if(W) for(let i=1;i<W.length;i++){ if(W[i]<W[i-1]-1e-9){ turn=i-1; break; } }
      Object.assign(out,{kind:"scan",W,P:P.map(v=>v*1e3),R,raw,turn,sweep:st.power_wheel_sweep_type||"",reps:st.repetitions});
    } else if(meas.includes("ni_daq_pulse_resp")){
      const b="measurement/ni_daq_pulse_resp/"; const t=get(b+"times"), c=get(b+"counts"), v=get(b+"set_voltages");
      if(!t||!c) throw new Error("ni_daq_pulse_resp 에 times / counts 가 없습니다");
      const st=attrs(b+"settings"); const pm=attrs("hardware/power_meter_2/settings");
      let on=-1,off=-1; if(v){ on=v.findIndex(x=>x>0); for(let i=v.length-1;i>=0;i--) if(v[i]>0){ off=i; break; } }
      const med=[...c].sort((a,b)=>a-b)[c.length>>1]; const bad0=c[0]>100*Math.max(med,1);
      Object.assign(out,{kind:"trace",t,c,on,off,bad0,fs:st.sampling_frequency,tHigh:st.t_high,tAcq:st.t_acquisition,PmW:pm.power!=null?pm.power*1e3:null});
    } else throw new Error("지원하는 측정이 아닙니다 (measurement/"+meas.join(", ")+")");
  } finally { try{ f.close(); }catch(e){} try{ h5wasm.FS.unlink(vp); }catch(e){} }
  return out;
}
function h5UI(el,chs,H){
  const m=$q(el,"#dreMap"); const pre=H5PRESET.map(([n,v],i)=>`<option value="${i}">${esc(n)}${v?` — ${v} kW/cm² per mW`:""}</option>`).join("");
  const stg=H.stage&&H.stage.x_position!=null?` · stage (${(+H.stage.x_position).toFixed(1)}, ${(+H.stage.y_position).toFixed(1)}) µm`:"";
  const conv=`<label class="dre-f">Power conversion <select class="dre-sel" id="hPre">${pre}</select></label><label class="dre-f">kW/cm² per mW <input type="number" id="hK" step="any" value="${H5PRESET[0][1]}"></label>`;
  if(H.kind==="scan"){
    const n=H.P.length, nu=H.turn+1;
    /* 파워 비례 배경 추정: up branch 아래쪽 30 % (가장 큰 점프 이전) 직선 */
    const lr=H.R.slice(0,nu).map((r,i)=>i?Math.abs(Math.log(Math.max(r,1)/Math.max(H.R[i-1],1))):0); let jk=nu; lr.forEach((v,i)=>{ if(v>Math.log(3)&&i<jk) jk=i; });
    const idx=[...Array(Math.max(3,Math.floor(jk*0.9))).keys()].filter(i=>i<jk);
    let a=0,b=0,sd=0; if(idx.length>=3){ const x=idx.map(i=>H.P[i]), y=idx.map(i=>H.R[i]); const mx=x.reduce((p,q)=>p+q,0)/x.length, my=y.reduce((p,q)=>p+q,0)/y.length; let sxy=0,sxx=0; x.forEach((v,i)=>{ sxy+=(v-mx)*(y[i]-my); sxx+=(v-mx)**2; }); b=sxx>0?sxy/sxx:0; a=my-b*mx; sd=Math.sqrt(y.reduce((p,q,i)=>p+(q-a-b*x[i])**2,0)/Math.max(1,y.length-2)); }
    const cvHi=H.raw?H.raw.filter(r=>{ const m=r.reduce((p,q)=>p+q,0)/r.length; const s=Math.sqrt(r.reduce((p,q)=>p+(q-m)**2,0)/Math.max(1,r.length-1)); return m>0&&s/m>0.3; }).length:0;
    m.innerHTML=`<div class="dre-map"><div><b>${esc(H.name)}</b> <span class="muted">— power scan ${n} pts (${esc(H.sweep||"")}${H.turn<n-1?`, up ${nu} · down ${n-nu}`:""}, repeats ${H.reps||"?"})${stg}</span></div>
      <div class="row">${conv}</div>
      <div class="row"><label class="dre-f">Channel <select class="dre-sel" id="hCh">${chs.map(c=>`<option>${esc(c)}</option>`).join("")}</select></label>
        <label class="dre-f">Background = <input type="number" id="hA" step="any" value="${a.toFixed(1)}"> + <input type="number" id="hB" step="any" value="${b.toFixed(1)}"> × P[mW] cps</label>
        <label class="dre-f">Background σ <input type="number" id="hS" step="any" value="${Math.max(sd,1).toFixed(1)}"></label><label class="dre-f">multiplicative noise <input type="number" id="hRel" step="0.01" value="0.05"></label></div>
      <div class="row"><label class="dre-tog"><input type="checkbox" id="hCv" ${cvHi?"checked":""}> Exclude points that switched within repeats (repeat CV &gt; 0.3: ${cvHi} pts)</label></div>
      <div class="muted">배경은 up branch의 첫 큰 점프 이전 ${idx.length}점에 맞춘 직선이다(파워에 비례하는 누설·산란 포함). 뺀 값이 0 이하이면 censored 점(상한)으로 들어가고 버리지 않는다.</div>
      <div class="row"><button class="btn sm primary" id="hOK">Add up/down S-curve</button></div></div>`;
  } else {
    const dt=H.t.length>1?H.t[1]-H.t[0]:0; const pre0=H.on>1?H.c.slice(H.bad0?1:0,H.on):[]; const bg=pre0.length?pre0.reduce((p,q)=>p+q,0)/pre0.length:0;
    m.innerHTML=`<div class="dre-map"><div><b>${esc(H.name)}</b> <span class="muted">— pulse response ${H.t.length} pts, ${dt?(dt*1e6).toFixed(0)+" µs spacing":""}, on ${H.on>=0?(H.t[H.on]*1e3).toFixed(2)+" ms":"?"} · off ${H.off>=0?(H.t[H.off]*1e3).toFixed(2)+" ms":"?"}${H.tAcq?` · accumulated ${H.tAcq} s`:""}${stg}</span></div>
      ${H.bad0?'<div class="note warn">첫 샘플이 나머지의 100배 이상 — 카운터 누적값으로 보고 버립니다.</div>':""}
      <div class="row">${conv}<label class="dre-f">Power <input type="number" id="hP" step="any" value="${H.PmW!=null?H.PmW.toPrecision(4):""}"> mW</label></div>
      <div class="note warn">파일의 파워는 측정 시점 파워미터 값이다. 펄스 중에 읽었다면 듀티 평균일 수 있으니 on 구간 peak 파워와 같은지 확인하세요.</div>
      <div class="row"><label class="dre-f">Channel <select class="dre-sel" id="hCh">${chs.map(c=>`<option>${esc(c)}</option>`).join("")}</select></label><label class="dre-f">Background (pre-on mean) <input type="number" id="hBg" step="any" value="${bg.toPrecision(3)}"></label></div>
      <div class="row"><button class="btn sm primary" id="hOK">Add time trace (on-window, diagnostic only)</button></div></div>`;
  }
  const pr=$q(m,"#hPre"), kk=$q(m,"#hK"); pr.onchange=()=>{ const v=H5PRESET[+pr.value][1]; if(v!=null) kk.value=v; };
  $q(m,"#hOK").onclick=()=>{
    const si=+$q(el,"#dreDSmp").value, ch=$q(m,"#hCh").value, K=+kk.value*1e3; /* W/cm² per mW */
    if(!(K>0)){ note("환산 계수를 넣으세요","warn"); return; }
    const smp=S.samples[si];
    if(H.kind==="scan"){ const A=+$q(m,"#hA").value||0, Bb=+$q(m,"#hB").value||0, sbg=+$q(m,"#hS").value||0, rel=+$q(m,"#hRel").value||0.05, cv=$q(m,"#hCv").checked;
      const skip=i=>{ if(!cv||!H.raw) return false; const r=H.raw[i]; const mm=r.reduce((p,q)=>p+q,0)/r.length; const s=Math.sqrt(r.reduce((p,q)=>p+(q-mm)**2,0)/Math.max(1,r.length-1)); return mm>0&&s/mm>0.3; };
      const pv=H5PRESET[+$q(m,"#hPre").value]; const cal=(pv&&pv[1]!=null&&Math.abs(+kk.value-pv[1])<1e-9)?(()=>{ const q=CALPRE[+$q(m,"#hPre").value]; return {name:q.name,C:q.C,A:q.A}; })():{kPerMW:+kk.value};
      const mk=(a,b,dir)=>{ const I=[],y=[],Pr=[],Wh=[]; for(let i=a;i<=b;i++){ if(skip(i)||!(H.P[i]>0)) continue; I.push(H.P[i]*K); Pr.push(H.P[i]); if(H.W) Wh.push(H.W[i]); y.push(H.R[i]-(A+Bb*H.P[i])); }
        const o={type:"scurve",ch,label:H.name+(dir<0?" down":" up"),I,y,sbg,k:3,relErr:rel,dir,raw:{kind:"power",unit:"mW",P:Pr,wheel:H.W?Wh:undefined}}; applyCal(o,cal,true); return o; };
      smp.data.push(mk(0,H.turn,1)); if(H.turn<H.P.length-1) smp.data.push(mk(H.turn+1,H.P.length-1,-1));
      note("S-curve 를 추가했습니다"+(H.turn<H.P.length-1?" (up·down)":""));
    } else { const P=+$q(m,"#hP").value; if(!(P>0)){ note("파워(mW)를 넣으세요","warn"); return; } const bg=+$q(m,"#hBg").value||0;
      const a=Math.max(H.on,H.bad0?1:0), b=H.off>=0?H.off:H.t.length-1; const t=[],y=[]; for(let i=a;i<=b;i++){ t.push(H.t[i]-H.t[a]); y.push(H.c[i]-bg); }
      const pv=H5PRESET[+$q(m,"#hPre").value]; const cal=(pv&&pv[1]!=null&&Math.abs(+kk.value-pv[1])<1e-9)?{name:CALPRE[+$q(m,"#hPre").value].name,C:CALPRE[+$q(m,"#hPre").value].C,A:CALPRE[+$q(m,"#hPre").value].A}:{kPerMW:+kk.value};
      const o={type:"trace",ch,label:H.name,I:[P*K],t,y,raw:{kind:"power",unit:"mW",P:[P]}}; applyCal(o,cal,true); smp.data.push(o); note("시간 트레이스를 추가했습니다 — 적합성 A6 에 쓰입니다"); }
    PENDING=null; S.result=null; save(); rerender();
  };
}

/* ---------- 데이터 점검 ----------
   모델에 넣기 전에 잡을 수 있는 입력 실수를 데이터마다 검사한다. level: err(피팅 막음) · warn(표시) · info.
   실제로 겪은 사례(2026-09-30): up·down 을 한 데이터셋에 넣음, pulse CSV 첫 행의 카운터 누적값이 rise 시간으로 들어감,
   배경을 빼지 않은 채 배경σ=0, 보간 행 포함. */
function turnIdx(I){ const n=I.length; if(n<8) return -1; let k=0; for(let i=1;i<n;i++) if(I[i]>I[k]) k=i;
  if(k<3||k>n-4) return -1; let ub=0,db=0; for(let i=1;i<=k;i++) if(I[i]<I[i-1]*0.98) ub++; for(let i=k+1;i<n;i++) if(I[i]>I[i-1]*1.02) db++;
  if(ub>0.1*k||db>0.1*(n-k)) return -1; if(Math.min(...I.slice(k+1))>0.7*I[k]) return -1; return k; }
/* 배경 추정: 파워 순으로 낮은 쪽 점들의 중앙값·MAD 로 바닥을 잡고, 바닥 안쪽 점에 (상수 + 파워 비례) 직선을 맞춘다 */
function bgEstimate(I,y){ const idx=[...I.keys()].filter(i=>isFinite(y[i])&&I[i]>0).sort((a,b)=>I[a]-I[b]); if(idx.length<8) return null;
  const m0=Math.max(5,Math.floor(idx.length*0.15)), base=idx.slice(0,m0).map(i=>y[i]).sort((a,b)=>a-b); const med=base[base.length>>1];
  const mad=[...base].map(v=>Math.abs(v-med)).sort((a,b)=>a-b)[base.length>>1]*1.4826||Math.max(1,Math.sqrt(Math.abs(med)));
  const use=[]; for(const i of idx){ if(y[i]>med+6*mad&&use.length>=m0) break; use.push(i); } if(use.length<5) return null;
  const x=use.map(i=>I[i]), v=use.map(i=>y[i]); const mx=x.reduce((p,q)=>p+q,0)/x.length, my=v.reduce((p,q)=>p+q,0)/v.length; let sxy=0,sxx=0; x.forEach((q,k)=>{ sxy+=(q-mx)*(v[k]-my); sxx+=(q-mx)**2; });
  let b=sxx>0?sxy/sxx:0; if(b<0) b=0; const a=my-b*mx; const sd=Math.sqrt(v.reduce((p,q,k)=>p+(q-a-b*x[k])**2,0)/Math.max(1,v.length-2));
  return {a,b,mean:my,sd:Math.max(sd,1e-9),n:use.length,Imax:x[x.length-1]}; }
function dataIssues(d){ const out=[];
  if(d.type==="scurve"){ const k=turnIdx(d.I); if(k>=0&&!d.split) out.push({lvl:"err",code:"mix",txt:`파워가 ${fmtIu(d.I[k])} 까지 올랐다가 다시 내려간다 — up(${k+1}점)·down(${d.I.length-k-1}점)이 한 데이터셋에 섞여 있다. 정상상태 모델은 한 방향 branch 만 설명하므로 나눠야 한다.`,fix:"split",fixTxt:"Split up/down"});
    const e=d.bgSub?null:bgEstimate(d.I,d.y);
    if(e&&e.n>=5&&e.mean>e.sd&&e.mean/(e.sd/Math.sqrt(e.n))>5) out.push({lvl:"warn",code:"bg",txt:`낮은 파워 ${e.n}점의 바닥 평균이 ${sciTxt(e.mean,3)} (점 간 SD ${sciTxt(e.sd,2)}) 로 0 이 아니다 — 배경을 빼지 않은 것으로 보인다${!(d.sbg>0)&&!d.sy?`. 배경σ 도 0 이라 바닥 근처 점의 σ 가 multiplicative noise(${Math.round((d.relErr||0.05)*100)} %)만으로 과소평가된다`:""}. 모델은 배경까지 신호로 설명하려 한다.`,fix:"bg",fixTxt:"Estimate &amp; subtract background"});
    if(d.bgSub) out.push({lvl:"info",code:"bgdone",txt:`배경 ${sciTxt(d.bgSub.a,3)} + ${sciTxt(d.bgSub.b,2)}·I 를 뺐다 (바닥 ${d.bgSub.n}점, σ ${sciTxt(d.sbg,2)}).`,fix:"bgundo",fixTxt:"Undo"});
    if(d.nInterp) out.push({lvl:"info",code:"interp",txt:`보간 행 ${d.nInterp}개를 뺐다.`}); }
  if(d.type==="rise"||d.type==="decay"){ const bad=d.t.map((t,i)=>i).filter(i=>!(d.t[i]>1e-9&&d.t[i]<1e4));
    if(bad.length) out.push({lvl:"err",code:"t",txt:`시간 값 ${bad.length}개가 1 ns–10⁴ s 밖이다 (예: ${sciTxt(d.t[bad[0]],3)} s). 카운터 누적값이나 단위가 잘못 들어온 것으로 보인다 — pulse 응답 파일은 '시간 트레이스'로 넣는다.`,fix:"del",fixTxt:"Delete this data"}); }
  return out; }
function fmtIu(W){ return sciTxt(W/PU[S.pu],3)+" "+S.pu; }
function issuesRow(d,si,di){ const L=dataIssues(d); if(!L.length) return ""; return `<tr class="dre-iss"><td colspan="6">${L.map(q=>`<div class="note ${q.lvl==="err"?"warn":q.lvl==="warn"?"warn":""}" style="margin:2px 0;padding:6px 8px">${q.lvl==="err"?"<b>Fix before fitting</b> — ":q.lvl==="warn"?"<b>Check</b> — ":""}${esc(q.txt)} ${q.fix?`<button class="btn xs" data-fix="${si}|${di}|${q.fix}">${esc(q.fixTxt)}</button>`:""}</div>`).join("")}</td></tr>`; }
const sliceD=(d,a,b,dir)=>{ const o=JSON.parse(JSON.stringify(d)); const cut=v=>Array.isArray(v)?v.slice(a,b):v; o.I=cut(d.I); o.y=cut(d.y); if(d.sy) o.sy=cut(d.sy); if(d.yRaw) o.yRaw=cut(d.yRaw);
  if(o.raw){ o.raw.P=cut(d.raw.P); if(d.raw.wheel) o.raw.wheel=cut(d.raw.wheel); } o.dir=dir; o.split=true; o.label=(d.label||d.ch)+(dir<0?" down":" up"); delete o._fit; return o; };
function splitUpDown(d){ const k=turnIdx(d.I); if(k<0) return null; return [sliceD(d,0,k+1,1),sliceD(d,k+1,d.I.length,-1)]; }
function subtractBg(d,e){ if(!e) return false; if(!d.yRaw) d.yRaw=d.y.slice(); d.y=d.yRaw.map((v,i)=>v-(e.a+e.b*d.I[i])); d.bgSub={a:e.a,b:e.b,unit:"W/cm²",n:e.n}; d.sbg=e.sd; return true; }
function fixData(el,si,di,act){ const smp=S.samples[si], d=smp.data[di]; if(!d) return;
  if(act==="split"){ const r=splitUpDown(d); if(!r){ note("전환점을 찾지 못했습니다","warn"); return; } smp.data.splice(di,1,r[0],r[1]); note(`up ${r[0].I.length}점 · down ${r[1].I.length}점으로 나눴습니다`); }
  else if(act==="bg"){ /* up branch(또는 전체)의 바닥으로 추정해 같은 시료의 같은 원본에서 온 데이터에 함께 쓴다 */
    const e=bgEstimate(d.I,d.yRaw||d.y); if(!subtractBg(d,e)){ note("배경을 추정하지 못했습니다","warn"); return; }
    note(`배경 ${sciTxt(e.a,3)} + ${sciTxt(e.b,2)}·I 를 뺐습니다 (σ ${sciTxt(e.sd,2)}, ${e.n}점). 0 이하가 된 점은 censored 로 들어갑니다.`); }
  else if(act==="del"){ smp.data.splice(di,1); note("삭제했습니다"); }
  else if(act==="bgundo"){ if(d.yRaw){ d.y=d.yRaw.slice(); delete d.yRaw; delete d.bgSub; d.sbg=0; note("배경 빼기를 되돌렸습니다"); } }
  S.result=null; save(); rerender(); }
function allIssues(){ const L=[]; S.samples.forEach((s,si)=>s.data.forEach((d,di)=>dataIssues(d).forEach(q=>L.push(Object.assign({si,di,sid:s.id,ch:d.ch},q))))); return L; }

/* ---------- pulse 응답 CSV (Time, Counts, SetV, Power) ----------
   ScopeFoundry ni_daq_pulse_resp 를 CSV 로 내보낸 형식. 첫 행 Counts 는 카운터 누적값, Power 는 첫 행에만 있다. */
/* 가져오기 형식별 마지막 단위 선택 (이 브라우저) */
function unitMem(kind,v){ const K="innoanalytics_dre_unit_v1"; try{ const o=JSON.parse(localStorage.getItem(K)||"{}"); if(v===undefined) return o[kind]||null; o[kind]=v; localStorage.setItem(K,JSON.stringify(o)); }catch(e){} return null; }
/* 파일 이름의 파워 표기(예: 8uW, 1.5mW)와 파일 안 값이 다르면 알린다 — 둘 중 하나는 보정된 값일 가능성 */
function fileUnitHint(name,v){ const q=/(\d+(?:\.\d+)?)\s*(uW|µW|mW|nW|W)(?![a-z])/i.exec(name||""); if(!q||!isFinite(v)) return "";
  return ` 파일 이름에는 <b>${esc(q[1]+" "+q[2].replace(/^u/i,"µ"))}</b>, 파일 안 첫 값은 <b>${sciTxt(v,3)}</b> — 이름은 측정 설정값, 파일 값은 보정된 값일 수 있다.`; }
function isPulseTable(T){ const h=T.header.map(x=>String(x).toLowerCase()); return h.some(x=>/^time/.test(x))&&h.some(x=>/setv|set_v|set.?volt/.test(x)); }
function pulseStats(t,c,on,off,bg){ const seg=[]; for(let i=on;i<=off;i++) seg.push(c[i]-bg); const n=seg.length; if(n<20) return null; const dt=(t[off]-t[on])/Math.max(1,n-1);
  const k=Math.max(1,Math.round(n/300)); const sm=seg.map((_,i)=>{ let s2=0,q=0; for(let j=Math.max(0,i-k);j<=Math.min(n-1,i+k);j++){ s2+=seg[j]; q++; } return s2/q; });
  const tail=Math.max(3,Math.floor(n*0.05)); const fin=sm.slice(-tail).reduce((a,b)=>a+b,0)/tail; const tf=f=>{ const i=sm.findIndex(v=>v>=f*fin); return i<0?null:i*dt; };
  return {fin,t10:tf(0.1),t50:tf(0.5),t90:tf(0.9),dt,dur:(n-1)*dt}; }
function binTrace(t,y,maxN){ const n=t.length; if(n<=maxN) return {t,y}; const k=Math.ceil(n/maxN), T=[],Y=[]; for(let i=0;i<n;i+=k){ let st=0,sy=0,q=0; for(let j=i;j<Math.min(n,i+k);j++){ st+=t[j]; sy+=y[j]; q++; } T.push(st/q); Y.push(sy/q); } return {t:T,y:Y}; }
function pulseCsvUI(el,chs,T){ const m=$q(el,"#dreMap"); const h=T.header.map(x=>String(x).toLowerCase());
  const jt=h.findIndex(x=>/^time/.test(x)), jv=h.findIndex(x=>/setv|set_v|set.?volt/.test(x)), jc=h.findIndex(x=>/count|signal|apd/.test(x)), jp=h.findIndex(x=>/power|intens/.test(x));
  const t=T.cols[jt], c=T.cols[jc>=0?jc:1], v=T.cols[jv]; const med=[...c.slice(1)].filter(isFinite).sort((a,b)=>a-b); const md=med[med.length>>1]||0;
  const bad0=c[0]>100*Math.max(md,1); const s0=bad0?1:0; let on=-1,off=-1; for(let i=s0+1;i<v.length;i++) if(v[i]>0.5&&!(v[i-1]>0.5)){ on=i; break; } for(let i=v.length-1;i>on;i--) if(v[i]>0.5){ off=i; break; }
  const P0=jp>=0?T.cols[jp].find(x=>isFinite(x)):NaN; const pre=c.slice(s0,Math.max(s0,on)); const bg=pre.length?pre.reduce((a,b)=>a+b,0)/pre.length:0;
  const st=on>0&&off>on?pulseStats(t,c,on,off,bg):null; const post=off>0?c.slice(off+1):[];
  let tOffHalf=null; if(st&&post.length){ const i=post.findIndex(x=>x-bg<=0.5*st.fin); if(i>=0) tOffHalf=(i+1)*(st.dt||0); }
  const nm=T.name||"pulse CSV"; const lastU=unitMem("pulse")||("I|"+S.pu);
  const uopt=`<optgroup label="Intensity (calibrated · sample plane)">${Object.keys(PU).map(u=>`<option value="I|${u}" ${"I|"+u===lastU?"selected":""}>${u}</option>`).join("")}</optgroup><optgroup label="Power-meter power (needs calibration)">${Object.keys(PWU).map(u=>`<option value="P|${u}" ${"P|"+u===lastU?"selected":""}>${u} (PM reading)</option>`).join("")}</optgroup>`;
  const hint=fileUnitHint(nm,P0);
  const slow=st&&st.t90!=null&&st.t90>10;
  m.innerHTML=`<div class="dre-map"><div><b>${esc(nm)}</b> <span class="muted">— pulse-response CSV ${T.n} rows · Δt ${st?sciTxt(st.dt,2)+" s":"?"} · on ${on>=0?sciTxt(t[on],3)+" s":"?"} · off ${off>=0?sciTxt(t[off],3)+" s":"?"}</span></div>
    ${bad0?`<div class="note warn">첫 행 Counts(${sciTxt(c[0],3)})는 나머지의 100배 이상 — 카운터 누적값으로 보고 버립니다. 이 값을 rise 시간으로 쓰면 안 됩니다.</div>`:""}
    ${st?`<div class="muted">Background (pre-on) ${sciTxt(bg,3)} · plateau ${sciTxt(st.fin,3)} (per bin) · rise t10 ${st.t10!=null?sciTxt(st.t10,2):"—"} s, t50 ${st.t50!=null?sciTxt(st.t50,2):"—"} s, t90 ${st.t90!=null?sciTxt(st.t90,2):"—"} s · time to half after off ${tOffHalf!=null?sciTxt(tOffHalf,2)+" s":"—"}</div>`:""}
    ${slow?`<div class="note warn">rise 가 ${sciTxt(st.t90,2)} s 에 걸쳐 오르고 끈 뒤에는 ${tOffHalf!=null?sciTxt(tOffHalf,2)+" s 안에":"빠르게"} 떨어진다 — 켤 때와 끌 때가 비대칭이다. 준위 수명(ms)으로 정해지는 정상상태 DRE 밖의 느린 변수일 수 있으므로 rise 제약으로 쓰지 말고 진단용 트레이스로만 넣는다.</div>`:""}
    <div class="note" style="margin:6px 0">Power 열이 <b>보정이 끝난 시료면 intensity</b>인지 <b>파워미터 값</b>인지 데이터마다 다르다 — 단위에서 고른다 (마지막 선택을 기억한다).${hint}</div>
    <div class="row"><label class="dre-f">Excitation power <input type="number" id="pP" step="any" value="${isFinite(P0)?P0:""}"></label><label class="dre-f">Unit <select class="dre-sel" id="pU">${uopt}</select></label>
      <label class="dre-f">Channel <select class="dre-sel" id="pCh">${chs.map(q=>`<option>${esc(q)}</option>`).join("")}</select></label><label class="dre-f">Background <input type="number" id="pBg" step="any" value="${bg.toPrecision(3)}"></label></div>
    <div class="row"><button class="btn sm primary" id="pOK">Add time trace (on-window, diagnostic only)</button><span class="muted">${off>on&&on>=0?`on 구간 ${off-on+1}점 → 최대 2000점으로 평균 묶음`:"SetV 로 on/off 를 찾지 못했습니다"}</span></div></div>`;
  $q(m,"#pOK").onclick=()=>{ const si=+$q(el,"#dreDSmp").value, P=+$q(m,"#pP").value, [uk,uu]=$q(m,"#pU").value.split("|"), b=+$q(m,"#pBg").value||0; unitMem("pulse",uk+"|"+uu);
    if(!(P>0)){ note("여기 파워를 넣으세요","warn"); return; } if(!(on>=0&&off>on)){ note("on/off 구간이 없습니다","warn"); return; }
    const tt=[],yy=[]; for(let i=on;i<=off;i++){ if(!isFinite(t[i])||!isFinite(c[i])) continue; tt.push(t[i]-t[on]); yy.push(c[i]-b); } const B=binTrace(tt,yy,2000);
    const o={type:"trace",ch:$q(m,"#pCh").value,label:nm,I:[uk==="I"?P*PU[uu]:P],t:B.t,y:B.y,raw:{kind:uk==="I"?"intensity":"power",unit:uu,P:[P]}};
    if(uk==="I"){ o.I=[P*PU[uu]]; o.cal={}; o.calOK=true; } else applyCal(o,S.lastCal||{name:CALPRE[0].name,C:CALPRE[0].C,A:CALPRE[0].A},false);
    if(st) o.stats={t50:st.t50,t90:st.t90,offHalf:tOffHalf,plateau:st.fin,bg};
    S.samples[si].data.push(o); PENDING=null; S.result=null; save(); rerender(); note("시간 트레이스를 추가했습니다 — 적합성 A6 에 쓰입니다"); };
}
function mapUI(el,chs){
  const T=PENDING, m=$q(el,"#dreMap"); const type=$q(el,"#dreDType").value;
  const opt=T.header.map((h,j)=>`<option value="${j}">${esc(h)}</option>`).join("");
  if(type==="trace"){
    m.innerHTML=`<div class="dre-map"><div class="muted">레이저를 켠 순간부터의 신호(예: 50 kHz 계수). 피팅에는 쓰지 않고, 모델이 같은 파워에서 같은 모양(점등 여부·과도 상승)을 내는지만 검사한다.</div>
      <div class="row"><label class="dre-f">Time column <select class="dre-sel" id="mI">${opt}</select></label><label class="dre-f">Time unit <select class="dre-sel" id="mTU"><option value="1">s</option><option value="1e-3">ms</option><option value="1e-6">µs</option></select></label>
      <label class="dre-f">Signal column <select class="dre-sel" id="mY">${T.header.map((h,j)=>`<option value="${j}" ${j===1?"selected":""}>${esc(h)}</option>`).join("")}</select></label>
      <label class="dre-f">Channel <select class="dre-sel" id="mCh">${chs.map(c=>`<option>${esc(c)}</option>`).join("")}</select></label></div>
      <div class="row"><label class="dre-f">Excitation power <input type="number" id="mTP" step="any"> <select class="dre-sel" id="mPU">${Object.keys(PU).map(u=>`<option ${u===S.pu?"selected":""}>${u}</option>`).join("")}</select></label>
      <label class="dre-f">On time <input type="number" id="mT0" step="any" value="0"></label><label class="dre-f">Off time <input type="number" id="mT1" step="any" placeholder="to end"></label>
      <label class="dre-f">Background (dark) <input type="number" id="mBg" step="any" value="0"></label></div>
      <div class="row"><button class="btn sm primary" id="mOK">Add</button><span class="muted">${T.n} rows · 켠 시각 이전과 끈 시각 이후는 버림</span></div></div>`;
    $q(m,"#mOK").onclick=()=>{
      const si=+$q(el,"#dreDSmp").value, jt=+$q(m,"#mI").value, jy=+$q(m,"#mY").value, tu=+$q(m,"#mTU").value, P0=+$q(m,"#mTP").value*PU[$q(m,"#mPU").value];
      const t0=+$q(m,"#mT0").value||0, t1v=$q(m,"#mT1").value, t1=t1v===""?Infinity:+t1v, bg=+$q(m,"#mBg").value||0;
      if(!(P0>0)){ note("여기 파워를 넣으세요","warn"); return; } if(jt===jy){ note("시간 열과 신호 열이 같습니다","warn"); return; }
      const t=[],y=[]; T.cols[jt].forEach((tv,r)=>{ const v=T.cols[jy][r]; if(!isFinite(tv)||!isFinite(v)) return; if(tv<t0||tv>t1) return; t.push((tv-t0)*tu); y.push(v-bg); });
      if(t.length<20){ note("점이 너무 적습니다 (20 개 이상)","warn"); return; }
      const pu=$q(m,"#mPU").value; S.samples[si].data.push({type:"trace",ch:$q(m,"#mCh").value,label:T.header[jy],I:[P0],t,y,raw:{kind:"intensity",unit:pu,P:[+$q(m,"#mTP").value]},cal:{},calOK:true});
      PENDING=null; save(); rerender(); note("시간 트레이스를 추가했습니다 — 적합성 A6 에 쓰입니다"); };
    return; }
  const H=T.header.map(x=>String(x)); const jSig=H.findIndex(x=>/(^|_|\b)(std|sd|sigma|σ|err|error|unc)/i.test(x)); const jInt=H.findIndex(x=>/interp/i.test(x));
  const jPow=Math.max(0,H.findIndex(x=>/power|intens|^p$|^p[_ ]/i.test(x))); const repN=(()=>{ const q=/rep(?:eat)?[_ ]?(\d+)/i.exec(T.name||""); return q?+q[1]:1; })();
  const interpRows=jInt>=0?T.txt[jInt].map(v=>/^(true|1|yes)$/i.test(v.trim())):null; const nInt=interpRows?interpRows.filter(Boolean).length:0;
  const kTurn=type==="scurve"?turnIdx(T.cols[jPow].filter((p,r)=>isFinite(p)&&p>0&&!(interpRows&&interpRows[r]))):-1;
  m.innerHTML=`<div class="dre-map">${T.name?`<div><b>${esc(T.name)}</b> <span class="muted">— ${T.n} rows</span></div>`:""}<div class="row"><label class="dre-f">Power column <select class="dre-sel" id="mI">${T.header.map((h,j)=>`<option value="${j}" ${j===jPow?"selected":""}>${esc(h)}</option>`).join("")}</select></label>
    <label class="dre-f">Unit <select class="dre-sel" id="mPU">${(()=>{ const lu=unitMem("table")||"P|mW"; return `<optgroup label="Intensity (calibrated · sample plane)">${Object.keys(PU).map(u=>`<option value="I|${u}" ${"I|"+u===lu?"selected":""}>${u}</option>`).join("")}</optgroup><optgroup label="Power-meter power (needs calibration)">${Object.keys(PWU).map(u=>`<option value="P|${u}" ${"P|"+u===lu?"selected":""}>${u} (PM reading)</option>`).join("")}</optgroup>`; })()}</select></label></div>
    <div class="muted">파워 열이 보정이 끝난 시료면 intensity 인지 파워미터 값인지는 데이터마다 다르다 — 확인해서 고른다 (마지막 선택을 기억한다).${fileUnitHint(T.name||"",T.cols[jPow].find(v=>isFinite(v)))}</div>
    <div class="muted">파워미터 값을 넣으면 원본을 그대로 두고, 마지막으로 쓴 보정(없으면 NIR 950 DMSP·Oil)을 임시로 걸어 둔다. 목록의 ‘파워 보정’에서 확인·수정한다.</div>
    ${T.header.map((h,j)=>j===jPow||(type==="scurve"&&(j===jSig||j===jInt))?"":`<div class="row"><label class="dre-f">Column “${esc(h)}” → channel <select class="dre-sel" data-col="${j}"><option value="">(unused)</option>${chs.map(c=>`<option>${esc(c)}</option>`).join("")}</select></label></div>`).join("")}
    ${type==="scurve"?`<div class="row"><label class="dre-f" title="점마다 다른 σ. 반복 측정 SD 라면 반복 수 n 으로 나눠 평균의 σ 로 쓴다">Per-point σ column <select class="dre-sel" id="mSy"><option value="-1">(none)</option>${T.header.map((h,j)=>`<option value="${j}" ${j===jSig?"selected":""}>${esc(h)}</option>`).join("")}</select></label>
      <label class="dre-f">Repeats n <input type="number" id="mRep" min="1" step="1" value="${repN}" title="파일 이름의 rep 숫자에서 읽었다"></label>
      ${jInt>=0?`<label class="dre-tog"><input type="checkbox" id="mInt" checked> Exclude interpolated rows (${esc(H[jInt])} = True: ${nInt} rows)</label>`:""}</div>
      <div class="row"><label class="dre-tog"><input type="checkbox" id="mSplit" ${kTurn>=0?"checked":""} ${kTurn>=0?"":"disabled"}> Auto-split up/down ${kTurn>=0?`(파워 최대 지점에서 전환: up ${kTurn+1} · down 나머지)`:"(전환 없음)"}</label>
      <label class="dre-tog" title="낮은 파워 바닥(중앙값 ± 6·MAD 안쪽 점)에 상수 + 파워 비례 직선을 맞춰 빼고, 그 산포를 배경σ로 쓴다. 0 이하가 된 점은 censored"><input type="checkbox" id="mBgAuto" checked> Auto-estimate &amp; subtract background</label></div>
      <div class="row"><label class="dre-f">Background (dark) σ <input type="number" id="mSbg" step="any" value="0" title="배경을 뺀 신호의 잡음 표준편차. 이 값의 k배 아래 점은 버리지 않고 '상한'으로 쓴다(censored). 배경 자동 추정을 켜면 추정값으로 바뀐다"></label>
      <label class="dre-f">k <input type="number" id="mK" step="0.5" value="3"></label><label class="dre-f">multiplicative noise <input type="number" id="mRel" step="0.01" value="0.05"></label>
      <label class="dre-f">Direction <select class="dre-sel" id="mDir"><option value="1">up</option><option value="-1">down</option></select></label></div>
      <div class="muted">신호는 <b>배경을 뺀 값</b>이어야 합니다(자동 추정을 켜면 앱이 뺍니다). 0 이나 음수는 censored 점(상한)으로 처리되고 버리지 않습니다 — 0 으로 바꿔 log 에 넣지 않습니다.</div>`:
      `<div class="row"><label class="dre-f">Relative error <input type="number" id="mRel" step="0.01" value="0.05"></label>${type==="rise"?`<label class="dre-f">Reach fraction <input type="number" id="mFrac" step="0.05" value="0.95"></label>`:""}</div>`}
    <div class="row"><button class="btn sm primary" id="mOK">Add</button><span class="muted">${T.n} rows</span></div></div>`;
  const firstCh=m.querySelector('select[data-col="1"]'); if(firstCh&&chs.length) firstCh.value=chs[0];
  $q(m,"#mOK").onclick=()=>{
    const si=+$q(el,"#dreDSmp").value, jI=+$q(m,"#mI").value; const [uk,uu]=$q(m,"#mPU").value.split("|"); const f=uk==="I"?PU[uu]:1; unitMem("table",uk+"|"+uu);
    const defCal=S.lastCal||Object.assign({name:CALPRE[0].name},{C:CALPRE[0].C,A:CALPRE[0].A});
    const wrap=(o,Praw)=>{ o.raw={kind:uk==="I"?"intensity":"power",unit:uu,P:Praw}; if(uk==="I"){ o.cal={}; o.calOK=true; } else { applyCal(o,defCal,false); } return o; };
    let added=0;
    const jSy=type==="scurve"?+$q(m,"#mSy").value:-1, nRep=Math.max(1,+(type==="scurve"?$q(m,"#mRep").value:1)||1), skipI=type==="scurve"&&$q(m,"#mInt")&&$q(m,"#mInt").checked&&interpRows;
    const doSplit=type==="scurve"&&$q(m,"#mSplit").checked, doBg=type==="scurve"&&$q(m,"#mBgAuto").checked; const msgs=[];
    m.querySelectorAll("select[data-col]").forEach(sel=>{
      if(!sel.value) return; const j=+sel.dataset.col; if(j===jI) return;
      const I=[],y=[],Pr=[],sy=[]; let nI=0; T.cols[jI].forEach((p,r)=>{ const v=T.cols[j][r]; if(skipI&&interpRows[r]){ nI++; return; } if(isFinite(p)&&p>0&&!isNaN(v)){ I.push(p*f); Pr.push(p); y.push(v); if(jSy>=0) sy.push(T.cols[jSy][r]/Math.sqrt(nRep)); } });
      if(!I.length) return;
      if(type==="scurve"){ let d=wrap({type,ch:sel.value,label:T.header[j],I,y,sbg:+$q(m,"#mSbg").value||0,k:+$q(m,"#mK").value||3,relErr:+$q(m,"#mRel").value||0.05,dir:+$q(m,"#mDir").value},Pr);
        if(jSy>=0){ d.sy=sy; d.syNote=`${T.header[jSy]}/√${nRep}`; } if(nI) d.nInterp=nI;
        let parts=[d]; if(doSplit){ const r=splitUpDown(d); if(r) parts=r; }
        if(doBg){ const up=parts.find(q=>(q.dir||1)>0)||parts[0]; const e=bgEstimate(up.I,up.y); if(e){ parts.forEach(q=>subtractBg(q,e)); msgs.push(`배경 ${sciTxt(e.a,3)} + ${sciTxt(e.b,2)}·I (σ ${sciTxt(e.sd,2)})`); } else msgs.push("배경 추정 실패 — 입력값 그대로"); }
        parts.forEach(q=>S.samples[si].data.push(q)); added+=parts.length-1; }
      else S.samples[si].data.push(wrap({type,ch:sel.value,label:T.header[j],I,t:y,relErr:+$q(m,"#mRel").value||0.05,frac:type==="rise"?(+$q(m,"#mFrac").value||0.95):undefined},Pr));
      added++;
    });
    if(!added){ note("채널을 하나 이상 고르세요","warn"); return; }
    PENDING=null; S.result=null; save(); rerender(); note(`데이터 ${added}개 추가${msgs.length?" · "+[...new Set(msgs)].join(" · "):""}`);
  };
}
function synthUI(el,chs){
  const m=$q(el,"#dreMap");
  m.innerHTML=`<div class="dre-map"><div class="muted">현재 매개변수를 참값으로 삼아 잡음 있는 S-curve 를 만든다. 피팅 전략·identifiability를 참값과 비교해 시험하는 용도.</div>
   <div class="row"><label class="dre-f">Channel <select class="dre-sel" id="sCh">${chs.map(c=>`<option>${esc(c)}</option>`).join("")}</select></label>
   <label class="dre-f">Points <input type="number" id="sN" value="40" min="5" max="200"></label>
   <label class="dre-f">Max counts <input type="number" id="sMax" value="1e6" step="any"></label>
   <label class="dre-f">Background <input type="number" id="sBg" value="40" step="any"></label>
   <label class="dre-f">multiplicative noise <input type="number" id="sRel" value="0.05" step="0.01"></label></div>
   <div class="row"><button class="btn sm primary" id="sOK">Generate</button></div></div>`;
  $q(m,"#sOK").onclick=()=>{
    const ch=$q(m,"#sCh").value, k=MODEL.channels.findIndex(c=>c.id===ch), n=Math.max(5,+$q(m,"#sN").value|0);
    const f=PU[S.pu], I=[...Array(n).keys()].map(i=>S.simRange[0]*f*Math.pow(S.simRange[1]/S.simRange[0],i/(n-1)));
    const si=+$q(el,"#dreDSmp").value; const env=envOf(si);
    const o=D.observe(MODEL,env,I,1,S.nm,k,{beam:S.beam}); const mx=Math.max(...o.y); const Sc=(+$q(m,"#sMax").value||1e6)/mx, bg=+$q(m,"#sBg").value||0, rel=+$q(m,"#sRel").value||0;
    const rnd=D.mulberry(Date.now()%100000); const gauss=()=>{ let u=0,v=0; while(!u) u=rnd(); v=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); };
    const y=o.y.map(v=>{ const y0=Sc*v; return Math.max(y0+bg+gauss()*Math.sqrt(y0+bg),0)*Math.exp(rel*gauss())-bg; });
    S.samples[si].data.push({type:"scurve",ch,label:"Synthetic",I,y,sbg:Math.sqrt(bg),k:3,relErr:Math.max(rel,0.01),dir:1,truth:Object.fromEntries(Object.keys(S.params).map(n=>[n,S.params[n].v]))});
    S.result=null; save(); rerender(); note("합성 데이터를 추가했습니다 — 참값은 데이터에 함께 저장됩니다");
  };
}
/* ---------- 피팅 ---------- */
function problemOf(){
  const params={}; Object.entries(S.params).forEach(([n,p])=>{ params[n]={v:+p.v,lo:+p.lo,hi:+p.hi,fixed:!!p.fixed,per:p.per||"global"}; });
  return {spec:S.spec,nm:S.nm,beam:S.beam,params,samples:S.samples.map(s=>({id:s.id,vars:s.vars||{},data:s.data.map(d=>Object.assign({},d,{truth:undefined}))}))};
}
function panelFit(el){
  const nData=S.samples.reduce((a,s)=>a+s.data.length,0);
  const free=Object.entries(S.params).filter(([n,p])=>!p.fixed);
  const nFree=free.reduce((a,[n,p])=>a+(p.per==="sample"?S.samples.length:1),0);
  const hasKin=S.samples.some(s=>s.data.some(d=>d.type==="rise"||d.type==="decay"));
  const running=!!RUN;
  el.innerHTML=`${!nData?'<div class="note warn">데이터가 없습니다 — ‘시료·데이터’에서 넣거나 합성 데이터를 만드세요.</div>':""}
   <div class="muted">Free parameters ${nFree} · datasets ${nData} (배율은 데이터셋마다 profiled out)</div>
   ${(()=>{ const L=allIssues().filter(q=>q.lvl!=="info"); if(!L.length) return ""; const ne=L.filter(q=>q.lvl==="err").length;
      return `<div class="note warn" style="margin-top:6px"><b>Data check: ${L.length} issue${L.length>1?"s":""}${ne?` (${ne} to fix before fitting)`:""}</b> — ${L.map(q=>`${esc(q.sid)} · ${esc(q.ch)}: ${esc(q.code==="mix"?"up·down 이 한 데이터셋에 섞임":q.code==="t"?"rise/decay 시간 값 비정상":q.code==="bg"?"배경을 빼지 않은 것으로 보임":q.code)}`).join(" / ")}. ‘시료·데이터’ 목록에서 바로 고칠 수 있다.</div>`; })()}
   ${(()=>{ const u=S.samples.reduce((a,s)=>a+s.data.filter(d=>d.raw&&!d.calOK).length,0); return u?`<div class="note warn" style="margin-top:6px">파워 보정을 확인하지 않은 데이터 ${u}개 — 임시 보정(마지막 사용값 또는 NIR 950 DMSP·Oil)이 걸려 있다. σ·I_th 같은 절대값은 보정을 확인한 뒤에만 의미가 있다. ⑦의 ‘환산 무관’ 양은 영향이 없다. ‘시료·데이터’ 목록의 ‘파워 보정’에서 확인하세요.</div>`:""; })()}
   ${!hasKin&&nData?`<div class="note warn" style="margin-top:6px">시간 데이터(수명·상승 시간)가 없습니다. 검출 배율을 모르는 S-curve 만으로는 <b>모든 속도·단면적을 같은 배수로 바꾸는 방향</b>이 정해지지 않습니다 — 절대 속도 하나(예: 저농도 수명)를 고정하거나 감쇠·상승 시간을 넣으세요.</div>`:""}
   <div class="row" style="margin-top:8px"><label class="dre-f">Starts <input type="number" id="fS" value="${S.fit.starts}" min="1" max="1000"></label>
     <label class="dre-f">Max iterations <input type="number" id="fM" value="${S.fit.maxit}" min="10" max="2000"></label>
     <label class="dre-f">seed <input type="number" id="fSeed" value="${S.fit.seed}"></label>
     <label class="dre-f" title="논리 스레드 중 피팅에 쓸 비율. 피팅 중에 바꿔도 곧바로 반영된다">CPU limit <select id="fCpu">${CPU_OPTS.map(([v,l])=>`<option value="${v}" ${v===cpuPct()?"selected":""}>${l}</option>`).join("")}</select></label>
     <label class="dre-f" title="0 = CPU 사용 한도와 측정 결과로 자동 결정">Workers <input type="number" id="fW" value="${S.fit.workers||0}" min="0" max="32"> <span class="muted">(0 = auto)</span></label>
     <span class="muted" id="fCpuN" style="flex-basis:100%">${cpuNote()}</span>
     <label class="dre-tog"><input type="checkbox" id="fP" ${S.fit.profile?"checked":""}> profile likelihood</label>
     <label class="dre-f">±<input type="number" id="fSpan" value="${S.fit.span}" step="0.5" min="0.5" max="6"> decade · <input type="number" id="fN" value="${S.fit.npts}" min="5" max="41"> pts</label></div>
   <div class="row"><button class="btn primary" id="fRun" ${running||!nData||!nFree?"disabled":""}>Run fit</button><button class="btn" id="fStop" ${running?"":"disabled"}>Stop</button><button class="btn" id="fAdq" ${running||!nData?"disabled":""} title="피팅 없이 지금 매개변수 값으로 모델이 데이터와 구조적으로 어긋나는지 검사">Check adequacy at current values</button></div>
   <div class="dre-prog" id="fProg" ${running?"":"hidden"}><div class="bar"><i id="fBar"></i></div><span id="fTxt"></span></div>
   <div id="fHw" style="margin-top:6px">${hwHTML()}</div>
   <details class="dre-src"><summary>Methods</summary><ul>
    <li>log₁₀ 매개변수 공간의 Latin hypercube multi-start → Levenberg–Marquardt(box constraints, gain ratio 갱신). 기울기는 implicit function theorem로 구한 정확한 감도.</li>
    <li>residual = (ln y − ln ŷ − ln s)/σ, σ² = multiplicative noise² + (배경σ/y)². 배율 s 는 데이터셋마다 가중 평균으로 profile out. 배경 k·σ 아래 점은 censored(상한) — 모델이 상한을 넘을 때만 벌점.</li>
    <li>global minimum 도달 판정: 최적 χ² 에서 0.5 이내로 끝난 출발점의 비율(waterfall).</li>
    <li>95 % 구간: profile likelihood(나머지를 다시 최적화) Δχ² ≤ 3.84·max(1, χ²/자유도). 구간 끝이 탐색 범위에 닿으면 ‘열림’ = 그 방향으로는 데이터가 값을 정하지 못함.</li></ul></details>`;
  const upd=()=>{ S.fit={cpu:+$q(el,"#fCpu").value||75,workers:Math.max(0,+$q(el,"#fW").value|0),starts:+$q(el,"#fS").value||40,maxit:+$q(el,"#fM").value||150,seed:+$q(el,"#fSeed").value||1,profile:$q(el,"#fP").checked,span:+$q(el,"#fSpan").value||2,npts:+$q(el,"#fN").value||13}; save(); };
  el.querySelectorAll("input").forEach(i=>i.onchange=upd);
  const cpuCh=()=>{ upd(); const n=$q(el,"#fCpuN"); if(n) n.textContent=cpuNote(); if(RUN&&LIVE.pump){ LIVE.K=nWorkers(); LIVE.pump(); note(`CPU 한도 변경 → 워커 ${LIVE.K}개 (진행 중인 계산이 끝나는 대로 반영)`); } };
  $q(el,"#fCpu").onchange=cpuCh; $q(el,"#fW").onchange=cpuCh;
  $q(el,"#fRun").onclick=()=>{ upd(); startFit(); };
  $q(el,"#fStop").onclick=()=>stopFit();
  $q(el,"#fAdq").onclick=()=>adequacyNow();
  const hw=$q(el,"#fHwRe"); if(hw) hw.onclick=async()=>{ if(RUN) return; const src=workerSource(); if(!src) return; hw.disabled=true; hw.textContent="Measuring…"; try{ await calibrateHW(src); note("parallel scaling을 다시 측정했습니다"); }catch(e){ note("측정 실패: "+e.message,"warn"); } renderPanel(); };
}
function workerSource(){
  const core=document.getElementById("dre-core-src");
  const src=core?core.textContent:null;
  if(!src) return null;
  return src+`
;(function(){ const D=self.DRE_CORE; let stop=false; const CACHE={pid:null,B:null};
const getB=m=>{ if(m.pid&&CACHE.pid===m.pid&&CACHE.B) return CACHE.B; const B=D.buildProblem(m.problem); CACHE.pid=m.pid||null; CACHE.B=B; return B; };
let BENCH=null;
self.onmessage=ev=>{ const m=ev.data;
  if(m.cmd==="start1"){ try{ const B=getB(m); const r=D.lm(B,Float64Array.from(m.theta),{maxit:m.maxit}); self.postMessage({type:"start1",r:{k:m.k,chi2:r.chi2,theta:Array.from(r.theta),it:r.it,nev:r.nev,conv:r.conv,why:r.why}}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="slicerow"){ try{ const B=getB(m); const th=Float64Array.from(m.theta); const row=m.xs.map(x=>{ const t=Float64Array.from(th); t[m.j]=x; t[m.k]=m.y; try{ return D.chi2(D.residuals(B,t,false).r); }catch(e){ return NaN; } }); self.postMessage({type:"slicerow",b:m.b,row}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="bench"){ /* 실제 DRE 연산(정상상태 스윕 + 감도)을 정해진 시간 동안 반복 — throughput 측정 */
    try{ if(!BENCH){ const T=D.templateLee(); const M=D.compile(T.spec); const env=Object.assign({},...Object.entries(T.params).map(([k,p])=>({[k]:p.v})),{c:8}); BENCH={M,env,I:D.logGrid(1e3,2e5,12)}; }
      const t0=performance.now(); let n=0; while(performance.now()-t0<m.ms){ D.observe(BENCH.M,BENCH.env,BENCH.I,1,1064,0,{pnames:["sE","k31"]}); n++; }
      self.postMessage({type:"bench",n,dt:performance.now()-t0}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="slice"){ try{ const B=D.buildProblem(m.problem); const sl=D.slice2D(B,m.theta,m.j,m.k,{n:17,spanJ:m.span,spanK:m.span,progress:(d,t)=>{ if(d%17===0) self.postMessage({type:"progress",done:d,total:t}); }}); self.postMessage({type:"slice",slice:sl}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="starts"){ try{ const B=D.buildProblem(m.problem); const out=[]; m.list.forEach(q=>{ const r=D.lm(B,Float64Array.from(q.theta),{maxit:m.maxit}); out.push({k:q.k,chi2:r.chi2,theta:Array.from(r.theta),it:r.it,nev:r.nev,conv:r.conv,why:r.why}); self.postMessage({type:"tick",chi2:r.chi2}); }); self.postMessage({type:"starts",out}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="profile"){ try{ const B=D.buildProblem(m.problem); const pr=D.profile(B,m.best,m.j,m.opt); self.postMessage({type:"profile",j:m.j,pr}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd==="final"){ try{ self.postMessage({type:"result",...finalize(D.buildProblem(m.problem),m.ms,m.prof,m.opts,m.t0)}); }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); } return; }
  if(m.cmd!=="fit") return;
  try{
    const B=D.buildProblem(m.problem); const o=m.opts; const t0=Date.now();
    const ms=D.multistart(B,{starts:o.starts,seed:o.seed,maxit:o.maxit},p=>self.postMessage({type:"progress",phase:"multi-start",done:p.done,total:p.total,best:p.best}));
    const best=ms[0]; const info=D.dofInfo(B,best); const thr=3.84*info.s2;
    const prof=[];
    if(o.profile){ for(let j=0;j<B.free.length;j++){ self.postMessage({type:"progress",phase:"Profile "+B.free[j].key,done:j,total:B.free.length,best:best.chi2});
      prof.push(D.profile(B,best,j,{span:o.span,npts:o.npts,thr,maxit:Math.min(o.maxit,120)})); } }
    self.postMessage({type:"result",...finalize(B,ms,prof,o,t0)});
  }catch(e){ self.postMessage({type:"error",msg:String(e&&e.message||e)}); }
};
function finalize(B,ms,prof,o,t0){
    const best=ms[0]; const info=D.dofInfo(B,best); const thr=3.84*info.s2;
    const ev2=D.fim(B,best.theta);
    const R=D.residuals(B,Float64Array.from(best.theta),false);
    const scales={}; B.problem.samples.forEach((s,si)=>s.data.forEach((d,di)=>{ if(d._fit&&d._fit.logS!=null) scales[si+"|"+di]=Math.exp(d._fit.logS); }));
    const bestBySample=B.problem.samples.map((s,si)=>{ const e=D.envFor(B,Float64Array.from(best.theta),si); const o2={}; B.free.forEach(f=>{ if(f.si===-1||f.si===si) o2[f.nm]=e[f.nm]; }); return o2; });
    const nTot=R.r.length, k=B.free.length+info.nScale;
    self.postMessage({type:"progress",phase:"Adequacy · detailed analysis",done:0,total:1,best:best.chi2});
    let adequacy=null, insight=null; const fI=v=>(v/o.puF).toPrecision(3)+" "+o.pu;
    try{ adequacy=D.adequacy(B,best.theta,{fmtI:fI}); }catch(e){ adequacy={items:[{key:"A0",title:"Adequacy check",level:"warn",summary:"계산 실패: "+(e&&e.message||e)}],worst:"warn",verdict:"적합성 검사를 끝내지 못했습니다."}; }
    try{ insight=D.insight(B,best,ms,{thr}); }catch(e){ insight=null; }
    let derived=null; try{ self.postMessage({type:"progress",phase:"Predictions",done:0,total:1,best:best.chi2}); derived=D.derivedIntervals(B,best,prof,thr,{}); }catch(e){ derived=null; }
    return {adequacy,insight,derived,free:B.free,ms:ms.map(r=>({chi2:r.chi2,theta:r.theta,conv:r.conv})),best,info,thr,prof,fim:ev2,scales,bestBySample,aicc:D.aicc(best.chi2,nTot,k),n:nTot,k,time:(Date.now()-t0)/1000,workers:o.workers||1};
}
})();`;
}
let POOL=[];
const HWKEY="innoanalytics_dre_hw_v1";
function hwInfo(){ const nav=typeof navigator!=="undefined"?navigator:{}; return {hc:nav.hardwareConcurrency||2, mem:nav.deviceMemory||null, ua:(nav.userAgentData&&nav.userAgentData.platform)||nav.platform||"", cores:null}; }
function hwCal(){ try{ const c=JSON.parse(localStorage.getItem(HWKEY)||"null"); const h=hwInfo(); if(c&&c.hc===h.hc&&c.ua===h.ua) return c; }catch(e){} return null; }
/* CPU 사용 한도: 논리 스레드의 cpu % 까지만 워커를 쓴다. 브라우저는 스레드 우선순위를 낮출 수 없으므로
   '남겨 두는 스레드 수'로 다른 작업의 응답성을 지킨다. 한도 안에서는 측정된 최적 워커 수를 넘지 않는다. */
const CPU_OPTS=[[100,"Max (100 %)"],[90,"90 %"],[75,"Balanced (75 %)"],[50,"Relaxed (50 %)"],[25,"Background (25 %)"]];
function cpuPct(){ const v=+(S.fit&&S.fit.cpu); return v>0&&v<=100?v:75; }
function capK(pct){ const hc=hwInfo().hc; return Math.max(1,Math.min(32,Math.floor(hc*pct/100+1e-9))); }
function nWorkers(pct){ const h=hwInfo(); const want=+S.fit.workers||0; if(want>0) return Math.min(32,want); const cap=capK(pct==null?cpuPct():pct); const c=hwCal(); if(c&&c.best) return Math.min(cap,c.best); return Math.min(cap,Math.max(1,h.hc-1)); }
/* 측정 곡선에서 K 워커의 throughput / 최고 throughput (선형 보간) */
function speedAt(K){ const c=hwCal(); if(!c||!c.rows||!c.rows.length) return null; const r=c.rows, top=Math.max(...r.map(q=>q.thr)); if(K<=r[0].K) return r[0].thr*K/r[0].K/top;
  for(let i=1;i<r.length;i++) if(K<=r[i].K){ const a=r[i-1],b=r[i],f=(K-a.K)/(b.K-a.K); return (a.thr+f*(b.thr-a.thr))/top; } return r[r.length-1].thr/top; }
function cpuNote(){ const hc=hwInfo().hc, K=nWorkers(), sp=speedAt(K); return `${K} worker${K>1?"s":""} / ${hc} logical threads → CPU ≈ ${Math.round(100*Math.min(1,K/hc))} %${sp!=null?` · speed ${Math.round(100*Math.min(1,sp))} % of max`:""}${+S.fit.workers>0?" (workers set manually)":""}`; }
const LIVE={K:1,pump:null};
/* scaling 측정: 워커 수 K 마다 같은 시간 동안 실제 DRE 연산을 돌려 throughput을 잰다.
   논리 스레드 수만으로는 SMT·성능/효율 코어·전원 모드·다른 부하를 알 수 없어서 실측으로 고른다. */
async function calibrateHW(src,onStep){
  const h=hwInfo(); const hc=h.hc; const Ks=[...new Set([1,2,3,4,6,8,12,16,24,32].filter(k=>k<=hc).concat([Math.max(1,hc-1),hc]))].sort((a,b)=>a-b);
  const url=URL.createObjectURL(new Blob([src],{type:"text/javascript"})); const ws=[...Array(hc)].map(()=>new Worker(url)); URL.revokeObjectURL(url);
  const bench=(w,ms)=>new Promise((res,rej)=>{ w.onmessage=ev=>{ const m=ev.data; if(m.type==="bench") res(m); else if(m.type==="error") rej(new Error(m.msg)); }; w.onerror=e=>rej(new Error(e.message||"")); w.postMessage({cmd:"bench",ms}); });
  try{ await Promise.all(ws.map(w=>bench(w,60)));      /* warm-up: 컴파일·JIT */
    const rows=[]; for(const K of Ks){ const r=await Promise.all(ws.slice(0,K).map(w=>bench(w,350))); const thr=r.reduce((a,q)=>a+q.n/q.dt*1000,0); rows.push({K,thr}); onStep&&onStep(rows.length,Ks.length); }
    const t1=rows[0].thr; rows.forEach(q=>{ q.speedup=q.thr/t1; q.eff=q.speedup/q.K; });
    const top=Math.max(...rows.map(q=>q.thr)); const pick=rows.find(q=>q.thr>=0.95*top);       /* 최고의 95 % 안에서 가장 적은 워커 */
    /* 물리 코어 추정: 효율(속도비/K)이 0.8 아래로 처음 떨어지기 직전의 K */
    let phys=null; for(const q of rows){ if(q.eff>=0.8) phys=q.K; else break; }
    const smt=phys&&hc>phys?(rows.find(q=>q.K===hc)||rows[rows.length-1]).thr/(rows.find(q=>q.K===phys)||rows[0]).thr:null;
    const cal=Object.assign(h,{rows,best:Math.max(1,Math.min(pick.K,hc>1?hc-(pick.K===hc&&hc>2?1:0):1)),top,phys,smt,date:new Date().toISOString()});
    try{ localStorage.setItem(HWKEY,JSON.stringify(cal)); }catch(e){}
    return cal;
  } finally { ws.forEach(w=>w.terminate()); }
}
function hwHTML(){
  const h=hwInfo(), c=hwCal();
  const base=`Logical threads ${h.hc}${h.mem?` · memory ≥${h.mem} GB`:""}`;
  if(!c) return `<div class="muted">${base} · scaling not measured (첫 피팅 때 약 3–4초 측정)</div>`;
  const P=palette(); const W=300,H=110, mx=Math.max(...c.rows.map(q=>q.speedup),1), kx=Math.max(...c.rows.map(q=>q.K));
  const X=k=>36+(k-1)/(Math.max(1,kx-1))*(W-46), Y=v=>H-18-(v/mx)*(H-28);
  const svg=`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="${W}" role="img" aria-label="parallel scaling"><rect width="${W}" height="${H}" fill="${P.bg}"/>
    <line x1="${X(1)}" y1="${Y(1)}" x2="${X(Math.min(kx,mx))}" y2="${Y(Math.min(kx,mx))}" stroke="${P.muted}" stroke-dasharray="3 3"/>
    <path d="${c.rows.map((q,i)=>(i?"L":"M")+X(q.K).toFixed(1)+" "+Y(q.speedup).toFixed(1)).join("")}" fill="none" stroke="${P.accent}" stroke-width="2"/>
    ${c.rows.map(q=>`<circle cx="${X(q.K)}" cy="${Y(q.speedup)}" r="3" fill="${q.K===c.best?P.accent:P.bg}" stroke="${P.accent}"><title>Workers ${q.K}: ×${q.speedup.toFixed(2)} (efficiency ${(100*q.eff).toFixed(0)} %)</title></circle>`).join("")}
    <text x="${W/2}" y="${H-3}" font-size="9.5" text-anchor="middle" fill="${P.muted}">Workers (dashed = ideal linear)</text><text x="4" y="12" font-size="9.5" fill="${P.muted}">Speed ×${mx.toFixed(1)}</text></svg>`;
  return `<div class="dre-hw"><div>${base} · measured optimal workers <b>${c.best}</b> (×${(c.rows.find(q=>q.K===c.best)||{speedup:1}).speedup.toFixed(1)})${c.phys?` · est. physical cores ≈ ${c.phys}`:""}${c.smt?` · extra-thread (SMT / E-core) gain ×${c.smt.toFixed(2)}`:""} <button class="btn xs" id="fHwRe">Re-measure</button></div>${svg}
    <div class="muted">측정 ${esc(c.date?new Date(c.date).toLocaleString():"")}. 브라우저는 논리 스레드 수만 알려 주므로, 물리 코어·SMT·성능/효율 코어·전원 모드는 throughput 곡선으로 추정한다(효율 ≥ 80 % 까지를 물리 코어로 간주). 다른 프로그램이 무거울 때 측정했다면 다시 측정.</div></div>`;
}
function startFit(){
  if(RUN) return;
  const src=workerSource(); if(!src){ note("Worker 원본을 찾지 못했습니다","warn"); return; }
  { const E=allIssues().filter(q=>q.lvl==="err"); if(E.length){ note(`피팅 전에 고칠 데이터 ${E.length}건 — ‘시료·데이터’ 목록의 표시를 확인하세요: ${E.map(q=>`${q.sid} · ${q.ch}: ${q.code==="mix"?"up·down 섞임":q.code==="t"?"시간 값 비정상":q.code}`).join(" / ")}`,"warn"); return; } }
  let prob, B; try{ prob=problemOf(); B=D.buildProblem(JSON.parse(JSON.stringify(prob))); }catch(e){ note("피팅 준비 실패: "+e.message,"warn"); return; }
  RUN={t0:Date.now()}; renderPanel();
  (async()=>{
    if(!(+S.fit.workers>0)&&!hwCal()){ setProg("Measuring parallel scaling",0,1,null); try{ await calibrateHW(src,(d,t)=>RUN&&setProg("Measuring parallel scaling",d,t,null)); }catch(e){} if(!RUN) return; renderPanel(); }
    runFit(src,prob,B);
  })();
}
function runFit(src,prob,B){
  const o=Object.assign({},S.fit,{pu:S.pu,puF:PU[S.pu]}); const t0=Date.now(); RUN.t0=t0;
  const n=Math.max(1,o.starts|0), rnd=D.mulberry(o.seed||12345);
  const starts=[Array.from(D.thetaFromParams(B))].concat(D.lhs(n-1,B.free.length,rnd).map(u=>u.map((v,j)=>B.free[j].lo+v*(B.free[j].hi-B.free[j].lo))));
  LIVE.K=Math.min(nWorkers(),starts.length); let K=LIVE.K; let Kmax=K; const pid="p"+t0;
  const url=URL.createObjectURL(new Blob([src],{type:"text/javascript"})); RUN.url=url;
  POOL=[...Array(K)].map(()=>new Worker(url)); WORKER=POOL[0];
  const fail=msg=>{ stopFit(true); note("피팅 오류: "+msg,"warn"); };
  const call=(w,msg,kind)=>new Promise((res,rej)=>{ w.onmessage=ev=>{ const m=ev.data; if(m.type===kind) res(m); else if(m.type==="error") rej(new Error(m.msg)); }; w.onerror=e=>rej(new Error(e.message||"Worker 오류")); w.postMessage(msg); });
  /* 동적 분배: 먼저 끝난 워커가 다음 일을 가져간다 (성능·효율 코어 혼재, 출발점마다 다른 계산량에 강함) */
  /* LIVE.K(현재 CPU 한도)만큼만 동시에 돌린다. 한도를 낮추면 진행 중인 계산이 끝난 워커는 쉬고, 높이면 워커를 더 만든다.
     결과는 출발점 번호·매개변수 번호로 정렬되므로 워커 수와 무관하게 같다. */
  const queue=(items,mk,kind,onRes)=>new Promise((res,rej)=>{ let next=0,busy=0,done=0; const idle=POOL.slice(); let dead=false;
    const pump=()=>{ if(dead) return; if(!RUN){ dead=true; return res(); }
      while(busy<LIVE.K&&next<items.length){ let w=idle.pop(); if(!w){ if(POOL.length>=32) break; w=new Worker(url); POOL.push(w); }
        const it=items[next++]; busy++; K=Math.max(busy,K); Kmax=Math.max(Kmax,busy);
        call(w,mk(it),kind).then(m=>{ busy--; idle.push(w); if(dead) return; onRes(m,it); done++; if(done===items.length){ dead=true; res(); } else pump(); }).catch(e=>{ dead=true; rej(e); }); }
      if(next>=items.length&&busy===0&&!dead){ dead=true; res(); } };
    LIVE.pump=pump; pump(); });
  const ms=[]; let best=Infinity; setProg(`multi-start · ${LIVE.K} worker${LIVE.K>1?"s":""}`,0,starts.length,null);
  queue(starts.map((t,k)=>({k,t})),it=>({cmd:"start1",pid,problem:prob,k:it.k,theta:it.t,maxit:o.maxit}),"start1",m=>{ ms.push(m.r); best=Math.min(best,m.r.chi2); if(RUN) setProg(`multi-start · ${LIVE.K} worker${LIVE.K>1?"s":""}`,ms.length,starts.length,best); })
  .then(async()=>{
    if(!RUN) return;
    ms.sort((a,b)=>a.chi2-b.chi2||a.k-b.k);
    const bst=ms[0]; const info=D.dofInfo(B,bst); const thr=3.84*info.s2; const prof=[];
    if(o.profile){ const P=B.free.length; let fin=0; setProg(`Profile · ${LIVE.K} worker${LIVE.K>1?"s":""}`,0,P,bst.chi2);
      await queue([...Array(P).keys()],j=>({cmd:"profile",problem:prob,best:bst,j,opt:{span:o.span,npts:o.npts,thr,maxit:Math.min(o.maxit,120)}}),"profile",m=>{ prof[m.j]=m.pr; fin++; if(RUN) setProg(`Profile · ${LIVE.K} worker${LIVE.K>1?"s":""}`,fin,P,bst.chi2); }); }
    if(!RUN) return;
    setProg("Adequacy · detailed analysis · predictions",0,1,bst.chi2);
    o.workers=Kmax; const m=await call(POOL[0],{cmd:"final",problem:prob,ms,prof,opts:o,t0},"result");
    m.at=Date.now(); S.result=m; S.useFit=true; URL.revokeObjectURL(url); RUN=null; LIVE.pump=null; POOL.forEach(w=>w.terminate()); POOL=[]; WORKER=null; save(); rerender(); note(`피팅 완료 (${m.time.toFixed(1)} s · 최대 워커 ${Kmax}개)`);
  }).catch(e=>{ if(RUN) fail(e.message||String(e)); });
}
function stopFit(silent){ LIVE.pump=null; if(RUN&&RUN.url) try{ URL.revokeObjectURL(RUN.url); }catch(e){} POOL.forEach(w=>{ try{ w.terminate(); }catch(e){} }); POOL=[]; if(WORKER){ try{ WORKER.terminate(); }catch(e){} WORKER=null; } RUN=null; renderPanel(); if(!silent) note("피팅을 중단했습니다"); }
function setProg(phase,done,total,best){
  const p=$q(ROOT,"#fProg"); if(!p) return; p.hidden=false;
  $q(ROOT,"#fBar").style.width=Math.round(100*done/Math.max(1,total))+"%";
  $q(ROOT,"#fTxt").textContent=`${phase} ${done}/${total}${best!=null&&isFinite(best)?` · min χ² ${best.toFixed(2)}`:""} · ${((Date.now()-RUN.t0)/1000).toFixed(0)} s`;
}

/* ---------------- 진단 ---------------- */
const SER10=["#7a3f6e","#1f6fa3","#b4580f","#2a7141","#6a4bb0","#b23a3a","#4a7c8c","#a0772b","#5c5c9e","#8a4f3d"], SER10D=["#e39ad0","#7fc6e8","#f0b26a","#8fd19e","#c9a8f0","#f08c86","#8cc7d4","#e0c27a","#aeb0ea","#d9a08c"];
const S10=P=>P.dark?SER10D:SER10;
const legendHTML=items=>`<div class="dre-legend">${items.map(([n,c])=>`<span><i style="background:${c}"></i>${esc(n)}</span>`).join("")}</div>`;
const LV={ok:["g0","Pass"],warn:["g3","Warn"],fail:["g4","Fail"],na:["g2","N/A"]};
const lvBadge=l=>`<span class="gr ${LV[l][0]}">${LV[l][1]}</span>`;
function normInv(p){ /* Acklam */ const a=[-39.69683028665376,220.9460984245205,-275.9285104469687,138.357751867269,-30.66479806614716,2.506628277459239],b=[-54.47609879822406,161.5858368580409,-155.6989798598866,66.80131188771972,-13.28068155288572],c=[-0.007784894002430293,-0.3223964580411365,-2.400758277161838,-2.549732539343734,4.374664141464968,2.938163982698783],d=[0.007784695709041462,0.3224671290700398,2.445134137142996,3.754408661907416];
  const pl=0.02425; if(p<pl){ const q=Math.sqrt(-2*Math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  if(p>1-pl){ const q=Math.sqrt(-2*Math.log(1-p)); return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  const q=p-0.5, r=q*q; return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1); }
/* 행렬 히트맵 (상관: 발산, fraction: 순차) */
function heat(o){
  const P=palette(), R=o.rows.length, C=o.cols.length, cw=o.cw||34, ch=o.ch||24, lw=o.lw||92, th=o.th||70;
  const W=Math.max(lw+C*cw+8,o.minw||0), H=th+R*ch+6+(o.groups&&o.groups.length?20:0); const g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="${Math.min(W,o.maxw||W)}" role="img" aria-label="${esc(o.label||"")}"><rect width="${W}" height="${H}" fill="${P.bg}"/>`];
  o.cols.forEach((c,j)=>g.push(`<text transform="translate(${lw+j*cw+cw/2+4},${th-6}) rotate(-50)" font-size="10" fill="${P.ink2}">${esc(c)}</text>`));
  (o.groups||[]).forEach(gp=>{ const x0=lw+gp.a*cw+2, x1=lw+(gp.b+1)*cw-2; g.push(`<line x1="${x0}" x2="${x1}" y1="${th+R*ch+6}" y2="${th+R*ch+6}" stroke="${P.ink2}"/><text x="${(x0+x1)/2}" y="${th+R*ch+19}" font-size="10" text-anchor="middle" fill="${P.ink2}">${esc(gp.label)}</text>`); });
  o.rows.forEach((r,i)=>{ g.push(`<text x="${lw-6}" y="${th+i*ch+ch/2+4}" font-size="10.5" text-anchor="end" fill="${P.ink2}">${esc(r)}</text>`);
    o.cols.forEach((c,j)=>{ const v=o.M[i][j]; const f=o.color(v,P); g.push(`<rect x="${lw+j*cw+1}" y="${th+i*ch+1}" width="${cw-2}" height="${ch-2}" rx="3" fill="${f.c}" fill-opacity="${f.a}" stroke="${P.line}" stroke-width=".5"><title>${esc(r)} × ${esc(c)}: ${esc(o.fmt(v))}</title></rect>`);
      if(o.text) g.push(`<text x="${lw+j*cw+cw/2}" y="${th+i*ch+ch/2+3.5}" font-size="9.5" text-anchor="middle" fill="${f.a>0.55?P.bg:P.ink2}">${esc(o.text(v))}</text>`); }); });
  g.push("</svg>"); return g.join(""); }
/* 누적 막대 (과정별 flux fraction) */
function stackBars(o){
  const P=palette(), n=o.rows.length, W=o.w||620, lw=150, bh=18, gap=7, H=n*(bh+gap)+34; const pw=W-lw-10;
  const g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="${esc(o.label||"")}"><rect width="${W}" height="${H}" fill="${P.bg}"/>`];
  o.rows.forEach((r,i)=>{ const y=6+i*(bh+gap); g.push(`<text x="${lw-6}" y="${y+bh/2+4}" font-size="10.5" text-anchor="end" fill="${P.ink2}">${esc(r.label)}</text>`); let x=lw;
    r.v.forEach((v,k)=>{ const w=Math.max(0,v*pw); if(w<0.3) return; g.push(`<rect x="${x.toFixed(1)}" y="${y}" width="${w.toFixed(1)}" height="${bh}" fill="${o.colors[k]}"><title>${esc(o.names[k])}: ${(100*v).toFixed(1)} %</title></rect>`); if(w>26) g.push(`<text x="${(x+w/2).toFixed(1)}" y="${y+bh/2+3.5}" font-size="9" text-anchor="middle" fill="${P.dark?"#1b2422":"#fff"}">${Math.round(100*v)}</text>`); x+=w; }); });
  g.push(`<text x="${lw}" y="${H-8}" font-size="10" fill="${P.muted}">0 %</text><text x="${lw+pw}" y="${H-8}" font-size="10" text-anchor="end" fill="${P.muted}">100 % (all transitions at that state)</text></svg>`); return g.join(""); }
/* 신뢰구간 비교 (forest): 프로파일 대 곡률(Wald), log10(값/최적) */
function forest(R,X){
  const P=palette(), n=R.free.length, W=620, lw=110, rh=26, H=n*rh+44, span=Math.max(1,...R.prof.map(p=>p?Math.max(Math.abs(p.lo-R.best.theta[R.prof.indexOf(p)]),Math.abs(p.hi-R.best.theta[R.prof.indexOf(p)])):1),...(X?X.se.map(s=>isFinite(s)?Math.min(3,1.96*s):1):[1]));
  const sp=Math.min(4,Math.ceil(span*2)/2); const x0=lw, pw=W-lw-20; const Xp=v=>x0+(v+sp)/(2*sp)*pw;
  const g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Confidence interval comparison"><rect width="${W}" height="${H}" fill="${P.bg}"/>`];
  for(let t=-sp;t<=sp+1e-9;t+=sp>2?1:0.5){ g.push(`<line x1="${Xp(t)}" x2="${Xp(t)}" y1="8" y2="${H-30}" stroke="${t===0?P.muted:P.grid}"/><text x="${Xp(t)}" y="${H-16}" font-size="10" text-anchor="middle" fill="${P.muted}">${t>0?"+":""}${t}</text>`); }
  g.push(`<text x="${x0+pw/2}" y="${H-3}" font-size="10.5" text-anchor="middle" fill="${P.ink2}">log₁₀(value / best) — thick bar: profile 95 %, thin line: curvature (Wald) 95 %</text>`);
  R.free.forEach((f,j)=>{ const y=12+j*rh, b=R.best.theta[j]; const pr=R.prof[j];
    g.push(`<text x="${lw-8}" y="${y+5}" font-size="11" text-anchor="end" fill="${P.ink}">${esc(f.key)}</text>`);
    if(X){ const w=isFinite(X.se[j])?1.96*X.se[j]:sp*1.01; const a=Math.max(-sp,-w), c=Math.min(sp,w); g.push(`<line x1="${Xp(a)}" x2="${Xp(c)}" y1="${y+7}" y2="${y+7}" stroke="${P.ink2}" stroke-width="1.4"/>${w>sp?`<text x="${Xp(sp)-2}" y="${y+11}" font-size="9" text-anchor="end" fill="${P.muted}">▸</text>`:""}`); }
    if(pr){ const a=Math.max(-sp,pr.lo-b), c=Math.min(sp,pr.hi-b); const col=pr.openLo&&pr.openHi?"#b23a3a":(pr.openLo||pr.openHi)?"#b4580f":P.accent;
      g.push(`<rect x="${Xp(a)}" y="${y-2}" width="${Math.max(2,Xp(c)-Xp(a))}" height="7" rx="3" fill="${col}" fill-opacity=".75"><title>${esc(f.key)} profile 95 %: ${sciTxt(Math.pow(10,pr.lo),3)} – ${sciTxt(Math.pow(10,pr.hi),3)}${pr.openLo?" (open below)":""}${pr.openHi?" (open above)":""}</title></rect>`);
      if(pr.openLo) g.push(`<text x="${Xp(a)-3}" y="${y+5}" font-size="10" text-anchor="end" fill="${col}">◂</text>`); if(pr.openHi) g.push(`<text x="${Xp(c)+3}" y="${y+5}" font-size="10" fill="${col}">▸</text>`); }
    g.push(`<circle cx="${Xp(0)}" cy="${y+1.5}" r="3" fill="${P.ink}"/>`); });
  g.push("</svg>"); return g.join(""); }
function insOpen(k,def){ S.insOpen=S.insOpen||{}; return (k in S.insOpen?S.insOpen[k]:def)?"open":""; }
function adequacyHTML(A,src){
  if(!A) return "";
  return `<div class="dre-adq"><div class="note ${A.worst==="fail"?"danger":A.worst==="warn"?"warn":""}"><b>Model adequacy:</b> ${esc(A.verdict)} <span class="muted">(${esc(src)})</span></div>
   <table class="dre-tbl"><tbody>${A.items.map(q=>`<tr><td style="width:34px" class="mono">${q.key}</td><td style="width:74px">${lvBadge(q.level)}</td><td><b>${esc(adqTitle(q))}</b> — ${esc(q.summary)}${q.detail?`<div class="muted" style="margin-top:2px">${esc(q.detail)}</div>`:""}</td></tr>`).join("")}</tbody></table>
   <div class="muted" style="margin-top:4px">모든 항목은 같은 파워·시간에서 실측과 DRE 모델을 비교한 residual 진단이다 (INNO Prism SegPA 의 s·I<sub>th</sub> 추정과 다른 것). ‘통과’는 모델이 옳다는 증명이 아니라, 찾을 수 있는 구조적 모순이 없다는 뜻이다.</div></div>`;
}
function renderDiag(){
  const el=$q(ROOT,"#dreDiag"); if(!el) return;
  const R=S.result;
  if(!R){ if(S.adeqNow){ el.hidden=false; el.innerHTML=`<h3>Model adequacy <span class="hint">at current parameter values · before fit</span></h3>${adequacyHTML(S.adeqNow,"current parameters")}`; } else el.hidden=true; return; }
  el.hidden=false;
  const P=palette(), X=R.insight||null, A=R.adequacy||null;
  const pre=S.adeqNow&&S.adeqNow.at>(R.at||0)?`<div style="margin-bottom:10px"><h4 class="dre-h4">Model adequacy — at current parameters (values changed after fit)</h4>${adequacyHTML(S.adeqNow,"current parameters")}</div>`:"";
  const succ=R.ms.filter(r=>r.chi2-R.best.chi2<=0.5).length;
  const truth=(()=>{ for(const s of S.samples) for(const d of s.data) if(d.truth) return d.truth; return null; })();
  const idf=R.free.map((f,j)=>{ const pr=R.prof[j]; if(!pr) return "na"; if(pr.openLo&&pr.openHi) return "g4"; if(pr.openLo||pr.openHi) return "g3"; const st=pr.pts.length>1?Math.abs(pr.pts[1].v-pr.pts[0].v):0; if(pr.hi-pr.lo<1e-9) return st<=0.5?"g0":"g2"; return (pr.hi-pr.lo)<=0.5?"g0":"g2"; });
  const cnt=k=>idf.filter(v=>v===k).length;
  const rows=R.free.map((f,j)=>{
    const v=Math.pow(10,R.best.theta[j]); const pr=R.prof[j];
    const stp=pr&&pr.pts.length>1?Math.abs(pr.pts[1].v-pr.pts[0].v):0;
    const ci=!pr?"—":(pr.hi-pr.lo<1e-9&&!pr.openLo&&!pr.openHi)?`<span title="프로파일 격자 간격보다 좁다 — 점 수를 늘리거나 범위를 줄이면 실제 폭이 보인다">narrower than grid (±${stp.toFixed(2)} decade)</span>`:`${pr.openLo?"≤ ":""}${sci(Math.pow(10,pr.lo),3)} – ${sci(Math.pow(10,pr.hi),3)}${pr.openHi?" ≥":""}`;
    const wald=!X?"—":(isFinite(X.se[j])&&1.96*X.se[j]>6)?'<span class="muted" title="곡률이 거의 0 — 선형 근사로는 6 decade 이상 퍼짐">≫ 10<sup>6</sup></span>':!isFinite(X.se[j])?(X.atBound&&X.atBound[j]?'<span class="muted" title="탐색 경계에 붙어 곡률 구간이 무의미">at bound</span>':'<span class="muted" title="곡률이 0 인 방향 — 선형 근사로는 정해지지 않음">∞</span>'):`×/÷ ${Math.pow(10,1.96*X.se[j])<1e3?Math.pow(10,1.96*X.se[j]).toFixed(2):"10<sup>"+(1.96*X.se[j]).toFixed(1)+"</sup>"}`;
    const flag={g4:'<span class="gr g4">non-identifiable</span>',g3:'<span class="gr g3">One-sided</span>',g0:'<span class="gr g0">Determined</span>',g2:'<span class="gr g2">Broad</span>',na:""}[idf[j]];
    const tv=truth&&truth[f.nm]!=null?Math.log10(v/truth[f.nm]):null;
    return `<tr><td><b>${esc(f.key)}</b></td><td class="num">${sci(v,4)}</td><td class="num">${ci}</td><td class="num">${wald}</td><td>${flag}</td>${truth?`<td class="num">${tv==null?"":(tv>0?"+":"")+tv.toFixed(2)}</td>`:""}</tr>`; }).join("");
  const glob=succ/R.ms.length>=0.3?"ok":succ>=2?"warn":"fail";
  const ncl=X?X.clusters.filter(c=>c.chi2-R.best.chi2>R.thr).length:null;
  const strip=`<div class="dre-strip">
    <div><span class="k">Global minimum</span>${lvBadge(glob)}<span>${succ}/${R.ms.length} starts reached the lowest χ²${ncl!=null?` · other local minima: ${ncl}`:""}</span></div>
    <div><span class="k">identifiability</span>${lvBadge(cnt("g4")?"fail":cnt("g3")?"warn":"ok")}<span>Determined ${cnt("g0")} · Broad ${cnt("g2")} · One-sided ${cnt("g3")} · Non-identifiable ${cnt("g4")}</span></div>
    <div><span class="k">Model adequacy</span>${A?lvBadge(A.worst):lvBadge("na")}<span>${A?esc(A.items.filter(q=>q.level==="fail").map(q=>q.key).join(", ")||"No structural contradiction"):"—"}</span></div>
    <div><span class="k">Residual size</span>${lvBadge(R.info.s2>5?"fail":R.info.s2>2?"warn":"ok")}<span>χ²/dof ${(R.best.chi2/R.info.dof).toFixed(2)}</span></div></div>`;
  let html=pre+`<h3>Fit results <span class="hint">χ² ${R.best.chi2.toFixed(2)} · dof ${R.info.dof} · AICc ${R.aicc.toFixed(1)} · ${R.time.toFixed(1)} s</span></h3>${strip}
   ${A&&A.worst==="fail"?`<div class="note danger" style="margin-top:6px"><b>${esc(A.verdict)}</b> 아래 ‘Model adequacy’에서 이유를 확인하세요.</div>`:""}
   <table class="dre-tbl" style="margin-top:8px"><thead><tr><th>Parameter</th><th class="num">Best</th><th class="num">Profile 95 %</th><th class="num" title="최적점 곡률(FIM⁻¹·s²)로 본 95 % 배수 구간 — 선형 근사">Curvature 95 %</th><th>Status</th>${truth?'<th class="num">log₁₀(best/true)</th>':""}</tr></thead><tbody>${rows}</tbody></table>
   <h4 class="dre-h4" style="margin-top:10px">Detailed analysis — evidence</h4>`;
  const sec=(k,title,sub,body,def)=>`<details class="dre-ins" data-k="${k}" ${insOpen(k,def)}><summary>${title} <span class="muted">${sub}</span></summary><div class="dre-insb">${body}</div></details>`;
  /* 1 적합성 */
  html+=sec("adq","① Model adequacy",A?`${LV[A.worst][1]} · ${A.items.filter(q=>q.level==="fail"||q.level==="warn").length} item(s) flagged`:"",A?adequacyHTML(A,"fit optimum"):'<div class="muted">이 결과에는 적합성 검사가 없습니다 — 다시 피팅하세요.</div>',A&&A.worst!=="ok");
  /* 2 매개변수별 근거 */
  { let b=`<div class="dre-h4">Confidence intervals: profile likelihood vs. curvature (Wald)</div>${forest(R,X)}
      <div class="muted">두 구간이 비슷하면 최적점 근처가 2차 곡면(선형 근사 유효)이다. 프로파일이 더 넓거나 한쪽으로 치우치면 비선형 — 곡률 구간만 보고하면 과신이 된다. 막대 끝 ◂▸ 는 탐색 범위까지 열림.</div>`;
    if(X){ /* 정보 기여 히트맵: 매개변수 × (데이터셋·파워 구간) */
      const cols=[], key=[], groups=[]; X.sets.forEach(st=>{ const pts=X.pts.filter(p=>p.si===st.si&&p.di===st.di); const a0=cols.length; if(st.type==="scurve"){ const Is=pts.map(p=>p.I).filter(v=>v>0); if(!Is.length) return; const lo=Math.log(Math.min(...Is)), hi=Math.log(Math.max(...Is))+1e-9; const nb=Math.min(6,Math.max(1,Math.round(pts.length/6)));
          for(let q=0;q<nb;q++){ const a=Math.exp(lo+(hi-lo)*q/nb), c=Math.exp(lo+(hi-lo)*(q+1)/nb); cols.push(`${sciTxt(a/PU[S.pu],2)}–${sciTxt(c/PU[S.pu],2)}`); key.push(p=>p.si===st.si&&p.di===st.di&&p.I>=a&&p.I<c); } }
        else { cols.push(st.type==="rise"?"Rise":"Decay"); key.push(p=>p.si===st.si&&p.di===st.di); }
        groups.push({a:a0,b:cols.length-1,label:enL(st.label)}); });
      const M=R.free.map((f,j)=>{ const tot=X.pts.reduce((a,p)=>a+p.g[j],0)||1; return key.map(fn=>X.pts.filter(fn).reduce((a,p)=>a+p.g[j],0)/tot); });
      b+=`<div class="dre-h4" style="margin-top:10px">Which data determined each parameter — information contribution (rows sum to 100 %)</div><div class="dre-tw">${heat({rows:R.free.map(f=>f.key),cols,M,groups,cw:40,th:74,lw:84,label:"Information contribution",color:(v,P)=>({c:P.accent,a:Math.min(1,v*2.2)}),fmt:v=>(100*v).toFixed(1)+" %",text:v=>v>=0.1?Math.round(100*v):""})}</div>
        <div class="muted">열 = 데이터셋별 파워 구간 (${esc(S.pu)}). 점 i 가 매개변수 j 에 주는 Fisher 정보 J<sub>ij</sub>² 의 fraction (파워 단위 ${esc(S.pu)}). 한 구간에 몰려 있으면 그 구간의 점 몇 개가 값을 좌우한다 — 그 구간을 더 촘촘히 재거나 전환 중인 점이 섞였는지 확인. 상관이 강한 매개변수끼리는 이 fraction이 겹쳐 보인다(③ 참조).</div>`; }
    const profs=R.prof.length?R.free.map((f,j)=>{ const pr=R.prof[j]; if(!pr) return ""; return `<div class="dre-pf">${chart({id:"p"+j,xlog:false,ylog:false,series:[{x:pr.pts.map(p=>p.v),y:pr.pts.map(p=>p.chi2-R.best.chi2),points:true,label:f.key,color:S10(P)[j%10]}],xlabel:`log₁₀ ${f.key}`,ylabel:"Δχ²",w:260,h:170,nt:5,hlines:[{y:R.thr,label:"95 %"}],ymin:0,ymax:Math.max(R.thr*3,...pr.pts.map(p=>Math.min(p.chi2-R.best.chi2,R.thr*6))),legend:false})}</div>`; }).join(""):"";
    if(profs) b+=`<div class="dre-h4" style="margin-top:10px">Profile likelihood (Δχ² with the other parameters re-optimized)</div><div class="dre-pfs">${profs}</div>`;
    html+=sec("par","② Parameter evidence","Interval comparison · information contribution · profiles",b,false); }
  /* 3 종속성 */
  { let b="";
    if(X){ const names=R.free.map(f=>f.key);
      b+=`<div class="dre-grid2"><div><div class="dre-h4">Correlation matrix (curvature-based)</div><div class="dre-tw">${heat({rows:names,cols:names,M:X.corr,cw:38,ch:26,lw:84,th:74,label:"Correlation matrix",color:(v,P)=>({c:v>=0?P.accent:"#1f6fa3",a:Math.min(1,Math.abs(v))}),fmt:v=>v.toFixed(3),text:v=>Math.abs(v)>=0.1?v.toFixed(1):""})}</div>
        <div class="muted">|ρ| → 1 이면 두 매개변수는 따로 정해지지 않고 한 조합만 정해진다. 보라 = 같이 움직임, 파랑 = 서로 보상.</div></div>`;
      const ev=R.fim; const bar=(e,title)=>{ const W=300, lw=84, H=names.length*18+28, g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="100%"><rect width="${W}" height="${H}" fill="${P.bg}"/><line x1="${lw+(W-lw-10)/2}" x2="${lw+(W-lw-10)/2}" y1="4" y2="${H-20}" stroke="${P.muted}"/>`];
        e.vec.forEach((v,j)=>{ const y=6+j*18, x0=lw+(W-lw-10)/2, w=v*(W-lw-10)/2; g.push(`<text x="${lw-6}" y="${y+10}" font-size="10" text-anchor="end" fill="${P.ink2}">${esc(names[j])}</text><rect x="${Math.min(x0,x0+w)}" y="${y+2}" width="${Math.abs(w)}" height="11" fill="${v>=0?P.accent:"#1f6fa3"}" fill-opacity=".8"><title>${v.toFixed(3)}</title></rect>`); });
        g.push(`<text x="${W/2}" y="${H-5}" font-size="10" text-anchor="middle" fill="${P.muted}">${esc(title)}</text></svg>`); return g.join(""); };
      const nd=Math.min(2,ev.length);
      b+=`<div><div class="dre-h4">Eigenvectors — least / best determined combinations</div>${ev.slice(0,nd).map((e,i)=>bar(e,`Least determined #${i+1} · eigenvalue 10^${Math.log10(Math.max(Math.abs(e.val),1e-30)).toFixed(1)}`)).join("")}${bar(ev[ev.length-1],`Best determined · eigenvalue 10^${Math.log10(Math.max(Math.abs(ev[ev.length-1].val),1e-30)).toFixed(1)}`)}
        <div class="muted">덜 정해지는 방향의 막대가 여러 매개변수에 같은 부호로 퍼져 있으면 ‘모든 속도·단면적 동시 배율’ 같은 scale degeneracy다 — 절대 시간 데이터가 깬다.</div></div></div>`; }
    const pj=S.insProfJ!=null&&S.insProfJ<R.free.length?S.insProfJ:0;
    if(R.prof.length&&R.prof[pj]){ const pr=R.prof[pj]; const b0=R.best.theta;
      const series=R.free.map((f,k)=>k===pj?null:{x:pr.pts.map(p=>p.v-b0[pj]),y:pr.pts.map(p=>p.theta[k]-b0[k]),points:true,label:f.key,color:S10(P)[k%10]}).filter(Boolean);
      b+=`<div class="dre-h4" style="margin-top:10px">Compensation path — shifting <select class="dre-sel" id="insProfJ">${R.free.map((f,k)=>`<option value="${k}" ${k===pj?"selected":""}>${esc(f.key)}</option>`).join("")}</select> while the others re-optimize</div>
        ${chart({id:"ct",xlog:false,ylog:false,series,xlabel:`Δlog₁₀ ${R.free[pj].key} (imposed shift)`,ylabel:"Δlog₁₀ (re-optimized value)",w:640,h:260,legend:false,vlines:[{x:pr.lo-b0[pj]},{x:pr.hi-b0[pj]}]})}${legendHTML(series.map(q=>[q.label,q.color]))}
        <div class="muted">기울기 −1 로 따라가는 매개변수는 곱(또는 비)만 정해진다. 수평선이면 서로 독립. 점선 = 이 매개변수의 95 % 구간.</div>`; }
    html+=sec("dep","③ Dependence &amp; independence","Correlation · eigenvectors · compensation path",b||'<div class="muted">No data</div>',false); }
  /* 4 χ² landscape */
  { const wf=chart({id:"wf",xlog:false,ylog:true,series:[{x:R.ms.map((_,i)=>i+1),y:R.ms.map(r=>Math.max(r.chi2,1e-12)),points:true,line:true,width:1,label:"Final χ² per start",color:P.accent}],xlabel:"Start (sorted by χ²)",ylabel:"χ²",w:420,h:230,hlines:[{y:R.best.chi2+R.thr,label:"95 % threshold",color:P.muted}],legend:false});
    let b=`<div class="dre-grid2"><div><div class="dre-h4">waterfall</div>${wf}<div class="muted">평평한 첫 계단의 길이 = global minimum에 도달한 출발점 수. 계단이 여럿이면 local minimum가 여럿.</div></div>`;
    if(X){ b+=`<div><div class="dre-h4">Solution clusters (by χ² level)</div><div class="dre-tw"><table class="dre-tbl"><thead><tr><th class="num">Δχ²</th><th class="num">Starts</th>${R.free.map(f=>`<th class="num">${esc(f.key)} range (log₁₀)</th>`).join("")}</tr></thead><tbody>${X.clusters.slice(0,8).map(c=>`<tr><td class="num">${(c.chi2-R.best.chi2).toFixed(2)}</td><td class="num">${c.count}</td>${R.free.map((f,j)=>`<td class="num">${c.lo[j].toFixed(2)}${c.hi[j]-c.lo[j]>0.05?` – ${c.hi[j].toFixed(2)}`:""}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
      <div class="muted">같은 χ² 수준 안에서 값이 넓게 퍼지면 ‘local minimum’가 아니라 평평한 골짜기(non-identifiable)다.</div></div>`; }
    b+=`</div>`;
    /* 쌍 산점도 */
    const Pn=R.free.length; if(Pn>=2){ const fin=R.ms.filter(r=>isFinite(r.chi2)); const cls=r=>r.chi2-R.best.chi2<=R.thr?0:r.chi2-R.best.chi2<=10*R.thr?1:2; const cc=[P.accent,"#b4580f",P.muted], cl=["Global (Δχ² ≤ threshold)","Near (≤ 10×)","Distant local minimum"];
      const maxP=Math.min(Pn,6); const pairs=[]; for(let a=0;a<maxP;a++) for(let c=a+1;c<maxP;c++) pairs.push([a,c]);
      b+=`<div class="dre-h4" style="margin-top:10px">Endpoint distribution — parameter pairs (log₁₀)${Pn>6?" · first 6":""}</div><div class="muted">${cl.map((t,i)=>`<span style="color:${cc[i]}">●</span> ${t}`).join(" · ")} · line = best</div><div class="dre-pairs">${pairs.map(([a,c])=>`<div>${chart({id:`pr${a}_${c}`,xlog:false,ylog:false,series:[0,1,2].map(k=>({x:fin.filter(r=>cls(r)===k).map(r=>r.theta[a]),y:fin.filter(r=>cls(r)===k).map(r=>r.theta[c]),points:true,line:false,r:2.6,color:cc[k],label:cl[k]})).filter(s=>s.x.length),xlabel:R.free[a].key,ylabel:R.free[c].key,w:220,h:190,nt:4,legend:false,vlines:[{x:R.best.theta[a]}],hlines:[{y:R.best.theta[c]}]})}</div>`).join("")}</div>`; }
    /* χ² slice */
    const sl=S.insSlice&&S.insSlice.ver===R.time?S.insSlice:null;
    b+=`<div class="dre-h4" style="margin-top:10px">χ² slice (two-parameter grid, others fixed at best)</div><div class="row"><select class="dre-sel" id="slA">${R.free.map((f,k)=>`<option value="${k}" ${sl&&sl.j===k||(!sl&&k===0)?"selected":""}>${esc(f.key)}</option>`).join("")}</select> × <select class="dre-sel" id="slB">${R.free.map((f,k)=>`<option value="${k}" ${sl&&sl.k===k||(!sl&&k===1)?"selected":""}>${esc(f.key)}</option>`).join("")}</select>
       <label class="dre-f">± <input type="number" id="slS" value="${sl?sl.span:1.5}" step="0.5" min="0.25" max="4"> decade</label><button class="btn sm" id="slGo">Compute slice (17×17)</button><span class="muted" id="slMsg"></span></div>`;
    if(sl) b+=sliceSVG(sl,R);
    b+=`<div class="muted">slice는 ‘나머지 고정’이라 프로파일보다 좁게 보인다 — 골짜기의 방향(상관)과 여러 우물(local minimum)을 눈으로 확인하는 용도. 점 = multi-start 끝점의 투영.</div>`;
    html+=sec("land","④ χ² landscape","waterfall · solution clusters · endpoint distribution · χ² slice",b,false); }
  /* 5 메커니즘 */
  if(X){ const names=X.procs.map(p=>enL(p.name)); const cols=names.map((_,k)=>P.ser[k%6]+(k>=6?"":"")), dash=names.map((_,k)=>k>=6);
    const colors=names.map((_,k)=>S10(P)[k%10]);
    let b=`<div class="dre-h4">Per-process flux decomposition — fraction of each process among the transitions at that state</div>`;
    X.flux.forEach(fx=>{ if(!fx) return; b+=`<div class="muted" style="margin:4px 0">Sample ${esc(S.samples[fx.si]?S.samples[fx.si].id:"")}</div>`+stackBars({rows:fx.rows.map(r=>({label:`${sciTxt(r.I/PU[S.pu],3)} ${S.pu} · ${r.branch}`,v:r.share})),names,colors,label:"Per-process flux"}); });
    b+=`<div class="dre-legend">${names.map((n,k)=>`<span><i style="background:${colors[k]}"></i>${esc(n)} <span class="mono">${esc(X.procs[k].expr)}</span></span>`).join("")}</div>`;
    const maxSh=names.map((_,k)=>Math.max(0,...X.flux.filter(Boolean).flatMap(f=>f.rows.map(r=>r.share[k]))));
    const weak=names.map((n,k)=>({n,k,m:maxSh[k]})).filter(q=>q.m<0.01);
    if(weak.length){ const pn=new Set(weak.flatMap(q=>X.procs[q.k].names)); const free=R.free.filter(f=>pn.has(f.nm)).map(f=>f.key);
      b+=`<div class="note warn" style="margin-top:6px">데이터 범위 전체에서 fraction이 1 % 미만인 과정: ${weak.map(q=>esc(q.n)).join(", ")}. 이 과정의 계수${free.length?`(${free.map(esc).join(", ")})`:""}는 이 데이터로 정해질 수 없다 — ‘non-identifiable’ 판정이 나와도 이상한 것이 아니며, 물성으로 보고하면 안 된다.</div>`; }
    b+=`<div class="muted">fraction = 과정 속도(population × 계수 × photon flux) / 모든 과정 속도의 합. threshold 부근에서 CR·ESA fraction이 커지고 GSA fraction이 작아지는 것이 광자 사태의 전형적 서명이다.</div>`;
    html+=sec("mech","⑤ Mechanism","Per-process flux decomposition",b,false); }
  /* 예측값 */
  if(R.derived&&R.derived.length){ const f=PU[S.pu]; const hasProf=R.prof&&R.prof.length;
    const rowsD=R.derived.map(q=>{ const conv=q.unit==="W/cm²"; const sc=conv?1/f:1; const u=conv?S.pu:q.unit; const w=q.lo&&q.hi?Math.log10(q.hi/q.lo):null;
      const gr=w==null||!hasProf?"":w<=0.1?'<span class="gr g0">Determined</span>':w<=0.5?'<span class="gr g2">Broad</span>':'<span class="gr g3">Spread</span>';
      return `<tr><td>${esc(enL(q.label))}${q.calFree?' <span class="gr g2" title="nominal intensity 기준 속도 — intensity 환산 오차가 σ 와 I 에서 상쇄된다">Cal-independent</span>':""}${q.calNote?' <span class="muted" title="intensity 환산 오차가 그대로 들어간다">Cal-dependent</span>':""}</td><td class="num">${sci(q.v*sc,4)} ${esc(u)}</td><td class="num">${q.lo!=null&&hasProf?`${sci(q.lo*sc,3)} – ${sci(q.hi*sc,3)}`:"—"}</td><td class="num">${w!=null&&hasProf?w.toFixed(3):"—"}</td><td>${gr}</td></tr>`; }).join("");
    const b=`<div class="dre-tw"><table class="dre-tbl"><thead><tr><th>Quantity</th><th class="num">Best</th><th class="num">95 % range (profile ensemble)</th><th class="num">Width (decade)</th><th>Status</th></tr></thead><tbody>${rowsD}</tbody></table></div>
      <div class="muted">범위 = 프로파일 경로 위에서 Δχ² ≤ threshold인 매개변수 집합 전부로 계산한 최소–최대${hasProf?"":" (프로파일을 켜야 계산된다)"}. I_ref = 데이터 파워 범위의 기하 평균(${sciTxt(R.derived[0].Iref/f,3)} ${esc(S.pu)}).
      매개변수가 ‘non-identifiable’이어도 이 양들은 좁게 정해지는 경우가 많다 — 합성 참값 시험(Lee 3준위, intensity 환산 오차 ±0.3 decade 포함)에서 S-curve + 저파워 감쇠 1점만으로 ESA 속도·CR fraction의 오차 중앙값이 0.01 decade 이하, σ_ESA 자체는 0.26 decade(환산 오차 그대로)였다. 보고할 때는 σ 대신 ‘nominal intensity 기준 속도’를 쓰는 편이 정확하다.</div>`;
    html+=sec("pred","⑦ Predictions — quantities better determined than the parameters","Absorption rate · transfer fraction · threshold · population",b,false); }
  /* 6 residual */
  if(X){ const sets=X.sets.filter(s=>s.type==="scurve"); const f=PU[S.pu];
    const series=sets.map((st,k)=>{ const pts=X.pts.filter(p=>p.si===st.si&&p.di===st.di&&!p.cens); return {x:pts.map(p=>p.I/f),y:pts.map(p=>p.r),points:true,line:false,r:2.6,label:enL(st.label),color:P.ser[k%6]}; }).filter(s=>s.x.length);
    const rr=X.pts.filter(p=>!p.cens&&p.r!==0).map(p=>p.r).sort((a,b)=>a-b); const q=rr.map((_,i)=>normInv((i+0.5)/rr.length));
    let b=`<div class="dre-grid2"><div>${series.length?chart({id:"rs",xlog:true,ylog:false,series,xlabel:`excitation intensity (${S.pu})`,ylabel:"Residual (σ units)",w:420,h:250,legend:false,hlines:[{y:0},{y:3,label:"+3σ"},{y:-3,label:"−3σ"}]})+legendHTML(series.map(q=>[q.label,q.color])):""}<div class="muted">무리 지은 같은 부호 = 모양 차이(적합성 A2). 한 점만 튀면 전환 중인 점일 수 있다.</div></div>
      <div>${rr.length>3?chart({id:"qq",xlog:false,ylog:false,series:[{x:q,y:rr,points:true,line:false,r:2.4,label:"residual",color:P.accent},{x:[q[0],q[q.length-1]],y:[q[0],q[q.length-1]],label:"Normal",color:P.muted,dash:"4 3"}],xlabel:"normal quantile",ylabel:"residual quantile",w:320,h:250,legend:false}):""}<div class="muted">정규 Q-Q: 직선을 따르면 잡음 모델이 맞다. 기울기 &gt; 1 이면 잡음 과소 입력 또는 모델 부족, 꼬리만 벗어나면 이상점.</div></div></div>`;
    html+=sec("res","⑥ Residuals","Residual vs. power · normal Q–Q",b,false); }
  el.innerHTML=html;
  el.querySelectorAll("details.dre-ins").forEach(d=>d.addEventListener("toggle",()=>{ S.insOpen=S.insOpen||{}; S.insOpen[d.dataset.k]=d.open; save(); }));
  const pj=$q(el,"#insProfJ"); if(pj) pj.onchange=()=>{ S.insProfJ=+pj.value; save(); renderDiag(); };
  const go=$q(el,"#slGo"); if(go) go.onclick=()=>runSlice(+$q(el,"#slA").value,+$q(el,"#slB").value,+$q(el,"#slS").value||1.5);
}
function sliceSVG(sl,R){
  const P=palette(), n=sl.xs.length, W=440, H=380, m={l:62,r:70,t:12,b:44}, pw=W-m.l-m.r, ph=H-m.t-m.b; const c0=R.best.chi2;
  const dz=sl.z.map(r=>r.map(v=>Math.max(0,v-c0))); const mx=Math.max(...dz.flat().filter(isFinite));
  const col=v=>{ if(!isFinite(v)) return {c:P.muted,a:.2}; const t=Math.log10(1+v)/Math.log10(1+mx); return {c:v<=R.thr?P.accent:"#1f6fa3",a:v<=R.thr?0.85:0.1+0.75*(1-t)}; };
  const cw=pw/n, chh=ph/n; const g=[`<svg class="dre-chart" viewBox="0 0 ${W} ${H}" width="100%" style="max-width:520px" role="img" aria-label="χ² slice"><rect width="${W}" height="${H}" fill="${P.bg}"/>`];
  dz.forEach((row,b)=>row.forEach((v,a)=>{ const f=col(v); g.push(`<rect x="${(m.l+a*cw).toFixed(1)}" y="${(m.t+ph-(b+1)*chh).toFixed(1)}" width="${(cw+0.5).toFixed(1)}" height="${(chh+0.5).toFixed(1)}" fill="${f.c}" fill-opacity="${f.a}"><title>${R.free[sl.j].key} ${sl.xs[a].toFixed(2)}, ${R.free[sl.k].key} ${sl.ys[b].toFixed(2)}: Δχ² ${isFinite(v)?v.toFixed(2):"failed"}</title></rect>`); }));
  const X=v=>m.l+(v-sl.xs[0])/(sl.xs[n-1]-sl.xs[0])*pw, Y=v=>m.t+ph-(v-sl.ys[0])/(sl.ys[n-1]-sl.ys[0])*ph;
  R.ms.forEach(r=>{ const x=r.theta[sl.j], y=r.theta[sl.k]; if(x<sl.xs[0]||x>sl.xs[n-1]||y<sl.ys[0]||y>sl.ys[n-1]) return; g.push(`<circle cx="${X(x)}" cy="${Y(y)}" r="2.4" fill="none" stroke="${P.ink}" stroke-width=".8"/>`); });
  g.push(`<path d="M${X(R.best.theta[sl.j])-6} ${Y(R.best.theta[sl.k])}h12M${X(R.best.theta[sl.j])} ${Y(R.best.theta[sl.k])-6}v12" stroke="${P.ink}" stroke-width="1.6"/>`);
  g.push(`<rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}" fill="none" stroke="${P.line}"/>`);
  [0,Math.floor(n/2),n-1].forEach(i=>{ g.push(`<text x="${X(sl.xs[i])}" y="${m.t+ph+14}" font-size="10" text-anchor="middle" fill="${P.muted}">${sl.xs[i].toFixed(2)}</text><text x="${m.l-5}" y="${Y(sl.ys[i])+3}" font-size="10" text-anchor="end" fill="${P.muted}">${sl.ys[i].toFixed(2)}</text>`); });
  g.push(`<text x="${m.l+pw/2}" y="${H-8}" font-size="11" text-anchor="middle" fill="${P.ink2}">log₁₀ ${esc(R.free[sl.j].key)}</text><text x="14" y="${m.t+ph/2}" font-size="11" text-anchor="middle" fill="${P.ink2}" transform="rotate(-90 14 ${m.t+ph/2})">log₁₀ ${esc(R.free[sl.k].key)}</text>`);
  g.push(`<rect x="${W-m.r+12}" y="${m.t+4}" width="12" height="12" fill="${P.accent}" fill-opacity=".85"/><text x="${W-m.r+28}" y="${m.t+14}" font-size="9.5" fill="${P.ink2}">Δχ²≤threshold</text><rect x="${W-m.r+12}" y="${m.t+24}" width="12" height="12" fill="#1f6fa3" fill-opacity=".8"/><text x="${W-m.r+28}" y="${m.t+34}" font-size="9.5" fill="${P.ink2}">Low</text><rect x="${W-m.r+12}" y="${m.t+44}" width="12" height="12" fill="#1f6fa3" fill-opacity=".15"/><text x="${W-m.r+28}" y="${m.t+54}" font-size="9.5" fill="${P.ink2}">High</text>`);
  g.push("</svg>"); return g.join(""); }
function runSlice(j,k,span){
  if(j===k){ note("서로 다른 두 매개변수를 고르세요","warn"); return; } const R=S.result; if(!R) return;
  const src=workerSource(); if(!src) return; const n=17; const th=R.best.theta;
  const xs=[...Array(n).keys()].map(i=>th[j]-span+2*span*i/(n-1)), ys=[...Array(n).keys()].map(i=>th[k]-span+2*span*i/(n-1));
  const K=Math.min(nWorkers(),n); const url=URL.createObjectURL(new Blob([src],{type:"text/javascript"})); const ws=[...Array(K)].map(()=>new Worker(url)); URL.revokeObjectURL(url);
  const prob=problemOf(), pid="s"+Date.now(); const z=new Array(n); let next=0, done=0;
  const msg=()=>$q(ROOT,"#slMsg"); if(msg()) msg().textContent=`계산 중… (워커 ${K}개)`;
  const call=(w,m)=>new Promise((res,rej)=>{ w.onmessage=ev=>{ const d=ev.data; if(d.type==="slicerow") res(d); else if(d.type==="error") rej(new Error(d.msg)); }; w.onerror=e=>rej(new Error(e.message||"")); w.postMessage(m); });
  Promise.all(ws.map(async w=>{ while(next<n){ const b=next++; const d=await call(w,{cmd:"slicerow",pid,problem:prob,theta:th,j,k,xs,y:ys[b],b}); z[d.b]=d.row; done++; const e=msg(); if(e) e.textContent=`${done}/${n} 행 · 워커 ${K}개`; } }))
    .then(()=>{ ws.forEach(w=>w.terminate()); S.insSlice={xs,ys,z,j,k,span,ver:R.time}; S.insOpen=Object.assign(S.insOpen||{},{land:true}); save(); renderDiag(); })
    .catch(e=>{ ws.forEach(w=>w.terminate()); note("slice 계산 오류: "+e.message,"warn"); });
}
function adequacyNow(){
  try{ const B=D.buildProblem(JSON.parse(JSON.stringify(problemOf()))); const f=PU[S.pu];
    S.adeqNow=Object.assign(D.adequacy(B,D.thetaFromParams(B),{fmtI:v=>sciTxt(v/f,3)+" "+S.pu}),{at:Date.now()}); save(); renderDiag(); const el=$q(ROOT,"#dreDiag"); if(el&&el.scrollIntoView) el.scrollIntoView({behavior:"smooth",block:"start"}); note("적합성 검사 완료 — 현재 매개변수 기준"); }
  catch(e){ note("적합성 검사 실패: "+e.message,"warn"); }
}

/* ---------------- 내보내기 ---------------- */
function download(name,blob){ const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },400); }
function exportSimCSV(){
  if(!SIM||SIM.err) return; const f=PU[S.pu];
  const hdr=["I_"+S.pu].concat(MODEL.channels.map(c=>c.id+" (up)"),MODEL.channels.map(c=>c.id+" (down)"),MODEL.lv.map(l=>`n_${l.sp}_${l.id}`));
  const rows=SIM.I.map((I,i)=>[I/f].concat(SIM.up.map(a=>a[i]),SIM.dn.map(a=>a[i]),Array.from(SIM.pop[i])));
  const csv="﻿# "+(S.name||"")+" | lambda "+S.nm+" nm | beam "+S.beam+"\n"+[hdr].concat(rows).map(r=>r.map(v=>typeof v==="number"?v.toPrecision(8):'"'+String(v).replace(/"/g,'""')+'"').join(",")).join("\n");
  download("DRE_simulation.csv",new Blob([csv],{type:"text/csv"}));
}
function exportJSON(){ download("DRE_model.json",new Blob([JSON.stringify(S,null,1)],{type:"application/json"})); }
function importJSON(f){ const r=new FileReader(); r.onload=()=>{ try{ const o=JSON.parse(r.result); if(!o.spec||!o.params) throw new Error("형식이 다릅니다"); D.compile(o.spec); S=Object.assign(defState(),o); save(); rerender(); note("불러왔습니다"); }catch(e){ note("불러오기 실패: "+e.message,"warn"); } }; r.readAsText(f); }

window.addEventListener("resize",()=>{ clearTimeout(window._dreRz); window._dreRz=setTimeout(()=>{ if(ROOT&&ROOT.isConnected) drawSim(); },200); });
window.DRE_UI={render,rerender,_state:()=>S,_model:()=>MODEL,_sim:()=>SIM,_setState:o=>{ S=Object.assign(defState(),o); save(); rerender(); }};
})();
