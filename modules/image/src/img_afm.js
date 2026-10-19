/* ============================================================
   INNO Analytics — AFM·SPM 형식 읽기 (IMG_IO 에 덧붙임)
   · Gwyddion GWY (GWYP 직렬화): GwyContainer 안 GwyDataField(/N/data)·mask(/N/mask)·GwyBrick(/brick/N), 제목·SI 단위·meta.
     그래프(1D)·XYZ 점 자료·curve map·spectra 는 이미지가 아니므로 읽지 않고 알린다. GWYO(Gwyddion 1.x)는 읽지 않는다.
   · Gwyddion Simple Field (GSF): 머리(Key = Value) + NUL 채움(4 바이트 경계) + float32 LE.
   · Bruker/Veeco Nanoscope (.spm, .000–.999, .pfc): ASCII 머리(\*File list … \*File list end) + 16/32-bit 정수 자료.
     값 = 원시 정수 / 256^(Bytes/pixel) × @2:Z scale 의 hard value(V) × soft scale(@Sens. …, 예 nm/V).
     해석 규칙(머리 길이 찾기, 9.2 이후 32-bit, Ciao scan list 의 전역 해상도, 종횡비, 위아래 뒤집기)은 Gwyddion 2.64
     modules/file/nanoscope.c 의 동작을 따랐다 — 시험도 Gwyddion 이 같은 파일을 읽은 결과와 대조한다.
     힘 곡선 파일(\*Force file list)·텍스트형(?*File list)은 읽지 않는다(분명한 오류). Force Volume 은 이미지 채널만.
   · Asylum Research ARDF (force map·fast force mapping): FTOC → IMAG(TTOC·IDEF·IBOX·IDAT·GAMI) 이미지 채널.
     IDEF 의 간격·단위·이름을 쓴다(Gwyddion 은 FastScanSize 와 제목으로 추정 — 정사각이 아니면 틀림).
     이미지는 Igor wave 처럼 첫 줄이 아래(IBW 와 같은 방향) — magic-afm 도 origin="lower" 로 그린다. 힘 곡선(VOLM)은 알리기만.
   값 단위: 길이(m)는 nm 로, N·Pa·A·F·S 는 값 크기에 맞는 접두어로(예 MPa, nN) 바꾸고 바꾼 배율을 정보에 적는다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IO=G.IMG_IO; if(!IO) throw new Error("IMG_IO 가 먼저 있어야 합니다");
const tdl=new TextDecoder("latin1"), tdu=new TextDecoder("utf-8");
let tdw; try{ tdw=new TextDecoder("windows-1252"); }catch(e){ tdw=tdl; }

/* ---------------- SI 단위 (Gwyddion gwysiunit 의 관례: ~·u·µ·μ = 마이크로, 접두어와 겹치는 단위 deg·Pa·cd·mol 은 통째로) ---------------- */
const PFX={y:-24,z:-21,a:-18,f:-15,p:-12,n:-9,u:-6,"~":-6,"µ":-6,"μ":-6,m:-3,c:-2,k:3,M:6,G:9,T:12,P:15,E:18,Z:21,Y:24};
const KNOWN=new Set(["m","g","s","A","K","mol","cd","Hz","N","Pa","J","W","C","V","F","Ω","S","Wb","T","H","deg","rad","px","cts","cps","eV","L","Gy","Gauss","pt","cal"]);
function unitParse(str){
  const out={p10:0,u:{}}; let s=String(str==null?"":str).trim(); if(!s) return out;
  s=s.replace(/^microns?$/i,"µm").replace(/^[°\u00ba]$/,"deg").replace(/µ/g,"µ");
  const parts=s.split("/");
  parts.forEach((part,pi)=>{ const sign=pi===0?1:-1;
    for(let tok of part.split(/[\s*·]+/).filter(Boolean)){
      let pw=1; const m=tok.match(/^(.*?)\^(-?\d+)$/); if(m){ tok=m[1]; pw=+m[2]; }
      if(/^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?$/.test(tok)){ const q=Math.round(Math.log10(+tok)); if(+tok>0&&Math.abs(+tok/Math.pow(10,q)-1)<1e-12) out.p10+=sign*pw*q; continue; } // "1e-9" 같은 수 배율
      let sym=tok, p=0;
      if(tok==="Å"||tok==="Ang"||tok==="AA"){ sym="m"; p=-10; }
      else if(tok==="°"||tok==="\u00ba"||tok==="Deg"||/^degrees?$/.test(tok)) sym="deg";
      else if(!KNOWN.has(tok)&&tok.length>1&&PFX[tok[0]]!=null&&KNOWN.has(tok.slice(1))){ p=PFX[tok[0]]; sym=tok.slice(1); }
      else if(!KNOWN.has(tok)&&tok.length>2&&tok.startsWith("da")&&KNOWN.has(tok.slice(2))){ p=1; sym=tok.slice(2); }
      out.p10+=sign*pw*p; out.u[sym]=(out.u[sym]||0)+sign*pw; if(out.u[sym]===0) delete out.u[sym];
    } });
  return out;
}
function unitMul(a,b){ const o={p10:a.p10+b.p10,u:Object.assign({},a.u)}; for(const k in b.u){ o.u[k]=(o.u[k]||0)+b.u[k]; if(o.u[k]===0) delete o.u[k]; } return o; }
function unitStr(o){ const pos=[],neg=[]; for(const k in o.u){ const p=o.u[k]; if(p>0) pos.push(k+(p>1?"^"+p:"")); else neg.push(k+(p<-1?"^"+(-p):"")); } return (pos.join(" ")||(neg.length?"1":""))+(neg.length?"/"+neg.join(" "):""); }
const isOnly=(o,k)=>{ const ks=Object.keys(o.u); return ks.length===1&&ks[0]===k&&o.u[k]===1; };
/* 값(기본 SI 단위) → 표시 단위. m 은 늘 nm, N·Pa·A·F·S 는 중앙값 크기로 접두어. 그 밖은 그대로 */
const PFXN={"-15":"f","-12":"p","-9":"n","-6":"µ","-3":"m","0":"","3":"k","6":"M","9":"G","12":"T"};
function dispUnit(o,data){
  const base=unitStr(o);
  if(isOnly(o,"m")) return {unit:"nm",f:1e9,base};
  if(["N","Pa","A","F","S"].some(k=>isOnly(o,k))){
    const n=data.length, st=Math.max(1,Math.floor(n/20000)); const a=[]; for(let i=0;i<n;i+=st){ const v=Math.abs(data[i]); if(v>0&&isFinite(v)) a.push(v); }
    if(!a.length) return {unit:base,f:1,base}; a.sort((x,y)=>x-y); const med=a[a.length>>1];
    let e=Math.floor(Math.log10(med)/3)*3; e=Math.max(-15,Math.min(12,e)); return {unit:PFXN[e]+base,f:Math.pow(10,-e),base};
  }
  return {unit:base,f:1,base};
}
/* 가로·세로 보정: 길이면 calFromMeters, 그 밖 단위는 그 단위 그대로(예 s), 단위가 없으면 보정 없음 */
function calOf(dx,dy,uo,src){
  if(!(dx>0)||!isFinite(dx)) return null;
  const k=Math.pow(10,uo.p10);
  if(isOnly(uo,"m")) return IO.calFromMeters(dx*k,(dy>0?dy:dx)*k,src);
  if(!Object.keys(uo.u).length) return null;
  return {pw:dx*k,ph:(dy>0?dy:dx)*k,unit:unitStr(uo),src:src+" (길이 단위 아님)"};
}
const g6=v=>(+v).toPrecision(6).replace(/\.?0+(e|$)/,"$1");

/* ============================================================ Gwyddion GWY ============================================================ */
function gwyParse(u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength), N=u8.length;
  const bad=(o,what)=>{ throw new Error(`GWY: 손상된 파일 — ${what} (바이트 ${o})`); };
  function cstr(o,end){ let e=o; while(e<end&&u8[e]!==0) e++; if(e>=end) bad(o,"끝나지 않은 문자열"); return [IO._txt(u8.subarray(o,e)),e+1]; } /* UTF-8, 아니면 latin1 (GWYO 의 µ = 0xB5) */
  function obj(o,end){
    const [name,a]=cstr(o,end); if(a+4>end) bad(a,"객체 크기"); const size=dv.getUint32(a,true); const s=a+4, e=s+size; if(e>end) bad(a,`객체 '${name}' 가 파일보다 김`);
    const c={}; let p=s;
    while(p<e){ const [key,q]=cstr(p,e); if(q>=e) bad(q,"성분 형식"); const t=String.fromCharCode(u8[q]); let r=q+1, v;
      const need=n=>{ if(r+n>e) bad(r,`성분 '${key}' 가 객체보다 김`); };
      switch(t){
        case "b": need(1); v=u8[r]!==0; r+=1; break;
        case "c": need(1); v=u8[r]; r+=1; break;
        case "i": need(4); v=dv.getInt32(r,true); r+=4; break;
        case "q": need(8); v=dv.getUint32(r,true)+dv.getInt32(r+4,true)*4294967296; r+=8; break;
        case "d": need(8); v=dv.getFloat64(r,true); r+=8; break;
        case "s": [v,r]=cstr(r,e); break;
        case "o": { const ob=obj(r,e); v=ob.v; r=ob.end; break; }
        case "C": case "I": case "Q": case "D": case "S": case "O": {
          need(4); const n=dv.getUint32(r,true); r+=4; const sz={C:1,I:4,Q:8,D:8}[t];
          if(sz){ need(n*sz); v={arr:t,off:r,n}; r+=n*sz; }
          else if(t==="S"){ v=[]; for(let i=0;i<n;i++){ let x; [x,r]=cstr(r,e); v.push(x); } }
          else { v=[]; for(let i=0;i<n;i++){ const ob=obj(r,e); v.push(ob.v); r=ob.end; } }
          break; }
        default: bad(q,`알 수 없는 성분 형식 '${t}' (성분 '${key}')`);
      }
      c[key]=v; p=r; }
    return {v:{_t:name,c},end:e};
  }
  /* Gwyddion 1.x (GWYO): 맨 위 GwyContainer 만 옛 직렬화 — 항목마다 int32 GType, 이름, 값 (GLib 기본형 번호). 안의 객체는 새 형식과 같다 */
  function objOld(o,end){
    const [name,a]=cstr(o,end); if(name!=="GwyContainer") bad(o,"GWYO 맨 위가 GwyContainer 가 아님"); if(a+4>end) bad(a,"객체 크기"); const size=dv.getUint32(a,true); const s=a+4, e=s+size; if(e>end) bad(a,"GwyContainer 가 파일보다 김");
    const c={}; let p=s;
    while(p<e){ if(p+4>e) bad(p,"GType"); const ty=dv.getInt32(p,true); let key; [key,p]=cstr(p+4,e); let v;
      const need=n=>{ if(p+n>e) bad(p,`항목 '${key}' 가 객체보다 김`); };
      switch(ty){
        case 80: { const ob=obj(p,e); v=ob.v; p=ob.end; break; }   /* G_TYPE_OBJECT */
        case 20: case 16: case 12: need(1); v=u8[p]; p+=1; if(ty===20) v=v!==0; break;   /* BOOLEAN·UCHAR·CHAR */
        case 24: need(4); v=dv.getInt32(p,true); p+=4; break;   /* INT */
        case 40: need(8); v=dv.getUint32(p,true)+dv.getInt32(p+4,true)*4294967296; p+=8; break;   /* INT64 */
        case 60: need(8); v=dv.getFloat64(p,true); p+=8; break;   /* DOUBLE */
        case 64: [v,p]=cstr(p,e); break;   /* STRING */
        default: bad(p,`GWYO 항목 '${key}' 의 모르는 GType ${ty}`);
      }
      c[key]=v; }
    return {_t:"GwyContainer",c,old:true};
  }
  const magic=tdl.decode(u8.subarray(0,4));
  if(magic!=="GWYP"&&magic!=="GWYO") throw new Error("GWY 머리(GWYP·GWYO)가 아닙니다");
  const top=magic==="GWYO"?objOld(4,N):obj(4,N).v;
  return {top,f64:(a,i)=>dv.getFloat64(a.off+8*i,true)};
}
const gstr=o=>o&&o._t==="GwySIUnit"?(o.c.unitstr||""):"";
/* Gwyddion 은 GwySIUnit·GSF 단위 문자열의 접두어를 버린다(gwy_si_unit_new: power10 무시) — 값·크기는 늘 기본 SI 단위로 저장돼 있다 */
const unitBase=s=>({p10:0,u:unitParse(s).u});
function gwyOldMeta(C,id){ /* GWYO: meta 가 '/meta/…'·'/N/data/meta/…' 문자열 열쇠로 흩어져 있음 (Gwyddion gwyfile_pack_metadata 와 같은 자리) */
  const pre=[`/${id}/data/meta/`,"/meta/"]; return Object.keys(C).filter(k=>pre.some(p=>k.startsWith(p))&&typeof C[k]==="string").map(k=>`${k.replace(/^.*\/meta\//,"")}: ${C[k]}`).join("\n"); }
function gwyMetaText(m){ if(!m||m._t!=="GwyContainer") return ""; return Object.keys(m.c).map(k=>`${k}: ${typeof m.c[k]==="object"?"(…)":m.c[k]}`).join("\n"); }
function readGWY(name,u8){
  const {top,f64}=gwyParse(u8); const fileWarn=top.old?["Gwyddion 1.x 형식(GWYO) 파일입니다"]:[]; if(top._t!=="GwyContainer") throw new Error("GWY 최상위가 GwyContainer 가 아닙니다 ("+top._t+")");
  const C=top.c, keys=Object.keys(C), imgs=[], warn=fileWarn; const ids=re=>keys.map(k=>{ const m=k.match(re); return m?+m[1]:null; }).filter(v=>v!=null).sort((a,b)=>a-b);
  const field=(df,label,info,meta,isMask)=>{
    const c=df.c, W=c.xres, H=c.yres; const D=c.data;
    if(!(W>0&&H>0)||!D||D.arr!=="D"||D.n!==W*H) throw new Error(`GWY: '${label}' 의 자료 크기가 맞지 않습니다 (${W}×${H}, 값 ${D&&D.n})`);
    if(W*H*4>IO.MAXBYTES) throw new Error(`GWY: '${label}' 가 너무 큼 (${W}×${H})`);
    const uz=unitBase(gstr(c.si_unit_z)), uxy=unitBase(gstr(c.si_unit_xy));
    const raw=new Float64Array(W*H); for(let i=0;i<W*H;i++) raw[i]=f64(D,i);
    const du=isMask?{unit:"",f:1,base:""}:dispUnit({p10:0,u:uz.u},raw); const p=new Float32Array(W*H);
    if(isMask){ for(let i=0;i<W*H;i++) p[i]=raw[i]>0.5?255:0; }
    else for(let i=0;i<W*H;i++) p[i]=raw[i]*du.f;
    const xr=c.xreal!=null?c.xreal:1, yr=c.yreal!=null?c.yreal:1;
    const cal=calOf(xr/W,yr/H,uxy,"Gwyddion xreal/xres");
    const inf=[["Channel",label],...info]; if(!isMask) inf.push(["Value unit",du.unit?(du.f!==1?`${du.unit} (파일 ${du.base} × ${g6(du.f)})`:du.unit):"(없음)"]);
    if(c.xoff||c.yoff) inf.push(["Offset (x, y)",`${g6(c.xoff||0)}, ${g6(c.yoff||0)} ${unitStr(uxy)}`]);
    return {name:`${name} · ${label}`,w:W,h:H,bit:isMask?8:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal,calCandidates:[],info:inf,warn:[],valueUnit:isMask?"":du.unit,siScale:isMask?1:du.f,siUnit:isMask?"":du.base,metaText:meta};
  };
  for(const id of ids(/^\/(\d+)\/data$/)){
    const df=C[`/${id}/data`]; if(!df||df._t!=="GwyDataField"){ warn.push(`/${id}/data 가 GwyDataField 가 아닙니다 — 건너뜀`); continue; }
    const title=C[`/${id}/data/title`]||`Channel ${id}`; const meta=top.old?gwyOldMeta(C,id):gwyMetaText(C[`/${id}/meta`]);
    imgs.push(field(df,title,[["Gwyddion id",String(id)]],meta,false));
    const mk=C[`/${id}/mask`]; if(mk&&mk._t==="GwyDataField") imgs.push(field(mk,title+" · mask",[["Gwyddion id",id+" (mask)"]],"",true));
  }
  for(const id of ids(/^\/brick\/(\d+)$/)){
    const b=C[`/brick/${id}`]; if(!b||b._t!=="GwyBrick"){ warn.push(`/brick/${id} 가 GwyBrick 이 아닙니다 — 건너뜀`); continue; }
    const c=b.c, W=c.xres, H=c.yres, Z=c.zres, D=c.data; const title=C[`/brick/${id}/title`]||`Volume ${id}`;
    if(!(W>0&&H>0&&Z>0)||!D||D.arr!=="D"||D.n!==W*H*Z) throw new Error(`GWY: '${title}' 의 자료 크기가 맞지 않습니다`);
    if(W*H*Z*4>IO.MAXBYTES) throw new Error(`GWY: volume '${title}' 가 너무 큼 (${W}×${H}×${Z})`);
    const uw=unitBase(gstr(c.si_unit_w)), uxy=unitBase(gstr(c.si_unit_x)), uzz=unitBase(gstr(c.si_unit_z));
    const raw=new Float64Array(W*H*Z); for(let i=0;i<raw.length;i++) raw[i]=f64(D,i);
    const du=dispUnit({p10:0,u:uw.u},raw); const planes=[]; for(let z=0;z<Z;z++){ const p=new Float32Array(W*H); for(let i=0;i<W*H;i++) p[i]=raw[z*W*H+i]*du.f; planes.push(p); }
    const cal=calOf((c.xreal!=null?c.xreal:1)/W,(c.yreal!=null?c.yreal:1)/H,uxy,"Gwyddion brick xreal/xres");
    const zr=(c.zreal!=null?c.zreal:1)/Z; const inf=[["Volume",title],["Value unit",du.f!==1?`${du.unit} (파일 ${du.base} × ${g6(du.f)})`:du.unit||"(없음)"],["Z axis",`${Z} 단계, 간격 ${g6(zr)} ${unitStr(uzz)||"(단위 없음)"}`]];
    if(c.calibration) inf.push(["Z calibration","고르지 않은 z 눈금(calibration)이 있습니다 — z 간격은 평균값"]);
    const img={name:`${name} · ${title}`,w:W,h:H,bit:32,rgb:false,nc:1,nz:Z,nt:1,order:"czt",planes,cal,calCandidates:[],info:inf,warn:[],valueUnit:du.unit,siScale:du.f,siUnit:du.base,metaText:gwyMetaText(C[`/brick/${id}/meta`])};
    if(isOnly(uzz,"m")) img.zStep=zr*1e6;
    imgs.push(img);
  }
  const cnt=re=>keys.filter(k=>re.test(k)).length;
  const ng=cnt(/^\/\d+\/graph\/graph\/\d+$/)+cnt(/^\/graph\/\d+$/), nx=cnt(/^\/xyz\/\d+$/), nl=cnt(/^\/lawn\/\d+$/), ns=cnt(/^\/sps\/\d+$/)+cnt(/^\/\d+\/spectra$/);
  if(ng) warn.push(`그래프(1D) ${ng} 개는 이미지가 아니므로 읽지 않았습니다`); if(nx) warn.push(`XYZ 점 자료 ${nx} 개는 격자가 아니므로 읽지 않았습니다`);
  if(nl) warn.push(`curve map(lawn) ${nl} 개는 읽지 않았습니다`); if(ns) warn.push(`분광(spectra) ${ns} 묶음은 읽지 않았습니다`);
  if(!imgs.length) throw new Error("GWY 안에 이미지(GwyDataField·GwyBrick)가 없습니다"+(warn.length?" — "+warn.join(", "):""));
  return {images:imgs,warn};
}

/* ============================================================ Gwyddion Simple Field ============================================================ */
const GSF_MAGIC="Gwyddion Simple Field 1.0\n";
function readGSF(name,u8){
  if(tdl.decode(u8.subarray(0,GSF_MAGIC.length))!==GSF_MAGIC) throw new Error("GSF 머리가 아닙니다");
  let z=GSF_MAGIC.length; while(z<u8.length&&u8[z]!==0) z++; if(z>=u8.length) throw new Error("GSF: 머리 끝(NUL)을 찾지 못했습니다");
  const hdr=tdu.decode(u8.subarray(GSF_MAGIC.length,z)); const H={}; hdr.split(/\r?\n/).forEach(l=>{ const i=l.indexOf("="); if(i>0) H[l.slice(0,i).trim()]=l.slice(i+1).trim(); });
  const W=parseInt(H.XRes,10), Hh=parseInt(H.YRes,10); if(!(W>0&&Hh>0)) throw new Error("GSF: XRes·YRes 가 없습니다");
  const off=z+4-(z%4); if(off+W*Hh*4!==u8.length) throw new Error(`GSF: 자료 크기가 맞지 않습니다 (기대 ${off+W*Hh*4} 바이트, 파일 ${u8.length})`);
  const dv=new DataView(u8.buffer,u8.byteOffset+off,W*Hh*4); const uz=unitBase(H.ZUnits), uxy=unitBase(H.XYUnits);
  const raw=new Float64Array(W*Hh); for(let i=0;i<W*Hh;i++) raw[i]=dv.getFloat32(4*i,true);
  const du=dispUnit({p10:0,u:uz.u},raw); const p=new Float32Array(W*Hh); for(let i=0;i<W*Hh;i++) p[i]=raw[i]*du.f;
  const num=k=>{ const v=parseFloat(H[k]); return isFinite(v)&&v>0?v:1; };
  const cal=calOf(num("XReal")/W,num("YReal")/Hh,uxy,"GSF XReal/XRes");
  const nn=p.reduce((a,v)=>a+(v!==v?1:0),0); const info=[["Title",H.Title||""],["Value unit",du.f!==1?`${du.unit} (파일 ${du.base} × ${g6(du.f)})`:du.unit||"(없음)"]];
  return {images:[{name:H.Title?`${name} · ${H.Title}`:name,w:W,h:Hh,bit:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal,calCandidates:[],info,warn:nn?[`NaN ${nn} 화소 (Gwyddion 은 이 화소를 mask 로 봄)`]:[],valueUnit:du.unit,siScale:du.f,siUnit:du.base,metaText:hdr}],warn:[]};
}

/* ============================================================ Bruker Nanoscope ============================================================ */
// 값 한 줄 해석 — Gwyddion nanoscope.c parse_value 와 같은 규칙
function nsValue(key,line){
  const v={type:"old",soft:null,hscale:0,hscaleU:null,hv:NaN,hvStr:line,hvU:null};
  const strtod=s=>{ const m=s.match(/^\s*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/); return m?[parseFloat(m[0]),m[0].length]:[NaN,0]; };
  if(key[0]!=="@"){ const [x,n]=strtod(line); v.hv=x; if(n>0&&line[n]===" "){ let p=line.slice(n).replace(/^\s+/,""); const q=p.indexOf("("); if(q>=0){ const r=p.slice(q+1); const [hs,hn]=strtod(r); v.hscale=r[hn]===")"?hs:0; p=p.slice(0,q); } v.hvU=p; } return v; }
  const t=line[0]; v.type=t==="V"?"value":t==="S"?"select":t==="C"?"scale":null; if(!v.type||line[1]!==" ") return null;
  let s=line.slice(1).replace(/^\s+/,"");
  if(s[0]==="["){ const e=s.indexOf("]"); if(e<0) return null; if(e>1) v.soft=s.slice(1,e); s=s.slice(e+1); if(s[0]!==" ") return null; s=s.replace(/^\s+/,""); }
  if(s[0]==="("){ s=s.slice(1).replace(/^\s+/,""); let lv=1,i=0; for(;i<s.length&&lv;i++){ if(s[i]===")") lv--; else if(s[i]==="(") lv++; } if(lv) return null; const inner=s.slice(0,i-1); const [hs,hn]=strtod(inner); v.hscale=hs; let u=inner.slice(hn).trim(); if(u.endsWith("/LSB")) u=u.slice(0,-4); v.hscaleU=u||null; s=s.slice(i); if(s[0]!==" ") return null; s=s.replace(/^\s+/,""); }
  if(v.type==="select"){ v.hvStr=s.length>1&&s[0]==='"'&&s[s.length-1]==='"'?s.slice(1,-1):s; }
  else { const [x,n]=strtod(s); v.hv=x; v.hvStr=s; if(v.type==="value"&&n>0&&s[n]===" "&&s.indexOf(" ",n+1)<0){ v.hvU=s.slice(n).replace(/^\s+/,""); } }
  return v;
}
function readNanoscope(name,u8){
  const head=tdl.decode(u8.subarray(0,Math.min(u8.length,40)));
  const forceMagic=head.startsWith("\\*Force file list\r\n");   /* Gwyddion: 처음 종류일 뿐 — Operating mode 가 Force Volume·Image 면 바뀐다 */
  if(head.startsWith("?*File list\r\n")) throw new Error("텍스트형 Nanoscope 파일(?*File list)은 읽지 않습니다 — 이진(binary) 형식으로 내보내 주세요");
  if(!forceMagic&&!head.startsWith("\\*File list\r\n")&&!head.startsWith("\\*EC File list\r\n")){
    if(/^[\\?]\*File list/.test(head)) throw new Error("줄 끝이 바뀌어(CR LF → LF) 이진 자료가 망가진 Nanoscope 파일입니다 — 메일·FTP(ascii)·압축 과정에서 생깁니다. 원본을 다시 받아 주세요");
    throw new Error("Nanoscope 머리(\\*File list)가 아닙니다");
  }
  // 머리 길이: 앞쪽 몇 필드 안의 "\Data length: " (Gwyddion extract_header)
  let p=0, hlen=-1; const pre=tdl.decode(u8.subarray(0,Math.min(u8.length,4096)));
  for(let i=0;i<8;i++){ p=pre.indexOf("\\",p+1); if(p<0) break; if(pre.startsWith("\\Data length: ",p)){ const m=pre.slice(p+14).match(/^\d+/); hlen=m?+m[0]:-1; break; } }
  if(!(hlen>0)||hlen>u8.length) throw new Error("Nanoscope: 머리 길이(\\Data length)를 찾지 못했거나 파일보다 깁니다");
  const header=tdl.decode(u8.subarray(0,hlen)).replace(/\0[\s\S]*$/,"");
  // 절(\*…) 나누기 — 열쇠는 대소문자 무시
  const lists=[]; let cur=null;
  for(const ln of header.split("\n")){ const l=ln.replace(/\r$/,""); if(!l) continue;
    if(l.startsWith("\\*")){ const nm=l.slice(2); if(nm==="File list end") break; cur={self:nm,h:{},order:[]}; lists.push(cur); continue; }
    if(!cur||l[0]!=="\\") continue;
    const line=l.slice(1); let c0=0; if(line[0]==="@"&&/\d/.test(line[1])&&line[2]===":") c0=3; const ci=line.indexOf(":",c0);
    if(ci<0||!/\s/.test(line[ci+1]||"")) continue;
    const key=line.slice(0,ci).trim(), val=line.slice(ci+1).replace(/^\s+/,""); const v=nsValue(key,val); if(v){ cur.h[key.toLowerCase()]=v; cur.order.push(key); } }
  const get=(L,k)=>L?L.h[k.toLowerCase()]:undefined;
  let mode=null, version=0, startCtx=null;
  for(const L of lists){ const om=get(L,"Operating mode"); if(om) mode=om.hvStr; if(/^(File list|EC File list|Force file list)$/.test(L.self)){ const v=get(L,"Version"); if(v) version=parseInt(v.hvStr,16)||0; const sc=get(L,"Start context"); if(sc) startCtx=sc.hvStr; } }
  let kind=mode==="Force Volume"?"fvol":mode==="Force"?"force":mode==="Image"?"image":(forceMagic?"force":"image");
  if(startCtx==="FVOL"&&kind==="image") kind="fvol"; if(startCtx==="FOL") kind="force"; else if(startCtx&&/VAR$/.test(startCtx)&&kind==="image") kind="profiles";
  if(kind==="force") throw new Error("Nanoscope 힘 곡선 파일(Force)입니다 — 이미지가 아니어서 읽지 않습니다");
  if(kind==="profiles") throw new Error("Nanoscope 단면(profile) 묶음 파일(Deep Trench 등)입니다 — 이미지가 아니어서 읽지 않습니다");
  const bpp=version>=0x09200000?4:2;
  let scanner=null, scan=null, contr=null, gx=0, gy=0, gNonsq=false; const imgs=[], warn=[]; let nForceImg=0;
  const nonsq=L=>{ const v=get(L,"Aspect ratio"); if(!v||v.hvStr==="1:1") return false; const a=parseFloat(v.hvStr); return a>0&&a!==1; };
  const IMAGE_LISTS=/^(AFM image list|Ciao image list|STM image list|NCAFM image list|Ciao force image list|Image list)$/;
  for(const L of lists){
    if(L.self==="Scanner list"||L.self==="Microscope list"){ scanner=L; continue; }
    if(L.self==="Controller list"){ contr=L; continue; }
    if(/^(Ciao scan list|Afm list|Stm list|NC Afm list|Ciao force list)$/.test(L.self)){ const a=get(L,"Samps/line"), b=get(L,"Lines"); if(a) gx=Math.trunc(a.hv); if(b) gy=Math.trunc(b.hv); gNonsq=nonsq(L); if(L.self!=="Ciao force list") scan=L; continue; }
    if(!IMAGE_LISTS.test(L.self)) continue;
    if(L.self==="Ciao force image list"&&kind==="fvol"){ nForceImg++; continue; }
    imgs.push(nsImage(L));
  }
  void contr;
  function nsImage(L){
    { const sl=get(L,"Samps/line"); if(sl&&!get(L,"Number of lines")&&sl.hvU&&/^\d/.test(sl.hvU)){ const n=parseFloat(sl.hvU); L.h["number of lines"]={type:"old",hv:n,hvStr:sl.hvU}; } } /* Nanoscope E: "Samps/line: 512 512" (Gwyddion read_hash) */
    for(const k of ["Samps/line","Number of lines","Scan size","Data offset","Data length"]) if(!get(L,k)) throw new Error(`Nanoscope: '${k}' 가 없습니다 (${L.self})`);
    let W=Math.trunc(get(L,"Samps/line").hv), H=Math.trunc(get(L,"Number of lines").hv); const bq=get(L,"Bytes/pixel"); const qbpp=bq?Math.trunc(bq.hv):2;
    // Scan size: "7 7 ~m" 또는 옛 "5000 nm"
    const ss=get(L,"Scan size").hvStr; const m2=ss.match(/^\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?) \s*([+-]?[\d.]+(?:[eE][+-]?\d+)?) \s*(\S+)/), m1=ss.match(/^\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s+(\S+)/);
    let xr,yr,uxy; if(m2){ xr=+m2[1]; yr=+m2[2]; uxy=m2[3]; } else if(m1){ xr=yr=+m1[1]; uxy=m1[2]; } else throw new Error("Nanoscope: 'Scan size' 를 해석하지 못했습니다 ("+ss+")");
    uxy=uxy.slice(0,7); const uo=unitParse(uxy);   /* 접두어(~m → µm)는 calOf 가 곱한다 */
    const gxx=gx||W, gyy=gy||H; const off=Math.trunc(get(L,"Data offset").hv), size=Math.trunc(get(L,"Data length").hv);
    if(!(Number.isInteger(off)&&off>=0&&Number.isInteger(size)&&size>=0)) throw new Error(`Nanoscope: Data offset·Data length 가 수가 아닙니다 ('${get(L,"Data offset").hvStr}', '${get(L,"Data length").hvStr}')`);
    if(off>u8.length||size>u8.length-off) throw new Error("Nanoscope: 파일이 잘렸습니다 (자료가 파일 끝을 넘음)");
    let ok=false, glob=false; const ns=nonsq(L);
    if(size===bpp*W*H) ok=true; else if(size===bpp*gxx*gyy){ ok=true; glob=true; }
    else if(size>bpp*Math.max(W*H,gxx*gyy)){ ok=true; glob=W*H<gxx*gyy; } else if(size>bpp*Math.min(W*H,gxx*gyy)){ ok=true; glob=W*H>gxx*gyy; }
    if(!ok) throw new Error(`Nanoscope: 자료 크기가 맞지 않습니다 (기대 ${bpp*W*H} 바이트, 머리 ${size})`);
    if(glob){ if(gxx){ xr*=gxx/W; W=gxx; } if(gyy){ yr*=gyy/H; H=gyy; } }
    else if(ns){ if(gNonsq){ yr*=H; yr/=W; } else { yr*=H; yr/=gyy; } }
    if(!(W>0&&H>0)) throw new Error("Nanoscope: 이미지 크기가 0 입니다");
    if(W*H*4>IO.MAXBYTES) throw new Error("Nanoscope: 이미지가 너무 큽니다");
    const sc=nsScale(L); const k=sc.q/Math.pow(256,qbpp);
    const dv=new DataView(u8.buffer,u8.byteOffset+off,bpp*W*H); const raw=new Float64Array(W*H);
    for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const i=y*W+x; const r=bpp===4?dv.getInt32(4*i,true):dv.getInt16(2*i,true); raw[(H-1-y)*W+x]=r*k; } // 첫 줄이 아래 (Gwyddion: invert)
    const du=dispUnit({p10:0,u:sc.unit.u},raw); const pl=new Float32Array(W*H); for(let i=0;i<W*H;i++) pl[i]=raw[i]*du.f;
    const idn=get(L,"@2:Image Data")||get(L,"@3:Image Data")||get(L,"@4:Image Data"); const iname=idn?(idn.hvStr||idn.soft):(get(L,"Image data")||{}).hvStr;
    const dir=(get(L,"Line Direction")||{}).hvStr; const label=(iname||"Channel")+(dir?` (${dir})`:"");
    const cal=calOf(xr/W,yr/H,uo,"Nanoscope Scan Size / Samps·Lines");
    const info=[["Channel",iname||""],["Soft scale",idn&&idn.soft||""],["Line Direction",dir||""],["Version","0x"+version.toString(16)],["Raw data",`${bpp*8}-bit 정수 (눈금 Bytes/pixel=${qbpp})`],
      ["Value = raw ×",`${g6(k)} ${sc.base}/LSB  (${sc.how})`],["Value unit",du.f!==1?`${du.unit} (× ${g6(du.f)})`:du.unit||"(없음)"],["Scan Size",ss]];
    return {name:`${name} · ${label}`,w:W,h:H,bit:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[pl],cal,calCandidates:[],info,warn:[],valueUnit:du.unit,siScale:du.f,siUnit:du.base,metaText:header};
  }
  function nsScale(L){
    let v=get(L,"Z scale");
    if(v){ const u=unitParse(v.hvU); return {q:v.hv*Math.pow(10,u.p10),unit:u,base:unitStr(u),how:"Z scale (옛 형식)"}; }
    v=get(L,"@4:Z scale")||get(L,"@2:Z scale");
    if(!v) throw new Error("Nanoscope: Z 눈금(@2:Z scale)을 찾지 못했습니다 — 아주 오래된 형식은 읽지 않습니다");
    if(v.soft){ const sk="@"+v.soft; const sv=get(scanner,sk)||get(scan,sk);
      if(!sv){ warn.push(`soft scale '${v.soft}' 를 찾지 못해 hard value 만 곱했습니다 (단위 없음 — Gwyddion 과 같음)`); return {q:v.hv,unit:{p10:0,u:{}},base:"",how:`soft scale ${v.soft} 없음`}; }
      if(!sv.hvU){ const u=unitParse(v.soft==="Sens. Phase"?"deg":"V"); return {q:v.hv*sv.hv,unit:u,base:unitStr(u),how:`${g6(v.hv)} × ${v.soft} ${g6(sv.hv)}`}; }
      const u2=unitParse(sv.hvU), u1=v.hvU?unitParse(v.hvU):unitParse("V"); const u=unitMul(u2,u1); const q=v.hv*sv.hv*Math.pow(10,u.p10);
      return {q,unit:{p10:0,u:u.u},base:unitStr(u),how:`${g6(v.hv)} ${v.hvU||"V"} × ${v.soft} ${g6(sv.hv)} ${sv.hvU}`}; }
    const u=unitParse(v.hvU); return {q:v.hv*Math.pow(10,u.p10),unit:{p10:0,u:u.u},base:unitStr(u),how:`${g6(v.hv)} ${v.hvU||""} (soft scale 없음)`};
  }
  if(nForceImg) warn.push(`Force Volume: 힘 곡선 묶음 ${nForceImg} 개는 읽지 않았습니다 (이미지 채널만)`);
  if(!imgs.length) throw new Error("Nanoscope 파일에 이미지 채널이 없습니다");
  return {images:imgs,warn};
}

/* ============================================================ Asylum ARDF ============================================================ */
const CRC_T=(()=>{ const t=new Uint32Array(256); for(let n=0;n<256;n++){ let c=n; for(let k=0;k<8;k++) c=c&1?0xEDB88320^(c>>>1):c>>>1; t[n]=c>>>0; } return t; })();
function crc32(u8,a,b){ let c=0xFFFFFFFF; for(let i=a;i<b;i++) c=CRC_T[(c^u8[i])&255]^(c>>>8); return (c^0xFFFFFFFF)>>>0; }
function readARDF(name,u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength), N=u8.length;
  const u64=o=>dv.getUint32(o,true)+dv.getUint32(o+4,true)*4294967296;
  const cs=(o,n)=>{ let e=o; while(e<o+n&&u8[e]) e++; return tdw.decode(u8.subarray(o,e)); };
  function hdr(o,want,size){
    if(!(o>=0&&o+16<=N)) throw new Error(`ARDF: 블록 위치 ${o} 가 파일 밖입니다 (파일이 잘림)`);
    const h={o,crc:dv.getUint32(o,true),size:dv.getUint32(o+4,true),name:tdl.decode(u8.subarray(o+8,o+12)),flags:dv.getUint32(o+12,true)};
    if(want&&h.name!==want) throw new Error(`ARDF: ${want} 블록 자리(바이트 ${o})에 '${h.name.replace(/\0/g,"·")}' 가 있습니다 — 손상됐거나 모르는 변형`);
    if(size!=null&&h.size!==size) throw new Error(`ARDF: ${h.name} 블록 크기 ${h.size} (기대 ${size})`);
    if(h.size<16||o+h.size>N) throw new Error(`ARDF: ${h.name} 블록(바이트 ${o})이 파일 끝을 넘습니다 (파일이 잘림)`);
    if(crc32(u8,o+4,o+h.size)!==h.crc) throw new Error(`ARDF: ${h.name} 블록(바이트 ${o})의 CRC32 가 맞지 않습니다 — 파일 손상`);
    return h;
  }
  function toc(o,want,stride){ const h=hdr(o,want,32); const tsize=u64(o+16), n=dv.getUint32(o+24,true), st=dv.getUint32(o+28,true);
    if(st!==stride) throw new Error(`ARDF: ${want} 목차 항목 크기 ${st} (기대 ${stride})`); if(tsize!==32+n*st||o+tsize>N) throw new Error(`ARDF: ${want} 목차 크기가 맞지 않습니다`);
    const ents=[]; for(let e=o+32;e<o+tsize;e+=st){ if(!dv.getUint32(e,true)&&!dv.getUint32(e+4,true)) break; const eh=hdr(e,null,null); if(eh.size>st) throw new Error(`ARDF: ${want} 목차 항목(바이트 ${e}) 크기 ${eh.size} 가 항목 간격 ${st} 보다 큽니다`); ents.push({name:eh.name,p:st===24?u64(e+16):u64(e+24),idx:st===32?u64(e+16):0}); }
    return {h,tsize,ents,end:o+tsize}; }
  function text(p){ const h=hdr(p,"TEXT"); if(h.size<24) throw new Error("ARDF: TEXT 블록이 너무 짧습니다"); const len=dv.getUint32(p+20,true); if(24+len>h.size) throw new Error("ARDF: TEXT 길이가 블록보다 깁니다"); return tdw.decode(u8.subarray(p+24,p+24+len)).replace(/\r\n?/g,"\n"); }
  hdr(0,"ARDF",16);
  const ftoc=toc(16,"FTOC",24); const ttoc=toc(ftoc.end,"TTOC",32);
  const notes=ttoc.ents.filter(e=>e.name==="TOFF").map(e=>text(e.p)).join("\n");
  const nv={}; notes.split("\n").forEach(l=>{ const i=l.indexOf(":"); if(i>0) nv[l.slice(0,i).trim()]=l.slice(i+1).trim(); });
  const imgs=[], warn=[];
  for(const e of ftoc.ents){
    if(e.name==="IMAG"){
      const it=toc(e.p,"IMAG",24); const tt=toc(it.end,"TTOC",32); const inote=tt.ents.filter(x=>x.name==="TOFF").map(x=>text(x.p)).join("\n");
      const d=hdr(tt.end,"IDEF",184); const o=d.o+16; const W=dv.getUint32(o,true), H=dv.getUint32(o+4,true); const dx=dv.getFloat64(o+24,true), dy=dv.getFloat64(o+32,true);
      const xu=cs(o+40,32), yu=cs(o+72,32), title=cs(o+104,32), units=cs(o+136,32);
      const ib=hdr(d.o+184,"IBOX",32); const isz=u64(ib.o+16), nl=dv.getUint32(ib.o+24,true), st=dv.getUint32(ib.o+28,true);
      if(nl!==H||st!==16+4*W) throw new Error(`ARDF: '${title}' 의 IBOX(${nl} 줄, 줄 ${st} 바이트)가 IDEF(${W}×${H})와 맞지 않습니다`);
      if(ib.o+32+H*st>N||ib.o+isz>N) throw new Error(`ARDF: '${title}' 이미지 자료가 파일 끝을 넘습니다 (파일이 잘림)`);
      hdr(ib.o+isz,"GAMI",16);
      if(W*H*4>IO.MAXBYTES) throw new Error(`ARDF: '${title}' 이 너무 큽니다`);
      const raw=new Float64Array(W*H); let empty=0;
      for(let y=0;y<H;y++){ const lo=ib.o+32+y*st; const nm=tdl.decode(u8.subarray(lo+8,lo+12));
        if(nm==="\0\0\0\0"&&!dv.getUint32(lo,true)&&!dv.getUint32(lo+4,true)){ empty++; for(let x=0;x<W;x++) raw[(H-1-y)*W+x]=NaN; continue; } /* 비어 있는 줄(중단된 주사) — Gwyddion 도 자료 줄로 셈 */
        if(nm!=="IDAT") throw new Error(`ARDF: '${title}' ${y} 째 줄에 IDAT 가 없습니다 ('${nm.replace(/\0/g,"·")}')`);
        for(let x=0;x<W;x++) raw[(H-1-y)*W+x]=dv.getFloat32(lo+16+4*x,true); } // 첫 줄이 아래 (Igor wave 와 같은 방향)
      const uz=unitParse(units); const kz=Math.pow(10,uz.p10); if(kz!==1) for(let i=0;i<raw.length;i++) raw[i]*=kz;
      const du=dispUnit({p10:0,u:uz.u},raw); const p=new Float32Array(W*H); for(let i=0;i<W*H;i++) p[i]=raw[i]*du.f;
      const cal=calOf(dx,dy,unitParse(xu),"ARDF IDEF 간격 (x, y step)");
      const info=[["Channel",title],["IDEF units",`x ${xu}, y ${yu}, 값 ${units||"(없음)"}`],["Value unit",du.f!==1?`${du.unit} (파일 ${du.base} × ${g6(du.f)})`:du.unit||"(없음)"]];
      ["ImagingMode","FastScanSize","SlowScanSize","ScanRate","SpringConstant","InvOLS","FMapScanPoints","FMapScanLines","Date","Time"].forEach(k=>{ if(nv[k]!=null) info.push([k,nv[k]]); });
      imgs.push({name:`${name} · ${title}`,w:W,h:H,bit:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal,calCandidates:[],info,warn:empty?[`비어 있는 줄 ${empty} 개(중단된 주사?)는 NaN 으로 두었습니다`]:[],valueUnit:du.unit,siScale:du.f,siUnit:du.base,metaText:(inote?inote+"\n\n":"")+notes});
    } else if(e.name==="VOLM"){
      try{ const vt=toc(e.p,"VOLM",24); const tt=toc(vt.end,"TTOC",32); const vd=hdr(tt.end,"VDEF"); const P=dv.getUint32(vd.o+16,true), Lh=dv.getUint32(vd.o+20,true); const ch=[]; let q=vd.o+vd.size;
        for(let i=0;i<16;i++){ const h=hdr(q); if(h.name!=="VCHN") break; ch.push(cs(q+16,32)); q+=h.size; }
        warn.push(`힘 곡선 묶음(VOLM) ${P}×${Lh} 점${ch.length?" — 채널 "+ch.join("·"):""} 은 읽지 않았습니다 (이미지 채널만 읽음)`); }
      catch(err){ warn.push("힘 곡선 묶음(VOLM)이 있으나 구조를 해석하지 못했습니다 — 이미지 채널은 정상 ("+err.message+")"); }
    } else warn.push(`ARDF 목차 항목 '${e.name}' 은 읽지 않았습니다`);
  }
  if(!imgs.length) throw new Error("ARDF 안에 이미지 채널(IMAG)이 없습니다"+(warn.length?" — "+warn.join(" "):""));
  return {images:imgs,warn};
}

IO.readGWY=readGWY; IO.readGSF=readGSF; IO.readNanoscope=readNanoscope; IO.readARDF=readARDF;
IO._unitParse=unitParse; IO._unitStr=unitStr; IO._nsValue=nsValue; IO._crc32=crc32; IO._dispUnit=dispUnit;
})();
