/* ============================================================
   INNO Analytics — XPS 파일 읽기 (DOM 없음: 브라우저·Worker·node 공용)
   지원: Thermo Avantage .vgd (OLE2 원본, counts·dwell·scans·transmission 계수 포함)
         VAMAS .vms/.npl (ISO 14976, REGULAR/IRREGULAR, 여러 블록)
         Avantage 내보내기 .xls (BIFF8) / .xlsx, Avantage .avg (DataSpace 텍스트)
         2~3열 텍스트 .txt/.csv/.dat/.xy (BE 또는 KE, CPS)
   모든 결과는 같은 형태의 스펙트럼 객체로 돌려준다(아래 mkSpec).
   ============================================================ */
(function(root){
"use strict";

/* ---------------- 공통 ---------------- */
const TD8=(typeof TextDecoder!=="undefined")?new TextDecoder("latin1"):null;
const TD16=(typeof TextDecoder!=="undefined")?new TextDecoder("utf-16le"):null;
const TDU=(typeof TextDecoder!=="undefined")?new TextDecoder("utf-8"):null;
function latin1(u8){ if(TD8) return TD8.decode(u8); let s=""; for(let i=0;i<u8.length;i++) s+=String.fromCharCode(u8[i]); return s; }
function utf16(u8){ if(TD16) return TD16.decode(u8); let s=""; for(let i=0;i+1<u8.length;i+=2) s+=String.fromCharCode(u8[i]|(u8[i+1]<<8)); return s; }
function utf8(u8){ return TDU?TDU.decode(u8):latin1(u8); }
function dv(u8){ return new DataView(u8.buffer,u8.byteOffset,u8.byteLength); }
const NULLV=1e37;                          // VAMAS 의 '값 없음'
const isNullV=v=>!isFinite(v)||Math.abs(v)>=0.99e37;

/* 영역 이름 → 원소·궤도 ("Co2p", "Ta 4f7/2", "O1s Scan", "fe2p", "Survey")
   원소 기호는 주기율표(Z = 1–92 + 몇 개)에 있는 것만 인정한다 — "Peak 2p" 같은 이름을 원소로 오인하지 않게. */
const SYMS="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm".split(" ");
const SYMSET=new Set(SYMS), SYMLC=new Map(SYMS.map(e=>[e.toLowerCase(),e]));
function parseRegionName(s){
  const t=String(s||"").replace(/\.(txt|vgd|vms|npl|avg|csv|dat|xy)$/i,"").trim();
  if(/survey|wide|overview|와이드|서베이/i.test(t)) return {kind:"survey",label:"Survey"};
  if(/valence|\bVB\b|가전자/i.test(t)) return {kind:"vb",label:"VB"};
  // 1) 대소문자 그대로 (Co2p, Ta 4f7/2)
  const re=/(?:^|[^A-Za-z])([A-Z][a-z]?)[\s_\-]*([1-7])\s*([spdf])\s*(1\/2|3\/2|5\/2|7\/2)?(?![a-z])/g; let m;
  while((m=re.exec(" "+t))){ if(SYMSET.has(m[1])) return {kind:"region",el:m[1],orb:m[2]+m[3],j:m[4]||null,label:m[1]+" "+m[2]+m[3]+(m[4]||"")}; }
  // 2) 대소문자 무시 (fe2p, CO2P) — 기호표에 있는 것만
  const re2=/(?:^|[^A-Za-z])([A-Za-z]{1,2})[\s_\-]*([1-7])\s*([spdfSPDF])\s*(1\/2|3\/2|5\/2|7\/2)?(?![a-z])/g;
  while((m=re2.exec(" "+t))){ const e=SYMLC.get(m[1].toLowerCase()); if(e) return {kind:"region",el:e,orb:m[2]+m[3].toLowerCase(),j:m[4]||null,label:e+" "+m[2]+m[3].toLowerCase()+(m[4]||""),caseFixed:true}; }
  const a=/(?:^|[^A-Za-z])([A-Z][a-z]?)\s*(KLL|LMM|MNN|KVV|LVV|MVV|NOO)/.exec(" "+t);
  if(a&&SYMSET.has(a[1])) return {kind:"auger",el:a[1],label:a[1]+" "+a[2]};
  return {kind:"region",label:t||"?",unknown:true};
}
/* 파일 이름에서 시료 이름: "[D.I. on BCT (600'C 5h)] Ta4f.txt" → "D.I. on BCT (600'C 5h)" */
function sampleFromName(fn){
  const b=String(fn||"").replace(/^.*[\\/]/,"");
  const m=/^\[(.+?)\]/.exec(b); if(m) return m[1].trim();
  const r=parseRegionName(b);
  let s=b.replace(/(\.(txt|vms|npl|csv|dat|xy|vgd|avg|xlsx?))+$/i,"");
  s=s.replace(/(Survey|Scan|Snap|Region)/ig,"");
  if(r.el) s=s.replace(new RegExp("[_\\s\\-]*"+r.el+"\\s*"+(r.orb||"")+"(1\\/2|3\\/2|5\\/2|7\\/2)?$"),"");
  s=s.replace(/[_\s\-]+$/,"").replace(/^[_\s\-]+/,"");
  return s||"";
}

/* Avantage 경로 "…\\Experiment2\\X-Ray018 400um - FG  ON\\W1\\Group\\C1s Scan.VGD" → "W1" (Group 위 폴더 = 측정 지점/시료) */
/* Avantage 원래 경로 → 측정 지점(점) 폴더 이름. 예: …\\W1\\Group\\C1s Scan.VGD → "W1" */
function pathPoint(p){
  const a=String(p||"").split(/[\\/]/).filter(Boolean); if(a.length<2) return "";
  a.pop(); if(/^group$/i.test(a[a.length-1])) a.pop();
  return a.length?a[a.length-1].replace(/\s+/g," "):"";
}
/* VGD 의 시료 키: 측정 지점 폴더 + 위로 세 단계(자료·실험·측정 조건) — 다른 실험의 같은 지점 이름(W1 등)이 한 시료로 합쳐지지 않게.
   예: C:\\DATA\\K231213C\\Experiment2\\X-Ray018 400um - FG  ON\\W1\\Group\\C1s Scan.VGD → "K231213C / Experiment2 / X-Ray018 400um - FG ON / W1"
   (표시: 이름 = 마지막 부분, 맥락 = 앞부분) */
function pathSample(p){
  const a=String(p||"").split(/[\\/]/).filter(Boolean); if(a.length<2) return "";
  a.pop(); if(/^group$/i.test(a[a.length-1])) a.pop();
  let seg=a.filter(s=>!/^[a-z]:$/i.test(s)).map(s=>s.replace(/\s+/g," ").trim()).filter(Boolean);
  seg=seg.slice(-4); if(seg.length>1&&/^data$/i.test(seg[0])) seg.shift();
  return seg.join(" / ");
}
/* 스펙트럼 객체: be 는 오름차순 Float64Array, y 는 CPS. c = 1 CPS 당 count (모르면 null) */
let _uid=0;
function mkSpec(o){
  const n=o.be.length; let be=Float64Array.from(o.be), y=Float64Array.from(o.y), yr=o.yRaw?Float64Array.from(o.yRaw):null, tf=o.tfArr?Float64Array.from(o.tfArr):null;
  if(n>1&&be[0]>be[n-1]){ be.reverse(); y.reverse(); if(yr) yr.reverse(); if(tf) tf.reverse(); }
  for(let i=1;i<n;i++) if(!(be[i]>be[i-1])) throw new Error("에너지 축이 단조롭지 않습니다 ("+o.name+")");
  const step=n>1?(be[n-1]-be[0])/(n-1):0;
  const reg=parseRegionName(o.region||o.name);
  const region=reg.el||reg.kind==="survey"?reg.label:(o.region||reg.label);
  return {uid:"x"+(++_uid)+"_"+Date.now().toString(36), name:o.name, region, reg, sample:o.sample||"", file:o.file||"", format:o.format,
    hv:o.hv||1486.6, hvSrc:o.hvSrc||"기본값(Al Kα 1486.6)", be, y, yRaw:yr, c:o.c||null, cSrc:o.cSrc||null, step, n,
    tfArr:tf, meta:o.meta||{}, warn:o.warn||[], extra:o.extra||null};
}

/* ============================================================
   1) Microsoft Compound File (OLE2/CFB) — VGD·XLS 공용
   ============================================================ */
function cfbIs(u8){ return u8.length>=512&&u8[0]===0xD0&&u8[1]===0xCF&&u8[2]===0x11&&u8[3]===0xE0&&u8[4]===0xA1&&u8[5]===0xB1&&u8[6]===0x1A&&u8[7]===0xE1; }
function cfbOpen(u8){
  if(!cfbIs(u8)) throw new Error("OLE2(Compound File) 형식이 아닙니다");
  const D=dv(u8);
  const ss=1<<D.getUint16(30,true), mss=1<<D.getUint16(32,true);
  const nFat=D.getUint32(44,true), dir0=D.getUint32(48,true), cutoff=D.getUint32(56,true)||4096;
  const mf0=D.getUint32(60,true), nMf=D.getUint32(64,true), df0=D.getUint32(68,true), nDf=D.getUint32(72,true);
  const off=s=>(s+1)*ss;
  const ENDC=0xFFFFFFFE, FREE=0xFFFFFFFF;
  // DIFAT → FAT 섹터 목록
  const fatSecs=[];
  for(let i=0;i<109&&fatSecs.length<nFat;i++){ const s=D.getUint32(76+4*i,true); if(s!==FREE) fatSecs.push(s); }
  let ds=df0, guard=0;
  while(fatSecs.length<nFat&&ds!==ENDC&&ds!==FREE&&guard++<nDf+2){
    const o=off(ds); for(let i=0;i<ss/4-1&&fatSecs.length<nFat;i++){ const s=D.getUint32(o+4*i,true); if(s!==FREE) fatSecs.push(s); }
    ds=D.getUint32(o+ss-4,true);
  }
  const per=ss/4, fat=new Uint32Array(fatSecs.length*per);
  fatSecs.forEach((s,k)=>{ const o=off(s); for(let i=0;i<per;i++) fat[k*per+i]=o+4*i+4<=u8.length?D.getUint32(o+4*i,true):FREE; });
  const chain=(s,maxN)=>{ const out=[]; let g=0; while(s!==ENDC&&s!==FREE&&s<fat.length&&g++<(maxN||1e7)){ out.push(s); s=fat[s]; } return out; };
  const readChain=(s,size)=>{ const cs=chain(s); const buf=new Uint8Array(cs.length*ss); cs.forEach((c,i)=>{ const o=off(c); buf.set(u8.subarray(o,Math.min(o+ss,u8.length)),i*ss); }); return size!=null?buf.subarray(0,size):buf; };
  // 디렉터리
  const dbuf=readChain(dir0); const DD=dv(dbuf); const ents=[];
  for(let o=0;o+128<=dbuf.length;o+=128){
    const nl=DD.getUint16(o+64,true); const type=dbuf[o+66];
    const name=nl>=2?utf16(dbuf.subarray(o,o+nl-2)):"";
    ents.push({name,type,left:DD.getUint32(o+68,true),right:DD.getUint32(o+72,true),child:DD.getUint32(o+76,true),start:DD.getUint32(o+116,true),size:DD.getUint32(o+120,true)+(ss===4096?DD.getUint32(o+124,true)*4294967296:0)});
  }
  const rootE=ents[0];
  let mini=null, mfat=null;
  const getMini=()=>{ if(mini) return; mini=readChain(rootE.start,rootE.size); const mb=readChain(mf0); const MD=dv(mb); mfat=new Uint32Array(Math.min(nMf*per,mb.length/4)); for(let i=0;i<mfat.length;i++) mfat[i]=MD.getUint32(4*i,true); };
  const readMini=(s,size)=>{ getMini(); const out=new Uint8Array(size); let p=0,g=0; while(s!==ENDC&&s!==FREE&&p<size&&g++<1e7){ const o=s*mss; const n=Math.min(mss,size-p); out.set(mini.subarray(o,o+n),p); p+=n; s=mfat[s]; } return out; };
  // 경로 목록 (트리 순회)
  const paths={};
  const walk=(id,pre,depth)=>{ if(id===FREE||id>=ents.length||depth>64) return; const e=ents[id]; walk(e.left,pre,depth+1); const p=pre?pre+"/"+e.name:e.name; paths[p]=id; if(e.type===1||e.type===5) walk(e.child,p,depth+1); walk(e.right,pre,depth+1); };
  walk(rootE.child,"",0);
  return {
    list:()=>Object.keys(paths),
    has:p=>p in paths,
    read:p=>{ const id=paths[p]; if(id==null) throw new Error("스트림 없음: "+p); const e=ents[id]; if(e.type!==2) throw new Error("스트림이 아님: "+p); return e.size<cutoff?readMini(e.start,e.size):readChain(e.start,e.size); },
  };
}

/* ---------------- OLE property set (MS-OLEPS) ---------------- */
function propsetParse(u8){
  const D=dv(u8); if(u8.length<48) return [];
  const nsec=D.getUint32(24,true), out=[];
  for(let s=0;s<nsec&&28+20*s+20<=u8.length;s++){
    const fmt=[...u8.subarray(28+20*s,44+20*s)].map(b=>b.toString(16).padStart(2,"0")).join("");
    const so=D.getUint32(44+20*s,true); if(so+8>u8.length) continue;
    const np=D.getUint32(so+4,true), props={};
    for(let i=0;i<np&&so+8+8*i+8<=u8.length;i++){
      const pid=D.getUint32(so+8+8*i,true), po=D.getUint32(so+12+8*i,true);
      if(pid===0||pid===1) continue;
      try{ props[pid]=propVal(u8,D,so+po); }catch(e){ props[pid]={t:-1,v:null,err:String(e.message||e)}; }
    }
    out.push({fmt,props});
  }
  return out;
}
function propVal(u8,D,o){
  const t=D.getUint16(o,true); o+=4;
  const one=(t,o,vec)=>{ switch(t){
    case 0: case 1: return [null,o];
    case 2: return [D.getInt16(o,true),o+(vec?2:4)];
    case 0x12: return [D.getUint16(o,true),o+(vec?2:4)];
    case 3: case 0x16: return [D.getInt32(o,true),o+4];
    case 0x13: case 0x17: return [D.getUint32(o,true),o+4];
    case 4: return [D.getFloat32(o,true),o+4];
    case 5: return [D.getFloat64(o,true),o+8];
    case 0x14: return [Number(D.getBigInt64(o,true)),o+8];
    case 0xB: return [D.getUint16(o,true)!==0,o+(vec?2:4)];
    case 0x40: { const lo=D.getUint32(o,true), hi=D.getUint32(o+4,true); const ms=(hi*4294967296+lo)/1e4-11644473600000; return [new Date(ms).toISOString(),o+8]; }
    case 0x1F: { const n=D.getUint32(o,true); const s=utf16(u8.subarray(o+4,o+4+2*n)).replace(/\0+$/,""); return [s,(o+4+2*n+3)&~3]; }
    case 0x1E: case 0x41: { const n=D.getUint32(o,true); const b=u8.subarray(o+4,o+4+n); return [t===0x1E?latin1(b).replace(/\0+$/,""):"(blob "+n+" B)",(o+4+n+3)&~3]; }
    case 0x48: return [[...u8.subarray(o,o+16)].map(b=>b.toString(16).padStart(2,"0")).join(""),o+16];
    case 8: { const n=D.getUint32(o,true); const s=utf16(u8.subarray(o+4,o+4+n)).replace(/\0+$/,""); return [s,(o+4+n+3)&~3]; } // BSTR (바이트 수)
    default: throw new Error("형식 0x"+t.toString(16));
  } };
  if(t&0x1000){ const bt=t&0xFFF; const n=D.getUint32(o,true); o+=4; const arr=[]; for(let i=0;i<n;i++){ if(bt===0xC){ const tt=D.getUint16(o,true); let v; [v,o]=one(tt,o+4,false); arr.push(v); } else { let v; [v,o]=one(bt,o,true); arr.push(v); } } return {t,v:arr}; }
  const [v]=one(t,o,false); return {t,v};
}

/* ============================================================
   2) Thermo Avantage .VGD
   property id 의미(같은 장비 5개 파일 대조·Avantage .avg 의 이름표와 맞춰 확인):
     105 광원 에너지 · 108/109 spot(µm) · 305 dwell(s) · 306 scan 수 · 404 pass energy · 406 분석기 work function
     411 lens mode · 605 flood gun 모드 · 303/304 측정 시작·끝 · 6 장비 · 17 Avantage 판 · 16 원래 경로
     100000–100003 transmission 함수 계수(DS_ANPROPID_TXFN_COEFF[0..3])
     415 Acq_TXFN_Correction (bool: 측정 때 transmission 을 이미 보정·제거) · 416 TXFN 보정 모드 (3 = 'TXF Shape Correction, retaining Intensity' 에서 봄)
       → meta.tfRemoved = (415 === true), meta.tfMode = 416. tfRemoved 이면 정량에서 TF 를 다시 나누면 안 된다 (이중 보정, test/ref/avantage_a2_patch.md §5)
   VGData = scan 합 counts(double). CPS = counts/(dwell·scans).
   ============================================================ */
function vgdParse(u8,fileName){
  const C=cfbOpen(u8);
  if(!C.has("VGData")) throw new Error("VGData 스트림이 없습니다 (VGD 가 아님)");
  const raw=C.read("VGData"); const RD=dv(raw); const total=raw.length>>3;
  const vals=new Float64Array(total); for(let i=0;i<total;i++) vals[i]=RD.getFloat64(8*i,true);
  // 에너지 축: VGSpaceAxes 안의 (u32 개수, f64 시작, f64 간격) — 개수가 데이터 길이와 맞는 자리를 찾는다
  let ax=null;
  if(C.has("VGSpaceAxes")){
    const sa=C.read("VGSpaceAxes"), SD=dv(sa);
    for(let o=0;o+20<=sa.length;o++){
      const n=SD.getUint32(o,true); if(n<2||n>total||total%n) continue;
      const st=SD.getFloat64(o+4,true), w=SD.getFloat64(o+12,true);
      if(isFinite(st)&&isFinite(w)&&Math.abs(w)>1e-6&&Math.abs(w)<20&&st>-50&&st<20000){ ax={n,start:st,step:w,off:o}; break; }
    }
  }
  if(!ax) throw new Error("VGD 에너지 축(VGSpaceAxes)을 해석하지 못했습니다");
  const nsp=total/ax.n;
  // 속성
  let P={}; const streams=C.list().filter(p=>/^\u0005/.test(p.split("/").pop())&&!/SummaryInformation/.test(p));
  let best=null; streams.forEach(p=>{ try{ const ps=propsetParse(C.read(p)); ps.forEach(sec=>{ const k=Object.keys(sec.props).length; if(!best||k>best.k) best={k,props:sec.props,stream:p}; }); }catch(e){} });
  if(best) P=best.props;
  const gv=id=>P[id]&&P[id].v!=null?P[id].v:null;
  const r6=v=>typeof v==="number"?+v.toPrecision(6):v;   // float32 → 사람이 쓴 값(0.05000000074 → 0.05)
  const hv=r6(gv(105)), dwell=r6(gv(305)), scans=gv(306)!=null?gv(306):gv(328), pe=r6(gv(404)), wf=r6(gv(406));
  const tfc=[100000,100001,100002,100003].map(gv); const tfOK=tfc.every(v=>typeof v==="number"&&isFinite(v));
  const txfnCorr=gv(415), txfnMode=gv(416), tfRemoved=txfnCorr===true;   // Acq_TXFN_Correction (bool) · 모드 (3 = TXF Shape Correction, retaining Intensity)
  let title="", summ=null; try{ if(C.has("\u0005SummaryInformation")){ summ=propsetParse(C.read("\u0005SummaryInformation")); const pr=summ[0]&&summ[0].props; if(pr&&pr[2]) title=String(pr[2].v||""); } }catch(e){}
  const warn=[];
  const hvUse=(typeof hv==="number"&&hv>100)?hv:1486.68;
  if(!(typeof hv==="number"&&hv>100)) warn.push("광원 에너지를 파일에서 찾지 못해 1486.68 eV 로 둠");
  const c=(dwell>0&&scans>0)?+(dwell*scans).toPrecision(10):null;
  if(!c) warn.push("dwell·scan 수를 찾지 못함 → counts 환산 불가(잡음 척도는 데이터에서 추정)");
  { let fr=0; for(let i=0;i<total;i++) fr+=Math.abs(vals[i]-Math.round(vals[i])); fr/=Math.max(1,total);
    if(c&&fr>0.05) warn.push("VGData 값이 정수가 아님(검출기 채널→에너지 보간으로 추정) — 순수 Poisson counts 가 아니므로 잡음의 절대 척도는 피팅의 χ²ν(φ)로 다시 정한다. counts 환산(dwell·scans)은 분산의 모양(∝ 세기)에만 쓴다"); }
  const ke=Float64Array.from({length:ax.n},(_,i)=>ax.start+i*ax.step);
  // 축이 KE 인지 BE 인지: Avantage 는 KE 로 저장한다. 시작값이 광원 에너지보다 크면 해석 불가.
  if(ke[0]>hvUse+5) throw new Error("VGD 에너지 축이 광원 에너지보다 큽니다 — 축 해석 실패");
  const be=ke.map(k=>hvUse-k);
  const base=(title||String(fileName||"").replace(/^.*[\\/]/,"")).replace(/\.vgd$/i,"").replace(/\s*(Scan|Snap|Region)$/i,"").trim();
  const meta={instrument:gv(6)||"", software:gv(17)||"", origPath:gv(16)||"", passEnergy:pe, dwell, scans, workFunction:wf, lensMode:gv(411)||"", floodGun:gv(605)||"",
    spot:[gv(108),gv(109)], start:gv(303), end:gv(304), tf:tfOK?{model:"thermo",coef:tfc,pe,src:"VGD 파일의 transmission 계수 (형태 ln T = Σaₖ(ln(KE/PE))ᵏ 은 추정)"+(tfRemoved?" (removed at acquisition, prop 415)":"")}:null,
    tfRemoved, tfMode:txfnMode,
    propIds:Object.keys(P).length, propStream:best?best.stream.replace(/\u0005/g,"\\x05"):""};
  const sample=pathSample(meta.origPath);
  const out=[];
  for(let s=0;s<nsp;s++){
    const cnt=vals.subarray(s*ax.n,(s+1)*ax.n);
    const y=c?Array.from(cnt,v=>v/c):Array.from(cnt);
    out.push(mkSpec({name:base+(nsp>1?" #"+(s+1):""),region:base,sample,file:fileName,format:"VGD",hv:hvUse,hvSrc:"VGD 파일",be,y,yRaw:cnt,c,cSrc:c?`VGD: dwell ${+dwell.toPrecision(4)} s × ${scans} scans`:null,meta,warn:warn.slice()}));
  }
  return out;
}

/* ============================================================
   3) VAMAS (ISO 14976)
   ============================================================ */
function vamasParse(text,fileName){
  const L=text.split(/\r?\n/); let p=0;
  const nx=()=>{ if(p>=L.length) throw new Error("VAMAS 파일이 중간에 끝났습니다 (줄 "+p+")"); return L[p++]; };
  const nxt=()=>nx().trim(), nxi=()=>parseInt(nxt(),10), nxf=()=>parseFloat(nxt());
  const H={};
  H.format=nxt(); if(!/VAMAS/i.test(H.format)) throw new Error("VAMAS 머리줄이 아닙니다");
  H.inst=nxt(); H.model=nxt(); H.operator=nxt(); H.exp=nxt();
  const nc=nxi(); H.comment=[]; for(let i=0;i<nc;i++) H.comment.push(nx());
  H.mode=nxt(); H.scan=nxt();
  const MODES1=["MAP","MAPDP","NORM","SDP"];
  if(MODES1.includes(H.mode)) H.nReg=nxi();
  if(H.mode==="MAP"||H.mode==="MAPDP"){ H.nPos=nxi(); H.nX=nxi(); H.nY=nxi(); }
  H.nExp=nxi(); H.expVars=[]; for(let i=0;i<H.nExp;i++) H.expVars.push({label:nxt(),unit:nxt()});
  const nInc=nxi(); if(nInc!==0) throw new Error("VAMAS 매개변수 포함/제외 목록(값 "+nInc+")은 아직 지원하지 않습니다");
  const nMan=nxi(); for(let i=0;i<nMan;i++) nxt();
  const nFutE=nxi(), nFutB=nxi(); for(let i=0;i<nFutE;i++) nxt();
  const nb=nxi(); const out=[];
  const SPUT_MODES=["MAPDP","MAPSVDP","SDP","SDPSV"], SPUT_TECH=["FABMS","FABMS energy spec","ISS","SIMS","SIMS energy spec","SNMS","SNMS energy spec"];
  for(let b=0;b<nb;b++){
    const B={}; B.id=nxt(); B.sample=nxt();
    B.date=[nxi(),nxi(),nxi(),nxi(),nxi(),parseInt(nxt(),10)]; B.gmt=nxt();
    const ncl=nxi(); B.comment=[]; for(let i=0;i<ncl;i++) B.comment.push(nx());
    B.tech=nxt();
    if(H.mode==="MAP"||H.mode==="MAPDP"){ B.x=nxt(); B.y=nxt(); }
    B.expVals=[]; for(let i=0;i<H.nExp;i++) B.expVals.push(nxt());
    B.src=nxt();
    if(SPUT_MODES.includes(H.mode)||SPUT_TECH.includes(B.tech)){ B.sputIon=[nxt(),nxt(),nxt()]; }
    B.hv=nxf(); B.power=nxt(); B.bwx=nxt(); B.bwy=nxt();
    if(["MAP","MAPDP","MAPSV","MAPSVDP","SEM"].includes(H.mode)){ B.fovX=nxf(); B.fovY=nxf(); }
    if(["MAPSV","MAPSVDP","SEM"].includes(H.mode)){ B.linescan=[]; for(let i=0;i<6;i++) B.linescan.push(nxf()); }   // ISO 14976: first/last linescan start·finish x,y (MAPSV·MAPSVDP·SEM; v1.7 fix — was MAP)
    B.polar=nxf(); B.azim=nxf(); B.anMode=nxt(); B.pe=nxf();
    if(B.tech==="AES diff") nxt();
    B.mag=nxf(); B.wf=nxf(); B.bias=nxf(); B.awx=nxf(); B.awy=nxf(); B.toPolar=nxf(); B.toAz=nxf();
    B.species=nxt(); B.trans=nxt(); B.charge=nxt();
    let nvar;
    if(H.scan==="REGULAR"){ B.absLabel=nxt(); B.absUnit=nxt(); B.absStart=nxf(); B.absStep=nxf(); nvar=nxi(); }
    else { nvar=nxi()-1; B.absLabel=nxt(); B.absUnit=nxt(); }
    B.vars=[]; for(let i=0;i<nvar;i++) B.vars.push({label:nxt(),unit:nxt()});
    B.sigMode=nxt(); B.dwell=nxf(); B.scans=nxi(); B.tcorr=nxf();
    if(SPUT_MODES.includes(H.mode)&&["AES","AES diff","AES dir","EDX","ELS","UPS","XPS","XRF"].includes(B.tech)){ B.sput={energy:nxf(),current:nxf(),widthX:nxf(),widthY:nxf(),polar:nxf(),azimuth:nxf(),mode:nxt()}; }
    B.tilt=nxf(); B.tiltAz=nxf(); B.rot=nxf();
    const nadd=nxi(); B.add=[]; for(let i=0;i<nadd;i++) B.add.push({label:nxt(),unit:nxt(),value:nxt()});
    for(let i=0;i<nFutB;i++) nxt();
    const nord=nxi();
    const nv=H.scan==="REGULAR"?nvar:nvar+1;
    for(let i=0;i<nv;i++){ nxt(); nxt(); }      // 각 변수 최소·최대
    const raw=new Float64Array(nord); for(let i=0;i<nord;i++) raw[i]=parseFloat(nxt());
    const np=Math.floor(nord/nv);
    let x, cols=[];
    if(H.scan==="REGULAR"){ x=Float64Array.from({length:np},(_,i)=>B.absStart+i*B.absStep); for(let v=0;v<nv;v++) cols.push(Float64Array.from({length:np},(_,i)=>raw[i*nv+v])); }
    else { x=Float64Array.from({length:np},(_,i)=>raw[i*nv]); for(let v=1;v<nv;v++) cols.push(Float64Array.from({length:np},(_,i)=>raw[i*nv+v])); }
    if(B.tech!=="XPS"&&B.tech!=="UPS"&&B.tech!=="AES"&&B.tech!=="AES dir") { /* 그대로 읽되 표시 */ }
    const warn=[];
    const hvOK=!isNullV(B.hv)&&B.hv>0;
    const hv=hvOK?B.hv:1486.6; if(!hvOK) warn.push("광원 에너지가 비어 있어 1486.6 eV 로 둠");
    const isKE=/kinetic/i.test(B.absLabel), isBE=/binding/i.test(B.absLabel);
    if(!isKE&&!isBE) warn.push("가로축 이름 '"+B.absLabel+"' — 운동에너지로 간주");
    const be=isBE?Array.from(x):Array.from(x,k=>hv-k);
    // 세기 열: 이름에 Transmission 이 없는 첫 열
    let yi=B.vars.findIndex(v=>!/transmission/i.test(v.label)); if(yi<0) yi=0;
    const ti=B.vars.findIndex(v=>/transmission/i.test(v.label));
    const yr=cols[yi];
    const c=(!isNullV(B.dwell)&&B.dwell>0&&B.scans>0)?B.dwell*B.scans:null;
    const y=Array.from(yr,v=>c?v/c:v);
    const regRaw=(B.species&&B.species!=="Survey"&&B.species!=="survey"?B.species+" ":"")+(B.trans&&B.trans.trim()?B.trans.trim():"");
    let region=regRaw.trim()||B.id;
    if(/^survey$/i.test(B.species)&&!/survey/i.test(B.id)) region=B.id.replace(/\.(txt|vgd)$/i,"").replace(/^\[.*?\]\s*/,"")||"Survey";
    if(/^survey$/i.test(B.species)&&/survey/i.test(B.id)) region="Survey";
    const meta={instrument:[H.inst,H.model].filter(s=>s&&s!=="Not Specified").join(" "),sourceLabel:B.src,passEnergy:isNullV(B.pe)?null:B.pe,anMode:B.anMode,
      workFunction:isNullV(B.wf)?null:B.wf,dwell:isNullV(B.dwell)?null:B.dwell,scans:B.scans,sigMode:B.sigMode,ordUnit:B.vars[yi]&&B.vars[yi].unit,
      date:B.date.slice(0,3).some(v=>v>0)?B.date:null,blockId:B.id,sampleId:B.sample,casa:H.comment.some(s=>/casa/i.test(s)),technique:B.tech,
      tf:null,vamas:vamasBlockMeta(H,B,b),casaModel:casaCommentParse(B.comment)};
    meta.vamas.tag=seriesTag(H,B);   // 표시용 "Etch Time 30 s" (이름은 바꾸지 않음 — CasaXPS 는 모든 파일에 'Data Set' 변수를 쓴다)
    if(meta.casa&&c===1) warn.push("CasaXPS 가 만든 VAMAS(dwell·scan = 1) — 세기가 이미 CPS 일 수 있어 counts 환산은 하지 않음");
    out.push(mkSpec({name:region,region,file:fileName,format:"VAMAS",hv,hvSrc:hvOK?"VAMAS 블록":"기본값",be,y,yRaw:c&&c!==1?yr:null,c:c&&!(meta.casa&&c===1)?c:null,
      cSrc:c&&!(meta.casa&&c===1)?`VAMAS: collection time ${B.dwell} s × ${B.scans} scans`:null,tfArr:ti>=0?cols[ti]:null,
      sample:(B.id.match(/^\[(.+?)\]/)||[])[1]||(B.sample&&B.sample!==B.id?B.sample:""),meta,warn}));
  }
  return out;
}

/* CasaXPS 가 VAMAS 블록 주석에 남기는 모형 (v1.8, 실측 CasaXPS 2.3.25 파일에서 읽은 문법 — CasaXPS 코드·문서 원문은 쓰지 않음):
   'Calib M = 284.52 A = 284.8 BE ADD'  → 표시 BE 이동 (A − M; 여러 줄이면 합)
   'CASA region (*name*) (*BG*) KEstart KEend RSF AvWidth StartOffset EndOffset c1 c2 c3 c4 (*tag*) mass …'  (c1…c4 = Cross Section 칸)
   'CASA comp (*name*) (*shape*) Area v lo hi link k MFWHM v lo hi link k Position v(KE) lo hi link k RSF … MASS … INDEX …'
   link = 연결 대상 성분 번호(0 부터, −1 = 연결 없음), k = 배율 (CasaXPS 화면 'A*0.667' 등). 에너지는 KE(파일 축 그대로). */
function casaCommentParse(lines){
  const L=(lines||[]).map(s=>String(s)); const out={calib:[],shift:0,regions:[],comps:[]};
  const strs=s=>{ const a=[]; s.replace(/\(\*(.*?)\*\)/g,(m,t)=>{ a.push(t); return ""; }); return a; };
  const rest=s=>s.replace(/\(\*.*?\*\)/g," ").trim().split(/\s+/).filter(Boolean);
  L.forEach(s=>{ let m;
    if((m=s.match(/^\s*Calib\s+M\s*=\s*([-\d.eE+]+)\s+A\s*=\s*([-\d.eE+]+)\s+(BE|KE)\s+ADD/i))){ const M=+m[1], A=+m[2]; out.calib.push({M,A,axis:m[3].toUpperCase()}); out.shift+=(m[3].toUpperCase()==="BE"?1:-1)*(A-M); return; }
    if(/^\s*CASA region\b/i.test(s)){ const t=strs(s), n=rest(s.replace(/^\s*CASA region/i,"")).map(Number);
      out.regions.push({name:t[0]||"",bg:t[1]||"",keStart:n[0],keEnd:n[1],rsf:n[2],aw:n[3],startOffset:n[4],endOffset:n[5],cross:n.slice(6,10),tag:t[2]||"",mass:n[10]}); return; }
    if(/^\s*CASA comp\b/i.test(s)){ const t=strs(s), body=s.replace(/^\s*CASA comp\s*\(\*.*?\*\)\s*\(\*.*?\*\)/i,""); const c={name:t[0]||"",shape:t[1]||""};
      const grab=key=>{ const r=new RegExp("\\b"+key+"\\s+([-\\d.eE+]+)\\s+([-\\d.eE+]+)\\s+([-\\d.eE+]+)\\s+([-\\d]+)\\s+([-\\d.eE+]+)","i").exec(body); return r?{v:+r[1],lo:+r[2],hi:+r[3],link:parseInt(r[4],10),k:+r[5]}:null; };
      c.area=grab("Area"); c.fwhm=grab("MFWHM"); c.pos=grab("Position"); const rs=/\bRSF\s+([-\d.eE+]+)/i.exec(body); c.rsf=rs?+rs[1]:null; const ms=/\bMASS\s+([-\d.eE+]+)/i.exec(body); c.mass=ms?+ms[1]:null;
      out.comps.push(c); } });
  return out.regions.length||out.comps.length||out.calib.length?out:null;
}
/* VAMAS 실험 방식·블록 변수 (v1.7 — 깊이 profile·각도·map 묶음을 만들 수 있게).
   expVals: 실험 변수(예: Etch Time s, Etch Level, Angle deg)의 값(문자열 그대로, expNum 은 숫자), x·y: MAP/MAPDP 위치,
   takeOffPolar/Az: 분석기 축 take-off 각(ARXPS), tilt: 시료 기울기. 숫자가 VAMAS '값 없음'(1e37)이면 null. */
function vamasBlockMeta(H,B,index){
  const num=v=>{ const x=typeof v==="number"?v:parseFloat(v); return isFinite(x)&&!isNullV(x)?x:null; };
  return {mode:H.mode,scan:H.scan,nReg:H.nReg!=null?H.nReg:null,nPos:H.nPos!=null?H.nPos:null,nX:H.nX!=null?H.nX:null,nY:H.nY!=null?H.nY:null,
    expVars:H.expVars.slice(),expVals:B.expVals.slice(),expNum:B.expVals.map(num),blockIndex:index,
    x:B.x!=null?num(B.x):null,y:B.y!=null?num(B.y):null,fov:B.fovX!=null?[num(B.fovX),num(B.fovY)]:null,linescan:B.linescan||null,
    takeOffPolar:num(B.toPolar),takeOffAz:num(B.toAz),srcPolar:num(B.polar),srcAz:num(B.azim),tilt:num(B.tilt),tiltAz:num(B.tiltAz),rotation:num(B.rot),
    sputter:B.sput||null,sputIon:B.sputIon||null,add:B.add.slice()};
}
function seriesTag(H,B){
  const a=[]; H.expVars.forEach((v,i)=>a.push(v.label+" "+B.expVals[i]+(v.unit&&v.unit!=="n"?" "+v.unit:"")));
  if(H.mode==="MAP"||H.mode==="MAPDP") a.push("x "+B.x+" y "+B.y);
  return a.join(", ");
}
/* VAMAS 머리 요약: readFile 결과의 vamas 칸 */
function vamasHeader(text){
  const L=text.split(/\r?\n/); let p=0; const nx=()=>(L[p++]||"").trim();
  const H={}; nx(); H.inst=nx(); H.model=nx(); H.operator=nx(); H.exp=nx(); const nc=parseInt(nx(),10); for(let i=0;i<nc;i++) nx();
  H.mode=nx(); H.scan=nx(); if(["MAP","MAPDP","NORM","SDP"].includes(H.mode)) H.nReg=parseInt(nx(),10);
  if(H.mode==="MAP"||H.mode==="MAPDP"){ H.nPos=parseInt(nx(),10); H.nX=parseInt(nx(),10); H.nY=parseInt(nx(),10); }
  const ne=parseInt(nx(),10); H.expVars=[]; for(let i=0;i<ne;i++) H.expVars.push({label:nx(),unit:nx()});
  return H;
}

/* ============================================================
   6b) 여러 열 텍스트: 에너지 1열 + 스펙트럼 여러 열 (깊이 profile·각도 series 내보내기)
   머리줄(숫자 줄 앞의 마지막 글 줄)을 열 이름으로 쓰고, 이름 안의 첫 숫자를 좌표(식각 시간·각도 등)로 읽는다.
   단일 스펙트럼 읽기(textParse)는 바꾸지 않는다 — 이 함수는 series 화면이 직접 부른다.
   ============================================================ */
function textMultiParse(text,fileName,opt){
  opt=opt||{}; const L=text.split(/\r?\n/); const cd=commaDecimal(L); const rows=[]; let hdr=[], lastHdr="";
  const split=s=>cd?s.split(/[\t;]+|\s{2,}/).map(t=>t.trim()).filter(t=>t!=="").map(t=>t.replace(",",".")):s.split(/[\t,;]+|\s+/).filter(t=>t!=="");
  for(const raw of L){ const s=raw.trim(); if(!s) continue; const parts=split(s); const nums=parts.map(Number);
    if(parts.length>=2&&nums.every(isFinite)) rows.push(nums); else if(!rows.length){ hdr.push(s); lastHdr=raw; } }
  if(rows.length<5) throw new Error("숫자 행이 5개 미만입니다");
  const nc=Math.min(...rows.map(r=>r.length)); if(nc<2) throw new Error("열이 2개 미만입니다");
  let names=lastHdr?lastHdr.split(/\t|;|,(?!\d)/).map(t=>t.trim()).filter(t=>t!==""):[]; if(names.length!==nc) names=lastHdr?lastHdr.trim().split(/\s{2,}|\t/).map(t=>t.trim()):[];
  if(names.length!==nc) names=Array.from({length:nc},(_,i)=>i?"col"+i:"E");
  const h=hdr.join(" ").toLowerCase(); const isKE=opt.axis==="KE"?true:opt.axis==="BE"?false:/kinetic|\bke\b/.test(h)&&!/binding|\bbe\b/.test(h);
  const hv=opt.hv||1486.6; let E=rows.map(r=>isKE?hv-r[0]:r[0]); let Y=[]; for(let c=1;c<nc;c++) Y.push(rows.map(r=>r[c]));
  if(E.length>1&&E[0]>E[E.length-1]){ E=E.reverse(); Y=Y.map(y=>y.reverse()); }
  for(let i=1;i<E.length;i++) if(!(E[i]>E[i-1])) throw new Error("에너지 축이 단조롭지 않습니다");
  const labels=names.slice(1); const coord=labels.map(s=>{ const m=/[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?/.exec(s); return m?parseFloat(m[0]):null; });
  return {E:Float64Array.from(E),Y:Y.map(y=>Float64Array.from(y)),labels,coord,axis:isKE?"KE":"BE",hv,file:fileName,header:hdr.slice(0,5),commaDecimal:cd,
    warn:coord.some(v=>v==null)?["열 이름에서 좌표(숫자)를 읽지 못한 열이 있습니다 — 순서 번호를 쓰세요"]:[]};
}

/* ============================================================
   4) BIFF8 (.xls) — Avantage 'Excel 로 내보내기'
   ============================================================ */
function xlsWorkbook(u8){
  const C=cfbOpen(u8);
  const nm=C.has("Workbook")?"Workbook":C.has("Book")?"Book":null;
  if(!nm) throw new Error("Workbook 스트림이 없습니다 (.xls 아님)");
  const w=C.read(nm), D=dv(w);
  const recs=[]; for(let o=0;o+4<=w.length;){ const t=D.getUint16(o,true), n=D.getUint16(o+2,true); recs.push({t,o:o+4,n}); o+=4+n; }
  // SST (CONTINUE 경계에서 문자열이 갈라지는 경우 처리)
  let sst=[]; const sheets=[];
  for(let i=0;i<recs.length;i++){
    const r=recs[i];
    if(r.t===0x0085){ const pos=D.getUint32(r.o,true); const cch=w[r.o+6], fl=w[r.o+7]; const s=fl&1?utf16(w.subarray(r.o+8,r.o+8+2*cch)):latin1(w.subarray(r.o+8,r.o+8+cch)); sheets.push({name:s,pos,type:w[r.o+5]}); }
    if(r.t===0x00FC){
      const segs=[w.subarray(r.o,r.o+r.n)]; let j=i+1; while(j<recs.length&&recs[j].t===0x003C){ segs.push(w.subarray(recs[j].o,recs[j].o+recs[j].n)); j++; }
      sst=readSST(segs);
    }
  }
  const out=[];
  for(const sh of sheets){
    if(sh.type!==0) continue;          // 워크시트만
    const cells=new Map(); let maxR=-1,maxC=-1; let o=sh.pos; let pendingStr=null;
    const put=(r,c,v)=>{ cells.set(r*1024+c,v); if(r>maxR) maxR=r; if(c>maxC) maxC=c; };
    let first=true;
    while(o+4<=w.length){
      const t=D.getUint16(o,true), n=D.getUint16(o+2,true), q=o+4; o=q+n;
      if(first){ first=false; if(t!==0x0809) break; continue; }
      if(t===0x000A) break;
      if(t===0x0203) put(D.getUint16(q,true),D.getUint16(q+2,true),D.getFloat64(q+6,true));
      else if(t===0x027E) put(D.getUint16(q,true),D.getUint16(q+2,true),rk(D.getUint32(q+6,true)));
      else if(t===0x00BD){ const r=D.getUint16(q,true), c0=D.getUint16(q+2,true), c1=D.getUint16(q+n-2,true); for(let c=c0;c<=c1;c++) put(r,c,rk(D.getUint32(q+4+(c-c0)*6+2,true))); }
      else if(t===0x00FD) put(D.getUint16(q,true),D.getUint16(q+2,true),sst[D.getUint32(q+6,true)]);
      else if(t===0x0204){ const cch=D.getUint16(q+6,true), fl=w[q+8]; put(D.getUint16(q,true),D.getUint16(q+2,true),fl&1?utf16(w.subarray(q+9,q+9+2*cch)):latin1(w.subarray(q+9,q+9+cch))); }
      else if(t===0x0006){ const r=D.getUint16(q,true), c=D.getUint16(q+2,true); if(D.getUint16(q+12,true)===0xFFFF){ if(w[q+6]===0) pendingStr=[r,c]; else if(w[q+6]===1) put(r,c,!!w[q+8]); } else put(r,c,D.getFloat64(q+6,true)); }
      else if(t===0x0207&&pendingStr){ const cch=D.getUint16(q,true), fl=w[q+2]; put(pendingStr[0],pendingStr[1],fl&1?utf16(w.subarray(q+3,q+3+2*cch)):latin1(w.subarray(q+3,q+3+cch))); pendingStr=null; }
      else if(t===0x0205){ /* BOOLERR */ }
    }
    out.push({name:sh.name,get:(r,c)=>cells.get(r*1024+c),nRows:maxR+1,nCols:maxC+1});
  }
  return out;
}
function rk(v){ let x; if(v&2){ x=(v|0)>>2; } else { const b=new DataView(new ArrayBuffer(8)); b.setUint32(4,v&0xFFFFFFFC,true); b.setUint32(0,0,true); x=b.getFloat64(0,true); } return v&1?x/100:x; }
function readSST(segs){
  let si=0, o=0; const cur=()=>segs[si];
  const need=()=>{ if(o>=cur().length&&si+1<segs.length){ si++; o=0; } };
  const u8=()=>{ need(); return cur()[o++]; };
  const u16=()=>{ const a=u8(), b=u8(); return a|(b<<8); };
  const u32=()=>{ const a=u16(), b=u16(); return (a|(b<<16))>>>0; };
  const tot=u32(), uniq=u32(); const out=[];
  for(let k=0;k<uniq;k++){
    if(si>=segs.length) break;
    need(); const cch=u16(); let fl=u8(); let runs=0, ext=0;
    if(fl&8) runs=u16(); if(fl&4) ext=u32();
    let s="", left=cch;
    while(left>0){
      if(o>=cur().length){ if(si+1>=segs.length) break; si++; o=0; fl=cur()[o++]; }  // 이어지는 조각의 첫 바이트 = 인코딩 표지
      const wide=fl&1; const avail=wide?Math.floor((cur().length-o)/2):cur().length-o; const n=Math.min(left,avail);
      if(n<=0){ if(si+1>=segs.length) break; si++; o=0; fl=cur()[o++]; continue; }
      s+=wide?utf16(cur().subarray(o,o+2*n)):latin1(cur().subarray(o,o+n)); o+=wide?2*n:n; left-=n;
    }
    let skip=runs*4+ext; while(skip>0){ need(); const n=Math.min(skip,cur().length-o); o+=n; skip-=n; if(n===0) break; }
    out.push(s);
  }
  return out;
}

/* ---------------- .xlsx (zip + XML) ---------------- */
async function inflateRaw(u8){
  if(typeof DecompressionStream==="undefined") throw new Error("이 브라우저는 압축 해제(DecompressionStream)를 지원하지 않습니다");
  const ds=new DecompressionStream("deflate-raw"); const st=new Blob([u8]).stream().pipeThrough(ds);
  return new Uint8Array(await new Response(st).arrayBuffer());
}
async function zipRead(u8){
  const D=dv(u8); let e=-1; for(let i=u8.length-22;i>=Math.max(0,u8.length-66000);i--) if(D.getUint32(i,true)===0x06054b50){ e=i; break; }
  if(e<0) throw new Error("zip 끝 레코드를 찾지 못했습니다");
  const n=D.getUint16(e+10,true); let cd=D.getUint32(e+16,true); const files={};
  for(let k=0;k<n;k++){
    if(D.getUint32(cd,true)!==0x02014b50) break;
    const meth=D.getUint16(cd+10,true), csz=D.getUint32(cd+20,true), nl=D.getUint16(cd+28,true), xl=D.getUint16(cd+30,true), cl=D.getUint16(cd+32,true), lo=D.getUint32(cd+42,true);
    const name=utf8(u8.subarray(cd+46,cd+46+nl));
    files[name]={meth,csz,lo}; cd+=46+nl+xl+cl;
  }
  return {list:()=>Object.keys(files),has:p=>p in files,read:async p=>{ const f=files[p]; if(!f) return null; const nl=D.getUint16(f.lo+26,true), xl=D.getUint16(f.lo+28,true); const d=u8.subarray(f.lo+30+nl+xl,f.lo+30+nl+xl+f.csz); return f.meth===0?d:await inflateRaw(d); }};
}
const xmlUnesc=s=>s.replace(/&lt;/g,"<").replace(/&gt;/g,">").replace(/&quot;/g,'"').replace(/&apos;/g,"'").replace(/&#(\d+);/g,(m,d)=>String.fromCharCode(+d)).replace(/&#x([0-9a-f]+);/ig,(m,h)=>String.fromCharCode(parseInt(h,16))).replace(/&amp;/g,"&");
async function xlsxWorkbook(u8){
  const Z=await zipRead(u8);
  const txt=async p=>{ const b=await Z.read(p); return b?utf8(b):null; };
  const wb=await txt("xl/workbook.xml"); if(!wb) throw new Error("xl/workbook.xml 이 없습니다");
  const rels=(await txt("xl/_rels/workbook.xml.rels"))||"";
  const relMap={}; rels.replace(/<Relationship\b[^>]*>/g,t=>{ const id=/Id="([^"]+)"/.exec(t), tg=/Target="([^"]+)"/.exec(t); if(id&&tg) relMap[id[1]]=tg[1]; return t; });
  const ss=[]; const sx=await txt("xl/sharedStrings.xml");
  if(sx) sx.replace(/<si>([\s\S]*?)<\/si>/g,(m,b)=>{ let s=""; b.replace(/<t[^>]*>([\s\S]*?)<\/t>/g,(m2,t)=>{ s+=xmlUnesc(t); return m2; }); ss.push(s); return m; });
  const out=[];
  const sheetTags=wb.match(/<sheet\b[^>]*>/g)||[];
  for(const t of sheetTags){
    const nm=xmlUnesc((/name="([^"]*)"/.exec(t)||[])[1]||""); const rid=(/r:id="([^"]+)"/.exec(t)||/\bid="([^"]+)"/.exec(t)||[])[1];
    let tg=relMap[rid]||""; if(!tg) continue; tg=tg.replace(/^\/?xl\//,""); const x=await txt("xl/"+tg); if(!x) continue;
    const cells=new Map(); let maxR=-1,maxC=-1;
    x.replace(/<c\b([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g,(m,a,body)=>{
      const ref=(/r="([A-Z]+)(\d+)"/.exec(a)); if(!ref) return m;
      let c=0; for(const ch of ref[1]) c=c*26+(ch.charCodeAt(0)-64); c--; const r=+ref[2]-1;
      const ty=(/t="([^"]+)"/.exec(a)||[])[1]||"n"; let v=null;
      if(body){ const vv=/<v>([\s\S]*?)<\/v>/.exec(body); if(ty==="s"&&vv) v=ss[+vv[1]]; else if(ty==="inlineStr"){ const tt=/<t[^>]*>([\s\S]*?)<\/t>/.exec(body); v=tt?xmlUnesc(tt[1]):""; } else if(ty==="str"&&vv) v=xmlUnesc(vv[1]); else if(ty==="b"&&vv) v=vv[1]==="1"; else if(vv) v=parseFloat(vv[1]); }
      if(v!=null){ cells.set(r*1024+c,v); if(r>maxR) maxR=r; if(c>maxC) maxC=c; }
      return m;
    });
    out.push({name:nm,get:(r,c)=>cells.get(r*1024+c),nRows:maxR+1,nCols:maxC+1});
  }
  return out;
}

/* Avantage 시트 → 스펙트럼. 'Binding Energy'/'Kinetic Energy' 머리칸 아래가 데이터. 피팅 성분 열이 있으면 함께 보관 */
function avantageSheets(sheets,fileName){
  const out=[], tables=[], warn=[]; let when="";
  for(const sh of sheets){
    const isNum=v=>typeof v==="number"&&isFinite(v);
    let hr=-1, hc=0, isKE=false;
    // 머리칸 = 'Binding/Kinetic Energy' 이고 그 열 아래 몇 줄 안에 숫자가 있는 칸 (v1.8: Kratos→Avantage 내보내기는 위쪽에 'Axis | Kinetic energy | Elements=' 설명 줄이 있어 그 칸을 머리로 잡으면 읽지 못했다)
    const numBelow=(r,c)=>{ for(let k=r+1;k<=Math.min(sh.nRows-1,r+6);k++) if(isNum(sh.get(k,c))) return true; return false; };
    for(let r=0;r<Math.min(sh.nRows,80)&&hr<0;r++) for(let c=0;c<Math.min(sh.nCols,6);c++){ const v=sh.get(r,c); if(typeof v==="string"&&/^(binding|kinetic)\s+energy/i.test(v.trim())&&numBelow(r,c)){ hr=r; hc=c; isKE=/kinetic/i.test(v); break; } }
    if(/^titles$/i.test(sh.name)){ const v=sh.get(0,0); if(typeof v==="string") when=v.trim(); }
    if(hr<0){
      // 표 형태(Peak Table 등)는 그대로 보관
      if(sh.nRows>1&&sh.nCols>1){ const rows=[]; for(let r=0;r<sh.nRows;r++){ const row=[]; for(let c=0;c<sh.nCols;c++){ const v=sh.get(r,c); row.push(v==null?"":v); } if(row.some(v=>v!=="")) rows.push(row); } if(rows.length) tables.push({name:sh.name,rows}); }
      continue;
    }
    // 데이터 시작 행(첫 숫자 행)과 열 이름(머리 아래 줄들)
    let r0=hr+1; while(r0<sh.nRows&&!isNum(sh.get(r0,hc))) r0++;
    const colName=c=>{ const a=[]; for(let r=hr;r<r0;r++){ const v=sh.get(r,c); if(v!=null&&v!=="") a.push(String(v).trim()); } return a.join(" "); };
    const dataCols=[]; for(let c=hc+1;c<Math.max(sh.nCols,hc+2);c++){ if(isNum(sh.get(r0,c))) dataCols.push(c); }
    if(!dataCols.length) continue;
    const ic=dataCols.find(c=>/counts|cps|intensity/i.test(colName(c)))!=null?dataCols.find(c=>/counts|cps|intensity/i.test(colName(c))):dataCols[0];
    const x=[], y=[], others={}; dataCols.filter(c=>c!==ic).forEach(c=>others[colName(c)||("열"+c)]=[]);
    for(let r=r0;r<sh.nRows;r++){ const e=sh.get(r,hc), v=sh.get(r,ic); if(!isNum(e)||!isNum(v)) break; x.push(e); y.push(v); dataCols.filter(c=>c!==ic).forEach(c=>{ const w=sh.get(r,c); others[colName(c)||("열"+c)].push(isNum(w)?w:NaN); }); }
    if(x.length<5) continue;
    let path=""; for(let r=0;r<hr;r++){ const v=sh.get(r,0); if(typeof v==="string"&&/\.vgd/i.test(v)) path=v; }
    let desc=""; for(let r=0;r<hr;r++) for(let c=0;c<sh.nCols;c++){ const v=sh.get(r,c); if(typeof v==="string"&&/experiment\\/i.test(v)) desc=v; }
    const region=sh.name.replace(/\s*(Scan|Snap|Region)$/i,"").trim();
    const hv=1486.68;  // Avantage 내보내기에는 광원 에너지가 없다 — BE 축 그대로 쓰므로 KE 계산용으로만 쓰인다
    const be=isKE?x.map(k=>hv-k):x;
    const unit=colName(ic);
    const w2=[]; if(!/counts\s*\/\s*s|cps/i.test(unit)) w2.push("세기 단위 '"+unit+"' — CPS 로 간주");
    if(be.length>1&&be[0]>be[be.length-1]) Object.keys(others).forEach(k=>others[k].reverse());   // v1.8: mkSpec 가 BE 오름차순으로 뒤집으므로 추가 열도 같은 순서로 (전에는 뒤집히지 않아 채널이 어긋났다)
    const ext=Object.keys(others).length?{avantageCols:others}:null;
    const pt=pathPoint(path);
    out.push(mkSpec({name:region,region,file:fileName,format:"Avantage XLS",hv,hvSrc:"Avantage 기본값 1486.68 (파일에 없음)",be,y,
      sample:String(fileName||"").replace(/^.*[\\/]/,"").replace(/\.xlsx?$/i,"")+(pt?" · "+pt:""),
      meta:{origPath:path,expDesc:desc,yUnit:unit,when,point:pt},warn:w2,extra:ext}));
  }
  // v1.8: 시트 이름이 'Sheet5' 처럼 일반 이름이면 Peak Table 의 Start/End BE 가 Avantage 배경 열(0 아닌 구간)을 감싸는 줄의 Name 으로 region 이름을 붙인다
  const pt=tables.find(t=>/peak\s*table/i.test(t.name));
  if(pt){ const hi=pt.rows.findIndex(r=>r.some(v=>/^name$/i.test(String(v).trim()))&&r.some(v=>/^start\s*be$/i.test(String(v).trim())));
    if(hi>=0){ const H=pt.rows[hi].map(v=>String(v).trim()), cN=H.findIndex(v=>/^name$/i.test(v)), cS=H.findIndex(v=>/^start\s*be$/i.test(v)), cE=H.findIndex(v=>/^end\s*be$/i.test(v));
      const peaks=pt.rows.slice(hi+1).filter(r=>typeof r[cS]==="number"&&typeof r[cE]==="number"&&String(r[cN]).trim()).map(r=>({name:String(r[cN]).trim(),lo:Math.min(r[cS],r[cE]),hi:Math.max(r[cS],r[cE])}));
      out.forEach(sp=>{ if(!/^sheet\s*\d+$/i.test(sp.region)) return; const cols=sp.extra&&sp.extra.avantageCols; const k=cols&&Object.keys(cols).find(n=>/backgnd/i.test(n)); if(!k) return;
        const b=cols[k]; let lo=Infinity, hi2=-Infinity; for(let i=0;i<b.length;i++) if(b[i]&&isFinite(b[i])){ const e=sp.be[i]; if(e<lo) lo=e; if(e>hi2) hi2=e; }
        if(!(hi2>lo)) return; const st=Math.abs(sp.step)||0.05; const m=peaks.filter(p=>Math.abs(lo-p.lo)<=st/2+1e-6&&Math.abs(p.hi-hi2)<=st/2+1e-6);   // Avantage 배경 끝 = Start/End BE 에 가장 가까운 채널
        if(m.length===1){ sp.meta.sheet=sp.region; const rg=parseRegionName(m[0].name); sp.reg=rg; sp.name=m[0].name; sp.region=rg.el||rg.kind==="survey"?rg.label:m[0].name; } }); } }
  return {spectra:out,tables,warn};
}

/* ============================================================
   5) Avantage .avg (DataSpace_BatchDump 텍스트)
   ============================================================ */
function avgParse(text,fileName){
  const L=text.split(/\r?\n/); const P={}; let axes=[], data=[], cur=null, inData=false;
  for(const raw of L){
    const s=raw.trim();
    const pm=/^([A-Z_]+(?:\[\d+\])?)\s*:\s*VT_(\w+)\s*=\s*(.*)$/.exec(s);
    if(pm){ let v=pm[3].trim(); if(/^'.*'$/.test(v)) v=v.slice(1,-1); else if(/^(R4|R8|I4|I2|UI4)$/.test(pm[2])) v=parseFloat(v); P[pm[1]]=v; continue; }
    const am=/^(\d+)=\s*([-\d.eE+]+),\s*([-\d.eE+]+),\s*(\d+),\s*(\w+)\s*,\s*(\w+)?\s*,?\s*'([^']*)'\s*,\s*'([^']*)'\s*,\s*'([^']*)'/.exec(s);
    if(am&&!inData){ axes.push({i:+am[1],start:+am[2],step:+am[3],n:+am[4],type:am[5],sym:am[7],unit:am[8],label:am[9]}); continue; }
    if(/^\$DATA\s*=/.test(s)){ if(cur&&cur.length) data.push(cur); cur=[]; inData=true; continue; }
    if(inData&&/^LIST@/.test(s)){ const q=s.split("=").slice(1).join("="); q.split(",").forEach(t=>{ t=t.trim(); if(t==="#empty#") cur.push(NaN); else if(t!=="") { const v=parseFloat(t); if(isFinite(v)||t==="NaN") cur.push(v); } }); continue; }
    if(inData&&/^\$/.test(s)&&!/^\$DATA/.test(s)) inData=false;
  }
  if(cur&&cur.length) data.push(cur);
  const eax=axes.find(a=>/ENERGY/i.test(a.type))||axes[0];
  if(!eax) throw new Error(".avg 에서 에너지 축($SPACEAXES)을 찾지 못했습니다");
  const flat=[].concat(...data); const nsp=Math.max(1,Math.floor(flat.length/eax.n));
  const hv=+P.DS_SOPROPID_ENERGY>100?+P.DS_SOPROPID_ENERGY:1486.68;
  const dwell=+P.DS_ACPROPID_ACQ_TIME||null, scans=+P.DS_ACPROPID_PERIODS||null;
  const c=dwell>0&&scans>0?dwell*scans:null;
  const tfc=[0,1,2,3].map(k=>P["DS_ANPROPID_TXFN_COEFF["+k+"]"]); const tfOK=tfc.every(v=>typeof v==="number"&&isFinite(v));
  const pe=+P.DS_ANPROPID_PASS||null;
  const isBE=/binding|BE/i.test(eax.label+" "+eax.sym);
  const ke=Array.from({length:eax.n},(_,i)=>eax.start+i*eax.step);
  const be=isBE?ke:ke.map(k=>hv-k);
  const base=(P.DS_EXT_SUPROPID_TITLE||String(fileName).replace(/^.*[\\/]/,"").replace(/\.avg$/i,"")).replace(/\s*(Scan|Snap|Region)$/i,"").trim();
  const meta={instrument:P.DS_GEPROPID_INSTRUMENT||"",passEnergy:pe,dwell,scans,workFunction:+P.DS_ANPROPID_WORK_FTN||null,lensMode:P.DS_ANPROPID_LENS_MODE_NAME||"",
    tf:tfOK?{model:"thermo",coef:tfc,pe,src:".avg 의 DS_ANPROPID_TXFN_COEFF (형태는 추정)"}:null,props:P};
  const out=[];
  for(let s=0;s<nsp;s++){ const cnt=flat.slice(s*eax.n,(s+1)*eax.n); out.push(mkSpec({name:base+(nsp>1?" #"+(s+1):""),region:base,file:fileName,format:"Avantage AVG",hv,hvSrc:".avg",be,y:cnt.map(v=>c?v/c:v),yRaw:c?cnt:null,c,cSrc:c?`.avg: ${dwell} s × ${scans}`:null,meta})); }
  return out;
}

/* ============================================================
   6) 텍스트 2~3열
   ============================================================ */
/* 소수점이 쉼표인 표("284,50;1234,5" · "284,50<TAB>1234,5")인가: 공백·세미콜론으로 나눈 칸이 2 개 이상이고, 칸이 모두
   쉼표 소수(점 없음)·정수이며 쉼표 소수가 하나 이상인 줄이 숫자 줄의 대부분. "1200,345" 처럼 칸이 하나뿐인 줄은 보통 CSV 라 제외. */
function commaDecimal(lines){
  const NUMC=/^[-+]?\d+(,\d+)?([eE][-+]?\d+)?$/; let yes=0, num=0;
  for(const raw of lines){ const s=raw.trim(); if(!s||!/\d/.test(s)) continue; const parts=s.split(/[\s;]+/).filter(t=>t!=="");
    if(parts.length<2||!parts.slice(0,2).every(t=>/^[-+]?[\d.,]+([eE][-+]?\d+)?$/.test(t))) continue; num++;
    if(parts.every(t=>NUMC.test(t))&&parts.some(t=>t.includes(","))) yes++; }
  return num>=5&&yes>=0.8*num;
}
function textParse(text,fileName,opt){
  opt=opt||{};
  const L=text.split(/\r?\n/); const rows=[]; let hdr=[];
  const cd=commaDecimal(L);
  for(const raw of L){
    const s=raw.trim(); if(!s) continue;
    const parts=cd?s.split(/[\s;]+/).filter(t=>t!=="").map(t=>t.replace(",",".")):s.split(/[\s,;]+/).filter(t=>t!=="");
    const nums=parts.map(Number);
    if(parts.length>=2&&nums.slice(0,2).every(isFinite)) rows.push(nums.filter(isFinite));
    else if(!rows.length) hdr.push(s);
  }
  if(rows.length<5) throw new Error("숫자 2열 이상인 행이 5개 미만입니다");
  const h=hdr.join(" ").toLowerCase();
  let isKE=opt.axis==="KE"?true:opt.axis==="BE"?false:/kinetic|\bke\b/.test(h)&&!/binding|\bbe\b/.test(h);
  const hv=opt.hv||1486.6;
  const x=rows.map(r=>r[0]), y=rows.map(r=>r[1]);
  const warn=[];
  const axisGuess=!/binding|kinetic|\bbe\b|\bke\b/.test(h)&&!opt.axis;
  if(axisGuess) warn.push("머리줄에 축 이름이 없어 결합에너지(BE)로 간주 — 운동에너지면 Info 탭에서 KE 로 바꾸세요");
  if(cd) warn.push("소수점이 쉼표(comma decimals)인 표로 읽었습니다");
  const bs=bgSubtracted(y); if(bs) warn.push(bs);
  const be=isKE?x.map(k=>hv-k):x;
  const base=String(fileName||"").replace(/^.*[\\/]/,"");
  const reg=parseRegionName(base.replace(/^\[.*?\]\s*/,""));
  try{
    return [mkSpec({name:reg.label,region:reg.label,file:fileName,format:"텍스트",hv,hvSrc:opt.hv?"직접 입력":"기본값(Al Kα 1486.6) — 텍스트에는 없음",be,y,sample:sampleFromName(base),meta:{header:hdr.slice(0,5),axis:isKE?"KE":"BE",axisGuess,commaDecimal:cd,bgSub:!!bs},warn})];
  }catch(e){
    // 단조롭지 않은 축: 쉼표 소수를 쉼표 구분으로 잘못 나눈 경우가 흔하다 — 그 가능성을 알려 준다
    if(/단조/.test(e.message)&&!cd&&/\d,\d/.test(text.slice(0,4000))) throw new Error(e.message+" — 소수점이 쉼표(comma decimals)인 파일이면 칸을 탭·세미콜론으로 나누어 저장하세요");
    throw e;
  }
}
/* 배경을 이미 뺀 자료인가: 음수 세기가 있거나(측정 CPS 는 음수가 될 수 없음), 10 % 넘는 점이 최댓값의 0.5 % 이하 */
function bgSubtracted(y){
  const n=y.length; if(n<8) return null; let mx=-Infinity, neg=0, low=0; for(const v of y) mx=Math.max(mx,v); if(!(mx>0)) return null;
  for(const v of y){ if(v<0) neg++; if(v<=0.005*mx) low++; }
  if(neg>0||low>0.1*n) return "음수·0 근처 세기 — 배경을 이미 뺀 자료로 보입니다. Background ‘None’ + Weighting ‘Constant’ 로 두었습니다";
  return null;
}
/* 텍스트 자료의 축 해석을 바꾼다(다시 읽기와 같음): 지금 be 에서 파일의 x 를 되살리고 새 축(BE/KE)·hν 로 다시 계산.
   sp 를 직접 바꾸고 true/false(바뀜) 를 돌려준다. */
function reaxis(sp,axis,hv){
  const old=(sp.meta&&sp.meta.axis)||"BE", hv0=sp.hv, hv1=hv>0?+hv:hv0, n=sp.n;
  if(old===axis&&Math.abs(hv1-hv0)<1e-9) return false;
  const x=Array.from(sp.be,b=>old==="KE"?hv0-b:b);
  let be=x.map(v=>axis==="KE"?hv1-v:v), y=Array.from(sp.y), yr=sp.yRaw?Array.from(sp.yRaw):null, tf=sp.tfArr?Array.from(sp.tfArr):null;
  if(n>1&&be[0]>be[n-1]){ be.reverse(); y.reverse(); if(yr) yr.reverse(); if(tf) tf.reverse(); }
  sp.be=Float64Array.from(be); sp.y=Float64Array.from(y); if(yr) sp.yRaw=Float64Array.from(yr); if(tf) sp.tfArr=Float64Array.from(tf);
  sp.step=n>1?(sp.be[n-1]-sp.be[0])/(n-1):0; sp.hv=hv1; if(Math.abs(hv1-hv0)>1e-9) sp.hvSrc="직접 입력";
  sp.meta=Object.assign({},sp.meta||{},{axis,axisGuess:false});
  sp.warn=(sp.warn||[]).filter(w=>!/머리줄에 축 이름/.test(w));
  return true;
}

/* ============================================================
   진입점: 파일 이름 + 바이트 → {spectra, tables, warn}
   ============================================================ */
async function readFile(name,u8,opt){
  const r=await readFile0(name,u8,opt);
  (r.spectra||[]).forEach(sp=>{ if(sp.meta&&!sp.meta.bgSub&&!(sp.reg&&sp.reg.kind==="survey")){ const w=bgSubtracted(sp.y); if(w){ sp.meta.bgSub=true; sp.warn.push(w); } } });
  return r;
}
async function readFile0(name,u8,opt){
  const low=String(name).toLowerCase();
  const head=latin1(u8.subarray(0,Math.min(u8.length,200)));
  if(cfbIs(u8)){
    let C; try{ C=cfbOpen(u8); }catch(e){ throw e; }
    if(C.has("VGData")) return {spectra:vgdParse(u8,name),tables:[],warn:[]};
    if(C.has("Workbook")||C.has("Book")){ const r=avantageSheets(xlsWorkbook(u8),name); if(!r.spectra.length) throw new Error("엑셀 파일에서 에너지·세기 열을 찾지 못했습니다 ('Binding Energy' 머리칸)"); return r; }
    throw new Error("알 수 없는 OLE2 파일 (VGData·Workbook 스트림 없음)");
  }
  if(u8[0]===0x50&&u8[1]===0x4B){ const r=avantageSheets(await xlsxWorkbook(u8),name); if(!r.spectra.length) throw new Error("xlsx 에서 에너지·세기 열을 찾지 못했습니다"); return r; }
  const text=utf8(u8);
  if(/^VAMAS/i.test(head.trim())) return {spectra:vamasParse(text,name),tables:[],warn:[],vamas:vamasHeader(text)};
  if(/\.avg$/i.test(low)||/\$SPACEAXES|DS_EXT_SUPROPID/.test(text.slice(0,20000))) return {spectra:avgParse(text,name),tables:[],warn:[]};
  return {spectra:textParse(text,name,opt),tables:[],warn:[]};
}

const XPS_IO={SYMS,readFile,casaCommentParse,textMultiParse,vamasHeader,vamasBlockMeta,reaxis,bgSubtracted,commaDecimal,pathSample,pathPoint,vgdParse,vamasParse,xlsWorkbook,xlsxWorkbook,avantageSheets,avgParse,textParse,cfbOpen,cfbIs,propsetParse,parseRegionName,sampleFromName,mkSpec};
if(typeof module!=="undefined"&&module.exports) module.exports=XPS_IO; else root.XPS_IO=XPS_IO;
})(typeof self!=="undefined"?self:this);
