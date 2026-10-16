/* ============================================================
   INNO Analytics — 광학 현미경 형식 읽기 (IMG_IO 에 덧붙임)
   · Nikon ND2 (NIS-Elements, 형식 v2.x·v3.x): chunk map, CLx Lite Variant(이진)·CLxVariant(XML) 메타데이터,
     비압축·lossless(zlib) 프레임, 실험 loop(T·NETime·XY·Z) → 위치별 영상, 채널 이름·색, 보정(µm/px), z step.
     해석 규칙은 Python `nd2` 패키지(T. Lambert, BSD-3) 의 modern reader 를 따랐다 — 시험도 그 패키지가 읽은 배열과 대조.
     v1.0(JPEG 2000 로 저장한 NIS-Elements 3 이전 형식)은 읽지 않는다.
   · Zeiss CZI (ZISRAW): subblock directory, Gray8/16/32f/Bgr24/48/96f, 비압축·zstd0·zstd1(HiLoByteUnpack)·JPG,
     pyramid 층 0 만, mosaic tile 합성(M 순서, 뒤 tile 이 위 — libCZI 와 같음), scene 별 영상, Scaling·채널 이름·색.
     JPEG XR 압축(slide scanner 의 손실 압축)은 읽지 않는다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IO=G.IMG_IO; if(!IO) throw new Error("IMG_IO 가 먼저 있어야 합니다");
const tdu=new TextDecoder("utf-8"), td16=new TextDecoder("utf-16le");
const LUT_BY_RGB=[["Red",[255,0,0]],["Green",[0,255,0]],["Blue",[0,0,255]],["Cyan",[0,255,255]],["Magenta",[255,0,255]],["Yellow",[255,255,0]],["Grays",[255,255,255]]];
function nearestLut(r,g,b){ const m=Math.max(r,g,b); if(!(m>0)) return null; const n=[r/m*255,g/m*255,b/m*255]; let best=null,bd=1e9;
  for(const [nm,c] of LUT_BY_RGB){ const d=(c[0]-n[0])**2+(c[1]-n[1])**2+(c[2]-n[2])**2; if(d<bd){ bd=d; best=nm; } } return best; }
const u64=(dv,o)=>dv.getUint32(o,true)+dv.getUint32(o+4,true)*4294967296;
const i64=(dv,o)=>dv.getUint32(o,true)+dv.getInt32(o+4,true)*4294967296;
// 광학 현미경은 µm 가 관례 (10 nm 보다 작은 화소 = SEM 등은 nm)
function calOpt(pwm,phm,src){ const u=pwm>=1e-8?["µm",1e6]:["nm",1e9]; return {pw:pwm*u[1],ph:(phm||pwm)*u[1],unit:u[0],src}; }
function memCheck(n,W,H,ctx,warn){ const bytes=n*W*H*4; const lim=(ctx&&ctx.maxBytes)||1.6e9; return bytes<=lim?null:`${n} 장 × ${W}×${H} = ${(bytes/1e9).toFixed(2)} GB (32-bit 메모리) — 한도 ${(lim/1e9).toFixed(1)} GB`; }

/* ============================================================ ND2 ============================================================ */
const ND2_MAGIC=0x0ABECEDA;
// CLx Lite Variant (이진) — nd2/_parse/_clx_lite.py 의 json_from_clx_lite_variant (lists_to_indexed_dicts=True, strip_prefix=False)
function lvDecode(u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength);
  function str(p){ let e=p; while(e+1<u8.length&&!(u8[e]===0&&u8[e+1]===0)) e+=2; return [td16.decode(u8.subarray(p,e)),e+2]; }
  function parse(p,end,count){ const out={}; const emp=[];
    for(let k=0;k<count&&p<end;k++){ const curs=p; const type=u8[p], nl=u8[p+1]; p+=2; if(type===0||type===10) throw new Error("ND2 메타데이터: 알 수 없는 형 "+type);
      if(type===76) throw new Error("ND2 메타데이터: 안쪽 압축 블록");
      const name=td16.decode(u8.subarray(p,p+2*nl)).replace(/\0$/,""); p+=2*nl; let v;
      switch(type){
        case 1: v=!!u8[p]; p+=1; break; case 2: v=dv.getInt32(p,true); p+=4; break; case 3: v=dv.getUint32(p,true); p+=4; break;
        case 4: v=i64(dv,p); p+=8; break; case 5: v=u64(dv,p); p+=8; break; case 6: v=dv.getFloat64(p,true); p+=8; break; case 7: v=u64(dv,p); p+=8; break;
        case 8: { const [s,q]=str(p); v=s; p=q; break; }
        case 9: { const n=u64(dv,p); p+=8; v=Array.from(u8.subarray(p,p+n)); p+=n; break; }
        case 11: { const cnt=dv.getUint32(p,true), len=u64(dv,p+4); p+=12; const cend=curs+len; const val=parse(p,cend,cnt); p=cend+cnt*8;
          const ks=Object.keys(val); if(ks.length===1&&ks[0]===""){ const x=val[""]; const L=Array.isArray(x)?x:[x]; v={}; L.forEach((y,i)=>v["i"+String(i).padStart(10,"0")]=y); } else v=val; break; }
        default: throw new Error("ND2 메타데이터: 형 "+type);
      }
      if(name==="") emp.push(v); else out[name]=v; }
    if(emp.length) out[""]=emp.length===1?emp[0]:emp; // nd2: 이름 없는 항목 = 목록
    return out; }
  return parse(0,u8.length,1); // nd2 패키지와 같이 맨 위 항목 하나
}
// CLxVariant (XML, v2.x) — nd2/_parse/_clx_xml.py
function xmlVariant(s){
  s=s.replace(/^﻿/,"").replace(/^<\?xml[^>]*\?>/,""); let p=0;
  const attrRe=/([\w:.-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g;
  const unesc=t=>t.replace(/&#x([0-9a-fA-F]+);/g,(_,h)=>String.fromCodePoint(parseInt(h,16))).replace(/&#(\d+);/g,(_,d)=>String.fromCodePoint(+d)).replace(/&lt;/g,"<").replace(/&gt;/g,">").replace(/&quot;/g,'"').replace(/&apos;/g,"'").replace(/&amp;/g,"&");
  function tagEnd(q){ let quote=0; for(let i=q;i<s.length;i++){ const ch=s.charCodeAt(i); if(quote){ if(ch===quote) quote=0; } else if(ch===34||ch===39) quote=ch; else if(ch===62) return i; } return -1; } // 따옴표 안의 '>' 는 건너뜀, 길이 제한 없음
  function node(){ // p 는 '<' 위치
    const nm=/^<([\w:.-]+)/.exec(s.substr(p,256)); const e=nm?tagEnd(p):-1; if(!nm||e<0) throw new Error("ND2 XML 메타데이터 해석 실패 @"+p);
    let body=s.slice(p+nm[0].length,e); const self=body.endsWith("/"); if(self) body=body.slice(0,-1); p=e+1;
    const tag=nm[1], at={}; let a; attrRe.lastIndex=0; while((a=attrRe.exec(body))) at[a[1]]=unesc(a[2]!=null?a[2]:a[3]);
    const kids=[]; if(!self){ for(;;){ const q=s.indexOf("<",p); if(q<0) throw new Error("ND2 XML: 끝 태그 없음"); p=q; if(s.startsWith("</",p)){ p=tagEnd(p)+1; break; } if(s.startsWith("<!--",p)){ p=s.indexOf("-->",p)+3; continue; } kids.push(node()); } }
    return {tag,at,kids}; }
  const CAST={bool:x=>/^(true|1)$/i.test(x),CLxByteArray:x=>x,CLxStringW:x=>x,double:x=>parseFloat(x),float:x=>parseFloat(x),lx_int32:x=>parseInt(x),lx_int64:x=>parseInt(x),lx_uint32:x=>parseInt(x),lx_uint64:x=>parseInt(x),unknown:x=>x};
  function val(n){ const rt=n.at.runtype; if(rt in CAST) return CAST[rt](n.at.value);
    const o={}; n.kids.forEach((c,i)=>{ let k=c.tag; const v=val(c); if(["no_name","","BinaryItem","TextInfoItem","Wavelength","MinSrc","MaxSrc"].includes(k)){ if(v==null||(typeof v==="object"&&!Object.keys(v).length)||v===""||v===0||v===false) return; k="i"+String(i).padStart(10,"0"); } o[k]=v; }); return o; }
  p=s.indexOf("<"); const root=node(); let v=val(root); const ks=v&&typeof v==="object"?Object.keys(v):[]; if(ks.length===1&&ks[0]==="i0000000000") v=v.i0000000000;
  return root.tag==="variant"?v:{[root.tag]:v};
}
function nd2ChunkMap(u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const L=u8.length;
  const sig=tdu.decode(u8.subarray(L-40,L-8)); if(sig!=="ND2 CHUNK MAP SIGNATURE 0000001!") throw new Error("ND2 chunk map 서명이 없습니다 (파일이 잘렸거나 저장 중 끊김)");
  const loc=u64(dv,L-8); if(dv.getUint32(loc,true)!==ND2_MAGIC) throw new Error("ND2 chunk map 위치가 올바르지 않습니다");
  const nl=dv.getUint32(loc+4,true), dl=u64(dv,loc+8); const d0=loc+16+nl; const data=u8.subarray(d0,d0+dl);
  const map={}; let p=0; for(;;){ let q=p; while(q<data.length&&data[q]!==0x21) q++; if(q>=data.length) break; const name=tdu.decode(data.subarray(p,q+1)); if(name==="ND2 CHUNK MAP SIGNATURE 0000001!") break;
    const dd=new DataView(data.buffer,data.byteOffset+q+1,16); map[name]=u64(dd,0); p=q+1+16; }
  return map;
}
function nd2Chunk(u8,off){ const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); if(dv.getUint32(off,true)!==ND2_MAGIC) throw new Error("ND2 chunk 머리가 올바르지 않습니다 @"+off);
  const nl=dv.getUint32(off+4,true), dl=u64(dv,off+8); const s=off+16+nl; return u8.subarray(s,s+dl); }
async function nd2Decode(u8){ if(u8[0]===0x3C) return xmlVariant(tdu.decode(u8)); if(u8[0]===76){ return lvDecode(await IO.inflate(u8.subarray(12))); } return lvDecode(u8); }
const LOOPAX={1:"T",2:"P",4:"Z",8:"T",7:"U"};
function nd2Loops(src,level,dest){ // nd2 _parse.load_experiment
  dest=dest||[]; level=level||0; const loop=nd2Loop(src); if(!loop||loop.count===0) return dest;
  if(loop.type===6) level-=1;
  else if(!dest.length||dest[dest.length-1].nest<level){ loop.nest=level; dest.push(loop); }
  else { const prev=dest[dest.length-1]; if(prev.nest===level&&prev.type===loop.type){ loop.nest=level; if(prev.count<loop.count) dest[dest.length-1]=loop; } }
  const lp=src.uLoopPars||{}; if(lp.pSubLoops&&lp.pSubLoops.i0000000000&&lp.pSubLoops.i0000000000.ppNextLevelEx&&lp.pSubLoops.i0000000000.ppNextLevelEx.i0000000000) dest.push(...nd2Loops(lp.pSubLoops.i0000000000.ppNextLevelEx.i0000000000.SLxExperiment,0,[])); // nd2: 새 목록을 만들어 이어 붙임
  const nx=src.ppNextLevelEx; if(nx&&typeof nx==="object") for(const it of Object.values(nx)) dest=nd2Loops(it,level+1,dest);
  return dest; }
function nd2Loop(e){ if(!e) return null; const t=e.eType||0; let lp=e.uLoopPars; if(!lp||!Object.keys(lp).length||t>8) return null; const ks=Object.keys(lp); if(ks.length===1&&ks[0]==="i0000000000") lp=lp.i0000000000;
  if(t===1){ const n=lp.uiCount||0; return n?{type:1,count:n,period:lp.dPeriod}:null; }
  if(t===2){ let valid=e.pItemValid; if(valid&&!Array.isArray(valid)&&typeof valid==="object") valid=Object.keys(valid).sort().map(k=>valid[k]); if(typeof valid==="string") valid=null;
    let pts=lp.Points?Object.values(lp.Points):(lp.dPosX?Object.keys(lp.dPosX).map(k=>({dPosX:lp.dPosX[k]})):[]); const names=pts.map(q=>q.dPosName||q.pPosName||"");
    let keep=pts.map((_,i)=>i); if(valid&&valid.length) keep=keep.filter(i=>valid[i]); return {type:2,count:keep.length,names:keep.map(i=>names[i])}; }
  if(t===4){ const n=lp.uiCount||0; let st=lp.dZStep||0; if(st===0&&n>1) st=Math.abs((lp.dZHigh||0)-(lp.dZLow||0))/(n-1); return {type:4,count:n,step:Math.abs(st)}; }
  if(t===8){ let pv=lp.pPeriodValid||[]; if(!Array.isArray(pv)&&typeof pv==="object") pv=Object.keys(pv).sort().map(k=>pv[k]); let n=0; const per=Object.values(lp.pPeriod||{}); per.forEach((q,i)=>{ if(pv[i]&&q.uiCount) n+=q.uiCount; }); return {type:8,count:n}; }
  if(t===6){ return {type:6,count:(lp.pPlanes&&lp.pPlanes.uiCount)||lp.uiCount||0}; }
  if(t===7){ return {type:7,count:lp.uiCount||0}; }
  return null; }

async function readND2(name,u8,ctx){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const warn=[], info=[];
  if(dv.getUint32(0,true)===0x0C000000) throw new Error("ND2 v1.0 (JPEG 2000 로 저장한 NIS-Elements 3 이전 형식)은 아직 읽지 않습니다 — NIS-Elements 에서 다시 저장하거나 OME-TIFF 로 내보내세요");
  if(dv.getUint32(0,true)!==ND2_MAGIC) throw new Error("ND2 머리가 올바르지 않습니다");
  const ver=tdu.decode(u8.subarray(48,56)).replace(/\0.*$/,""); const major=+ver[3]||0, minor=+ver[5]||0; info.push(["Format",`Nikon ND2 ${ver}`]);
  const cm=nd2ChunkMap(u8); const v3=major>=3;
  const dec=async k=>{ if(cm[k]==null) return null; return nd2Decode(nd2Chunk(u8,cm[k])); };
  let A=await dec(v3?"ImageAttributesLV!":"ImageAttributes!"); if(!A) throw new Error("ND2: ImageAttributes 가 없습니다"); A=A.SLxImageAttributes||A;
  let M=await dec(v3?"ImageMetadataSeqLV|0!":"ImageMetadataSeq|0!")||{}; M=M.SLxPictureMetadata||M;
  let X=await dec(v3?"ImageMetadataLV!":"ImageMetadata!"); X=X?(X.SLxExperiment||X):null;
  let TI=await dec(v3?"ImageTextInfoLV!":"ImageTextInfo!"); TI=TI?(TI.SLxImageTextInfo||TI):null;
  const W=A.uiWidth, H=A.uiHeight, ncomp=A.uiComp||1, bpc=A.uiBpcInMemory, sig=A.uiBpcSignificant, wb=A.uiWidthBytes||W*ncomp*bpc/8; const bypc=bpc/8;
  if(![8,16,32].includes(bpc)) throw new Error("ND2: 화소당 "+bpc+" bit 는 읽지 않습니다");
  const ecomp=A.eCompression==null?2:A.eCompression; if(ecomp===1) throw new Error("ND2: 손실 압축(lossy) 프레임은 읽지 않습니다");
  const PP=M.sPicturePlanes||{}; const nCh=PP.uiCount||1; const cpc=Math.floor(ncomp/nCh); const rgb=cpc===3;
  const nOut=rgb?nCh:ncomp; // RGB 가 아니면 모든 component 를 채널로 (nd2 패키지 sizes 규칙)
  const rawPlanes=PP.sPlaneNew||PP.sPlane||{}; const pl=Object.values(rawPlanes);
  const chNames=[], chLut=[];
  if(!rgb){ for(let c=0;c<nOut;c++){ const q=pl[Math.floor(c/Math.max(1,cpc))]; chNames.push(q&&q.sDescription||`C${c+1}`); const col=q&&q.uiColor; chLut.push(col!=null?nearestLut(col&255,(col>>8)&255,(col>>16)&255):null); } }
  const loops=X?nd2Loops(X):[]; const seq=loops.reduce((a,l)=>a*l.count,1);
  if(A.uiSequenceCount>seq) warn.push(`ND2: 프레임 ${A.uiSequenceCount} 장 중 실험 loop 이 정의한 ${seq} 장만 읽습니다`);
  const axes=loops.map(l=>LOOPAX[l.type]||"U");
  const zi=axes.indexOf("Z"), pi=axes.indexOf("P"); const tIdx=axes.map((a,i)=>(a==="T"||a==="U")?i:-1).filter(i=>i>=0);
  const nz=zi>=0?loops[zi].count:1, np=pi>=0?loops[pi].count:1, nt=tIdx.reduce((a,i)=>a*loops[i].count,1);
  if(tIdx.length>1) warn.push(`ND2: 시간형 loop ${tIdx.length} 개(${tIdx.map(i=>axes[i]).join("·")})를 하나의 t 축으로 이었습니다`);
  const cal=(M.bCalibrated&&M.dCalibration>0)?calOpt(M.dCalibration*1e-6,M.dCalibration*1e-6,"ND2 dCalibration (µm/px)"):null;
  const zStep=zi>=0?loops[zi].step:0;
  if(M.wsObjectiveName) info.push(["Objective",`${M.wsObjectiveName}${M.dObjectiveNA>0?`, NA ${M.dObjectiveNA}`:""}`]);
  if(sig&&sig!==bpc) info.push(["Bit depth",`${sig}-bit significant (${bpc}-bit 저장)`]);
  if(tIdx.length&&loops[tIdx[0]].period>0) info.push(["Frame interval",`${+(loops[tIdx[0]].period/1000).toPrecision(6)} s`]);
  if(TI){ for(const [k,l] of [["TextInfoItem_5","Description"],["TextInfoItem_13","Optics"],["TextInfoItem_6","Capturing"],["TextInfoItem_9","Date"]]) if(TI[k]) info.push([l,String(TI[k]).slice(0,300)]); }
  info.push(["Dimensions",`${axes.map((a,i)=>a+loops[i].count).join(" ")}${axes.length?" ":""}C${nCh}${rgb?" RGB":""} · ${W}×${H} · ${bpc===32?"32-bit float":bpc+"-bit"}${ecomp===0?" · lossless":""}`]);
  // 위치·Z·T 고르기 (큰 파일) / 메모리
  const all=n=>[...Array(n).keys()]; let sel={p:all(np),z:all(nz),t:all(nt)}; const per1=(rgb?3*nCh:nOut);
  if(ctx&&ctx.choose){ const r=await ctx.choose({kind:"nd2",name,positions:np,posNames:pi>=0?loops[pi].names:[],nz,nt,nc:nOut,w:W,h:H,bytesPerPlane:W*H*4,planesPer:per1}); if(!r) return {images:[],warn:[],cancel:true}; sel=Object.assign(sel,r); }
  const zs=sel.z, ts=sel.t, nzo=zs.length, nto=ts.length; const zMap=new Map(zs.map((z,i)=>[z,i])), tMap=new Map(ts.map((t,i)=>[t,i]));
  const mem=memCheck(per1*nzo*nto*sel.p.length,W,H,ctx); if(mem) throw new Error("ND2: 너무 큼 — "+mem+". 위치·Z·T 일부만 고르세요");
  if(nzo<nz) warn.push(`ND2: Z ${nz} 장 중 ${nzo} 장(${zs[0]+1}–${zs[nzo-1]+1}${nzo>1&&zs[1]-zs[0]>1?`, ${zs[1]-zs[0]} 장마다`:""})만 읽었습니다`); if(nto<nt) warn.push(`ND2: T ${nt} 장 중 ${nto} 장만 읽었습니다`);
  const zUni=nzo>1&&zs.every((z,i)=>i===0||z-zs[i-1]===zs[1]-zs[0]); const zStepOut=zStep>0?(nzo>1?(zUni?zStep*(zs[1]-zs[0]):0):zStep):0;
  const getFrame=async s=>{ const k="ImageDataSeq|"+s+"!"; if(cm[k]==null) return null; let d=nd2Chunk(u8,cm[k]).subarray(8); if(ecomp===0) d=await IO.inflate(d); return d; };
  const out=new Map(sel.p.map(p=>[p,new Array(per1*nzo*nto)])); let missing=0;
  for(let s=0;s<seq;s++){ // s 를 loop 색인으로 (바깥 loop 먼저 — nd2 패키지 loop_indices)
    let r=s; const ix=new Array(loops.length); for(let i=loops.length-1;i>=0;i--){ ix[i]=r%loops[i].count; r=Math.floor(r/loops[i].count); }
    const p=pi>=0?ix[pi]:0; if(!out.has(p)) continue; const z0=zi>=0?ix[zi]:0; let t0=0; for(const i of tIdx) t0=t0*loops[i].count+ix[i];
    if(!zMap.has(z0)||!tMap.has(t0)) continue; const z=zMap.get(z0), t=tMap.get(t0);
    const d=await getFrame(s); const fdv=d?new DataView(d.buffer,d.byteOffset,d.byteLength):null; if(!d) missing++;
    if(d&&d.length<wb*(H-1)+W*ncomp*bypc) throw new Error(`ND2: 프레임 ${s} 자료가 짧습니다`);
    const arr=out.get(p);
    for(let k=0;k<per1;k++){ const P=new Float32Array(W*H);
      if(d){ // 화소 (y,x) 의 component k: y*widthBytes + (x*nComp + k)*bytes  (nd2: (H, W, C, comp) 와 같은 순서)
        for(let y=0;y<H;y++){ const ro=y*wb; for(let x=0;x<W;x++){ const o=ro+(x*ncomp+k)*bypc; P[y*W+x]=bpc===8?d[o]:bpc===16?fdv.getUint16(o,true):fdv.getFloat32(o,true); } } }
      if(rgb){ const c=Math.floor(k/3), comp=k%3; arr[comp+3*(c+nCh*(z+nzo*t))]=P; } else arr[k+nOut*(z+nzo*t)]=P; }
    if(ctx&&ctx.onProgress&&(s&15)===0) ctx.onProgress(s/seq); await null; }
  if(missing) warn.push(`ND2: 프레임 ${missing} 장이 파일에 없어 0 으로 채웠습니다`);
  const images=[];
  for(const p of sel.p){ const arr=out.get(p); const pn=pi>=0&&np>1?` · ${loops[pi].names[p]||"XY"+(p+1)}`:"";
    if(rgb){ for(let c=0;c<nCh;c++){ const planes=[]; for(let zt=0;zt<nzo*nto;zt++) for(let comp=0;comp<3;comp++) planes.push(arr[comp+3*(c+nCh*zt)]);
        images.push({name:`${name}${pn}${nCh>1?` · ${pl[c]&&pl[c].sDescription||"C"+(c+1)}`:""}`,w:W,h:H,...(bpc===8?{bit:8,rgb:true}:{bit:bpc,rgb:false,chNames:["R","G","B"],chLut:["Red","Green","Blue"]}),nc:3,nz:nzo,nt:nto,order:"czt",planes,cal,calCandidates:[],info:info.slice(),warn:[],zStep:zStepOut}); } }
    else images.push({name:name+pn,w:W,h:H,bit:bpc===32?32:bpc,rgb:false,nc:nOut,nz:nzo,nt:nto,order:"czt",planes:arr,cal,calCandidates:[],info:info.slice(),warn:[],chNames,chLut,zStep:zStepOut});
  }
  return {images,warn};
}
IO.readND2=readND2; IO._lvDecode=lvDecode; IO._xmlVariant=xmlVariant;

/* ============================================================ CZI ============================================================ */
const CZI_PT={0:[8,1,"u8"],1:[16,1,"u16"],2:[32,1,"f32"],3:[8,3,"u8"],4:[16,3,"u16"],8:[32,3,"f32"],9:[8,4,"u8"],12:[32,1,"i32"],13:[32,1,"f64"]};
const CZI_PTN={0:"Gray8",1:"Gray16",2:"Gray32Float",3:"Bgr24",4:"Bgr48",8:"Bgr96Float",9:"Bgra32",10:"Gray64ComplexFloat",11:"Bgr192ComplexFloat",12:"Gray32",13:"Gray64"};
const CZI_COMP={0:"비압축",1:"JPG",2:"LZW",4:"JPEG XR",5:"zstd0",6:"zstd1"};
function cziEntry(dv,u8,o){ // DirectoryEntryDV
  if(u8[o]!==0x44||u8[o+1]!==0x56) throw new Error("CZI directory 항목이 올바르지 않습니다 @"+o);
  const e={pt:dv.getInt32(o+2,true),pos:i64(dv,o+6),part:dv.getInt32(o+14,true),comp:dv.getInt32(o+18,true),pyr:u8[o+22],dims:{}}; const nd=dv.getInt32(o+28,true);
  for(let k=0;k<nd;k++){ const q=o+32+20*k; const dn=String.fromCharCode(...u8.subarray(q,q+4)).replace(/\0.*$/,""); e.dims[dn]={start:dv.getInt32(q+4,true),size:dv.getInt32(q+8,true),stored:dv.getInt32(q+16,true)}; }
  e.len=32+20*nd; return e; }
async function cziPixels(u8,e,ctx){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const o=e.pos;
  if(String.fromCharCode(...u8.subarray(o,o+14))!=="ZISRAWSUBBLOCK") throw new Error("CZI subblock 머리가 올바르지 않습니다 @"+o);
  const d=o+32; const msz=dv.getInt32(d,true), dsz=i64(dv,d+8); const ent=cziEntry(dv,u8,d+16); const hl=Math.max(256,16+ent.len); const p0=d+hl+msz; let px=u8.subarray(p0,p0+dsz);
  const W=e.dims.X.stored, H=e.dims.Y.stored; const [bit,nComp,ty]=CZI_PT[e.pt]; const bpp=(ty==="f64"?8:bit/8)*nComp; const need=W*H*bpp;
  if(e.comp===5||e.comp===6){ if(!G.IMG_ZSTD) throw new Error("zstd 해독기가 없습니다"); let hilo=false;
    if(e.comp===6){ const hs=px[0]; let q=1; while(q<hs){ const typ=px[q]; if(typ===1){ hilo=!!(px[q+1]&1); q+=2; } else throw new Error("CZI zstd1: 알 수 없는 머리 항목 "+typ); } px=px.subarray(hs); }
    let raw=G.IMG_ZSTD.decompress(px); if(hilo){ const n=raw.length>>1; const r2=new Uint8Array(raw.length); for(let i=0;i<n;i++){ r2[2*i]=raw[i]; r2[2*i+1]=raw[n+i]; } raw=r2; } px=raw; }
  else if(e.comp===1){ const bmp=await createImageBitmap(new Blob([px],{type:"image/jpeg"}),{colorSpaceConversion:"none",premultiplyAlpha:"none"}); const cv=new OffscreenCanvas(bmp.width,bmp.height); const cx=cv.getContext("2d"); cx.drawImage(bmp,0,0); const rgba=cx.getImageData(0,0,bmp.width,bmp.height).data;
    const r=new Uint8Array(W*H*nComp); for(let i=0;i<W*H;i++){ if(nComp===1) r[i]=rgba[4*i]; else { r[3*i]=rgba[4*i+2]; r[3*i+1]=rgba[4*i+1]; r[3*i+2]=rgba[4*i]; } } px=r; }
  else if(e.comp!==0) throw new Error(`CZI: ${CZI_COMP[e.comp]||"압축 "+e.comp} subblock 은 아직 읽지 않습니다`);
  if(px.length<need) throw new Error(`CZI subblock 자료가 짧습니다 (${px.length} < ${need})`);
  const pdv=new DataView(px.buffer,px.byteOffset,px.byteLength); const out=[...Array(Math.min(nComp,3))].map(()=>new Float32Array(W*H));
  const bs=bpp/nComp; // 한 component 바이트
  for(let i=0;i<W*H;i++){ for(let c=0;c<out.length;c++){ const oo=i*bpp+c*bs; const v=ty==="u8"?px[oo]:ty==="u16"?pdv.getUint16(oo,true):ty==="f32"?pdv.getFloat32(oo,true):ty==="i32"?pdv.getInt32(oo,true):pdv.getFloat64(oo,true);
    if(nComp>=3) out[2-c][i]=v; else out[c][i]=v; } } // B,G,R → R,G,B
  return {W,H,comps:out}; }
// XML 경로(ImageDocument/Metadata/...)의 요소 내용 — libCZI 와 같이 정해진 자리만 본다 (다른 곳의 같은 이름 요소는 무시)
function xmlPath(xml,pathArr,from,to){ from=from||0; to=to==null?xml.length:to; let depth=0, want=0, start=-1; const re=/<(\/?)([\w:.-]+)((?:[^>"']|"[^"]*"|'[^']*')*?)(\/?)>|<!--[\s\S]*?-->|<\?[\s\S]*?\?>|<!\[CDATA\[[\s\S]*?\]\]>/g; re.lastIndex=from; let m;
  while((m=re.exec(xml))&&m.index<to){ if(!m[2]) continue; const close=m[1]==="/", self=m[4]==="/", tag=m[2];
    if(close){ depth--; if(start>=0&&depth===want-1&&want===pathArr.length) return {inner:xml.slice(start,m.index),open:xml.slice(openAt,start)}; if(want>depth) want=depth; continue; } // 맞은 조상이 닫히면 경로를 되돌림
    if(depth===want&&tag===pathArr[want]){ if(want===pathArr.length-1){ if(self) return {inner:"",open:m[0]}; var openAt=m.index; start=m.index+m[0].length; want++; depth++; continue; } want++; }
    if(!self) depth++; }
  return null; }
function cziMeta(xml,info){
  const g=(re)=>{ const m=xml.match(re); return m?m[1]:null; };
  const sc={}; const S=xmlPath(xml,["ImageDocument","Metadata","Scaling"]);
  if(S){ const re=/<Distance\s+Id="([XYZ])"\s*>([\s\S]*?)<\/Distance>/g; let m; while((m=re.exec(S.inner))){ const v=(m[2].match(/<Value>\s*([^<\s]+)\s*<\/Value>/)||[])[1]; if(v&&+v>0) sc[m[1]]=+v; } }
  const C=xmlPath(xml,["ImageDocument","Metadata","Information","Image","Dimensions","Channels"]); const ch=[];
  if(C){ const re=/<Channel\b([^>]*?)(\/>|>([\s\S]*?)<\/Channel>)/g; let m; while((m=re.exec(C.inner))){ const nm=(m[1].match(/\bName="([^"]*)"/)||[])[1]||(m[1].match(/\bId="([^"]*)"/)||[])[1]||""; const inner=m[3]||"";
      const col=(inner.match(/<Color>\s*#?([0-9A-Fa-f]{6,8})\s*<\/Color>/)||[])[1]; ch.push({name:nm,color:col}); } }
  const DS=xmlPath(xml,["ImageDocument","Metadata","DisplaySetting","Channels"]); if(DS){ const re=/<Channel\b[^>]*>([\s\S]*?)<\/Channel>/g; let m,i=0; while((m=re.exec(DS.inner))){ const col=(m[1].match(/<Color>\s*#?([0-9A-Fa-f]{6,8})\s*<\/Color>/)||[])[1]; if(col&&ch[i]&&!ch[i].color) ch[i].color=col; i++; } }
  const obj=g(/<Objective\b[^>]*Name="([^"]*)"/), na=g(/<LensNA>\s*([\d.]+)/), mag=g(/<NominalMagnification>\s*([\d.]+)/), mic=g(/<Microscope\b[^>]*Name="([^"]*)"/), sys=g(/<System>\s*([^<]+)<\/System>/), date=g(/<AcquisitionDateAndTime>\s*([^<]+)</), app=g(/<Application>[\s\S]*?<Name>([^<]+)<\/Name>[\s\S]*?<Version>([^<]+)</);
  if(mic||sys) info.push(["Instrument",[mic,sys].filter(Boolean).join(" · ")]); if(obj) info.push(["Objective",`${obj}${mag?` (${mag}×)`:""}${na?`, NA ${na}`:""}`]); if(date) info.push(["Acquired",date]); if(app) info.push(["Software",app]);
  return {sc,ch}; }
async function readCZI(name,u8,ctx){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const warn=[], info=[];
  if(String.fromCharCode(...u8.subarray(0,10))!=="ZISRAWFILE") throw new Error("CZI 머리가 올바르지 않습니다");
  const major=dv.getInt32(32,true), minor=dv.getInt32(36,true), fpart=dv.getInt32(80,true); let dirPos=i64(dv,84); const metaPos=i64(dv,92);
  info.push(["Format",`Zeiss CZI ${major}.${minor}`]); if(fpart>0) warn.push("CZI: 여러 파일로 나뉜 파일의 일부입니다 — 이 파일 안의 subblock 만 읽습니다");
  const ents=[];
  if(dirPos>0&&dirPos<u8.length&&String.fromCharCode(...u8.subarray(dirPos,dirPos+15))==="ZISRAWDIRECTORY"){ const n=dv.getInt32(dirPos+32,true); let o=dirPos+32+128; for(let i=0;i<n;i++){ const e=cziEntry(dv,u8,o); ents.push(e); o+=e.len; } }
  else { // directory 가 없거나 잘못됨 → segment 를 처음부터 훑는다 (알림)
    warn.push("CZI: subblock directory 가 없어 파일 전체를 훑어 subblock 을 찾았습니다"); let o=0;
    while(o+32<=u8.length){ const id=String.fromCharCode(...u8.subarray(o,o+16)).replace(/\0.*$/,""); const alloc=i64(dv,o+16); if(id==="ZISRAWSUBBLOCK"){ const e=cziEntry(dv,u8,o+32+16); e.pos=o; ents.push(e); } if(!(alloc>0)) break; o+=32+alloc; } }
  if(!ents.length) throw new Error("CZI: subblock 이 없습니다");
  let xml=""; if(metaPos>0&&String.fromCharCode(...u8.subarray(metaPos,metaPos+14))==="ZISRAWMETADATA"){ const xs=dv.getInt32(metaPos+32,true); xml=tdu.decode(u8.subarray(metaPos+32+256,metaPos+32+256+xs)); }
  const meta=cziMeta(xml,info);
  // 층 0 만 (pyramid 층 = 저장 크기 ≠ 논리 크기)
  const L0=ents.filter(e=>e.pyr===0&&e.dims.X&&e.dims.Y&&e.dims.X.stored===e.dims.X.size&&e.dims.Y.stored===e.dims.Y.size);
  if(L0.length<ents.length) info.push(["Pyramid",`축소 층 subblock ${ents.length-L0.length} 개는 건너뜀`]);
  const bad=[...new Set(L0.map(e=>e.pt))].filter(t=>!CZI_PT[t]); if(bad.length) throw new Error("CZI: 화소 형식 "+bad.map(t=>CZI_PTN[t]||t).join(", ")+" 은 읽지 않습니다");
  const comps=[...new Set(L0.map(e=>e.comp))]; info.push(["Pixel",[...new Set(L0.map(e=>CZI_PTN[e.pt]))].join(", ")+" · "+comps.map(c=>CZI_COMP[c]||c).join(", ")]);
  // 차원 값 → 색인 (start 값 정렬)
  const val=(e,d)=>e.dims[d]?e.dims[d].start:0;
  const DIMS=[...new Set(L0.flatMap(e=>Object.keys(e.dims)))].filter(d=>!["X","Y","M"].includes(d));
  const vals={}; for(const d of DIMS){ vals[d]=[...new Set(L0.map(e=>val(e,d)))].sort((a,b)=>a-b); }
  const n=d=>vals[d]?vals[d].length:1; const idx=(e,d)=>vals[d]?vals[d].indexOf(val(e,d)):0;
  const extra=DIMS.filter(d=>!["C","Z","T","S"].includes(d)&&n(d)>1); // H·R·I·V·B 등 → 따로 영상
  const nc=n("C"), nz=n("Z"), nt=n("T"), ns=n("S");
  // 채널마다 화소 형식: RGB(Bgr*) 채널 = 채널마다 따로 RGB 영상, 회색 채널 = 하나의 다채널 영상 (밝은 시야 RGB + 형광 회색이 섞인 파일)
  const chPt=[...Array(nc)].map(()=>new Set()); for(const e of L0) chPt[idx(e,"C")].add(e.pt);
  const isRgbPt=t=>CZI_PT[t][1]>=3, isFPt=t=>[2,8,12,13].includes(t);
  const chKind=chPt.map(S=>{ const L=[...S]; if(!L.length) return {rgb:false,bit:8}; const rgb=L.some(isRgbPt); if(rgb&&!L.every(isRgbPt)) throw new Error("CZI: 한 채널 안에 RGB 와 회색 subblock 이 섞여 있습니다"); return {rgb,bit:L.some(isFPt)?32:Math.max(...L.map(t=>CZI_PT[t][0]))}; });
  const grayC=[...Array(nc).keys()].filter(c=>!chKind[c].rgb), rgbC=[...Array(nc).keys()].filter(c=>chKind[c].rgb);
  const grayBit=grayC.length?(grayC.some(c=>chKind[c].bit===32)?32:Math.max(...grayC.map(c=>chKind[c].bit))):8;
  const pts=[...new Set(L0.map(e=>e.pt))];
  const gPts=[...new Set(grayC.flatMap(c=>[...chPt[c]]))]; if(gPts.length>1) warn.push(`CZI: 채널마다 화소 형식이 다릅니다(${gPts.map(t=>CZI_PTN[t]).join(", ")}) — ${grayBit}-bit 로 맞췄습니다(값은 그대로)`);
  if(grayC.length&&rgbC.length) warn.push(`CZI: RGB 채널(${rgbC.length})과 회색 채널(${grayC.length})이 섞여 있어 따로 엽니다`);
  if(pts.includes(13)) warn.push("CZI: Gray64(double) 를 32-bit float 로 읽었습니다(유효숫자 약 7 자리)");
  if(pts.includes(12)) warn.push("CZI: Gray32(정수) 를 32-bit float 로 읽었습니다(2²⁴ 넘는 값은 반올림)");
  // 영상 묶음 = scene × extra
  const groups=new Map(); for(const e of L0){ const key=[idx(e,"S"),...extra.map(d=>idx(e,d))].join(","); if(!groups.has(key)) groups.set(key,[]); groups.get(key).push(e); }
  const cal=meta.sc.X?calOpt(meta.sc.X,meta.sc.Y||meta.sc.X,"CZI Scaling (Distance X·Y)"):null; const zStep=meta.sc.Z?meta.sc.Z*1e6:0;
  info.push(["Dimensions",`${DIMS.filter(d=>n(d)>1).map(d=>d+n(d)).join(" ")}${L0.some(e=>e.dims.M)?` M${new Set(L0.map(e=>e.dims.M?e.dims.M.start:0)).size}`:""}`]);
  const chNames=[...Array(nc)].map((_,c)=>meta.ch[c]&&meta.ch[c].name||`C${c+1}`);
  const chLut=[...Array(nc)].map((_,c)=>{ const q=meta.ch[c]; if(!q||!q.color) return null; const h=q.color.length===8?q.color.slice(2):q.color; return nearestLut(parseInt(h.slice(0,2),16),parseInt(h.slice(2,4),16),parseInt(h.slice(4,6),16)); });
  const keys=[...groups.keys()].sort((a,b)=>{ const A=a.split(",").map(Number), B=b.split(",").map(Number); for(let i=0;i<A.length;i++) if(A[i]!==B[i]) return A[i]-B[i]; return 0; });
  const bbox=list=>{ let x0=Infinity,y0=Infinity,x1=-Infinity,y1=-Infinity; for(const e of list){ x0=Math.min(x0,e.dims.X.start); y0=Math.min(y0,e.dims.Y.start); x1=Math.max(x1,e.dims.X.start+e.dims.X.size); y1=Math.max(y1,e.dims.Y.start+e.dims.Y.size); } return {x0,y0,W:x1-x0,H:y1-y0}; };
  const planesPer=grayC.length+3*rgbC.length; const posBytes=keys.map(k=>{ const b=bbox(groups.get(k)); return b.W*b.H*4*planesPer*nz*nt; });
  const all=k=>[...Array(k).keys()]; let sel={p:all(keys.length),z:all(nz),t:all(nt)};
  const keyName=k=>{ const q=k.split(",").map(Number); return [ns>1?`S${q[0]+1}`:"",...extra.map((d,i)=>`${d}${q[1+i]+1}`)].filter(Boolean).join(" ")||"영상"; };
  if(ctx&&ctx.choose){ const b0=bbox(groups.get(keys[0]));
    const r=await ctx.choose({kind:"czi",name,positions:keys.length,posNames:keys.map(keyName),posBytes,nz,nt,nc,w:b0.W,h:b0.H,bytesPerPlane:b0.W*b0.H*4,planesPer}); if(!r) return {images:[],warn:[],cancel:true}; sel=Object.assign(sel,r); }
  const pick=sel.p.map(i=>keys[i]); const zMap=new Map(sel.z.map((z,i)=>[z,i])), tMap=new Map(sel.t.map((t,i)=>[t,i])); const nzo=sel.z.length, nto=sel.t.length;
  { const tot=sel.p.reduce((a,i)=>a+posBytes[i],0)*nzo*nto/(nz*nt); const lim=(ctx&&ctx.maxBytes)||1.6e9; if(tot>lim) throw new Error(`CZI: 너무 큼 — 고른 ${sel.p.length} 개 ${(tot/1e9).toFixed(2)} GB (32-bit 메모리, 한도 ${(lim/1e9).toFixed(1)} GB). scene·Z·T 일부만 고르세요`); }
  if(nzo<nz) warn.push(`CZI: Z ${nz} 장 중 ${nzo} 장만 읽었습니다`); if(nto<nt) warn.push(`CZI: T ${nt} 장 중 ${nto} 장만 읽었습니다`);
  const zUni=nzo>1&&sel.z.every((z,i)=>i===0||z-sel.z[i-1]===sel.z[1]-sel.z[0]); const zStepOut=zStep>0?(nzo>1?(zUni?zStep*(sel.z[1]-sel.z[0]):0):zStep):0;
  const gIdx=new Map(grayC.map((c,i)=>[c,i])), ng=grayC.length;
  const images=[]; let done=0;
  for(const key of pick){ const list=groups.get(key); const {x0,y0,W,H}=bbox(list);
    const Pg=[...Array(ng*nzo*nto)].map(()=>new Float32Array(W*H)); const Pr=new Map(rgbC.map(c=>[c,[...Array(3*nzo*nto)].map(()=>new Float32Array(W*H))]));
    const ord=list.map((e,i)=>({e,i,m:e.dims.M?e.dims.M.start:-1})).sort((a,b)=>a.m-b.m||a.i-b.i); // libCZI: M 순서로 그림, 뒤가 위
    for(const {e} of ord){ const z0=idx(e,"Z"), t0=idx(e,"T"); if(!zMap.has(z0)||!tMap.has(t0)) continue; const s=await cziPixels(u8,e,ctx); const c=idx(e,"C"), z=zMap.get(z0), t=tMap.get(t0); const ox=e.dims.X.start-x0, oy=e.dims.Y.start-y0;
      for(let k=0;k<s.comps.length;k++){ const dst=chKind[c].rgb?Pr.get(c)[k+3*(z+nzo*t)]:Pg[gIdx.get(c)+ng*(z+nzo*t)]; const src=s.comps[k];
        for(let y=0;y<s.H;y++) dst.set(src.subarray(y*s.W,(y+1)*s.W),(oy+y)*W+ox); }
      if(ctx&&ctx.onProgress&&((++done)&7)===0) ctx.onProgress(done/L0.length); await null; }
    const sIdx=+key.split(",")[0]; const exN=extra.map((d,i)=>`${d}${+key.split(",")[1+i]+1}`).join(" ");
    const tag=[ns>1||pick.length>1?`S${sIdx+1}`:"",exN].filter(Boolean).join(" ");
    const base=`${name}${tag?` · ${tag}`:""}`; const common={w:W,h:H,nz:nzo,nt:nto,order:"czt",cal,calCandidates:[],info:info.slice(),warn:[],zStep:zStepOut};
    if(ng) images.push({...common,name:rgbC.length?`${base} · ${grayC.map(c=>chNames[c]).join("·")}`:base,bit:grayBit,rgb:false,nc:ng,planes:Pg,chNames:grayC.map(c=>chNames[c]),chLut:grayC.map(c=>chLut[c])});
    for(const c of rgbC){ const P=Pr.get(c); const nm=nc>1?`${base} · ${chNames[c]}`:base;
      if(chKind[c].bit===8) images.push({...common,name:nm,bit:8,rgb:true,nc:3,planes:P});
      else images.push({...common,name:nm,bit:chKind[c].bit,rgb:false,nc:3,planes:P,chNames:["R","G","B"],chLut:["Red","Green","Blue"]}); }
  }
  if(images.length) images[0].metaText=xml.length>200000?xml.slice(0,200000)+"…":xml; // 원문은 첫 영상에만 (저장 공간)
  return {images,warn};
}
IO.readCZI=readCZI;
})();
