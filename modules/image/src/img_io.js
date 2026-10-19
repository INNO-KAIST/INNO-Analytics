/* ============================================================
   INNO Analytics — Image 탭 파일 읽기·쓰기 (IMG_IO)
   읽기: TIFF/BigTIFF(8·16·32-bit, float, RGB, 다중 쪽, 타일, LZW·Deflate·PackBits, predictor 2·3,
         ImageJ hyperstack·OME-TIFF 메타데이터, Zeiss·FEI/Thermo·TESCAN SEM 태그, Hitachi·JEOL 옆 파일(.txt)),
         PNG(1–16-bit, 직접 해독), JPEG·BMP·GIF·WebP(브라우저), DM3/DM4·MRC·EMD/HDF5·BCF(앱의 기존 파서 재사용),
         Igor IBW(Asylum AFM), Nanonis SXM, NumPy NPY, 텍스트 행렬, Raw(사용자 지정)
   쓰기: ImageJ 호환 TIFF(보정·스택·hyperstack 메타데이터 포함)
   원칙: 화소 크기는 출처가 확실한 값만 적용한다. 출처가 불확실한 값(TIFF 해상도 태그 등)은 '후보'로만 보인다.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IO={};
const td=new TextDecoder("latin1"), tdu=new TextDecoder("utf-8"), tdf=new TextDecoder("utf-8",{fatal:true});
/* 문자열 태그: UTF-8 로 읽히면 UTF-8(µ 등), 아니면 latin1 */
const txt=b=>{ try{ return tdf.decode(b); }catch(e){ return td.decode(b); } }; IO._txt=txt;

/* ---------- 압축 풀기 ---------- */
IO.MAXBYTES=1.6e9; // 32-bit 메모리 한도 (nd2·czi 와 같음)
IO.inflate=async function(u8,raw){ // zlib(deflate) — 브라우저 DecompressionStream, node 시험은 IO.inflate 를 바꿔 끼운다
  if(typeof DecompressionStream==="undefined") throw new Error("이 브라우저는 deflate 압축을 풀 수 없습니다");
  const rs=new Blob([u8]).stream().pipeThrough(new DecompressionStream(raw?"deflate-raw":"deflate"));
  return new Uint8Array(await new Response(rs).arrayBuffer());
};
function lzw(src,expected){ // TIFF LZW (MSB 우선, early change)
  let out=new Uint8Array(Math.max(expected||0,src.length*3)), op=0;
  const grow=n=>{ if(op+n>out.length){ const o=new Uint8Array(Math.max(out.length*2,op+n)); o.set(out); out=o; } };
  const pre=new Int32Array(4096), suf=new Uint8Array(4096), len=new Int32Array(4096), first=new Uint8Array(4096);
  for(let i=0;i<256;i++){ pre[i]=-1; suf[i]=i; len[i]=1; first[i]=i; }
  let next=258, bits=9, old=-1, bp=0; const nbits=src.length*8;
  const read=()=>{ if(bp+bits>nbits) return 257; let v=0; for(let k=0;k<bits;k++){ const b=bp+k; v=(v<<1)|((src[b>>3]>>(7-(b&7)))&1); } bp+=bits; return v; };
  const emit=c=>{ const L=len[c]; grow(L); let p=op+L-1, k=c; while(k>=0&&p>=op){ out[p--]=suf[k]; k=pre[k]; } op+=L; };
  for(;;){ const code=read(); if(code===257) break;
    if(code===256){ next=258; bits=9; old=-1; continue; }
    if(code>next||(old===-1&&code>255)) throw new Error("LZW 자료가 손상되었습니다");
    if(old===-1){ emit(code); old=code; continue; }
    if(code<next){ emit(code); if(next<4096){ pre[next]=old; suf[next]=first[code]; len[next]=len[old]+1; first[next]=first[old]; next++; } }
    else { // code === next (KwKwK)
      if(next<4096){ pre[next]=old; suf[next]=first[old]; len[next]=len[old]+1; first[next]=first[old]; next++; } emit(code); }
    old=code; if(next+1>=(1<<bits)&&bits<12) bits++; }
  return out.subarray(0,op);
}
function packbits(src,expected){ const out=new Uint8Array(expected); let i=0,o=0; while(i<src.length&&o<expected){ const n=(src[i]<<24>>24); i++; if(n>=0){ for(let k=0;k<=n&&o<expected;k++) out[o++]=src[i++]; } else if(n!==-128){ const b=src[i++]; for(let k=0;k<1-n&&o<expected;k++) out[o++]=b; } } return out; }
IO._lzw=lzw; IO._packbits=packbits;

/* ============================================================ TIFF ============================================================ */
const TSZ={1:1,2:1,3:2,4:4,5:8,6:1,7:1,8:2,9:4,10:8,11:4,12:8,13:4,16:8,17:8,18:8};
function tiffIFDs(u8){
  const le=u8[0]===0x49; const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength);
  const magic=dv.getUint16(2,le); const big=magic===43; if(magic!==42&&!big) throw new Error("TIFF 머리가 올바르지 않습니다");
  const u64=o=>{ const lo=dv.getUint32(o,le), hi=dv.getUint32(o+4,le); return le?hi*4294967296+lo:lo*4294967296+hi; };
  const hiLo=o=>{ const a=dv.getUint32(o,le), b=dv.getUint32(o+4,le); return le?b*4294967296+a:a*4294967296+b; };
  let ifd=big?hiLo(8):dv.getUint32(4,le); const pages=[]; const seen=new Set();
  while(ifd>0&&ifd<u8.byteLength&&!seen.has(ifd)){
    seen.add(ifd);
    const n=big?hiLo(ifd):dv.getUint16(ifd,le); const es=big?20:12, base=ifd+(big?8:2); const tags={};
    for(let i=0;i<n;i++){ const o=base+i*es; const code=dv.getUint16(o,le), type=dv.getUint16(o+2,le); const cnt=big?hiLo(o+4):dv.getUint32(o+4,le);
      const sz=(TSZ[type]||1)*cnt; const inl=big?8:4; let vo=o+(big?12:8); if(sz>inl) vo=big?hiLo(o+12):dv.getUint32(o+8,le);
      if(vo+sz>u8.byteLength){ continue; }
      const rd=k=>{ const q=vo+k*(TSZ[type]||1); switch(type){ case 1: case 7: case 2: return dv.getUint8(q); case 3: return dv.getUint16(q,le); case 4: case 13: return dv.getUint32(q,le);
        case 5: return dv.getUint32(q,le)/(dv.getUint32(q+4,le)||1); case 6: return dv.getInt8(q); case 8: return dv.getInt16(q,le); case 9: return dv.getInt32(q,le);
        case 10: return dv.getInt32(q,le)/(dv.getInt32(q+4,le)||1); case 11: return dv.getFloat32(q,le); case 12: return dv.getFloat64(q,le); case 16: case 18: return hiLo(q); case 17: return hiLo(q); default: return dv.getUint8(q); } };
      const t={type,cnt,off:vo,len:sz};
      if(type===2) t.str=txt(u8.subarray(vo,vo+sz)).replace(/\0+$/,"");
      else if(cnt<=1) t.v=rd(0);
      else if(cnt<=1e6&&type!==7&&type!==1) t.a=Array.from({length:cnt},(_,k)=>rd(k));
      t.raw=u8.subarray(vo,vo+sz);
      tags[code]=t; }
    pages.push(tags);
    const nx=base+n*es; ifd=big?hiLo(nx):dv.getUint32(nx,le);
    if(pages.length>100000) break;
  }
  void u64;
  return {le,big,dv,pages};
}
const tv=(t,d)=>t?(t.v!=null?t.v:(t.a?t.a[0]:d)):d;
const ta=(t,d)=>t?(t.a?t.a:(t.v!=null?[t.v]:d)):d;
async function tiffDecode(T,u8,pi,limitPx){
  const tg=T.pages[pi], le=T.le;
  const W=tv(tg[256]), H=tv(tg[257]); if(!(W>0&&H>0)) throw new Error("쪽 "+(pi+1)+": 크기 태그 없음");
  const bpsA=ta(tg[258],[1]), spp=tv(tg[277],1), bits=bpsA[0], fmt=tv(tg[339],1), comp=tv(tg[259],1), photo=tv(tg[262],spp>=3?2:1), pred=tv(tg[317],1), planar=tv(tg[284],1);
  if(bpsA.some(b=>b!==bits)) throw new Error("채널마다 비트 수가 다른 TIFF 는 읽지 않습니다");
  if(![1,2,4,8,16,32,64].includes(bits)) throw new Error(bits+"-bit TIFF 는 읽지 않습니다");
  if(comp===7||comp===6) throw new Error("JPEG 압축 TIFF 는 아직 읽지 않습니다 (무손실 TIFF 로 다시 저장하세요)");
  if(![1,5,8,32946,32773].includes(comp)) throw new Error("압축 방식 "+comp+" 은 지원하지 않습니다");
  if(limitPx&&W*H>limitPx) throw new Error(`이미지가 너무 큽니다 (${W}×${H})`);
  if(W*H*Math.min(spp,4)*4>IO.MAXBYTES) throw new Error(`이미지가 너무 큽니다: ${W}×${H}×${spp} = ${(W*H*Math.min(spp,4)*4/1e9).toFixed(2)} GB (32-bit 메모리, 한도 ${(IO.MAXBYTES/1e9).toFixed(1)} GB)`);
  const tiled=!!tg[322]; const tw=tiled?tv(tg[322]):W, th=tiled?tv(tg[323]):Math.min(H,tv(tg[278],H)); // RowsPerStrip 기본값 2³²−1
  const offs=ta(tiled?tg[324]:tg[273]), cnts=ta(tiled?tg[325]:tg[279]);
  if(!offs) throw new Error("화소 자료 위치(StripOffsets) 태그가 없습니다");
  const sppChunk=planar===2?1:spp; const bytesPS=bits/8;
  const rowBytes=Math.ceil(tw*sppChunk*bits/8);
  const chunkBytes=rowBytes*th;
  const nPlanes=planar===2?spp:1;
  const out=[...Array(spp)].map(()=>new Float32Array(W*H));
  const across=Math.ceil(W/tw), down=Math.ceil(H/th), perPlane=across*down;
  const dv0=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); let short=0;
  for(let pl=0;pl<nPlanes;pl++){
    for(let ci=0;ci<perPlane;ci++){
      const k=pl*perPlane+ci; if(k>=offs.length) break;
      const nb=cnts?cnts[k]:(comp===1?chunkBytes:((offs[k+1]||u8.length)-offs[k])); // StripByteCounts 가 없으면 다음 strip 까지
      if(offs[k]>=u8.length||offs[k]+nb>u8.length) short++;
      let src=u8.subarray(offs[k],offs[k]+nb);
      let buf;
      if(comp===1) buf=src; else if(comp===5) buf=lzw(src,chunkBytes); else if(comp===32773) buf=packbits(src,chunkBytes); else buf=await IO.inflate(src);
      const y0c=Math.floor((ci%perPlane)/across)*th, need=rowBytes*Math.min(th,H-y0c); if(buf.length<(tiled?chunkBytes:need)) short++; // 마지막 strip 은 짧아도 됨
      if(buf.length<chunkBytes){ const b=new Uint8Array(chunkBytes); b.set(buf); buf=b; } else if(buf===src) buf=Uint8Array.from(src.subarray(0,chunkBytes));
      const dv=new DataView(buf.buffer,buf.byteOffset,buf.byteLength);
      // predictor
      if(pred===2){ for(let y=0;y<th;y++){ const o=y*rowBytes; const ns=tw*sppChunk;
          if(bits===8) for(let i=sppChunk;i<ns;i++) buf[o+i]=(buf[o+i]+buf[o+i-sppChunk])&255;
          else if(bits===16) for(let i=sppChunk;i<ns;i++) dv.setUint16(o+2*i,(dv.getUint16(o+2*i,le)+dv.getUint16(o+2*(i-sppChunk),le))&0xffff,le);
          else if(bits===32) for(let i=sppChunk;i<ns;i++) dv.setUint32(o+4*i,(dv.getUint32(o+4*i,le)+dv.getUint32(o+4*(i-sppChunk),le))>>>0,le); } }
      else if(pred===3){ // floating point predictor (Adobe TN3)
        const bps=bits/8, ns=tw*sppChunk, tmp=new Uint8Array(ns*bps);
        for(let y=0;y<th;y++){ const o=y*rowBytes; for(let i=sppChunk;i<ns*bps;i++) buf[o+i]=(buf[o+i]+buf[o+i-sppChunk])&255; tmp.set(buf.subarray(o,o+ns*bps));
          for(let i=0;i<ns;i++) for(let b=0;b<bps;b++){ const v=tmp[b*ns+i]; buf[o+i*bps+(le?bps-1-b:b)]=v; } } }
      const x0=(ci%across)*tw, y0=Math.floor(ci/across)*th;
      for(let y=0;y<th;y++){ const Y=y0+y; if(Y>=H) break; const o=y*rowBytes;
        for(let x=0;x<tw;x++){ const X=x0+x; if(X>=W) break;
          for(let s=0;s<sppChunk;s++){ const si=x*sppChunk+s; let v;
            if(bits===8) v=fmt===2?(buf[o+si]<<24>>24):buf[o+si];
            else if(bits===16) v=fmt===2?dv.getInt16(o+2*si,le):fmt===3?f16(dv.getUint16(o+2*si,le)):dv.getUint16(o+2*si,le);
            else if(bits===32) v=fmt===3?dv.getFloat32(o+4*si,le):fmt===2?dv.getInt32(o+4*si,le):dv.getUint32(o+4*si,le);
            else if(bits===64) v=fmt===3?dv.getFloat64(o+8*si,le):fmt===2?Number(dv.getBigInt64(o+8*si,le)):Number(dv.getBigUint64(o+8*si,le));
            else { const bitPos=si*bits; const byte=buf[o+(bitPos>>3)]; v=(byte>>(8-bits-(bitPos&7)))&((1<<bits)-1); }
            out[planar===2?pl:s][Y*W+X]=v; } } }
    }
  }
  void dv0;
  if(short) throw new Error(`쪽 ${pi+1}: 화소 자료가 잘렸거나 손상되었습니다 (strip/tile ${short} 개가 짧음) — 파일이 끝까지 저장·복사되었는지 확인하세요`);
  // 형식
  let bit=bits<=8?8:bits===16&&fmt!==2&&fmt!==3?16:32;
  if(bits===16&&fmt===2){ bit=32; }
  if(bits<8&&photo!==3){ const mx=(1<<bits)-1; for(const p of out) for(let i=0;i<p.length;i++) p[i]=Math.round(p[i]*255/mx); } // 팔레트 색인은 그대로
  if(bits===8&&fmt===2){ bit=32; }
  const notes=[];
  if(photo===0&&spp===1){ notes.push("WhiteIsZero(0=흰색) TIFF — 화소 값은 그대로 두고 LUT 를 뒤집어 보입니다"); }
  let rgb=false, palette=null, photo0inv=false;
  if(photo===3&&tg[320]){ const cm=ta(tg[320]); const n=cm.length/3; let gray=true; for(let i=0;i<n;i++) if(cm[i]!==cm[n+i]||cm[i]!==cm[2*n+i]){ gray=false; break; }
    if(!gray){ const r=new Float32Array(W*H),g=new Float32Array(W*H),b=new Float32Array(W*H); const s=out[0]; for(let i=0;i<W*H;i++){ const k=Math.min(n-1,s[i]|0); r[i]=cm[k]>>8; g[i]=cm[n+k]>>8; b[i]=cm[2*n+k]>>8; } return {W,H,bit:8,planes:[r,g,b],rgb:true,inv:false,notes:["색 팔레트 TIFF 를 RGB 로 바꿨습니다"],bits,fmt,comp}; }
    palette="gray"; if(n>1&&cm[0]>cm[n-1]) { photo0inv=true; } }
  if(photo===2&&spp>=3&&bit===8){ rgb=true; }
  const planes=rgb?out.slice(0,3):out.slice(0,photo===2?Math.min(spp,3):spp);
  if(spp===4&&photo===2) notes.push("알파 채널은 버렸습니다");
  if(photo===5) notes.push("CMYK TIFF — 4 채널(C·M·Y·K)로 열었습니다 (RGB 로 바꾸지 않음)"); if(photo===6) notes.push("YCbCr TIFF — 3 채널(Y·Cb·Cr)로 열었습니다 (RGB 로 바꾸지 않음)");
  return {W,H,bit,planes,rgb,inv:photo===0||photo0inv,notes,bits,fmt,comp,palette};
}
function f16(h){ const s=(h&0x8000)?-1:1, e=(h>>10)&31, m=h&1023; if(e===0) return s*Math.pow(2,-14)*(m/1024); if(e===31) return m?NaN:s*Infinity; return s*Math.pow(2,e-15)*(1+m/1024); }

/* ---------- 메타데이터 ---------- */
const UNIT={ "m":1,"meter":1,"mm":1e-3,"millimeter":1e-3,"cm":1e-2,"um":1e-6,"µm":1e-6,"μm":1e-6,"micron":1e-6,"microns":1e-6,"nm":1e-9,"nanometer":1e-9,"pm":1e-12,"Å":1e-10,"A":1e-10,"angstrom":1e-10,"inch":0.0254,"in":0.0254 };
IO.unitToM=u=>{ if(u==null) return null; const s=String(u).trim().replace(/^\\u00B5/,"µ"); if(UNIT[s]!=null) return UNIT[s]; if(UNIT[s.toLowerCase()]!=null) return UNIT[s.toLowerCase()]; if(/^\\?mu?m$/i.test(s)||/^micro/i.test(s)) return 1e-6; return null; };
/* 미터 단위 화소 크기 → 사람이 읽기 좋은 단위 */
IO.niceUnit=function(m){ if(!(m>0)) return {unit:"px",f:1}; if(m>=1e-4) return {unit:"mm",f:1e3}; if(m>=1e-7) return {unit:"µm",f:1e6}; if(m>=1e-11) return {unit:"nm",f:1e9}; return {unit:"Å",f:1e10}; }; // pm 은 쓰지 않음(사용자 요청 10-19: nm·Å 가 편함) — 0.01 nm 미만만 Å. 화소 값이 0.1 이상이 되는 가장 큰 단위 — TEM 0.55 nm/px 는 nm (예전엔 552 pm 이라 Size 칸에 nm² 로 넣은 값이 10⁶ 배 틀림, 10-19). pm 은 0.01 nm 미만에서만
IO.calFromMeters=function(pwm,phm,src,unit){ const k=unit&&IO.unitToM(unit); const nu=k&&k!==0.0254&&k!==1e-2&&k!==1?{unit:k===1e-6?"µm":k===1e-9?"nm":k===1e-12?null:k===1e-10?"Å":k===1e-3?"mm":null,f:1/k}:null; const u=nu&&nu.unit?nu:IO.niceUnit(pwm); return {pw:pwm*u.f,ph:(phm||pwm)*u.f,unit:u.unit,src}; }; // unit: 파일이 밝힌 단위가 있으면 그 단위로 (micron → µm)
function parseImageJDesc(s){ const o={}; s.split(/\r?\n/).forEach(l=>{ const i=l.indexOf("="); if(i>0) o[l.slice(0,i).trim()]=l.slice(i+1).trim(); }); return o; }
function parseOME(s){ const L=parseOMEAll(s); return L.length?L[0]:null; }
function parseOMEAll(s){ // <Image><Pixels .../> 마다 하나 (multi-series OME-TIFF)
  const out=[]; const re=/<(?:\w+:)?Pixels\b[^>]*>/g; let m;
  while((m=re.exec(s))){ const t=m[0]; const a=k=>{ const q=t.match(new RegExp("\\b"+k+'="([^"]*)"')); return q?q[1]:null; };
    out.push({sizeX:+a("SizeX"),sizeY:+a("SizeY"),sizeC:+a("SizeC")||1,sizeZ:+a("SizeZ")||1,sizeT:+a("SizeT")||1,order:a("DimensionOrder")||"XYCZT",type:a("Type"),
      psx:a("PhysicalSizeX")?+a("PhysicalSizeX"):null,psy:a("PhysicalSizeY")?+a("PhysicalSizeY"):null,psxu:a("PhysicalSizeXUnit")||"µm",psyu:a("PhysicalSizeYUnit")||"µm",psz:a("PhysicalSizeZ")?+a("PhysicalSizeZ"):null,pszu:a("PhysicalSizeZUnit")||"µm"}); }
  return out; }
function iniSections(s){ const o={}; let cur=""; s.split(/\r?\n/).forEach(l=>{ const m=l.match(/^\s*\[(.+)\]\s*$/); if(m){ cur=m[1]; o[cur]=o[cur]||{}; return; } const i=l.indexOf("="); if(i>0){ (o[cur]=o[cur]||{})[l.slice(0,i).trim()]=l.slice(i+1).trim(); } }); return o; }
/* SEM·TEM 장비 태그에서 화소 크기 */
function semCal(tg,info,warn){
  // FEI / Thermo Fisher (34680 FEI_SFEG, 34682 FEI_HELIOS, 34683)
  for(const code of [34682,34680,34683]){ const t=tg[code]; if(!t) continue; const s=t.str!=null?t.str:txt(t.raw); const ini=iniSections(s);
    const pw=(ini.Scan&&+ini.Scan.PixelWidth)||(ini.EScan&&+ini.EScan.PixelWidth)||(ini.IScan&&+ini.IScan.PixelWidth);
    const ph=(ini.Scan&&+ini.Scan.PixelHeight)||(ini.EScan&&+ini.EScan.PixelHeight)||pw;
    info.push(["Instrument",(ini.System&&(ini.System.SystemType||ini.System.Type))||"FEI/Thermo"]);
    if(ini.Beam&&ini.Beam.HV) info.push(["HV (V)",ini.Beam.HV]); if(ini.EBeam&&ini.EBeam.HV) info.push(["HV (V)",ini.EBeam.HV]);
    if(ini.Image&&ini.Image.MagnificationMode) info.push(["Mag mode",ini.Image.MagnificationMode]);
    const ry=ini.Image&&+ini.Image.ResolutionY; const dbh=ini.PrivateFei&&+ini.PrivateFei.DatabarHeight;
    const r={raw:s}; if(pw>0) r.cal=IO.calFromMeters(pw,ph,"FEI/Thermo 태그 "+code+" [Scan] PixelWidth"); if(ry>0) r.imageRows=ry; if(dbh>0) r.databar=dbh; return r; }
  // Zeiss (34118 CZ_SEM)
  if(tg[34118]){ const s=txt(tg[34118].raw); const m=s.match(/Image Pixel Size\s*=\s*([\d.eE+-]+)\s*(pm|nm|µm|um|mm)/)||s.match(/Pixel Size\s*=\s*([\d.eE+-]+)\s*(pm|nm|µm|um|mm)/);
    const hv=s.match(/EHT\s*=\s*([\d.]+\s*\w+)/); if(hv) info.push(["EHT",hv[1]]); const mag=s.match(/\bMag\s*=\s*([\d.,]+\s*\w+)/); if(mag) info.push(["Mag",mag[1]]);
    info.push(["Instrument","Zeiss SEM"]);
    if(m){ const u=m[2]==="um"?"µm":m[2]; const k=IO.unitToM(u); return {cal:IO.calFromMeters(+m[1]*k,+m[1]*k,"Zeiss 태그 34118 'Image Pixel Size'"),raw:s}; }
    warn.push("Zeiss 태그에서 Pixel Size 를 찾지 못했습니다"); return {raw:s}; }
  // TESCAN (50431)
  if(tg[50431]){ const s=txt(tg[50431].raw); const m=s.match(/PixelSizeX\s*=\s*([\d.eE+-]+)/), n=s.match(/PixelSizeY\s*=\s*([\d.eE+-]+)/); info.push(["Instrument","TESCAN"]);
    if(m) return {cal:IO.calFromMeters(+m[1],n?+n[1]:+m[1],"TESCAN 태그 50431 PixelSizeX"),raw:s}; }
  return null;
}
function sidecarCal(name,ctx,info){
  const sc=ctx&&ctx.sidecars; if(!sc) return null; const base=name.replace(/\.[^.]+$/,"").toLowerCase(); const txt=sc[base+".txt"]||sc[base+".txt".toUpperCase()]; if(!txt) return null;
  // Hitachi: PixelSize=…(nm)
  let m=txt.match(/^\s*PixelSize\s*=\s*([\d.eE+-]+)/mi); if(m){ info.push(["Instrument",(txt.match(/InstructName\s*=\s*(.+)/i)||[])[1]||"Hitachi SEM"]); const mg=txt.match(/^\s*Magnification\s*=\s*(\S+)/mi); if(mg) info.push(["Mag",mg[1]]);
    return {cal:IO.calFromMeters(+m[1]*1e-9,+m[1]*1e-9,"Hitachi 옆 파일 "+base+".txt 'PixelSize' (nm)")}; }
  // JEOL: $$SM_MICRON_BAR (px), $$SM_MICRON_MARKER (길이)
  const bar=txt.match(/\$\$?SM_MICRON_BAR\s+([\d.]+)/), mk=txt.match(/\$\$?SM_MICRON_MARKER\s+([\d.]+)\s*(nm|um|µm|mm)/i);
  if(bar&&mk){ info.push(["Instrument","JEOL SEM"]); const L=+mk[1]*IO.unitToM(mk[2].toLowerCase()==="um"?"µm":mk[2]); const px=L/(+bar[1]); return {cal:IO.calFromMeters(px,px,`JEOL 옆 파일 ${base}.txt (축척 막대 ${mk[1]} ${mk[2]} = ${bar[1]} px)`)}; }
  return null;
}

async function readTIFF(name,u8,ctx){
  const T=tiffIFDs(u8); if(!T.pages.length) throw new Error("TIFF: 읽을 수 있는 IFD 가 없습니다 (파일이 손상되었거나 잘림)"); const p0=T.pages[0]; const warn=[], info=[];
  const desc=p0[270]&&p0[270].str||"";
  const ij=desc.startsWith("ImageJ=")?parseImageJDesc(desc):null;
  const omeAll=/<OME[\s>]/.test(desc)?parseOMEAll(desc):[]; const ome=omeAll[0]||null;
  const sem=semCal(p0,info,warn);
  // 쪽 묶기: 같은 크기·형식인 쪽들 = 한 스택 (썸네일 등 다른 크기 쪽은 따로)
  const decoded=[]; const nP=T.pages.length;
  { let tot=0; for(let i=0;i<nP;i++){ const tg=T.pages[i]; if(i>0&&(tv(tg[254],0)&1)) continue; tot+=(tv(tg[256],0)*tv(tg[257],0))*Math.min(3,tv(tg[277],1))*4; } if(ij&&+ij.images>nP) tot=Math.max(tot,(tv(p0[256],0)*tv(p0[257],0))*Math.min(3,tv(p0[277],1))*4*(+ij.images));
    if(tot>IO.MAXBYTES) throw new Error(`TIFF: 너무 큼 — ${(tot/1e9).toFixed(2)} GB (32-bit 메모리, 한도 ${(IO.MAXBYTES/1e9).toFixed(1)} GB). ImageJ 에서 일부(substack)만 저장해 여세요`); }
  let nImg=ij&&+ij.images>1?+ij.images:nP;
  for(let i=0;i<nP;i++){ const tg=T.pages[i]; const sub=tv(tg[254],0); if(i>0&&(sub&1)) continue; // reduced-resolution 쪽 건너뜀
    try{ decoded.push(await tiffDecode(T,u8,i,ctx&&ctx.limitPx)); }catch(e){ if(i===0) throw e; warn.push(`쪽 ${i+1}: ${e.message}`); }
    if(ctx&&ctx.onProgress) ctx.onProgress(i/nP); }
  // ImageJ 큰 스택: IFD 가 하나뿐이고 images=N 이면 이어진 비압축 자료
  if(ij&&+ij.images>decoded.length&&decoded.length===1&&tv(p0[259],1)===1){
    const W=tv(p0[256]),H=tv(p0[257]), bits=ta(p0[258],[8])[0], fmt=tv(p0[339],1), spp=tv(p0[277],1), off=ta(p0[273])[0]; const per=W*H*spp*bits/8; const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const rgb0=decoded[0].rgb;
    const rd=o=>bits===8?(fmt===2?(u8[o]<<24>>24):u8[o]):bits===16?(fmt===2?dv.getInt16(o,T.le):dv.getUint16(o,T.le)):fmt===3?dv.getFloat32(o,T.le):fmt===2?dv.getInt32(o,T.le):dv.getUint32(o,T.le);
    for(let k=1;k<+ij.images;k++){ const o=off+k*per; if(o+per>u8.byteLength){ warn.push(`ImageJ 스택 ${ij.images} 장 중 ${k} 장만 파일 안에 있습니다`); break; }
      const ch=[...Array(rgb0?3:1)].map(()=>new Float32Array(W*H)); const bs=bits/8;
      for(let i=0;i<W*H;i++) for(let c=0;c<ch.length;c++) ch[c][i]=rd(o+(i*spp+c)*bs);
      decoded.push({W,H,bit:decoded[0].bit,planes:ch,rgb:rgb0,notes:[]}); }
    nImg=decoded.length; }
  void nImg;
  let groups=[];
  if(omeAll.length>1&&omeAll.reduce((a,o)=>a+o.sizeC*o.sizeZ*o.sizeT,0)===decoded.reduce((a,d)=>a+(d.rgb?1:d.planes.length),0)){ // multi-series OME-TIFF: series 마다 따로
    let k=0; omeAll.forEach((o,si)=>{ const n=o.sizeC*o.sizeZ*o.sizeT; const list=[]; let got=0; while(got<n&&k<decoded.length){ const d=decoded[k++]; list.push(d); got+=d.rgb?1:d.planes.length; } if(list.length) groups.push({W:list[0].W,H:list[0].H,bit:list[0].bit,rgb:list[0].rgb,np:list[0].planes.length,list,planes:list[0].planes,ome:o,series:si}); });
  } else for(const d of decoded){ const g=groups.find(g=>g.W===d.W&&g.H===d.H&&g.bit===d.bit&&g.rgb===d.rgb&&g.planes[0]&&g.np===d.planes.length); if(g){ g.list.push(d); } else groups.push({W:d.W,H:d.H,bit:d.bit,rgb:d.rgb,np:d.planes.length,list:[d],planes:d.planes}); }
  const images=[];
  groups.forEach((g,gi)=>{
    const planes=[]; g.list.forEach(d=>d.planes.forEach(p=>planes.push(p)));
    let nc=g.rgb?3:g.np, nz=g.list.length, nt=1, order="czt";
    if(gi===0&&ij&&!g.rgb){ const C=+ij.channels||1, Z=+ij.slices||1, Tt=+ij.frames||1; if(C*Z*Tt===planes.length){ nc=C; nz=Z; nt=Tt; } }
    if(gi===0&&ij&&g.rgb){ const Z=+ij.slices||1, Tt=+ij.frames||1; if(Z*Tt===g.list.length){ nz=Z; nt=Tt; } } /* RGB hyperstack (z·t) — ImageJ 와 같게 */
    const om=g.ome||(gi===0?ome:null);
    if(om&&!(gi===0&&ij&&!g.rgb)&&!g.rgb&&om.sizeC*om.sizeZ*om.sizeT===planes.length){ const ome=om; nc=ome.sizeC; nz=ome.sizeZ; nt=ome.sizeT;
      // OME DimensionOrder(XYZCT 등)의 파일 순서를 앱 순서 c + nc·(z + nz·t) 로 다시 놓는다
      const ord=ome.order.slice(2).toLowerCase(); const size={c:nc,z:nz,t:nt}; const re=new Array(planes.length);
      for(let k=0;k<planes.length;k++){ let r=k; const ix={}; for(const ch of ord){ ix[ch]=r%size[ch]; r=Math.floor(r/size[ch]); } re[ix.c+nc*(ix.z+nz*ix.t)]=planes[k]; }
      for(let k=0;k<planes.length;k++) planes[k]=re[k]; }
    let cal=null; const cands=[];
    if(gi===0&&sem&&sem.cal) cal=sem.cal;
    if(!cal&&gi===0&&ij&&ij.unit&&p0[282]){ const xr=tv(p0[282]), yr=tv(p0[283],xr); const u=ij.unit.replace(/\\u00B5/g,"µ"); const k=IO.unitToM(u);
      if(xr>0&&k){ cal=IO.calFromMeters(k/xr,k/yr,"ImageJ 메타데이터 (unit="+u+", XResolution)",u); } else if(xr>0&&u!=="pixel"){ cal={pw:1/xr,ph:1/yr,unit:u,src:"ImageJ 메타데이터"}; } }
    if(!cal&&om&&om.psx){ const k=IO.unitToM(om.psxu)||1e-6; cal=IO.calFromMeters(om.psx*k,(om.psy||om.psx)*(IO.unitToM(om.psyu)||k),"OME-XML PhysicalSizeX"+(omeAll.length>1?` (series ${(g.series||0)+1})`:""),om.psxu); }
    if(!cal&&gi===0){ const sc=sidecarCal(name,ctx,info); if(sc) cal=sc.cal; }
    if(!cal&&gi===0&&typeof G.tiffGatanScale==="function"&&p0[65027]){ try{ const gs=G.tiffGatanScale({65027:{raw:p0[65027].raw}}); if(gs&&gs.ok) cal=IO.calFromMeters(gs.nmPerPx*1e-9,gs.nmPerPx*1e-9,"Gatan 태그 65027 → "+gs.how); else if(gs&&gs.why) warn.push("Gatan TIFF: "+gs.why); }catch(e){} }
    if(!cal&&gi===0&&p0[282]&&!ij){ const xr=tv(p0[282]), ru=tv(p0[296],2); if(xr>0&&ru===3) cands.push({pw:1e4/xr,unit:"µm",src:`TIFF XResolution ${xr.toFixed(3)} /cm`}); else if(xr>0&&ru===2&&![72,96,150,300,600,1200].includes(Math.round(xr))) cands.push({pw:25400/xr,unit:"µm",src:`TIFF XResolution ${xr.toFixed(3)} /inch`}); }
    const im={name:groups.length>1?`${name} #${gi+1}`:name,w:g.W,h:g.H,bit:g.bit,rgb:g.rgb,nc,nz,nt,order,planes,cal,calCandidates:cands,info:info.slice(),warn:(g.list[0]&&g.list[0].notes||[]).slice(),invertLut:!!(g.list[0]&&g.list[0].inv)}; // 파일 전체 알림은 r.warn 으로 한 번만
    if(om&&om.psz>0&&nz>1){ const kz=IO.unitToM(om.pszu)||1e-6; im.zStep=om.psz*kz/1e-6; }
    // ImageJ 'signed 16-bit': 값 = 저장값 − 32768 (Calibration function 0, c0=−32768, c1=1)
    if(gi===0&&ij&&g.bit===16&&ij.cf==="0"&&+ij.c0===-32768&&+ij.c1===1){ for(const P of planes) for(let i=0;i<P.length;i++) P[i]-=32768; im.bit=32; im.info.push(["Signed 16-bit","ImageJ calibration c0=−32768 적용 → 32-bit 로 열었습니다"]); if(ij.min!=null) ij.min=+ij.min-32768; if(ij.max!=null) ij.max=+ij.max-32768; }
    if(gi===0&&sem&&sem.raw) im.metaText=sem.raw;
    if(gi===0&&desc) im.desc=desc.length>20000?desc.slice(0,20000)+"…":desc;
    if(gi===0&&sem&&sem.imageRows&&sem.imageRows<g.H) im.databar={rows:g.H-sem.imageRows,src:"FEI [Image] ResolutionY"};
    if(gi===0&&ij){ if(ij.spacing){ const ku=IO.unitToM(String(ij.unit||"").replace(/\\u00B5/g,"µ")); if(ku) im.zStep=+ij.spacing*ku/1e-6; } /* spacing 은 ImageJ unit 단위 → µm */ if(ij.finterval) im.frameInterval=+ij.finterval; if(ij.min!=null&&ij.max!=null&&!(+ij.channels>1)) im.display=[{min:+ij.min,max:+ij.max}]; }
    images.push(im);
  });
  return {images,warn};
}

/* ============================================================ PNG ============================================================ */
async function readPNG(name,u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); let p=8; let W,H,bd,ct,il; const idat=[]; let plte=null, phys=null;
  while(p<u8.length){ const len=dv.getUint32(p); const type=td.decode(u8.subarray(p+4,p+8)); const data=u8.subarray(p+8,p+8+len);
    if(type==="IHDR"){ W=dv.getUint32(p+8); H=dv.getUint32(p+12); bd=u8[p+16]; ct=u8[p+17]; il=u8[p+20]; }
    else if(type==="PLTE") plte=data; else if(type==="IDAT") idat.push(data); else if(type==="pHYs"){ const q=new DataView(data.buffer,data.byteOffset,data.byteLength); phys={x:q.getUint32(0),y:q.getUint32(4),unit:data[8]}; } else if(type==="IEND") break;
    p+=12+len; }
  if(!(W>0&&H>0)) throw new Error("PNG: IHDR 가 없거나 크기가 0 입니다");
  const nch={0:1,2:3,3:1,4:2,6:4}[ct]; if(!nch) throw new Error("PNG: 알 수 없는 color type "+ct); const bpp=Math.max(1,(nch*bd)>>3);
  if(W*H*Math.min(nch,3)*4>IO.MAXBYTES) throw new Error(`PNG: 너무 큼 — ${W}×${H} (32-bit 메모리 ${(W*H*Math.min(nch,3)*4/1e9).toFixed(2)} GB, 한도 ${(IO.MAXBYTES/1e9).toFixed(1)} GB)`);
  let tot=0; for(const d of idat) tot+=d.length; const z=new Uint8Array(tot); let o=0; for(const d of idat){ z.set(d,o); o+=d.length; }
  const raw=await IO.inflate(z);
  { let need=0; const passes=il?[[0,0,8,8],[4,0,8,8],[0,4,4,8],[2,0,4,4],[0,2,2,4],[1,0,2,2],[0,1,1,2]]:[[0,0,1,1]]; for(const [x0,y0,dx,dy] of passes){ const pw=Math.ceil((W-x0)/dx), ph=Math.ceil((H-y0)/dy); if(pw>0&&ph>0) need+=ph*(1+Math.ceil(pw*nch*bd/8)); }
    if(raw.length<need) throw new Error(`PNG: 화소 자료가 짧습니다 (${raw.length} < ${need} 바이트) — 파일이 잘렸거나 손상됨`); }
  const outC=[...Array(nch)].map(()=>new Float32Array(W*H));
  let rp=0;
  const pass=(x0,y0,dx,dy)=>{ const pw=Math.ceil((W-x0)/dx), ph=Math.ceil((H-y0)/dy); if(pw<=0||ph<=0) return; const rb=Math.ceil(pw*nch*bd/8); let prev=new Uint8Array(rb);
    for(let y=0;y<ph;y++){ const f=raw[rp++]; const line=raw.slice(rp,rp+rb); rp+=rb;
      for(let i=0;i<rb;i++){ const a=i>=bpp?line[i-bpp]:0, b=prev[i], c=i>=bpp?prev[i-bpp]:0; let v=line[i];
        if(f===1) v+=a; else if(f===2) v+=b; else if(f===3) v+=(a+b)>>1; else if(f===4){ const pp=a+b-c, pa=Math.abs(pp-a), pb=Math.abs(pp-b), pc=Math.abs(pp-c); v+=(pa<=pb&&pa<=pc)?a:(pb<=pc?b:c); }
        line[i]=v&255; }
      prev=line;
      for(let x=0;x<pw;x++) for(let c=0;c<nch;c++){ let v; const si=x*nch+c;
        if(bd===8) v=line[si]; else if(bd===16) v=(line[2*si]<<8)|line[2*si+1]; else { const bp=si*bd; v=(line[bp>>3]>>(8-bd-(bp&7)))&((1<<bd)-1); }
        outC[c][(y0+y*dy)*W+x0+x*dx]=v; } } };
  if(il) [[0,0,8,8],[4,0,8,8],[0,4,4,8],[2,0,4,4],[0,2,2,4],[1,0,2,2],[0,1,1,2]].forEach(a=>pass(...a)); else pass(0,0,1,1);
  const warn=[]; let planes, rgb=false, bit=bd===16?16:8;
  if(ct===3){ const r=new Float32Array(W*H),g=new Float32Array(W*H),b=new Float32Array(W*H); let gray=true; for(let i=0;i<plte.length/3;i++) if(plte[3*i]!==plte[3*i+1]||plte[3*i]!==plte[3*i+2]){ gray=false; break; }
    for(let i=0;i<W*H;i++){ const k=outC[0][i]; r[i]=plte[3*k]; g[i]=plte[3*k+1]; b[i]=plte[3*k+2]; } planes=gray?[r]:[r,g,b]; rgb=!gray; bit=8; }
  else if(ct===0||ct===4){ planes=[outC[0]]; if(bd<8){ const mx=(1<<bd)-1; for(let i=0;i<W*H;i++) planes[0][i]=Math.round(planes[0][i]*255/mx); } if(ct===4) warn.push("알파 채널은 버렸습니다"); }
  else { planes=outC.slice(0,3); if(ct===6) warn.push("알파 채널은 버렸습니다");
    let gray=true; for(let i=0;i<W*H;i+=7) if(planes[0][i]!==planes[1][i]||planes[0][i]!==planes[2][i]){ gray=false; break; }
    if(gray){ for(let i=0;i<W*H;i++) if(planes[0][i]!==planes[1][i]||planes[0][i]!==planes[2][i]){ gray=false; break; } }
    if(gray){ planes=[planes[0]]; warn.push("R=G=B 인 회색 PNG 를 한 채널로 읽었습니다"); } else rgb=bit===8; }
  const cands=[]; if(phys&&phys.unit===1&&phys.x>0&&![2835,3780,3779,11811].includes(phys.x)) cands.push({pw:1e6/phys.x,unit:"µm",src:`PNG pHYs ${phys.x} px/m`});
  return {images:[{name,w:W,h:H,bit,rgb,nc:planes.length,nz:1,nt:1,order:"czt",planes,cal:null,calCandidates:cands,info:[["PNG",`${bd}-bit, color type ${ct}${il?", interlaced":""}`]],warn}],warn:[]};
}
/* 브라우저 해독 (JPEG·BMP·GIF·WebP·AVIF) */
async function readBitmap(name,u8,mime){
  if(typeof createImageBitmap==="undefined") throw new Error("이 환경에서는 "+name+" 형식을 해독할 수 없습니다");
  const bmp=await createImageBitmap(new Blob([u8],{type:mime||""}));
  const W=bmp.width,H=bmp.height; const cv=(typeof OffscreenCanvas!=="undefined")?new OffscreenCanvas(W,H):Object.assign(document.createElement("canvas"),{width:W,height:H});
  const cx=cv.getContext("2d"); cx.drawImage(bmp,0,0); const d=cx.getImageData(0,0,W,H).data;
  const r=new Float32Array(W*H),g=new Float32Array(W*H),b=new Float32Array(W*H); let gray=true;
  for(let i=0;i<W*H;i++){ r[i]=d[4*i]; g[i]=d[4*i+1]; b[i]=d[4*i+2]; if(gray&&(r[i]!==g[i]||r[i]!==b[i])) gray=false; }
  const warn=[]; if(/\.jpe?g$/i.test(name)) warn.push("JPEG 는 손실 압축입니다 — 세기 측정에는 원본(무손실) 파일을 쓰세요");
  return {images:[{name,w:W,h:H,bit:8,rgb:!gray,nc:gray?1:3,nz:1,nt:1,order:"czt",planes:gray?[r]:[r,g,b],cal:null,calCandidates:[],info:[["Decoder","browser"]],warn}],warn:[]};
}

/* ============================================================ 앱의 기존 파서 재사용 (DM·MRC·EMD·BCF) ============================================================ */
function fromRec(rec,name,kind){
  const out=[]; const warn=[];
  for(const im of (rec.images||[])){
    const ds=rec.h5&&rec.h5.get(im.base+"/Data"); if(!ds) continue;
    const sh=ds.shape||[]; const H=sh[0], W=sh[1]; let N=sh.length>2?sh[2]:1;
    let val=ds.value; if(val&&val.length>W*H*N) N=Math.floor(val.length/(W*H)); // DM 스택(앱 파서는 첫 장 모양만 적음)
    let planes=[];
    const per=W*H;
    if(kind==="emd"&&N>1){ // Velox: [H,W,N] 프레임이 마지막 축
      const maxFrames=Math.max(1,Math.floor(4e8/(per*4)));
      if(N>maxFrames){ warn.push(`${im.title}: 프레임 ${N}장은 메모리에 다 올릴 수 없어 합친 이미지 1장으로 읽었습니다`); const acc=new Float32Array(per); const chunk=Math.max(1,Math.min(N,Math.floor(8e6/per)));
        for(let f0=0;f0<N;f0+=chunk){ const f1=Math.min(N,f0+chunk), blk=ds.slice([[0,H],[0,W],[f0,f1]]), nf=f1-f0; for(let i=0;i<per;i++){ let s=0; for(let k=0;k<nf;k++) s+=blk[i*nf+k]; acc[i]+=s; } } planes=[acc]; N=1; }
      else { const blk=ds.slice([[0,H],[0,W],[0,N]]); for(let k=0;k<N;k++){ const p=new Float32Array(per); for(let i=0;i<per;i++) p[i]=blk[i*N+k]; planes.push(p); } }
    } else {
      if(!val) val=ds.slice?ds.slice():null; if(!val) continue;
      for(let k=0;k<N;k++){ const p=new Float32Array(per); for(let i=0;i<per;i++) p[i]=val[k*per+i]; planes.push(p); }
    }
    // 형식: 정수이고 범위가 맞으면 8/16-bit
    let isInt=true, mn=Infinity, mx=-Infinity; for(const p of planes) for(let i=0;i<p.length;i+=Math.max(1,(p.length/20000)|0)){ const v=p[i]; if(v!==Math.floor(v)) isInt=false; if(v<mn) mn=v; if(v>mx) mx=v; }
    if(isInt){ for(const p of planes) for(let i=0;i<p.length;i++){ const v=p[i]; if(v!==Math.floor(v)){ isInt=false; break; } if(v<mn) mn=v; if(v>mx) mx=v; } }
    const dt=String(im.dtype||""); let bit=32; if(isInt&&mn>=0&&mx<=255&&/u?int8|u1|i1/.test(dt)) bit=8; else if(isInt&&mn>=0&&mx<=65535) bit=16; else if(isInt&&mn>=0&&mx<=255) bit=16;
    let cal=null; if(im.px&&im.px.x>0) cal=IO.calFromMeters(im.px.x,im.px.y||im.px.x,kind.toUpperCase()+" 메타데이터");
    else if(kind==="dm"&&im.dmScale>0&&im.dmUnit){ cal={pw:im.dmScale,ph:im.dmScale,unit:im.dmUnit,src:"DM 보정 (길이 단위 아님)"}; }
    const info=[]; if(im.det) info.push(["Detector",im.det]); const md=im.md||{};
    if(md.Instrument&&md.Instrument.InstrumentModel) info.push(["Instrument",md.Instrument.InstrumentModel]); if(md.Optics&&md.Optics.AccelerationVoltage) info.push(["HV (V)",md.Optics.AccelerationVoltage]);
    if(im.dmMag) info.push(["Mag (actual)",im.dmMag]);
    out.push({name:(rec.images.length>1?`${name} · ${im.title||im.det||""}`:name),w:W,h:H,bit,rgb:false,nc:1,nz:planes.length,nt:1,order:"czt",planes,cal,calCandidates:[],info,warn:[]});
  }
  if(!out.length) throw new Error("이미지가 들어 있지 않습니다");
  return {images:out,warn};
}

/* ============================================================ Igor IBW (v5) ============================================================ */
/* Asylum Research(Igor IBW) 채널 이름 → 물리 단위. 이름 = [Nap]<종류>[숫자](Trace|Retrace).
   Gwyddion 2.64 igorfile 이 같은 파일에 붙이는 단위와 같다(실측 Cypher: Height·ZSensor·Amplitude·Deflection = m, Phase = deg, Potential = V, Current = A, Frequency = Hz).
   모르는 종류는 null — 값은 그대로 두고 알린다(wave 하나에 단위가 하나라 추측하지 않음). */
const AR_UNIT=[[/^(Height|ZSensor|Z|Amplitude|Deflection|Amp)\d*$/i,"m"],[/^Phase\d*$/i,"deg"],[/^(Potential|SurfacePotential|Voltage|Bias|UserIn|UserCalc|Input)\w*$/i,"V"],
  [/^Current\d*$/i,"A"],[/^(Frequency|Freq)\d*$/i,"Hz"],[/^Capacitance\d*$/i,"F"]];
IO.arChannelUnit=function(label){ const b=String(label||"").replace(/^Nap/,"").replace(/(Trace|Retrace)$/,""); for(const [re,u] of AR_UNIT) if(re.test(b)) return u; return null; };

function readIBW(name,u8){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); let le=true; let ver=dv.getUint16(0,true); if(ver>10){ le=false; ver=dv.getUint16(0,false); }
  if(ver!==5) throw new Error("Igor IBW 버전 "+ver+" 은 아직 읽지 않습니다 (버전 5 만)");
  const i32=o=>dv.getInt32(o,le), u16=o=>dv.getUint16(o,le), f64=o=>dv.getFloat64(o,le);
  const formulaSize=i32(8), noteSize=i32(12), dataEUnitsSize=i32(16), dimEU=[0,1,2,3].map(k=>i32(20+4*k)), dimLab=[0,1,2,3].map(k=>i32(36+4*k));
  const wh=64; const npnts=i32(wh+12), type=u16(wh+16); const nDim=[0,1,2,3].map(k=>i32(wh+68+4*k)); const sfA=[0,1,2,3].map(k=>f64(wh+84+8*k)); const sfB=[0,1,2,3].map(k=>f64(wh+116+8*k));
  const dimUnits=[0,1,2,3].map(k=>td.decode(u8.subarray(wh+152+4*k,wh+156+4*k)).replace(/\0.*$/,"")); const bname=td.decode(u8.subarray(wh+28,wh+60)).replace(/\0.*$/,""); // WaveHeader5.bname 은 28 바이트 뒤
  const dataOff=wh+320; const tsz={2:4,4:8,8:1,16:2,32:4}; const base=type&~0x41; const uns=!!(type&0x40); if(type&1) throw new Error("복소수 wave 는 읽지 않습니다");
  const bs=tsz[base]; if(!bs) throw new Error("wave 자료형 "+type+" 은 읽지 않습니다");
  if(!(npnts>0)||dataOff+npnts*bs>u8.length) throw new Error(`IBW: 자료가 짧습니다 (${npnts} 점 × ${bs} 바이트) — 파일이 잘렸거나 손상됨`);
  const rd=i=>{ const o=dataOff+i*bs; switch(base){ case 2: return dv.getFloat32(o,le); case 4: return dv.getFloat64(o,le); case 8: return uns?dv.getUint8(o):dv.getInt8(o); case 16: return uns?dv.getUint16(o,le):dv.getInt16(o,le); case 32: return uns?dv.getUint32(o,le):dv.getInt32(o,le); } };
  let o=dataOff+npnts*bs+formulaSize; const note=tdu.decode(u8.subarray(o,o+noteSize)).replace(/\r/g,"\n"); o+=noteSize+dataEUnitsSize;
  const dimEUstr=[]; for(let k=0;k<4;k++){ dimEUstr.push(td.decode(u8.subarray(o,o+dimEU[k]))); o+=dimEU[k]; }
  const labels=[]; for(let k=0;k<4;k++){ const L=[]; for(let j=0;j*32<dimLab[k];j++) L.push(td.decode(u8.subarray(o+32*j,o+32*j+32)).replace(/\0.*$/,"")); labels.push(L); o+=dimLab[k]; }
  const nx=nDim[0]||npnts, ny=nDim[1]||1, nl=nDim[2]||1; if(ny<2) throw new Error("1D wave 는 이미지가 아닙니다");
  const iwarn=[]; if(nDim[3]>1) iwarn.push(`4D wave (${nDim.join("×")}) — 넷째 차원의 첫 묶음만 읽었습니다`);
  // Igor: 첫 차원(행)이 가로(x), 둘째 차원이 세로(y, 위로 증가) → 화면은 위가 큰 y
  const planes=[]; for(let l=0;l<nl;l++){ const p=new Float32Array(nx*ny); for(let y=0;y<ny;y++) for(let x=0;x<nx;x++) p[(ny-1-y)*nx+x]=rd(x+nx*y+nx*ny*l); planes.push(p); }
  const u0=dimUnits[0]||dimEUstr[0]; const k=IO.unitToM(u0); let cal=null; if(k&&sfA[0]!==0&&isFinite(sfA[0])) cal=IO.calFromMeters(Math.abs(sfA[0])*k,Math.abs(sfA[1]||sfA[0])*(IO.unitToM(dimUnits[1]||dimEUstr[1])||k),"Igor wave scaling (sfA, "+u0+")",u0); // 음의 간격 = 축 방향만 반대
  const info=[["Wave",bname]]; const nv={}; note.split(/\n/).forEach(l=>{ const i=l.indexOf(":"); if(i>0) nv[l.slice(0,i).trim()]=l.slice(i+1).trim(); });
  ["ImagingMode","ScanSize","ScanRate","SetPoint","DriveAmplitude","DriveFrequency","SpringConstant","InvOLS","Date","Time"].forEach(kk=>{ if(nv[kk]!=null) info.push([kk,nv[kk]]); });
  const chLab=labels[2]&&labels[2].length>1?labels[2].slice(1):null;
  // 값 단위: Asylum 은 한 wave 에 여러 채널(높이·위상·전위·전류…)을 담고 wave 단위는 하나(m) — 채널 이름으로 정한다 (Gwyddion igorfile 과 같은 결과)
  const dataUnits=td.decode(u8.subarray(wh+148,wh+152)).replace(/\0.*$/,"")||"";
  const isAR=nv.ImagingMode!=null||nv.MicroscopeModel!=null;
  const imgs=planes.map((p,l)=>{ const lab=chLab?chLab[l]:null; let u=null, uwarn=null;
    if(isAR&&lab){ u=IO.arChannelUnit(lab); if(u==null) uwarn=`채널 '${lab}' 의 물리 단위를 모릅니다 — 값은 파일의 수 그대로(단위 없음)`; }
    else if(nl===1) u=dataUnits||null;
    let valueUnit="", siScale=1, siUnit="", iinfo=info;
    if(u&&IO._unitParse&&IO._dispUnit){ const uo=IO._unitParse(u); const du=IO._dispUnit({p10:0,u:uo.u},p); const f=du.f*Math.pow(10,uo.p10);
      if(f!==1) for(let i=0;i<p.length;i++) p[i]*=f; valueUnit=du.unit; siScale=f; siUnit=du.base;
      iinfo=info.concat([["Value unit",f!==1?`${du.unit} (파일 ${u} × ${+f.toPrecision(6)})`:du.unit]]); }
    return {name:nl>1?`${name} · ${lab||"layer "+(l+1)}`:name,w:nx,h:ny,bit:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal,calCandidates:[],info:iinfo,warn:uwarn?[uwarn]:[],metaText:note,
      channelLabel:lab,valueUnit,siScale,siUnit}; });
  void sfB;
  return {images:imgs,warn:iwarn};
}
/* ============================================================ Nanonis SXM ============================================================ */
function readSXM(name,u8){
  let end=-1; for(let i=0;i<Math.min(u8.length-1,2e6);i++) if(u8[i]===0x1A&&u8[i+1]===0x04){ end=i; break; }
  if(end<0) throw new Error("SXM 머리 끝(\\x1A\\x04)을 찾지 못했습니다");
  const hdr=td.decode(u8.subarray(0,end)); const H={}; let cur=null;
  hdr.split(/\r?\n/).forEach(l=>{ const m=l.match(/^:(.+):$/); if(m){ cur=m[1]; H[cur]=[]; } else if(cur) H[cur].push(l); });
  const g=k=>(H[k]||[]).join("\n").trim();
  const [nx,ny]=g("SCAN_PIXELS").split(/\s+/).map(Number); const rng=g("SCAN_RANGE").split(/\s+/).map(Number); const dir=g("SCAN_DIR");
  const rows=(H["DATA_INFO"]||[]).map(s=>s.trim()).filter(Boolean); const head=rows[0]?rows[0].split(/\t/):[]; const ci=head.indexOf("Name"), di=head.indexOf("Direction"), ui=head.indexOf("Unit");
  const chans=rows.slice(1).map(r=>{ const c=r.split(/\t/); return {name:c[ci],unit:c[ui],both:/both/i.test(c[di])}; });
  const dv=new DataView(u8.buffer,u8.byteOffset+end+2,u8.byteLength-end-2); let off=0; const imgs=[];
  const cal=rng[0]>0?IO.calFromMeters(rng[0]/nx,(rng[1]||rng[0])/ny,"Nanonis SCAN_RANGE / SCAN_PIXELS"):null;
  for(const c of chans){ for(const d of (c.both?["fwd","bwd"]:["fwd"])){
    if(off+nx*ny*4>dv.byteLength) break; const p=new Float32Array(nx*ny);
    for(let y=0;y<ny;y++) for(let x=0;x<nx;x++){ const v=dv.getFloat32(off+4*(y*nx+x),false); const yy=dir==="up"?ny-1-y:y; const xx=d==="bwd"?nx-1-x:x; p[yy*nx+xx]=v; }
    off+=nx*ny*4; imgs.push({name:`${name} · ${c.name} (${d})`,w:nx,h:ny,bit:32,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal,calCandidates:[],info:[["Channel",c.name],["Unit",c.unit],["Direction",d],["SCAN_DIR",dir],["Bias (V)",g("BIAS")]],warn:[],valueUnit:c.unit}); } }
  if(!imgs.length) throw new Error("SXM 에 채널 자료가 없습니다");
  return {images:imgs,warn:[]};
}
/* ============================================================ NumPy NPY ============================================================ */
function readNPY(name,u8){
  if(!(u8[0]===0x93&&td.decode(u8.subarray(1,6))==="NUMPY")) throw new Error("NPY 머리가 아닙니다");
  const major=u8[6]; const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const hl=major>=2?dv.getUint32(8,true):dv.getUint16(8,true); const ho=major>=2?12:10;
  const h=td.decode(u8.subarray(ho,ho+hl)); const descr=(h.match(/'descr':\s*'([^']+)'/)||[])[1], fo=/'fortran_order':\s*True/.test(h), shape=((h.match(/'shape':\s*\(([^)]*)\)/)||[])[1]||"").split(",").map(s=>s.trim()).filter(Boolean).map(Number);
  if(fo) throw new Error("Fortran 순서 NPY 는 읽지 않습니다 (np.ascontiguousarray 로 저장)");
  if(!descr) throw new Error("NPY: 구조체(structured) dtype 은 이미지로 읽지 않습니다");
  const le=descr[0]!==">"; const t=descr.replace(/^[<>|=]/,""); const sz=+t.slice(1); const kind=t[0];
  if(!"fiub".includes(kind)||!(sz>0)||(kind==="f"&&![2,4,8].includes(sz))) throw new Error("NPY dtype '"+descr+"' 은 이미지로 읽지 않습니다 (실수·정수·bool 만)");
  const rd=(o)=>kind==="f"?(sz===4?dv.getFloat32(o,le):sz===8?dv.getFloat64(o,le):f16(dv.getUint16(o,le))):kind==="u"?(sz===1?dv.getUint8(o):sz===2?dv.getUint16(o,le):sz===4?dv.getUint32(o,le):Number(dv.getBigUint64(o,le))):kind==="i"?(sz===1?dv.getInt8(o):sz===2?dv.getInt16(o,le):sz===4?dv.getInt32(o,le):Number(dv.getBigInt64(o,le))):kind==="b"?dv.getUint8(o):NaN;
  let H,W,N=1,rgb=false,chLast=0; if(shape.length===2){ [H,W]=shape; }
  else if(shape.length===3&&(shape[2]===3||shape[2]===4)&&shape[0]>4&&shape[1]>4){ [H,W]=shape; chLast=shape[2]; rgb=kind==="u"&&sz===1; } // (H, W, 3|4) = 색 채널
  else if(shape.length===3){ [N,H,W]=shape; } else throw new Error("NPY 모양 "+JSON.stringify(shape)+" 은 이미지로 읽지 않습니다");
  const off=ho+hl, per=W*H; const planes=[]; if(off+per*N*Math.max(1,chLast)*sz>u8.length) throw new Error("NPY: 자료가 짧습니다 — 파일이 잘렸거나 손상됨");
  if(per*N*Math.min(3,Math.max(1,chLast))*4>IO.MAXBYTES) throw new Error("NPY: 너무 큼 (32-bit 메모리 한도 1.6 GB)");
  if(chLast){ for(let c=0;c<3;c++){ const p=new Float32Array(per); for(let i=0;i<per;i++) p[i]=rd(off+(i*chLast+c)*sz); planes.push(p); } }
  else for(let k=0;k<N;k++){ const p=new Float32Array(per); for(let i=0;i<per;i++) p[i]=rd(off+(k*per+i)*sz); planes.push(p); }
  const bit=kind==="u"&&sz===1||kind==="b"?8:(kind==="u"&&sz===2)?16:32;
  const w2=[]; if(chLast===4) w2.push("넷째 채널(알파)은 버렸습니다"); if(kind==="i"&&sz===8||kind==="u"&&sz===8) w2.push("64-bit 정수를 32-bit float 로 읽었습니다 (2²⁴ 넘는 값은 반올림)");
  return {images:[{name,w:W,h:H,bit,rgb,nc:chLast?3:1,nz:chLast?1:N,nt:1,order:"czt",planes,cal:null,calCandidates:[],info:[["dtype",descr],["shape",shape.join("×")]],warn:w2,...(chLast&&!rgb?{chNames:["R","G","B"],chLut:["Red","Green","Blue"]}:{})}],warn:[]};
}
/* ============================================================ 텍스트 행렬 ============================================================ */
function readText(name,u8){
  const s=tdu.decode(u8); const lines=s.split(/\r?\n/).filter(l=>l.trim()&&!/^\s*[#%]/.test(l)); if(!lines.length) throw new Error("빈 텍스트 파일입니다");
  // 구분자: 탭 > ; (이때 , 는 소수점) > , > 공백
  const smp=lines.slice(0,50).join("\n"); const delim=smp.includes("\t")?"\t":smp.includes(";")?";":smp.includes(",")?",":null; const commaDec=delim===";"&&/\d,\d/.test(smp);
  const num=t=>{ t=t.trim(); if(t==="") return NaN; if(commaDec) t=t.replace(",","."); if(/^[+-]?(nan)$/i.test(t)) return NaN; if(/^\+?inf(inity)?$/i.test(t)) return Infinity; if(/^-inf(inity)?$/i.test(t)) return -Infinity; const v=Number(t); return t!==""&&!isNaN(v)?v:null; };
  const rows=[]; let header=0, bad=0;
  for(const l of lines){ let tk=delim?l.split(delim):l.trim().split(/\s+/); if(delim&&tk.length>1&&tk[tk.length-1].trim()==="") tk=tk.slice(0,-1); // 줄 끝 구분자
    const v=tk.map(num); if(v.some(x=>x===null)){ if(!rows.length) header++; else bad++; continue; } rows.push(v); }
  const W=Math.max(0,...rows.map(r=>r.length)); const good=rows.filter(r=>r.length===W); const short=rows.length-good.length;
  if(good.length<2||W<2) throw new Error("숫자 행렬(텍스트 이미지)로 읽을 수 없습니다");
  const H=good.length, p=new Float32Array(W*H); good.forEach((r,y)=>r.forEach((v,x)=>p[y*W+x]=v));
  let isInt=true,mn=Infinity,mx=-Infinity; for(const v of p){ if(!isFinite(v)){ isInt=false; continue; } if(v!==Math.floor(v)) isInt=false; if(v<mn) mn=v; if(v>mx) mx=v; }
  const bit=isInt&&mn>=0&&mx<=255?8:isInt&&mn>=0&&mx<=65535?16:32;
  const w=[]; if(header) w.push(`머리 줄 ${header}개를 건너뛰었습니다`); if(bad) w.push(`숫자가 아닌 칸이 있는 줄 ${bad}개를 건너뛰었습니다 — 행 위치가 어긋날 수 있습니다`); if(short) w.push(`칸 수가 다른 줄 ${short}개를 건너뛰었습니다`); if(commaDec) w.push("';' 구분, ',' 소수점으로 읽었습니다");
  return {images:[{name,w:W,h:H,bit,rgb:false,nc:1,nz:1,nt:1,order:"czt",planes:[p],cal:null,calCandidates:[],info:[["Text image",`${W}×${H}${delim?` (구분자 '${delim==="\t"?"tab":delim}')`:""}`]],warn:w}],warn:[]};
}
/* ============================================================ Raw (사용자 지정) ============================================================ */
IO.readRaw=function(name,u8,o){
  const {w,h,n=1,type="u8",offset=0,gap=0,le=true}=o; const sz={u8:1,i8:1,u16:2,i16:2,u32:4,i32:4,f32:4,f64:8}[type]; if(!sz) throw new Error("자료형: "+type);
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const per=w*h; const planes=[];
  for(let k=0;k<n;k++){ const base=offset+k*(per*sz+gap); if(base+per*sz>u8.byteLength) break; const p=new Float32Array(per);
    for(let i=0;i<per;i++){ const q=base+i*sz; p[i]=type==="u8"?dv.getUint8(q):type==="i8"?dv.getInt8(q):type==="u16"?dv.getUint16(q,le):type==="i16"?dv.getInt16(q,le):type==="u32"?dv.getUint32(q,le):type==="i32"?dv.getInt32(q,le):type==="f32"?dv.getFloat32(q,le):dv.getFloat64(q,le); }
    planes.push(p); }
  if(!planes.length) throw new Error("파일이 지정한 크기보다 작습니다");
  const bit=type==="u8"?8:type==="u16"?16:32;
  return {images:[{name,w,h,bit,rgb:false,nc:1,nz:planes.length,nt:1,order:"czt",planes,cal:null,calCandidates:[],info:[["Raw",`${type} ${le?"little":"big"}-endian, offset ${offset}`]],warn:[]}],warn:[]};
};

/* ============================================================ 판별·읽기 ============================================================ */
IO.FORMATS=[
  ["TIFF / BigTIFF","tif,tiff,ome.tif,btf,tf8,gel","8·16·32-bit, float, RGB, 스택, ImageJ·OME 메타데이터, Zeiss·FEI·TESCAN SEM 태그"],
  ["PNG","png","1–16-bit (16-bit 그대로)"],["JPEG · BMP · GIF · WebP","jpg,jpeg,bmp,gif,webp,avif","브라우저 해독 (8-bit)"],
  ["Gatan DM3 / DM4","dm3,dm4","보정·스택"],["MRC / CCP4","mrc,mrcs,st,rec,ali,map","CELLA 보정"],["Velox EMD / HDF5","emd,h5,hdf5","프레임 스택"],["Bruker BCF","bcf","HeaderData 안 이미지"],
  ["Nikon ND2","nd2","NIS-Elements v2·v3: 채널·Z·T·XY 위치, lossless, 채널 이름·색, µm/px·z step (v1 JPEG 2000 제외)"],["Leica LIF","lif","LAS X·LAS AF: series·tile, 채널·Z·T(λ 포함), XZ 주사, LUT 색, µm/px·z step"],["Olympus/Evident OIR · POIR","oir,poir","FluoView FV3000·FV4000: 채널·Z·T·λ(T 에 이음), 8/12/16-bit·float, RGB, 선 주사, 참조 이미지, 채널 이름·파장, µm/px·z step (보조 파일 _00001… 은 함께 읽지 않음)"],["Zeiss CZI","czi","ZEN: 채널·Z·T·scene·mosaic tile 합성, 비압축·zstd·JPG, Scaling·채널 이름·색 (JPEG XR 제외)"],
  ["Gwyddion GWY · GSF","gwy,gsf","채널(DataField)·mask·volume(Brick), 제목·SI 단위·meta"],["Bruker Nanoscope","spm,000–999,pfc","채널·Trace/Retrace, Z scale × soft scale(nm/V 등) → 물리값, 16/32-bit (힘 곡선 파일 제외)"],
  ["Asylum ARDF","ardf","force map 이미지 채널(Height·Adhesion·Modulus 등), IDEF 간격·단위 (힘 곡선은 알리기만)"],
  ["Igor IBW (Asylum AFM)","ibw","채널별 이미지, wave scaling"],["Nanonis SXM","sxm","채널·방향별"],["NumPy","npy","2D·3D"],["텍스트 행렬","txt,csv,tsv,dat,asc","숫자 표"],["Raw","raw,bin","크기·형식 직접 지정"]];
IO.ACCEPT=".nd2,.czi,.lif,.oir,.poir,.gwy,.gsf,.spm,.ardf,.pfc,"+Array.from({length:1000},(_,i)=>"."+String(i).padStart(3,"0")).join(",")+",.tif,.tiff,.btf,.png,.jpg,.jpeg,.bmp,.gif,.webp,.avif,.dm3,.dm4,.mrc,.mrcs,.st,.rec,.ali,.map,.emd,.h5,.hdf5,.bcf,.ibw,.sxm,.npy,.txt,.csv,.tsv,.dat,.asc,.raw,.bin";
IO.sniff=function(name,u8){
  const n=name.toLowerCase(); const b=u8;
  if(b.length>=4&&((b[0]===0x49&&b[1]===0x49&&(b[2]===42||b[2]===43))||(b[0]===0x4D&&b[1]===0x4D&&b[2]===0&&(b[3]===42||b[3]===43)))) return "tiff";
  if(b[0]===0x89&&b[1]===0x50&&b[2]===0x4E&&b[3]===0x47) return "png";
  if(b[0]===0xFF&&b[1]===0xD8) return "jpeg"; if(b[0]===0x42&&b[1]===0x4D) return "bmp"; if(b[0]===0x47&&b[1]===0x49&&b[2]===0x46) return "gif";
  if(b[0]===0x52&&b[1]===0x49&&b[2]===0x46&&b[3]===0x46&&b[8]===0x57) return "webp";
  if(b[0]===0x93&&b[1]===0x4E&&b[2]===0x55) return "npy";
  if(b[0]===0x89&&b[1]===0x48&&b[2]===0x44&&b[3]===0x46) return "hdf5";
  if(b.length>=4&&b[0]===0xDA&&b[1]===0xCE&&b[2]===0xBE&&b[3]===0x0A) return "nd2"; if(/\.nd2$/.test(n)&&b[0]===0&&b[1]===0&&b[2]===0&&b[3]===0x0C) return "nd2";
  if(b.length>=10&&String.fromCharCode(...b.subarray(0,10))==="ZISRAWFILE") return "czi";
  if(b.length>=16&&String.fromCharCode(...b.subarray(0,16))==="OLYMPUSRAWFORMAT") return "oir"; /* 보조 파일(이름_00001)은 확장자가 없음 */
  if(/\.poir$/.test(n)&&b[0]===0x50&&b[1]===0x4B&&b[2]===3&&b[3]===4) return "poir";
  { const h=String.fromCharCode(...b.subarray(0,Math.min(b.length,32)));
    if(h.startsWith("GWYP")||h.startsWith("GWYO")) return "gwy"; if(h.startsWith("Gwyddion Simple Field 1.0\n")) return "gsf";
    if(/^[\\?]\*(File list|EC File list|Force file list)/.test(h)) return "nanoscope";
    if(b.length>=16&&h.slice(8,12)==="ARDF") return "ardf";
    if(b.length>=21&&b[0]===0x70&&b[1]===0&&b[2]===0&&b[3]===0&&b[8]===0x2A&&b[13]===0x3C&&b[14]===0&&b[15]===0x4C&&b[17]===0x4D) return "lif"; } /* 0x70 … 0x2A "<LM" (UTF-16) */
  if(/\.dm[34]$/.test(n)) return "dm"; if(/\.(mrc|mrcs|st|rec|ali|map)$/.test(n)) return "mrc"; if(/\.bcf$/.test(n)) return "bcf"; if(/\.(emd|h5|hdf5)$/.test(n)) return "hdf5";
  if(/\.ibw$/.test(n)) return "ibw"; if(/\.sxm$/.test(n)) return "sxm"; if(/\.(avif)$/.test(n)) return "bitmap";
  if(/\.(txt|csv|tsv|dat|asc)$/.test(n)) return "text"; if(/\.(raw|bin)$/.test(n)) return "raw";
  return null;
};
IO.readFile=async function(name,u8,ctx){
  ctx=ctx||{}; const k=IO.sniff(name,u8);
  switch(k){
    case "tiff": return readTIFF(name,u8,ctx);
    case "png": return readPNG(name,u8);
    case "jpeg": case "bmp": case "gif": case "webp": case "bitmap": return readBitmap(name,u8,{jpeg:"image/jpeg",bmp:"image/bmp",gif:"image/gif",webp:"image/webp"}[k]);
    case "npy": return readNPY(name,u8);
    case "nd2": return IO.readND2(name,u8,ctx);
    case "czi": return IO.readCZI(name,u8,ctx);
    case "gwy": return IO.readGWY(name,u8); case "gsf": return IO.readGSF(name,u8);
    case "nanoscope": return IO.readNanoscope(name,u8); case "ardf": return IO.readARDF(name,u8);
    case "lif": return IO.readLIF(name,u8,ctx);
    case "oir": return IO.readOIR(name,u8,ctx); case "poir": return IO.readPOIR(name,u8,ctx);
    case "ibw": return readIBW(name,u8);
    case "sxm": return readSXM(name,u8);
    case "text": return readText(name,u8);
    case "dm": if(typeof G.parseDM!=="function") throw new Error("DM 파서를 찾지 못했습니다"); return fromRec(G.parseDM(u8,name,u8.length),name,"dm");
    case "mrc": if(typeof G.parseMRC!=="function") throw new Error("MRC 파서를 찾지 못했습니다"); return fromRec(G.parseMRC(u8,name,u8.length),name,"mrc");
    case "bcf": if(typeof G.parseBcf!=="function") throw new Error("BCF 파서를 찾지 못했습니다"); return fromRec(await G.parseBcf(u8,name,u8.length),name,"bcf");
    case "hdf5": {
      if(typeof G.ensureH5!=="function"||typeof G.parseEmd!=="function") throw new Error("HDF5 엔진을 찾지 못했습니다");
      await G.ensureH5(); const vp="/imgtab/"+Date.now()+"_"+Math.random().toString(36).slice(2)+".emd"; try{ G.h5wasm.FS.mkdir("/imgtab"); }catch(e){}
      G.h5wasm.FS.writeFile(vp,u8); const rec=G.parseEmd(vp,name,u8.length); try{ return fromRec(rec,name,"emd"); } finally { try{ rec.h5&&rec.h5.close&&rec.h5.close(); G.h5wasm.FS.unlink(vp); }catch(e){} } }
    case "raw": { const e=new Error("Raw 파일은 크기와 자료형을 지정해야 합니다"); e.needRaw=true; throw e; }
    default: throw new Error("알 수 없는 형식입니다 ("+name+")");
  }
};

/* ============================================================ ImageJ 호환 TIFF 쓰기 ============================================================
   img: {w,h,bit,rgb,nc,nz,nt,planes,cal}  — little-endian, 비압축, ImageJ ImageDescription(보정·hyperstack) */
IO.writeTIFF=function(img,o){
  o=o||{}; const W=img.w,H=img.h, rgb=!!img.rgb, bit=rgb?8:img.bit, bps=bit/8, spp=rgb?3:1;
  const nImg=rgb?Math.max(1,img.planes.length/3):img.planes.length;
  const descL=["ImageJ=1.54u"]; if(nImg>1) descL.push("images="+nImg);
  if(rgb&&(img.nz>1||img.nt>1)){ if(img.nz>1) descL.push("slices="+img.nz); if(img.nt>1) descL.push("frames="+img.nt); if(img.nz>1&&img.nt>1) descL.push("hyperstack=true"); }
  if(!rgb&&(img.nc>1||img.nz>1||img.nt>1)){ if(img.nc>1) descL.push("channels="+img.nc); if(img.nz>1) descL.push("slices="+img.nz); if(img.nt>1) descL.push("frames="+img.nt); if(img.nc>1||img.nt>1) descL.push("hyperstack=true"); if(img.nc>1) descL.push("mode=composite"); }
  const cal=img.cal; const uOut=cal?(cal.unit==="µm"?"micron":cal.unit):null;
  if(cal&&cal.unit&&cal.unit!=="px") descL.push("unit="+(cal.unit==="µm"?"\\u00B5m":cal.unit));
  { const ku=cal&&IO.unitToM(cal.unit); if(img.zStep>0&&img.nz>1&&ku) descL.push("spacing="+(img.zStep*1e-6/ku)); } // zStep(µm) → ImageJ unit 단위
  if(o.min!=null&&o.max!=null) { descL.push("min="+o.min); descL.push("max="+o.max); }
  descL.push("loop=false"); const desc=descL.join("\n")+"\n\0"; void uOut;
  const descB=new TextEncoder().encode(desc);
  const nTags=rgb?12:12; const ifdSize=2+nTags*12+4;
  const pixBytes=W*H*bps*spp;
  // 배치: [머리 8][desc][res 16][bps 6][IFD...][자료...]
  const head=8, dOff=head, resOff=dOff+descB.length+(descB.length&1), bpsOff=resOff+16, ifd0=bpsOff+8;
  const ifdTotal=ifdSize*nImg; const data0=ifd0+ifdTotal; const total=data0+pixBytes*nImg;
  const buf=new ArrayBuffer(total), dv=new DataView(buf), u8=new Uint8Array(buf);
  u8[0]=0x49; u8[1]=0x49; dv.setUint16(2,42,true); dv.setUint32(4,ifd0,true);
  u8.set(descB,dOff);
  const xr=cal&&cal.pw>0?1/cal.pw:1, yr=cal&&cal.ph>0?1/cal.ph:1; const rat=(o,v)=>{ // RATIONAL: 분자가 2³² 를 넘지 않게 분모를 고른다
    const den=Math.max(1,Math.floor(4294967295/Math.max(1,v))); dv.setUint32(o,Math.min(4294967295,Math.round(v*den)),true); dv.setUint32(o+4,den,true); }; // 분모를 가능한 크게 (상대 오차 < 1e-9)
  rat(resOff,xr); rat(resOff+8,yr); if(rgb){ dv.setUint16(bpsOff,8,true); dv.setUint16(bpsOff+2,8,true); dv.setUint16(bpsOff+4,8,true); }
  for(let k=0;k<nImg;k++){
    const io=ifd0+k*ifdSize; let p=io; dv.setUint16(p,nTags,true); p+=2;
    const tag=(c,t,n,v)=>{ dv.setUint16(p,c,true); dv.setUint16(p+2,t,true); dv.setUint32(p+4,n,true); if(t===3&&n===1) dv.setUint16(p+8,v,true); else dv.setUint32(p+8,v,true); p+=12; };
    tag(254,4,1,0); tag(256,4,1,W); tag(257,4,1,H);
    if(rgb) tag(258,3,3,bpsOff); else tag(258,3,1,bit);
    tag(259,3,1,1); tag(262,3,1,rgb?2:1);
    if(k===0) tag(270,2,descB.length,dOff); else tag(270,2,descB.length,dOff);
    tag(273,4,1,data0+k*pixBytes); tag(277,3,1,spp); tag(278,4,1,H); tag(279,4,1,pixBytes);
    tag(339,3,1,bit===32?3:1);
    dv.setUint32(p,k<nImg-1?io+ifdSize:0,true);
    const o0=data0+k*pixBytes;
    if(rgb){ const r=img.planes[3*k],g=img.planes[3*k+1],b=img.planes[3*k+2]; for(let i=0;i<W*H;i++){ u8[o0+3*i]=r[i]; u8[o0+3*i+1]=g[i]; u8[o0+3*i+2]=b[i]; } }
    else { const d=img.planes[k]; const cl=(v,m)=>v!==v?0:v<0?0:v>m?m:Math.floor(v+0.5); // 8/16-bit 은 ImageJ 처럼 반올림·자르기
      if(bit===8) for(let i=0;i<W*H;i++) u8[o0+i]=cl(d[i],255); else if(bit===16) for(let i=0;i<W*H;i++) dv.setUint16(o0+2*i,cl(d[i],65535),true); else for(let i=0;i<W*H;i++) dv.setFloat32(o0+4*i,d[i],true); }
  }
  // 해상도 태그는 뒤에 붙이지 않고 ImageJ 가 쓰는 282/283/296 을 IFD 에 넣지 못했으므로(12 개 고정) 첫 IFD 에만 덧붙인다
  return addResolution(u8,cal,resOff);
};
function addResolution(u8,cal,resOff){
  if(!cal||!(cal.pw>0)) return u8;
  // 첫 IFD 를 해상도 태그 3 개를 더한 새 IFD 로 다시 쓴다(파일 끝)
  const dv=new DataView(u8.buffer); const ifd0=dv.getUint32(4,true); const n=dv.getUint16(ifd0,true); const ents=[];
  for(let i=0;i<n;i++){ const o=ifd0+2+12*i; ents.push(u8.slice(o,o+12)); } const next=dv.getUint32(ifd0+2+12*n,true);
  const add=(c,t,cnt,v)=>{ const e=new Uint8Array(12), d=new DataView(e.buffer); d.setUint16(0,c,true); d.setUint16(2,t,true); d.setUint32(4,cnt,true); if(t===3) d.setUint16(8,v,true); else d.setUint32(8,v,true); return e; };
  ents.push(add(282,5,1,resOff),add(283,5,1,resOff+8),add(296,3,1,1));
  ents.sort((a,b)=>new DataView(a.buffer).getUint16(0,true)-new DataView(b.buffer).getUint16(0,true));
  const newOff=u8.length+(u8.length&1); const size=2+12*ents.length+4; const out=new Uint8Array(newOff+size); out.set(u8);
  const dvo=new DataView(out.buffer); dvo.setUint16(newOff,ents.length,true); ents.forEach((e,i)=>out.set(e,newOff+2+12*i)); dvo.setUint32(newOff+2+12*ents.length,next,true); dvo.setUint32(4,newOff,true);
  return out;
}

IO._tiffIFDs=tiffIFDs;
G.IMG_IO=IO;
if(typeof module!=="undefined") module.exports=IO;
})();
