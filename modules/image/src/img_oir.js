/* ============================================================
   INNO Analytics — Olympus/Evident OIR 읽기 (FluoView FV3000·FV4000 ".oir", Olympus Image Format Raw) + POIR(ZIP 묶음)
   · 구조: 머리 48 바이트 "OLYMPUSRAWFORMAT" … u64 파일 크기(@32), u64 block index 위치(@40). index = 0xFFFFFFFF 뒤 u64 블록 위치 목록(0 에서 끝).
     블록 = [u32 길이][u32 종류] + 내용: 0 METADATA(XML 여러 개 — fileinfo·lsmimage·lut·cameraimage…), 1 FRAMEPROPERTIES(면 XML),
     2 BMP(썸네일), 3 UID(u32 ×3 뒤 ASCII UID — 바로 다음 index 항목의 PIXEL 블록 이름), 4 PIXEL(화소).
     UID = [l###][z###][t###][_a_b]_<채널 GUID>_<블록 번호> — 한 면이 PIXEL 블록 여러 개로 나뉘면 번호 차례로 이어 붙인다. REF_… = 참조(개관) 이미지.
   · 크기: FRAMEPROPERTIES 의 width·height·depth(바이트/화소: 1 → 8-bit, 2 → 16-bit, 4 → 32-bit float)·bitCounts·colorType,
     채널 = lsmimage 의 channel(id·order·name·파장), 화소 크기 = lsmimage 의 첫 length 요소(x·y, µm), Z 간격 = 면 위치(axisValue ZSTACK)
     또는 축 설정 step.
   · 앱 배치: c + nc·(z + nz·t). λ(파장 주사) 축은 T 에 잇는다 (t' = λ + nλ·t, LIF 읽기와 같은 방식).
   · 시험 (test/test_oir.js): 공개 실측 FV4000 시료 13 개(ref/real/oir, CC BY 4.0 — ref/gen_oir_real.py) 의 모든 면이 oirfile = Bio-Formats 7.3.1
     (= 같은 출처의 .ome.tif 3 개) 와 같음 (λ 는 oirfile — Bio-Formats 는 λ1 면을 되풀이), 화소 크기·z step·frame interval·채널 이름·참조 이미지까지.
     경계 사례는 실측 배치를 본뜬 합성 파일(ref/gen_oir_synth.py). 실측으로 확인 못 한 것: 선 주사·RGB·깊이 4 바이트·보조 파일(_00001)·POIR.

   Format knowledge from oirfile (BSD-3-Clause, Copyright (c) 2025-2026 Christoph Gohlke)
   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
   following conditions are met:
   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided with the distribution.
   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
   BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IO=G.IMG_IO; if(!IO) throw new Error("IMG_IO 가 먼저 있어야 합니다");
const tdu=new TextDecoder("utf-8"), tdl=new TextDecoder("latin1"), tdf=new TextDecoder("utf-8",{fatal:true});
const MAGIC="OLYMPUSRAWFORMAT";
const B_META=0, B_FRAME=1, B_BMP=2, B_UID=3, B_PIXEL=4;
const UID_RE=/^(?:l(\d+))?(?:z(\d+))?(?:t(\d+))?(?:_(\d+)_(\d+))?_([0-9a-f-]+)_(\d+)$/i; // GUID 는 대소문자 무시 (소문자로 맞춰 비교)
const AXN={ZSTACK:"Z",TIMELAPSE:"T",LAMBDA:"λ"};
const tick=()=>new Promise(r=>typeof setTimeout==="function"?setTimeout(r,0):r()); // 진행 막대가 다시 그려지도록 macrotask 로 양보

/* XML: img_lif.js 의 작은 XML 나무(IO._xmlTree)를 쓰고, 요소는 접두어를 뺀 이름(local name)으로 찾는다 (oirfile 과 같음) */
const local=t=>{ const i=t.lastIndexOf(":"); return i<0?t:t.slice(i+1); };
function tree(s){ if(typeof IO._xmlTree!=="function") throw new Error("XML 해석기(img_lif.js)가 먼저 있어야 합니다"); return IO._xmlTree(s); }
function* iter(e){ yield e; for(const c of e.c) yield* iter(c); }   // 문서 순서 (자기 자신 포함)
const firstText=(root,nm)=>{ for(const e of iter(root)) if(local(e.tag)===nm&&e.t) return e.t; return ""; };
const kidsL=(e,nm)=>e.c.filter(c=>local(c.tag)===nm);
const num=v=>{ const x=parseFloat(v); return Number.isFinite(x)?x:null; };
const u64=(dv,o)=>dv.getUint32(o,true)+dv.getUint32(o+4,true)*4294967296;
function calOpt(pwm,phm,src){ const u=pwm>=1e-8?["µm",1e6]:["nm",1e9]; return {pw:pwm*u[1],ph:(phm||pwm)*u[1],unit:u[0],src}; }
const g6=v=>+(+v).toPrecision(6);

/* METADATA 블록 안의 XML 문서들 — '<?xml' 마다 하나, 끝의 XML 아닌 바이트는 마지막 '>' 뒤에서 자름 (oirfile 과 같음) */
function splitXml(bytes){ const out=[]; const n=bytes.length; const dv=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength);
  const starts=[]; for(let i=0;i+5<=n;i++) if(bytes[i]===60&&bytes[i+1]===63&&bytes[i+2]===120&&bytes[i+3]===109&&bytes[i+4]===108) starts.push(i); // "<?xml"
  starts.forEach((p,k)=>{ const q=k+1<starts.length?starts[k+1]:n; let x=tdu.decode(bytes.subarray(p,q)); const e=x.lastIndexOf(">"); if(e>=0) x=x.slice(0,e+1);
    // 실측 파일: LUT 등 채널별 XML 앞에 [u32 36][채널 GUID 36 글자][u32 XML 길이] 가 있다
    let guid=null; if(p>=44&&dv.getUint32(p-44,true)===36){ const g=tdl.decode(bytes.subarray(p-40,p-4)); if(/^[0-9a-f-]{36}$/i.test(g)) guid=g.toLowerCase(); }
    out.push({x,guid}); });
  return out; }
function xmlKind(x){ const h=x.slice(0,200); for(const k of ["fileinfo","lsmimage","annotation","overlay","lut","base","event","cameraimage"]) if(h.includes(k+":")) return k; return null; }

/* 블록 목록: index 를 먼저 (0xFFFFFFFF 표지가 없어도 위치들이 모두 블록 머리를 가리키면 씀), 없거나 망가졌으면 앞에서부터 차례로 (잘린 파일) */
function blockList(u8,dv,warn){
  const N=u8.length; const fsize=u64(dv,32), ipos=u64(dv,40);
  if(fsize>N) warn.push(`파일이 머리에 적힌 크기보다 짧습니다 (${N.toLocaleString()} / ${fsize.toLocaleString()} 바이트) — 끝이 잘렸습니다`);
  if(ipos>=48&&ipos+12<=N){ const cand=[];
    for(let p=ipos+4;p+8<=N;p+=8){ const o=u64(dv,p); if(o===0||o>fsize) break; cand.push(o); }
    const marker=dv.getUint32(ipos,true)===0xFFFFFFFF;
    const valid=cand.length>0&&cand.every(o=>o>=48&&o+8<=N&&dv.getUint32(o+4,true)<=5);
    if(cand.length&&(marker||valid)) return {offs:cand,scanned:false,noMarker:!marker}; }
  // 차례 읽기: [u32 길이][u32 종류] 가 그럴듯하면 블록, 아니면 4 바이트씩 넘겨 다음 머리를 찾는다 (블록 사이 틈·-1 표지)
  const offs=[]; let p=48, gap=0;
  while(p+8<=N&&gap<=64){ const len=dv.getUint32(p,true), typ=dv.getUint32(p+4,true); let next=null;
    if(typ===B_UID){ if(p+20<=N){ const ul=dv.getUint32(p+16,true); if(ul>0&&ul<512) next=p+20+ul; } }
    else if(typ===B_FRAME||typ===B_META){ if(len<=N&&tdl.decode(u8.subarray(p+8,Math.min(N,p+8+64))).includes("<?xml")) next=p+8+len; }
    else if(typ===B_PIXEL){ if(len<=N) next=p+8+len; }
    else if(typ===B_BMP){ if(len<=N&&tdl.decode(u8.subarray(p+8,Math.min(N,p+12)))==="BMP ") next=p+8+len; } // 실측: 'BMP ' 로 시작
    else if(typ===5){ if(len===0) next=p+8; }
    if(next!=null&&typ===B_UID&&!/^[\x21-\x7e]+$/.test(tdl.decode(u8.subarray(p+20,Math.min(N,p+20+Math.min(64,next-p-20)))))) next=null; // UID 는 ASCII
    if(next==null){ p+=4; gap+=4; continue; }
    offs.push(p); if(next>=N) break; p=next; gap=0; }
  if(offs.length) warn.push("블록 목록(index)이 없습니다 — 파일 끝이 잘렸을 수 있습니다. 앞에서부터 블록을 차례로 읽었습니다");
  return {offs,scanned:true,noMarker:false};
}

function parseOIR(u8){
  const N=u8.length; if(N<48||tdl.decode(u8.subarray(0,16))!==MAGIC) throw new Error("OIR 머리(OLYMPUSRAWFORMAT)가 아닙니다");
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength); const warn=[];
  const {offs,scanned,noMarker}=blockList(u8,dv,warn);
  const pairs=[], frames=[], xml={}, xmlAll=[], lutBy=new Map(); let bmp=null, cutPix=0;
  for(let i=0;i<offs.length;i++){ const off=offs[i]; if(off+8>N) break; const len=dv.getUint32(off,true), typ=dv.getUint32(off+4,true);
    if(typ===B_UID){ if(off+20>N) break; const ul=dv.getUint32(off+16,true); if(off+20+ul>N) break; const uid=tdl.decode(u8.subarray(off+20,off+20+ul));
      if(i+1<offs.length){ const no=offs[i+1]; if(no+8<=N&&dv.getUint32(no+4,true)===B_PIXEL){ const pl=dv.getUint32(no,true); const avail=Math.max(0,Math.min(pl,N-(no+8))); if(avail<pl) cutPix++;
          pairs.push({uid,off:no+8,len:pl,avail,at:dv.getUint32(off+8,true),cnt:dv.getUint32(off+12,true)}); i++; continue; } } }
    else if(typ===B_FRAME){ const c=u8.subarray(off+8,Math.min(N,off+8+len)); const s=tdu.decode(c); const k=s.indexOf("<?xml"); if(k>=0) frames.push(s.slice(k).replace(/[\0\s]+$/,"")); }
    else if(typ===B_META){ for(const {x,guid} of splitXml(u8.subarray(off+8,Math.min(N,off+8+len)))){ const k=xmlKind(x); xmlAll.push(x); if(k) (xml[k]=xml[k]||[]).push(x); if(k==="lut"&&guid&&!lutBy.has(guid)) lutBy.set(guid,x); } }
    else if(typ===B_BMP){ bmp={off:off+8,len}; } }
  return {dv,warn,pairs,frames,xml,xmlAll,lutBy,bmp,cutPix,scanned,noMarker,truncated:u64(dv,32)>N,nBlocks:offs.length};
}

/* lsmimage: 채널·축·화소 크기
   · 채널 이름: channel 의 바로 아래 name 요소 (없으면 안쪽 첫 name). oirfile 은 안쪽 name 중 마지막(염료 이름 등)을 쓰고, Bio-Formats 는 바로 아래 name 을 쓴다 → Bio-Formats 와 같게.
   · 같은 축(ZSTACK 등)의 enable 정의가 여럿이면 paramEnable="true" 인 것 — 실측 FV4000 파일은 쓰인 정의(paramEnable="true", 예 step 2.67·maxSize 4) 뒤에
     쓰이지 않은 "Range" 정의(paramEnable="false", step 1·maxSize 21)를 함께 둔다 (oirfile 은 마지막 것을 써서 maxSize 21). paramEnable 이 없으면 마지막 것 (oirfile 과 같음).
   · 화소 크기 = 문서에서 첫 length 요소의 x·y (oirfile 과 같음), Z 설정 = length 요소의 z (Bio-Formats 가 z step 으로 쓰는 값). */
function parseLSM(x){
  const R={channels:[],axes:{},pxX:0,pxY:0,pxZ:0,date:"",scanner:null,imgW:0,imgH:0,chDepth:0,chBits:0,objective:""};
  if(!x) return R; let T; try{ T=tree(x); }catch(e){ return R; }
  const seen=new Set();
  for(const e of iter(T)){ const n=local(e.tag);
    if(n==="channel"&&e.a.id&&!seen.has(e.a.id.toLowerCase())){ const id=e.a.id.toLowerCase(); seen.add(id); let sw=null, ew=null;
      const dn=e.c.find(c=>local(c.tag)==="name"&&c.t); let name=dn?dn.t:"";
      for(const c of iter(e)){ const m=local(c.tag); if(m==="name"&&c.t&&!name) name=c.t; else if(m==="startWavelength"&&c.t&&sw==null) sw=num(c.t); else if(m==="endWavelength"&&c.t&&ew==null) ew=num(c.t);
        else if(m==="depth"&&c.t&&!R.chDepth) R.chDepth=+c.t; else if(m==="bitCounts"&&c.t&&!R.chBits) R.chBits=+c.t; }
      R.channels.push({id,order:+(e.a.order||0)||0,name,sw,ew}); }
    else if(n==="axis"&&e.a.enable!=null){ let type="", info={enable:e.a.enable==="true",param:e.a.paramEnable};
      for(const c of e.c){ const m=local(c.tag); if(m==="axis"&&c.t) type=c.t; else if(m==="startPosition"&&c.t) info.start=num(c.t); else if(m==="endPosition"&&c.t) info.end=num(c.t); else if(m==="step"&&c.t) info.step=num(c.t); else if(m==="maxSize"&&c.t) info.maxSize=parseInt(c.t,10); }
      if(type&&info.enable){ const old=R.axes[type]; if(!old||(old.param==="false"&&info.param!=="false")||(info.param==="true"&&old.param!=="true")||(old.param==null&&info.param==null)) R.axes[type]=info; } }
    else if((n==="objectiveLens"||n==="objective")&&!R.objective){ const nm=[...iter(e)].find(c=>c!==e&&/^(name|model|objectiveName)$/.test(local(c.tag))&&c.t); R.objective=nm?nm.t:e.t||""; } }
  R.channels.sort((a,b)=>a.order-b.order); // 같은 order 는 문서 순서 (안정 정렬)
  for(const e of iter(T)) if(local(e.tag)==="length"){ for(const c of e.c){ const m=local(c.tag); if(m==="x"&&c.t) R.pxX=num(c.t)||0; else if(m==="y"&&c.t) R.pxY=num(c.t)||0; } break; }
  for(const e of iter(T)) if(local(e.tag)==="length"){ const z=e.c.find(c=>local(c.tag)==="z"&&c.t); if(z){ R.pxZ=num(z.t)||0; break; } }
  R.date=firstText(T,"creationDateTime");
  for(const e of iter(T)) if(local(e.tag)==="imageInfo"){ for(const c of e.c){ const m=local(c.tag); if(m==="width"&&c.t) R.imgW=+c.t; else if(m==="height"&&c.t) R.imgH=+c.t; } break; }
  const st=firstText(T,"scannerType"); if(st){ for(const e of iter(T)) if(local(e.tag)==="scannerSettings"&&e.a.type===st){ const g=nm=>{ for(const c of iter(e)) if(local(c.tag)===nm&&c.t) return c.t; return ""; };
      R.scanner={type:st,seriesInterval:num(g("seriesInterval")),lineSpeed:num(g("lineSpeed")),frameSpeed:num(g("frameSpeed"))}; break; } }
  return R; }
function parseFileInfo(x){ if(!x) return null; let T; try{ T=tree(x); }catch(e){ return null; } const o={version:"",system:""};
  for(const e of iter(T)){ const m=local(e.tag); if(!e.t) continue; if(m==="version"&&!o.version) o.version=e.t; else if(/^system/i.test(m)&&!o.system) o.system=e.t; }
  return o.version||o.system?o:null; }
function parseFrame(x){ let T; try{ T=tree(x); }catch(e){ return null; }
  const f={w:+firstText(T,"width")||0,h:+firstText(T,"height")||0,depth:+firstText(T,"depth")||0,bits:+firstText(T,"bitCounts")||0,color:firstText(T,"colorType")};
  // 면 이름 (실측: 루트 id·base:name = "l001z002_0_1" — UID 앞부분과 같아 위치를 면 번호에 바로 붙일 수 있다)
  const id=(T.c[0]&&T.c[0].a.id)||firstText(T,"name"); const m=/^(?:l(\d+))?(?:z(\d+))?(?:t(\d+))?_\d+_\d+$/.exec(id||"");
  f.idx=m&&(m[1]||m[2]||m[3])?{l:m[1]?+m[1]:0,z:m[2]?+m[2]:0,t:m[3]?+m[3]:0}:null;
  f.pos=[]; for(const e of iter(T)) if(local(e.tag)==="axisValue"){ let t="", p=0; for(const c of e.c){ const m=local(c.tag); if(m==="axisType"&&c.t) t=c.t; else if(m==="position"&&c.t) p=num(c.t)||0; } if(t) f.pos.push([t,p]); }
  return f; }
const median=a=>{ if(!a.length) return 0; const b=a.slice().sort((x,y)=>x-y), k=b.length>>1; return b.length&1?b[k]:(b[k-1]+b[k])/2; };
/* Z 간격: 면 위치(axisValue ZSTACK)를 같은 T·λ 안에서만 비교한 간격의 중앙값 (T 사이 Z 표류에 흔들리지 않게).
   설정 값(length z → 없으면 ZSTACK 축 step)과 1 % 넘게 다르면 설정 값을 쓰고 알린다. 위치가 하나뿐이면 설정 값 (oirfile 은 이때 간격을 버림). */
function zStepOf(fps,lsm,warn){
  const d=[]; let first=null; const withIdx=fps.filter(q=>q&&q.idx&&q.pos.some(p=>p[0]==="ZSTACK"));
  if(withIdx.length>1){ // 면 이름으로: 같은 (t, λ) 안에서 z 번호 차례로 (위치 차 / 번호 차)
    const groups=new Map(); for(const q of withIdx){ const k=q.idx.t+"|"+q.idx.l; if(!groups.has(k)) groups.set(k,new Map()); groups.get(k).set(q.idx.z,q.pos.find(p=>p[0]==="ZSTACK")[1]); }
    for(const M of groups.values()){ const Z=[...M.keys()].sort((a,b)=>a-b); for(const z of Z){ const v=M.get(z); if(first==null||v<first) first=v; } for(let i=1;i<Z.length;i++) d.push((M.get(Z[i])-M.get(Z[i-1]))/(Z[i]-Z[i-1])); } }
  else { const groups=new Map(); for(const q of fps){ if(!q) continue; let z=null, key=""; for(const [a,p] of q.pos){ if(a==="ZSTACK") z=p; else key+=a+"="+p+";"; } if(z==null) continue; if(!groups.has(key)) groups.set(key,new Set()); groups.get(key).add(z); }
    for(const S of groups.values()){ const P=[...S].sort((a,b)=>a-b); if(first==null||P[0]<first) first=P[0]; for(let i=1;i<P.length;i++) d.push(P[i]-P[i-1]); } }
  const pos=Math.abs(median(d)); const ax=lsm.axes.ZSTACK; const set=lsm.pxZ>0?{v:lsm.pxZ,src:"lsmimage length z"}:ax&&ax.step?{v:Math.abs(ax.step),src:"ZSTACK 축 step"}:null;
  if(pos>0&&set&&Math.abs(pos-set.v)>0.01*set.v){ warn.push(`면 위치로 본 Z 간격 ${g6(pos)} µm 이 설정 값 ${g6(set.v)} µm (${set.src}) 과 달라 설정 값을 썼습니다`); return {v:set.v,src:set.src+" (면 위치와 다름)",first}; }
  if(pos>0) return {v:pos,src:"면 위치(axisValue)",first};
  if(set) return {v:set.v,src:set.src,first};
  return {v:0,src:"",first}; }
function axisScale(type,positions,ax){ const P=[...new Set(positions[type]||[])].sort((a,b)=>a-b);
  if(P.length>1) return {v:P[1]-P[0],pos:P}; if(ax&&ax.step) return {v:ax.step,pos:P}; return {v:0,pos:P}; }
/* 방출 파장(가운데) → 표시 LUT 이름 (채널의 LUT XML 이 없을 때) */
function lutOfWl(c){ if(!c||(c.sw==null&&c.ew==null)) return null; const w=c.sw!=null&&c.ew!=null?(c.sw+c.ew)/2:(c.sw!=null?c.sw:c.ew);
  return w<470?"Blue":w<505?"Cyan":w<560?"Green":w<590?"Yellow":w<700?"Red":"Magenta"; }
/* 깊이 4 바이트: float32 로 보면 대부분 비정규(subnormal)·NaN·Inf·아주 큰 값이면 uint32 로 읽는다 (bitCounts < 32 도 uint32 쪽 근거) */
function looksU32(bytes,bits){ const a=new Uint32Array(bytes.buffer,bytes.byteOffset,bytes.byteLength>>2); let nz=0, bad=0;
  for(let i=0;i<a.length;i++){ const w=a[i]; if(!(w&0x7fffffff)) continue; nz++; const e=(w>>>23)&255; if(e===0||e===255||e>=0xE0) bad++; }
  return nz>0&&(bad>0.5*nz||(bits>0&&bits<32&&bad>0)); }

async function readOIR(name,u8,ctx,opt){
  opt=opt||{}; const R=parseOIR(u8); const warn=R.warn; const dname=opt.inner?`${name} · ${opt.inner}`:name;
  const lsm=parseLSM(R.xml.lsmimage&&R.xml.lsmimage[0]); const finfo=parseFileInfo(R.xml.fileinfo&&R.xml.fileinfo[0]);
  const fr=R.frames.length?parseFrame(R.frames[0]):null;
  let W=fr&&fr.w||lsm.imgW, H=fr&&fr.h||lsm.imgH, depth=fr&&fr.depth||lsm.chDepth, bits=fr&&fr.bits||lsm.chBits; const color=(fr&&fr.color)||"";
  if(!fr&&W>0) warn.push("FRAMEPROPERTIES 가 없어 크기·형식을 lsmimage 설정에서 읽었습니다");
  // UID → 면 지도. 한 면의 PIXEL 블록은 block index 차례대로 잇는다 (oirfile 과 같음). UID 끝 번호는 파일 전체에서 세는 경우도 있어
  // (채널끼리 엇갈려 0,2,4 / 1,3,5) 차례·빈 번호로 자르지 않고, 같은 번호가 되풀이되면 첫 것만 쓰고, 번호가 줄면 알리기만 한다.
  const map=new Map(), refs=new Map(); let badUid=0, dup=0, back=0;
  for(const p of R.pairs){
    if(/^REF_/i.test(p.uid)){ const q=p.uid.split("_"); if(q.length>=4){ const g=q.slice(2,-1).join("_").toLowerCase(); if(!refs.has(g)) refs.set(g,[]); refs.get(g).push(p); } continue; }
    const m=UID_RE.exec(p.uid); if(!m){ badUid++; continue; }
    const l=m[1]?+m[1]:0, z=m[2]?+m[2]:0, t=m[3]?+m[3]:0, g=m[6].toLowerCase(), b=+m[7]; const k=t+"|"+l+"|"+z+"|"+g;
    if(!map.has(k)) map.set(k,{t,l,z,g,blocks:[],nums:new Set()}); const v=map.get(k);
    if(v.nums.has(b)){ dup++; continue; } if(v.blocks.length&&b<v.blocks[v.blocks.length-1].b) back++;
    v.nums.add(b); v.blocks.push({...p,b}); }
  if(badUid) warn.push(`이름을 해석할 수 없는 화소 블록 ${badUid} 개는 건너뛰었습니다`);
  if(dup) warn.push(`같은 UID 의 화소 블록이 ${dup} 번 되풀이되어 첫 것만 썼습니다`);
  if(back) warn.push(`면 안의 블록 번호가 ${back} 번 거꾸로 갑니다 — 파일에 적힌 차례대로 이었습니다`);
  const hasXml=R.xmlAll.length>0;
  if((!hasXml&&!R.scanned&&!R.truncated&&R.pairs.length)||(/_\d{5}$/.test(name)&&!R.xml.lsmimage))
    throw new Error("이 파일에는 메타데이터가 없습니다 — 큰 OIR 의 보조 파일(이름_00001 등)로 보입니다. 본 .oir 파일을 여세요");
  if(!map.size&&!refs.size){
    if(R.scanned) throw new Error("OIR: 파일이 잘려 화소 블록을 찾지 못했습니다 (block index 없음) — 파일이 끝까지 복사되었는지 확인하세요");
    if(!R.pairs.length&&fr&&fr.w*fr.h>1e8) throw new Error(`OIR: 지도(Map) 파일입니다 — ${fr.w}×${fr.h} 화소 크기의 개관 설정만 있고 화소 자료가 없습니다. 같은 폴더의 타일 파일(…_A01_G001_0001.oir 등)을 여세요`);
    throw new Error(R.pairs.length?"OIR: 화소 블록 이름(UID)을 해석하지 못했습니다":"OIR 안에 화소 자료가 없습니다 (메타데이터만 있는 파일)"); }
  if(!(W>0&&H>0)) throw new Error("OIR: 이미지 크기(width·height)를 찾지 못했습니다"+(R.scanned?" — 파일이 잘려 메타데이터가 없습니다":""));
  if(!R.xml.lsmimage) warn.push(`메타데이터(lsmimage)가 없습니다${R.truncated||R.scanned?" (파일이 잘림)":""} — 채널 이름·화소 크기·Z 간격을 알 수 없어 화소만 읽었습니다`);
  if(![1,2,4].includes(depth)) throw new Error(`OIR: 화소 깊이 ${depth||"?"} 바이트는 읽지 않습니다 (1·2·4 만)`);
  if(!map.size){ // 참조(개관) 이미지만
    const ri=refImage(R,lsm,refs,{W,H,depth,name:dname,u8}); if(!ri) throw new Error("OIR: 참조 이미지 자료가 비었습니다");
    ri.warn.unshift("이 파일에는 본 이미지(화소 블록)가 없고 참조(개관) 이미지만 있습니다 — 그것만 열었습니다"); return {images:[ri],warn}; }
  // 축 범위: UID 번호의 최소~최대 (가운데 빠진 면 = 0 으로 두고 알림; oirfile 은 있는 번호만 모아 당기고, Bio-Formats 는 번호 자리를 지킨다 → Bio-Formats 와 같게)
  const rng=key=>{ let lo=Infinity, hi=-Infinity; for(const v of map.values()){ lo=Math.min(lo,v[key]); hi=Math.max(hi,v[key]); } return [lo,hi]; };
  const [t0,t1]=rng("t"), [l0,l1]=rng("l"), [z0,z1]=rng("z"); const nt0=t1-t0+1, nl=l1-l0+1, nz=z1-z0+1;
  // 채널: lsmimage 의 order 순 (GUID 대소문자 무시), 메타데이터에 없는 GUID 는 뒤에 정렬해서
  const gset=new Set([...map.values()].map(v=>v.g)); let guids=[]; for(const c of lsm.channels) if(gset.has(c.id)){ guids.push(c.id); gset.delete(c.id); } guids.push(...[...gset].sort());
  const nc=guids.length;
  const rgb=/^rgb$/i.test(color)&&nc===3;
  let rgbNote="";
  if(rgb){ const nm=g=>{ const c=lsm.channels.find(x=>x.id===g); return c?c.name.toLowerCase():""; }; const idx=k=>guids.findIndex(g=>[["r","red"],["g","green"],["b","blue"]][k].includes(nm(g)));
    const ord=[0,1,2].map(idx); if(ord.every(i=>i>=0)&&new Set(ord).size===3&&ord.join()!=="0,1,2"){ guids=ord.map(i=>guids[i]); rgbNote=`RGB 성분을 채널 이름으로 맞췄습니다 (파일 차례 ${[0,1,2].map(i=>"RGB"[ord.indexOf(i)]).join("")} → RGB)`; } }
  // 세로: 면 자료가 width·height 보다 길면 — 선 주사(T 축 없이 TIMELAPSE maxSize = 줄 수, 면 높이 1 이거나 더 짧음)면 실제 줄 수 (oirfile 과 같음),
  // 아니면 FRAMEPROPERTIES height 를 지키고 남는 바이트는 버린다 (Bio-Formats 와 같음; oirfile 은 모든 면의 높이를 늘림)
  const fh=H; let Hreal=H, extra=0; for(const v of map.values()){ const tot=v.blocks.reduce((a,b)=>a+b.len,0); const h=Math.floor(tot/(W*depth)); if(h>fh) extra++; Hreal=Math.max(Hreal,h); }
  const tAx=lsm.axes.TIMELAPSE; const lineScan=nt0===1&&!!tAx&&tAx.maxSize===Hreal&&Hreal>1&&(fh===1||Hreal>fh);
  if(lineScan) H=Hreal; else if(extra) warn.push(`면 자료가 ${W}×${fh} 보다 긴 면이 ${extra} 개 있습니다 — FRAMEPROPERTIES 높이 ${fh} 줄만 읽고 남는 바이트는 버렸습니다`);
  const nt=nt0*nl; // λ 는 T 에 잇는다 (λ 가 빠르게)
  // 고르기·메모리
  let sel={z:null,t:null};
  if(ctx&&ctx.choose){ const r=await ctx.choose({kind:"oir",name:dname,positions:1,posNames:[],nz,nt,nc,w:W,h:H,bytesPerPlane:W*H*4,planesPer:nc}); if(!r) return {images:[],warn:[],cancel:true}; sel=Object.assign(sel,r); }
  const pick=(list,n)=>n<=1?[0]:list?list.filter(v=>v>=0&&v<n):[...Array(n).keys()];
  const zs=pick(sel.z,nz), ts=pick(sel.t,nt); if(!zs.length||!ts.length) throw new Error("OIR: 고른 Z·T 범위가 비었습니다");
  const tot=W*H*4*nc*zs.length*ts.length; const lim=(ctx&&ctx.maxBytes)||IO.MAXBYTES;
  if(tot>lim) throw new Error(`OIR: 고른 이미지가 너무 큼 (${(tot/1e9).toFixed(2)} GB, 32-bit 메모리 한도 ${(lim/1e9).toFixed(1)} GB) — Z·T 일부만 고르세요`);
  // 화소
  const per=W*H*depth; const buf=new Uint8Array(per);
  // 실측 파일의 UID 블록 머리 [u32 면 안 위치][u32 PIXEL 길이] 가 맞으면 그 위치에 놓는다 (빠진 블록이 있어도 뒤 블록이 밀리지 않음); 아니면 차례로 잇는다 (oirfile 과 같음)
  const placed=v=>v.blocks.length>0&&v.blocks.every(B=>B.cnt===B.len&&B.at+B.len<=per)&&v.blocks.some(B=>B.at>0)&&new Set(v.blocks.map(B=>B.at)).size===v.blocks.length;
  const fill=v=>{ buf.fill(0); let got=0, short=false;
    if(placed(v)){ for(const B of v.blocks){ buf.set(u8.subarray(B.off,B.off+B.avail),B.at); got+=B.avail; if(B.avail<B.len) short=true; } return short||got<per; }
    let o=0; for(const B of v.blocks){ if(o>=per) break; const n=Math.min(B.avail,per-o); buf.set(u8.subarray(B.off,B.off+n),o); o+=n; if(B.avail<B.len){ short=true; break; } } return short||o<per; };
  let u32=false;
  if(depth===4){ for(const v of map.values()){ fill(v); if(buf.some(x=>x)){ u32=looksU32(buf,bits); break; } } }
  const planes=new Array(nc*zs.length*ts.length); let missing=0, partial=0, done=0; const all=planes.length||1;
  for(let ti=0;ti<ts.length;ti++){ const tt=ts[ti], li=tt%nl, tI=Math.floor(tt/nl);
    for(let zi=0;zi<zs.length;zi++){
      for(let c=0;c<nc;c++){ const P=new Float32Array(W*H); planes[c+nc*(zi+zs.length*ti)]=P;
        const v=map.get((tI+t0)+"|"+(li+l0)+"|"+(zs[zi]+z0)+"|"+guids[c]);
        if(!v){ missing++; continue; }
        if(fill(v)) partial++;
        if(depth===1){ for(let i=0;i<W*H;i++) P[i]=buf[i]; }
        else if(depth===2) P.set(new Uint16Array(buf.buffer,0,W*H));
        else P.set(u32?new Uint32Array(buf.buffer,0,W*H):new Float32Array(buf.buffer,0,W*H)); }
      done++; if(ctx&&ctx.onProgress) ctx.onProgress(done*nc/all); if((done&15)===0) await tick(); } }
  // 보정
  const cal=lsm.pxX>0?calOpt(lsm.pxX*1e-6,(lineScan?lsm.pxX:(lsm.pxY||lsm.pxX))*1e-6,"OIR lsmimage length x·y (µm)"):null;
  const fps=R.frames.map(parseFrame); const positions={}; for(const q of fps) if(q) for(const [a,p] of q.pos) (positions[a]=positions[a]||[]).push(p);
  const w2=[];
  const zS=nz>1?zStepOf(fps,lsm,w2):{v:0}; let zStep=0;
  if(nz>1&&zS.v){ zStep=zS.v; if(zs.length>1){ const st=zs[1]-zs[0]; zStep=zs.every((x,i)=>i===0||x-zs[i-1]===st)?zStep*st:0; } }
  const bit=depth===1?8:depth===2?16:32;
  const info=[["Format",`Olympus OIR${opt.inner?` (${opt.inner})`:""}`]];
  if(finfo&&finfo.system) info.push(["System",finfo.system]); if(finfo&&finfo.version) info.push(["Software version",finfo.version]);
  if(lsm.objective) info.push(["Objective",lsm.objective]);
  if(lsm.date) info.push(["Acquired",lsm.date]);
  info.push(["Dimensions",`X${W} Y${H}${nl>1?` λ${nl}`:""}${nz>1?` Z${nz}`:""}${nt0>1?` T${nt0}`:""} · C${nc}`],["Bit depth",`${bits||depth*8}-bit${depth===4?(u32?" 정수":" float"):""} (${depth*8}-bit 저장)`]);
  if(depth===4){ const s=u32?"uint32 로 읽었습니다 (float 로 보면 값이 비정규수라)":"32-bit float 로 읽었습니다"; info.push(["깊이 4 바이트",`${s} — Bio-Formats 는 uint32, oirfile 은 float 로 읽음`]); w2.push(`화소 깊이 4 바이트를 ${u32?"uint32":"float"} 로 읽었습니다 (Bio-Formats 는 uint32, oirfile 은 float — 실측으로 확인 못 함)`); }
  const chMeta=guids.map(g=>lsm.channels.find(c=>c.id===g)||null);
  const chNames=guids.map((g,i)=>chMeta[i]&&chMeta[i].name?chMeta[i].name:`CH${i+1}`);
  chMeta.forEach((c,i)=>{ if(c&&(c.sw!=null||c.ew!=null)) info.push([`${chNames[i]} emission`,`${c.sw!=null?g6(c.sw):"?"}–${c.ew!=null?g6(c.ew):"?"} nm`]); });
  // 채널 색: 채널 GUID 가 붙은 LUT XML 의 lut:name (실측 FV4000: Green·Red·Cyan·Blue·Gray …), 없으면 방출 파장으로 추정
  const LUTN={red:"Red",green:"Green",blue:"Blue",cyan:"Cyan",magenta:"Magenta",yellow:"Yellow",gray:"Grays",grey:"Grays",grays:"Grays"};
  let fromLut=0, fromWl=0; const luts=guids.map((g,i)=>{ const x=R.lutBy.get(g); if(x){ let T=null; try{ T=tree(x); }catch(e){} const nm=T?firstText(T,"name").toLowerCase():""; if(LUTN[nm]){ fromLut++; return LUTN[nm]; } }
    const w=lutOfWl(chMeta[i]); if(w) fromWl++; return w; });
  if(fromLut||fromWl) info.push(["채널 색",fromLut?`파일의 LUT 이름${fromWl?" + 방출 파장 추정":""}`:"방출 파장으로 추정 (LUT 없음)"]);
  if(nz>1) info.push(["Z",`${nz} 장${zS.v?` · 간격 ${g6(zS.v)} µm (${zS.src})`:""}${zS.first!=null?` · 처음 ${g6(zS.first)} µm`:""}`]);
  if(nt0>1){ const P=positions.TIMELAPSE||[]; const U=[...new Set(P)].sort((a,b)=>a-b);
    const byT=new Map(); for(const q of fps) if(q&&q.idx){ const p=q.pos.find(x=>x[0]==="TIMELAPSE"); if(p&&(!byT.has(q.idx.t)||p[1]<byT.get(q.idx.t))) byT.set(q.idx.t,p[1]); } // 면 이름으로: 시점마다 첫 면의 시각
    const TT=[...byT.keys()].sort((a,b)=>a-b); const dts=[]; for(let i=1;i<TT.length;i++) dts.push((byT.get(TT[i])-byT.get(TT[i-1]))/(TT[i]-TT[i-1]));
    if(TT.length>1&&TT.length>=nt0-1) info.push(["Frame interval",`${g6(median(dts)/1000)} s (면 위치, 시점마다 첫 면)`]);
    else if(U.length===nt0) info.push(["Frame interval",`${g6((U[nt0-1]-U[0])/(nt0-1)/1000)} s (면 위치)`]);
    else if(P.length&&P.length%nt0===0&&P.length/nt0<P.length) info.push(["Frame interval",`${g6((P[P.length/nt0]-P[0])/1000)} s (면 위치)`]);
    else if(tAx&&tAx.step) info.push(["Frame interval",`${g6(tAx.step)} s (축 설정 step)`]); }
  if(nl>1){ const S=axisScale("LAMBDA",positions,lsm.axes.LAMBDA); info.push(["λ",`${nl} 단계${S.pos.length?` · ${g6(S.pos[0])} nm 부터`:""}${S.v?` · ${g6(S.v)} nm 간격`:""} — T 에 이어 놓음 (T = λ${nt0>1?" + "+nl+"·시간":""})`]); }
  if(lsm.scanner) info.push(["Scanner",lsm.scanner.type]);
  const nBlk=[...map.values()].reduce((a,v)=>Math.max(a,v.blocks.length),0); if(nBlk>1) info.push(["Pixel blocks",`면마다 ${nBlk} 개`]);
  if(R.noMarker) info.push(["Block index","0xFFFFFFFF 표지 없음 — 위치들이 블록 머리를 가리켜 그대로 씀"]);
  if(R.bmp) info.push(["Thumbnail","BMP 썸네일 있음 (읽지 않음)"]);
  if(lineScan){ w2.push("선 주사(line scan) — 세로 축 = 시간 (줄마다 한 번 주사). 화소 높이는 화소 너비와 같게 두었습니다 (실측 선 주사 파일로는 확인하지 못한 읽기)"); info.push(["이미지 면","X–시간 (kymograph)"]); if(lsm.scanner&&lsm.scanner.seriesInterval>0) info.push(["Line interval",`${g6(lsm.scanner.seriesInterval)} ms`]); }
  if(/^rgb$/i.test(color)&&nc!==3) w2.push(`colorType 이 RGB 인데 채널이 ${nc} 개라 채널별로 열었습니다`);
  if(rgbNote) w2.push(rgbNote);
  if(rgb) w2.push("RGB(colorType RGB) OIR 읽기는 실측 파일로 확인하지 못했습니다 — 색 순서를 FluoView 와 비교해 주세요");
  if(missing) w2.push(`화소 블록이 없는 면 ${missing} 개는 0 으로 채웠습니다`);
  if(partial) w2.push(`자료가 모자란 면 ${partial} 개의 빈 부분은 0 으로 채웠습니다${R.cutPix?" (파일이 잘림)":""}`);
  // 획득 설정보다 면이 적음: 중단된 획득이거나 나머지가 보조 파일(이름_00001 …)에 있음
  { const short=[]; for(const [a,n] of [["ZSTACK",nz],["TIMELAPSE",nt0],["LAMBDA",nl]]){ const ax=lsm.axes[a]; if(ax&&ax.maxSize>n&&!(a==="TIMELAPSE"&&lineScan)) short.push(`${AXN[a]} ${n}/${ax.maxSize}`); }
    if(short.length){ const base=name.replace(/\.[^.]+$/,""); w2.push(`획득 설정보다 면이 적습니다 (${short.join(", ")}) — 획득이 중간에 멈췄거나, 나머지 면이 보조 파일(${base}_00001 …)에 있습니다. 보조 파일은 함께 읽지 않습니다`); } }
  if(zs.length<nz) w2.push(`Z ${nz} 장 중 ${zs.length} 장만 읽었습니다`); if(ts.length<nt) w2.push(`T ${nt} 장 중 ${ts.length} 장만 읽었습니다`);
  const meta=[...(R.xml.lsmimage||[]),...R.xmlAll.filter(x=>!(R.xml.lsmimage||[]).includes(x)),...R.frames.slice(0,1)].join("\n");
  const im={name:dname,w:W,h:H,nc,nz:zs.length,nt:ts.length,order:"czt",planes,cal,calCandidates:[],info,warn:w2,zStep,
    metaText:meta.length>200000?meta.slice(0,200000)+"…":meta};
  if(rgb&&bit===8) Object.assign(im,{bit:8,rgb:true});
  else if(rgb) Object.assign(im,{bit,rgb:false,chNames:["R","G","B"],chLut:["Red","Green","Blue"]});
  else Object.assign(im,{bit,rgb:false,chNames,chLut:luts.some(Boolean)?luts:null});
  const images=[im];
  // 참조(개관) 이미지: REF_… 블록 — 크기는 imagedefinition → 정사각 → 본 이미지 높이·너비 순으로 맞춤 (oirfile 과 같음)
  if(refs.size&&!opt.noRef){ try{ const ri=refImage(R,lsm,refs,{W,H:fh,depth,name:dname,u8}); if(ri) images.push(ri); }catch(e){ warn.push("참조 이미지: "+e.message); } }
  return {images,warn};
}
function refImage(R,lsm,refs,o){
  let rdep=0; for(const x of R.xml.cameraimage||[]){ let T; try{ T=tree(x); }catch(e){ continue; } for(const e of iter(T)) if(local(e.tag)==="elementChannel"){ for(const c of e.c) if(local(c.tag)==="depth"&&+c.t>0){ rdep=+c.t; break; } if(rdep) break; } if(rdep) break; }
  const dep=rdep||o.depth; if(![1,2,4].includes(dep)) throw new Error(`깊이 ${dep} 바이트`);
  let rw=0, rh=0; const ids=R.xml.base&&R.xml.base[0]; if(ids){ const T=tree(ids); rw=+firstText(T,"width")||0; rh=+firstText(T,"height")||0; }
  const gs=[]; for(const c of lsm.channels) if(refs.has(c.id)) gs.push(c.id); gs.push(...[...refs.keys()].filter(g=>!gs.includes(g)).sort());
  const first=refs.get([...refs.keys()][0]); const tp=Math.floor(first.reduce((a,b)=>a+b.len,0)/dep); if(!tp) return null;
  const sq=Math.floor(Math.sqrt(tp)); let w,h; const rw2=[];
  if(rw>0&&rh>0&&rw*rh===tp){ w=rw; h=rh; } else if(sq*sq===tp){ w=h=sq; } else if(o.H>0&&tp%o.H===0){ h=o.H; w=tp/o.H; } else if(o.W>0&&tp%o.W===0){ w=o.W; h=tp/o.W; } else { w=tp; h=1; rw2.push("참조 이미지의 크기를 알 수 없어 한 줄로 열었습니다"); }
  let short=0;
  const planes=gs.map(g=>{ const P=new Float32Array(w*h), b=new Uint8Array(w*h*dep); let q=0; for(const B of refs.get(g)){ const n=Math.min(B.avail,b.length-q); if(n<=0) break; b.set(o.u8.subarray(B.off,B.off+n),q); q+=n; }
    if(q<b.length) short++;
    if(dep===1) for(let i=0;i<w*h;i++) P[i]=b[i]; else if(dep===2) P.set(new Uint16Array(b.buffer,0,w*h)); else P.set(new Float32Array(b.buffer,0,w*h)); return P; });
  if(short) rw2.push(`참조 이미지 채널 ${short} 개의 자료가 모자라 빈 부분을 0 으로 채웠습니다`);
  const names=gs.map((g,i)=>{ const c=lsm.channels.find(x=>x.id===g); return c&&c.name?c.name:`CH${i+1}`; });
  const cal=lsm.pxX>0?calOpt(lsm.pxX*1e-6,(lsm.pxY||lsm.pxX)*1e-6,"OIR lsmimage length x·y (µm) — 본 이미지의 값"):null;
  if(cal) rw2.push("참조 이미지(REF, FluoView 가 함께 저장한 주사 이미지)의 화소 크기는 파일에 따로 없어 본 이미지의 값을 붙였습니다 — 확대(zoom)를 바꿨다면 다를 수 있습니다");
  return {name:o.name+" · reference",w,h,nc:planes.length,nz:1,nt:1,order:"czt",planes,bit:dep===1?8:dep===2?16:32,rgb:false,chNames:names,chLut:null,cal,calCandidates:[],
    info:[["Format","Olympus OIR 참조(개관) 이미지"],["Dimensions",`X${w} Y${h} · C${planes.length}`]],warn:rw2,zStep:0}; }

/* ============================================================ POIR (ZIP 안의 OIR 여러 개) ============================================================
   메모리 한도는 묶음 전체에 (앞 파일이 쓴 만큼 줄여 넘김), 진행 막대는 항목마다 나눠, 고르기 창을 도중에 취소하면 이미 읽은 이미지는 남긴다. */
async function readPOIR(name,u8,ctx){
  ctx=ctx||{}; const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength), N=u8.length; let e=-1;
  for(let p=N-22;p>=Math.max(0,N-65557);p--) if(dv.getUint32(p,true)===0x06054b50){ e=p; break; }
  if(e<0) throw new Error("POIR: ZIP 끝 목록을 찾지 못했습니다 (잘린 파일?)");
  const Z64="ZIP64 묶음(4 GB 넘음)은 읽지 않습니다 — FluoView 에서 .oir 로 내보내 여세요";
  const cnt=dv.getUint16(e+10,true); let cd=dv.getUint32(e+16,true);
  if(cd===0xFFFFFFFF||cnt===0xFFFF) throw new Error("POIR: "+Z64);
  const ents=[]; for(let k=0;k<cnt&&cd+46<=N;k++){ if(dv.getUint32(cd,true)!==0x02014b50) break;
    const flag=dv.getUint16(cd+8,true), meth=dv.getUint16(cd+10,true), csz=dv.getUint32(cd+20,true), usz=dv.getUint32(cd+24,true), nl=dv.getUint16(cd+28,true), xl=dv.getUint16(cd+30,true), cl=dv.getUint16(cd+32,true), lo=dv.getUint32(cd+42,true);
    const nb=u8.subarray(cd+46,cd+46+nl); let fn; try{ fn=tdf.decode(nb); }catch(err){ fn=`항목 ${k+1} (${tdl.decode(nb)})`; } // bit 11 이 없어도 UTF-8 이면 그대로, 아니면 차례 번호
    ents.push({fn,flag,meth,csz,usz,lo,raw:tdl.decode(nb)}); cd+=46+nl+xl+cl; }
  const oirs=ents.filter(x=>/\.oir$/i.test(x.raw)); if(!oirs.length) throw new Error("POIR 안에 .oir 파일이 없습니다");
  const images=[], warn=[]; const lim=ctx.maxBytes||IO.MAXBYTES; let used=0;
  for(let k=0;k<oirs.length;k++){ const x=oirs[k]; try{
      if(x.flag&1) throw new Error("암호화된 항목은 읽지 않습니다");
      if(x.csz===0xFFFFFFFF||x.usz===0xFFFFFFFF||x.lo===0xFFFFFFFF) throw new Error(Z64);
      if(x.lo+30>N) throw new Error("자료가 파일 밖 (잘림)"); const p=x.lo+30+dv.getUint16(x.lo+26,true)+dv.getUint16(x.lo+28,true); if(p+x.csz>N) throw new Error("자료가 잘렸습니다");
      if(lim-used<=0) throw new Error("메모리 한도에 닿아 읽지 않았습니다");
      let d; if(x.meth===0) d=u8.subarray(p,p+x.csz); else if(x.meth===8) d=await IO.inflate(u8.subarray(p,p+x.csz),true); else throw new Error("압축 방식 "+x.meth);
      const sub=Object.assign({},ctx,{maxBytes:lim-used,onProgress:ctx.onProgress?(q=>ctx.onProgress((k+q)/oirs.length)):null});
      const r=await readOIR(name,d,sub,{inner:x.fn});
      if(r.cancel){ if(!images.length) return r; warn.push(`${x.fn} 부터는 열기를 취소했습니다 (앞의 ${images.length} 개만 열었습니다)`); break; }
      for(const im of r.images) used+=im.planes.length*im.w*im.h*4;
      images.push(...r.images); r.warn.forEach(w=>warn.push(`${x.fn}: ${w}`)); }
    catch(err){ warn.push(`${x.fn}: ${err.message} — 건너뜀`); } }
  if(!images.length) throw new Error("POIR: 읽은 이미지가 없습니다"+(warn.length?" — "+warn.join(", "):""));
  return {images,warn};
}

IO.readOIR=readOIR; IO.readPOIR=readPOIR; IO._parseOIR=parseOIR;
})();
