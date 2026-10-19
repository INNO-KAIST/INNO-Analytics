/* ============================================================
   INNO Analytics — Leica LIF 읽기 (IMG_IO 에 덧붙임)
   · 구조: [0x70, 크기, 0x2A, 글자 수, XML(UTF-16LE)] 다음에 memory block 이 잇달아 — [0x70, ?, 0x2A, 크기(v1 int32 · v2 int64), 0x2A, 글자 수, ID(UTF-16LE), 자료].
   · XML 의 Element 마다 Data/Image/ImageDescription 이 이미지 하나: Channels(ChannelDescription: Resolution·DataType·BytesInc·LUTName),
     Dimensions(DimensionDescription: DimID 1 X·2 Y·3 Z·4 T·5 λ·9 Λ·10 tile(M) …, NumberOfElements·Length·Unit·BytesInc), Memory(MemoryBlockID·Size).
   · 앞의 두 차원이 이미지 면(XZ 주사면 세로 = Z), 남은 공간 차원 하나 = Z, 시간·파장 차원(T·λ·Λ·그 밖)은 하나의 T 로 잇는다, tile(M)마다 이미지 따로 —
     Bio-Formats(LIFReader) 와 같은 배치. 화소 크기 = Length/(N−1) (Bio-Formats·liffile 과 같음).
   · 시험: 실측 LIF 9 파일(readlif·bioio-lif 시험 자료)의 모든 series 를 Bio-Formats 7.3.1 과 liffile(C. Gohlke) 이 읽은 값과 대조.
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IO=G.IMG_IO; if(!IO) throw new Error("IMG_IO 가 먼저 있어야 합니다");
const td16=new TextDecoder("utf-16le");
/* 작은 XML 나무: 요소·속성·글자(t, 앞뒤 공백 뺌)·원문 위치(s0 = '<' 위치, s1 = 닫는 태그 끝) */
function xmlTree(s){
  const root={tag:"#root",a:{},c:[],t:"",s0:0,s1:s.length}; const st=[root]; let i=0; const n=s.length;
  const cp=k=>k>=0&&k<=0x10FFFF&&!(k>=0xD800&&k<=0xDFFF)?String.fromCodePoint(k):"\uFFFD";
  const ent=v=>v.indexOf("&")<0?v:v.replace(/&(#x[0-9a-fA-F]+|#\d+|lt|gt|amp|quot|apos);/g,(m,e)=>e[0]==="#"?cp(e[1]==="x"?parseInt(e.slice(2),16):+e.slice(1)):{lt:"<",gt:">",amp:"&",quot:'"',apos:"'"}[e]);
  while(i<n){ const lt=s.indexOf("<",i); if(lt<0) break;
    if(lt>i){ const tx=s.slice(i,lt).trim(); if(tx){ const cur=st[st.length-1]; if(cur.t.length<4096) cur.t+=ent(tx); } }
    if(s.startsWith("<!--",lt)){ const e=s.indexOf("-->",lt+4); i=e<0?n:e+3; continue; }
    if(s.startsWith("<![CDATA[",lt)){ const e=s.indexOf("]]>",lt+9); const cur=st[st.length-1]; if(cur.t.length<4096) cur.t+=s.slice(lt+9,e<0?n:e).trim(); i=e<0?n:e+3; continue; }
    if(s[lt+1]==="?"||s[lt+1]==="!"){ const e=s.indexOf(">",lt); i=e<0?n:e+1; continue; }
    if(s[lt+1]==="/"){ const e0=s.indexOf(">",lt), e=e0<0?n:e0; const tag=s.slice(lt+2,e).trim(); for(let k=st.length-1;k>0;k--){ if(st[k].tag===tag){ st[k].s1=e+1; st.length=k; break; } } i=e+1; continue; }
    // 여는 태그: 따옴표 안의 '>' 를 건너뛰며 끝 찾기
    let j=lt+1, q=null; while(j<n){ const ch=s[j]; if(q){ if(ch===q) q=null; } else if(ch==='"'||ch==="'") q=ch; else if(ch===">") break; j++; }
    const body=s.slice(lt+1,j), self=body.endsWith("/"); const m=body.match(/^([^\s/>]+)/); if(!m){ i=j+1; continue; }
    const node={tag:m[1],a:{},c:[],t:"",s0:lt,s1:j+1}; const re=/([^\s=\/]+)\s*=\s*("([^"]*)"|'([^']*)')/g; let r; while((r=re.exec(body))) node.a[r[1]]=ent(r[3]!=null?r[3]:r[4]);
    st[st.length-1].c.push(node); if(!self) st.push(node); i=j+1; }
  return root;
}
const kid=(e,t)=>e&&e.c.find(x=>x.tag===t)||null;
const kids=(e,t)=>e?e.c.filter(x=>x.tag===t):[];
const path=(e,p)=>{ for(const t of p.split("/")){ e=kid(e,t); if(!e) return null; } return e; };
const LUTN={red:"Red",green:"Green",blue:"Blue",cyan:"Cyan",magenta:"Magenta",yellow:"Yellow",gray:"Grays",grey:"Grays"};
const DIMN={1:"X",2:"Y",3:"Z",4:"T",5:"λ",6:"rotation",7:"?",8:"?",9:"Λ",10:"M"};
function calOpt(pwm,phm,src){ const u=pwm>=1e-8?["µm",1e6]:["nm",1e9]; return {pw:pwm*u[1],ph:(phm||pwm)*u[1],unit:u[0],src}; }
const toM={m:1,mm:1e-3,um:1e-6,"µm":1e-6,nm:1e-9};

function lifStep(img){ let v=0; const go=e=>{ if(v) return; if(e.tag==="ScannerSettingRecord"&&e.a.Identifier==="dblStepSize"){ const x=parseFloat(e.a.Variant); if(x>0&&(e.a.Unit||"m")==="m") v=x; return; } for(const c of e.c) go(c); }; go(img); return v; }
async function readLIF(name,u8,ctx){
  const dv=new DataView(u8.buffer,u8.byteOffset,u8.byteLength), N=u8.length;
  if(N<13||dv.getUint32(0,true)!==0x70||u8[8]!==0x2A) throw new Error("LIF 머리(0x70 … 0x2A)가 아닙니다");
  const nch=dv.getUint32(9,true); if(13+2*nch>N) throw new Error("LIF: XML 머리가 파일보다 깁니다 (잘린 파일)");
  const xml=td16.decode(u8.subarray(13,13+2*nch)); const T=xmlTree(xml); const hdr=kid(T,"LMSDataContainerHeader");
  if(!hdr) throw new Error(xml.startsWith("<LMSDataContainerEnhancedHeader")?"LIF 확장 형식(LMSDataContainerEnhancedHeader)은 아직 읽지 않습니다":"LIF: LMSDataContainerHeader 가 없습니다");
  const ver=+hdr.a.Version||0; if(ver!==1&&ver!==2) throw new Error("LIF 판 "+hdr.a.Version+" 은 읽지 않습니다");
  // memory block 목록 (머리가 잘렸으면 거기서 멈추고 알림)
  const blocks=new Map(); let o=13+2*nch; const warn=[];
  while(o+9<=N){ if(dv.getUint32(o,true)!==0x70) break; if(u8[o+8]!==0x2A) throw new Error(`LIF: memory block 머리 손상 (바이트 ${o})`);
    let p=o+9, size; const hb=ver===2?8:4; if(p+hb+5>N){ warn.push(`memory block 머리가 파일 끝에서 잘렸습니다 (바이트 ${o})`); break; }
    if(ver===2){ size=dv.getUint32(p,true)+dv.getUint32(p+4,true)*4294967296; p+=8; } else { size=dv.getUint32(p,true); p+=4; }
    if(u8[p]!==0x2A) throw new Error(`LIF: memory block 머리 손상 (바이트 ${p})`); const sl=dv.getUint32(p+1,true); p+=5;
    if(p+2*sl>N){ warn.push(`memory block 이름이 파일 끝에서 잘렸습니다 (바이트 ${o})`); break; }
    const id=td16.decode(u8.subarray(p,p+2*sl)); p+=2*sl;
    blocks.set(id,{off:p,size,avail:Math.max(0,Math.min(size,N-p))}); o=p+size; }
  // 이미지 요소 (문서 순서)
  const series=[]; const walk=(el,pre)=>{ for(const e of kids(el,"Element")){ const nm=e.a.Name||""; const full=pre?pre+"/"+nm:nm;
      const img=path(e,"Data/Image"), mem=kid(e,"Memory");
      if(img&&kid(img,"ImageDescription")&&mem&&+mem.a.Size>0) series.push({el:e,img,mem,name:nm,path:full});
      const ch=kid(e,"Children"); if(ch) walk(ch,full); } };
  walk(hdr,"");
  if(!series.length) throw new Error("LIF 안에 이미지가 없습니다");
  // 이름: 잎 이름이 겹치면 경로(맨 위 요소 뺌)로
  { const cnt=new Map(); for(const s of series) cnt.set(s.name,(cnt.get(s.name)||0)+1);
    for(const s of series) if(cnt.get(s.name)>1){ const q=s.path.split("/"); s.name=q.length>1?q.slice(1).join("/"):s.path; } }
  const plan=[]; const skip=(s,why)=>warn.push(`'${s.name}': ${why} — 건너뜀`);
  const bpsOf=c=>c.dt===1?4:c.res<=8?1:c.res<=16?2:4;
  for(const s of series){
    const desc=kid(s.img,"ImageDescription"); const chEls=kids(kid(desc,"Channels"),"ChannelDescription"), dims=kids(kid(desc,"Dimensions"),"DimensionDescription");
    if(!chEls.length||dims.length<2){ skip(s,"채널·차원 정보가 없음"); continue; }
    const D=dims.map(d=>({id:+d.a.DimID,n:+d.a.NumberOfElements,len:+d.a.Length,unit:d.a.Unit||"",inc:+d.a.BytesInc}));
    const C=chEls.map((c,k)=>({el:c,k,res:+c.a.Resolution,dt:+c.a.DataType||0,inc:+c.a.BytesInc,lut:c.a.LUTName||"",tag:c.a.ChannelTag!=null?+c.a.ChannelTag:NaN}));
    const badD=D.find(d=>!(Number.isInteger(d.n)&&d.n>=1)||!(Number.isFinite(d.inc)&&d.inc>0)); if(badD){ skip(s,`차원 ${DIMN[badD.id]||badD.id} 의 NumberOfElements·BytesInc 가 맞지 않음 (${badD.n}, ${badD.inc})`); continue; }
    if(C.some(c=>!(Number.isFinite(c.inc)&&c.inc>=0))){ skip(s,"채널 BytesInc 가 맞지 않음"); continue; }
    const bps=bpsOf(C[0]);
    if(C.some(c=>c.dt!==C[0].dt||bpsOf(c)!==bps)){ skip(s,"채널마다 자료형이 다름"); continue; }
    const fx=D[0], fy=D[1]; const rest=D.slice(2);
    const zD=rest.find(d=>d.id===3)||(fy.id===3?rest.find(d=>d.id===2):null)||null;
    const mD=rest.find(d=>d.id===10)||null; const tDs=rest.filter(d=>d!==zD&&d!==mD); /* T·λ·Λ·그 밖 → 하나의 T */
    const nz=zD?zD.n:1, nM=mD?mD.n:1; let nt=1; for(const d of tDs) nt*=d.n;
    const blk=blocks.get(s.mem.a.MemoryBlockID); if(!blk){ skip(s,`memory block ${s.mem.a.MemoryBlockID} 이 파일에 없음 (잘린 파일?)`); continue; }
    const need=(fx.n-1)*fx.inc+(fy.n-1)*fy.inc+rest.reduce((a,d)=>a+(d.n-1)*d.inc,0)+Math.max(...C.map(c=>c.inc))+bps;
    if(need>blk.size){ skip(s,`차원·채널 크기(${need} 바이트)가 memory block(${blk.size})보다 큼`); continue; }
    // RGB (색 카메라): 화소 안에 세 채널이 붙어 있음 (X BytesInc = 3·bps). Leica 는 B·G·R 순서로 저장 (ChannelTag 3·2·1) → R·G·B 로
    let rgb=false, order=C.map((_,k)=>k);
    if(C.length===3&&fx.inc===3*bps&&new Set(C.map(c=>c.inc)).size===3&&C.every(c=>c.inc<3*bps)){ rgb=true;
      const tags=C.map(c=>c.tag); order=tags.every(t=>t>=1&&t<=3)&&new Set(tags).size===3?[1,2,3].map(t=>tags.indexOf(t)):C.map((_,k)=>k).sort((a,b)=>C[b].inc-C[a].inc); }
    for(let m=0;m<nM;m++) plan.push({s,D,C,bps,fx,fy,zD,mD,tDs,nz,nt,m,nM,blk,rgb,order});
  }
  if(!plan.length) throw new Error("LIF 안에 읽을 이미지가 없습니다"+(warn.length?" — "+warn.join(", "):""));
  // 고르기 (위치 = series·tile). series 마다 Z·T 가 다르면 가장 큰 값으로 묻고 series 마다 그 안의 것만
  const nc0=plan[0].C.length; const same=plan.every(p=>p.fx.n===plan[0].fx.n&&p.fy.n===plan[0].fy.n&&p.nz===plan[0].nz&&p.nt===plan[0].nt&&p.C.length===nc0);
  const posBytes=plan.map(p=>p.fx.n*p.fy.n*p.C.length*p.nz*p.nt*4);
  const NZ=Math.max(...plan.map(p=>p.nz)), NT=Math.max(...plan.map(p=>p.nt));
  let sel={p:plan.map((_,i)=>i),z:null,t:null};
  if(ctx&&ctx.choose){ const r=await ctx.choose({kind:"lif",name,positions:plan.length,posNames:plan.map(p=>p.s.name+(p.nM>1?` · tile ${p.m+1}`:"")),posBytes:same?null:posBytes,posNz:same?null:plan.map(p=>p.nz),posNt:same?null:plan.map(p=>p.nt),nz:NZ,nt:NT,nc:nc0,w:plan[0].fx.n,h:plan[0].fy.n,bytesPerPlane:plan[0].fx.n*plan[0].fy.n*4,planesPer:nc0});
    if(!r) return {images:[],warn:[],cancel:true}; sel=Object.assign(sel,r); }
  const pick=(list,n)=>n<=1?[0]:list?list.filter(v=>v<n):[...Array(n).keys()];
  const tot=sel.p.reduce((a,i)=>{ const P=plan[i]; return a+posBytes[i]*pick(sel.z,P.nz).length/P.nz*pick(sel.t,P.nt).length/P.nt; },0);
  const lim=(ctx&&ctx.maxBytes)||IO.MAXBYTES; if(tot>lim) throw new Error(`LIF: 고른 이미지가 너무 큼 (${(tot/1e9).toFixed(2)} GB, 32-bit 메모리 한도 ${(lim/1e9).toFixed(1)} GB) — 위치·Z·T 일부만 고르세요`);
  const images=[]; let done=0; const allPl=sel.p.reduce((a,i)=>a+pick(sel.z,plan[i].nz).length*pick(sel.t,plan[i].nt).length,0)||1;
  for(const pi of sel.p){ const P=plan[pi]; const {D,C,bps,fx,fy,zD,mD,tDs,nz,nt,blk,rgb,order}=P; const W=fx.n, H=fy.n, nc=C.length;
    const zs=pick(sel.z,nz), ts=pick(sel.t,nt);
    if(!zs.length||!ts.length){ warn.push(`'${P.s.name}': 고른 Z·T 범위가 이 series(Z ${nz} · T ${nt}) 밖이라 건너뜀`); continue; }
    const base=blk.off+(mD?P.m*mD.inc:0), end=blk.off+blk.avail; const planes=new Array(nc*zs.length*ts.length); let cut=0;
    const rd=bps===1?(q=>u8[q]):bps===2?(q=>dv.getUint16(q,true)):C[0].dt===1?(q=>dv.getFloat32(q,true)):(q=>dv.getUint32(q,true));
    for(let ti=0;ti<ts.length;ti++){ let r=ts[ti], toff=0; for(const d of tDs){ toff+=(r%d.n)*d.inc; r=Math.floor(r/d.n); } /* 앞 T 차원이 빠르게 (Bio-Formats 와 같은 순서) */
      for(let zi=0;zi<zs.length;zi++){ const zoff=zD?zs[zi]*zD.inc:0;
        for(let c=0;c<nc;c++){ const cc=order[c]; const pl=new Float32Array(W*H); const o0=base+toff+zoff+C[cc].inc; const last=o0+(H-1)*fy.inc+(W-1)*fx.inc+bps;
          if(last<=end){ if(bps===1&&fx.inc===1){ for(let y=0;y<H;y++){ const ro=o0+y*fy.inc; for(let x=0;x<W;x++) pl[y*W+x]=u8[ro+x]; } }
            else for(let y=0;y<H;y++){ const ro=o0+y*fy.inc; for(let x=0;x<W;x++) pl[y*W+x]=rd(ro+x*fx.inc); } }
          else { cut++; for(let y=0;y<H;y++){ const ro=o0+y*fy.inc; for(let x=0;x<W;x++){ const q=ro+x*fx.inc; pl[y*W+x]=q+bps<=end?rd(q):0; } } } /* 파일 끝 넘는 화소 = 0 (Bio-Formats 처럼) */
          planes[c+nc*(zi+zs.length*ti)]=pl; }
        if(ctx&&ctx.onProgress) ctx.onProgress(++done/allPl); if((done&7)===0) await null; } }
    const lenM=d=>d.n>1&&d.len>0?Math.abs(d.len)*(toM[d.unit]||NaN)/(d.n-1):NaN;
    const spatial=d=>d.id===1||d.id===2||d.id===3;
    const px=lenM(fx), py=spatial(fy)?lenM(fy):NaN; const cal=px>0?calOpt(px,py>0?py:px,`LIF DimensionDescription Length/(N−1) (${DIMN[fx.id]||fx.id}·${DIMN[fy.id]||fy.id})`):null;
    let zStep=0; if(zD){ let z=lenM(zD); const st=lifStep(P.s.img); if(st>0&&zD.id===3) z=st; /* ScannerSetting dblStepSize (Bio-Formats 와 같음 — Length/(N−1) 보다 자릿수가 많음) */ if(z>0){ zStep=z*1e6; if(zs.length>1){ const st=zs[1]-zs[0]; zStep=zs.every((v,i)=>i===0||v-zs[i-1]===st)?zStep*st:0; } } }
    const bit=C[0].dt===1?32:bps===1?8:bps===2?16:32;
    const info=[["Series",P.s.path+(P.nM>1?` (tile ${P.m+1}/${P.nM})`:"")],["Dimensions",D.map(d=>`${DIMN[d.id]||"dim"+d.id}${d.n}`).join(" ")+` · C${nc}`],["Bit depth",`${C[0].res}-bit${C[0].dt===1?" float":""} (${bps*8}-bit 저장)`]];
    const w2=[];
    if(fy.id!==2){ info.push(["Image plane",`${DIMN[fx.id]||fx.id}–${DIMN[fy.id]||fy.id} (세로 = ${DIMN[fy.id]||fy.id})`]);
      if(!spatial(fy)){ const st=fy.n>1&&fy.len>0?Math.abs(fy.len)/(fy.n-1):NaN; info.push([`세로 간격 (${DIMN[fy.id]||fy.id})`,st>0?`${+st.toPrecision(6)} ${fy.unit}`:"?"]);
        w2.push(`세로 축이 공간이 아닙니다 (${DIMN[fy.id]||fy.id}) — 화소 높이는 화소 너비와 같게 두었습니다 (세로 ${fy.unit||"단위"} 간격은 Image Info)`); } }
    if(tDs.length>1||tDs.some(d=>d.id!==4)) info.push(["T",`${tDs.map(d=>`${DIMN[d.id]||"dim"+d.id} ${d.n}`).join(" × ")} 을 T 하나로`]);
    if(zD&&zD.n===1&&tDs.some(d=>d.id!==4&&d.n>1)) info.push(["Dimension layout","Z 가 1 장이라 Bio-Formats 는 λ 등을 Z 로 둡니다 — 여기서는 T (화소는 같음)"]);
    const tD=tDs.find(d=>d.id===4); if(tD&&tD.n>1&&tD.unit==="s"&&tD.len>0) info.push(["Frame interval",`${+(tD.len/(tD.n-1)).toPrecision(6)} s`]);
    if(zs.length<nz) w2.push(`Z ${nz} 장 중 ${zs.length} 장만 읽었습니다`); if(ts.length<nt) w2.push(`T ${nt} 장 중 ${ts.length} 장만 읽었습니다`);
    if(cut) w2.push(`파일이 잘려 ${cut} 면의 끝부분을 0 으로 채웠습니다 (자료 ${blk.avail}/${blk.size} 바이트)`);
    const nm=plan.length>1?`${name} · ${P.s.name}${P.nM>1?` · tile ${P.m+1}`:""}`:name;
    const common={name:nm,w:W,h:H,nc,nz:zs.length,nt:ts.length,order:"czt",planes,cal,calCandidates:[],info,warn:w2,zStep};
    if(rgb){ if(bit===8) images.push({...common,bit:8,rgb:true}); else images.push({...common,bit,rgb:false,chNames:["R","G","B"],chLut:["Red","Green","Blue"]}); }
    else { const names=chanNames(P.s,C); const luts=C.map(c=>LUTN[(c.lut||"").toLowerCase()]||null);
      images.push({...common,bit,rgb:false,chNames:names,chLut:luts.some(Boolean)?luts:null}); }
  }
  if(!images.length) throw new Error("LIF: 읽은 이미지가 없습니다"+(warn.length?" — "+warn.join(", "):""));
  { const fst=plan[sel.p[0]].s.el; const t=xml.slice(fst.s0,fst.s1); images[0].metaText=t.length>200000?t.slice(0,200000)+"…":t; } // 첫 이미지 요소의 XML 원문 (저장 공간 — 첫 이미지에만)
  return {images,warn};
}
/* 채널 이름: ChannelProperty DyeName → (모두 없으면) HardwareSetting 의 MultiBand DyeName(채널 수와 같을 때) → 순차 주사 단계마다 하나씩 → C1… (liffile 과 같은 차례) */
function chanNames(s,C){ const n=C.length; const strip=v=>v.replace(/^Leica\//,"");
  const nm=C.map(c=>{ for(const p of kids(c.el,"ChannelProperty")){ const k=kid(p,"Key"), v=kid(p,"Value"); if(k&&k.t==="DyeName"&&v&&v.t) return strip(v.t); } return ""; });
  if(nm.some(Boolean)) return nm.map((v,k)=>v||`C${k+1}`);
  const hs=kids(s.img,"Attachment").find(a=>a.a.Name==="HardwareSetting");
  if(hs){ const byInc=C.map((c,k)=>k).sort((a,b)=>C[a].inc-C[b].inc); const put=list=>{ const r=new Array(n); byInc.forEach((k,i)=>r[k]=list[i]); return r; };
    const all=(e,t,out)=>{ for(const c of e.c){ if(c.tag===t) out.push(c); all(c,t,out); } return out; };
    const def=kid(hs,"ATLConfocalSettingDefinition"); if(def){ const mb=all(def,"MultiBand",[]).map(m=>m.a.DyeName||"").filter(Boolean); if(mb.length===n) return put(mb.map(strip)); }
    const seq=kid(hs,"LDM_Block_Sequential"); if(seq){ const steps=[]; for(const L of kids(seq,"LDM_Block_Sequential_List")) for(const d of kids(L,"ATLConfocalSettingDefinition")){ const mb=all(d,"MultiBand",[]).map(m=>m.a.DyeName||"").filter(Boolean); if(mb.length===1) steps.push(strip(mb[0])); }
      if(steps.length===n) return put(steps); } }
  return C.map((_,k)=>`C${k+1}`); }
IO.readLIF=readLIF; IO._xmlTree=xmlTree;
})();
