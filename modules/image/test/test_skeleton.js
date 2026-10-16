// Analyze Skeleton — 앱 구현(src/core/85_skeleton.js: 동작 명세 cleanroom/SPEC_analyze_skeleton.md 와 블랙박스 대조만으로 새로 쓴 독립 구현)을
// Fiji AnalyzeSkeleton_ (GPL-3, 앱에 넣지 않음 — 시험에서 실행만) 의 결과와 대조
// 기준: ref/askel/SkelDriver.java + github.com/fiji/AnalyzeSkeleton 의 6 파일을 ImageJ 1.54u jar 와 컴파일 (/tmp/askel_out)
//   git clone --depth 1 https://github.com/fiji/AnalyzeSkeleton /tmp/AnalyzeSkeleton
//   S=/tmp/AnalyzeSkeleton/src/main/java/sc/fiji/analyzeSkeleton; javac -cp /tmp/ij.jar -d /tmp/askel_out ref/askel/SkelDriver.java $S/{AnalyzeSkeleton_,Edge,Vertex,Graph,SkeletonResult,Point}.java
// 사용: node test/test_skeleton.js [ij.jar] [askel_out]
//   SKEL_IMPL=<모듈 경로>  : 그 모듈의 analyzeSkeleton(vol,cal,opt) 을 시험 (없으면 앱 core 의 IC.analyzeSkeleton). 모듈이 PRUNE=true 를 내보내면 가지치기 사례도
//   FUZZ=<n>               : 무작위 사례 n 개 더 (2D·3D·골격·잡음)   SEED=<수>: 무작위 씨앗
const fs=require("fs"), cp=require("child_process"), path=require("path"), os=require("os");
const IC=require("./load_core.js"); const IO=require("../src/img_io.js");
const IMPL=process.env.SKEL_IMPL?require(path.resolve(process.env.SKEL_IMPL)):{analyzeSkeleton:IC.analyzeSkeleton,PRUNE:!!IC.SKELETON_PRUNE};
const RUN=IMPL.analyzeSkeleton;
const JAR=process.argv[2]||"/tmp/ij.jar", OUT=process.argv[3]||"/tmp/askel_out";
if(!fs.existsSync(JAR)||!fs.existsSync(path.join(OUT,"SkelDriver.class"))){ console.log("ImageJ jar 또는 AnalyzeSkeleton 기준 없음 — 건너뜀"); process.exit(0); }
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(c||!d?"":" — "+d)); };
const TMP=fs.mkdtempSync(path.join(os.tmpdir(),"askel-"));
let seed=(+(process.env.SEED||97531))>>>0; const rnd=()=>{ seed=(seed+0x6D2B79F5)>>>0; let t=seed; t=Math.imul(t^(t>>>15),t|1); t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; }; // mulberry32 (이전 LCG 는 double 정밀도를 잃어 짧게 되풀이)
const cases=[]; // {name,w,h,d,data(Uint8 0/255),cal:[pw,ph,pd],lsp}
const blank=(w,h,d)=>({w,h,d,data:new Uint8Array(w*h*d)});
const put=(v,x,y,z)=>{ z=z||0; if(x>=0&&x<v.w&&y>=0&&y<v.h&&z>=0&&z<v.d) v.data[x+v.w*(y+v.h*(z||0))]=255; };
function line(v,x0,y0,x1,y1,z){ const n=Math.max(Math.abs(x1-x0),Math.abs(y1-y0),1); for(let k=0;k<=n;k++) put(v,Math.round(x0+(x1-x0)*k/n),Math.round(y0+(y1-y0)*k/n),z||0); }
function add(name,v,cal,lsp,pr){ cases.push(Object.assign({name,cal:cal||[1,1,1],lsp:lsp!==false,prune:0,pruneEnds:false,orig:null},v,pr||{})); }
const PRUNE_NAMES=["none","shortest branch","lowest intensity voxel","lowest intensity branch"];
// 손으로 만든 경우
{ const v=blank(20,15,1); line(v,2,7,16,7); line(v,9,1,9,13); add("cross",v); }
{ const v=blank(9,9,1); put(v,4,4); add("one pixel",v); }
{ const v=blank(9,9,1); put(v,4,4); put(v,5,5); add("two pixels",v); }
{ const v=blank(30,30,1); for(let a=0;a<360;a+=1) put(v,Math.round(15+10*Math.cos(a*Math.PI/180)),Math.round(15+10*Math.sin(a*Math.PI/180))); add("ring (no junction)",v); }
{ const v=blank(22,22,1); line(v,10,2,18,10); line(v,18,10,10,18); line(v,10,18,2,10); line(v,2,10,10,2); add("diamond loop (all slab)",v,[0.7,0.7,1]); }
{ const v=blank(30,22,1); line(v,10,2,18,10); line(v,18,10,10,18); line(v,10,18,2,10); line(v,2,10,10,2); line(v,18,10,28,10); add("loop + tail",v); }
{ const v=blank(40,30,1); for(let a=0;a<360;a+=1){ put(v,Math.round(12+8*Math.cos(a*Math.PI/180)),Math.round(15+8*Math.sin(a*Math.PI/180))); put(v,Math.round(28+8*Math.cos(a*Math.PI/180)),Math.round(15+8*Math.sin(a*Math.PI/180))); } add("figure 8",v); }
{ const v=blank(30,20,1); line(v,0,10,29,10); line(v,15,0,15,10); line(v,5,19,15,10); add("T + edges at border",v); }
{ const v=blank(25,25,1); for(let y=2;y<23;y+=5) line(v,2,y,22,y); for(let x=2;x<23;x+=5) line(v,x,2,x,22); add("grid (loops)",v); }
{ const v=blank(25,25,1); line(v,3,3,20,20); line(v,3,20,20,3); line(v,11,0,11,24); add("star (clustered junction)",v,[0.3,0.5,1]); }
// ImageJ Skeletonize 로 만든 골격 (임의 원·선 덩어리)
for(let k=0;k<14;k++){ const w=60+(k%3)*30, h=50+(k%4)*20; const p={w,h,bit:8,data:new Float32Array(w*h)};
  const nb=3+Math.floor(rnd()*6); for(let b=0;b<nb;b++){ const cx=rnd()*w, cy=rnd()*h, rx=3+rnd()*12, ry=3+rnd()*12; for(let y=0;y<h;y++) for(let x=0;x<w;x++) if(((x-cx)/rx)**2+((y-cy)/ry)**2<=1) p.data[y*w+x]=255; }
  for(let b=0;b<3;b++){ const v={w,h,d:1,data:new Uint8Array(w*h)}; const x0=rnd()*w,y0=rnd()*h,x1=rnd()*w,y1=rnd()*h,th=1+rnd()*3; for(let y=0;y<h;y++) for(let x=0;x<w;x++){ const t=Math.max(0,Math.min(1,((x-x0)*(x1-x0)+(y-y0)*(y1-y0))/((x1-x0)**2+(y1-y0)**2||1))); if(Math.hypot(x-x0-t*(x1-x0),y-y0-t*(y1-y0))<th) p.data[y*w+x]=255; } }
  const sk=IC.skeletonize(p); const v={w,h,d:1,data:Uint8Array.from(sk.data,x=>x?255:0)}; add(`skeletonize #${k}`,v,k%2?[0.25,0.25,1]:[1,1,1]); }
// 아무 이진 잡음 (골격이 아님 — 접점 덩어리·이상한 이웃이 많음)
for(let k=0;k<12;k++){ const w=30+k*3, h=25+k*2, dens=0.08+0.03*(k%6); const v=blank(w,h,1); for(let i=0;i<w*h;i++) if(rnd()<dens) v.data[i]=255; add(`noise 2D #${k} (${dens.toFixed(2)})`,v,[1,1,1],k%2===0); }
// 3D: 무작위 걸음·잡음, z 간격 다름
for(let k=0;k<8;k++){ const w=24,h=20,d=6+k%4; const v=blank(w,h,d); for(let s=0;s<3+k%3;s++){ let x=Math.floor(rnd()*w),y=Math.floor(rnd()*h),z=Math.floor(rnd()*d); for(let n=0;n<30;n++){ put(v,x,y,z); x+=Math.floor(rnd()*3)-1; y+=Math.floor(rnd()*3)-1; z+=Math.floor(rnd()*3)-1; x=Math.max(0,Math.min(w-1,x)); y=Math.max(0,Math.min(h-1,y)); z=Math.max(0,Math.min(d-1,z)); } }
  add(`3D walk #${k}`,v,[0.5,0.5,1.7]); }
for(let k=0;k<6;k++){ const w=18,h=16,d=5; const v=blank(w,h,d); for(let i=0;i<w*h*d;i++) if(rnd()<0.05+0.02*k) v.data[i]=255; add(`noise 3D #${k}`,v,[1,1,2]); }

// 무작위 사례 (FUZZ)
const NF=+(process.env.FUZZ||0);
for(let k=0;k<NF;k++){ const kind=k%5;
  if(kind<=1){ const w=20+Math.floor(rnd()*50), h=15+Math.floor(rnd()*40); const p={w,h,bit:8,data:new Float32Array(w*h)}; const nb=1+Math.floor(rnd()*5);
    for(let b=0;b<nb;b++){ const cx=rnd()*w, cy=rnd()*h, rx=2+rnd()*9, ry=2+rnd()*9; for(let y=0;y<h;y++) for(let x=0;x<w;x++) if(((x-cx)/rx)**2+((y-cy)/ry)**2<=1) p.data[y*w+x]=255; }
    if(kind===1) for(let i=0;i<w*h;i++) if(rnd()<0.03) p.data[i]=p.data[i]?0:255;
    const sk=IC.skeletonize(p); add(`fuzz skel #${k}`,{w,h,d:1,data:Uint8Array.from(sk.data,x=>x?255:0)},[0.3+rnd(),0.3+rnd(),1]); }
  else if(kind===2){ const w=8+Math.floor(rnd()*30), h=6+Math.floor(rnd()*25); const v=blank(w,h,1); const dens=0.05+rnd()*0.4; for(let i=0;i<w*h;i++) if(rnd()<dens) v.data[i]=255; add(`fuzz noise2D #${k}`,v,[1,1+rnd(),1],rnd()<0.7); }
  else if(kind===3){ const w=6+Math.floor(rnd()*14),h=6+Math.floor(rnd()*12),d=2+Math.floor(rnd()*6); const v=blank(w,h,d); const dens=0.03+rnd()*0.15; for(let i=0;i<w*h*d;i++) if(rnd()<dens) v.data[i]=255; add(`fuzz noise3D #${k}`,v,[0.5+rnd(),0.5+rnd(),0.5+2*rnd()]); }
  else { const w=16+Math.floor(rnd()*14),h=14+Math.floor(rnd()*12),d=3+Math.floor(rnd()*6); const v=blank(w,h,d); for(let s2=0;s2<2+Math.floor(rnd()*4);s2++){ let x=Math.floor(rnd()*w),y=Math.floor(rnd()*h),z=Math.floor(rnd()*d); for(let n=0;n<15+rnd()*40;n++){ put(v,x,y,z); x=Math.max(0,Math.min(w-1,x+Math.floor(rnd()*3)-1)); y=Math.max(0,Math.min(h-1,y+Math.floor(rnd()*3)-1)); z=Math.max(0,Math.min(d-1,z+Math.floor(rnd()*3)-1)); } } add(`fuzz walk3D #${k}`,v,[1,1,1+rnd()]); } }
// 가지치기 (구현이 PRUNE 을 내보낼 때): 고리 있는 골격·잡음 + 원본 세기(매끈한 무작위)
if(IMPL.PRUNE){ const base=cases.slice(); let j=0;
  for(const c of base){ if(c.d>1&&j%3) { j++; continue; } const orig=new Uint8Array(c.w*c.h*c.d); const fx=rnd()*0.7, fy=rnd()*0.7, ph=rnd()*6; for(let i=0;i<orig.length;i++){ const x=i%c.w, y=((i/c.w)|0)%c.h, z=(i/(c.w*c.h))|0; orig[i]=Math.max(0,Math.min(255,Math.round(128+60*Math.sin(fx*x+ph)+50*Math.cos(fy*y+z)+20*(rnd()-0.5)))); }
    const pr=1+(j%3), pe=j%2===0; add(c.name+` · prune ${PRUNE_NAMES[pr]}${pe?" + ends":""}`,{w:c.w,h:c.h,d:c.d,data:c.data},c.cal,c.lsp,{prune:pr,pruneEnds:pe,orig});
    if(j%4===0) add(c.name+" · prune ends only",{w:c.w,h:c.h,d:c.d,data:c.data},c.cal,c.lsp,{prune:0,pruneEnds:true,orig:null}); j++; } }

// Fiji 실행
const list=[]; cases.forEach((c,i)=>{ const f=path.join(TMP,`c${i}.tif`); const planes=[]; for(let z=0;z<c.d;z++) planes.push(Float32Array.from(c.data.subarray(z*c.w*c.h,(z+1)*c.w*c.h)));
  fs.writeFileSync(f,IO.writeTIFF({w:c.w,h:c.h,bit:8,rgb:false,nc:1,nz:c.d,nt:1,planes,cal:null}));
  let of=""; if(c.orig){ of=path.join(TMP,`o${i}.tif`); const op=[]; for(let z=0;z<c.d;z++) op.push(Float32Array.from(c.orig.subarray(z*c.w*c.h,(z+1)*c.w*c.h))); fs.writeFileSync(of,IO.writeTIFF({w:c.w,h:c.h,bit:8,rgb:false,nc:1,nz:c.d,nt:1,planes:op,cal:null})); }
  list.push([f,c.prune,c.pruneEnds?1:0,c.lsp?1:0,...c.cal].concat(of?[of]:[]).join("\t")); });
fs.writeFileSync(path.join(TMP,"list.txt"),list.join("\n"));
const outTxt=cp.execSync(`JAVA_TOOL_OPTIONS= java -Xmx2g -cp ${JAR}:${OUT} SkelDriver -list ${path.join(TMP,"list.txt")} 2>/dev/null`,{encoding:"utf8",maxBuffer:1<<28});
const blocks=outTxt.split(/^CASE \d+\n/m).slice(1);
console.log(`== Analyze Skeleton: ${cases.length} 경우 (Fiji AnalyzeSkeleton 실행 결과와 대조)`);
const close=(a,b)=>a===b||Math.abs(a-b)<=1e-12*Math.max(1,Math.abs(b));
cases.forEach((c,i)=>{ const L=(blocks[i]||"").trim().split("\n"); const T=L.filter(l=>l.startsWith("T ")).map(l=>l.split(" ").slice(2).map(Number)); const E=L.filter(l=>l.startsWith("E ")).map(l=>l.split(" ").slice(1).map(Number));
  const TAG=new Map(((L.find(l=>l.startsWith("TAG"))||"TAG").split(" ").slice(1)).map(s=>s.split(":").map(Number)));
  const LSPI=new Set(((L.find(l=>l.startsWith("LSPIMG"))||"LSPIMG").split(" ").slice(1)).map(Number));
  let r, ex=""; try{ r=RUN({w:c.w,h:c.h,d:c.d,data:c.data},{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]},Object.assign({lsp:c.lsp},c.prune||c.pruneEnds?{prune:PRUNE_NAMES[c.prune],pruneEnds:c.pruneEnds,orig:c.orig?{data:c.orig}:null}:{})); }catch(e){ ok(false,c.name,"오류: "+e.message); return; }
  if(r.trees.length!==T.length) ex+=` 나무 수 앱 ${r.trees.length} Fiji ${T.length}`;
  else r.trees.forEach((t,k)=>{ const a=[t.branches,t.junctions,t.endPoints,t.junctionVoxels,t.slabs,t.avgLength,t.triples,t.quadruples,t.maxLength].concat(c.lsp?[t.lsp,t.spx,t.spy,t.spz]:[]); const b=T[k];
    const names=["branches","junctions","endPoints","junctionVoxels","slabs","avgLength","triples","quadruples","maxLength","lsp","spx","spy","spz"]; a.forEach((v,j)=>{ if(!close(v,b[j])&&!ex) ex+=` 나무 ${k} ${names[j]} 앱 ${v} Fiji ${b[j]}`; }); });
  if(!ex&&r.edges.length!==E.length) ex+=` 가지 수 앱 ${r.edges.length} Fiji ${E.length}`;
  if(!ex) r.edges.forEach((e,k)=>{ const a=[e.tree,e.length,...e.v1,...e.v2,e.lengthRa,e.nSlabs], b=E[k]; a.forEach((v,j)=>{ if(!close(v,b[j])&&!ex) ex+=` 가지 ${k} 항목 ${j} 앱 ${v} Fiji ${b[j]}`; }); });
  if(!ex){ let bad=0; for(let p=0;p<r.tag.length;p++){ const t=TAG.get(p)||0; if(r.tag[p]!==t) bad++; } if(bad) ex+=` 태그 다른 화소 ${bad}`; }
  if(!ex&&c.lsp){ let bad=0; const li=r.lspImage; for(let p=0;p<c.w*c.h*c.d;p++){ if(((li&&li[p]===96)?1:0)!==(LSPI.has(p)?1:0)) bad++; } if(bad) ex+=` LSP 영상 다른 화소 ${bad} (Fiji 경로 ${LSPI.size})`; }
  ok(!ex,`${c.name} (${c.w}×${c.h}${c.d>1?"×"+c.d:""}, 나무 ${T.length}, 가지 ${E.length}${c.lsp?", LSP":""})`,ex); });
console.log(`\n${pass} 통과, ${fail} 실패`); try{ fs.rmSync(TMP,{recursive:true}); }catch(e){} if(fail){ console.log(fails.join("\n")); process.exit(1); }
