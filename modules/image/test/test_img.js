/* Image 모듈 검증 — node test/test_img.js  (modules/image 에서)
   1) ImageJ 1.54u 실행 결과(ref/data, ref/gen_ref.py 가 만듦)와 화소·수치 대조
   2) 정답을 아는 합성 영상으로 성질 시험 (EDM, 입자 크기, 방향, 이동 추정 등)
   3) 파일 읽기 회귀 (img_io)                                                        */
const fs=require("fs"), path=require("path");
const IC=require("./load_core.js");
const REF=path.join(__dirname,"..","ref","data");
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!cond) console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
const rdF32=(p)=>{ const b=fs.readFileSync(p); return new Float32Array(b.buffer,b.byteOffset,b.byteLength/4).slice(); };

function compareImg(name,got,ref,bit){
  const n=ref.length; if(got.length!==n) return ok(false,name,`크기 ${got.length} ≠ ${n}`);
  let eq=0, mx=0, mi=-1; let lo=Infinity, hi=-Infinity;
  for(let i=0;i<n;i++){ const a=got[i], b=ref[i]; if(b<lo) lo=b; if(b>hi) hi=b; const d=(a!==a&&b!==b)?0:Math.abs(a-b); if(d===0) eq++; if(d>mx||d!==d){ mx=d; mi=i; } }
  if(bit===32){ const rg=Math.max(1e-12,hi-lo); const rel=mx/rg; ok(rel<2e-5,name,`최대차/범위 ${rel.toExponential(2)} (일치 ${(100*eq/n).toFixed(2)} %)`+(rel>=2e-5?` @${mi} got ${got[mi]} ref ${ref[mi]}`:"")); }
  else { const fr=eq/n; const need=bit===16?0.99:0.999; // 16-bit: ImageJ 의 float32 중간 계산과 우리 float64 의 차이로 ±1 (전체 눈금의 1.5e-5)
    ok(fr>=need&&mx<=1,name,`일치 ${(100*fr).toFixed(3)} %, 최대차 ${mx}`+(mx>1?` @${mi} (${mi%W0},${Math.floor(mi/W0)}) got ${got[mi]} ref ${ref[mi]}`:"")); }
}
let W0=0;

(async()=>{
const HAVE_REF=fs.existsSync(path.join(REF,"manifest.json"));
if(!HAVE_REF) console.log("!! ref/data 가 없어 1 (ImageJ 대조) 을 건너뜁니다 — ImageJ 를 컴파일한 ij.jar 로 python3 -I ref/gen_ref.py ij.jar ref/data (용량 59 MB 라 저장소에 두지 않음)");
const MAN=HAVE_REF?JSON.parse(fs.readFileSync(path.join(REF,"manifest.json"),"utf8")):{meta:{},cases:[]};
const IN={}; for(const k in MAN.meta){ const m=MAN.meta[k]; const d=rdF32(path.join(REF,"in",k+".f32")); IN[k]=m.n?{...m,data:d}:{w:m.w,h:m.h,bit:m.bit,data:d}; }
const P=k=>({w:IN[k].w,h:IN[k].h,bit:IN[k].bit,data:Float32Array.from(IN[k].data)});
const planesOf=k=>{ const m=IN[k], N=m.w*m.h; return [...Array(m.n)].map((_,i)=>({w:m.w,h:m.h,bit:m.bit,data:m.data.slice(i*N,(i+1)*N)})); };
const RANGE={g8:null,g16:null,g32:null};
const T={ g8:{thrO:20,proms:[10,30]}, g16:{thrO:3000,proms:[2000,6000]}, g32:{thrO:0.3,proms:[0.2,0.6]} };

console.log("== 1 ImageJ 1.54u 대조 ("+MAN.cases.length+" 사례)");
for(const c of MAN.cases){
  const nm=c.name; let s=c.src; W0=MAN.meta[s]?MAN.meta[s].w:0;
  try{
    // ----- 영상 결과 -----
    if(c.kind==="img"||c.kind==="imgnew"){
      const ref=rdF32(path.join(REF,"out",nm+".f32")); let got=null; const p=s==="stk"?null:P(s); let m;
      if((m=nm.match(/^gauss_(g\d+)_([\d.]+)$/))) got=IC.gaussian(p,{sigma:+m[2]});
      else if((m=nm.match(/^(median|mean|minimum|maximum|variance)_(g\d+)_([\d.]+)$/))) got=IC.rank(p,{median:"median",mean:"mean",minimum:"min",maximum:"max",variance:"variance"}[m[1]],{radius:+m[3]});
      else if((m=nm.match(/^out([BD])_(g\d+)$/))) got=IC.rank(p,"outliers",{radius:2,threshold:T[s].thrO,which:m[1]==="B"?"bright":"dark"});
      else if((m=nm.match(/^tophat(L?)_/))) got=IC.rank(p,"tophat",{radius:5,light:!!m[1]});
      else if(nm.startsWith("despeckle_")) got=IC.rank(p,"despeckle",{});
      else if(nm.startsWith("smooth_")) got=IC.filter3(p,"smooth"); else if(nm.startsWith("sharpen_")) got=IC.filter3(p,"sharpen"); else if(nm.startsWith("findedges_")) got=IC.filter3(p,"edges");
      else if(nm.startsWith("unsharp_")) got=IC.unsharp(p,{sigma:2,weight:0.6});
      else if((m=nm.match(/^rb_(g\d+)_(\d+)$/))) got=IC.subtractBackground(p,{radius:+m[2]});
      else if(nm.startsWith("rbL_")) got=IC.subtractBackground(p,{radius:15,light:true});
      else if(nm.startsWith("rbNS_")) got=IC.subtractBackground(p,{radius:8,smooth:false});
      else if(nm.startsWith("rbC_")) got=IC.subtractBackground(p,{radius:25,createBackground:true});
      else if(nm.startsWith("conv_")) got=IC.convolve(p,{kw:3,kh:3,k:[1,2,1,2,4,2,1,2,1],normalize:true});
      else if(nm.startsWith("bpH_")) got=IC.bandpass(p,{large:25,small:2,stripes:"horizontal"});
      else if(nm.startsWith("bpN_")) got=IC.bandpass(p,{large:30,small:4,autoscale:false,saturate:false});
      else if(nm.startsWith("bp_")) got=IC.bandpass(p,{large:40,small:3});
      else if((m=nm.match(/^math_(\w+)_(g\d+)$/))){ const op={log:["log"],exp:["exp"],squareroot:["sqrt"],square:["sqr"],gamma:["gamma",0.6],multiply:["mul",1.3],add:["add",25],subtract:["sub",17.6],divide:["div",2.5],min:["min",60],max:["max",150]}[m[1]]; got=IC.math(p,op[0],op[1]); }
      else if(nm.startsWith("enh_")) got=IC.enhanceContrast(p,{saturated:0.35,normalize:true}).plane;
      else if(nm.startsWith("eqC_")) got=IC.equalize(p,true); else if(nm.startsWith("eq_")) got=IC.equalize(p,false);
      else if(nm.startsWith("inv_")) got=IC.invert(p);
      else if((m=nm.match(/^max_(g\d+)_([\d.]+)_(single|within|seg)$/))) got=IC.findMaximaImage(p,+m[2],{output:{single:"single",within:"within",seg:"segmented"}[m[3]]});
      else if(nm.startsWith("to8_")){ const [a,b]=IC.minMax(p.data); got=IC.convertBit(p,8,a,b,true); }
      else if(nm.startsWith("to16_")){ const [a,b]=IC.minMax(p.data); got=IC.convertBit(p,16,a,b,true); }
      else if((m=nm.match(/^bin_(erode|dilate|open|close)_(\d)_(\d)$/))) got=IC.binaryMorph(p,m[1],{iterations:+m[2],count:+m[3]});
      else if(nm==="bin_outlineW") got=IC.outline(p); else if(nm==="bin_fillholes") got=IC.fillHoles(p); else if(nm==="bin_skeletonize") got=IC.skeletonize(p);
      else if(nm==="bin_watershed") got=IC.watershed(p); else if(nm==="bin_ultimatepoints") got=IC.ultimatePoints(p);
      else if(nm==="bin_edm8") got=IC.distanceMap(p,"8-bit"); else if(nm==="bin_edm32") got=IC.distanceMap(p,"32-bit");
      else if(nm==="bin_g8_mask"){ const t=IC.autoThreshold(p,"Otsu",true); got=IC.applyThreshold(p,t.lower,t.upper); }
      else if((m=nm.match(/^calc_(\w+)$/))){ const b=P("bin8"); got=IC.calc(p,b,{add:"add",sub32:"sub",diff:"difference",avg:"average",div:"div",mul:"mul"}[m[1]],m[1]==="sub32"); }
      else if((m=nm.match(/^zp_(\w+)$/))) got=IC.zProject(planesOf("stk"),m[1]);
      if(!got){ ok(false,nm,"대응 함수 없음"); continue; }
      const bit=c.bit||got.bit;
      if(got.bit!==c.bit&&!(c.bit===16&&got.bit===16)) { ok(false,nm,`bit ${got.bit} ≠ ImageJ ${c.bit}`); continue; }
      W0=got.w; compareImg(nm,got.data,ref,bit);
      continue;
    }
    if(c.kind==="txt"){ // 문턱
      const p=P(s); let worst=0, bad=[];
      for(const k in c.thr){ const [meth,mode]=k.split("|"); const r=IC.autoThreshold(p,meth,mode==="dark"); const [l,u]=c.thr[k];
        const tol=p.bit===32?Math.abs(r.max-r.min)*1e-6+1e-9:1e-6; const d=Math.max(Math.abs(r.lower-l),Math.abs(r.upper-u)); if(d>tol) bad.push(`${k}: [${r.lower},${r.upper}] vs [${l},${u}]`); worst=Math.max(worst,d); }
      ok(bad.length===0,`${nm}: 16 방법 × dark/light 문턱 일치`,bad.slice(0,3).join("; "));
      continue;
    }
    if(c.kind==="csv"){
      let m;
      if((m=nm.match(/^max_(g\d+)_([\d.]+)_(list|listX|listS|listL|count)$/))){
        if(m[3]==="count") continue;
        const p=P(s); const r=IC.findMaxima(p,+m[2],{output:"points",excludeEdges:m[3]==="listX",strict:m[3]==="listS",light:m[3]==="listL"});
        const A=new Set(r.points.map(q=>q[0]+","+q[1])), B=new Set(c.rows.map(q=>q.X+","+q.Y)); let miss=0; for(const k of B) if(!A.has(k)) miss++; let extra=0; for(const k of A) if(!B.has(k)) extra++;
        ok(miss===0&&extra===0,`${nm}: 극대 위치 일치`,`ImageJ ${B.size}, 우리 ${A.size}, 빠짐 ${miss}, 더함 ${extra}`); continue; }
      if(nm.startsWith("pa")){
        const b=P("bin8"); const bin=Uint8Array.from(b.data,v=>v===255?1:0); let o={}; let cal={pw:1,ph:1}; let inten=b;
        if(nm==="paX") o={excludeEdges:true}; else if(nm==="paH") o={includeHoles:true}; else if(nm==="paS") o={minSize:20,maxSize:400}; else if(nm==="paC") o={minCirc:0.6,maxCirc:1}; else if(nm==="paXH") o={minSize:5,excludeEdges:true,includeHoles:true};
        else if(nm==="paCal") cal={pw:2.5,ph:2.5}; else if(nm==="paRed") inten=P("g8");
        if(nm==="paCal"){ o={}; }
        const r=IC.analyzeParticles(bin,b.w,b.h,inten,cal,o);
        comparePA(nm,r.particles.map(q=>q.meas),c.rows); continue; }
      if((m=nm.match(/^(roi\w+)_(g\d+)$/))){
        const p=P(m[2]); const roi=makeRoi(m[1],p.w,p.h); const got=IC.measureRoi(p,roi,{pw:1,ph:1});
        const skip=m[1].startsWith("roiRectOut")?["X","Y","XM","YM","Major","Minor","Angle","AR","Round"]:m[1].startsWith("roiOval")?["Feret","FeretX","FeretY","FeretAngle","MinFeret","Solidity"]:m[1]==="roiFree"?["X","Y","Feret","FeretX","FeretY","FeretAngle","MinFeret","Solidity","BX","BY","Width","Height"]:[];
        comparePA(nm,[got],c.rows,skip); continue; }
      ok(false,nm,"대응 시험 없음"); continue;
    }
  }catch(e){ ok(false,nm,"예외: "+(e.stack||e).toString().split("\n").slice(0,3).join(" | ")); }
}
function makeRoi(kind,W,H){
  if(kind==="roiRect"){ const r={x:31,y:22,w:57,h:41}; return {mask:null,rect:r,poly:null,kind:"rect"}; }
  if(kind==="roiOval"){ const m=IC.ovalMask(40,30,61,45,W,H); return {mask:m.mask,rect:m.rect,poly:null,kind:"oval"}; }
  if(kind==="roiOval2"){ const m=IC.ovalMask(20,20,50,50,W,H); return {mask:m.mask,rect:m.rect,poly:null,kind:"oval"}; }
  const geo=(type,X0,Y0,BW,BH)=>{ if(type==="rect"){ const x=Math.max(0,X0), y=Math.max(0,Y0), w=Math.min(W,X0+BW)-x, h=Math.min(H,Y0+BH)-y; return {mask:null,rect:{x,y,w,h},bounds:{x:X0,y:Y0,w:BW,h:BH},poly:{xs:[X0,X0+BW,X0+BW,X0],ys:[Y0,Y0,Y0+BH,Y0+BH]},kind:"rect"}; }
    const m=IC.ovalMask(X0,Y0,BW,BH,W,H); return {mask:m.mask,rect:m.rect,bounds:{x:X0,y:Y0,w:BW,h:BH},poly:null,kind:"oval"}; };
  if(kind==="roiOvalOut") return geo("oval",-10,5,30,20); if(kind==="roiRectOut") return geo("rect",-10,5,30,20); if(kind==="roiOvalOut2") return geo("oval",220,150,40,50);
  if(kind==="roiPoly"){ const xs=[30,120,140,60,20], ys=[20,35,110,150,90]; const m=IC.polygonMask(xs,ys,W,H); return {mask:m.mask,rect:m.rect,poly:{xs,ys},kind:"polygon"}; }
  if(kind==="roiFree"){ const xs=[30,40,55,70,90,95,85,60,40,28], ys=[30,22,25,20,35,60,90,100,85,60]; const m=IC.polygonMask(xs,ys,W,H); return {mask:m.mask,rect:m.rect,poly:{xs,ys},kind:"freehand"}; }
}
function comparePA(nm,got,rows,skip){
  skip=new Set(skip||[]);
  if(got.length!==rows.length){ ok(false,nm+": 입자 수",`우리 ${got.length} vs ImageJ ${rows.length}`); return; }
  const cols=Object.keys(rows[0]||{}).filter(k=>k&&k!==""&&!skip.has(k)); let worst=0, wk="", wi=-1;
  for(let i=0;i<rows.length;i++) for(const k of cols){ const a=got[i][k], b=rows[i][k]; if(b==null&&(a==null||!isFinite(a))) continue; if(b==null||a==null){ if(!(b==null&&a!=null&&!isFinite(a))){ worst=Infinity; wk=k; wi=i; } continue; }
    const d=Math.abs(a-b)/Math.max(1,Math.abs(b)); if(d>worst){ worst=d; wk=k; wi=i; } }
  ok(worst<1e-6,`${nm}: ${rows.length}개 × ${cols.length}열 일치`,`최대 상대차 ${worst.toExponential(2)} (${wk}, #${wi+1}: ${wi>=0?got[wi][wk]:""} vs ${wi>=0?rows[wi][wk]:""})`);
}

console.log("== 2 정답을 아는 합성 영상 (성질 시험)");
{ const W=200,H=160; const mk=f=>{ const d=new Float32Array(W*H); for(let y=0;y<H;y++) for(let x=0;x<W;x++) d[y*W+x]=f(x,y); return {w:W,h:H,bit:32,data:d}; };
  // Orientation: 화면에서 반시계 θ 방향으로 뻗은 줄무늬
  for(const th of [0,30,-45,72,-12.5]){ const t=th*Math.PI/180; const p=mk((x,y)=>100+50*Math.sin(2*Math.PI*(x*Math.sin(t)+y*Math.cos(t))/9));
    const r=IC.orientation(p,{sigma:3}); let d=r.dominant-th; d=((d+90)%180+180)%180-90; ok(Math.abs(d)<0.3,`Orientation: 줄무늬 ${th}° → ${r.dominant.toFixed(2)}°`,`coherency ${r.meanCoherency.toFixed(3)}`); }
  // 위상 상관: 알려진 부화소 이동을 되돌린다
  { const F=(x,y)=>Math.exp(-((x-80)**2+(y-70)**2)/200)*200+Math.exp(-((x-130)**2+(y-50)**2)/90)*150+30*Math.sin(x/7)*Math.cos(y/11)+100; const a=mk(F);
    // b = 연속 함수를 (dx,dy) 만큼 옮겨 다시 표본화 (보간 인공물 없는 참 이동)
    for(const [dx,dy] of [[3,2],[-5.4,2.3],[0.6,-7.2],[0.25,0.5]]){ const b=mk((x,y)=>F(x-dx,y-dy)); const s=IC.phaseCorr(a,b); ok(Math.abs(s.dx+dx)<0.03&&Math.abs(s.dy+dy)<0.03,`Phase correlation: (${dx}, ${dy}) 이동 → 추정 (${s.dx.toFixed(2)}, ${s.dy.toFixed(2)})`); } }
  // Radial profile: 반지름 30 의 고리
  { const p=mk((x,y)=>Math.exp(-((Math.hypot(x+0.5-100,y+0.5-80)-30)**2)/4)); const r=IC.radialProfile(p,100,80,70,1); let k=0; r.mean.forEach((v,i)=>{ if(v>r.mean[k]) k=i; }); ok(Math.abs(r.r[k]-30)<=1,`Radial profile: 고리 r=30 → 극대 ${r.r[k]}`); }
  // EDM: 완전 탐색과 비교 (ImageJ 알고리즘은 근사 — 최대 오차 확인)
  { const bw=60,bh=50; const b=new Float32Array(bw*bh); const rnd=(()=>{ let z=7; return ()=>(z=(z*1103515245+12345)%2147483648)/2147483648; })(); for(let i=0;i<bw*bh;i++) b[i]=rnd()<0.93?255:0;
    const e=IC.makeFloatEDM({w:bw,h:bh,bit:8,data:b},0,false).data; let mx=0; const bg=[]; for(let i=0;i<bw*bh;i++) if(!b[i]) bg.push([i%bw,Math.floor(i/bw)]);
    for(let i=0;i<bw*bh;i++){ if(!b[i]) continue; const x=i%bw,y=Math.floor(i/bw); let d=Infinity; for(const [u,v] of bg) d=Math.min(d,Math.hypot(x-u,y-v)); mx=Math.max(mx,Math.abs(e[i]-d)); }
    ok(mx<0.05,"EDM (ImageJ 알고리즘) ≈ 정확한 유클리드 거리",`최대 오차 ${mx.toFixed(4)} px`); }
  // 입자: 알려진 반지름의 원 → ECD, 원형도
  { const bw=300,bh=200, b=new Uint8Array(bw*bh); const C=[[50,50,10],[150,60,20],[240,140,30.5],[80,150,6]]; for(let y=0;y<bh;y++) for(let x=0;x<bw;x++) for(const [cx,cy,r] of C) if((x+0.5-cx)**2+(y+0.5-cy)**2<=r*r) b[y*bw+x]=1;
    const r=IC.analyzeParticles(b,bw,bh,null,{pw:2,ph:2},{}); const ecd=r.particles.map(q=>q.meas.ECD); const exp=C.map(c=>4*c[2]);
    ok(r.count===4&&ecd.every((v,i)=>Math.abs(v-exp[i])/exp[i]<0.02),"Analyze Particles: 원 4 개 ECD (보정 2 nm/px) 오차 < 2 %",ecd.map(v=>v.toFixed(2)).join(", ")+" vs "+exp.join(", "));
    const circ=r.particles.map(q=>q.meas["Circ."]); ok(circ.slice(1).every(c=>c>0.88),"원의 Circularity > 0.88 (ImageJ 꺾임 보정 둘레)",circ.map(c=>c.toFixed(3)).join(", ")); }
  // 로그정규 크기 분포 회복
  { let z=12345; const rn=()=>{ z=(z*16807)%2147483647; return z/2147483647; }; const v=[]; for(let i=0;i<4000;i++){ const u1=rn(), u2=rn(); v.push(Math.exp(Math.log(20)+0.3*Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2))); }
    const s=IC.sizeSummary(v); ok(Math.abs(s.gMean-20)<0.5&&Math.abs(s.lnSigma-0.3)<0.015,`Log-normal (x̄g=20, σ=0.3, N=4000) → x̄g ${s.gMean.toFixed(3)}, σ ${s.lnSigma.toFixed(4)}`); }
  // Local Mean threshold = 원형 창 평균을 따로 계산한 것
  { const bw=50,bh=40, d=new Float32Array(bw*bh); let z=3; for(let i=0;i<d.length;i++){ z=(z*69069+1)%4294967296; d[i]=Math.floor(z/4294967296*256); } const p={w:bw,h:bh,bit:8,data:d};
    const m=IC.localThreshold(p,"Mean",{radius:4,p1:0}).data; const L=IC.lineRadii(4); let bad=0;
    for(let y=0;y<bh;y++) for(let x=0;x<bw;x++){ let s=0,n=0; for(let k=-L.kR;k<=L.kR;k++) for(let u=-L.dx[k+L.kR];u<=L.dx[k+L.kR];u++){ const xx=Math.min(bw-1,Math.max(0,x+u)), yy=Math.min(bh-1,Math.max(0,y+k)); s+=d[yy*bw+xx]; n++; } const e=d[y*bw+x]>s/n?255:0; if(e!==m[y*bw+x]) bad++; }
    ok(bad===0,"Local threshold (Mean) = 원형 창 평균 직접 계산",`다른 화소 ${bad}`); }
}

console.log("== 3 파일 읽기 (img_io) — tifffile·pypng·numpy·igor2 로 확인한 정답 파일");
{ const IO=require("../src/img_io.js"); require("../src/zstd_vendor.js"); require("../src/img_micro.js"); const zlib=require("zlib"); IO.inflate=async(u8,raw)=>new Uint8Array(raw?zlib.inflateRawSync(u8):zlib.inflateSync(u8));
  const D=path.join(__dirname,"io"); const truth=JSON.parse(fs.readFileSync(path.join(D,"truth.json"),"utf8"));
  const sidecars={}; for(const f of fs.readdirSync(D)) if(/\.txt$/i.test(f)) sidecars[f.toLowerCase()]=fs.readFileSync(path.join(D,f),"utf8");
  const cache={};
  for(const key of Object.keys(truth)){ const t=truth[key]; if(!t) continue; const [file,idx]=key.split("#");
    try{
      if(!cache[file]){ const u8=new Uint8Array(fs.readFileSync(path.join(D,file))); cache[file]=await IO.readFile(file,u8,{sidecars}); }
      const r=cache[file]; const im=r.images[idx?+idx:0]; const tr=rdF32(path.join(D,key+".truth")); const sh=t.shape;
      let planes=im.planes; let flat;
      if(sh.length===3&&sh[2]===3&&!(t.nc)){ flat=new Float32Array(tr.length); const N=sh[0]*sh[1]; for(let i=0;i<N;i++) for(let c=0;c<3;c++) flat[i*3+c]=planes[c][i]; }
      else { flat=new Float32Array(planes.length*planes[0].length); planes.forEach((p,k)=>flat.set(p,k*p.length)); }
      let mx=0; for(let i=0;i<tr.length;i++) mx=Math.max(mx,Math.abs(flat[i]-tr[i])); const same=flat.length===tr.length&&mx<=1e-6*Math.max(1,...[...tr.slice(0,50)].map(Math.abs));
      let extra="";
      if(t.nc&&im.nc!==t.nc) extra+=` nc ${im.nc}≠${t.nc}`; if(t.nz&&im.nz!==t.nz) extra+=` nz ${im.nz}≠${t.nz}`; if(t.nt&&im.nt!==t.nt) extra+=` nt ${im.nt}≠${t.nt}`;
      if(t.pw){ if(!im.cal) extra+=" 보정 없음"; else { const k=IO.unitToM(im.cal.unit)/IO.unitToM(t.unit); if(Math.abs(im.cal.pw*k-t.pw)>1e-6*t.pw) extra+=` 화소 ${im.cal.pw} ${im.cal.unit} ≠ ${t.pw} ${t.unit}`; } }
      if(t.cand){ if(im.cal) extra+=" 해상도 태그를 보정으로 적용함(후보여야 함)"; else if(!im.calCandidates.length||Math.abs(im.calCandidates[0].pw-t.cand)>1e-9) extra+=" 후보 없음"; }
      if(t.databar&&(!im.databar||im.databar.rows!==t.databar)) extra+=" 데이터 막대 미검출";
      if(t.names&&t.names.length&&!t.rgb&&JSON.stringify(im.chNames)!==JSON.stringify(t.names)) extra+=` 채널 이름 ${JSON.stringify(im.chNames)}≠${JSON.stringify(t.names)}`;
      if(t.zstep&&!(Math.abs((im.zStep||0)-t.zstep)<1e-9*t.zstep)) extra+=` z step ${im.zStep}≠${t.zstep}`; if(t.rgb&&!im.rgb) extra+=" RGB 아님";
      ok(same&&!extra,`${key}: 화소${t.pw?" · 보정":""}${t.nc?" · 차원":""} 일치`,`${im.w}×${im.h} ${im.bit}-bit${im.rgb?" RGB":""}, 최대차 ${mx}${extra}`+(im.cal?`, ${im.cal.pw.toPrecision(5)} ${im.cal.unit} (${im.cal.src})`:""));
    }catch(e){ ok(false,key,"예외: "+(e.stack||e).toString().split("\n").slice(0,2).join(" | ")); } }
  // 쓰기 → 다시 읽기 (ImageJ 메타데이터 보정·hyperstack)
  const src={w:20,h:10,bit:16,rgb:false,nc:2,nz:3,nt:1,planes:[...Array(6)].map((_,k)=>Float32Array.from({length:200},(_,i)=>(i*7+k*1000)%65536)),cal:{pw:0.25,ph:0.25,unit:"µm"}};
  const tif=IO.writeTIFF(src); const back=await IO.readFile("w.tif",tif,{});
  const b=back.images[0]; let same=b.planes.length===6&&b.planes.every((p,k)=>p.every((v,i)=>v===src.planes[k][i]));
  ok(same&&b.nc===2&&b.nz===3&&b.cal&&Math.abs(b.cal.pw*IO.unitToM(b.cal.unit)-0.25e-6)<1e-15,"TIFF 쓰기 → 읽기: 화소·hyperstack·보정 보존",`nc ${b.nc} nz ${b.nz} cal ${b.cal&&b.cal.pw} ${b.cal&&b.cal.unit}`);
  fs.writeFileSync(path.join(__dirname,"io","written_by_inno.tif"),tif);
  { const t2=IO.writeTIFF({w:8,h:6,bit:8,nc:1,nz:1,nt:1,planes:[Float32Array.from({length:48},(_,i)=>i)],cal:{pw:0.0001,ph:0.0001,unit:"µm"}}); const r2=(await IO.readFile("s.tif",t2,{})).images[0];
    ok(r2.cal&&Math.abs(r2.cal.pw*IO.unitToM(r2.cal.unit)-1e-10)<1e-14,"TIFF 쓰기: 아주 작은 화소(0.0001 µm) 도 XResolution 넘침 없이",r2.cal&&`${r2.cal.pw} ${r2.cal.unit}`);
    // RowsPerStrip = 2³²−1 (규격 기본값)
    const T=IO._tiffIFDs(t2); const e=T.pages[0][278]; const dv=new DataView(t2.buffer,t2.byteOffset); const ifd=dv.getUint32(4,true); const n=dv.getUint16(ifd,true); for(let i=0;i<n;i++){ const o=ifd+2+12*i; if(dv.getUint16(o,true)===278){ dv.setUint32(o+8,4294967295,true); } } void e;
    let okR=false; try{ const r3=(await IO.readFile("s.tif",t2,{})).images[0]; okR=r3.planes[0][47]===47; }catch(err){ okR=false; } ok(okR,"TIFF RowsPerStrip = 4294967295 (기본값) 읽기"); }
  { // ND2 v2 XML: 긴 태그(>4096 자)·속성 값 안의 '>'·숫자 문자 참조·작은따옴표
    const long="x".repeat(9000)+"a&gt;b&#x000d;&#x000a;c&#65;";
    const xml=`<?xml version="1.0"?><variant version="1.0"><no_name runtype="CLxListVariant"><TextInfoItem_5 runtype="CLxStringW" value="${long}"/><sOdd runtype='CLxStringW' value='p>q'/><uiWidth runtype="lx_uint32" value="12"/></no_name></variant>`;
    let v=null, err=""; try{ v=IO._xmlVariant(xml); }catch(e){ err=e.message; }
    ok(v&&v.TextInfoItem_5===("x".repeat(9000)+"a>b\r\ncA")&&v.sOdd==="p>q"&&v.uiWidth===12,"ND2 XML 메타데이터: 긴 태그·따옴표 안 '>'·&#..; 참조·작은따옴표",err||JSON.stringify(v&&{len:v.TextInfoItem_5&&v.TextInfoItem_5.length,odd:v.sOdd,w:v.uiWidth})); }
  { // TIFF z step 왕복: 보정 단위가 nm 여도 zStep(µm) 보존
    const src={w:4,h:3,bit:16,nc:1,nz:3,nt:1,planes:[0,1,2].map(k=>new Float32Array(12).fill(k)),cal:{pw:50,ph:50,unit:"nm"},zStep:0.05};
    const b=(await IO.readFile("z.tif",IO.writeTIFF(src),{})).images[0]; const b2=(await IO.readFile("z2.tif",IO.writeTIFF({...src,cal:b.cal,zStep:b.zStep}),{})).images[0];
    ok(Math.abs(b.zStep-0.05)<1e-12&&Math.abs(b2.zStep-0.05)<1e-12,"TIFF spacing ↔ z step(µm) 왕복 (nm 보정, 두 번 저장)",`${b.zStep}, ${b2.zStep}`); }
  { // 전수 검수 회귀 (ref/gen_audit.py)
    const DA=path.join(D,"audit"); const TA=JSON.parse(fs.readFileSync(path.join(DA,"truth.json"),"utf8")); const C2={};
    for(const key of Object.keys(TA)){ const t=TA[key]; const [file,idx]=key.split("#"); const u8=new Uint8Array(fs.readFileSync(path.join(DA,file)));
      if(t.expect_error){ let err=null; try{ await IO.readFile(file,u8,{}); }catch(e){ err=e; } ok(err&&err.message.includes(t.expect_error),`${key}: 잘린 파일 → 분명한 오류`,err?err.message:"오류 없이 열림"); continue; }
      try{ if(!C2[file]) C2[file]=await IO.readFile(file,u8,{}); const r=C2[file]; const im=r.images[idx?+idx:0]; const tr=rdF32(path.join(DA,key+".truth"));
        const flat=new Float32Array(im.planes.length*im.planes[0].length); im.planes.forEach((q,k)=>flat.set(q,k*q.length)); let bad=flat.length!==tr.length?1e9:0;
        if(!bad) for(let i=0;i<tr.length;i++){ const a=flat[i], b=tr[i]; if(!(a===b||(a!==a&&b!==b))) bad++; }
        let ex=""; if(t.nc&&im.nc!==t.nc) ex+=` nc ${im.nc}≠${t.nc}`; if(t.nz&&im.nz!==t.nz) ex+=` nz ${im.nz}≠${t.nz}`; if(t.rgb&&!im.rgb) ex+=" RGB 아님"; if(!t.rgb&&im.rgb) ex+=" RGB 로 읽힘";
        if(t.pw){ const k=im.cal?IO.unitToM(im.cal.unit)/1e-6:0; if(!im.cal||Math.abs(im.cal.pw*k-t.pw)>1e-6*t.pw) ex+=` 보정 ${im.cal&&im.cal.pw} ${im.cal&&im.cal.unit}`; }
        if(t.zstep&&!(Math.abs((im.zStep||0)-t.zstep)<1e-9)) ex+=` z step ${im.zStep}`; if(t.display&&!(im.display&&im.display[0].min===t.display[0]&&im.display[0].max===t.display[1])) ex+=` 표시 범위 ${JSON.stringify(im.display)}`;
        if(idx==null&&r.images.length!==Object.keys(TA).filter(k=>k.split("#")[0]===file).length) ex+=` 영상 수 ${r.images.length}`;
        ok(!bad&&!ex,`검수 회귀 ${key}`,`다른 값 ${bad}${ex}`); }catch(e){ ok(false,`검수 회귀 ${key}`,"예외: "+e.message); } }
    // 쓰기: 해상도 유리수 정밀도, 8/16-bit 반올림·자르기
    const src={w:3,h:2,bit:16,nc:1,nz:1,nt:1,planes:[Float32Array.from([70000,-5,100.7,100.4,0,65535])],cal:{pw:2512.3,ph:2512.3,unit:"nm"}};
    const b=(await IO.readFile("w.tif",IO.writeTIFF(src),{})).images[0]; const pw=b.cal.pw*IO.unitToM(b.cal.unit)/1e-9;
    ok(Math.abs(pw-2512.3)<1e-6*2512.3&&JSON.stringify(Array.from(b.planes[0]))==="[65535,0,101,100,0,65535]","TIFF 쓰기: 큰 화소 크기 정밀도·16-bit 반올림·자르기",`${pw} nm, ${Array.from(b.planes[0])}`); }
  { const s=IC.sizeSummary([-2,-1,0,1,2,3]); ok(s.n===6&&Math.abs(s.mean-0.5)<1e-12&&s.nPos===3,"크기 요약: 선형 통계는 모든 값, 로그정규는 양수만",JSON.stringify({n:s.n,mean:s.mean,nPos:s.nPos})); }
}

console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.slice(0,80).join("\n")); process.exit(1); }
})();
