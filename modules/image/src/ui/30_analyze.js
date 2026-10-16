/* ============================================================
   Analyze · Measure 패널, 측정 연산, 아래쪽 표(Results·Particles·Summary·Distribution·Plots), 그래프
   ============================================================ */
const ALLCOLS=[["Area","Area"],["Mean","Mean"],["StdDev","StdDev"],["Mode","Modal"],["Min","Min"],["Max","Max"],["X","Centroid X"],["Y","Centroid Y"],["XM","Center of mass X"],["YM","Center of mass Y"],["Perim.","Perimeter"],
  ["BX","Bounding BX"],["BY","BY"],["Width","Width"],["Height","Height"],["Major","Ellipse Major"],["Minor","Minor"],["Angle","Angle"],["Circ.","Circularity"],["AR","Aspect ratio"],["Round","Roundness"],["Solidity","Solidity"],
  ["Feret","Feret"],["MinFeret","MinFeret"],["FeretAngle","FeretAngle"],["FeretX","FeretX"],["FeretY","FeretY"],["IntDen","IntDen"],["RawIntDen","RawIntDen"],["Median","Median"],["Skew","Skewness"],["Kurt","Kurtosis"],["ECD","ECD (등가 원 지름)"]];
const LENCOL=new Set(["X","Y","XM","YM","Perim.","BX","BY","Width","Height","Major","Minor","Feret","MinFeret","ECD"]);
/* 열 단위: 측정할 때의 보정 단위(행·분석 결과에 저장한 값)로 — 지금 고른 영상의 보정이 아니라 */
const colUnit=(k,u)=>{ if(u==null) return ""; if(k==="Area") return u==="px"?"px²":u+"²"; if(LENCOL.has(k)||k==="Length") return u; if(k==="Angle"||k==="FeretAngle") return "°"; return ""; };

/* ---------- Measure ---------- */
function measureRoiRow(im,r,label){
  const p=curPlaneOf(im); const c=calOf(im); const base={Label:label,Image:im.name,Slice:isStack(im)?`${im.c+1}:${im.z+1}:${im.t+1}`:"",Unit:c.unit};
  if(!r||AREA_T.has(r.type)){ const g=r?roiGeom(im,r):{mask:null,rect:{x:0,y:0,w:im.w,h:im.h},poly:null,kind:"rect"}; const m=IC.measureRoi(p,g,c); delete m._fe; delete m._el; delete m._n; return Object.assign(base,m,{Type:r?r.type:"image"}); }
  if(r.type==="line"){ const pr=profileOf(im,r); const st=IC.stats({w:pr.y.length,h:1,bit:32,data:Float32Array.from(pr.y)}); return Object.assign(base,{Type:"line",Length:lineLen(im,r),Angle:roiAngle(r),Mean:st.mean,StdDev:st.stdDev,Min:st.min,Max:st.max,X:(r.pts[0][0]+r.pts[1][0])/2*c.pw,Y:(r.pts[0][1]+r.pts[1][1])/2*c.ph}); }
  if(r.type==="angle") return Object.assign(base,{Type:"angle",Angle:roiAngle(r)});
  if(r.type==="point") return r.pts.map((q,k)=>{ const x=Math.floor(q[0]), y=Math.floor(q[1]); return Object.assign({},base,{Label:`${label}:${k+1}`,Type:"point",X:q[0]*c.pw,Y:q[1]*c.ph,Mean:p.data[y*im.w+x],Class:r.cls?r.cls[k]+1:1}); });
  return base;
}
function addResults(rows){ rows=[].concat(rows); rows.forEach(r=>S.results.push(r)); S.ui.bottom="results"; S.ui.bottomOpen=true; save(); drawBottom(); }
op("measure",{label:"Measure",group:"meas",hidden:true,kind:"measure",params:[R("all","ROI Manager 전체",false,{type:"bool"})],
  run:(im,P)=>{ const rec=P.roi&&REPLAY?Object.assign({id:"rec"},P.roi):null; const L=P.all?roiList(im):[rec||activeRoi(im)].filter(Boolean); const rows=[]; // 재생: 기록한 ROI 로 잼
    if(!L.length) rows.push(measureRoiRow(im,null,im.name)); else L.forEach((r,i)=>{ const x=measureRoiRow(im,r,`${im.name}:${r.name||(P.all?i+1:r.type)}`); [].concat(x).forEach(q=>rows.push(q)); });
    const a=activeRoi(im); addResults(rows); record({op:"measure",params:Object.assign(clean(P),!P.all&&a?{roi:{type:a.type,pts:a.pts}}:{}),label:P.all?"Measure (all ROIs)":`Measure${a?` (${a.type} ROI)`:""}`},im); },
  ij:P=>P.all?`roiManager("Measure");`:`${roiIJM(P.roi)}run("Measure");`});
function roiIJM(r){ if(!r) return "run(\"Select None\");\n"; const f=v=>+(+v).toFixed(4); const xs=r.pts.map(p=>f(p[0])), ys=r.pts.map(p=>f(p[1]));
  if(r.type==="rect"){ const x=Math.min(...xs), y=Math.min(...ys); return `makeRectangle(${x}, ${y}, ${Math.max(...xs)-x}, ${Math.max(...ys)-y});\n`; }
  if(r.type==="oval"){ const x=Math.min(...xs), y=Math.min(...ys); return `makeOval(${x}, ${y}, ${Math.max(...xs)-x}, ${Math.max(...ys)-y});\n`; }
  if(r.type==="line") return `makeLine(${xs[0]}, ${ys[0]}, ${xs[1]}, ${ys[1]});\n`;
  if(r.type==="point") return `makeSelection("point", newArray(${xs.join(",")}), newArray(${ys.join(",")}));\n`;
  const kind=r.type==="free"?"freehand":r.type==="traced"?"traced":r.type==="angle"?"angle":"polygon"; return `makeSelection("${kind}", newArray(${xs.join(",")}), newArray(${ys.join(",")}));\n`; }

/* ---------- Analyze Particles ---------- */
async function binaryFor(im,P){ // 이진이면 255, 아니면 Threshold 설정(방법 또는 수동 범위)
  const p=curPlaneOf(im); const isB=IC.isBinary(p); let lo=255, hi=255, how="binary";
  if(!isB&&REPLAY&&"thr" in P&&!P.thr) throw new Error("기록할 때는 이진(mask) 영상이었는데 이 영상은 이진이 아닙니다 — 매크로에 Convert to Mask 단계를 넣으세요");
  if(!isB){ const T=P.thr||S.ui.thr; if(T.manual){ lo=+T.lower; hi=+T.upper; how=`threshold ${g4(lo)}–${g4(hi)}`; } else { const t=await call("autoThreshold",[p,T.method,!!T.dark]); lo=t.lower; hi=t.upper; how=`${T.method}${T.dark?" dark":""} [${g4(lo)}, ${g4(hi)}]`; } }
  const bin=new Uint8Array(p.data.length); for(let i=0;i<bin.length;i++){ const v=p.data[i]; bin[i]=(v>=lo&&v<=hi)?1:0; }
  return {bin,how,thr:isB?null:{lo,hi}};
}
op("analyze",{label:"Analyze Particles",group:"meas",hidden:true,kind:"measure",params:[],
  run:async(im,P)=>{
    if(im.rgb) throw new Error("RGB 영상은 먼저 Type ▸ RGB → 8-bit 또는 Threshold ▸ Convert to Mask");
    const {bin,how,thr}=await binaryFor(im,P); let red=im;
    const rootOf=x=>{ let r=x; while(r.parent&&byId(r.parent)) r=byId(r.parent); return r; };
    if(P.redirect==="__parent") red=byId(im.parent)||im; else if(P.redirect==="__root") red=rootOf(im); else if(P.redirect) red=byId(P.redirect)||S.images.find(x=>x.name===P.redirectName)||im;
    if(red.w!==im.w||red.h!==im.h) throw new Error("측정용(redirect) 영상의 크기가 다릅니다");
    const inten=red===im&&IC.isBinary(curPlaneOf(im))?curPlaneOf(im):(red.rgb?curPlaneOf(red):planeAt(red,red.c,im.z<red.nz?im.z:0,im.t<red.nt?im.t:0));
    const c=calOf(im), a1=c.pw*c.ph; const toPx=v=>v===""||v==null||!isFinite(+v)?null:(P.unitPx?+v:+v/a1);
    const o={minSize:toPx(P.minSize)||0,maxSize:toPx(P.maxSize)==null?Infinity:toPx(P.maxSize),minCirc:+P.minCirc||0,maxCirc:P.maxCirc===""||P.maxCirc==null?1:+P.maxCirc,excludeEdges:!!P.exclude,includeHoles:!!P.holes};
    setProg("Analyze Particles…",0.1);
    const r=await call("analyzeParticles",[bin,im.w,im.h,inten,{pw:c.pw,ph:c.ph},o],f=>setProg("Analyze Particles…",f)); hideProg();
    const rows=r.particles.map(q=>{ const m=q.meas; delete m._fe; delete m._el; delete m._n; return Object.assign({"#":q.n},m); });
    const polys=r.particles.map(q=>[q.poly.xs,q.poly.ys]); const cent=r.particles.map(q=>[q.meas.X/c.pw,q.meas.Y/c.ph]);
    const totPx=bin.reduce((a,v)=>a+v,0)||1; const areas=rows.map(x=>x.Area); const tot=areas.reduce((a,b)=>a+b,0);
    const fieldA=im.w*im.h*a1;
    const sum={Image:im.name,Slice:isStack(im)?(im.nt>1?`z${im.z+1} t${im.t+1}`:`${im.z+1}`):"",Count:rows.length,"Total Area":tot,"Average Size":rows.length?tot/rows.length:NaN,"%Area":100*tot/fieldA,
      "ECD mean":mean(rows.map(x=>x.ECD)),"ECD SD":sd(rows.map(x=>x.ECD)),"Feret mean":mean(rows.map(x=>x.Feret)),"Circ. mean":mean(rows.map(x=>x["Circ."])),"Density (/area)":rows.length/fieldA,Unit:c.unit,Threshold:how};
    S.particles[im.id]={params:clean(P),rows,polys,cent,meanArea:rows.length?mean(rows.map(x=>x.Area))/a1:25,summary:sum,sel:[],thr,redirect:red.name,unit:c.unit};
    S.summary.push(sum);
    if(P.addRoi){ const L=roiList(im); r.particles.slice(0,5000).forEach(q=>L.push({id:uid(),type:"traced",pts:q.poly.xs.map((x,k)=>[x,q.poly.ys[k]]),name:String(q.n)})); }
    // 기록: 세기 영상이 이 영상의 원본 쪽이면 상대 지정(__parent/__root)으로 — batch 에서 파일마다 자기 원본을 쓴다
    let rdRec=P.redirect||""; if(red!==im){ if(byId(im.parent)===red) rdRec="__parent"; else if(rootOf(im)===red) rdRec="__root"; }
    record({op:"analyze",params:Object.assign(clean(P),{redirect:rdRec,thr:P.thr||(IC.isBinary(curPlaneOf(im))?null:Object.assign({},S.ui.thr)),redirectName:red!==im&&!rdRec.startsWith("__")?red.name:""}),label:`Analyze Particles (${rows.length})`},im);
    S.ui.bottom="particles"; S.ui.bottomOpen=true; save();
    note(`입자 ${rows.length}개 — ${how}${red!==im?`, 세기는 '${red.name}' 에서`:""}${totPx&&rows.length===0?" (조건에 맞는 입자 없음)":""}`,"ok"); },
  ij:P=>{ const s=`size=${P.minSize||0}-${P.maxSize===""||P.maxSize==null?"Infinity":P.maxSize}${P.unitPx?" pixel":""} circularity=${(+P.minCirc||0).toFixed(2)}-${(P.maxCirc===""||P.maxCirc==null?1:+P.maxCirc).toFixed(2)} show=${P.show==="masks"?"Masks":P.show==="none"?"Nothing":"Overlay"} display${P.exclude?" exclude":""}${P.holes?" include":""}${P.addRoi?" add":""} summarize`;
    const T=P.thr; const thr=!T?"":T.manual?`setThreshold(${T.lower}, ${T.upper});\n`:`setAutoThreshold("${T.method}${T.dark?" dark":""}");\n`;
    const rel=P.redirect==="__parent"||P.redirect==="__root"; const rd=rel?"None":(P.redirectName?`[${ijq(P.redirectName)}]`:"None");
    return `${rel?"// redirect: 원본(세기) 영상 이름을 아래 redirect=[...] 에 넣으세요 (앱에서는 이 영상의 원본)\n":""}${thr}run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis redirect=${rd} decimal=4");\nrun("Analyze Particles...", "${s}");`; }});
const ijq=s=>String(s).replace(/\\/g,"\\\\").replace(/"/g,'\\"'); // ImageJ 매크로 문자열 안전하게

const mean=a=>{ const v=a.filter(x=>x!=null&&isFinite(x)); return v.length?v.reduce((s,x)=>s+x,0)/v.length:NaN; };
const sd=a=>{ const v=a.filter(x=>x!=null&&isFinite(x)); if(v.length<2) return NaN; const m=mean(v); return Math.sqrt(v.reduce((s,x)=>s+(x-m)*(x-m),0)/(v.length-1)); };

/* ---------- Find Maxima ---------- */
op("findMaxima",{label:"Find Maxima",group:"meas",hidden:true,kind:"measure",params:[R("prominence","Prominence >",10,{min:0,step:"any"}),R("light","Light background (어두운 점)",false,{type:"bool"}),R("exclude","Exclude edge maxima",false,{type:"bool"}),R("strict","Strict",false,{type:"bool"}),
  R("output","Output","points",{type:"sel",opts:[["points","Point Selection + Count"],["single","Single Points (영상)"],["within","Maxima Within Tolerance (영상)"],["segmented","Segmented Particles (영상)"]]})],
  run:async(im,P)=>{ const p=curPlaneOf(im); const o={output:P.output,light:!!P.light,excludeEdges:!!P.exclude,strict:!!P.strict};
    if(P.output==="points"){ const r=await call("findMaxima",[p,+P.prominence,o]); const pts=r.points.map(([x,y])=>[x+0.5,y+0.5]); setActive(im,{id:uid(),type:"point",pts,cls:pts.map(()=>0),name:"maxima"});
      S.summary.push({Image:im.name,Slice:isStack(im)?(im.nt>1?`z${im.z+1} t${im.t+1}`:`${im.z+1}`):"",Count:pts.length,"Total Area":"","Average Size":"","%Area":"",Threshold:`Find Maxima prominence>${P.prominence}${P.light?" light":""}`,Unit:calOf(im).unit,"Density (/area)":pts.length/(im.w*im.h*calOf(im).pw*calOf(im).ph)});
      record({op:"findMaxima",params:clean(P),label:`Find Maxima (${pts.length})`},im); note(`극대 ${pts.length}개 (prominence > ${P.prominence})`,"ok"); S.ui.bottom="summary"; S.ui.bottomOpen=true; return; }
    const r=await call("findMaximaImage",[p,+P.prominence,o]); const nm={single:"Maxima",within:"Maxima",segmented:"Segmented"}[P.output];
    const n=addImage({name:`${im.name} ${nm}`,w:im.w,h:im.h,bit:8,planes:[r.data],cal:im.cal},{parent:im.id,src:"derived",kind:"mask"}); record({op:"findMaxima",params:clean(P),label:`Find Maxima → ${nm}`},im); return n; },
  ij:P=>`run("Find Maxima...", "prominence=${P.prominence}${P.strict?" strict":""}${P.exclude?" exclude":""}${P.light?" light":""} output=${{points:"[Point Selection]",single:"[Single Points]",within:"[Maxima Within Tolerance]",segmented:"[Segmented Particles]"}[P.output]}");`,
  help:"MaximumFinder 이식: 더 높은 극대로 이어지는 능선보다 prominence 이상 솟은 점만. 원자 열·나노입자 세기·점 계수에."});
/* ---------- Radial profile · Orientation · Z profile ---------- */
op("radial",{label:"Radial Profile",group:"meas",hidden:true,kind:"measure",params:[R("rmax","최대 반지름 (0 = 자동)",0,{min:0,unit:"px"}),R("step","칸 폭",1,{min:0.25,step:0.25,unit:"px"})],
  run:async(im,P)=>{ const a=activeRoi(im); let cx=im.w/2, cy=im.h/2; let how="영상 중심";
    if(im.kind==="fft"){ cx=Math.floor(im.w/2)+0.5; cy=Math.floor(im.h/2)+0.5; how="FFT 원점(DC 화소 중심)"; }
    if(a&&a.type==="point"&&a.pts.length){ [cx,cy]=a.pts[0]; how="점 ROI"; } else if(a&&(a.type==="oval"||a.type==="rect")){ const b=roiBounds(a); cx=b.x+b.w/2; cy=b.y+b.h/2; how=a.type+" ROI 중심"; }
    const rmax=+P.rmax>0?+P.rmax:Math.min(cx,cy,im.w-cx,im.h-cy); const r=await call("radialProfile",[curPlaneOf(im),cx,cy,rmax,+P.step||1]); const c=calOf(im);
    const recip=/^1\//.test(c.unit); S.plots.push({id:uid(),kind:"xy",title:`Radial profile — ${im.name} (${how} ${cx.toFixed(1)}, ${cy.toFixed(1)})`,xl:`r (${c.unit})`,yl:"mean",series:[{name:"radial",x:r.r.map(v=>v*c.pw),y:r.mean}],recip});
    record({op:"radial",params:clean(P),label:"Radial Profile"},im); S.ui.bottom="plots"; S.ui.bottomOpen=true; save(); note(`Radial profile (${how}) — Plots 탭${recip?". x 축이 역공간이라 각 피크의 d = 1/r 를 함께 보입니다":""}`,"ok"); },
  ij:()=>`// Radial Profile (plugin): Analyze ▸ Radial Profile`,help:"방위 평균 세기. 회절(SAED/FFT) 고리의 반지름, 둥근 입자 단면 등. 중심 = 점 ROI · 타원/사각형 ROI 중심 · 없으면 영상 중심."});
op("orientation",{label:"Orientation (구조 텐서)",group:"meas",hidden:true,kind:"measure",params:[R("sigma","Local window σ",2,{min:0.5,step:0.5,unit:"px"}),R("minC","최소 coherency",0.1,{min:0,max:1,step:0.05}),R("map","방향 색 지도 영상 만들기",true,{type:"bool"})],
  run:async(im,P)=>{ const r=await call("orientation",[curPlaneOf(im),{sigma:+P.sigma,minCoherency:+P.minC}]);
    S.plots.push({id:uid(),kind:"xy",title:`Orientation — ${im.name}: 주 방향 ${r.dominant===r.dominant?r.dominant.toFixed(1)+"°":"없음"}, order ${r.orderParam.toFixed(3)}, coherency ${r.meanCoherency.toFixed(3)}`,xl:"angle (°, 가로축 기준 반시계 +)",yl:"weighted count",series:[{name:"histogram",x:[...Array(180).keys()].map(k=>k-89.5),y:[...r.hist]}]});
    if(P.map){ const N=im.w*im.h, R1=new Float32Array(N),G1=new Float32Array(N),B1=new Float32Array(N); let emax=0; for(let i=0;i<N;i++) if(r.energy[i]>emax) emax=r.energy[i];
      for(let i=0;i<N;i++){ const h=((r.theta[i]+90)/180)%1; const [rr,gg,bb]=hsb(h); const s=r.coherency[i], b=Math.min(1,Math.sqrt(r.energy[i]/(emax||1))*1.5); R1[i]=Math.round(255*(1-s+s*rr/255)*b); G1[i]=Math.round(255*(1-s+s*gg/255)*b); B1[i]=Math.round(255*(1-s+s*bb/255)*b); }
      var map=addImage({name:"Orientation of "+im.name,w:im.w,h:im.h,bit:8,rgb:true,planes:[R1,G1,B1],cal:im.cal},{parent:im.id,src:"derived"}); }
    record({op:"orientation",params:clean(P),label:"Orientation"},im); S.ui.bottom="plots"; S.ui.bottomOpen=true; save(); note(r.dominant===r.dominant?`주 방향 ${r.dominant.toFixed(1)}° · order parameter ${r.orderParam.toFixed(3)} (0 무질서 – 1 완전 정렬)`:"방향이 없습니다 (일정한 영상)","ok"); return P.map?map:undefined; }, // 색 지도가 다음 단계의 대상 (재생도 같게)
  ij:()=>`// OrientationJ (plugin) Analysis/Distribution 과 같은 정의`,help:"OrientationJ 와 같은 구조 텐서 정의. 색 지도: 색상 = 방향, 채도 = coherency, 밝기 = energy. 나노선·섬유·결 방향 분석."});
op("zprofile",{label:"Plot Z-axis Profile",group:"meas",hidden:true,kind:"measure",params:[],
  run:(im)=>{ if(!isStack(im)) throw new Error("스택이 아닙니다"); const a=activeRoi(im); const g=a&&AREA_T.has(a.type)?roiGeom(im,a):null; const ys=[];
    for(let t=0;t<im.nt;t++) for(let z=0;z<im.nz;z++){ const p=planeAt(im,im.c,z,t); const full=g?(()=>{ const m=new Uint8Array(im.w*im.h); const Rr=g.rect; for(let y=0;y<Rr.h;y++) for(let x=0;x<Rr.w;x++) if(!g.mask||g.mask[y*Rr.w+x]) m[(y+Rr.y)*im.w+x+Rr.x]=1; return m; })():null; ys.push(IC.stats(p,full).mean); }
    S.plots.push({id:uid(),kind:"xy",title:`Z-axis profile — ${im.name}${a?" ("+a.type+" ROI)":""}`,xl:"slice",yl:"mean",series:[{name:"mean",x:ys.map((_,i)=>i+1),y:ys}]}); S.ui.bottom="plots"; S.ui.bottomOpen=true; save(); },
  ij:()=>`run("Plot Z-axis Profile");`});

/* ===== Analyze 패널 ===== */
/* ---------- Analyze Skeleton (2D/3D) — Fiji AnalyzeSkeleton 과 같은 결과 (core/85_skeleton.js: 동작 명세·블랙박스 대조로 새로 쓴 독립 구현) ---------- */
const SK_PRUNE=[["none","none"],["shortest branch","shortest branch"],["lowest intensity voxel","lowest intensity voxel"],["lowest intensity branch","lowest intensity branch"]];
op("skeleton",{label:"Analyze Skeleton (2D/3D)",group:"analyze",kind:"analyze",hidden:true,params:[R("skel","먼저 Skeletonize (2D, ImageJ)",false,{type:"bool"}),
  R("prune","Prune cycle method","none",{type:"sel",opts:SK_PRUNE}),R("orig","Original image (세기)","",{type:"img",showIf:P=>/intensity/.test(P.prune||"")}),R("pruneEnds","Prune ends",false,{type:"bool"}),R("lsp","Longest shortest path",true,{type:"bool"}),R("tagged","Tagged skeleton 영상",true,{type:"bool"}),R("lspImg","Longest shortest path 영상",false,{type:"bool"})],
  run:async(im,P)=>{ if(im.rgb) throw new Error("이진(0/255) 영상에서 — RGB 는 Color Threshold 로 mask 를 만든 뒤");
    const d=im.nz, w=im.w, h=im.h, P0=planesOf(im); const data=new Uint8Array(w*h*d); let nonBin=0;
    for(let z=0;z<d;z++){ let pl=P0[pIndex(im,im.c,z,im.t)];
      if(P.skel){ if(d>1) throw new Error("3D 스택의 Skeletonize 는 없습니다 (Fiji Skeletonize3D) — 이미 골격인 스택에서 쓰거나, 장마다 2D 로는 Binary ▸ Skeletonize");
        if(im.bit!==8||!IC.isBinary({w,h,bit:im.bit,data:pl})) throw new Error("Skeletonize 는 이진(0/255) 8-bit 영상에서 — Threshold ▸ Convert to Mask 먼저"); pl=(await call("skeletonize",[{w,h,bit:8,data:pl}])).data; }
      for(let i=0;i<w*h;i++){ const v=pl[i]; if(v!==0&&v!==255) nonBin++; data[z*w*h+i]=v!==0&&v===v?255:0; } }
    if(nonBin) throw new Error(`이진(0/255) 영상이 아닙니다 (0·255 아닌 화소 ${nonBin}) — Threshold ▸ Convert to Mask, Binary ▸ Skeletonize 뒤에`);
    const c=calOf(im); let pd=1; const later=[];
    if(d>1){ const k=im.cal?IO.unitToM(c.unit):0; if(im.zStep>0&&k) pd=im.zStep*1e-6/k;
      else later.push(`z 간격을 모릅니다 (${!(im.zStep>0)?"z step 없음":`단위 '${c.unit}' 를 길이로 못 바꿈`}) — z 간격 = 1 화소로 셈. Image ▸ Properties 에서 지정하면 3D 길이가 맞습니다`); }
    if(im.nc>1||im.nt>1) later.push(`hyperstack: 현재 ${im.nc>1?`채널 ${im.c+1}`:""}${im.nc>1&&im.nt>1?" · ":""}${im.nt>1?`시점 ${im.t+1}`:""} 의 z 스택만 분석했습니다 (Fiji 는 hyperstack 의 모든 장을 한 3D 로 봄 — 같은 결과는 Duplicate 로 떼어 낸 뒤)`);
    const wantLsp=!!(P.lsp||P.lspImg); const prune=P.prune||"none"; let orig=null, origName="";
    if(/intensity/.test(prune)){ const o=byId(P.orig)||S.images.find(x=>x.name===P.origName); if(!o) throw new Error("세기 기준 가지치기: 원본(회색) 영상을 고르세요");
      if(o.w!==w||o.h!==h||o.nz!==d) throw new Error(`원본 영상의 크기가 골격과 같아야 합니다 (${w}×${h}${d>1?"×"+d:""})`); if(o.rgb) throw new Error("원본 영상은 회색(8/16/32-bit)이어야 합니다");
      const OP=planesOf(o); orig=new Float32Array(w*h*d); for(let z=0;z<d;z++) orig.set(OP[pIndex(o,o.c,z,o.t)],z*w*h); origName=o.name; }
    setProg("Analyze Skeleton …",0.3); let r; try{ r=await call("analyzeSkeleton",[{w,h,d,data},{pw:c.pw,ph:c.ph,pd},{lsp:wantLsp,prune,pruneEnds:!!P.pruneEnds,orig:orig?{data:orig}:null}],f=>setProg("Analyze Skeleton …",0.3+0.6*f)); } finally { hideProg(); }
    const slice=[im.nc>1?`c${im.c+1}`:"",im.nt>1?`t${im.t+1}`:""].filter(Boolean).join(" "); const u=c.unit;
    const CAP_T=20000, CAP_B=100000; const room=n=>Math.max(0,n);
    const nT=Math.min(r.trees.length,room(CAP_T-S.skelTrees.length)); if(nT<r.trees.length) later.push(`Skeleton 표가 ${CAP_T} 줄을 넘어 골격 ${r.trees.length-nT} 개를 표에 넣지 않았습니다 — Clear 뒤 다시, 또는 작은 영역으로`);
    const lspTxt=v=>v===v?v:"건너뜀";
    r.trees.slice(0,nT).forEach((t,k)=>S.skelTrees.push({Image:im.name,Slice:slice,Skeleton:k+1,"# Branches":t.branches,"# Junctions":t.junctions,"# End-point voxels":t.endPoints,"# Junction voxels":t.junctionVoxels,"# Slab voxels":t.slabs,
      "Average Branch Length":t.avgLength,"# Triple points":t.triples,"# Quadruple points":t.quadruples,"Maximum Branch Length":t.maxLength,...(wantLsp?{"Longest Shortest Path":lspTxt(t.lsp),spx:lspTxt(t.spx),spy:lspTxt(t.spy),spz:lspTxt(t.spz)}:{}),Unit:u}));
    // Branch information: 나무마다 길이 내림차순 (Fiji 'Show detailed info' 와 같은 순서 — 같은 길이는 원래 순서)
    const byT=r.trees.map(()=>[]); r.edges.forEach((e,i)=>byT[e.tree].push([e,i])); let bLeft=room(CAP_B-S.skelBranches.length), bDrop=0;
    for(let k=0;k<nT;k++){ const E=byT[k].sort((a,b)=>b[0].length-a[0].length||a[1]-b[1]).map(x=>x[0]);
      for(const e of E){ if(bLeft<=0){ bDrop++; continue; } bLeft--; S.skelBranches.push({Image:im.name,Slice:slice,"Skeleton ID":k+1,"Branch length":e.length,"V1 x":e.v1[0]*c.pw,"V1 y":e.v1[1]*c.ph,"V1 z":e.v1[2]*pd,"V2 x":e.v2[0]*c.pw,"V2 y":e.v2[1]*c.ph,"V2 z":e.v2[2]*pd,"Euclidean distance":e.euclid,"running average length":e.lengthRa,Unit:u}); } }
    for(let k=nT;k<r.trees.length;k++) bDrop+=byT[k].length;
    if(bDrop) later.push(`Branch information 이 ${CAP_B} 줄을 넘어 가지 ${bDrop} 개를 표에 넣지 않았습니다 (Skeleton 표의 수는 전체)`);
    const tb=r.trees.reduce((a,t)=>a+t.branches,0), tj=r.trees.reduce((a,t)=>a+t.junctions,0);
    S.summary.push({Image:im.name,Slice:slice,Count:r.trees.length,Threshold:`Analyze Skeleton: 골격 ${r.trees.length} · 가지 ${tb} · 접점 ${tj}`,Unit:u});
    let out;
    if(P.tagged){ const pl=[]; for(let z=0;z<d;z++) pl.push(Float32Array.from(r.tag.subarray(z*w*h,(z+1)*w*h)));
      out=addImage({name:"Tagged skeleton of "+im.name,w,h,bit:8,nc:1,nz:d,nt:1,planes:pl,cal:im.cal,zStep:im.zStep,info:[["Tag","끝점 30 · 줄기 127 · 접점 70 (Fiji 와 같음)"]]},{parent:im.id,src:"derived"}); out.disp[0].lut="Fire"; out.disp[0].min=0; out.disp[0].max=255; }
    if(P.lspImg&&r.lspImage){ const pl=[]; for(let z=0;z<d;z++) pl.push(Float32Array.from(r.lspImage.subarray(z*w*h,(z+1)*w*h)));
      const o2=addImage({name:"Longest shortest paths of "+im.name,w,h,bit:8,nc:1,nz:d,nt:1,planes:pl,cal:im.cal,zStep:im.zStep,info:[["Value","골격 255 · 가장 긴 최단 경로의 줄기 96"]]},{parent:im.id,src:"derived"}); o2.disp[0].lut="Fire"; o2.disp[0].min=0; o2.disp[0].max=255; if(!out) out=o2; }
    record({op:"skeleton",params:{...clean(P),origName},label:`Analyze Skeleton (${r.trees.length} 골격, 가지 ${tb}${prune!=="none"||P.pruneEnds?" · pruned":""})`},im);
    S.ui.bottom="skeleton"; S.ui.bottomOpen=true; save(); note(`골격 ${r.trees.length} · 가지 ${tb} · 접점 ${tj}${d>1?` (3D, z 간격 ${g4(pd)} ${u})`:""}`,"ok");
    (r.warn||[]).concat(later).forEach(x=>note(x,"warn")); /* 경고는 ok 뒤에 — 묻히지 않게 */
    return undefined; },
  ij:P=>{ const pr=P.prune||"none", dup=P.skel||pr!=="none"||P.pruneEnds; // Skeletonize·가지치기는 Fiji 에서 입력을 바꿈 → 사본에서 (앱은 원본을 바꾸지 않음)
    return (dup?`orig=getImageID(); run("Duplicate...", "title=[Skeleton of "+getTitle()+"]");${P.skel?' run("Skeletonize");':""} // 앱은 원본을 바꾸지 않음 — 사본에서\n`:"")
    +`run("Analyze Skeleton (2D/3D)", "prune=${pr==="none"?"none":"["+pr+"]"}${P.pruneEnds?" prune_0":""}${P.lsp||P.lspImg?" calculate":""} show${/intensity/.test(pr)?` original_image=[${ijq(P.origName||"original")}]`:""}"); // Fiji 는 Tagged skeleton${P.lsp||P.lspImg?"·Longest shortest paths":""} 영상을 늘 만듦`+(dup?"\nselectImage(orig);":""); },
  help:"Fiji Analyze Skeleton 과 같은 정의·결과: 끝점(<2 이웃)·접점(>2)·줄기(=2), 가지 길이 = 화소 중심 거리의 합. 가지치기: 고리(가장 짧은 가지·가장 어두운 화소·가장 어두운 가지 — 세기는 원본 영상)·끝 가지. Longest shortest path 는 크기 제한 없음."});
function panelAnalyze(el,im){
  const A=S.ui.pa; const bin=!im.rgb&&IC.isBinary(curPlaneOf(im)); const c=calOf(im); const us=A.unitPx?"px²":unitSq(im); const fm=formVals(OPS.findMaxima);
  const reds=S.images.filter(x=>x.id!==im.id&&x.w===im.w&&x.h===im.h);
  el.innerHTML=`<div class="im-h">Analyze Particles</div>
  <div class="xs muted">${bin?"이진 영상: 255 인 화소를 입자로 셉니다.":`문턱: <b>${S.ui.thr.manual?`수동 ${g4(S.ui.thr.lower)}–${g4(S.ui.thr.upper)}`:`${S.ui.thr.method}${S.ui.thr.dark?" dark":""}`}</b> (Threshold 탭 설정) 범위의 화소를 입자로 셉니다.`}</div>
  <div class="im-row"><label class="im-f"><span>Size</span><input type="number" id="paMin" value="${esc(A.minSize)}" step="any" min="0"><span>–</span><input type="number" id="paMax" value="${esc(A.maxSize)}" step="any" placeholder="∞"><span class="u">${esc(us)}</span></label><label class="im-chk"><input type="checkbox" id="paPx" ${A.unitPx?"checked":""}> pixel 단위</label></div>
  <div class="im-row"><label class="im-f"><span>Circularity</span><input type="number" id="paC0" value="${A.minCirc}" min="0" max="1" step="0.05"><span>–</span><input type="number" id="paC1" value="${A.maxCirc}" min="0" max="1" step="0.05"></label></div>
  <div class="im-row"><label class="im-chk"><input type="checkbox" id="paX" ${A.exclude?"checked":""}> Exclude on edges</label><label class="im-chk"><input type="checkbox" id="paH" ${A.holes?"checked":""}> Include holes</label><label class="im-chk"><input type="checkbox" id="paR" ${A.addRoi?"checked":""}> Add to ROI Manager</label></div>
  <div class="im-row"><label class="im-f"><span>세기 측정 영상 (Redirect)</span><select class="im-sel" id="paRed"><option value="">이 영상</option>${im.parent&&byId(im.parent)?`<option value="__parent" ${A.redirect==="__parent"?"selected":""}>원본: ${esc(byId(im.parent).name)}</option>`:""}${im.parent&&byId(im.parent)&&byId(im.parent).parent?`<option value="__root" ${A.redirect==="__root"?"selected":""}>맨 처음 원본</option>`:""}${reds.map(x=>`<option value="${x.id}" ${A.redirect===x.id?"selected":""}>${esc(x.name)}</option>`).join("")}</select></label>
    <label class="im-f"><span>Show</span><select class="im-sel" id="paShow">${[["outlines","Outlines + 번호"],["masks","Masks"],["none","없음"]].map(([a,b])=>`<option value="${a}" ${A.show===a?"selected":""}>${b}</option>`).join("")}</select></label></div>
  ${!im.cal?`<div class="note warn sm">보정이 없어 크기가 px 단위입니다 — Adjust ▸ Set Scale.</div>`:""}
  <div class="im-row"><span class="xs muted">Circularity = 4π·Area/Perimeter² · ECD = 2√(Area/π)</span><span class="sp"></span><button class="btn primary sm" id="paRun">Analyze</button></div>
  ${S.particles[im.id]?particleCard(im):""}
  <hr class="im-hr"><div class="im-h">Find Maxima <span class="muted xs">점 세기 · 원자 열</span></div>
  <div id="fmForm">${formHTML(OPS.findMaxima,im)}</div><div class="im-row"><button class="btn sm" id="fmPrev">Preview points</button><span class="sp"></span><button class="btn primary sm" id="fmRun">Find</button></div>
  <hr class="im-hr"><div class="im-h">Tools</div>
  <div class="im-tool"><b>Radial Profile</b><span class="xs muted">중심 = 점·타원 ROI</span><div id="rpForm">${formHTML(OPS.radial,im)}</div><button class="btn sm" id="rpRun">Plot</button></div>
  <div class="im-tool"><b>Orientation</b><span class="xs muted">구조 텐서 (OrientationJ 정의)</span><div id="orForm">${formHTML(OPS.orientation,im)}</div><button class="btn sm" id="orRun">Analyze</button></div>
  <div class="im-tool"><b>Analyze Skeleton (2D/3D)</b><span class="xs muted">가지·접점·길이 (Fiji 와 같은 결과)${im.nz>1?" · 스택 = 3D":""}</span><div id="skForm">${formHTML(OPS.skeleton,im)}</div><button class="btn sm" id="skRun">Analyze</button></div>
  ${isStack(im)?`<div class="im-tool"><b>Z-axis Profile</b><span class="xs muted">ROI 평균의 장별 변화</span><button class="btn sm" id="zpRun">Plot</button></div>`:""}`;
  const q=s=>$q(el,s); const up=()=>save();
  q("#paMin").onchange=e=>{ A.minSize=e.target.value===""?0:+e.target.value; up(); }; q("#paMax").onchange=e=>{ A.maxSize=e.target.value; up(); }; q("#paPx").onchange=e=>{ A.unitPx=e.target.checked; up(); renderPanel(); };
  q("#paC0").onchange=e=>{ A.minCirc=+e.target.value; up(); }; q("#paC1").onchange=e=>{ A.maxCirc=+e.target.value; up(); };
  q("#paX").onchange=e=>{ A.exclude=e.target.checked; up(); }; q("#paH").onchange=e=>{ A.holes=e.target.checked; up(); }; q("#paR").onchange=e=>{ A.addRoi=e.target.checked; up(); };
  q("#paRed").onchange=e=>{ A.redirect=e.target.value; up(); }; q("#paShow").onchange=e=>{ A.show=e.target.value; up(); draw(); };
  q("#paRun").onclick=()=>runOp("analyze",Object.assign({},A,{thr:bin?null:Object.assign({},S.ui.thr)}));
  bindForm(q("#fmForm"),OPS.findMaxima,re=>{ if(re) renderPanel(); }); bindForm(q("#rpForm"),OPS.radial,()=>{}); bindForm(q("#orForm"),OPS.orientation,()=>{}); bindForm(q("#skForm"),OPS.skeleton,re=>{ if(re) renderPanel(); }); q("#skRun").onclick=()=>runOp("skeleton");
  q("#fmRun").onclick=()=>runOp("findMaxima"); q("#fmPrev").onclick=()=>runOp("findMaxima",Object.assign({},fm,{output:"points"}));
  q("#rpRun").onclick=()=>runOp("radial"); q("#orRun").onclick=()=>runOp("orientation"); if(q("#zpRun")) q("#zpRun").onclick=()=>runOp("zprofile",{});
  $a(el,"[data-dist]").forEach(b=>b.onclick=()=>{ S.ui.bottom="dist"; S.ui.bottomOpen=true; drawBottom(); });
}
function particleCard(im){ const pa=S.particles[im.id]; const s=pa.summary; const ss=IC.sizeSummary(pa.rows.map(r=>r.ECD)); const u=pa.unit;
  return `<div class="im-card"><div class="im-row"><b>${s.Count} particles</b><span class="xs muted">${esc(s.Threshold)}</span><span class="sp"></span><button class="btn xs" data-dist="1">Size distribution</button></div>
  <div class="im-stats"><div><span>Area fraction</span><b>${fmt(s["%Area"],2)} %</b></div><div><span>Mean area</span><b>${g4(s["Average Size"])} ${esc(u==="px"?"px²":u+"²")}</b></div>
  <div><span>ECD mean ± SD</span><b>${g4(ss.mean)} ± ${g4(ss.sd)} ${esc(u)}</b></div><div><span>ECD median (D50)</span><b>${g4(ss.median)} ${esc(u)}</b></div>
  <div><span>Log-normal</span><b>x̄g ${g4(ss.gMean)}, σg ${g4(ss.gSD)}</b></div><div><span>Number density</span><b>${g4(s["Density (/area)"])} /${esc(u==="px"?"px²":u+"²")}</b></div></div>
  ${s.Count<100?`<div class="xs muted">입자 수가 ${s.Count}개라 평균의 표준오차가 ${g4(ss.se)} ${esc(u)} 입니다 — 분포를 말하려면 여러 영상을 모으세요 (Macro ▸ Batch).</div>`:""}</div>`; }

/* ===== Measure 패널 (ROI Manager · Set Measurements) ===== */
function panelMeasure(el,im){
  const L=roiList(im); const a=activeRoi(im);
  el.innerHTML=`<div class="im-h">ROI Manager <span class="muted xs">${L.length}</span><span class="sp"></span><label class="im-chk"><input type="checkbox" id="rmAll" ${S.ui.showAll!==false?"checked":""}> Show all</label></div>
  <div class="im-row"><button class="btn sm" id="rmAdd" ${a?"":"disabled"} title="T">Add [T]</button><button class="btn sm" id="rmUpd" ${a&&L.includes(a)?"":"disabled"}>Update</button><button class="btn sm" id="rmDel" ${a&&L.includes(a)?"":"disabled"}>Delete</button><button class="btn sm" id="rmClr" ${L.length?"":"disabled"}>Clear</button><span class="sp"></span><button class="btn sm" id="rmSave" ${L.length?"":"disabled"}>Save</button><label class="btn sm" style="cursor:pointer">Open<input type="file" id="rmOpen" accept=".json" hidden></label></div>
  <div class="im-roilist">${L.length?L.map((r,i)=>`<div class="it ${r===a?"sel":""}" data-i="${i}"><span class="ty">${r.type}</span><input class="nm" value="${esc(r.name||String(i+1))}" data-n="${i}"><span class="xs muted">${roiSummary(im,r)}</span></div>`).join(""):`<div class="xs muted">ROI 를 그리고 Add (T). 입자 분석에서 'Add to ROI Manager' 를 켜면 입자가 들어옵니다.</div>`}</div>
  <div class="im-row"><button class="btn primary sm" id="msOne" title="Ctrl+M">Measure [Ctrl+M]</button><button class="btn sm" id="msAll" ${L.length?"":"disabled"}>Measure all ROIs</button><span class="xs muted">${a?`현재 ROI: ${a.type}`:"ROI 가 없으면 영상 전체"}</span></div>
  ${a&&a.type==="point"?pointCounter(a):""}
  <div class="im-row"><label class="im-f"><span>Line width</span><input type="number" id="msLw" value="${S.ui.lw||1}" min="1" step="1"><span class="u">px</span></label><label class="im-f"><span>Wand tolerance</span><input type="number" id="msWt" value="${S.ui.wandTol||0}" min="0" step="any"></label></div>
  <hr class="im-hr"><div class="im-h">Set Measurements <span class="muted xs">표에 보일 열</span></div>
  <div class="im-cols">${ALLCOLS.map(([k,l])=>`<label class="im-chk"><input type="checkbox" data-col="${esc(k)}" ${S.measure.cols.includes(k)?"checked":""}> ${esc(l)}</label>`).join("")}</div>`;
  const q=s=>$q(el,s);
  q("#rmAll").onchange=e=>{ S.ui.showAll=e.target.checked; save(); draw(); };
  q("#rmAdd").onclick=()=>roiAdd(im); q("#rmUpd").onclick=()=>{ save(); renderPanel(); draw(); }; q("#rmDel").onclick=()=>{ const i=L.indexOf(a); if(i>=0) L.splice(i,1); setActive(im,null); renderPanel(); };
  q("#rmClr").onclick=()=>{ if(!confirm2("rmclr","ROI 를 모두 지울까요? 한 번 더 누르면 지웁니다")) return; S.rois[im.id]=[]; setActive(im,null); renderPanel(); };
  q("#rmSave").onclick=()=>downloadText(`${im.name}_RoiSet.json`,JSON.stringify({image:im.name,w:im.w,h:im.h,rois:L},null,1));
  q("#rmOpen").onchange=e=>{ const f=e.target.files[0]; if(!f) return; f.text().then(t=>{ try{ const o=JSON.parse(t); (o.rois||[]).forEach(r=>{ r.id=uid(); L.push(r); }); save(); renderPanel(); draw(); note(`ROI ${(o.rois||[]).length}개를 읽었습니다`,"ok"); }catch(err){ note("ROI 파일을 읽지 못했습니다","warn"); } }); e.target.value=""; };
  $a(el,".im-roilist .it").forEach(d=>d.onclick=e=>{ if(e.target.tagName==="INPUT") return; setActive(im,L[+d.dataset.i]); renderPanel(); });
  $a(el,".im-roilist .nm").forEach(n=>n.onchange=e=>{ L[+n.dataset.n].name=e.target.value; save(); draw(); });
  q("#msOne").onclick=()=>runOp("measure",{all:false}); q("#msAll").onclick=()=>runOp("measure",{all:true});
  q("#msLw").onchange=e=>{ S.ui.lw=Math.max(1,+e.target.value||1); save(); draw(); }; q("#msWt").onchange=e=>{ S.ui.wandTol=Math.max(0,+e.target.value||0); save(); };
  $a(el,"[data-col]").forEach(cb=>cb.onchange=()=>{ const k=cb.dataset.col; const s=new Set(S.measure.cols); if(cb.checked) s.add(k); else s.delete(k); S.measure.cols=ALLCOLS.map(x=>x[0]).filter(x=>s.has(x)); save(); drawBottom(); });
  $a(el,"[data-cls]").forEach(b=>b.onclick=()=>{ S.ui.ptCls=+b.dataset.cls; save(); renderPanel(); });
}
function pointCounter(r){ const n=new Array(8).fill(0); (r.cls||[]).forEach(c=>n[c]++); const cur=S.ui.ptCls||0;
  return `<div class="im-card"><div class="im-row"><b>Point counter</b><span class="xs muted">클릭 = 추가, Alt+클릭 = 지우기</span></div><div class="im-row">${n.map((v,k)=>`<button class="btn xs cls ${k===cur?"on":""}" data-cls="${k}" style="--c:${CLS_COL[k]}"><i></i>Type ${k+1}: <b>${v}</b></button>`).join("")}</div><div class="xs muted">합계 ${r.pts.length}</div></div>`; }
function roiSummary(im,r){ const c=calOf(im); if(r.type==="line") return `${g4(lineLen(im,r))} ${c.unit}`; if(r.type==="point") return `${r.pts.length} pts`; if(r.type==="angle") return `${roiAngle(r).toFixed(1)}°`; const b=roiBounds(r); return `${Math.round(b.w)}×${Math.round(b.h)}`; }
function roiAdd(im){ const a=activeRoi(im); if(!a){ note("ROI 가 없습니다"); return; } const L=roiList(im); if(!L.includes(a)){ a.name=a.name||String(L.length+1); L.push(a); } const n={...a}; void n; save(); renderPanel(); draw(); note(`ROI Manager 에 추가 (${L.length})`); }

/* ============================================================ 아래쪽 표 ============================================================ */
function drawBottom(){
  const el=ROOT&&$q(ROOT,"#imBottom"); if(!el) return; const tabs=[["results",`Results (${S.results.length})`],["particles","Particles"],["summary",`Summary (${S.summary.length})`],["dist","Distribution"],["skeleton",`Skeleton (${S.skelTrees.length})`],["plots",`Plots (${S.plots.length})`]];
  el.innerHTML=`<div class="im-btabs">${tabs.map(([k,l])=>`<button data-b="${k}" class="${S.ui.bottom===k?"active":""}">${l}</button>`).join("")}<span class="sp"></span><button class="btn xs" id="imBtog">${S.ui.bottomOpen?"접기":"펼치기"}</button></div><div class="im-bbody" ${S.ui.bottomOpen?"":"hidden"} id="imBbody"></div>`;
  $a(el,"[data-b]").forEach(b=>b.onclick=()=>{ S.ui.bottom=b.dataset.b; S.ui.bottomOpen=true; save(); drawBottom(); });
  $q(el,"#imBtog").onclick=()=>{ S.ui.bottomOpen=!S.ui.bottomOpen; save(); drawBottom(); };
  if(!S.ui.bottomOpen) return; const B=$q(el,"#imBbody");
  ({results:bResults,particles:bParticles,summary:bSummary,dist:bDist,skeleton:bSkeleton,plots:bPlots}[S.ui.bottom]||bResults)(B);
}
function tableHTML(rows,cols,o){ o=o||{}; if(!rows.length) return `<div class="xs muted im-empty">${o.empty||"표가 비었습니다"}</div>`;
  return `<div class="im-tbl"><table><thead><tr>${cols.map(c=>`<th>${esc(c.label||c.k)}${c.u?`<span class="u">${esc(c.u)}</span>`:""}</th>`).join("")}</tr></thead><tbody>${rows.map((r,i)=>`<tr data-i="${i}" class="${o.sel&&o.sel.has(i)?"sel":""}">${cols.map(c=>{ const v=r[c.k]; return `<td>${typeof v==="number"?g4(v):esc(v==null?"":v)}</td>`; }).join("")}</tr>`).join("")}</tbody></table></div>`; }
const csvQ=v=>{ const s=v==null?"":String(v); return /[",\n\r]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s; };
function csvOf(rows,cols){ return cols.map(c=>csvQ((c.label||c.k)+(c.u?` (${c.u})`:""))).join(",")+"\n"+rows.map(r=>cols.map(c=>csvQ(r[c.k])).join(",")).join("\n"); }
// 복사(Excel·Origin 붙여 넣기): 탭 구분, 값 안의 탭·줄바꿈은 공백 (쉼표는 그대로)
function tsvOf(rows,cols){ const t=v=>(v==null?"":String(v)).replace(/[\t\r\n]+/g," "); return cols.map(c=>t((c.label||c.k)+(c.u?` (${c.u})`:""))).join("\t")+"\n"+rows.map(r=>cols.map(c=>t(r[c.k])).join("\t")).join("\n"); }
function bResults(B){ const extra=["Length","Class"].filter(k=>S.results.some(r=>r[k]!=null));
  const units=[...new Set(S.results.map(r=>r.Unit).filter(Boolean))]; const one=units.length===1?units[0]:null; // 단위가 섞이면 Unit 열로
  const cols=[{k:"Label"},...(S.results.some(r=>r.Slice)?[{k:"Slice"}]:[]),...S.measure.cols.filter(k=>S.results.some(r=>r[k]!=null)).map(k=>({k,u:one?colUnit(k,one):""})),...extra.map(k=>({k,u:one?colUnit(k,one):""})),...(one?[]:[{k:"Unit"}])];
  B.innerHTML=`<div class="im-row"><button class="btn xs" id="rsCsv" ${S.results.length?"":"disabled"}>CSV</button><button class="btn xs" id="rsCopy" ${S.results.length?"":"disabled"}>복사 (표 붙여넣기)</button><button class="btn xs" id="rsClr" ${S.results.length?"":"disabled"}>Clear</button><span class="xs muted">열은 Measure ▸ Set Measurements 에서</span></div>${tableHTML(S.results,cols,{empty:"Measure (Ctrl+M) 결과가 여기에 쌓입니다."})}`;
  $q(B,"#rsCsv").onclick=()=>downloadText("Results.csv",csvOf(S.results,cols)); $q(B,"#rsCopy").onclick=()=>copyText(tsvOf(S.results,cols));
  $q(B,"#rsClr").onclick=()=>{ if(confirm2("rsclr","Results 를 비울까요? 한 번 더 누르면 비웁니다")){ S.results=[]; save(); drawBottom(); } }; }
function bParticles(B){ const im=cur(); const pa=im&&S.particles[im.id]; if(!pa){ B.innerHTML=`<div class="xs muted im-empty">Analyze ▸ Analyze Particles 를 실행하면 입자별 측정값이 여기에 나옵니다.</div>`; return; }
  const cols=[{k:"#"},...S.measure.cols.filter(k=>pa.rows.length&&pa.rows[0][k]!==undefined).map(k=>({k,u:colUnit(k,pa.unit)}))]; const sel=new Set(pa.sel||[]);
  B.innerHTML=`<div class="im-row"><b>${esc(im.name)}</b><span class="xs muted">${pa.rows.length} particles · ${esc(pa.summary.Threshold)}${pa.redirect&&pa.redirect!==im.name?` · 세기: ${esc(pa.redirect)}`:""}</span><span class="sp"></span><button class="btn xs" id="ptCsv">CSV</button><button class="btn xs" id="ptCopy">복사</button></div>${tableHTML(pa.rows,cols,{sel})}`;
  $q(B,"#ptCsv").onclick=()=>downloadText(`${im.name}_particles.csv`,csvOf(pa.rows,cols)); $q(B,"#ptCopy").onclick=()=>copyText(tsvOf(pa.rows,cols));
  $a(B,"tbody tr").forEach(tr=>tr.onclick=e=>{ const i=+tr.dataset.i; if(e.shiftKey||e.ctrlKey||e.metaKey){ const s=new Set(pa.sel); s.has(i)?s.delete(i):s.add(i); pa.sel=[...s]; } else pa.sel=[i];
    const c=pa.cent[i]; const v=viewOf(im); const {W,H}=canvasSize(); if(c){ v.x=c[0]-W/v.zoom/2; v.y=c[1]-H/v.zoom/2; } draw(); $a(B,"tbody tr").forEach(t=>t.classList.toggle("sel",pa.sel.includes(+t.dataset.i))); }); }
function bSkeleton(B){ const T=S.skelTrees, Bb=S.skelBranches; const units=[...new Set(T.map(r=>r.Unit))]; const one=units.length===1?units[0]:null;
  const tc=["Image","Slice","Skeleton","# Branches","# Junctions","# End-point voxels","# Junction voxels","# Slab voxels","Average Branch Length","# Triple points","# Quadruple points","Maximum Branch Length","Longest Shortest Path","spx","spy","spz"].filter(k=>T.some(r=>r[k]!==undefined&&r[k]!=="")).map(k=>({k,u:one&&/Length|Path|^sp/.test(k)?one:""})).concat(one?[]:[{k:"Unit"}]);
  const bc=["Image","Slice","Skeleton ID","Branch length","V1 x","V1 y","V1 z","V2 x","V2 y","V2 z","Euclidean distance","running average length"].filter(k=>Bb.some(r=>r[k]!==undefined&&r[k]!=="")).map(k=>({k,u:one&&/length|distance|^V/.test(k)?one:""})).concat(one?[]:[{k:"Unit"}]);
  B.innerHTML=`<div class="im-row"><button class="btn xs" id="skCsv" ${T.length?"":"disabled"}>골격 CSV</button><button class="btn xs" id="skBcsv" ${Bb.length?"":"disabled"}>가지 CSV</button><button class="btn xs" id="skClr" ${T.length?"":"disabled"}>Clear</button><span class="xs muted">열 이름 = Fiji Analyze Skeleton 의 Results · Branch information (Image·Slice·Unit 열은 앱에서 더함)</span></div>
    ${tableHTML(T,tc,{empty:"Analyze ▸ Analyze Skeleton 결과가 쌓입니다 (골격마다 한 줄)."})}${Bb.length?`<div class="im-h" style="margin-top:8px">Branch information <span class="muted xs">${Bb.length} 가지</span></div>${tableHTML(Bb.slice(0,2000),bc)}${Bb.length>2000?`<div class="xs muted">처음 2000 줄만 — 전부는 CSV</div>`:""}`:""}`;
  $q(B,"#skCsv").onclick=()=>downloadText("Skeleton.csv",csvOf(T,tc)); $q(B,"#skBcsv").onclick=()=>downloadText("Branch information.csv",csvOf(Bb,bc));
  $q(B,"#skClr").onclick=()=>{ if(confirm2("skclr","Skeleton 표를 비울까요? 한 번 더")){ S.skelTrees=[]; S.skelBranches=[]; save(); drawBottom(); } }; }
function bSummary(B){ const cols=["Image","Slice","Count","Total Area","Average Size","%Area","ECD mean","ECD SD","Feret mean","Circ. mean","Density (/area)","Unit","Threshold"].filter(k=>S.summary.some(r=>r[k]!==undefined&&r[k]!=="")).map(k=>({k}));
  B.innerHTML=`<div class="im-row"><button class="btn xs" id="smCsv" ${S.summary.length?"":"disabled"}>CSV</button><button class="btn xs" id="smClr" ${S.summary.length?"":"disabled"}>Clear</button><span class="xs muted">영상(장)마다 한 줄 — batch 결과도 여기에</span></div>${tableHTML(S.summary,cols,{empty:"Analyze Particles · Find Maxima 요약이 쌓입니다."})}`;
  $q(B,"#smCsv").onclick=()=>downloadText("Summary.csv",csvOf(S.summary,cols)); $q(B,"#smClr").onclick=()=>{ if(confirm2("smclr","Summary 를 비울까요? 한 번 더")){ S.summary=[]; save(); drawBottom(); } }; }
/* 크기 분포 */
function bDist(B){
  const im=cur(); const pa=im&&S.particles[im.id]; const D=S.ui.dist;
  const pool=S.ui.dist.pool==="all"?Object.values(S.particles).filter(p=>p.unit===(pa?pa.unit:"")):(pa?[pa]:[]);
  if(!pool.length){ B.innerHTML=`<div class="xs muted im-empty">입자 분석 결과가 없습니다.</div>`; return; }
  const col=D.col; const vals=[].concat(...pool.map(p=>p.rows.map(r=>r[col]))).filter(v=>v!=null&&isFinite(v)); const u=pa?((col==="Area")?(pa.unit==="px"?"px²":pa.unit+"²"):(LENCOL.has(col)?pa.unit:"")):"";
  const ss=IC.sizeSummary(vals);
  const bins=IC.binsFD(D.log?vals.filter(v=>v>0).map(Math.log10):vals,D.bins||0);
  B.innerHTML=`<div class="im-row"><label class="im-f"><span>Variable</span><select class="im-sel" id="dsCol">${ALLCOLS.filter(([k])=>pa&&pa.rows[0]&&pa.rows[0][k]!==undefined).map(([k,l])=>`<option value="${esc(k)}" ${k===col?"selected":""}>${esc(l)}</option>`).join("")}</select></label>
    <label class="im-f"><span>Bins</span><input type="number" id="dsBins" value="${D.bins||""}" placeholder="FD" min="3" step="1"></label><label class="im-chk"><input type="checkbox" id="dsLog" ${D.log?"checked":""}> log x</label>
    <label class="im-chk"><input type="checkbox" id="dsPool" ${D.pool==="all"?"checked":""}> 열린 모든 영상 합치기 (${Object.keys(S.particles).length})</label><span class="sp"></span><button class="btn xs" id="dsCsv">CSV</button><button class="btn xs" id="dsSvg">SVG</button></div>
  <div class="im-dist">${distSVG(vals,bins,ss,{log:D.log,xl:`${col}${u?" ("+u+")":""}`})}<table class="im-kv sm">${[["N",ss.n,""],["Mean ± SD",`${g4(ss.mean)} ± ${g4(ss.sd)}`,u],["SE of mean",g4(ss.se),u],["CV",fmt(ss.cv*100,1)+" %",""],["Median (D50)",g4(ss.median),u],["D10 / D90 (개수 기준)",`${g4(ss.d10)} / ${g4(ss.d90)}`,u],["Span (D90−D10)/D50",fmt(ss.span,3),""],["Min – Max",`${g4(ss.min)} – ${g4(ss.max)}`,u],
    ["Log-normal x̄g (기하평균)",g4(ss.gMean),u],["Log-normal σg (무차원)",g4(ss.gSD),""],["Log-normal mode",g4(ss.lnMode),u],["μ, σ of ln",`${fmt(ss.lnMu,4)}, ${fmt(ss.lnSigma,4)}`,""]].map(([k,v,uu])=>`<tr><td>${k}</td><td>${esc(v)}${uu?" "+esc(uu):""}</td></tr>`).join("")}</table></div>
  <div class="xs muted">Log-normal 매개변수는 ln(값)의 표본 평균·표준편차(n−1)로 구한 값${ss.nPos<ss.n?` — 0 이하 값 ${ss.n-ss.nPos}개는 로그정규에서 뺐습니다`:""}. 곡선 = 같은 N 과 칸 폭으로 맞춘 확률밀도. D10·D50·D90 은 입자 개수 기준(부피 기준 아님). 입자 수가 적으면(< 100) 분포 모양 판단을 삼가세요.</div>`;
  const q=s=>$q(B,s); q("#dsCol").onchange=e=>{ D.col=e.target.value; save(); drawBottom(); }; q("#dsBins").onchange=e=>{ D.bins=+e.target.value||0; save(); drawBottom(); }; q("#dsLog").onchange=e=>{ D.log=e.target.checked; save(); drawBottom(); };
  q("#dsPool").onchange=e=>{ D.pool=e.target.checked?"all":""; save(); drawBottom(); };
  q("#dsCsv").onclick=()=>downloadText(`distribution_${col}.csv`,"bin_lo,bin_hi,count\n"+bins.counts.map((c,i)=>{ const a=bins.lo+i*bins.w, b=a+bins.w; return (D.log?Math.pow(10,a):a)+","+(D.log?Math.pow(10,b):b)+","+c; }).join("\n"));
  q("#dsSvg").onclick=()=>{ const s=$q(B,".im-dist svg"); downloadText(`distribution_${col}.svg`,s.outerHTML,"image/svg+xml"); };
}
function distSVG(vals,bins,ss,o){ const W=640,H=260,ml=48,mr=12,mt=10,mb=38, pw=W-ml-mr, ph=H-mt-mb; if(!bins) return "";
  const cmax=Math.max(...bins.counts,1); const xOf=x=>ml+(x-bins.lo)/(bins.hi-bins.lo||1)*pw; const yOf=y=>mt+ph-(y/cmax/1.12)*ph;
  let bars=""; bins.counts.forEach((c,i)=>{ const x0=xOf(bins.lo+i*bins.w), x1=xOf(bins.lo+(i+1)*bins.w); bars+=`<rect x="${x0+0.5}" y="${yOf(c)}" width="${Math.max(0.5,x1-x0-1)}" height="${mt+ph-yOf(c)}" class="bar"><title>${c}</title></rect>`; });
  // 로그정규 pdf (칸 폭 × N)
  let curve=""; if(ss.n>2&&ss.lnSigma>0){ const pts=[]; for(let k=0;k<=160;k++){ const t=bins.lo+(bins.hi-bins.lo)*k/160; let dens;
      if(o.log){ const x=Math.pow(10,t); const lnx=Math.log(x); const pdfLn=Math.exp(-0.5*Math.pow((lnx-ss.lnMu)/ss.lnSigma,2))/(ss.lnSigma*Math.sqrt(2*Math.PI)); dens=pdfLn*Math.LN10*bins.w*ss.n; }
      else { const x=t; if(x<=0){ pts.push([xOf(t),yOf(0)]); continue; } const pdf=Math.exp(-0.5*Math.pow((Math.log(x)-ss.lnMu)/ss.lnSigma,2))/(x*ss.lnSigma*Math.sqrt(2*Math.PI)); dens=pdf*bins.w*ss.n; }
      pts.push([xOf(t),yOf(dens)]); } curve=`<polyline points="${pts.map(p=>p[0].toFixed(1)+","+Math.max(mt,p[1]).toFixed(1)).join(" ")}" class="fit"/>`; }
  const tk=niceTicksI(bins.lo,bins.hi,6).map(t=>`<line x1="${xOf(t)}" x2="${xOf(t)}" y1="${mt+ph}" y2="${mt+ph+4}" class="ax"/><text x="${xOf(t)}" y="${mt+ph+16}" text-anchor="middle">${o.log?g4(Math.pow(10,t)):g4(t)}</text>`).join("");
  const ty=niceTicksI(0,cmax*1.12,4).map(t=>`<line x1="${ml}" x2="${ml+pw}" y1="${yOf(t)}" y2="${yOf(t)}" class="gd"/><text x="${ml-6}" y="${yOf(t)+4}" text-anchor="end">${t}</text>`).join("");
  return `<svg class="im-chart" viewBox="0 0 ${W} ${H}" width="100%" style="max-width:${W*1.25}px" xmlns="http://www.w3.org/2000/svg">${ty}${bars}${curve}<line x1="${ml}" x2="${ml+pw}" y1="${mt+ph}" y2="${mt+ph}" class="ax"/>${tk}<text x="${ml+pw/2}" y="${H-4}" text-anchor="middle" class="al">${esc(o.xl)}</text><text transform="translate(12,${mt+ph/2}) rotate(-90)" text-anchor="middle" class="al">count</text></svg>`; }
function niceTicksI(a,b,n){ const span=b-a; if(!(span>0)) return [a]; const st0=span/(n||5); const p=Math.pow(10,Math.floor(Math.log10(st0))); const m=st0/p; const st=(m<1.5?1:m<3?2:m<7?5:10)*p; const out=[]; for(let t=Math.ceil(a/st)*st;t<=b+1e-9*st;t+=st) out.push(+t.toFixed(10)); return out; }
const PAL=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#6250d6"];
function svgXY(series,o){ o=o||{}; const W=o.w||640,H=o.h||220,ml=54,mr=12,mt=10,mb=36,pw=W-ml-mr,ph=H-mt-mb;
  let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity; series.forEach(s=>s.x.forEach((x,i)=>{ const y=s.y[i]; if(!isFinite(x)||!isFinite(y)) return; x0=Math.min(x0,x); x1=Math.max(x1,x); y0=Math.min(y0,y); y1=Math.max(y1,y); }));
  if(!isFinite(x0)){ x0=0; x1=1; y0=0; y1=1; } if(x1===x0) x1=x0+1; if(y1===y0){ y1=y0+1; } const pad=(y1-y0)*0.06; y0-=pad; y1+=pad;
  const X=x=>ml+(x-x0)/(x1-x0)*pw, Y=y=>mt+ph-(y-y0)/(y1-y0)*ph;
  const lines=series.map((s,k)=>`<polyline class="ln" style="stroke:${PAL[k%PAL.length]}" points="${s.x.map((x,i)=>isFinite(s.y[i])?X(x).toFixed(1)+","+Y(s.y[i]).toFixed(1):"").filter(Boolean).join(" ")}"/>`).join("");
  const tx=niceTicksI(x0,x1,6).map(t=>`<line x1="${X(t)}" x2="${X(t)}" y1="${mt}" y2="${mt+ph}" class="gd"/><text x="${X(t)}" y="${mt+ph+15}" text-anchor="middle">${g4(t)}</text>`).join("");
  const ty=niceTicksI(y0,y1,4).map(t=>`<line x1="${ml}" x2="${ml+pw}" y1="${Y(t)}" y2="${Y(t)}" class="gd"/><text x="${ml-6}" y="${Y(t)+4}" text-anchor="end">${g4(t)}</text>`).join("");
  const leg=series.length>1?series.map((s,k)=>`<tspan style="fill:${PAL[k%PAL.length]}">■ </tspan><tspan>${esc(s.name)}   </tspan>`).join(""):"";
  return `<svg class="im-chart" viewBox="0 0 ${W} ${H}" width="100%" style="max-width:${W*1.25}px" xmlns="http://www.w3.org/2000/svg">${tx}${ty}<rect x="${ml}" y="${mt}" width="${pw}" height="${ph}" class="fr"/>${lines}<text x="${ml+pw/2}" y="${H-3}" text-anchor="middle" class="al">${esc(o.xl||"")}</text><text transform="translate(12,${mt+ph/2}) rotate(-90)" text-anchor="middle" class="al">${esc(o.yl||"")}</text>${leg?`<text x="${ml+6}" y="${mt+13}" class="lg">${leg}</text>`:""}</svg>`; }
function peaksOf(x,y){ const out=[]; const n=y.length; let mx=-Infinity; for(const v of y) if(isFinite(v)&&v>mx) mx=v; let mn=Infinity; for(const v of y) if(isFinite(v)&&v<mn) mn=v;
  for(let i=2;i<n-2;i++){ if(y[i]>y[i-1]&&y[i]>=y[i+1]&&y[i]>y[i-2]&&y[i]>=y[i+2]){ const lo=Math.min(...y.slice(Math.max(0,i-6),i+7)); if(y[i]-lo>0.04*(mx-mn)) out.push(i); } } return out; }
function bPlots(B){ if(!S.plots.length){ B.innerHTML=`<div class="xs muted im-empty">Line profile · Radial profile · Orientation · Z-axis profile · Drift 그래프가 여기에 쌓입니다.</div>`; return; }
  B.innerHTML=S.plots.slice().reverse().map(p=>{ let extra=""; if(p.recip){ const s=p.series[0]; const pk=peaksOf(s.x,s.y).slice(0,12); if(pk.length) extra=`<div class="xs">피크: ${pk.map(i=>`r = ${g4(s.x[i])} → d = <b>${g4(1/s.x[i])}</b>`).join(" · ")}</div>`; }
    return `<div class="im-plot" data-id="${p.id}"><div class="im-row"><b class="sm">${esc(p.title)}</b><span class="sp"></span><button class="btn xs" data-csv="${p.id}">CSV</button><button class="btn xs" data-svg="${p.id}">SVG</button><button class="btn xs" data-x="${p.id}">×</button></div>${svgXY(p.series,{xl:p.xl,yl:p.yl})}${extra}</div>`; }).join("");
  $a(B,"[data-x]").forEach(b=>b.onclick=()=>{ S.plots=S.plots.filter(p=>p.id!==b.dataset.x); save(); drawBottom(); });
  $a(B,"[data-csv]").forEach(b=>b.onclick=()=>{ const p=S.plots.find(q=>q.id===b.dataset.csv); const s=p.series; downloadText(p.title.replace(/[^\w.-]+/g,"_")+".csv",[p.xl,...s.map(x=>x.name)].map(csvQ).join(",")+"\n"+s[0].x.map((x,i)=>[x,...s.map(q=>q.y[i])].join(",")).join("\n")); });
  $a(B,"[data-svg]").forEach(b=>b.onclick=()=>{ const s=$q(B,`.im-plot[data-id="${b.dataset.svg}"] svg`); downloadText("plot.svg",s.outerHTML,"image/svg+xml"); }); }

/* ---------- 내려받기·복사 ---------- */
function downloadBlob(name,blob){ const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },1500); }
function downloadText(name,text,mime){ downloadBlob(name,new Blob([mime&&mime.includes("svg")?text:"﻿"+text],{type:(mime||"text/csv")+";charset=utf-8"})); }
function copyText(t){ try{ navigator.clipboard.writeText(t).then(()=>note("복사했습니다 — Excel·Origin 에 붙여 넣으세요","ok")); }catch(e){ note("복사하지 못했습니다","warn"); } }
