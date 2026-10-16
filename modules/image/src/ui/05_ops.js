/* ============================================================
   연산 등록부 (한 곳): 이름·매개변수·실행·ImageJ 매크로 대응
   패널의 입력 칸, 매크로 기록·재생, batch, .ijm 내보내기가 모두 이 표를 쓴다.
   p: {k, label, type:"num"|"int"|"bool"|"sel"|"text", def, min, max, step, opts:[[v,label]], unit, hint}
   ============================================================ */
const OPS={};
const GROUPS=[["smooth","Smooth & denoise"],["rank","Rank & morphology (gray)"],["bg","Background"],["edge","Edges & sharpen"],["freq","Frequency (FFT)"],["contrast","Contrast"],["math","Math"],["geom","Geometry"],["type","Type"],["stack","Stack"],["calc","Image Calculator"]];
function op(id,o){ OPS[id]=Object.assign({id,params:[],kind:"plane"},o); }
const R=(k,label,def,o)=>Object.assign({k,label,type:"num",def},o||{});
/* 값 단위(valueUnit·siScale — AFM 높이 nm 등): 같은 물리량을 남기는 연산만 이어받는다. 화소를 바꾸는 나머지(형 변환의 눈금 바꾸기·문턱·FFT·곱셈 등)는 지운다.
   제자리 연산은 im.ver 가 바뀌었는지로, 새 영상은 돌려받은 영상에 붙인다. 되돌리기는 pushUndo 의 meta 가 되살린다. */
const KEEP_UNIT={gaussian:1,median:1,mean:1,despeckle:1,outliers:1,smooth3:1,minimum:1,maximum:1,tophat:1,subtractBg:1,sharpen:1,unsharp:1,
  convolve:P=>{ if(!P.normalize) return false; const v=String(P.kernel).trim().split(/[\s,]+/).map(Number); return Math.abs(v.reduce((a,b)=>a+b,0))>1e-12; }, /* 합으로 나눈 커널 = 가중 평균 */
  bandpass:P=>!P.autoscale&&!P.saturate, math:P=>["add","sub","min","max","set","abs"].includes(P.fn), nanbg:1, invert:1,
  crop:1, cropbar:1, rot90:1, flip:1, rotate:1, resize:1, duplicate:1, register:1,
  type:P=>P.to==="32"||P.to==="split"||((P.to==="8"||P.to==="16")&&!P.scale), zproject:P=>["avg","max","min","median","sd"].includes(P.method)};
async function runWithUnit(o,im,P){
  const u={valueUnit:im.valueUnit||"",siScale:im.siScale||null,siUnit:im.siUnit||""}, ver0=im.ver; const r=await o.run(im,P);
  if(!u.valueUnit&&!u.siScale) return r;
  const k=KEEP_UNIT[o.id]; const keep=typeof k==="function"?!!k(P,im):!!k;
  if(!keep&&im.ver!==ver0){ im.valueUnit=""; im.siScale=null; im.siUnit=""; }
  if(keep&&r&&r.id&&r.id!==im.id&&!r.valueUnit) Object.assign(r,u);
  return r; }

/* 선택한 장(현재 / 전체)에 core 함수 적용 → 새 planes */
async function mapPlanes(im,fn,argsOf,all,label){
  const P=planesOf(im).slice(); const idx=[];
  if(all){ for(let k=0;k<P.length;k++) idx.push(k); }
  else if(im.rgb){ for(let c=0;c<3;c++) idx.push(pIndex(im,c,im.z,im.t)); }
  else idx.push(pIndex(im,im.c,im.z,im.t));
  let done=0, outBit=im.bit; const plBit=im.rgb?8:im.bit;
  setProg(label+" …",0.02);
  const jobs=idx.map(k=>call(fn,argsOf({w:im.w,h:im.h,bit:plBit,data:P[k]}),f=>{ if(idx.length===1) setProg(label+" …",f); }).then(r=>{ const pl=r.plane||r; P[k]=pl.data; outBit=pl.bit; done++; if(idx.length>1) setProg(`${label} ${done}/${idx.length}`,done/idx.length); }));
  await Promise.all(jobs); hideProg();
  return {planes:P,bit:im.rgb?im.bit:outBit};
}
async function runPlaneOp(im,P,fn,argsOf,label){
  if(P.__pv) return await mapPlanes(im,fn,argsOf,false,label); // 미리 보기: 현재 장만 계산해 돌려주고 영상은 그대로
  const all=!!P.stack; const r=await mapPlanes(im,fn,argsOf,all||(!isStack(im)),label);
  setPlanes(im,r.planes,r.bit!==im.bit?{bit:r.bit}:null,label,{op:P.__op,params:clean(P),label});
}
const clean=P=>{ const o={}; for(const k in P) if(k!=="__op"&&k!=="__pv") o[k]=P[k]; return o; };
const stackParam={k:"stack",label:"Process all slices",type:"bool",def:false,stackOnly:true};
const ijStack=P=>P.stack?" stack":"";

/* ---------- Smooth & denoise ---------- */
op("gaussian",{label:"Gaussian Blur",group:"smooth",params:[R("sigma","Sigma (radius)",2,{min:0.1,step:0.1,unit:"px"}),R("scaled","Scaled units",false,{type:"bool"}),stackParam],
  run:(im,P)=>{ const sx=P.scaled&&im.cal?+P.sigma/calOf(im).pw:+P.sigma, sy=P.scaled&&im.cal?+P.sigma/calOf(im).ph:+P.sigma; return runPlaneOp(im,P,"gaussian",p=>[p,{sigma:sx,sigmaY:sy,rgb:im.rgb}],`Gaussian Blur σ=${P.sigma}`); },
  ij:P=>`run("Gaussian Blur...", "sigma=${P.sigma}${P.scaled?" scaled":""}${ijStack(P)}");`,
  help:"ImageJ 와 같은 커널(정확도 8-bit 0.002, 16/32-bit 0.0002)과 가장자리 처리(가장자리 화소 복제). σ > 8.5 px 는 축소-합성곱-확대 경로. Scaled units 를 켜면 Sigma 를 보정 단위(예: µm)로 읽습니다."});
op("median",{label:"Median",group:"smooth",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"median",{radius:+P.radius}],`Median r=${P.radius}`),ij:P=>`run("Median...", "radius=${P.radius}${ijStack(P)}");`,help:"원형 커널 (ImageJ RankFilters 와 같은 모양)."});
op("mean",{label:"Mean",group:"smooth",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"mean",{radius:+P.radius}],`Mean r=${P.radius}`),ij:P=>`run("Mean...", "radius=${P.radius}${ijStack(P)}");`});
op("despeckle",{label:"Despeckle",group:"smooth",params:[stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"despeckle",{}],"Despeckle"),ij:P=>`run("Despeckle"${P.stack?', "stack"':', " "'});`,help:"Median r=1."});
op("outliers",{label:"Remove Outliers",group:"smooth",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),R("threshold","Threshold",50,{step:1}),R("which","Which",["bright"],{type:"sel",opts:[["bright","Bright"],["dark","Dark"]],def:"bright"}),stackParam],
  run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"outliers",{radius:+P.radius,threshold:+P.threshold,which:P.which}],`Remove Outliers (${P.which})`),ij:P=>`run("Remove Outliers...", "radius=${P.radius} threshold=${P.threshold} which=${P.which==="dark"?"Dark":"Bright"}${ijStack(P)}");`,
  help:"중앙값과의 차이가 threshold 보다 큰 화소를 중앙값으로 바꿉니다."});
op("smooth3",{label:"Smooth (3×3)",group:"smooth",params:[stackParam],run:(im,P)=>runPlaneOp(im,P,"filter3",p=>[p,"smooth"],"Smooth"),ij:P=>`run("Smooth"${P.stack?', "stack"':', " "'});`});
/* ---------- Rank & morphology ---------- */
op("minimum",{label:"Minimum",group:"rank",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"min",{radius:+P.radius}],`Minimum r=${P.radius}`),ij:P=>`run("Minimum...", "radius=${P.radius}${ijStack(P)}");`,help:"회색조 erosion."});
op("maximum",{label:"Maximum",group:"rank",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"max",{radius:+P.radius}],`Maximum r=${P.radius}`),ij:P=>`run("Maximum...", "radius=${P.radius}${ijStack(P)}");`,help:"회색조 dilation."});
op("variance",{label:"Variance",group:"rank",params:[R("radius","Radius",2,{min:0.5,step:0.5,unit:"px"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"variance",{radius:+P.radius}],`Variance r=${P.radius}`),ij:P=>`run("Variance...", "radius=${P.radius}${ijStack(P)}");`,help:"지역 분산(모집단). 8/16-bit 는 최댓값에서 잘립니다 — 32-bit 로 바꾼 뒤 쓰는 것이 좋습니다."});
op("tophat",{label:"Top Hat",group:"rank",params:[R("radius","Radius",5,{min:0.5,step:0.5,unit:"px"}),R("light","Light background",false,{type:"bool"}),stackParam],run:(im,P)=>runPlaneOp(im,P,"rank",p=>[p,"tophat",{radius:+P.radius,light:!!P.light}],`Top Hat r=${P.radius}`),ij:P=>`run("Top Hat...", "radius=${P.radius}${P.light?" light":""}${ijStack(P)}");`,help:"원본 − opening. 반지름보다 작은 밝은 구조만 남깁니다."});
/* ---------- Background ---------- */
op("subtractBg",{label:"Subtract Background",group:"bg",params:[R("radius","Rolling ball radius",50,{min:1,step:1,unit:"px"}),R("light","Light background",false,{type:"bool"}),R("create","Create background (don't subtract)",false,{type:"bool"}),R("nosmooth","Disable smoothing",false,{type:"bool"}),stackParam],
  run:(im,P)=>runPlaneOp(im,P,"subtractBackground",p=>[p,{radius:+P.radius,light:!!P.light,createBackground:!!P.create,smooth:!P.nosmooth}],`Subtract Background r=${P.radius}`),
  ij:P=>`run("Subtract Background...", "rolling=${P.radius}${P.light?" light":""}${P.create?" create":""}${P.nosmooth?" disable":""}${ijStack(P)}");`,help:"반지름은 배경에 속하지 않는 가장 큰 물체의 반지름보다 크게. 반지름 > 10 px 는 ImageJ 처럼 영상을 줄여 계산합니다."});
/* ---------- Edges ---------- */
op("edges",{label:"Find Edges (Sobel)",group:"edge",params:[stackParam],run:(im,P)=>runPlaneOp(im,P,"filter3",p=>[p,"edges"],"Find Edges"),ij:P=>`run("Find Edges"${P.stack?', "stack"':', " "'});`});
op("sharpen",{label:"Sharpen (3×3)",group:"edge",params:[stackParam],run:(im,P)=>runPlaneOp(im,P,"filter3",p=>[p,"sharpen"],"Sharpen"),ij:P=>`run("Sharpen"${P.stack?', "stack"':', " "'});`});
op("unsharp",{label:"Unsharp Mask",group:"edge",params:[R("sigma","Radius (sigma)",2,{min:0.1,step:0.1,unit:"px"}),R("weight","Mask weight",0.6,{min:0.1,max:0.9,step:0.05}),stackParam],run:(im,P)=>runPlaneOp(im,P,"unsharp",p=>[p,{sigma:+P.sigma,weight:+P.weight}],`Unsharp Mask σ=${P.sigma}`),ij:P=>`run("Unsharp Mask...", "radius=${P.sigma} mask=${(+P.weight).toFixed(2)}${ijStack(P)}");`});
op("convolve",{label:"Convolve",group:"edge",params:[R("kernel","Kernel (one row per line)","-1 -1 -1\n-1 8 -1\n-1 -1 -1",{type:"text"}),R("normalize","Normalize kernel",true,{type:"bool"}),stackParam],
  run:(im,P)=>{ const rows=String(P.kernel).trim().split(/\n+/).map(l=>l.trim().split(/[\s,]+/).map(Number)); const kh=rows.length, kw=rows[0].length; if(rows.some(r=>r.length!==kw||r.some(v=>!isFinite(v)))) throw new Error("커널의 각 행은 같은 개수의 숫자여야 합니다");
    return runPlaneOp(im,P,"convolve",p=>[p,{kw,kh,k:[].concat(...rows),normalize:!!P.normalize}],"Convolve"); },
  ij:P=>`run("Convolve...", "text1=[${String(P.kernel).trim().replace(/\n/g,"\\n")}\\n]${P.normalize?" normalize":""}${ijStack(P)}");`});
/* ---------- Frequency ---------- */
op("bandpass",{label:"Bandpass Filter",group:"freq",params:[R("large","Filter large structures down to",40,{min:1,unit:"px"}),R("small","Filter small structures up to",3,{min:0,unit:"px"}),R("stripes","Suppress stripes","none",{type:"sel",opts:[["none","None"],["horizontal","Horizontal"],["vertical","Vertical"]]}),R("tolerance","Tolerance of direction",5,{unit:"%"}),R("autoscale","Autoscale after filtering",true,{type:"bool"}),R("saturate","Saturate image when autoscaling",true,{type:"bool"}),stackParam],
  run:(im,P)=>runPlaneOp(im,P,"bandpass",p=>[p,{large:+P.large,small:+P.small,stripes:P.stripes,tolerance:+P.tolerance,autoscale:!!P.autoscale,saturate:!!P.saturate}],`Bandpass ${P.small}–${P.large} px`),
  ij:P=>`run("Bandpass Filter...", "filter_large=${P.large} filter_small=${P.small} suppress=${P.stripes==="none"?"None":P.stripes==="horizontal"?"Horizontal":"Vertical"} tolerance=${P.tolerance}${P.autoscale?" autoscale":""}${P.saturate?" saturate":""}${P.stack?" process":""}");`,
  help:"거울상으로 2ⁿ 크기에 채운 뒤 주파수 영역 Gaussian 대역 통과 (ImageJ FFTFilter 정의). DC(평균)는 그대로."});
op("fft",{label:"FFT",group:"freq",kind:"new",params:[R("window","Window","none",{type:"sel",opts:[["none","None"],["hann","Hann"]]})],
  run:async(im,P)=>{ const p=curPlaneOf(im); const r=await call("powerSpectrum",[p,{window:P.window}]); const cs=calOf(im); const n=addImage({name:"FFT of "+im.name,w:r.w,h:r.h,bit:8,planes:[r.data],cal:im.cal?{pw:1/(r.w*cs.pw),ph:1/(r.h*cs.ph),unit:"1/"+cs.unit,src:`FFT: 1/(${r.w}·화소 크기) — 역공간`}:null},{parent:im.id,src:"derived",kind:"fft",fft:{of:im.id,S:r.w,pad:r.pad,pw:calOf(im).pw,unit:calOf(im).unit,c:im.c,z:im.z,t:im.t}});
    n.info=[["Description","ln(1+|F|) 를 0–255 로. 가운데 = 0 주파수"],["Frequency step",im.cal?`1/(${r.w}·${g4(calOf(im).pw)} ${calOf(im).unit})`:`1/${r.w} px⁻¹`]]; record({op:"fft",params:clean(P),label:"FFT"},im); return n; },
  ij:()=>`run("FFT");`,help:"표시용 power spectrum. 그 위에 ROI 를 그리고 Inverse FFT 로 주기 잡음 제거·격자 필터를 할 수 있습니다 (Process ▸ Frequency ▸ Inverse FFT)."});
op("ifft",{label:"Inverse FFT",group:"freq",kind:"new",params:[R("mode","Inside the ROI","keep",{type:"sel",opts:[["keep","Pass (clear the rest)"],["remove","Filter out (pass the rest)"]]}),R("sym","Add symmetric ROI about origin",true,{type:"bool"})],
  run:async(im,P)=>{ if(im.kind!=="fft"||!im.fft) throw new Error("FFT power spectrum 영상에서 실행하세요"); const src=byId(im.fft.of); if(!src) throw new Error("원본 영상이 닫혔습니다");
    const Sx=im.fft.S; const m=roisMask(im,true); if(!m) throw new Error("FFT 영상 위에 ROI(사각형·타원·다각형)를 그리세요");
    if(P.sym){ for(let y=0;y<Sx;y++) for(let x=0;x<Sx;x++){ const xs=(Sx-x)%Sx, ys=(Sx-y)%Sx; if(m[y*Sx+x]) m[ys*Sx+xs]=1; } }
    const p=planeAt(src,im.fft.c,im.fft.z,im.fft.t); const r=await call("fftMaskFilter",[p,m,Sx,P.mode==="keep"]);
    const n=addImage({name:"Inverse FFT of "+src.name,w:r.w,h:r.h,bit:32,planes:[r.data],cal:src.cal},{parent:src.id,src:"derived"}); return n; },
  ij:P=>`run("Inverse FFT");`,help:"FFT 영상 위에 그린 ROI 를 마스크로 역변환합니다 (주기 잡음 제거·격자 필터)."});
/* ---------- Contrast ---------- */
op("enhance",{label:"Enhance Contrast",group:"contrast",params:[R("saturated","Saturated pixels",0.35,{min:0,max:50,step:0.05,unit:"%"}),R("normalize","Normalize",false,{type:"bool"}),R("equalize","Equalize histogram",false,{type:"bool"}),R("classic","Classic equalization",false,{type:"bool"}),stackParam],
  run:async(im,P)=>{ if(P.equalize){ if(im.bit===32&&!im.rgb) throw new Error("32-bit 은 Equalize 하지 않습니다 (ImageJ 와 같음)"); return runPlaneOp(im,P,"equalize",p=>[p,!!P.classic],"Equalize"); }
    if(!P.normalize){ if(P.__pv) return null; const p=curPlaneOf(im); const r=await call("enhanceContrast",[p,{saturated:+P.saturated}]); if(!r.range){ note("영상이 거의 일정해 표시 범위를 바꾸지 않았습니다 (ImageJ 와 같음)","warn"); return; } const c=im.rgb?0:im.c; im.disp[c].min=r.range[0]; im.disp[c].max=r.range[1]; invalidate(); save(); rerender(); /* 표시만 바뀜 — ver 를 올리면 저장된 화소가 '오래됨' 으로 버려짐 */ note(`표시 범위 ${g4(r.range[0])} – ${g4(r.range[1])} (화소 값은 그대로)`); return; }
    return runPlaneOp(im,P,"enhanceContrast",p=>[p,{saturated:+P.saturated,normalize:true}],`Normalize ${P.saturated}%`); },
  ij:P=>P.equalize?`${P.classic?'setKeyDown("alt"); ':""}run("Enhance Contrast...", "saturated=${P.saturated} equalize${P.stack?" process_all":""}");${P.classic?' setKeyDown("none");':""}`:`run("Enhance Contrast...", "saturated=${P.saturated}${P.normalize?" normalize":""}${P.stack&&P.normalize?" process_all":""}");`,
  help:"Normalize 를 끄면 표시 범위만 바꿉니다(측정값 그대로). 세기 비교에는 Normalize·Equalize 를 쓰지 마세요 (둘 다 화소 값을 바꿉니다). Equalize 는 ImageJ 기본이 제곱근 방식이고, Classic 을 켜면 원래 방식입니다."});
op("applyLut",{label:"Apply LUT",group:"contrast",params:[stackParam],
  run:(im,P)=>{ if(im.rgb||im.bit===32) throw new Error("8/16-bit 영상에서만 (ImageJ 와 같음)"); const d=im.disp[im.c]; if(P.min==null||P.max==null){ P.min=d.min; P.max=d.max; } // 기록 = 그때의 표시 범위 (재생도 같게)
    const mn=+P.min, mx=+P.max; return runPlaneOp(im,P,"applyLut",p=>[p,mn,mx],`Apply LUT ${g4(mn)}–${g4(mx)}`); },
  ij:P=>`${P.min!=null?`setMinAndMax(${P.min}, ${P.max});\n`:""}run("Apply LUT"${P.stack?', "stack"':', " "'});`,hidden:true});
/* ---------- Math ---------- */
const MATHS=[["add","Add",25],["sub","Subtract",25],["mul","Multiply",1.25],["div","Divide",1.25],["min","Min",0],["max","Max",255],["set","Set",0],["gamma","Gamma",0.5],["log","Log"],["exp","Exp"],["sqr","Square"],["sqrt","Square Root"],["abs","Abs"],["reciprocal","Reciprocal"]];
op("math",{label:"Math",group:"math",params:[R("fn","Operation","add",{type:"sel",opts:MATHS.map(m=>[m[0],m[1]])}),R("value","Value",25,{step:"any",showIf:P=>(MATHS.find(m=>m[0]===P.fn)||[])[2]!==undefined}),stackParam],
  run:(im,P)=>runPlaneOp(im,P,"math",p=>[p,P.fn,+P.value],`${(MATHS.find(m=>m[0]===P.fn)||[0,P.fn])[1]}${MATHS.find(m=>m[0]===P.fn)[2]!==undefined?" "+P.value:""}`),
  ij:P=>{ const m=MATHS.find(x=>x[0]===P.fn); const nm=m[1]; return m[2]!==undefined?`run("${nm}...", "value=${P.value}${ijStack(P)}");`:`run("${nm}"${P.stack?', "stack"':', " "'});`; },
  help:"8/16-bit 은 ImageJ 정수 규칙(반올림·범위 자르기). 소수 결과가 필요하면 먼저 32-bit 로."});
op("nanbg",{label:"NaN Background",group:"math",params:[R("lower","Lower",0),R("upper","Upper",255)],help:"32-bit 영상에서 Lower–Upper 범위 밖 화소를 NaN 으로 바꿉니다.",run:(im,P)=>{ if(im.bit!==32) throw new Error("32-bit 영상에서만"); return runPlaneOp(im,P,"math",p=>[p,"nan",[+P.lower,+P.upper]],"NaN Background"); },ij:P=>`setThreshold(${P.lower}, ${P.upper}); run("NaN Background");`});
op("invert",{label:"Invert",group:"math",params:[stackParam],run:(im,P)=>runPlaneOp(im,P,"invert",p=>[p],"Invert"),ij:P=>`run("Invert"${P.stack?', "stack"':', " "'});`,help:"8-bit: 255−v. 16/32-bit: (min+max)−v (화소 범위)."});
/* ---------- Geometry ---------- */
op("crop",{label:"Crop",group:"geom",kind:"image",params:[],help:"사각형 ROI 를 먼저 그린 뒤 그 영역으로 자릅니다.",
  run:(im,P)=>{ const r=P.rect||roiRect(im); if(!r) throw new Error("사각형 ROI 를 먼저 그리세요"); const rr={x:clamp(Math.round(r.x),0,im.w-1),y:clamp(Math.round(r.y),0,im.h-1)}; rr.w=clamp(Math.round(r.w),1,im.w-rr.x); rr.h=clamp(Math.round(r.h),1,im.h-rr.y);
    const P2=planesOf(im).map(d=>IC.crop({w:im.w,h:im.h,bit:im.bit,data:d},rr).data); clearRois(im); setPlanes(im,P2,{w:rr.w,h:rr.h},"Crop",{op:"crop",params:{rect:rr},label:`Crop ${rr.w}×${rr.h}`}); resetView(im); },
  ij:P=>`makeRectangle(${P.rect.x}, ${P.rect.y}, ${P.rect.w}, ${P.rect.h}); run("Crop");`});
op("cropbar",{label:"Crop data bar",group:"geom",kind:"image",hidden:true,params:[R("rows","Rows at bottom",0,{type:"int"})],
  run:(im,P)=>{ const n=+P.rows; if(!(n>0&&n<im.h)) throw new Error("줄 수가 올바르지 않습니다"); return OPS.crop.run(im,{rect:{x:0,y:0,w:im.w,h:im.h-n}}); },ij:P=>`makeRectangle(0, 0, getWidth(), getHeight()-${P.rows}); run("Crop");`});
op("rot90",{label:"Rotate 90 Degrees",group:"geom",kind:"image",params:[R("dir","Direction","cw",{type:"sel",opts:[["cw","Right"],["ccw","Left"]]})],
  run:(im,P)=>{ const cw=P.dir==="cw"; const P2=planesOf(im).map(d=>IC.rotate90({w:im.w,h:im.h,bit:im.bit,data:d},cw).data); const c=im.cal; clearRois(im); setPlanes(im,P2,{w:im.h,h:im.w,cal:c?{...c,pw:c.ph,ph:c.pw}:c},"Rotate 90",{op:"rot90",params:clean(P),label:`Rotate 90 ${cw?"Right":"Left"}`}); resetView(im); },
  ij:P=>`run("Rotate 90 Degrees ${P.dir==="cw"?"Right":"Left"}");`});
op("flip",{label:"Flip",group:"geom",kind:"image",params:[R("dir","Direction","h",{type:"sel",opts:[["h","Horizontally"],["v","Vertically"]]})],
  run:(im,P)=>{ const P2=planesOf(im).map(d=>IC.flip({w:im.w,h:im.h,bit:im.bit,data:d},P.dir).data); clearRois(im); setPlanes(im,P2,{keepDisplay:true},"Flip",{op:"flip",params:clean(P),label:"Flip "+(P.dir==="h"?"Horizontally":"Vertically")}); },
  ij:P=>`run("Flip ${P.dir==="h"?"Horizontally":"Vertically"}", "stack");`}); // 앱은 모든 장을 뒤집음
op("rotate",{label:"Rotate",group:"geom",kind:"image",params:[R("angle","Angle (clockwise +)",15,{unit:"°",step:0.1}),R("enlarge","Enlarge to fit",false,{type:"bool"}),R("interp","Interpolation","bilinear",{type:"sel",opts:[["bilinear","Bilinear"],["nearest","None"]]})],
  run:async(im,P)=>{ const outs=await Promise.all(planesOf(im).map(d=>call("rotate",[{w:im.w,h:im.h,bit:im.bit,data:d},+P.angle,{enlarge:!!P.enlarge,nearest:P.interp==="nearest"}]))); clearRois(im); setPlanes(im,outs.map(o=>o.data),{w:outs[0].w,h:outs[0].h},"Rotate",{op:"rotate",params:clean(P),label:`Rotate ${P.angle}°`}); resetView(im); },
  ij:P=>`run("Rotate... ", "angle=${P.angle} grid=1 interpolation=${P.interp==="nearest"?"None":"Bilinear"}${P.enlarge?" enlarge":""} stack");`});
op("resize",{label:"Scale",group:"geom",kind:"image",params:[R("factor","Scale factor",0.5,{min:0.01,step:0.05}),R("interp","Interpolation","bilinear",{type:"sel",opts:[["bilinear","Bilinear (average when downsizing)"],["nearest","None"]]})],
  run:async(im,P)=>{ const f=+P.factor; const W=Math.max(1,Math.round(im.w*f)), H=Math.max(1,Math.round(im.h*f)); const outs=await Promise.all(planesOf(im).map(d=>call("resize",[{w:im.w,h:im.h,bit:im.bit,data:d},W,H,{method:P.interp}])));
    const c=im.cal; clearRois(im); setPlanes(im,outs.map(o=>o.data),{w:W,h:H,cal:c?{...c,pw:c.pw*im.w/W,ph:c.ph*im.h/H}:c},"Scale",{op:"resize",params:clean(P),label:`Scale ×${f}`}); resetView(im); },
  ij:P=>`run("Scale...", "x=${P.factor} y=${P.factor} interpolation=${P.interp==="nearest"?"None":"Bilinear"} average process");`,help:"보정(화소 크기)은 배율만큼 함께 바뀝니다."});
op("duplicate",{label:"Duplicate",group:"geom",kind:"new",params:[R("slice","Current slice only",false,{type:"bool"})],
  run:(im,P)=>{ let planes, nz=im.nz, nt=im.nt, nc=im.nc; if(P.slice&&isStack(im)){ planes=[]; for(let c=0;c<im.nc;c++) planes.push(Float32Array.from(planesOf(im)[pIndex(im,c,im.z,im.t)])); nz=1; nt=1; } else planes=planesOf(im).map(p=>Float32Array.from(p));
    const n=addImage({name:im.name+"-1",w:im.w,h:im.h,bit:im.bit,rgb:im.rgb,nc:im.rgb?1:nc,nz,nt,planes,cal:im.cal,chNames:im.chNames,zStep:nz>1?im.zStep:0},{parent:im.id,src:"derived",ops:im.ops.slice()}); n.disp=JSON.parse(JSON.stringify(im.disp));
    record({op:"duplicate",params:clean(P),label:"Duplicate"+(P.slice?" (slice)":"")},im); return n; },
  ij:P=>P.slice?`run("Duplicate...", " ");`:`run("Duplicate...", "duplicate");`});
/* ---------- Type ---------- */
op("type",{label:"Type",group:"type",kind:"image",params:[R("to","Convert to","8",{type:"sel",opts:[["8","8-bit"],["16","16-bit"],["32","32-bit"],["rgb2gray","RGB → 8-bit gray"],["split","Split channels"]]}),R("weighted","Weighted RGB (0.299/0.587/0.114)",false,{type:"bool",showIf:P=>P.to==="rgb2gray"}),R("scale","Scale when converting",true,{type:"bool",showIf:P=>P.to==="8"||P.to==="16"})],
  run:(im,P)=>{
    if(P.to==="split"){ if(im.nc<2) throw new Error("채널이 하나뿐입니다"); const names=im.rgb?["red","green","blue"]:im.chNames; let last;
      for(let c=0;c<im.nc;c++){ const planes=[]; for(let t=0;t<im.nt;t++) for(let z=0;z<im.nz;z++) planes.push(Float32Array.from(planesOf(im)[pIndex(im,c,z,t)]));
        last=addImage({name:`C${c+1}-${im.name}${names&&names[c]?" ("+names[c]+")":""}`,w:im.w,h:im.h,bit:im.rgb?8:im.bit,nc:1,nz:im.nz,nt:im.nt,planes,cal:im.cal,zStep:im.zStep,valueUnit:im.valueUnit,siScale:im.siScale,siUnit:im.siUnit},{parent:im.id,src:"derived"}); if(!im.rgb&&im.disp[c]) last.disp[0].lut=im.disp[c].lut; }
      record({op:"type",params:clean(P),label:"Split Channels"},im); return last; }
    if(P.to==="rgb2gray"){ if(!im.rgb) throw new Error("RGB 영상이 아닙니다"); const out=[]; for(let t=0;t<im.nt;t++) for(let z=0;z<im.nz;z++){ const r=planeAt(im,0,z,t), g=planeAt(im,1,z,t), b=planeAt(im,2,z,t); out.push(IC.rgbToGray(r,g,b,!!P.weighted).data); }
      setPlanes(im,out,{rgb:false,nc:1,bit:8,c:0,composite:false},"RGB → 8-bit",{op:"type",params:clean(P),label:"RGB → 8-bit"}); im.disp=[{min:0,max:255,lut:"Grays"}]; return; }
    if(im.rgb) throw new Error("RGB 는 먼저 'RGB → 8-bit gray' 또는 'Split channels'");
    const to=+P.to; if(to===im.bit){ note("이미 "+to+"-bit 입니다"); return; }
    // 비례 변환 범위: 기록된 값이 있으면 그대로(매크로 재생), 없으면 지금 표시 범위 — 화소 min–max 와 같으면 "minmax" 로 기록해 영상마다 다시 구함
    const mm=c=>{ let a=Infinity,b=-Infinity; for(let t=0;t<im.nt;t++) for(let z=0;z<im.nz;z++){ const r=IC.minMax(planesOf(im)[pIndex(im,c,z,t)]); if(r[0]<a) a=r[0]; if(r[1]>b) b=r[1]; } return [a,b]; };
    if(P.range===undefined){ const disp=im.disp.map(d=>[d.min,d.max]); P.range=disp.every((r,c)=>{ const m=mm(c); return Math.abs(r[0]-m[0])<1e-9*(1+Math.abs(m[0]))&&Math.abs(r[1]-m[1])<1e-9*(1+Math.abs(m[1])); })?"minmax":disp; }
    const rng=c=>P.range==="minmax"||!Array.isArray(P.range)?mm(c):(P.range[c]||P.range[0]);
    const R0=[...Array(im.nc).keys()].map(rng);
    const out=planesOf(im).map((d,k)=>{ const c=k%im.nc; return IC.convertBit({w:im.w,h:im.h,bit:im.bit,data:d},to,R0[c][0],R0[c][1],P.scale!==false).data; });
    setPlanes(im,out,{bit:to},`${to}-bit`,{op:"type",params:clean(P),label:`${to}-bit`}); },
  ij:P=>P.to==="split"?`run("Split Channels");`:P.to==="rgb2gray"?`run("Conversions...", "scale${P.weighted?" weighted":""}"); run("8-bit");`:`${Array.isArray(P.range)?`setMinAndMax(${P.range[0][0]}, ${P.range[0][1]}); `:"resetMinAndMax(); "}run("Conversions...", "${P.scale!==false?"scale":" "}"); run("${P.to}-bit");`,
  help:"16/32 → 8-bit 은 현재 표시 범위로 비례 변환합니다(ImageJ 와 같음). 측정 전에 정보 손실을 확인하세요."});
/* ---------- Stack ---------- */
op("zproject",{label:"Z Project",group:"stack",kind:"new",params:[R("method","Projection","avg",{type:"sel",opts:[["avg","Average Intensity"],["max","Max Intensity"],["min","Min Intensity"],["sum","Sum Slices"],["sd","Standard Deviation"],["median","Median"]]}),R("float","32-bit result",false,{type:"bool"})],
  run:async(im,P)=>{ if(!isStack(im)) throw new Error("스택이 아닙니다");
    // hyperstack: z 가 있으면 시간(t)마다 z 를 투영(ImageJ 와 같음), z 가 하나면 t 를 투영
    const overZ=im.nz>1, nOut=overZ?im.nt:1, out=[]; let bit=im.bit;
    const nch=im.rgb?3:im.nc;
    for(let to=0;to<nOut;to++) for(let c=0;c<nch;c++){ const pl=[];
      if(overZ) for(let z=0;z<im.nz;z++) pl.push({w:im.w,h:im.h,bit:im.rgb?8:im.bit,data:planesOf(im)[pIndex(im,c,z,to)]});
      else for(let t=0;t<im.nt;t++) pl.push({w:im.w,h:im.h,bit:im.rgb?8:im.bit,data:planesOf(im)[pIndex(im,c,0,t)]});
      const src=P.float?pl.map(p=>({...p,bit:32})):pl; const r=await call("zProject",[src,P.method]); out[c+nch*to]=r.data; bit=r.bit; }
    const nm={avg:"AVG",max:"MAX",min:"MIN",sum:"SUM",sd:"STD",median:"MED"}[P.method];
    const rgbOut=im.rgb&&bit===8;
    const n=addImage({name:`${nm}_${im.name}`,w:im.w,h:im.h,bit,rgb:rgbOut,nc:rgbOut?1:nch,nz:1,nt:nOut,planes:out,cal:im.cal,chNames:im.rgb?null:im.chNames},{parent:im.id,src:"derived"});
    if(!im.rgb) im.disp.forEach((d,c)=>{ if(n.disp[c]) n.disp[c].lut=d.lut; });
    if(im.rgb&&!rgbOut){ ["Red","Green","Blue"].forEach((l,c)=>n.disp[c].lut=l); n.composite=true; } // RGB 의 합·SD 는 R·G·B 3 채널(32-bit) 로
    record({op:"zproject",params:clean(P),label:"Z Project "+nm},im); return n; },
  ij:P=>`run("Z Project...", "projection=[${{avg:"Average Intensity",max:"Max Intensity",min:"Min Intensity",sum:"Sum Slices",sd:"Standard Deviation",median:"Median"}[P.method]}] all");`,help:"Hyperstack 은 시간마다 z 를 투영합니다(z 가 하나면 시간을 투영). ImageJ 1.54u 와 같이 평균·최대·최소는 원래 형식으로 (int)(v+0.5) — 정밀도가 필요하면 32-bit 를 켜세요. RGB 의 Sum·SD 는 R·G·B 3 채널 32-bit 영상이 됩니다."});
op("register",{label:"StackReg",group:"stack",kind:"image",params:[
    R("method","Transformation","rigid",{type:"sel",opts:[["translation","Translation"],["rigid","Rigid Body"],["scaled","Scaled Rotation"],["affine","Affine"],["phase","Translation (phase corr.)"]]}),
    R("ref","Reference","previous",{type:"sel",opts:[["previous","Previous slice"],["first","First slice(s)"],["mean","Mean image"]]}),
    R("nFirst","Slices to average",1,{type:"int",min:1,showIf:P=>P.ref==="first"&&P.method!=="phase"}),
    R("channel","Reference channel",1,{type:"int",min:1})],
  run:async(im,P)=>{ if(im.nz*im.nt<2) throw new Error("스택이 아닙니다"); if(im.rgb) throw new Error("RGB 는 먼저 8-bit 또는 Split channels");
    const c0=clamp((+P.channel||1)-1,0,im.nc-1); const both=im.nz>1&&im.nt>1; const along=im.nt>1?"t":"z"; const N=along==="t"?im.nt:im.nz; const all=planesOf(im);
    // 추정용 장: t 축(또는 z 축). z·t 가 모두 있으면 t 마다 z 최대 투영으로 추정해 모든 z 에 같은 변환
    const est=[]; for(let k=0;k<N;k++){ if(both){ const pl=[]; for(let z=0;z<im.nz;z++) pl.push({w:im.w,h:im.h,bit:im.bit,data:all[pIndex(im,c0,z,k)]}); est.push(await call("zProject",[pl,"max"])); }
      else est.push({w:im.w,h:im.h,bit:im.bit,data:all[along==="t"?pIndex(im,c0,0,k):pIndex(im,c0,k,0)]}); }
    // '이전 장' 기준: ImageJ StackReg 와 같이 지금 보이는 장(anchor)에서 앞뒤로 이웃끼리 맞추며 누적
    const anc=P.ref==="previous"?clamp(P.anchor!=null?(+P.anchor-1):(along==="t"?im.t:im.z),0,N-1):0; // 매크로 재생 = 기록한 기준 장
    const anchored=async(a,ref,fn)=>{ if(ref!=="previous"||a===0) return fn(est,0,1); const fw=est.slice(a), bw=est.slice(0,a+1).reverse();
      const F=await fn(fw,0,fw.length/(N+1)), B=await fn(bw,fw.length/(N+1),1); const r=new Array(N); for(let k=0;k<F.length;k++) r[a+k]=F[k]; for(let k=1;k<B.length;k++) r[a-k]=B[k]; return r; };
    const out=all.slice(); let rows;
    if(P.method==="phase"){ if(P.ref==="mean") throw new Error("phase correlation 은 '첫 장' 또는 '이전 장' 기준만");
      setProg("이동량 추정…",0.05); const sh=await anchored(anc,P.ref,(L,f0,f1)=>call("registerShifts",[L,P.ref],f=>setProg("이동량 추정…",f0+(f1-f0)*f)),{dx:0,dy:0});
      for(let k=0;k<N;k++) for(let c=0;c<im.nc;c++) for(let z=0;z<(both?im.nz:1);z++){ const i=along==="t"?pIndex(im,c,both?z:0,k):pIndex(im,c,k,0); out[i]=IC.translate({w:im.w,h:im.h,bit:im.bit,data:out[i]},sh[k].dx,sh[k].dy).data; }
      rows=sh.map(q=>({dx:-q.dx,dy:-q.dy})); } // 표시 = drift(내용이 기준 대비 움직인 양) = −보정 이동
    else { setProg("StackReg 변환 추정…",0.05); const mats=await anchored(anc,P.ref,(L,f0,f1)=>call("stackRegMatrices",[L,P.method,P.ref,Math.max(1,+P.nFirst||1)],f=>setProg("StackReg 변환 추정…",(f0+(f1-f0)*f)*0.7)));
      const jobs=[]; for(let k=0;k<N;k++) for(let c=0;c<im.nc;c++) for(let z=0;z<(both?im.nz:1);z++){ const i=along==="t"?pIndex(im,c,both?z:0,k):pIndex(im,c,k,0); jobs.push([i,k]); }
      let done=0; await Promise.all(jobs.map(([i,k])=>call("turboTransform",[{w:im.w,h:im.h,bit:im.bit,data:out[i]},mats[k],P.method]).then(r=>{ out[i]=IC.roundTo(Float32Array.from(r),im.bit); setProg("변환 적용…",0.7+0.3*(++done)/jobs.length); })));
      // M: 기준 좌표 → 움직인 장 좌표 (x=열, y=행, y 아래로). drift = 영상 중심이 옮겨진 양, 각도 = 화면에서 반시계 + (Orientation 과 같은 규약)
      const cx=(im.w-1)/2, cy=(im.h-1)/2; rows=mats.map(m=>{ const a=m[0][0],b=m[0][1],c=m[1][0],d=m[1][1]; const X=a*cx+b*cy+m[0][2], Y=c*cx+d*cy+m[1][2];
        return {dx:X-cx,dy:Y-cy,angle:-Math.atan2(c-b,a+d)*180/Math.PI,scale:Math.sqrt(Math.abs(a*d-b*c))}; }); }
    hideProg(); delete S.particles[im.id]; clearRois(im);
    setPlanes(im,out,{keepDisplay:true},"Register",{op:"register",params:{...clean(P),...(P.ref==="previous"?{anchor:anc+1}:{})},label:`Register (${P.method==="phase"?"phase correlation":P.method})`});
    const xs=rows.map((_,i)=>i+1); const ser=[{name:"drift x",x:xs,y:rows.map(r=>r.dx)},{name:"drift y",x:xs,y:rows.map(r=>r.dy)}];
    S.plots.push({id:uid(),kind:"xy",title:`Drift (${along}) — ${im.name}`,xl:along==="t"?"frame":"slice",yl:"drift of image centre (px)",series:ser}); S.lastReg={id:im.id,method:P.method,rows};
    if(P.method!=="phase"&&P.method!=="translation") S.plots.push({id:uid(),kind:"xy",title:`Rotation — ${im.name}`,xl:along==="t"?"frame":"slice",yl:"rotation (°, CCW +)",series:[{name:"θ",x:xs,y:rows.map(r=>r.angle)}].concat(P.method!=="rigid"?[{name:"scale×100",x:xs,y:rows.map(r=>100*r.scale)}]:[])});
    const mx=Math.max(...rows.map(r=>Math.hypot(r.dx,r.dy))); note(`정렬 완료 — ${N} 장${P.ref==="previous"?`, 기준 장 ${anc+1}`:""}, 최대 이동 ${mx.toFixed(2)} px${both?" (t 마다 z 최대 투영으로 추정, 모든 z 에 적용)":""} · Plots 탭`,"ok"); S.ui.bottom="plots"; },
  ij:P=>P.method==="phase"?`// translation 정렬(phase correlation): Fiji 'Correct 3D Drift' 또는 StackReg Translation 과 같은 목적`:(P.ref!=="previous"?`// 주의: ImageJ StackReg 는 '이전 장' 기준만 있습니다 (여기서는 ${P.ref==="first"?"첫 장":"평균"} 기준)\n`:"")+(P.ref==="previous"&&P.anchor?`setSlice(${P.anchor});\n`:"")+`run("StackReg ", "transformation=[${{translation:"Translation",rigid:"Rigid Body",scaled:"Scaled Rotation",affine:"Affine"}[P.method]}]");`,
  help:"StackReg/TurboReg (Thévenaz, Ruttimann, Unser, IEEE Trans. Image Process. 7, 27–41, 1998) — B-spline 피라미드에서 세기 차 제곱합 최소화, pystackreg 와 같은 계산. Rigid Body = 이동+회전, Scaled Rotation = +배율, Affine = +전단. Previous slice = 지금 보이는 장에서 앞뒤로 이웃끼리 맞춰 누적(ImageJ StackReg 와 같음). 기준 채널의 변환을 모든 채널에, z·t 가 모두 있으면 t 마다 z 최대 투영으로. 영상 밖은 0."});
/* ---------- Image Calculator ---------- */
op("calc",{label:"Image Calculator",group:"calc",kind:"new",params:[R("b","Image 2","",{type:"img"}),R("fn","Operation","sub",{type:"sel",opts:[["add","Add"],["sub","Subtract"],["mul","Multiply"],["div","Divide"],["and","AND"],["or","OR"],["xor","XOR"],["min","Min"],["max","Max"],["average","Average"],["difference","Difference"],["copy","Copy"],["transparent-zero","Transparent-zero"]]}),R("float","32-bit (float) result",true,{type:"bool"})],
  run:async(im,P)=>{ const b=byId(P.b)||S.images.find(x=>x.name===P.bName); if(!b) throw new Error("두 번째 영상을 고르세요"); if(b.w!==im.w||b.h!==im.h) throw new Error("두 영상의 크기가 같아야 합니다");
    const r=await call("calc",[curPlaneOf(im),curPlaneOf(b),P.fn,!!P.float]); const n=addImage({name:`Result of ${im.name}`,w:im.w,h:im.h,bit:r.bit,planes:[r.data],cal:im.cal},{parent:im.id,src:"derived"}); record({op:"calc",params:{...clean(P),bName:b.name},label:`Calculator ${P.fn}`},im); return n; },
  ij:P=>`imageCalculator("${{add:"Add",sub:"Subtract",mul:"Multiply",div:"Divide",and:"AND",or:"OR",xor:"XOR",min:"Min",max:"Max",average:"Average",difference:"Difference",copy:"Copy","transparent-zero":"Transparent-zero"}[P.fn]} create${P.float?" 32-bit":""}", getTitle(), "${ijq(P.bName||"image2")}");`});

/* ---------- Threshold · Binary (패널이 따로 있음) ---------- */
op("threshold",{label:"Convert to Mask",group:"thr",kind:"plane",hidden:true,params:[],
  run:async(im,P)=>{ // P: {method, dark, manual, lower, upper, stack, newImage, local:{method,radius,p1,p2}}
    const planes=planesOf(im); const idx=[]; for(let k=0;k<planes.length;k++) if(k%im.nc===(im.rgb?0:im.c)||im.nc===1) idx.push(k); // 스택은 언제나 모든 장 (ImageJ 처럼 스택 유지)
    let fixed=null; if(isStack(im)&&!P.stack&&!P.manual&&!P.local){ const t=await call("autoThreshold",[{w:im.w,h:im.h,bit:im.bit,data:planes[pIndex(im,im.rgb?0:im.c,im.z,im.t)]},P.method,!!P.dark]); fixed=[t.lower,t.upper]; } // '모든 장' 끔 = 지금 장에서 구한 문턱을 모든 장에
    if(im.rgb) throw new Error("RGB 는 먼저 Type ▸ RGB → 8-bit gray");
    const out=[]; let msg="";
    for(const k of idx){ const p={w:im.w,h:im.h,bit:im.bit,data:planes[k]}; let m;
      if(P.local){ m=await call("localThreshold",[p,P.local.method,{radius:+P.local.radius,p1:P.local.p1===""?null:+P.local.p1,p2:P.local.p2===""?null:+P.local.p2,dark:P.dark}]); msg=`Local ${P.local.method} r=${P.local.radius}`; }
      else { let lo=+P.lower, hi=+P.upper; if(fixed){ [lo,hi]=fixed; } else if(!P.manual){ const t=await call("autoThreshold",[p,P.method,!!P.dark]); lo=t.lower; hi=t.upper; } m=IC.applyThreshold(p,lo,hi); msg=P.manual?`Threshold ${g4(lo)}–${g4(hi)}`:`${P.method}${P.dark?" dark":""} → [${g4(lo)}, ${g4(hi)}]`; }
      out.push(m.data); }
    const step={op:"threshold",params:clean(P),label:"Convert to Mask ("+msg+")"};
    if(P.newImage){ const n=addImage({name:"Mask of "+im.name,w:im.w,h:im.h,bit:8,nc:1,nz:idx.length>1?im.nz:1,nt:idx.length>1?im.nt:1,planes:out,cal:im.cal},{parent:im.id,src:"derived",kind:"mask"}); record(step,im); note(msg,"ok"); return n; }
    setPlanes(im,out,{bit:8,nc:1,c:0,composite:false,nz:idx.length>1?im.nz:1,nt:idx.length>1?im.nt:1},"Convert to Mask",step); im.disp=[{min:0,max:255,lut:"Grays"}]; im.kind="mask"; note(msg,"ok"); },
  ij:P=>P.local?`run("Auto Local Threshold", "method=${P.local.method} radius=${P.local.radius} parameter_1=${P.local.p1||0} parameter_2=${P.local.p2||0}${P.dark?" white":""}");`:(P.manual?`setThreshold(${P.lower}, ${P.upper}); setOption("BlackBackground", true); run("Convert to Mask", "method=Default background=${P.dark?"Dark":"Light"} black");`:`setAutoThreshold("${P.method}${P.dark?" dark":""}"); setOption("BlackBackground", true); run("Convert to Mask", "method=${P.method} background=${P.dark?"Dark":"Light"}${P.stack?" calculate":""} black");`)}); // 스택: calculate = 장마다 문턱, 아니면 지금 장의 문턱을 모든 장에
/* Color Threshold (ImageJ ColorThresholder): RGB 영상 → 8-bit mask. P: {space, min[3], max[3], pass[3], stack, newImage} */
const ctFilt=P=>({space:P.space,min:P.min.map(Number),max:P.max.map(Number),pass:P.pass.map(Boolean)});
op("colorThreshold",{label:"Color Threshold",group:"thr",kind:"plane",hidden:true,params:[],
  run:async(im,P)=>{ if(!im.rgb) throw new Error("Color Threshold 는 RGB 영상에서 — 회색조는 Threshold");
    if(!P.newImage&&isStack(im)&&!P.stack) throw new Error("스택을 제자리에서 바꾸려면 '모든 장' 을 켜세요 (아니면 '새 영상으로')");
    if(!IC.COLOR_SPACES.includes(P.space)) throw new Error("색 공간: "+P.space); const o=ctFilt(P); const planes=planesOf(im); const out=[];
    const zt=[]; if(P.stack){ for(let t=0;t<im.nt;t++) for(let z=0;z<im.nz;z++) zt.push([z,t]); } else zt.push([im.z,im.t]);
    let k=0; for(const [z,t] of zt){ setProg(`Color Threshold ${++k}/${zt.length}`,k/zt.length);
      if(z===im.z&&t===im.t&&o.space===S.ui.cthr.space&&ctReady(im)){ out.push(IC.colorMask(CT.ch,o)); continue; } /* 패널이 이미 변환한 지금 장은 다시 변환하지 않음 */
      const u=IC.rgbU8(planes[pIndex(im,0,z,t)],planes[pIndex(im,1,z,t)],planes[pIndex(im,2,z,t)]); const r=await call("colorThreshold",[u[0],u[1],u[2],im.w,im.h,o]); out.push(r.data); } hideProg();
    const nm=IC.COLOR_CHANNEL_NAMES[o.space]; const msg=`${o.space}: `+[0,1,2].map(i=>`${nm[i]} ${o.pass[i]?"":"stop "}${o.min[i]}–${o.max[i]}`).join(", ");
    const step={op:"colorThreshold",params:clean(P),label:"Color Threshold ("+msg+")"}; const nz=P.stack?im.nz:1, nt=P.stack?im.nt:1;
    if(P.newImage){ const n=addImage({name:"Mask of "+im.name,w:im.w,h:im.h,bit:8,nc:1,nz,nt,planes:out,cal:im.cal},{parent:im.id,src:"derived",kind:"mask"}); record(step,im); note(msg,"ok"); return n; }
    setPlanes(im,out,{rgb:false,bit:8,nc:1,c:0,composite:false,nz,nt},"Color Threshold",step); im.disp=[{min:0,max:255,lut:"Grays"}]; im.kind="mask"; note(msg,"ok"); },
  ij:P=>{ const o=ctFilt(P); const L=[];
    // 한 장(지금 선택된 RGB 영상 하나)에 대한 ImageJ ColorThresholder 매크로 단계 → 결과 mask 창이 선택된 채로 끝남
    const pipe=[...(o.space==="HSB"?['run("HSB Stack");','run("Convert Stack to Images");','selectWindow("Hue"); rename("0");','selectWindow("Saturation"); rename("1");','selectWindow("Brightness"); rename("2");']
        :[...(o.space==="Lab"?['call("ij.plugin.frame.ColorThresholder.RGBtoLab");']:o.space==="YUV"?['call("ij.plugin.frame.ColorThresholder.RGBtoYUV");']:[]),'run("RGB Stack");','run("Convert Stack to Images");','selectWindow("Red"); rename("0");','selectWindow("Green"); rename("1");','selectWindow("Blue"); rename("2");']),
      `ct_min=newArray(${o.min.join(",")}); ct_max=newArray(${o.max.join(",")}); ct_f=newArray("${o.pass.map(p=>p?"pass":"stop").join('","')}");`,
      'for (ct_i=0; ct_i<3; ct_i++) { selectWindow(""+ct_i); setThreshold(ct_min[ct_i], ct_max[ct_i]); run("Convert to Mask"); if (ct_f[ct_i]=="stop") run("Invert"); }',
      'imageCalculator("AND create", "0","1");','imageCalculator("AND create", "Result of 0","2");','for (ct_i=0; ct_i<3; ct_i++) { selectWindow(""+ct_i); close(); }','selectWindow("Result of 0"); close();','selectWindow("Result of Result of 0");'];
    L.push("// ImageJ ColorThresholder 가 만드는 매크로와 같은 단계 (변수 ct_ — batch 고리의 i 를 덮지 않음)","ct_a=getTitle();");
    if(!P.stack){ if(P.newImage) L.push('run("Duplicate...", "title=[ct_src]");'); L.push(...pipe, `rename(${P.newImage?'"Mask of "+ct_a':"ct_a"});`); return L.join("\n"); }
    // 스택: 장마다 한 장씩 복제해 같은 단계 → Images to Stack (hyperstack 이면 z·t 되살림)
    L.push("ct_id=getImageID(); getDimensions(ct_w, ct_h, ct_c, ct_z, ct_t); ct_n=nSlices;","for (ct_s=1; ct_s<=ct_n; ct_s++) {","  selectImage(ct_id); setSlice(ct_s); run(\"Duplicate...\", \"title=[ct_src]\");",
      ...pipe.map(x=>"  "+x),'  rename("ct_m"+IJ.pad(ct_s,5));',"}",'run("Images to Stack", "name=[ct_stack] title=ct_m use");',
      'if (ct_z>1 && ct_t>1) run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+ct_z+" frames="+ct_t+" display=Grayscale");');
    if(P.newImage) L.push('rename("Mask of "+ct_a);'); else L.push("ct_m=getImageID(); selectImage(ct_id); close(); selectImage(ct_m); rename(ct_a);");
    return L.join("\n"); },
  help:"ImageJ Image ▸ Adjust ▸ Color Threshold 와 같은 색 변환(HSB·RGB·Lab·YUV)과 pass/stop 규칙."});
const BIN=[["erode","Erode"],["dilate","Dilate"],["open","Open"],["close","Close-"],["fill","Fill Holes"],["outline","Outline"],["skeleton","Skeletonize"],["watershed","Watershed"],["edm","Distance Map"],["uep","Ultimate Points"],["voronoiX","—"]].slice(0,10);
op("binary",{label:"Binary",group:"bin",hidden:true,params:[R("fn","Operation","erode",{type:"sel",opts:BIN}),R("iterations","Iterations",1,{type:"int",min:1}),R("count","Count",1,{type:"int",min:1,max:8}),R("edm32","Distance Map as new 32-bit image",false,{type:"bool"}),stackParam],
  run:async(im,P)=>{ const f=P.fn; const all=!!P.stack||!isStack(im); const ix=all?planesOf(im).map((_,k)=>k):[pIndex(im,0,im.z,im.t)];
    if(f==="edm"&&P.edm32){ const r=await call("distanceMap",[curPlaneOf(im),"32-bit"]); const n=addImage({name:"EDM of "+im.name,w:im.w,h:im.h,bit:32,planes:[r.data],cal:im.cal},{parent:im.id,src:"derived"}); record({op:"binary",params:clean(P),label:"Distance Map (32-bit)"},im); return n; }
    const P2=planesOf(im).slice(); const lab=(BIN.find(b=>b[0]===f)||[0,f])[1];
    for(const k of ix){ const p={w:im.w,h:im.h,bit:8,data:P2[k]}; let r;
      if(f==="erode"||f==="dilate"||f==="open"||f==="close") r=await call("binaryMorph",[p,f,{iterations:+P.iterations,count:+P.count}]);
      else if(f==="fill") r=await call("fillHoles",[p]); else if(f==="outline") r=await call("outline",[p]); else if(f==="skeleton") r=await call("skeletonize",[p]);
      else if(f==="watershed") r=await call("watershed",[p]); else if(f==="edm") r=await call("distanceMap",[p,"8-bit"]); else if(f==="uep") r=await call("ultimatePoints",[p]);
      P2[k]=r.data; }
    setPlanes(im,P2,{keepDisplay:f!=="edm"&&f!=="uep"},lab,{op:"binary",params:clean(P),label:lab+((f==="erode"||f==="dilate"||f==="open"||f==="close")?` (×${P.iterations}, count ${P.count})`:"")});
    if(f==="edm"||f==="uep"){ im.kind=""; } },
  ij:P=>{ const nm=(BIN.find(b=>b[0]===P.fn)||[0,P.fn])[1]; if(["erode","dilate","open","close"].includes(P.fn)) return `run("Options...", "iterations=${P.iterations} count=${P.count} black do=Nothing"); run("${nm==="Close-"?"Close-":nm}"${P.stack?', "stack"':', " "'});`; if(P.fn==="edm"&&P.edm32) return `run("Duplicate...", "title=edm_src"); run("Options...", "iterations=1 count=1 black edm=32-bit"); run("Distance Map"); run("Options...", "iterations=1 count=1 black edm=Overwrite");`; /* 앱: 지금 장에서 32-bit 새 영상 */ return `run("${nm}"${P.stack?', "stack"':', " "'});`; }});

/* 현재 장 plane (RGB 면 밝기) */
function curPlaneOf(im){ if(im.rgb){ return IC.rgbToGray(planeAt(im,0),planeAt(im,1),planeAt(im,2),false); } return planeAt(im); }
