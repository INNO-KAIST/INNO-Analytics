// Color Threshold — ImageJ 1.54u 를 직접 실행해 대조
//  1) 모든 24-bit 색(2²⁴ = 4096×4096 RGB) 의 HSB·Lab·YUV 변환 = ImageJ (HSB Stack / ColorThresholder.RGBtoLab·RGBtoYUV)
//  2) 문턱(pass/stop·세 채널 AND) = ImageJ ColorThresholder 가 만드는 매크로(setThreshold → Convert to Mask → Invert → AND)
//  3) Auto(autoSetThreshold) 의 채널 선택·문턱 = ImageJ AutoThresholder
// 사용: node test/test_color.js [ij.jar]
const fs=require("fs"), cp=require("child_process"), path=require("path"), os=require("os");
const IC=require("./load_core.js"); const IO=require("../src/img_io.js");
const JAR=process.argv[2]||"/tmp/ij.jar";
let pass=0, fail=0; const fails=[]; const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(c||!d?"":" — "+d)); };
if(!fs.existsSync(JAR)){ console.log("ImageJ jar 없음 — 건너뜀"); process.exit(0); }
const TMP=fs.mkdtempSync(path.join(os.tmpdir(),"ijcolor-"));
const ij=(macro,mem)=>{ const f=path.join(TMP,"m"+Math.random().toString(36).slice(2)+".ijm"); fs.writeFileSync(f,macro); return cp.execSync(`JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx${mem||"3g"} -cp ${JAR} ij.ImageJ -batch ${f} 2>&1`,{encoding:"utf8",maxBuffer:1<<26}); };
function rgbTiff(file,w,h,R,G,B){ fs.writeFileSync(file,IO.writeTIFF({w,h,bit:8,rgb:true,nc:3,nz:1,nt:1,planes:[R,G,B],cal:null})); }

// ---------- 1) 모든 색 ----------
{ const W=4096,H=4096,N=W*H; const R=new Float32Array(N),G=new Float32Array(N),B=new Float32Array(N);
  for(let i=0;i<N;i++){ R[i]=(i>>16)&255; G[i]=(i>>8)&255; B[i]=i&255; }
  const tif=path.join(TMP,"all.tif"); rgbTiff(tif,W,H,R,G,B);
  const out={HSB:path.join(TMP,"hsb.raw"),Lab:path.join(TMP,"lab.raw"),YUV:path.join(TMP,"yuv.raw")};
  ij(`setBatchMode(true);
open("${tif}"); run("HSB Stack"); saveAs("Raw Data", "${out.HSB}"); close();
open("${tif}"); call("ij.plugin.frame.ColorThresholder.RGBtoLab"); run("RGB Stack"); saveAs("Raw Data", "${out.Lab}"); close();
open("${tif}"); call("ij.plugin.frame.ColorThresholder.RGBtoYUV"); run("RGB Stack"); saveAs("Raw Data", "${out.YUV}"); close();
`,"6g");
  for(const sp of ["HSB","Lab","YUV","RGB"]){ const ch=IC.colorChannels(R,G,B,sp);
    if(sp==="RGB"){ ok(ch[0][0x123456]===0x12&&ch[1][0x123456]===0x34&&ch[2][0x123456]===0x56,"RGB 채널 = R,G,B"); continue; }
    const ref=fs.readFileSync(out[sp]); let bad=[0,0,0], ex="";
    for(let k=0;k<3;k++){ const a=ch[k]; for(let i=0;i<N;i++){ if(a[i]!==ref[k*N+i]){ if(!bad[k]&&!ex) ex=`채널 ${k} 색 #${i.toString(16).padStart(6,"0")}: 앱 ${a[i]} ImageJ ${ref[k*N+i]}`; bad[k]++; } } }
    ok(!bad[0]&&!bad[1]&&!bad[2],`${sp}: 2²⁴ 색 모두 ImageJ 와 같음`,`다른 수 ${bad.join("/")} ${ex}`); }
}
// ---------- 2) 문턱 = ImageJ 생성 매크로 ----------
const W=97,H=73,N=W*H; const rnd=(()=>{ let s=2468; return ()=>((s=(s*1103515245+12345)&0x7fffffff)/0x7fffffff); })();
const R=new Float32Array(N),G=new Float32Array(N),B=new Float32Array(N);
for(let y=0;y<H;y++) for(let x=0;x<W;x++){ const i=y*W+x; R[i]=Math.round(255*Math.abs(Math.sin(x/9)))^((rnd()*16)|0); G[i]=Math.round(255*(y/H))^((rnd()*32)|0); B[i]=Math.round(128+120*Math.cos((x+y)/11))&255; }
const tif=path.join(TMP,"img.tif"); rgbTiff(tif,W,H,R,G,B);
function ijMask(space,min,max,pass){
  const conv=space==="HSB"?`run("HSB Stack"); run("Convert Stack to Images"); selectWindow("Hue"); rename("0"); selectWindow("Saturation"); rename("1"); selectWindow("Brightness"); rename("2");`
    :`${space==="Lab"?'call("ij.plugin.frame.ColorThresholder.RGBtoLab");':space==="YUV"?'call("ij.plugin.frame.ColorThresholder.RGBtoYUV");':""} run("RGB Stack"); run("Convert Stack to Images"); selectWindow("Red"); rename("0"); selectWindow("Green"); rename("1"); selectWindow("Blue"); rename("2");`;
  const o=path.join(TMP,"mask_"+Math.random().toString(36).slice(2)+".raw");
  ij(`setBatchMode(true); setOption("BlackBackground", true); open("${tif}"); ${conv}
min=newArray(${min.join(",")}); max=newArray(${max.join(",")}); filter=newArray("${pass.map(p=>p?"pass":"stop").join('","')}");
for (i=0;i<3;i++){ selectWindow(""+i); setThreshold(min[i], max[i]); run("Convert to Mask"); if (filter[i]=="stop") run("Invert"); }
imageCalculator("AND create", "0","1"); imageCalculator("AND create", "Result of 0","2"); saveAs("Raw Data", "${o}");`);
  return fs.readFileSync(o); }
const cases=[["HSB",[20,40,30],[120,255,230],[true,true,true]],["HSB",[200,0,0],[40,255,255],[false,true,true]],["HSB",[30,60,50],[90,200,180],[false,false,true]],
  ["RGB",[50,0,100],[200,128,255],[true,false,true]],["Lab",[40,100,90],[200,160,170],[true,true,false]],["YUV",[60,110,100],[190,150,160],[true,true,true]],["RGB",[0,0,0],[255,255,255],[true,true,true]]];
for(const [sp,mn,mx,ps] of cases){ const a=IC.colorThreshold(R,G,B,W,H,{space:sp,min:mn,max:mx,pass:ps}).data; const r=ijMask(sp,mn,mx,ps); let bad=0; for(let i=0;i<N;i++) if(a[i]!==r[i]) bad++;
  const nsel=a.reduce((s,v)=>s+(v>0),0); ok(!bad,`${sp} [${mn}]–[${mx}] ${ps.map(p=>p?"pass":"stop").join("/")} (${nsel} px) = ImageJ`,`다른 화소 ${bad}`); }
// ---------- 3) Auto ----------
{ const {hist}=IC.colorHist(R,G,B,"HSB"); const a=IC.colorAuto(hist,"HSB","Otsu",true); const keys=Object.keys(a);
  const L=ij(`setBatchMode(true); open("${tif}"); run("HSB Stack"); setSlice(3); setAutoThreshold("Otsu dark"); getThreshold(lo,hi); print("T "+lo+" "+hi);`).split("\n").find(l=>l.startsWith("T "));
  const [lo,hi]=L.slice(2).split(" ").map(Number); ok(keys.join()==="2"&&a[2][0]===lo&&a[2][1]===hi,"Auto HSB: Brightness 만, Otsu dark = ImageJ",`앱 ${JSON.stringify(a)} ImageJ ${lo}..${hi}`);
  const r=IC.colorAuto(IC.colorHist(R,G,B,"RGB").hist,"RGB","Default",false); ok(Object.keys(r).join()==="0,1,2"&&r[0][0]===0,"Auto RGB: 세 채널 모두, 밝은 배경 → [0,t]",JSON.stringify(r)); }
// ---------- 4) Sample (Hue 빈 구간) ----------
{ const ch=[Uint8Array.from([250,252,3,5,128]),Uint8Array.from([10,20,30,40,50]),Uint8Array.from([100,90,80,70,60])];
  const s=IC.colorSample(ch,[0,1,2,3],"HSB"); ok(!s.pass0&&s.min[0]===6&&s.max[0]===249&&s.min[1]===10&&s.max[1]===40,"Sample: 빨강(0 근처)을 감싸는 Hue → stop 6–249",JSON.stringify(s));
  const t=IC.colorSample(ch,[4],"HSB"); ok(t.pass0&&t.min[0]===128&&t.max[0]===128,"Sample: 한 색 → pass",JSON.stringify(t)); }
// ---------- 5) RGB hyperstack TIFF 읽기·쓰기 (z·t 유지, ImageJ 로 확인) ----------
(async()=>{ const u8=new Uint8Array(fs.readFileSync(path.join(__dirname,"io","rgb_hyper.tif"))); const r=await IO.readFile("rgb_hyper.tif",u8,{}); const im=r.images[0];
  ok(im.rgb&&im.nz===2&&im.nt===2,"RGB hyperstack TIFF 읽기 → z 2 · t 2 (ImageJ 와 같음)",`${im.nz}×${im.nt}`);
  const f=path.join(TMP,"rh.tif"); fs.writeFileSync(f,IO.writeTIFF(im)); const L=ij(`open("${f}"); getDimensions(w,h,c,z,t); print("D "+z+" "+t+" "+bitDepth());`);
  ok(/D 2 2 24/.test(L),"RGB hyperstack TIFF 쓰기 → ImageJ 가 z 2 · t 2 RGB 로 읽음",L.split("\n").filter(x=>x.startsWith("D ")).join(""));
  console.log(`\n${pass} 통과, ${fail} 실패`); try{ fs.rmSync(TMP,{recursive:true}); }catch(e){} if(fail){ console.log(fails.join("\n")); process.exit(1); } })();
