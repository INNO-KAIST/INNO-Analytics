// 경계 사례(NaN·±∞·일정한 영상·퇴화) — ImageJ 1.54u 를 직접 실행해 대조 (2026-10-16 전수 검수 회귀)
// 사용: node test/test_ij_edge.js [ij.jar]   (ImageJ 소스를 javac 로 컴파일한 jar; 없으면 건너뜀)
const fs=require("fs"), cp=require("child_process"), path=require("path"), os=require("os");
const IC=require("./load_core.js");
const JAR=process.argv[2]||"/tmp/ij.jar";
if(!fs.existsSync(JAR)){ console.log("ImageJ jar 없음 — 건너뜀 ("+JAR+")"); process.exit(0); }
const TMP=fs.mkdtempSync(path.join(os.tmpdir(),"ijedge-")); let K=0, pass=0, fail=0; const fails=[];
const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d?" — "+d:"")); } console.log((c?"  ✓ ":"  ✗ ")+m+(c||!d?"":" — "+d)); };
function ij(macro){ const f=path.join(TMP,"m"+(K++)+".ijm"); fs.writeFileSync(f,macro); return cp.execSync(`JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx1g -cp ${JAR} ij.ImageJ -batch ${f} 2>/dev/null`,{encoding:"utf8",maxBuffer:1<<28}).split("\n"); }
const jv=v=>v!==v?"NaN":v===Infinity?"1/0":v===-Infinity?"-1/0":String(v);
function mk(name,p){ const t=p.bit===8?"8-bit":p.bit===16?"16-bit":"32-bit"; let s=`newImage("${name}","${t} black",${p.w},${p.h},1);\n`;
  for(let y=0;y<p.h;y++) for(let x=0;x<p.w;x++){ const v=p.data[y*p.w+x]; if(v!==0) s+=`setPixel(${x},${y},${jv(v)});\n`; } return s+"resetMinAndMax();\n"; }
const dump=`function dump(){ s=""; for(y=0;y<getHeight();y++){ for(x=0;x<getWidth();x++) s=s+d2s(getPixel(x,y),6)+" "; } print("PIX "+s); }\n`;
const pix=lines=>lines.filter(l=>l.startsWith("PIX ")).map(l=>l.trim().split(/\s+/).slice(1).map(Number));
const P=(w,h,bit,arr)=>({w,h,bit,data:Float32Array.from(arr)});
function same(a,b,tol){ if(a.length!==b.length) return `길이 ${a.length}≠${b.length}`; let n=0,ex=""; for(let i=0;i<a.length;i++){ const x=a[i],y=b[i]; if(x!==x&&y!==y) continue; if(x===y) continue; if(!(Math.abs(x-y)<=(tol||1e-5)*Math.max(1,Math.abs(y)))){ n++; if(!ex) ex=` 예: #${i} 앱 ${x} ImageJ ${y}`; } } return n?`${n} 곳 다름${ex}`:""; }
const rnd=(()=>{ let s=12345; return ()=>((s=(s*1103515245+12345)&0x7fffffff)/0x7fffffff); })();
// 1) rank filter 에 NaN
{ const W=30,H=20; const d=[]; for(let i=0;i<W*H;i++) d.push(Math.round(rnd()*1000)/10); d[7*W+11]=NaN; d[3*W+2]=NaN; const p=P(W,H,32,d);
  for(const [f,ijn,o] of [["max","Maximum...",{radius:2}],["min","Minimum...",{radius:2}],["mean","Mean...",{radius:2}],["variance","Variance...",{radius:1.5}],["median","Median...",{radius:2}]]){
    const r=pix(ij(dump+mk("a",p)+`run("${ijn}", "radius=${o.radius}");\ndump();\n`))[0]; const a=Array.from(IC.rank(p,f,o).data);
    if(f==="mean"||f==="variance"){ // ImageJ 는 이동 합이 NaN 창을 벗어난 다음 화소 하나도 NaN 으로 남김(줄마다) — 여기서는 커널 안에 NaN 이 있는 화소만 NaN (SPEC §3)
      const {kR,dx}=IC.lineRadii(o.radius); let foot=0, appNaN=0; for(let y=0;y<H;y++) for(let x=0;x<W;x++){ let hit=false; for(let k=-kR;k<=kR&&!hit;k++){ const yy=Math.min(H-1,Math.max(0,y+k)); for(let xx=x-dx[k+kR];xx<=x+dx[k+kR];xx++){ const X=Math.min(W-1,Math.max(0,xx)); if(d[yy*W+X]!==d[yy*W+X]){ hit=true; break; } } } if(hit) foot++; if(a[y*W+x]!==a[y*W+x]) appNaN++; }
      const e=same(a.map((v,i)=>r[i]!==r[i]?NaN:v),r,1e-5); ok(!e&&foot===appNaN,`${f} (32-bit, NaN 2 개): 유한값 = ImageJ, NaN = 커널이 닿는 화소만`,`${e} NaN 앱 ${appNaN} / 커널 ${foot} / ImageJ ${r.filter(v=>v!==v).length}`); continue; }
    const e=same(a,r,1e-5); ok(!e,`${f} (32-bit, NaN 2 개) = ImageJ`,e); } }
// 2) Find Maxima 32-bit + NaN 덩어리
{ const W=12,H=9, d=[]; for(let y=0;y<H;y++) for(let x=0;x<W;x++) d.push(10*Math.exp(-((x-9)**2+(y-4)**2)/3)+8*Math.exp(-((x-2)**2+(y-2)**2)/3)+6*Math.exp(-((x-7)**2+(y-7)**2)/2));
  for(let y=1;y<6;y++) for(let x=4;x<7;x++) d[y*W+x]=NaN; const p=P(W,H,32,d);
  const L=ij(mk("a",p)+`run("Find Maxima...", "prominence=3 output=[Point Selection]");\ngetSelectionCoordinates(xs,ys); s=""; for(i=0;i<xs.length;i++) s=s+xs[i]+","+ys[i]+" "; print("MAX "+s);\n`).find(l=>l.startsWith("MAX "))||"MAX ";
  const ijp=L.slice(4).trim().split(/\s+/).filter(Boolean).sort().join(" "); const r=IC.findMaxima(p,3,{output:"points"}); const pts=(r.points||r.xy||r).map?((r.points||r.xy||r)).map(q=>Array.isArray(q)?q.join(","):`${q.x},${q.y}`).sort().join(" "):"";
  ok(pts===ijp,"Find Maxima (32-bit, NaN 덩어리) = ImageJ",`앱 [${pts}] ImageJ [${ijp}]`); }
// 3) Image Calculator 32-bit Divide (정수 결과 체크 꺼도 실수)
{ const a=P(3,1,32,[7.5,0,-5]), b=P(3,1,32,[2,0,0]); const r=IC.calc(a,b,"div",false).data; const e=same(Array.from(r),[3.75,NaN,-Infinity]); ok(!e,"Image Calculator Divide 32-bit = 실수 나눗셈 (0/0=NaN, −5/0=−∞)",e); }
// 4) Z Project 32-bit NaN
{ const s=[P(2,1,32,[NaN,1]),P(2,1,32,[2,5]),P(2,1,32,[3,NaN])];
  for(const [m,e] of [["avg",[2.5,3]],["max",[3,5]],["min",[2,1]]]){ const r=Array.from(IC.zProject(s,m).data); const er=same(r,e); ok(!er,`Z Project ${m} 은 NaN 을 뺌 (ImageJ)`,er); } }
// 5) Apply LUT 8-bit 반올림
{ const p=P(256,1,8,[...Array(256).keys()]); const r=pix(ij(dump+mk("a",p)+`setMinAndMax(50,100);\nrun("Apply LUT");\ndump();\n`))[0]; const e=same(Array.from(IC.applyLut(p,50,100).data),r,0); ok(!e,"Apply LUT (8-bit, 50–100) = ImageJ",e); }
// 6) 임의 각도 회전 가장자리
{ const W=40,H=30, p=P(W,H,8,new Array(W*H).fill(100)); const r=pix(ij(dump+mk("a",p)+`run("Rotate... ", "angle=1 grid=1 interpolation=Bilinear");\ndump();\n`))[0]; const a=Array.from(IC.rotate(p,1,{}).data);
  const e=same(a,r,0); ok(!e,"Rotate 1° (일정한 영상) 가장자리 = ImageJ",e); }
// 7) 정수·축 평행 다각형 ROI 둘레 = traced
{ const W=20,H=20, p=P(W,H,8,new Array(W*H).fill(5)); const L=ij(mk("a",p)+`makePolygon(2,2,12,2,12,12,2,12);\nrun("Set Measurements...", "area perimeter redirect=None decimal=6");\nrun("Measure");\nprint("PER "+getResult("Perim.",0));\n`).find(l=>l.startsWith("PER "));
  const mask=IC.polygonMask([2,12,12,2],[2,2,12,12],W,H); const m=IC.measureRoi(p,{mask:mask.mask,rect:mask.rect,poly:{xs:[2,12,12,2],ys:[2,2,12,12]},kind:"poly"},{pw:1,ph:1}); ok(Math.abs(m["Perim."]-(+L.slice(4)))<1e-4,"정수 축 평행 다각형 둘레 = ImageJ (traced)",`앱 ${m["Perim."]} ImageJ ${L}`); }
// 8) ±∞ 측정
{ const p=P(4,3,32,[1,2,3,NaN,Infinity,-Infinity,4,5,1,2,3,0]); const L=ij(mk("a",p)+`run("Set Measurements...", "area mean min redirect=None decimal=6");\nrun("Measure");\nprint("M "+d2s(getResult("Area",0),9)+" "+d2s(getResult("Mean",0),9)+" "+d2s(getResult("Min",0),9)+" "+d2s(getResult("Max",0),9));\n`).find(l=>l.startsWith("M "));
  const s=IC.stats(p); const a=[s.area,s.mean,s.min,s.max], b=L.slice(2).trim().split(/\s+/).map(Number); const e=same(a,b,1e-5); ok(!e,"측정: ±∞ 는 빼고 (ImageJ FloatStatistics)",e+` 앱 ${a} ImageJ ${b}`); }
// 9) 일정한 영상: Enhance Contrast normalize 는 그대로, 8-bit 변환은 0
{ const p=P(5,4,16,new Array(20).fill(1234)); const r=IC.enhanceContrast(p,{saturated:0.35,normalize:true}); ok(r.unchanged&&r.plane.data.every(v=>v===1234),"Enhance Contrast (일정한 영상) → 바꾸지 않음 (ImageJ)");
  const q=P(5,4,32,new Array(20).fill(3.5)); const c=IC.to8bitScaled(q,3.5,3.5); ok(c.every(v=>v===0),"32→8-bit (min=max) → 0, NaN 없음 (ImageJ)"); }
// 10) 빈 histogram 문턱: 8-bit 범위 안
{ const p=P(6,4,8,new Array(24).fill(0)); let okAll=true, ex=""; for(const m of IC.THRESHOLD_METHODS){ for(const dark of [true,false]){ const t=IC.autoThreshold(p,m,dark); if(!(t.lower>=0&&t.lower<=255&&t.upper>=0&&t.upper<=255)){ okAll=false; ex=`${m} ${dark}: ${t.lower}..${t.upper}`; } } } ok(okAll,"모두 0 인 8-bit: 문턱이 0..255 안",ex);
  const q=P(6,4,16,new Array(24).fill(700)); const t=IC.autoThreshold(q,"Otsu",true); ok(t.none===true,"일정한 16-bit: 문턱 없음 (ImageJ −1,−1)",JSON.stringify(t.lower)); }
// 11) Subtract Background 얇은 영상 → 분명한 오류
{ let err=null; try{ IC.subtractBackground(P(3,5,8,new Array(15).fill(10)),{radius:50}); }catch(e){ err=e; } ok(err&&/얇/.test(err.message),"Subtract Background (3×5, r=50) → 분명한 오류",err&&err.message); }
// 13) 큰 창 Median·Remove Outliers (16/32-bit, n ≥ 25 → 순위 이동 히스토그램) — NaN·±∞·같은 값 많음, 가장자리
{ const W=37,H=29; const mk16=(f)=>{ const d=[]; for(let i=0;i<W*H;i++) d.push(f(i)); return P(W,H,16,d); };
  const g16=[["16-bit 넓은 범위",mk16(()=>Math.floor(rnd()*65536))],["16-bit 같은 값 많음",mk16(()=>30000+Math.floor(rnd()*4))]];
  const d32=[]; for(let i=0;i<W*H;i++) d32.push(Math.round(rnd()*2000-1000)/10); for(let y=8;y<12;y++) for(let x=14;x<19;x++) d32[y*W+x]=NaN; d32[3*W+30]=Infinity; d32[20*W+5]=-Infinity; d32[0]=NaN; d32[W*H-1]=Infinity;
  const d32b=d32.map(v=>v===v?Math.round(v/50)*5:v); // 같은 값 많은 float + NaN
  const all=g16.concat([["32-bit NaN·±∞",P(W,H,32,d32)],["32-bit 같은 값 많음 + NaN",P(W,H,32,d32b)]]);
  for(const [nm,p] of all) for(const r of [2.5,5,8.5]){
    const L=pix(ij(dump+mk("a",p)+`run("Median...", "radius=${r}");\ndump();\n`))[0]; const a=Array.from(IC.rank(p,"median",{radius:r}).data);
    const e=same(a,L,0); ok(!e,`Median r=${r} (${nm}) = ImageJ`,e); }
  for(const [nm,p] of all) for(const [r,thr] of [[6,p.bit===16?(nm.includes("같은")?1:5000):20],[1,p.bit===16?0.5:12.3]]) for(const wh of ["Bright","Dark"]){
    const L=pix(ij(dump+mk("a",p)+`run("Remove Outliers...", "radius=${r} threshold=${thr} which=${wh}");\ndump();\n`))[0]; const a=Array.from(IC.rank(p,"outliers",{radius:r,threshold:thr,which:wh.toLowerCase()}).data);
    const e=same(a,L,0); ok(!e,`Remove Outliers r=${r} thr=${thr} ${wh} (${nm}) = ImageJ`,e); } }
console.log(`\n${pass} 통과, ${fail} 실패`); try{ fs.rmSync(TMP,{recursive:true}); }catch(e){} if(fail){ console.log(fails.join("\n")); process.exit(1); }
