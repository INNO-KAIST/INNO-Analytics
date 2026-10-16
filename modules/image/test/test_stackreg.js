/* StackReg/TurboReg 검증 — pystackreg 0.2.8 결과와 대조
     python3 -I test/stackreg_parity.py [폴더]   (사례 생성, 기본 <tmpdir>/stackreg_parity)
     node test/test_stackreg.js [폴더]
   기준: 행렬·점 |차| ≤ 1e-9, 변환 영상 |차| ≤ 1e-6 × (max−min).
   정답 행렬(T)이 있는 사례는 우리 행렬의 정답 복원 오차도 따로 본다 (pystackreg 와 무관한 시험). */
const fs=require("fs"), path=require("path"), os=require("os");
const IC=require("./load_core.js");
const DIR=process.argv[2]||path.join(os.tmpdir(),"stackreg_parity");
const MAT_TOL=1e-9, IMG_TOL=1e-6;
let pass=0, fail=0; const fails=[];
function ok(cond,name,info){ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } if(process.env.V||!cond) console.log((cond?"  ok  ":"  FAIL")+" "+name+(info?"  ("+info+")":"")); }
if(!fs.existsSync(path.join(DIR,"manifest.json"))){ console.log("사례가 없습니다: python3 -I test/stackreg_parity.py 를 먼저 실행하세요 ("+DIR+")"); process.exit(1); }
const MAN=JSON.parse(fs.readFileSync(path.join(DIR,"manifest.json"),"utf8"));
const rd=f=>{ const b=fs.readFileSync(path.join(DIR,f)); return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength)); };
function toPlane(d,w,h,bit,off){ const N=w*h, a=new Float32Array(N); off=off||0; for(let i=0;i<N;i++){ a[i]=d[off+i]; if(a[i]!==d[off+i]) throw new Error("float32 로 정확히 담기지 않는 입력값"); } return {w,h,bit,data:a}; }
const matDiff=(A,B)=>{ let m=0; for(let i=0;i<B.length;i++) for(let j=0;j<B[i].length;j++){ const d=Math.abs(A[i][j]-B[i][j]); if(!(d<=m)) m=d; } return m; };
function imgDiff(got,ref){ let lo=Infinity,hi=-Infinity,mx=0; for(let i=0;i<ref.length;i++){ const r=ref[i]; if(r<lo) lo=r; if(r>hi) hi=r; const d=Math.abs(got[i]-r); if(!(d<=mx)) mx=d; } return {abs:mx,rel:mx/Math.max(hi-lo,1e-300)}; }
function gtErr(M,T,W,H){ const P=[[0,0],[W-1,0],[0,H-1],[W-1,H-1],[(W-1)/2,(H-1)/2]]; let dmax=0, dc=0;
  for(const [x,y] of P){ const d=Math.hypot(M[0][0]*x+M[0][1]*y+M[0][2]-(T[0][0]*x+T[0][1]*y+T[0][2]),M[1][0]*x+M[1][1]*y+M[1][2]-(T[1][0]*x+T[1][1]*y+T[1][2])); dmax=Math.max(dmax,d); dc=d; }
  const ang=Math.abs(Math.atan2(M[1][0],M[0][0])-Math.atan2(T[1][0],T[0][0]))*180/Math.PI; return {dmax,dcenter:dc,angle:ang}; }

const perMode={}; const note=(mode,k,v)=>{ const o=perMode[mode]||(perMode[mode]={mat:0,pts:0,img:0,n:0}); if(v>o[k]) o[k]=v; };
const gtRows=[]; const timing=[];
for(const c of MAN.cases){
  if(c.kind==="pair"){
    const ref=toPlane(rd(c.ref),c.w,c.h,c.bit), mov=toPlane(rd(c.mov),c.w,c.h,c.bit);
    const t0=process.hrtime.bigint(); const r=IC.turboReg(ref,mov,c.mode); const ms=Number(process.hrtime.bigint()-t0)/1e6;
    if(c.id.startsWith("big_")) timing.push([c.id,ms,c.py_ms]);
    const dm=matDiff(r.mat,c.mat), dp=Math.max(matDiff(r.refPts,c.refPts),matDiff(r.movPts,c.movPts));
    note(c.mode,"mat",dm); note(c.mode,"pts",dp); perMode[c.mode].n++;
    ok(dm<=MAT_TOL,c.id+" 행렬",`최대차 ${dm.toExponential(2)}`); ok(dp<=MAT_TOL,c.id+" 점",`최대차 ${dp.toExponential(2)}`);
    const out=IC.turboTransform(mov,c.mat,c.mode), di=imgDiff(out,rd(c.out)); note(c.mode,"img",di.rel);
    ok(out.length===c.w*c.h&&di.rel<=IMG_TOL,c.id+" 변환영상",`최대차/범위 ${di.rel.toExponential(2)}`);
    if(c.T){ const g=gtErr(r.mat,c.T,c.w,c.h); gtRows.push([c.id,c.mode,g]);
      const good=c.mode==="affine"?g.dmax<0.05:(g.dcenter<0.05&&g.angle<0.05);
      ok(good,c.id+" 정답복원",`dmax ${g.dmax.toFixed(4)} px, 중심 ${g.dcenter.toFixed(4)} px, 각 ${g.angle.toFixed(4)}°`); }
  } else {
    const d=rd(c.frames), N=c.w*c.h, planes=[]; for(let i=0;i<c.n;i++) planes.push(toPlane(d,c.w,c.h,c.bit,i*N));
    const mats=IC.stackRegMatrices(planes,c.mode,c.reference,c.nFrames);
    let dm=0; for(let i=0;i<c.n;i++) dm=Math.max(dm,matDiff(mats[i],c.mats[i]));
    note("stack:"+c.reference,"mat",dm); (perMode["stack:"+c.reference]).n++;
    ok(dm<=MAT_TOL,c.id+" 행렬",`최대차 ${dm.toExponential(2)}`);
    const outs=IC.stackRegApply(planes,c.mats,c.mode), ref=rd(c.out); let rel=0;
    for(let i=0;i<c.n;i++) rel=Math.max(rel,imgDiff(outs[i],ref.subarray(i*N,(i+1)*N)).rel);
    note("stack:"+c.reference,"img",rel);
    ok(rel<=IMG_TOL,c.id+" 변환스택",`최대차/범위 ${rel.toExponential(2)}`);
    // 우리 행렬로 변환해도 같은 결과
    const outs2=IC.stackRegApply(planes,mats,c.mode); let rel2=0; for(let i=0;i<c.n;i++) rel2=Math.max(rel2,imgDiff(outs2[i],ref.subarray(i*N,(i+1)*N)).rel);
    ok(rel2<=IMG_TOL,c.id+" 변환스택(우리 행렬)",`최대차/범위 ${rel2.toExponential(2)}`);
  }
}
// 오류 처리: bilinear + previous 는 pystackreg 처럼 거부
{ const p={w:30,h:30,bit:32,data:new Float32Array(900)}; let threw=false; try{ IC.stackRegMatrices([p,p],"bilinear","previous"); }catch(e){ threw=true; } ok(threw,"bilinear+previous 거부"); }

console.log("\n변환별 최대 차 (행렬 abs, 점 abs, 영상 차/범위):");
for(const k in perMode){ const o=perMode[k]; console.log(`  ${k.padEnd(18)} 사례 ${String(o.n).padStart(2)}  행렬 ${o.mat.toExponential(2)}  점 ${o.pts.toExponential(2)}  영상 ${o.img.toExponential(2)}`); }
console.log("정답 복원 오차 (우리 행렬):");
for(const [id,mode,g] of gtRows) console.log(`  ${id.padEnd(34)} dmax ${g.dmax.toFixed(4)} px  중심 ${g.dcenter.toFixed(4)} px  각 ${g.angle.toFixed(4)}°`);
// 속도: 512×512 다시 재기 (첫 실행의 JIT 예열 제외)
for(const c of MAN.cases.filter(c=>c.id.startsWith("big_"))){
  const ref=toPlane(rd(c.ref),c.w,c.h,c.bit), mov=toPlane(rd(c.mov),c.w,c.h,c.bit); const ts=[];
  for(let k=0;k<3;k++){ const t0=process.hrtime.bigint(); IC.turboReg(ref,mov,c.mode); ts.push(Number(process.hrtime.bigint()-t0)/1e6); }
  const first=timing.find(t=>t[0]===c.id);
  console.log(`속도 ${c.id}: 첫 실행 ${first[1].toFixed(0)} ms, 이후 ${ts.map(t=>t.toFixed(0)).join("/")} ms (pystackreg ${c.py_ms.toFixed(0)} ms)`);
}
console.log(`\n${pass} 통과, ${fail} 실패`);
if(fail){ for(const f of fails) console.log("  - "+f); process.exit(1); }
