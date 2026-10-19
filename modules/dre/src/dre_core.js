/* ============================================================
   INNO Analytics — DRE(rate equation) 엔진
   브라우저(메인·Worker)와 node 에서 같은 코드가 돈다. DOM 의존 없음.

   설계 요지 (기획서 v0.1)
   - 정상상태는 적분이 아니라 근 찾기: pseudo-transient continuation(ΨTC, implicit Euler + 증가하는 dt).
     앞 파워의 해에서 출발하므로 up/down 스윕이 실제 측정처럼 안정 branch를 따라가고, branch가 사라지면 점프한다.
   - 감도는 implicit function theorem(IFT): dn/dθ = −J⁻¹ ∂f/∂θ. f 가 과정 계수에 선형이라 ∂f/∂θ 는 정확하다.
   - 동역학은 L-안정 TR-BDF2 + 단계 배가 오차 제어.
   - 피팅은 log10 매개변수, log residual, 데이터셋별 배율 profiled out, 검출 하한 아래 점은 censored residual.
   ============================================================ */
(function(root){
"use strict";
const H_PLANCK=6.62607015e-34, C_LIGHT=2.99792458e8;
const photonEnergyJ = nm => H_PLANCK*C_LIGHT/(nm*1e-9);

/* ---------------- 식(expression) ---------------- */
/* 허용: 숫자, 이름, + - * / ^ ( ), 단항 −. 이름은 매개변수 또는 시료 변수(c 등). */
function compileExpr(src){
  const s=String(src).trim(); let i=0;
  const toks=[];
  while(i<s.length){
    const ch=s[i];
    if(/\s/.test(ch)){ i++; continue; }
    const m=/^(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)/.exec(s.slice(i));
    if(m){ toks.push({t:"n",v:parseFloat(m[1])}); i+=m[1].length; continue; }
    const id=/^[A-Za-z_][A-Za-z_0-9]*/.exec(s.slice(i));
    if(id){ toks.push({t:"id",v:id[0]}); i+=id[0].length; continue; }
    if("+-*/^()".includes(ch)){ toks.push({t:ch}); i++; continue; }
    throw new Error(`식 해석 불가: "${src}" (${i}번째 문자 '${ch}')`);
  }
  let p=0; const peek=()=>toks[p], eat=t=>{ if(!toks[p]||toks[p].t!==t) throw new Error(`식 오류: "${src}"`); p++; };
  const names=new Set();
  function prim(){
    const k=peek(); if(!k) throw new Error(`식이 끝남: "${src}"`);
    if(k.t==="n"){ p++; const v=k.v; return ()=>v; }
    if(k.t==="id"){ p++; names.add(k.v); const nm=k.v; return env=>{ const v=env[nm]; if(v===undefined) throw new Error(`정의되지 않은 이름: ${nm}`); return v; }; }
    if(k.t==="("){ p++; const e=add(); eat(")"); return e; }
    if(k.t==="-"){ p++; const e=pow(); return env=>-e(env); }
    if(k.t==="+"){ p++; return pow(); }
    throw new Error(`식 오류: "${src}"`);
  }
  function pow(){ const b=prim(); if(peek()&&peek().t==="^"){ p++; const e=pow(); return env=>Math.pow(b(env),e(env)); } return b; }
  function mul(){ let a=pow(); while(peek()&&(peek().t==="*"||peek().t==="/")){ const op=peek().t; p++; const b=pow(); const x=a; a= op==="*"?env=>x(env)*b(env):env=>x(env)/b(env); } return a; }
  function add(){ let a=mul(); while(peek()&&(peek().t==="+"||peek().t==="-")){ const op=peek().t; p++; const b=mul(); const x=a; a= op==="+"?env=>x(env)+b(env):env=>x(env)-b(env); } return a; }
  const fn=add(); if(p!==toks.length) throw new Error(`식 끝에 남은 부분: "${src}"`);
  fn.names=[...names]; fn.src=s; return fn;
}

/* ---------------- 모델 컴파일 ----------------
 spec = {
   species:[{id:"Tm", levels:["3H6","3F4","3H4"], cvar:"c"}],   // 첫 준위 = 바닥
   procs:[{type:"abs", sp:"Tm", from:"3H6", to:"3F4", sigma:"sG"},        // σ [cm²], 속도 = σ·φ
          {type:"dec", sp:"Tm", from:"3F4", to:"3H6", k:"w2r+w2nr"},      // k [s⁻¹]
          {type:"et", d:["Tm","3H4","3F4"], a:["Tm","3H6","3F4"], W:"s31", scale:"none"|"x"}],
   channels:[{id:"800", terms:[["Tm","3H4",1]]}],
   laser_nm:1064
 }
 et 의 scale:"x" — donor 감소 W·x_A, acceptor 증가 W·x_D (x = mole fraction, 시료 변수 cvar[%]/100). "none" 이면 둘 다 W. */
function compile(spec){
  const lv=[], idx={}, spOf=[], ground=[];
  spec.species.forEach((s,si)=>{
    s.levels.forEach((l,k)=>{ const key=s.id+":"+l; if(idx[key]!=null) throw new Error("중복 준위 "+key); idx[key]=lv.length; lv.push({sp:s.id,id:l,si,ground:k===0}); spOf.push(si); });
    ground[si]=idx[s.id+":"+s.levels[0]];
  });
  const L=(sp,l)=>{ const k=idx[sp+":"+l]; if(k==null) throw new Error(`준위 없음: ${sp} ${l}`); return k; };
  const procs=spec.procs.map((p,pi)=>{
    const o={type:p.type,src:p,pi};
    if(p.type==="abs"){ o.from=L(p.sp,p.from); o.to=L(p.sp,p.to); o.expr=compileExpr(p.sigma); }
    else if(p.type==="dec"){ o.from=L(p.sp,p.from); o.to=L(p.sp,p.to); o.expr=compileExpr(p.k); }
    else if(p.type==="et"){ o.di=L(p.d[0],p.d[1]); o.df=L(p.d[0],p.d[2]); o.ai=L(p.a[0],p.a[1]); o.af=L(p.a[0],p.a[2]);
      o.spD=p.d[0]; o.spA=p.a[0]; o.scale=p.scale||"none"; o.expr=compileExpr(p.W); }
    else throw new Error("알 수 없는 과정 "+p.type);
    return o;
  });
  const N=lv.length;
  const red=[]; for(let i=0;i<N;i++) if(!lv[i].ground) red.push(i);     /* 축약 좌표: 바닥 제외 */
  const redPos=new Array(N).fill(-1); red.forEach((g,k)=>redPos[g]=k);
  const channels=(spec.channels||[]).map(c=>({id:c.id,terms:c.terms.map(t=>[L(t[0],t[1]),t[2]==null?1:+t[2]])}));
  const names=new Set(); procs.forEach(p=>p.expr.names.forEach(n=>names.add(n)));
  const cvars={}; spec.species.forEach(s=>{ if(s.cvar) cvars[s.id]=s.cvar; });
  return {spec,lv,idx,N,red,redPos,spOf,ground,procs,channels,names:[...names],cvars,nSp:spec.species.length,L};
}

/* 과정 계수 (선형 값) — env = 매개변수 + 시료 변수 */
function coeffs(M,env){
  const c=new Float64Array(M.procs.length);
  M.procs.forEach((p,k)=>{ const v=p.expr(env); if(!(v>=0)||!isFinite(v)) throw new Error(`계수가 음수/비정상: ${p.expr.src} = ${v}`); c[k]=v; });
  return c;
}
/* mole fraction (et scale:x 용) */
function molFrac(M,env){
  const x=new Float64Array(M.nSp);
  M.spec.species.forEach((s,si)=>{ x[si]= s.cvar && env[s.cvar]!=null ? env[s.cvar]/100 : 1; });
  return x;
}

/* 전체 n → dn/dt, 각 준위의 flux 크기(수렴 척도) */
function rhs(M,c,x,n,phi,out,mag){
  const N=M.N; out.fill(0); if(mag) mag.fill(0);
  const spIdx=s=>M.spec.species.findIndex(q=>q.id===s);
  for(let k=0;k<M.procs.length;k++){
    const p=M.procs[k];
    if(p.type==="abs"||p.type==="dec"){
      const R=(p.type==="abs"?c[k]*phi:c[k])*n[p.from];
      out[p.from]-=R; out[p.to]+=R; if(mag){ mag[p.from]+=R; mag[p.to]+=R; }
    } else {
      const R=c[k]*n[p.di]*n[p.ai];
      let sD=1,sA=1;
      if(p.scale==="x"){ sD=x[p._sA]; sA=x[p._sD]; }
      out[p.di]-=sD*R; out[p.df]+=sD*R; out[p.ai]-=sA*R; out[p.af]+=sA*R;
      if(mag){ mag[p.di]+=sD*R; mag[p.df]+=sD*R; mag[p.ai]+=sA*R; mag[p.af]+=sA*R; }
    }
  }
}
function prep(M){
  M.procs.forEach(p=>{ if(p.type==="et"){ p._sD=M.spec.species.findIndex(q=>q.id===p.spD); p._sA=M.spec.species.findIndex(q=>q.id===p.spA); } });
  return M;
}
/* 축약 Jacobian (red×red): ∂f_i/∂n_j − ∂f_i/∂n_ground(j) */
function jacRed(M,c,x,n,phi,Jr){
  const N=M.N, F=new Float64Array(N*N);
  const add=(i,j,v)=>{ F[i*N+j]+=v; };
  for(let k=0;k<M.procs.length;k++){
    const p=M.procs[k];
    if(p.type==="abs"||p.type==="dec"){ const r=p.type==="abs"?c[k]*phi:c[k]; add(p.from,p.from,-r); add(p.to,p.from,r); }
    else {
      let sD=1,sA=1; if(p.scale==="x"){ sD=x[p._sA]; sA=x[p._sD]; }
      const dRdDi=c[k]*n[p.ai], dRdAi=c[k]*n[p.di];
      [[p.di,-sD],[p.df,sD],[p.ai,-sA],[p.af,sA]].forEach(([i,s])=>{ add(i,p.di,s*dRdDi); add(i,p.ai,s*dRdAi); });
    }
  }
  const R=M.red, m=R.length;
  for(let a=0;a<m;a++){ const i=R[a]; for(let b=0;b<m;b++){ const j=R[b]; const g=M.ground[M.spOf[j]]; Jr[a*m+b]=F[i*N+j]-F[i*N+g]; } }
  return Jr;
}
function fullFromRed(M,xr,n){
  n.fill(0); const sums=new Float64Array(M.nSp);
  M.red.forEach((g,k)=>{ n[g]=xr[k]; sums[M.spOf[g]]+=xr[k]; });
  for(let s=0;s<M.nSp;s++) n[M.ground[s]]=1-sums[s];
  return n;
}
/* 작은 조밀 선형계: 부분 피벗 LU */
function solve(A,b,m){
  const a=Float64Array.from(A), x=Float64Array.from(b);
  for(let k=0;k<m;k++){
    let p=k, mx=Math.abs(a[k*m+k]); for(let i=k+1;i<m;i++){ const v=Math.abs(a[i*m+k]); if(v>mx){ mx=v; p=i; } }
    if(mx===0) return null;
    if(p!==k){ for(let j=0;j<m;j++){ const t=a[k*m+j]; a[k*m+j]=a[p*m+j]; a[p*m+j]=t; } const t=x[k]; x[k]=x[p]; x[p]=t; }
    const d=a[k*m+k];
    for(let i=k+1;i<m;i++){ const f=a[i*m+k]/d; if(f===0) continue; for(let j=k;j<m;j++) a[i*m+j]-=f*a[k*m+j]; x[i]-=f*x[k]; }
  }
  for(let k=m-1;k>=0;k--){ let s=x[k]; for(let j=k+1;j<m;j++) s-=a[k*m+j]*x[j]; x[k]=s/a[k*m+k]; }
  return x;
}

/* ---------------- 정상상태 ----------------
   1) 앞 해(또는 바닥상태)에서 Newton + 선 탐색. 같은 branch가 아직 있으면 여기서 끝난다.
   2) 실패하면(branch가 사라짐 = threshold를 넘음) implicit Euler 로 pseudo-time 진행: dt 를 키워 가며 실제 동역학처럼
      다음 안정 상태로 흘러가게 한 뒤 Newton 으로 마무리. critical slowing down에도 dt 가 기하급수로 커져 안전하다. */
function steady(M,c,x,phi,n0,opt){
  opt=opt||{};
  const m=M.red.length, N=M.N, tol=opt.tol||1e-11;
  const n=Float64Array.from(n0||initGround(M)), f=new Float64Array(N), mag=new Float64Array(N);
  let xr=new Float64Array(m); M.red.forEach((g,k)=>xr[k]=n[g]);
  const Jr=new Float64Array(m*m), fr=new Float64Array(m);
  const evalF=v=>{ fullFromRed(M,v,n); rhs(M,c,x,n,phi,f,mag); M.red.forEach((g,k)=>fr[k]=f[g]); let r=0; for(const g of M.red){ const q=Math.abs(f[g])/(mag[g]+1e-300); if(q>r) r=q; } return r; };
  const valid=v=>{ const sums=new Float64Array(M.nSp); for(let k=0;k<m;k++){ if(v[k]<0) return false; sums[M.spOf[M.red[k]]]+=v[k]; } for(let s=0;s<M.nSp;s++) if(sums[s]>1+1e-12) return false; return true; };
  let rmax=1e-30; for(let k=0;k<M.procs.length;k++){ const p=M.procs[k]; rmax=Math.max(rmax,p.type==="abs"?c[k]*phi:c[k]); }
  /* 가중 2-노름 척도: 시작 시점의 flux 크기로 고정 */
  const newton=(y0,maxit)=>{
    let y=Float64Array.from(y0), r=evalF(y); const sc=new Float64Array(m); M.red.forEach((g,k)=>sc[k]=mag[g]+1e-300);
    const merit=()=>{ let s=0; for(let k=0;k<m;k++){ const q=fr[k]/sc[k]; s+=q*q; } return s; };
    let mf=merit();
    for(let it=0;it<maxit;it++){
      if(r<tol) return {y,r,ok:true,it};
      jacRed(M,c,x,n,phi,Jr); const A=new Float64Array(m*m); for(let a=0;a<m*m;a++) A[a]=-Jr[a];
      const d=solve(A,fr,m); if(!d) return {y,r,ok:false};
      let t=1, acc=false;
      for(let ls=0;ls<30;ls++){
        const yn=new Float64Array(m); for(let k=0;k<m;k++) yn[k]=y[k]+t*d[k];
        if(valid(yn)){ const rn=evalF(yn), mn=merit(); if(mn<mf*(1-1e-4*t)||rn<tol){ y=yn; r=rn; mf=mn; acc=true; break; } }
        t*=0.5;
      }
      if(!acc){ evalF(y); return {y,r,ok:r<tol,stalled:true}; }
    }
    return {y,r,ok:r<tol};
  };
  let res=newton(xr,opt.newtonIt||40);
  let path="newton";
  if(!res.ok){
    /* pseudo-time: y − dt·f(y) = y_prev */
    path="transient";
    let y=Float64Array.from(res.stalled?xr:res.y); if(!valid(y)) y=Float64Array.from(xr);
    let dt=opt.dt0||(1/rmax), T=0; const maxSteps=opt.maxSteps||600, dl=Math.min(opt.deadline||Infinity,opt.maxMs?Date.now()+opt.maxMs:Infinity)<Infinity?Math.min(opt.deadline||Infinity,opt.maxMs?Date.now()+opt.maxMs:Infinity):0;
    for(let step=0; step<maxSteps; step++){
      if(dl&&(step&7)===7&&Date.now()>dl) break;   /* 화면용 계산: 점 하나·스윕 전체 시간 상한 */
      /* implicit Euler 한 단계: Newton */
      let z=Float64Array.from(y), okStep=false;
      for(let it=0;it<25;it++){
        evalF(z); jacRed(M,c,x,n,phi,Jr);
        const A=new Float64Array(m*m); for(let a=0;a<m*m;a++) A[a]=-dt*Jr[a]; for(let a=0;a<m;a++) A[a*m+a]+=1;
        const rr=new Float64Array(m); for(let k=0;k<m;k++) rr[k]=y[k]-z[k]+dt*fr[k];
        const d=solve(A,rr,m); if(!d) break;
        let t=1; for(let k=0;k<m;k++) if(z[k]+d[k]<0) t=Math.min(t,0.9*z[k]/(-d[k]));
        let mx=0; for(let k=0;k<m;k++){ z[k]=Math.max(0,z[k]+t*d[k]); mx=Math.max(mx,Math.abs(d[k])/(z[k]+1e-18)); }
        if(mx<1e-9){ okStep=valid(z); break; }
      }
      if(!okStep){ dt*=0.25; if(dt<1e-18) break; continue; }
      y=z; T+=dt; dt*=2;
      if(dt>1e3*(opt.tmax||10)){ const nr=newton(y,40); if(nr.ok){ res=nr; break; } }
      if(step%20===19){ const nr=newton(y,20); if(nr.ok){ res=nr; break; } }
    }
    if(!res.ok){ const nr=newton(y,dl&&Date.now()>dl?8:60); res=nr; }
  }
  fullFromRed(M,res.y,n);
  return {n,ok:!!res.ok,res:res.r,path};
}
function initGround(M){ const n=new Float64Array(M.N); M.ground.forEach(g=>n[g]=1); return n; }

/* 파워 스윕: dir=+1 up, −1 down(up 끝 해에서 시작). I [W/cm²]
   opt(화면용 상한, 피팅은 쓰지 않음): maxFails — 수렴 실패가 이만큼 쌓이면 스윕 중단(남은 점 n=null, aborted),
   deadline — Date.now() 기준 끝 시각(넘으면 실패로 보고 중단), stOpt — steady 옵션(maxSteps 등). */
function sweep(M,env,I,dir,nm,opt){
  opt=opt||{};
  const c=coeffs(M,env), x=molFrac(M,env), Eph=photonEnergyJ(nm);
  const order=[...I.keys()].sort((a,b)=>I[a]-I[b]); if(dir<0) order.reverse();
  const out=new Array(I.length).fill(null); let n=null, fails=0, aborted=false; const failIdx=[];
  const maxF=opt.maxFails||Infinity, dl=opt.deadline||0, so=Object.assign({},opt.stOpt||{},dl?{deadline:dl}:{});
  if(dir<0){ /* 가장 높은 파워에서 먼저 upper branch로 */ const top=I[order[0]]; const hi=steady(M,c,x,top/Eph,null,so); n=hi.n; }
  for(const k of order){
    if(dl&&Date.now()>dl){ aborted=true; break; }
    const r=steady(M,c,x,I[k]/Eph,n,so); out[k]=r.n; n=r.n; if(!r.ok){ fails++; failIdx.push(k); if(fails>=maxF){ aborted=true; break; } }
  }
  const done=out.filter(v=>v).length;
  return {n:out,fails,failIdx,aborted,done};
}

/* ---------------- 감도 (IFT) ----------------
   θ = log10 매개변수(이름 목록). 반환: dn_full/dθ (N × P) */
function sensitivity(M,env,n,phi,pnames){
  const m=M.red.length, N=M.N, P=pnames.length;
  const c=coeffs(M,env), x=molFrac(M,env);
  const Jr=jacRed(M,c,x,n,phi,new Float64Array(m*m));
  /* ∂c_k/∂θ_j: 식의 수치 미분 (log10 공간, 중앙 차분) */
  const dc=[]; const h=1e-6;
  for(let j=0;j<P;j++){
    const nm=pnames[j], v=env[nm]; const e1=Object.assign({},env), e2=Object.assign({},env);
    e1[nm]=v*Math.pow(10,h); e2[nm]=v*Math.pow(10,-h);
    const c1=coeffs(M,e1), c2=coeffs(M,e2); const d=new Float64Array(c.length); for(let k=0;k<c.length;k++) d[k]=(c1[k]-c2[k])/(2*h); dc.push(d);
  }
  /* 과정별 g_k(n) 을 이용해 ∂f/∂θ_j */
  const S=new Float64Array(N*P);
  const g=new Float64Array(N);
  for(let j=0;j<P;j++){
    g.fill(0);
    for(let k=0;k<M.procs.length;k++){
      const d=dc[j][k]; if(d===0) continue; const p=M.procs[k];
      if(p.type==="abs"||p.type==="dec"){ const R=d*(p.type==="abs"?phi:1)*n[p.from]; g[p.from]-=R; g[p.to]+=R; }
      else { let sD=1,sA=1; if(p.scale==="x"){ sD=x[p._sA]; sA=x[p._sD]; } const R=d*n[p.di]*n[p.ai]; g[p.di]-=sD*R; g[p.df]+=sD*R; g[p.ai]-=sA*R; g[p.af]+=sA*R; }
    }
    const rhsv=new Float64Array(m); M.red.forEach((gi,k)=>rhsv[k]=-g[gi]);
    const dx=solve(Jr,rhsv,m); if(!dx) return null;
    const sums=new Float64Array(M.nSp);
    M.red.forEach((gi,k)=>{ S[gi*P+j]=dx[k]; sums[M.spOf[gi]]+=dx[k]; });
    for(let s=0;s<M.nSp;s++) S[M.ground[s]*P+j]=-sums[s];
  }
  return S;
}

/* ---------------- 동역학: TR-BDF2 ---------------- */
const GAM=2-Math.SQRT2, CG=GAM/2, D2=1/(GAM*(2-GAM)), E2=(1-GAM)*(1-GAM)/(GAM*(2-GAM));
function trbdf2Step(M,c,x,phi,xr,h){
  const m=M.red.length, N=M.N; const n=new Float64Array(N), f=new Float64Array(N), Jr=new Float64Array(m*m);
  const F=(v)=>{ fullFromRed(M,v,n); rhs(M,c,x,n,phi,f); const o=new Float64Array(m); M.red.forEach((g,k)=>o[k]=f[g]); return o; };
  const newton=(y0,rhsConst,coef)=>{ /* y − coef·h·F(y) = rhsConst */
    let y=Float64Array.from(y0);
    for(let it=0;it<12;it++){
      fullFromRed(M,y,n); rhs(M,c,x,n,phi,f); jacRed(M,c,x,n,phi,Jr);
      const A=new Float64Array(m*m); for(let a=0;a<m*m;a++) A[a]=-coef*h*Jr[a]; for(let a=0;a<m;a++) A[a*m+a]+=1;
      const r=new Float64Array(m); let nr=0;
      M.red.forEach((g,k)=>{ r[k]=rhsConst[k]-(y[k]-coef*h*f[g]); });
      const d=solve(A,r,m); if(!d) return null;
      let mx=0; for(let k=0;k<m;k++){ y[k]+=d[k]; mx=Math.max(mx,Math.abs(d[k])/(Math.abs(y[k])+1e-16)); }
      if(mx<1e-10) return y;
    }
    return y;
  };
  const f0=F(xr);
  const r1=new Float64Array(m); for(let k=0;k<m;k++) r1[k]=xr[k]+CG*h*f0[k];
  const yg=newton(xr,r1,CG); if(!yg) return null;
  const r2=new Float64Array(m); for(let k=0;k<m;k++) r2[k]=D2*yg[k]-E2*xr[k];
  const y1=newton(yg,r2,(1-GAM)/(2-GAM)); return y1;
}
/* t_end 까지 적분. onStep(t, n_full) 이 true 를 돌려주면 멈춤 */
function integrate(M,env,I,nm,tEnd,n0,opt){
  opt=opt||{};
  const c=coeffs(M,env), x=molFrac(M,env), phi=I/photonEnergyJ(nm);
  const m=M.red.length, N=M.N;
  let xr=new Float64Array(m); const n0f=n0||initGround(M); M.red.forEach((g,k)=>xr[k]=n0f[g]);
  let rmax=1e-30; for(let k=0;k<M.procs.length;k++){ const p=M.procs[k]; rmax=Math.max(rmax,p.type==="abs"?c[k]*phi:c[k]); }
  let t=0, h=opt.h0||Math.min(tEnd/1e4,1e-4/rmax);
  const rtol=opt.rtol||1e-6, atol=opt.atol||1e-15;
  const T=[0], Y=[fullFromRed(M,xr,new Float64Array(N))];
  let steps=0;
  const tOut=opt.tOut||null; let io=0; const outT=[], outN=[];
  while(t<tEnd && steps<(opt.maxSteps||20000)){
    if(t+h>tEnd) h=tEnd-t;
    let hitOut=false; if(tOut && io<tOut.length && t+h>=tOut[io]){ h=tOut[io]-t; hitOut=true; if(h<=0){ outT.push(t); outN.push(fullFromRed(M,xr,new Float64Array(N))); io++; continue; } }
    const yA=trbdf2Step(M,c,x,phi,xr,h);
    const yh=trbdf2Step(M,c,x,phi,xr,h/2); const yB=yh?trbdf2Step(M,c,x,phi,yh,h/2):null;
    if(!yA||!yB){ h*=0.25; continue; }
    let err=0; for(let k=0;k<m;k++){ const sc=atol+rtol*Math.max(Math.abs(yB[k]),Math.abs(xr[k])); err=Math.max(err,Math.abs(yB[k]-yA[k])/3/sc); }
    if(err<=1||h<1e-14){
      for(let k=0;k<m;k++) xr[k]=Math.max(0,yB[k]+(yB[k]-yA[k])/3);
      t+=h; steps++;
      const nf=fullFromRed(M,xr,new Float64Array(N)); T.push(t); Y.push(nf);
      if(hitOut){ outT.push(t); outN.push(nf); io++; }
      if(opt.onStep && opt.onStep(t,nf)) break;
    }
    h*=Math.min(4,Math.max(0.2,0.9*Math.pow(Math.max(err,1e-10),-1/3)));
  }
  return {t:T,n:Y,steps,outT,outN};
}
/* 상승 시간: 바닥 상태에서 켠 뒤 채널 신호가 frac×정상값에 닿는 시각 */
function channelValue(M,ch,n){ let s=0; for(const [i,w] of ch.terms) s+=w*n[i]; return s; }
function riseTime(M,env,I,nm,chIdx,frac){
  frac=frac||0.95; const ch=M.channels[chIdx];
  const c=coeffs(M,env), x=molFrac(M,env);
  const ss=steady(M,c,x,I/photonEnergyJ(nm),null); const target=frac*channelValue(M,ch,ss.n);
  let prevT=0, prevV=0, hit=null;
  /* 적분 끝은 '시간'이 아니라 목표 도달로 판정 — threshold 부근의 느려짐(critical slowing down)에 안전 */
  integrate(M,env,I,nm,1e3,null,{onStep:(t,n)=>{ const v=channelValue(M,ch,n); if(v>=target){ hit=prevT+(t-prevT)*(target-prevV)/((v-prevV)||1e-300); return true; } prevT=t; prevV=v; return false; }});
  return {t:hit,ss:ss.n};
}
/* 감쇠: 파워 I 정상상태에서 끈 뒤 신호가 1/e 로 떨어지는 시각 */
function decayTime(M,env,I,nm,chIdx){
  const ch=M.channels[chIdx]; const c=coeffs(M,env), x=molFrac(M,env);
  const ss=steady(M,c,x,I/photonEnergyJ(nm),null); const v0=channelValue(M,ch,ss.n), target=v0/Math.E;
  let prevT=0, prevV=v0, hit=null;
  integrate(M,env,0,nm,1e3,ss.n,{onStep:(t,n)=>{ const v=channelValue(M,ch,n); if(v<=target){ hit=prevT+(t-prevT)*(prevV-target)/((prevV-v)||1e-300); return true; } prevT=t; prevV=v; return false; }});
  return {t:hit};
}

/* ---------------- 관측 ----------------
   beam: "peak" (단일입자·peak intensity) | "gauss" (가우시안 빔 공간 적분: ∫₀^{I0} y(I) dI/I, Lee 2021) */
function logGrid(a,b,perDec){ const n=Math.max(2,Math.ceil(Math.log10(b/a)*perDec)+1); const o=[]; for(let i=0;i<n;i++) o.push(a*Math.pow(b/a,i/(n-1))); return o; }
function observe(M,env,I,dir,nm,chIdx,opt){
  const r=observeMulti(M,env,I,dir,nm,[chIdx],opt);
  const o={y:r.Y[0],dy:r.dY[0],fails:r.fails}; if((opt||{}).beam!=="gauss") o.n=r.n; return o;
}
/* 여러 채널을 스윕 한 번으로 (화면 시뮬레이션: 채널마다 스윕을 되풀이하지 않는다).
   반환 Y[채널][점], dY, n(I 위치의 population, 못 구한 점은 null), fails, aborted, done/total(스윕 격자) */
const cv=(M,ch,n)=>n?channelValue(M,ch,n):NaN;
function observeMulti(M,env,I,dir,nm,chList,opt){
  opt=opt||{}; const beam=opt.beam||"peak"; const P=opt.pnames||null;
  const swOpt={maxFails:opt.maxFails,deadline:opt.deadline,stOpt:opt.stOpt};
  if(beam==="peak"){
    const sw=sweep(M,env,I,dir,nm,swOpt);
    const Y=[], dY=[];
    let Sk=null; if(P) Sk=sw.n.map((n,k)=>n?sensitivity(M,env,n,I[k]/photonEnergyJ(nm),P):null);
    chList.forEach(chIdx=>{ const ch=M.channels[chIdx]; Y.push(sw.n.map(n=>cv(M,ch,n)));
      dY.push(P?Sk.map(S=>{ if(!S) return null; const o=new Float64Array(P.length); for(const [i,w] of ch.terms) for(let j=0;j<P.length;j++) o[j]+=w*S[i*P.length+j]; return o; }):null); });
    return {Y,dY,fails:sw.fails,n:sw.n,aborted:sw.aborted,done:sw.done,total:I.length,failIdx:sw.failIdx};
  }
  /* gauss: ∫₀^{I0} y(I) d ln I (Lee 2021 의 공간 적분과 같은 양).
     구간마다 거듭제곱 법칙 y ∝ I^s 로 적분하고, 격자 간격 h 와 2h 의 결과로 Richardson 외삽(오차 ∝ h²) —
     threshold의 급한 굽이에서도 정확하다(시험: 1920 점/decade 기준 오차 < 10⁻³). */
  const Imin=Math.min(...I), Imax=Math.max(...I);
  let Gf=logGrid(Imin*1e-3,Imax,opt.perDec||80); if(Gf.length%2===0) Gf=logGrid(Imin*1e-3,Imax*1.0000001,(opt.perDec||80)+0.37);
  const all=[...new Set(Gf.concat(I))].sort((p,q)=>p-q);
  const sw=sweep(M,env,all,dir,nm,swOpt);
  const Sall=P?all.map((g,k)=>sw.n[k]?sensitivity(M,env,sw.n[k],g/photonEnergyJ(nm),P):null):null;
  const idx=new Map(); all.forEach((g,k)=>idx.set(g,k));
  const Y=[], dY=[];
  chList.forEach(chIdx=>{ const r=gaussInt(M,M.channels[chIdx],I,Gf,all,sw.n,Sall,P); Y.push(r.y); dY.push(r.dy); });
  return {Y,dY,fails:sw.fails,n:I.map(v=>sw.n[idx.get(v)]),aborted:sw.aborted,done:sw.done,total:all.length,failIdx:sw.failIdx};
}
function gaussInt(M,ch,I,Gf,all,nAll,Sall,P){
  const val=new Map(), sen=new Map();
  all.forEach((g,k)=>{ val.set(g,cv(M,ch,nAll[k]));
    if(P){ const S=Sall[k]; const o=new Float64Array(P.length); if(S) for(const [i,w] of ch.terms) for(let j=0;j<P.length;j++) o[j]+=w*S[i*P.length+j]; sen.set(g,o); } });
  const seg=(ya,yb,L)=>{
    if(!(ya>0)||!(yb>0)) return [0.5*(ya+yb)*L,0.5*L,0.5*L];
    const lr=Math.log(yb/ya);
    if(Math.abs(lr)<1e-6){ const m=0.5*(ya+yb); return [m*L,0.5*L,0.5*L]; }
    const F=(yb-ya)*L/lr; return [F,-L/lr+F/(lr*ya),L/lr-F/(lr*yb)];
  };
  const integ=G=>{
    const cum=[0], cd=[P?new Float64Array(P.length):null];
    for(let k=1;k<G.length;k++){ const [F,da,db]=seg(val.get(G[k-1]),val.get(G[k]),Math.log(G[k]/G[k-1])); cum.push(cum[k-1]+F);
      if(P){ const o=new Float64Array(P.length), sa=sen.get(G[k-1]), sb=sen.get(G[k]); for(let j=0;j<P.length;j++) o[j]=cd[k-1][j]+da*sa[j]+db*sb[j]; cd.push(o); } }
    const y0=val.get(G[0]), y1=val.get(G[1]); const s0=Math.max(0.5,Math.log(y1/y0)/Math.log(G[1]/G[0])||1);
    return I.map(I0=>{ let k=0; while(k<G.length-1&&G[k+1]<=I0) k++;
      let v=cum[k]+y0/s0; let d=null;
      if(P){ d=new Float64Array(P.length); const s00=sen.get(G[0]); for(let j=0;j<P.length;j++) d[j]=cd[k][j]+s00[j]/s0; }
      if(I0>G[k]){ const [F,da,db]=seg(val.get(G[k]),val.get(I0),Math.log(I0/G[k])); v+=F; if(P){ const sa=sen.get(G[k]), sb=sen.get(I0); for(let j=0;j<P.length;j++) d[j]+=da*sa[j]+db*sb[j]; } }
      return [v,d]; });
  };
  const fine=integ(Gf), coarse=integ(Gf.filter((_,k)=>k%2===0));
  const y=fine.map((f,i)=>f[0]+(f[0]-coarse[i][0])/3);
  const dy=P?fine.map((f,i)=>{ const o=new Float64Array(P.length); for(let j=0;j<P.length;j++) o[j]=f[1][j]+(f[1][j]-coarse[i][1][j])/3; return o; }):null;
  return {y,dy};
}

/* ---------------- 피팅 문제 ----------------
 problem = {
   spec, nm, beam,
   params:{name:{v, lo, hi, fixed, per:"global"|"sample"}},   // v·lo·hi 는 선형값
   samples:[{id, vars:{c:8}, data:[
      {type:"scurve", ch:"800", dir:+1, I:[W/cm²], y:[배경 뺀 신호], sbg:배경 표준편차, k:3, relErr:0.05},
      {type:"rise", ch:"800", I:[], t:[s], frac:0.95, relErr:0.05},
      {type:"decay", ch:"800", I:[], t:[s], relErr:0.05} ]}]
 }
 자유 매개변수 벡터 θ = log10 값. 'sample' 매개변수는 시료마다 하나. */
function buildProblem(problem){
  const M=prep(compile(problem.spec));
  const free=[];
  Object.entries(problem.params).forEach(([nm,p])=>{
    if(p.fixed) return;
    if(p.per==="sample") problem.samples.forEach((s,si)=>free.push({nm,si,key:nm+"@"+s.id,lo:Math.log10(p.lo),hi:Math.log10(p.hi)}));
    else free.push({nm,si:-1,key:nm,lo:Math.log10(p.lo),hi:Math.log10(p.hi)});
  });
  M.names.forEach(n=>{ if(!(n in problem.params) && !problem.samples.every(s=>s.vars&&n in s.vars)) throw new Error(`식에 쓰였지만 값이 없는 이름: ${n}`); });
  return {M,free,problem};
}
function envFor(B,theta,si){
  const pr=B.problem, env={};
  Object.entries(pr.params).forEach(([nm,p])=>{ env[nm]=p.v; });
  Object.assign(env,pr.samples[si].vars||{});
  B.free.forEach((f,j)=>{ if(f.si===-1||f.si===si) env[f.nm]=Math.pow(10,theta[j]); });
  return env;
}
/* residual과 Jacobian. 배율은 데이터셋마다 profiled out(가중 평균) */
function residuals(B,theta,wantJ){
  const pr=B.problem, P=B.free.length, r=[], J=wantJ?[]:null, info=[];
  let fails=0;
  pr.samples.forEach((smp,si)=>{
    const env=envFor(B,theta,si);
    const pn=B.free.map((f,j)=>(f.si===-1||f.si===si)?f.nm:null);
    const active=pn.map((n,j)=>n?j:-1).filter(j=>j>=0), names=active.map(j=>pn[j]);
    smp.data.forEach((d,di)=>{
      const chIdx=B.M.channels.findIndex(c=>c.id===d.ch); if(chIdx<0) throw new Error("채널 없음: "+d.ch);
      if(d.type==="scurve"){
        const o=observe(B.M,env,d.I,d.dir||1,pr.nm,chIdx,{beam:pr.beam,pnames:wantJ?names:null}); fails+=o.fails;
        const rel=d.relErr||0.05, kc=d.k==null?3:d.k;
        /* 가중치: σ_ln² = rel² + (s_i/y)², s_i = √(배경σ² + 점별σ²) ; censored: y ≤ k·s_i */
        const sA=d.y.map((v,i)=>sAbs(d,i));
        const cens=d.y.map((v,i)=>!(v>kc*sA[i])||!(v>0));
        const sig=d.y.map((v,i)=>cens[i]?Math.sqrt(rel*rel+0.25):Math.sqrt(rel*rel+Math.pow(sA[i]/v,2)));
        const lm=o.y.map(v=>Math.log(Math.max(v,1e-300)));
        let sw=0, sl=0; d.y.forEach((v,i)=>{ if(cens[i]) return; const w=1/(sig[i]*sig[i]); sw+=w; sl+=w*(Math.log(v)-lm[i]); });
        const logS= sw>0 ? sl/sw : 0;
        /* d(ln m)/dθ 와 가중 평균 */
        let dlm=null, mean=null;
        if(wantJ){ dlm=o.y.map((v,i)=>{ const g=new Float64Array(P); if(o.dy&&o.dy[i]) active.forEach((j,a)=>g[j]=o.dy[i][a]/Math.max(v,1e-300)); return g; });
          mean=new Float64Array(P); d.y.forEach((v,i)=>{ if(cens[i]) return; const w=1/(sig[i]*sig[i]); for(let j=0;j<P;j++) mean[j]+=w*dlm[i][j]/sw; }); }
        d.y.forEach((v,i)=>{
          if(cens[i]){
            const L=Math.max(kc*sA[i],1e-300), pred=lm[i]+logS;
            const rr=pred>Math.log(L)?(pred-Math.log(L))/sig[i]:0; r.push(rr); info.push({si,di,i,cens:true});
            if(wantJ){ const g=new Float64Array(P); if(rr!==0) for(let j=0;j<P;j++) g[j]=(dlm[i][j]-mean[j])/sig[i]; J.push(g); }
          } else {
            r.push((Math.log(v)-lm[i]-logS)/sig[i]); info.push({si,di,i,cens:false});
            if(wantJ){ const g=new Float64Array(P); for(let j=0;j<P;j++) g[j]=-(dlm[i][j]-mean[j])/sig[i]; J.push(g); }
          }
        });
        d._fit={model:o.y.map(v=>v*Math.exp(logS)),logS};
      } else if(d.type==="rise"||d.type==="decay"){
        const f=(th)=>{ const e=envFor(B,th,si); return d.I.map(I=> d.type==="rise"? riseTime(B.M,e,I,pr.nm,chIdx,d.frac||0.95).t : decayTime(B.M,e,I,pr.nm,chIdx).t); };
        const tm=f(theta); const rel=d.relErr||0.05;
        const rr=d.t.map((t,i)=>tm[i]?Math.log(t/tm[i])/rel:10);
        let Jd=null;
        if(wantJ){ Jd=d.t.map(()=>new Float64Array(P)); const h=0.01;
          active.forEach(j=>{ const tp=Float64Array.from(theta), tq=Float64Array.from(theta); tp[j]+=h; tq[j]-=h; const a=f(tp), b=f(tq);
            d.t.forEach((t,i)=>{ Jd[i][j]= (a[i]&&b[i]) ? -(Math.log(a[i])-Math.log(b[i]))/(2*h)/rel : 0; }); }); }
        rr.forEach((v,i)=>{ r.push(v); info.push({si,di,i,kin:true}); if(wantJ) J.push(Jd[i]); });
        d._fit={model:tm};
      }
    });
  });
  return {r:Float64Array.from(r),J,info,fails};
}
const chi2=r=>{ let s=0; for(const v of r) s+=v*v; return s; };

/* ---------------- Levenberg–Marquardt (box constraints, log10 공간) ----------------
   F = ½‖r‖². (JᵀJ + μ·diag(JᵀJ)) h = −Jᵀr. gain ratio ρ 로 μ 갱신 (Nielsen 1999). 경계에 붙어 밖을 향하는 성분은 고정. */
function lm(B,theta0,opt){
  opt=opt||{}; const P=B.free.length; const lo=B.free.map(f=>f.lo), hi=B.free.map(f=>f.hi);
  let th=Float64Array.from(theta0).map((v,j)=>Math.min(hi[j],Math.max(lo[j],v)));
  let R; try{ R=residuals(B,th,true); }catch(e){ return {theta:th,chi2:Infinity,it:0,nev:1,err:String(e.message||e),conv:false}; }
  let F=0.5*chi2(R.r), mu=null, nu=2, it=0, nev=1, conv=false, why="maxit";
  const maxit=opt.maxit||100;
  for(it=0; it<maxit; it++){
    if(opt.cancel&&opt.cancel()){ why="cancel"; break; }
    const n=R.r.length, A=new Float64Array(P*P), g=new Float64Array(P);
    for(let i=0;i<n;i++){ const Ji=R.J[i]; for(let a=0;a<P;a++){ g[a]+=Ji[a]*R.r[i]; for(let b=0;b<=a;b++) A[a*P+b]+=Ji[a]*Ji[b]; } }
    for(let a=0;a<P;a++) for(let b=0;b<a;b++) A[b*P+a]=A[a*P+b];
    const act=[]; for(let a=0;a<P;a++){ const atLo=th[a]<=lo[a]+1e-12&&g[a]>0, atHi=th[a]>=hi[a]-1e-12&&g[a]<0; act.push(!(atLo||atHi)); }
    let gmax=0; for(let a=0;a<P;a++) if(act[a]) gmax=Math.max(gmax,Math.abs(g[a]));
    if(gmax<(opt.gtol||1e-8)*(1+F)){ conv=true; why="gradient"; break; }
    if(mu==null){ let d=0; for(let a=0;a<P;a++) d=Math.max(d,A[a*P+a]); mu=1e-3; }
    let accepted=false, stepNorm=0;
    for(let tr=0;tr<40;tr++){
      const Ad=new Float64Array(P*P), rhsv=new Float64Array(P);
      for(let a=0;a<P;a++){ for(let b=0;b<P;b++) Ad[a*P+b]=(act[a]&&act[b])?A[a*P+b]:0; Ad[a*P+a]= act[a]? A[a*P+a]*(1+mu)+1e-14 : 1; rhsv[a]= act[a]? -g[a] : 0; }
      const h=solve(Ad,rhsv,P); if(!h){ mu*=nu; nu*=2; continue; }
      const tn=th.map((v,a)=>Math.min(hi[a],Math.max(lo[a],v+h[a])));
      const hh=tn.map((v,a)=>v-th[a]);
      let pred=0; for(let a=0;a<P;a++){ pred-=hh[a]*g[a]; let s=0; for(let b=0;b<P;b++) s+=A[a*P+b]*hh[b]; pred-=0.5*hh[a]*s; }
      let Rn; try{ Rn=residuals(B,tn,true); nev++; }catch(e){ mu*=nu; nu*=2; continue; }
      const Fn=0.5*chi2(Rn.r);
      const rho=(F-Fn)/Math.max(pred,1e-300);
      if(isFinite(Fn) && Fn<F){
        stepNorm=Math.max(...hh.map(Math.abs));
        const dF=F-Fn; th=tn; R=Rn; F=Fn; mu*=Math.max(1/3,1-Math.pow(2*Math.min(rho,1)-1,3)); nu=2; accepted=true;
        if(stepNorm<(opt.xtol||1e-8) || dF<(opt.ftol||1e-12)*(1+F)){ conv=true; why="step"; }
        break;
      } else { mu*=nu; nu*=2; if(mu>1e16) break; }
    }
    if(conv) { it++; break; }
    if(!accepted){ why="stalled"; break; }
  }
  return {theta:th,chi2:2*F,it,nev,R,conv,why};
}

/* ---------------- multi-start ---------------- */
function mulberry(seed){ return ()=>{ seed|=0; seed=seed+0x6D2B79F5|0; let t=Math.imul(seed^seed>>>15,1|seed); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
function lhs(n,d,rnd){ const X=[]; for(let j=0;j<d;j++){ const perm=[...Array(n).keys()]; for(let i=n-1;i>0;i--){ const k=Math.floor(rnd()*(i+1)); [perm[i],perm[k]]=[perm[k],perm[i]]; } X.push(perm.map(p=>(p+rnd())/n)); } return [...Array(n).keys()].map(i=>X.map(col=>col[i])); }
function thetaFromParams(B){ return Float64Array.from(B.free.map(f=>Math.log10(B.problem.params[f.nm].v))); }
function multistart(B,opt,progress){
  opt=opt||{}; const n=opt.starts||60, rnd=mulberry(opt.seed||12345);
  const starts=[thetaFromParams(B)].concat(lhs(n-1,B.free.length,rnd).map(u=>Float64Array.from(u.map((v,j)=>B.free[j].lo+v*(B.free[j].hi-B.free[j].lo)))));
  const res=[];
  for(let k=0;k<starts.length;k++){
    if(opt.cancel&&opt.cancel()) break;
    const r=lm(B,starts[k],{maxit:opt.maxit||80,cancel:opt.cancel});
    res.push({k,chi2:r.chi2,theta:Array.from(r.theta),it:r.it,nev:r.nev,conv:r.conv,why:r.why});
    if(progress) progress({done:k+1,total:starts.length,best:Math.min(...res.map(x=>x.chi2))});
  }
  res.sort((a,b)=>a.chi2-b.chi2);
  return res;
}
/* 신뢰 threshold: Δχ² ≤ 3.84·s² (s² = max(1, χ²/자유도)) */
function dofInfo(B,best){ const R=residuals(B,Float64Array.from(best.theta),false); const n=R.r.length; const nScale=B.problem.samples.reduce((a,s)=>a+s.data.filter(d=>d.type==="scurve").length,0); const dof=Math.max(1,n-B.free.length-nScale); return {n,dof,s2:Math.max(1,best.chi2/dof),nScale}; }

/* ---------------- profile likelihood ---------------- */
function profile(B,best,j,opt,progress){
  opt=opt||{}; const span=opt.span||2, npt=opt.npts||17;
  const f=B.free[j], th0=Float64Array.from(best.theta);
  const grid=[]; for(let i=0;i<npt;i++){ const v=th0[j]-span+2*span*i/(npt-1); if(v>=f.lo-1e-9&&v<=f.hi+1e-9) grid.push(v); }
  if(!grid.some(v=>Math.abs(v-th0[j])<1e-9)) grid.push(th0[j]);
  grid.sort((a,b)=>a-b);
  const out=new Map();
  /* 최적점에서 양쪽으로 따뜻한 시작 */
  const sub=(thFix,warm)=>{ const Bs=Object.assign({},B,{free:B.free.map((q,k)=>k===j?Object.assign({},q,{lo:thFix,hi:thFix}):q)}); const t=Float64Array.from(warm); t[j]=thFix; return lm(Bs,t,{maxit:opt.maxit||60,cancel:opt.cancel}); };
  const i0=grid.findIndex(v=>Math.abs(v-th0[j])<1e-9);
  out.set(i0,{v:grid[i0],chi2:best.chi2,theta:Array.from(th0)});
  let warm=th0; for(let i=i0+1;i<grid.length;i++){ if(opt.cancel&&opt.cancel()) break; const r=sub(grid[i],warm); out.set(i,{v:grid[i],chi2:r.chi2,theta:Array.from(r.theta)}); warm=r.theta; if(progress) progress(); }
  warm=th0; for(let i=i0-1;i>=0;i--){ if(opt.cancel&&opt.cancel()) break; const r=sub(grid[i],warm); out.set(i,{v:grid[i],chi2:r.chi2,theta:Array.from(r.theta)}); warm=r.theta; if(progress) progress(); }
  const pts=[...out.entries()].sort((a,b)=>a[0]-b[0]).map(e=>e[1]);
  const thr=opt.thr||3.84;
  const inside=pts.filter(p=>p.chi2-best.chi2<=thr);
  return {pts,lo:Math.min(...inside.map(p=>p.v)),hi:Math.max(...inside.map(p=>p.v)),
    openLo:pts[0].chi2-best.chi2<=thr, openHi:pts[pts.length-1].chi2-best.chi2<=thr};
}

/* ---------------- FIM eigen-decomposition (Jacobi) ---------------- */
function eigSym(A,n){
  const a=Float64Array.from(A), V=new Float64Array(n*n); for(let i=0;i<n;i++) V[i*n+i]=1;
  for(let sweep=0;sweep<100;sweep++){
    let off=0; for(let p=0;p<n;p++) for(let q=p+1;q<n;q++) off+=a[p*n+q]*a[p*n+q];
    if(off<1e-30) break;
    for(let p=0;p<n;p++) for(let q=p+1;q<n;q++){
      const apq=a[p*n+q]; if(Math.abs(apq)<1e-300) continue;
      const th=(a[q*n+q]-a[p*n+p])/(2*apq), t=Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1)), c=1/Math.sqrt(t*t+1), s=t*c;
      for(let k=0;k<n;k++){ const akp=a[k*n+p], akq=a[k*n+q]; a[k*n+p]=c*akp-s*akq; a[k*n+q]=s*akp+c*akq; }
      for(let k=0;k<n;k++){ const apk=a[p*n+k], aqk=a[q*n+k]; a[p*n+k]=c*apk-s*aqk; a[q*n+k]=s*apk+c*aqk; }
      for(let k=0;k<n;k++){ const vkp=V[k*n+p], vkq=V[k*n+q]; V[k*n+p]=c*vkp-s*vkq; V[k*n+q]=s*vkp+c*vkq; }
    }
  }
  const ev=[...Array(n).keys()].map(i=>({val:a[i*n+i],vec:[...Array(n).keys()].map(k=>V[k*n+i])}));
  ev.sort((x,y)=>x.val-y.val); return ev;
}
function fim(B,theta){
  const R=residuals(B,Float64Array.from(theta),true), P=B.free.length, A=new Float64Array(P*P);
  R.J.forEach(Ji=>{ for(let a=0;a<P;a++) for(let b=0;b<P;b++) A[a*P+b]+=Ji[a]*Ji[b]; });
  return eigSym(A,P);
}
function aicc(chi,n,k){ return chi+2*k+(n-k-1>0?2*k*(k+1)/(n-k-1):Infinity); }

/* ---------------- 편의: 곡선 기울기 ---------------- */
function slopes(I,y){ return I.map((_,k)=>{ const a=Math.max(0,k-1), b=Math.min(I.length-1,k+1); return Math.log(y[b]/y[a])/Math.log(I[b]/I[a]); }); }

/* ---------------- 기본 모델: Lee et al. 2021 Tm³⁺ 3준위 ----------------
   공개 코드(github.com/nawhgnahc/Photon_Avalanche_DRE_calculation)의 식과 값을 그대로 옮김.
   코드 주석은 σ 단위를 cm² 로 적었으나 실제 계산은 m² 이므로(photon flux m⁻²s⁻¹) cm² 로 바꿔 ×10⁴. */
/* ---------------- 모델 적합성 검사 (V4 교훈) ----------------
   피팅이 끝나면 숫자는 언제나 나온다. 모델 구조가 데이터를 설명하지 못하면 그 숫자는 물성이 아니다.
   여기서는 '맞는다는 증명'이 아니라 '구조적으로 어긋나는 흔적'만 찾는다. 각 항목: level ok|warn|fail|na.
   (INNO Prism SegPA 의 s·I_th 추정과 무관 — 모든 비교는 같은 파워에서의 실측 대 DRE 모델이다.) */
function linfit(x,y){ const n=x.length; let sx=0,sy=0; for(let i=0;i<n;i++){ sx+=x[i]; sy+=y[i]; } const mx=sx/n,my=sy/n; let sxx=0,sxy=0; for(let i=0;i<n;i++){ sxx+=(x[i]-mx)**2; sxy+=(x[i]-mx)*(y[i]-my); }
  const b=sxy/sxx, a=my-b*mx; let ss=0; for(let i=0;i<n;i++) ss+=(y[i]-a-b*x[i])**2; return {b,a,sxx,se:n>2?Math.sqrt(ss/(n-2)/sxx):Infinity}; }
function interpLL(X,Y,x){ /* log-log 선형 보간, 범위 밖·비양수 null */ if(x<X[0]||x>X[X.length-1]) return null; let k=0; while(k<X.length-2&&X[k+1]<x) k++;
  const a=Y[k],b=Y[k+1]; if(!(a>0)||!(b>0)) return null; const f=Math.log(x/X[k])/Math.log(X[k+1]/X[k]||1); return {v:Math.exp(Math.log(a)+f*Math.log(b/a)),k}; }
/* 점별 절대 σ: 배경 σ(sbg)와 점별 σ(sy, 예: 반복 SD/√n)를 합친다 */
function sAbs(d,i){ const a=d.sbg||0, b=d.sy&&isFinite(d.sy[i])?d.sy[i]:0; return Math.sqrt(a*a+b*b); }
function adequacy(B,theta,opt){
  opt=opt||{}; const fmtI=opt.fmtI||(v=>v.toExponential(2)+" W/cm²");
  const th=Float64Array.from(theta); const pr=B.problem; const items=[];
  const R=residuals(B,th,false);
  /* A1 residual 수준 */
  const dof=Math.max(1,R.r.length-B.free.length-pr.samples.reduce((a,s)=>a+s.data.filter(d=>d.type==="scurve").length,0));
  const s2=chi2(R.r)/dof;
  items.push({key:"A1",title:"residual 크기",level:s2>5?"fail":s2>2?"warn":"ok",
    summary:`χ²/자유도 = ${s2.toFixed(2)}`,
    detail:s2>2?"잡음 모델(multiplicative noise·배경σ)이 맞다면 모델이 데이터를 설명하지 못한다는 뜻이다. 잡음을 과소 입력했을 가능성도 함께 확인할 것.":"입력한 잡음 수준과 맞는다."});
  /* 데이터셋별 정리 */
  const sets=[]; pr.samples.forEach((smp,si)=>smp.data.forEach((d,di)=>{ if(d.type!=="scurve"||!d._fit) return;
    const rel=d.relErr||0.05, kc=d.k==null?3:d.k;
    const idx=[...d.I.keys()].sort((a,b)=>d.I[a]-d.I[b]);
    const pts=idx.map(i=>{ const sa=sAbs(d,i); return {I:d.I[i],y:d.y[i],m:d._fit.model[i],cens:!(d.y[i]>kc*sa)||!(d.y[i]>0),sig:Math.sqrt(rel*rel+(d.y[i]>0?Math.pow(sa/d.y[i],2):0.25))}; });
    sets.push({si,di,d,pts,rel,chIdx:B.M.channels.findIndex(c=>c.id===d.ch),logS:d._fit.logS,label:`${smp.id} · ${d.ch} · ${d.dir<0?"down":"up"}`}); }));
  /* A2 residual 구조: 부호 연속(runs) 시험 */
  { const res=[]; let worst=null, lvl="ok";
    sets.forEach(S=>{ const u=S.pts.filter(p=>!p.cens); if(u.length<10) return;
      const r=u.map(p=>(Math.log(p.y)-Math.log(p.m))/p.sig); const sg=r.map(v=>v>=0);
      const np=sg.filter(Boolean).length, nm=sg.length-np; if(!np||!nm){ res.push({S,z:-Infinity}); lvl="fail"; return; }
      let runs=1; for(let i=1;i<sg.length;i++) if(sg[i]!==sg[i-1]) runs++;
      const n=np+nm, mu=2*np*nm/n+1, va=2*np*nm*(2*np*nm-n)/(n*n*(n-1)); const z=(runs-mu)/Math.sqrt(va);
      /* 가장 어긋난 5점 구간 */
      let best=null; for(let i=0;i+5<=r.length;i++){ const m=r.slice(i,i+5).reduce((a,b)=>a+b,0)/5; if(!best||Math.abs(m)>Math.abs(best.m)) best={m,a:u[i].I,b:u[i+4].I}; }
      res.push({S,z,runs,mu,best}); if(z<-3) lvl="fail"; else if((z<-2||z>3)&&lvl==="ok") lvl="warn";
      if(best&&(!worst||Math.abs(best.m)>Math.abs(worst.best.m))) worst={S,best}; });
    if(res.length){ const zmin=Math.min(...res.map(q=>q.z)), zmax=Math.max(...res.map(q=>q.z)), hiZ=res.filter(q=>q.z>3);
      items.push({key:"A2",title:"residual의 계통적 구조",level:lvl,
        summary:isFinite(zmin)?`부호 연속 시험 z = ${zmin.toFixed(1)}${zmax>zmin+0.05?` … ${zmax.toFixed(1)}`:""}`:"한 데이터셋의 residual이 모두 같은 부호",
        detail:(zmin<-2?"residual이 파워를 따라 무리 지어 같은 부호를 낸다 — 잡음이 아니라 곡선 모양이 다르다는 흔적. ":"")+
          (hiZ.length?`부호가 우연보다 훨씬 자주 바뀐다(z > 3: ${hiZ.map(q=>q.S.label).join(", ")}) — 같은 파워에서 값이 다른 두 집단(예: up·down 을 한 데이터셋에 넣음, 반복 측정이 어긋남)이 번갈아 섞였을 때 나타난다. `:"")+
          (lvl==="ok"?"residual이 무작위에 가깝다.":"")+(worst?` 가장 어긋난 구간: ${worst.S.label}, ${fmtI(worst.best.a)} – ${fmtI(worst.best.b)} (평균 residual ${worst.best.m.toFixed(1)}σ).`:""),
        rows:res.map(q=>({set:q.S.label,z:q.z}))}); }
    else items.push({key:"A2",title:"residual의 계통적 구조",level:"na",summary:"점 10개 이상인 S-curve 없음"}); }
  /* A3 구간별 모양 비교: 같은 파워 구간에서 실측 국소 log-기울기 vs 모델 */
  { const W=opt.win||5, flags=[]; let nWin=0;
    sets.forEach(S=>{ const u=S.pts.filter(p=>!p.cens);
      const X=u.map(p=>Math.log(p.I)), span=opt.span||0.15;
      for(let i=0;i<u.length;){
        let j=i; while(j<u.length-1&&(j-i+1<W||X[j]-X[i]<span)) j++;
        if(j-i+1<W||X[j]-X[i]<span*0.8) break;
        const w=u.slice(i,j+1); const x=w.map(p=>Math.log(p.I)); const nw=w.length;
        { let nx=i+1; while(nx<u.length&&X[nx]<X[i]+span/2) nx++; i=nx; }
        let jump=false; for(let k=1;k<nw;k++){ if(Math.abs(Math.log(w[k].m/w[k-1].m))>Math.LN2||Math.abs(Math.log(w[k].y/w[k-1].y))>Math.log(3)) jump=true; } if(jump) continue;
        const fd=linfit(x,w.map(p=>Math.log(p.y))), fm=linfit(x,w.map(p=>Math.log(p.m)));
        const se=Math.max(fd.se,S.rel/Math.sqrt(fd.sxx)); const dz=(fd.b-fm.b)/se; nWin++;
        if(Math.abs(fd.b-fm.b)>0.15&&Math.abs(dz)>3) flags.push({S,a:w[0].I,b:w[W-1].I,sd:fd.b,se,sm:fm.b,z:dz}); } });
    flags.sort((p,q)=>Math.abs(q.z)-Math.abs(p.z)); const f=flags[0];
    items.push({key:"A3",title:"구간별 모양 비교 (모델 vs 실측)",level:!nWin?"na":flags.length>=2?"fail":flags.length?"warn":"ok",
      summary:!nWin?"비교할 연속 구간 없음":`${nWin}개 구간 중 ${flags.length}개 어긋남`,
      detail:f?`가장 큰 차이: ${f.S.label}, ${fmtI(f.a)} – ${fmtI(f.b)} 에서 실측 국소 기울기 ${f.sd.toFixed(2)} ± ${f.se.toFixed(2)}, 모델 ${f.sm.toFixed(2)}. model class이 이 구간의 모양을 만들 수 없는지(예: bistable Tm 3준위의 upper branch 기울기 상한 ≈ 1, V4b) 확인할 것.`:"같은 파워 구간에서 모델과 실측의 기울기가 오차 안에서 같다. (PA 기울기 추정이 아니라 residual 진단이다.)",
      rows:flags.slice(0,6).map(q=>({set:q.S.label,a:q.a,b:q.b,sd:q.sd,se:q.se,sm:q.sm}))}); }
  /* A4 up/down 정합 + A5 두 branch 사이의 중간값 */
  { const pairs=[]; sets.forEach(U=>{ if(U.d.dir<0) return; sets.forEach(Dn=>{ if(Dn.d.dir>=0||Dn.si!==U.si||Dn.d.ch!==U.d.ch) return; pairs.push([U,Dn]); }); });
    if(!pairs.length) items.push({key:"A4",title:"up·down 정합",level:"na",summary:"같은 채널의 up·down 쌍 없음"});
    const jumpAt=pts=>{ let b=0,k=-1; for(let i=1;i<pts.length;i++){ if(pts[i-1].cens&&pts[i].cens) continue; const a=Math.abs(Math.log(Math.max(pts[i].y,1e-300)/Math.max(pts[i-1].y,1e-300))); if(a>b){b=a;k=i;} } return b>Math.log(5)?[pts[k-1].I,pts[k].I]:null; };
    const jumpM=pts=>{ let b=0,k=-1; for(let i=1;i<pts.length;i++){ const a=Math.abs(Math.log(pts[i].m/pts[i-1].m)); if(a>b){b=a;k=i;} } return b>Math.log(5)?[pts[k-1].I,pts[k].I]:null; };
    pairs.forEach(([U,Dn])=>{
      const eU=Math.exp(U.logS), eD=Math.exp(Dn.logS);
      const XU=U.pts.map(p=>p.I), YU=U.pts.map(p=>p.cens?0:p.y), MU=U.pts.map(p=>p.m/eU);
      let n=0,bad=0,zmax=0,wz=null;
      const L3=Math.log(3), sgp=(rel,sbg,y)=>Math.sqrt(rel*rel+Math.pow(sbg/y,2));
      Dn.pts.forEach((p,ip)=>{ if(p.cens) return; const q=interpLL(XU,YU,p.I); if(!q) return; const qm=interpLL(XU,MU,p.I); if(!qm) return;
        const k=q.k; if(Math.abs(Math.log(YU[k+1]/YU[k]))>L3||Math.abs(Math.log(MU[k+1]/MU[k]))>L3) return;   /* 점프(실측·모델)를 가로지르는 보간 제외 */
        const nb=[Dn.pts[ip-1],Dn.pts[ip+1]].filter(Boolean);                                                /* 접힘 바로 옆 점 제외 */
        if(nb.some(o=>o.cens||Math.abs(Math.log(o.y/p.y))>L3||Math.abs(Math.log(o.m/p.m))>L3)) return;
        const dd=Math.log(q.v/p.y), dm=Math.log(qm.v/(p.m/eD));
        const z=(dd-dm)/Math.hypot(sgp(U.rel,U.d.sbg||0,q.v),sgp(Dn.rel,Dn.d.sbg||0,p.y)); n++;
        if(Math.abs(z)>3){ bad++; if(Math.abs(z)>Math.abs(zmax)){ zmax=z; wz={I:p.I,dd,dm}; } } });
      const ju=jumpAt(U.pts), jd=jumpAt(Dn.pts), mu=jumpM(U.pts), md=jumpM(Dn.pts);
      const hd=ju&&jd?[ju[0]/jd[1],ju[1]/jd[0]]:null, hm=mu&&md?[mu[0]/md[1],mu[1]/md[0]]:null;
      let lvl=bad>=3?"fail":bad?"warn":"ok";
      const hTxt=`실측 P↑/P↓ ${hd?`${Math.max(1,hd[0]).toFixed(2)}–${hd[1].toFixed(2)}`:"(점프 없음)"} · 모델 ${hm?`${Math.max(1,hm[0]).toFixed(2)}–${hm[1].toFixed(2)}`:"(점프 없음)"}`;
      if((hd&&!hm)||(!hd&&hm)){ lvl="fail"; }
      else if(hd&&hm&&(hd[1]<hm[0]*0.97||hm[1]<hd[0]*0.97)) lvl="fail";
      items.push({key:"A4",title:"up·down 정합",level:lvl,summary:`${U.si===Dn.si?pr.samples[U.si].id:""} ${U.d.ch}: 겹치는 ${n}점 중 ${bad}점 불일치 · ${hTxt}`,
        detail:(bad?`up/down 비(배율 무관)가 모델과 가장 크게 다른 곳 ${fmtI(wz.I)}: 실측 ln(y↑/y↓) = ${wz.dd.toFixed(2)}, 모델 ${wz.dm.toFixed(2)}. `:"")+
          ((hd&&!hm)?"실측에는 히스테리시스가 있는데 모델에는 없다. ":(!hd&&hm)?"모델에는 히스테리시스가 있는데 실측에는 없다. ":"")+
          "점등 직후 느린 접근이나 초 단위 변화는 정상상태 모델 밖이다 — 전환 중인 점은 데이터에서 빼거나 표시할 것."});
    });
    /* A5 중간값 */
    let nInt=0, where=[];
    sets.forEach(S=>{ const env=envFor(B,th,S.si); const other=observe(B.M,env,S.d.I,-(S.d.dir||1),pr.nm,S.chIdx,{beam:pr.beam}).y; const e=Math.exp(S.logS);
      const kc=S.d.k==null?3:S.d.k;
      S.d.I.forEach((I,i)=>{ const y=S.d.y[i], sbg=sAbs(S.d,i); if(!(y>kc*sbg)||!(y>0)) return; const a=S.d._fit.model[i], b=other[i]*e; if(!(a>0&&b>0)) return;
        const lo=Math.max(Math.min(a,b),kc*sbg), hi=Math.max(a,b); if(hi/lo<3) return; const sg=3*Math.sqrt(S.rel*S.rel+Math.pow(sbg/y,2));
        if(Math.log(y/lo)>sg&&Math.log(hi/y)>sg){ nInt++; where.push(I); } }); });
    items.push({key:"A5",title:"두 branch 사이의 중간값",level:nInt>=2?"fail":nInt?"warn":"ok",summary:`${nInt}점`,
      detail:nInt?`모델이 두 안정 상태를 갖는 파워에서 실측이 어느 쪽에도 없다 (${where.slice(0,6).map(fmtI).join(", ")}). 적분 시간 안의 전환·깜박임(작은 계의 요동), 여러 방출체의 부분 점등, 또는 느린 변수를 뜻한다 — mean-field 한 방출체 모델 밖.`:"모델이 bistable인 구간에서 실측은 둘 중 한 branch 위에 있다."}); }
  /* A6 시간 트레이스: 과도 상승·점등 여부 (진단 전용 데이터) */
  pr.samples.forEach((smp,si)=>smp.data.forEach(d=>{ if(d.type!=="trace") return;
    const chIdx=B.M.channels.findIndex(c=>c.id===d.ch); const env=envFor(B,th,si);
    const n=d.t.length; const idx=[...d.t.keys()].filter(i=>d.t[i]>=0&&isFinite(d.y[i])).sort((a,b)=>d.t[a]-d.t[b]); const T=idx.map(i=>d.t[i]), Y=idx.map(i=>d.y[i]);
    const tail=Math.max(3,Math.floor(Y.length*0.2)); const fin=Y.slice(-tail).reduce((a,b)=>a+b,0)/tail;
    const sd=Math.sqrt(Y.slice(-tail).reduce((a,b)=>a+(b-fin)**2,0)/(tail-1))/Math.abs(fin||1);
    const k=Math.max(1,Math.round(Y.length/200)); const sm=Y.map((_,i)=>{ let s=0,c=0; for(let j=Math.max(0,i-(k>>1));j<=Math.min(Y.length-1,i+(k>>1));j++){ s+=Y[j]; c++; } return s/c; });
    const i90=sm.findIndex(v=>v>=0.9*fin); let pk=-Infinity,tpk=null; if(i90>=0) for(let i=i90;i<sm.length-tail;i++) if(sm[i]>pk){ pk=sm[i]; tpk=T[i]; }
    const ovD=i90>=0&&isFinite(pk)?pk/fin-1:0, thr=Math.max(0.03,4*sd/Math.sqrt(k));
    const head=Math.max(3,Math.floor(Y.length*0.02)); const y0=Y.slice(0,head).reduce((a,b)=>a+b,0)/head; const riseD=fin/Math.max(y0,1e-300);
    /* 모델: 바닥 상태에서 on */
    const ns=Math.min(400,T.length), tOut=[...new Set(Array.from({length:ns},(_,i)=>T[Math.floor(i*(T.length-1)/Math.max(1,ns-1))]).filter(t=>t>0))];
    const I=d.I[0]; const run=integrate(B.M,env,I,pr.nm,tOut[tOut.length-1],null,{tOut,maxSteps:200000});
    const ch=B.M.channels[chIdx]; const ym=run.outN.map(nv=>channelValue(B.M,ch,nv)); const fm=ym[ym.length-1]||0;
    const c0=coeffs(B.M,env), x0=molFrac(B.M,env); const ss=steady(B.M,c0,x0,I/photonEnergyJ(pr.nm),null); const yss=channelValue(B.M,ch,ss.n);
    const ovM=ym.length?Math.max(...ym)/Math.max(fm,1e-300)-1:0; const riseM=fm/Math.max(ym[0]||1e-300,1e-300);
    let lvl="ok", msg=[];
    if(riseD>5&&riseM<5){ lvl="fail"; msg.push("실측은 이 파워에서 켜지는데(최종/초기 "+riseD.toFixed(0)+"배) 모델은 lower branch에 머문다 — 이 트레이스 파워에서 모델의 점등 threshold이 더 높다."); }
    if(ovD>thr&&ovM<0.5*ovD){ lvl="fail"; msg.push(`실측 과도 상승 ${(100*ovD).toFixed(1)} % (${tpk!=null?(tpk*1e3).toFixed(1)+" ms":""}, 판정 threshold ${(100*thr).toFixed(1)} %), 모델 ${(100*ovM).toFixed(1)} %. 이 모델 구조로는 재현되지 않는다(예: Tm 3준위처럼 방출 준위가 빨라 저장 준위에 종속되면 동역학이 1차원이라 원리상 과도 상승이 없다). 모델 밖의 느린 변수(레이저 출력 과도, 가역 photobleaching 등)를 의심할 것 — 레이저 파형을 함께 기록하면 구분된다.`); }
    else if(ovD>thr) msg.push(`실측 과도 상승 ${(100*ovD).toFixed(1)} %, 모델 ${(100*ovM).toFixed(1)} % — 모델이 재현한다.`);
    if(!msg.length) msg.push(`과도 상승 없음 (실측 ${(100*ovD).toFixed(1)} %, threshold ${(100*thr).toFixed(1)} %).`);
    if(lvl==="ok"&&fm>0&&yss>0&&Math.abs(Math.log(fm/yss))>0.1) msg.push(`트레이스 끝에서 모델이 아직 정상상태의 ${(100*fm/yss).toFixed(0)} % — 트레이스가 짧다.`);
    items.push({key:"A6",title:"시간 트레이스 ("+smp.id+" · "+d.ch+" · "+fmtI(I)+")",level:lvl,summary:`과도 상승 실측 ${(100*ovD).toFixed(1)} % · 모델 ${(100*ovM).toFixed(1)} %`,detail:msg.join(" "),
      curve:{t:run.outT,y:ym.map(v=>v/Math.max(fm,1e-300))}});
  }));
  if(!items.some(q=>q.key==="A6")) items.push({key:"A6",title:"시간 트레이스",level:"na",summary:"진단용 트레이스 없음 — ‘시료·데이터’에서 ‘시간 트레이스’로 넣을 수 있다"});
  const order={fail:3,warn:2,ok:1,na:0}; const worst=items.reduce((a,q)=>order[q.level]>order[a]?q.level:a,"na");
  const verdict=worst==="fail"?"모델 구조가 데이터를 설명하지 못합니다 — 피팅값을 물성으로 보고하지 마세요.":worst==="warn"?"구조적 어긋남의 징후가 있습니다 — 표시된 항목을 확인한 뒤 해석하세요.":"구조적 모순은 찾지 못했습니다 (모델이 옳다는 증명은 아닙니다).";
  return {items,worst,verdict,s2};
}

/* ---------------- 피팅 상세 분석 재료 (Worker 에서 한 번에 계산) ----------------
   - Wald 공분산(곡률) · 상관행렬 · 매개변수별 정보 기여(데이터셋 × 파워)
   - multi-start 끝점의 해 집단(χ² 수준 → 같은 수준 안의 매개변수 퍼짐)
   - 과정별 flux 분해(대표 파워에서 각 과정이 담당하는 전이율 fraction) */
function procRates(M,c,x,n,phi){
  return M.procs.map((p,k)=>{
    if(p.type==="abs") return c[k]*phi*n[p.from];
    if(p.type==="dec") return c[k]*n[p.from];
    const f=p.scale==="x"?x[M.spOf[p.ai]]:1; return c[k]*f*n[p.di]*n[p.ai];
  });
}
function pinv(A,P){ /* 대칭 semi-definite pseudo-inverse + eigen-decomposition */
  const e=eigSym(A,P); const inv=new Float64Array(P*P); const lmax=Math.max(...e.map(q=>Math.abs(q.val)),1e-300);
  e.forEach(q=>{ if(q.val>lmax*1e-15) for(let a=0;a<P;a++) for(let b=0;b<P;b++) inv[a*P+b]+=q.vec[a]*q.vec[b]/q.val; });
  return {inv,eig:e};
}
function insight(B,best,ms,opt){
  opt=opt||{}; const P=B.free.length, th=Float64Array.from(best.theta);
  const R=residuals(B,th,true); const n=R.r.length;
  const info=dofInfo(B,best); const s2=info.s2;
  const A=new Float64Array(P*P); for(let i=0;i<n;i++){ const g=R.J[i]; for(let a=0;a<P;a++) for(let b=0;b<P;b++) A[a*P+b]+=g[a]*g[b]; }
  const {inv,eig}=pinv(A,P);
  /* 곡률이 사실상 0 인 방향(eigenvalue < 최대×1e-12)에 걸친 매개변수, 탐색 경계에 붙은 매개변수는 곡률 구간이 무의미 → ∞ */
  const lmax=Math.max(...eig.map(q=>Math.abs(q.val)),1e-300); const nullW=[...Array(P).keys()].map(a=>eig.filter(q=>q.val<=lmax*1e-15).reduce((t,q)=>t+q.vec[a]*q.vec[a],0));
  const atB=B.free.map((f,a)=>th[a]<=f.lo+1e-6||th[a]>=f.hi-1e-6);
  const se=[...Array(P).keys()].map(a=>(nullW[a]>0.01||atB[a])?Infinity:Math.sqrt(Math.max(inv[a*P+a],0)*s2));   /* log10 단위 표준오차 (θ 가 log10) */
  const corr=[]; for(let a=0;a<P;a++){ const row=[]; for(let b=0;b<P;b++){ const d=Math.sqrt(inv[a*P+a]*inv[b*P+b]); row.push(d>0?inv[a*P+b]/d:0); } corr.push(row); }
  /* 점별 정보 기여: J_ij² (대각 FIM fraction) — 데이터셋·파워 위치와 함께 */
  const pts=R.info.map((q,i)=>{ const smp=B.problem.samples[q.si], d=smp.data[q.di]; const I=d.I?d.I[q.i]:null;
    return {si:q.si,di:q.di,I,cens:!!q.cens,kin:!!q.kin,r:R.r[i],g:Array.from(R.J[i]).map(v=>v*v)}; });
  const setsInfo=[]; B.problem.samples.forEach((smp,si)=>smp.data.forEach((d,di)=>{ if(d.type==="trace") return;
    setsInfo.push({si,di,label:`${smp.id} · ${d.type==="scurve"?d.ch+" "+(d.dir<0?"down":"up"):d.type==="rise"?"상승 "+d.ch:"감쇠 "+d.ch}`,type:d.type}); }));
  /* 해 집단 */
  const thr=opt.thr||3.84*s2; const sorted=ms.slice().sort((a,b)=>a.chi2-b.chi2).filter(r=>isFinite(r.chi2));
  const clusters=[]; sorted.forEach(r=>{ let c=clusters.find(c=>Math.abs(r.chi2-c.chi2)<=Math.max(0.5,thr/4)); if(!c){ c={chi2:r.chi2,members:[]}; clusters.push(c); } c.members.push(r); });
  clusters.forEach(c=>{ c.count=c.members.length; c.lo=[...Array(P).keys()].map(j=>Math.min(...c.members.map(m=>m.theta[j]))); c.hi=[...Array(P).keys()].map(j=>Math.max(...c.members.map(m=>m.theta[j]))); c.rep=c.members[0].theta; delete c.members; });
  /* 과정별 flux 분해: 시료마다 S-curve 파워 범위의 대표점 5개, up branch */
  const flux=B.problem.samples.map((smp,si)=>{
    const Is=[]; smp.data.forEach(d=>{ if(d.type==="scurve") Is.push(...d.I); else if(d.I) Is.push(...d.I); }); if(!Is.length) return null;
    const lo=Math.min(...Is), hi=Math.max(...Is); const G=logGrid(lo,hi,40); const env=envFor(B,th,si);
    const c=coeffs(B.M,env), x=molFrac(B.M,env);
    const up=sweep(B.M,env,G,1,B.problem.nm).n, dn=sweep(B.M,env,G,-1,B.problem.nm).n;
    const pick=[0,Math.floor(G.length*0.35),Math.floor(G.length*0.6),Math.floor(G.length*0.8),G.length-1];
    const rows=[]; pick.forEach(k=>{ [["up",up[k]],["down",dn[k]]].forEach(([br,nv])=>{ if(br==="down"&&up[k].every((v,q)=>Math.abs(v-nv[q])<=1e-9*(1+Math.abs(v)))) return;
      const r=procRates(B.M,c,x,nv,G[k]/photonEnergyJ(B.problem.nm)); const tot=r.reduce((a,b)=>a+b,0)||1; rows.push({I:G[k],branch:br,share:r.map(v=>v/tot),pop:Array.from(nv)}); }); });
    return {si,rows}; });
  const procs=B.M.procs.map(p=>{ const s=p.src; const lv=i=>B.M.lv[i].id;
    const name=p.type==="abs"?`흡수 ${lv(p.from)}→${lv(p.to)}`:p.type==="dec"?`감쇠 ${lv(p.from)}→${lv(p.to)}`:`전달 ${lv(p.di)}→${lv(p.df)} ⇒ ${lv(p.ai)}→${lv(p.af)}`;
    return {name,expr:(s.sigma||s.k||s.W||""),note:s.note||"",names:p.expr.names}; });
  return {se,atBound:atB,corr,pts,sets:setsInfo,clusters,flux,procs,s2,n,thr,levelNames:B.M.lv.map(l=>l.sp+" "+l.id)};
}
/* χ² slice: 두 매개변수 격자, 나머지는 최적값 고정 */
function slice2D(B,theta,j,k,opt){
  opt=opt||{}; const n=opt.n||17, sj=opt.spanJ||1.5, sk=opt.spanK||1.5; const th=Float64Array.from(theta);
  const xs=[...Array(n).keys()].map(i=>th[j]-sj+2*sj*i/(n-1)), ys=[...Array(n).keys()].map(i=>th[k]-sk+2*sk*i/(n-1));
  const z=[]; for(let b=0;b<n;b++){ const row=[]; for(let a=0;a<n;a++){ const t=Float64Array.from(th); t[j]=xs[a]; t[k]=ys[b]; let v; try{ v=chi2(residuals(B,t,false).r); }catch(e){ v=NaN; } row.push(v); if(opt.progress) opt.progress(b*n+a+1,n*n); } z.push(row); }
  return {xs,ys,z};
}

function templateLee(){
  return {
    name:"Tm³⁺ 3준위 PA (Lee et al. 2021)",
    spec:{species:[{id:"Tm",levels:["3H6","3F4","3H4"],cvar:"c"}],
      procs:[
        {type:"abs",sp:"Tm",from:"3H6",to:"3F4",sigma:"sG",note:"non-resonant GSA (³H₆→³H₅→³F₄)"},
        {type:"abs",sp:"Tm",from:"3F4",to:"3H4",sigma:"sE",note:"ESA ³F₄→³F₂,₃ (→³H₄)"},
        {type:"dec",sp:"Tm",from:"3F4",to:"3H6",k:"w2r+w2nr",note:"³F₄ 감쇠"},
        {type:"dec",sp:"Tm",from:"3H4",to:"3H6",k:"(1-b32)*w3r",note:"³H₄→³H₆ 복사"},
        {type:"dec",sp:"Tm",from:"3H4",to:"3F4",k:"b32*w3r+w3nr",note:"³H₄→³F₄ (복사 분기 + 비복사)"},
        {type:"et",d:["Tm","3H4","3F4"],a:["Tm","3H6","3F4"],W:"k31*c^2",scale:"none",note:"CR (³H₄,³H₆)→(³F₄,³F₄)"},
        {type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3F4"],W:"kQ22*c^3/(c^2+4.3^2)",scale:"none",note:"ETU 손실형 (³F₄,³F₄)→(³H₆,³H₅→³F₄)"},
        {type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3H4"],W:"kQ23*c^3/(c^2+4.3^2)",scale:"none",note:"ETU (³F₄,³F₄)→(³H₆,³H₄)"}],
      channels:[{id:"800 nm (³H₄)",terms:[["Tm","3H4",1]]},{id:"1.8 µm (³F₄)",terms:[["Tm","3F4",1]]}]},
    nm:1064,
    params:{
      sG:{v:6.0e-25,lo:1e-28,hi:1e-21,unit:"cm²",src:"Lee 2021 코드 f_ab_1·0.12×10⁻²⁵ m²"},
      sE:{v:6.4e-21,lo:1e-24,hi:1e-18,unit:"cm²",src:"Lee 2021 코드 f_ab_2·3.2×10⁻²⁵ m²"},
      w2r:{v:162.6,lo:10,hi:1e4,unit:"s⁻¹",fixed:true,src:"Lee 2021 코드 (Kwock 2021·Villanueva-Delgado 2016 은 83.3 s⁻¹ = τ 12 ms)"},
      w2nr:{v:40.65,lo:0.01,hi:1e5,unit:"s⁻¹",per:"sample",src:"Lee 2021 코드 Sample No. 4 (162.6×0.25+0.00239)"},
      w3r:{v:636.01,lo:10,hi:1e5,unit:"s⁻¹",fixed:true,src:"Lee 2021 코드"},
      w3nr:{v:87.32,lo:0.1,hi:1e6,unit:"s⁻¹",per:"sample",src:"Lee 2021 코드 Sample No. 4 (6.02+162.6×0.5)"},
      b32:{v:0.144,lo:0.01,hi:0.9,unit:"",fixed:true,src:"Lee 2021 코드"},
      k31:{v:160,lo:0.1,hi:1e6,unit:"s⁻¹·%⁻²",src:"Lee 2021 코드 s31 = 0.1·1.6·c²·1000"},
      kQ22:{v:25.6,lo:1e-3,hi:1e5,unit:"s⁻¹·%⁻¹",src:"Lee 2021 코드 0.08·320·c³/(c²+4.3²)"},
      kQ23:{v:9.0,lo:1e-3,hi:1e5,unit:"s⁻¹·%⁻¹",src:"Lee 2021 코드 0.1·90·c³/(c²+4.3²)"}},
    samples:[{id:"S1",vars:{c:8}}]
  };
}

/* Tm³⁺ 5준위(³H₆ ³F₄ ³H₅ ³H₄ ¹G₄) / 6준위(+¹D₂) — 1064 nm 광자 사태.
   값의 근거: Chen 2025(NaLuF₄:Tm 8준위 모델, SI Table 3; 프로젝트 검증노트 V4), Lee 2021(JO 복사율·Q 계수, SI Table S4–S7).
   에너지 불일치(Carnall 무게중심·Stark 범위, Energy Diagram 데이터): GSA ³H₆→³H₅ 는 광자가 830–1250 cm⁻¹ 남는 포논 방출형,
   ESA1 ³F₄→³F₂ resonant(창 안), ESA2 ³H₄→¹G₄ 는 440–1250 cm⁻¹ 남음, CR1 ³H₄+³H₆→³F₄+³F₄ 는 ≥670 cm⁻¹ 남음(포논 보조),
   c2 ³H₄→³F₄ ⇒ ¹G₄→¹D₂ 는 창이 겹쳐 resonant, ETU ³F₄+³F₄→³H₆+³H₄ 는 ≥670 cm⁻¹ endothermic. */
function templateTm(level){
  const six=level===6;
  const lv=six?["3H6","3F4","3H5","3H4","1G4","1D2"]:["3H6","3F4","3H5","3H4","1G4"];
  const N="NRE*c/100";   /* Tm number density [cm⁻³] — 전달 계수(cm³ s⁻¹) × number density = 정규화 population 기준 속도 [s⁻¹] */
  const procs=[
    {type:"abs",sp:"Tm",from:"3H6",to:"3H5",sigma:"sG",note:"GSA ³H₆→³H₅ (포논 방출형, non-resonant)"},
    {type:"abs",sp:"Tm",from:"3F4",to:"3H4",sigma:"sE1",note:"ESA1 ³F₄→³F₂,₃ (resonant) → ³H₄ 로 빠른 이완을 가정"},
    {type:"abs",sp:"Tm",from:"3H4",to:"1G4",sigma:"sE2",note:"ESA2 ³H₄→¹G₄ (포논 방출형)"},
    {type:"dec",sp:"Tm",from:"3F4",to:"3H6",k:"w2r+w2nr",note:"³F₄ 감쇠 (복사 + 시료별 비복사)"},
    {type:"dec",sp:"Tm",from:"3H5",to:"3F4",k:"w35*(1-b51)",note:"³H₅ → ³F₄ (다포논, 빠름)"},
    {type:"dec",sp:"Tm",from:"3H5",to:"3H6",k:"w35*b51",note:"³H₅ → ³H₆"},
    {type:"dec",sp:"Tm",from:"3H4",to:"3H6",k:"(1-b32)*w3r",note:"³H₄→³H₆ 복사 (800 nm)"},
    {type:"dec",sp:"Tm",from:"3H4",to:"3F4",k:"b32*w3r",note:"³H₄→³F₄ 복사 (1.47 µm)"},
    {type:"dec",sp:"Tm",from:"3H4",to:"3H5",k:"w3nr",note:"³H₄ 비복사 (시료별) → ³H₅"},
    {type:"dec",sp:"Tm",from:"1G4",to:"3H6",k:"0.30*w6",note:"¹G₄→³H₆ (475 nm)"},
    {type:"dec",sp:"Tm",from:"1G4",to:"3F4",k:"0.10*w6",note:"¹G₄→³F₄ (650 nm)"},
    {type:"dec",sp:"Tm",from:"1G4",to:"3H5",k:"0.45*w6",note:"¹G₄→³H₅ (780 nm)"},
    {type:"dec",sp:"Tm",from:"1G4",to:"3H4",k:"0.15*w6",note:"¹G₄→³H₄ (+³F₂,₃ 경유 0.03)"},
    {type:"et",d:["Tm","3H4","3F4"],a:["Tm","3H6","3F4"],W:`cCR*${N}`,scale:"none",note:"CR1 (³H₄,³H₆)→(³F₄,³F₄) — 사태 루프"},
    {type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3H4"],W:`cUC23*${N}`,scale:"none",note:"ETU (³F₄,³F₄)→(³H₆,³H₄) — endothermic"},
    {type:"et",d:["Tm","3F4","3H6"],a:["Tm","3F4","3H5"],W:`cUC22*${N}`,scale:"none",note:"ETU (³F₄,³F₄)→(³H₆,³H₅→³F₄) — 손실형"}];
  if(six) procs.push(
    {type:"et",d:["Tm","3H4","3F4"],a:["Tm","1G4","1D2"],W:`cCR2*${N}`,scale:"none",note:"CR2 (³H₄,¹G₄)→(³F₄,¹D₂) — resonant, 450 nm 경로"},
    {type:"dec",sp:"Tm",from:"1D2",to:"3H6",k:"0.36*w7",note:"¹D₂→³H₆ (360 nm)"},
    {type:"dec",sp:"Tm",from:"1D2",to:"3F4",k:"0.46*w7",note:"¹D₂→³F₄ (450 nm)"},
    {type:"dec",sp:"Tm",from:"1D2",to:"3H5",k:"0.02*w7",note:"¹D₂→³H₅"},
    {type:"dec",sp:"Tm",from:"1D2",to:"3H4",k:"0.15*w7",note:"¹D₂→³H₄ (+³F₂,₃ 경유 0.07)"},
    {type:"dec",sp:"Tm",from:"1D2",to:"1G4",k:"0.01*w7",note:"¹D₂→¹G₄"});
  /* 채널 가중치 = 기본값에서의 복사(방출) 속도 [s⁻¹] — 채널끼리 비율에만 의미가 있다 */
  const ch=[{id:"800 nm (³H₄→³H₆)",terms:[["Tm","3H4",0.856*636.01]]},
    {id:"475 nm (¹G₄→³H₆)",terms:[["Tm","1G4",0.30*2e4]]},
    {id:"NIR 필터 (³H₄ 800 + ¹G₄→³H₅ 780)",terms:[["Tm","3H4",0.856*636.01],["Tm","1G4",0.45*2e4]]},
    {id:"1.8 µm (³F₄)",terms:[["Tm","3F4",162.6]]}];
  if(six) ch.splice(2,0,{id:"450 nm (¹D₂→³F₄)",terms:[["Tm","1D2",0.46*2e4]]},{id:"청색 450+475",terms:[["Tm","1D2",0.46*2e4],["Tm","1G4",0.30*2e4]]});
  const params={
    sG:{v:7e-25,lo:1e-28,hi:1e-21,unit:"cm²",src:"Chen 2025 SI Table 3 σ_GSA (Lee 2021: 6.0×10⁻²⁵)"},
    sE1:{v:7e-21,lo:1e-24,hi:1e-18,unit:"cm²",src:"Chen 2025 σ_ESA1 (E/G=1; Lee 2021: 6.4×10⁻²¹)"},
    sE2:{v:7e-22,lo:1e-25,hi:1e-19,unit:"cm²",src:"Chen 2025 σ_ESA2 = σ_ESA1/10"},
    w2r:{v:162.6,lo:10,hi:1e4,unit:"s⁻¹",fixed:true,src:"Lee 2021 JO (SI Table S6)"},
    w2nr:{v:40.65,lo:0.01,hi:1e6,unit:"s⁻¹",per:"sample",src:"placeholder: Lee 2021 Sample 4. ³F₄ 수명으로 시료마다 정할 것 (Chen 2025 NaLuF₄ 15 %: τ³F₄ 26 µs)"},
    w35:{v:1e5,lo:1e3,hi:1e8,unit:"s⁻¹",fixed:true,src:"Chen 2025 w3 (³H₅ 전체 이완)"},
    b51:{v:0.1,lo:0.01,hi:0.9,unit:"",fixed:true,src:"Chen 2025 b31 (³H₅→³H₆)"},
    w3r:{v:636.01,lo:10,hi:1e5,unit:"s⁻¹",fixed:true,src:"Lee 2021 JO"},
    b32:{v:0.144,lo:0.01,hi:0.9,unit:"",fixed:true,src:"Lee 2021 JO (³H₄→³F₄ 분기; ³H₅ fraction 포함)"},
    w3nr:{v:87.32,lo:0.1,hi:1e7,unit:"s⁻¹",per:"sample",src:"placeholder: Lee 2021 Sample 4. ³H₄ 수명(저농도)으로 정할 것"},
    w6:{v:2e4,lo:1e2,hi:1e7,unit:"s⁻¹",src:"Chen 2025 w6 (W=1); 분기 0.30/0.10/0.45/0.15 = Chen b6x"},
    cCR:{v:5.4e-17,lo:1e-21,hi:1e-13,unit:"cm³ s⁻¹",src:"Chen 2025 c1 (CR=1, NaYF₄). NaLuF₄ 맞춤은 CR=6 → 3.2×10⁻¹⁶"},
    cUC23:{v:5.1e-20,lo:1e-24,hi:1e-15,unit:"cm³ s⁻¹",src:"Lee 2021 a_inv=9.0 을 c=8 % 에서 환산 (kQ·c³/(c²+4.3²)/N_Tm)"},
    cUC22:{v:1.45e-19,lo:1e-24,hi:1e-15,unit:"cm³ s⁻¹",src:"Lee 2021 a_uc=25.6 을 c=8 % 에서 환산"},
    NRE:{v:1.42e22,lo:1e22,hi:2e22,unit:"cm⁻³",fixed:true,src:"β-NaLuF₄ 희토류 자리 밀도 1.5 f.u./105.9 Å³ (Chen 2025 Rietveld). NaYF₄ 는 1.37×10²²"}};
  if(six){ params.cCR2={v:8e-18,lo:1e-22,hi:1e-14,unit:"cm³ s⁻¹",src:"Chen 2025 c2 (CR=1)"};
    params.w7={v:2e4,lo:1e2,hi:1e7,unit:"s⁻¹",src:"Chen 2025 w7 (W=1); 분기 = Chen b7x"}; }
  return {name:six?"Tm³⁺ 6준위 PA (³H₆–¹D₂, 1064 nm)":"Tm³⁺ 5준위 PA (³H₆–¹G₄, 1064 nm)",
    spec:{species:[{id:"Tm",levels:lv,cvar:"c"}],procs,channels:ch},nm:1064,params,samples:[{id:"S1",vars:{c:15}}],
    note:"Chen 2025 8준위 모델을 줄인 것: ³F₂,₃ 는 ³H₄ 로 즉시 이완(ESA1 도착점), "+(six?"¹I₆ 경로(c3)는 제외.":"¹D₂·¹I₆ 경로는 제외 — 450 nm 가 필요하면 6준위.")+" 전달 계수는 number density mean-field(W = c·N_Tm)."};
}

/* 유도량(예측값): 개별 매개변수보다 데이터가 잘 정하는 조합들.
   - 흡수 속도 σ·φ(I_ref): intensity 환산 오차와 무관(nominal intensity 기준)
   - 전달 fraction η = W/(W + donor 준위의 감쇠 합) (acceptor가 바닥 준위, 저여기 극한)
   - threshold I_th(첫 채널 최대 기울기 위치), 2·I_th 에서 각 준위 population */
function derived(B,theta,opt){
  opt=opt||{}; const si=opt.si||0; const th=Float64Array.from(theta); const pr=B.problem; const M=B.M;
  const env=envFor(B,th,si); const c=coeffs(M,env), x=molFrac(M,env);
  const Is=[]; pr.samples[si].data.forEach(d=>{ if(d.I) Is.push(...d.I.filter(v=>v>0)); });
  if(!Is.length) return [];
  const lo=Math.min(...Is), hi=Math.max(...Is), Iref=opt.Iref||Math.sqrt(lo*hi), phi=Iref/photonEnergyJ(pr.nm);
  const lvn=i=>M.lv[i].id; const out=[];
  M.procs.forEach((p,k)=>{ if(p.type==="abs") out.push({key:"abs"+k,label:`흡수 속도 ${lvn(p.from)}→${lvn(p.to)} @ I_ref`,unit:"s⁻¹",v:c[k]*phi,calFree:true}); });
  const decOut=i=>M.procs.reduce((a,p,k)=>a+(p.type==="dec"&&p.from===i?c[k]:0),0);
  M.procs.forEach((p,k)=>{ if(p.type!=="et"||!M.lv[p.ai].ground) return; const f=p.scale==="x"?x[M.spOf[p.ai]]:1; const W=c[k]*f; const dsum=decOut(p.di);
    out.push({key:"eta"+k,label:`전달 fraction ${lvn(p.di)}→${lvn(p.df)} ⇒ ${lvn(p.ai)}→${lvn(p.af)} (저여기)`,unit:"",v:W/(W+dsum)}); });
  const chIdx=Math.max(0,opt.ch||0), ch=M.channels[chIdx];
  if(ch){ const G=logGrid(lo,hi,80); const sw=sweep(M,env,G,1,pr.nm).n; const y=sw.map(n=>channelValue(M,ch,n)); const sl=slopes(G,y); let k=0; sl.forEach((v,i)=>{ if(v>sl[k]) k=i; });
    const Ith=G[k]; out.push({key:"Ith",label:`threshold I_th (${ch.id} 최대 기울기 위치)`,unit:"W/cm²",v:Ith,calNote:true});
    out.push({key:"smax",label:`최대 기울기 s_max (${ch.id})`,unit:"",v:sl[k]});
    const n2=sweep(M,env,[2*Ith],1,pr.nm).n[0]; M.lv.forEach((l,i)=>{ if(!l.ground&&n2[i]>1e-12) out.push({key:"pop"+i,label:`population ${l.sp} ${l.id} @ 2·I_th`,unit:"",v:n2[i]}); }); }
  return out.map(q=>Object.assign(q,{Iref}));
}
function derivedIntervals(B,best,prof,thr,opt){
  const b=derived(B,best.theta,opt); const ens=[];
  (prof||[]).forEach(pf=>{ if(!pf) return; pf.pts.forEach(p=>{ if(p.chi2-best.chi2<=thr&&p.theta) ens.push(p.theta); }); });
  const vals=ens.map(t=>{ try{ return derived(B,t,opt); }catch(e){ return null; } }).filter(Boolean);
  return b.map(q=>{ const vs=vals.map(v=>{ const o=v.find(z=>z.key===q.key); return o?o.v:null; }).filter(v=>v>0&&isFinite(v)).concat(q.v>0?[q.v]:[]);
    return Object.assign({},q,{lo:vs.length?Math.min(...vs):null,hi:vs.length?Math.max(...vs):null,n:vs.length}); });
}

const API={compileExpr,compile:m=>prep(compile(m)),coeffs,molFrac,steady,sweep,sensitivity,integrate,riseTime,decayTime,observe,observeMulti,
  buildProblem,residuals,lm,multistart,profile,fim,eigSym,aicc,dofInfo,slopes,logGrid,photonEnergyJ,templateLee,templateTm,envFor,chi2,thetaFromParams,mulberry,lhs,channelValue,adequacy,insight,slice2D,procRates,derived,derivedIntervals};
root.DRE_CORE=API;
if(typeof module!=="undefined") module.exports=API;
})(typeof self!=="undefined"?self:(typeof window!=="undefined"?window:globalThis));
