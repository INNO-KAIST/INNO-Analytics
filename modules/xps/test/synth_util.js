/* shared helpers for test_multi / test_arxps / test_image: seeded RNG, exact Poisson sampler, check/report */
"use strict";
function rng(seed){ let a=(seed>>>0)||1; const u=()=>{ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; };
  let sp=null; const g=()=>{ if(sp!=null){ const v=sp; sp=null; return v; } let u1=0; while(u1<1e-300) u1=u(); const u2=u(); const r=Math.sqrt(-2*Math.log(u1)); sp=r*Math.sin(2*Math.PI*u2); return r*Math.cos(2*Math.PI*u2); };
  const LGC=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  const lg=x=>{ x-=1; let s=LGC[0]; const t=x+7.5; for(let i=1;i<9;i++) s+=LGC[i]/(x+i); return 0.5*Math.log(2*Math.PI)+(x+0.5)*Math.log(t)-t+Math.log(s); };
  /* exact Poisson: Knuth (λ < 10), PTRS transformed rejection (Hörmann 1993, as in numpy) otherwise */
  const pois=lam=>{ if(!(lam>0)) return 0; if(lam<10){ const L=Math.exp(-lam); let k=0, p=1; do{ k++; p*=u(); }while(p>L); return k-1; }
    const sl=Math.sqrt(lam), ll=Math.log(lam), b=0.931+2.53*sl, a=-0.059+0.02483*b, ia=1.1239+1.1328/(b-3.4), vr=0.9277-3.6224/(b-2);
    for(;;){ const U=u()-0.5, V=u(), us=0.5-Math.abs(U); const k=Math.floor((2*a/us+b)*U+lam+0.43); if(us>=0.07&&V<=vr) return k; if(k<0||(us<0.013&&V>us)) continue;
      if(Math.log(V)+Math.log(ia)-Math.log(a/(us*us)+b)<=-lam+k*ll-lg(k+1)) return k; } };
  return {u,g,pois}; }
function reporter(title){
  let pass=0, fail=0; const fails=[]; console.log("### "+title);
  const ok=(cond,name,info)=>{ if(cond) pass++; else { fail++; fails.push(name+(info?" — "+info:"")); } console.log((cond?"  ok   ":"  FAIL ")+name+(info?"  ("+info+")":"")); };
  const done=()=>{ console.log(`\n${title}: ${pass} passed, ${fail} failed`); if(fail){ console.log(fails.map(f=>"  - "+f).join("\n")); process.exitCode=1; } return {pass,fail}; };
  return {ok,done};
}
const maxRel=(a,b,scale)=>{ let m=0; for(let i=0;i<a.length;i++){ const d=Math.abs(a[i]-b[i])/(scale!=null?scale:Math.max(Math.abs(b[i]),1e-300)); if(d>m) m=d; } return m; };
const F64=a=>Float64Array.from(a);
const gl=(E,c,w,m)=>{ m=m==null?0.3:m; return F64(Array.from(E,e=>{ const x=(e-c)/w; return (1-m)*Math.exp(-4*Math.LN2*x*x)+m/(1+4*x*x); })); };
const trapz=(E,y)=>{ let s=0; for(let i=1;i<E.length;i++) s+=0.5*(y[i]+y[i-1])*(E[i]-E[i-1]); return s; };
const stat=a=>{ const n=a.length, m=a.reduce((s,v)=>s+v,0)/n; return {m,sd:Math.sqrt(a.reduce((s,v)=>s+(v-m)*(v-m),0)/Math.max(n-1,1))}; };
module.exports={rng,reporter,maxRel,F64,gl,trapz,stat};
