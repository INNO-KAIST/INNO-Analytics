// internal consistency: lspLarge (no matrix) vs lspMatrix (exact Fiji semantics) on random multigraphs
const CR=require("../askel_cr.js"); const {lspMatrix,lspLarge}=CR._t;
let seed=12345; const rnd=()=>((seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff);
let bad=0, warned=0, tot=0, warnedOK=0;
for(let it=0;it<4000;it++){ const nv=2+Math.floor(rnd()*60); const E=[]; const integer=rnd()<0.3;
  const L=()=>integer?1+Math.floor(rnd()*3):0.5+rnd()*5;
  for(let v=1;v<nv;v++) if(rnd()<0.97) E.push({a:Math.floor(rnd()*v),b:v,len:L()});
  const extra=Math.floor(rnd()*nv*(rnd()<0.5?0.1:0.6)); for(let k=0;k<extra;k++){ const a=Math.floor(rnd()*nv), b=rnd()<0.15?a:Math.floor(rnd()*nv); E.push({a,b,len:L()}); }
  // shuffle creation order
  for(let i=E.length-1;i>0;i--){ const j=Math.floor(rnd()*(i+1)); [E[i],E[j]]=[E[j],E[i]]; }
  const m=lspMatrix(nv,E), g=lspLarge(nv,E); tot++;
  const okBest=Math.abs(m.best-g.best)<=1e-12*Math.max(1,m.best);
  let pm=[], pg=[]; if(m.a>=0&&m.a!==m.b){ let c=m.b; while(c!==m.a){ pm.push(c); c=m.pred(c);} }
  if(g.a>=0&&g.a!==g.b){ let c=g.b; while(c!==g.a){ pg.push(c); c=g.pred(c);} }
  const same=okBest&&m.a===g.a&&m.b===g.b&&pm.join()===pg.join();
  if(g.warn) warned++;
  if(!same&&!g.warn){ bad++; if(bad<5) console.log("MISMATCH",nv,JSON.stringify({m:[m.best,m.a,m.b,pm],g:[g.best,g.a,g.b,pg]}),integer); }
  if(!okBest){ bad++; console.log("BEST MISMATCH",m.best,g.best); }
  if(same&&g.warn) warnedOK++;
}
console.log({tot,bad,warned,warnedOK});
