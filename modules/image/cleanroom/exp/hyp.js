// cut-selection hypotheses for cycle pruning, evaluated against oracle removed-voxel sets
const O=require("./oracle.js"), CR=require("../askel_cr.js"), gen=require("./gen.js");
function predict(c,V){ const g=CR._t.debugGraph(c,{pw:c.cal[0],ph:c.cal[1],pd:c.cal[2]}); const rem=new Set(); const key=p=>p.join();
  const raw=(x,y,z)=>c.orig[x+c.w*(y+c.h*z)];
  const I=(p)=>{ if(!V.nb) return raw(p[0],p[1],p[2]);
    if(V.nb==="cal"){ const [pw,ph,pd]=c.cal; const R=V.R2d&&c.d===1?Math.max(pw,ph):Math.max(pw,ph,pd); const rx=Math.floor(R/pw+0.5), ry=Math.floor(R/ph+0.5), rz=c.d===1&&!V.zIn2d?0:Math.floor(R/pd+0.5);
      let s=0,n=0; for(let z=p[2]-rz;z<=p[2]+rz;z++) for(let y=p[1]-ry;y<=p[1]+ry;y++) for(let x=p[0]-rx;x<=p[0]+rx;x++){ n++; if(x<0||y<0||z<0||x>=c.w||y>=c.h||z>=c.d){ if(V.oob==="skip") n--; continue; } s+=raw(x,y,z); } return s/n; }
    if(V.oob==="shift"){ const x0=Math.max(0,Math.min(c.w-3,p[0]-1)), y0=Math.max(0,Math.min(c.h-3,p[1]-1)); const z0=c.d===1?0:Math.max(0,Math.min(c.d-3,p[2]-1)), nz=c.d===1?1:3; let s=0,n=0; for(let z=z0;z<z0+nz;z++) for(let y=y0;y<y0+3;y++) for(let x=x0;x<x0+3;x++){ if(x>=c.w||y>=c.h||z>=c.d) continue; s+=raw(x,y,z); n++; } return s/n; } let s=0, n=0; const dz=V.nb==="3x3"||c.d===1?0:1;
    for(let z=p[2]-dz;z<=p[2]+dz;z++) for(let y=p[1]-1;y<=p[1]+1;y++) for(let x=p[0]-1;x<=p[0]+1;x++){ if(x<0||y<0||z<0||x>=c.w||y>=c.h||z>=c.d){ if(V.oob==="zero") n++; else if(V.oob==="clamp"){ s+=raw(Math.min(c.w-1,Math.max(0,x)),Math.min(c.h-1,Math.max(0,y)),Math.min(c.d-1,Math.max(0,z))); n++; } else if(V.oob==="mirror"){ const m=(a,n)=>a<0?-a:a>=n?2*n-2-a:a; s+=raw(m(x,c.w),m(y,c.h),m(z,c.d)); n++; } continue; } s+=raw(x,y,z); n++; }
    return V.sumOnly?s:s/n; };
  g.trees.forEach((T,t)=>{ const E=g.edges.filter(e=>e.tree===t); const nv=T.nv; const adj=[]; for(let v=0;v<nv;v++) adj.push([]);
    E.forEach((e,k)=>{ adj[e.v1].push(k); if(e.v2!==e.v1||V.selfTwice) adj[e.v2].push(k); });
    const vis=new Uint8Array(nv), parent=new Int32Array(nv).fill(-1), pe=new Int32Array(nv).fill(-1), cls=new Int8Array(E.length); const backs=[];
    const other=(e,v)=>E[e].v1===v?E[e].v2:E[e].v1;
    const root=V.root==="first"?0:V.root(T,E);
    const visit=(v)=>{ vis[v]=1; for(const e of adj[v]){ if(cls[e]) continue; const u=other(e,v); if(!vis[u]){ cls[e]=1; parent[u]=v; pe[u]=e; visit(u); } else { cls[e]=2; backs.push([v,u,e]); } } };
    const visitS=(r)=>{ const st=[[r,-1,-1]]; while(st.length){ const [v,e,from]=st.pop();
        if(e>=0){ if(cls[e]) continue; if(vis[v]){ cls[e]=2; backs.push([from,v,e]); continue; } cls[e]=1; parent[v]=from; pe[v]=e; } else if(vis[v]) continue;
        vis[v]=1; for(const k of adj[v]){ if(cls[k]) continue; st.push([other(k,v),k,v]); } } };
    const go=V.dfs==="stackE"?visitS:visit;
    if(root>=0&&root<nv) go(root); if(V.allRoots) for(let v=0;v<nv;v++) if(!vis[v]) go(v);
    for(const [v,u,e] of backs){ const cyc=[e]; let x=v; let guard=0; while(x!==u&&x>=0&&guard++<1e5){ cyc.push(pe[x]); x=parent[x]; } if(x!==u) continue;
      const order=V.order==="rev"?cyc.slice().reverse():cyc;
      if(T.cycle&&V.cycleStartAll){ rem.add(key(E[cyc[0]].slabs[0])); continue; }
      if(V.method===1){ let best=-1; for(const k of order){ const m=V.metric==="slabs"?E[k].slabs.length:E[k].len; if(best<0||(V.le?m<=(V.metric==="slabs"?E[best].slabs.length:E[best].len):m<(V.metric==="slabs"?E[best].slabs.length:E[best].len))) best=k; }
        const e=E[best]; if(T.cycle&&V.cycleStart) rem.add(key(e.slabs[0])); else if(e.slabs.length) rem.add(key(e.slabs[Math.floor(e.slabs.length/2)])); }
      else if(V.method===2){ let bp=null, bi=Infinity; for(const k of order){ for(const p of V.vox(E[k],T)) { const i=I(p); if(i<bi||(V.le&&i<=bi)){ bi=i; bp=p; } } } if(bp) rem.add(key(bp)); }
      else if(V.method===3){ let best=-1, ba=Infinity; for(const k of order){ const vs=V.vox(E[k],T); if(!vs.length) continue; let s=0; for(const p of vs) s+=(V.avgRaw?raw(p[0],p[1],p[2]):I(p)); const a=V.sumOnly3?s:s/vs.length; if(a<ba||(V.le&&a<=ba)){ ba=a; best=k; } }
        if(best>=0){ let bp=null, bi=Infinity; for(const p of V.vox2(E[best],T)){ const i=I(p); if(i<bi){ bi=i; bp=p; } } if(bp) rem.add(key(bp)); } } } });
  return rem; }
function evaluate(P,variants,N,seed,show){ const cases=gen(seed,N); const base=cases.map(c=>Object.assign({},c,{prune:0,pruneEnds:false})); const pr=cases.map(c=>Object.assign({},c,{prune:P,pruneEnds:false}));
  const os=O.run(base.concat(pr)); const res={};
  for(const [name,V] of Object.entries(variants)){ let bad=0; const badList=[];
    pr.forEach((c,i)=>{ const o0=os[i], o1=os[N+i]; const want=new Set(); for(const [p] of o0.tag) if(!o1.tag.has(p)) want.add([p%c.w,((p/c.w)|0)%c.h,(p/(c.w*c.h))|0].join());
      const got=predict(c,V); const same=got.size===want.size&&[...got].every(k=>want.has(k)); if(!same){ bad++; badList.push([i,[...want].sort().join(" "),[...got].sort().join(" ")]); } });
    res[name]=bad; if(show) badList.slice(0,show).forEach(b=>console.log(name,"case",b[0],"want",b[1],"got",b[2])); }
  return {res,cases,os}; }
module.exports={predict,evaluate};
