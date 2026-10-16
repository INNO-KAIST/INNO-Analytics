/* askel_cr.js — skeleton graph analysis (2D/3D) compatible with Fiji's "Analyze Skeleton (2D/3D)".
 *
 * Independent clean-room implementation written only from a behavioural specification (cleanroom/SPEC_analyze_skeleton.md)
 * and black-box testing of the compiled plugin; no Fiji/AnalyzeSkeleton source code was consulted. Evidence for every rule
 * that was not in the specification is logged in cleanroom/NOTES.md.
 * Method: I. Arganda-Carreras, R. Fernández-González, A. Muñoz-Barrutia, C. Ortiz-De-Solorzano,
 *   "3D reconstruction of histological sections: Application to mammary gland tissue", Microsc. Res. Tech. 73:1019–1029 (2010).
 *
 * Design: foreground voxels are given compact ids in scan order (z, x, y); a padded id map gives O(1) neighbour access.
 * Labelling (components, junction clusters) → explicit per-tree graph construction (walks) → per-tree metrics →
 * optional end / cycle pruning (cycle cuts trigger a re-analysis of the cut skeleton) → longest shortest path: the exact
 * matrix semantics for graphs up to FW_MAX vertices; above that the same quantities without the n^3 matrix
 * (block-cut-tree eccentricities, or eccentricity bounding for big biconnected blocks, + Dijkstra), with a `warn` entry
 * when ties could make the start point / path depend on Fiji's relaxation order.
 * API: analyzeSkeleton(vol,cal,opt) → {trees, edges, tag, lspImage, warn}  (see SPEC_analyze_skeleton.md §1, §10).
 */
(function(root){ "use strict";
const END=30, JUNC=70, SLAB=127, SPATH=96;
const C_END=1, C_JUNC=2, C_SLAB=3;
let FW_MAX=600; let BLOCK_BUDGET=4e6; // exact matrix procedure up to this many vertices

// ------------------------------------------------------------------ volume context
function makeCtx(vol,cal){
  const w=vol.w|0, h=vol.h|0, d=(vol.d|0)||1, N=w*h*d, src=vol.data;
  const pw=cal&&cal.pw>0?+cal.pw:1, ph=cal&&cal.ph>0?+cal.ph:1, pd=cal&&cal.pd>0?+cal.pd:1;
  const is3=d>1, W=w+2, H=h+2, D=is3?d+2:1, pz=is3?1:0, sz=W*H;
  // compact ids in scan order: z outer, x, y inner
  let n=0; for(let i=0;i<N;i++) if(src[i]) n++;
  const X=new Int32Array(n), Y=new Int32Array(n), Z=new Int32Array(n), orig=new Int32Array(n);
  const idm=new Int32Array(W*H*D).fill(-1);
  let k=0;
  for(let z=0;z<d;z++) for(let x=0;x<w;x++){ let i=x+w*h*z, p=(x+1)+W*(1+H*(z+pz));
    for(let y=0;y<h;y++,i+=w,p+=W) if(src[i]){ X[k]=x; Y[k]=y; Z[k]=z; orig[k]=i; idm[p]=k; k++; } }
  // neighbour offsets in pick priority order: dx = 1,0,-1 ; dy = 1,0,-1 ; dz = -1,0,1
  const offs=[]; for(const dx of [1,0,-1]) for(const dy of [1,0,-1]) for(const dz of (is3?[-1,0,1]:[0])){ if(!dx&&!dy&&!dz) continue; offs.push(dx+W*dy+sz*dz); }
  const OFF=Int32Array.from(offs), NO=OFF.length;
  const pad=new Int32Array(n); for(let i=0;i<n;i++) pad[i]=(X[i]+1)+W*((Y[i]+1)+H*(Z[i]+pz));
  return {w,h,d,N,n,pw,ph,pd,X,Y,Z,orig,idm,OFF,NO,pad};
}
function dist(c,a,b){ const dx=(c.X[a]-c.X[b])*c.pw, dy=(c.Y[a]-c.Y[b])*c.ph, dz=(c.Z[a]-c.Z[b])*c.pd; return Math.sqrt(dx*dx+dy*dy+dz*dz); }

// classify: active[i] (1 = foreground) → cls (C_END/C_JUNC/C_SLAB, 0 for inactive)
function classify(c,active){ const n=c.n, cls=new Uint8Array(n), idm=c.idm, OFF=c.OFF, NO=c.NO, pad=c.pad;
  for(let i=0;i<n;i++){ if(!active[i]) continue; let m=0; const p=pad[i];
    for(let o=0;o<NO;o++){ const q=idm[p+OFF[o]]; if(q>=0&&active[q]) m++; }
    cls[i]=m<2?C_END:m>2?C_JUNC:C_SLAB; }
  return cls; }

// ------------------------------------------------------------------ graph construction (one pass over all trees)
function buildGraph(c,cls){
  const n=c.n, idm=c.idm, OFF=c.OFF, NO=c.NO, pad=c.pad;
  // --- components (26-connected over active voxels) and tree numbering
  const comp=new Int32Array(n).fill(-1), queue=new Int32Array(Math.max(1,n)); let nc=0;
  for(let s=0;s<n;s++){ if(!cls[s]||comp[s]>=0) continue; let qh=0, qt=0; queue[qt++]=s; comp[s]=nc;
    while(qh<qt){ const v=queue[qh++], p=pad[v]; for(let o=0;o<NO;o++){ const q=idm[p+OFF[o]]; if(q>=0&&cls[q]&&comp[q]<0){ comp[q]=nc; queue[qt++]=q; } } }
    nc++; }
  const treeOfComp=new Int32Array(nc).fill(-1); let nt=0; const cycleStart=[];
  for(const want of [C_END,C_JUNC,C_SLAB]) for(let i=0;i<n;i++) if(cls[i]===want&&treeOfComp[comp[i]]<0){ treeOfComp[comp[i]]=nt++; cycleStart.push(want===C_SLAB?i:-1); }
  const tree=new Int32Array(n); for(let i=0;i<n;i++) tree[i]=cls[i]?treeOfComp[comp[i]]:-1;
  // --- junction clusters (exploration order as specified)
  const clusterOf=new Int32Array(n).fill(-1); const cl=[]; // cl[k] = {tree, vox:Int32Array}
  { const mark=new Uint8Array(n); const buf=new Int32Array(Math.max(1,n));
    const pickJ=(v)=>{ const p=pad[v]; for(let o=0;o<NO;o++){ const q=idm[p+OFF[o]]; if(q>=0&&cls[q]===C_JUNC&&!mark[q]) return q; } return -1; };
    for(let s=0;s<n;s++){ if(cls[s]!==C_JUNC||mark[s]) continue; let len=0, head=0; buf[len++]=s; mark[s]=1; let nx=pickJ(s);
      while(nx>=0||head<len){ if(nx>=0){ buf[len++]=nx; mark[nx]=1; nx=pickJ(nx); } else { nx=pickJ(buf[head]); if(nx<0) head++; } }
      const k=cl.length, vox=buf.slice(0,len); for(let i=0;i<len;i++) clusterOf[vox[i]]=k; cl.push({tree:tree[s],vox}); } }
  // --- per tree lists
  const T=[]; for(let t=0;t<nt;t++) T.push({ends:[],clusters:[],endPoints:0,junctionVoxels:0,slabVox:0,cycle:cycleStart[t]});
  for(let i=0;i<n;i++){ if(!cls[i]) continue; const t=T[tree[i]]; if(cls[i]===C_END){ t.ends.push(i); t.endPoints++; } else if(cls[i]===C_JUNC) t.junctionVoxels++; else t.slabVox++; }
  for(let k=0;k<cl.length;k++) T[cl[k].tree].clusters.push(k);
  // --- walks
  const visited=new Uint8Array(n); const slabBuf=new Int32Array(Math.max(1,n)); let sb=0;
  const edges=[]; // {tree,v1,v2,len,ra,s0,s1}  (v1/v2 = local vertex index)
  const pickU=(v)=>{ const p=pad[v]; for(let o=0;o<NO;o++){ const q=idm[p+OFF[o]]; if(q>=0&&cls[q]&&!visited[q]) return q; } return -1; };
  const pickVJ=(v,exCl)=>{ const p=pad[v]; for(let o=0;o<NO;o++){ const q=idm[p+OFF[o]]; if(q>=0&&cls[q]===C_JUNC&&visited[q]&&clusterOf[q]!==exCl) return q; } return -1; };
  for(let t=0;t<nt;t++){ const TT=T[t];
    // vertex list: clusters first (local index = position), then added vertices; vFirst = first voxel, vCl = cluster id or -1
    const vFirst=[], vCl=[]; const clLocal=new Map();
    for(const k of TT.clusters){ clLocal.set(k,vFirst.length); vFirst.push(cl[k].vox[0]); vCl.push(k); }
    const newV=(vox)=>{ vFirst.push(vox); vCl.push(-1); return vFirst.length-1; };
    let slabCount=0, branches=0, total=0, maxLen=0;
    // walk state
    let wL=0, wR=0, wTerm=-1, wFar=-1; const hist=new Int32Array(6); let hk=0;
    const push=(v)=>{ hk++; hist[hk%6]=v; const m=hk<5?hk:5; wR+=dist(c,v,hist[(hk-m)%6])/m; };
    const walk=(s)=>{ visited[s]=1; hk=0; hist[0]=s; wL=0; wR=0; wTerm=-1; wFar=-1; let prev=s, nx=pickU(s); if(nx<0) return;
      while(nx>=0&&cls[nx]===C_SLAB){ slabCount++; slabBuf[sb++]=nx; wL+=dist(c,prev,nx); push(nx); visited[nx]=1; prev=nx; nx=pickU(prev); }
      if(nx>=0){ wL+=dist(c,prev,nx); push(nx); visited[nx]=1;
        if(cls[nx]===C_END) wFar=newV(nx); else { const k=clusterOf[nx]; wFar=clLocal.get(k); const dd=dist(c,cl[k].vox[0],nx); wL+=dd; wR+=dd; }
        wTerm=nx; }
      else wTerm=prev; };
    TT.edgeIdx=[]; const addEdge=(v1,v2,len,ra,s0,s1)=>{ TT.edgeIdx.push(edges.length); edges.push({tree:t,v1,v2,len,ra,s0,s1,removed:false}); };
    // phase A
    for(const e of TT.ends){ if(visited[e]) continue; const A=newV(e); const s0=sb; walk(e);
      if(wL===0){ const q=pickVJ(e,-1); if(q>=0){ const k=clusterOf[q]; const len=dist(c,e,q)+dist(c,cl[k].vox[0],e); addEdge(A,clLocal.get(k),len,0,s0,s0); branches++; total+=len; } continue; }
      let L=wL, far=wFar;
      if(cls[wTerm]===C_SLAB){ const q=pickVJ(wTerm,-1); if(q>=0){ const k=clusterOf[q]; far=clLocal.get(k); L+=dist(c,q,wTerm)+dist(c,cl[k].vox[0],q); } else { far=A; L+=dist(c,e,wTerm); } }
      addEdge(A,far,L,wR,s0,sb); branches++; total+=L; if(L>maxLen) maxLen=L; }
    // phase B
    for(const k of TT.clusters){ const V=clLocal.get(k), vox=cl[k].vox, f0=vox[0];
      for(let ii=0;ii<vox.length;ii++){ const j=vox[ii]; visited[j]=1; let nb;
        while((nb=pickU(j))>=0){
          if(cls[nb]===C_JUNC){ visited[nb]=1; continue; }
          const s0=sb; slabBuf[sb++]=nb; slabCount++; let L=dist(c,j,nb); walk(nb); L+=wL; const R=wR; let term=wTerm, far=wFar;
          if(term<0) term=nb;
          if(cls[term]===C_SLAB){ const q=pickVJ(term,k); if(q>=0){ far=clLocal.get(clusterOf[q]); L+=dist(c,q,term); } else far=V; }
          L+=dist(c,f0,j); addEdge(V,far,L,R,s0,sb); branches++; total+=L; if(L>maxLen) maxLen=L; } } }
    // phase C
    if(TT.cycle>=0){ const s=TT.cycle; const C=newV(s); const s0=sb; slabBuf[sb++]=s; slabCount++; walk(s); const L=wL;
      if(L!==0){ branches++; total+=L; if(L>maxLen) maxLen=L; } addEdge(C,C,L,wR,s0,sb); }
    TT.vFirst=vFirst; TT.vCl=vCl; TT.slabs=slabCount; TT.branches=branches; TT.total=total; TT.maxLength=maxLen;
  }
  return {nt,T,cl,clusterOf,tree,edges,slabBuf,comp};
}

// ------------------------------------------------------------------ longest shortest path
// g: {nv, E:[{a,b,len,idx}]}  returns {best, a, b, pred(function or array for a), warn}
function lspMatrix(nv,E){
  const n=nv, A=new Float64Array(n*n).fill(Infinity), P=new Int32Array(n*n).fill(-1);
  for(const e of E){ const r=e.a, c=e.b; A[r*n+r]=0; A[c*n+c]=0; A[r*n+c]=e.len; A[c*n+r]=e.len; P[r*n+r]=-1; P[c*n+c]=-1; P[r*n+c]=r; P[c*n+r]=c; }
  for(let k=0;k<n;k++){ const kn=k*n; for(let i=0;i<n;i++){ const in_=i*n, aik=A[in_+k]; if(aik===Infinity) continue;
      for(let j=0;j<n;j++){ const s=aik+A[kn+j]; if(s<A[in_+j]){ A[in_+j]=s; P[in_+j]=P[kn+j]; } } } }
  let best=0, a=-1, b=-1; for(let i=0;i<n;i++) for(let j=0;j<n;j++){ const v=A[i*n+j]; if(v!==Infinity&&v>best){ best=v; a=i; b=j; } }
  return {best,a,b,pred:(cur)=>P[a*n+cur],warn:null};
}

// large graphs: same quantities without the n^3 matrix (see NOTES.md)
function lspLarge(nv,E){
  // simple graph with "last edge wins" weights; diagonal values
  const diagInit=new Float64Array(nv).fill(Infinity); const wmap=new Map();
  for(const e of E){ if(e.a===e.b){ diagInit[e.a]=e.len; } else { diagInit[e.a]=0; diagInit[e.b]=0; const lo=Math.min(e.a,e.b), hi=Math.max(e.a,e.b); wmap.set(lo*nv+hi,e.len); } }
  const deg=new Int32Array(nv+1); for(const key of wmap.keys()){ const lo=Math.floor(key/nv), hi=key-lo*nv; deg[lo]++; deg[hi]++; }
  const st=new Int32Array(nv+1); for(let i=0;i<nv;i++) st[i+1]=st[i]+deg[i];
  const adj=new Int32Array(st[nv]), aw=new Float64Array(st[nv]), fill=st.slice(0,nv);
  for(const [key,len] of wmap){ const lo=Math.floor(key/nv), hi=key-lo*nv; adj[fill[lo]]=hi; aw[fill[lo]++]=len; adj[fill[hi]]=lo; aw[fill[hi]++]=len; }
  const diag=new Float64Array(nv); for(let v=0;v<nv;v++){ let m=Infinity; for(let p=st[v];p<st[v+1];p++) if(aw[p]<m) m=aw[p]; diag[v]=Math.min(diagInit[v],2*m); }
  const dg=(v)=>diag[v]===Infinity?0:diag[v];
  let M=0, a=-1, near=0;
  const ecc=eccentricities(nv,st,adj,aw);
  if(ecc){ for(let v=0;v<nv;v++){ if(ecc[v]>M) M=ecc[v]; if(dg(v)>M) M=dg(v); }
    if(!(M>0)) return {best:0,a:-1,b:-1,pred:null,warn:null};
    const tol=1e-12*M; for(let v=0;v<nv;v++){ if(Math.max(ecc[v],dg(v))>=M-tol){ near++; if(a<0) a=v; } } }
  else { const r=boundingEcc(nv,st,adj,aw,dg); M=r.M; a=r.a; near=r.near; if(!(M>0)) return {best:0,a:-1,b:-1,pred:null,warn:null}; }
  const tol=1e-12*M;
  const sp=dijkstra(nv,st,adj,aw,a);
  let b=-1; for(let j=0;j<nv;j++){ const v=j===a?(diag[a]===Infinity?0:diag[a]):sp.dist[j]; if(v!==Infinity&&v>=M-tol){ b=j; break; } }
  let warn=null; if(near>2) warn="Longest shortest path: 길이가 (거의) 같은 꼭짓점 쌍이 여럿 — 시작점이 Fiji 와 다를 수 있음 (길이는 같음)";
  if(b!==a){ let c=b; while(c!==a){ if(sp.tie[c]){ warn=warn||"Longest shortest path: 길이가 같은 다른 경로가 있음 — 경로 화소가 Fiji 와 다를 수 있음 (길이는 같음)"; break; } c=sp.pred[c]; } }
  return {best:b===a?(diag[a]):sp.dist[b],a,b,pred:(cur)=>sp.pred[cur],warn};
}

// exact eccentricity bounds (Takes & Kosters style) for graphs with large biconnected blocks:
// returns M = max over vertices of max(ecc, diag), a = smallest vertex index reaching M (within 1e-12 rel.), near = #vertices possibly reaching M
function boundingEcc(nv,st,adj,aw,dg){
  const eL=new Float64Array(nv), eU=new Float64Array(nv).fill(Infinity), ex=new Uint8Array(nv);
  let DL=0; for(let v=0;v<nv;v++){ const x=dg(v); if(x>DL) DL=x; if(st[v]===st[v+1]){ eU[v]=0; ex[v]=1; } }
  const slack=(x)=>x*1e-9; // rounding slack for the bound arithmetic
  const sweep=(v)=>{ const sp=dijkstra(nv,st,adj,aw,v); let e=0; for(let w=0;w<nv;w++){ const dw=sp.dist[w]; if(dw!==Infinity&&dw>e) e=dw; }
    for(let w=0;w<nv;w++){ const dw=sp.dist[w]; if(dw===Infinity) continue; const lo=Math.max(dw,e-dw), hi=e+dw; if(lo>eL[w]) eL[w]=lo; if(hi<eU[w]) eU[w]=hi; }
    eL[v]=eU[v]=e; ex[v]=1; if(e>DL) DL=e; return e; };
  // phase 1: settle the maximum
  let turn=0;
  for(;;){ let best=-1, bv=-1, low=Infinity, lv=-1, open=0; const lim=DL+slack(DL);
    for(let v=0;v<nv;v++){ if(ex[v]||eU[v]<=lim) continue; open++; if(eU[v]>best||(eU[v]===best&&eL[v]>eL[bv])){ best=eU[v]; bv=v; } if(eL[v]<low){ low=eL[v]; lv=v; } }
    if(!open) break; sweep((turn++&1)&&lv>=0?lv:bv); }
  const M=DL; if(!(M>0)) return {M:0,a:-1,near:0};
  const tol=1e-12*M, cut=M-tol-slack(M);
  const reach=(v)=>Math.max(dg(v),ex[v]?eL[v]:-1);
  let a=-1; for(let v=0;v<nv&&a<0;v++){ if(dg(v)>=M-tol){ a=v; break; } if(eU[v]<cut) continue; if(!ex[v]) sweep(v); if(reach(v)>=M-tol) a=v; }
  // count other vertices that may reach the maximum (resolve a limited number; unresolved ones count as possible ties)
  let near=0, budget=48; for(let v=0;v<nv;v++){ if(dg(v)>=M-tol){ near++; continue; } if(eU[v]<cut) continue; if(!ex[v]&&budget>0){ budget--; sweep(v); } if(!ex[v]||reach(v)>=M-tol) near++; }
  return {M,a,near};
}
function dijkstra(nv,st,adj,aw,s){
  const dist=new Float64Array(nv).fill(Infinity), pred=new Int32Array(nv).fill(-1), tie=new Uint8Array(nv), done=new Uint8Array(nv);
  const hv=new Int32Array(Math.max(16,st[nv]+1)), hd=new Float64Array(hv.length); let hn=0;
  const hpush=(v,dv)=>{ let i=hn++; while(i>0){ const p=(i-1)>>1; if(hd[p]<=dv) break; hv[i]=hv[p]; hd[i]=hd[p]; i=p; } hv[i]=v; hd[i]=dv; };
  const hpop=()=>{ const v=hv[0], lv=hv[--hn], ld=hd[hn]; let i=0; for(;;){ let l=2*i+1; if(l>=hn) break; if(l+1<hn&&hd[l+1]<hd[l]) l++; if(hd[l]>=ld) break; hv[i]=hv[l]; hd[i]=hd[l]; i=l; } hv[i]=lv; hd[i]=ld; return v; };
  dist[s]=0; hpush(s,0);
  while(hn>0){ const v=hpop(); if(done[v]) continue; done[v]=1; const dv=dist[v];
    for(let p=st[v];p<st[v+1];p++){ const u=adj[p]; if(done[u]) continue; const nd=dv+aw[p];
      if(nd<dist[u]){ const old=dist[u]; dist[u]=nd; pred[u]=v; tie[u]=(old-nd<=1e-12*nd)?1:0; hpush(u,nd); }
      else if(nd-dist[u]<=1e-12*nd) tie[u]=1; } }
  return {dist,pred,tie};
}
// eccentricity (max finite shortest-path distance to another vertex) of every vertex: block-cut tree rerooting
function eccentricities(nv,st,adj,aw){
  const ecc=new Float64Array(nv);
  // biconnected components (iterative Tarjan on edges)
  const disc=new Int32Array(nv).fill(-1), low=new Int32Array(nv), it=new Int32Array(nv), par=new Int32Array(nv).fill(-1);
  const blocks=[]; const estack=[]; let time=0; const stack=new Int32Array(nv);
  for(let r=0;r<nv;r++){ if(disc[r]>=0) continue; if(st[r]===st[r+1]){ disc[r]=time++; blocks.push([r]); continue; }
    let sp=0; stack[sp++]=r; disc[r]=low[r]=time++; it[r]=st[r];
    while(sp>0){ const v=stack[sp-1];
      if(it[v]<st[v+1]){ const p=it[v]++, u=adj[p];
        if(disc[u]<0){ par[u]=v; disc[u]=low[u]=time++; it[u]=st[u]; estack.push(v,u); stack[sp++]=u; }
        else if(u!==par[v]&&disc[u]<disc[v]){ estack.push(v,u); if(disc[u]<low[v]) low[v]=disc[u]; } }
      else { sp--; const p=par[v]; if(p>=0){ if(low[v]<low[p]) low[p]=low[v];
          if(low[v]>=disc[p]){ const set=new Set(); for(;;){ const b=estack.pop(), a=estack.pop(); set.add(a); set.add(b); if(a===p&&b===v) break; } blocks.push(Array.from(set)); } } } } }
  // per-block distance matrices (Dijkstra restricted to the block); too costly for big blocks → caller uses bounding
  { let cost=0; for(const b of blocks) cost+=b.length*b.length; if(cost>BLOCK_BUDGET) return null; }
  const nb=blocks.length; const bD=new Array(nb); const vb=new Array(nv); for(let v=0;v<nv;v++) vb[v]=[];
  const loc=new Int32Array(nv).fill(-1);
  for(let B=0;B<nb;B++){ const vs=blocks[B], m=vs.length; for(let i=0;i<m;i++){ loc[vs[i]]=i; vb[vs[i]].push(B); }
    const D=new Float64Array(m*m);
    if(m===2){ let wv=Infinity; const a=vs[0], b=vs[1]; for(let p=st[a];p<st[a+1];p++) if(adj[p]===b) wv=aw[p]; D[1]=D[2]=wv; }
    else if(m>2) blockAPSP(vs,loc,st,adj,aw,D);
    bD[B]=D; for(let i=0;i<m;i++) loc[vs[i]]=-1; }
  // block-cut tree rooted per connected component; order blocks by BFS from a root block
  const bparV=new Int32Array(nb).fill(-1), seen=new Uint8Array(nb), order=[];
  const vpar=new Int32Array(nv).fill(-1); // parent block of vertex v (block in which v is the "child side")
  for(let R=0;R<nb;R++){ if(seen[R]) continue; seen[R]=1; order.push(R); for(let qi=order.length-1;qi<order.length;qi++){ const B=order[qi];
      for(const v of blocks[B]){ if(v===bparV[B]) continue; vpar[v]=B; for(const B2 of vb[v]) if(!seen[B2]){ seen[B2]=1; bparV[B2]=v; order.push(B2); } } } }
  // down[v] = max distance from v into the sub-blocks hanging below v (excluding vpar[v]); downFar[B] = from bparV[B] into B and below
  const down=new Float64Array(nv), downFar=new Float64Array(nb).fill(0);
  for(let qi=order.length-1;qi>=0;qi--){ const B=order[qi], vs=blocks[B], m=vs.length, D=bD[B], p=bparV[B]; if(p<0) continue;
    const ip=vs.indexOf(p); let mx=0; for(let i=0;i<m;i++){ if(i===ip) continue; const val=D[ip*m+i]+down[vs[i]]; if(val!==Infinity&&val>mx) mx=val; }
    downFar[B]=mx; if(mx>down[p]) down[p]=mx; }
  // up pass
  const best1=new Float64Array(nv), best2=new Float64Array(nv), arg1=new Int32Array(nv).fill(-1);
  for(let B=0;B<nb;B++){ const p=bparV[B]; if(p<0) continue; const v=downFar[B]; if(v>best1[p]){ best2[p]=best1[p]; best1[p]=v; arg1[p]=B; } else if(v>best2[p]) best2[p]=v; }
  const upV=new Float64Array(nv); // upV[v] = max distance from v through its parent block vpar[v] (and everything above)
  for(const B of order){ const vs=blocks[B], m=vs.length, D=bD[B], p=bparV[B];
    // h[i] = height hanging at vs[i] excluding B
    const hgt=new Float64Array(m); for(let i=0;i<m;i++){ const v=vs[i]; if(v===p){ const other=arg1[p]===B?best2[p]:best1[p]; hgt[i]=Math.max(upV[p],other); } else hgt[i]=down[v]; }
    for(let i=0;i<m;i++){ const v=vs[i]; if(v===p) continue; let mx=0; for(let j=0;j<m;j++){ if(j===i) continue; const val=D[i*m+j]+hgt[j]; if(val!==Infinity&&val>mx) mx=val; } upV[v]=mx; } }
  for(let v=0;v<nv;v++) ecc[v]=Math.max(down[v],upV[v]);
  return ecc;
}
function blockAPSP(vs,loc,st,adj,aw,D){ const m=vs.length; const dist=new Float64Array(m), done=new Uint8Array(m);
  let ne=0; for(let i=0;i<m;i++){ const g=vs[i]; for(let p=st[g];p<st[g+1];p++) if(loc[adj[p]]>=0) ne++; }
  const hv=new Int32Array(ne+16), hd=new Float64Array(ne+16); let hn=0; // lazy heap: at most one push per directed edge + source
  const hpush=(v,dv)=>{ let i=hn++; while(i>0){ const p=(i-1)>>1; if(hd[p]<=dv) break; hv[i]=hv[p]; hd[i]=hd[p]; i=p; } hv[i]=v; hd[i]=dv; };
  const hpop=()=>{ const v=hv[0], lv=hv[--hn], ld=hd[hn]; let i=0; for(;;){ let l=2*i+1; if(l>=hn) break; if(l+1<hn&&hd[l+1]<hd[l]) l++; if(hd[l]>=ld) break; hv[i]=hv[l]; hd[i]=hd[l]; i=l; } hv[i]=lv; hd[i]=ld; return v; };
  for(let s=0;s<m;s++){ dist.fill(Infinity); done.fill(0); hn=0; dist[s]=0; hpush(s,0);
    while(hn>0){ const v=hpop(); if(done[v]) continue; done[v]=1; const g=vs[v];
      for(let p=st[g];p<st[g+1];p++){ const u=loc[adj[p]]; if(u<0||done[u]) continue; const nd=dist[v]+aw[p]; if(nd<dist[u]){ dist[u]=nd; hpush(u,nd); } } }
    D.set(dist,s*m); } }

// ------------------------------------------------------------------ pruning (behaviour derived by black-box experiments, see NOTES.md)
const PRUNE_NAMES=["none","shortest branch","lowest intensity voxel","lowest intensity branch"];
function pruneMode(p){ if(typeof p==="number") return p>=0&&p<=3?p|0:0; const k=PRUNE_NAMES.indexOf(String(p||"none").toLowerCase()); return k<0?0:k; }
// end pruning on the first analysis: for each end-point vertex (vertex order) that still has a branch, remove that branch
// (its slab voxels) and the end-point voxel; the vertex and edge are dropped from the graph
function pruneEndsStep(c,cls,g,rem){
  for(let t=0;t<g.nt;t++){ const T=g.T[t], nv=T.vFirst.length; const inc=new Array(nv); for(let v=0;v<nv;v++) inc[v]=[];
    T.edgeIdx.forEach((k)=>{ const e=g.edges[k]; inc[e.v1].push(k); if(e.v2!==e.v1) inc[e.v2].push(k); });
    T.vRemoved=new Uint8Array(nv);
    for(let v=0;v<nv;v++){ if(T.vCl[v]>=0||cls[T.vFirst[v]]!==C_END) continue; const live=inc[v].filter(k=>!g.edges[k].removed); if(!live.length) continue;
      for(const k of live){ const e=g.edges[k]; e.removed=true; for(let s=e.s0;s<e.s1;s++) rem[g.slabBuf[s]]=1; }
      T.vRemoved[v]=1; rem[T.vFirst[v]]=1; } } }
// window mean of the original image for cycle-voxel intensities: box of half-size round(R/p) per axis, R = max voxel size,
// out-of-bounds counted as 0 (the divisor is always the full box size)
function makeIntensity(c,orig){ const R=Math.max(c.pw,c.ph,c.pd), rx=Math.floor(R/c.pw+0.5), ry=Math.floor(R/c.ph+0.5), rz=c.d>1?Math.floor(R/c.pd+0.5):0;
  const cnt=(2*rx+1)*(2*ry+1)*(2*rz+1), w=c.w, h=c.h, d=c.d, cache=new Map();
  return (v)=>{ let m=cache.get(v); if(m!==undefined) return m; const x0=c.X[v], y0=c.Y[v], z0=c.Z[v]; let sum=0;
    for(let z=Math.max(0,z0-rz);z<=Math.min(d-1,z0+rz);z++) for(let y=Math.max(0,y0-ry);y<=Math.min(h-1,y0+ry);y++){ const row=w*(y+h*z); for(let x=Math.max(0,x0-rx);x<=Math.min(w-1,x0+rx);x++) sum+=orig[x+row]; }
    m=sum/cnt; cache.set(v,m); return m; }; }
// cycle detection: depth-first search from the tree's root (first non-cluster vertex, else the first cluster) with an explicit
// stack of (vertex, via edge); an edge popped towards an already visited vertex is a back edge; its cycle = back edge + tree path
function cycleCuts(c,cls,g,mode,orig){
  const cuts=[]; const I=mode>=2?makeIntensity(c,orig):null;
  const minC=new Map(), avgC=new Map(); // per-edge darkest slab (first of equals) and mean intensity over vertex voxels + slabs
  const edgeMin=(q)=>{ let r=minC.get(q); if(r) return r; const e=g.edges[q]; let bv=-1, bi=Infinity; for(let s=e.s0;s<e.s1;s++){ const v=g.slabBuf[s], i=I(v); if(i<bi){ bi=i; bv=v; } } r={v:bv,i:bi}; minC.set(q,r); return r; };
  const edgeAvg=(q,T)=>{ let r=avgC.get(q); if(r!==undefined) return r; const e=g.edges[q]; let sum=0, m=0;
    const vox=(vv)=>{ const ci=T.vCl[vv]; if(ci>=0){ for(const p of g.cl[ci].vox){ sum+=I(p); m++; } } else { sum+=I(T.vFirst[vv]); m++; } };
    vox(e.v1); for(let s=e.s0;s<e.s1;s++){ sum+=I(g.slabBuf[s]); m++; } if(e.v2!==e.v1) vox(e.v2); r=sum/m; avgC.set(q,r); return r; };
  for(let t=0;t<g.nt;t++){ const T=g.T[t], nv=T.vFirst.length; const E=T.edgeIdx.filter(k=>!g.edges[k].removed); if(!E.length) continue;
    const vr=T.vRemoved; let root=-1; for(let v=0;v<nv;v++) if(T.vCl[v]<0){ root=v; break; } if(root<0) root=0; if(vr&&vr[root]) continue;
    const adj=new Array(nv); for(let v=0;v<nv;v++) adj[v]=[]; for(const k of E){ const e=g.edges[k]; adj[e.v1].push(k); if(e.v2!==e.v1) adj[e.v2].push(k); }
    const vis=new Uint8Array(nv), par=new Int32Array(nv).fill(-1), pe=new Int32Array(nv).fill(-1), lab=new Map(); const backs=[];
    const st=[[root,-1,-1]];
    while(st.length){ const [v,k,from]=st.pop();
      if(k>=0){ if(lab.has(k)) continue; if(vis[v]){ lab.set(k,2); backs.push([from,v,k]); continue; } lab.set(k,1); par[v]=from; pe[v]=k; } else if(vis[v]) continue;
      vis[v]=1; for(const q of adj[v]){ if(lab.has(q)) continue; const e=g.edges[q]; st.push([e.v1===v?e.v2:e.v1,q,v]); } }
    for(const [from,u,k] of backs){ const cyc=[k]; let x=from; while(x!==u&&x>=0){ cyc.push(pe[x]); x=par[x]; } if(x!==u) continue;
      const slabs=(q)=>{ const e=g.edges[q]; return g.slabBuf.subarray(e.s0,e.s1); };
      if(T.cycle>=0){ cuts.push(T.cycle); continue; }
      if(mode===1){ let best=-1, bn=0; for(const q of cyc){ const m=g.edges[q].s1-g.edges[q].s0; if(best<0||m<bn){ best=q; bn=m; } } const sl=slabs(best); if(sl.length) cuts.push(sl[sl.length>>1]); }
      else if(mode===2){ let bv=-1, bi=Infinity; for(const q of cyc){ const m=edgeMin(q); if(m.v>=0&&m.i<bi){ bi=m.i; bv=m.v; } } if(bv>=0) cuts.push(bv); }
      else { let best=-1, ba=Infinity; for(const q of cyc){ const avg=edgeAvg(q,T); if(avg<ba){ ba=avg; best=q; } }
        if(best>=0){ const m=edgeMin(best); if(m.v>=0) cuts.push(m.v); } } } }
  return cuts; }

// ------------------------------------------------------------------ main
function analyzeSkeleton(vol,cal,opt){
  opt=opt||{}; const IC0=root&&root.IMG_CORE; const prog=typeof opt.progress==="function"?opt.progress:(IC0&&typeof IC0.progress==="function"?(f=>IC0.progress(f)):null); const warn=[];
  { const NV=(vol.w|0)*(vol.h|0)*((vol.d|0)||1); if(NV>2.5e8) throw new Error(`Analyze Skeleton: 부피가 너무 큽니다 (${(NV/1e6).toFixed(0)} M 화소 > 250 M) — 잘라서 분석하세요`); }
  if(!vol||!(vol.w>0)||!(vol.h>0)||!vol.data) throw new Error("analyzeSkeleton: invalid volume");
  const c=makeCtx(vol,cal); if(prog) prog(0.05);
  const active=new Uint8Array(c.n).fill(1);
  let cls=classify(c,active); if(prog) prog(0.15);
  let g=buildGraph(c,cls); if(prog) prog(0.45);
  const mode=pruneMode(opt.prune); const rem=new Uint8Array(c.n); let pruned=false;
  if(opt.pruneEnds){ pruneEndsStep(c,cls,g,rem); pruned=true; }
  if(mode>0){ let orig=opt.orig&&opt.orig.data;
    if(mode>=2&&(!orig||orig.length<c.N)){ warn.push("세기(intensity) 기준 고리 가지치기에는 같은 크기의 원본 영상이 필요합니다 — 고리는 자르지 않았습니다"); orig=null; }
    if(mode===1||orig){ const cuts=cycleCuts(c,cls,g,mode,orig);
      if(cuts.length){ for(const v of cuts) rem[v]=1; for(let i=0;i<c.n;i++) if(rem[i]) active[i]=0; // re-analyse the modified skeleton
        cls=classify(c,active); g=buildGraph(c,cls); pruned=false; if(prog) prog(0.7); } } }
  return finish(c,cls,g,opt,warn,prog,pruned?rem:null);
}
// rem (only after end pruning without a re-analysis): removed voxels / edges / vertices; counts that Fiji recomputes are recomputed
function finish(c,cls,g,opt,warn,prog,rem){
  const n=c.n, N=c.N, nt=g.nt; const alive=(v)=>cls[v]&&!(rem&&rem[v]);
  const tag=new Uint8Array(N); for(let i=0;i<n;i++) if(alive(i)) tag[c.orig[i]]=cls[i]===C_END?END:cls[i]===C_JUNC?JUNC:SLAB;
  const trees=[]; const edgesOut=[];
  const ft=(v)=>[c.X[v],c.Y[v],c.Z[v]];
  // triples / quadruples
  const tri=new Int32Array(nt), quad=new Int32Array(nt);
  for(const k of g.cl.keys()){ const C=g.cl[k]; let nb=0; for(const v of C.vox){ const p=c.pad[v]; for(let o=0;o<c.NO;o++){ const q=c.idm[p+c.OFF[o]]; if(q>=0&&(cls[q]===C_END||cls[q]===C_SLAB)&&!(rem&&rem[q])) nb++; } }
    if(nb===3) tri[C.tree]++; else if(nb===4) quad[C.tree]++; }
  let remEnd=null, remSlab=null; if(rem){ remEnd=new Int32Array(nt); remSlab=new Int32Array(nt); for(let i=0;i<n;i++) if(rem[i]&&cls[i]){ if(cls[i]===C_END) remEnd[g.tree[i]]++; else if(cls[i]===C_SLAB) remSlab[g.tree[i]]++; } }
  for(let t=0;t<nt;t++){ const T=g.T[t]; trees.push({branches:T.branches,junctions:T.clusters.length,endPoints:T.endPoints-(rem?remEnd[t]:0),junctionVoxels:T.junctionVoxels,slabs:T.slabs-(rem?remSlab[t]:0),
    avgLength:T.branches?T.total/T.branches:0,triples:tri[t],quadruples:quad[t],maxLength:T.maxLength}); }
  for(const e of g.edges){ if(e.removed) continue; const T=g.T[e.tree]; const f1=T.vFirst[e.v1], f2=T.vFirst[e.v2];
    edgesOut.push({tree:e.tree,length:e.len,v1:ft(f1),v2:ft(f2),euclid:dist(c,f1,f2),lengthRa:e.ra,nSlabs:e.s1-e.s0}); }
  let lspImage=null;
  if(opt.lsp){ lspImage=new Uint8Array(N); for(let i=0;i<n;i++) if(alive(i)) lspImage[c.orig[i]]=255;
    let sx=0, sy=0, sz=0;
    for(let t=0;t<nt;t++){ const T=g.T[t]; let nv=T.vFirst.length; let map=null, vF=T.vFirst;
      if(rem&&T.vRemoved){ map=new Int32Array(nv).fill(-1); vF=[]; for(let v=0;v<nv;v++) if(!T.vRemoved[v]){ map[v]=vF.length; vF.push(T.vFirst[v]); } nv=vF.length; }
      const E=[]; for(const k of T.edgeIdx){ const e=g.edges[k]; if(e.removed) continue; E.push({a:map?map[e.v1]:e.v1,b:map?map[e.v2]:e.v2,len:e.len,e}); }
      let best=0;
      if(nv===1){ const f=vF[0]; sx=c.X[f]*c.pw; sy=c.Y[f]*c.ph; sz=c.Z[f]*c.pd; }
      else if(nv>1){ const r=nv<=FW_MAX?lspMatrix(nv,E):lspLarge(nv,E); if(r.warn) warn.push(`골격 ${t+1}: ${r.warn}`); best=r.best;
        if(r.a>=0&&r.a!==r.b){ // path from b back to a
          const pairBest=new Map(); for(const x of E){ if(x.a===x.b) continue; const key=Math.min(x.a,x.b)*nv+Math.max(x.a,x.b); const o=pairBest.get(key); if(!o||x.len<o.len) pairBest.set(key,x); }
          let cur=r.b; while(cur!==r.a){ const pre=r.pred(cur); const x=pairBest.get(Math.min(pre,cur)*nv+Math.max(pre,cur));
            for(let s=x.e.s0;s<x.e.s1;s++) lspImage[c.orig[g.slabBuf[s]]]=SPATH; cur=pre; }
          const f=vF[r.b]; sx=c.X[f]*c.pw; sy=c.Y[f]*c.ph; sz=c.Z[f]*c.pd; } }
      Object.assign(trees[t],{lsp:best,spx:sx,spy:sy,spz:sz}); } }
  if(prog) prog(1);
  return {trees,edges:edgesOut,tag,lspImage,warn};
}

function debugGraph(vol,cal){ const c=makeCtx(vol,cal); const cls=classify(c,new Uint8Array(c.n).fill(1)); const g=buildGraph(c,cls);
  const P=(v)=>[c.X[v],c.Y[v],c.Z[v]];
  return {edges:g.edges.map(e=>({tree:e.tree,v1:e.v1,v2:e.v2,p1:P(g.T[e.tree].vFirst[e.v1]),p2:P(g.T[e.tree].vFirst[e.v2]),len:e.len,slabs:Array.from(g.slabBuf.subarray(e.s0,e.s1),P)})),
    trees:g.T.map(T=>({nv:T.vFirst.length,vFirst:T.vFirst.map(P),vCl:T.vCl,vVox:T.vFirst.map((f,i)=>T.vCl[i]>=0?Array.from(g.cl[T.vCl[i]].vox,P):[P(f)]),cycle:T.cycle>=0}))}; }
const api={analyzeSkeleton,PRUNE:true,TAGS:{END,JUNC,SLAB,SPATH},_t:{debugGraph,lspMatrix,lspLarge,setFwMax:(v)=>{FW_MAX=v;},setBlockBudget:(v)=>{BLOCK_BUDGET=v;}}};
if(typeof module!=="undefined"&&module.exports) module.exports=api;
if(root&&root.IMG_CORE){ root.IMG_CORE.analyzeSkeleton=analyzeSkeleton; root.IMG_CORE.SKELETON_TAGS=api.TAGS; root.IMG_CORE.SKELETON_PRUNE=api.PRUNE; }
})(typeof self!=="undefined"?self:globalThis);
