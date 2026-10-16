const D=require('../src/dre_core.js'); const T=D.templateLee(); const M=D.compile(T.spec);
const base={c:8}; Object.entries(T.params).forEach(([n,p])=>base[n]=p.v);
const I=[...Array(60).keys()].map(i=>10*Math.pow(1e5,i/59));
const rnd=D.mulberry(11); let found=0;
for(let t=0;t<400&&found<3;t++){
  const L=a=>Math.pow(10,(rnd()*2-1)*a);
  const e=Object.assign({},base,{k31:base.k31*L(3),sG:base.sG*L(3),sE:base.sE*L(1),w2nr:base.w2nr*L(2),w3nr:base.w3nr*L(2),kQ22:base.kQ22*L(3),kQ23:base.kQ23*L(3),b32:Math.min(0.9,base.b32*L(0.7))});
  const up=D.sweep(M,e,I,1,1064).n.map(n=>n[2]), dn=D.sweep(M,e,I,-1,1064).n.map(n=>n[2]);
  const w=I.filter((_,k)=>Math.abs(up[k]/dn[k]-1)>0.05);
  if(w.length){ found++; console.log(JSON.stringify(Object.fromEntries(['k31','sG','sE','w2nr','w3nr','kQ22','kQ23','b32'].map(k=>[k,+e[k].toPrecision(3)]))),'hyst',Math.min(...w).toFixed(0),'-',Math.max(...w).toFixed(0)); }
}
console.log('found',found);
