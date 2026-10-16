// sensitivity on the full image: which orig pixels change the set of removed voxels
const gen=require("./gen.js"), O=require("./oracle.js");
const [seed,N,idx]=process.argv.slice(2).map(Number); const c=gen(seed,N)[idx]; const val=+(process.env.VAL||0); const P=+(process.env.P||2);
const cases=[Object.assign({},c,{prune:0,lsp:false}),Object.assign({},c,{prune:P,lsp:false})];
for(let i=0;i<c.orig.length;i++){ const o=c.orig.slice(); o[i]=val; cases.push(Object.assign({},c,{orig:o,prune:P,lsp:false,i})); }
const os=O.run(cases); const W=c.w,H=c.h; const xy=(p)=>[p%W,((p/W)|0)%H,(p/(W*H))|0].join(",");
const remS=(o)=>new Set([...os[0].tag.keys()].filter(p=>!o.tag.has(p)));
const r0=remS(os[1]); console.log("removed:",[...r0].map(xy).join(" "));
for(let k=2;k<cases.length;k++){ const r=remS(os[k]); const a=[...r].filter(p=>!r0.has(p)), b=[...r0].filter(p=>!r.has(p)); if(a.length||b.length) console.log("pixel",xy(cases[k].i),"→ +",a.map(xy).join(" "),"−",b.map(xy).join(" ")); }
