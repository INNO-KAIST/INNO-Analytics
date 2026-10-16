const H=require("./hyp.js");
const root=(T,E)=>{ const i=T.vCl.findIndex(c=>c<0); return i<0?0:i; };
const V={m:{dfs:"stackE",order:"fwd",le:false,allRoots:true,root,method:+(process.env.P||2),vox:(e)=>e.slabs,vox2:(e)=>e.slabs,nb:process.env.NB||"3x3",oob:process.env.OOB||"skip",cycleStartAll:true}};
const r=H.evaluate(+(process.env.P||2),V,+(process.env.N||100),+(process.env.SEED||21),0);
// print mismatches with symmetric difference only
const gen=require("./gen.js"); 
const cases=r.cases; r.cases.forEach((c,i)=>{ const o0=r.os[i], o1=r.os[cases.length+i]; const want=new Set(); for(const [p] of o0.tag) if(!o1.tag.has(p)) want.add([p%c.w,((p/c.w)|0)%c.h,(p/(c.w*c.h))|0].join());
  const got=H.predict(c,V.m); const miss=[...want].filter(k=>!got.has(k)), extra=[...got].filter(k=>!want.has(k)); const tg=(k)=>{ const [x,y,z]=k.split(",").map(Number); return o0.tag.get(x+c.w*(y+c.h*z)); }; if(miss.length||extra.length) console.log(i,c.name,c.w+"x"+c.h+"x"+c.d,c.cal.map(x=>x.toFixed(3)).join(","),"want-only",miss.map(k=>k+"["+tg(k)+"]").join(" "),"| got-only",extra.map(k=>k+"["+tg(k)+"]").join(" ")); });
console.log(r.res);
