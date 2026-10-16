// full-output comparison with the oracle for all prune modes × prune ends, on diverse random volumes
const C=require("./cmp.js"), gen=require("./gen.js");
const N=+(process.env.N||600), seed=+(process.env.SEED||5);
const base=gen(seed,N); const cases=[]; base.forEach((c,i)=>{ const P=i%4, pe=((i>>2)&1)===1||process.env.PE==="1"; cases.push(Object.assign({},c,{prune:P,pruneEnds:pe,name:c.name+` P${P} pe${pe?1:0}`})); });
const {nbad}=C.cmp(cases); console.log("seed",seed,"cases",cases.length,"mismatching",nbad);
