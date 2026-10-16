// show oracle raw output (+ascii tag map) for named shapes, for each option set
const O=require("./oracle.js");
module.exports=function show(shapes,optsList,extra){ const cases=[]; for(const [k,rows] of Object.entries(shapes)) for(const op of optsList) cases.push(Object.assign(O.fromRows(rows),{lsp:true},op,{name:k+" "+JSON.stringify(op)}));
  const os=O.run(cases);
  cases.forEach((c,i)=>{ const o=os[i]; console.log("==",c.name); console.log(o.raw.split("\n").filter(l=>!l.startsWith("TAG")&&!l.startsWith("CASE")).join("\n"));
    const lsp=o.lsp; console.log(O.ascii(c.w,c.h,c.data,new Map([...o.tag].map(([p,v])=>[p,lsp.has(p)?"*":v===30?"E":v===70?"J":v===127?"s":"?"])))); });
  return {cases,os}; };
