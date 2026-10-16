const O=require("./oracle.js"), C=require("./cmp.js");
const shapes={
cross:[
"...........",
".....#.....",
".....#.....",
".....#.....",
".#########.",
".....#.....",
".....#.....",
"..........."],
T:[
"...........",
".#########.",
".....#.....",
".....#.....",
".....#.....",
"..........."],
line:[
"........",
".######.",
"........"],
looptail:[
"..........",
"...###....",
"..#...#...",
"..#...####",
"..#...#...",
"...###....",
".........."],
};
const cases=[]; for(const [k,rows] of Object.entries(shapes)) for(const pe of [0,1]) cases.push(Object.assign(O.fromRows(rows),{lsp:true,pruneEnds:!!pe,name:k+" pe"+pe}));
const os=O.run(cases);
cases.forEach((c,i)=>{ const o=os[i]; console.log("==",c.name); console.log(o.raw.split("\n").filter(l=>!l.startsWith("TAG")).join("\n")); console.log(O.ascii(c.w,c.h,c.data,new Map([...o.tag].map(([p,v])=>[p,v===30?"E":v===70?"J":v===127?"s":"?"])))); });
