const O=require("./oracle.js"), C=require("./cmp.js");
// loop with a tail: phase A edge (tail) then phase B loop (self-loop on the cluster)
const rows=[
"..........",
"...###....",
"..#...#...",
"..#...####",
"..#...#...",
"...###....",
".........."];
const c=Object.assign(O.fromRows(rows),{lsp:true});
const o=O.run([c])[0]; console.log(o.raw); const r=C.mine(c); console.log(JSON.stringify(r.trees),JSON.stringify(r.edges),r.warn);
