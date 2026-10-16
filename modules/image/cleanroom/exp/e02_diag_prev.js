const O=require("./oracle.js"), C=require("./cmp.js");
const rows=[
"#............",
"#.....###....",
"#....#...#...",
"#....#...####",
"#....#...#...",
"......###....",
"............."];
const c=Object.assign(O.fromRows(rows),{lsp:true,cal:[0.5,0.7,1]});
const o=O.run([c])[0]; console.log(o.raw); const r=C.mine(c); console.log(r.trees.map(t=>[t.lsp,t.spx,t.spy,t.spz].join(" ")).join("\n"), C.diff(c,o,r));
