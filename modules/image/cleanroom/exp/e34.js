// window extent probing vs calibration: rectangle loop, uniform 200, one dark background pixel at (20, 8+dy)
const O=require("./oracle.js"), CR=require("../askel_cr.js");
const w=40,h=26, v={w,h,d:1,data:new Uint8Array(w*h)}; const put=(x,y)=>v.data[x+w*y]=255;
for(let x=8;x<=31;x++){ put(x,8); put(x,17); } for(let y=9;y<=16;y++){ put(7,y); put(32,y); } for(let x=2;x<=6;x++) put(x,12);
const g=CR._t.debugGraph(v,{pw:1,ph:1,pd:1}); const L=g.edges.find(e=>e.slabs.length>20).slabs; const idx=(s)=>L.findIndex(q=>q[0]+","+q[1]===s);
console.log("list head",L.slice(0,4).map(p=>p.join()).join(" "),"...");
const cals=(process.env.CALS||"1,1;1,1.93;1,2.6;1,3.2;0.5,1;0.3,1;2,1;1,0.4;0.7,1.3").split(";").map(s=>s.split(",").map(Number));
const offs=[]; for(let dy=-4;dy<=4;dy++) offs.push([0,dy]);
const cases=[]; for(const cal of cals) for(const [dx,dy] of offs){ const o=new Uint8Array(w*h).fill(200); o[(20+dx)+w*(8+dy)]=0; cases.push(Object.assign({},v,{orig:o,prune:2,lsp:false,cal:[cal[0],cal[1],1],tag:`cal ${cal} dy ${dy}`})); }
for(const cal of cals) cases.push(Object.assign({},v,{prune:0,lsp:false,cal:[cal[0],cal[1],1],tag:"base"}));
const os=O.run(cases); const nb=cals.length*offs.length;
cals.forEach((cal,ci)=>{ const base=os[nb+ci]; const line=[]; offs.forEach(([dx,dy],k)=>{ const o=os[ci*offs.length+k]; const rem=[]; for(const [p] of base.tag) if(!o.tag.has(p)) rem.push(p%w+","+((p/w)|0)); line.push(`dy${dy}:${rem.join("/")}`); }); console.log("cal",cal.join(","),line.join("  ")); });
