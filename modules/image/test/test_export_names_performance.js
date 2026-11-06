const assert=require('node:assert/strict'),IO=require('../src/img_io.js');
const make=(sourceName,channelName='Height')=>({sourceName,channelName,c:0,z:0,t:0});
const old=IO.exportPlaneNames([make('same.ibw'),make('same.ibw'),make('same.ibw')],'png');
assert.deepEqual(old,['same_Height_C1.png','same_Height_C1_2.png','same_Height_C1_3.png']);
const planes=Array.from({length:4000},()=>make('same.ibw'));
const t=performance.now(),names=IO.exportPlaneNames(planes,'csv'),elapsed=performance.now()-t;
assert.equal(new Set(names.map(s=>s.normalize('NFC').toLowerCase())).size,4000);
assert(elapsed<1500,`4000 identical stems took ${Math.round(elapsed)} ms`);
const traps=[make('same.ibw'),make('same.ibw'),make('same_Height_C1_2',''),make('same.ibw'),make('SAME.ibw')];
const result=IO.exportPlaneNames(traps,'png');assert.equal(new Set(result.map(s=>s.normalize('NFC').toLowerCase())).size,result.length);
function previousAlgorithm(items){const bases=items.map(IO.exportPlaneStem),counts=new Map(),used=new Set();for(const b of bases){const k=b.normalize('NFC').toLowerCase();counts.set(k,(counts.get(k)||0)+1)}return items.map((p,i)=>{let b=bases[i],k=b.normalize('NFC').toLowerCase();if(counts.get(k)>1)b+='_C'+(p.c+1);const stem=b;let number=2;while(used.has((b+'.png').normalize('NFC').toLowerCase()))b=stem+'_'+number++;used.add((b+'.png').normalize('NFC').toLowerCase());return b+'.png'})}
const source=['x','X','x_Height_C1','x_Height_C1_2','é','e\u0301'];for(let round=0;round<40;round++){const a=Array.from({length:80},(_,i)=>make(source[(i*7+round*13+((i*i+round)%7))%source.length],i%4?'Height':''));assert.deepEqual(IO.exportPlaneNames(a,'png'),previousAlgorithm(a));}
console.log(`4000 duplicate planes: ${Math.round(elapsed)} ms, all names unique`);
