const assert=require('node:assert/strict'),IC=require('./load_core.js');
function plane(w,h,fn){return {w,h,bit:32,data:Float64Array.from({length:w*h},(_,i)=>fn(i%w,Math.floor(i/w)))}}function near(a,b=0,tol=1e-9){assert(Math.abs(a-b)<=tol,`${a} vs ${b}`)}
const p=plane(5,4,(x,y)=>5+2*x-3*y),r=IC.afmLevel(p,{method:'plane'});assert(r.data instanceof Float64Array);r.data.forEach(x=>near(x));assert.equal(p.data[0],5);
const quad=plane(3,3,(x,y)=>x*x+y*y);IC.afmLevel(quad,{method:'poly2'}).data.forEach(x=>near(x));
const row=plane(3,2,(x,y)=>1+x+3*y);assert.deepEqual(Array.from(IC.afmLevel(row,{method:'rows'}).data),[-1,0,1,-1,0,1]);
const noisy=plane(4,3,(x,y)=>5+2*x-3*y);noisy.data[11]+=50;const fit=new Uint8Array(12).fill(1);fit[11]=0;const apply=new Uint8Array(12);apply[11]=1;const masked=IC.afmLevel(noisy,{method:'plane',fitMask:fit,applyMask:apply});assert.equal(masked.data[0],noisy.data[0]);near(masked.data[11],50);
const special=plane(4,3,(x,y)=>5+2*x-3*y);special.data[0]=NaN;special.data[1]=Infinity;special.data[2]=-Infinity;const out=IC.afmLevel(special,{method:'plane'});assert(Number.isNaN(out.data[0])&&out.data[1]===Infinity&&out.data[2]===-Infinity);out.data.slice(3).forEach(x=>near(x));
assert.throws(()=>IC.afmLevel(plane(2,2,(x,y)=>x+y),{method:'poly2'}),/sample|rank|singular/i);const line=new Uint8Array(20);line.fill(1,0,5);assert.throws(()=>IC.afmLevel(p,{method:'plane',fitMask:line}),/sample|rank|singular/i);
const empty=new Uint8Array(6);assert.throws(()=>IC.afmLevel(row,{method:'rows',fitMask:empty}),/sample|row/i);
const big=plane(6,5,(x,y)=>1e12+2*x-3*y);big.cal={pw:0.1,ph:300,unit:'nm'};IC.afmLevel(big,{method:'plane'}).data.forEach(x=>near(x,0,0.001));
assert.throws(()=>IC.afmLevel(p,{method:'plane',applyMask:new Uint8Array(20)}),/apply|region/i);
const extreme=plane(2,2,(x,y)=>(x+y)%2?1e308:-1e308);assert.throws(()=>IC.afmLevel(extreme,{method:'plane'}),/overflow|range/i);assert.throws(()=>IC.afmLevel(extreme,{method:'rows'}),/overflow|range/i);
console.log('AFM plane, row, polynomial, masks, special values and stable fits passed');

