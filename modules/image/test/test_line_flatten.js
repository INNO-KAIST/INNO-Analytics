const assert=require('assert');
require('../src/core/45_afm_level.js');
const level=globalThis.IMG_CORE.afmLevel;
const img=(w,h,a)=>({w,h,data:Float64Array.from(a)});
const run=(p,degree,axis='x',extra={})=>level(p,{method:'lineFlatten',degree,axis,...extra});
const near=(a,b,t=1e-10)=>assert.ok(Math.abs(a-b)<=t,`${a} != ${b}`);
const p=img(4,2,[10,12,14,16,20,21,24,29]);
assert.deepStrictEqual(Array.from(run(p,0).data),Array.from(level(p,{method:'rows'}).data));
{ const q=img(3,1,[0,Number.MIN_VALUE,Number.MIN_VALUE]); assert.deepStrictEqual(Array.from(run(q,0).data),Array.from(level(q,{method:'rows'}).data)); }
for(const axis of ['x','y']) for(const degree of [0,1,2]) {
  const w=5,h=4,src=[];for(let y=0;y<h;y++)for(let x=0;x<w;x++){const t=axis==='x'?x:y, l=axis==='x'?y:x;src.push(1e9+l*3+(degree>=1?2*t:0)+(degree>=2?t*t:0));}
  const q=img(w,h,src), out=run(q,degree,axis); assert.ok(out.data instanceof Float64Array); out.data.forEach(v=>near(v,0,1e-6)); assert.deepStrictEqual(Array.from(q.data),src);
}
{
 const q=img(5,2,[0,1,52,3,4,10,12,14,16,18]),fit=Uint8Array.from([1,1,0,1,1,1,1,1,1,1]),apply=Uint8Array.from([1,1,1,1,0,0,0,0,0,0]);
 const o=run(q,1,'x',{fitMask:fit,applyMask:apply}).data; near(o[2],50); assert.strictEqual(o[4],4); assert.strictEqual(o[5],10);
}
{
 const q=img(4,2,[NaN,Infinity,-Infinity,3,1,2,3,4]); const o=run(q,0).data; assert.ok(Number.isNaN(o[0]));assert.strictEqual(o[1],Infinity);assert.strictEqual(o[2],-Infinity);near(o[3],0);
}
assert.throws(()=>run(img(3,1,[1,2,3]),2,'x',{fitMask:[1,0,1]}),/row 1.*insufficient/i);
assert.throws(()=>run(img(2,2,[1,2,3,4]),2,'y'),/column 1.*insufficient/i);
assert.throws(()=>run(p,1,'x',{applyMask:new Uint8Array(8)}),/no finite/i);
assert.throws(()=>run(img(2,1,[Number.MAX_VALUE,-Number.MAX_VALUE]),0),/overflow/i);
assert.throws(()=>run(p,3),/degree/i);assert.throws(()=>run(p,1,'z'),/axis/i);
assert.throws(()=>run(p,1,'x',{fitMask:[1]}),/mask/i);
assert.throws(()=>run(p,1,'x',{fitMask:[1,1,1,1,1,1,1,NaN]}),/mask/i);
{
 const n=800000, a=new Float64Array(n), fit=new Uint8Array(n), apply=new Uint8Array(n);
 for(let i=0;i<n;i++)a[i]=(i/(n-1))**2;
 fit.fill(1,0,3);apply[n-1]=1;
 assert.throws(()=>run({w:n,h:1,data:a},2,'x',{fitMask:fit,applyMask:apply}),/row 1.*(ill.condition|singular)/i);
}
{
 const a=[0].concat(Array(15).fill(1e308)), o=run(img(16,1,a),0).data;
 near(o[0]/1e308,-15/16,1e-14);near(o[15]/1e308,1/16,1e-14);
}
console.log('line flatten numerical tests passed');

