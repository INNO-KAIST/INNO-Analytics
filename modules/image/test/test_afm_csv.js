const assert=require('node:assert/strict');
const IO=require('../src/img_io.js');
const {TextDecoder}=require('node:util');
const decoder=new TextDecoder('utf-8');
function entries(bytes){
  const d=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength), out={};let p=0;
  while(d.getUint32(p,true)===0x04034b50){const flags=d.getUint16(p+6,true),method=d.getUint16(p+8,true),n=d.getUint16(p+26,true),extra=d.getUint16(p+28,true),size=d.getUint32(p+18,true);assert.equal(method,0);assert(flags&0x0800);const name=decoder.decode(bytes.subarray(p+30,p+30+n));assert(!name.includes('/')&&!name.includes('\\'));assert(!(name in out));out[name]=bytes.slice(p+30+n+extra,p+30+n+extra+size);p+=30+n+extra+size;}
  assert.equal(d.getUint32(p,true),0x02014b50);return out;
}
const vals=new Float64Array([0,-0,1+2**-52,-3.25,NaN,Infinity,-Infinity,1e-300,Number.MAX_VALUE]);
const original=vals.slice();const plane={name:'../실측:sample',w:3,h:3,data:vals,cal:{pw:0.2,ph:0.3,unit:'µm'},valueUnit:'nm',c:0,z:2,t:1,historyIncomplete:true,ops:[{op:'x'}]};
assert.equal(typeof IO.exportAfmCsvBundle,'function');
for(const format of ['matrix','coordinates']){
 const result=IO.exportAfmCsvBundle({format,planes:[plane]});const zip=entries(result.bytes);const names=Object.keys(zip);assert.equal(names.length,2);const m=JSON.parse(decoder.decode(zip['metadata.json']));assert.equal(m.format,format);assert.equal(m.planes.length,1);assert.equal(m.planes[0].valueUnit,'nm');assert.equal(m.planes[0].xy.unit,'µm');assert.equal(m.planes[0].xy.pw,0.2);assert.equal(m.planes[0].xy.ph,0.3);assert.equal(m.planes[0].c,0);assert.equal(m.planes[0].z,2);assert.equal(m.planes[0].t,1);
 const rows=decoder.decode(zip[m.planes[0].file]).trimEnd().split('\n').map(s=>s.replace(/\r$/,'').split(','));let cells;
 if(format==='matrix'){assert.equal(rows.length,3);assert(rows.every(r=>r.length===3));cells=rows.flat();}
 else {assert.deepEqual(rows.shift(),['x','y','value']);assert.deepEqual(rows[0].slice(0,2),['0','0']);assert.deepEqual(rows[4].slice(0,2),['0.2','0.3']);cells=rows.map(row=>row[2]);}
 for(let i=0;i<vals.length;i++){const got=Number(cells[i]);assert(Object.is(got,vals[i])||(Number.isNaN(got)&&Number.isNaN(vals[i])),`value ${i}: ${cells[i]}`)}
}
assert.deepEqual(vals,original);
const two=entries(IO.exportAfmCsvBundle({format:'matrix',planes:[plane,plane]}).bytes);assert.equal(Object.keys(two).length,3);
const giantHistory={op:'x',comment:'a'.repeat(12000),toJSON(){throw new Error('JSON was invoked before metadata budget check')}};
assert.throws(()=>IO.exportAfmCsvBundle({format:'matrix',planes:Array(200).fill({...plane,w:1,h:1,data:new Float64Array(1),ops:[giantHistory]})}),/metadata exceeds the memory limit/);
const pixel=JSON.parse(decoder.decode(entries(IO.exportAfmCsvBundle({format:'coordinates',planes:[{...plane,cal:null}]}).bytes)['metadata.json']));assert.equal(pixel.planes[0].xy.unit,'px');
const invalidCal=JSON.parse(decoder.decode(entries(IO.exportAfmCsvBundle({format:'coordinates',planes:[{...plane,cal:{pw:NaN,ph:3,unit:'nm'},ops:[{threshold:NaN,offset:-0,limit:Infinity}]}]}).bytes)['metadata.json']));assert.equal(invalidCal.planes[0].xy.unit,'px');assert.deepEqual(invalidCal.planes[0].ops,[{threshold:'NaN',offset:'-0',limit:'Infinity'}]);
const unknown=IO.exportAfmCsvBundle({format:'coordinates',planes:[{...plane,cal:{pw:2,ph:3,unit:'unknown'}}]});assert.equal(unknown.manifest.planes[0].xy.unit,'px');
const bogus=IO.exportAfmCsvBundle({format:'coordinates',planes:[{...plane,cal:{pw:2,ph:3,unit:'bogus'}}]});assert.equal(bogus.manifest.planes[0].xy.unit,'px');
if(typeof SharedArrayBuffer!=='undefined')assert.throws(()=>IO.exportAfmCsvBundle({format:'matrix',planes:[{...plane,data:new Float64Array(new SharedArrayBuffer(9*8))}]}));
for(const f of [()=>IO.exportAfmCsvBundle({format:'matrix',planes:[]}),()=>IO.exportAfmCsvBundle({format:'bad',planes:[plane]}),()=>IO.exportAfmCsvBundle({format:'matrix',planes:[{...plane,rgb:true}]}),()=>IO.exportAfmCsvBundle({format:'matrix',planes:[{...plane,w:4}]}),()=>IO.exportAfmCsvBundle({format:'matrix',planes:[plane],maxBytes:10}),()=>IO.exportAfmCsvBundle({format:'coordinates',planes:[{...plane,cal:{pw:Number.MAX_VALUE,ph:1,unit:'nm'}}]})])assert.throws(f);
console.log('AFM CSV matrix/coordinates roundtrip, ZIP and budget validation passed');
