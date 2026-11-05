const assert = require('node:assert/strict');
const IO = require('../src/img_io.js'); require('../src/img_afm.js');
function fixture({type, values, little=true, nx=2, ny=2, channels=1}) {
  const base=type & ~0x40, sizes={2:4,4:8,32:4}, size=sizes[base];
  const bytes=new Uint8Array(384+values.length*size), d=new DataView(bytes.buffer);
  d.setUint16(0,5,little); d.setInt32(76,values.length,little); d.setUint16(80,type,little);
  [nx,ny,channels,0].forEach((v,i)=>d.setInt32(132+4*i,v,little));
  values.forEach((v,i)=>{const o=384+i*size; if(base===4)d.setFloat64(o,v,little); else if(base===2)d.setFloat32(o,v,little); else d.setUint32(o,v,little);});
  return bytes;
}
(async()=>{
  for(const little of [true,false]) {
    let bytes=fixture({type:4,values:[1,1+2**-40,NaN,3,4,5,6,7],little,channels:2});
    const r=await IO.readFile('x.ibw',bytes); assert.equal(r.images.length,2);
    assert.deepEqual([r.images[0].w,r.images[0].h],[2,2]);
    assert.deepEqual(Array.from(r.images[0].sourceSamples),[NaN,3,1,1+2**-40]);
    assert.deepEqual(Array.from(r.images[1].sourceSamples),[6,7,4,5]);
    assert.deepEqual(r.images[0].sourceMeta.dimensions,[2,2,2,0]);
    assert.equal(r.images[0].sourceMeta.channelIndex,0);
    assert.equal(r.images[0].sourceMeta.orientation,'y-reversed');
    assert.equal(r.images[0].valueUnit,'');
    assert.equal(r.images[0].planes[0][3],1+2**-40);
    r.images[0].planes[0][3]=9; assert.equal(r.images[0].sourceSamples[3],1+2**-40);
    const scaled=fixture({type:4,values:[1e-9,2e-9,3e-9,4e-9],little}); scaled[212]=109;
    const scaledImage=(await IO.readFile('x.ibw',scaled)).images[0];
    assert.equal(scaledImage.valueUnit,'nm'); assert.equal(scaledImage.sourceSamples[2],1e-9);
    assert.notEqual(scaledImage.planes[0][2],scaledImage.sourceSamples[2]);
    bytes=fixture({type:96,values:[16777216,16777217,0,4294967295],little});
    const integer=(await IO.readFile('x.ibw',bytes)).images[0];
    assert.deepEqual(Array.from(integer.sourceSamples),[0,4294967295,16777216,16777217]);
    const f32=(await IO.readFile('x.ibw',fixture({type:2,values:[1.5,2.5,3.5,NaN],little}))).images[0];
    assert.deepEqual(Array.from(f32.sourceSamples),[3.5,NaN,1.5,2.5]);
    await assert.rejects(IO.readFile('x.ibw',bytes.slice(0,-1)),/자료가 짧습니다/);
  }
  const malformed=fixture({type:4,values:[1,2,3,4]}); new DataView(malformed.buffer).setInt32(140,2,true); await assert.rejects(IO.readFile('x.ibw',malformed),/차원/);
  const oldLimit=IO.MAXBYTES; IO.MAXBYTES=32; try{await assert.rejects(IO.readFile('x.ibw',fixture({type:4,values:[1,2,3,4]})),/너무 큽니다/);}finally{IO.MAXBYTES=oldLimit;}
  console.log('IBW precision: 2 endian variants passed');
})().catch(e=>{console.error(e);process.exitCode=1;});
