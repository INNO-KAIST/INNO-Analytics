const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),zlib=require('node:zlib');
const IO=require('../src/img_io.js');require('../src/img_lif.js');require('../src/img_oir.js');
IO.inflate=async(u,raw)=>new Uint8Array(raw?zlib.inflateRawSync(u):zlib.inflateSync(u));
function npy(shape){const d=`{'descr': '<u2', 'fortran_order': False, 'shape': (${shape.join(', ')},), }`,pad=(16-((10+d.length+1)%16))%16,s=d+' '.repeat(pad)+'\n',b=Buffer.alloc(10+Buffer.byteLength(s)+shape.reduce((a,b)=>a*b,1)*2);b.set([0x93,78,85,77,80,89,1,0]);b.writeUInt16LE(Buffer.byteLength(s),8);b.write(s,10);return new Uint8Array(b);}
(async()=>{
 const a=(await IO.readFile('color.npy',npy([5,6,3]))).images[0];assert.equal(a.signalOrigin,'color-component');
 const b=(await IO.readFile('numeric.npy',npy([5,6]))).images[0];assert.equal(b.signalOrigin,'numeric');
 for(const [f,want] of [['rgb8.oir','color-component'],['gray16.oir','numeric']]){const q=new Uint8Array(fs.readFileSync(path.join(__dirname,'io/oir',f))),im=(await IO.readFile(f,q)).images[0];assert.equal(im.signalOrigin,want);}
 // Missing RGB pixel block must not become a numerical multi-channel measurement.
 for(const count of [1,2]){const q=Buffer.from(fs.readFileSync(path.join(__dirname,'io/oir/rgb8.oir'))),uids=['_0_1_1a2b3c4d-5e6f-4a1b-9c2d-3e4f5a6b7c8d_0','_0_1_1a2b3d5e-5e6f-4a1b-9c2d-3e4f5a6b7c8e_0'];for(const uid of uids.slice(0,count)){const i=q.indexOf(Buffer.from(uid));assert(i>=0);q[i+5]=0x78;}const im=(await IO.readFile('missing-rgb.oir',new Uint8Array(q))).images[0];assert.equal(im.nc,3-count);assert(['unknown','color-component'].includes(im.signalOrigin),'incomplete RGB must never be labeled numeric');}
 let snapshot=IO.writeNumericSnapshot({name:'a',w:1,h:1,nc:1,nz:1,nt:1,bit:32,rgb:false,signalOrigin:'numeric',planes:[new Float64Array([1])]});
 let im=(await IO.readFile('a.innoimg',snapshot)).images[0];assert.equal(im.signalOrigin,'numeric');
 const n=new DataView(snapshot.buffer,snapshot.byteOffset,snapshot.byteLength).getUint32(8,true), meta=JSON.parse(new TextDecoder().decode(snapshot.subarray(12,12+n)));delete meta.signalOrigin;
 const txt=new TextEncoder().encode(JSON.stringify(meta)), old=new Uint8Array(12+txt.length+8);old.set(snapshot.subarray(0,8));new DataView(old.buffer).setUint32(8,txt.length,true);old.set(txt,12);old.set(snapshot.subarray(12+n),12+txt.length);
 assert.equal((await IO.readFile('old.innoimg',old)).images[0].signalOrigin,'unknown');
 console.log('NPY/OIR, snapshot and legacy signal origin verified');
})().catch(e=>{console.error(e);process.exitCode=1});
