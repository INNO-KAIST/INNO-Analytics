const assert=require('node:assert/strict'),fs=require('fs'),path=require('path'),IO=require('../src/img_io.js');require('../src/img_micro.js');require('../src/img_lif.js');
(async()=>{async function read(file){const f=path.join(__dirname,'io',file);return (await IO.readFile(file,new Uint8Array(fs.readFileSync(f)))).images;}
const nd=await read('n_rgb16.nd2');assert.equal(nd[0].colorChannels,true);assert.equal(nd[0].rgb,false);assert.equal(nd[0].bit,16);
const ndGray=await read('n_z4t3c2.nd2');assert.equal(ndGray[0].colorChannels,false);assert.equal(ndGray[0].bit,16);
const lif=await read('lif/rgb16.lif');assert.equal(lif[0].colorChannels,true);assert.equal(lif[0].bit,16);
const lifGray=await read('lif/dyes.lif');assert.equal(lifGray[0].colorChannels,false);
const czi=await read('c_mix_rgb_gray.czi');assert.equal(czi[0].colorChannels,false);assert.equal(czi[1].rgb,true);
console.log('ND2, LIF, CZI RGB and numeric provenance preserved');})().catch(e=>{console.error(e);process.exitCode=1});
