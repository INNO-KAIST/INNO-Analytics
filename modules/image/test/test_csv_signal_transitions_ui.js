const {chromium}=require('playwright'),assert=require('node:assert/strict'),path=require('node:path'),{pathToFileURL}=require('node:url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH});try{const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');
const outcome=await p.evaluate(async()=>{
 const u=IMG_UI, add=(name,signalOrigin,rgb=false)=>u._addImage({name,w:2,h:2,bit:8,rgb,nc:rgb?3:1,signalOrigin,planes:(rgb?[1,2,3]:[1]).map(x=>new Float32Array([x,x+1,x+2,x+3]))});
 const rgb=add('rgb','color-component',true);await u._run('type',{to:'rgb2gray'});const gray=u._cur().signalOrigin,unit=u._cur().valueUnit;u._undo();const restored=u._cur().signalOrigin;
 const color=add('component','color-component'),num=add('numeric','numeric');await u._run('calc',{b:color.id,fn:'add',float:true});const combined=u._cur().signalOrigin;
 u._state().sel=color.id;await u._run('calc',{b:num.id,fn:'add',float:true});const reversed=u._cur().signalOrigin;
 const unknown=add('unknown','unknown');u._state().sel=num.id;await u._run('calc',{b:unknown.id,fn:'add',float:true});const mixedUnknown=u._cur().signalOrigin;
 return {gray,unit,restored,combined,reversed,mixedUnknown};
});
assert.deepEqual(outcome,{gray:'display-gray',unit:'',restored:'color-component',combined:'color-component',reversed:'color-component',mixedUnknown:'unknown'});
console.log('Gray conversion/undo, color composition in both orders and unknown propagation verified');
}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
