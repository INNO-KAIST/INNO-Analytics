const {chromium}=require('playwright'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),{pathToFileURL}=require('node:url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH});try{const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');
await p.evaluate(()=>{IMG_UI._addImage({name:'component',w:4,h:3,bit:16,signalOrigin:'color-component',planes:[new Float32Array([0,1,2,3,4,5,6,7,8,9,10,11])]});IMG_UI.rerender()});
const event=p.waitForEvent('download');await p.click('#imSaveP');const raw=fs.readFileSync(await (await event).path());await p.evaluate(()=>{IMG_UI._close(IMG_UI._cur().id,true);IMG_UI.rerender()});await p.setInputFiles('#imFin',{name:'view.png',mimeType:'image/png',buffer:raw});await p.waitForFunction(()=>IMG_UI._state().images.length===1);
const origin=await p.evaluate(()=>IMG_UI._cur().signalOrigin);assert(['display-gray','color-component'].includes(origin));assert.notEqual(origin,'numeric');
console.log('PNG view reopens as display or color, never claimed as numeric measurement');
}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
