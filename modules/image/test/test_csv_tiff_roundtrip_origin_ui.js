const {chromium}=require('playwright'),assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),{pathToFileURL}=require('node:url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH});try{const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');
for(const origin of ['color-component','numeric']){
 await p.evaluate(origin=>{IMG_UI._addImage({name:'plane-'+origin,w:2,h:2,bit:16,signalOrigin:origin,planes:[new Float32Array([1,32769,65535,42])]});IMG_UI.rerender()},origin);
 const exp=p.waitForEvent('download');await p.click('#imSaveT');const raw=fs.readFileSync(await (await exp).path());await p.evaluate(()=>{IMG_UI._close(IMG_UI._cur().id,true);IMG_UI.rerender()});
 await p.setInputFiles('#imFin',{name:'roundtrip.tif',mimeType:'image/tiff',buffer:raw});await p.waitForFunction(()=>IMG_UI._state().images.length===1);const state=await p.evaluate(()=>({origin:IMG_UI._cur().signalOrigin,v:Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0])}));assert.equal(state.origin,origin);assert.deepEqual(state.v,[1,32769,65535,42]);
 await p.click('#imSaveCsv');if(origin==='numeric'){const exp2=p.waitForEvent('download',{timeout:3000});await p.click('#imCsvExport');await exp2;}else{await p.click('#imCsvExport');await p.locator('#imCsvError').waitFor({state:'visible',timeout:3000});await p.click('#imModal [data-close]')}
 await p.evaluate(()=>{IMG_UI._close(IMG_UI._cur().id,true);IMG_UI.rerender()});
}console.log('Own TIFF roundtrip preserves component rejection, numeric export, and 16-bit values');
}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
