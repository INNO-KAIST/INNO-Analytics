const {chromium}=require('playwright'),assert=require('node:assert/strict'),path=require('path'),{pathToFileURL}=require('url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH});try{const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');
await p.evaluate(()=>{IMG_UI._addImage({name:'many.tif',signalOrigin:'numeric',w:1,h:1,bit:32,nz:4000,planes:Array.from({length:4000},(_,i)=>new Float32Array([i]))},{file:'many.tif'});IMG_UI.rerender()});
const t=performance.now(),pending=p.waitForEvent('download',{timeout:10000});await p.click('#imSaveP');const dl=await pending,elapsed=performance.now()-t;assert(dl.suggestedFilename().includes('Z1'));assert(elapsed<5000,`4000-plane PNG export took ${Math.round(elapsed)}ms`);
console.log(`Browser 4000-plane PNG export: ${Math.round(elapsed)}ms`);
}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
