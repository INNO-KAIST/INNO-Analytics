const {chromium}=require('playwright'), assert=require('node:assert/strict'), path=require('path'),{pathToFileURL}=require('url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH});try{
const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');
await p.evaluate(()=>{IMG_UI._addImage({name:'height',signalOrigin:'numeric',w:3,h:2,bit:32,planes:[new Float64Array([1,2,3,4,5,6])],cal:{pw:2,ph:3,unit:'nm'},valueUnit:'nm'});IMG_UI.rerender()});
const before=await p.evaluate(()=>({data:Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0]),ops:IMG_UI._cur().ops.length}));
await p.click('#imSaveCsv',{timeout:3000});await p.selectOption('#imCsvFormat','coordinates');await p.selectOption('#imCsvSource','current');await p.selectOption('#imCsvScope','current');const event=p.waitForEvent('download');await p.click('#imCsvExport');const d=await event;assert(d.suggestedFilename().endsWith('.zip'));
assert.deepEqual(await p.evaluate(()=>({data:Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0]),ops:IMG_UI._cur().ops.length})),before);
console.log('CSV UI selected-plane coordinate download is non-mutating');
}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
