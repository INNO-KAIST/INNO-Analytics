const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const path=require('node:path');
const {pathToFileURL}=require('node:url');
(async()=>{
  const browser=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH||undefined,args:['--allow-file-access-from-files']});
  try {
    const page=await browser.newPage();
    const errors=[]; page.on('pageerror',e=>errors.push(String(e)));
    await page.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);
    await page.click('#nav button[data-v="imgx"]'); await page.waitForSelector('.imx');
    await page.evaluate(async()=>{const b=new Uint8Array(416),d=new DataView(b.buffer); d.setUint16(0,5,true);d.setInt32(76,4,true);d.setUint16(80,4,true);d.setInt32(132,2,true);d.setInt32(136,2,true);[1,1+2**-40,3,4].forEach((x,i)=>d.setFloat64(384+8*i,x,true)); const r=await IMG_IO.readFile('precision.ibw',b); IMG_UI._addImage(r.images[0]);});
    async function verify(){return page.evaluate(()=>{const i=IMG_UI._cur(), s=IMG_UI._source.get(i.id); return {samples:Array.from(s.samples),type:s.type,little:s.littleEndian};});}
    const expected={samples:[3,4,1,1+2**-40],type:4,little:true};
    assert.deepEqual(await verify(),expected);
    const before=await page.evaluate(()=>Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0]));
    await page.evaluate(async()=>{await IMG_UI._run('invert',{});});
    const changed=await page.evaluate(()=>Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0]));
    assert.notDeepEqual(changed,before); assert.deepEqual(await verify(),expected);
    await page.evaluate(async()=>{IMG_UI._undo(); await IMG_UI._flush();});
    assert.deepEqual(await page.evaluate(()=>Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0])),before);
    assert.deepEqual(await verify(),expected);
    await page.reload(); await page.click('#nav button[data-v="imgx"]'); await page.waitForFunction('window.IMG_UI && IMG_UI._state().images.length===1');
    assert.deepEqual(await verify(),expected);
    assert.deepEqual(errors,[]);
    console.log('IBW source survives add, invert/undo and IndexedDB reload');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
