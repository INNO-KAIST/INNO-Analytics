const {chromium}=require('playwright'),a=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),{pathToFileURL}=require('node:url');
(async()=>{const b=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH||undefined});try{const p=await b.newPage({acceptDownloads:true});await p.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await p.click('#nav button[data-v="imgx"]');await p.waitForSelector('.imx');await p.evaluate(()=>{IMG_UI._addImage({name:'nan profile',w:3,h:1,bit:32,planes:[new Float64Array([1,NaN,3])],valueUnit:'nm'});IMG_UI._setActive({id:'l',type:'line',pts:[[0.5,0.5],[2.5,0.5]]});IMG_UI.rerender()});a(await p.locator('.im-profile .afm-boundary').isVisible());a((await p.locator('.im-profile .afm-boundary').innerText()).includes('clamp')); const dl=p.waitForEvent('download');await p.click('#imProfCsv');const rows=fs.readFileSync(await (await dl).path(),'utf8').trim().split(/\r?\n/);a.equal(Number(rows[1].split(',')[1]),1);a.equal(Number(rows[3].split(',')[1]),3);await p.click('#imProfStep');a.equal(await p.evaluate(()=>IMG_UI._state().results.at(-1).Delta),2);console.log('profile CSV endpoints agree with endpoint height and boundary policy');}finally{await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});






