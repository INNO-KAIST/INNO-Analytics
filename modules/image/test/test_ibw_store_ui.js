const {chromium}=require('playwright'),assert=require('node:assert/strict'),path=require('node:path'),{pathToFileURL}=require('node:url');
(async()=>{const browser=await chromium.launch({headless:true,executablePath:process.env.INNO_CHROMIUM_PATH||undefined,args:['--allow-file-access-from-files']});try{
const page=await browser.newPage();await page.goto(pathToFileURL(path.resolve('INNO Analytics.html')).href);await page.click('#nav button[data-v="imgx"]');await page.waitForSelector('.imx');
await page.evaluate(async()=>{IMG_UI._addImage({name:'legacy',w:2,h:2,bit:8,planes:[new Float32Array([1,2,3,4])]});await IMG_UI._flush()});
const get=()=>page.evaluate(()=>new Promise((ok,fail)=>{const q=indexedDB.open('innoanalytics_image_v2',1);q.onsuccess=()=>{const r=q.result.transaction('kv').objectStore('kv').get('state');r.onsuccess=()=>ok(r.result);r.onerror=()=>fail(r.error)}}));
assert.equal((await get()).images.length,1);
await page.evaluate(async()=>{IMG_UI._addImage({name:'fail',w:2,h:2,bit:8,planes:[new Float32Array(4)]});const orig=IDBObjectStore.prototype.put;IDBObjectStore.prototype.put=function(v,k){if(String(k).startsWith('px:'))throw new DOMException('synthetic','DataCloneError');return orig.call(this,v,k)};try{await IMG_UI._flush()}finally{IDBObjectStore.prototype.put=orig}});
await page.waitForTimeout(200);assert.equal((await get()).images.length,1,'failed transaction must not commit state');
await page.evaluate(()=>new Promise((resolve,reject)=>{const q=indexedDB.open('innoanalytics_image_v2',1);q.onsuccess=()=>{const tx=q.result.transaction('kv','readwrite'),st=tx.objectStore('kv');const r=st.get('state');r.onsuccess=()=>{const state=r.result,im=state.images[0];delete im.hasSource;delete im.originalId;st.put(state,'state');st.put([new Float32Array([1,2,3,4])],'px:'+im.id)};tx.oncomplete=resolve;tx.onerror=()=>reject(tx.error)}}));
await page.reload();await page.click('#nav button[data-v="imgx"]');await page.waitForFunction('IMG_UI._state().images.length===1'); assert.equal(await page.evaluate(()=>IMG_UI._cur().name),'legacy');
assert.deepEqual(await page.evaluate(()=>Array.from(IMG_UI._px.get(IMG_UI._cur().id)[0])),[1,2,3,4]);
await page.evaluate(async()=>{
 const desc=Object.getOwnPropertyDescriptor(IDBTransaction.prototype,'oncomplete');
 Object.defineProperty(IDBTransaction.prototype,'oncomplete',{configurable:true,get:desc.get,set(fn){return desc.set.call(this,fn?function(e){setTimeout(()=>fn.call(this,e),150)}:fn)}});
 try{const im=IMG_UI._addImage({name:'transient',w:2,h:2,bit:8,planes:[new Float32Array(4)]});const a=IMG_UI._flush();await new Promise(r=>setTimeout(r,30));IMG_UI._close(im.id,true);const b=IMG_UI._flush();await Promise.all([a,b]);}finally{Object.defineProperty(IDBTransaction.prototype,'oncomplete',desc)}
});
await page.reload();await page.click('#nav button[data-v="imgx"]');await page.waitForFunction('IMG_UI._state().images.length===1');assert.equal(await page.evaluate(()=>IMG_UI._cur().name),'legacy');
await page.evaluate(async()=>{const orig=IDBObjectStore.prototype.put;let fail=true;IDBObjectStore.prototype.put=function(v,k){if(fail&&String(k).startsWith('px:')){fail=false;throw new DOMException('synthetic','DataCloneError')}return orig.call(this,v,k)};try{IMG_UI._addImage({name:'retry',w:2,h:2,bit:8,planes:[new Float32Array(4)]});await Promise.all([IMG_UI._flush(),IMG_UI._flush()]);}finally{IDBObjectStore.prototype.put=orig}});
await page.reload();await page.click('#nav button[data-v="imgx"]');await page.waitForFunction('IMG_UI._state().images.length===2');assert.equal(await page.evaluate(()=>IMG_UI._state().images.some(i=>i.name==='retry')),true);
console.log('DataCloneError rollback, delayed flush/close, retry and legacy restore passed');
}finally{await browser.close()}})().catch(e=>{console.error(e);process.exitCode=1});
