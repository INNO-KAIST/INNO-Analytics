/* core/*.js 를 이름 순으로 읽어 globalThis.IMG_CORE 를 만든다 (Worker 와 같은 순서) */
const fs=require("fs"), path=require("path");
const dir=path.join(__dirname,"..","src","core");
for(const f of fs.readdirSync(dir).filter(f=>f.endsWith(".js")).sort()) (0,eval)(fs.readFileSync(path.join(dir,f),"utf8"));
module.exports=globalThis.IMG_CORE;
