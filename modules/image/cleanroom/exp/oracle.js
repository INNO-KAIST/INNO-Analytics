// Black-box oracle runner for experiments (runs the compiled SkelDriver; never inspects it).
// cases: [{w,h,d,data:Uint8Array(0/255), cal:[pw,ph,pd], prune:0..3, pruneEnds:bool, lsp:bool, orig:Uint8Array|null}]
const fs=require("fs"), cp=require("child_process"), path=require("path"), os=require("os");
const IO=require("../../src/img_io.js");
function writeVol(f,w,h,d,data){ const planes=[]; for(let z=0;z<d;z++) planes.push(Float32Array.from(data.subarray(z*w*h,(z+1)*w*h)));
  fs.writeFileSync(f,IO.writeTIFF({w,h,bit:8,rgb:false,nc:1,nz:d,nt:1,planes,cal:null})); }
function run(cases,opts){ opts=opts||{}; const TMP=fs.mkdtempSync(path.join(os.tmpdir(),"askcr-")); const list=[];
  cases.forEach((c,i)=>{ const d=c.d||1; const f=path.join(TMP,`c${i}.tif`); writeVol(f,c.w,c.h,d,c.data); let of="";
    if(c.orig){ of=path.join(TMP,`o${i}.tif`); writeVol(of,c.w,c.h,d,c.orig); }
    const cal=c.cal||[1,1,1]; list.push([f,c.prune||0,c.pruneEnds?1:0,c.lsp===false?0:1,...cal].concat(of?[of]:[]).join("\t")); });
  fs.writeFileSync(path.join(TMP,"list.txt"),list.join("\n"));
  const out=cp.execSync(`JAVA_TOOL_OPTIONS= java -Xmx${opts.mem||"2g"} -cp /tmp/ij.jar:/tmp/askel_out SkelDriver -list ${path.join(TMP,"list.txt")} 2>/dev/null`,{encoding:"utf8",maxBuffer:1<<30});
  fs.rmSync(TMP,{recursive:true});
  return out.split(/^CASE \d+\n/m).slice(1).map(parse); }
function parse(block){ const L=block.trim().split("\n"); const r={trees:[],edges:[],tag:new Map(),lsp:new Set(),raw:block};
  for(const l of L){ const p=l.split(" ");
    if(p[0]==="T") r.trees.push(p.slice(2).map(Number));
    else if(p[0]==="E") r.edges.push(p.slice(1).map(Number));
    else if(p[0]==="TAG") for(const s of p.slice(1)){ const [a,b]=s.split(":").map(Number); r.tag.set(a,b); }
    else if(p[0]==="LSPIMG") for(const s of p.slice(1)) r.lsp.add(+s); }
  return r; }
// pretty ASCII grid of an image (2D) with marks from a map idx->char
function ascii(w,h,data,mark){ let s=""; for(let y=0;y<h;y++){ for(let x=0;x<w;x++){ const i=x+w*y; s+=mark&&mark.has(i)?mark.get(i):data[i]?"#":"."; } s+="\n"; } return s; }
// build a 2D image from strings ('#' or any non-'.' char = foreground)
function fromRows(rows){ const h=rows.length,w=rows[0].length,data=new Uint8Array(w*h); rows.forEach((r,y)=>{ for(let x=0;x<w;x++) if(r[x]!=="."&&r[x]!==" ") data[x+w*y]=255; }); return {w,h,d:1,data}; }
module.exports={run,parse,ascii,fromRows,writeVol};
