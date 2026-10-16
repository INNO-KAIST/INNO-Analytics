const gen=require("./gen.js"); const [seed,N,idx,x0,y0,x1,y1,z]=process.argv.slice(2).map(Number); const c=gen(seed,N)[idx]; const zz=z||0;
for(let y=y0;y<=y1;y++){ let l=String(y).padStart(3)+": "; for(let x=x0;x<=x1;x++){ const i=x+c.w*(y+c.h*zz); l+=(x<0||y<0||x>=c.w||y>=c.h)?"  OOB   ":String(c.orig[i]).padStart(4)+(c.data[i]?"#":" ")+"   "; } console.log(l); }
