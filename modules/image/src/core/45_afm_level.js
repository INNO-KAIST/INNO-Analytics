/* AFM background correction. Masks use 1=include; fit and apply are independent. */
(function(){
"use strict";
const G=typeof self!=="undefined"?self:globalThis, IC=G.IMG_CORE||(G.IMG_CORE={});
IC.afmLevel=function(p,o){
  o=o||{};const w=p.w,h=p.h,N=w*h,method=o.method;
  if(!Number.isSafeInteger(w)||!Number.isSafeInteger(h)||w<1||h<1||!Number.isSafeInteger(N)||!p.data||p.data.length!==N)throw new Error("AFM correction: invalid image dimensions");
  if(!["plane","rows","poly2","lineFlatten"].includes(method))throw new Error("AFM correction: invalid method");
  for(const m of [o.fitMask,o.applyMask])if(m&&m.length!==N)throw new Error("AFM correction: invalid mask dimensions");
  if(method==="lineFlatten"){
    if(![0,1,2].includes(o.degree))throw new Error("AFM correction: invalid degree");
    if(!["x","y"].includes(o.axis))throw new Error("AFM correction: invalid axis");
    for(const mask of [o.fitMask,o.applyMask])if(mask!=null){
      if(mask.length!==N)throw new Error("AFM correction: invalid mask dimensions");
      for(let i=0;i<N;i++)if(mask[i]!==0&&mask[i]!==1)throw new Error("AFM correction: invalid mask value (expected 0 or 1)");
    }
    const degree=o.degree,axis=o.axis,fit=o.fitMask,apply=o.applyMask,src=p.data;
    const length=axis==="x"?w:h,lines=axis==="x"?h:w,k=degree+1;
    if(N*8+length*(k*8+16)>1.2e9)throw new Error("AFM correction: image exceeds memory limit");
    const idx=(line,t)=>axis==="x"?line*w+t:t*w+line;
    let applicable=0;for(let i=0;i<N;i++)if((!apply||apply[i]===1)&&Number.isFinite(src[i]))applicable++;
    if(!applicable)throw new Error("AFM correction: apply region has no finite pixels");
    const out=Float64Array.from(src);
    for(let line=0;line<lines;line++){
      let hasApply=false,n=0,base=0;
      for(let t=0;t<length;t++){const i=idx(line,t);if((!apply||apply[i]===1)&&Number.isFinite(src[i]))hasApply=true;}
      if(!hasApply)continue;
      for(let t=0;t<length;t++){const i=idx(line,t);if((!fit||fit[i]===1)&&Number.isFinite(src[i])){if(!n)base=src[i];n++;}}
      const label=(axis==="x"?"row ":"column ")+(line+1);
      if(n<k)throw new Error(`AFM correction: ${label} has insufficient fit samples`);
      // Match rows' summation order when possible; scale only on overflow.
      if(degree===0){
        let sum=0,overflow=false;for(let t=0;t<length;t++){const i=idx(line,t);if((!fit||fit[i]===1)&&Number.isFinite(src[i])){const delta=src[i]-base;if(!Number.isFinite(delta))throw new Error(`AFM correction: ${label} numeric range overflow`);sum+=delta;if(!Number.isFinite(sum)){overflow=true;break;}}}
        let mean=sum/n;if(overflow){mean=0;for(let t=0;t<length;t++){const i=idx(line,t);if((!fit||fit[i]===1)&&Number.isFinite(src[i])){const delta=src[i]-base;if(!Number.isFinite(delta))throw new Error(`AFM correction: ${label} numeric range overflow`);mean+=delta/n;if(!Number.isFinite(mean))throw new Error(`AFM correction: ${label} numeric range overflow`);}}}
        for(let t=0;t<length;t++){const i=idx(line,t);if((!apply||apply[i]===1)&&Number.isFinite(src[i])){const result=(src[i]-base)-mean;if(!Number.isFinite(result))throw new Error(`AFM correction: ${label} numeric range overflow`);out[i]=result;}}
        continue;
      }
      const cols=Array.from({length:k},()=>new Float64Array(n)),target=new Float64Array(n);
      let z=0;for(let t=0;t<length;t++){const i=idx(line,t);if((!fit||fit[i]===1)&&Number.isFinite(src[i])){
        const x=length>1?2*t/(length-1)-1:0;cols[0][z]=1;if(k>1)cols[1][z]=x;if(k>2)cols[2][z]=x*x;
        target[z]=src[i]-base;if(!Number.isFinite(target[z]))throw new Error(`AFM correction: ${label} numeric range overflow`);z++;
      }}
      const R=Array.from({length:k},()=>new Float64Array(k)),qTy=new Float64Array(k);
      for(let j=0;j<k;j++){const v=cols[j];for(let pass=0;pass<2;pass++)for(let i=0;i<j;i++){const q=cols[i];let dot=0;for(let a=0;a<n;a++)dot+=q[a]*v[a];R[i][j]+=dot;for(let a=0;a<n;a++)v[a]-=dot*q[a];}
        let norm=0;for(let a=0;a<n;a++)norm+=v[a]*v[a];norm=Math.sqrt(norm);
        if(!(norm>1e-12*Math.sqrt(n)))throw new Error(`AFM correction: ${label} singular fit`);
        R[j][j]=norm;let dot=0;for(let a=0;a<n;a++){v[a]/=norm;dot+=v[a]*target[a];}if(!Number.isFinite(dot))throw new Error(`AFM correction: ${label} numeric range overflow`);qTy[j]=dot;
      }
      // The 1-norm condition number of R equals that of the fitted QR
      // system up to norm choice; reject fits that amplify roundoff excessively.
      let rNorm=0,invNorm=0;
      for(let j=0;j<k;j++){
        let colSum=0;for(let i=0;i<=j;i++)colSum+=Math.abs(R[i][j]);rNorm=Math.max(rNorm,colSum);
        const v=new Float64Array(k);for(let i=k-1;i>=0;i--){let x=i===j?1:0;for(let m=i+1;m<k;m++)x-=R[i][m]*v[m];v[i]=x/R[i][i];}
        let sum=0;for(let i=0;i<k;i++)sum+=Math.abs(v[i]);invNorm=Math.max(invNorm,sum);
      }
      if(!Number.isFinite(rNorm*invNorm)||rNorm*invNorm>1e8)throw new Error(`AFM correction: ${label} ill-conditioned fit`);
      const coef=new Float64Array(k);for(let j=k-1;j>=0;j--){let v=qTy[j];for(let i=j+1;i<k;i++)v-=R[j][i]*coef[i];coef[j]=v/R[j][j];if(!Number.isFinite(coef[j]))throw new Error(`AFM correction: ${label} numeric range overflow`);}
      for(let t=0;t<length;t++){const i=idx(line,t);if((!apply||apply[i]===1)&&Number.isFinite(src[i])){const x=length>1?2*t/(length-1)-1:0;let v=coef[0];if(k>1)v+=coef[1]*x;if(k>2)v+=coef[2]*x*x;const result=(src[i]-base)-v;if(!Number.isFinite(result))throw new Error(`AFM correction: ${label} numeric range overflow`);out[i]=result;}}
    }
    return {w,h,bit:32,data:out};
  }
  const fit=o.fitMask,apply=o.applyMask,src=p.data,include=i=>(!fit||fit[i])&&Number.isFinite(src[i]),use=i=>(!apply||apply[i])&&Number.isFinite(src[i]);
  let applicable=0;for(let i=0;i<N;i++)if(use(i))applicable++;if(!applicable)throw new Error("AFM correction: apply region has no finite pixels");
  const k=method==="poly2"?6:method==="plane"?3:1;if(N*(8*(k+3)+16)>1.2e9)throw new Error("AFM correction: image exceeds memory limit");
  const out=Float64Array.from(src);
  if(method==="rows"){
    for(let y=0;y<h;y++){let n=0,sum=0,base=0,hasApply=false;for(let x=0;x<w;x++){const i=y*w+x;if(use(i))hasApply=true;if(include(i)){if(!n)base=src[i];const delta=src[i]-base;if(!Number.isFinite(delta))throw new Error("AFM correction: numeric range overflow");sum+=delta;if(!Number.isFinite(sum))throw new Error("AFM correction: numeric range overflow");n++;}}
      if(hasApply&&!n)throw new Error(`AFM correction: row ${y+1} has insufficient fit samples`);if(!n)continue;const meanOffset=sum/n;for(let x=0;x<w;x++){const i=y*w+x;if(use(i)){const v=(src[i]-base)-meanOffset;if(!Number.isFinite(v))throw new Error("AFM correction: numeric range overflow");out[i]=v;}}}
    return {w,h,bit:32,data:out};
  }
  const ids=[];let base=0;for(let i=0;i<N;i++)if(include(i)){if(!ids.length)base=src[i];ids.push(i);}const n=ids.length;if(n<k)throw new Error("AFM correction: insufficient fit samples");
  const coords=i=>{const x=i%w,y=Math.floor(i/w);return [w>1?2*x/(w-1)-1:0,h>1?2*y/(h-1)-1:0];};
  const terms=i=>{const [x,y]=coords(i);return k===3?[1,x,y]:[1,x,y,x*x,x*y,y*y];};
  const cols=Array.from({length:k},()=>new Float64Array(n)),target=new Float64Array(n);for(let j=0;j<n;j++){const t=terms(ids[j]);for(let c=0;c<k;c++)cols[c][j]=t[c];target[j]=src[ids[j]]-base;if(!Number.isFinite(target[j]))throw new Error("AFM correction: numeric range overflow");}
  const R=Array.from({length:k},()=>new Float64Array(k)),qTy=new Float64Array(k);
  for(let j=0;j<k;j++){const v=cols[j];for(let pass=0;pass<2;pass++)for(let i=0;i<j;i++){const q=cols[i];let dot=0;for(let z=0;z<n;z++)dot+=q[z]*v[z];R[i][j]+=dot;for(let z=0;z<n;z++)v[z]-=dot*q[z];}
    let norm=0;for(let z=0;z<n;z++)norm+=v[z]*v[z];norm=Math.sqrt(norm);if(!(norm>1e-12*Math.sqrt(n)))throw new Error("AFM correction: singular fit or insufficient independent samples");R[j][j]=norm;let dot=0;for(let z=0;z<n;z++){v[z]/=norm;dot+=v[z]*target[z];}if(!Number.isFinite(dot))throw new Error("AFM correction: numeric range overflow");qTy[j]=dot;}
  const coef=new Float64Array(k);for(let i=k-1;i>=0;i--){let v=qTy[i];for(let j=i+1;j<k;j++)v-=R[i][j]*coef[j];coef[i]=v/R[i][i];if(!Number.isFinite(coef[i]))throw new Error("AFM correction: numeric range overflow");}
  for(let i=0;i<N;i++)if(use(i)){const t=terms(i);let v=0;for(let j=0;j<k;j++)v+=coef[j]*t[j];const result=(src[i]-base)-v;if(!Number.isFinite(result))throw new Error("AFM correction: numeric range overflow");out[i]=result;}
  return {w,h,bit:32,data:out};
};
})();
