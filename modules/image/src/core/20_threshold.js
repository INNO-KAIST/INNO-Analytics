/* ============================================================
   Image ▸ Adjust ▸ Threshold  —  ij/process/AutoThresholder.java (G. Landini 이식본) 그대로
   + ImageProcessor.setAutoThreshold 의 16/32-bit 처리(8-bit 로 비례 변환한 256 칸 히스토그램에서 구하고 되돌림)
   + Auto Local Threshold (Bernsen·Contrast·Mean·Median·MidGrey·Niblack·Phansalkar·Sauvola; 원형 창)
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;
const EPS=2.220446049250313e-16;
const idiv=(a,b)=>Math.trunc(a/b);

function bilevel(d){ let a=-1,b=-1,n=0; for(let i=0;i<d.length;i++) if(d[i]>0){ n++; if(n>2) return -1; if(a===-1) a=i; else b=i; }
  if(n===1) return a-1; if(n===2) return b-1; return -1; }
function IJIsoData(data){
  const maxValue=data.length-1; let result=0; const count0=data[0]; data[0]=0; const countMax=data[maxValue]; data[maxValue]=0;
  let min=0; while(data[min]===0&&min<maxValue) min++; let max=maxValue; while(data[max]===0&&max>0) max--;
  if(min>=max){ data[0]=count0; data[maxValue]=countMax; return data.length>>1; }
  let mi=min;
  do{ let s1=0,s2=0,s3=0,s4=0; for(let i=min;i<=mi;i++){ s1+=i*data[i]; s2+=data[i]; } for(let i=mi+1;i<=max;i++){ s3+=i*data[i]; s4+=data[i]; }
    result=(s1/s2+s3/s4)/2; mi++; }while((mi+1)<=result&&mi<max-1);
  data[0]=count0; data[maxValue]=countMax; return Math.round(result);
}
function defaultIsoData(data){
  const n=data.length, d2=Array.from(data); let mode=0,maxCount=0;
  for(let i=0;i<n;i++) if(d2[i]>maxCount){ maxCount=d2[i]; mode=i; }
  let maxCount2=0; for(let i=0;i<n;i++) if(d2[i]>maxCount2&&i!==mode) maxCount2=d2[i];
  if(maxCount>maxCount2*2&&maxCount2!==0){ d2[mode]=Math.trunc(maxCount2*1.5); }
  return IJIsoData(d2);
}
function Huang(data){
  let first=0; for(let ih=0;ih<data.length;ih++) if(data[ih]!==0){ first=ih; break; }
  let last=data.length-1; for(let ih=data.length-1;ih>=first;ih--) if(data[ih]!==0){ last=ih; break; }
  const term=1/(last-first); const mu0=new Float64Array(data.length), mu1=new Float64Array(data.length);
  let sp=0,np=0; for(let ih=first;ih<data.length;ih++){ sp+=ih*data[ih]; np+=data[ih]; mu0[ih]=sp/np; }
  sp=0; np=0; for(let ih=last;ih>0;ih--){ sp+=ih*data[ih]; np+=data[ih]; mu1[ih-1]=sp/np; }
  let thr=-1, minEnt=Number.MAX_VALUE;
  for(let it=0;it<data.length;it++){ let ent=0;
    for(let ih=0;ih<=it;ih++){ const mx=1/(1+term*Math.abs(ih-mu0[it])); if(!(mx<1e-6||mx>0.999999)) ent+=data[ih]*(-mx*Math.log(mx)-(1-mx)*Math.log(1-mx)); }
    for(let ih=it+1;ih<data.length;ih++){ const mx=1/(1+term*Math.abs(ih-mu1[it])); if(!(mx<1e-6||mx>0.999999)) ent+=data[ih]*(-mx*Math.log(mx)-(1-mx)*Math.log(1-mx)); }
    if(ent<minEnt){ minEnt=ent; thr=it; } }
  return thr;
}
function bimodalTest(y){ let modes=0; for(let k=1;k<y.length-1;k++) if(y[k-1]<y[k]&&y[k+1]<y[k]){ modes++; if(modes>2) return false; } return modes===2; }
function Intermodes(data){
  const n=data.length, h=Float64Array.from(data); let iter=0;
  while(!bimodalTest(h)){ let prev=0,cur=0,next=h[0]; for(let i=0;i<n-1;i++){ prev=cur; cur=next; next=h[i+1]; h[i]=(prev+cur+next)/3; } h[n-1]=(cur+next)/3; iter++; if(iter>10000) return -1; }
  let tt=0; for(let i=1;i<n-1;i++) if(h[i-1]<h[i]&&h[i+1]<h[i]) tt+=i; return Math.floor(tt/2);
}
function IsoData(data){
  let g=0; for(let i=1;i<data.length;i++) if(data[i]>0){ g=i+1; break; }
  for(;;){ let l=0,totl=0; for(let i=0;i<g;i++){ totl+=data[i]; l+=data[i]*i; } let h=0,toth=0; for(let i=g+1;i<data.length;i++){ toth+=data[i]; h+=data[i]*i; }
    if(totl>0&&toth>0){ l=idiv(l,totl); h=idiv(h,toth); if(g===Math.round((l+h)/2)) break; }
    g++; if(g>data.length-2) return -1; }
  return g;
}
function Li(data){
  let num=0; for(let ih=0;ih<data.length;ih++) num+=data[ih];
  let mean=0; for(let ih=0;ih<data.length;ih++) mean+=ih*data[ih]; mean/=num;
  let newT=mean, oldT, thr;
  do{ oldT=newT; thr=Math.trunc(oldT+0.5);
    let sb=0,nb=0; for(let ih=0;ih<=thr;ih++){ sb+=ih*data[ih]; nb+=data[ih]; } const mb=nb===0?0:sb/nb;
    let so=0,no=0; for(let ih=thr+1;ih<data.length;ih++){ so+=ih*data[ih]; no+=data[ih]; } const mo=no===0?0:so/no;
    const temp=(mb-mo)/(Math.log(mb)-Math.log(mo));
    newT=temp< -EPS?Math.trunc(temp-0.5):Math.trunc(temp+0.5);
  }while(Math.abs(newT-oldT)>0.5);
  return thr;
}
function normHist(data){ let total=0; for(const v of data) total+=v; const nh=new Float64Array(data.length), P1=new Float64Array(data.length), P2=new Float64Array(data.length);
  for(let i=0;i<data.length;i++) nh[i]=data[i]/total; P1[0]=nh[0]; P2[0]=1-P1[0]; for(let i=1;i<data.length;i++){ P1[i]=P1[i-1]+nh[i]; P2[i]=1-P1[i]; }
  let first=0; for(let i=0;i<data.length;i++) if(!(Math.abs(P1[i])<EPS)){ first=i; break; }
  let last=data.length-1; for(let i=data.length-1;i>=first;i--) if(!(Math.abs(P2[i])<EPS)){ last=i; break; }
  return {nh,P1,P2,first,last}; }
function MaxEntropy(data){
  const {nh,P1,P2,first,last}=normHist(data); let thr=-1, maxEnt=Number.MIN_VALUE;
  for(let it=first;it<=last;it++){ let eb=0; for(let ih=0;ih<=it;ih++) if(data[ih]!==0) eb-=(nh[ih]/P1[it])*Math.log(nh[ih]/P1[it]);
    let eo=0; for(let ih=it+1;ih<data.length;ih++) if(data[ih]!==0) eo-=(nh[ih]/P2[it])*Math.log(nh[ih]/P2[it]);
    const t=eb+eo; if(maxEnt<t){ maxEnt=t; thr=it; } }
  return thr;
}
function Mean(data){ let tot=0,sum=0; for(let i=0;i<data.length;i++){ tot+=data[i]; sum+=i*data[i]; } return Math.floor(idiv(sum,tot)); }
function MinErrorI(data){
  const A=j=>{ let x=0; for(let i=0;i<=j;i++) x+=data[i]; return x; }, B=j=>{ let x=0; for(let i=0;i<=j;i++) x+=i*data[i]; return x; }, C=j=>{ let x=0; for(let i=0;i<=j;i++) x+=i*i*data[i]; return x; };
  const L=data.length-1; let thr=Mean(data), Tprev=-2;
  while(thr!==Tprev){
    const mu=B(thr)/A(thr), nu=(B(L)-B(thr))/(A(L)-A(thr)), p=A(thr)/A(L), q=(A(L)-A(thr))/A(L);
    const s2=C(thr)/A(thr)-mu*mu, t2=(C(L)-C(thr))/(A(L)-A(thr))-nu*nu;
    const w0=1/s2-1/t2, w1=mu/s2-nu/t2, w2=mu*mu/s2-nu*nu/t2+Math.log10((s2*(q*q))/(t2*(p*p)));
    const sq=w1*w1-w0*w2; if(sq<0) return thr;
    Tprev=thr; const temp=(w1+Math.sqrt(sq))/w0; if(temp!==temp) thr=Tprev; else thr=Math.floor(temp);
  }
  return thr;
}
function Minimum(data){
  const n=data.length; if(n<2) return 0; let iter=0, thr=-1, max=-1;
  const h=new Float64Array(n), t=new Float64Array(n); for(let i=0;i<n;i++){ h[i]=data[i]; if(data[i]>0) max=i; }
  while(!bimodalTest(h)){ for(let i=1;i<n-1;i++) t[i]=(h[i-1]+h[i]+h[i+1])/3; t[0]=(h[0]+h[1])/3; t[n-1]=(h[n-2]+h[n-1])/3; h.set(t); iter++; if(iter>10000) return -1; }
  for(let i=1;i<max;i++) if(h[i-1]>h[i]&&h[i+1]>=h[i]){ thr=i; break; }
  return thr;
}
function Moments(data){
  let total=0; for(const v of data) total+=v; const h=new Float64Array(data.length); for(let i=0;i<data.length;i++) h[i]=data[i]/total;
  let m0=1,m1=0,m2=0,m3=0; for(let i=0;i<data.length;i++){ m1+=i*h[i]; m2+=i*i*h[i]; m3+=i*i*i*h[i]; }
  const cd=m0*m2-m1*m1, c0=(-m2*m2+m1*m3)/cd, c1=(m0*-m3+m2*m1)/cd;
  const z0=0.5*(-c1-Math.sqrt(c1*c1-4*c0)), z1=0.5*(-c1+Math.sqrt(c1*c1-4*c0)); const p0=(z1-m1)/(z1-z0);
  let sum=0; for(let i=0;i<data.length;i++){ sum+=h[i]; if(sum>p0) return i; } return -1;
}
function Otsu(data){
  let N=0; for(const v of data) N+=v; const term=1/N, L=data.length;
  const h=new Float64Array(L), cnh=new Float64Array(L), mean=new Float64Array(L);
  for(let i=0;i<L;i++) h[i]=term*data[i]; cnh[0]=h[0]; for(let i=1;i<L;i++) cnh[i]=cnh[i-1]+h[i];
  mean[0]=0; for(let i=1;i<L;i++) mean[i]=mean[i-1]+i*h[i]; const tm=mean[L-1];
  let thr=-2147483648, maxB=0;
  for(let i=0;i<L;i++){ let b=tm*cnh[i]-mean[i]; b*=b/(cnh[i]*(1-cnh[i])); if(maxB<b){ maxB=b; thr=i; } }
  return thr;
}
function Percentile(data){
  let thr=-1; const ptile=0.5; let total=0; for(const v of data) total+=v; let temp=1, ps=0;
  for(let i=0;i<data.length;i++){ ps+=data[i]; const a=Math.abs(ps/total-ptile); if(a<temp){ temp=a; thr=i; } } return thr;
}
function RenyiEntropy(data){
  const {nh,P1,P2,first,last}=normHist(data); const L=data.length;
  let thr=0, maxEnt=0;
  for(let it=first;it<=last;it++){ let eb=0; for(let ih=0;ih<=it;ih++) if(data[ih]!==0) eb-=(nh[ih]/P1[it])*Math.log(nh[ih]/P1[it]);
    let eo=0; for(let ih=it+1;ih<L;ih++) if(data[ih]!==0) eo-=(nh[ih]/P2[it])*Math.log(nh[ih]/P2[it]); const t=eb+eo; if(maxEnt<t){ maxEnt=t; thr=it; } }
  let ts2=thr;
  thr=0; maxEnt=0; let term=1/(1-0.5);
  for(let it=first;it<=last;it++){ let eb=0; for(let ih=0;ih<=it;ih++) eb+=Math.sqrt(nh[ih]/P1[it]); let eo=0; for(let ih=it+1;ih<L;ih++) eo+=Math.sqrt(nh[ih]/P2[it]);
    const t=term*((eb*eo)>0?Math.log(eb*eo):0); if(t>maxEnt){ maxEnt=t; thr=it; } }
  let ts1=thr;
  thr=0; maxEnt=0; term=1/(1-2);
  for(let it=first;it<=last;it++){ let eb=0; for(let ih=0;ih<=it;ih++) eb+=(nh[ih]*nh[ih])/(P1[it]*P1[it]); let eo=0; for(let ih=it+1;ih<L;ih++) eo+=(nh[ih]*nh[ih])/(P2[it]*P2[it]);
    const t=term*((eb*eo)>0?Math.log(eb*eo):0); if(t>maxEnt){ maxEnt=t; thr=it; } }
  let ts3=thr, tmp;
  if(ts2<ts1){ tmp=ts1; ts1=ts2; ts2=tmp; } if(ts3<ts2){ tmp=ts2; ts2=ts3; ts3=tmp; } if(ts2<ts1){ tmp=ts1; ts1=ts2; ts2=tmp; }
  let b1,b2,b3;
  if(Math.abs(ts1-ts2)<=5){ if(Math.abs(ts2-ts3)<=5){ b1=1;b2=2;b3=1; } else { b1=0;b2=1;b3=3; } }
  else { if(Math.abs(ts2-ts3)<=5){ b1=3;b2=1;b3=0; } else { b1=1;b2=2;b3=1; } }
  const omega=P1[ts3]-P1[ts1];
  return Math.trunc(ts1*(P1[ts1]+0.25*omega*b1)+0.25*ts2*omega*b2+ts3*(P2[ts3]+0.25*omega*b3));
}
function Shanbhag(data){
  const {nh,P1,P2,first,last}=normHist(data); let thr=-1, minEnt=Number.MAX_VALUE;
  for(let it=first;it<=last;it++){ let eb=0, term=0.5/P1[it]; for(let ih=1;ih<=it;ih++) eb-=nh[ih]*Math.log(1-term*P1[ih-1]); eb*=term;
    let eo=0; term=0.5/P2[it]; for(let ih=it+1;ih<data.length;ih++) eo-=nh[ih]*Math.log(1-term*P2[ih]); eo*=term;
    const t=Math.abs(eb-eo); if(t<minEnt){ minEnt=t; thr=it; } }
  return thr;
}
function Triangle(data0){
  const data=Array.from(data0), L=data.length; let min=0,dmax=0,max=0,min2=0;
  for(let i=0;i<L;i++) if(data[i]>0){ min=i; break; } if(min>0) min--;
  for(let i=L-1;i>0;i--) if(data[i]>0){ min2=i; break; } if(min2<L-1) min2++;
  for(let i=0;i<L;i++) if(data[i]>dmax){ max=i; dmax=data[i]; }
  let inverted=false;
  if((max-min)<(min2-max)){ inverted=true; data.reverse(); min=L-1-min2; max=L-1-max; }
  if(min===max) return min;
  let nx=data[max], ny=min-max, d=Math.sqrt(nx*nx+ny*ny); nx/=d; ny/=d; d=nx*min+ny*data[min];
  let split=min, sd=0; for(let i=min+1;i<=max;i++){ const nd=nx*i+ny*data[i]-d; if(nd>sd){ split=i; sd=nd; } }
  split--;
  return inverted?(L-1-split):split;
}
function Yen(data){
  let total=0; for(const v of data) total+=v; const L=data.length, nh=new Float64Array(L), P1=new Float64Array(L), P1s=new Float64Array(L), P2s=new Float64Array(L);
  for(let i=0;i<L;i++) nh[i]=data[i]/total; P1[0]=nh[0]; for(let i=1;i<L;i++) P1[i]=P1[i-1]+nh[i];
  P1s[0]=nh[0]*nh[0]; for(let i=1;i<L;i++) P1s[i]=P1s[i-1]+nh[i]*nh[i];
  P2s[L-1]=0; for(let i=L-2;i>=0;i--) P2s[i]=P2s[i+1]+nh[i+1]*nh[i+1];
  let thr=-1, maxC=Number.MIN_VALUE;
  for(let it=0;it<L;it++){ const c=-1*((P1s[it]*P2s[it])>0?Math.log(P1s[it]*P2s[it]):0)+2*((P1[it]*(1-P1[it]))>0?Math.log(P1[it]*(1-P1[it])):0); if(c>maxC){ maxC=c; thr=it; } }
  return thr;
}
const METHODS={Default:defaultIsoData,Huang,Intermodes,IsoData,Li,MaxEntropy,Mean,MinError:MinErrorI,Minimum,Moments,Otsu,Percentile,RenyiEntropy,Shanbhag,Triangle,Yen,IJ_IsoData:d=>IJIsoData(Array.from(d))};
IC.THRESHOLD_METHODS=["Default","Huang","Intermodes","IsoData","Li","MaxEntropy","Mean","MinError","Minimum","Moments","Otsu","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"];

/* AutoThresholder.getThreshold(method, histogram) — 히스토그램 → 정수 문턱 */
IC.thresholdLevel=function(method,hist){
  const h=Array.from(hist,v=>Math.round(v));
  const b=bilevel(h); if(b>=0) return b;
  const fn=METHODS[method]; if(!fn) throw new Error("알 수 없는 threshold 방법: "+method);
  let t=fn(h); if(t===-1) t=0; return t;
};

/* ImageProcessor.setAutoThreshold(method, darkBackground) 와 같은 결과의 [lower, upper] (원래 값 단위).
   dark=true: 밝은 물체(어두운 배경) → [t+1, 상한]. dark=false: 어두운 물체 → [0(하한), t].
   16/32-bit: 화소 범위 [min,max] 로 8-bit 변환 후 256 칸에서 t 를 구해 lower=min+(t/255)(max−min) 로 되돌린다.
   반환 {lower, upper, level, hist(8-bit), min, max} */
IC.autoThreshold=function(p,method,dark,mask){
  let h8, mn=0, mx=255;
  if(p.bit===8){ h8=IC.histogram(p,mask).h; }
  else{ [mn,mx]=IC.minMax(p.data,mask); const d8=IC.to8bitScaled(p,mn,mx); h8=IC.histogram({w:p.w,h:p.h,bit:8,data:d8},mask).h; }
  const t=IC.thresholdLevel(method,h8);
  let lower,upper; const tmax=255;
  if(dark){ lower=t+1; upper=tmax; } else { lower=0; upper=t; }
  if(p.bit===8){ lower=Math.min(255,Math.max(0,lower)); upper=Math.min(255,Math.max(0,upper)); } // 빈 histogram 등에서 방법이 −2³¹ 을 돌려줘도 0..255 (ImageJ)
  if(p.bit!==8){
    if(lower>255) lower=255;
    if(mx>mn){
      if(lower===0){ lower=p.bit===32?Math.min(mn,-1e30):0; } else lower=mn+(lower/255)*(mx-mn);
      if(upper===255){ upper=p.bit===16?65535:Math.max(mx,1e30); } else upper=mn+(upper/255)*(mx-mn);
    } else { return {lower:NaN,upper:NaN,level:t,hist:h8,min:mn,max:mx,none:true}; } // 일정한 16/32-bit 이미지: ImageJ 는 문턱 없음(−1,−1)
    if(p.bit===16){ lower=Math.round(lower); upper=Math.round(upper); } // ShortProcessor.setThreshold: 정수로 반올림
  }
  return {lower,upper,level:t,hist:h8,min:mn,max:mx};
};
/* 같은 histogram 으로 모든 방법을 한 번에 (Try all) */
IC.tryAllThresholds=function(p,dark,mask){ const out={}; for(const m of IC.THRESHOLD_METHODS){ try{ out[m]=IC.autoThreshold(p,m,dark,mask); }catch(e){ out[m]={err:String(e.message||e)}; } } return out; };

/* Convert to Mask: [lower, upper] 안 = 255 (Black background 규칙), 밖 = 0 → 8-bit */
IC.applyThreshold=function(p,lower,upper){
  const d=p.data, out=new Float32Array(d.length);
  for(let i=0;i<d.length;i++){ const v=d[i]; out[i]=(v>=lower&&v<=upper)?255:0; }
  return {w:p.w,h:p.h,bit:8,data:out};
};

/* ---------- Auto Local Threshold (8-bit 전용; 원형 창 = RankFilters 커널) ----------
   식: Landini 의 Fiji 'Auto Local Threshold' 문서 정의.
   Niblack: t = mean + k·sd − c (k=0.2, c=0)        Sauvola: t = mean·(1 + k·(sd/r − 1)) (k=0.5, r=128)
   Phansalkar: 0..1 로 정규화, t = mean·(1 + p·exp(−q·mean) + k·(sd/r − 1)) (k=0.25, r=0.5, p=2, q=10)
   Bernsen: 대비(max−min) < c(15) 이면 중간값 ≥ 128 을 물체로, 아니면 화소 ≥ (max+min)/2
   Mean / Median / MidGrey: 화소 > 지역값 − c (c=0)     Contrast: max·min 중 가까운 쪽
   sd 는 모집단 표준편차(RankFilters VARIANCE 의 제곱근). 물체 = 255 (밝은 물체 기준; dark=false 면 반전 이미지에서) */
IC.localThreshold=function(p,method,o){
  if(p.bit!==8) throw new Error("Local threshold 는 8-bit 이미지에서만 씁니다 (Image ▸ Type ▸ 8-bit)");
  o=o||{}; const r=+o.radius||15, W=p.w,H=p.h, N=W*H; let src=p;
  if(o.dark===false){ const d=new Float32Array(N); for(let i=0;i<N;i++) d[i]=255-p.data[i]; src={w:W,h:H,bit:8,data:d}; }
  const v=src.data, out=new Float32Array(N);
  const P1=o.p1, P2=o.p2;
  const meanVar=()=>{ const m=IC.rank(Object.assign({},src,{bit:32}),"mean",{radius:r}).data, q=IC.rank(Object.assign({},src,{bit:32}),"variance",{radius:r}).data; return [m,q]; };
  if(method==="Niblack"||method==="Sauvola"){
    const [m,q]=meanVar(); const k=P1!=null&&P1!==0?+P1:(method==="Niblack"?0.2:0.5);
    const c2=P2!=null&&P2!==0?+P2:(method==="Niblack"?0:128);
    for(let i=0;i<N;i++){ const sd=Math.sqrt(q[i]); const t=method==="Niblack"?m[i]+k*sd-c2:m[i]*(1+k*(sd/c2-1)); out[i]=v[i]>t?255:0; }
  }else if(method==="Phansalkar"){
    const s01={w:W,h:H,bit:32,data:Float32Array.from(v,x=>x/255)};
    const m=IC.rank(s01,"mean",{radius:r}).data, q=IC.rank(s01,"variance",{radius:r}).data;
    const k=P1!=null&&P1!==0?+P1:0.25, rr=P2!=null&&P2!==0?+P2:0.5, pp=2, qq=10;
    for(let i=0;i<N;i++){ const sd=Math.sqrt(q[i]); const t=m[i]*(1+pp*Math.exp(-qq*m[i])+k*(sd/rr-1)); out[i]=(v[i]/255)>t?255:0; }
  }else if(method==="Mean"||method==="Median"||method==="MidGrey"||method==="Bernsen"||method==="Contrast"){
    const c=P1!=null?+P1:0;
    if(method==="Mean"){ const m=IC.rank(Object.assign({},src,{bit:32}),"mean",{radius:r}).data; for(let i=0;i<N;i++) out[i]=v[i]>(m[i]-c)?255:0; }
    else if(method==="Median"){ const m=IC.rank(src,"median",{radius:r}).data; for(let i=0;i<N;i++) out[i]=v[i]>(m[i]-c)?255:0; }
    else { const mn=IC.rank(src,"min",{radius:r}).data, mx=IC.rank(src,"max",{radius:r}).data;
      if(method==="MidGrey") for(let i=0;i<N;i++) out[i]=v[i]>((mx[i]+mn[i])/2-c)?255:0;
      else if(method==="Contrast") for(let i=0;i<N;i++) out[i]=(Math.abs(mx[i]-v[i])<=Math.abs(v[i]-mn[i]))?255:0;
      else { const ct=P1!=null&&P1!==0?+P1:15; for(let i=0;i<N;i++){ const con=mx[i]-mn[i], mid=(mx[i]+mn[i])/2; out[i]=con<ct?(mid>=128?255:0):(v[i]>=mid?255:0); } } }
  }else throw new Error("알 수 없는 local threshold: "+method);
  return {w:W,h:H,bit:8,data:out};
};
IC.LOCAL_METHODS=["Bernsen","Contrast","Mean","Median","MidGrey","Niblack","Phansalkar","Sauvola"];

})();
