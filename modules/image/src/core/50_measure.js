/* ============================================================
   Analyze ▸ Measure / Analyze Particles
   원본: ij/gui/Wand.java (traceEdge, 8-연결 legacy), ij/gui/PolygonRoi.getTracedPerimeter,
         ij/process/FloatPolygon.getConvexHull, ij/gui/Roi.getFeretValues, ij/process/EllipseFitter,
         ij/plugin/filter/Analyzer.saveResults (Circ./AR/Round/Solidity), ij/plugin/filter/ParticleAnalyzer
   ============================================================ */
(function(){
"use strict";
const G=(typeof self!=="undefined")?self:globalThis;
const IC=G.IMG_CORE;

/* Wand.traceEdge (8-연결). inside(x,y) → bool. 시작: 물체 화소 (x,y) 에서 오른쪽 끝까지 간 경계.
   반환: 꼭짓점(화소 모서리 좌표) 배열 xs, ys — 방향이 바뀌는 점만 */
IC.traceOutline=function(inside,startX,startY,fourConnected){
  let x=startX, y=startY;
  if(inside(x,y)){ do{ x++; }while(inside(x,y)); }
  else { do{ x++; }while(!inside(x,y)); }
  const ins=(x,y,d)=>{ switch(d&3){ case 0: return inside(x,y); case 1: return inside(x,y-1); case 2: return inside(x-1,y-1); case 3: return inside(x-1,y); } return false; };
  const xs=[], ys=[];
  let sDir; let sx=x, sy=y;
  if(inside(sx,sy)) sDir=1; else { sDir=3; sy++; }
  x=sx; y=sy; let dir=sDir;
  do{
    let nd;
    if(fourConnected){ nd=dir; do{ if(!ins(x,y,nd)) break; nd++; }while(nd<dir+2); nd--; }
    else { nd=dir+1; do{ if(ins(x,y,nd)) break; nd--; }while(nd>=dir); }
    if(nd!==dir){ xs.push(x); ys.push(y); }
    switch(nd&3){ case 0: x++; break; case 1: y--; break; case 2: x--; break; case 3: y++; break; }
    dir=nd;
  }while(x!==sx||y!==sy||(dir&3)!==sDir);
  if(xs[0]!==x){ xs.push(x); ys.push(y); }
  return {xs,ys};
};
/* Wand (Worker 에서): 화소값이 [lo, hi] 인 8-연결 성분의 바깥 외곽선. 성분의 가장 위·왼쪽 화소에서 따라가야 늘 바깥 경계 (구멍 경계 방지).
   반환 {xs, ys, n(성분 화소 수)} · 클릭한 화소가 범위 밖이면 null */
IC.wand=function(p,x,y,lo,hi){ const W=p.w,H=p.h,d=p.data; x=Math.floor(x); y=Math.floor(y);
  const ins=(a,b)=>a>=0&&b>=0&&a<W&&b<H&&d[b*W+a]>=lo&&d[b*W+a]<=hi;
  if(!ins(x,y)) return null;
  const seen=new Uint8Array(W*H), st=new Int32Array(W*H); let sp=0, n=0; st[sp++]=y*W+x; seen[y*W+x]=1; let top=y*W+x;
  while(sp){ const i=st[--sp]; n++; if(i<top) top=i; const xx=i%W, yy=(i-xx)/W;
    for(let dy=-1;dy<=1;dy++){ const b2=yy+dy; if(b2<0||b2>=H) continue; for(let dx=-1;dx<=1;dx++){ const a2=xx+dx; if(a2<0||a2>=W) continue; const j=b2*W+a2; if(!seen[j]){ const v=d[j]; if(v>=lo&&v<=hi){ seen[j]=1; st[sp++]=j; } } } }
    if((n&0xfffff)===0) IC.progress(Math.min(0.9,n/(W*H))); }
  const ol=IC.traceOutline(ins,top%W,Math.floor(top/W),false); ol.n=n; return ol; };
/* PolygonRoi.getTracedPerimeter */
IC.tracedPerimeter=function(xs,ys,pw,ph){
  const n=xs.length; if(n<4) return 0; pw=pw||1; ph=ph||1;
  let sdx=0,sdy=0,nc=0; let dx1=xs[0]-xs[n-1], dy1=ys[0]-ys[n-1]; let side1=Math.abs(dx1)+Math.abs(dy1); let corner=false;
  for(let i=0;i<n;i++){ const ni=i+1===n?0:i+1; const dx2=xs[ni]-xs[i], dy2=ys[ni]-ys[i]; sdx+=Math.abs(dx1); sdy+=Math.abs(dy1); const side2=Math.abs(dx2)+Math.abs(dy2);
    if(side1>1||!corner){ corner=true; nc++; } else corner=false; dx1=dx2; dy1=dy2; side1=side2; }
  return sdx*pw+sdy*ph-(nc*((pw+ph)-Math.sqrt(pw*pw+ph*ph)));
};
/* 다각형 둘레(닫힘) */
IC.polyLength=function(xs,ys,pw,ph,closed){ pw=pw||1; ph=ph||1; let L=0; const n=xs.length; for(let i=1;i<n;i++) L+=Math.hypot((xs[i]-xs[i-1])*pw,(ys[i]-ys[i-1])*ph); if(closed&&n>1) L+=Math.hypot((xs[0]-xs[n-1])*pw,(ys[0]-ys[n-1])*ph); return L; };
/* FREEROI 둘레 (PolygonRoi.getSmoothedPerimeter = 3 점 이동평균 선 길이 + 끝→처음 직선) */
IC.smoothedPerimeter=function(xs,ys,pw,ph){ const n=xs.length; pw=pw||1; ph=ph||1; if(n<3) return IC.polyLength(xs,ys,pw,ph,true);
  const w2=pw*pw, h2=ph*ph; let L=0, dx=(xs[0]+xs[1]+xs[2])/3-xs[0], dy=(ys[0]+ys[1]+ys[2])/3-ys[0]; L+=Math.sqrt(dx*dx*w2+dy*dy*h2);
  for(let i=1;i<n-2;i++){ dx=(xs[i+2]-xs[i-1])/3; dy=(ys[i+2]-ys[i-1])/3; L+=Math.sqrt(dx*dx*w2+dy*dy*h2); }
  dx=xs[n-1]-(xs[n-3]+xs[n-2]+xs[n-1])/3; dy=ys[n-1]-(ys[n-3]+ys[n-2]+ys[n-1])/3; L+=Math.sqrt(dx*dx*w2+dy*dy*h2);
  dx=xs[n-1]-xs[0]; dy=ys[n-1]-ys[0]; L+=Math.sqrt(dx*dx*w2+dy*dy*h2); return L; };
/* FloatPolygon.getConvexHull (gift wrapping, 반시계, 공선점 제외) */
IC.convexHull=function(xp,yp){
  const n=xp.length; const xx=[], yy=[]; if(!n) return {xs:xx,ys:yy};
  const fx=xp.map(Math.fround), fy=yp.map(Math.fround);
  let sY=Infinity; for(let i=0;i<n;i++) if(fy[i]<sY) sY=fy[i];
  let sX=Infinity, p1=0; for(let i=0;i<n;i++) if(fy[i]===sY&&fx[i]<sX){ sX=fx[i]; p1=i; }
  const pstart=p1; const sq=v=>v*v;
  do{
    const x1=fx[p1], y1=fy[p1]; let p2=p1, x2=0, y2=0;
    do{ p2++; if(p2===n) p2=0; if(p2===p1) break; x2=fx[p2]; y2=fy[p2]; }while(x2===x1&&y2===y1);
    let p3=p2+1; if(p3===n) p3=0;
    if(p2!==p1) do{ const x3=fx[p3], y3=fy[p3]; const det=x1*(y2-y3)-y1*(x2-x3)+(y3*x2-y2*x3); let cf=false;
      if(det===0&&p3!==p2){ cf=(sq(x3-x1)+sq(y3-y1))>(sq(x2-x1)+sq(y2-y1)); }
      if(det>0||cf){ x2=x3; y2=y3; p2=p3; } p3++; if(p3===n) p3=0; }while(p3!==p1);
    xx.push(x1); yy.push(y1);
    if(p2===p1) break; p1=p2;
    if(xx.length>1&&fx[p1]===xx[0]&&fy[p1]===yy[0]) break;
  }while(p1!==pstart);
  return {xs:xx,ys:yy};
};
IC.polyArea=function(xs,ys){ let a=0; const n=xs.length; for(let i=0;i<n;i++){ const j=(i+1)%n; a+=xs[i]*ys[j]-xs[j]*ys[i]; } return Math.abs(a)/2; };
/* Roi.getFeretValues: [Feret, FeretAngle, MinFeret, FeretX, FeretY] (+ 끝점) */
IC.feret=function(xs,ys,pw,ph){
  pw=pw||1; ph=ph||1; const H=IC.convexHull(xs,ys); const X=H.xs, Y=H.ys, n=X.length; if(!n) return null;
  let p1=0,p2=0,dmax=0;
  for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ const dx=(X[i]-X[j])*pw, dy=(Y[i]-Y[j])*ph, d=dx*dx+dy*dy; if(d>dmax){ dmax=d; p1=i; p2=j; } }
  if(X[p1]>X[p2]){ const t=p1; p1=p2; p2=t; }
  const xf1=X[p1],yf1=Y[p1],xf2=X[p2],yf2=Y[p2];
  let ang=(180/Math.PI)*Math.atan2((yf1-yf2)*ph,(xf2-xf1)*pw); if(ang<0) ang+=180;
  let x0=X[n-1], y0=Y[n-1], minF=Number.MAX_VALUE; let mEnd=[0,0,0,0];
  for(let i=0;i<n;i++){ const xp=x0, yp=y0; x0=X[i]; y0=Y[i]; let xn=(y0-yp)*ph, yn=(xp-x0)*pw; const nf=1/Math.sqrt(xn*xn+yn*yn); xn*=nf*pw; yn*=nf*ph;
    let md=0; let e=[0,0,0,0]; for(let j=0;j<n;j++){ const dx=X[j]-x0, dy=Y[j]-y0; const dist=dx*xn+dy*yn; if(dist>md){ md=dist; e=[X[j],Y[j],X[j]-(xn/pw*dist)/pw,Y[j]-(yn/ph*dist)/ph]; } }
    if(md<minF){ minF=md; mEnd=e; } }
  return {feret:Math.sqrt(dmax),angle:ang,minFeret:minF,x:xf1,y:yf1,ends:[xf1,yf1,xf2,yf2],minEnds:mEnd,hull:H};
};
/* EllipseFitter (마스크의 2차 모멘트, 넓이 보존 보정). mask: rect 크기 Uint8Array 또는 null */
IC.ellipse=function(mask,rect){
  const sqrtPi=1.772453851, HALFPI=1.5707963267949; const {x:left,y:top,w:width,h:height}=rect;
  if(!mask){ let major=(width*2)/sqrtPi, minor=(height*2)/sqrtPi, angle=0, theta=0; if(major<minor){ const t=major; major=minor; minor=t; angle=90; theta=Math.PI/2; } return {major,minor,angle,theta,xc:left+width/2,yc:top+height/2}; }
  let xsum=0,ysum=0,x2sum=0,y2sum=0,xysum=0,bc=0;
  for(let y=0;y<height;y++){ let bl=0, xl=0; const o=y*width; for(let x=0;x<width;x++) if(mask[o+x]){ bl++; xl+=x; x2sum+=x*x; } xsum+=xl; ysum+=bl*y; xysum+=xl*y; y2sum+=y*y*bl; bc+=bl; }
  if(!bc) return {major:0,minor:0,angle:0,theta:0,xc:NaN,yc:NaN};
  x2sum+=0.08333333*bc; y2sum+=0.08333333*bc; const n=bc;
  const x1=xsum/n, y1=ysum/n, x2=x2sum/n, y2=y2sum/n, xy=xysum/n; const u20=x2-x1*x1, u02=y2-y1*y1, u11=xy-x1*y1;
  let m4=4*Math.abs(u02*u20-u11*u11); if(m4<0.000001) m4=0.000001;
  let a11=u02/m4, a12=u11/m4, a22=u20/m4; let tmp=a11-a22; if(tmp===0) tmp=0.000001;
  let theta=0.5*Math.atan(2*a12/tmp); if(theta<0) theta+=HALFPI;
  if(a12>0) theta+=HALFPI; else if(a12===0){ if(a22>a11){ theta=0; tmp=a22; a22=a11; a11=tmp; } else if(a11!==a22) theta=HALFPI; }
  tmp=Math.sin(theta); if(tmp===0) tmp=0.000001; const z=a12*Math.cos(theta)/tmp;
  let major=Math.sqrt(1/Math.abs(a22+z)), minor=Math.sqrt(1/Math.abs(a11-z)); const scale=Math.sqrt(bc/(Math.PI*major*minor));
  major=major*scale*2; minor=minor*scale*2; let angle=180*theta/Math.PI; if(angle===180) angle=0;
  if(major<minor){ const t=major; major=minor; minor=t; }
  return {major,minor,angle,theta,xc:left+x1+0.5,yc:top+y1+0.5};
};

/* ROI 하나를 ImageJ 'Measure' 처럼 잰다.
   roi = {mask:Uint8Array(rect.w*rect.h)|null, rect:{x,y,w,h}, poly:{xs,ys}|null, kind:"traced"|"rect"|"oval"|"polygon"|"freehand"}
   cal = {pw, ph, unit}. 반환은 ImageJ 열 이름 */
function isAxisInt(P){ const n=P.xs.length; if(n<3) return false; for(let i=0;i<n;i++){ const j=(i+1)%n; if(P.xs[i]!==Math.trunc(P.xs[i])||P.ys[i]!==Math.trunc(P.ys[i])) return false; if(!(P.xs[j]-P.xs[i]===0||P.ys[j]-P.ys[i]===0)) return false; } return true; }
IC.measureRoi=function(p,roi,cal){
  cal=cal||{pw:1,ph:1}; const pw=cal.pw||1, ph=cal.ph||1; const W=p.w;
  const r=roi.rect; void W;
  const s=IC.stats(p,roi.mask||null,r,{pw,ph},true);
  const out={Area:s.area,Mean:s.mean,StdDev:s.stdDev,Mode:s.mode,Min:s.min,Max:s.max,X:s.xCentroid,Y:s.yCentroid,XM:s.xCM,YM:s.yCM,
    IntDen:s.intDen,RawIntDen:s.rawIntDen,Median:s.median,Skew:s.skewness,Kurt:s.kurtosis,
    BX:r.x*pw,BY:r.y*ph,Width:r.w*pw,Height:r.h*ph,_n:s.n}; // 경계 상자는 영상 안으로 자른 것 (ImageJ 와 같음)
  const bb=roi.bounds||r; // 둘레: 영상 밖으로 나간 부분도 포함한 ROI 모양 (ImageJ 와 같음)
  // 둘레
  let per;
  if(roi.kind==="traced") per=IC.tracedPerimeter(roi.poly.xs,roi.poly.ys,pw,ph);
  else if(roi.kind==="rect") per=2*bb.w*pw+2*bb.h*ph;
  else if(roi.kind==="oval"){ const a=bb.w*pw/2, b=bb.h*ph/2; per=Math.PI*(3*(a+b)-Math.sqrt((3*a+b)*(a+3*b))); }
  else if(roi.kind==="freehand") per=IC.smoothedPerimeter(roi.poly.xs,roi.poly.ys,pw,ph);
  else if(roi.kind==="poly"&&roi.poly&&isAxisInt(roi.poly)) per=IC.tracedPerimeter(roi.poly.xs,roi.poly.ys,pw,ph); // ImageJ PolygonRoi.isTraced: 변이 모두 가로·세로이고 꼭짓점이 정수면 traced 둘레
  else if(roi.poly) per=IC.polyLength(roi.poly.xs,roi.poly.ys,pw,ph,true);
  else per=2*bb.w*pw+2*bb.h*ph;
  out["Perim."]=per;
  const el=IC.ellipse(roi.mask,r); const k=pw; // ImageJ: pw==ph 일 때 major·minor × pw
  out.Major=el.major*k; out.Minor=el.minor*k; out.Angle=el.angle;
  let circ=per===0?0:4*Math.PI*(s.area/(per*per)); if(circ>1) circ=1; out["Circ."]=circ;
  const poly=roi.poly||{xs:[r.x,r.x+r.w,r.x+r.w,r.x],ys:[r.y,r.y,r.y+r.h,r.y+r.h]};
  const fe=IC.feret(poly.xs,poly.ys,pw,ph);
  out.Feret=fe?fe.feret:NaN; out.FeretAngle=fe?fe.angle:NaN; out.MinFeret=fe?fe.minFeret:NaN; out.FeretX=fe?fe.x:NaN; out.FeretY=fe?fe.y:NaN;
  out.AR=el.major/el.minor; out.Round=4*s.area/(Math.PI*out.Major*out.Major);
  // Solidity = 화소 수 / 볼록 껍질 넓이(반올림 정수 다각형, 화소 단위)
  if(fe){ const hx=fe.hull.xs.map(Math.round), hy=fe.hull.ys.map(Math.round); const ca=IC.polyArea(hx,hy); out.Solidity=ca>0?s.n/ca:NaN; } else out.Solidity=NaN;
  out.ECD=2*Math.sqrt(s.area/Math.PI); // 같은 넓이의 원 지름 (ImageJ 기본 열 아님 — 입자 크기 분포용)
  out._fe=fe; out._el=el;
  return out;
};

/* ---------- Analyze Particles ----------
   bin: 문턱 판정 함수 대신 0/1 Uint8Array(w*h) (물체=1). intensity: 측정용 plane(redirect; 없으면 bin 영상)
   o: {minSize, maxSize (화소 수), minCirc, maxCirc, excludeEdges, includeHoles, fourConnected(false)}
   반환: {particles:[{n, rect, poly, mask, meas}], count, totalArea, labels:Int32Array} */
IC.analyzeParticles=function(bin,W,H,intensity,cal,o){
  o=o||{}; cal=cal||{pw:1,ph:1};
  const minS=o.minSize!=null?o.minSize:0, maxS=o.maxSize!=null?o.maxSize:Infinity, minC=o.minCirc||0, maxC=o.maxCirc!=null?o.maxCirc:1;
  const done=new Uint8Array(W*H); const labels=new Int32Array(W*H); const parts=[]; let count=0;
  const cur=new Uint8Array(W*H);
  for(let y=0;y<H;y++){ for(let x=0;x<W;x++){ const i0=y*W+x; if(!bin[i0]||done[i0]) continue;
    // 8-연결 성분 (처리 안 된 물체 화소)
    const st=[i0]; const pix=[]; cur[i0]=1; let minx=x,maxx=x,miny=y,maxy=y;
    while(st.length){ const i=st.pop(); pix.push(i); const xx=i%W, yy=(i-xx)/W; if(xx<minx)minx=xx; if(xx>maxx)maxx=xx; if(yy<miny)miny=yy; if(yy>maxy)maxy=yy;
      for(let dy=-1;dy<=1;dy++){ const y2=yy+dy; if(y2<0||y2>=H) continue; for(let dx=-1;dx<=1;dx++){ const x2=xx+dx; if(x2<0||x2>=W) continue; const j=y2*W+x2; if(!cur[j]&&bin[j]&&!done[j]){ cur[j]=1; st.push(j); } } } }
    const inside=(a,b)=>a>=0&&b>=0&&a<W&&b<H&&cur[b*W+a]===1;
    const ol=IC.traceOutline(inside,x,y,!!o.fourConnected);
    const rect={x:minx,y:miny,w:maxx-minx+1,h:maxy-miny+1};
    let mask=null;
    let maskPix=pix;
    if(o.includeHoles){ const b=IC.componentFilledBox(cur,W,H,x,y,v=>v===1); maskPix=[]; for(let yy=0;yy<b.h;yy++) for(let xx=0;xx<b.w;xx++) if(b.m[xx+yy*b.w]) maskPix.push((b.y0+yy)*W+b.x0+xx); }
    if(rect.w>1&&rect.h>1){ mask=new Uint8Array(rect.w*rect.h); for(const i of maskPix){ const xx=i%W, yy=(i-xx)/W; mask[(yy-miny)*rect.w+(xx-minx)]=1; } }
    for(const i of pix) cur[i]=0;
    for(const i of maskPix) done[i]=1;
    const np=mask?maskPix.length:rect.w*rect.h;
    let include=true;
    if(o.excludeEdges&&(rect.x===0||rect.y===0||rect.x+rect.w===W||rect.y+rect.h===H)) include=false;
    if(minC>0||maxC!==1){ const per=IC.tracedPerimeter(ol.xs,ol.ys,1,1); let c=per===0?0:4*Math.PI*(np/(per*per)); if(c>1&&maxC<=1) c=1; if(c<minC||c>maxC) include=false; }
    if(np>=minS&&np<=maxS&&include){
      count++;
      const roi={mask,rect,poly:ol,kind:"traced"};
      if(!intensity) intensity={w:W,h:H,bit:8,data:Float32Array.from(bin,v=>v?255:0)}; // 한 번만 만듦
      const meas=IC.measureRoi(intensity,roi,cal);
      for(const i of maskPix) labels[i]=count;
      parts.push({n:count,rect,poly:ol,mask,meas,start:[x,y]});
    }
  }
  if((y&31)===0) IC.progress(y/H); }
  let tot=0; for(const q of parts) tot+=q.meas.Area;
  return {particles:parts,count,totalArea:tot,labels};
};

/* 크기 분포 요약: 값 배열 → n, mean, sd, median, D10/D50/D90, 로그정규 MLE (μ, σ of ln), 기하평균·기하표준편차 */
IC.sizeSummary=function(vals){
  // 선형 통계 = 유한한 모든 값, 로그정규 = 양수만 (nPos 로 따로 알림). 백분위(D10…)는 개수 기준
  const v=vals.filter(x=>x!=null&&isFinite(x)).sort((a,b)=>a-b); const n=v.length; if(!n) return {n:0};
  const mean=v.reduce((a,b)=>a+b,0)/n; const sd=n>1?Math.sqrt(v.reduce((a,b)=>a+(b-mean)*(b-mean),0)/(n-1)):NaN; /* n=1 이면 SD·SE·CV 정의 없음 */
  const q=pct=>{ const h=(n-1)*pct; const lo=Math.floor(h), hi=Math.ceil(h); return v[lo]+(v[hi]-v[lo])*(h-lo); };
  const pos=v.filter(x=>x>0); const np=pos.length; const ln=pos.map(Math.log); const mu=np?ln.reduce((a,b)=>a+b,0)/np:NaN; const sg=np>1?Math.sqrt(ln.reduce((a,b)=>a+(b-mu)*(b-mu),0)/(np-1)):NaN;
  return {n,nPos:np,mean,sd,se:sd/Math.sqrt(n),cv:mean?sd/mean:NaN,min:v[0],max:v[n-1],median:q(0.5),d10:q(0.1),d50:q(0.5),d90:q(0.9),
    span:(q(0.9)-q(0.1))/q(0.5),lnMu:mu,lnSigma:sg,gMean:Math.exp(mu),gSD:Math.exp(sg),lnMuSE:sg/Math.sqrt(np),
    lnMode:Math.exp(mu-sg*sg),lnMean:Math.exp(mu+sg*sg/2)};
};
/* 히스토그램 칸 (Freedman–Diaconis, 최소 5 최대 60) */
IC.binsFD=function(vals,nb){ const v=vals.filter(isFinite).sort((a,b)=>a-b); const n=v.length; if(!n) return null;
  const lo=v[0], hi=v[n-1]; if(!nb){ const iqr=v[Math.floor(0.75*(n-1))]-v[Math.floor(0.25*(n-1))]; const w=2*iqr/Math.cbrt(n); nb=w>0?Math.ceil((hi-lo)/w):Math.ceil(Math.sqrt(n)); nb=Math.max(5,Math.min(60,nb)); }
  const w=(hi-lo)/nb||1; const c=new Array(nb).fill(0); for(const x of v){ let k=Math.floor((x-lo)/w); if(k>=nb) k=nb-1; c[k]++; }
  return {lo,hi,w,counts:c,n}; };

/* 다각형(꼭짓점, 화소 모서리 좌표) → 마스크 — ij/process/PolygonFiller.fillByteProcessorMask 이식.
   정수 꼭짓점이면 int 경로, 아니면 float 경로(꼭짓점 y 반올림). 마스크 범위 = 경계 상자 [floor(min), ceil(max)) */
IC.polygonMask=function(xs,ys,W,H){
  xs=Array.from(xs,Math.fround); ys=Array.from(ys,Math.fround); // ImageJ FloatPolygon 은 float32 좌표
  const n=xs.length; let minx=Infinity,maxx=-Infinity,miny=Infinity,maxy=-Infinity;
  for(let i=0;i<n;i++){ minx=Math.min(minx,xs[i]); maxx=Math.max(maxx,xs[i]); miny=Math.min(miny,ys[i]); maxy=Math.max(maxy,ys[i]); }
  const bx=Math.floor(minx), by=Math.floor(miny), bw=Math.ceil(maxx)-bx, bh=Math.ceil(maxy)-by;
  const isInt=xs.every(Number.isInteger)&&ys.every(Number.isInteger);
  const ex=[], ey1=[], ey2=[], esl=[]; let yMin=Infinity, yMax=-Infinity;
  for(let i=0;i<n;i++){ const j=i===n-1?0:i+1;
    if(isInt){ let y1=ys[i]-by, y2=ys[j]-by, x1=xs[i]-bx, x2=xs[j]-bx; if(y1===y2) continue; if(y1>y2){ [y1,y2]=[y2,y1]; [x1,x2]=[x2,x1]; }
      const sl=(x2-x1)/(y2-y1); ex.push(x1+0.5*sl+1e-8); ey1.push(y1); ey2.push(y2); esl.push(sl); if(y1<yMin) yMin=y1; if(y2>yMax) yMax=y2; }
    else { let y1f=ys[i]-by, y2f=ys[j]-by, x1f=xs[i]-bx, x2f=xs[j]-bx; let y1=Math.round(y1f), y2=Math.round(y2f);
      if(y1===y2||(y1<=0&&y2<=0)) continue; if(y1>y2){ [y1,y2]=[y2,y1]; [y1f,y2f]=[y2f,y1f]; [x1f,x2f]=[x2f,x1f]; }
      const sl=(x2f-x1f)/(y2f-y1f); ex.push(x1f+(y1-y1f+0.5)*sl+1e-8); ey1.push(y1); ey2.push(y2); esl.push(sl); if(y1<yMin) yMin=y1; if(y2>yMax) yMax=y2; } }
  const E=ex.length, mask=new Uint8Array(Math.max(0,bw*bh)); let act=[];
  const sortAct=()=>{ for(let i=0;i<act.length;i++){ let m=i; for(let k=i;k<act.length;k++) if(ex[act[k]]<ex[act[m]]) m=k; const t=act[m]; act[m]=act[i]; act[i]=t; } };
  const yStart=yMin>0?yMin:0;
  if(yMin!==0){ for(let i=0;i<E;i++) if(ey1[i]<yStart&&ey2[i]>=yStart){ ex[i]+=esl[i]*(yStart-ey1[i]); act.push(i); } sortAct(); }
  for(let y=yStart;y<Math.min(bh,yMax+1);y++){
    act=act.filter(i=>!(y<ey1[i]||y>=ey2[i]));
    for(let e=0;e<E;e++) if(y===ey1[e]){ let idx=0; while(idx<act.length&&ex[e]>ex[act[idx]]) idx++; act.splice(idx,0,e); }
    for(let i=0;i+1<act.length;i+=2){ let x1=Math.trunc(ex[act[i]]+0.5), x2=Math.trunc(ex[act[i+1]]+0.5); if(x1<0) x1=0; if(x1>bw) x1=bw; if(x2<0) x2=0; if(x2>bw) x2=bw; for(let x=x1;x<x2;x++) mask[y*bw+x]=1; }
    let prev=-Number.MAX_VALUE, sorted=true; for(const i of act){ ex[i]+=esl[i]; if(ex[i]<prev) sorted=false; prev=ex[i]; } if(!sorted) sortAct();
  }
  // 영상 밖으로 나간 부분 자르기
  const rx=Math.max(0,bx), ry=Math.max(0,by), rx2=Math.min(W,bx+bw), ry2=Math.min(H,by+bh), rw=Math.max(0,rx2-rx), rh=Math.max(0,ry2-ry);
  if(rx===bx&&ry===by&&rw===bw&&rh===bh) return {mask,rect:{x:bx,y:by,w:bw,h:bh}};
  const m2=new Uint8Array(rw*rh); for(let y=0;y<rh;y++) for(let x=0;x<rw;x++) m2[y*rw+x]=mask[(y+ry-by)*bw+(x+rx-bx)];
  return {mask:m2,rect:{x:rx,y:ry,w:rw,h:rh}};
};
/* 타원 마스크: 화소 중심이 타원 안 */
IC.ovalMask=function(x,y,w,h,W,H){
  const rx=Math.max(0,Math.floor(x)), ry=Math.max(0,Math.floor(y)), rx2=Math.min(W,Math.ceil(x+w)), ry2=Math.min(H,Math.ceil(y+h)); const rw=rx2-rx, rh=ry2-ry;
  const m=new Uint8Array(Math.max(0,rw*rh)); const cx=x+w/2, cy=y+h/2, a=w/2, b=h/2;
  for(let j=0;j<rh;j++) for(let i=0;i<rw;i++){ const dx=(rx+i+0.5-cx)/a, dy=(ry+j+0.5-cy)/b; if(dx*dx+dy*dy<=1) m[j*rw+i]=1; }
  return {mask:m,rect:{x:rx,y:ry,w:rw,h:rh}};
};

/* 선 프로파일: (x1,y1)-(x2,y2) 화소 좌표(화소 중심 = 정수+0.5). 쌍선형 보간, 간격 1 px, width>1 이면 수직 방향 평균 */
IC.bilinear=function(p,x,y){ // x,y: 화소 중심이 정수+0.5 인 좌표
  const W=p.w,H=p.h,d=p.data; let fx=x-0.5, fy=y-0.5; if(fx<0) fx=0; if(fy<0) fy=0; if(fx>W-1) fx=W-1; if(fy>H-1) fy=H-1;
  const x0=Math.floor(fx), y0=Math.floor(fy), x1=Math.min(W-1,x0+1), y1=Math.min(H-1,y0+1), ax=fx-x0, ay=fy-y0;
  return (d[y0*W+x0]*(1-ax)+d[y0*W+x1]*ax)*(1-ay)+(d[y1*W+x0]*(1-ax)+d[y1*W+x1]*ax)*ay;
};
IC.lineProfile=function(p,x1,y1,x2,y2,width){
  const dx=x2-x1, dy=y2-y1, L=Math.hypot(dx,dy); const n=Math.round(L); const N=n+1; const out=new Float64Array(N);
  const ux=L?dx/L:0, uy=L?dy/L:0, nx=-uy, ny=ux; const wd=Math.max(1,Math.round(width||1));
  for(let i=0;i<N;i++){ const t=n?i/n:0; const cx=x1+dx*t, cy=y1+dy*t; let s=0; for(let k=0;k<wd;k++){ const o=k-(wd-1)/2; s+=IC.bilinear(p,cx+nx*o,cy+ny*o); } out[i]=s/wd; }
  return {values:out,step:n?L/n:0,length:L};
};

})();
