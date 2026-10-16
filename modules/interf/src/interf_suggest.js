/* INNO Analytics — 간섭 검사: 원소별 분석 라인 자동 제안
   목표: 선택한 원소마다 intensity 는 크고 interference 는 작은 라인(= 적분 window)을 고른다.
   판정량: 그 window 에서 원소 신호 S 를 얻을 때의 예상 상대오차
     err = sqrt( S + ΣB_int + B_bg + Σ(δ_f·B_f)² + (δ_bg·B_bg)² ) / S
     - 계수 통계(Poisson) + 겹친 선을 빼낼 때의 계통오차(line ratio 불확도 δ) + 배경 차감 오차
     - interference 가 없으면 √S/S — intensity 가 클수록 작다. 겹치면 B 가 커져 커진다.
   intensity 출처:
     1) 스펙트럼 파일이 열려 있으면 실측 — 선택 원소의 계열별 Gaussian 모형을 NNLS 로 맞춘 계수(counts)
     2) 없으면 모형 — 조성 × 이온화 단면적(Bethe, ln U / E_c) × 형광수율·branching(xraydb) × 검출효율(windowless Si)
        계열 사이(K/L/M) 비교는 모형 오차가 커서, 차이가 3배 안쪽이면 '실측 확인'으로 표시한다. */
(function(){
"use strict";
const YT={"Be":{"K":[111.5,2.2e-05,"K","Ka1"]},"B":{"K":[188.0,0.000466,"K","Ka1"]},"C":{"K":[284.2,0.000932,"K","Ka1"]},"N":{"K":[409.9,0.002064,"K","Ka1"]},"O":{"K":[543.1,0.003861,"K","Ka1"]},"F":{"K":[696.7,0.006124,"K","Ka1"]},"Ne":{"K":[870.2,0.010651,"K","Ka1"]},"Na":{"K":[1070.8,0.013984,"K","Ka1"]},"Mg":{"K":[1303.0,0.015203,"K","Ka1"],"L":[88.6,1.2e-05,"L1","Lb3"]},"Al":{"K":[1559.0,0.019502,"K","Ka1"],"L":[117.8,1.5e-05,"L1","Lb3"]},"Si":{"K":[1839.0,0.025347,"K","Ka1"],"L":[99.2,0.00038,"L3","Ll"]},"P":{"K":[2145.5,0.03273,"K","Ka1"],"L":[135.0,3.1e-05,"L3","Ll"]},"S":{"K":[2472.0,0.041648,"K","Ka1"],"L":[162.5,0.00026,"L3","Ll"]},"Cl":{"K":[2822.0,0.052093,"K","Ka1"],"L":[200.0,0.00024,"L3","Ll"]},"Ar":{"K":[3205.9,0.06406,"K","Ka1"],"L":[248.4,0.00022,"L3","Ll"]},"K":{"K":[3608.4,0.077545,"K","Ka1"],"L":[294.6,0.00027,"L3","Ll"]},"Ca":{"K":[4038.5,0.086083,"K","Ka1"],"L":[346.2,0.00033,"L3","Ll"]},"Sc":{"K":[4492.0,0.107224,"K","Ka1"],"L":[398.7,0.00084,"L3","Ll"]},"Ti":{"K":[4966.0,0.127681,"K","Ka1"],"L":[453.8,0.001217,"L3","La1"]},"V":{"K":[5465.0,0.147563,"K","Ka1"],"L":[512.1,0.00211,"L3","La1"]},"Cr":{"K":[5989.0,0.166852,"K","Ka1"],"L":[574.1,0.003003,"L3","La1"]},"Mn":{"K":[6539.0,0.185553,"K","Ka1"],"L":[638.7,0.004058,"L3","La1"]},"Fe":{"K":[7112.0,0.203669,"K","Ka1"],"L":[706.8,0.005114,"L3","La1"]},"Co":{"K":[7709.0,0.221205,"K","Ka1"],"L":[778.1,0.006345,"L3","La1"]},"Ni":{"K":[8333.0,0.238166,"K","Ka1"],"L":[852.7,0.00775,"L3","La1"]},"Cu":{"K":[8979.0,0.254557,"K","Ka1"],"L":[932.7,0.009243,"L3","La1"]},"Zn":{"K":[9659.0,0.270384,"K","Ka1"],"L":[1021.8,0.010144,"L3","La1"]},"Ga":{"K":[10367.0,0.28565,"K","Ka1"],"L":[1116.4,0.011037,"L3","La1"]},"Ge":{"K":[11103.0,0.300363,"K","Ka1"],"L":[1217.0,0.012776,"L3","La1"]},"As":{"K":[11867.0,0.314077,"K","Ka1"],"L":[1323.6,0.01366,"L3","La1"]},"Se":{"K":[12658.0,0.326655,"K","Ka1"],"L":[1433.9,0.015399,"L3","La1"]},"Br":{"K":[13474.0,0.338624,"K","Ka1"],"L":[1550.0,0.017143,"L3","La1"]},"Kr":{"K":[14326.0,0.35,"K","Ka1"],"L":[1678.4,0.018892,"L3","La1"],"M":[95.0,0.0027,"M4","Mz"]},"Rb":{"K":[15200.0,0.360794,"K","Ka1"],"L":[1804.0,0.020645,"L3","La1"],"M":[113.0,0.0027,"M4","Mz"]},"Sr":{"K":[16105.0,0.371017,"K","Ka1"],"L":[1940.0,0.0224,"L3","La1"],"M":[136.0,0.0027,"M4","Mz"]},"Y":{"K":[17038.0,0.381106,"K","Ka1"],"L":[2080.0,0.024146,"L3","La1"],"M":[157.7,0.0027,"M4","Mz"]},"Zr":{"K":[17998.0,0.390667,"K","Ka1"],"L":[2223.0,0.026732,"L3","La1"],"M":[181.1,0.0027,"M4","Mz"]},"Nb":{"K":[18986.0,0.40004,"K","Ka1"],"L":[2371.0,0.029332,"L3","La1"],"M":[205.0,0.00275,"M4","Mz"]},"Mo":{"K":[20000.0,0.408918,"K","Ka1"],"L":[2520.0,0.031933,"L3","La1"],"M":[231.1,0.0028,"M4","Mz"]},"Tc":{"K":[21044.0,0.41731,"K","Ka1"],"L":[2677.0,0.034534,"L3","La1"],"M":[257.6,0.00285,"M4","Mz"]},"Ru":{"K":[22117.0,0.425224,"K","Ka1"],"L":[2838.0,0.037134,"L3","La1"],"M":[284.2,0.0029,"M4","Mz"]},"Rh":{"K":[23220.0,0.432648,"K","Ka1"],"L":[3004.0,0.036469,"L3","La1"],"M":[311.9,0.002833,"M4","Mz"]},"Pd":{"K":[24350.0,0.439753,"K","Ka1"],"L":[3173.0,0.038405,"L3","La1"],"M":[340.5,0.002767,"M4","Mz"]},"Ag":{"K":[25514.0,0.44617,"K","Ka1"],"L":[3351.0,0.040334,"L3","La1"],"M":[374.0,0.0027,"M4","Mz"]},"Cd":{"K":[26711.0,0.452071,"K","Ka1"],"L":[3538.0,0.04303,"L3","La1"],"M":[411.9,0.0027,"M4","Mz"]},"In":{"K":[27940.0,0.457999,"K","Ka1"],"L":[3730.0,0.045716,"L3","La1"],"M":[451.4,0.0027,"M4","Mz"]},"Sn":{"K":[29200.0,0.463513,"K","Ka1"],"L":[3929.0,0.048395,"L3","La1"],"M":[493.2,0.0027,"M4","Mz"]},"Sb":{"K":[30491.0,0.468496,"K","Ka1"],"L":[4132.0,0.051826,"L3","La1"],"M":[528.2,0.0,"M5","Ma"]},"Te":{"K":[31814.0,0.473082,"K","Ka1"],"L":[4341.0,0.055251,"L3","La1"],"M":[573.0,0.0,"M5","Ma"]},"I":{"K":[33169.0,0.476899,"K","Ka1"],"L":[4557.0,0.058678,"L3","La1"],"M":[619.3,0.0,"M5","Ma"]},"Xe":{"K":[34561.0,0.480336,"K","Ka1"],"L":[4786.0,0.06285,"L3","La1"],"M":[676.4,0.0,"M5","Ma"]},"Cs":{"K":[35985.0,0.483405,"K","Ka1"],"L":[5012.0,0.067027,"L3","La1"],"M":[726.6,0.0,"M5","Ma"]},"Ba":{"K":[37441.0,0.486117,"K","Ka1"],"L":[5247.0,0.071214,"L3","La1"],"M":[780.5,0.0,"M5","Ma"]},"La":{"K":[38925.0,0.488483,"K","Ka1"],"L":[5483.0,0.07615,"L3","La1"],"M":[836.0,0.0,"M5","Ma"]},"Ce":{"K":[40443.0,0.490514,"K","Ka1"],"L":[5723.0,0.081102,"L3","La1"],"M":[883.8,0.001067,"M5","Ma"]},"Pr":{"K":[41991.0,0.492095,"K","Ka1"],"L":[5964.0,0.086077,"L3","La1"],"M":[928.8,0.002133,"M5","Ma"]},"Nd":{"K":[43569.0,0.493619,"K","Ka1"],"L":[6208.0,0.091093,"L3","La1"],"M":[980.4,0.0032,"M5","Ma"]},"Pm":{"K":[45184.0,0.495107,"K","Ka1"],"L":[6459.0,0.096125,"L3","La1"],"M":[1027.0,0.0041,"M5","Ma"]},"Sm":{"K":[46834.0,0.274111,"K","Ka2"],"L":[6716.0,0.101191,"L3","La1"],"M":[1083.4,0.005,"M5","Ma"]},"Eu":{"L":[6977.0,0.107022,"L3","La1"],"M":[1127.5,0.0059,"M5","Ma"]},"Gd":{"L":[7243.0,0.112892,"L3","La1"],"M":[1189.6,0.007075,"M5","Ma"]},"Tb":{"L":[7514.0,0.119533,"L3","La1"],"M":[1241.1,0.00825,"M5","Ma"]},"Dy":{"L":[7790.0,0.126948,"L3","La1"],"M":[1292.0,0.009425,"M5","Ma"]},"Ho":{"L":[8071.0,0.13295,"L3","La1"],"M":[1351.0,0.0106,"M5","Ma"]},"Er":{"L":[8358.0,0.14046,"L3","La1"],"M":[1409.0,0.012033,"M5","Ma"]},"Tm":{"L":[8648.0,0.147288,"L3","La1"],"M":[1468.0,0.013467,"M5","Ma"]},"Yb":{"L":[8944.0,0.159957,"L3","La1"],"M":[1528.0,0.0149,"M5","Ma"]},"Lu":{"L":[9244.0,0.160135,"L3","La1"],"M":[1589.0,0.016767,"M5","Ma"]},"Hf":{"L":[9561.0,0.166994,"L3","La1"],"M":[1662.0,0.018633,"M5","Ma"]},"Ta":{"L":[9881.0,0.174646,"L3","La1"],"M":[1735.0,0.0205,"M5","Ma"]},"W":{"L":[10207.0,0.182358,"L3","La1"],"M":[1809.0,0.0214,"M5","Ma"]},"Re":{"L":[10535.0,0.190837,"L3","La1"],"M":[1883.0,0.0223,"M5","Ma"]},"Os":{"L":[10871.0,0.199357,"L3","La1"],"M":[1960.0,0.0232,"M5","Ma"]},"Ir":{"L":[11215.0,0.20791,"L3","La1"],"M":[2040.0,0.024,"M5","Ma"]},"Pt":{"L":[11564.0,0.215783,"L3","La1"],"M":[2122.0,0.0248,"M5","Ma"]},"Au":{"L":[11919.0,0.221468,"L3","La1"],"M":[2206.0,0.0256,"M5","Ma"]},"Hg":{"L":[12284.0,0.22947,"L3","La1"],"M":[2295.0,0.027325,"M5","Ma"]},"Tl":{"L":[12658.0,0.238135,"L3","La1"],"M":[2389.0,0.02905,"M5","Ma"]},"Pb":{"L":[13035.0,0.246076,"L3","La1"],"M":[2484.0,0.030775,"M5","Ma"]},"Bi":{"L":[13419.0,0.253974,"L3","La1"],"M":[2580.0,0.0325,"M5","Ma"]},"Po":{"L":[13814.0,0.261825,"L3","La1"],"M":[2683.0,0.033733,"M5","Ma"]},"At":{"L":[14214.0,0.269621,"L3","La1"],"M":[2787.0,0.034967,"M5","Ma"]},"Rn":{"L":[14619.0,0.276688,"L3","La1"],"M":[2892.0,0.0362,"M5","Ma"]},"Fr":{"L":[15031.0,0.284372,"L3","La1"],"M":[3000.0,0.039575,"M5","Ma"]},"Ra":{"L":[15444.0,0.291996,"L3","La1"],"M":[3105.0,0.04295,"M5","Ma"]},"Ac":{"L":[15871.0,0.299567,"L3","La1"],"M":[3219.0,0.046325,"M5","Ma"]},"Th":{"L":[16300.0,0.307091,"L3","La1"],"M":[3332.0,0.0497,"M5","Ma"]},"Pa":{"L":[16733.0,0.314576,"L3","La1"],"M":[3442.0,0.053075,"M5","Ma"]},"U":{"L":[17166.0,0.32204,"L3","La1"],"M":[3552.0,0.05645,"M5","Ma"]}};   /* el → fam → [E_edge eV, ω·branching(주선), 초기준위, 주선] — xraydb 4.5.8 (Elam 2002) */
const NS={K:2,L1:2,L2:2,L3:4,M1:2,M2:2,M3:4,M4:4,M5:6};
const DELTA={K:0.03,L:0.05,M:0.15};        /* 겹친 선을 빼낼 때 쓰는 line ratio 의 상대 불확도 (M 은 DB 자체가 부정확) */
const DBG=0.03;                            /* 배경 차감 상대 불확도 */
const NOM=1e5;                             /* 모형 모드: 가장 밝은 선 = 1e5 counts 로 둔 명목 척도 */
const WIN=0.6;                             /* 적분 window 반폭 = 0.6·FWHM(E) (±1.41σ, 선의 84 %) */

function erf(x){ const s=x<0?-1:1; x=Math.abs(x); const t=1/(1+0.3275911*x);
  const y=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-x*x); return s*y; }
const gfrac=(mu,sig,a,b)=>0.5*(erf((b-mu)/(sig*Math.SQRT2))-erf((a-mu)/(sig*Math.SQRT2)));
function siMu(EeV){ const t=(typeof MACDB!=="undefined")&&MACDB.mac&&MACDB.mac.Si; if(!t) return null;
  const X=t[0],Y=t[1]; if(EeV<=X[0]) return Y[0]; if(EeV>=X[X.length-1]) return Y[Y.length-1];
  let i=1; while(X[i]<EeV) i++; const f=Math.log(EeV/X[i-1])/Math.log(X[i]/X[i-1]); return Math.exp(Math.log(Y[i-1])+f*Math.log(Y[i]/Y[i-1])); }
/* windowless SDD: Si dead layer 0.1 µm 투과 × Si 450 µm 흡수 */
function detEff(EkeV){ const m=siMu(EkeV*1000); if(m==null) return 1; const mr=m*2.33; return Math.exp(-mr*1e-5)*(1-Math.exp(-mr*0.045)); }
function famYield(el,fam,kV){ const t=YT[el]&&YT[el][fam]; if(!t) return 0; const U=kV*1000/t[0]; if(!(U>1)) return 0; return (NS[t[2]]||2)*Math.log(U)/t[0]*t[1]; }
const sigAt=(EkeV,fw)=>fwhmAt(EkeV*1000,fw)/2355;      /* keV */

function linesOf(el,Emin,Emax){ const r=XRL[el], out=[]; if(!r) return out;
  for(const fam in r.lines) for(const [nm,e,rel] of r.lines[fam]) if(e>=Emin&&e<=Emax&&rel>=0.01) out.push({el,fam,nm,e,rel});
  return out; }
function parseComp(txt,els){ const c={}; els.forEach(e=>c[e]=1); let any=false;
  (txt||"").split(/[,;]+/).forEach(p=>{ const m=p.trim().match(/^([A-Za-z]{1,2})\s*[:=]?\s*([0-9.eE+-]+)$/); if(!m) return;
    const el=m[1][0].toUpperCase()+m[1].slice(1).toLowerCase(); const v=+m[2]; if(el in c&&v>=0&&isFinite(v)){ c[el]=v; any=true; } });
  return {c,any}; }

/* 작은 NNLS (Lawson–Hanson, 정규방정식) */
function nnls(G,g){ const n=g.length, x=new Float64Array(n), P=new Array(n).fill(false);
  const solveP=()=>{ const idx=[...Array(n).keys()].filter(j=>P[j]); const m=idx.length; const A=idx.map(i=>idx.map(j=>G[i][j])), b=idx.map(i=>g[i]);
    for(let k=0;k<m;k++){ let p=k; for(let r=k+1;r<m;r++) if(Math.abs(A[r][k])>Math.abs(A[p][k])) p=r; [A[k],A[p]]=[A[p],A[k]]; [b[k],b[p]]=[b[p],b[k]];
      const d=A[k][k]||1e-300; for(let r=k+1;r<m;r++){ const f=A[r][k]/d; for(let c=k;c<m;c++) A[r][c]-=f*A[k][c]; b[r]-=f*b[k]; } }
    const z=new Array(m).fill(0); for(let k=m-1;k>=0;k--){ let s=b[k]; for(let c=k+1;c<m;c++) s-=A[k][c]*z[c]; z[k]=s/(A[k][k]||1e-300); }
    const out=new Float64Array(n); idx.forEach((j,k)=>out[j]=z[k]); return out; };
  const grad=()=>{ const w=new Float64Array(n); for(let i=0;i<n;i++){ let s=g[i]; for(let j=0;j<n;j++) s-=G[i][j]*x[j]; w[i]=s; } return w; };
  const scale=Math.max(1e-300,...g.map(Math.abs));
  for(let it=0;it<3*n+10;it++){ const w=grad(); let jm=-1,wm=1e-10*scale; for(let j=0;j<n;j++) if(!P[j]&&w[j]>wm){ wm=w[j]; jm=j; } if(jm<0) break; P[jm]=true;
    for(let inner=0;inner<3*n+10;inner++){ const z=solveP(); if([...Array(n).keys()].every(j=>!P[j]||z[j]>0)){ x.set(z); break; }
      let a=1; for(let j=0;j<n;j++) if(P[j]&&z[j]<=0) a=Math.min(a,x[j]/(x[j]-z[j])); for(let j=0;j<n;j++){ x[j]+=a*(z[j]-x[j]); if(P[j]&&x[j]<=1e-12*scale){ P[j]=false; x[j]=0; } } } }
  return x; }

/* 쓰는 계열: 본체 규칙(famAllowed)을 바탕으로 한다 — K 외 계열은 주선 ≥ 0.8 keV, M 은 Z≥57.
   그 아래 계열은 형광수율이 사실상 0이라 실측 봉우리를 만들지 못하고 모형 오차만 흡수한다(본체 주석의 Ar L·V L 사례). */
function famOK(el,fam,Emin,Emax){ const r=XRL[el]; if(!r||!r.lines[fam]) return false;
  const pr=r.lines[fam].reduce((m,l)=>l[2]>m[2]?l:m,r.lines[fam][0]); if(pr[1]<Emin||pr[1]>Emax) return false;
  /* 사용자가 직접 고른 원소이므로 본체 자동 식별보다 한 단계 넓힌다: L 은 Z≥20 이고 주선 ≥ 0.8 keV 이면 쓴다.
     실측 파일(Ni 그리드)에서 Ni Lα(0.851)를 빼면 0.8–1.2 keV 에 설명되지 않는 counts 가 남았다. */
  const z=r.z||0; if(fam==="K") return true; if(pr[1]<0.8) return false; return fam==="L"?z>=20:z>=57; }

/* 실측: 본체와 같은 선 모양(elementProfile — Gaussian + 저에너지 불완전 전하수집 꼬리)으로
   선택 원소의 계열별 열을 만들어, SNIP 배경을 뺀 스펙트럼에 NNLS 로 맞춘다.
   꼬리 크기 tau 는 residual에서 재고(본체 estimateTau 와 같은 방식) 한 번 더 푼다. */
function fitMeasured(S,els,fw,Emin,Emax){
  const snip=specBackground(S.y,S,fw);
  const ch=[]; for(let i=0;i<S.nch;i++) if(S.E[i]>=Emin&&S.E[i]<=Emax&&i>=(S.i0|0)) ch.push(i);
  const cols=[]; els.forEach(el=>{ for(const fam in (XRL[el]||{lines:{}}).lines){ if(!famOK(el,fam,Emin,Emax)) continue;
    const pr=elementProfile(S,el,fam,fw,[Emin,Emax]); if(pr) cols.push({el,fam,pr}); } });
  const nEl=cols.length;
  /* SNIP 은 연속선을 낮게 잡는 경향이 있다(실측 파일: 9–11 keV 에서 0, 6–7 keV 에서 필요량의 2 %).
     남은 연속선을 원소 열이 흡수하지 않도록 넓은 삼각(hat) 기저를 음이 아닌 배경 보정 열로 함께 맞춘다.
     마디 간격 ≥ 4·FWHM 이라 피크와 구별된다. */
  const kn=[]; { let e=Emin; while(e<Emax){ kn.push(e); e+=Math.max(0.2,4*fwhmAt(e*1000,fw)/1000,0.15*e); } kn.push(Emax); }
  const hats=kn.map((k,j)=>{ const lo=kn[j-1]!=null?kn[j-1]:k, hi=kn[j+1]!=null?kn[j+1]:k;
    return ch.map(i=>{ const e=S.E[i]; if(e<=k) return lo===k?1:(e<=lo?0:(e-lo)/(k-lo)); return hi===k?1:(e>=hi?0:(hi-e)/(hi-k)); }); });
  const w=ch.map(i=>1/Math.max(1,S.y[i])), net=ch.map(i=>S.y[i]-snip[i]);
  const solve=tau=>{ const X=cols.map(c=>ch.map(i=>c.pr.prof[i]+(tau&&c.pr.dtail?tau*c.pr.dtail[i]:0))).concat(hats);
    const n=X.length, G=[...Array(n)].map(()=>new Array(n).fill(0)), g=new Array(n).fill(0);
    for(let a=0;a<n;a++){ const xa=X[a]; for(let b=a;b<n;b++){ const xb=X[b]; let t=0; for(let k=0;k<ch.length;k++) t+=xa[k]*xb[k]*w[k]; G[a][b]=G[b][a]=t; }
      let t=0; for(let k=0;k<ch.length;k++) t+=xa[k]*w[k]*net[k]; g[a]=t; }
    return {X,A:nnls(G,g)}; };
  let tau=0, R=solve(0);
  for(let pass=0;pass<2;pass++){      /* residual이 요구하는 꼬리 증분 */
    const m=ch.map((i,k)=>cols.reduce((s,c,a)=>s+(c.pr.dtail?R.A[a]*c.pr.dtail[i]:0),0));
    const f=ch.map((i,k)=>R.X.reduce((s,x,a)=>s+R.A[a]*x[k],0));
    let num=0,den=0; ch.forEach((i,k)=>{ num+=w[k]*(net[k]-f[k])*m[k]; den+=w[k]*m[k]*m[k]; });
    if(!(den>0)||num*num/den<9) break;
    const nt=Math.max(0,Math.min(2,tau+num/den)); if(Math.abs(nt-tau)<0.02) break; tau=nt; R=solve(tau); }
  const bg=Float64Array.from(snip), comps=[], model=new Float64Array(S.nch);
  cols.forEach((c,a)=>{ const arr=new Float64Array(S.nch); ch.forEach((i,k)=>{ arr[i]=R.A[a]*R.X[a][k]; model[i]+=arr[i]; }); comps.push({el:c.el,fam:c.fam,A:R.A[a],arr}); });
  hats.forEach((hx,h)=>{ const a=nEl+h; ch.forEach((i,k)=>bg[i]+=R.A[a]*hx[k]); });
  return {bg,snip,comps,cols:comps,model,ch,tau};
}

/* 모형: 가상 채널(5 eV) 위에 계열별 Gaussian 을 그린다 */
function modelComps(els,fw,Emin,Emax,kV,comp){
  const disp=0.005, n=Math.ceil((Emax+0.6)/disp), E=Float64Array.from({length:n},(_,i)=>i*disp);
  const comps=[]; els.forEach(el=>{ for(const fam in (XRL[el]||{lines:{}}).lines){ if(!famOK(el,fam,Emin,Emax)) continue;
    const Y=comp.c[el]*famYield(el,fam,kV); if(!(Y>0)) continue; const arr=new Float64Array(n);
    XRL[el].lines[fam].forEach(([nm,e,rel])=>{ if(rel<0.01) return; const I=Y*rel*detEff(e), sg=sigAt(e,fw); const c0=Math.round(e/disp), wv=Math.ceil(6*sg/disp);
      for(let i=Math.max(0,c0-wv);i<=Math.min(n-1,c0+wv);i++){ const d=(E[i]-e)/sg; arr[i]+=I*Math.exp(-0.5*d*d)/(sg*Math.sqrt(2*Math.PI))*disp; } });
    comps.push({el,fam,A:Y,arr}); } });
  const mx=Math.max(1e-300,...comps.map(c=>Math.max(...c.arr))); const k=NOM/(mx*Math.sqrt(2*Math.PI)*sigAt(1,fw)/disp);
  comps.forEach(c=>{ for(let i=0;i<n;i++) c.arr[i]*=k; c.A*=k; });
  return {E,n,comps,disp};
}

function compute(opt){
  const {els,fw,emax,kV,compTxt,S}=opt; const Emin=Math.max(0.15,S?(S.E[S.i0|0]||0.15):0.15,S&&S.cal&&S.cal.beginE?S.cal.beginE/1000:0);
  let mode, comp=null, fit=null, grid;
  if(S){ mode="measured"; fit=fitMeasured(S,els,fw,Emin,emax); grid={E:S.E,n:S.nch,comps:fit.comps,bg:fit.bg,y:S.y}; }
  else { mode="model"; comp=parseComp(compTxt,els); const M=modelComps(els,fw,Emin,emax,kV,comp); grid={E:M.E,n:M.n,comps:M.comps,bg:null,y:null}; }
  const sumWin=(arr,a,b)=>{ let s=0; const E=grid.E; for(let i=0;i<grid.n;i++) if(E[i]>=a&&E[i]<b) s+=arr[i]; return s; };
  /* 선 목록(이름·에너지·세기) — 후보 window 를 만들고 이름 붙이는 데 쓴다 */
  const L=[]; grid.comps.forEach(c=>XRL[c.el].lines[c.fam].forEach(([nm,e,rel])=>{ if(e>=Emin-0.05&&e<=emax+0.05&&rel>=0.01) L.push({el:c.el,fam:c.fam,nm,e,rel,I:c.A*rel}); }));
  const out=els.map(el=>{
    const own=L.filter(l=>l.el===el&&l.I>0&&l.e>=Emin&&l.e<=emax).sort((x,y)=>x.e-y.e), groups=[];
    own.forEach(l=>{ const g=groups[groups.length-1]; if(g&&g.fam===l.fam&&(l.e-g.ls[g.ls.length-1].e)*1000<0.5*fwhmAt(l.e*1000,fw)) g.ls.push(l); else groups.push({fam:l.fam,ls:[l]}); });
    const cands=groups.map(g=>{ const I=g.ls.reduce((a,m)=>a+m.I,0), rel=g.ls.reduce((a,m)=>a+m.rel,0); const main=g.ls.slice().sort((x,y)=>y.I-x.I)[0];
      const major=g.ls.filter(m=>m.I>=0.1*main.I).sort((x,y)=>y.I-x.I);
      return {el,fam:g.fam,nm:major.map(m=>m.nm).join("+"),e:g.ls.reduce((a,m)=>a+m.e*m.I,0)/I,rel,I}; })
    .filter(l=>l.rel>=0.1).map(l=>{
      const h=WIN*fwhmAt(l.e*1000,fw)/1000, a=l.e-h, b=l.e+h; let Sg=0; const Bf={};
      grid.comps.forEach(c=>{ const v=sumWin(c.arr,a,b); if(!(v>0)) return; if(c.el===el) Sg+=v; else Bf[c.el+" "+c.fam]=(Bf[c.el+" "+c.fam]||0)+v; });
      const Bint=Object.values(Bf).reduce((x,y)=>x+y,0), Bbg=grid.bg?Math.max(0,sumWin(grid.bg,a,b)):0;
      let sys=0; for(const k in Bf) sys+=Math.pow(DELTA[k.split(" ")[1]]*Bf[k],2);
      const err=Sg>0?Math.sqrt(Sg+Bint+Bbg+sys+Math.pow(DBG*Bbg,2))/Sg:Infinity;
      const top=Object.entries(Bf).sort((x,y)=>y[1]-x[1])[0];
      let topLine=null; if(top){ const [te,tf]=top[0].split(" "); let best=0; L.forEach(m=>{ if(m.el===te&&m.fam===tf){ const f=m.I*gfrac(m.e,sigAt(m.e,fw),a,b); if(f>best){ best=f; topLine=m; } } }); }
      let unexpl=0; if(grid.y){ const r=sumWin(grid.y,a,b)-Bbg-Sg-Bint; if(r>3*Math.sqrt(Math.max(1,Sg+Bint+Bbg))&&r>0.1*Sg) unexpl=r; }
      return {line:l,a,b,S:Sg,Bint,Bbg,Bf,err,purity:Sg/(Sg+Bint),pb:Bbg>0?Sg/Bbg:null,top:top?{key:top[0],frac:top[1]/(Sg+Bint),line:topLine}:null,unexpl};
    }).filter(c=>c.S>0);
    cands.sort((x,y)=>x.err-y.err);
    const best=cands[0]||null; let alt=null;
    if(best){ alt=cands.find(c=>c.line.fam!==best.line.fam)||cands.find(c=>c!==best)||null; }
    const famUnsure=mode==="model"&&best&&alt&&alt.line.fam!==best.line.fam&&alt.err/best.err<3;
    const close=mode==="measured"&&best&&alt&&alt.err/best.err<1.25;
    return {el,best,alt,cands,famUnsure,close,noSignal:mode==="measured"&&!best};
  });
  return {mode,Emin,emax,fw,kV,comp,out,fit,L,S,grid};
}

/* ---------------- 화면 ---------------- */
const pct=v=>!isFinite(v)?"—":v<0.001?"<0.1 %":v<0.1?(100*v).toFixed(1)+" %":(100*v).toFixed(0)+" %";
const badge=c=>c.purity>=0.9?'<span class="badge ok">clean</span>':c.purity>=0.6?'<span class="badge warn">overlap-corrected</span>':'<span class="badge risk">needs deconvolution</span>';
const G1=n=>n.replace(/^(K|L|M)([a-z])(\d?)$/,(m,f,g,d)=>f+({a:"α",b:"β",g:"γ",l:"ℓ",n:"η"}[g]||g)+d);
const lname=l=>{ const ps=l.nm.split("+").map(G1); if(ps.length>1){ const pre=ps[0].replace(/\d+$/,""); if(ps.every(q=>q.replace(/\d+$/,"")===pre)) return `${l.el} ${pre}${ps.map(q=>q.slice(pre.length)||"").join(",")}`; } return `${l.el} ${ps.join("+")}`; };
function spark(R,c){
  const fw=R.fw, g=R.grid, l=c.line; const H=fwhmAt(l.e*1000,fw)/1000, lo=l.e-3*H, hi=l.e+3*H, W=120, Hh=40;
  const idx=[]; for(let i=0;i<g.n;i++) if(g.E[i]>=lo&&g.E[i]<=hi) idx.push(i); if(idx.length<3) return "";
  const own=idx.map(i=>g.comps.reduce((s,q)=>q.el===l.el?s+q.arr[i]:s,0)), oth=idx.map(i=>g.comps.reduce((s,q)=>q.el!==l.el?s+q.arr[i]:s,0)), bgv=idx.map(i=>g.bg?Math.max(0,g.bg[i]):0), dat=g.y?idx.map(i=>g.y[i]):null;
  const top=Math.max(1e-9,...idx.map((_,k)=>bgv[k]+oth[k]+own[k]),...(dat||[0]));
  const xs=idx.map(i=>g.E[i]), N=xs.length, X=x=>((x-lo)/(hi-lo)*W).toFixed(1), Y=v=>(Hh-2-(v/top)*(Hh-6)).toFixed(1);
  const area=(base,add)=>"M"+xs.map((x,k)=>X(x)+" "+Y(base[k]+add[k])).join("L")+"L"+xs.slice().reverse().map((x,k)=>X(x)+" "+Y(base[N-1-k])).join("L")+"Z";
  const zero=xs.map(()=>0), b2=bgv.map((v,k)=>v+oth[k]);
  return `<svg viewBox="0 0 ${W} ${Hh}" width="${W}" height="${Hh}" role="img" aria-label="${esc(lname(l))} window" style="display:block">
    <rect x="${X(c.a)}" y="0" width="${(X(c.b)-X(c.a)).toFixed(1)}" height="${Hh}" style="fill:var(--accent-soft)"/>
    ${g.bg?`<path d="${area(zero,bgv)}" style="fill:var(--grid)"/>`:""}
    <path d="${area(bgv,oth)}" style="fill:var(--warn);opacity:.55"/>
    <path d="${area(b2,own)}" style="fill:var(--accent);opacity:.75"/>
    ${dat?`<path d="M${xs.map((x,k)=>X(x)+" "+Y(dat[k])).join("L")}" style="fill:none;stroke:var(--ink-2);stroke-width:1"/>`:""}
    <line x1="${X(l.e)}" x2="${X(l.e)}" y1="0" y2="${Hh}" style="stroke:var(--accent-strong);stroke-width:1;stroke-dasharray:2 2"/></svg>`;
}
function candRow(c,R,rank){ const l=c.line; return `<tr>
  <td class="l"><b>${esc(lname(l))}</b> <span class="muted">${l.fam}</span></td><td class="num">${l.e.toFixed(3)}</td>
  <td class="num">${R.mode==="measured"?Math.round(c.S).toLocaleString():pct(c.S/Math.max(1e-300,...R.out.filter(o=>o.el===l.el).flatMap(o=>o.cands.map(q=>q.S))))}</td>
  <td class="num">${pct(c.purity)}</td><td class="num">${pct(c.err)}</td>
  <td class="l">${c.top&&c.top.frac>=0.01?`${esc(c.top.line?lname(c.top.line):c.top.key)} <span class="muted">${pct(c.top.frac)}</span>`:'<span class="muted">—</span>'}</td>
  <td>${badge(c)}</td></tr>`; }
function html(R){
  const meas=R.mode==="measured";
  const head=meas?`<b>실측 스펙트럼 기반</b> — 선택 원소의 계열별 선 모양(Gaussian + 저에너지 꼬리, 본체와 같은 모형)을 배경(SNIP + 보정) 뺀 스펙트럼에 NNLS 로 맞춘 counts 를 씁니다.${R.fit&&R.fit.tau?` 꼬리 크기 τ = ${R.fit.tau.toFixed(2)}.`:""}`
    :`<b>모형 기반</b> — 스펙트럼 파일이 없어 이론 intensity 로 계산했습니다. 조성 ${R.comp.any?"입력값":"<b>미입력(모두 같은 원자비)</b>"}, ${R.kV} kV.
      파일을 열면 실측 counts 로 자동 대체됩니다.`;
  const rows=R.out.map(o=>{
    if(!o.best) return `<tr><td class="l"><b>${esc(o.el)}</b></td><td class="l" colspan="${meas?7:6}"><span class="muted">${o.noSignal?"맞춘 신호가 없습니다 — 이 스펙트럼에서 검출되지 않거나 범위 밖입니다":"범위 안에 쓸 만한 라인이 없습니다"}</span></td></tr>`;
    const c=o.best, l=c.line, a=o.alt;
    const notes=[]; if(o.famUnsure) notes.push(`<span class="badge info" title="K/L/M 계열 사이 intensity 비는 모형 오차가 커서, 차이가 3배 안쪽이면 실측으로 확인해야 합니다">Verify series</span>`);
    if(o.close) notes.push(`<span class="badge info" title="alternative의 예상 상대오차가 25 % 안쪽으로 가깝습니다 — 둘 다 쓸 만합니다">Similar to alternative</span>`);
    if(c.err>0.5) notes.push(`<span class="badge risk" title="예상 상대오차가 50 %를 넘습니다 — 계수 시간을 늘리거나 ROI 를 넓히세요">${meas?"Low signal in this spectrum":"Low signal"}</span>`);
    if(c.unexpl) notes.push(`<span class="badge warn" title="모형(선택 원소 + 배경)이 설명하지 못한 counts 가 window 안에 3σ 이상, 신호의 10 % 이상 남았습니다. 원인 후보: 목록에 없는 원소, 선 모양 차이(특히 M 계열 상대세기), escape·sum peak">Large window residual</span>`);
    if(l.fam==="M") notes.push(`<span class="badge gray" title="M 계열 상대세기 DB 가 부정확합니다 — 정량이면 L 계열과 교차 확인하세요${o.alt&&o.alt.line.fam==="L"?` (alternative ${esc(lname(o.alt.line))})`:""}">M series</span>`);
    return `<tr><td class="l"><b>${esc(o.el)}</b></td>
      <td class="l" style="max-width:220px;white-space:normal"><b>${esc(lname(l))}</b> <span class="muted">${l.e.toFixed(3)} keV</span><div style="display:flex;flex-wrap:wrap;gap:3px;margin-top:3px">${badge(c)} ${notes.join(" ")}</div></td>
      <td class="num">${pct(c.err)}</td><td class="num">${pct(c.purity)}</td>
      ${meas?`<td class="num">${Math.round(c.S).toLocaleString()}${c.pb!=null?`<br><span class="muted">P/B ${c.pb.toFixed(c.pb<10?1:0)}</span>`:""}</td>`:""}
      <td class="l">${c.top&&c.top.frac>=0.01?`${esc(c.top.line?lname(c.top.line):c.top.key)}<br><span class="muted">${pct(c.top.frac)} of window</span>`:'<span class="muted">None</span>'}</td>
      <td class="l">${a?`${esc(lname(a.line))} <span class="muted">${a.line.e.toFixed(3)}</span><br><span class="muted">Error ${pct(a.err)} · purity ${pct(a.purity)}</span>`:'<span class="muted">—</span>'}</td>
      <td>${spark(R,c)}</td></tr>`; }).join("");
  const all=R.out.filter(o=>o.cands.length).map(o=>`<tr><th class="l" colspan="7" style="padding-top:10px">${esc(o.el)} — ${o.cands.length} candidates</th></tr>${o.cands.map((c,k)=>candRow(c,R,k)).join("")}`).join("");
  return `<div class="card" id="ifSug"><h3>Recommended Line per Element <span class="hint">예상 상대오차가 가장 작은 라인 = intensity 최대 · interference 최소의 절충</span>
      <span style="flex:1"></span><button class="btn sm" id="ifSugCsv">CSV</button></h3>
    <div class="note ${meas?"okay":"warn"}" style="margin-bottom:8px">${head}</div>
    <div class="tblwrap"><table class="tbl"><thead><tr><th class="l">Element</th><th class="l">Recommended line</th><th title="계수 통계 + 겹친 선 차감 오차 + 배경 차감 오차${meas?"":" — 모형 모드는 가장 밝은 선 1e5 counts 로 둔 명목값"}">Expected rel. error${meas?"":" (nominal)"}</th>
      <th title="window 안 counts 중 이 원소 fraction">purity</th>${meas?"<th>counts (window)</th>":""}<th class="l">Main interference</th><th class="l">Alternative</th><th class="l">window</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div class="muted" style="font-size:12px;margin-top:6px;display:flex;gap:14px;flex-wrap:wrap;align-items:center">
      <span><i style="display:inline-block;width:10px;height:10px;background:var(--accent);opacity:.75;border-radius:2px;vertical-align:-1px"></i> This element</span>
      <span><i style="display:inline-block;width:10px;height:10px;background:var(--warn);opacity:.55;border-radius:2px;vertical-align:-1px"></i> Other elements (interference)</span>
      ${meas?`<span><i style="display:inline-block;width:10px;height:10px;background:var(--grid);border-radius:2px;vertical-align:-1px"></i> Background</span><span>— measured</span>`:""}
      <span><i style="display:inline-block;width:10px;height:10px;background:var(--accent-soft);border-radius:2px;vertical-align:-1px"></i> Integration window (±${WIN}·FWHM)</span><span>x-axis = line ± 3·FWHM</span></div>
    <details class="if-more" style="margin-top:8px"><summary>All Candidates per Element · Method</summary>
      <div class="tblwrap" style="margin-top:6px"><table class="tbl"><thead><tr><th class="l">Line</th><th>E [keV]</th><th>${meas?"counts":"Relative intensity"}</th><th>purity</th><th>Expected rel. error</th><th class="l">Main interference</th><th>Verdict</th></tr></thead><tbody>${all}</tbody></table></div>
      <ul class="muted" style="margin:8px 0 0 18px;font-size:12px;line-height:1.55">
        <li>window = 라인 에너지 ± ${WIN}·FWHM(E) (선의 84 %). 그 안에서 원소 신호 S, 다른 원소 선 B<sub>int</sub>, 배경 B<sub>bg</sub> 를 셉니다.</li>
        <li>예상 상대오차 = √(S + B<sub>int</sub> + B<sub>bg</sub> + Σ(δ·B<sub>f</sub>)² + (${DBG}·B<sub>bg</sub>)²) / S. δ(겹친 선을 빼낼 때의 line ratio 불확도) = K ${DELTA.K}, L ${DELTA.L}, M ${DELTA.M}.</li>
        <li>purity ≥ 90 % clean · 60–90 % overlap 보정(다른 선으로 interferer 를 재서 빼기) · &lt; 60 % deconvolution(MLLS) 필요.</li>
        <li>${meas?"실측 모드: counts 는 이 스펙트럼의 NNLS 계수입니다. 선택하지 않은 원소(그리드 Cu, 극편 Fe·Co, Si escape 등)가 있으면 목록에 넣어야 정확합니다.":`모형 모드: intensity ∝ 조성 × n<sub>s</sub>·ln(U)/E<sub>c</sub> (Bethe) × ω·branching (xraydb) × windowless Si 검출효율, 가장 밝은 선 = ${NOM.toExponential(0)} counts 로 둔 명목 척도. 같은 계열 안의 선 비교는 믿을 만하지만, K/L/M 계열 사이 비는 수 배까지 틀릴 수 있습니다.`}</li>
        <li>escape·sum peak, 시료 자체 흡수는 넣지 않았습니다.${meas?" 약한 선과 0.6 keV 아래 선의 counts 는 배경 추정에 민감합니다(합성 시험에서 주성분 ±5 %, 2 % 미량 L 선 ±40 %).":""}</li></ul>
    </details></div>`;
}
let LAST=null;
function readOpts(){ const q=id=>document.getElementById(id); const P=state.spec;
  const els=(q("ifEls").value||"").split(/[,\s]+/).map(s=>s.trim()).filter(Boolean).map(s=>s[0].toUpperCase()+s.slice(1).toLowerCase()).filter(s=>XRL[s]);
  const fw=+q("ifFwhm").value||P.fwhmMn||130, emax=+q("ifEmax").value||20, kV=+((q("ifKV")||{}).value)||200, compTxt=(q("ifComp")||{}).value||"";
  P.ifComp=compTxt; P.ifKV=kV;
  let S=null; try{ const B=interfBasis(); if(B.fromFile){ const rec=curRec(); S=(LASTSPEC&&LASTSPEC.sp&&rec.spectra.indexOf(LASTSPEC.sp)>=0)?LASTSPEC:loadSpectrum(); } }catch(e){ S=null; }
  return {els,fw,emax,kV,compTxt,S}; }
function setup(root){ const row=root.querySelector("#ifRun"); if(!row) return; const card=row.closest(".card"); const P=state.spec;
  const d=document.createElement("div"); d.className="row"; d.style.cssText="gap:10px;align-items:flex-end;margin-top:8px";
  d.innerHTML=`<label class="f" style="flex:1;min-width:240px" title="스펙트럼 파일이 없을 때만 씁니다. 파일이 열려 있으면 실측 counts 를 씁니다">Composition (optional · relative atomic ratio, used only without a file)
      <input type="text" id="ifComp" value="${esc(P.ifComp||"")}" placeholder="예: Na 1, Yb 0.78, Tm 0.02, Lu 0.2, F 4 — 비우면 같은 비율"></label>
    <label class="f">Accel. voltage [kV]<input type="number" id="ifKV" class="mono" style="width:80px" value="${P.ifKV||200}" min="30" max="300" step="10"></label>`;
  row.closest(".row").after(d); }
function after(){ const out=document.getElementById("ifOut"); if(!out) return; const o=readOpts();
  const old=[...out.childNodes]; out.innerHTML="";
  if(!o.els.length){ old.forEach(n=>out.appendChild(n)); return; }
  let R=null; try{ R=compute(o); }catch(e){ console.error(e); }
  if(R){ out.insertAdjacentHTML("beforeend",html(R)); LAST=R;
    const b=document.getElementById("ifSugCsv"); if(b) b.onclick=()=>{ const rows=[["element","line","family","E_keV","rel_err","purity","S","B_int","B_bg","top_interferer","alt_line","alt_rel_err","mode"]];
      R.out.forEach(q=>{ const c=q.best; if(!c) return; rows.push([q.el,c.line.nm,c.line.fam,c.line.e,c.err.toFixed(4),c.purity.toFixed(4),c.S.toFixed(1),c.Bint.toFixed(1),c.Bbg.toFixed(1),c.top&&c.top.line?c.top.line.el+" "+c.top.line.nm:"",q.alt?q.alt.line.nm:"",q.alt?q.alt.err.toFixed(4):"",R.mode]); });
      downloadText(rows.map(r=>r.join(",")).join("\n"),"eds_line_recommendation.csv","text/csv"); }; }
  if(o.els.length<2) return;             /* 원소 하나면 겹침 목록이 없다 */
  const det=document.createElement("details"); det.className="card if-all";
  const h3=old.map(n=>n.querySelector&&n.querySelector("h3")).find(Boolean); const m=h3&&h3.textContent.match(/(?:겹침\s*(\d+)건|Overlaps\s*\((\d+)\))/);
  det.innerHTML=`<summary style="cursor:pointer;font-weight:600">All overlaps ${m?`(${m[1]||m[2]})`:"(none)"} <span class="muted" style="font-weight:400">— 모든 라인쌍의 ΔE/FWHM 판정</span></summary>`;
  const inner=document.createElement("div"); inner.style.marginTop="8px"; old.forEach(n=>{ if(n.classList&&n.classList.contains("card")){ n.classList.remove("card"); n.style.padding="0"; n.style.border="0"; n.style.boxShadow="none"; } inner.appendChild(n); }); det.appendChild(inner); out.appendChild(det); }
window.IF_SUG={setup,after,compute,_last:()=>LAST,_detEff:detEff,_famYield:famYield};
})();
