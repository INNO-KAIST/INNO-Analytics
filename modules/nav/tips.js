/* INNO tips (2026-10-19) — 화면에 늘어놓던 설명 문장을 '오래 대고 있으면 뜨는' 말풍선으로.
   ① [data-tip] 요소: 마우스를 TIP_DELAY(2 s) 대고 있으면 말풍선. 손가락은 0.7 s 길게 누르기. 움직여 벗어나거나 누르거나 굴리면 사라짐.
   ② 자동 옮김: 회색 설명 문장(Image Analysis 의 div.xs.muted·.im-help·긴 span, 다른 탭의 .hint·60 자 넘는 회색 문단)을 숨기고
      (입력·단추·링크·표가 든 것, 경고 .note 안, .keep 은 그대로)
      그 글을 바로 앞 요소(없으면 부모)의 data-tip 에 붙임 — 렌더할 때마다 MutationObserver 가 처리.
      남겨야 할 문장은 class "keep" (빈 화면 안내 im-empty 는 남김). 경고(.note.warn 등)는 건드리지 않음. */
(function(){
"use strict";
const TIP_DELAY=2000, TOUCH_HOLD=700;
const TIP_SCOPES=["body"]; // 앱 전체 — 아래 규칙에 맞는 '설명 문장'만
// Image Analysis(.imx): 회색 설명 div 전부 · 다른 탭: .hint 와 60 자 넘는 회색 문단(div/p.muted)만 (짧은 회색 글은 값·단위·이름이라 둠)
const SEL_BLOCK=".imx div.xs.muted:not(.im-empty):not(.keep), .imx p.xs.muted:not(.keep), .im-help:not(.keep), .hint:not(.keep), div.muted:not(.keep), p.muted:not(.keep), div.muted.xs:not(.keep), div.muted.sm:not(.keep)";
const SEL_SPAN=".imx span.xs.muted:not(.keep)";
const MIN_OUT=60; // .imx 밖은 이보다 긴 문장만
function eligible(el){ if(el.closest(".keep,table,button,label,select,summary,.inno-tip,.note,[contenteditable]")) return false; if(el.querySelector("input,button,select,textarea,a,canvas,svg,table")) return false;
  if(/empty|placeholder|drop/i.test(el.className)) return false; // 빈 화면 안내는 남김
  const sib=[...(el.parentElement?el.parentElement.children:[])].filter(c=>c!==el&&!c.classList.contains("inno-tip-moved")); if(!sib.length) return false; // 혼자 있는 글(빈 칸 안내)은 남김
  const t=(el.textContent||"").trim(); if(!t) return false; if(!el.closest(".imx")&&!el.classList.contains("hint")&&t.length<MIN_OUT) return false; if(el.closest(".imx")&&el.classList.contains("im-empty")) return false; return true; }
let box=null, timer=0, cur=null;
function ensure(){ if(box) return box; box=document.createElement("div"); box.className="inno-tip"; box.setAttribute("role","tooltip"); document.body.appendChild(box); return box; }
function show(el,x,y){ const t=el.getAttribute("data-tip"); if(!t) return; const b=ensure(); b.textContent=t; b.style.display="block";
  const vw=innerWidth, vh=innerHeight, r=b.getBoundingClientRect(); let L=x+14, T=y+18; if(L+r.width>vw-8) L=Math.max(8,vw-8-r.width); if(T+r.height>vh-8) T=Math.max(8,y-12-r.height); b.style.left=L+"px"; b.style.top=T+"px"; }
function hide(){ clearTimeout(timer); timer=0; cur=null; if(box) box.style.display="none"; }
let lx=0, ly=0;
document.addEventListener("pointermove",e=>{ lx=e.clientX; ly=e.clientY; if(e.pointerType!=="mouse") return; const el=e.target&&e.target.closest?e.target.closest("[data-tip]"):null;
  if(el!==cur){ hide(); if(el){ cur=el; timer=setTimeout(()=>{ if(cur===el) show(el,lx,ly); },TIP_DELAY); } } },{passive:true});
document.addEventListener("pointerdown",e=>{ hide(); if(e.pointerType==="mouse") return; const el=e.target&&e.target.closest?e.target.closest("[data-tip]"):null; if(!el) return; cur=el; const x=e.clientX, y=e.clientY; timer=setTimeout(()=>{ if(cur===el) show(el,x,y); },TOUCH_HOLD); },{passive:true});
for(const ev of ["pointerup","pointercancel"]) document.addEventListener(ev,e=>{ if(e.pointerType!=="mouse"&&!(box&&box.style.display==="block")) hide(); },{passive:true});
for(const ev of ["scroll","wheel","keydown","blur"]) addEventListener(ev,()=>hide(),{passive:true,capture:true});
document.addEventListener("mouseleave",hide);
function addTip(host,text){ text=(text||"").replace(/\s+/g," ").trim(); if(!host||!text) return; const old=host.getAttribute("data-tip"); if(old&&old.includes(text)) return; host.setAttribute("data-tip",old?old+"\n\n"+text:text); }
function convert(root){ for(const sc of TIP_SCOPES){ const scopes=root.matches&&root.matches(sc)?[root]:[...(root.querySelectorAll?root.querySelectorAll(sc):[])].concat(root.closest&&root.closest(sc)?[root]:[]);
    for(const s of scopes){ s.querySelectorAll(SEL_BLOCK).forEach(el=>{ if(el.dataset.tipMoved) return; if(!eligible(el)) return; el.dataset.tipMoved="1"; const host=el.previousElementSibling||el.parentElement; addTip(host,el.innerText||el.textContent); el.classList.add("inno-tip-moved"); });
      s.querySelectorAll(SEL_SPAN).forEach(el=>{ if(el.dataset.tipMoved) return; const t=(el.textContent||"").trim(); if(t.length<=28) return; el.dataset.tipMoved="1"; addTip(el.parentElement,t); el.classList.add("inno-tip-moved"); }); } } }
let pend=false; const mo=new MutationObserver(()=>{ if(pend) return; pend=true; requestAnimationFrame(()=>{ pend=false; convert(document.body); }); });
function start(){ convert(document.body); mo.observe(document.body,{childList:true,subtree:true}); }
if(document.readyState==="loading") document.addEventListener("DOMContentLoaded",start); else start();
window.INNO_TIPS={convert,addTip,hide,DELAY:TIP_DELAY};
})();
