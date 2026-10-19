/* ============================================================
   INNO_TOUCH — 공용 터치 제스처 (pointer events) · 2026-10-17
   한 구현을 모든 탭(EM Image·EDS Maps·EDS Spectrum·Lattice FFT·Image Analysis·XPS, 나중의 Diffraction)이 같이 쓴다.

   원칙
   - pointerType==="touch" 인 포인터만 다룬다. 마우스·펜은 손대지 않는다 → 데스크톱 동작은 그대로.
   - 붙인 요소(el)에 capture 단계 수신기를 단다. 이 도우미가 맡은 손가락의 사건은 stopImmediatePropagation 으로
     앱의 기존 pointer 처리기(같은 요소 또는 자식)에 가지 않게 한다. 'pass' 로 넘긴 손가락은 앱이 그대로 받는다.
   - 화면 갱신은 requestAnimationFrame 한 번에 한 번(setView). 손가락 사건마다 다시 그리지 않는다.
   - touch-action:none 은 붙인 요소에만(opts.touchAction===false 이면 앱 CSS 에 맡김). 나머지 쪽은 평소처럼 스크롤된다.

   API
   const detach = INNO_TOUCH.attach(el, {
     getView()            → 제스처 시작 때의 보기 상태(불투명 값, 앱이 복사해서 돌려준다)
     setView(v0, g)       → 시작 상태 v0 기준으로 제스처 g 를 적용(한 프레임에 한 번)
                            g = {kind:"pan"|"pinch", x0,y0, x,y, dx,dy, scale, scaleX, scaleY,
                                 clientX0,clientY0, clientX,clientY}
                            x0,y0 = 시작 중심(손가락 둘이면 가운데), x,y = 지금 중심 — el 기준 CSS px(opts.box 로 바꿀 수 있음)
                            scale = 손가락 거리 비, scaleX/scaleY = 가로/세로 간격 비(간격이 24 px 보다 좁게 시작하면 scale)
     oneFinger            "pan"(기본) | "pass" | "none" | (pointerdown 사건)=>그중 하나
                            pass: 앱의 처리기가 그 손가락을 받는다(그리기 도구). 둘째 손가락이 오면 onPassCancel 후 핀치로 넘어간다.
                            none: 한 손가락은 탭·길게 누르기만(끌어도 보기 안 움직임).
     onPassCancel()       pass 로 넘긴 손가락을 도우미가 가져갈 때(둘째 손가락·길게 누르기) 앱이 시작한 일을 되돌린다.
                            없으면 합성 pointercancel 을 el 에 보낸다.
     onTap(p)             p = {x,y,clientX,clientY,target,pointerType:"touch"}  (pan/none 모드의 한 번 탭)
     onDoubleTap(p,info)  두 번 탭(보통 맞춤·초기화). 있으면 한 번 탭은 doubleTapMs 만큼 늦게 onTap 으로 간다.
     doubleTap            true(기본) | false | ()=>bool — 지금 두 번 탭을 기다릴지(예: 점 찍기 도구에서는 끔)
     doubleTapInPass      true 면 pass 모드의 빠른 두 탭도 onDoubleTap(p,{pass:true}) 로 알린다(앱의 탭 처리는 그대로 일어남)
     onLongPress(p)       길게 누르기(값 읽기 등). false 를 돌려주면 무시. 이어서 onLongPressMove(p), onLongPressEnd(p)
     longPressInPass      true 면 pass 모드에서도 길게 누르기를 잡는다(그 손가락은 onPassCancel 후 도우미가 가져감)
     onGestureStart(), onGestureEnd()
     filter(ev)           false 면 이 pointerdown 을 아예 다루지 않음(앱 그대로)
     box()                지역 좌표의 원점 {left,top} (client 좌표). 기본 el.getBoundingClientRect()
     touchAction          false 면 el 의 touch-action 을 건드리지 않음
     longPressMs=500, tapSlop=10, doubleTapMs=300, doubleTapSlop=30
   })  → detach()

   보조
   INNO_TOUCH.affine(v0,g,{min,max})   보기 {zoom,x,y} (화면 = (이미지 − x)·zoom) 에 핀치·이동을 적용한 새 값 — 핀치 중심 고정
   INNO_TOUCH.cssZoom(target,{frame,view,min,max,clip,onTap,onDoubleTap,onChange,...})
                                       확대 기능이 없는 캔버스에 CSS transform 으로 확대·이동(다시 그리지 않음).
                                       view = {s,tx,ty} 객체(앱이 다시 그려도 그대로 쓰도록 앱이 들고 있음)
                                       → {view, reset(), apply(), detach()}
   INNO_TOUCH.touchUI()                마지막 포인터가 손가락이었거나 거친 포인터(pointer:coarse) 기기인가 — 손잡이 크게 그리기 등
   INNO_TOUCH.lastPointer()            "touch" | "mouse" | "pen"
   ============================================================ */
(function(){
"use strict";
if(window.INNO_TOUCH) return;
const DEF={oneFinger:"pan",longPressMs:500,tapSlop:10,doubleTapMs:300,doubleTapSlop:30,doubleTap:true};
const now=()=>performance.now();
let LAST={type:"mouse",t:0};
const MQ=(()=>{ try{ return window.matchMedia("(pointer: coarse)"); }catch(e){ return null; } })();
try{ document.addEventListener("pointerdown",e=>{ LAST={type:e.pointerType||"mouse",t:now()};
  const r=document.documentElement; const on=e.pointerType==="touch"; if(r.classList.contains("inno-touch-last")!==on) r.classList.toggle("inno-touch-last",on); },true); }catch(e){}
const touchUI=()=>LAST.type==="touch"||!!(MQ&&MQ.matches);
const SWALLOW=["pointerover","pointerenter","pointerout","pointerleave"];

function attach(el,opts){
  if(!el) return ()=>{};
  const o=Object.assign({},DEF,opts||{});
  const P=new Map();            // 내려 있는 손가락 id → {x,y,x0,y0} (client)
  const owned=new Set();        // 도우미가 맡은(소비하는) 손가락 id
  let passId=null;              // 앱에 넘긴 손가락
  let gest=null, started=false, pendG=null, raf=0, lastG=null;
  let lp=null, lpOn=false, seq=null;   // seq: 이번 손가락 묶음 {t0,moved,multi,mode}
  let pendTap=null, passTap=null, suppressUntil=0;
  const prevTA=el.style.touchAction; if(o.touchAction!==false) el.style.touchAction="none";
  el.classList.add("inno-touch");
  const box=()=>o.box?o.box():el.getBoundingClientRect();
  const loc=(cx,cy)=>{ const b=box(); return {x:cx-b.left,y:cy-b.top}; };
  const mk=(cx,cy,target)=>Object.assign(loc(cx,cy),{clientX:cx,clientY:cy,target:target||el,pointerType:"touch"});
  const dtOn=()=>!!o.onDoubleTap&&(typeof o.doubleTap==="function"?o.doubleTap():o.doubleTap!==false);
  const modeOf=e=>{ const m=typeof o.oneFinger==="function"?o.oneFinger(e):o.oneFinger; return m==="pass"||m==="none"?m:"pan"; };
  const stop=e=>{ e.stopImmediatePropagation(); if(e.cancelable) e.preventDefault(); };
  const cap=id=>{ try{ el.setPointerCapture(id); }catch(e){} };
  const clearLP=()=>{ if(lp){ clearTimeout(lp); lp=null; } };
  function firePending(){ if(!pendTap) return; clearTimeout(pendTap.timer); const p=pendTap.p; pendTap=null; if(o.onTap) o.onTap(p); }
  function applyNow(){ if(raf){ cancelAnimationFrame(raf); raf=0; } if(pendG&&gest&&o.setView){ const g=pendG; pendG=null; o.setView(gest.v0,g); } pendG=null; }
  function schedule(g){ pendG=g; lastG=g; if(!raf) raf=requestAnimationFrame(()=>{ raf=0; const gg=pendG; pendG=null; if(gg&&gest&&o.setView) o.setView(gest.v0,gg); }); }
  function centre(fromStart){ const ids=[...P.keys()].filter(id=>owned.has(id)).slice(0,2); const a=P.get(ids[0]), b=ids.length>1?P.get(ids[1]):null; if(!a) return null;
    if(!b&&fromStart) return {ids,x:a.x0,y:a.y0,d:0,sx:0,sy:0};
    return b?{ids,x:(a.x+b.x)/2,y:(a.y+b.y)/2,d:Math.hypot(a.x-b.x,a.y-b.y),sx:Math.abs(a.x-b.x),sy:Math.abs(a.y-b.y)}:{ids,x:a.x,y:a.y,d:0,sx:0,sy:0}; }
  function begin(fromStart){ applyNow(); const c=centre(fromStart); if(!c){ gest=null; return; }
    if(!o.getView||!o.setView){ gest=null; return; }
    if(!started){ started=true; if(o.onGestureStart) o.onGestureStart(); }
    const L=loc(c.x,c.y); gest={ids:c.ids,v0:o.getView(),cx0:c.x,cy0:c.y,x0:L.x,y0:L.y,d0:c.d,sx0:c.sx,sy0:c.sy,pinch:c.ids.length>1}; }
  function update(){ if(!gest) return; const c=centre(); if(!c||c.ids.length!==gest.ids.length) return;
    const L=loc(c.x,c.y); const s=gest.pinch&&gest.d0>0?Math.max(0.02,c.d/gest.d0):1;
    const sX=gest.pinch&&gest.sx0>=24?Math.max(0.02,c.sx/gest.sx0):s, sY=gest.pinch&&gest.sy0>=24?Math.max(0.02,c.sy/gest.sy0):s;
    schedule({kind:gest.pinch?"pinch":"pan",x0:gest.x0,y0:gest.y0,x:L.x,y:L.y,dx:L.x-gest.x0,dy:L.y-gest.y0,scale:s,scaleX:sX,scaleY:sY,clientX0:gest.cx0,clientY0:gest.cy0,clientX:c.x,clientY:c.y}); }
  function endAll(){ applyNow(); const was=started; gest=null; started=false; if(was&&o.onGestureEnd) o.onGestureEnd(lastG); lastG=null; }
  function takePass(){ if(passId==null) return; const id=passId; passId=null; owned.add(id);
    if(o.onPassCancel){ try{ o.onPassCancel(); }catch(e){ console.error(e); } }
    else { try{ el.dispatchEvent(new PointerEvent("pointercancel",{pointerId:id,pointerType:"touch",bubbles:true})); }catch(e){} }
    cap(id); }

  function down(e){
    if(e.pointerType!=="touch"||!e.isTrusted) return;
    if(e.isPrimary&&P.size) reset();                         // 앞 묶음이 끝을 못 받았으면(요소가 바뀜 등) 버림
    if(!P.size&&o.filter&&o.filter(e)===false) return;
    P.set(e.pointerId,{x:e.clientX,y:e.clientY,x0:e.clientX,y0:e.clientY});
    if(P.size===1){
      if(pendTap&&Math.hypot(e.clientX-pendTap.p.clientX,e.clientY-pendTap.p.clientY)>o.doubleTapSlop) firePending();
      const m=modeOf(e); seq={t0:now(),moved:false,multi:false,mode:m,x:e.clientX,y:e.clientY,target:e.target}; lpOn=false;
      if(m==="pass"){ passId=e.pointerId; suppressUntil=0; if(o.onLongPress&&o.longPressInPass) armLP(e); return; }
      owned.add(e.pointerId); stop(e); cap(e.pointerId);
      if(o.onLongPress) armLP(e);
      return;
    }
    // 둘째(이상) 손가락: 언제나 핀치·이동
    stop(e); owned.add(e.pointerId); cap(e.pointerId);
    clearLP(); firePending(); if(seq) seq.multi=true;
    if(lpOn){ lpOn=false; if(o.onLongPressEnd) o.onLongPressEnd(null); }
    if(passId!=null) takePass();
    begin();
  }
  function armLP(e){ clearLP(); const t=e.target; lp=setTimeout(()=>{ lp=null; if(P.size!==1||!seq||seq.moved||seq.multi) return;
      const [id,q]=[...P.entries()][0]; const p=mk(q.x,q.y,t);
      if(passId!=null) takePass();
      if(o.onLongPress(p)===false) return; lpOn=true; firePending(); if(gest){ applyNow(); gest=null; } },o.longPressMs); }
  function move(e){
    if(e.pointerType!=="touch") return; const q=P.get(e.pointerId); if(!q) return;
    q.x=e.clientX; q.y=e.clientY;
    if(seq&&!seq.moved&&Math.hypot(q.x-q.x0,q.y-q.y0)>o.tapSlop){ seq.moved=true; clearLP(); firePending(); }
    if(e.pointerId===passId) return;                         // 앱이 받는다
    if(!owned.has(e.pointerId)) return;
    stop(e);
    if(lpOn){ if(o.onLongPressMove) o.onLongPressMove(mk(q.x,q.y,e.target)); return; }
    if(!gest&&seq&&seq.moved&&seq.mode==="pan"&&P.size===1) begin(true);   // 끌기 시작: 누른 자리 기준이라 보기가 손가락을 정확히 따라감
    if(gest) update();
  }
  function up(e){
    if(e.pointerType!=="touch") return; const q=P.get(e.pointerId); if(!q) return;
    P.delete(e.pointerId);
    if(e.pointerId===passId){ passId=null; clearLP();
      const tap=e.type==="pointerup"&&seq&&!seq.moved&&!seq.multi&&now()-seq.t0<o.longPressMs;
      if(tap&&o.doubleTapInPass&&o.onDoubleTap){ const p=mk(e.clientX,e.clientY,e.target);
        if(passTap&&now()-passTap.t<o.doubleTapMs&&Math.hypot(p.clientX-passTap.x,p.clientY-passTap.y)<o.doubleTapSlop){ passTap=null; o.onDoubleTap(p,{pass:true}); }
        else passTap={t:now(),x:p.clientX,y:p.clientY}; }
      if(!P.size) seq=null; return; }
    if(!owned.has(e.pointerId)) return;
    stop(e); suppressUntil=now()+700; const pid=e.pointerId; setTimeout(()=>owned.delete(pid),500);
    if(lpOn&&!P.size){ lpOn=false; clearLP(); if(o.onLongPressEnd) o.onLongPressEnd(mk(q.x,q.y,e.target)); seq=null; endAll(); return; }
    if(P.size){ if(gest) begin(); return; }                  // 손가락 하나 남음 → 그 손가락으로 이동을 이어 감
    clearLP();
    const s=seq; seq=null; endAll();
    const tap=e.type==="pointerup"&&s&&!s.moved&&!s.multi&&!lpOn&&now()-s.t0<o.longPressMs;
    if(!tap) return;
    const p=mk(e.clientX,e.clientY,s.target);
    if(dtOn()){
      if(pendTap&&now()-pendTap.t<o.doubleTapMs+120&&Math.hypot(p.clientX-pendTap.p.clientX,p.clientY-pendTap.p.clientY)<o.doubleTapSlop){ clearTimeout(pendTap.timer); pendTap=null; o.onDoubleTap(p,{pass:false}); }
      else { firePending(); pendTap={p,t:now(),timer:setTimeout(()=>{ const t=pendTap; pendTap=null; if(t&&o.onTap) o.onTap(t.p); },o.doubleTapMs)}; }
    } else if(o.onTap) o.onTap(p);
  }
  function swallowCompat(e){ if(now()<suppressUntil){ e.stopImmediatePropagation(); e.preventDefault(); } }
  function swallowHover(e){ if(e.pointerType==="touch"&&owned.has(e.pointerId)) e.stopImmediatePropagation(); }
  function lost(e){ if(e.pointerType==="touch"&&owned.has(e.pointerId)&&P.has(e.pointerId)&&!el.isConnected) reset(); }
  function reset(){ clearLP(); if(raf){ cancelAnimationFrame(raf); raf=0; } pendG=null; const was=started; gest=null; started=false; P.clear(); passId=null; lpOn=false; seq=null; if(was&&o.onGestureEnd) o.onGestureEnd(null); }
  function ctx(e){ if(LAST.type==="touch"&&now()-LAST.t<2000) e.preventDefault(); }
  const L=[["pointerdown",down],["pointermove",move],["pointerup",up],["pointercancel",up],["click",swallowCompat],["dblclick",swallowCompat],["contextmenu",ctx],["lostpointercapture",lost]].concat(SWALLOW.map(t=>[t,swallowHover]));
  L.forEach(([t,f])=>el.addEventListener(t,f,{capture:true,passive:false}));
  const detach=()=>{ L.forEach(([t,f])=>el.removeEventListener(t,f,{capture:true})); clearLP(); if(pendTap) clearTimeout(pendTap.timer); if(raf) cancelAnimationFrame(raf);
    if(o.touchAction!==false) el.style.touchAction=prevTA; el.classList.remove("inno-touch"); P.clear(); owned.clear(); };
  detach.el=el; detach.opts=o;
  return detach;
}

/* 보기 {zoom,x,y}: 화면 = (이미지 − x)·zoom. 시작 중심 아래의 이미지 점이 지금 중심 아래에 오도록 */
function affine(v0,g,lim){ lim=lim||{}; const mn=lim.min!=null?lim.min:1e-3, mx=lim.max!=null?lim.max:1e3;
  const ix=v0.x+g.x0/v0.zoom, iy=v0.y+g.y0/v0.zoom; const z=Math.min(mx,Math.max(mn,v0.zoom*g.scale));
  return {zoom:z,x:ix-g.x/z,y:iy-g.y/z}; }

/* CSS transform 확대 (그 자체로 확대 기능이 없는 캔버스용). target 의 배치 상자를 기준으로 translate·scale.
   frame = 손가락을 받는 요소(기본 target). frame===target 이면 clip-path 로 원래 상자 밖을 잘라 레이아웃·이웃을 덮지 않는다. */
function cssZoom(target,opt){
  opt=opt||{}; const frame=opt.frame||target, self=frame===target; const view=opt.view||{s:1,tx:0,ty:0};
  if(!(view.s>0)){ view.s=1; view.tx=0; view.ty=0; }
  const mn=opt.min!=null?opt.min:1, mx=opt.max!=null?opt.max:8, clip=opt.clip!=null?opt.clip:self;
  const size=()=>({W:target.offsetWidth||target.clientWidth||1,H:target.offsetHeight||target.clientHeight||1});
  const layoutLeftTop=()=>{ const r=target.getBoundingClientRect(); return {left:r.left-view.tx,top:r.top-view.ty}; };
  const box=()=>self?layoutLeftTop():frame.getBoundingClientRect();
  const lt=()=>{ if(self) return {x:0,y:0}; const a=layoutLeftTop(), f=frame.getBoundingClientRect(); return {x:a.left-f.left,y:a.top-f.top}; };
  function bound(){ const {W,H}=size(); if(view.s<=1){ view.s=Math.max(mn,view.s); if(view.s===1){ view.tx=0; view.ty=0; } return; }
    view.tx=Math.min(0,Math.max(W*(1-view.s),view.tx)); view.ty=Math.min(0,Math.max(H*(1-view.s),view.ty)); }
  function apply(){ bound(); const id=view.s===1&&!view.tx&&!view.ty;
    target.style.transformOrigin="0 0"; target.style.transform=id?"":`translate(${view.tx}px,${view.ty}px) scale(${view.s})`;
    if(clip){ if(id) target.style.clipPath=""; else { const {W,H}=size(), s=view.s; const l=-view.tx/s, t=-view.ty/s, r=W-(W-view.tx)/s, b=H-(H-view.ty)/s;
      target.style.clipPath=`inset(${Math.max(0,t)}px ${Math.max(0,r)}px ${Math.max(0,b)}px ${Math.max(0,l)}px)`; } }
    target.classList.toggle("inno-zoomed",!id); if(opt.onChange) opt.onChange(view); }
  function reset(){ view.s=1; view.tx=0; view.ty=0; apply(); }
  const det=attach(frame,Object.assign({},opt,{
    box,
    getView:()=>{ const L=lt(); return {s:view.s,tx:view.tx,ty:view.ty,lx:L.x,ly:L.y}; },
    setView:(v0,g)=>{ const ux=(g.x0-v0.lx-v0.tx)/v0.s, uy=(g.y0-v0.ly-v0.ty)/v0.s; const s=Math.min(mx,Math.max(mn,v0.s*g.scale));
      view.s=s; view.tx=g.x-v0.lx-s*ux; view.ty=g.y-v0.ly-s*uy; apply(); },
    onDoubleTap:opt.onDoubleTap||(()=>reset())
  }));
  if(view.s!==1||view.tx||view.ty) apply();
  return {view,reset,apply,detach:det};
}

window.INNO_TOUCH={version:"1.0",attach,affine,cssZoom,touchUI,lastPointer:()=>LAST.type};
})();
