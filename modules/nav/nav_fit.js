/* 탭 막대 맞춤 (데스크톱 ≥ 721 px) — 10-19
   한 줄에 다 들어가면 전체 이름 → 안 들어가면 휴대폰과 같은 짧은 이름(data-short) → 그래도 넘치면 뒤쪽 탭부터 'More ▾' 메뉴로.
   지금 보고 있는 탭은 접지 않는다. 휴대폰(≤ 720 px)은 아래 고정 막대(CSS 격자)라 손대지 않는다.
   renderView 가 탭을 다시 그릴 때마다, 그리고 막대 폭이 바뀔 때마다 부른다. */
(function(){
  "use strict";
  let menu=null, moreBtn=null, raf=0, ro=null;
  const mobile=()=>window.matchMedia&&matchMedia("(max-width:720px)").matches;
  function closeMenu(){ if(menu){ menu.remove(); menu=null; } if(moreBtn) moreBtn.setAttribute("aria-expanded","false"); }
  function openMenu(btn,hid){ closeMenu(); moreBtn=btn; btn.setAttribute("aria-expanded","true");
    menu=document.createElement("div"); menu.className="nav-menu"; menu.setAttribute("role","menu");
    hid.forEach(b=>{ const it=document.createElement("button"); it.type="button"; it.setAttribute("role","menuitem"); it.disabled=b.disabled;
      const c=b.querySelector(".cnt"); it.innerHTML=""; it.textContent=b.getAttribute("aria-label")||b.textContent; if(c&&c.textContent){ const s=document.createElement("span"); s.className="cnt"; s.textContent=c.textContent; it.appendChild(s); }
      it.dataset.v=b.dataset.v; it.onclick=()=>{ closeMenu(); b.click(); }; menu.appendChild(it); });
    document.body.appendChild(menu); const r=btn.getBoundingClientRect(); const w=menu.offsetWidth;
    menu.style.top=Math.round(r.bottom+4)+"px"; menu.style.left=Math.round(Math.max(8,Math.min(window.innerWidth-w-8,r.right-w)))+"px";
    const f=menu.querySelector("button:not(:disabled)"); if(f) f.focus(); }
  function fit(){
    const nav=document.getElementById("nav"); if(!nav) return;
    closeMenu(); nav.classList.remove("nav-short");
    const old=nav.querySelector(".nav-more"); if(old) old.remove();
    const bs=[...nav.querySelectorAll("button[data-v]")]; bs.forEach(b=>b.classList.remove("nav-hid"));
    if(!ro&&window.ResizeObserver){ ro=new ResizeObserver(()=>{ cancelAnimationFrame(raf); raf=requestAnimationFrame(fit); }); ro.observe(nav); }
    if(mobile()) return;
    const over=()=>nav.scrollWidth>nav.clientWidth+1;
    if(!over()) return;
    nav.classList.add("nav-short"); if(!over()) return;
    const m=document.createElement("button"); m.type="button"; m.className="nav-more"; m.setAttribute("aria-haspopup","menu"); m.setAttribute("aria-expanded","false"); m.textContent="More ▾"; nav.appendChild(m);
    const hid=[]; for(let i=bs.length-1;i>=0&&over();i--){ if(bs[i].classList.contains("active")) continue; bs[i].classList.add("nav-hid"); hid.unshift(bs[i]); }
    if(!hid.length){ m.remove(); return; }
    if(hid.some(b=>b.querySelector(".cnt")&&+b.querySelector(".cnt").textContent>0)) m.classList.add("has");
    m.title=hid.map(b=>b.getAttribute("aria-label")).join(" · ");
    m.onclick=e=>{ e.stopPropagation(); if(menu) closeMenu(); else openMenu(m,hid); };
  }
  document.addEventListener("click",e=>{ if(menu&&!menu.contains(e.target)) closeMenu(); });
  document.addEventListener("keydown",e=>{ if(e.key==="Escape"&&menu){ const b=moreBtn; closeMenu(); if(b) b.focus(); } });
  window.addEventListener("resize",()=>{ cancelAnimationFrame(raf); raf=requestAnimationFrame(fit); });
  window.NAV_FIT=fit;
})();
