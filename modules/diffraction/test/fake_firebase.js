/* 시험용 Firebase (app·auth·firestore) — 메모리 위에, 게시된 Firestore 규칙을 그대로 집행한다.
   diff_db.js 가 쓰는 모듈러 SDK 함수만: initializeApp/getApps, getAuth/GoogleAuthProvider/signInWithPopup/signInWithRedirect/
   getRedirectResult/onAuthStateChanged/signOut, getFirestore/doc/collection/getDoc(FromServer)/getDocs(FromServer)/setDoc/updateDoc/
   deleteDoc/runTransaction/serverTimestamp.
   규칙의 '오류 = 거부' 의미: 없는 필드끼리 비교(== 양쪽 undefined)·없는 문서 get() 은 false 로 본다. */
"use strict";
class FakeTs{ constructor(ms){ this.seconds=Math.floor(ms/1000); this.nanoseconds=(ms%1000)*1e6; } toMillis(){ return this.seconds*1000+this.nanoseconds/1e6; } }
const SENT={__serverTs:true};
const fe=(code,message)=>{ const e=new Error(message||code); e.code=code; e.name="FirebaseError"; return e; };
const dcopy=v=>{ if(v instanceof FakeTs) return new FakeTs(v.toMillis()); if(v instanceof Date) return new Date(v.getTime()); if(Array.isArray(v)) return v.map(dcopy); if(v&&typeof v==="object"){ const o={}; for(const k of Object.keys(v)) o[k]=dcopy(v[k]); return o; } return v; };
function checkVal(v,p){ if(v===undefined) throw fe("invalid-argument","Function setDoc() called with invalid data. Unsupported field value: undefined (found in field "+p+")"); if(v&&typeof v==="object"&&!(v instanceof FakeTs)&&!(v instanceof Date)&&v!==SENT) for(const k of Object.keys(v)) checkVal(v[k],p+"."+k); }

function makeFake(opt){ opt=opt||{};
  const store=new Map();        // "a/b/c/d" → data
  const F={store,offline:false,log:[],clock:1.7e12,nextUser:null,popupError:null,redirectPending:null,denials:0};
  const tick=()=>(F.clock+=1000);
  const resolveTs=d=>{ const o={}; for(const k of Object.keys(d)){ const v=d[k]; o[k]=v===SENT?new FakeTs(tick()):(v&&typeof v==="object"&&!(v instanceof FakeTs)&&!(v instanceof Date)&&!Array.isArray(v))?resolveTs(v):dcopy(v); } return o; };

  /* ---------- 규칙 ---------- */
  const eq=(a,b)=>a!==undefined&&b!==undefined&&JSON.stringify(a)===JSON.stringify(b);
  function R(auth,state){ const st=state||store;
    const signedIn=()=>!!(auth&&auth.token&&auth.token.email_verified===true);
    const email=()=>auth.token.email;
    const isAdmin=()=>{ if(!signedIn()) return false; const a=st.get("config/admins"); return !!(a&&Array.isArray(a.emails)&&a.emails.includes(email())); };
    const isMember=()=>{ if(!signedIn()) return false; if(isAdmin()) return true; const m=st.get("members/"+email()); return !!(m&&m.active===true); };
    return {signedIn,isAdmin,isMember}; }
  /* op: get | list | create | update | delete ; path = 문서 경로 (list 는 컬렉션 경로) */
  function allowed(op,path,auth,newData,oldData){ const seg=path.split("/"); const r=R(auth); const uid=auth&&auth.uid;
    const read=op==="get"||op==="list";
    if(seg[0]==="config"&&(seg.length===2||(op==="list"&&seg.length===1))) return read?r.isMember():r.isAdmin();
    if(seg[0]==="members"&&(seg.length===2||(op==="list"&&seg.length===1))) return read?r.isMember():r.isAdmin();
    if(seg[0]==="phases"&&(seg.length===2||(op==="list"&&seg.length===1))){
      if(read) return r.isMember();
      if(op==="create") return r.isMember()&&eq(newData.createdBy,uid)&&newData.version===1&&newData.hidden===false;
      if(op==="update") return r.isMember()&&(r.isAdmin()||eq(oldData.createdBy,uid))&&eq(newData.createdBy,oldData.createdBy)&&typeof oldData.version==="number"&&newData.version===oldData.version+1&&(r.isAdmin()||eq(newData.verified,oldData.verified));
      return false; }
    if(seg[0]==="phases"&&seg[2]==="history"&&(seg.length===4||(op==="list"&&seg.length===3))){ if(read||op==="create") return r.isMember(); return false; }
    if(seg[0]==="calibrations"&&(seg.length===2||(op==="list"&&seg.length===1))){
      if(read) return r.isMember();
      if(op==="create") return r.isMember()&&eq(newData.createdBy,uid);
      if(op==="update") return r.isMember()&&(r.isAdmin()||eq(oldData.createdBy,uid));
      return false; }
    return false; }
  function deny(op,path){ F.denials++; F.log.push({deny:op,path}); return fe("permission-denied","Missing or insufficient permissions."); }

  /* ---------- app ---------- */
  const apps=[];
  const app={ initializeApp:(cfg,name)=>{ const a={name:name||"[DEFAULT]",options:Object.assign({},cfg)}; apps.push(a); F.lastConfig=a.options; return a; }, getApps:()=>apps.slice() };

  /* ---------- auth ---------- */
  const AU={currentUser:null,listeners:[]};
  const setUser=u=>{ AU.currentUser=u?Object.assign({},u):null; for(const cb of AU.listeners.slice()) setTimeout(()=>cb(AU.currentUser),0); };
  const tokenOf=()=>AU.currentUser?{uid:AU.currentUser.uid,token:{email:AU.currentUser.email,email_verified:AU.currentUser.emailVerified!==false}}:null;
  class GoogleAuthProvider{ setCustomParameters(p){ this.params=p; } }
  const auth={ getAuth:()=>AU, GoogleAuthProvider,
    signInWithPopup:async()=>{ if(F.popupError){ const c=F.popupError; F.popupError=null; throw fe(c); } if(!F.nextUser) throw fe("auth/popup-closed-by-user"); setUser(F.nextUser); return {user:AU.currentUser}; },
    signInWithRedirect:async()=>{ F.redirectPending=F.nextUser; F.redirectCalls=(F.redirectCalls||0)+1; },
    getRedirectResult:async()=>{ if(F.redirectPending){ const u=F.redirectPending; F.redirectPending=null; setUser(u); return {user:AU.currentUser}; } return null; },
    onAuthStateChanged:(a,cb)=>{ AU.listeners.push(cb); setTimeout(()=>cb(AU.currentUser),0); return ()=>{ AU.listeners=AU.listeners.filter(x=>x!==cb); }; },
    signOut:async()=>setUser(null) };

  /* ---------- firestore ---------- */
  const DBO={};
  const join=(base,segs)=>[base].concat(segs).filter(Boolean).join("/");
  const mkDoc=p=>({type:"document",path:p,id:p.split("/").pop()});
  const mkCol=p=>({type:"collection",path:p,id:p.split("/").pop()});
  const doc=(b,...s)=>{ const p=b===DBO?join("",s):join(b.path,s); if(p.split("/").length%2) throw fe("invalid-argument","odd doc path "+p); return mkDoc(p); };
  const collection=(b,...s)=>{ const p=b===DBO?join("",s):join(b.path,s); if(!(p.split("/").length%2)) throw fe("invalid-argument","even col path "+p); return mkCol(p); };
  const snap=(p,d)=>({id:p.split("/").pop(),ref:mkDoc(p),exists:()=>d!==undefined,data:()=>d===undefined?undefined:dcopy(d),metadata:{fromCache:false}});
  const online=()=>{ if(F.offline) throw fe("unavailable","Failed to get document because the client is offline."); };
  async function getDoc(ref){ online(); if(!allowed("get",ref.path,tokenOf())) throw deny("get",ref.path); return snap(ref.path,store.get(ref.path)); }
  async function getDocs(col){ online(); if(!allowed("list",col.path,tokenOf())) throw deny("list",col.path); const n=col.path.split("/").length+1; const docs=[];
    for(const [p,d] of store) if(p.startsWith(col.path+"/")&&p.split("/").length===n) docs.push(snap(p,d));
    docs.sort((a,b)=>a.id<b.id?-1:1); return {docs,size:docs.length,empty:!docs.length,forEach:f=>docs.forEach(f),metadata:{fromCache:false}}; }
  /* 쓰기 하나를 규칙으로 판정 (state = 쓰기 전 상태) */
  function judge(w,auth,state){ const old=state.get(w.path); let next;
    if(w.kind==="set"){ checkVal(w.data,"data"); next=resolveTs(w.data); }
    else if(w.kind==="update"){ if(old===undefined) throw fe("not-found","No document to update: "+w.path); checkVal(w.data,"data"); next=Object.assign(dcopy(old),resolveTs(w.data)); }
    else if(w.kind==="delete"){ next=undefined; }
    const op=w.kind==="delete"?"delete":old===undefined?"create":"update";
    if(!allowed(op,w.path,auth,next,old)) throw deny(op,w.path);
    return next; }
  async function commit(writes){ online(); const auth=tokenOf(); const state=new Map(store); const out=[];
    for(const w of writes){ const n=judge(w,auth,state); out.push([w.path,n]); }      // 규칙의 get() 은 쓰기 전 상태를 본다
    for(const [p,n] of out){ if(n===undefined) store.delete(p); else store.set(p,n); F.log.push({write:p}); } }
  const firestore={ getFirestore:()=>DBO, doc, collection, serverTimestamp:()=>SENT,
    getDoc, getDocFromServer:getDoc, getDocs, getDocsFromServer:getDocs,
    setDoc:async(ref,data)=>commit([{kind:"set",path:ref.path,data}]),
    updateDoc:async(ref,data)=>commit([{kind:"update",path:ref.path,data}]),
    deleteDoc:async ref=>commit([{kind:"delete",path:ref.path}]),
    runTransaction:async(db,fn)=>{ online(); const writes=[]; const t={
        get:async ref=>{ if(writes.length) throw fe("invalid-argument","Firestore transactions require all reads to be executed before all writes."); return getDoc(ref); },
        set:(ref,data)=>{ writes.push({kind:"set",path:ref.path,data}); return t; },
        update:(ref,data)=>{ writes.push({kind:"update",path:ref.path,data}); return t; },
        delete:ref=>{ writes.push({kind:"delete",path:ref.path}); return t; } };
      const r=await fn(t); await commit(writes); return r; } };

  F.sdk={app,auth,firestore};
  F.user=(uid,email,o)=>Object.assign({uid,email,emailVerified:true,displayName:uid},o||{});
  F.loginAs=u=>{ setUser(u); };
  F.seed=(path,data)=>store.set(path,dcopy(data));
  F.get=path=>store.has(path)?dcopy(store.get(path)):undefined;
  F.Ts=FakeTs; F.SENT=SENT;
  return F; }
module.exports={makeFake,FakeTs};
