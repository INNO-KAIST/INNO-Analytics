// Phase DB 시험: Local CRUD · CIF 가져오기 · 검사 · 중복 · 원소 거르기 · 내보내기/가져오기 · Lab(규칙 집행 가짜 Firestore)
// 사용: node test/test_db.js   (fake-indexeddb 는 modules/diffraction/node_modules 에 — npm i --no-save fake-indexeddb)
"use strict";
const fs=require("fs"), path=require("path");
const {IDBFactory}=require("fake-indexeddb");
require("../src/diff_sg.js"); require("../src/diff_scat.js"); const DC=require("../src/diff_core.js");
const DB=require("../src/diff_db.js");
const {makeFake}=require("./fake_firebase.js");
let pass=0, fail=0; const fails=[];
const ok=(c,m,d)=>{ if(c) pass++; else { fail++; fails.push(m+(d!==undefined?" — "+(typeof d==="string"?d:JSON.stringify(d)):"")); } if(!c||process.env.V) console.log((c?"  ✓ ":"  ✗ ")+m+(c||d===undefined?"":" — "+(typeof d==="string"?d:JSON.stringify(d)))); };
async function rejects(p,code,m){ try{ await p; ok(false,m,"예외 없음"); return null; }catch(e){ ok(!code||e.code===code,m,`${e.code}: ${e.message}`); return e; } }
const CIF=n=>fs.readFileSync(path.join(__dirname,"cif",n),"utf8");
const CIFS=fs.readdirSync(path.join(__dirname,"cif")).filter(f=>/\.cif$/.test(f)).sort();
const content=r=>JSON.stringify(DB.contentOf(r));
const flush=()=>new Promise(r=>setTimeout(r,5));

(async()=>{
let idb=new IDBFactory(); DB._setIDB(idb); DB._setEnv({protocol:"https:"});
await DB.init();
ok(DB.storage()==="indexeddb","IndexedDB 열림");
const packBefore=DB.packInfo().count;
console.log(`== 기본 상(pack): 지금 ${packBefore} 개 (diff_pack.js ${fs.existsSync(path.join(__dirname,"../src/diff_pack.js"))?"있음":"없음"})`);
{ const save=globalThis.DIFF_PACK; globalThis.DIFF_PACK={version:"test-1",phases:[{id:"pack:Si",name:"Silicon",formula:"Si",cell:{a:0.5431,b:0.5431,c:0.5431,al:90,be:90,ga:90},sg:{number:227,symbol:"Fd-3m"},source:{kind:"COD",codId:"9008565"}}]};
  ok(DB.list({layers:["pack"]}).length===1&&DB.packInfo().version==="test-1","나중에 들어온 DIFF_PACK 을 집어 든다");
  ok(DB.get("pack:Si").layer==="pack"&&DB.refOf("pack:Si").packVersion==="test-1","pack refOf = packVersion");
  ok(DB.phase("pack:Si").ops.length===192,"pack 상 → DC.phaseFromRecord (Fd-3m 192 연산)");
  // 아래 시험은 빈 pack 으로 (실제 diff_pack.js 가 생겨도 결과가 같게)
  globalThis.DIFF_PACK={version:"empty",phases:[]}; ok(DB.packInfo().count===0,"빈 pack 으로 바꿈"); void save; }

console.log("== CIF 가져오기 (test/cif 10 개)");
const imported={};
for(const f of CIFS){ const r=DB.importCIF(CIF(f),{fileName:f}); const it=r[0];
  ok(r.length===1&&it.rec&&it.validation.ok,`${f}: 블록 1 · 검사 통과`,it.error||(it.validation&&it.validation.errors));
  imported[f]=it.rec; }
{ const multi=DB.importCIF(CIF("Si_Fd-3m.cif")+"\n"+CIF("NaCl_Fm-3m.cif")+"\ndata_broken\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n",{fileName:"multi.cif"});
  ok(multi.length===3&&multi[0].rec&&multi[1].rec&&multi[2].error,"여러 블록 CIF: 2 개 읽고 깨진 블록은 오류로",multi.map(x=>x.error||x.rec.name));
  const twice=DB.importCIF(CIF("Si_Fd-3m.cif")+"\n"+CIF("Si_Fd-3m.cif").replace(/^data_\S+/m,"data_again"));
  ok(twice[1].duplicate&&twice[1].duplicate.layer==="import","같은 파일 안 중복 블록 표시");
  ok(DB.importCIF("hello")[0].error,"구조 없는 글 → 오류");
  const si=imported["Si_Fd-3m.cif"]; ok(si.symops&&si.symops.length>=96&&si.source.kind&&/CIF Si_Fd-3m.cif/.test(si.source.note),"CIF 기록: symops·출처 note",si.source); }

console.log("== Local CRUD");
const ids={};
for(const f of CIFS){ const res=await DB.saveLocal(imported[f]); ok(res.ok&&/^local:/.test(res.rec.id)&&res.rec.version===1,`${f} → Local 저장`,res.reason); ids[f]=res.rec&&res.rec.id; }
ok(DB.list({layers:["local"]}).length===CIFS.length,"Local 목록 수");
{ const dup=await DB.saveLocal(imported["Au_Fm-3m.cif"]); ok(!dup.ok&&dup.reason==="duplicate"&&dup.duplicate.id===ids["Au_Fm-3m.cif"],"같은 상 다시 → 중복으로 돌려줌 (저장 안 함)");
  ok(DB.list({layers:["local"]}).length===CIFS.length,"중복은 쓰지 않음");
  const forced=await DB.saveLocal(imported["Au_Fm-3m.cif"],{force:true}); ok(forced.ok,"force → 저장"); await DB.removeLocal(forced.rec.id);
  ok(!DB.get(forced.rec.id),"removeLocal"); }
{ const id=ids["NaCl_Fm-3m.cif"]; const r0=DB.get(id); const ref1=DB.refOf(id);
  const ed=await DB.saveLocal(Object.assign({},r0,{note:"고친 메모",cell:Object.assign({},r0.cell,{a:0.5650,b:0.5650,c:0.5650})}));
  ok(ed.ok&&ed.rec.version===2&&ed.rec.id===id,"편집 → version 2, id 그대로");
  const ref2=DB.refOf(id); ok(ref2.version===2&&ref2.hash!==ref1.hash&&ref2.layer==="local","refOf: version·hash 가 바뀜");
  const back=await DB.getVersion(ref1); ok(back&&back.exact&&Math.abs(back.rec.cell.a-r0.cell.a)<1e-12,"getVersion(옛 ref) → 옛 내용 (Local history)");
  const same=await DB.saveLocal(DB.get(id)); ok(same.unchanged&&DB.get(id).version===2,"내용이 같으면 version 그대로"); }

console.log("== 검사 결과 전달");
{ const bad=await DB.saveLocal({name:"Si Å",formula:"Si",cell:{a:5.431},sg:{number:227}}); ok(bad.ok&&bad.validation.warnings.some(w=>/Å/.test(w)),"Å 로 넣은 격자 → 경고와 함께 저장",bad.validation&&bad.validation.warnings); if(bad.ok) await DB.removeLocal(bad.rec.id);
  const bad2=await DB.saveLocal({name:"x",formula:"Fe",cell:{a:0.3,b:0.31,c:0.3},sg:{number:229}}); ok(!bad2.ok&&bad2.reason==="invalid"&&bad2.validation.errors.some(e=>/입방정/.test(e)),"입방정에 a≠b → 오류, 저장 안 함",bad2.validation.errors);
  const bad3=await DB.saveLocal({name:"x",formula:"Fe",cell:{a:0.3},sg:{symbol:"Qq"}}); ok(!bad3.ok&&bad3.validation.errors.length,"모르는 공간군 기호 → 오류"); }

console.log("== 원소 거르기·검색");
{ const L=o=>DB.list(Object.assign({layers:["local"]},o)).map(r=>r.formula.replace(/\s+/g,""));
  const sub=L({elements:"Na Y F"}); ok(sub.length===1&&/Na/.test(sub[0])&&/Y/.test(sub[0]),"⊆ {Na,Y,F} → β-NaYF4 만",sub);
  const sub2=L({elements:["Y","F","Li","Na"]}); ok(sub2.length===2,"⊆ {Y,F,Li,Na} → LiYF4·NaYF4",sub2);
  const all=L({elements:"Y",elementMode:"all"}); ok(all.length===3,"Y 를 포함 → Y2O3·LiYF4·NaYF4",all);
  ok(L({elements:"Y",elementMode:"contains"}).length===3,"elementMode 'contains' = 'all' (diff_ui.js 값)");
  const ox=L({elements:"o",elementMode:"all"}); ok(ox.length===5,"O 포함 (소문자 입력) 5 개",ox);
  ok(DB.list({query:"TiO2"}).length===1,"검색: 화학식 비 (TiO2 = 'Ti2 O4')"); ok(DB.list({query:"Na2Y2F12"}).length===1,"검색: 비가 같은 식 (Na2Y2F12 = 'Na1 Y1 F6')"); ok(DB.list({query:"periclase"}).length===0&&DB.list({query:"ti2"}).length===1,"검색: 소문자 부분 문자열");
  const cod=imported["Au_Fm-3m.cif"].source.codId; if(cod){ ok(DB.list({query:cod}).length>=1,"검색: COD id "+cod); } else ok(true,"(Au CIF 에 COD id 없음)"); }

console.log("== 내보내기·가져오기 왕복");
{ const before=DB.list({layers:["local"]}); const js=DB.exportJSON(); const csv=DB.exportCSV(null,{layers:["local"]});
  ok(csv.split("\r\n").filter(Boolean).length===before.length+1&&/^id,layer,name,formula/.test(csv),"CSV 요약: 머리 + 행 수");
  await DB._reset(); idb=new IDBFactory(); DB._setIDB(idb); await DB.init(); ok(DB.list({layers:["local"]}).length===0,"새 저장소는 비어 있음");
  const rep=await DB.importJSON(js); ok(rep.added===before.length&&!rep.errors.length,"JSON 가져오기: 전부 추가",rep);
  const after=DB.list({layers:["local"]}); const m=new Map(after.map(r=>[r.id,r]));
  ok(before.every(r=>m.has(r.id)&&content(m.get(r.id))===content(r)&&m.get(r.id).version===r.version),"왕복: id·내용·version 같음");
  const nacl=ids["NaCl_Fm-3m.cif"]; ok((m.get(nacl).history||[]).length===1,"왕복: Local history 도 보존");
  const rep2=await DB.importJSON(js); ok(rep2.unchanged===before.length&&rep2.added===0,"같은 백업 다시 → 바뀐 것 없음",rep2);
  await DB._reset(); DB._setIDB(idb); await DB.init(); ok(DB.list({layers:["local"]}).length===before.length,"다시 열어도(IndexedDB) 남아 있음");
  await rejects(DB.importJSON("{nope"),"bad-json","깨진 JSON → 오류"); }

/* ============================================================ Lab */
console.log("== Lab: 가짜 Firestore (게시된 규칙 집행)");
const FK=makeFake(); DB._setSDK(FK.sdk);
const ADMIN=FK.user("u-admin","innolab@example.com"), BOB=FK.user("u-bob","ALLOWED_EMAIL_1@example.com"), CAROL=FK.user("u-carol","ALLOWED_EMAIL_2@example.com"), EVE=FK.user("u-eve","ALLOWED_EMAIL_3@example.com"), UNV=FK.user("u-unv","ALLOWED_EMAIL_5@example.com",{emailVerified:false});
FK.seed("config/admins",{emails:["innolab@example.com"]});
async function as(u){ await DB._reset(); DB._setIDB(idb); DB._setSDK(FK.sdk); DB._setEnv({protocol:"https:"}); FK.loginAs(u); await DB.init(); await DB.connect(); await flush(); return DB.firebaseStatus(); }

{ await DB._reset(); DB._setIDB(idb); DB._setEnv({protocol:"file:"}); await DB.init();
  const st=await DB.connect(); ok(st.fileProtocol&&st.error&&st.error.code==="file-protocol"&&/웹 주소/.test(st.error.message),"file:// → 로그인 불가 안내",st.error);
  await rejects(DB.signIn(),"file-protocol","file:// 에서 signIn → 오류");
  ok(DB.list({layers:["local"]}).length===CIFS.length,"file:// 에서도 Local 은 그대로"); }
{ await DB._reset(); DB._setIDB(idb); DB._setSDK(FK.sdk); DB._setEnv({protocol:"https:"}); FK.loginAs(null); await DB.init();
  const st0=DB.firebaseStatus(); ok(!st0.sdkLoaded&&!st0.connected,"init 만으로는 SDK 를 부르지 않음 (지연 로딩)");
  await DB.connect(); const st=DB.firebaseStatus(); ok(st.connected&&!st.user&&st.projectId==="inno-analytics-base"&&st.configOk,"로그인 전: projectId 표시, user 없음");
  ok(FK.lastConfig&&FK.lastConfig.projectId==="inno-analytics-base","initializeApp 에 기본 config");
  await rejects(DB.addLab(imported["Si_Fd-3m.cif"]),"not-signed-in","로그인 전 addLab → 거부");
  FK.popupError="auth/popup-blocked"; FK.nextUser=BOB; const r=await DB.signIn(); ok(r.redirect&&FK.redirectCalls===1,"팝업 막힘 → redirect 로 전환");
  await DB._reset(); DB._setIDB(idb); DB._setSDK(FK.sdk); await DB.init(); await DB.connect(); await flush(); const st2=DB.firebaseStatus();
  ok(st2.user&&st2.user.email==="ALLOWED_EMAIL_1@example.com"&&st2.wasSignedIn,"redirect 뒤 다시 열면 getRedirectResult 로 로그인됨");
  FK.popupError=null; FK.nextUser=CAROL; await DB.signOut(); const r3=await DB.signIn(); ok(r3.user&&r3.user.email==="ALLOWED_EMAIL_2@example.com"&&!r3.redirect,"팝업 로그인"); await DB.signOut(); }
{ const st=await as(EVE); ok(st.user&&st.user.email==="ALLOWED_EMAIL_3@example.com"&&st.isMember===false&&/구성원으로 등록되지 않았습니다/.test(st.error.message)&&st.error.message.includes("ALLOWED_EMAIL_3@example.com"),"비구성원: 안내 + 이메일",st.error);
  await rejects(DB.sync(),"not-member","비구성원 sync → 거부");
  let denied=false; try{ await FK.sdk.firestore.getDocs(FK.sdk.firestore.collection(FK.sdk.firestore.getFirestore(),"phases")); }catch(e){ denied=e.code==="permission-denied"; } ok(denied,"(가짜 규칙) 비구성원 phases 읽기 = permission-denied"); }
{ const st=await as(UNV); ok(st.isMember===false&&st.error.code==="email-unverified","이메일 미확인 계정 → 구성원 아님"); }
{ const st=await as(ADMIN); ok(st.isMember&&st.isAdmin&&st.admins.includes("innolab@example.com"),"관리자 로그인: isAdmin");
  await DB.addMember("ALLOWED_EMAIL_1@example.com","Bob"); await DB.addMember("ALLOWED_EMAIL_2@example.com","Carol");
  const ms=await DB.listMembers(); ok(ms.length===2&&ms[0].email==="ALLOWED_EMAIL_1@example.com"&&ms[0].active===true&&ms[0].addedBy==="innolab@example.com","구성원 추가 (이메일 소문자) · 목록",ms);
  await rejects(DB.addMember("not-an-email"),"invalid","이메일 형식 검사");
  await rejects(DB.removeAdmin("innolab@example.com"),"last-admin","마지막 관리자 지우기 거부");
  ok(FK.get("config/admins").emails.length===1,"관리자 목록 그대로");
  await DB.addAdmin("second@example.com"); ok(FK.get("config/admins").emails.length===2,"관리자 추가"); await DB.removeAdmin("second@example.com"); ok(FK.get("config/admins").emails.length===1,"관리자 제거 (남은 사람 있음)"); }
let labSi, labNaCl;
{ const st=await as(BOB); ok(st.isMember===true&&!st.isAdmin,"Bob: 구성원, 관리자 아님");
  await rejects(DB.addMember("x@y.z","x"),"not-admin","비관리자 addMember → 거부");
  const ws=FK.sdk.firestore; let d=false; try{ await ws.setDoc(ws.doc(ws.getFirestore(),"members","x@y.z"),{name:"x",active:true}); }catch(e){ d=e.code==="permission-denied"; } ok(d,"(가짜 규칙) 비관리자 members 쓰기 거부");
  // 로컬 Si 와 중복 → 돌려줌
  const dup=await DB.addLab(imported["Si_Fd-3m.cif"]); ok(!dup.ok&&dup.reason==="duplicate"&&dup.duplicate.layer==="local","addLab: Local 과 중복 → 쓰지 않고 중복 반환");
  const n0=[...FK.store.keys()].filter(k=>/^phases\/[^/]+$/.test(k)).length; ok(n0===0,"중복이면 Firestore 에 쓰지 않음");
  const a=await DB.promoteLocal(ids["Si_Fd-3m.cif"]); ok(a.ok&&/^lab:/.test(a.rec.id)&&a.rec.version===1&&a.rec.verified===false&&a.rec.hidden===false,"promoteLocal: Local → Lab (version 1, verified false)",a.reason);
  labSi=a.rec.id; const sv=FK.get("phases/"+labSi); ok(sv&&sv.createdBy==="u-bob"&&sv.createdByEmail==="ALLOWED_EMAIL_1@example.com"&&sv.promotedFrom===ids["Si_Fd-3m.cif"]&&sv.createdAt instanceof FK.Ts,"서버 문서: createdBy=uid, serverTimestamp");
  ok(DB.get(ids["Si_Fd-3m.cif"]).promotedTo===labSi,"Local 항목에 promotedTo");
  const b=await DB.addLab(Object.assign({},imported["NaCl_Fm-3m.cif"],{name:"Halite (lab)"}),{force:true}); ok(b.ok,"addLab force (Local 과 중복이어도)"); labNaCl=b.rec.id;
  const dup2=await DB.addLab(imported["NaCl_Fm-3m.cif"]); ok(!dup2.ok&&dup2.reason==="duplicate","addLab 중복 확인 (Lab·Local)");
  const inv=await DB.addLab({name:"bad",cell:{a:-1},sg:{number:225}}); ok(!inv.ok&&inv.reason==="invalid","addLab 검사 오류 → 쓰지 않음");
  // 고치기 → history/1 + version 2
  const cur=DB.get(labNaCl); const up=await DB.updateLab(labNaCl,Object.assign({},cur,{note:"a 다시 잼",cell:Object.assign({},cur.cell,{a:0.5641,b:0.5641,c:0.5641})}));
  ok(up.ok&&up.rec.version===2&&FK.get("phases/"+labNaCl).version===2,"updateLab: version 2");
  const h1=FK.get("phases/"+labNaCl+"/history/1"); ok(h1&&h1.version===1&&!h1.note&&h1.archivedReason==="edit"&&h1.archivedBy==="u-bob","history/1 = 옛 내용");
  // 원자 지우기가 실제로 반영되는지 (update 가 아니라 전체 바꾸기)
  const noAt=Object.assign({},DB.get(labNaCl)); delete noAt.atoms; const up2=await DB.updateLab(labNaCl,noAt); ok(up2.ok&&!FK.get("phases/"+labNaCl).atoms&&FK.get("phases/"+labNaCl+"/history/2").atoms.length>0,"필드 지우기 반영 (원자 삭제), history/2 에 옛 원자");
  await rejects(DB.setVerified(labNaCl,true),"not-admin","비관리자 setVerified → 거부");
  // 규칙 직접: verified 바꾸기·version 건너뛰기·createdBy 바꾸기·지우기·history 덮어쓰기
  const F_=FK.sdk.firestore, D=F_.getFirestore(), ref=F_.doc(D,"phases",labNaCl), cur3=FK.get("phases/"+labNaCl);
  const tryW=async(fn,m)=>{ let c=null; try{ await fn(); }catch(e){ c=e.code; } ok(c==="permission-denied",m,c); };
  await tryW(()=>F_.setDoc(ref,Object.assign({},cur3,{version:cur3.version+1,verified:true})),"(규칙) 비관리자가 verified 바꾸기 거부");
  await tryW(()=>F_.setDoc(ref,Object.assign({},cur3,{version:cur3.version+2})),"(규칙) version+2 거부");
  await tryW(()=>F_.setDoc(ref,Object.assign({},cur3,{version:cur3.version})),"(규칙) version 그대로 거부");
  await tryW(()=>F_.setDoc(ref,Object.assign({},cur3,{version:cur3.version+1,createdBy:"u-carol"})),"(규칙) createdBy 바꾸기 거부");
  await tryW(()=>F_.deleteDoc(ref),"(규칙) 지우기 거부");
  await tryW(()=>F_.setDoc(F_.doc(D,"phases",labNaCl,"history","1"),{x:1}),"(규칙) history 덮어쓰기 거부");
  await tryW(()=>F_.setDoc(F_.doc(D,"phases","lab:x"),Object.assign({},cur3,{version:2})),"(규칙) 새 문서 version≠1 거부");
  await tryW(()=>F_.setDoc(F_.doc(D,"phases","lab:y"),Object.assign({},cur3,{version:1,createdBy:"u-carol"})),"(규칙) 새 문서 createdBy≠uid 거부");
  // 보정
  const cl=await DB.saveCalibration({name:"Talos 200kV CL 300",microscope:"Talos F200X",kV:200,cameraLength:300,mode:"SAED",scale:0.01234,ellip:1.003,phi:12,sigmaRel:0.002,standard:labSi,sourceFile:"Au_std.dm3"},{lab:true});
  ok(/^lab:/.test(cl.id)&&FK.get("calibrations/"+cl.id).createdBy==="u-bob"&&cl.standardRef&&cl.standardRef.id===labSi&&cl.standardRef.version===1,"Lab 보정 저장 (+ 표준 상 refOf)");
  const lc=await DB.saveCalibration({name:"local cal",scale:0.02}); ok(lc.layer==="local"&&DB.listCalibrations().length===2,"Local 보정 + 목록");
  await rejects(DB.saveCalibration({name:"x",scale:0}),"invalid","보정 scale 검사");
  window_dummy=0; }
var window_dummy;
let carolFirst;
{ const st=await as(CAROL); ok(st.isMember&&!st.isAdmin,"Carol: 구성원");
  const s=await DB.sync(); ok(s.phases===2&&s.calibrations===1&&DB.firebaseStatus().lastSync>0,"sync: phases 2 · calibrations 1 · lastSync",s);
  ok(DB.list({layers:["lab"]}).length===2&&DB.listCalibrations({layers:["lab"]}).length===1,"Lab 캐시 목록");
  await rejects(DB.updateLab(labNaCl,DB.get(labNaCl)),"not-owner","남의 Lab 항목 고치기 → 거부 (클라이언트)");
  await rejects(DB.hideLab(labNaCl),"not-owner","남의 항목 숨기기 → 거부");
  const F_=FK.sdk.firestore, D=F_.getFirestore(), cur=FK.get("phases/"+labNaCl); let c=null; try{ await F_.setDoc(F_.doc(D,"phases",labNaCl),Object.assign({},cur,{version:cur.version+1,note:"carol"})); }catch(e){ c=e.code; } ok(c==="permission-denied","(규칙) 남의 항목 update 거부");
  const cal=DB.listCalibrations({layers:["lab"]})[0]; await rejects(DB.saveCalibration(cal,{lab:true}),"not-owner","남의 보정 고치기 → 거부");
  c=null; try{ await F_.setDoc(F_.doc(D,"calibrations",cal.id),{name:"hack",createdBy:"u-bob"}); }catch(e){ c=e.code; } ok(c==="permission-denied","(규칙) 남의 보정 update 거부");
  const h=await DB.labHistory(labNaCl); ok(h.length===2&&h[0].version===2&&h[1].version===1,"labHistory: 판본 2·1");
  carolFirst=DB.refOf(labNaCl); ok(carolFirst.layer==="lab"&&carolFirst.version===3&&carolFirst.projectId==="inno-analytics-base","refOf (lab): version 3, projectId");
  FK.offline=true; const e=await rejects(DB.sync(),"unavailable","오프라인 sync → 오류"); ok(e&&/캐시/.test(e.message)&&DB.list({layers:["lab"]}).length===2,"오프라인: 캐시 그대로"); FK.offline=false; }
{ const st=await as(ADMIN); await DB.sync();
  const v=await DB.setVerified(labNaCl,true); ok(v.verified&&v.version===4&&v.verifiedVersion===4&&FK.get("phases/"+labNaCl).verified===true,"관리자 verify → version 4, verifiedVersion 4");
  ok(DB.verifyState(labNaCl).state==="verified","verifyState: verified");
  await DB.setMemberActive("ALLOWED_EMAIL_2@example.com",false); ok(FK.get("members/carol@kaist.ac.kr").active===false,"관리자: Carol 비활성"); }
{ const st=await as(BOB); await DB.sync();
  const cur=DB.get(labNaCl); const up=await DB.updateLab(labNaCl,Object.assign({},cur,{note:"verified 뒤 편집"}));
  ok(up.ok&&up.rec.version===5&&up.rec.verified===true&&up.rec.verifiedVersion===4,"주인 편집: verified 유지(규칙), verifiedVersion 4");
  ok(DB.verifyState(labNaCl).state==="edited","verifyState: 검증 뒤 편집됨");
  const hd=await DB.hideLab(labNaCl); ok(hd.hidden===true&&hd.version===6&&FK.get("phases/"+labNaCl).hidden===true&&FK.get("phases/"+labNaCl+"/history/5"),"hideLab → hidden, version 6, history/5");
  ok(!DB.list({layers:["lab"]}).some(r=>r.id===labNaCl)&&DB.list({layers:["lab"],includeHidden:true}).some(r=>r.id===labNaCl),"숨긴 항목: 기본 목록에서 빠짐, includeHidden 이면 보임");
  ok(DB.get(labNaCl)&&DB.phase(labNaCl),"숨겨도 get/phase 는 됨 (재현)");
  const old=await DB.getVersion(carolFirst); ok(old&&old.exact&&old.rec.version===3,"getVersion(lab v3) → history 에서",old&&old.rec.version);
  const un=await DB.hideLab(labNaCl,false); ok(un.hidden===false&&un.version===7,"숨김 해제"); }
{ const st=await as(CAROL); ok(st.isMember===false&&/구성원으로 등록되지 않았습니다/.test(st.error.message),"비활성 구성원 → 안내");
  await rejects(DB.sync(),"not-member","비활성 구성원 sync 거부");
  ok(DB.list({layers:["lab"]}).length===2,"그래도 캐시는 남음 (오프라인 사용)"); }
{ await DB._reset(); DB._setIDB(idb); DB._setSDK(FK.sdk); DB._setEnv({protocol:"https:"}); await DB.init();
  await DB.setConfig({apiKey:"x",authDomain:"other.firebaseapp.com",projectId:"some-other-project",appId:"1:2:web:3"});
  FK.loginAs(BOB); await DB.connect(); await flush(); const st=DB.firebaseStatus();
  ok(!st.configOk&&st.projectId==="some-other-project"&&st.configOverridden,"다른 projectId → configOk false");
  ok(DB.list({layers:["lab"]}).length===0&&st.ignoredCache&&st.ignoredCache.cacheProject==="inno-analytics-base","다른 프로젝트의 Lab 캐시를 섞지 않음");
  await rejects(DB.addLab({name:"Fe",formula:"Fe",cell:{a:0.2866},sg:{number:229},source:{kind:"manual",note:"t"}},{force:true}),"wrong-project","다른 projectId 로는 쓰기 거부");
  await DB.setConfig(null); await DB._reset(); DB._setIDB(idb); DB._setSDK(FK.sdk); await DB.init(); ok(DB.list({layers:["lab"]}).length===2,"기본 config 로 돌아오면 캐시 다시 보임"); }
{ // 이벤트
  let evs=[]; const off=DB.onChange(e=>evs.push(e.type)); await DB.saveLocal({name:"Fe bcc",formula:"Fe",cell:{a:0.2866},sg:{number:229},source:{kind:"manual",note:"t"}}); off();
  ok(evs.includes("local"),"onChange: local 이벤트"); }
console.log("== Firebase 10.12.2 실제 SDK 에 쓰는 이름이 있는지 (npm firebase@10.12.2)");
try{ const [a,au,f]=await Promise.all([import("firebase/app"),import("firebase/auth"),import("firebase/firestore")]);
  const need={app:["initializeApp","getApps"],auth:["getAuth","GoogleAuthProvider","signInWithPopup","signInWithRedirect","getRedirectResult","onAuthStateChanged","signOut"],firestore:["getFirestore","doc","collection","getDoc","getDocFromServer","getDocs","getDocsFromServer","setDoc","updateDoc","runTransaction","serverTimestamp"]};
  const mods={app:a,auth:au,firestore:f}; for(const k of Object.keys(need)){ const miss=need[k].filter(n=>typeof mods[k][n]!=="function"); ok(!miss.length,`firebase/${k}: ${need[k].length} 함수`,miss); }
  const fbv=require("firebase/package.json").version; ok(fbv==="10.12.2","npm firebase 판 = "+fbv);
}catch(e){ ok(false,"firebase npm 패키지 불러오기",e.message); }

console.log(`\n${pass} 통과, ${fail} 실패`); if(fail){ console.log(fails.map(f=>"  ✗ "+f).join("\n")); process.exit(1); }
process.exit(0);
})().catch(e=>{ console.error(e); process.exit(2); });
