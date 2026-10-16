"""Phase DB 화면 브라우저 실측 (Playwright · Chromium, 네트워크 없음)
사용: python3 test/ui_db.py [그림 폴더] [앱 빌드 html(선택, 겹침 시험용)]
- 시험 페이지를 만든다: 앱 CSS(빌드 html 의 기본·터치 style) + diff_db.css + diff_sg/scat/core/db/db_ui.js
- Firebase 는 test/fake_firebase.js (규칙 집행 가짜) 로 바꿔 끼우고, 모든 외부 요청은 막는다 (firebasejs 요청 0 개 확인)
- 데스크톱 1400×900 · 휴대폰 390×844 / 360×780 (pointer:coarse): 가로 스크롤 없음, 단추 ≥ 40 px,
  추가 대화상자 검사(오류·경고·중복), 목록·원소 거르기, 링 표, 오프라인 상태 줄, 콘솔 오류 0."""
import sys, os, re, json, glob, tempfile
from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src')
OUT = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix='ui_db_')
BUILD = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, '..', '..', '..', 'INNO Analytics.html')   # 앱 CSS 를 가져올 html
os.makedirs(OUT, exist_ok=True)
rd = lambda p: open(p, encoding='utf-8').read()

host_css = ''
if BUILD and os.path.exists(BUILD):
    h = rd(BUILD)
    m = re.search(r'<style>(.*?)</style>', h, re.S); host_css += m.group(1) if m else ''
    m = re.search(r'<style id="inno-touch-css">(.*?)</style>', h, re.S); host_css += m.group(1) if m else ''
fake = rd(os.path.join(HERE, 'fake_firebase.js'))
pack_js = os.path.join(SRC, 'diff_pack.js')
mini_pack = """window.DIFF_PACK={version:"ui-test",phases:[
 {id:"pack:Au",name:"Gold",formula:"Au",cell:{a:0.40782,b:0.40782,c:0.40782,al:90,be:90,ga:90},sg:{number:225,symbol:"Fm-3m"},atoms:[{label:"Au",el:"Au",x:0,y:0,z:0,occ:1,B:0.6}],source:{kind:"COD",codId:"9008463"},verified:true},
 {id:"pack:Si",name:"Silicon",formula:"Si",cell:{a:0.54309,b:0.54309,c:0.54309,al:90,be:90,ga:90},sg:{number:227,symbol:"Fd-3m:2"},source:{kind:"manual",note:"no atoms"}}]};"""
scripts = ''.join(f'<script>\n{rd(os.path.join(SRC, n))}\n</script>\n' for n in ['diff_sg.js', 'diff_scat.js', 'diff_core.js'])
scripts += (f'<script>\n{rd(pack_js)}\n</script>\n' if os.path.exists(pack_js) else f'<script>{mini_pack}</script>\n')
scripts += f'<script>(function(){{ const module={{exports:{{}}}}; {fake}\n window.FAKEFB=module.exports; }})();</script>\n'
scripts += ''.join(f'<script>\n{rd(os.path.join(SRC, n))}\n</script>\n' for n in ['diff_db.js', 'diff_db_ui.js'])
page_html = f"""<!doctype html><html lang="ko" data-theme="light"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Phase DB test</title><style>{host_css}</style><style>{rd(os.path.join(SRC, 'diff_db.css'))}</style></head>
<body><div id="app"><div id="root"></div></div><div id="toasts"></div>
<script>window.__toasts=[]; function toast(m,k){{ window.__toasts.push([k||"",m]); const t=document.createElement("div"); t.className="toast "+(k||""); t.textContent=m; document.getElementById("toasts").appendChild(t); setTimeout(()=>t.remove(),2500); }}</script>
{scripts}
<script>
/* 시험 설정: ?env=https → 가짜 Firebase 로 웹 주소처럼, 아니면 file:// 그대로 (오프라인 상태 줄) */
const P=new URLSearchParams(location.search);
window.FK=FAKEFB.makeFake(); DIFF_DB._setSDK(FK.sdk);
FK.seed("config/admins",{{emails:["innolab@example.com"]}});
if(P.get("env")==="https") DIFF_DB._setEnv({{protocol:"https:"}});
DIFF_DB_UI.render(document.getElementById("root"),{{onUse:(ids,refs)=>{{ window.__used={{ids,refs}}; }}}});
</script></body></html>"""
test_html = os.path.join(OUT, 'ui_db_test.html')
open(test_html, 'w', encoding='utf-8').write(page_html)
CIFS = sorted(glob.glob(os.path.join(HERE, 'cif', '*.cif')))

res = {}; errs = []; fails = []; netreq = []
def ok(c, m, d=None):
    print(('  ✓ ' if c else '  ✗ ') + m + ('' if c or d is None else ' — ' + str(d)))
    if not c: fails.append(m)

def guard(pg, tag):
    pg.on('console', lambda m: errs.append(f'{tag} console.{m.type}: {m.text}') if m.type == 'error' else None)
    pg.on('pageerror', lambda e: errs.append(f'{tag} pageerror: {e}'))
    pg.on('dialog', lambda d: d.accept())
    def route(r):
        u = r.request.url
        if u.startswith('file:') or u.startswith('data:') or u.startswith('blob:'): return r.continue_()
        netreq.append(u)
        if 'fonts.g' in u: return r.fulfill(status=200, body='', content_type='text/css')   # 앱의 Google Fonts — 빈 응답
        return r.abort()
    pg.route('**/*', route)

def layout(pg, tag, coarse):
    """가로 스크롤·단추 크기"""
    m = pg.evaluate("""(coarse)=>{ const de=document.documentElement; const vw=window.innerWidth;
      const vis=e=>{ const r=e.getBoundingClientRect(); const s=getComputedStyle(e); return r.width>0&&r.height>0&&s.visibility!=="hidden"; };
      // 표 안(.tblwrap)은 자기 상자 안에서 스크롤 — 그 밖의 요소가 화면 밖으로 나가면 안 됨
      const over=[...document.querySelectorAll('.ddb *, .ddb-modal *')].filter(e=>vis(e)&&!e.closest('.tblwrap')).filter(e=>{ const r=e.getBoundingClientRect(); return r.right>vw+1; }).slice(0,5).map(e=>e.tagName+'.'+e.className+'#'+e.id+' '+Math.round(e.getBoundingClientRect().right));
      const small=[...document.querySelectorAll('.ddb button, .ddb .btn, .ddb-modal button, .ddb-modal .btn')].filter(vis).filter(e=>{ const r=e.getBoundingClientRect(); return r.height<40-0.5||(e.closest('.chip')&&r.width<40-0.5); }).map(e=>(e.textContent||e.getAttribute('aria-label')||'').trim().slice(0,24)+' '+Math.round(e.getBoundingClientRect().height));
      const wraps=[...document.querySelectorAll('.ddb .tblwrap')].filter(vis).map(w=>({sw:w.scrollWidth,cw:w.clientWidth,ov:getComputedStyle(w).overflowX}));
      return {sw:de.scrollWidth,vw,coarse:matchMedia('(pointer:coarse)').matches,over,small,wraps}; }""", coarse)
    ok(m['sw'] <= m['vw'], f'{tag}: 가로 스크롤 없음 (scrollWidth {m["sw"]} ≤ {m["vw"]})', m['over'])
    ok(not m['over'], f'{tag}: 화면 밖으로 나간 요소 없음', m['over'])
    if coarse:
        ok(m['coarse'], f'{tag}: pointer:coarse 적용')
        ok(not m['small'], f'{tag}: 보이는 단추 모두 ≥ 40 px', m['small'][:8])
    ok(all(w['ov'] in ('auto', 'scroll') for w in m['wraps']), f'{tag}: 표는 자기 상자 안에서 스크롤', m['wraps'])
    return m

with sync_playwright() as p:
    b = p.chromium.launch()
    # ================= 데스크톱 · file:// (오프라인, Firebase 없음) =================
    pg = b.new_page(viewport={'width': 1400, 'height': 900}); guard(pg, 'desk')
    pg.goto('file://' + test_html); pg.wait_for_selector('.ddb #ddbStLine'); pg.wait_for_timeout(300)
    st = pg.inner_text('#ddbStatus')
    ok('로그인하지 않음' in st and 'inno-analytics-base' in st, '오프라인 상태 줄: 로그인하지 않음 · projectId', st[:200])
    ok('웹 주소' in st and pg.is_disabled('#ddbStatus button[data-a=signin]'), 'file:// 안내 + Sign in 비활성')
    ok(pg.evaluate('!DIFF_DB.firebaseStatus().sdkLoaded'), 'Firebase SDK 를 부르지 않음 (오프라인)')
    n_pack = pg.evaluate('DIFF_DB.counts().pack'); ok(n_pack >= 1, f'Built-in {n_pack} 개 목록')
    pg.screenshot(path=os.path.join(OUT, 'd01_empty.png'))
    # --- CIF 로 추가 ---
    pg.click('#ddbAdd'); pg.wait_for_selector('.ddb-modal #ddbCifF', state='attached')
    pg.set_input_files('#ddbCifF', CIFS); pg.wait_for_function('document.querySelectorAll(".ddb-modal .blk").length>=10')
    nb = pg.evaluate('document.querySelectorAll(".ddb-modal .blk").length'); nsel = pg.evaluate('document.querySelectorAll(".ddb-modal .blk input[data-u]:checked").length')
    ok(nb == len(CIFS), f'CIF {len(CIFS)} 개 → 블록 {nb} 개, 선택 {nsel}')
    # 기대 중복 수: DB 를 거치지 않고 DIFF_PACK 과 직접 DC.sameRecord (pack 내용이 바뀌어도 시험이 맞게)
    exp_dup = pg.evaluate('''(texts)=>texts.filter(t=>{ const r=DIFF_CORE.recordFromPhase(DIFF_CORE.parseCIF(t)[0]); return (window.DIFF_PACK?DIFF_PACK.phases:[]).some(p=>DIFF_CORE.sameRecord(r,p)); }).length''', [rd(f) for f in CIFS])
    n_dupnote = pg.evaluate('[...document.querySelectorAll(".ddb-modal .blk .note.warn")].filter(n=>/중복 후보/.test(n.textContent)).length')
    ok(nb - nsel == exp_dup and n_dupnote == exp_dup, f'Built-in 과 중복 {exp_dup} 개 → 체크 해제 + 중복 안내 ({nb - nsel} · {n_dupnote})')
    # 중복도 체크해서 그래도 저장 (Local 에 10 개 모두 — 이후 시험이 pack 내용과 무관하게)
    for i in pg.evaluate('[...document.querySelectorAll(".ddb-modal .blk input[data-u]:not(:checked):not(:disabled)")].map(x=>x.dataset.u)'):
        pg.check(f'.ddb-modal .blk input[data-u="{i}"]')
    ok(pg.inner_text('#ddbSave').startswith(f'Save {nb}'), f'모두 체크 → {pg.inner_text("#ddbSave")}')
    pg.screenshot(path=os.path.join(OUT, 'd02_cif_blocks.png'))
    pg.click('#ddbSave'); pg.wait_for_timeout(400)
    n_loc = pg.evaluate('DIFF_DB.counts().local'); ok(n_loc == nb, f'Local 에 {n_loc} 개 저장 (중복 {exp_dup} 개는 체크해서 그래도 저장)')
    ok(pg.evaluate('!document.querySelector(".ddb-modal")'), '저장 뒤 대화상자 닫힘')
    if pg.evaluate('!!document.querySelector(".ddb-modal")'): pg.click('.ddb-modal .closex')
    # --- 매개변수로: 검사 ---
    pg.click('#ddbAdd'); pg.click('#ddbMode button[data-m=params]')
    pg.fill('.ddb-modal input[data-k=name]', 'Test cubic'); pg.fill('.ddb-modal input[data-k=formula]', 'Fe')
    pg.fill('.ddb-modal input[data-k=sgn]', '229'); pg.fill('.ddb-modal input[data-k=a]', '2.866')
    pg.wait_for_timeout(250)
    sysv = pg.input_value('#ddbSys'); sym = pg.input_value('.ddb-modal input[data-k=sgs]')
    ok(sysv == 'Cubic' and sym.replace(' ', '') in ('Im-3m',), f'SG 229 → System {sysv}, 기호 {sym}')
    ok(pg.is_disabled('.ddb-modal input[data-k=b]') and pg.is_disabled('.ddb-modal input[data-k=ga]') and pg.input_value('.ddb-modal input[data-k=c]') == '2.866', '입방정: b·c·각 잠김, c = a')
    warn = pg.inner_text('#ddbVal'); ok('Å' in warn and not pg.is_disabled('#ddbSave'), 'a = 2.866 nm → Å 경고 (저장은 가능)', warn[:160])
    pg.fill('.ddb-modal input[data-k=a]', '0.2866'); pg.fill('.ddb-modal input[data-k=sgn]', '300'); pg.wait_for_timeout(250)
    ok(pg.is_visible('#ddbErr') and pg.is_disabled('#ddbSave'), 'SG 300 → 오류, Save 비활성', pg.inner_text('#ddbVal')[:160])
    pg.fill('.ddb-modal input[data-k=sgn]', '14'); pg.wait_for_timeout(200)
    ok(pg.input_value('#ddbSys') == 'Monoclinic' and not pg.is_disabled('.ddb-modal input[data-k=be]') and pg.is_disabled('.ddb-modal input[data-k=al]'), '단사정: β 만 풀림')
    # 중복: NaCl 과 같은 값
    nacl = pg.evaluate('DIFF_DB.list({query:"NaCl"}).concat(DIFF_DB.list({query:"Cl Na"}))[0]')
    pg.fill('.ddb-modal input[data-k=formula]', 'NaCl'); pg.fill('.ddb-modal input[data-k=sgn]', '225'); pg.fill('.ddb-modal input[data-k=a]', str(nacl['cell']['a']))
    pg.wait_for_timeout(350)
    ok(pg.is_visible('#ddbDup') and pg.is_disabled('#ddbSave'), '같은 상(NaCl) → 중복 경고, Save 비활성', pg.inner_text('#ddbVal')[:200])
    pg.check('#ddbDupOk'); pg.wait_for_timeout(50)
    ok(not pg.is_disabled('#ddbSave') and 'anyway' in pg.inner_text('#ddbSave'), 'Save anyway 체크 → 저장 가능')
    pg.click('.ddb-modal details.sec summary'); pg.click('#ddbAddAt'); pg.fill('.ddb-modal input[data-at="0"][data-f=label]', 'Na1'); pg.fill('.ddb-modal input[data-at="0"][data-f=el]', 'Na')
    pg.screenshot(path=os.path.join(OUT, 'd03_params_dup.png'))
    pg.click('#ddbSave'); pg.wait_for_timeout(300)
    ok(pg.evaluate('DIFF_DB.counts().local') == n_loc + 1 and not pg.evaluate('!!document.querySelector(".ddb-modal")'), '중복 확인 후 저장 · 대화상자 닫힘')
    # --- 원소 거르기 ---
    pg.fill('#ddbEl', 'Na Y F'); pg.press('#ddbEl', 'Enter'); pg.wait_for_timeout(200)
    rows = pg.evaluate('[...document.querySelectorAll(".ddb-tbl tbody tr")].map(t=>t.children[2].textContent)')
    ok(len(rows) >= 1 and all(set(re.findall(r'[A-Z][a-z]?', f)) <= {'Na', 'Y', 'F'} for f in rows), f'원소 ⊆ {{Na,Y,F}} → {rows}')
    pg.select_option('#ddbElMode', 'all'); pg.wait_for_timeout(150)
    rows2 = pg.evaluate('[...document.querySelectorAll(".ddb-tbl tbody tr")].map(t=>t.children[2].textContent)')
    ok(all({'Na', 'Y', 'F'} <= set(re.findall(r'[A-Z][a-z]?', f)) for f in rows2), f'Contains all → {rows2}')
    pg.click('#ddbElClr'); pg.select_option('#ddbElMode', 'subset')
    pg.fill('#ddbQ', 'TiO2'); pg.wait_for_timeout(300)
    exp_ti = 1 + pg.evaluate('''(window.DIFF_PACK?DIFF_PACK.phases:[]).filter(p=>{ const f=DIFF_CORE.parseFormula(p.formula); const k=Object.keys(f).sort(); return k.join()==="O,Ti"&&Math.abs(f.O/f.Ti-2)<1e-6; }).length''')
    n_ti = pg.evaluate('document.querySelectorAll(".ddb-tbl tbody tr").length')
    ok(n_ti == exp_ti, f'검색 TiO2 → {n_ti} 줄 (Local 1 + Built-in {exp_ti - 1})')
    pg.click('.ddb-layers button[data-l=pack]'); pg.wait_for_timeout(150)
    ok(pg.evaluate('document.querySelectorAll(".ddb-tbl tbody tr").length') == 1, '검색 TiO2, Built-in 끔 → Local 1 줄 (Ti2 O4)')
    pg.click('.ddb-layers button[data-l=pack]')
    pg.fill('#ddbQ', ''); pg.wait_for_timeout(300)
    # --- 층 거르기 ---
    pg.click('.ddb-layers button[data-l=pack]'); pg.wait_for_timeout(100)
    lays = pg.evaluate('[...new Set([...document.querySelectorAll(".ddb-tbl tbody .lay")].map(x=>x.textContent))]')
    ok(lays == ['Local'], f'Built-in 끄기 → {lays}'); pg.click('.ddb-layers button[data-l=pack]')
    # --- 상세 · 링 표 ---
    nid = pg.evaluate('DIFF_DB.list({layers:["local"]}).find(r=>r.formula.replace(/\\s/g,"")==="Na4Cl4").id')
    pg.click(f'.ddb-tbl tbody tr[data-id="{nid}"] td.nm', timeout=3000)
    pg.wait_for_selector('#ddbRings tbody tr')
    r1 = pg.evaluate('[...document.querySelectorAll("#ddbRings tbody tr")].slice(0,3).map(t=>[...t.children].map(c=>c.textContent))')
    nr = pg.evaluate('document.querySelectorAll("#ddbRings tbody tr").length')
    ok(r1[0][3] == '{111}' and abs(float(r1[0][1]) - 0.3256) < 0.002 and r1[1][3] == '{200}' and r1[0][4] == '8', f'NaCl 링: {r1[0]} · {r1[1]} ({nr} 줄)')
    ok(all(float(r[2]) <= 10.0001 for r in pg.evaluate('[...document.querySelectorAll("#ddbRings tbody tr")].map(t=>[...t.children].map(c=>c.textContent))')), 'g ≤ 10 1/nm')
    det = pg.inner_text('#ddbDetail'); ok('Fm-3m' in det.replace(' ', '') and 'Cubic' in det, '상세: 공간군·결정계')
    pg.screenshot(path=os.path.join(OUT, 'd04_detail.png'), full_page=True)
    # Built-in(원자 없음) → 세기 없음 안내
    # 원자 없는 상 (시험에서 직접 만듦 — pack 내용과 무관)
    wid = pg.evaluate('DIFF_DB.saveLocal({name:"No-atom W",formula:"W",cell:{a:0.31652},sg:{number:229},source:{kind:"manual",note:"ui test"}}).then(r=>r.rec.id)')
    pg.evaluate(f'DIFF_DB_UI._state().sel="{wid}"; DIFF_DB_UI.rerender()'); pg.wait_for_timeout(200)
    wr = pg.evaluate('[...document.querySelectorAll("#ddbRings tbody tr")].slice(0,2).map(t=>t.children[3].textContent)')
    ok('원자 좌표가 없어' in pg.inner_text('#ddbDetail') and wr == ['{110}', '{200}'], f'원자 없는 상: 세기 없음 안내 + 링(소멸칙, bcc → {wr})')
    # 편집 → version 2, history
    lid = pg.evaluate('DIFF_DB.list({layers:["local"]})[0].id')
    pg.evaluate(f'DIFF_DB_UI._state().sel="{lid}"; DIFF_DB_UI.rerender()')
    pg.click('#ddbDetail button[data-a=edit]'); pg.fill('.ddb-modal input[data-k=note]', 'edited in UI test'); pg.wait_for_timeout(200); pg.click('#ddbSave'); pg.wait_for_timeout(300)
    ok(pg.evaluate(f'DIFF_DB.get("{lid}").version') == 2 and pg.evaluate('document.querySelectorAll("#ddbHist tbody tr").length') == 1, 'Edit → v2, History 1 줄')
    # 선택 → onUse(refs)
    pg.check('.ddb-tbl tbody tr:nth-child(1) input[data-c]'); pg.check('.ddb-tbl tbody tr:nth-child(2) input[data-c]'); pg.click('#ddbUse')
    used = pg.evaluate('window.__used'); ok(used and len(used['refs']) == 2 and all('version' in r and 'layer' in r and 'hash' in r for r in used['refs']), 'Use selected → refOf 2 개', used)
    layout(pg, 'desktop', False)
    pg.evaluate('document.documentElement.setAttribute("data-theme","dark")'); pg.wait_for_timeout(100)
    pg.screenshot(path=os.path.join(OUT, 'd05_dark.png'))
    pg.evaluate('document.documentElement.setAttribute("data-theme","light")')
    backup = pg.evaluate('DIFF_DB.exportJSON()')
    pg.close()

    # ================= Lab 흐름 (가짜 Firebase, https 처럼) =================
    pg = b.new_page(viewport={'width': 1400, 'height': 900}); guard(pg, 'lab')
    pg.goto('file://' + test_html + '?env=https'); pg.wait_for_selector('.ddb #ddbStLine'); pg.wait_for_timeout(300)
    ok(not pg.is_disabled('#ddbStatus button[data-a=signin]'), 'https: Sign in 활성')
    pg.evaluate('FK.nextUser=FK.user("u-admin","innolab@example.com")'); pg.click('#ddbStatus button[data-a=signin]')
    pg.wait_for_function('DIFF_DB.firebaseStatus().isAdmin===true'); pg.wait_for_timeout(200)
    st = pg.inner_text('#ddbStatus')
    ok('innolab@example.com' in st and 'Admin' in st and 'inno-analytics-base' in st, '로그인: 이메일 · Admin · projectId', st[:200])
    ok(pg.is_visible('#ddbMem'), '관리자: Members 패널')
    pg.click('#ddbMemD summary'); pg.fill('#ddbMemE', 'ALLOWED_EMAIL_1@example.com'); pg.fill('#ddbMemN', 'Bob'); pg.click('#ddbMem button[data-m=add]')
    pg.wait_for_selector('#ddbMemT tbody tr:has-text("ALLOWED_EMAIL_1@example.com")'); ok(True, '구성원 추가 → 목록')
    ok(pg.is_disabled('#ddbMem button[data-rmadm]'), '마지막 관리자 × 비활성')
    pg.screenshot(path=os.path.join(OUT, 'l01_admin.png'), full_page=True)
    exp_bk = sum(1 for r in json.loads(backup)['phases'] if r.get('layer') == 'local')
    rep = pg.evaluate(f'DIFF_DB.importJSON({json.dumps(backup)})'); ok(rep['added'] == exp_bk and exp_bk >= len(CIFS), f'백업 가져오기 {rep["added"]} / {exp_bk} 개')
    lid = pg.evaluate('DIFF_DB.list({layers:["local"]}).find(r=>/Cl/.test(r.formula)).id')
    pg.evaluate(f'DIFF_DB_UI._state().sel="{lid}"; DIFF_DB_UI.rerender()')
    pg.click('#ddbDetail button[data-a=promote]'); pg.wait_for_timeout(400)   # 다른 Local NaCl 과 중복 → confirm(자동 수락) → 올림
    lab = pg.evaluate('DIFF_DB.list({layers:["lab"]})'); ok(len(lab) == 1 and lab[0]['version'] == 1, 'Copy to Lab → Lab v1', lab and lab[0].get('id'))
    pg.click('#ddbDetail button[data-a=verify]'); pg.wait_for_timeout(300)
    ok(pg.evaluate(f'DIFF_DB.get("{lab[0]["id"]}").verified') and 'Verified' in pg.inner_text('#ddbDetail'), 'Verify (관리자) → Verified 배지')
    pg.click('#ddbDetail button[data-a=hist]'); pg.wait_for_selector('#ddbHist tbody tr'); ok(pg.evaluate('document.querySelectorAll("#ddbHist tbody tr").length') == 1, 'Load history → v1 한 줄')
    pg.click('#ddbStatus button[data-a=sync]'); pg.wait_for_timeout(300); ok('just now' in pg.inner_text('#ddbStatus'), 'Sync → Last sync just now')
    pg.screenshot(path=os.path.join(OUT, 'l02_lab_detail.png'), full_page=True)
    # 비구성원
    pg.click('#ddbStatus button[data-a=signout]'); pg.evaluate('FK.nextUser=FK.user("u-eve","ALLOWED_EMAIL_3@example.com")'); pg.click('#ddbStatus button[data-a=signin]'); pg.wait_for_timeout(400)
    st = pg.inner_text('#ddbStatus'); ok('구성원으로 등록되지 않았습니다' in st and 'ALLOWED_EMAIL_3@example.com' in st and 'Not a member' in st, '비구성원 안내 + 이메일', st[:240])
    # 다른 projectId
    pg.evaluate('DIFF_DB.setConfig({apiKey:"x",projectId:"other-proj",authDomain:"x",appId:"y"})'); pg.wait_for_timeout(200)
    ok(pg.is_visible('#ddbStatus .note.danger') and 'other-proj' in pg.inner_text('#ddbStatus'), '다른 projectId → 위험 경고')
    pg.screenshot(path=os.path.join(OUT, 'l03_wrong_project.png'))
    pg.click('#ddbStatus button[data-a=resetcfg]'); pg.wait_for_timeout(200); ok(not pg.is_visible('#ddbStatus .note.danger:has-text("프로젝트가 다릅니다")'), 'Reset config → 경고 사라짐')
    pg.close()

    # ================= 휴대폰 (pointer:coarse) =================
    for (w, h) in [(390, 844), (360, 780)]:
        ctx = b.new_context(viewport={'width': w, 'height': h}, is_mobile=True, has_touch=True, device_scale_factor=2)
        pg = ctx.new_page(); guard(pg, f'm{w}')
        pg.goto('file://' + test_html + '?env=https'); pg.wait_for_selector('.ddb #ddbStLine'); pg.wait_for_timeout(300)
        pg.evaluate(f'DIFF_DB.importJSON({json.dumps(backup)})'); pg.wait_for_timeout(300)
        pg.fill('#ddbEl', 'O'); pg.press('#ddbEl', 'Enter'); pg.select_option('#ddbElMode', 'all'); pg.wait_for_timeout(200)
        pg.tap('.ddb-tbl tbody tr:nth-child(1) td.nm'); pg.wait_for_selector('#ddbRings tbody tr')
        layout(pg, f'{w}px list+detail', True)
        pg.screenshot(path=os.path.join(OUT, f'm{w}_01_list.png'))   # full_page 는 pointer:coarse 흉내를 풀어 버림
        pg.tap('#ddbAdd'); pg.tap('#ddbMode button[data-m=params]'); pg.fill('.ddb-modal input[data-k=sgn]', '194'); pg.fill('.ddb-modal input[data-k=a]', '0.3'); pg.tap('.ddb-modal details.sec summary'); pg.tap('#ddbAddAt'); pg.wait_for_timeout(200)
        ok(pg.input_value('#ddbSys') == 'Hexagonal' and pg.input_value('.ddb-modal input[data-k=ga]') == '120', f'{w}px: 육방정 γ = 120 자동')
        layout(pg, f'{w}px add dialog', True)
        mw = pg.evaluate('(()=>{ const m=document.querySelector(".ddb-modal"); return [m.scrollWidth,m.clientWidth,m.getBoundingClientRect().right,innerWidth]; })()')
        ok(mw[0] <= mw[1] + 1 and mw[2] <= mw[3], f'{w}px: 대화상자 가로 넘침 없음 {mw}')
        pg.screenshot(path=os.path.join(OUT, f'm{w}_02_dialog.png'))
        pg.tap('.ddb-modal .closex')
        pg.evaluate('FK.nextUser=FK.user("u-admin","innolab@example.com")'); pg.tap('#ddbStatus button[data-a=signin]'); pg.wait_for_function('DIFF_DB.firebaseStatus().isAdmin===true'); pg.wait_for_timeout(200)
        pg.tap('#ddbMemD summary'); pg.tap('#ddbMem button[data-m=load]'); pg.wait_for_timeout(300)
        layout(pg, f'{w}px signed in + members', True)
        pg.screenshot(path=os.path.join(OUT, f'm{w}_03_signed.png'))
        pg.evaluate('document.querySelector("#ddbDetail").scrollIntoView()'); pg.screenshot(path=os.path.join(OUT, f'm{w}_04_detail.png'))
        ctx.close()

    # ================= 앱 빌드에 넣어 겹침 확인 (선택) =================
    if BUILD and os.path.exists(BUILD):
        pg = b.new_page(viewport={'width': 1400, 'height': 900}); guard(pg, 'build')
        pg.goto('file://' + os.path.abspath(BUILD)); pg.wait_for_timeout(800)
        for n in ['diff_sg.js', 'diff_scat.js', 'diff_core.js', 'diff_db.js', 'diff_db_ui.js']: pg.add_script_tag(path=os.path.join(SRC, n))
        pg.add_style_tag(path=os.path.join(SRC, 'diff_db.css'))
        pg.evaluate('(()=>{ const d=document.createElement("div"); d.id="ddbHost"; document.getElementById("app").prepend(d); DIFF_DB_UI.render(d); })()'); pg.wait_for_timeout(500)
        ok('로그인하지 않음' in pg.inner_text('#ddbHost'), '앱 빌드 안에서 렌더 (전역 충돌 없음)')
        ok(pg.evaluate('typeof toast==="function"'), '앱 toast 사용')
        pg.screenshot(path=os.path.join(OUT, 'b01_in_app.png'))
        pg.close()
    b.close()

fb = [u for u in netreq if 'firebasejs' in u or 'firestore' in u or 'googleapis.com/identitytoolkit' in u]
ok(not fb, 'Firebase 네트워크 요청 0 개 (가짜 SDK)', fb[:3])
ext = [u for u in netreq if 'fonts.g' not in u]
ok(not ext, '외부 요청 없음 (앱의 Google Fonts 제외)', ext[:3])
ok(not errs, '콘솔 오류 0', errs[:6])
print(f'\n그림: {OUT}\n{"실패 " + str(len(fails)) + " 개" if fails else "모두 통과"}')
for f in fails: print('  ✗ ' + f)
sys.exit(1 if fails else 0)
