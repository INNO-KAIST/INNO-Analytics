"""전수 검수(2026-10-16) 회귀 — 계산 중 되돌리기·닫기 막기, 중단, .ijm 을 ImageJ 가 끝까지 실행, Measure ROI 재생, batch 옆 파일, Convert to Mask 스택 유지
사용: python3 test/ui_audit.py <주입된 html> <그림 폴더> [ij.jar]"""
import sys, os, json, subprocess, tempfile, shutil, numpy as np
from playwright.sync_api import sync_playwright
html, out = sys.argv[1], sys.argv[2]; JAR = sys.argv[3] if len(sys.argv) > 3 else '/tmp/ij.jar'; os.makedirs(out, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__)); IOD = os.path.join(H, 'io'); REFIN = os.path.join(H, '..', 'ref', 'data', 'in')
res = []; errs = []
def chk(c, name, info=''): res.append((bool(c), name, str(info)))
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1500, 'height': 950})
    pg.on('pageerror', lambda e: errs.append(str(e))); pg.on('console', lambda m: errs.append(m.text) if m.type == 'error' and 'ERR_' not in m.text else None)
    pg.goto('file://' + os.path.abspath(html)); pg.wait_for_timeout(700); pg.click('#nav button[data-v="imgx"]'); pg.wait_for_selector('.imx')
    # 1) 계산 중 되돌리기·닫기 막기
    pg.evaluate("(()=>{ const W=1600,H=1600,d=new Float32Array(W*H); for(let i=0;i<d.length;i++) d[i]=(i*2654435761>>>24)&255; IMG_UI._addImage({name:'big',w:W,h:H,bit:8,planes:[d]}); IMG_UI.rerender(); })()")
    pg.evaluate("IMG_UI._run('math',{fn:'add',value:1})"); pg.wait_for_timeout(500)
    pg.evaluate("window._p=IMG_UI._run('type',{to:'32',scale:false})"); pg.wait_for_function('!IMG_UI._busy()', timeout=60000); pg.evaluate("(()=>{ window._p=IMG_UI._run('median',{radius:10}); })()"); pg.wait_for_timeout(150) # 약속을 기다리지 않게
    busy = pg.evaluate('IMG_UI._busy()'); pg.evaluate('IMG_UI._undo()'); pg.evaluate('IMG_UI._close(IMG_UI._cur().id)')
    pg.wait_for_function('!IMG_UI._busy()', timeout=120000); pg.wait_for_timeout(300)
    st = pg.evaluate("(()=>{ const im=IMG_UI._cur(); return {n:IMG_UI._state().images.length, ops:im?im.ops.map(o=>o.op):null}; })()")
    chk(busy and st['n'] == 1 and st['ops'] == ['math', 'type', 'median'], '계산 중 되돌리기·닫기는 막힘 (기록·화소 일치)', st)
    # 2) 중단: 기다리던 장까지 멈추고 영상은 그대로
    pg.evaluate("(()=>{ const W=700,H=700,P=[]; for(let k=0;k<24;k++){ const d=new Float32Array(W*H); for(let i=0;i<d.length;i++) d[i]=((i*(k+1)*40503>>>20)&255)+0.5; P.push(d); } IMG_UI._addImage({name:'stk24',w:W,h:H,bit:32,nz:24,planes:P}); IMG_UI.rerender(); })()")
    before = pg.evaluate("Array.from(IMG_UI._px.get(IMG_UI._cur().id)[20].slice(0,50))")
    pg.evaluate("(()=>{ window._p2=IMG_UI._run('median',{radius:8,stack:true}); })()"); pg.wait_for_timeout(400); pg.evaluate('IMG_UI._cancel()')
    pg.wait_for_function('!IMG_UI._busy()', timeout=60000); pg.wait_for_timeout(1500)
    st2 = pg.evaluate("(()=>{ const im=IMG_UI._cur(); return {ops:im.ops.length, prog:document.querySelector('#imProg').hidden, px:Array.from(IMG_UI._px.get(im.id)[20].slice(0,50))}; })()")
    chk(st2['ops'] == 0 and st2['prog'] and st2['px'] == before, '중단 → 작업 안 됨·진행 막대 사라짐 (기다리던 장 포함)', {k: st2[k] for k in ('ops', 'prog')})
    # 3) Convert to Mask: 스택에서 '모든 장' 끔 → 스택 유지 (지금 장의 문턱을 모든 장에)
    pg.evaluate("IMG_UI._run('threshold',{method:'Default',dark:true,stack:false,newImage:false})"); pg.wait_for_function('!IMG_UI._busy()', timeout=60000); pg.wait_for_timeout(300)
    nz = pg.evaluate("(()=>{ const im=IMG_UI._cur(); return [im.nz, IMG_UI._px.get(im.id).length, im.bit]; })()"); chk(nz == [24, 24, 8], 'Convert to Mask (스택, 모든 장 끔) → 24 장 유지', nz)
    # 4) Measure 를 ROI 와 함께 기록 → 다른 영상에서 같은 ROI 로 재생
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.macro.steps=[]; S.macro.rec=true; const d=new Float32Array(100*80).fill(5); IMG_UI._addImage({name:'m1',w:100,h:80,bit:8,planes:[d]}); IMG_UI.rerender(); IMG_UI._setActive({id:'r1',type:'rect',pts:[[10,10],[30,20]]}); })()")
    pg.evaluate("IMG_UI._run('measure',{all:false})"); pg.wait_for_timeout(200)
    pg.evaluate("(async()=>{ const d=new Float32Array(100*80).fill(9); const im=IMG_UI._addImage({name:'m2',w:100,h:80,bit:8,planes:[d]}); await IMG_UI._runSteps(IMG_UI._steps().slice(-1),im); })()"); pg.wait_for_timeout(400)
    last = pg.evaluate("IMG_UI._state().results.slice(-1)[0]"); chk(last and last['Area'] == 200 and last['Mean'] == 9, 'Measure 재생 = 기록한 ROI (20×10)', last and {'Area': last['Area'], 'Mean': last['Mean']})
    # 5) batch 파일: 옆 .txt 보정 쓰고, .txt 자체는 영상으로 열지 않음
    pg.evaluate("(()=>{ const S=IMG_UI._state(); S.macro.steps=[{op:'measure',params:{all:false},label:'Measure'}]; S.ui.panel='macro'; IMG_UI.rerender(); })()"); pg.wait_for_timeout(300)
    n0 = pg.evaluate('IMG_UI._state().results.length'); pg.set_input_files('#mcFiles', [os.path.join(IOD, 'thitachi.tif'), os.path.join(IOD, 'thitachi.txt')])
    pg.wait_for_function(f'IMG_UI._state().results.length>{n0} && !IMG_UI._busy()', timeout=30000); pg.wait_for_timeout(300)
    rr = pg.evaluate(f'IMG_UI._state().results.slice({n0})'); chk(len(rr) == 1 and rr[0].get('Unit') == 'nm', 'batch 파일: Hitachi 옆 .txt 보정 사용 (nm), .txt 는 영상 아님', [(r.get('Image'), r.get('Unit')) for r in rr])
    # 6) .ijm 을 ImageJ 1.54u 가 대화 상자 없이 끝까지 실행 (batch mode)
    if os.path.exists(JAR) and os.path.exists(os.path.join(REFIN, 'g8.tif')):
        pg.set_input_files('#imFin', [os.path.join(REFIN, 'g8.tif'), os.path.join(REFIN, 'stk.tif')]); pg.wait_for_function("IMG_UI._state().images.some(i=>i.name==='stk.tif')", timeout=20000)
        steps = [("gaussian", {'sigma': 1.5, 'stack': True}), ("flip", {'dir': 'h'}), ("rotate", {'angle': 5, 'interp': 'bilinear', 'enlarge': False}), ("resize", {'factor': 0.8, 'interp': 'bilinear'}),
                 ("enhance", {'saturated': 0.35, 'normalize': True, 'stack': True}), ("bandpass", {'large': 40, 'small': 3, 'stripes': 'none', 'tolerance': 5, 'autoscale': True, 'saturate': True, 'stack': True}),
                 ("duplicate", {'slice': False}), ("type", {'to': '16', 'scale': False}), ("type", {'to': '8', 'scale': True}),
                 ("threshold", {'method': 'Otsu', 'dark': True, 'stack': False, 'newImage': False}), ("binary", {'fn': 'fill', 'stack': True}), ("binary", {'fn': 'edm', 'edm32': True})]
        pg.evaluate("(()=>{ const S=IMG_UI._state(); S.sel=S.images.find(i=>i.name==='stk.tif').id; S.macro.steps=[]; S.macro.rec=true; IMG_UI.rerender(); })()")
        for op, prm in steps:
            pg.evaluate("([o,p])=>IMG_UI._run(o,p)", [op, prm]); pg.wait_for_function('!IMG_UI._busy()', timeout=60000); pg.wait_for_timeout(150)
        pg.evaluate("IMG_UI._setActive({id:'rr',type:'rect',pts:[[5,5],[40,30]]})"); pg.evaluate("IMG_UI._run('measure',{all:false})"); pg.wait_for_timeout(200)
        pg.evaluate("IMG_UI._run('analyze',{minSize:0,maxSize:'',minCirc:0,maxCirc:1,exclude:true,holes:false,show:'outlines',unitPx:true,redirect:''})"); pg.wait_for_function('!IMG_UI._busy()', timeout=60000)
        ijm = pg.evaluate('IMG_UI._toIJM()'); tmp = tempfile.mkdtemp(); shutil.copy(os.path.join(REFIN, 'stk.tif'), tmp)
        ijm = ijm.replace('dir = getDirectory("영상 폴더");', f'dir = "{tmp}/";') + '\nprint("IJM_DONE");\n'
        open(os.path.join(out, 'exported.ijm'), 'w').write(ijm); mf = os.path.join(out, 'run_exported.ijm'); open(mf, 'w').write(ijm)
        try:
            o = subprocess.run(f'JAVA_TOOL_OPTIONS= xvfb-run -a java -Xmx1g -cp {JAR} ij.ImageJ -batch {mf}', shell=True, capture_output=True, text=True, timeout=180)
            txt = o.stdout + o.stderr; done = 'IJM_DONE' in txt; bad = [l for l in txt.splitlines() if 'Macro Error' in l or 'Unrecognized' in l or 'not found' in l.lower()]
            chk(done and not bad, '.ijm (기록 13 단계) → ImageJ 1.54u 가 대화 상자 없이 끝까지 실행', (bad or txt.strip().splitlines()[-3:]))
        except subprocess.TimeoutExpired:
            chk(False, '.ijm → ImageJ 실행', '180 s 시간 초과 (대화 상자에서 멈춤)')
        shutil.rmtree(tmp, ignore_errors=True)
    b.close()
chk(not errs, '콘솔 오류 0', errs[:5])
print('\n'.join(('ok  ' if c else 'FAIL') + ' ' + n + ('' if c else '  (' + i + ')') for c, n, i in res)); print('FAILS', sum(1 for c, _, _ in res if not c))
