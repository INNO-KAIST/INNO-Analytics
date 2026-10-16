"""실측 AFM 파일만 정답 만들기 (합성은 만들지 않음) — gen_afm.py 의 Gwyddion 변환·저장 함수를 그대로 쓴다.
사용: python3 -I ref/gen_afm_real.py <실측 폴더>   (예 ref/real/afm, ref/real/ibw)
※ ref/real/ibw 는 사용자 연구실 측정 자료 — 저장소·배포에 올리지 않는다(inno_deploy 의 analytics_src 는 'real' 폴더 제외)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); REAL = sys.argv[1]
src = open(os.path.join(HERE, 'gen_afm.py'), encoding='utf-8').read()
head = src.split('# =============================================================== 합성 GWY')[0]
tail = src.split('# =============================================================== 실측')[1]
sys.argv = ['gen_afm.py', os.path.join(REAL, '_unused'), REAL]
exec(compile(head, 'gen_afm.py(head)', 'exec'))
exec(compile('# 실측' + tail, 'gen_afm.py(real)', 'exec'))
