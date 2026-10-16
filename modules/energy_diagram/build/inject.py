"""더 쓰지 않음 — Energy Diagram·DRE 를 함께 주입하는 ../../inject_all.py 를 쓸 것.
(두 모듈이 같은 본 앱 코드 줄(탭 목록 등)을 고치므로 주입기는 하나여야 멱등이 보장된다.)"""
import sys, os, subprocess
sys.exit(subprocess.call([sys.executable, os.path.join(os.path.dirname(__file__),'..','..','inject_all.py')]+sys.argv[1:]))
