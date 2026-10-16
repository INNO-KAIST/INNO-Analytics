"""비직관적 번역어 → 원어 치환표 (2026-09-30, 연구원 요청). 사용: python3 terms_map.py 파일...
MAP 먼저, MAP2 는 남은 up/down·branch·population 문맥 보정."""
import sys
MAP=[
("절단(censored)","censored"),("절단 점(상한)","censored 점(상한)"),("절단","censored"),
("첨두","peak"),
("올림·내림 가지가","up·down branch가"),("올림·내림","up·down"),("안정 가지","안정 branch"),("올림 가지","up branch"),("내림 가지","down branch"),
("위 가지","upper branch"),("아래 가지","lower branch"),("가지가 사라지면","branch가 사라지면"),("가지가 사라짐","branch가 사라짐"),("가지를 따라","branch를 따라"),
("내림 스윕","down sweep"),("올림 스윕","up sweep"),("(올림)","(up)"),("(내림)","(down)"),("· 내림","· down"),("· 올림","· up"),("올림=","up="),("내림=","down="),
('"올림"','"up"'),('"내림"','"down"'),
("프로파일 가능도","profile likelihood"),("가능도","likelihood"),
("비식별","non-identifiable"),("식별성","identifiability"),
("먼 국소해","먼 local minimum"),("국소해","local minimum"),("전역해","global minimum"),
("다중 출발","multi-start"),("χ² 단면","χ² slice"),("단면 계산","slice 계산"),
("고유방향","eigenvector"),("고유분해","eigen-decomposition"),("고유 분해","eigen-decomposition"),("고유값","eigenvalue"),
("유사역행렬","pseudo-inverse"),("준정부호","semi-definite"),("음함수 정리","implicit function theorem"),
("광자속","photon flux"),("유량","flux"),("몫","fraction"),
("임계 부근 느려짐","critical slowing down"),("95 % 임계","95 % threshold"),("임계","threshold"),("문턱","threshold"),
("계층 소거","profiled out"),("평균으로 소거","평균으로 profile out"),
("곱잡음","multiplicative noise"),
("다중포논","multiphonon"),("교차이완","cross-relaxation"),("상향변환","upconversion"),("광자 눈사태","photon avalanche"),("눈사태","avalanche"),
("쌍안정","bistable"),("준안정","metastable"),("몰분율","mole fraction"),("흡열","endothermic"),("배위장","crystal field"),
("다중항 바리센터","multiplet barycenter"),("다중항","multiplet"),("자리표시","placeholder"),("예열","warm-up"),("처리량","throughput"),
("병렬 확장성","parallel scaling"),("확장성","scaling"),("변동계수","CV"),("정규 분위수","normal quantile"),("잔차 분위수","residual quantile"),("분위수","quantile"),
("이득비","gain ratio"),("상자 경계","box constraints"),("multi-start의 상자","multi-start box"),
("암시적 Euler","implicit Euler"),("의사시간","pseudo-time"),("평균장","mean-field"),("척도 퇴화","scale degeneracy"),("퇴화","degeneracy"),("광퇴색","photobleaching"),
("감광 이온","sensitizer 이온"),("받개","acceptor"),("공여 준위","donor 준위"),("최소값 지형","χ² landscape"),("모델 계급","model class"),("명목 강도","nominal 강도"),("중적외","mid-IR"),
("켬 구간","on 구간"),("켬 ","on "),("· 끔","· off"),("읽음값","reading"),("흘러내림","cascade"),("수밀도","number density"),
("비공명","non-resonant"),("공명","resonant"),
]
MAP2=[("dir=+1 올림, −1 내림(올림","dir=+1 up, −1 down(up"),("올림/내림","up/down"),('" 내림"','" down"'),('" 올림"','" up"'),("`, 올림 ${nu}","`, up ${nu}"),
(">올림</option>",">up</option>"),(">내림</option>",">down</option>"),
("같은 가지가","같은 branch가"),("두 가지 사이","두 branch 사이"),("한 가지 위에","한 branch 위에"),("바닥 가지","lower branch"),
("실제로 점유될 수 있는","실제로 population이 생길 수 있는"),("점유되는 준위","population이 쌓이는 준위"),("점유 가능한","population 가능한"),("점유 가능성","population 가능성"),
("점유율","population"),("점유","population"),
("단면은 ‘","slice는 ‘"),("가장자리 resonant은","가장자리 resonance는"),("resonant이 자동","resonance가 자동"),
]
FIX=[("잔차가","residual이"),("잔차는","residual은"),("잔차와","residual과"),("잔차를","residual을"),("잔차","residual"),("차선과","alternative와"),("차선의","alternative의"),("차선은","alternative는"),("차선이","alternative가"),("화살표 차선","화살표 CHASEON_KEEP"),("차선","alternative"),("화살표 CHASEON_KEEP","화살표 차선"),("여기 강도","excitation intensity"),("피크 강도","peak intensity"),("peak 강도","peak intensity"),("강도","intensity"),("population가","population이"),("population는","population은"),("threshold을","threshold를"),("identifiability을","identifiability를"),("mean-field(mean-field)","mean-field")]
for p in sys.argv[1:]:
    s=open(p,encoding='utf-8').read(); n0=s
    for a,b in MAP+MAP2+FIX: s=s.replace(a,b)
    s=s.replace('/^placeholder/.test(','/^(placeholder|자리표시)/.test(')
    if s!=n0: open(p,'w',encoding='utf-8').write(s); print('changed',p)
