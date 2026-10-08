# M11. Maxwell / CPW — 판단 노트

- `cap_xmon.py`(커패시턴스 행렬), `maxwell_xmon.py`(Maxwell eigenmode, inductor port/EPR),
  `maxwell_meandered_resonator.py`(적응형 메시 공진기) 모두 Reference에서 끝까지 수렴.
- 커패시턴스 행렬: adaptive mesh 8회 반복(445,991 노드) 후 수렴, xmon_cross 자기-커패시턴스
  ≈ 9.85e-14 F (≈ 98.5 fF) — 트랜스몬 패드의 전형적인 범위.
- Maxwell eigenmode: xmon 공진 주파수(mode 0) ≈ 5.70 GHz (21회 반복, 96,486 노드),
  meandered resonator mode 0/1 ≈ 4.75 GHz / 14.4 GHz — 모두 실제 CPW 큐비트/공진기에서
  보는 수 GHz대 값과 정성적으로 일치함.
- 즉 QTCAD의 정전용량·Maxwell eigenmode 결과는 물리적으로 타당한 범위에 들어오지만,
  기존(비-QTCAD) CPW 설계 워크플로의 같은 구조에 대한 결과와 **직접 비교는 아직 하지 않음**.

**판단**: 결과 자체는 신뢰할 만한 범위이나, "기존 워크플로 대체/보조 가치"를 확정하려면
같은 레이아웃 하나를 기존 도구로도 돌려 두 결과를 나란히 비교하는 단계가 필요함 —
그 비교 전까지는 보조/교차검증 도구로만 사용 권장.
