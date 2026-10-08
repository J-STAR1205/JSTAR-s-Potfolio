# M10. Quantum transport — 판단 노트

- Transport 4 (`fdsoi_negf.py`)와 Transport 5 (`poisson_negf_master.py`) 모두 Reference에서 끝까지 실행 완료.
  Transport 4는 band diagram(`band_0.png`)과 LDOS(`ldos_0.png`, `spectral_current_0.png`)를 생성했고,
  Transport 5는 Poisson-NEGF-Master equation으로 charge stability diagram(`PNEGF_Charge_diagram.png`)과
  전류 곡선(`PNEGF_current_featureless.png`, `PNEGF_Current_scale.png`)을 생성함.
- M06(Master equation만 사용)과의 차이: M10은 실제 전류(current)를 NEGF로 직접 계산하므로,
  Coulomb peak의 **높이·선폭**까지 원칙적으로 비교 가능 (M06은 점유 상태 기반 CSD만 제공).
- 다만 `PNEGF_current_scale.txt`가 전부 0으로 나오는 등 기본 튜토리얼 파라미터로는 전류가
  featureless(특징 없는 평탄한 곡선)에 가까움 — 실제 측정 Coulomb peak와 정량 비교하려면
  튜토리얼 기본 바이어스/온도/결합 세기를 실험값에 맞게 재조정하는 선행 작업이 필요해 보임.

**판단**: 전류 기반 정량 비교 틀 자체는 유효하나, 기본 설정 그대로는 바로 쓸 수 없음 —
측정 Coulomb peak와 비교할 실제 필요가 생기면(파라미터 보정 여력이 있을 때) 재검토.
