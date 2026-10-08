# QTCAD M03 / M06 / M08 물리 모델 및 결과 신뢰성 검토 보고서

**검토 범위**: 코드 수정 없음 (추적·검증만 수행). 기존 tutorial 재현과 분석 notebook에서 추가된 해석을 구분하여 표시.
**검토 방법**: 설치된 QTCAD 2.2.6 실제 소스/docstring 확인, 저장된 output 파일 직접 로드, 일부 항목은 독립 재계산/재실행으로 교차검증.

---

## 1. Executive Summary

### M03 — Source/Drain 접합부 Potential Wiggle
Source/channel(y=-45nm), channel/drain(y=+45nm) 경계 바로 바깥쪽에서 작은 potential wiggle(좌측 ~3.6mV, 우측 ~2.7mV)이 관찰된다. **원본 스크립트가 쓰는 `linecut(method="default")` 결과는 비교적 완만한 wiggle**이지만, **`method="pyvista"`로 바꾸면 서브-옹스트롬 간격에서 최대 26mV까지 튀는 물리적으로 불가능한 jump**가 나타난다 — 이는 pyvista 보간기가 이 비정형 사면체 메시에서 셀 탐색 아티팩트를 일으킨다는 강력한 증거다. 좌우 전기장(Ey)의 거동이 비대칭(한쪽은 부호반전, 한쪽은 부호 유지)인 점은 메시가 완벽히 대칭이 아니라는 신호다. **결론: Mixed/Unresolved에 가까움** — 작은 wiggle 자체는 "default" 방식에서도 완전히 사라지지 않아 순수 plotting artifact로 단정할 수 없지만, 메시 수렴 테스트(coarse vs refined 비교)가 수행되지 않아 numerical vs physical 여부를 확정할 수 없다.

### M06 — Charge Stability Diagram이 지나치게 직사각형인 이유
**가장 중요한 발견**: 의심했던 두 항목(bias reference 처리, 10K 온도) 중 하나는 **버그가 아니고(API 문서로 확인됨)**, 다른 하나는 **원본 저자가 의도적으로 명시한 설정**이다. 실제 원인은 다음 세 가지의 조합이다: (1) cross lever arm이 dominant 대비 1/100 수준(0.006~0.007 vs 0.765)이라 transition line이 거의 수직/수평이 되고, (2) lever-arm 행렬·에너지·Coulomb 행렬이 전체 90×90 스윕에서 **단 한 번만 계산되어 고정 재사용**되며, (3) `overlap=False`로 인해 Coulomb 행렬이 대각 성분만 남아 궤도 혼성화(anticrossing, triple-point curvature)를 원천적으로 표현할 수 없다. Noise는 전혀 없음(코드에 random 생성 자체가 없음, grep으로 확인).

### M08 — 82µs T2* 수치와 Rabi 주파수 불일치
**가장 중요한 발견**: 노트북이 "RWA 컨벤션 차이"로 설명한 1.6899MHz vs 1.5177MHz 불일치는 **실험적으로 반증됨** — `EDSR_noise.py`의 하드코딩된 `h0`,`u` 행렬을 동일한 `transition_2_levels()` 함수에 다시 넣으면 1.5177MHz가 재현되어, 두 스크립트 자체는 자기일관적이다. 즉 **`EDSR_noise.py`의 하드코딩된 행렬이 현재 `MOS_EDSR.py`의 최신 실행 결과와 맞지 않는 stale-data 문제**다. 82µs라는 T2* 값은 시뮬레이션 시간(13.2µs)의 6.2배에 달하는 외삽이며, 그 구간에서 envelope는 겨우 2.5%만 감소한 상태 — 사실상 "관찰되지 않은 decay"를 fit한 것이다. 게다가 Ramsey/Hahn-echo 시퀀스가 전혀 구현되어 있지 않으므로 이 수치를 "T2*"라 부르는 것 자체가 용어상 부적절하다(최대 "Rabi decay time").

---

## 2. Confirmed Findings (실제 코드/데이터로 확인된 사실)

| # | 모듈 | 확인 사실 |
|---|------|-----------|
| 1 | M03 | `dqdfdsoi.geo` 기준 source/channel 경계 y=-45nm, channel/drain 경계 y=+45nm. 게이트: B1=[-40,-30], P1=[-25,-10], B2=[-5,5], P2=[10,25], B3=[30,40]. 접합과 가장 가까운 게이트 사이 5nm의 **무게이트 갭**(gate_oxide만 존재) 존재. |
| 2 | M03 | Wiggle 정량화(default linecut): 좌측 min 632.20mV@y=-48.739nm, max 635.76mV@y=-47.850nm (진폭 3.56mV); 우측 max 634.73mV@y=48.240nm, min 632.01mV@y=49.554nm (진폭 2.72mV). |
| 3 | M03 | 접합부 근처 로컬 메시 요소 크기: min 0.526nm / mean 1.071nm / max 2.970nm. Dot 중앙 영역 mean 0.697nm. 접합부가 dot 영역보다 약 1.5배만 조대(粗大) — 극단적 차이 아님. |
| 4 | M03 | `linecut(method="pyvista", res=2000)`는 default와 **질적으로 다른 결과**(날카로운 불연속 jump, 최대 26mV/0.065nm)를 낸다. `an.linecut` docstring: `method="default"`는 QTCAD 자체 FEAcpp 선형 보간(min_interval=1e-13 병합), `res`는 pyvista 전용 — 두 방식이 완전히 다른 보간 엔진을 쓴다. |
| 5 | M03 | 좌우 게이트 바이어스는 거의 대칭(B1=B3=0.5V, P1=P2=0.59V±0.1mV detuning, B2=0.51V는 대칭축상 섭동)이지만 Ey(y) 거동은 비대칭 — 좌측은 부호반전(−2.89, −3.88 → +22 MV/m), 우측은 부호 유지하며 비단조 wobble. |
| 6 | M03 | `Ec(y) = -eφ(y) - χ_Si`이고 source/drain/channel 모두 동일 재질(Si, 도핑만 다름) — 즉 Ec는 -φ(y)의 **대수적으로 동일한** 상수이동이며 독립적 물리량이 아님 (노트북 자체도 이미 이 점을 명시). |
| 7 | M06 | 전체 계산 체인: Mesh→`get_double_dot_fdsoi()`(T=0.1K)→비선형 Poisson(adaptive, tol=1e-3)→Schrödinger(SubDevice, num_states=4)→`LeverArmSolver.solve(bias_increment=1e-3)` **1회**→`get_coulomb_matrix(overlap=False)` **1회**→`Junction` 생성(에너지·Coulomb 고정)→90×90 루프(`set_biases`+`add_spectrum`만 반복)→plot. |
| 8 | M06 | 저장된 `lever_arm_matrix.npy` 실측값: `[[0.00622386,0.76485381],[0.76499903,0.00669277],[0.00623200,0.76449898],[0.76457140,0.00658892]]` — 사용자 claim과 정확히 일치. |
| 9 | M06 | `LeverArmSolver`는 `plunger_gate_1_bias=0.6V, plunger_gate_2_bias=0.601V`에서 스윕 **이전에** 1회 계산되고, 루프 바디(90×90)는 이 행렬을 다시 계산하지 않음 — 코드로 직접 확인. |
| 10 | M06 | `many_body_solver_params.energies`, `.coulomb_mat`도 루프 **이전**에 고정 — 루프 안에서는 게이트 바이어스만 바뀜. |
| 11 | M06 | `get_coulomb_matrix(overlap=False)` 확인. 설치된 `qtcad/device/many_body.py` docstring: "overlap=True → 4D 배열 V_ijkl, False → 2D 배열 V_ij=V_ijij" — **문서로 확인된 사실**이며, "혼성화를 못 나타낸다"는 결론은 이 행렬 구조에서 **직접 따르는 추론**(docstring이 명시적으로 "no hybridization"이라 쓴 것은 아님). |
| 12 | M06 | `dvc.set_temperature(0.1)` (helper, 전기정역학용) vs `temperature_spec=10` (addition spectrum용) — **원본 스크립트에 "10K로 설정해 CSD 선을 더 두껍고 관찰하기 쉽게 만든다"는 명시적 주석이 바로 위에 존재** → 의도된 튜토리얼 설정, 버그 아님. |
| 13 | M06 | `Junction` 클래스의 `biases`/`set_bias_vector`/`set_biases`/`set_bias` **네 곳 모두 동일한 문구**로 "바이어스는 reference contact potential(초기 electrostatics를 푼 경계전위) 기준 상대값"이라고 명시. 전기정역학이 P1=0.6V, P2=0.601V로 풀렸으므로, 스윕의 30~140mV는 그 기준의 delta이고, plot에서 `+0.6V`/`+0.601V`를 더해 절대축을 복원하는 것은 **API 문서와 정합** — 버그 아님. |
| 14 | M06 | `interpolation="bilinear"`가 `double_dot_stability.py` 1곳 + notebook 2곳에서 확인됨. |
| 15 | M06 | `double_dot_stability.py`와 notebook 전체에서 `random`, `noise`, `seed`, `normal(` 등을 grep — **일치 항목 없음** (bilinear 3건 제외). CSD는 완전히 deterministic. |
| 16 | M08 | `Bfield(x,y,z)`: `B0=0.6, b=0.3e6, return [B0,0,b*x]` — 균일장 + 완전 선형 gradient, 다른 항 없음. 실측/COMSOL/magnetostatic 필드맵 사용 없음. |
| 17 | M08 | `mos_edsr_log.txt`: `f_rabi = 1689879.77` Hz. `EDSR_noise.py`의 하드코딩 `h0,u`로 **동일한 `transition_2_levels()` 함수**를 재호출하면 `f_rabi=1517665.99` Hz(1.5177MHz, 로그의 1.6899MHz가 아님) — **RWA 설명이 실험적으로 반증**됨. `h0` 대각성분은 현재 로그의 에너지와 거의 일치하지만 drive 행렬 `u`(또는 그것을 만든 게이트 바이어스 스윕)는 현재 실행과 불일치 — stale/mismatched hardcoded 값. |
| 18 | M08 | `qtcad.qubit.spectra`에 `lorentz`와 `power_law`(1/f^α) **두 모델 모두 존재** — "Lorentzian만 지원"이라는 전제는 틀림. `Noise.dynamics()` 시그니처는 단일 `H0`/`omega`+하나의 확률과정만 받음, 그리고 스크립트 주석("charge noise는 drive delta_V를 만드는 게이트에만 영향")과 결합하면 현재 구현은 **drive-amplitude noise만** 구조적으로 가능 — detuning/주파수 noise 경로는 이 호출에 없음. |
| 19 | M08 | `T_Rabi=658.91ns`, 시뮬레이션 전체시간=13178.2ns(13.18µs). Fit한 `T2star=81991.22ns`(82µs). 비율 fit/window≈6.2배. `exp(-(13178/81991)^2)=0.9745` — **창 끝에서 envelope가 겨우 2.5% 감소**. |
| 20 | M08 | M08 전체 스크립트/노트북에 Ramsey(π/2-자유진화-π/2) 또는 Hahn-echo 시퀀스 **없음** — 유일한 동역학은 연속 구동 Rabi oscillation. |

---

## 3. Suspected Issues (추가 실행/검증 필요)

| # | 모듈 | 의심 사항 | 추가로 필요한 것 |
|---|------|-----------|------------------|
| 1 | M03 | Wiggle이 numerical artifact인지 physical feature인지 — mesh convergence 테스트(coarse vs refined 비교)가 수행되지 않음. 현재는 refined mesh의 solution만 존재하고 coarse mesh의 독립 solution이 없음. | source/drain 접합 ±5nm만 국소 재정제한 새 Poisson 재계산 (원본 파일 수정 없이 별도 스크립트로). |
| 2 | M03 | `cond_band_edge()`가 저장되지 않아 phi/−phi/Ec 3중 비교가 Ec 계산 자체로는 안 됨(다만 Ec≡−φ+상수이므로 이 특정 소자에서는 비교가 애초에 무의미함을 확인). | 재질이 다른 소자에서 재검증하면 더 의미 있음 (이 소자는 전부 Si). |
| 3 | M06 | `overlap=True`로 바꾸면 실제로 anticrossing/triple-point curvature가 나타나는지 — 아직 실행 안 함. | `overlap=True` 재실행 비교 (계산 비용 확인 필요). |
| 4 | M06 | Cross lever arm이 게이트/배리어 geometry를 바꾸면 실제로 커지는지 — 아직 실행 안 함. | M06-C(geometry 변경) 실험. |
| 5 | M08 | `EDSR_noise.py`의 `h0,u`가 정확히 "몇 번째 과거 실행"에서 나온 값인지 — git 이력이나 별도 로그가 없어 확정 불가. | `MOS_EDSR.py`를 현재 파라미터로 재실행해 `h0,u`를 다시 뽑아 `EDSR_noise.py`에 갱신하면 해결되는지 확인. |
| 6 | M08 | Noise 유무 결과 차이가 왜 작은지(100-run 평균에서 noise effect가 제대로 나타나는지) — 정량 비교 미수행. | S0 sweep (0, 0.001...0.1) 실험. |

---

## 4. M03 Review — Potential Wiggle

### 4.1 Geometry and Doping Mapping

| y position (nm) | 영역 | 비고 |
|---|---|---|
| -65 ~ -45 | source (n+, 1e20 cm⁻³) | |
| -45 | **source/channel 접합** | wiggle 좌측 중심 (y≈-47.9 ~ -48.7) |
| -45 ~ -40 | 무게이트 갭 (gate_oxide only) | |
| -40 ~ -30 | barrier_gate_1 (B1) | |
| -25 ~ -10 | plunger_gate_1 (P1) | |
| -5 ~ 5 | barrier_gate_2 (B2) | |
| 10 ~ 25 | plunger_gate_2 (P2) | |
| 30 ~ 40 | barrier_gate_3 (B3) | |
| 40 ~ 45 | 무게이트 갭 | |
| 45 | **channel/drain 접합** | wiggle 우측 중심 (y≈48.2 ~ 49.6) |
| 45 ~ 65 | drain (n+, 1e20 cm⁻³) | |

### 4.2 Mesh Audit

| 영역 | min edge (nm) | mean edge (nm) | max edge (nm) |
|---|---|---|---|
| 접합부 ±5nm (양쪽) | 0.526 | 1.071 | 2.970 |
| Dot/채널 중앙 (±10nm) | — | 0.697 | — |

접합부가 dot 영역보다 **약 1.5배만** 조대 — `params_poisson.refined_region`에 접합부가 명시적으로 포함되지 않았음에도 전역 적응정제가 어느 정도 접합부까지 미쳤다는 뜻. 극단적 차이는 아니므로 "완전히 방치된 영역"은 아니다.

### 4.3 Wiggle Quantification

| 항목 | 좌측 (source측) | 우측 (drain측) |
|---|---|---|
| Local min (mV) | 632.20 | 632.01 |
| Local max (mV) | 635.76 | 634.73 |
| 진폭 (mV) | 3.56 | 2.72 |
| 위치 (접합에서 거리, nm) | 2.85 ~ 3.74 | 3.24 ~ 4.55 |

좌우 진폭이 ~30% 차이 — 완벽한 대칭은 아님 (대칭 바이어스에도 불구).

### 4.4 Linecut Interpolation Method Comparison (결정적 테스트)

| Method | 좌측 거동 | 우측 거동 |
|---|---|---|
| `default` (원본 스크립트 사용) | 완만한 overshoot/dip, 진폭 수 mV | 동일 |
| `pyvista, res=2000` | **불연속 jump 26mV / 0.065nm** (y=-48.100nm) | jump 5.7mV, 0.9mV |

`pyvista` 방식의 jump 크기(수십 mV를 서브-옹스트롬 거리에서)는 어떤 물리적 전기장으로도 설명 불가능 — 비정형 사면체 메시에서 VTK 셀 탐색 아티팩트로 강하게 의심됨. **이것이 "plotting artifact 가능성"에 대한 가장 강력한 증거**이지만, `default` 방식에서도 작은 wiggle(2.7~3.6mV)이 남아있어 전체 wiggle이 순수 plotting 문제라고 단정할 수 없다.

### 4.5 Physical Consistency (Ey = -dφ/dy)

- **좌측**: 멀리서 약한 양의 Ey(~0.1~1 MV/m) → **부호반전**(-2.89, -3.88 MV/m) → 다시 양으로 전환 → 접합에서 +22 MV/m까지 단조증가.
- **우측**: 같은 구간에서 Ey는 **부호반전 없이** 음의 값 유지, 크기만 비단조 진동(-0.78→-0.32→-1.21→-0.60→-0.31 MV/m).

좌측의 "부호반전"은 실제 depletion-region 물리(n+ → intrinsic 접합에서 band bending 방향이 바뀌는 지점이 있을 수 있음)로 설명 가능할 수도 있으나, 우측에 대칭적으로 나타나지 않는 점은 의심스럽다.

### 4.6 Final Judgment

**Case 3 (Mixed) 또는 Case 4 (Unresolved)에 가장 가깝다.**

근거:
1. `pyvista` 방식의 극단적 jump는 명백한 numerical/plotting artifact (근거: 물리적으로 불가능한 크기).
2. 그러나 원본이 실제로 쓰는 `default` 방식에서도 작은(2.7~3.6mV) wiggle이 남아있고, 이것이 사라지는지 확인할 **mesh convergence 테스트(coarse vs refined 재계산)가 수행되지 않았다** — 이 테스트 없이는 "default 방식의 작은 wiggle"조차 numerical인지 physical인지 결정할 수 없다.
3. 좌우 비대칭(진폭 30% 차이, Ey 부호거동 차이)은 완전한 물리적 대칭이 기대되는 상황에서 발생 — 메시 비대칭 또는 미세한 수치 오차의 신호일 수 있다.

### 4.7 Recommended Modification

**즉시 코드 수정은 권장하지 않음.** 대신:
1. 원본 `tunnel_coupling_1.py`/`tunnel_coupling_2.py`/notebook은 그대로 보존.
2. 별도 `M03_mesh_convergence.py` (접합부 ±5nm만 국소 재정제, 동일 bias로 재계산)를 작성해 wiggle amplitude가 메시 세분화에 따라 감소하는지 확인.
3. 확인 결과에 따라 `M03_Potential_Wiggle_Review.ipynb`에 "numerical artifact로 확인" 또는 "물리적 특징, 단 mesh-dependent 정밀도 주의" 중 하나로 결론 업데이트.

---

## 5. M06 Review — Charge Stability

| 항목 | 현재 구현 | 물리적 의미 | 이상화 수준 | 결과에 미치는 영향 | 수정 필요 | 우선순위 |
|---|---|---|---|---|---|---|
| Fixed lever-arm | 1회 계산(0.6V/0.601V 기준), 90×90 전체 재사용 | 국소 선형화(linearization) | 중간 — 110mV 스윕 범위에 비해 단일 선형화점 | transition line을 정확한 직선으로 만듦 (CSD가 "너무 반듯한" 1차 원인) | 가능하면 Yes (여러 bias point에서 재계산) | **High** |
| Fixed energies/Coulomb | 루프 이전 1회 계산, 고정 | wavefunction/dot 위치 불변 가정 | 중간~높음 | lever-arm 고정과 결합해 CSD를 완전히 선형적인 모델로 만듦 | 선택적 (대규모 재설계 필요) | Medium |
| `overlap=False` | Coulomb 행렬 대각성분만 사용 | 궤도 혼성화(tunnel coupling에 의한 anticrossing) 배제 | 높음 — triple-point curvature를 원천적으로 못 그림 | honeycomb 모서리가 날카로운 직선 교차로 보임 (CSD가 "직사각형"으로 보이는 **두 번째** 핵심 원인) | Yes, `overlap=True` 비교 실험 권장 | **High** |
| 10K 온도 (addition spectrum) | 전기정역학 0.1K, transport 10K | 의도된 튜토리얼 설정(주석 확인) | 낮음 (문서화된 선택) | 선 폭을 두껍게 — "반듯함"의 원인이 아니라 "보기 쉬움"의 원인 | No (의도된 설정, 필요시 sensitivity study로 비교만) | Low |
| Bias reference 처리 | Junction API 문서와 정합 확인됨 | - | 없음 (버그 아님) | 없음 | No | — |
| `interpolation="bilinear"` | 3곳에서 사용 | 시각화만 | 낮음 | 선을 부드럽게 보이게 하지만 물리 결과 자체는 아님 | 비교용으로 `nearest` 테스트 권장 | Low |
| Noise 부재 | 전혀 없음 (grep 확인) | 실험 CSD의 broadening/drift/jump 없음 | 높음 (가장 단순한 모델) | **"반듯함"의 세 번째 원인** — 단, cross-lever-arm/overlap 문제보다 후순위로 다뤄야 함 (문제 지시대로) | Yes, 별도 study로 | Medium |

### 핵심 질문에 대한 답

1. **CSD가 직사각형에 가까운 가장 큰 이유**: cross lever arm이 극단적으로 작음(0.006~0.007) **+** `overlap=False`로 혼성화 자체가 모델에 없음 — 두 요인의 조합. Noise 부재는 "날카로움"의 원인이지 "직사각형 모양" 자체의 원인은 아님.
2. **Noise 부재의 영향**: 선이 완벽하게 가늘고 깨끗하게 보이는 원인. 모양(직사각형 여부)에는 영향 없음.
3. **Cross lever arm 크기**: dominant 대비 ~1/100 — 매우 작음. 실제 FD-SOI 소자에서 이 정도로 작은 cross-capacitance가 전형적인지는 추가 문헌 비교 필요.
4. **Fixed-alpha approximation의 이상화 정도**: 110mV 스윕 전체에 단일 선형화 — transition line을 완전한 직선으로 만드는 직접적 원인.
5. **Fixed energies/Coulomb의 영향**: lever-arm 고정과 결합해 전체 모델을 "순수 선형" many-body 모델로 만듦 — 실제 게이트 바이어스에 따른 dot 재형성 효과를 누락.
6. **overlap=False의 영향**: anticrossing/triple-point curvature를 원천적으로 표현 불가 — High priority 이슈.
7. **10K 설정의 목적/부작용**: 목적은 명시된 대로 "관찰 용이성". 부작용은 실제 저온 큐비트 조건(0.1K)과의 정량적 차이를 감춤.
8. **Plot interpolation의 착시**: 가능성 있으나 bilinear는 일반적으로 데이터 포인트 사이만 부드럽게 하므로, 90×90 해상도에서 transition line 각도 자체를 바꾸지는 않을 것 — 영향은 작을 것으로 추정 (비교 실험 권장).
9. **Voltage-axis/reference 처리 오류**: **없음** — API 문서로 확인된 정상 동작.
10. **실제 Si DQD CSD에 가까워지기 위한 최우선 수정**: `overlap=True`로 바꿔 혼성화를 켜는 것(High), 그 다음 cross lever arm의 geometry 기반 민감도 분석(M06-C).

---

## 6. M08 Review — Micromagnet EDSR

| 항목 | 현재 구현 | 물리적 의미 | 이상화 수준 | 결과에 미치는 영향 | 수정 필요 | 우선순위 |
|---|---|---|---|---|---|---|
| Magnetic field 모델 | `B=[B0, 0, b*x]` (균일+완전선형) | 이상화된 micromagnet gradient | 매우 높음 — fringing/curvature/fabrication variation 전부 없음 | EDSR 구동 강도가 과도하게 균일/예측 가능하게 나옴 | 이번 단계 불필요 (명칭과 구현 수준 구분만 명시) | Low (문서화만) |
| `h0,u` 하드코딩 불일치 | `MOS_EDSR.py` 로그(1.69MHz)와 `EDSR_noise.py` 하드코딩(1.52MHz) 불일치, RWA 때문 아님(반증됨) | - | 버그 (stale data) | EDSR_noise.py 전체 결과가 현재 MOS_EDSR.py 파라미터와 대응하지 않음 | **Yes** | **Critical** |
| Noise → drive-amplitude만 | `Noise.dynamics()` 시그니처가 단일 H0/omega만 받음 | 실제 charge noise 중 detuning 경로 누락 | 높음 | 실제 스핀 큐비트의 주요 decoherence 경로(주파수 노이즈) 미반영 | 가능하면 QuTiP로 custom Hamiltonian 구현 | Medium |
| 82µs "T2*" fit | 시뮬레이션 창(13.2µs)의 6.2배 외삽, envelope 2.5%만 감소 | 사실상 "관찰되지 않은 decay"의 추정 | 높음 | 보고된 수치가 데이터로 constrain되지 않음 | **Yes** — lower bound로 재표기 또는 시뮬레이션 시간 연장 | **Critical** |
| "T2*" 용어 사용 | Ramsey 시퀀스 없음, 실제로는 Rabi decay time | 용어 오용 | 중간 | 보고서/미팅에서 오해 유발 가능 | Yes — "Rabi decay time"으로 재명명 | High |
| Noise 모델 범위 | Lorentzian + 1/f 지원 확인(API에 이미 존재) | - | - | "Lorentzian만 가능"이라는 전제가 틀렸음을 확인 — 1/f 비교 실험 가능 | 참고만 (제약 아님) | Low |

### 핵심 질문에 대한 답

1. **Noiseless Rabi가 이상적으로 보이는 원인**: 필드 모델 자체가 이상화(균일+선형)되어 있고, drive-amplitude noise 외 다른 noise 경로가 없기 때문.
2. **현재 필드가 실제 micromagnet 시뮬레이션인가**: **아니다** — 순수 analytic 선형 gradient 모델.
3. **MOS_EDSR.py와 EDSR_noise.py의 Hamiltonian/drive 일치 여부**: **불일치** — stale hardcoded 데이터 문제로 확인.
4. **1.6899 vs 1.5177 MHz 차이의 정확한 원인**: RWA 컨벤션이 아니라 `EDSR_noise.py`의 `h0,u`가 현재 `MOS_EDSR.py` 실행 결과와 맞지 않는 **데이터 불일치**.
5. **Lorentzian noise가 작용하는 Hamiltonian 항**: drive amplitude 항(Omega_R → Omega_R(1+beta(t))) — detuning 항은 이 API로 불가능.
6. **Noiseless/noisy 차이가 작은 이유**: (추가 실행 필요 — Suspected Issues 참고) 가능성 있는 원인은 S0=0.01이 작거나, drive-amplitude noise 자체가 Rabi oscillation의 envelope에 미치는 영향이 제한적이기 때문.
7. **100-run 평균에서 noise effect**: 미검증 (Suspected Issues).
8. **82µs fit이 데이터로 constrain되는가**: **아니다** — 창 끝에서 2.5%만 감소, 명백한 외삽.
9. **"T2*" 명명이 적절한가**: **아니다** — Ramsey 시퀀스가 없으므로 최대 "Rabi decay time".
10. **실제 T2*/T2를 얻기 위해 필요한 것**: 별도 Ramsey 시퀀스(T2*), Hahn-echo 시퀀스(T2) 구현.
11. **1/f noise 및 detuning noise 추가 가능성**: 1/f는 API에 이미 존재(`power_law`). Detuning noise는 `Noise.dynamics()` 시그니처로 불가 — custom QuTiP Hamiltonian 필요.
12. **QTCAD API vs custom QuTiP 구분**: Rabi 구동/필드/궤도 계산은 QTCAD, noise의 detuning 경로와 Ramsey/Hahn-echo 시퀀스는 QuTiP로 직접 구현해야 함.

---

## 7. Recommended Experiments

### M06-B: Cross-lever-arm sensitivity (parameter sensitivity study, device prediction 아님)
- **목적**: cross lever arm → CSD 기울기 관계를 정량 확인.
- **변경**: `alpha_cross` = 0.006, 0.02, 0.05, 0.10, 0.15, 0.20 (dominant alpha는 고정).
- **고정**: 그 외 모든 파라미터.
- **예상 변화**: cross alpha가 커질수록 transition line이 대각선에 가까워짐.
- **판정 기준**: honeycomb 셀의 기울기 각도가 cross/dominant 비율과 일치하는지.
- **계산 비용**: 낮음 (기존 energies/Coulomb 재사용, alpha만 수동 교체 — 새 Poisson/Schrödinger 불필요).

### M06-C: overlap=True 비교
- **목적**: 혼성화(anticrossing)가 실제로 복원되는지 확인.
- **변경**: `get_coulomb_matrix(overlap=True)`, `many_body_solver_params.overlap=True`.
- **판정 기준**: triple-point 부근에 곡률이 나타나는지.
- **계산 비용**: 중간 (Coulomb 행렬이 4D로 커짐 — 메모리 확인 필요).

### M08-B: Noise model sensitivity audit
- **목적**: 현재 noise 모델의 민감도 확인 (물리값 예측 아님).
- **변경**: S0 = 0, 0.001, 0.003, 0.01, 0.03, 0.10.
- **비교**: 1-run vs 100-run 평균.
- **계산 비용**: 중간 (100-run 평균을 6개 S0에 대해 — 기존 1회 실행 시간 x 6).

### M08-C: 시뮬레이션 시간 연장 (decay 수렴 확인)
- **목적**: 82µs fit이 실제 데이터로 constrain되는 시점을 찾음.
- **변경**: 20 → 50 → 100 Rabi cycles (필요한 해상도만 유지해 비용 제한).
- **판정 기준**: fit한 decay time이 시뮬레이션 창을 늘려도 안정되는지(진짜 수렴) 아니면 계속 커지는지(여전히 외삽).

### M03-X: Junction-local mesh convergence
- **목적**: wiggle이 메시에 의존하는지 확인.
- **변경**: source/channel, channel/drain 접합 ±5nm만 국소 재정제 (0.5nm 요소 크기 목표).
- **고정**: 나머지 메시, bias, Poisson tolerance.
- **판정 기준**: wiggle 진폭이 재정제 후 감소하면 numerical, 유지되면 physical.

---

## 8. Proposed File Structure (실제 수정 단계 진입 시)

```text
Simulation Example/
├── M03_M06_M08_Physics_Review_Report.md   (이 문서)
├── M03. Gate sweep/Reference/
│   ├── M03_mesh_convergence.py            (신규, 원본 불변)
│   └── M03_Potential_Wiggle_Review.ipynb  (신규)
├── M06. Charge stability/Reference/
│   ├── M06_cross_coupling_study.py        (신규, M06-B)
│   ├── M06_overlap_comparison.py          (신규, M06-C)
│   └── M06_Charge_stability_Review.ipynb  (신규)
└── M08. Micromagnet EDSR/Reference/
    ├── M08_noise_sensitivity.py           (신규, M08-B)
    ├── M08_sim_time_convergence.py        (신규, M08-C)
    └── M08_EDSR_Review.ipynb              (신규)
```

## 9. Modification Order

```text
Step 1 (M08, Critical, 독립적)
  MOS_EDSR.py를 현재 파라미터로 재실행 → h0,u를 다시 추출 →
  EDSR_noise.py의 하드코딩 값을 갱신 (또는 자동 전달 구조로 변경)
  필요 결과: 일치하는 Rabi frequency

Step 2 (M08, Critical, Step 1 이후)
  M08-C 실행 (시뮬레이션 시간 연장) → 82us fit의 타당성 재평가
  필요 결과: Step 1의 올바른 h0,u

Step 3 (M06, High, 독립적)
  M06-C 실행 (overlap=True) → 혼성화 복원 확인
  필요 결과: 없음 (기존 energies 재사용 가능)

Step 4 (M06, High, Step 3과 독립적으로 병행 가능)
  M06-B 실행 (cross lever arm sweep) → CSD 기울기 관계 정량화
  필요 결과: 없음 (기존 energies/Coulomb 재사용)

Step 5 (M03, Medium, 독립적)
  M03-X 실행 (junction-local mesh convergence) → wiggle 원인 확정
  필요 결과: 없음

Step 6 (전체)
  3개 Review notebook 작성 (원본 불변, 신규 파일에 결과 통합)
```
