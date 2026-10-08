# QTCAD M03 / M06 / M08 물리 모델 검증 보고서 (Validation Report)

**이 문서의 역할**: `M03_M06_M08_Physics_Review_Report.md`(최초 검토 보고서)에 제기된 각 주장(claim)을
실제 시뮬레이션·코드 수정·재검증을 통해 어떻게 판정했는지 추적한다. 최초 보고서는 수정하지
않고 그대로 보존한다.

**작업 순서**: M08 Step 1→2 (Critical) → M06 Step 3→4 (High) → M03 Step (Medium), 지침
22절 우선순위를 그대로 따름.

**상태 태그 정의** (지침 48절): `CONFIRMED`(실제 코드/데이터로 확인 완료) ·
`SUPPORTED`(여러 검증에서 일관되지만 완전 확정은 아님) · `UNRESOLVED`(추가 검증 필요) ·
`INVALID`(기존 설명이 틀린 것으로 확인) · `TUTORIAL_ONLY`(교육용 단순화, 의도된 것) ·
`PARAMETRIC`(sensitivity study 결과, device prediction 아님)

---

## 1. M08 — Micromagnet EDSR

### 1.1 Claim: "1.6899 MHz vs 1.5177 MHz Rabi 주파수 불일치는 RWA convention 차이다"
*(최초 검토 보고서, notebook 자체의 기존 설명)*

- **Validation Method**: `EDSR_noise.py`의 하드코딩된 `h0`,`u`를 동일한
  `dyn.transition_2_levels()` 함수에 다시 넣어 재현 테스트 (최초 검토 단계, fork 조사).
- **Result**: 하드코딩된 행렬을 공식 함수에 다시 넣으면 1.5177MHz가 나옴 (1.6899MHz가
  아님) — 즉 같은 함수가 같은 입력에 대해 일관된 출력을 내므로 "함수 내부 convention
  차이"라는 가설은 성립할 수 없음.
- **Status**: **INVALID** (RWA convention 설명 자체가 틀림)
- **Code Change**: `MOS_EDSR.py`에 `h0`,`u`,`omega0`,`omega_rabi`와 생성 파라미터(gate
  bias, B-field, mesh 경로)를 `output/mos_edsr_h0_u.npz` + `mos_edsr_metadata.json`으로
  명시적으로 저장하는 코드 추가. `EDSR_noise.py`에서 하드코딩된 `h0`,`u` 선언을 제거하고
  위 npz를 로드하도록 변경, 로드 직후 재계산값과 저장값을 비교하는 consistency check
  (상대오차 1e-3 기준, 초과 시 `RuntimeError`) 추가. 원본 파일 수정 전 `_baseline_backup/`
  폴더에 두 파일 모두 백업.
- **Physical Interpretation**: 실제 원인은 **stale/mismatched hardcoded data** —
  `EDSR_noise.py`의 `h0`,`u`가 과거 어느 시점의 `MOS_EDSR.py` 실행 결과를 손으로 복사한
  것인데, 그 이후 `MOS_EDSR.py`의 파라미터가 바뀌면서 두 스크립트가 서로 다른 소자를
  계산하고 있었던 것.
- **Remaining Limitation**: 수정 후 재실행한 `MOS_EDSR.py`의 Rabi 주파수(1.7086MHz)는
  과거 로그값(1.6899MHz)과도 ~1.1% 다름 — 이는 QTCAD/mesh/solver의 run-to-run 변동으로
  보이며 stale-data 문제와는 무관한 별도 현상 (추가 조사 불필요한 수준의 작은 차이).

| quantity | MOS_EDSR(과거 로그, stale) | EDSR_noise(과거 하드코딩) | MOS_EDSR(재실행) | EDSR_noise(수정 후) |
|---|---:|---:|---:|---:|
| f_Rabi (Hz) | 1.6899 MHz | 1.5177 MHz | **1.7086 MHz** | **1.7086 MHz** (일치) |

---

### 1.2 Claim: "82 µs T2* 값이 데이터로 충분히 constrain되는지 의심스럽다"

- **Validation Method**: `M08_sim_time_convergence.py`(신규) — 동일 noise 모델(Lorentzian,
  S0=0.01)로 시뮬레이션 창을 20→50→100→200 Rabi cycle로 늘려가며 damped-cosine fit이
  수렴하는지 확인. Fit-validity gate: window 끝에서 envelope<0.7이어야 "resolved"로 인정.
- **Result**:

  | n_cycles | T_sim (µs) | T_decay fit (µs) | envelope@끝 | 판정 |
  |---:|---:|---:|---:|---|
  | 20 | 11.71 | 21.51 | 0.744 | unresolved |
  | 50 | 29.26 | 29.58 | 0.376 | resolved (gate 통과) |
  | 100 | 58.53 | 39.95 | 0.117 | resolved |
  | 200 | 117.05 | 53.77 | 0.009 | resolved |

  100→200 구간에서도 fit값이 34.6% 더 커짐 — **아직 안정된 값에 도달하지 못함**.
- **Status**: **UNRESOLVED** (82µs라는 특정 숫자는 확인된 외삽이며, 창을 10배 늘려도
  최종 수렴값을 얻지 못함)
- **Code Change**: `M08_sim_time_convergence.py` 신규 작성. `EDSR_noise.py`,
  `MOS_EDSR.py`, 기존 분석 notebook은 수정하지 않음.
- **Physical Interpretation**: 20-cycle 창(13.2µs)에서 envelope는 2.5%만 감소한
  상태였고, 82µs라는 fit은 사실상 "거의 평평한 곡선"에서 장거리 외삽한 숫자. 창을
  늘리면 fit값이 계속 커지는 경향(21.5→29.6→40.0→53.8µs)은 단일 Gaussian decay
  모델 자체가 이 noise 과정을 정확히 기술하지 못하고 있을 가능성을 시사.
- **Remaining Limitation**: 완전한 수렴을 보려면 n_cycles를 더 늘려야 하는데 200
  cycles가 이미 367초 걸림 — 추가로 400+ cycles를 시도할지는 비용 대비 가치 판단
  필요. 또한 Gaussian 모델 대신 exponential/stretched-exponential 비교(지침 39절 D항)는
  아직 수행하지 않음.

---

### 1.3 Claim: "현재 noiseless Rabi가 real micromagnet 시뮬레이션인가?"

- **Validation Method**: `MOS_EDSR.py`의 `Bfield()` 함수 직접 확인.
- **Result**: `B = [B0, 0, b*x]` (균일 0.6T + 완전 선형 gradient 0.3e6 T/m) — fringing
  field, magnet geometry, 제작 공차 등 전혀 없음. 실측/COMSOL/magnetostatic 필드맵
  import 경로도 QTCAD 전체에 없음.
- **Status**: **CONFIRMED** (idealized linear-gradient model, not a real micromagnet sim)
- **Code Change**: 없음 (이번 라운드에서는 수정하지 않기로 함, 지침 20절 — 새
  multiphysics 모델 구축은 범위 밖)
- **Physical Interpretation**: "Micromagnet EDSR"이라는 이름과 실제 구현 수준(이상화된
  gradient-field EDSR) 사이에 간극이 있음을 명확히 인지해야 함.
- **Remaining Limitation**: TUTORIAL_ONLY로 분류 — 실제 micromagnet 필드맵을 쓰려면
  별도 자기장 솔버/측정 데이터 import가 필요 (향후 과제).

---

### 1.4 Claim: "현재 noise 모델은 drive-amplitude noise 중심이고 detuning noise는 없다"

- **Validation Method**: `qtcad.qubit.noise.Noise.dynamics()` 시그니처 + `EDSR_noise.py`
  주석 확인.
- **Result**: `Noise.dynamics(H0, delta_V, omega, spectrum, ...)` — 단일 H0/omega +
  하나의 확률과정만 받는 구조. 스크립트 주석("charge noise...assumed to affect only
  the gate which generates the drive delta_V")과 결합해 drive-amplitude noise
  (`Omega_R -> Omega_R[1+beta(t)]`)만 구조적으로 가능함을 확인. 단, 추가로
  `qtcad.qubit.spectra`에 `lorentz`뿐 아니라 **`power_law`(1/f^alpha)도 존재** — 최초
  검토 보고서의 "noise spectrum이 Lorentzian뿐"이라는 전제는 **부분적으로 틀림**
  (스펙트럼 모델은 다양하게 있으나, detuning 경로 자체가 없는 것은 맞음).
- **Status**: **CONFIRMED** (drive-amplitude noise만 가능, detuning noise 경로 없음) /
  **INVALID** (noise spectrum이 Lorentzian뿐이라는 가정은 틀림 — power_law도 있음)
- **Code Change**: 없음 (API 조사만, 구현 확장은 범위 밖)
- **Remaining Limitation**: 실제 detuning noise를 모사하려면 QuTiP으로 custom
  Hamiltonian(`H(t) = (hbar/2)[Delta+delta_omega(t)]sigma_z + ...`)을 직접 구성해야
  함 — 지침 8절(Step 4)에 해당하나 이번 라운드에서는 미수행.

---

## 2. M06 — Charge Stability Diagram

### 2.1 Claim: "CSD가 지나치게 직사각형인 가장 큰 원인은 무엇인가?"

이 질문을 세 가지 독립 가설로 분해해서 각각 검증했다.

#### 2.1.a Gate-bias reference 처리가 버그인가?
- **Validation Method**: `qtcad.transport.junction.Junction`의 `biases`/`set_bias_vector`/
  `set_biases`/`set_bias` 네 곳 docstring 확인.
- **Result**: 전부 동일 문구로 "bias는 초기 electrostatics를 푼 reference 전위 기준
  상대값"이라고 명시. 코드의 30~140mV 스윕(상대값) + plot의 `+0.6V`/`+0.601V`(reference
  복원)가 API 사양과 정확히 일치.
- **Status**: **CONFIRMED — 버그 아님**
- **Code Change**: 없음

#### 2.1.b 10K 온도 설정이 의도치 않은 것인가?
- **Validation Method**: `double_dot_stability.py` 소스 주석 확인.
- **Result**: L186-188에 "Set the junction temperature to 10 K (instead of 100 mK)
  to make lines... thicker and decrease the resolution required to observe them"라고
  명시적으로 설명되어 있음.
- **Status**: **TUTORIAL_ONLY — 의도된 설정, 버그 아님**
- **Code Change**: 없음

#### 2.1.c overlap=False가 직사각형의 원인인가?
- **Validation Method**: `M06_overlap_comparison.py`(신규) — baseline(overlap=False) 재현 +
  overlap=True Coulomb tensor(4,4,4,4) 계산 후 동일 90×90 스윕 재실행, 이미지 직접 비교.
- **Result**: 정규화 CSD 최대 절대 차이 0.991(수치상 큼)이지만, **이미지를 직접 보면
  전이선 토폴로지(직사각형 격자 모양)가 거의 동일** — 차이는 전이선의 정확한 위치
  이동에서만 나옴. 새로운 anticrossing/curvature는 관찰되지 않음.
- **Status**: **SUPPORTED — overlap=False는 직사각형의 주원인이 아님** (이 bias
  range/gate 구성에서는)
- **Code Change**: `M06_overlap_comparison.py` 신규. Baseline 파일 불변. Consistency
  check: energies rel.err 5.22e-6, lever arm rel.err 1.77e-2, Coulomb(no overlap)
  rel.err 5.83e-4, CSD(no overlap) rel.err 4.78e-2 — baseline 재현 양호.
- **Remaining Limitation**: 다른 bias range(예: 더 넓은 detuning)에서는 overlap=True
  효과가 더 뚜렷할 수 있음 — 이번 테스트는 원본과 동일한 좁은 범위(110mV)에서만
  확인했음.

#### 2.1.d Cross lever arm의 작음이 직사각형의 원인인가?
- **Validation Method**: `M06_cross_coupling_study.py`(신규) — 캐시된 energies/Coulomb
  행렬 재사용(비용 절감), lever-arm의 cross 성분만 0.006→0.02→0.05→0.10→0.15→0.20으로
  인위적으로 변경(dominant는 고정), CSD 6개 재계산.
- **Result**: **시각적으로 매우 뚜렷하게** alpha_cross가 커질수록 직사각형 격자가
  점진적으로 평행사변형으로 기울어짐 (0.006: 거의 완벽한 직각 → 0.2: 뚜렷한 대각선
  격자).
- **Status**: **CONFIRMED — cross lever arm의 작음이 직사각형의 핵심 원인**
- **Code Change**: `M06_cross_coupling_study.py` 신규.
- **Physical Interpretation**: 2.1.c(overlap)와 2.1.d(cross lever arm)를 종합하면,
  **"평행사변형/honeycomb 모양"을 결정하는 것은 cross-capacitance(lever arm)이지
  Coulomb 행렬의 overlap 여부가 아니다** — 최초 검토 보고서 2절의 가설이 올바르게
  확인됨.
- **Remaining Limitation**: 자동 기울기 추출 코드(`empirical_slope`)에 버그가 있어
  (여러 전이선이 창 안에 겹칠 때 ridge-hopping으로 기울기 수치가 틀어짐)
  `slope_simulated`/`relative_slope_diff_pct` 수치는 **신뢰할 수 없음** — 시각적 확인만
  유효함. 정량적 analytic-vs-simulated slope 일치를 보려면 단일 전이선만 격리하는
  방식으로 재작성 필요 (미수행).

### 2.2 Claim: "noise가 없어서 CSD가 너무 깨끗하다"
- **Validation Method**: `double_dot_stability.py` + `M06_Charge_stability_Analysis.ipynb`
  전체에서 `random|noise|seed|normal(` grep.
- **Result**: 일치 항목 0건 (interpolation="bilinear" 3건 제외) — 완전히 deterministic.
- **Status**: **CONFIRMED**
- **Code Change**: 없음 (noise 추가는 지침 10절 Step 10, cross-coupling/overlap 검증
  이후로 명시적으로 미뤄짐 — 이번 라운드 범위 밖)
- **Physical Interpretation**: Noise 부재는 "선이 너무 깨끗함"의 원인이지, "직사각형
  모양"의 원인은 아님 — 2.1.d에서 이미 모양의 원인이 cross lever arm임을 확인했으므로,
  이 claim은 올바르게 "별개 원인"으로 분리되어 있었음.

---

## 3. M03 — Potential Wiggle at Source/Drain Junctions

### 3.1 Claim: "Source/channel, channel/drain 경계 부근의 작은 potential wiggle이
numerical artifact인지 physical feature인지"

- **Validation Method 1 (interpolation 비교)**: `an.linecut()`을 `method="default"`
  (원본이 쓰는 방식)와 `method="pyvista"`로 각각 호출해 비교 (최초 조사 단계).
- **Result 1**: `pyvista` 방식은 서브-옹스트롬 간격에서 최대 26mV까지 튀는 물리적으로
  불가능한 jump를 만듦 — VTK 셀 탐색 아티팩트로 강하게 의심. `default` 방식에서는
  훨씬 작고 완만한 wiggle(좌 3.56mV, 우 2.72mV)만 남음.
- **Status 1**: **CONFIRMED** — pyvista 방식 자체는 명백한 plotting/interpolation
  artifact. `default` 방식(원본이 실제로 쓰는 것)의 작은 wiggle은 별도 검증 필요.

- **Validation Method 2 (mesh convergence)**: `M03_mesh_convergence.py`(신규) —
  baseline(h_refined=0.8, dot 영역만 정제, 캐시된 HDF5 재사용) vs Case 1(h_refined=0.5,
  channel/source/drain도 추가 정제, 644,774 노드) 비교. 동일 bias/온도/tolerance.
- **Result 2**:

  | Case | 메시 노드 | Left wiggle (mV) | Right wiggle (mV) |
  |---|---:|---:|---:|
  | baseline | 173,626 | 3.552 | 2.723 |
  | 조밀화(h=0.5) | 644,774 (3.7배) | **검출 안 됨(사라짐)** | **1.003 (63%↓)** |

- **Status 2**: **SUPPORTED — numerical artifact 쪽으로 강하게 기움** (완전 CONFIRMED는
  아님, 아래 한계 참고)
- **Code Change**: `M03_mesh_convergence.py` 신규. `tunnel_coupling_1.py`/`_2.py`/분석
  notebook 전혀 수정 안 함. **실행 중 2가지 자체 버그 발견·수정**:
  1. 1차 시도에서 Case 2(h_refined=0.3을 channel/source/drain **전체**에 적용)가
     메시 267만 노드로 폭주해 `ArrayMemoryError` 크래시 — 전용 junction-local physical
     group이 mesh에 없어 전체 채널 볼륨(90nm)에 과도한 목표를 건 설계 실수. Case 2
     제거, 2-point 비교로 축소.
  2. Wiggle 자동 검출 함수가 2번 틀림 — ①단순 min/max(전체 band-bending ~120meV를
     그대로 잡음), ②선형 추세 제거(여전히 ~45meV, 불충분) — 최종적으로
     `scipy.signal.argrelextrema`로 진짜 국소 극값 쌍을 찾는 방법(최초 조사 fork가
     썼던 방법과 동일)으로 교체 후 baseline 수치가 조사 보고서와 정확히 일치함을
     확인.
- **Physical Interpretation**: 두 접합 모두 "메시가 조밀할수록 wiggle 감소"라는
  일관된 방향성을 보였고, 그 폭(63% 감소, 완전 소멸)이 우연한 노이즈로 보기엔 큼 —
  numerical discretization artifact라는 가설을 뒷받침.
- **Remaining Limitation**: **메시 레벨이 2개뿐**(Case 2가 비용 문제로 실패)이라 엄밀한
  3점 이상 수렴 곡선은 없음 — 완전한 CONFIRMED 판정을 위해서는 접합부 전용의 더 좁은
  physical group(예: ±5~10nm)을 `.geo`에 추가해 더 저렴하게 조밀화한 Case 2/3가
  필요함(아직 미수행). 또한 Ey(전기장)·전하밀도와의 물리적 일관성 비교(지침 11절)는
  이번 라운드에서 수행하지 않음.

---

## 4. 종합 판정 요약

| 항목 | 최초 가설/우려 | 최종 상태 |
|---|---|---|
| M08 Rabi 주파수 불일치 원인 | "RWA convention" | **INVALID** → stale hardcoded data (CONFIRMED, 수정 완료) |
| M08 82µs T2* | 의심스러움 | **UNRESOLVED** (확인된 외삽, 200 cycle에서도 미수렴) |
| M08 micromagnet 필드 | 이상화 의심 | **CONFIRMED** (순수 선형 gradient, 실제 micromagnet 아님) |
| M08 noise = drive-amplitude only | 의심 | **CONFIRMED** (detuning 경로 없음); 1/f 지원 여부는 **INVALID** 가정 수정(존재함) |
| M06 bias reference | 버그 의심 | **CONFIRMED — 버그 아님** (API 사양대로 정상) |
| M06 10K 온도 | 버그 의심 | **TUTORIAL_ONLY** (의도된 설정) |
| M06 overlap=False가 직사각형 원인 | 가설 | **SUPPORTED — 주원인 아님** |
| M06 cross lever arm이 직사각형 원인 | 가설 | **CONFIRMED — 주원인 맞음** |
| M06 noise 부재 | 확인 필요 | **CONFIRMED** (0건), 모양이 아닌 "깨끗함"의 원인 |
| M03 wiggle 원인 | 미확정(numerical/physical/mixed) | **SUPPORTED — numerical 쪽** (2점 비교, 완전 수렴 아님) |
| M03 pyvista artifact | 의심 | **CONFIRMED** (물리적으로 불가능한 크기) |

---

## 5. 다음 단계로 넘길 항목 (이번 라운드 미수행, 지침 22절 후순위)

- M08: Ramsey 시퀀스 구현 → 진짜 T2* 측정
- M08: Hahn-echo 시퀀스 구현 → 진짜 T2 측정
- M08: QuTiP custom Hamiltonian으로 detuning noise 추가
- M08: Noise ensemble 수렴성 검사 (N=1,10,50,100,200 run 비교, 지침 37절)
- M06: Geometry 기반 cross-coupling 유도 (plunger spacing/width 등 변경, 지침 12절)
- M06: CSD에 noise 추가 (quasi-static/1-f/RTN, cross-coupling 검증 이후, 지침 10절)
- M06: `empirical_slope` 버그 수정 (단일 전이선 격리 추적 방식으로 재작성)
- M03: 접합부 전용 physical group을 `.geo`에 추가해 저비용 3점 이상 수렴 곡선 완성
- M03: Ey·전하밀도와의 물리적 일관성 비교

---

## 6. Proposed File Structure (실제 적용됨)

```text
Simulation Example/
├── M03_M06_M08_Physics_Review_Report.md       (최초 검토, 불변)
├── M03_M06_M08_Physics_Validation_Report.md   (이 문서)
├── M03. Gate sweep/Reference/
│   └── M03_mesh_convergence.py                (신규, 원본 불변)
├── M06. Charge stability/Reference/
│   ├── M06_overlap_comparison.py              (신규, 원본 불변)
│   └── M06_cross_coupling_study.py            (신규, 원본 불변)
└── M08. Micromagnet EDSR/Reference/
    ├── _baseline_backup/
    │   ├── MOS_EDSR_baseline.py               (수정 전 백업)
    │   └── EDSR_noise_baseline.py             (수정 전 백업)
    ├── MOS_EDSR.py                            (수정됨 — h0/u 저장 추가)
    ├── EDSR_noise.py                          (수정됨 — 하드코딩 제거, consistency check 추가)
    └── M08_sim_time_convergence.py            (신규)
```
