# QTCAD 체험 커리큘럼

## 폴더 구조
```
QTCAD Simulation/
├── CLAUDE.md
├── CURRICULUM.md
└── Simulation Example/
    ├── 01. SiGe One QD_part.1/      ← 기준 소자 (Philips 2022 적층)
    │   ├── Reference/
    │   └── Study/
    ├── M01. 3D Poisson/
    │   ├── Reference/
    │   │   └── output/
    │   └── Study/
    │       ├── output/
    │       └── notes.md
    ├── M02. Schrodinger/
    └── ...
```

## 진행 규칙
- 모든 주제 폴더는 `Reference/`와 `Study/`를 하나씩 가진다. 결과 파일은 각 폴더 안의 `output/`에 저장한다.
- **Reference/**: 공식 튜토리얼 원본을 복사해 Claude가 실행·디버깅까지 완료한다 (정답지).
- **Study/**: 학습 모드. 기준 소자에 적용한 코드를 Claude가 작성하되,
  - 파이프라인 단계별로 파일을 나눈다 (메시 → 재료·영역 → 경계조건 → Poisson → Schrödinger → 후처리)
  - 한 번에 한 단계씩 작성하고, 사용자가 이해했다고 하면 다음 단계로 넘어간다
  - 각 코드 블록에 물리적 의미와 공식 문서 위치를 주석으로 단다
  - 스크립트를 직접 실행하지 않는다. 실행은 사용자가 하고 결과나 에러를 전달한다
  - 바꿔볼 파라미터는 파일 상단에 모으고, 결과 변화를 예측하는 질문을 함께 준다
  - 같은 주제의 `Reference/` 코드는 사용자가 요청할 때만 보여준다
- Study 코드는 기준 소자 폴더(`01. SiGe One QD_part.1/Study/`)의 메시와 설정을 출발점으로 삼는다.
- 튜토리얼 스크립트와 메시는 설치된 qtcad 패키지의 `examples/tutorials/`에서 찾는다. 이름이 문서와 다르면 목록을 먼저 보여주고 확인받는다.
- API가 불확실하면 추측하지 말고 https://docs.nanoacademic.com/qtcad/ 를 확인한다.
- 주제마다 `Study/notes.md`에 기록한다: 결과 수치, 노드 수, 계산 시간, 미팅 질문에 대한 답.
- 주제를 마치면 해당 모듈의 `- [ ] 완료`를 `- [x] 완료`로 바꾼다.
- 진행 순서: M01 → M02 → M04 → M03 → M05 → M07 → M08 → M06 → (M09~M11 선택)

---

## Phase 1 — 핵심 전기정적 체인 (S 등급)

### M01. 3D Poisson
- [ ] 완료
- 폴더: `Simulation Example/M01. 3D Poisson/`
- 튜토리얼: Device 2 (Poisson solver with adaptive meshing), Device 6 (Band alignment in heterostructures), Atoms 1 (SiGe Part 1)
- 확인된 파일: `adaptive.py`, `band_alignment.py` (M01/Reference/에 복사됨), `multiscale_1.py` (01. SiGe One QD_part.1/Reference/에 있음)
- 실습: 게이트 전압 0 → 목표값 스윕, 스페이서 두께 2~3개 비교
- 확인: φ(r), E_c 라인컷, 양자우물에서 전도대 바닥 위치
- 미팅 질문: 극저온 electrostatics 신뢰 범위
- 준비할 결과: adaptive mesh 수렴 기준과 온도 설정을 바꿔가며 결과 변화량 기록

### M02. Schrödinger
- [ ] 완료
- 폴더: `Simulation Example/M02. Schrodinger/`
- 튜토리얼: Device 4 (Adaptive-mesh Poisson + Schrödinger), Device 5 (Schrödinger equation for a quantum dot), Atoms 1 (SiGe Part 1)
- 확인된 파일: `adaptive_schrodinger.py`, `creating_dot.py` (M02/Reference/에 복사됨), `multiscale_1.py` (01. SiGe One QD_part.1/Reference/에 있음)
- 실습: SubMesh(계산 영역) 크기와 메시 밀도를 각각 3단계로 바꿔 E₀, E₁ 비교
- 확인: 고유에너지, 파동함수, 궤도 에너지 간격
- 미팅 질문: 메시·계산 영역을 어떻게 정해야 정확하고 빠른가
- 준비할 결과: 정확도-계산시간 비교표

### M03. Gate sweep
- [ ] 완료
- 폴더: `Simulation Example/M03. Gate sweep/`
- 튜토리얼: Device 16, 17 (Tunnel coupling in a DQD in FD-SOI, Part 1 plunger / Part 2 barrier)
- 확인된 파일: `tunnel_coupling_1.py`, `tunnel_coupling_2.py` (M03/Reference/에 복사됨)
- 실습: 기준 소자에서 P, B 게이트를 각각 스윕
- 확인: E₀, E₁ 변화, 양자점 중심 위치 이동
- 미팅 질문: 대규모 design-space sweep/optimization 방법
- 준비할 결과: 스윕 1점당 계산시간

### M04. Lever arm
- [x] 완료
- 폴더: `Simulation Example/M04. Lever arm/`
- 튜토리얼: Lever arm 이론 문서 (theory_spin_fem/leverarm) + lever arm 모듈
- 확인된 파일: 별도 튜토리얼 스크립트 없음. `leverarm_fdsoi.py`(Reference에 직접 작성)가
  `qtcad.device.leverarm`(단일 게이트 선형 피팅), `qtcad.device.leverarm_matrix`(5게이트
  유한차분)를 M03의 `get_double_dot_fdsoi` 소자에 적용. plunger_gate_1 레버암 ≈ 0.207,
  cross-coupling 행렬은 공간적으로 타당한 패턴(각 state가 가까운 게이트에 강하게 반응)
- 실습: P 게이트 레버암과 인접 게이트 cross-coupling 행렬 계산. M03 결과의 기울기와 모듈 결과를 교차 확인
- 확인: ∂E/∂V_G, cross-coupling
- 미팅 질문: 실험 CSD로 어떻게 calibration하는가

---

## Phase 2 — 소자 동작 (A 등급)

### M05. Many-body
- [ ] 완료
- 폴더: `Simulation Example/M05. Many-body/`
- 튜토리얼: Device 14 (Many-body analysis of a nanowire QD), Device 19 (Exchange coupling in a DQD, Part 2: exact diagonalization)
- 확인된 파일: `manybody.py`, `exchange_2.py` (M05/Reference/에 복사됨)
- 실습: 기준 소자에서 N = 0~3 전자 addition energy 계산
- 확인: N = 0, 1, 2, 3 에너지, addition energy, Coulomb 상호작용
- 미팅 질문: 실제 single-electron operating window 예측 정확도

### M06. Charge stability
- [ ] 완료
- 폴더: `Simulation Example/M06. Charge stability/`
- 튜토리얼: Transport 1 (Master equation), Transport 3 (Charge stability diagram of a DQD)
- 확인된 파일: `nanowire_diamond.py`(Transport 1, 문서상 이름과 달리 나노와이어+다이아몬드 결함 소자), `double_dot_stability.py`(Transport 3) — M06/Reference/에 복사됨
- 실습: 튜토리얼 CSD 재현 후, 기준 소자를 이중 양자점으로 확장할지 판단
- 확인: single/double dot CSD
- 미팅 질문: tuning 자동화와 연결 가능한가

### M07. Si valley splitting
- [ ] 완료
- 폴더: `Simulation Example/M07. Valley splitting/`
- 튜토리얼: Device 12 (Valley splitting, MVEMT), Atoms 2 (SiGe Part 2, tight-binding)
- 확인된 파일: `valleysplitting_1D_MVEMT.py` (M07/Reference/에 복사됨), `multiscale_2.py` (01. SiGe One QD_part.1/Reference/에 있음)
- 실습: 같은 소자에서 MVEMT와 TB 결과 비교, 계면 폭과 random alloy 시드 변경
- 확인: FEM → atomistic TB 밸리 분리 값과 분포
- 미팅 질문: interface roughness / random alloy 처리 방법

### M08. Micromagnet EDSR
- [ ] 완료
- 폴더: `Simulation Example/M08. Micromagnet EDSR/`
- 튜토리얼: Qubit 1 (EDSR—Dynamics), Qubit 2 (EDSR—Noise), Atoms 3 (SiGe Part 3, quantum control)
- 확인된 파일: `MOS_EDSR.py`(Qubit 1), `EDSR_noise.py`(Qubit 2) — M08/Reference/에 복사됨, `multiscale_3.py` (01. SiGe One QD_part.1/Reference/에 있음)
- 실습: 튜토리얼의 자기장 대신 기존 마이크로마그넷 시뮬레이터의 B-gradient를 넣어 Rabi 주파수 비교
- 확인: B-gradient → Rabi 주파수, quantum dynamics
- 미팅 질문: 외부 micromagnet simulation과 연동 방법
- 준비할 결과: 실제로 연동해 본 결과

---

## Phase 3 — 필요성 판단 (B·C 등급)
목적은 "우리 팀에 필요한가" 판단이다. `Reference/`에서 튜토리얼 실행까지만 하고, `Study/`는 만들지 않는다. 판단 결과는 `Reference/notes.md`에 한 줄로 남긴다.

### M09. g-tensor (B)
- [x] 완료
- 폴더: `Simulation Example/M09. g-tensor/Reference/`
- 튜토리얼: Atoms 4 (Computing g-tensors in an FD-SOI device), Device 10 (Spin–orbit coupling and magnetic effects)
- 확인된 파일: `g_tensor.py`(Atoms 4), `holes_soc_B.py`(Device 10) — M09/Reference/에 복사됨
- 판단 기준: Si/SiGe 전자에서 g 변화가 큐비트 주파수 오차에 의미 있는 크기인가

### M10. Quantum transport (B)
- [ ] 완료
- 폴더: `Simulation Example/M10. Quantum transport/Reference/`
- 튜토리얼: Transport 4 (NEGF–Poisson, FD-SOI FET with a QD), Transport 5 (Poisson–NEGF–Master equation)
- 확인된 파일: `fdsoi_negf.py`(Transport 4), `poisson_negf_master.py`(Transport 5) — M10/Reference/에 복사됨
- 판단 기준: 측정 Coulomb peak와 정량 비교가 필요한가

### M11. Maxwell / CPW (C)
- [ ] 완료
- 폴더: `Simulation Example/M11. Maxwell CPW/Reference/`
- 튜토리얼: Capacitance matrix, Maxwell eigenmode 튜토리얼 (Superconductors 섹션)
- 확인된 파일: `cap_xmon.py`(용량 행렬), `maxwell_xmon.py`(Maxwell eigenmode, inductor ports/EPR), `maxwell_meandered_resonator.py`(적응형 메시 공진기) — M11/Reference/에 복사됨
- 판단 기준: 기존 CPW 워크플로 결과와 한 케이스를 비교해, QTCAD로 옮길 가치가 있는가


---

## Phase 4 — 우리 소자 확장 (S·A 등급)
기준 소자 작업과 바로 연결되는 주제다. Phase 1·2와 같이 `Reference/`와 `Study/`를 둘 다 만든다.
진행 순서 (M11 이후): M12 → M13 → M14 → M15 → M16 → (M17~M20 선택)
M12를 먼저 하는 이유: M13(Tunnel Falls)도 Builder로 메시를 만든다.

### M12. Builder (레이아웃 → 메시) (S)
- [x] 완료 (Ge_hole Part1만; GaAs_gated는 구버전 API라 보류 — notes 참고)
- 폴더: `Simulation Example/M12. Builder/`
- 튜토리얼: QTCAD Builder "Gated Quantum Dot", Practical App Ge hole Part 1 (Builder workflow: Ge/SiGe DQD)
- 확인된 파일: 설치 패키지의 `examples/tutorials/`가 아니라 별도 `examples/practical_application/` 폴더에 있음.
  `Ge_hole/1-builder_ge.py` + `Ge_hole/ge_dqd.oas`(마스크) — Practical App Ge hole Part 1.
  `GaAs_gated/1-devicegen.py` + `GaAs_gated/gated_qd.gds`(마스크) — QTCAD Builder "Gated Quantum Dot" 튜토리얼로 추정.
- 실습: KLayout 게이트 마스크(.gds/.oas)의 폴리곤에 이름을 붙여 기준 소자 메시를 생성 → 손으로 만든 .geo 메시와 E₀, E₁, 노드 수 비교
- 확인: 마스크 → 적층 압출 → dot region → 메시 품질
- 미팅 질문: 다중 큐비트(6-dot) 레이아웃을 그대로 메시로 만들 때 권장 설정

### M13. Tunnel Falls 디튜닝 스펙트럼 (A)
- [x] 완료
- 폴더: `Simulation Example/M13. Tunnel Falls DAPS/`
- 튜토리얼: Practical App "Tunnel Falls Detuning Spectrum" Part 1~4 (Builder 메시 → 소자 생성 → 위치 의존 2-valley k·p → 디튜닝별 에너지 스펙트럼)
- 확인된 파일: `examples/practical_application/Tunnel_Falls_DAPS/` 안에 `1-tunnel_falls_builder.py`(Part 1, Builder 메시),
  `2-energy_vs_detuning.py`(디튜닝별 에너지 스펙트럼), `double_dot_tunnel_falls.py`·`valley_kp_model.py`(소자 생성·위치 의존 2-valley k·p, 헬퍼 모듈로 추정) — 문서상 "Part 1~4"이지만 실제 파일은 2개 번호 스크립트 + 2개 헬퍼 모듈로 구성됨
- 실습: Si/SiGe 이중 양자점에서 valley coupling 값을 바꿔 디튜닝 스펙트럼 변화 확인 → M07 VS 결과를 입력값으로 사용
- 확인: valley-orbit 준위 교차, 디튜닝 축 스펙트럼
- 미팅 질문: 측정한 valley splitting을 시뮬레이션에 반영하는 방법

### M14. 전하 잡음·결함 (A)
- [ ] 완료
- 폴더: `Simulation Example/M14. Charge noise/`
- 튜토리얼: Practical App GaAs "Charge noise in quantum dots", Device 3 (Poisson solver with background charges), Device 20 (Point charges in a DQD in FD-SOI), Atoms 5 (Atomistic disorder in a Si spherical QD in a Ge matrix)
- 확인된 파일: `background_charges.py`(Device 3), `fdsoi_point_charges.py`(Device 20), `spherical_dot.py`(Atoms 5로 추정) — 설치된 `examples/tutorials/`에 있음.
  GaAs "Practical App Charge noise"는 `examples/practical_application/GaAs_gated/8-noise.py` (Part 8) — M12에서 쓰는 `GaAs_gated/` 폴더의 마지막 단계
- 실습: 산화막·계면 고정 전하 밀도와 위치를 바꿔 양자점 에너지·위치 변화 계산 → GHZ dephasing 시뮬레이션의 노이즈 입력값으로 사용
- 확인: 전하 하나당 에너지 이동량, 거리 의존성
- 미팅 질문: 공정 결함 밀도 → 큐비트 편차 예측 방법

### M15. 메시·솔버 정확도 도구 (A)
- [ ] 완료
- 폴더: `Simulation Example/M15. Mesh accuracy/`
- 튜토리얼: Device 22 (Symmetric meshing with Poisson and Schrödinger solvers), Device 21 (Periodic boundary condition in a QD in FD-SOI), Device 9 (Schrödinger simulation of a quantum well), Practical App GaAs Part 3 (quantum-well solver)
- 확인된 파일: `sym_dqdfdsoi.py`(Device 22, 대칭 메시), `periodic.py`(Device 21, 주기경계), `quantum_well_holes.py`(Device 9로 추정) — `examples/tutorials/`에 있음.
  GaAs Part 3 quantum-well solver는 `examples/practical_application/GaAs_gated/3-QW_schrodinger.py` (M12에서 쓰는 폴더의 세 번째 단계)
- 실습: 기준 소자를 대칭 메시로 다시 만들어 E₁/E₂ 축퇴와 단면 아티팩트 개선 확인, z 방향을 quantum-well solver로 처리해 계산 시간 비교
- 확인: 정확도-계산시간 표 (M02 표에 추가)
- 미팅 질문: 대규모 다중 dot 계산에서 메시 비용을 줄이는 방법

### M16. Exchange (섭동론) (A)
- [x] 완료
- 폴더: `Simulation Example/M16. Exchange perturbation/`
- 튜토리얼: Device 18 (Exchange coupling in a DQD in FD-SOI—Part 1: Perturbation theory)
- 확인된 파일: `exchange_1.py` (M05에서 쓴 `exchange_2.py`의 Part 1 짝) — `examples/tutorials/`에 있음
- 실습: M05의 Device 19(exact diagonalization) 결과와 같은 조건에서 J 비교
- 확인: J 대 barrier 전압, 두 방법의 오차와 계산 시간
- 미팅 질문: 스윕에서 어느 방법을 써야 하는가

---

## Phase 5 — 필요성 판단 (B·C 등급)
Phase 3과 같이 `Reference/`에서 튜토리얼 실행까지만 하고, `Study/`는 만들지 않는다. 판단 결과는 `Reference/notes.md`에 한 줄로 남긴다.

### M17. 전하 센서 SET (B)
- [ ] 완료
- 폴더: `Simulation Example/M17. Charge sensor SET/Reference/`
- 튜토리얼: Practical App "FD-SOI SET" Part 1~6 (Builder 메시, 선형 Poisson, 화학 퍼텐셜, Coulomb peak, 주변 dot 전하 검출)
- 확인된 파일: `examples/practical_application/FDSOI/` 안에 `1-fdsoi_builder.py`, `2-linear_poisson.py`, `3-coulomb_peaks.py`,
  `4-charge_detection.py`, `chemical_potential.py`, `double_dot_fdsoi.py`(헬퍼) — 문서상 "Part 1~6"이지만 실제 번호 스크립트는 4개 + 헬퍼 2개로 구성됨
- 주의: 문서 기준 약 1.5시간 (전자 약 12개 Schrödinger–Poisson 수렴)
- 판단 기준: 센서 dot 감도를 시뮬레이션으로 설계할 가치가 있는가

### M18. Ge 정공 DQD 전체 흐름 (B)
- [ ] 완료
- 폴더: `Simulation Example/M18. Ge hole DQD/Reference/`
- 튜토리얼: Practical App Ge hole Part 2~4 (Poisson·Schrödinger → lever-arm matrix → Coulomb 상호작용·CSD)
- 확인된 파일: `examples/practical_application/Ge_hole/` 안에 `2-poisson_schrod.py`, `3-leverarm.py`, `4-CSD.py`
  (+ 결과 대조용 `addition_spectrum.txt`) — Part 1 메시는 M12의 `1-builder_ge.py`+`ge_dqd.oas`에서 생성
- 판단 기준: 헤테로구조 DQD 전 과정을 우리 소자(전자)로 옮길 때 참고 흐름으로 쓸 만한가

### M19. Strain (C)
- [ ] 완료
- 폴더: `Simulation Example/M19. Strain/Reference/`
- 튜토리얼: Device 11 (Including strain in a simulation)
- 확인된 파일: `holes_strain.py` — `examples/tutorials/`에 있음
- 판단 기준: 인장 변형 Si 우물의 밴드 오프셋·valley 변화가 결과에 의미 있는 크기인가

### M20. WKB 터널링 (C)
- [ ] 완료
- 폴더: `Simulation Example/M20. WKB/Reference/`
- 튜토리얼: Transport 2 (Quantum transport—WKB approximation)
- 확인된 파일: `WKB_current.py` — `examples/tutorials/`에 있음
- 판단 기준: barrier 전압별 터널 레이트 추정이 CSD 해석(M06)에 필요한가

---

## 제외한 튜토리얼
Device 1(나노와이어 입문), Device 7(ParaView, 필요할 때만 참고), Device 8(MOS 커패시터), Device 13(도너 MVEMT), Device 23(Jackiw–Rebbi k·p), Qiskit Metal 렌더러(M11 판단이 "옮길 가치 있음"일 때만 추가)
