# M07. Valley splitting — notes

## 폴더 구성
```
M07. Valley splitting/
├─ notes.md                      ← 이 파일
├─ Reference/                    ← Device 12: 1D MVEMT (Si/SiO2), Gamble 2016 Fig. S1 재현
│  ├─ valleysplitting_1D_MVEMT.py   (재실행 시 결과 덮어쓰기 + 점별 시간 기록 추가)
│  ├─ analysis.ipynb                (Part A: MVEMT 해석 정정·계산 시간 / Part B: TB 원자 구조 분석)
│  └─ TB_multiscale_2/README.md     ← SiGe Part 2 (TB) 정답지 실행 안내
└─ Study/                        ← 기준 소자(Si/SiGe)로 옮긴 단계별 코드
   ├─ README.md
   ├─ atoms_builder.py
   ├─ 01_mvemt_sige_well.py
   ├─ atoms_plot.py
   ├─ 02_tb_atoms_distribution.py
   ├─ 02b_plot_atoms_from_xyz.py
   └─ 03_tb_seed_ensemble.py
```

---

## 1. Reference (Device 12, 1D MVEMT) 결과

구조: Si 영역 13 nm + SiO2 장벽 3 nm(+3 eV), 계면 1개(계단), 1D 메시 1,601 노드(h = 0.01 nm).
전기장 5~50 MV/m, 10점 × 근사 3종.

| 항목 | 값 |
|---|---|
| Full VS 범위 | 0.105 – 1.058 meV |
| ff 근사, full 대비 평균 오차 | 1.85 % |
| trivial / full 비 | 4.86배 (E에 거의 무관) |
| 전기장 지수 (VS ∝ Eᵖ) | p = 1.004 (full), 0.999 (trivial) |
| 1점당 풀이 시간 중앙값 (run3) | trivial 28 s / ff 50 s / full 138 s |
| ff의 full 대비 속도 | 약 2.8배 (중앙값 기준) |

- 계산 시간은 30점을 모두 끝낸 `valleysplitting_run3.log` 기준 (총 2,391 s).
  run1은 점당 420~760 s로 다른 조건(동시 작업 등)으로 보여 제외.
  같은 근사 안에서도 점마다 40~160 s로 흔들림 → 미팅용 수치는 중앙값.
- 결론: 파라미터 스윕은 ff 근사로 충분 (오차 2 % 미만, 약 3배 빠름).

### 정정한 해석
- wavefunction 그림에서 큰 성분은 **±x valley** (1D 메시 축 = 구속 방향 = x).
  ±y, ±z는 진폭 ~1 수준의 수치 잡음. 초판의 "±z valley" 해석은 틀렸음.
- VS ∝ E는 계면 1개인 삼각형 우물의 표준 결과. 초판의 "resonant 영역" 언급 삭제.
- 결과 txt가 덧쓰기 모드라 재실행 시 행이 쌓이던 문제 → 스크립트 시작 시 삭제하도록 수정.

---

## 2. Reference/TB_multiscale_2 (SiGe Part 2) 결과
- [ ] 실행 완료
- VS: ____ meV (정답 1.4855 meV)
- 원자 수: ____ / ____ (정답 896,337 / 87,189)
- 시간: 구조 생성 ____ s, Keating ____ s, TB ____ s
- 최대 메모리: ____ GB
- analysis.ipynb Part B 결과
  - 장벽 Ge 평균 / 층간 표준편차: ____ % / ____ % (이항분포 예상 ____ %)
  - 계면 위치(50 %)·10–90 % 폭: 상부 ____ nm / ____ nm, 하부 ____ nm / ____ nm
  - 결합 길이 Si–Si / Si–Ge / Ge–Ge: ____ / ____ / ____ pm
  - 원자층 간격 우물 / 스페이서: ____ / ____ pm
  - 스페이서 국소 조성 요동(1 nm 셀): ____ % (이항분포 ____ %)
- VESTA 관찰:

---

## 3. Study 계획과 예측 질문 (실행 전에 답부터 적기)

### 01 — 1D MVEMT, Si/SiGe 8 nm 우물
바꾸는 것: ΔEc 3 → 0.15 eV, 계면 1 → 2개, E 1~10 MV/m, 계면 폭 w = 0 / 0.25 / 0.5 / 1 nm.

- Q1. 계단 계면(w = 0), E = 5 MV/m에서 VS는 Reference(SiO2, 같은 E: 0.105 meV)보다 클까 작을까? 이유는?
  - 내 예측:
  - 결과:
- Q2. w를 0 → 1 nm로 늘리면 VS는 대략 몇 배 줄어들까? (힌트: 2k0 주기 약 0.32 nm)
  - 내 예측:
  - 결과:
- Q3. w가 클 때도 VS가 E에 비례할까?
  - 내 예측:
  - 결과:
- 검증 1점 (ΔEc = 3 eV, w = 0, E = 50 MV/m, ff): VS ____ meV (Reference 1.037 meV), ⟨x⟩ ____ nm

### 02 — 원자 배열·조성 분포
- Q4. 10 × 10 nm 원자층 하나(약 680개)의 Ge 비율은 30 %에서 얼마나 흔들릴까? 1 nm³ 셀(약 50개)은?
  - 내 예측:
  - 결과 (z 프로파일 / xy 지도 표준편차):
- Q5. rough 모드에서 우물 안 Ge 원자 수는? graded(w = 0.5 nm) 모드에서는?
  - 내 예측:
  - 결과:
- 배치도 관찰 (`*_side_zoom.png`, `*_top_layers.png`, `*_lattice_top_interface.png`):
- 튜토리얼 `multiscale.xyz` 배치도 (02b): Keating 이완 후 SiGe 층 부풂, 거친 계면 모양

### 03 — 시드별 TB VS 분포
- 양자점 위치 수직 전기장 |E| = ____ MV/m (`03_field.txt`)
- 같은 구조 MVEMT (01을 WELL_WIDTH_NM = 3, w = 0, E = 위 값으로): VS ____ meV
- Q6. 시드 5개에서 VS의 상대 표준편차는 몇 %쯤일까? TB 평균은 MVEMT보다 클까 작을까?
  - 내 예측:
  - 결과: 평균 ____ meV, 표준편차 ____ meV
- Q7. 바닥 상태가 Ge 위에 있는 확률 P(Ge)와 VS 사이에 경향이 있을까?
  - 내 예측:
  - 결과:
- 시드 1개당 시간: 생성 ____ / Keating ____ / TB ____ s

---

## 4. 확인이 필요한 API 가정 (처음 실행 때 체크)
- [ ] 1D `Device.set_V`가 tanh 함수(배열 입력)를 그대로 받는가 (01)
- [ ] `atom_species` 자료형(문자열/원자번호), `atom_pos` 단위(m) (02, 03)
- [ ] `Atoms.new_region`에 위치 인자 없이 농도 함수 배열을 넘기면 구조 전체에 적용되는가 (graded 모드)
- [ ] `save_vtu`의 `out_dict`에 원자별 0/1 배열을 넘길 수 있는가 (02)
- [ ] `num_states = 4`로 줄여도 seed 104에서 VS 1.4855 meV가 그대로 나오는가 (03)
- [ ] `PART1_DIR` 경로 (03)

---

## 5. 미팅 질문에 대한 내 답
**Q. interface roughness / random alloy를 어떻게 처리하나?**
- QTCAD 방식 (확인한 것):
  - 거칠기: `RoughSurface` (Hurst 지수, 평균 높이, rms) → 층 경계로 사용
  - 합금: `new_layer` / `new_region`에 조성(상수 또는 위치 함수) → 원자별 무작위 배치, `rng_seed`로 재현
  - 변형: Keating VFF 이완
  - 전자 구조: FEM 퍼텐셜을 원자에 옮긴 뒤 sp³d⁵s* TB
- 우리 소자 기준 결과 (Study 후 작성):
- 남은 질문 (Nanoacademic에 물어볼 것):
  - 계면 거칠기와 조성 기울기를 동시에 넣는 방법 (rough 계면 + graded 조성)
  - 시드 앙상블 계산 비용을 줄이는 권장 방법 (SubAtoms 크기, num_states, 병렬화)
  - Philips 적층(8 nm 우물, 30 nm 스페이서)에서 TB 상자 크기 권장값
