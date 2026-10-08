# M07 Study — 실행 순서

모든 스크립트는 파일 맨 위 "바꿔볼 파라미터" 블록만 수정한다.
실행 전에 `../notes.md`의 예측 질문에 먼저 답을 적는다.

| 순서 | 파일 | 내용 | 필요 조건 | 대략 시간 |
|---|---|---|---|---|
| 1 | `01_mvemt_sige_well.py` | 1D MVEMT, Si/SiGe 8 nm 우물, 계면 폭 × 전기장 스윕 | 없음 | 17점 × 수십 초~수 분 |
| 2 | `02_tb_atoms_distribution.py` | 무작위 합금 원자 배열·조성 분포 확인 | 없음 | 수십 초 (Keating 끄면) |
| 2b | `02b_plot_atoms_from_xyz.py` | 저장된 `.xyz`(튜토리얼 `multiscale.xyz` 포함)로 배치도만 그림 | `.xyz` 파일 (QTCAD 불필요) | 수십 초 |
| 3 | `03_tb_seed_ensemble.py` | 시드별 TB VS 분포, MVEMT 비교용 전기장 | Part 1 결과, Reference/TB_multiscale_2 성공 | 시드당 튜토리얼 1회분 |

`atoms_builder.py`(원자 구조 생성)와 `atoms_plot.py`(배치도 그리기)는 공통 모듈이다 (직접 실행하지 않음).

## 원자 배치·격자 구조 그림 (2, 2b 실행 후 `output/`)
| 파일 | 내용 |
|---|---|
| `*_side_full.png` | 측면 단면도 (x–z, 두께 0.3 nm 슬라이스) — 전체 상자 |
| `*_side_zoom.png` | 측면 단면도 — 우물 ± 3 nm 확대. 우물(파랑만)과 장벽(Ge 빨강 섞임), 계면 거칠기 |
| `*_top_layers.png` | 원자층 하나씩 평면도: 스페이서 / 상부 계면 / 우물 중앙 |
| `*_lattice_top_interface.png` | 상부 계면 2.5 × 2.5 nm 격자 확대도 — 다이아몬드 격자 결합선, Si/Ge 원자 |
| `*_ge_3d.png` | 계면 부근 Ge 원자만 3D |
| `*_z_profile.png`, `*_xy_map.png` | 조성 분포(원자층별 Ge 비율, xy 국소 조성 지도) |

## 처음 실행할 때 확인할 출력
- 01: `[검증]` 줄의 VS가 1.037 meV 근처이고 `<x>`가 0에 가까운 음수인가
- 02: `atom_species 예시`, `atom_pos 예시`가 예상한 자료형·단위(m)인가.
  `positions_nm()`이 단위 오류를 내면 그 출력을 보고 변환을 고친다.
- 02: graded 모드가 오류 없이 도는가 (`new_region`에 농도 함수를 넘기는 방식은
  Atoms 5 튜토리얼과 같지만, 층 구조에서 위치 인자 없이 쓰는 것은 문서에 직접 예시가 없다)
- 03: seed 104의 `[검증]` 차이가 0 %에 가까운가 (튜토리얼과 같은 상자·시드·생성 순서)

## 결과 파일 (`output/`)
- `01_mvemt_vs.csv`, `01_vs_vs_E.png`, `01_vs_vs_width.png`
- `02_<mode>_seed<n>_z_profile.png/.csv`, `_xy_map.png`, `_rough_top_2d/3d.png`, `.xyz`, `.vtu`
- `03_seed_ensemble.csv`, `03_vs_distribution.png`, `03_field.txt`
