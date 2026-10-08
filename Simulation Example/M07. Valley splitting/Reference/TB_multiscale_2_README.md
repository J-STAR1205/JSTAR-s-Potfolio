# Reference/TB_multiscale_2 — SiGe Part 2 (tight-binding valley splitting) 정답지

공식 튜토리얼 원본을 그대로 실행해 정답지로 둔다. 원본 코드는 수정하지 않는다.

## 준비
1. QTCAD 설치본의 `examples/tutorials/multiscale_2.py`를 이 폴더로 복사한다.
2. Part 1을 실행해 생긴 파일을 이 폴더 아래 같은 상대 경로로 복사한다.
   - `meshes/multiscale.msh`
   - `output/multiscale_phi.hdf5`
   - `output/multiscale_der.hdf5`
   (원본 스크립트가 자기 폴더 기준 `meshes/`, `output/`을 읽는다.)
3. 실행: `python multiscale_2.py > multiscale_2_run.log 2>&1`

## 정답지 값 (공식 문서 기준, rng_seed = 104, 거칠기 시드 0·1)

| 항목 | 값 |
|---|---|
| 전체 원자 수 / 양자점 상자 원자 수 | 896,337 / 87,189 |
| TB 해밀토니안 크기 | 871,890 × 871,890 |
| 고유에너지 (eV) | 0.01042060, 0.01190612, 0.03169008, 0.03306697, … |
| Valley splitting | 1.4855 meV |
| Valley phase | 2.2296 rad |
| 바닥 상태: Si / Ge 위 확률 | 96.6 % / 3.4 % |
| 바닥 상태 p_z 성분 | 40 % |
| 바닥 상태 크기 (x, y, z) | 111, 115, 30 Å |
| SOC 행렬원소 | (-1.036 + 1.102j)×10⁻⁵ eV |
| D_BG 행렬원소 | 6.12×10⁻⁶ eV/V |

## 분석
실행 후 `Reference/analysis.ipynb`의 **Part B** 셀(8~12절)을 실행한다.
`output/multiscale.xyz`를 읽어 원자 배치도·격자 구조·조성 분포를 그리고,
그림은 `output/analysis_tb_*.png`로 저장된다.

| 절 | 내용 |
|---|---|
| 8 | 원자층별 Ge 조성 프로파일, 실제 계면 위치와 10–90 % 계면 폭 |
| 9 | 측면 단면 배치도(전체·확대), 원자층 평면도 4장, 계면 부근 Ge 3D |
| 10 | 상부 계면 격자 확대도(결합선), 결합 길이 분포(Si–Si/Si–Ge/Ge–Ge), 원자층 간격(변형) |
| 11 | 스페이서 1 nm slab의 국소 Ge 지도와 이항분포 비교 |
| 12 | 양자점 상자 원자 구성, `multiscale_params.txt`(VS, C, D) 정답지 대조 |

## 실행 후 기록할 것 (`../../notes.md`)
- 위 값과의 일치 여부 (특히 VS)
- 단계별 시간: 원자 구조 생성, Keating 이완, TB 풀이
- 최대 메모리 사용량 (작업 관리자로 확인)
- `output/multiscale.xyz`를 VESTA로 열어 본 결과 (우물, 거친 계면, Ge 분포, SiGe 층 부풂)
