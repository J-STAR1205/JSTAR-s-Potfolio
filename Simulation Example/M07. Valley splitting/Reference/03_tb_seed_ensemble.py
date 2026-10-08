"""M07 Study 3단계 — 합금 시드별 TB valley splitting 분포 + MVEMT 비교용 전기장 추출

목적
------------------------------
같은 평균 조성·같은 게이트 전압이라도 Ge 원자가 놓인 위치(시드)가 다르면
VS가 달라진다. 시드만 바꿔 여러 번 풀어 VS의 평균과 산포를 구한다.
  - 첫 시드 104 + 거칠기 시드 (0, 1) + 튜토리얼 상자 = SiGe Part 2 재현
    → VS = 1.4855 meV가 나와야 한다 (검증).
  - 이후 시드: 합금 배치만 바뀐다 (VARY_ROUGHNESS=False일 때).

MVEMT와의 비교
------------------------------
1D MVEMT(01 스크립트)에 넣을 전기장을 FEM 퍼텐셜에서 직접 뽑는다:
양자점 중심(x, y ≈ 0) 근처 우물 안 노드에서 φ(z)를 직선 맞춤 → |E| = |dφ/dz|.
01 스크립트를 WELL_WIDTH_NM = 3(튜토리얼 우물), E_FIELDS_MVM = [이 값],
INTERFACE_WIDTHS_NM = [0] 으로 돌리면 같은 구조의 MVEMT VS가 나온다.
MVEMT는 거칠기·합금 무질서가 없으므로 TB 분포의 "평균 프로파일 기준값"이다.

필요 파일 (Part 1을 먼저 실행해 둘 것)
------------------------------
PART1_DIR/meshes/multiscale.msh
PART1_DIR/output/multiscale_phi.hdf5

실행: (qtcad) python 03_tb_seed_ensemble.py
출력: output/03_seed_ensemble.csv, output/03_vs_distribution.png, output/03_field.txt
예상 시간: 시드 1개당 튜토리얼 1회분(워크스테이션 기준 수 분, 원자 약 90만 개).
           처음에는 SEEDS를 2개로 시간부터 재고 늘릴 것.
"""

import pathlib
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device import io
from qtcad.atoms import SubAtoms
from qtcad.atoms.keating import Solver as KeatingSolver
from qtcad.atoms.keating import SolverParams as KeatingSolverParams
from qtcad.atoms.schrodinger import Solver as SchrodingerSolver
from qtcad.atoms.schrodinger import SolverParams as SchrodingerSolverParams

from atoms_builder import build_atoms, is_ge

# =============================================================================
# 바꿔볼 파라미터
# =============================================================================
# Part 1 튜토리얼 폴더 (meshes/, output/ 가 들어 있는 곳) — 실제 경로로 바꿀 것
PART1_DIR = pathlib.Path(
    r"C:\Users\norma\Desktop\연구자료\QTCAD Simulation\01. SiGe One QD_part.1\Reference"
)
SEEDS = [104, 1, 2, 3, 4]   # 첫 값 104 = 튜토리얼 재현
VARY_ROUGHNESS = False      # True면 거칠기 시드도 (seed, seed+1)로 바꿈
MODE = "rough"              # "rough" / "graded" (graded면 104도 튜토리얼과 다름)
INTERFACE_WIDTH_NM = 0.5    # graded 모드에서만 사용
X_GE = 0.30
WELL_TOP_NM, WELL_BOT_NM = -35.0, -38.0           # 튜토리얼
BOX_XY_NM, BOX_Z_NM = 10.0, (-50.0, -5.0)          # 튜토리얼 원자 상자
QD_XY_NM, QD_Z_NM = 8.0, (-40.0, -33.0)            # 튜토리얼 SubAtoms 상자
NUM_STATES = 4              # VS에는 2개면 충분 (튜토리얼은 10)
ENERGY_TARGET_EV = -0.1     # 튜토리얼 값
REF_VS_MEV = 1.4855196868505457                    # 튜토리얼 seed 104 결과

# =============================================================================
# 경로
# =============================================================================
script_dir = pathlib.Path(__file__).parent.resolve()
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)
path_csv = out_dir / "03_seed_ensemble.csv"
path_mesh = PART1_DIR / "meshes" / "multiscale.msh"
path_phi = PART1_DIR / "output" / "multiscale_phi.hdf5"
for p in (path_mesh, path_phi):
    assert p.exists(), f"Part 1 결과가 없습니다: {p}\nPART1_DIR을 확인하세요."

nm = 1e-9

# =============================================================================
# 1. FEM 소자 불러오기 + 양자점 위치 수직 전기장 추출
# =============================================================================
mesh = Mesh(1e-9, str(path_mesh))
d = Device(mesh, conf_carriers="e")
phi = np.asarray(io.load(str(path_phi), "phi")).ravel()
d.set_potential(phi)

nodes = np.asarray(mesh.glob_nodes, dtype=float) / nm   # nm
r_xy = np.hypot(nodes[:, 0], nodes[:, 1])
for r_max in (2.0, 5.0, 10.0):
    sel = (r_xy < r_max) & (nodes[:, 2] < WELL_TOP_NM) & (nodes[:, 2] > WELL_BOT_NM)
    if sel.sum() >= 5:
        break
slope, _ = np.polyfit(nodes[sel, 2] * nm, phi[sel], 1)   # V/m
E_mvm = abs(slope) / 1e6
msg = (f"양자점 중심(r < {r_max:g} nm, 노드 {sel.sum()}개) 우물 내 수직 전기장 "
       f"|E| = {E_mvm:.2f} MV/m\n"
       f"→ 01 스크립트: WELL_WIDTH_NM = {WELL_TOP_NM - WELL_BOT_NM:g}, "
       f"E_FIELDS_MVM = [{E_mvm:.2f}], INTERFACE_WIDTHS_NM = [0]")
print(msg)
(out_dir / "03_field.txt").write_text(msg + "\n", encoding="utf-8")

# =============================================================================
# 2. 시드 루프
# =============================================================================
x = np.array([-BOX_XY_NM, BOX_XY_NM]) * nm
z = np.array(BOX_Z_NM) * nm
x_qd = np.array([-QD_XY_NM, QD_XY_NM]) * nm
z_qd = np.array(QD_Z_NM) * nm

rows = []
for seed in SEEDS:
    rough_seeds = (seed, seed + 1) if VARY_ROUGHNESS else (0, 1)

    t0 = time.perf_counter()
    atoms, _ = build_atoms(x, x.copy(), z, WELL_TOP_NM * nm, WELL_BOT_NM * nm, X_GE,
                           mode=MODE, interface_width=INTERFACE_WIDTH_NM * nm,
                           seed=seed, rough_seeds=rough_seeds)
    t_build = time.perf_counter() - t0

    t0 = time.perf_counter()
    KeatingSolver(atoms=atoms, solver_params=KeatingSolverParams()).solve()
    t_keat = time.perf_counter() - t0

    atoms_qd = SubAtoms(parent=atoms, x=x_qd, y=x_qd.copy(), z=z_qd)
    atoms_qd.set_potential_from_device(d=d)
    sp = SchrodingerSolverParams()
    sp.num_states = NUM_STATES
    t0 = time.perf_counter()
    SchrodingerSolver(atoms=atoms_qd, solver_params=sp).solve(
        energy_target=ENERGY_TARGET_EV * ct.e)
    t_tb = time.perf_counter() - t0

    vs = (atoms_qd.energies[1] - atoms_qd.energies[0]) / ct.e * 1e3
    ge_qd = is_ge(atoms_qd.atom_species)
    prob = np.asarray(atoms_qd.eigenfunctions[0].get_prob_on_atoms()).ravel()
    p_ge = float(prob[ge_qd].sum() / prob.sum())   # 바닥 상태 확률 중 Ge 위에 있는 비율

    row = {"seed": seed, "rough_seeds": str(rough_seeds), "mode": MODE,
           "VS_meV": vs, "E0_meV": atoms_qd.energies[0] / ct.e * 1e3,
           "P_on_Ge": p_ge, "Ge_frac_QDbox": float(ge_qd.mean()),
           "n_atoms": len(atoms.atom_species), "n_atoms_QD": len(ge_qd),
           "t_build_s": t_build, "t_keating_s": t_keat, "t_TB_s": t_tb}
    rows.append(row)
    pd.DataFrame(rows).to_csv(path_csv, index=False)   # 시드마다 중간 저장
    print(f"seed {seed}: VS = {vs:.4f} meV, P(Ge) = {p_ge*100:.2f} %, "
          f"시간 build/Keating/TB = {t_build:.0f}/{t_keat:.0f}/{t_tb:.0f} s")
    if seed == 104 and MODE == "rough" and not VARY_ROUGHNESS:
        print(f"   [검증] 튜토리얼 VS = {REF_VS_MEV:.4f} meV, "
              f"차이 {abs(vs - REF_VS_MEV) / REF_VS_MEV * 100:.2f} %")

    del atoms, atoms_qd   # 메모리 회수

# =============================================================================
# 3. 분포 요약과 그림
# =============================================================================
df = pd.DataFrame(rows)
vs = df["VS_meV"].to_numpy()
print(f"\nVS 평균 {vs.mean():.4f} meV, 표준편차 {vs.std(ddof=1) if len(vs) > 1 else 0:.4f} meV, "
      f"범위 {vs.min():.4f}–{vs.max():.4f} meV (시드 {len(vs)}개)")

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
axs[0].plot(np.zeros_like(vs), vs * 1e3, "o", alpha=0.7)
for s_, v_ in zip(df["seed"], vs):
    axs[0].annotate(str(s_), (0, v_ * 1e3), xytext=(6, 0), textcoords="offset points", fontsize=8)
axs[0].axhline(vs.mean() * 1e3, color="k", ls="--", lw=1, label="mean")
axs[0].set_xticks([]); axs[0].set_ylabel("Valley splitting (µeV)")
axs[0].set_title(f"TB valley splitting by alloy seed ({MODE})"); axs[0].legend()
axs[1].plot(df["P_on_Ge"] * 100, vs * 1e3, "o")
axs[1].set_xlabel("Ground-state probability on Ge atoms (%)")
axs[1].set_ylabel("Valley splitting (µeV)")
axs[1].set_title("VS vs. wavefunction overlap with Ge"); axs[1].grid(True)
fig.tight_layout(); fig.savefig(out_dir / "03_vs_distribution.png", dpi=150)
plt.show()
