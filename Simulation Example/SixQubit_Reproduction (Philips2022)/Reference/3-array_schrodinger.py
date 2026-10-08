__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [EXTENSION] Six-qubit Si/SiGe array reproduction — Phase 1 (Schrödinger 검증)
#
# [목적] 2-array_poisson.py에서 확인한 이중 우물 포텐셜 위에서 Schrödinger를 풀어,
#   가장 낮은 N_dots개 고유상태가 실제로 "각 dot에 하나씩 국소화"되는지 확인한다
#   (완전한 many-body (1,1,...,1) 점유 확정은 아니지만, 단일입자 궤도가 dot당 하나씩
#   분리되어 있다는 것은 그 전제조건임 — M03 tunnel_coupling_1.py와 동일한 SubDevice
#   패턴 사용).
#
# [BASELINE] Poisson 설정은 2-array_poisson.py와 100% 동일(바이어스/리전/geo_file) —
#   이 스크립트는 재실행 시 Poisson부터 다시 풀고(아직 저장된 phi 재사용 파이프라인
#   미구축), 그 위에 Schrödinger만 추가한다.

import pathlib
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams

N_DOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2
NUM_STATES = max(2 * N_DOTS, 4)  # [NUMERICAL] dot당 1개씩 + 여유분

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / f"array_n{N_DOTS}.msh"
path_out = script_dir / "output"
path_out.mkdir(parents=True, exist_ok=True)
path_geo = pathlib.Path(r"C:\temp\sixqubit") / f"array_n{N_DOTS}.xao"

scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

dvc = Device(mesh, conf_carriers="e")
dvc.set_temperature(0.1)
dvc.new_region("oxide", mt.SiO2)
dvc.new_region("channel", mt.Si)
dvc.new_region("source", mt.Si, ndoping=1e20 * 1e6)
dvc.new_region("drain", mt.Si, ndoping=1e20 * 1e6)
for i in range(1, N_DOTS + 1):
    dvc.new_region(f"oxide.QD{i}", mt.SiO2)
    dvc.new_region(f"channel.QD{i}", mt.Si)

barrier_bias = 0.3
plunger_bias = 0.6
back_gate_bias = -0.5
Ew = mt.Si.Eg / 2 + mt.Si.chi
for i in range(N_DOTS + 1):
    dvc.new_gate_bnd(f"barrier_gate_{i}_bnd", barrier_bias, Ew)
for i in range(1, N_DOTS + 1):
    dvc.new_gate_bnd(f"plunger_gate_{i}_bnd", plunger_bias, Ew)
dvc.new_ohmic_bnd("source_bnd")
dvc.new_ohmic_bnd("drain_bnd")
dvc.new_frozen_bnd("back_gate_bnd", back_gate_bias, mt.Si, 1e15 * 1e6, "n", 46 * 1e-3 * ct.e)

dot_region_list = [f"oxide.QD{i}" for i in range(1, N_DOTS + 1)] + \
                   [f"channel.QD{i}" for i in range(1, N_DOTS + 1)]
dvc.set_dot_region(dot_region_list)

params = PoissonSolverParams()
params.tol = 1e-3
params.initial_ref_factor = 0.1
params.final_ref_factor = 0.75
params.min_nodes = 20000
params.max_nodes = 2e5
params.maxiter_adapt = 30
params.maxiter = 200
params.refined_region = dot_region_list
params.h_refined = 1.0

slv = PoissonSolver(dvc, solver_params=params, geo_file=str(path_geo))
slv.solve()
dvc.set_V_from_phi()
print(f"[Phase 1] Poisson 완료, 최종 노드 수: {dvc.mesh.node_number}")

# --- Schrödinger: dot 영역만 SubMesh로 잘라서 풀기 (M03 tunnel_coupling_1.py 패턴) ---
submesh = SubMesh(dvc.mesh, dot_region_list)
subdvc = SubDevice(dvc, submesh)

params_schrod = SchrodingerSolverParams()
params_schrod.tol = 1e-6
params_schrod.num_states = NUM_STATES

schrod_slv = SchrodingerSolver(subdvc, solver_params=params_schrod)
schrod_slv.solve()

energies_eV = np.array(subdvc.energies) / ct.e
print(f"\n[Schrödinger] 가장 낮은 {NUM_STATES}개 고유에너지 (eV):")
for i, E in enumerate(energies_eV):
    print(f"  state {i}: E = {E:.6f} eV")

# --- 각 상태의 |psi|^2 y-중심(centroid)으로 국소화 위치 확인 ---
x_sub, y_sub, z_sub = subdvc.mesh.glob_nodes.T
# [PHYSICS] plunger 중심 기대 위치 (1-array_builder.py와 동일 계산)
plunger_w, barrier_w, pitch = 45, 25, 10
QD_w = plunger_w + barrier_w + 2 * pitch
channel_len = N_DOTS * QD_w + barrier_w
y_cursor = -channel_len / 2
expected_plunger_y = []
for i in range(N_DOTS):
    y_cursor += barrier_w + pitch
    expected_plunger_y.append((y_cursor + plunger_w / 2) * 1e-9)
    y_cursor += plunger_w + pitch

centroids_y = []
for i in range(NUM_STATES):
    psi = subdvc.eigenfunctions[:, i]
    density = np.abs(psi) ** 2
    if density.ndim > 1:
        density = np.sum(density, axis=-1)
    density = density / np.sum(density)
    y_centroid = np.sum(density * y_sub)
    centroids_y.append(y_centroid)
    nearest_dot = int(np.argmin(np.abs(np.array(expected_plunger_y) - y_centroid))) + 1
    print(f"  state {i}: <y> = {y_centroid/1e-9:.2f} nm "
          f"(가장 가까운 plunger: QD{nearest_dot} at {expected_plunger_y[nearest_dot-1]/1e-9:.1f}nm)")

# --- 판정: 가장 낮은 N_dots개 상태가 서로 다른 dot에 하나씩 국소화됐는가? ---
lowest_centroids = centroids_y[:N_DOTS]
assigned_dots = [int(np.argmin(np.abs(np.array(expected_plunger_y) - c))) for c in lowest_centroids]
unique_dots = len(set(assigned_dots))
print(f"\n[판정] 가장 낮은 {N_DOTS}개 상태가 할당된 서로 다른 dot 개수: {unique_dots} "
      f"(목표: {N_DOTS})")
if unique_dots == N_DOTS:
    print(f"[결론] 가장 낮은 {N_DOTS}개 단일입자 궤도가 각 dot에 하나씩 분리 국소화됨 — "
          f"(1,1,...,1) 점유의 단일입자 전제조건 충족. (many-body 전자수 확정은 아직 아님)")
else:
    print(f"[결론] 낮은 상태들이 같은 dot에 쏠려 있음(bonding/antibonding 분리가 큰 "
          f"tunnel coupling 때문일 수 있음) — barrier 바이어스를 낮춰 재시도 필요.")

# --- Plot: 에너지 준위 + 상태별 y-중심 ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.hlines(energies_eV, 0, 1, color="tab:blue")
for i, E in enumerate(energies_eV):
    ax1.text(1.05, E, f"state {i}", va="center", fontsize=8)
ax1.set_xlim(0, 1.6); ax1.set_xticks([])
ax1.set_ylabel("Energy (eV)"); ax1.set_title("Lowest eigenstates")

ax2.scatter(range(NUM_STATES), np.array(centroids_y) / 1e-9, c="tab:red", zorder=5)
for ypos in expected_plunger_y:
    ax2.axhline(ypos / 1e-9, color="gray", ls="--", lw=1)
ax2.set_xlabel("state index"); ax2.set_ylabel("<y> centroid (nm)")
ax2.set_title("State localization vs. expected plunger positions (gray dashed)")
fig.tight_layout()
fig.savefig(path_out / f"array_n{N_DOTS}_schrodinger_states.png", dpi=150)
plt.close(fig)

metadata = {
    "N_dots": N_DOTS, "num_states": NUM_STATES, "energies_eV": energies_eV.tolist(),
    "centroids_y_nm": [c / 1e-9 for c in centroids_y],
    "expected_plunger_y_nm": [y / 1e-9 for y in expected_plunger_y],
    "unique_dots_in_lowest_N": unique_dots,
}
with open(path_out / f"array_n{N_DOTS}_schrodinger_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"\nExports: {path_out}")
