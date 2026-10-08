__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [EXTENSION] Six-qubit Si/SiGe array reproduction — Phase 1 (Poisson 검증)
#
# 레퍼런스: Philips et al., Nature 609, 919 (2022).
# 분류: C. 물리 확장 — 새 geometry, 공식 튜토리얼 아님.
#
# [목적] 1-array_builder.py 2 로 만든 N_dots=2 mesh에서, barrier/plunger 게이트
#   바이어스로 실제 "두 개의 독립된 quantum dot"(이중 우물 포텐셜)이 형성되는지
#   비선형 Poisson으로 확인한다. 이것이 Phase 1의 핵심 검증 목표이며, 이게
#   통과해야 Phase 2(N_dots=6)로 확장할 근거가 생긴다.
#
# [BASELINE] 리전/바이어스 구조는 M03/M17의 get_double_dot_fdsoi() 패턴을 그대로
#   따른다(oxide/channel 본체 + oxide.QD{i}/channel.QD{i} dot 서브리전, barrier_gate_
#   {0..N}_bnd/plunger_gate_{1..N}_bnd 게이트, source/drain ohmic, back_gate frozen).
#   N_dots=2 전용이 아니라 N_dots 파라미터로 일반화해 Phase 2에서 코드 변경 없이
#   재사용 가능하게 작성.
#
# [NUMERICAL] 아직 실제 Si/SiGe 이종구조(cap/spacer/QW) 스택이 아니라 M03/M17과
#   동일한 단층 Si channel 근사를 쓴다 (설계도안 Phase 3에서 실제 스택으로 교체 예정).

import pathlib
import sys
import numpy as np
import matplotlib.pyplot as plt
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.mesh3d import Mesh
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device import analysis as an

N_DOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / f"array_n{N_DOTS}.msh"
path_out = script_dir / "output"
path_out.mkdir(parents=True, exist_ok=True)

# [WARNING] 한글 경로에서 gmsh XAO reader 실패 (M03/M06/M08 등과 동일 버그) —
# adaptive Poisson(h_refined)이 내부적으로 geo_file을 다시 읽어야 하므로 ASCII
# 경로의 .xao(1-array_builder.py가 이미 C:\temp\sixqubit\에 저장해둔 것)를 사용.
path_geo = pathlib.Path(r"C:\temp\sixqubit") / f"array_n{N_DOTS}.xao"

scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# --- Device: 리전 ---------------------------------------------------------
dvc = Device(mesh, conf_carriers="e")
dvc.set_temperature(0.1)

dvc.new_region("oxide", mt.SiO2)
dvc.new_region("channel", mt.Si)
dvc.new_region("source", mt.Si, ndoping=1e20 * 1e6)
dvc.new_region("drain", mt.Si, ndoping=1e20 * 1e6)
for i in range(1, N_DOTS + 1):
    dvc.new_region(f"oxide.QD{i}", mt.SiO2)
    dvc.new_region(f"channel.QD{i}", mt.Si)

# --- Device: 바이어스 -------------------------------------------------------
# [PHYSICS] barrier를 낮게(채널 핀치오프), plunger를 높게(전자 축적) 걸어 두
#   개의 독립된 우물을 만든다 -- M03 tunnel_coupling_1.py의 바이어스 스케일
#   (barrier~0.5V, plunger~0.59V)을 그대로 참고.
barrier_bias = 0.3
plunger_bias = 0.6
back_gate_bias = -0.5

Ew = mt.Si.Eg / 2 + mt.Si.chi  # midgap
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

# --- Poisson 솔브 (adaptive, dot 영역만 국소 정제) ---------------------------
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
print(f"[Phase 1, N_dots={N_DOTS}] Poisson 수렴 완료. 최종 mesh 노드 수: {dvc.mesh.node_number}")

# --- 채널 축(y) 전도대 바닥 linecut: dot 형성(이중/다중 우물) 확인 -------------
x, y, z = dvc.mesh.glob_nodes.T
ymin, ymax = float(np.min(y)), float(np.max(y))
distance, Ec = an.linecut(dvc.mesh, dvc.cond_band_edge(), (0, ymin, -1e-9), (0, ymax, -1e-9))
Ec_eV = Ec / ct.e

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(distance / 1e-9, Ec_eV, "-b")
ax.set_xlabel("y (nm, along channel)")
ax.set_ylabel("$E_C$ (eV)")
ax.set_title(f"Conduction band edge along channel axis (N_dots={N_DOTS}, "
             f"barrier={barrier_bias}V, plunger={plunger_bias}V)")
ax.grid(True)
fig.tight_layout()
fig.savefig(path_out / f"array_n{N_DOTS}_bandedge_linecut.png", dpi=150)
plt.close(fig)

np.savez(path_out / f"array_n{N_DOTS}_bandedge_linecut.npz", distance_m=distance, Ec_eV=Ec_eV)

# --- 국소 극소점(우물) 탐지: plunger 개수만큼의 독립된 dot이 형성됐는지 확인 ----
from scipy.signal import argrelextrema
order = np.argsort(distance)
d_sorted, Ec_sorted = distance[order], Ec_eV[order]
minima_idx = argrelextrema(Ec_sorted, np.less_equal, order=5)[0]
# 평평한 구간에서 중복 검출 방지: 인접한 극소점은 5nm 이내면 하나로 묶음
minima_y = d_sorted[minima_idx] / 1e-9
minima_E = Ec_sorted[minima_idx]
merged_y, merged_E = [], []
for yv, Ev in zip(minima_y, minima_E):
    if merged_y and abs(yv - merged_y[-1]) < 5.0:
        if Ev < merged_E[-1]:
            merged_y[-1], merged_E[-1] = yv, Ev
        continue
    merged_y.append(yv)
    merged_E.append(Ev)

print(f"\n검출된 국소 극소(우물) 개수: {len(merged_y)} (목표: {N_DOTS}개)")
for i, (yv, Ev) in enumerate(zip(merged_y, merged_E)):
    print(f"  dot candidate {i+1}: y={yv:.2f}nm, E_C={Ev:.4f}eV")

if len(merged_y) == N_DOTS:
    print(f"[판정] N_dots={N_DOTS}개의 독립된 우물이 형성됨 — Phase 1 통과 가능성 있음 "
          f"(완전한 (1,1,...,1) 전하점유 확인은 Schrödinger/ManyBody 솔브가 필요, 아직 미수행).")
else:
    print(f"[판정] 검출된 우물 개수({len(merged_y)})가 목표(N_dots={N_DOTS})와 다름 — "
          f"바이어스 재조정 또는 치수 재검토 필요.")

import json
metadata = {
    "N_dots": N_DOTS, "barrier_bias_V": barrier_bias, "plunger_bias_V": plunger_bias,
    "back_gate_bias_V": back_gate_bias, "final_mesh_nodes": dvc.mesh.node_number,
    "n_wells_detected": len(merged_y), "well_positions_nm": merged_y,
}
with open(path_out / f"array_n{N_DOTS}_phase1_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"\nExports: {path_out}")
