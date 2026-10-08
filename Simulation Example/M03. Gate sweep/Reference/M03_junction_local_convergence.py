__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M03 2차 라운드 Step 5 (지침 25절) — 접합부 국소 3단계 mesh 수렴성
#
# 분류: B. Parameter sensitivity / numerical convergence study.
#   tunnel_coupling_1.py/2.py, dqdfdsoi.geo/.xao, M03_Gate_sweep_Analysis.ipynb는
#   전혀 수정하지 않는다. M03_generate_junction_local_mesh.py가 gmsh Field(Box)로
#   만든 새 .msh 3개(h_local=2.0/1.0/0.5nm, 접합 y=-45/+45nm 주변 ±10nm만 국소
#   정제, 원본 .xao는 읽기만 함)를 입력으로 사용한다.
#
# [BASELINE] 바이어스/detuning/Poisson tol/dot_region 적응 정제는 1차 라운드
#   baseline(Case 0)과 100% 동일하게 유지 — 이번에 바뀌는 변수는 "접합부 국소
#   시작 mesh 조밀도"뿐이다 (지침 10절: 한 번에 하나의 변수만 바꾼다).
#
# [NUMERICAL] 시작 mesh의 노드 수가 이미 커서(h_local=0.5 -> 267,832) QTCAD의
#   dot_region 적응 정제가 더해지면 총 노드수가 더 늘어날 수 있다 — max_nodes를
#   시작 mesh 크기보다 넉넉히 높게 잡는다 (과거 1차 라운드에서 max_nodes 기본값
#   때문에 암묵적으로 잘릴 위험을 본 적 있어 명시적으로 체크).

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from qtcad.device import constants as ct
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

script_dir = pathlib.Path(__file__).parent.resolve()
path_out = script_dir / "output"
export_dir = path_out / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

ascii_staging = pathlib.Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation"
    r"\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo"
)
path_geo = ascii_staging / "dqdfdsoi.xao"

scaling = 1e-9

# [BASELINE] tunnel_coupling_1.py / M03_mesh_convergence.py와 100% 동일
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias_high = 0.57
plunger_gate_2_bias = 0.59
barrier_gate_3_bias = 0.5

detuning_file = path_out / "detuning_and_coupling.txt"
detuning_optimum = float(np.loadtxt(detuning_file)[0])
print(f"[캐시된 최적 detuning 로드] {detuning_optimum:.6e} V")

dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]
Z_LINECUT_NM = -1.0
junction_left_nm, junction_right_nm = -45.0, 45.0


def run_case(tag, msh_path, max_nodes):
    mesh = Mesh(scaling, str(msh_path))
    n_start = mesh.node_number
    dvc = get_double_dot_fdsoi(
        mesh, back_gate_bias, barrier_gate_1_bias, plunger_gate_1_bias,
        barrier_gate_2_bias_high, plunger_gate_2_bias, barrier_gate_3_bias,
    )
    dvc.set_applied_potential("plunger_gate_2_bnd", plunger_gate_2_bias + detuning_optimum)

    params_poisson = PoissonSolverParams()
    params_poisson.tol = 1e-3  # [BASELINE] 동일
    params_poisson.initial_ref_factor = 0.1
    params_poisson.final_ref_factor = 0.75
    params_poisson.min_nodes = 50000
    params_poisson.max_nodes = max_nodes
    params_poisson.maxiter_adapt = 30
    params_poisson.maxiter = 200
    params_poisson.refined_region = dot_region_list  # [BASELINE] baseline과 동일(dot만)
    params_poisson.h_refined = 0.8  # [BASELINE] baseline과 동일
    params_poisson.refined_mesh_filename = str(ascii_staging / f"refined_{tag}.msh")

    t0 = time.time()
    slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=str(path_geo))
    slv.solve()
    elapsed = time.time() - t0
    n_final = dvc.mesh.node_number
    print(f"[{tag}] start_nodes={n_start}, final_nodes={n_final}, time={elapsed:.1f}s")

    x, y, z = dvc.mesh.glob_nodes.T
    ymin, ymax = np.min(y), np.max(y)
    distance, Ec = an.linecut(dvc.mesh, dvc.cond_band_edge(), (0, ymin, Z_LINECUT_NM * scaling),
                               (0, ymax, Z_LINECUT_NM * scaling))
    return distance, Ec / ct.e, elapsed, n_final, n_start


def quantify_wiggle(distance_m, Ec_eV, junction_y_nm, window_nm=5.0):
    """[BASELINE] M03_mesh_convergence.py(1차 라운드)와 동일한 방법
    (scipy.signal.argrelextrema 기반 국소 극값 탐지)을 그대로 재사용."""
    y_nm = distance_m / 1e-9 - 65.0
    mask = (y_nm > junction_y_nm - window_nm) & (y_nm < junction_y_nm + window_nm)
    if np.sum(mask) < 7:
        return np.nan, np.nan
    sub_y, sub_E = y_nm[mask], Ec_eV[mask]
    order = np.argsort(sub_y)
    sub_y, sub_E = sub_y[order], sub_E[order]
    i_min = argrelextrema(sub_E, np.less, order=3)[0]
    i_max = argrelextrema(sub_E, np.greater, order=3)[0]
    if len(i_min) == 0 or len(i_max) == 0:
        return np.nan, np.nan
    best = None
    for im in i_min:
        for ix in i_max:
            dist = abs(sub_y[im] - sub_y[ix])
            if best is None or dist < best[0]:
                best = (dist, abs(sub_E[ix] - sub_E[im]), sub_y[ix])
    _, amp, y_max = best
    return amp * 1e3, y_max


cases_spec = [
    ("junction_local_h2.0", ascii_staging / "dqdfdsoi_junction_local_h2.0.msh", 2.5e5),
    ("junction_local_h1.0", ascii_staging / "dqdfdsoi_junction_local_h1.0.msh", 3e5),
    ("junction_local_h0.5", ascii_staging / "dqdfdsoi_junction_local_h0.5.msh", 5e5),
]

cases_data = {}
for tag, msh_path, max_nodes in cases_spec:
    cases_data[tag] = run_case(tag, msh_path, max_nodes)
    distance, Ec_eV, elapsed, n_final, n_start = cases_data[tag]
    np.savez(export_dir / f"m03_raw_linecut_{tag}.npz", distance_m=distance, Ec_eV=Ec_eV)

rows = []
for tag, (distance, Ec_eV, elapsed, n_final, n_start) in cases_data.items():
    left_amp, left_pos = quantify_wiggle(distance, Ec_eV, junction_left_nm)
    right_amp, right_pos = quantify_wiggle(distance, Ec_eV, junction_right_nm)
    rows.append({
        "case": tag, "start_nodes": n_start, "final_nodes": n_final, "wall_time_s": elapsed,
        "left_wiggle_amplitude_meV": left_amp, "left_wiggle_position_nm": left_pos,
        "right_wiggle_amplitude_meV": right_amp, "right_wiggle_position_nm": right_pos,
    })

conv_df = pd.DataFrame(rows)
print("\n" + "=" * 100)
print("JUNCTION-LOCAL CONVERGENCE TABLE")
print("=" * 100)
print(conv_df.to_string(index=False))
conv_df.to_csv(export_dir / "m03_junction_local_convergence_results.csv", index=False)

# --- Plot ---
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
colors = ["tab:green", "tab:orange", "tab:purple"]
for (tag, (distance, Ec_eV, _, _, _)), c in zip(cases_data.items(), colors):
    y_nm = distance / 1e-9 - 65.0
    axs[0].plot(y_nm, Ec_eV, color=c, label=tag, lw=1.2)
    mask_l = (y_nm > junction_left_nm - 6) & (y_nm < junction_left_nm + 6)
    axs[1].plot(y_nm[mask_l], Ec_eV[mask_l], color=c, label=tag, lw=1.5, marker=".")
axs[0].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[0].axvline(junction_right_nm, color="k", ls=":", lw=1)
axs[0].set_xlabel("y (nm)"); axs[0].set_ylabel("$E_C$ (eV)")
axs[0].set_title("Junction-local refinement: full linecut"); axs[0].legend(fontsize=7)
axs[1].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[1].set_xlabel("y (nm)"); axs[1].set_ylabel("$E_C$ (eV)")
axs[1].set_title("Zoom: source/channel junction (±6nm)"); axs[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(export_dir / "m03_junction_local_convergence.png", dpi=150)
plt.close(fig)

# --- 수렴 판정 (project practical criterion) ---
if len(conv_df) >= 2:
    last_two = conv_df.iloc[-2:]
    for side in ["left", "right"]:
        a0 = last_two[f"{side}_wiggle_amplitude_meV"].iloc[0]
        a1 = last_two[f"{side}_wiggle_amplitude_meV"].iloc[1]
        p0 = last_two[f"{side}_wiggle_position_nm"].iloc[0]
        p1 = last_two[f"{side}_wiggle_position_nm"].iloc[1]
        if np.isnan(a0) or np.isnan(a1):
            print(f"[{side}] 마지막 두 case 중 하나에서 wiggle 미검출 (NaN) — amplitude->0 수렴 가능성")
            continue
        rel = abs(a1 - a0) / abs(a0) if a0 else np.nan
        pos_shift = abs(p1 - p0)
        converged = (rel < 0.05) and (pos_shift < 0.2)
        print(f"[{side}] rel_change={rel*100:.1f}%, pos_shift={pos_shift:.3f}nm -> "
              f"{'수렴(artifact 쪽)' if converged else '미수렴/유지(physical 가능성)'}")

metadata = {
    "classification": "B. Parameter sensitivity / numerical convergence study",
    "method": "gmsh Mesh.Field(Box) local refinement at y=-45/+45nm +-10nm "
              "(x: +-25nm, z: -22 to +3nm), VOut=100nm, Thickness=5nm",
    "junction_left_nm": junction_left_nm, "junction_right_nm": junction_right_nm,
    "h_local_values_nm": [2.0, 1.0, 0.5],
    "convergence_criterion": "rel_amplitude_change<5% AND position_shift<0.2nm (project practical)",
}
with open(export_dir / "m03_junction_local_convergence_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
