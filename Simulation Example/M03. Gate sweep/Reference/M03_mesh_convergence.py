__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M03 Step (지침 15절) — Junction-local mesh convergence
#
# 분류: B. Parameter sensitivity / numerical convergence study.
#   tunnel_coupling_1.py/2.py, M03_Gate_sweep_Analysis.ipynb는 전혀 수정하지 않는다.
#
# [BASELINE] baseline(Case 0)은 이미 저장된 tunnel_coupling_final_potential.hdf5를
#   그대로 사용한다 (재계산 없음). 이 스크립트는 Case 1/2만 새로 계산해서, source/
#   channel/drain 영역을 추가로 적응 정제했을 때 source/channel, channel/drain
#   접합부의 potential wiggle amplitude/position이 어떻게 바뀌는지 비교한다.
#
# [NUMERICAL] 원본은 dot_region_list(oxide_dot 등 4개)만 refined_region으로 지정하고
#   h_refined=0.8을 쓴다. source/channel/drain 자체는 적응 정제 대상이 아니었다
#   (M03 조사 보고서 5절). 여기서는 이 세 영역을 refined_region에 추가해 h_refined를
#   0.5(Case 1), 0.3(Case 2)로 낮춰 접합부 근처 메시를 강제로 더 조밀하게 만든다.
#   "source/channel 접합 ±5nm만" 별도로 떼어 정제하는 전용 physical group은 mesh에
#   없으므로, 이 영역들을 포함하는 전체 channel/source/drain 볼륨을 정제 대상으로
#   쓴다(지침이 요구한 이상적인 "접합부만" 국소 정제보다는 범위가 넓지만, 동일한
#   정성적 질문 — "접합부 메시를 조밀화하면 wiggle이 줄어드는가" — 에 답할 수 있다).
#
# [WARNING] 바이어스/온도/Poisson tolerance는 baseline과 100% 동일하게 유지한다.
#   detuning 최적화(minimize_scalar)는 재실행하지 않고, 이미 저장된 최적값
#   (detuning_and_coupling.txt)을 그대로 사용한다 — mesh sensitivity와 solver
#   sensitivity를 동시에 바꾸지 않기 위함 (지침 10절).

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from qtcad.device import constants as ct
from qtcad.device import io
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"
export_dir = path_out / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

ascii_staging = pathlib.Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation"
    r"\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo"
)
path_geo = ascii_staging / "dqdfdsoi.xao"

scaling = 1e-9

# [BASELINE] 게이트 바이어스 — tunnel_coupling_1.py와 100% 동일 (최종/최적 detuning 지점)
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias_high = 0.57
plunger_gate_2_bias = 0.59
barrier_gate_3_bias = 0.5

detuning_file = path_out / "detuning_and_coupling.txt"
detuning_optimum = float(np.loadtxt(detuning_file)[0])
print(f"[캐시된 최적 detuning 로드] {detuning_optimum:.6e} V (재최적화 안 함)")

dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]

Z_LINECUT_NM = -1.0  # [BASELINE] tunnel_coupling_1.py L86과 동일


def run_case(tag, extra_refined_regions, h_refined, mesh_out_tag):
    mesh = Mesh(scaling, str(path_mesh / "dqdfdsoi.msh"))
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
    params_poisson.max_nodes = 1.5e5  # [NUMERICAL] 추가 영역 정제로 더 필요할 수 있어 상향
    params_poisson.maxiter_adapt = 30
    params_poisson.maxiter = 200
    params_poisson.refined_region = dot_region_list + extra_refined_regions
    params_poisson.h_refined = h_refined
    params_poisson.refined_mesh_filename = str(ascii_staging / f"refined_dqdfdsoi_{mesh_out_tag}.msh")

    t0 = time.time()
    slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=str(path_geo))
    slv.solve()
    elapsed = time.time() - t0
    print(f"[{tag}] Poisson solve: {elapsed:.1f} s, nodes={dvc.mesh.node_number}")

    x, y, z = dvc.mesh.glob_nodes.T
    ymin, ymax = np.min(y), np.max(y)
    distance, Ec = an.linecut(dvc.mesh, dvc.cond_band_edge(), (0, ymin, Z_LINECUT_NM * scaling),
                               (0, ymax, Z_LINECUT_NM * scaling))
    return distance, Ec / ct.e, elapsed, dvc.mesh.node_number


def quantify_wiggle(distance_m, Ec_eV, junction_y_nm, window_nm=5.0):
    """junction_y_nm 주변 ±window_nm 안에서 scipy.signal.argrelextrema로 진짜
    국소 극값(local min/max)을 찾고, 서로 가장 가까이 붙어 있는 min-max 쌍을
    '작은 wiggle'로 식별해 진폭을 계산한다.

    [NUMERICAL] [수정 이력] 1차 시도(window 내 단순 min/max)는 접합부 전체
    band-bending(~120meV)을 그대로 잡았고, 2차 시도(선형 추세 제거)도 band-bending
    자체가 강하게 휘어 있어 충분히 분리하지 못했다(~45meV). M03 조사 보고서가 쓴
    원래 방법(scipy.signal.argrelextrema로 실제 국소 극값을 찾는 것)으로 교체해
    baseline 값(left 3.55mV, right 2.72mV)이 조사 보고서와 정확히 일치함을
    확인했다.
    """
    from scipy.signal import argrelextrema
    y_nm = distance_m / 1e-9 - 65.0  # ymin(-65nm) 기준 절대좌표로 변환
    mask = (y_nm > junction_y_nm - window_nm) & (y_nm < junction_y_nm + window_nm)
    if np.sum(mask) < 7:
        return np.nan, np.nan, (np.nan, np.nan)
    sub_y, sub_E = y_nm[mask], Ec_eV[mask]
    order = np.argsort(sub_y)
    sub_y, sub_E = sub_y[order], sub_E[order]
    i_min = argrelextrema(sub_E, np.less, order=3)[0]
    i_max = argrelextrema(sub_E, np.greater, order=3)[0]
    if len(i_min) == 0 or len(i_max) == 0:
        return np.nan, np.nan, (np.nan, np.nan)  # wiggle 검출 안 됨 (사라짐)
    best = None
    for im in i_min:
        for ix in i_max:
            dist = abs(sub_y[im] - sub_y[ix])
            if best is None or dist < best[0]:
                best = (dist, abs(sub_E[ix] - sub_E[im]), sub_y[im], sub_y[ix])
    _, amp, y_min, y_max = best
    amplitude_meV = amp * 1e3
    position_nm = y_max
    return amplitude_meV, position_nm, (y_min, y_max)


# --- Case 0: baseline (저장된 HDF5 재사용, 재계산 없음) ---
print("=" * 80)
print("Case 0 (baseline): 저장된 tunnel_coupling_final_potential.hdf5 재사용")
print("=" * 80)
mesh0 = Mesh(scaling, str(path_mesh / "dqdfdsoi.msh"))
# baseline은 refined_dqdfdsoi.msh 위에서 풀렸으므로 그 메시를 로드해야 phi와 매칭됨
mesh0_refined = Mesh(scaling, str(path_mesh / "refined_dqdfdsoi.msh"))
dvc0 = get_double_dot_fdsoi(
    mesh0_refined, back_gate_bias, barrier_gate_1_bias, plunger_gate_1_bias,
    barrier_gate_2_bias_high, plunger_gate_2_bias, barrier_gate_3_bias,
)
phi0 = io.load(str(path_out / "tunnel_coupling_final_potential.hdf5"), var_name="phi")
dvc0.set_potential(phi0)
x0, y0_arr, z0 = dvc0.mesh.glob_nodes.T
distance0, Ec0 = an.linecut(dvc0.mesh, dvc0.cond_band_edge(), (0, np.min(y0_arr), Z_LINECUT_NM * scaling),
                             (0, np.max(y0_arr), Z_LINECUT_NM * scaling))
Ec0_eV = Ec0 / ct.e
node_count0 = dvc0.mesh.node_number
print(f"[Case 0] baseline 메시 nodes={node_count0}")

cases_data = {
    "Case 0 (baseline, h_refined=0.8 dot-only)": (distance0, Ec0_eV, 0.0, node_count0),
}

# --- Case 1: 새로 계산 (channel/source/drain 추가 정제) ---
# [WARNING] Case 2(h_refined=0.3, channel/source/drain 전체)는 최초 시도에서
# 메시가 2,673,597 노드(16.4M 요소)까지 폭주해 "numpy._core._exceptions.
# ArrayMemoryError: Unable to allocate 620. MiB for an array with shape
# (162568510,)"로 크래시했다. channel/source/drain은 "접합부 ±5nm"보다 훨씬 큰
# 전체 볼륨(channel만 90nm 길이)이므로, h_refined=0.3을 그 전체에 적용하는 것은
# 애초에 감당 불가능한 설계였다 — 전용 junction-local physical group이 mesh에
# 없어 이렇게 넓은 영역을 쓸 수밖에 없었던 것이 근본 원인이다(M03 조사 보고서
# 5절 참고). 따라서 이번 비교는 baseline(h_refined=0.8, dot-only) vs
# Case 1(h_refined=0.5, +channel/source/drain) 2점 비교로 축소한다.
cases_data["Case 1 (h_refined=0.5, +channel/source/drain)"] = run_case(
    "Case 1", ["channel", "source", "drain"], 0.5, "case1")

# [NUMERICAL] raw linecut 저장 — 향후 wiggle 정의를 다시 조정해야 할 경우
# Poisson을 또 풀지 않고 재처리할 수 있도록.
for i, (case_name, (distance, Ec_eV, _, _)) in enumerate(cases_data.items()):
    np.savez(export_dir / f"m03_raw_linecut_case{i}.npz", distance_m=distance, Ec_eV=Ec_eV,
              case_name=case_name)

# --- 접합 위치 (M03 조사 보고서에서 확정된 exact 값) ---
junction_left_nm = -45.0   # source/channel boundary
junction_right_nm = 45.0   # channel/drain boundary

rows = []
for case_name, (distance, Ec_eV, elapsed, nodes) in cases_data.items():
    left_amp, left_pos, left_range = quantify_wiggle(distance, Ec_eV, junction_left_nm)
    right_amp, right_pos, right_range = quantify_wiggle(distance, Ec_eV, junction_right_nm)
    rows.append({
        "case": case_name, "mesh_nodes": nodes, "wall_time_s": elapsed,
        "left_wiggle_amplitude_meV": left_amp, "left_wiggle_position_nm": left_pos,
        "right_wiggle_amplitude_meV": right_amp, "right_wiggle_position_nm": right_pos,
    })

conv_df = pd.DataFrame(rows)
print("\n" + "=" * 80)
print("CONVERGENCE TABLE")
print("=" * 80)
print(conv_df.to_string(index=False))
conv_df.to_csv(export_dir / "m03_mesh_convergence_results.csv", index=False)

# --- Plot: 전체 linecut + 접합부 zoom ---
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
colors = ["tab:gray", "tab:blue", "tab:red"]
for (case_name, (distance, Ec_eV, _, _)), c in zip(cases_data.items(), colors):
    y_nm = distance / 1e-9 - 65.0
    axs[0].plot(y_nm, Ec_eV, color=c, label=case_name, lw=1.2)
    mask_l = (y_nm > junction_left_nm - 6) & (y_nm < junction_left_nm + 6)
    axs[1].plot(y_nm[mask_l], Ec_eV[mask_l], color=c, label=case_name, lw=1.5, marker=".")
axs[0].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[0].axvline(junction_right_nm, color="k", ls=":", lw=1)
axs[0].set_xlabel("y (nm)"); axs[0].set_ylabel("$E_C$ (eV)")
axs[0].set_title("Full linecut + junction markers"); axs[0].legend(fontsize=7)
axs[1].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[1].set_xlabel("y (nm)"); axs[1].set_ylabel("$E_C$ (eV)")
axs[1].set_title("Zoom: source/channel junction (±6nm)"); axs[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(export_dir / "m03_mesh_convergence.png", dpi=150)
plt.close(fig)

# --- 수렴 판정 (project practical criterion, 지침 31절) ---
if len(conv_df) >= 2:
    last_two = conv_df.iloc[-2:]
    left_change = abs(last_two["left_wiggle_amplitude_meV"].iloc[1]
                       - last_two["left_wiggle_amplitude_meV"].iloc[0])
    left_rel = left_change / abs(last_two["left_wiggle_amplitude_meV"].iloc[0]) if last_two["left_wiggle_amplitude_meV"].iloc[0] else np.nan
    right_change = abs(last_two["right_wiggle_amplitude_meV"].iloc[1]
                        - last_two["right_wiggle_amplitude_meV"].iloc[0])
    right_rel = right_change / abs(last_two["right_wiggle_amplitude_meV"].iloc[0]) if last_two["right_wiggle_amplitude_meV"].iloc[0] else np.nan
    pos_shift_l = abs(last_two["left_wiggle_position_nm"].iloc[1] - last_two["left_wiggle_position_nm"].iloc[0])
    pos_shift_r = abs(last_two["right_wiggle_position_nm"].iloc[1] - last_two["right_wiggle_position_nm"].iloc[0])
    print(f"\n마지막 두 case 간 변화: left_amp rel={left_rel*100:.1f}%, right_amp rel={right_rel*100:.1f}%, "
          f"pos_shift left={pos_shift_l:.3f}nm, right={pos_shift_r:.3f}nm")
    converged = (left_rel < 0.05 and right_rel < 0.05 and pos_shift_l < 0.2 and pos_shift_r < 0.2)
    print("판정 (project practical criterion, rel<5% AND pos_shift<0.2nm):",
          "수렴 가능성 있음(numerical artifact 쪽)" if converged else
          "아직 수렴하지 않음 또는 amplitude가 유지됨(physical feature 가능성)")

metadata = {
    "classification": "B. Parameter sensitivity / numerical convergence study",
    "junction_left_nm": junction_left_nm, "junction_right_nm": junction_right_nm,
    "detuning_used_V": detuning_optimum,
    "convergence_criterion": "rel_amplitude_change<5% AND position_shift<0.2nm (project practical, not QTCAD spec)",
}
with open(export_dir / "m03_mesh_convergence_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
