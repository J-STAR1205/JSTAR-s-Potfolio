__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M08 2차 라운드 Step 1 — 동일 입력 반복 실행 재현성 검증
#
# 분류: B. Parameter sensitivity / numerical reproducibility study.
#   MOS_EDSR.py, EDSR_noise.py는 전혀 수정하지 않는다 (이미 Step 1에서 수정된
#   h0/u 파이프라인은 그대로 보존).
#
# [BASELINE] 동일 geometry/mesh/gate bias/magnetic field/solver 설정으로
#   MOS_EDSR.py의 h0/u 계산 부분(mesh -> Poisson -> Bfield -> Schrodinger ->
#   Gate operator sweep -> transition_2_levels)을 3회 반복 실행한다.
#
# [NUMERICAL] MOS_EDSR.py 뒷부분(500,000-point qutip.mesolve 전체 Rabi 시뮬레이션,
#   plot_slices 호출 6개)은 h0/u 값과 무관하고 극도로 느리므로(~30분+) 포함하지
#   않는다 — 이 스크립트는 h0/u 계산까지만 재현한다.
#
# [WARNING] 결과가 run마다 달라지더라도 있는 그대로 보고한다. 원하는 결과를 얻기
#   위해 파라미터를 조정하지 않는다.

import json
import pathlib
import time

import numpy as np
import pandas as pd
from qtcad.device.mesh3d import Mesh
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.schrodinger import Solver as ssolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.poisson import Solver as psolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.operators import Gate
from qtcad.qubit import dynamics

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "MOS_EDSR_example.msh"
output_dir = script_dir / "output"
export_dir = output_dir / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

scaling = 1e-9
num_states = 6
phi_tg = 1.0
phi_sg = 0.0


def run_once(run_idx):
    """[BASELINE] MOS_EDSR.py 1-155행과 100% 동일한 계산 (plot_slices/plt.show() 제외)."""
    t0 = time.time()
    mesh = Mesh(scaling, str(path_mesh))
    d = Device(mesh, conf_carriers="e")
    d.statistics = "FD_approx"
    d.new_region("barrier", mt.SiO2)
    d.new_region("confined", mt.Si, pdoping=0, ndoping=0)
    d.new_region("undoped", mt.Si, pdoping=0, ndoping=0)
    d.new_region("doped", mt.Si, pdoping=5e18 * 1e6, ndoping=0)
    d.new_ohmic_bnd("bottom")
    Ew = mt.Si.chi + mt.Si.Eg / 2
    d.new_gate_bnd("top_gate_bnd", phi_tg, Ew)
    d.new_gate_bnd("side_gate_bnd", phi_sg, Ew)

    ps = psolver(d)
    ps.solve()
    d.set_V_from_phi()

    def Bfield(x, y, z):
        B0 = 0.6
        b = 0.3 * 1e6
        return np.array([B0, 0, b * x])

    d.set_Bfield(Bfield)

    params_schrod = SchrodingerSolverParams()
    params_schrod.num_states = num_states
    params_schrod.tol = 1e-12
    ss = ssolver(d, solver_params=params_schrod)
    ss.solve()

    eigenenergies_eV = np.array(d.energies) / ct.e
    mesh_nodes = mesh.node_number

    gate_biases = np.linspace(0, 2, num=6)
    gate = "side_gate_bnd"
    UU = np.zeros((len(gate_biases), num_states, num_states), dtype=complex)
    gate_params = PoissonSolverParams()
    gate_params.tol = 1e-3
    for i, V in enumerate(gate_biases):
        G = Gate(d, gate, V, params=gate_params)
        UU[i, :, :] = G.get_operator_matrix()
    delta_V = UU[-1]

    E = d.energies
    E = E - E[0]
    H0 = np.diagflat(E)

    dyn = dynamics.Dynamics()
    result, h0, u, omega0, omega_rabi = dyn.transition_2_levels(
        1, 0, H0, delta_V, plot=False, npts=4000
    )
    elapsed = time.time() - t0

    return {
        "run": run_idx, "mesh_nodes": mesh_nodes, "wall_time_s": elapsed,
        "f_qubit_Hz": omega0 / (2 * np.pi), "f_Rabi_Hz": omega_rabi / (2 * np.pi),
        "omega0_rad_s": omega0, "omega_rabi_rad_s": omega_rabi,
        "eigenenergies_eV": eigenenergies_eV, "h0": h0, "u": u,
    }


N_RUNS = 3
results = []
for i in range(1, N_RUNS + 1):
    print("=" * 80)
    print(f"RUN {i} / {N_RUNS}")
    print("=" * 80)
    r = run_once(i)
    results.append(r)
    print(f"[Run {i}] f_qubit={r['f_qubit_Hz']:.6e} Hz, f_Rabi={r['f_Rabi_Hz']:.6e} Hz, "
          f"nodes={r['mesh_nodes']}, time={r['wall_time_s']:.1f}s")

# --- 요약 표 ---
summary_rows = []
for r in results:
    summary_rows.append({
        "run": r["run"], "mesh_nodes": r["mesh_nodes"], "wall_time_s": r["wall_time_s"],
        "f_qubit_Hz": r["f_qubit_Hz"], "f_Rabi_Hz": r["f_Rabi_Hz"],
    })
summary_df = pd.DataFrame(summary_rows)

h0_ref = results[0]["h0"]
u_ref = results[0]["u"]
for i, r in enumerate(results):
    max_diff_h0 = np.max(np.abs(r["h0"] - h0_ref)) if i > 0 else 0.0
    max_diff_u = np.max(np.abs(r["u"] - u_ref)) if i > 0 else 0.0
    summary_df.loc[summary_df["run"] == r["run"], "max_abs_diff_h0_vs_run1"] = max_diff_h0
    summary_df.loc[summary_df["run"] == r["run"], "max_abs_diff_u_vs_run1"] = max_diff_u

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(summary_df.to_string(index=False))

f_rabi_vals = summary_df["f_Rabi_Hz"].values
f_qubit_vals = summary_df["f_qubit_Hz"].values
stats = {
    "f_Rabi_mean_Hz": float(np.mean(f_rabi_vals)),
    "f_Rabi_std_Hz": float(np.std(f_rabi_vals)),
    "f_Rabi_max_minus_min_Hz": float(np.max(f_rabi_vals) - np.min(f_rabi_vals)),
    "f_Rabi_relative_spread_pct": float((np.max(f_rabi_vals) - np.min(f_rabi_vals))
                                         / np.mean(f_rabi_vals) * 100),
    "f_qubit_mean_Hz": float(np.mean(f_qubit_vals)),
    "f_qubit_std_Hz": float(np.std(f_qubit_vals)),
    "f_qubit_relative_spread_pct": float((np.max(f_qubit_vals) - np.min(f_qubit_vals))
                                          / np.mean(f_qubit_vals) * 100),
    "mesh_nodes_all_runs": [int(r["mesh_nodes"]) for r in results],
}
print("\n" + "=" * 80)
print("REPRODUCIBILITY STATISTICS")
print("=" * 80)
for k, v in stats.items():
    print(f"  {k}: {v}")

summary_df.to_csv(export_dir / "m08_reproducibility_results.csv", index=False)
with open(export_dir / "m08_reproducibility_stats.json", "w", encoding="utf-8") as f:
    json.dump(stats, f, indent=2)

# 비교 참고: 과거 로그값(stale) vs 현재 npz(수정된 파이프라인)
print(f"\n참고 — 과거 로그(수정 전 1회성 실행): f_Rabi = 1.6899e6 Hz")
print(f"참고 — 이전 재실행(M08 1차 라운드 Step 1): f_Rabi = 1.7086e6 Hz")
print(f"이번 3회 반복: mean={stats['f_Rabi_mean_Hz']:.6e} Hz, "
      f"spread={stats['f_Rabi_relative_spread_pct']:.4f}%")

print(f"\nExports: {export_dir}")
