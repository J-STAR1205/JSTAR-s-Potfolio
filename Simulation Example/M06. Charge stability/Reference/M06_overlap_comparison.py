__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M06 Step 3 (지침 10절) — overlap=False vs overlap=True 비교
#
# 분류: C. Physics-oriented extension / validation study.
#   Official tutorial(double_dot_stability.py)의 baseline(overlap=False)은 보존하고,
#   별도로 overlap=True Coulomb tensor를 계산해 같은 90x90 gate-bias 스윕으로 CSD를
#   재생성한 뒤 비교한다.
#
# [BASELINE] Mesh, device, Poisson/Schrodinger 파라미터, LeverArmSolver, 게이트 바이어스
#   스윕 범위는 double_dot_stability.py와 100% 동일하게 유지한다. 바뀌는 것은
#   get_coulomb_matrix()의 overlap 인자와 그로부터 만들어지는 Junction/add_spectrum
#   결과뿐이다.
#
# [WARNING] overlap=True가 "더 현실적인 결과"를 자동으로 만든다고 가정하지 않는다
#   (지침 35절). Coulomb tensor 차원, 계산 시간, addition spectrum의 실제 차이를
#   수치로 확인한 뒤에만 해석한다.

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from progress.bar import ChargingBar as Bar
from qtcad.device import constants as ct
from qtcad.device import analysis as an
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from qtcad.device.leverarm_matrix import Solver as LeverArmSolver
from qtcad.device.leverarm_matrix import SolverParams as LeverArmSolverParams
from qtcad.transport.junction import Junction
from qtcad.transport.mastereq import add_spectrum
from helper.double_dot_fdsoi import get_double_dot_fdsoi

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes"
path_out = script_dir / "output"
export_dir = path_out / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh / "dqdfdsoi.msh"))

# [BASELINE] 게이트 바이어스 — double_dot_stability.py와 100% 동일
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.6
barrier_gate_2_bias = 0.5
plunger_gate_2_bias = 0.6 + 1e-3
barrier_gate_3_bias = 0.5

dvc = get_double_dot_fdsoi(
    mesh, back_gate_bias, barrier_gate_1_bias, plunger_gate_1_bias,
    barrier_gate_2_bias, plunger_gate_2_bias, barrier_gate_3_bias,
)

dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]

params_poisson = PoissonSolverParams()
params_poisson.tol = 1e-3
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 50000
params_poisson.max_nodes = 1e5
params_poisson.maxiter_adapt = 30
params_poisson.maxiter = 200
params_poisson.refined_region = dot_region_list
params_poisson.h_refined = 0.8
# [NUMERICAL] 동일한 Korean-path XAO 우회 (double_dot_stability.py와 동일 staging 경로 재사용)
params_poisson.size_map_filename = r"C:\temp\qtcad_dqdfdsoi\refined_dqdfdsoi_ovl.pos"
params_poisson.refined_mesh_filename = r"C:\temp\qtcad_dqdfdsoi\refined_dqdfdsoi_ovl.msh"

num_states = 4
params_schrod = SchrodingerSolverParams()
params_schrod.num_states = num_states

geo_file_ascii = r"C:\temp\qtcad_dqdfdsoi\dqdfdsoi.xao"
t0 = time.time()
poisson_slv = PoissonSolver(dvc, solver_params=params_poisson, geo_file=geo_file_ascii)
poisson_slv.solve()
t_poisson = time.time() - t0
print(f"[시간] Poisson solve: {t_poisson:.1f} s")

dvc.set_V_from_phi()
submesh = SubMesh(dvc.mesh, dot_region_list)
subdvc = SubDevice(dvc, submesh)

t0 = time.time()
schrod_solver = SchrodingerSolver(subdvc, solver_params=params_schrod)
schrod_solver.solve()
t_schrod = time.time() - t0
print(f"[시간] Schrodinger solve: {t_schrod:.1f} s")

energies = subdvc.energies
print("Energies (eV):", energies / ct.e)

# [VALIDATION] baseline(dqdfdsoi_energies.npy)과 교차검증
baseline_energies_file = path_out / "dqdfdsoi_energies.npy"
if baseline_energies_file.exists():
    baseline_energies = np.load(baseline_energies_file)
    rel_err_E = np.max(np.abs(energies - baseline_energies) / np.abs(baseline_energies))
    print(f"[consistency check] energies vs baseline(dqdfdsoi_energies.npy): "
          f"max rel.err = {rel_err_E:.2e}")
else:
    print("[주의] baseline energies 파일이 없어 교차검증을 건너뜁니다.")

bias_vector = np.array([plunger_gate_1_bias, plunger_gate_2_bias])
gate_labels_la = ["plunger_gate_1_bnd", "plunger_gate_2_bnd"]
lam_params = LeverArmSolverParams()
lam_params.pot_solver_params = params_poisson
lam_params.schrod_solver_params = params_schrod

t0 = time.time()
slv_la = LeverArmSolver(dvc, gate_labels_la, bias_vector, dot_region=dot_region_list,
                         solver_params=lam_params)
lever_arm_matrix = slv_la.solve(bias_increment=1e-3)
t_la = time.time() - t0
print(f"[시간] LeverArmSolver: {t_la:.1f} s")
print("Lever arm matrix:\n", lever_arm_matrix)

baseline_la_file = path_out / "lever_arm_matrix.npy"
if baseline_la_file.exists():
    baseline_la = np.load(baseline_la_file)
    rel_err_la = np.max(np.abs(lever_arm_matrix - baseline_la) / np.abs(baseline_la))
    print(f"[consistency check] lever_arm_matrix vs baseline: max rel.err = {rel_err_la:.2e}")

# --- Case A: overlap=False (baseline 재현) vs Case B: overlap=True ---
many_body_params = ManyBodySolverParams()
many_body_params.n_degen = 2
many_body_params.num_states = num_states
mb_slv = ManyBodySolver(subdvc, solver_params=many_body_params)

t0 = time.time()
coulomb_no_overlap = mb_slv.get_coulomb_matrix(overlap=False, verbose=True)
t_coulomb_false = time.time() - t0
print(f"[시간] Coulomb matrix (overlap=False): {t_coulomb_false:.1f} s, "
      f"shape={coulomb_no_overlap.shape}, size={coulomb_no_overlap.nbytes} bytes")

t0 = time.time()
coulomb_overlap = mb_slv.get_coulomb_matrix(overlap=True, verbose=True)
t_coulomb_true = time.time() - t0
print(f"[시간] Coulomb matrix (overlap=True): {t_coulomb_true:.1f} s, "
      f"shape={coulomb_overlap.shape}, size={coulomb_overlap.nbytes} bytes")

np.save(path_out / "coulomb_mat_no_overlap_rerun.npy", coulomb_no_overlap)
np.save(path_out / "coulomb_mat_overlap.npy", coulomb_overlap)

# baseline Coulomb(no overlap) 교차검증
baseline_coulomb_file = path_out / "coulomb_mat_no_overlap.npy"
if baseline_coulomb_file.exists():
    baseline_coulomb = np.load(baseline_coulomb_file)
    rel_err_c = np.max(np.abs(coulomb_no_overlap - baseline_coulomb)
                        / (np.abs(baseline_coulomb) + 1e-300))
    print(f"[consistency check] coulomb_no_overlap vs baseline: max rel.err = {rel_err_c:.2e}")

temperature_spec = 10  # [BASELINE] double_dot_stability.py와 동일
gate_labels = ["source_bnd", "plunger_gate_1_bnd", "plunger_gate_2_bnd", "drain_bnd"]

min_gate_1_bias, max_gate_1_bias = 30e-3, 140e-3
min_gate_2_bias, max_gate_2_bias = 30e-3, 140e-3
gate_1_biases = np.linspace(min_gate_1_bias, max_gate_1_bias, 90)
gate_2_biases = np.linspace(min_gate_2_bias, max_gate_2_bias, 90)


def get_add_spectrum(junc, gate_biases_full, temperature):
    junc.set_biases(gate_labels, gate_biases_full, verbose=False)
    return add_spectrum(junc, temperature=temperature)


def run_sweep(coulomb_mat, overlap_flag, tag):
    many_body_params_local = ManyBodySolverParams()
    many_body_params_local.n_degen = 2
    many_body_params_local.num_states = num_states
    many_body_params_local.energies = energies
    many_body_params_local.overlap = overlap_flag
    many_body_params_local.coulomb_mat = coulomb_mat
    many_body_params_local.alpha = lever_arm_matrix
    junc = Junction(many_body_solver_params=many_body_params_local,
                     temperature=temperature_spec, contact_labels=gate_labels)

    mat = np.zeros((len(gate_1_biases), len(gate_2_biases)))
    bar = Bar(f"CSD sweep ({tag})", max=len(gate_1_biases) * len(gate_2_biases))
    t0 = time.time()
    for i1, v1 in enumerate(gate_1_biases):
        for i2, v2 in enumerate(gate_2_biases):
            mat[i1, i2] = get_add_spectrum(junc, np.array([0, v1, v2, 0]), temperature_spec)
            bar.next()
    bar.finish()
    elapsed = time.time() - t0
    print(f"[시간] CSD sweep ({tag}): {elapsed:.1f} s")
    return mat, elapsed


mat_false, t_sweep_false = run_sweep(coulomb_no_overlap, False, "overlap=False")
np.savetxt(path_out / "addition_spectrum_overlapFalse_rerun.txt", mat_false)

mat_true, t_sweep_true = run_sweep(coulomb_overlap, True, "overlap=True")
np.savetxt(path_out / "addition_spectrum_overlapTrue.txt", mat_true)

# baseline CSD(addition_spectrum.txt, overlap=False) 교차검증
baseline_csd_file = path_out / "addition_spectrum.txt"
if baseline_csd_file.exists():
    baseline_csd = np.loadtxt(baseline_csd_file)
    rel_err_csd = np.max(np.abs(mat_false - baseline_csd) / (np.abs(baseline_csd) + 1e-300))
    print(f"[consistency check] CSD(overlap=False) vs baseline addition_spectrum.txt: "
          f"max rel.err = {rel_err_csd:.2e}")

# --- 비교 plot ---
def to_image(mat):
    return np.flip(np.transpose(mat), axis=0)

extent = [min_gate_1_bias + plunger_gate_1_bias, max_gate_1_bias + plunger_gate_1_bias,
          min_gate_2_bias + plunger_gate_2_bias, max_gate_2_bias + plunger_gate_2_bias]

fig, axs = plt.subplots(1, 3, figsize=(16, 5))
im0 = axs[0].imshow(to_image(mat_false) / np.max(mat_false), cmap="jet",
                     interpolation="bilinear", extent=extent, aspect="auto")
axs[0].set_title("overlap=False (baseline)"); fig.colorbar(im0, ax=axs[0])
im1 = axs[1].imshow(to_image(mat_true) / np.max(mat_true), cmap="jet",
                     interpolation="bilinear", extent=extent, aspect="auto")
axs[1].set_title("overlap=True"); fig.colorbar(im1, ax=axs[1])
diff = to_image(mat_true) / np.max(mat_true) - to_image(mat_false) / np.max(mat_false)
im2 = axs[2].imshow(diff, cmap="RdBu_r", extent=extent, aspect="auto")
axs[2].set_title("Difference (True - False, normalized)"); fig.colorbar(im2, ax=axs[2])
for ax in axs:
    ax.set_xlabel("$V_{g1}$ (V)")
axs[0].set_ylabel("$V_{g2}$ (V)")
fig.tight_layout()
fig.savefig(export_dir / "m06_overlap_comparison.png", dpi=150)
plt.close(fig)

max_abs_diff = np.max(np.abs(diff))
print(f"\n정규화된 CSD의 최대 절대 차이 (overlap True vs False): {max_abs_diff:.4f}")

metadata = {
    "classification": "C. Physics-oriented extension / validation study",
    "num_states": num_states,
    "coulomb_shape_overlap_false": list(coulomb_no_overlap.shape),
    "coulomb_shape_overlap_true": list(coulomb_overlap.shape),
    "coulomb_bytes_overlap_false": int(coulomb_no_overlap.nbytes),
    "coulomb_bytes_overlap_true": int(coulomb_overlap.nbytes),
    "t_coulomb_overlap_false_s": t_coulomb_false,
    "t_coulomb_overlap_true_s": t_coulomb_true,
    "t_sweep_overlap_false_s": t_sweep_false,
    "t_sweep_overlap_true_s": t_sweep_true,
    "max_abs_normalized_csd_diff": float(max_abs_diff),
    "temperature_K": temperature_spec,
    "gate_bias_range_V": [min_gate_1_bias, max_gate_1_bias],
}
with open(export_dir / "m06_overlap_comparison_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
