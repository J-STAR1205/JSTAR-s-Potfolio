__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M06 Step 4 (지침 11절) — Cross-lever-arm sensitivity study
#
# 분류: B. Parameter sensitivity study (device prediction 아님, geometry 변경 아님).
#   baseline(dominant~0.765, cross~0.006~0.007)에서 dominant는 고정하고 cross만
#   인위적으로 키워가며, CSD transition-line 기울기가 cross/dominant 비율과
#   analytic하게 일치하는지, 그리고 직사각형 모양이 실제로 평행사변형/honeycomb
#   쪽으로 바뀌는지 확인한다.
#
# [BASELINE] 재사용(지침 43절 Caching Rule): Poisson/Schrodinger/LeverArmSolver/
#   Coulomb matrix는 M06_overlap_comparison.py(Step 3)에서 이미 baseline과 교차
#   검증된 값을 그대로 로드해서 재사용한다 — 바뀌는 것은 lever-arm 행렬의 cross
#   성분뿐이므로 비싼 Poisson/Schrodinger 재계산이 필요 없다.
#
# [WARNING] 여기서 바뀌는 alpha_cross는 실제 geometry에서 유도된 값이 아니라
#   임의로 지정한 sensitivity-study 파라미터다. "실제 소자의 cross capacitance"를
#   의미한다고 해석하지 않는다 (지침 33절).

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from qtcad.device import constants as ct
from qtcad.transport.junction import Junction
from qtcad.transport.mastereq import add_spectrum
from qtcad.device.many_body import SolverParams as ManyBodySolverParams

script_dir = pathlib.Path(__file__).parent.resolve()
path_out = script_dir / "output"
export_dir = path_out / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

# --- [BASELINE] 캐시된 결과 로드 (Poisson/Schrodinger/Coulomb 재계산 없음) ---
energies = np.load(path_out / "dqdfdsoi_energies.npy")
baseline_la = np.load(path_out / "lever_arm_matrix.npy")  # shape (4, 2)
coulomb_no_overlap = np.load(path_out / "coulomb_mat_no_overlap.npy")
num_states = len(energies)
print(f"로드됨: energies({energies.shape}), lever_arm({baseline_la.shape}), "
      f"coulomb({coulomb_no_overlap.shape})")
print("Baseline lever arm matrix:\n", baseline_la)

# 각 행에서 "dominant"(큰 값) vs "cross"(작은 값) 성분 식별
dominant_mask = np.abs(baseline_la) > 0.1  # [NUMERICAL] baseline 기준 0.1V/V 임계값
print("Dominant mask:\n", dominant_mask)

temperature_spec = 10  # [BASELINE] double_dot_stability.py와 동일
gate_labels = ["source_bnd", "plunger_gate_1_bnd", "plunger_gate_2_bnd", "drain_bnd"]
plunger_gate_1_bias = 0.6
plunger_gate_2_bias = 0.6 + 1e-3

min_gate_1_bias, max_gate_1_bias = 30e-3, 140e-3
min_gate_2_bias, max_gate_2_bias = 30e-3, 140e-3
gate_1_biases = np.linspace(min_gate_1_bias, max_gate_1_bias, 90)
gate_2_biases = np.linspace(min_gate_2_bias, max_gate_2_bias, 90)


def build_modified_lever_arm(alpha_cross):
    """Dominant 성분은 baseline 그대로, cross 성분만 alpha_cross 크기로 치환
    (원래 부호는 유지)."""
    new_la = baseline_la.copy()
    for i in range(new_la.shape[0]):
        for j in range(new_la.shape[1]):
            if not dominant_mask[i, j]:
                sign = np.sign(baseline_la[i, j]) if baseline_la[i, j] != 0 else 1.0
                new_la[i, j] = sign * alpha_cross
    return new_la


def run_sweep(lever_arm_matrix, tag):
    mb_params = ManyBodySolverParams()
    mb_params.n_degen = 2
    mb_params.num_states = num_states
    mb_params.energies = energies
    mb_params.overlap = False
    mb_params.coulomb_mat = coulomb_no_overlap
    mb_params.alpha = lever_arm_matrix
    junc = Junction(many_body_solver_params=mb_params, temperature=temperature_spec,
                     contact_labels=gate_labels)
    mat = np.zeros((len(gate_1_biases), len(gate_2_biases)))
    t0 = time.time()
    for i1, v1 in enumerate(gate_1_biases):
        for i2, v2 in enumerate(gate_2_biases):
            junc.set_biases(gate_labels, np.array([0, v1, v2, 0]), verbose=False)
            mat[i1, i2] = add_spectrum(junc, temperature=temperature_spec)
    elapsed = time.time() - t0
    print(f"[{tag}] CSD sweep: {elapsed:.1f} s")
    return mat, elapsed


def empirical_slope(mat, window_frac=0.3):
    """국소 영역(중앙 윈도우)에서 가장 강한 ridge의 기울기를 선형회귀로 추정."""
    n1, n2 = mat.shape
    i0, i1 = int(n1 * (0.5 - window_frac / 2)), int(n1 * (0.5 + window_frac / 2))
    j0, j1 = int(n2 * (0.5 - window_frac / 2)), int(n2 * (0.5 + window_frac / 2))
    sub = mat[i0:i1, j0:j1]
    # 각 열(gate_1_bias)에서 최대 응답을 갖는 행(gate_2_bias) 인덱스를 ridge로 간주
    ridge_j = np.arange(sub.shape[1])
    ridge_i = np.argmax(sub, axis=0)
    # 선형회귀: i = slope_pixel * j + intercept
    A = np.vstack([ridge_j, np.ones_like(ridge_j)]).T
    slope_pixel, _ = np.linalg.lstsq(A, ridge_i, rcond=None)[0]
    dv1 = (max_gate_1_bias - min_gate_1_bias) / len(gate_1_biases)
    dv2 = (max_gate_2_bias - min_gate_2_bias) / len(gate_2_biases)
    slope_physical = slope_pixel * dv2 / dv1  # dV2/dV1
    return slope_physical


alpha_cross_values = [0.006, 0.02, 0.05, 0.10, 0.15, 0.20]
alpha_dominant_repr = np.abs(baseline_la[dominant_mask]).mean()
print(f"\n대표 dominant lever arm (평균): {alpha_dominant_repr:.4f}")

results = []
mats = {}
for alpha_cross in alpha_cross_values:
    la = build_modified_lever_arm(alpha_cross)
    mat, elapsed = run_sweep(la, f"alpha_cross={alpha_cross}")
    mats[alpha_cross] = mat
    slope_analytic = -alpha_cross / alpha_dominant_repr  # dV2/dV1 = -alpha1/alpha2 (지침 33절)
    slope_sim = empirical_slope(mat)
    rel_diff = (abs(slope_sim - slope_analytic) / abs(slope_analytic)
                if slope_analytic != 0 else np.nan)
    results.append({
        "alpha_cross": alpha_cross,
        "alpha_dominant": alpha_dominant_repr,
        "cross_to_dominant_ratio": alpha_cross / alpha_dominant_repr,
        "slope_analytic": slope_analytic,
        "slope_simulated": slope_sim,
        "relative_slope_diff_pct": rel_diff * 100 if not np.isnan(rel_diff) else np.nan,
        "wall_time_s": elapsed,
    })
    print(f"  alpha_cross={alpha_cross}: slope_analytic={slope_analytic:.3f}, "
          f"slope_sim={slope_sim:.3f}, rel.diff={rel_diff*100:.1f}%")

df = pd.DataFrame(results)
print("\n" + "=" * 80)
print(df.to_string(index=False))
df.to_csv(export_dir / "m06_cross_coupling_results.csv", index=False)

# --- 6-panel CSD 그리드 (직사각형 -> honeycomb 변화 시각화) ---
fig, axs = plt.subplots(2, 3, figsize=(15, 9))
extent = [min_gate_1_bias + plunger_gate_1_bias, max_gate_1_bias + plunger_gate_1_bias,
          min_gate_2_bias + plunger_gate_2_bias, max_gate_2_bias + plunger_gate_2_bias]
for ax, alpha_cross in zip(axs.flat, alpha_cross_values):
    mat = mats[alpha_cross]
    img = np.flip(np.transpose(mat), axis=0)
    ax.imshow(img / np.max(mat), cmap="jet", interpolation="bilinear", extent=extent,
              aspect="auto")
    ax.set_title(f"alpha_cross={alpha_cross} (ratio={alpha_cross/alpha_dominant_repr:.3f})",
                 fontsize=10)
    ax.set_xlabel("$V_{g1}$ (V)")
axs[0, 0].set_ylabel("$V_{g2}$ (V)")
axs[1, 0].set_ylabel("$V_{g2}$ (V)")
fig.suptitle("M06 cross-lever-arm sensitivity study (PARAMETER SENSITIVITY, not device prediction)",
             y=1.02)
fig.tight_layout()
fig.savefig(export_dir / "m06_cross_coupling_grid.png", dpi=150, bbox_inches="tight")
plt.close(fig)

metadata = {
    "classification": "B. Parameter sensitivity study",
    "alpha_dominant_fixed": float(alpha_dominant_repr),
    "alpha_cross_values_tested": alpha_cross_values,
    "baseline_alpha_cross": "0.006-0.007 (see lever_arm_matrix.npy)",
    "analytic_formula": "slope dV2/dV1 = -alpha_cross/alpha_dominant (constant-mu line)",
    "warning": "alpha_cross values are NOT derived from geometry; do not interpret as "
               "actual device predictions.",
}
with open(export_dir / "m06_cross_coupling_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
