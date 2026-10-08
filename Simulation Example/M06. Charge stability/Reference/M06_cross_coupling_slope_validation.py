__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M06 2차 라운드 Step 4 — Cross-lever-arm slope 정량 검증 (버그 수정판)
#
# 분류: B. Parameter sensitivity study (device prediction 아님).
#   double_dot_stability.py, M06_cross_coupling_study.py(1차)는 수정하지 않는다.
#
# [WARNING] [이전 버그] M06_cross_coupling_study.py의 `empirical_slope()`는 ROI
#   안에서 열(column)마다 단순 argmax를 썼는데, 두 개의 서로 다른 전이선 계열
#   (steep: dominant=V1,cross=V2 / shallow: dominant=V2,cross=V1)이 ROI 안에
#   동시에 존재해 ridge-hopping이 발생, slope_simulated가 수백~수천% 오차를 냈다.
#   이번에는 "이전 column에서 찾은 row 근처(±search_radius)에서만 다음 column을
#   찾는" continuity-constrained tracking(Method A)으로 단일 전이선만 격리한다.
#
# [PHYSICS] [축 정의 유도, 암기 아님] 이 소자의 lever-arm 행렬은 baseline에서
#   두 가지 row 유형을 번갈아 갖는다: (a) dominant=V2,cross=V1 (row0,2유형),
#   (b) dominant=V1,cross=V2 (row1,3유형). 등-화학퍼텐셜선(constant-mu line)의
#   정의 alpha_1*V1+alpha_2*V2=const로부터 dV2/dV1=-alpha_1/alpha_2이므로,
#   (a)유형 상태는 slope=-alpha_cross/alpha_dominant(완만, |slope|<1),
#   (b)유형 상태는 slope=-alpha_dominant/alpha_cross(가파름, |slope|>1)를 낸다.
#   이번 스크립트는 **완만한(shallow) 계열**을 추적 대상으로 삼아
#   slope_analytic = -alpha_cross/alpha_dominant 와 직접 비교한다.

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

energies = np.load(path_out / "dqdfdsoi_energies.npy")
baseline_la = np.load(path_out / "lever_arm_matrix.npy")
coulomb_no_overlap = np.load(path_out / "coulomb_mat_no_overlap.npy")
num_states = len(energies)
dominant_mask = np.abs(baseline_la) > 0.1

temperature_spec = 10
gate_labels = ["source_bnd", "plunger_gate_1_bnd", "plunger_gate_2_bnd", "drain_bnd"]
plunger_gate_1_bias = 0.6
plunger_gate_2_bias = 0.6 + 1e-3

min_gate_1_bias, max_gate_1_bias = 30e-3, 140e-3
min_gate_2_bias, max_gate_2_bias = 30e-3, 140e-3
NPTS = 90
gate_1_biases = np.linspace(min_gate_1_bias, max_gate_1_bias, NPTS)
gate_2_biases = np.linspace(min_gate_2_bias, max_gate_2_bias, NPTS)
dv1 = gate_1_biases[1] - gate_1_biases[0]
dv2 = gate_2_biases[1] - gate_2_biases[0]


def build_modified_lever_arm(alpha_cross):
    new_la = baseline_la.copy()
    for i in range(new_la.shape[0]):
        for j in range(new_la.shape[1]):
            if not dominant_mask[i, j]:
                sign = np.sign(baseline_la[i, j]) if baseline_la[i, j] != 0 else 1.0
                new_la[i, j] = sign * alpha_cross
    return new_la


def run_sweep(lever_arm_matrix):
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
    for i1, v1 in enumerate(gate_1_biases):
        for i2, v2 in enumerate(gate_2_biases):
            junc.set_biases(gate_labels, np.array([0, v1, v2, 0]), verbose=False)
            mat[i1, i2] = add_spectrum(junc, temperature=temperature_spec)
    return mat


def track_shallow_ridge(mat, search_radius=4, roi_j_frac=(0.1, 0.9)):
    """[Method A: ROI tracking / continuity-constrained]
    mat[i1, i2]: i1 indexes gate_1 (V1, columns to scan along), i2 indexes gate_2 (V2,
    row to track). Shallow family -> for each V1 column, there's a well-defined V2 row
    of max response that moves slowly as V1 increases. Track it with continuity.
    """
    n1, n2 = mat.shape
    j0, j1 = int(n1 * roi_j_frac[0]), int(n1 * roi_j_frac[1])
    # seed: global argmax in the FIRST column of the ROI
    col0 = mat[j0, :]
    prev_row = int(np.argmax(col0))
    ridge_cols, ridge_rows = [], []
    for j in range(j0, j1):
        lo = max(0, prev_row - search_radius)
        hi = min(n2, prev_row + search_radius + 1)
        local = mat[j, lo:hi]
        row = lo + int(np.argmax(local))
        ridge_cols.append(j)
        ridge_rows.append(row)
        prev_row = row
    return np.array(ridge_cols), np.array(ridge_rows)


def robust_linear_fit(x, y, n_iter=3, sigma_clip=2.0):
    """Simple iterative sigma-clipping linear fit (outlier-robust)."""
    mask = np.ones(len(x), dtype=bool)
    for _ in range(n_iter):
        p = np.polyfit(x[mask], y[mask], 1)
        resid = y - np.polyval(p, x)
        std = np.std(resid[mask])
        mask = np.abs(resid) < sigma_clip * std if std > 0 else mask
    p = np.polyfit(x[mask], y[mask], 1)
    return p[0], p[1], mask


alpha_cross_values = [0.006, 0.02, 0.05, 0.10, 0.15, 0.20]
alpha_dominant_repr = np.abs(baseline_la[dominant_mask]).mean()

results = []
ridge_plots = {}
for alpha_cross in alpha_cross_values:
    t0 = time.time()
    la = build_modified_lever_arm(alpha_cross)
    mat = run_sweep(la)
    elapsed = time.time() - t0

    ridge_cols, ridge_rows = track_shallow_ridge(mat)
    # pixel -> physical units
    x_phys = gate_1_biases[ridge_cols]
    y_phys = gate_2_biases[ridge_rows]
    slope_pixel, intercept_pixel, inlier_mask = robust_linear_fit(
        np.arange(len(ridge_cols), dtype=float), ridge_rows.astype(float))
    # convert pixel slope (d(row)/d(col)) to physical dV2/dV1
    slope_sim = slope_pixel * dv2 / dv1

    slope_analytic = -alpha_cross / alpha_dominant_repr
    rel_err = abs(slope_sim - slope_analytic) / abs(slope_analytic) * 100 if slope_analytic != 0 else np.nan

    results.append({
        "alpha_cross": alpha_cross, "alpha_dominant": alpha_dominant_repr,
        "slope_analytic": slope_analytic, "slope_simulated": slope_sim,
        "relative_error_pct": rel_err, "n_ridge_points": len(ridge_cols),
        "n_inliers": int(np.sum(inlier_mask)), "wall_time_s": elapsed,
    })
    ridge_plots[alpha_cross] = (mat, ridge_cols, ridge_rows, inlier_mask)
    np.savez(export_dir / f"m06_raw_csd_alpha{alpha_cross}.npz", mat=mat,
             ridge_cols=ridge_cols, ridge_rows=ridge_rows, inlier_mask=inlier_mask)
    print(f"[alpha_cross={alpha_cross}] slope_analytic={slope_analytic:.4f}, "
          f"slope_sim={slope_sim:.4f}, rel.err={rel_err:.1f}%, "
          f"inliers={np.sum(inlier_mask)}/{len(ridge_cols)}, time={elapsed:.1f}s")

results_df = pd.DataFrame(results)
print("\n" + "=" * 100)
print(results_df.to_string(index=False))
results_df.to_csv(export_dir / "m06_cross_coupling_slope_results.csv", index=False)

# --- Plot: ROI + extracted ridge for each case ---
# [NUMERICAL] [버그 수정] 이전 버전은 배경 이미지 extent에 plunger_gate_bias를 더했지만
#   ridge 좌표(x_phys/y_phys)는 gate_1_biases/gate_2_biases 원값 그대로 써서 둘이
#   서로 다른 좌표계로 그려짐(배경이 우상단 구석에 작게 표시되는 버그) — extent를
#   ridge 좌표와 동일한 gate_1_biases/gate_2_biases 범위로 맞춤.
fig, axs = plt.subplots(2, 3, figsize=(16, 10))
extent = [min_gate_1_bias, max_gate_1_bias, min_gate_2_bias, max_gate_2_bias]
for ax, alpha_cross in zip(axs.flat, alpha_cross_values):
    mat, ridge_cols, ridge_rows, inlier_mask = ridge_plots[alpha_cross]
    img = np.flip(np.transpose(mat), axis=0)
    ax.imshow(img / np.max(mat), cmap="gray", extent=extent, aspect="auto", alpha=0.7)
    x_phys = gate_1_biases[ridge_cols]
    y_phys = gate_2_biases[ridge_rows]
    ax.plot(x_phys[inlier_mask], y_phys[inlier_mask], "g.", ms=4, label="ridge (inlier)")
    ax.plot(x_phys[~inlier_mask], y_phys[~inlier_mask], "rx", ms=4, label="ridge (outlier)")
    row = results_df[results_df["alpha_cross"] == alpha_cross].iloc[0]
    ax.set_title(f"alpha_cross={alpha_cross}: analytic={row['slope_analytic']:.3f}, "
                 f"sim={row['slope_simulated']:.3f}", fontsize=9)
    ax.set_xlabel("$V_{g1}$ (V)"); ax.legend(fontsize=6)
axs[0, 0].set_ylabel("$V_{g2}$ (V)"); axs[1, 0].set_ylabel("$V_{g2}$ (V)")
fig.suptitle("M06 cross-coupling: extracted shallow-family ridge vs. analytic line", y=1.01)
fig.tight_layout()
fig.savefig(export_dir / "m06_ridge_tracking_validation.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- Plot: relative error vs alpha_cross ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(results_df["alpha_cross"], results_df["relative_error_pct"], "o-")
ax.axhline(5, color="r", ls="--", lw=1, label="5% (project practical criterion)")
ax.set_xlabel("alpha_cross"); ax.set_ylabel("relative slope error (%)")
ax.set_title("Analytic vs. simulated slope: relative error"); ax.legend(); ax.grid(True)
fig.tight_layout()
fig.savefig(export_dir / "m06_slope_relative_error.png", dpi=150)
plt.close(fig)

metadata = {
    "classification": "B. Parameter sensitivity study",
    "method": "Method A: continuity-constrained ridge tracking (search_radius=4 px), "
              "shallow-family only (dominant=V2,cross=V1 states)",
    "analytic_formula": "slope dV2/dV1 = -alpha_cross/alpha_dominant, derived from "
                        "constant-mu line alpha_1*V1+alpha_2*V2=const",
    "acceptance_criterion": "relative_error<5% (project-level, NOT official QTCAD spec)",
}
with open(export_dir / "m06_slope_validation_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
