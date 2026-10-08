__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M03 2차 라운드 Step 5 부속 — Ey = -dEc/dy 물리적 일관성 체크
#
# 분류: B. 수치 진단(diagnostic), 저장된 phi/Ec 원본은 전혀 건드리지 않는다.
#
# [WARNING] 지침 11절: "raw phi는 smoothing 금지, smoothing은 derivative 진단용으로만
#   별도 적용하고 저장/보고되는 phi 자체에는 적용하지 않는다." 이 스크립트는 이미 저장된
#   Ec(y) linecut(raw, 수정 전)을 읽기만 하고, 미분 계산 시에만 매끄러운 spline으로
#   국소적으로 미분해 Ey를 얻는다 — 저장 파일은 전혀 다시 쓰지 않는다.

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

script_dir = pathlib.Path(__file__).parent.resolve()
export_dir = script_dir / "output" / "analysis_exports"

cases = ["junction_local_h2.0", "junction_local_h1.0", "junction_local_h0.5"]
junction_left_nm, junction_right_nm = -45.0, 45.0

fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for tag, color in zip(cases, ["tab:green", "tab:orange", "tab:purple"]):
    d = np.load(export_dir / f"m03_raw_linecut_{tag}.npz")
    distance_m, Ec_eV = d["distance_m"], d["Ec_eV"]
    y_nm = distance_m / 1e-9 - 65.0
    order = np.argsort(y_nm)
    y_sorted, Ec_sorted = y_nm[order], Ec_eV[order]
    y_u, idx_u = np.unique(y_sorted, return_index=True)  # 중복 y 제거 (FEA linecutter 특성)
    Ec_u = Ec_sorted[idx_u]
    # [NUMERICAL] 진단 전용 Savitzky-Golay smoothing + gradient (저장 데이터는 변경하지 않음)
    win = min(31, len(Ec_u) - (1 - len(Ec_u) % 2))
    if win % 2 == 0:
        win -= 1
    win = max(win, 5)
    Ec_smooth = savgol_filter(Ec_u, window_length=win, polyorder=3)
    dEc_dy = np.gradient(Ec_smooth, y_u)  # eV/nm
    y_sorted, Ec_sorted = y_u, Ec_u
    Ey = -dEc_dy  # [PHYSICS] Ec = -e*phi + const -> dEc/dy = -e*dphi/dy = e*Ey -> Ey = -dEc/dy/e,
    # 여기서는 eV/nm 단위를 그대로 사용 (상대적 매끄러움/부호 일관성만 확인 목적)

    mask_l = (y_sorted > junction_left_nm - 8) & (y_sorted < junction_left_nm + 8)
    axs[0].plot(y_sorted[mask_l], Ec_sorted[mask_l], color=color, lw=1, label=f"{tag} (Ec)")
    axs[1].plot(y_sorted[mask_l], Ey[mask_l], color=color, lw=1.5, label=f"{tag} (Ey=-dEc/dy)")

    sign_changes = np.sum(np.diff(np.sign(Ey[mask_l])) != 0)
    print(f"[{tag}] junction 부근(±8nm) Ey 부호 변화 횟수: {sign_changes}, "
          f"Ey range=[{Ey[mask_l].min():.4f}, {Ey[mask_l].max():.4f}] eV/nm")

axs[0].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[0].set_xlabel("y (nm)"); axs[0].set_ylabel("$E_C$ (eV)")
axs[0].set_title("Ec(y) near source/channel junction"); axs[0].legend(fontsize=7)
axs[1].axvline(junction_left_nm, color="k", ls=":", lw=1)
axs[1].set_xlabel("y (nm)"); axs[1].set_ylabel("$-dE_C/dy$ (eV/nm, $\\propto E_y$)")
axs[1].set_title("Derivative diagnostic: sign/smoothness check"); axs[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(export_dir / "m03_Ey_consistency_check.png", dpi=150)
plt.close(fig)
print(f"\nExports: {export_dir / 'm03_Ey_consistency_check.png'}")
