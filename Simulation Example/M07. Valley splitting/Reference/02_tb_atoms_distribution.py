"""M07 Study 2단계 — SiGe 무작위 합금의 원자 배열과 조성 분포 확인

목적
------------------------------
파동함수가 아니라 "원자가 실제로 어떻게 놓였는지"를 본다.
  (1) z 방향 원자층별 Ge 비율   — 목표 조성(30 %)·계면 폭이 의도대로인가
  (2) 우물 바로 위 스페이서의 xy Ge 지도 — 국소 조성 요동이 얼마나 큰가
  (3) 국소 조성 히스토그램 vs 이항분포 — "무작위"가 정말 무작위인가
  (4) 원자 배치도 PNG — 측면 단면, 원자층 평면, Ge 3D, 격자 확대(결합선)
  (5) VESTA(.xyz)·ParaView(.vtu)용 파일 — 3D로 돌려 보고 싶을 때

이 단계는 슈뢰딩거 방정식을 풀지 않으므로 FEM 결과가 필요 없고 빠르다.
기본 상자(10 × 10 × 15 nm, 약 7만 원자)는 튜토리얼(약 90만 원자)보다 작다.

물리 메모
------------------------------
- 다이아몬드 격자에서 [001] 원자층 간격은 a/4 ≈ 0.136 nm(Si). 층 하나의
  면적 L×L 안 원자 수는 2·L²/a² 이므로 10 × 10 nm면 층당 약 680개.
- 각 원자가 독립적으로 확률 x로 Ge가 되면, N개 원자의 Ge 비율 표준편차는
  sqrt(x(1-x)/N). 층 전체(N≈680)면 ±1.8 %, 1 × 1 nm × 1 nm 셀(N≈50)이면 ±6.5 %.
  전자가 느끼는 것은 이런 국소 요동이고, 이것이 소자마다 VS가 달라지는 원인이다.

실행: (qtcad) python 02_tb_atoms_distribution.py
출력: output/02_<mode>_seed<seed>_*.png / .csv / .xyz / .vtu
"""

import pathlib
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qtcad.atoms.keating import Solver as KeatingSolver
from qtcad.atoms.keating import SolverParams as KeatingSolverParams
from qtcad.atoms.analysis import save_vtu

from atoms_builder import build_atoms, is_ge, positions_nm, ge_fraction_profile
from atoms_plot import plot_all

# =============================================================================
# 바꿔볼 파라미터
# =============================================================================
MODE = "rough"              # "rough"(튜토리얼 방식) / "graded"(조성 기울기)
SEED = 104                  # 합금 배치 시드 (튜토리얼과 같게 시작)
ROUGH_SEEDS = (0, 1)        # 상·하부 거친 계면 시드
X_GE = 0.30                 # 장벽 Ge 조성
WELL_TOP_NM = -35.0         # 튜토리얼 값. Philips 메시로 옮길 때 바꿀 것
WELL_BOT_NM = -38.0
INTERFACE_WIDTH_NM = 0.5    # graded 모드 tanh 특성 길이 w
HURST, RMS_NM = 0.3, 0.5    # rough 모드 계면 거칠기
BOX_XY_NM = 5.0             # 상자 x, y 범위 = [-BOX_XY_NM, +BOX_XY_NM]
BOX_Z_NM = (-44.0, -29.0)   # 상자 z 범위
RUN_KEATING = False         # True면 격자 이완(변형) 후 분석. 느려지지만 SiGe 층 부풂이 보임
SLAB_NM = 1.0               # xy 지도를 그릴 스페이서 두께 (우물 상부 계면 바로 위)
CELL_NM = 1.0               # xy 지도 셀 크기
SAVE_XYZ = True             # VESTA용
SAVE_VTU = True             # ParaView용

# =============================================================================
# 경로
# =============================================================================
script_dir = pathlib.Path(__file__).parent.resolve()
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)
tag = f"02_{MODE}_seed{SEED}"

nm = 1e-9
x = np.array([-BOX_XY_NM, BOX_XY_NM]) * nm
y = x.copy()
z = np.array(BOX_Z_NM) * nm

# =============================================================================
# 1. 원자 구조 생성 (+ 선택: Keating 이완)
# =============================================================================
t0 = time.perf_counter()
atoms, surfaces = build_atoms(
    x, y, z, WELL_TOP_NM * nm, WELL_BOT_NM * nm, X_GE, mode=MODE,
    interface_width=INTERFACE_WIDTH_NM * nm, hurst=HURST, rms=RMS_NM * nm,
    seed=SEED, rough_seeds=ROUGH_SEEDS)
t_build = time.perf_counter() - t0

t_keat = 0.0
if RUN_KEATING:
    kp = KeatingSolverParams()
    kp.verbose = True
    t0 = time.perf_counter()
    KeatingSolver(atoms=atoms, solver_params=kp).solve()
    t_keat = time.perf_counter() - t0

# 자료형·단위 확인 (처음 한 번은 출력을 꼭 볼 것)
species = np.asarray(atoms.atom_species)
print("atom_species 예시:", species[:6], "| dtype:", species.dtype)
print("atom_pos 예시 (원본 단위):\n", np.asarray(atoms.atom_pos)[:3])
pos = positions_nm(atoms)
ge = is_ge(species)
print(f"원자 수 {len(ge):,} | 전체 Ge 비율 {ge.mean()*100:.2f} % | "
      f"생성 {t_build:.1f} s, Keating {t_keat:.1f} s")

# =============================================================================
# 2. z 방향 원자층별 Ge 비율
# =============================================================================
dz = 0.5431 / 4  # 원자층 간격 (nm)
zc = pos[:, 2]
# 이완 전 원자는 층 위에 정확히 놓이므로, 빈 경계를 층 사이(반 층 비킴)에 둔다
edges = np.arange(zc.min() - dz / 2, zc.max() + dz, dz)
idx = np.digitize(zc, edges) - 1
nb = len(edges) - 1
n_layer = np.bincount(idx, minlength=nb)[:nb]
n_ge = np.bincount(idx, weights=ge.astype(float), minlength=nb)[:nb]
ok = n_layer > 0
z_mid = 0.5 * (edges[:-1] + edges[1:])[ok]
frac = (n_ge[ok] / n_layer[ok])
target = ge_fraction_profile(
    z_mid, WELL_TOP_NM, WELL_BOT_NM,
    INTERFACE_WIDTH_NM if MODE == "graded" else 0.0, X_GE)
sigma = np.sqrt(np.clip(target * (1 - target), 0, None) / n_layer[ok])

prof = pd.DataFrame({"z_nm": z_mid, "n_atoms": n_layer[ok], "n_Ge": n_ge[ok],
                     "Ge_frac": frac, "target": target, "binomial_sigma": sigma})
prof.to_csv(out_dir / f"{tag}_z_profile.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.fill_between(z_mid, (target - 2 * sigma) * 100, (target + 2 * sigma) * 100,
                color="0.85", label="target ± 2σ (binomial)")
ax.plot(z_mid, target * 100, "k--", lw=1, label="target profile")
ax.plot(z_mid, frac * 100, "o-", ms=3, lw=1, label="generated (per atomic layer)")
ax.axvspan(WELL_BOT_NM, WELL_TOP_NM, color="tab:blue", alpha=0.08, label="Si well (mean)")
ax.set_xlabel("z (nm)"); ax.set_ylabel("Ge fraction (%)")
ax.set_title(f"Ge profile per atomic layer — {MODE}, seed {SEED}")
ax.legend(fontsize=8); ax.grid(True)
fig.tight_layout(); fig.savefig(out_dir / f"{tag}_z_profile.png", dpi=150)

in_well = (zc > WELL_BOT_NM + 0.3) & (zc < WELL_TOP_NM - 0.3)
print(f"우물 내부(계면 0.3 nm 제외) Ge 원자 수: {ge[in_well].sum()} / {in_well.sum()}")

# =============================================================================
# 3. 우물 바로 위 스페이서 slab의 xy Ge 지도 + 국소 조성 히스토그램
# =============================================================================
slab = (zc > WELL_TOP_NM) & (zc <= WELL_TOP_NM + SLAB_NM)
xy_edges = np.arange(-BOX_XY_NM, BOX_XY_NM + 1e-9, CELL_NM)
n_cell, _, _ = np.histogram2d(pos[slab, 0], pos[slab, 1], bins=[xy_edges, xy_edges])
g_cell, _, _ = np.histogram2d(pos[slab & ge, 0], pos[slab & ge, 1], bins=[xy_edges, xy_edges])
with np.errstate(invalid="ignore", divide="ignore"):
    f_cell = g_cell / n_cell

fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
im = axs[0].imshow(f_cell.T * 100, origin="lower", cmap="viridis",
                   extent=[xy_edges[0], xy_edges[-1], xy_edges[0], xy_edges[-1]])
fig.colorbar(im, ax=axs[0], label="local Ge fraction (%)")
axs[0].set_xlabel("x (nm)"); axs[0].set_ylabel("y (nm)")
axs[0].set_title(f"Spacer slab {WELL_TOP_NM:g} → {WELL_TOP_NM + SLAB_NM:g} nm, "
                 f"{CELL_NM:g} nm cells")

vals = f_cell[np.isfinite(f_cell)].ravel()
n_mean = n_cell[n_cell > 0].mean()
x_slab = ge[slab].mean()
s_binom = np.sqrt(x_slab * (1 - x_slab) / n_mean)
axs[1].hist(vals * 100, bins=15, density=True, alpha=0.7, label="generated cells")
xx = np.linspace(vals.min(), vals.max(), 200)
axs[1].plot(xx * 100, np.exp(-(xx - x_slab) ** 2 / (2 * s_binom ** 2))
            / (np.sqrt(2 * np.pi) * s_binom * 100), "k--",
            label=f"binomial σ = {s_binom*100:.1f} %")
axs[1].set_xlabel("local Ge fraction (%)"); axs[1].set_ylabel("density")
axs[1].set_title(f"cell std {vals.std()*100:.1f} % vs binomial {s_binom*100:.1f} % "
                 f"(⟨N⟩ = {n_mean:.0f}/cell)")
axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(out_dir / f"{tag}_xy_map.png", dpi=150)
print(f"slab Ge 비율 {x_slab*100:.2f} %, 셀 간 표준편차 {vals.std()*100:.2f} % "
      f"(이항분포 예상 {s_binom*100:.2f} %)")

# =============================================================================
# 4. 원자 배치도 (Si 파랑, Ge 빨강) — output/<tag>_side_*.png, _top_layers.png,
#    _ge_3d.png, _lattice_top_interface.png
# =============================================================================
plot_all(pos, ge, WELL_TOP_NM, WELL_BOT_NM, str(out_dir / tag),
         title=f"[{MODE}, seed {SEED}]")

# =============================================================================
# 5. 계면 거칠기 지도 (rough 모드)
# =============================================================================
if MODE == "rough":
    xp = np.linspace(x[0], x[1], 100)
    yp = np.linspace(y[0], y[1], 100)
    surfaces["top"].plot(x=xp, y=yp, type="2D",
                         path=str(out_dir / f"{tag}_rough_top_2d.png"), show_figure=False)
    surfaces["top"].plot(x=xp, y=yp, type="3D",
                         path=str(out_dir / f"{tag}_rough_top_3d.png"), show_figure=False)

# =============================================================================
# 6. 3D 시각화 파일
# =============================================================================
if SAVE_XYZ:
    atoms.save_xyz(path=str(out_dir / f"{tag}.xyz"))   # VESTA: Si/Ge 색 구분
if SAVE_VTU:
    # ParaView: Representation = Point Gaussian, 'is_Ge'로 색칠, Threshold로 Ge만 보기
    save_vtu(atoms=atoms, out_dict={"is_Ge": ge.astype(float)},
             path=str(out_dir / f"{tag}.vtu"))

print(f"\n저장 위치: {out_dir}")
plt.show()
