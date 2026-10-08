"""M07 Study 2b — 저장된 .xyz 파일로 원자 배치도 그리기 (QTCAD 불필요)

용도
------------------------------
- 튜토리얼 SiGe Part 2가 저장한 `output/multiscale.xyz`(약 90만 원자, Keating 이완 후)를
  VESTA 없이 바로 그림으로 본다.
- 02 스크립트가 저장한 .xyz도 같은 방식으로 다시 그릴 수 있다.

그림 (02와 동일, Si 파랑 / Ge 빨강)
  <이름>_side_full.png / _side_zoom.png   측면 단면 (x–z)
  <이름>_top_layers.png                   원자층 평면도 3장 (스페이서 / 계면 / 우물 중앙)
  <이름>_ge_3d.png                        계면 부근 Ge 원자 3D
  <이름>_lattice_top_interface.png        상부 계면 격자 확대도 (결합선)
  <이름>_z_profile.png                    원자층별 Ge 비율

실행: python 02b_plot_atoms_from_xyz.py
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from atoms_plot import plot_all

# =============================================================================
# 바꿔볼 파라미터
# =============================================================================
XYZ_PATH = pathlib.Path(
    r"..\..\01. SiGe One QD_part.1\Reference\output\multiscale.xyz"   # 기준소자 multiscale_2.py 결과 (공유)
)
XYZ_UNIT = "auto"           # "auto" / "angstrom" / "nm" / "m"
WELL_TOP_NM = -35.0         # 튜토리얼 우물 (평균 계면 높이)
WELL_BOT_NM = -38.0

# =============================================================================
# 1. .xyz 읽기 (1행: 원자 수, 2행: 주석, 이후 "원소 x y z")
# =============================================================================
script_dir = pathlib.Path(__file__).parent.resolve()
path = (script_dir / XYZ_PATH).resolve() if not XYZ_PATH.is_absolute() else XYZ_PATH
assert path.exists(), f".xyz 파일이 없습니다: {path}"

with open(path, encoding="utf-8", errors="ignore") as f:
    n_atoms = int(f.readline().split()[0])
    comment = f.readline().strip()
data = pd.read_csv(path, sep=r"\s+", skiprows=2, header=None, nrows=n_atoms,
                   usecols=[0, 1, 2, 3], engine="c")
species = data[0].astype(str).str.strip().to_numpy()
xyz = data[[1, 2, 3]].to_numpy(dtype=float)
print(f"{path.name}: {n_atoms:,} atoms | 주석: {comment[:80]}")

# 단위 판정: 튜토리얼 z가 -50 ~ -5 nm이므로 |좌표| 최대가
# ~1e-8 → m, 수십 → nm, 수백 → Å
span = np.abs(xyz).max()
unit = XYZ_UNIT
if unit == "auto":
    unit = "m" if span < 1e-6 else ("nm" if span < 150 else "angstrom")
to_nm = {"m": 1e9, "nm": 1.0, "angstrom": 0.1}[unit]
pos = xyz * to_nm
print(f"좌표 최대 |값| {span:.3g} → 단위 '{unit}'로 판정 (틀리면 XYZ_UNIT 지정)")
print("범위 (nm): x {:.1f}~{:.1f}, y {:.1f}~{:.1f}, z {:.1f}~{:.1f}".format(
    pos[:, 0].min(), pos[:, 0].max(), pos[:, 1].min(), pos[:, 1].max(),
    pos[:, 2].min(), pos[:, 2].max()))

ge = species == "Ge"
print(f"Ge 비율 {ge.mean()*100:.2f} % ({ge.sum():,} / {len(ge):,})")

# =============================================================================
# 2. 배치도
# =============================================================================
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)
prefix = str(out_dir / f"02b_{path.stem}")
plot_all(pos, ge, WELL_TOP_NM, WELL_BOT_NM, prefix, title=f"[{path.name}]")

# 원자층별 Ge 비율. 이완 후에는 SiGe 층 간격이 늘어나 고정 폭 빈은 누적 오차가 생기므로,
# 정렬한 z에서 0.04 nm보다 큰 간격을 층 경계로 본다.
order = np.argsort(pos[:, 2], kind="stable")
zs, gs = pos[order, 2], ge[order]
brk = np.concatenate([[0], np.where(np.diff(zs) > 0.04)[0] + 1, [len(zs)]])
lz = np.array([zs[brk[i]:brk[i+1]].mean() for i in range(len(brk) - 1)])
lf = np.array([gs[brk[i]:brk[i+1]].mean() for i in range(len(brk) - 1)])
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(lz, lf * 100, "o-", ms=2, lw=0.8)
ax.axvspan(WELL_BOT_NM, WELL_TOP_NM, color="tab:blue", alpha=0.08)
ax.set_xlabel("z (nm)"); ax.set_ylabel("Ge fraction per layer (%)")
ax.set_title(f"Ge profile — {path.name}"); ax.grid(True)
fig.tight_layout(); fig.savefig(f"{prefix}_z_profile.png", dpi=150)

print(f"\n저장: {out_dir}  (파일 이름 02b_{path.stem}_*.png)")
plt.show()
