__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [EXTENSION] Phase 1 후처리 — 우물 검출 범위를 게이트 채널 구간으로 제한
#
# [NUMERICAL] [버그 수정] 2-array_poisson.py의 1차 검출은 전체 도메인(source+channel+
#   drain)에서 극값을 찾아, source/drain 접합부(heavily doped, E_C가 깊게 꺾이는
#   구간)에서 가짜 극값(edge artifact)이 다수 검출됨(7개, 목표 2개). source_drain_w=20nm
#   레드 구간을 제외하고 게이트 채널 구간(distance in [20, domain_l-20])만 검색하도록
#   수정 -- Poisson을 다시 풀 필요 없이 이미 저장된 linecut만 재사용.

import sys
import pathlib
import numpy as np
from scipy.signal import argrelextrema

N_DOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2
source_drain_w = 20  # 1-array_builder.py/2-array_poisson.py와 동일

script_dir = pathlib.Path(__file__).parent.resolve()
path_out = script_dir / "output"
d = np.load(path_out / f"array_n{N_DOTS}_bandedge_linecut.npz")
distance_m, Ec_eV = d["distance_m"], d["Ec_eV"]

total_len_nm = distance_m.max() / 1e-9
mask = (distance_m / 1e-9 > source_drain_w) & (distance_m / 1e-9 < total_len_nm - source_drain_w)
d_ch, E_ch = distance_m[mask], Ec_eV[mask]
order = np.argsort(d_ch)
d_ch, E_ch = d_ch[order], E_ch[order]

minima_idx = argrelextrema(E_ch, np.less_equal, order=5)[0]
minima_y = d_ch[minima_idx] / 1e-9
minima_E = E_ch[minima_idx]
merged_y, merged_E = [], []
for yv, Ev in zip(minima_y, minima_E):
    if merged_y and abs(yv - merged_y[-1]) < 10.0:
        if Ev < merged_E[-1]:
            merged_y[-1], merged_E[-1] = yv, Ev
        continue
    merged_y.append(yv)
    merged_E.append(Ev)

print(f"채널 구간(source/drain {source_drain_w}nm 레드 제외, distance in "
      f"[{source_drain_w},{total_len_nm-source_drain_w:.0f}]nm)에서 검출된 우물 개수: "
      f"{len(merged_y)} (목표: {N_DOTS}개)")
for i, (yv, Ev) in enumerate(zip(merged_y, merged_E)):
    print(f"  dot candidate {i+1}: linecut distance={yv:.2f}nm, E_C={Ev:.4f}eV")

if len(merged_y) == N_DOTS:
    print(f"\n[판정] N_dots={N_DOTS}개의 독립된 우물이 채널 구간 내에서 형성됨 확인 — "
          f"Phase 1 전기정전학적 형성 검증 통과.")
else:
    print(f"\n[판정] 여전히 목표와 다름 — 추가 조정 필요.")
