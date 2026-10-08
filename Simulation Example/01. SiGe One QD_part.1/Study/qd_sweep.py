import pathlib
from qtcad.device.mesh3d import Mesh

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = str(script_dir / "meshes" / "multiscale.msh")
scaling = 1e-9

mesh = Mesh(scaling, path_mesh) 

from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device

# 소자 정의: 전자에 대해 속박상태를 풀 것
d = Device(mesh, conf_carriers="e")
d.set_temperature(0.01) # 힌트 : 10 mK을 켈빈(K) 단위로

# tight-binding 계산과 맞추기 위한 밴드정렬 보정 (원본과 동일)
VBO= 0.68 * ct.e
mt.Ge.set_param("chi", mt.Si.chi - VBO + mt.Si.Eg - mt.Ge.Eg)

# SiGe 합금 조성
x = 0.3
SiGe =mt.SiGe
SiGe.set_alloy_composition(x)

# 영역별 재료 할당 (H단계에서 지정한 9개 이름과 정확히 일치해야 함)
d.new_region("oxide", mt.Al2O3)        # 힌트: Al2O3
d.new_region("cap_top", SiGe)
d.new_region("cap_bot", SiGe)
d.new_region("cap_qd", SiGe)
d.new_region("well", mt.Si)         # 힌트: Si
d.new_region("well_qd", mt.Si)
d.new_region("buffer_top", SiGe)
d.new_region("buffer_qd", SiGe)
d.new_region("buffer_bot", SiGe)

import numpy as np
from qtcad.device.mesh3d import SubMesh
from qtcad.device import SubDevice
from qtcad.device.poisson_linear import Solver as PoissonSolver
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.analysis import analyze_dot

Ew = 4.33 * ct.e
V_confinement = -1.9

qd = ["cap_qd", "well_qd", "buffer_qd"]
submesh = SubMesh(d.mesh, qd)          # 반복문 밖! 한 번만
poisson_solver = PoissonSolver(d=d)    # 이것도 반복문 밖! 한 번만

V_plunger_list = np.arange(-1.0, 1.0 + 1e-9, 0.1)
ground_energies = []
dot_sizes = []

for V_plunger in V_plunger_list:
    try:
        d.new_gate_bnd("gate_confinement", phi=V_confinement, Ew=Ew)
        d.new_gate_bnd("gate_plunger", phi=V_plunger, Ew=Ew)
        poisson_solver.solve()

        subdvc = SubDevice(d, submesh)
        subdvc.set_V_from_phi()

        schrod_solver = SchrodingerSolver(d=subdvc)
        schrod_solver.solve()

        E0 = subdvc.energies[0] / ct.e
        props = analyze_dot(subdvc, eigenstate=0, verbose=False)

        ground_energies.append(E0)
        dot_sizes.append(props["size"])
        print(f"V_plunger={V_plunger:.2f} V -> E0={E0:.4f} eV")
    except Exception as e:
        print(f"V_plunger={V_plunger:.2f} V -> 실패: {e}")
        ground_energies.append(np.nan)
        dot_sizes.append([np.nan, np.nan, np.nan])

import matplotlib.pyplot as plt

ground_energies = np.array(ground_energies)
dot_sizes = np.array(dot_sizes)   # shape: (전압 개수, 3) — 열마다 x,y,z 크기

# 결과 저장 (나중에 다시 불러와서 분석 가능하도록)
np.savez(str(script_dir / "output" / "sweep_results.npz"),
         V_plunger=V_plunger_list, ground_energy=ground_energies, dot_size=dot_sizes)

# 그래프 1: 전압에 따른 바닥상태 에너지
plt.figure()
plt.plot(V_plunger_list, ground_energies, "o-")
plt.xlabel("V_plunger (V)")
plt.ylabel("Ground state energy (eV)")
plt.savefig(str(script_dir / "output" / "sweep_energy.png"))

# 그래프 2: 전압에 따른 양자점 크기 (x, y, z 각각)
plt.figure()
plt.plot(V_plunger_list, dot_sizes[:, 0] / scaling, "o-", label="x")
plt.plot(V_plunger_list, dot_sizes[:, 1] / scaling, "o-", label="y")
plt.plot(V_plunger_list, dot_sizes[:, 2] / scaling, "o-", label="z")
plt.xlabel("V_plunger (V)")
plt.ylabel("Dot size (nm)")
plt.legend()
plt.savefig(str(script_dir / "output" / "sweep_size.png"))
