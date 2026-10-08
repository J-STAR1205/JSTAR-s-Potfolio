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

V_confinement = -1.9
V_plunger = -0.4

Ew = 4.33 * ct.e   # Ti 일함수

d.new_gate_bnd("gate_confinement", phi=V_confinement, Ew=Ew)
d.new_gate_bnd("gate_plunger", phi=V_plunger, Ew=Ew)

from qtcad.device.poisson_linear import Solver as PoissonSolver
from qtcad.device import io
from qtcad.device.analysis  import plot_slice

poisson_solver = PoissonSolver(d=d)
poisson_solver.solve()

# phi 저장
path_phi_hdf5 = str(script_dir / "output" / "qd_phi.hdf5")
arrays_dict = {"phi": d.phi}
io.save(path_phi_hdf5, arrays_dict)

# 전도대 바닥 계산 및 저장
path_CB_vtu = str(script_dir / "output" / "qd_CB.vtu")
CB = d.cond_band_edge() / ct.e
out_dict = {"EC (eV)": CB}
io.save(path_CB_vtu, out_dict, d.mesh)

# 슬라이스 이미지 저장
normal = (0.0, 0.0, 1.0)
origin = (0.0, 0.0, -46 * scaling)   # 힌트: 방금 계산한 -46
path_CB_slice = str(script_dir / "output" / "qd_CB_slice.png")
plot_slice(d.mesh, CB, normal=normal, origin=origin, cb_axis_label="Conduction band edge (eV)", path=path_CB_slice, show_figure=False)

from qtcad.device.mesh3d import SubMesh
from qtcad.device import SubDevice
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.analysis import analyze_dot

# Poisson 결과를 전위 에너지로 변환
d.set_V_from_phi()

# QD 영역만 서브메시로 추출
qd = ["cap_qd", "well_qd", "buffer_qd"]
submesh = SubMesh(d.mesh, qd)
subdvc = SubDevice(d, submesh)
subdvc.set_V_from_phi()
V_0 = subdvc.get_Vconf()

# Schrodinger 솔버
schrod_solver = SchrodingerSolver(d=subdvc)
schrod_solver.solve()

subdvc.print_energies()
analyze_dot(subdvc, eigenstate=0, verbose=True)

path_wf_vtu = str(script_dir / "output" / "qd_wf.vtu")
wf0 = subdvc.eigenfunctions[:, 0]     # 힌트 :  바닥상태
wf1 = subdvc.eigenfunctions[:, 1]
wf2 = subdvc.eigenfunctions[:, 2]
out_dict = {"Ground state": wf0, "1st excited": wf1, "2nd excited":wf2}
io.save(path_wf_vtu, out_dict, submesh)

path_wf0_slice = str(script_dir / "output" / "qd_wf0_slice.png")
plot_slice(
    submesh,           # 힌트: submesh
    wf0,
    normal=normal,
    origin=origin,
    cb_axis_label="Ground state wavefunction",
    title="Ground state",
    path=path_wf0_slice,
    show_figure=False,
)

V_epsilon = 0.01
V_confinement_2 = V_confinement + V_epsilon
d.new_gate_bnd("gate_confinement", phi=V_confinement_2, Ew=Ew)
d.new_gate_bnd("gate_plunger", phi=V_plunger, Ew=Ew)

poisson_solver.solve()

subdvc_2 = SubDevice(d, submesh)
subdvc_2.set_V_from_phi()
V_1 = subdvc_2.get_Vconf()

derivative = (V_1 - V_0) /V_epsilon

path_der_hdf5 = str(script_dir / "output" / "qd_der.hdf5")
io.save(path_der_hdf5, {"der": derivative})

path_der_slice = str(script_dir / "output" / "qd_der_slice.png")
plot_slice(submesh, derivative / ct.e, normal=normal, origin=origin, cb_axis_label="Derivative (eV/V)", path=path_der_slice, show_figure=False)

