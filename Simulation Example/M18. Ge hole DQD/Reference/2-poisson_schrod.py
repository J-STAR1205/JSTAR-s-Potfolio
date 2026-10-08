__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import io
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device, SubDevice
from qtcad.device.poisson import Solver as PoissonSolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device import analysis as an
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt

script_dir = Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
mesh_file_dir = mesh_dir / "ge_dqd.msh"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

# NOTE: gmsh's XAO reader fails on non-ASCII (Korean) paths with "Could not
# load XML file" (same bug hit in M12/M13/M14/M15/M17). The project
# directory is under a Korean-named folder, so geo_file is staged through an
# ASCII-only path in C:\temp.
import shutil as _shutil

mesh_dir_ascii = Path(r"C:\temp\m18_gedqd")
mesh_dir_ascii.mkdir(parents=True, exist_ok=True)
_shutil.copy(mesh_dir / "ge_dqd.xao", mesh_dir_ascii / "ge_dqd.xao")
xao_dir = mesh_dir_ascii / "ge_dqd.xao"

scaling = 1e-9  # nanometers
mesh = Mesh(scaling, mesh_file_dir)
# NOTE: Mesh.show()/Device.show() launch an interactive QTCAD/pyvista
# browser viewer that crashes outright in this headless Windows environment
# (confirmed via an OSError access-violation when Builder.view() was tried
# in M17 -- same underlying GUI machinery). Disabled here for the same
# reason; the device schematic is instead reconstructed independently from
# the Builder mask geometry already captured in M12's analysis notebook.
# mesh.show()

dvc = Device(mesh, conf_carriers="h", hole_kp_model="luttinger_kohn_foreman")

dvc.set_temperature(0.1)

# Assign the material to regions
dvc.new_region("substrate", mt.SiGe_DFT)
dvc.new_region("SiGe_barrier", mt.SiGe_DFT)
dvc.new_region("Ge_well", mt.Ge)
dvc.new_region("SiGe_cap", mt.SiGe_DFT)
dvc.new_region("oxide", mt.Al2O3)

# Dot-region volumes
dvc.new_region("Ge_well.dot_region", mt.Ge)
dvc.new_region("SiGe_barrier.dot_region", mt.SiGe_DFT)
dvc.new_region("SiGe_cap.dot_region", mt.SiGe_DFT)

# boundary conditions
Ew = mt.Ge.chi + 1.1 * mt.Ge.Eg
dvc.new_gate_bnd("P1", -0.6, Ew)
dvc.new_gate_bnd("P2", -0.6, Ew)
dvc.new_gate_bnd("BC", 0.9, Ew)
dvc.new_gate_bnd("BL", 0.5, Ew)
dvc.new_gate_bnd("BR", 0.5, Ew)

# visualize
# dvc.show()  # disabled: crashes headless (see note above)

dot_region_list = [
    "Ge_well.dot_region",
    "SiGe_barrier.dot_region",
    "SiGe_cap.dot_region",
]

dvc.set_dot_region(dot_region_list)

poisson_params = PoissonSolverParams()
poisson_params.tol = 1e-3
poisson_params.maxiter = 50
poisson_params.refined_region = dot_region_list
poisson_params.h_refined = 1.5
poisson_params.initial_ref_factor = 0.1
poisson_params.final_ref_factor = 0.75
poisson_params.refined_mesh_filename = mesh_dir_ascii / "refined_gedqd.msh"

poisson_solver = PoissonSolver(dvc, solver_params=poisson_params, geo_file=xao_dir)

poisson_solver.solve()

# Copy the adaptively-refined mesh back to the project's meshes/ folder for
# downstream scripts (3-leverarm.py, 4-CSD.py) that reference it.
_shutil.copy(poisson_params.refined_mesh_filename, mesh_dir / "refined_gedqd.msh")

io.save(out_dir / "electrostatic_potential.hdf5", dvc.phi)

begin = (-20e-9, 35e-9, 60e-9)
end = (125e-9, 35e-9, 60e-9)

an.plot_bands(
    dvc, begin, end, show_figure=True, path=out_dir / "band_edges_linecut.png"
)

x_cords, val_cut = an.linecut(dvc.mesh, dvc.vlnce_band_edge(), begin, end)
fig, ax = plt.subplots()
ax.plot(x_cords * 1e9, val_cut / ct.e, label="Valence band edge (eV)")
ax.axhline(0, color="k", linestyle="--", label="Fermi level")
ax.set_xlabel("x (nm)")
ax.set_ylabel("Energy (eV)")
ax.set_title("Valence band edge along line cut")
ax.legend()
fig.savefig(out_dir / "valence_band_linecut.png")
plt.show()

dvc.set_V_from_phi()

submesh = SubMesh(dvc.mesh, dot_region_list)
subdvc = SubDevice(dvc, submesh)

schrodinger_params = SchrodingerSolverParams()
schrodinger_params.num_states = 6

schrodinger_solver = SchrodingerSolver(subdvc, solver_params=schrodinger_params)
schrodinger_solver.solve()

subdvc.print_energies()
an.analyze_dot(subdvc, verbose=True)

psi = subdvc.eigenfunctions
psi2 = np.sum(np.abs(psi) ** 2, axis=2)

an.plot_slices(
    subdvc.mesh,
    psi2[:, 0],
    title="Ground state |psi|^2",
    show_figure=False,
    path=out_dir / "Ge_Ground_state.png",
)
an.plot_slices(
    subdvc.mesh,
    psi2[:, 1],
    title="First excited state |psi|^2",
    show_figure=False,
    path=out_dir / "Ge_First_excited_state.png",
)
an.plot_slices(
    subdvc.mesh,
    psi2[:, 2],
    title="Second excited state |psi|^2",
    show_figure=False,
    path=out_dir / "Ge_Second_excited_state.png",
)
