__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
from qtcad.device import constants as ct
from qtcad.device.poisson_linear import Solver as LinearPoissonSolver
from qtcad.device.poisson_linear import SolverParams as LinearPoissonSolverParams
from qtcad.device import io
from double_dot_fdsoi import get_double_dot_fdsoi, save_slice

# Files -----------------------------------------------------------------------
script_dir = pathlib.Path(__file__).parent.resolve()

path_mesh = script_dir / "meshes"
path_out = script_dir / "output"
(path_out / "phi").mkdir(parents=True, exist_ok=True)

# NOTE: gmsh's XAO reader/writer fails on non-ASCII (Korean) paths with
# "Could not load XML file" (same bug hit in M12/M13/M14/M15/M17). The
# project directory is under a Korean-named folder, so geo_file/refined-mesh
# I/O is staged through an ASCII-only path in C:\temp.
path_mesh_ascii = pathlib.Path(r"C:\temp\m17_dqdfdsoi")
path_mesh_ascii.mkdir(parents=True, exist_ok=True)
import shutil as _shutil

_shutil.copy(path_mesh / "dqdfdsoi.xao", path_mesh_ascii / "dqdfdsoi.xao")

mesh_file = str(path_mesh / "dqdfdsoi.msh")
geo_file = str(path_mesh_ascii / "dqdfdsoi.xao")
ref_file = str(path_mesh_ascii / "refined_dqdfdsoi.msh")
pos_file = str(path_mesh / "dqdfdsoi.pos")

# output files
phi_out_file = str(path_out / "phi" / "phi.hdf5")
CB_out_file = str(path_out / "phi" / "CB.vtu")
CB_slice_file = str(path_out / "phi" / "CB_slice.png")

# Create the device -----------------------------------------------------------
d, set_region, qd_region = get_double_dot_fdsoi(mesh_file_name="dqdfdsoi.msh")

# Solve Linear Poisson equation -----------------------------------------------

# Solver parameters
p_adapt_params = LinearPoissonSolverParams()

p_adapt_params.tol = 1e-6
p_adapt_params.eta0 = 0.10
p_adapt_params.refined_region = qd_region + set_region
p_adapt_params.h_refined = 0.8
p_adapt_params.refined_mesh_filename = ref_file

# Solve the linear Poisson equation
p_adapt_solver = LinearPoissonSolver(
    d=d, solver_params=p_adapt_params, geo_file=geo_file
)
p_adapt_solver.solve()

# Save results ----------------------------------------------------------------

# Save the Poisson equation solution in .hdf5 format
arrays_dict = {"phi": d.phi}
io.save(phi_out_file, arrays_dict)

# Save the conduction band edge in .vtu format
CB = d.cond_band_edge() / ct.e
io.save(CB_out_file, {"EC [eV]": CB}, d.mesh)
# Save a slice of the conduction band edge in .png format
label_CB = "Conduction-band edge [eV]"
save_slice(d, CB, CB_slice_file, label_CB)

# Copy the adaptively-refined mesh back from the ASCII staging path to the
# project's meshes/ folder, since downstream scripts (3-coulomb_peaks.py via
# get_double_dot_fdsoi()'s default mesh_file_name="refined_dqdfdsoi.msh")
# expect it there.
_shutil.copy(ref_file, path_mesh / "refined_dqdfdsoi.msh")
