__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device, SubDevice
from qtcad.device.poisson_linear import Solver as PoissonSolver
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device import io
from qtcad.device.analysis import plot_slice, analyze_dot

# Load the mesh
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = str(script_dir / "meshes" / "multiscale.msh")
scaling = 1e-9
mesh = Mesh(scaling, path_mesh)

# Create device
d = Device(mesh, conf_carriers="e")
d.set_temperature(1.6)

# Ensure band alignment between Ge and Si matches that assumed in TB
# calculations
VBO = 0.68 * ct.e  # valence band offset from Niquet et al.
mt.Ge.set_param("chi", mt.Si.chi - VBO + mt.Si.Eg - mt.Ge.Eg)

# SiGe alloy
x = 0.3  # Ge fraction in SiGe alloy
SiGe = mt.SiGe
SiGe.set_alloy_composition(x)

# Define material stack in heterostructure
d.new_region("oxide", mt.HfO2)
d.new_region("cap_top", SiGe)
d.new_region("cap_bot", SiGe)
d.new_region("cap_qd", SiGe)
d.new_region("well", mt.Si)
d.new_region("well_qd", mt.Si)
d.new_region("buffer_top", SiGe)
d.new_region("buffer_qd", SiGe)
d.new_region("buffer_bot", SiGe)

# Applied potentials
V_confinement = -1.9
V_plunger = -0.4

# Set up boundary conditions
Ew = 2.9 * ct.e  # workfunction of TiN (100)
d.new_gate_bnd("gate_confinement", phi=V_confinement, Ew=Ew)
d.new_gate_bnd("gate_plunger", phi=V_plunger, Ew=Ew)

# Set up and run adaptive Poisson solver
poisson_solver = PoissonSolver(d=d)
poisson_solver.solve()

# Save the Poisson equation solution in .hdf5 format
path_phi_hdf5 = str(script_dir / "output" / "multiscale_phi.hdf5")
arrays_dict = {"phi": d.phi}
io.save(path_phi_hdf5, arrays_dict)

# Save the conduction band edge and electron density in .vtu format
path_CB_vtu = str(script_dir / "output" / "multiscale_CB.vtu")
CB = d.cond_band_edge() / ct.e
out_dict = {"EC (eV)": CB}
io.save(path_CB_vtu, out_dict, d.mesh)

# Save slices of the conduction band edge and electron density
normal = (0.0, 0.0, 1.0)  # the plane of the slice is normal to the z axis
# the slice is just below the interface between the SiGe spacer and the Si QW
origin = (0.0, 0.0, -36.0 * scaling)
label_CB = "Conduction band edge (eV)"
path_CB_slice = str(script_dir / "output" / "multiscale_CB_slice.png")
plot_slice(
    d.mesh,
    CB,
    normal=normal,
    origin=origin,
    cb_axis_label=label_CB,
    path=path_CB_slice,
    show_figure=False,
)

# Get the potential energy from the band edge for usage in the Schrodinger
# solver
d.set_V_from_phi()

# Create a subdevice for the dot region
qd = ["cap_qd", "well_qd", "buffer_qd"]  # quantum dot regions
submesh = SubMesh(d.mesh, qd)
subdvc = SubDevice(d, submesh)

# Save confinement potential in subdevice for later use
subdvc.set_V_from_phi()
V_0 = subdvc.get_Vconf()

# Set up and run Schrodinger solver
schrod_solver = SchrodingerSolver(d=subdvc)
schrod_solver.solve()

# Print energies and analyze ground state
subdvc.print_energies()
analyze_dot(subdvc, eigenstate=0, verbose=True)

# Save the first few eigenstates in .vtu format
path_wf_vtu = str(script_dir / "output" / "multiscale_wf.vtu")
wf0 = subdvc.eigenfunctions[:, 0]
wf1 = subdvc.eigenfunctions[:, 1]
wf2 = subdvc.eigenfunctions[:, 2]
out_dict = {"Ground state": wf0, "1st excited": wf1, "2nd excited": wf2}
io.save(path_wf_vtu, out_dict, submesh)

# Save slice of the ground state
path_wf0_slice = str(script_dir / "output" / "multiscale_wf0_slice.png")
label_wf = "Ground state wavefunction"
plot_slice(
    submesh,
    wf0,
    normal=normal,
    origin=origin,
    cb_axis_label=label_wf,
    title="Ground state",
    path=path_wf0_slice,
    show_figure=False,
)

# Compute derivative of quantum dot confinement potential with respect to
# barrier gate voltage
V_epsilon = 0.01  # small increment in voltage to numerically compute derivative
V_confinement += V_epsilon
d.new_gate_bnd("gate_confinement", phi=V_confinement, Ew=Ew)
d.new_gate_bnd("gate_plunger", phi=V_plunger, Ew=Ew)
poisson_solver.solve()
subdvc = SubDevice(d, submesh)
subdvc.set_V_from_phi()
V_1 = subdvc.get_Vconf()
derivative = (V_1 - V_0) / V_epsilon
path_der_hdf5 = str(script_dir / "output" / "multiscale_der.hdf5")
arrays_dict = {"der": derivative}
io.save(path_der_hdf5, arrays_dict)  # save the derivative in .hdf5 format
path_der_slice = str(script_dir / "output" / "multiscale_der_slice.png")
plot_slice(
    submesh,
    derivative / ct.e,
    normal=normal,
    origin=origin,
    cb_axis_label="Derivative (meV/V)",
    path=path_der_slice,
    show_figure=False,
)  # save slice of the derivative
