__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.builder import Builder, MeshAlgorithm3D
from pathlib import Path

script_dir = Path(__file__).parent.resolve()
# NOTE: gdstk (used internally to read .oas layout files) fails to open
# files on non-ASCII (Korean) Windows paths with "Error opening input file" /
# "[GDSTK] Unable to open OASIS file for input." -- same class of bug seen
# earlier with gmsh's XAO loader. Worked around by keeping a copy of the
# layout on an ASCII-only path.
layout_dir = Path(r"C:\temp\m12_ge_hole")
mesh_dir = script_dir / "meshes"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

# Mesh characteristic lengths
char_len = 5
dot_char_len = char_len / 2

# Thickness of each material layer (in nm)
high_k_gate = 10
SiGe_cap = 5
Ge_well = 20
SiGe_barrier = 40
substrate_thick = 10

builder = Builder(name="Quantum Dot").load_layout(
    layout_dir / "ge_dqd.oas", cell_name="TOP"
)

builder.print_mask_tree()

(
    builder.set_mesh_size(char_len)
    .use_mask("substrate")
    # substrate
    .set_group_name("substrate")
    .extrude(substrate_thick)
    # SiGe barrier
    .set_group_name("SiGe_barrier")
    .extrude(SiGe_barrier)
    # ge well
    .set_group_name("Ge_well")
    .extrude(Ge_well)
    # sige cap
    .set_group_name("SiGe_cap")
    .extrude(SiGe_cap)
    # high-k gate
    .set_group_name("oxide")
    .extrude(high_k_gate)
)

# NOTE: builder.view(..., save="*.svg") on this geometry was taking 15+
# minutes / 1100+ CPU-s with no progress (hidden-line-removal rendering
# appears pathologically slow here, not a memory issue -- WorkingSet stayed
# ~200 MB throughout). Skipped for the first successful run; the mesh file
# is the actual deliverable needed for the M12 comparison against the
# hand-written .geo mesh. Re-enable selectively later if a quick picture is
# needed.
# builder.view(
#     surfaces=False,
#     volume_labels=True,
#     angles=(-90, 0, 85),
#     save=str(out_dir / "heterostructure.svg"),
#     zoom=1.1,
# )

dot_height = Ge_well + SiGe_cap / 2 + SiGe_barrier * (1 / 5)

starting_point = builder.get_z_from_group("Ge_well", bottom=True) - SiGe_barrier * (
    1 / 5
)

(
    builder.set_z(starting_point)
    .overlay_mode()
    .set_mesh_size(dot_char_len)
    .use_mask("dot-region")
    .group_from_shape()
    .extrude(dot_height)
)
# builder.view(
#     surfaces=False,
#     volume_labels=True,
#     save=str(out_dir / "model_dr.svg"),
#     zoom=1.1,
# )

(
    builder.use_mask("gates")
    .displace_mode()
    .set_mesh_size(dot_char_len)
    .set_z_from_group("oxide")
    .group_from_shape()
    .add_surface()
)

# builder.view(
#     surface_labels=True,
#     surfaces=True,
#     angles=(0, 0, 0),
#     save=str(out_dir / "surfs.svg"),
#     zoom=1.2,
#     groups=[
#         "BR",
#         "P1",
#         "P2",
#         "BC",
#         "BL",
#     ],
# )

builder.mesh(3, algorithm3d=MeshAlgorithm3D.HXT, show_gmsh_output=True).write(
    mesh_dir / "ge_dqd.msh"
)

# NOTE: same Korean-path bug as the .oas loader, but on the write side this
# time (gmsh's XAO writer). Write to an ASCII temp path, then copy back.
import shutil
xao_tmp = Path(r"C:\temp\m12_ge_hole\ge_dqd.xao")
builder.write(xao_tmp)
shutil.copy(xao_tmp, mesh_dir / "ge_dqd.xao")


# builder.view(
#     surfaces=False,
#     volume_labels=True,
#     angles=(-90, 0, 85),
#     save=str(out_dir / "mesh.png"),
# )
#
# builder.view(
#     surfaces=False,
#     volume_labels=True,
#     angles=(-90, 0, 85),
#     save=str(out_dir / "mesh_dr.png"),
#     groups=[
#         "SiGe_barrier.dot_region",
#         "Ge_well.dot_region",
#         "SiGe_cap.dot_region",
#     ],
# )
print("Done.")
