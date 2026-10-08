__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
from qtcad.builder import Builder, Polygon, Mask

script_dir = Path(__file__).parent.resolve()
pa_dir = script_dir / "figs"
pa_dir.mkdir(parents=True, exist_ok=True)

masks_file = pa_dir / "fdsoi_masks.svg"
oxide_mask_file = pa_dir / "fdsoi_oxide_mask.svg"
BOX_file = pa_dir / "fdsoi_BOX.svg"
BOX_rename_file = pa_dir / "fdsoi_BOX_renamed.svg"
channel_file = pa_dir / "fdsoi_channel.svg"
STI_file = pa_dir / "fdsoi_STI.svg"
STI_rename_file = pa_dir / "fdsoi_STI_renamed.svg"
sd_file = pa_dir / "fdsoi_sd.svg"
sd_rename_file = pa_dir / "fdsoi_sd_renamed.svg"
gates_file = pa_dir / "fdsoi_gates.svg"
dot_file = pa_dir / "fdsoi_dot.svg"
dot_rename_file = pa_dir / "fdsoi_dot_renamed.svg"
view_mesh_file = pa_dir / "fdsoi_mesh.png"

# view angles
v_angles = [-65, 0, -22.5]

# Length scales ---------------------------------------------------------------
char_len = 4  # characteristic mesh length

domain_w = 60  # width of the simulation domain
domain_l = 130  # length of the simulation domain
channel_w = 40  # width of the silicon channel
plunger_w = 15  # width of the plunger gates
barrier_w = 10  # width of the barrier gates
source_drain_w = 20  # width of the source and drain regions
pitch = 5  # distance between the gates
# Width of the quantum-dot regions
QD_w = plunger_w + barrier_w + 2 * pitch
# Thicknesses
box_thick = 10  # Buried oxide thickness
channel_thick = 10  # Channel thickness
EOT_thick = 2  # Equivalent oxide thickness (gate oxide)

# Masks -----------------------------------------------------------------------

# Silicon channel
channel = Polygon.box(channel_w, domain_l, name="channel").centered()
channel_mask = Mask("channel")
channel_mask.add_shape(channel)


# oxides - Full domain
oxide = Polygon.box(domain_w, domain_l, name="oxide").centered()
oxide_mask = Mask("oxide")
oxide_mask.add_shape(oxide)

# Source and drain
source = (
    Polygon.box(channel_w, source_drain_w, name="source")
    .centered()
    .translated(0, -(domain_l - source_drain_w) / 2)
)
drain = (
    Polygon.box(channel_w, source_drain_w, name="drain")
    .centered()
    .translated(0, (domain_l - source_drain_w) / 2)
)

sd_mask = Mask("source_drain")
sd_mask.add_shapes([source, drain])

# Gates
B1 = (
    Polygon.box(domain_w, barrier_w, name="B1")
    .centered()
    .translated(0, -(barrier_w + 2 * pitch + plunger_w))
)
P1 = (
    Polygon.box(domain_w, plunger_w, name="P1")
    .centered()
    .translated(0, -(barrier_w / 2 + pitch + plunger_w / 2))
)
B2 = Polygon.box(domain_w, barrier_w, name="B2").centered()
P2 = (
    Polygon.box(domain_w, plunger_w, name="P2")
    .centered()
    .translated(0, (barrier_w / 2 + pitch + plunger_w / 2))
)
B3 = (
    Polygon.box(domain_w, barrier_w, name="B3")
    .centered()
    .translated(0, barrier_w + 2 * pitch + plunger_w)
)

gate_mask = Mask("gates")
gate_mask.add_shapes([B1, P1, B2, P2, B3])

# Quantum dot regions
# (SET)
QD1 = (
    Polygon.box(channel_w + pitch, QD_w, name="QD1")
    .centered()
    .translated(0, -(barrier_w / 2 + pitch + plunger_w / 2))
)
# (Qubit)
QD2 = (
    Polygon.box(channel_w + pitch, QD_w, name="QD2")
    .centered()
    .translated(0, barrier_w / 2 + pitch + plunger_w / 2)
)

dot_mask = Mask("dots")
dot_mask.add_shapes([QD1, QD2])

# Setup builder ---------------------------------------------------------------
# NOTE: all intermediate builder.view()/.view_shapes() calls below are
# disabled. In M12/M13, Builder.view(..., save="*.svg"/"*.png") was found to
# be pathologically slow in this environment (15+ min, 1100+ CPU-seconds per
# call, confirmed NOT memory-related) because it launches a full VTK/pyvista
# off-screen render for very little benefit over a direct reconstruction from
# the mask polygon geometry. Here the masks are pure analytic Polygon.box()
# shapes defined above (domain_w, channel_w, plunger_w, barrier_w, pitch,
# etc.), so the schematic is reconstructed independently and exactly in the
# analysis notebook using those same parameters -- no file parsing needed.
# However, view() was found to have a load-bearing side effect: it forces
# Gmsh/OCC to synchronize pending boolean-fragment operations (e.g. the
# auto-named overlap group "B3.QD2_top" created when overlay_mode() dots
# intersect the B3 gate), which merge_groups()/dissolve_physical_group()
# downstream rely on. Removing view() entirely caused a
# "ValueError: Unknown physical group 'B3.QD2_top'" the first time this
# script was run. The cheap fix is to call get_groups(sync=True) (same
# synchronization, no rendering) wherever a view() call used to sit.
builder = Builder()
# Set some global parameters
builder.set_mesh_size(char_len).number_hull_surfaces(True)

# Add the masks to the builder
builder.add_mask(channel_mask).add_mask(oxide_mask).add_mask(sd_mask).add_mask(
    gate_mask
).add_mask(dot_mask)

# Visualize the masks
# builder.use_all_masks().view_shapes(save=masks_file, font_size=20)
builder.use_all_masks()

# Build the BOX ---------------------------------------------------------------
builder.use_mask("oxide")
# builder.view_shapes(save=oxide_mask_file, font_size=20)

builder.set_z(-box_thick - channel_thick).extrude(box_thick)  # BOX volume
builder.get_groups(sync=True)

# Rename the bottom oxide surface to back gate boundary
builder.rename_group("oxide_bottom", "back_gate_bnd")
# Remove unnecessary surfaces
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build the channel -----------------------------------------------------------
builder.use_mask("channel")
builder.set_z(-channel_thick).extrude(channel_thick)  # Channel volume
# Remove unnecessary surfaces
builder.dissolve_physical_group(lambda g: "channel" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build the STI regions -------------------------------------------------------
builder.fill_mode()

builder.use_mask("oxide").set_z_from_group("channel", bottom=True)
builder.extrude(channel_thick)  # STI volume
builder.get_groups(sync=True)

# Remove unnecessary surfaces
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)
builder.dissolve_physical_group(lambda g: "channel" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build source and drain regions ----------------------------------------------
builder.displace_mode()

# Build the volumes
builder.set_z_from_group("channel", bottom=True)
builder.use_mask("source_drain")
builder.extrude(channel_thick)  # Lead volumes
builder.get_groups(sync=True)

# Rename source and drain boundaries
builder.rename_group("source_side[3]", "source_bnd")
builder.rename_group("drain_side[1]", "drain_bnd")
# Remove unnecessary surfaces
builder.dissolve_physical_group(
    lambda g: "source" in g.name and g.dim == 2 and "bnd" not in g.name
)
builder.dissolve_physical_group(
    lambda g: "drain" in g.name and g.dim == 2 and "bnd" not in g.name
)

builder.get_groups(sync=True)

builder.number_hull_surfaces(False)

# Build gate-oxide layer ------------------------------------------------------
builder.use_mask("oxide").extrude(EOT_thick)  # gate-oxide volume
# Remove unnecessary surfaces
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)

# Deposit gates ---------------------------------------------------------------
builder.use_mask("gates")
builder.add_surface()

builder.get_groups(sync=True)

# Build dot regions -----------------------------------------------------------
builder.overlay_mode()
builder.use_mask("dots").set_z_from_group("channel", bottom=True)
builder.extrude(channel_thick + EOT_thick)

builder.get_groups(sync=True)

# Merge gate groups
builder.merge_groups(lambda g: "P1" in g.name and g.dim == 2, "plunger_gate_1_bnd")
builder.merge_groups(lambda g: "P2" in g.name and g.dim == 2, "plunger_gate_2_bnd")
builder.merge_groups(lambda g: "B1" in g.name and g.dim == 2, "barrier_gate_1_bnd")
builder.merge_groups(lambda g: "B2" in g.name and g.dim == 2, "barrier_gate_2_bnd")
builder.merge_groups(lambda g: "B3" in g.name and g.dim == 2, "barrier_gate_3_bnd")
# Remove unnecessary surfaces.
# NOTE: the naive predicate `"QD" in g.name` also matches compound overlap
# names like "B3.QD2_top" (the fragment boundary between gate B3 and dot
# region QD2) that the merge_groups() calls above already consumed into
# "barrier_gate_3_bnd" etc. Querying dissolve_physical_group() for a name
# that a prior merge_groups() call already renamed away raised
# "ValueError: Unknown physical group 'B3.QD2_top'" (observed when first
# running this script). The gate-prefixed overlap names are explicitly
# excluded here so only the genuinely leftover, non-gate "QD" surfaces
# (e.g. "QD1_top"/"QD2_top") are dissolved, matching the original intent.
builder.dissolve_physical_group(
    lambda g: "QD" in g.name
    and g.dim == 2
    and not any(gate in g.name for gate in ("P1", "P2", "B1", "B2", "B3"))
)

builder.get_groups(sync=True)

# Save the mesh -----------------------------------------------------------------
# NOTE: gmsh's XAO writer fails on non-ASCII (Korean) paths with "Could not
# open file" (same bug hit in M12/M13/M14/M15). The project directory is
# under a Korean-named folder, so we write through an ASCII-only staging path
# in C:\temp and copy the results back to the project folder afterward.
builder.mesh()
mesh_dir_ascii = Path(r"C:\temp\m17_fdsoi_set")
mesh_dir_ascii.mkdir(parents=True, exist_ok=True)
# NOTE: named dqdfdsoi.msh/.xao (not fdsoi_set.*) because downstream scripts
# (2-linear_poisson.py via double_dot_fdsoi.py, etc.) hardcode this filename.
# This is a separate mesh from M06's/M14's/M15's own dqdfdsoi.msh -- each
# lives in its own module's meshes/ folder, so there is no collision.
builder.write(mesh_dir_ascii / "dqdfdsoi.msh").write(mesh_dir_ascii / "dqdfdsoi.xao")
import shutil

shutil.copy(mesh_dir_ascii / "dqdfdsoi.msh", script_dir / "meshes" / "dqdfdsoi.msh")
shutil.copy(mesh_dir_ascii / "dqdfdsoi.xao", script_dir / "meshes" / "dqdfdsoi.xao")
