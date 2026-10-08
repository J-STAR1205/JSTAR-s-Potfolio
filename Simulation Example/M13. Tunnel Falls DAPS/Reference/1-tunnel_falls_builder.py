__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
from qtcad.builder import Builder, Mask, Polygon

script_dir = Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
mesh_dir.mkdir(parents=True, exist_ok=True)
path_out = script_dir / "output"
path_local_out = path_out / "tunnel_falls_builder"
path_local_out.mkdir(parents=True, exist_ok=True)

# Mesh controls ---------------------------------------------------------------
char_len = 20.0
dot_char_len = 10.0

# Builder figure controls -----------------------------------------------------
# NOTE: builder.view(..., save=...) was observed to take 15+ minutes / 1100+
# CPU-s with no progress on a comparable Builder geometry in M12 (not a
# memory issue). Disabled here pre-emptively so the mesh-generation pipeline
# (the actual M13 deliverable) isn't blocked by slow intermediate renders.
save_builder_figures = False
view_angles = (-58.0, 20.0, 38.0)
top_view_angles = (0.0, 0.0, 90.0)

# Lateral dimensions (nm) -----------------------------------------------------
# x follows the gate row / transport direction. y points from the
# center-screening side toward the finger-gate side.
#
# References:
# - Marcks et al., Nat. Commun. 16, 11381 (2025).
# - George et al., Nano Lett. 25, 793-799 (2025).
# - Neyens et al., Nature 629, 80-85 (2024).
#
# The pitch is taken from the references above. The remaining
# lateral widths and distances are modeling choices chosen to approximate
# the layout of Fig. 1a of Marcks et al.
gate_pitch = 60.0
finger_length = 30.0
gate_span = 90.0
center_screen_gap = 20.0  # Gap between finger gates and central screening gate.
center_screen_width = 30.0

domain_length = 5.0 * gate_pitch
domain_width = gate_span + center_screen_gap + center_screen_width

# The finger-gate row starts at the +y edge of the footprint.
qubit_row_y = 0.5 * (domain_width - gate_span)

# The center-screening gate ends at the opposite (-y) edge of the footprint.
center_screen_length = domain_length
center_screen_y = -0.5 * (domain_width - center_screen_width)

# The buried screening gate is one continuous gate below the finger gates.
screening_length = domain_length
screening_span = gate_span - 2.0 * finger_length
screening_row_y = qubit_row_y + 0.5 * (gate_span - screening_span)

# The dot regions are pinned between the SG-facing and CS-facing edges.
screening_inner_edge_y = screening_row_y - 0.5 * screening_span
center_screen_inner_edge_y = center_screen_y + 0.5 * center_screen_width
dot_bottom_edge_y = center_screen_inner_edge_y
dot_top_edge_y = screening_inner_edge_y
dot_width = dot_top_edge_y - dot_bottom_edge_y
dot_row_y = 0.5 * (dot_top_edge_y + dot_bottom_edge_y)

# Dot regions extend from the edge of an outer barrier gate to the x = 0 plane.
dot_outer_edge_x = 2.0 * gate_pitch - 0.5 * finger_length
dot_length = dot_outer_edge_x
dot_center_offset = 0.5 * dot_length

# Vertical dimensions (nm) ----------------------------------------------------
# References:
# - Marcks et al., Nat. Commun. 16, 11381 (2025), Fig. 1b and device-modeling
#   text: 4.6 nm Si0.972Ge0.028 quantum well in bulk Si0.7Ge0.3.
# - George et al., Nano Lett. 25, 793-799 (2025), Fig. 2 process flow:
#   30-75 nm SiGe barrier, 1-2 nm Si cap, and 5-10 nm SiO2 + 5-10 nm HfO2.
# - Neyens et al., Nature 629, 80-85 (2024): optimized 60 nm-pitch Tunnel
#   Falls devices including a 50 nm SiGe barrier variant.
#
# The relaxed buffer and interlayer dielectric are not specified in the
# references above.
relaxed_buffer_slice_thick = 50.0
quantum_well_thick = 4.6
upper_barrier_thick = 50.0
si_cap_thick = 1.0
gate_oxide_sio2_thick = 5.0
gate_oxide_hfo2_thick = 5.0
screening_ild_thick = 5.0

# The dot refinement volumes extend only slightly into the adjacent SiGe.
dot_sige_extension = 6.0


def remove_unwanted_surfaces(builder: Builder) -> None:
    """Remove intermediate 2D surface labels while preserving named gates.

    Args:
        builder: Builder containing the current geometry.
    """
    for layer_name in [
        "relaxed_buffer",
        "quantum_well",
        "upper_barrier",
        "si_cap",
        "gate_oxide_sio2",
        "gate_oxide_hfo2",
        "screening_ild",
        "dot_",
    ]:
        builder.dissolve_physical_group(
            lambda group, layer_name=layer_name: (
                group.dim == 2 and layer_name in group.name
            )
        )


def save_shapes_if_enabled(builder: Builder, file_name: str) -> None:
    """Save selected masks in a top-down view when figures are enabled.

    Args:
        builder: Builder containing the selected masks to save.
        file_name: Name of the image file to save in ``path_local_out``.
    """
    if save_builder_figures:
        builder.view_shapes(
            save=path_local_out / file_name,
            show=False,
            angles=top_view_angles,
            font_size=16,
        )


def save_model_if_enabled(
    builder: Builder,
    file_name: str,
    *,
    show_mesh_faces: bool = False,
) -> None:
    """Save the 3D model in Gmsh when figures are enabled.

    Args:
        builder: Builder containing the model to save.
        file_name: Name of the image file to save in ``path_local_out``.
        show_mesh_faces: Whether to show mesh faces in the saved image.
    """
    if save_builder_figures:
        builder.view(
            surface_labels=True,
            volume_labels=True,
            angles=view_angles,
            save=path_local_out / file_name,
            show=False,
            font_size=16,
            show_mesh_faces=show_mesh_faces,
        )


# Masks -----------------------------------------------------------------------
footprint_mask = Mask("footprint")
footprint_mask.add_shape(
    Polygon.box(domain_length, domain_width, name="footprint").centered()
)

screening_mask = Mask("screening_gate")
screening_mask.add_shape(
    Polygon.box(screening_length, screening_span, name="SG")
    .centered()
    .translated(0.0, screening_row_y)
)

center_screen_mask = Mask("center_screen_gate")
center_screen_mask.add_shape(
    Polygon.box(center_screen_length, center_screen_width, name="CS")
    .centered()
    .translated(0.0, center_screen_y)
)

finger_gate_mask = Mask("finger_gates")
finger_gate_mask.add_shapes(
    [
        Polygon.box(finger_length, gate_span, name="B4")
        .centered()
        .translated(-2.0 * gate_pitch, qubit_row_y),
        Polygon.box(finger_length, gate_span, name="P5")
        .centered()
        .translated(-1.0 * gate_pitch, qubit_row_y),
        Polygon.box(finger_length, gate_span, name="B5")
        .centered()
        .translated(0.0, qubit_row_y),
        Polygon.box(finger_length, gate_span, name="P6")
        .centered()
        .translated(gate_pitch, qubit_row_y),
        Polygon.box(finger_length, gate_span, name="B6")
        .centered()
        .translated(2.0 * gate_pitch, qubit_row_y),
    ]
)

dot_mask = Mask("dots")
dot_mask.add_shapes(
    [
        Polygon.box(dot_length, dot_width, name="dot_left")
        .centered()
        .translated(-dot_center_offset, dot_row_y),
        Polygon.box(dot_length, dot_width, name="dot_right")
        .centered()
        .translated(dot_center_offset, dot_row_y),
    ]
)

# Setup builder ---------------------------------------------------------------
builder = Builder(name="Tunnel Falls DAPS")
builder.set_mesh_size(char_len)

builder.add_mask(footprint_mask).add_mask(screening_mask).add_mask(
    center_screen_mask
).add_mask(finger_gate_mask).add_mask(dot_mask)

builder.use_all_masks()
save_shapes_if_enabled(builder, "builder_mask_layout.png")

# Build the heterostructure ---------------------------------------------------
builder.use_mask("footprint")
builder.set_group_name("relaxed_buffer").extrude(relaxed_buffer_slice_thick)

builder.set_group_name("quantum_well").extrude(quantum_well_thick)

builder.set_group_name("upper_barrier").extrude(upper_barrier_thick)
builder.set_group_name("si_cap").extrude(si_cap_thick)
builder.set_group_name("gate_oxide_sio2").extrude(gate_oxide_sio2_thick)
builder.set_group_name("gate_oxide_hfo2").extrude(gate_oxide_hfo2_thick)

remove_unwanted_surfaces(builder)
save_model_if_enabled(builder, "builder_heterostructure_stack.png")

# Add the buried screening gate -----------------------------------------------
builder.use_mask("screening_gate")
builder.set_group_name("SG").add_surface()
save_model_if_enabled(builder, "builder_buried_screening_gate.png")

# Build the interlayer dielectric ---------------------------------------------
builder.fill_mode()
builder.use_mask("footprint")
builder.set_group_name("screening_ild").extrude(screening_ild_thick)
builder.displace_mode()
remove_unwanted_surfaces(builder)
save_model_if_enabled(builder, "builder_screening_ild.png")

# Add the upper finger gates and the center-screening gate --------------------
builder.use_mask("finger_gates")
builder.group_from_shape().add_surface()

builder.use_mask("center_screen_gate")
builder.group_from_shape().add_surface()
save_model_if_enabled(builder, "builder_upper_gates.png")

# Build the left and right dot regions ----------------------------------------
dot_stack_height = quantum_well_thick + 2.0 * dot_sige_extension

builder.overlay_mode()
builder.set_mesh_size(dot_char_len).minimum_mesh_size()
builder.set_z_from_group(
    "quantum_well",
    bottom=True,
    offset=-dot_sige_extension,
)
builder.use_mask("dots").group_from_shape().extrude(dot_stack_height)
remove_unwanted_surfaces(builder)
save_model_if_enabled(builder, "builder_dot_regions.png")

# Save the mesh ---------------------------------------------------------------
builder.mesh()
save_model_if_enabled(
    builder,
    "builder_final_mesh.png",
    show_mesh_faces=True,
)
builder.write(mesh_dir / "tunnel_falls_double_dot.msh")
# NOTE: gmsh's XAO writer fails on non-ASCII (Korean) Windows paths with
# "Could not open file" -- same bug seen in M12. Write to an ASCII temp path
# and copy back.
import shutil
xao_tmp = Path(r"C:\temp\m13_tunnel_falls") / "tunnel_falls_double_dot.xao"
xao_tmp.parent.mkdir(parents=True, exist_ok=True)
builder.write(xao_tmp)
shutil.copy(xao_tmp, mesh_dir / "tunnel_falls_double_dot.xao")
