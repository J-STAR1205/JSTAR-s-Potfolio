// Author: Raphaël Prentki
// Copyright 2023, Nanoacademic Technologies Inc.

// Field-effect transistor (FET) with a (single) quantum dot (QD) in a
// fully-depleted silicon-on-insulator (FDSOI) geometry

// use built-in meshing kernel
SetFactory("OpenCASCADE");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// the y axis is the transport direction (source - channel - drain)
// the z axis is the heterostructure growth direction (back gate - bottom oxide - silicon - gate oxide - gate)
// the x axis is chosen such that (x, y, z) forms a Cartesian coordinate system

// geometrical parameters (in units of nm)
// along transport direction (y)
l_source = 20; // length of the n-doped silicon source (along y)
l_buffer_1 = 10; // length of the intrinsic silicon region between the doped source and the leftmost barrier gate (along y)
l_barrier_1 = 4; // length of the leftmost barrier gate (along y)
l_buffer_2 = 10; // length between the leftmost barrier gate and the plunger gate (along y)
l_plunger = 15; // length of the plunger gate (along y)
l_buffer_3 = 10; // length between the plunger gate and the rightmost barrier gate (along y)
l_barrier_2 = 4; // length of the rightmost barrier gate (along y)
l_buffer_4 = 10; // length of the intrinsic silicon region between the rightmost barrier gate and the doped drain
l_drain = 20; // length of the n-doped silicon drain (along y)
// along heterostructure growth direction (z)
t_top_ox = 2; // top oxide thickness (along z)
t_si = 3; // silicon thickness (along z)
t_bot_ox = 15; // bottom oxide thickness (along z)
// along x
w_left = 5; // width of the oxide on the left of silicon (along x)
w_si = 10; // width of the silicon (along x)
w_right = 5; // width of the oxide on the right of the silicon (along x)

// mesh granularity parameters
// number of principal layers along transport direction (y)
layers_source = 20; // number of layers of the n-doped silicon source (along y)
layers_buffer_1 = 10; // number of layers of the intrinsic silicon region between the doped source and the leftmost barrier gate (along y)
layers_barrier_1 = 4; // number of layers of the leftmost barrier gate (along y)
layers_buffer_2 = 10; // number of layers between the leftmost barrier gate and the plunger gate (along y)
layers_plunger = 15; // number of layers of the plunger gate (along y)
layers_buffer_3 = 10; // number of layers between the plunger gate and the rightmost barrier gate (along y)
layers_barrier_2 = 4; // number of layers of the rightmost barrier gate (along y)
layers_buffer_4 = 10; // number of layers of the intrinsic silicon region between the rightmost barrier gate and the doped drain
layers_drain = 20; // number of layers of the n-doped silicon drain (along y)
// characteristic lengths
coarse = 3;
medium = 2;
fine = 1;

// source surface
Point(1) = {-w_si/2, 0, -t_top_ox, medium};
Point(2) = {w_si/2, 0, -t_top_ox, medium};
Point(3) = {w_si/2, 0, -t_top_ox - t_si, medium};
Point(4) = {-w_si/2, 0, -t_top_ox - t_si, medium};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// top oxide surface normal to the y axis, on the source side
Point(5) = {-w_si/2 - w_left, 0, 0, fine};
Point(6) = {w_si/2 + w_right, 0, 0, fine};
Point(7) = {w_si/2 + w_right, 0, -t_top_ox, fine};
Point(8) = {-w_si/2 - w_left, 0, -t_top_ox, fine};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 2};
Line(8) = {1, 8};
Line(9) = {8, 5};
Curve Loop(2) = {5, 6, 7, -1, 8, 9};
Plane Surface(2) = {2};

// bottom and side oxide surface normal to the y axis, on the source side
Point(9) = {w_si/2 + w_right, 0, -t_top_ox - t_si - t_bot_ox, coarse};
Point(10) = {-w_si/2 - w_left, 0, -t_top_ox - t_si - t_bot_ox, coarse};
Line(10) = {7, 9};
Line(11) = {9, 10};
Line(12) = {10, 8};
Curve Loop (3) = {10, 11, 12, -8, -4, -3, -2, -7};
Plane Surface(3) = {3};

// extrusions along y axis
// lists containing extrusion lengths and their numbers of layers
l_list = {l_source, l_buffer_1, l_barrier_1, l_buffer_2, l_plunger, l_buffer_3, l_barrier_2, l_buffer_4, l_drain};
layers_list = {layers_source, layers_buffer_1, layers_barrier_1, layers_buffer_2, layers_plunger, layers_buffer_3, layers_barrier_2, layers_buffer_4, layers_drain};
// first extrusion
extrusion_si[] = Extrude {0, l_source, 0} {Surface{1}; Layers{layers_source};};
extrusion_top_ox[] = Extrude {0, l_source, 0} {Surface{2}; Layers{layers_source};};
extrusion_side_bot_ox[] = Extrude {0, l_source, 0} {Surface{3}; Layers{layers_source};};
// remaining 8 extrusions
For id In {1:8}
    extrusion_si[] = Extrude {0, l_list[id], 0} {Surface{extrusion_si[0]}; Layers{layers_list[id]};};
    extrusion_top_ox[] = Extrude {0, l_list[id], 0} {Surface{extrusion_top_ox[0]}; Layers{layers_list[id]};};
    extrusion_side_bot_ox[] = Extrude {0, l_list[id], 0} {Surface{extrusion_side_bot_ox[0]}; Layers{layers_list[id]};};
EndFor

// physical groups
Physical Surface("barrier_gate_1_bnd") = {51};
Physical Surface("plunger_gate_bnd") = {93};
Physical Surface("barrier_gate_2_bnd") = {135};
Physical Surface("source_bnd") = {1};
Physical Surface("drain_bnd") = {176};
Physical Surface("back_gate_bnd") = {17, 38, 59, 80, 101, 122, 143, 164, 185};
Physical Volume("top_ox") = {26, 23, 20, 17, 14, 11, 8, 5, 2};
Physical Volume("source") = {1};
Physical Volume("channel") = {22, 19, 16, 13, 10, 7, 4};
Physical Volume("drain") = {25};
Physical Volume("side_bot_ox") = {27, 24, 21, 18, 15, 12, 9, 6, 3};

// mesh
Mesh 1;
Mesh 2;
Mesh 3;

// save
Save "fdsoi_negf.msh";