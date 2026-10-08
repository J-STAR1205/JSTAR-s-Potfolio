// Author: Raphaël Prentki
// Copyright: 2025 Nanoacademic Technologies
// SiO2-SiGe-Si-SiGe heterostructure
// Remark: length units are nanometers

SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// characteristic lengths
alpha = 6;
fine = 10 / alpha;
medium = 40 / alpha;
coarse = 60 / alpha;

// number of layers in extrusions
layers_ox = alpha * 3;
layers_cap = alpha * 15;
layers_well_above = alpha * 5;
layers_well = alpha * 5;
layers_well_below = alpha * 5;
layers_buffer = alpha * 15;

// lengthscales
size = 200; // size of simulation domain along x and y axes
d = 100; // diameter of plunger gate
gap = 10; // gap between plunger gate and confinement gate
thick_ox = 5; // thickness of oxide
thick_well_above = 2;
thick_cap = 30 - thick_well_above; // thickness of cap layer
thick_well = 3; // thickness of well layer
thick_well_below = 2;
thick_buffer = 100 - thick_well_below; // thickness of buffer layer

// plunger gate
Point(1) = {0, 0, 0, fine};
Point(2) = {d/2, 0, 0, fine};
Point(3) = {0, d/2, 0, fine};
Point(4) = {-d/2, 0, 0, fine};
Point(5) = {0, -d/2, 0, fine};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// annular gap between plunger gate and confinement gate
Point(6) = {d/2+gap, 0, 0, medium};
Point(7) = {0, d/2+gap, 0, medium};
Point(8) = {-d/2-gap, 0, 0, medium};
Point(9) = {0, -d/2-gap, 0, medium};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2, 1};

// confinement gate
Point(10) = {-size/2, size/2, 0, coarse};
Point(11) = {size/2, size/2, 0, coarse};
Point(12) = {size/2, -size/2, 0, coarse};
Point(13) = {-size/2, -size/2, 0, coarse};
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3, 2};
 
// extrusions along z axis
thick_list = {thick_ox, thick_cap, thick_well_above, thick_well, thick_well_below, thick_buffer};
layers_list = {layers_ox, layers_cap, layers_well_above, layers_well, layers_well_below, layers_buffer};
// first extrusions
extrusion_01[] = Extrude {0, 0, -thick_list[0]} {Surface{1}; Layers{layers_list[0]};};
extrusion_02[] = Extrude {0, 0, -thick_list[0]} {Surface{2}; Layers{layers_list[0]};};
extrusion_03[] = Extrude {0, 0, -thick_list[0]} {Surface{3}; Layers{layers_list[0]};};
// remaining extrusions
For id In {1:5}
    extrusion_01[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_01[0]}; Layers{layers_list[id]};};
    extrusion_02[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_02[0]}; Layers{layers_list[id]};};
    extrusion_03[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_03[0]}; Layers{layers_list[id]};};
EndFor

Physical Surface("gate_plunger") = {1};
Physical Surface("gate_confinement") = {3};
Physical Volume("oxide") = {1, 2, 3};
Physical Volume("cap_top") = {4, 5, 6};
Physical Volume("cap_qd") = {7};
Physical Volume("cap_bot") = {8, 9};
Physical Volume("well_qd") = {10};
Physical Volume("well") = {11, 12};
Physical Volume("buffer_qd") = {13};
Physical Volume("buffer_top") = {14, 15};
Physical Volume("buffer_bot") = {16, 17, 18};

Mesh 1;
Mesh 2;
Mesh 3;
Save "multiscale.msh";