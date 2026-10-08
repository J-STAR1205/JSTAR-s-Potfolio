// Author: Pericles Philippopoulos
// Copyright 2020 NanoAcademic Technologies
// Remarks: length units are nanometers

// Use built-in meshing kernel
SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Geometry parameters
domain_thick = 40;     // Thickness of the simulation domain
domain_x = 120;        // Domain dimensions
d = 45;                // Distance between center of top gate and left domain edge
domain_y = 80;         // Domain dimensions
barrier_thick = 8;     // Barrier thickness
channel_thick = 20;    // Channel thickness
confined_thick = 6;

doped_thick = domain_thick - barrier_thick - channel_thick;  // Thickness of doped region
gate_radius = 10;      // Gate radius
side_gate_thick = 10;
side_gate_distance = 60;

// Characteristic element length
coarse_ms = 5;
gate_ms = 1;

// Define density of nodes along z
dens_barrier = 0.5;
dens_confined = 2;

// GATE (Plunger)
// The points forming the gate
Point(1) = {0, 0, domain_thick, coarse_ms};
Point(2) = {gate_radius, 0, domain_thick, coarse_ms};
Point(3) = {0, gate_radius, domain_thick, coarse_ms};
Point(4) = {-gate_radius, 0, domain_thick, coarse_ms};
Point(5) = {0, -gate_radius, domain_thick, coarse_ms};
// Define the gate
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
// Gate surface
Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};
// Create physical surface for the top gate
Physical Surface("top_gate_bnd") = {1};


// SIDE GATE
// The points forming the gate
side_gate_middle = domain_x/6 + side_gate_distance;
Point(6) = {side_gate_distance + side_gate_thick/2, -4*domain_y/10, domain_thick, coarse_ms}; 
Point(7) = {side_gate_distance - side_gate_thick/2, -4*domain_y/10, domain_thick, coarse_ms};
Point(8) = {side_gate_distance + side_gate_thick/2, 4*domain_y/10, domain_thick, coarse_ms};
Point(9) = {side_gate_distance - side_gate_thick/2, 4*domain_y/10, domain_thick, coarse_ms};
// Define the gate
Line(5) = {9, 8};
Line(6) = {8, 6};
Line(7) = {6, 7};
Line(8) = {7, 9};
// Gate surface
Curve Loop(2) = {8, 5, 6, 7};
Plane Surface(2) = {2};
// Create physical surface for the side gate
Physical Surface("side_gate_bnd") = {2};


//DOMAIN - barrier
Point(10) = {-d, -domain_y/2, domain_thick, coarse_ms};
Point(11) = {domain_x - d, -domain_y/2, domain_thick, coarse_ms};
Point(12) = {-d, domain_y/2, domain_thick, coarse_ms};
Point(13) = {domain_x - d, domain_y/2, domain_thick, coarse_ms};
Line(9) = {10, 12};
Line(10) = {12, 13};
Line(11) = {13, 11};
Line(12) = {11, 10};
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {1, 2, 3};

// Extrude to form the barrier
nlayers = dens_barrier * barrier_thick;

top_gate[] = Extrude {0, 0, -barrier_thick} {
    Surface{1}; Layers{nlayers}; 
};

side_gate[] = Extrude {0, 0, -barrier_thick} {
    Surface{2}; Layers{nlayers}; 
};

domain[] = Extrude {0, 0, -barrier_thick} {
    Surface{3}; Layers{nlayers}; 
};

Physical Volume("barrier") = {top_gate[1], side_gate[1], domain[1]};

//Extrude to form region under top gate forming a quantum dot.
nlayers = dens_confined * confined_thick;

confined[] = Extrude {0, 0, -confined_thick} {
    Surface{top_gate[0]}; Layers{nlayers}; 
};

Physical Volume("confined") = {confined[1]};

//DOMAIN - undoped
Point(110) = {-d, -domain_y/2, doped_thick, coarse_ms};
Point(111) = {domain_x - d, -domain_y/2, doped_thick, coarse_ms};
Point(112) = {-d, domain_y/2, doped_thick, coarse_ms};
Point(113) = {domain_x - d, domain_y/2, doped_thick, coarse_ms};

Line(135) = {110, 111};
Line(136) = {111, 113};
Line(137) = {113, 112};
Line(138) = {112, 110};
Line(139) = {112, 82};
Line(140) = {110, 73};
Line(141) = {113, 78};
Line(142) = {111, 74};
Curve Loop(4) = {137, 138, 135, 136};
Plane Surface(141) = {4};
Curve Loop(5) = {139, -68, -141, 137};
Plane Surface(142) = {5};
Curve Loop(6) = {138, 140, -69, -139};
Plane Surface(143) = {6};
Curve Loop(7) = {141, -67, -142, 136};
Plane Surface(144) = {7};
Curve Loop(8) = {142, -66, -140, 135};
Plane Surface(145) = {8};
Surface Loop(1) = {145, 144, 142, 143, 141, 140, 127, 131, 135, 139, 118, 56};
Volume(5) = {1};
Physical Volume("undoped") = {5};

//DOMAIN - doped
Point(210) = {-d, -domain_y/2, 0, coarse_ms};
Point(211) = {domain_x - d, -domain_y/2, 0, coarse_ms};
Point(212) = {-d, domain_y/2, 0, coarse_ms};
Point(213) = {domain_x - d, domain_y/2, 0, coarse_ms};//+
Line(143) = {210, 212};
Line(144) = {212, 213};
Line(145) = {213, 211};
Line(146) = {211, 210};
Line(147) = {210, 110};
Line(148) = {211, 111};
Line(149) = {213, 113};
Line(150) = {112, 212};
Curve Loop(9) = {144, 145, 146, 143};
Plane Surface(146) = {9};
Curve Loop(10) = {149, 137, 150, 144};
Plane Surface(147) = {10};
Curve Loop(11) = {143, -150, 138, -147};
Plane Surface(148) = {11};
Curve Loop(12) = {147, 135, -148, 146};
Plane Surface(149) = {12};
Curve Loop(13) = {148, 136, -149, 145};
Plane Surface(150) = {13};
Surface Loop(2) = {149, 148, 146, 147, 150, 141};
Volume(6) = {2};
Physical Volume("doped") = {6};
// Create physical surface for the bottom
Physical Surface("bottom") = {146};

// After creating the geometry, we now specify the mesh using size
// fields. For more information on size fields, see gmsh tutorial 10:
// https://gmsh.info/doc/texinfo/gmsh.html#t1
Field[1] = Cylinder;
Field[1].Radius = 1.5*gate_radius;
Field[1].VIn = gate_ms;     // Mesh size inside cylinder
Field[1].VOut = coarse_ms;  // Mesh size outside cylinder
Field[1].ZAxis = 2*(barrier_thick);
Field[1].ZCenter = domain_thick - (barrier_thick);
Background Field = 1;

Mesh 1;
Mesh 2;
Mesh 3;
Save "MOS_EDSR_example.msh";
