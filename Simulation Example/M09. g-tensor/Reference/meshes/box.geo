// Author: Pericles Philippopoulos
// Copyright 2022 NanoAcademic Technologies
// Box

// Remarks: length units are nanometers

// Use built-in meshing kernel
SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Characteristic length
h = 1.2;

// Dimensions of the box
Lxy = 60;
Lz = 3;

// Corners of the bottom face
Point(1) = {-Lxy / 2, -Lxy / 2, 0, h};
Point(2) = {-Lxy / 2, Lxy / 2, 0, h};
Point(3) = {Lxy / 2, Lxy / 2, 0, h};
Point(4) = {Lxy / 2, -Lxy / 2, 0, h};
// Lines joining the corners
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

// Extrude surface to create a box
Extrude {0, 0, Lz} {
    Surface{1}; Layers{2*Lz/h}; 
}

// Define physical groups
Physical Surface("bnd") = {1, 13, 17, 21, 25, 26};
Physical Volume("domain") = {1};

Mesh 1;
Mesh 2;
Mesh 3;
Save "box.msh";