// Author: Pericles Philippopoulos
// Copyright 2023 NanoAcademic Technologies

// Use built-in meshing kernel
SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Dimensions of the slab
Lx = 60;
Ly = 60;
Lz = 12;

// Resolution along each direction
h = 3;
nlayers_y = Ly/h;
nlayers_z = Lz/h;

// Define geometric entities
Point(1) = {-Lx/2, -Ly/2, -Lz/2, h};
Point(2) = {Lx/2, -Ly/2, -Lz/2, h};
Line(1) = {1, 2};

// Extrude line to form plane
Extrude {0, Ly, 0} {
    Curve{1}; Layers{nlayers_y}; 
  }

// Extrude plane to form box
Extrude {0, 0, Lz} {
    Surface{5}; Layers{nlayers_z}; 
  }

// Define physical groups
Physical Surface("bnd") = {27, 26, 14, 22, 18, 5};
Physical Volume("domain") = {1};

// Save mesh file
Mesh 1;
Mesh 2;
Mesh 3;
// Generate a second order mesh
Mesh.SecondOrderLinear=1;
SetOrder 2;

Save "box_order2.msh";