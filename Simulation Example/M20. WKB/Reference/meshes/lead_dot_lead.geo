// Author: Pericles Philppopoulos
// Copyright 2021 NanoAcademic Technologies
// Tunnelling test

// Remarks: length units are nanometers

// Use built-in meshing kernel
SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Dimesions
Lx = 300;
dot_length = 190;
barrier_length = 5;
lead_length = (Lx - dot_length - 2*barrier_length) / 2;

Ly = 100;
Lz = 20;

// Define density of nodes along z
dens_dot = 0.5;
dens_lead = 0.5;

// Characteristic element lengths
h = 2;

// Left-end (square)
// Channel
Point(1) = {0, -Ly/2, -Lz/2, h};
Point(2) = {0, -Ly/2, Lz/2, h};
Point(3) = {0, Ly/2, -Lz/2, h};
Point(4) = {0, Ly/2, Lz/2, h};

Line(1) = {3, 4};
Line(2) = {4, 2};
Line(3) = {2, 1};
Line(4) = {1, 3};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

// Extrude 
// Lead 
nlayers = dens_lead * lead_length;
lead[] = Extrude {lead_length, 0, 0} {
    Surface{1}; Layers{nlayers}; 
};

// Barrier 1 
nlayers = dens_dot * barrier_length;
barrier[] = Extrude {barrier_length, 0, 0} {
    Surface{lead[0]}; Layers{nlayers}; 
};

//dot
nlayers = dens_dot * dot_length;
dot[] = Extrude {dot_length, 0, 0} {
    Surface{barrier[0]}; Layers{nlayers}; 
};

// Barrier 2
nlayers = dens_dot * barrier_length;
barrier2[] = Extrude {barrier_length, 0, 0} {
    Surface{dot[0]}; Layers{nlayers}; 
};

// Lead 2
nlayers = dens_lead * lead_length;
lead2[] = Extrude {lead_length, 0, 0} {
    Surface{barrier2[0]}; Layers{nlayers}; 
};

// Physical Volumes
Physical Volume('lead') = {lead[1]};
Physical Volume('barrier') = {barrier[1]};
Physical Volume('dot') = {dot[1]};
Physical Volume('barrier2') = {barrier2[1]};
Physical Volume('lead2') = {lead2[1]};


Mesh 1;
Mesh 2;
Mesh 3;

Save "lead_dot_lead.msh";