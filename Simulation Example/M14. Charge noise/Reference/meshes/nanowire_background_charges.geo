// Author: Felix Beaudoin
// Copyright 2020 NanoAcademic Technologies
// Nanowire quantum dot

// Remarks: length units are nanometers

// Use built-in meshing kernel
SetFactory("Built-in");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Characteristic element lengths
h = 1;

// Number of layers along z axis
layers_sd = 5;
layers_channel = 25;

// Geometry parameters
L = 10;     // Channel length
r0 = 1.5;   // Channel radius
tox = 1;    // Oxide thickness
Lc = 5;     // Source/drain thickness

// Bottom drain surface
Point(1) = {0,0,0,h};
Point(2) = {r0,0,0,h};
Point(3) = {0,r0,0,h};
Point(4) = {-r0,0,0,h};
Point(5) = {0,-r0,0,h};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Bottom oxide surface
Point(6) = {(r0+tox),0,0,h};
Point(7) = {0,(r0+tox),0,h};
Point(8) = {-(r0+tox),0,0,h};
Point(9) = {0,-(r0+tox),0,h};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2,1};

// Extrude to form the drain
out[] = Extrude {0,0,Lc} {Surface{1}; Layers{layers_sd};};
top_drain = out[0];
volume_drain = out[1];
For i In {0:3}
    sides_drain[i] = out[i+2];
EndFor

out[] = Extrude {0,0,Lc} {Surface{2}; Layers{layers_sd};};
top_drain_ox = out[0];
volume_drain_ox = out[1];
For i In {0:7}
    sides_drain_ox[i] = out[i+2];
EndFor

// Physical groups of the drain
Physical Surface("drain_bnd") = {1};
Physical Volume("drain") = {volume_drain};
Physical Volume("drain_ox") = {volume_drain_ox};

// Extrude to form the channel
out[] = Extrude {0,0,L} {Surface{top_drain}; Layers{layers_channel};};
top_channel = out[0];
volume_channel = out[1];
For i In {0:3}
    sides_channel[i] = out[i+2];
EndFor

out[] = Extrude {0,0,L} {Surface{top_drain_ox}; Layers{layers_channel};};
top_channel_ox = out[0];
volume_channel_ox = out[1];
For i In {0:7}
    sides_channel_ox[i] = out[i+2];
EndFor

// Physical groups of the channel
Physical Volume("channel") = {volume_channel};
Physical Surface("gate_bnd") = {sides_channel_ox[0],
    sides_channel_ox[1], sides_channel_ox[2], sides_channel_ox[3]};
Physical Volume("channel_ox") = {volume_channel_ox};

// Extrude to form the source
out[] = Extrude {0,0,Lc} {Surface{top_channel}; Layers{layers_sd};};
top_source = out[0];
volume_source = out[1];
For i In {0:3}
    sides_source[i] = out[i+2];
EndFor

out[] = Extrude {0,0,Lc} {Surface{top_channel_ox}; Layers{layers_sd};};
top_source_ox = out[0];
volume_source_ox = out[1];
For i In {0:7}
    sides_source_ox[i] = out[i+2];
EndFor

// Physical groups of the source
Physical Surface("source_bnd") = {top_source};
Physical Volume("source") = {volume_source};
Physical Volume("source_ox") = {volume_source_ox};

// Physical groups of the charged interface
Physical Surface("interface", 201) = {sides_channel[0],
    sides_channel[1], sides_channel[2], sides_channel[3]};

Mesh 1;
Mesh 2;
Mesh 3;
Save "nanowire_background_charges.msh";
Save "nanowire_background_charges.geo_unrolled";