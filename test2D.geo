// Inputs
length = 2;
width  = 2;
meshSize = 0.5;
//+
Point(1) = {0,      0,      0, meshSize};
Point(2) = {width,  0,      0, meshSize};
Point(3) = {width,  length, 0, meshSize};
Point(4) = {0,      length, 0, meshSize};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
//Physical Line( "lower" )  = { 1 };
//Physical Line( "right" )  = { 2 };
//Physical Line( "upper" )  = { 3 };
//Physical Line( "left" )   = { 4 };
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface {1};
Physical Surface( "domain" ) = {1};
