// Gmsh project created on Wed Jan 22 11:26:45 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0.0, 0.0, 0, 1.0};
//+
Point(2) = {50.0, 0.0, 0, 1.0};
//+
Point(3) = {0, 50.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("sud", 4) = {1};
//+
Physical Curve("est", 5) = {2};
//+
Physical Curve("ouest", 6) = {3};
//+
Physical Surface("triangle", 100) = {1};
//+
Physical Point("encastre", 101) = {1};
