m1 = 1e-2; 

L = 1;
W = 0.25;


// Corners
Point(1) = {0, 0, 0, m1};
Point(2) = {L, 0, 0, m1};
Point(3) = {L, W, 0, m1};
Point(4) = {0, W, 0, m1};

// Interface 
Point(5) = {0, W/2, 0, m1};
Point(6) = {L, W/2, 0, m1};

Line(12) = {1,2};
Line(26) = {2,6};
Line(63) = {6,3};
Line(34) = {3,4};
Line(45) = {4,5};
Line(51) = {5,1};
Line(56) = {5,6};

Line Loop(1) = {12,26,-56,51};
Line Loop(2) = {56,63,34,45};

Plane Surface(0) = {1};
Plane Surface(1) = {2}; 


// Define boundaries
Physical Curve("lower_left") = {51};
Physical Curve("lower_right") = {26};
Physical Curve("top") = {34};
Physical Curve("bottom") = {12};
Physical Curve("upper_left") = {45};
Physical Curve("upper_right") = {63};
Physical Curve("middle") = {56};


Physical Surface(0) = {0};
Physical Surface(1) = {1};
Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
