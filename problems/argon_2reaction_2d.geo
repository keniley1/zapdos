dom0Mult = 1e3;

// X extent
Point(1) = {0, 0, 0, 400e-8 * dom0Mult};
Point(2) = {.5e-3 * dom0Mult, 0, 0, 100e-6 * dom0Mult};
Point(3) = {1e-3 * dom0Mult, 0, 0, 400e-8 * dom0Mult};
Point(4) = {0, 1e-3 * dom0Mult, 0, 400e-8 * dom0Mult}; 
Point(5) = {0.5e-3 * dom0Mult, 1e-3 * dom0Mult, 0, 100e-6 * dom0Mult};
Point(6) = {1e-3 * dom0Mult, 1e-3 * dom0Mult, 0, 400e-8 * dom0Mult};
Line(12) = {1,2};
Line(23) = {2,3};
Line(41) = {4,1};
Line(54) = {5,4};
Line(65) = {6,5};
Line(36) = {3,6};

Curve Loop(1) = {12,23,36,65,54,41};
Plane Surface(0) = {1}; // Define the physical surface of the domain encompassed by curve loop 1

Physical Surface(0) = {0};
Physical Curve("bottom") = {12,23};
Physical Curve("left") = {41};
Physical Curve("top") = {54,65};
Physical Curve("right") = {36};

//Point(20) = {0, 0.000015, 0.00101, lc/1000};
//Point{20} In Surface{1}; 
Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
