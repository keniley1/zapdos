xmin = -0.5;
xmax = 0.5;
ymin = -0.5;
ymax = 0.5;

lc =  0.03;

// Point locations are defined in meters
Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

// Upper I
Point(5) = {xmax * 2, ymax, 0, lc};
Point(6) = {xmax * 2, ymax * 2, 0, lc};
Point(7) = {xmin * 2, ymax * 2, 0, lc};
Point(8) = {xmin * 2, ymax, 0, lc};

// Lower I
Point(9) = {xmax * 2, ymin, 0, lc};
Point(10) = {xmax * 2, ymin * 2, 0, lc};
Point(11) = {xmin * 2, ymin * 2, 0, lc};
Point(12) = {xmin * 2, ymin, 0, lc};

Line(1) = {1,12};
Line(2) = {12,11}; // Bottom of domain (along liquid interface)
Line(3) = {11,10};
Line(4) = {10,9};
Line(5) = {9,2};
Line(6) = {2,3};
Line(7) = {3,5};
Line(8) = {5,6};
Line(9) = {6,7};
Line(10) = {7,8};
Line(11) = {8,4};
Line(12) = {4,1};

Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(0) = {1}; 

Physical Surface(0) = {0};
Physical Curve("bottom") = {3};
Physical Curve("top") = {9};
Physical Curve("top_left") = {10};
Physical Curve("top_right") = {8};
Physical Curve("bottom_left") = {2};
Physical Curve("bottom_right") = {4};
Physical Line("left_electrode") = {1,12,11};
Physical Line("right_electrode") = {5,6,7};


//Recombine Surface {0};

Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;

