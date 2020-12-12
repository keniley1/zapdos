xlim = 0.5e-3;
ylim = 0.1e-3;

Point(1) = {0, 0, 0, 1e-5};
Point(2) = {xlim, 0, 0, 1e-5};
Point(3) = {xlim, ylim, 0, 1e-6};
Point(4) = {0, ylim, 0, 2e-7};

// Interface Region
Point(5) = {xlim, ylim - 1e-5, 0, 2e-6};
Point(6) = {0, ylim - 1e-5, 0, 1e-6};

Line(11) = {1,2};
Line(12) = {2,5};
Line(13) = {5,3};
Line(14) = {3,4};
Line(15) = {4,6};
Line(16) = {6,1};

// Fake line, just for meshing
Line(10) = {5,6};
//Line(14) = {4,5};
//Line(15) = {5,6};

//Physical Line(0) = {11,12};
//Physical Line(1) = {13,14,15};

Line Loop(1) = {11,12,13,14,15,16};
Plane Surface(0) = {1};

//Physical Curve("electrode") = {1112,1218};
//Physical Curve("dielectric_gas") = {1617,112,1415};
//Physical Curve("ground") = {23};
Physical Curve("top") = {14};
Physical Curve("ground") = {13};

Physical Surface(0) = {0};

Line{10} In Surface{0};

Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
