Point(1) = {0, 0, 0, 2e-9};
Point(2) = {.5e-3, 0, 0, 12e-6};
Point(3) = {1e-3, 0, 0, 2e-9};
Point(4) = {1e-3 + 100e-9, 0, 0, 2e-9};
Point(5) = {1e-3 + 5e-6, 0, 0, 0.05e-6};
Point(6) = {1e-3 + 10e-6, 0, 0, 0.01e-9};

Line(11) = {1,2};
Line(12) = {2,3};
Line(13) = {3,4};
Line(14) = {4,5};
Line(15) = {5,6};

Physical Line(0) = {11,12};
Physical Line(1) = {13,14,15};
