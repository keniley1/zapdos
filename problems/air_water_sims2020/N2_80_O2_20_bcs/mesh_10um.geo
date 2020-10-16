dom0Mult = 1e3;
dom1Mult = 1e6;

Point(1) = {0, 0, 0, 2e-9 * dom0Mult};
Point(2) = {.5e-3 * dom0Mult, 0, 0, 5e-6 * dom0Mult};
Point(3) = {1e-3 * dom0Mult, 0, 0, 2e-9 * dom0Mult};
Point(4) = {(1e-3 * dom0Mult) + (100e-9 * dom1Mult), 0, 0, 2e-9 * dom1Mult};
Point(5) = {(1e-3 * dom0Mult) + (5e-6 * dom1Mult), 0, 0, 0.1e-6 * dom1Mult};
Point(6) = {(1e-3 * dom0Mult) + (10e-6 * dom1Mult), 0, 0, 2e-9 * dom1Mult};
//Point(5) = {1e-3 + 5e-6, 0, 0, 0.025e-6};
//Point(6) = {1e-3 + 10e-6, 0, 0, 0.05e-6};

Line(11) = {1,2};
Line(12) = {2,3};
Line(13) = {3,4};
Line(14) = {4,5};
Line(15) = {5,6};

Physical Line(0) = {11,12};
Physical Line(1) = {13,14,15};
