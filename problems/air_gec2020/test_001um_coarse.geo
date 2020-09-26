// Mesh scale factor
d0 = 1;
d1 = 1;
//d0 = 1e3;
//d1 = 5e4;

// Domain size
L_gas = 1e-3;
L_water = 1e-6;

// Gas Phase
Point(1) = {0, 0, 0, 2e-9 * d0};
Point(2) = {(L_gas / 2) * d0, 0, 0, 10e-6 * d0};
Point(3) = {L_gas * d0, 0, 0, 0.1e-9 * d0};
Line(1) = {1,2};
Line(2) = {2,3};

// Water Phase
Point(4) = {L_gas * d0 + 100e-9 * d1, 0, 0, 5e-9 * d1};
Point(5) = {L_gas * d0 + (L_water / 2) * d1, 0, 0, 0.05e-6 * d1};
Point(6) = {L_gas * d0 + L_water * d1, 0, 0, 1e-9 * d1};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Physical Line(0) = {1,2};
Physical Line(1) = {3,4,5};
