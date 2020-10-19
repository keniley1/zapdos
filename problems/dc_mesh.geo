dom0Mult = 1;

d0_length = 0.38;

Point(1) = {0, 0, 0, 150e-6 * dom0Mult};
Point(2) = {d0_length/2 * dom0Mult, 0, 0, 7500e-6 * dom0Mult};
Point(3) = {d0_length * dom0Mult, 0, 0, 150e-6 * dom0Mult};

Line(11) = {1,2};
Line(21) = {2,3};

Physical Line(0) = {11,21};
