clc, clear all, close all

syms s_epsi c_epsi s_theta c_theta s_phi c_phi mg
R1 = [...
    1       0     0;
    0   c_phi s_phi;
    0  -s_phi c_phi;
];
R2 = [...
    c_theta   0 -s_theta;
          0   1        0;
    s_theta   0  c_theta;
];
R3 = [...
    c_epsi s_epsi 0;
   -s_epsi c_epsi 0;
         0      0 1;
];
R_i2b = R1*R2*R3;
R_b2i = R_i2b';

gravity_inertia = [0;0;-mg];
gravity_body = R_i2b*gravity_inertia;


