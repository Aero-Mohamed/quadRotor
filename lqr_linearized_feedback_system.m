
%% Inputs

A = [
    0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0;
];
B = [
    0 0 0 0;
    1 0 0 0;
    0 0 0 0;
    0 1 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0;
    0 0 0 1;
];

Q = diag([50 50 5 1 5 1 5 1]);
R = diag([0.002 0.002 0.002 0.002]);

%% LQR
LQR_K = lqr(A, B, Q, R);

