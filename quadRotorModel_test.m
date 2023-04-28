clc, clear , close all
%% Inputs
g = 9.81;
Mass = 1.605;               % Mass Kg
aLen = 0.24;           % arm length meter
Ixy=0;  Iyz=0;  Ixz=0;      % inertia Kg m^2
Ixx = 0.0017; Iyy = 0.0166; Izz = 0.0218;
Inertia = [...
         Ixx, -Ixy, -Ixz;...
        -Ixy,  Iyy, -Iyz;...
        -Ixz, -Iyz,  Izz];

K_th = 2.7778e-08;
K_tor = 3.6559e-07;

K_matrix = [...
    K_th K_th K_th K_th;
    K_th*aLen 0 -K_th*aLen 0;
    0 -K_th*aLen 0 K_th*aLen;
    K_tor -K_tor K_tor -K_tor
];
K_matrix_inv = inv(K_matrix);
u_0 = [11903.8 11903.9 11904.0 11903.9]';
% [u; v; w; p; q; r; phi; theta; epsi; xe0; ye0; ze0]
initialState = [0; 0; 0; 10*pi/180; 5*pi/180; 0; 0; 0; 0; 0; 0; -10]; %% example of free fall0
gravity_Iframe = [0;0;Mass*g];


%% simulation options
final_time = 20;
Ts = 0.01;
timeSteps = final_time/Ts;
time_V = 0:Ts:final_time;
result = NaN(12, timeSteps);
result(:, 1) = initialState;

%% Desired values
T_d = 2*Mass*g;
L = 0;
M = 0;
N = 0;
z_d = -10;

%% inital rigid body solver
rigidBodySolver = RigidBodySolver(Mass, Inertia, Ts, g);

%% Solving

dForces = [0;0;0];
dMoments = [0;0;0];

for i =1:timeSteps
    
    result(:, i+1) = rigidBodySolver.nextStep(result(:, i), dForces, dMoments);
    
    omega_d = [5*pi/180;0;0]; % phi theta epsi
    
    
    % Angluar velocitis [p, q, r] from desired eular angles [phi, theta,
    % epsi]
    wr_d = 5*rigidBodySolver.currentWr2omegaDot^-1*(omega_d - result(7:9, i+1) );
    
    
    % desired Torque from anglur rates error
    Moments_d = cross(wr_d, Inertia*wr_d) + 2*Inertia*(wr_d - result(4:6, i+1));
    
    % desired Thrust from altitute error
    T_d = (Mass*g) + .5*tanh(result(12, i+1)-z_d);  % Ok but it will occilate!
    
    % Requred motor speed from desired torque and moements
    d0 = K_matrix_inv*[T_d Moments_d']';
    u_0 = sign(d0).*sqrt(abs(d0));
    
    thrust_moment = K_matrix*u_0.^2;
    dForces = [0;0;-thrust_moment(1)];
    dMoments = thrust_moment(2:4);
    % then: AirFrame generate the Forces and moments    
end

%% Simulink
%sim('quadRotor.slx');

%% Plotting
u = result(1,:);        p = result(4,:)*180/pi;
v = result(2,:);        q = result(5,:)*180/pi;
w = result(3,:);        r = result(6,:)*180/pi;

phi = result(7,:)*180/pi;      x = result(10,:);
theta = result(8,:)*180/pi;    y = result(11,:);
psi = result(9,:)*180/pi;      z = result(12,:);

Vtotal = (u.^2 + v.^2 + w.^2).^(0.5);

beta_deg = asind(v./Vtotal); alpha_deg = atand(w./u);


%fig1 = figure;
plot3(x,y,z);
%hold on;
%plot3(xyz.Data(:, 1), xyz.Data(:, 2), xyz.Data(:, 3));
zlim([-15 15]);
fig1.CurrentAxes.ZDir = 'Reverse';
fig1.CurrentAxes.YDir = 'Reverse';
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
%legend('my solver', 'simulink');

figure;
subplot(3, 1, 1);
plot(time_V, phi); title('phi in x direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 1));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 2);
plot(time_V, theta); title('theta in y direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 2));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 3);
plot(time_V, psi); title('epsi in z direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 3));
grid on; %legend('my solver', 'simulink');


figure;
subplot(3, 1, 1);
plot(time_V, p); title('p in x direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 1));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 2);
plot(time_V, q); title('q in y direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 2));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 3);
plot(time_V, r); title('r in z direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 3));
grid on; %legend('my solver', 'simulink');


figure;
subplot(3, 1, 1);
plot(time_V, x); title('x direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 1));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 2);
plot(time_V, -y); title('y direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 2));
grid on; %legend('my solver', 'simulink');
subplot(3, 1, 3);
plot(time_V, -z); title('z direction');
%hold on;
%plot(phi_theta_epsi.Data(:, 3));
grid on; %legend('my solver', 'simulink');
