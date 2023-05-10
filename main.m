clc, clear , close all
%% Inputs & Configurations
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

% [u; v; w; p; q; r; phi; theta; epsi; xe0; ye0; ze0]
initialState = [0; 0; 0; 5*pi/180; 5*pi/180; 0; 10*pi/180; 0*pi/180; 10*pi/180; 0; 0; -5]; %% example of free fall0
gravity_Iframe = [0;0;Mass*g];


%% simulation options
final_time = 10;
Ts = 0.01;
timeSteps = final_time/Ts;
time_V = 0:Ts:final_time;
result = NaN(12, timeSteps);
dForces = NaN(timeSteps, 3);
dMoments = NaN(timeSteps, 3);
result(:, 1) = initialState;

%% Desired values - ref signal
z_d = -10;
omega_d = [0;0;0]; % phi theta epsi

%% inital rigid body solver
rigidBodySolver = RigidBodySolver(Mass, Inertia, Ts, g);

dForces(1, :) = [0 0 -Mass*g];
dMoments(1, :) = [0 0 0];

%% LQR Controller Gains
run('lqr_linearized_feedback_system');

%% Simulation

for i =1:timeSteps
    
    % Rigid Body Solver
    result(:, i+1) = rigidBodySolver.nextStep(result(:, i), dForces(i, :)', dMoments(i, :)');
    
    % Controller 1
    u_0 = controller_lqr(Inertia, Mass, g, K_matrix_inv, rigidBodySolver, omega_d, z_d, result(:, i+1), LQR_K);    
    
    % Motor Rotional Speed to Force and Moments
    thrust_moment = K_matrix*u_0.^2;
    dForces(i+1, :) = [0 0 -thrust_moment(1)];
    dMoments(i+1, :) = thrust_moment(2:4)';
end


%% Plotting
u = result(1,:);        p = result(4,:)*180/pi;
v = result(2,:);        q = result(5,:)*180/pi;
w = result(3,:);        r = result(6,:)*180/pi;

phi = result(7,:)*180/pi;      x = result(10,:);
theta = result(8,:)*180/pi;    y = result(11,:);
psi = result(9,:)*180/pi;      z = result(12,:);

Vtotal = (u.^2 + v.^2 + w.^2).^(0.5);

beta_deg = asind(v./Vtotal); alpha_deg = atand(w./u);

plot3(x,-y,-z); title('Position');
text(x(1),-y(1), -z(1),'\leftarrow Start Position', 'Color', 'red', 'FontSize', 10);
text(x(end),-y(end), -z(end),'\leftarrow Final Position')
zlim([-15 15]);xlim([-15 15]); ylim([-15 15]);
grid on;
xlabel('x'); ylabel('y'); zlabel('z');


figure;
subplot(3, 1, 1);
plot(time_V, phi); title('Angles (degs)');
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('$\phi$', 'Interpreter', 'latex');
hold on;
yline(omega_d(1)*180/pi, 'Color', 'red');
grid on; 
subplot(3, 1, 2);
plot(time_V, theta); 
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('$\theta$', 'Interpreter', 'latex');
hold on;
yline(omega_d(2)*180/pi, 'Color', 'red');
grid on; 
subplot(3, 1, 3);
plot(time_V, psi);
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('$\psi$', 'Interpreter', 'latex');
hold on;
yline(omega_d(3)*180/pi, 'Color', 'red');
grid on;


figure;
subplot(3, 1, 1);
plot(time_V, p); title('Angler Velocities');
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('p');
grid on; 
subplot(3, 1, 2);
plot(time_V, q);
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('q');
grid on; 
subplot(3, 1, 3);
plot(time_V, r); 
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('r');
grid on;


figure;
subplot(3, 1, 1);
plot(time_V, x); title('Position');
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('x');

grid on;
subplot(3, 1, 2);
plot(time_V, -y);
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('y');

grid on;
subplot(3, 1, 3);
plot(time_V, -z);
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('z');
grid on;


%% Plotting Input Forces & Moments
figure;
subplot(1, 2, 1); 
plot(time_V, dForces(:, 3));title('Force In Body Frame');
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('Z-Force');
legend({'$F_{z}$'},'Interpreter', 'latex', 'FontSize', 14);

subplot(1, 2, 2);
plot(time_V, dMoments(:, 1)); hold on;
plot(time_V, dMoments(:, 2)); hold on;
plot(time_V, dMoments(:, 3)); hold on;
xlabel('time (sec)', 'Interpreter', 'latex');
ylabel('Moments'); title('Moments');
legend({'$\tau_{\phi}$', '$\tau_{\theta}$', '$\tau_{\psi}$'},'Interpreter', 'latex', 'FontSize', 14);


function u = controller_lqr(Inertia, Mass, g, K_matrix_inv, rigidBodySolver, omega_d, z_d, state, LQR_K)

    % LQR Gains
    K_z = LQR_K(1, 1); K_z_dot = LQR_K(1, 2);
    K_phi = LQR_K(2, 3); K_phi_dot = LQR_K(2, 4);
    K_theta = LQR_K(3, 5); K_theta_dot = LQR_K(3, 6);
    K_epsi = LQR_K(4, 7); K_epsi_dot = LQR_K(4, 8);

    % Inertia Values
    Ixx = Inertia(1, 1); Iyy = Inertia(2, 2); Izz = Inertia(3, 3);
    a1 = (Iyy - Izz) / Ixx;
    a2 = (Izz - Ixx) / Iyy;
    a3 = (Ixx - Iyy) / Izz;
    b1 = 1/Ixx; b2 = 1/Iyy; b3 = 1/Izz;
    
    % Calculate Angler Rates & velocity in Inertial Frame
    phi_theta_epsi_dot = rigidBodySolver.currentWr2omegaDot * state(4:6);
    x_y_z_dot = rigidBodySolver.inertial2Body^-1*state(1:3);
    
    % Error
    angles_err = state(7:9) - omega_d;
    wr_err = phi_theta_epsi_dot;
    z_err = z_d - state(12);
    z_dot_err = -x_y_z_dot(3);
    
    % Control Input To Feedbak Linearized System
    v1 = -K_z*tanh(z_err) - K_z_dot*tanh(z_dot_err);
    v2 = -K_phi*angles_err(1) - K_phi_dot*wr_err(1);
    v3 = -K_theta*angles_err(2) - K_theta_dot*wr_err(2);
    v4 = -K_epsi*angles_err(3) - K_epsi_dot*wr_err(3);
    
    % Force & Moments as an Input to Original System
    v1_sat = saturation([-4 4]);
    u1 = Mass*(v1_sat.evaluate(v1) + g);
    u2 = 1/b1*(v2 - a1*phi_theta_epsi_dot(2)*phi_theta_epsi_dot(3));
    u3 = 1/b2*(v3 - a2*phi_theta_epsi_dot(1)*phi_theta_epsi_dot(3));
    u4 = 1/b3*(v4 - a3*phi_theta_epsi_dot(1)*phi_theta_epsi_dot(2));
       
    % Forces & Moments to Motor Rotional Speeds
    d0 = K_matrix_inv*[u1 u2 u3 u4]';
    u = sign(d0).*sqrt(abs(d0));
end
