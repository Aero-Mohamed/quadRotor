clc, clear vars, close all
%% Inputs
g = 9.81;
Mass = 1.605;               % Mass Kg
armLength = 0.24;           % arm length meter
Ixy=0;  Iyz=0;  Ixz=0;      % inertia Kg m^2
Ixx = 0.0017; Iyy = 0.0166; Izz = 0.0218;
Inertia = [...
         Ixx, -Ixy, -Ixz;...
        -Ixy,  Iyy, -Iyz;...
        -Ixz, -Iyz,  Izz];
% [u; v; w; p; q; r; phi; theta; epsi; xe0; ye0; ze0]
initialState = [1; 0; -50; 0; 0; 0; 0; 0; 0; 10; 10; 10]; %% example of free fall0



%% simulation options
final_time = 20;
dt = 0.01;
timeSteps = final_time/dt;
result = NaN(12, timeSteps);
result(:, 1) = initialState;

%% inital rigid body solver
rigidBodySolver = RigidBodySolver(Mass, Inertia, dt, g);

%% Solving

dForces = [0 ; 0; 0];
dMoments = [0 ; 0; 0];

for i =1:timeSteps
    
    result(:, i+1) = rigidBodySolver.nextStep(result(:, i), dForces, dMoments);
    
    % then: AirFrame generate the Forces and moments    
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


fig1 = figure;
plot3(x,y,z);
fig1.CurrentAxes.ZDir = 'Reverse';
fig1.CurrentAxes.YDir = 'Reverse';
grid on;
xlabel('x'); ylabel('y'); zlabel('z');

figure;
subplot(3, 1, 1);
plot(u); title('u in x direction');
grid on; axis equal
subplot(3, 1, 2);
plot(v); title('v in x direction');
grid on; axis equal
subplot(3, 1, 3);
plot(w); title('w in x direction');
grid on; axis equal

