% Clear workspace
close all;
clear;
clc;

% Define Parameters 
u = input('Enter flywheel speed [rpm]: ');
d = input('Enter flywheel diameter [mm]: ');
d = d/1000;
w = input('Enter flywheel width [mm]: ');
w = w/1000;
beta = input('Enter skew angle [degrees]: ');
n = input('Enter number of gimbals: ');
max_gimbal_rate = input('Enter max gimbal rate [rad/s]: ');
rho = 7500; % Flywheel density [kg/m^3]
omega = u * 2 * pi / 60;
beta = beta * pi / 180; % Pyramid skew angle [rad]

% Flywheel Inertia & Momentum
radius = d / 2;
flywheelInertia = 0.5 * (rho * pi * radius^2 * w) * radius^2;
H = flywheelInertia * omega;

g_axes = zeros(3, n);
h_axes = zeros(3, n);

% Gimbal & Spin Axis Geometries
for i = 1:n
    psi = (i - 1) * (2 * pi / n);
    g_axes(:, i) = [sin(beta) * cos(psi); 
                    sin(beta) * sin(psi); 
                    cos(beta)];
    h_axes(:,i) = cross([cos(psi); sin(psi); 0],[0,0,1]);
end

% Monte Carlo Torque Envelope Sampling
num_samples = 500000;          % Number of random configurations to test

torque_cloud = zeros(3, num_samples);

for k = 1:num_samples
    % Randomize gimbal angles uniformly between -pi and pi
    gammas = 2 * pi * rand(n, 1) - pi;
    
    % Randomize gimbal rates within the velocity limits (-max to max)
    gamma_dot = 2 * max_gimbal_rate * rand(n, 1) - max_gimbal_rate;
    
    % Build Jacobian and compute instantaneous torque vector
    J = build_jacobian(g_axes, h_axes, gammas, H);
    torque_cloud(:, k) = J * gamma_dot;
end

% Plot the 3D Torque Envelope
figure('Color', 'w');
% Plot scatter point cloud of reachable torque vectors
scatter3(torque_cloud(1, :), torque_cloud(2, :), torque_cloud(3, :), ...
         5, vector_magnitude(torque_cloud), 'filled');
colormap(jet);
colorbar;
grid on;
axis equal;
title(n+"-CMG 3D Torque Capability Envelope");
xlabel('Torque X (N\cdot m)');
ylabel('Torque Y (N\cdot m)');
zlabel('Torque Z (N\cdot m)');
view(45, 30);

% Extract component-wise torque limits and flywheel mass
torque_x = torque_cloud(1, :);
torque_y = torque_cloud(2, :);
torque_z = torque_cloud(3, :);
flywheelMass = rho * pi * radius^2 * w;
fprintf('Max torque in x is %.2f Nm\n', max(torque_x));
fprintf('Max torque in y is %.2f Nm\n', max(torque_y));
fprintf('Max torque in z is %.2f Nm\n', max(torque_z));
fprintf('Mass of all flywheels is %d kg \n', flywheelMass*n)

% Helper Functions
function h = rotate_h(g, h_0, angle)
    g = g / norm(g);
    h = h_0 * cos(angle) + cross(g, h_0) * sin(angle);
end 

function J = build_jacobian(g_axes, h_axes, gammas, H)
    num_cmgs = size(g_axes, 2);
    J = zeros(3, num_cmgs);
    for i = 1:num_cmgs
        h = rotate_h(g_axes(:, i), h_axes(:, i), gammas(i));
        J(:, i) = H * cross(g_axes(:, i), h);
    end 
end

function mags = vector_magnitude(vecs)
    mags = sqrt(sum(vecs.^2, 1));
end