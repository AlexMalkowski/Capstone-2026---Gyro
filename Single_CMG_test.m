% Clear workspace
close all;
clear;
clc;

% Define Parameters
gamma_dot = 30; % Gimbal speed [rad/s]
u = 20*1000;    % Flywheel speed [rpm]
d = 50/1000;    % Flywheel diameter [meters]
w = 20/1000;    % Flywheel width [meters]
rho = 7500;     % Flywheel density [kg/m^3]
% Convert flywheel speed to angular velocity
omega = u * 2 * pi / 60;

% Calculate flywheel geometry and mass properties
radius = d / 2;
flywheelVolume = pi * radius^2 * w;
flywheelMass = rho * flywheelVolume;
flywheelInertia = 0.5 * flywheelMass * radius^2;

% Define Flywheel & Gimbal Orientations
g = [1; 0; 0];  % Gimbal axis (along x)
h_0 = [0; 0; 1];  % Unperturbed flywheel axis (along z)

gamma = linspace(0, 2*pi, 100);

function h = rotate_h(g,h_0,angle)
g = g / norm(g);
    h = h_0 * cos(angle) + cross(g, h_0) * sin(angle);
end 

torque = zeros(3, length(gamma));
angularMomentum = zeros(3, length(gamma));

for i = 1:length(gamma)
    h = rotate_h(g, h_0, gamma(i));
    angularMomentum(:,i) = flywheelInertia * omega * h;
    torque(:,i) = gamma_dot * cross(g, angularMomentum(:,i));
end

torque_y = torque(2,:);
torque_z = torque(3,:);
figure;
plot(gamma, torque_y, 'LineWidth', 1.5);
hold on;
plot(gamma, torque_z, 'LineWidth', 1.5);
grid on;
xlabel('Gimbal Angle [rad]');
ylabel('Torque [N·m]');
legend('Torque Y', 'Torque Z', 'Location', 'best');

fprintf('Max torque in z is %d Nm \n', max(torque_z))
fprintf('Max torque in y is %d Nm \n', max(torque_y))
fprintf('Flywheel mass is %d kg \n', flywheelMass)

%% Animation

% Animation Setup in the YZ Plane
figure('Color', 'white');
grid on;
axis equal;
xlabel('Torque Y [N\cdot{}m]');
ylabel('Torque Z [N\cdot{}m]');
title('Torque Vector Rotation in the YZ Plane');

% Set axis limits dynamically with a margin
max_t = max(sqrt(torque_y.^2 + torque_z.^2)) * 1.2;
if max_t == 0, max_t = 1; end
xlim([-max_t, max_t]);
ylim([-max_t, max_t]);

hold on;

% Plot the full trajectory path as a faint reference
plot(torque_y, torque_z, 'Color', [0.7, 0.7, 0.7], 'LineStyle', '--', 'LineWidth', 1.5);

% Initialize the animated vector (quiver) and tip marker
h_vector = quiver(0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
h_point = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

% Infinite Animation Loop
while ishghandle(h_vector) % Keeps running until the figure window is closed
    for i = 1:length(gamma)
        if ~ishghandle(h_vector), break; end % Exit if figure was closed mid-loop
        
        set(h_vector, 'UData', torque_y(i), 'VData', torque_z(i));
        set(h_point, 'XData', torque_y(i), 'YData', torque_z(i));
        
        drawnow;
        pause(0.03); % Control animation playback speed
    end
end