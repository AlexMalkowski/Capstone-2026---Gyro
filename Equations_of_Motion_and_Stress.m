% Flywheel 
% Define flywheel geometry and material properties
thickness = 0.01;       % m
density_brass = 8400;   % kg/m^3
density_steel = 7800;   % kg/m^3
nu_steel = 0.3;
density = density_steel;% kg/m^3

% Flywheel Mass
rest_of_system_mass = 0.5; % kg

mass_budget = 1.5; % kg

% Flywheel configs
num_of_flywheel = 2;
% num_of_flywheel = 4;

% Determine allowable mass and radius for each flywheel
flywheelMass = (mass_budget - rest_of_system_mass) / num_of_flywheel;
flywheelVolume = flywheelMass / density;
flywheelRadius = sqrt(flywheelVolume / (pi * thickness));

% Solid disc inertia
flywheelInertia = 0.5 * flywheelMass * flywheelRadius^2;


% Paramaters
angular_speed_rpm = 10000; % rpm
angular_speed = angular_speed_rpm*2*pi/60;

% Stress Eqaution
motor_shaft_radius = 0.005; % m
nu = nu_steel;

rotor_max_stress = ((3 + nu)/4)*density*angular_speed^2 ... 
 *flywheelRadius^2*(1+(1-nu)/(3+nu)*(motor_shaft_radius/flywheelRadius^2));

rotor_max_stress_MPa = rotor_max_stress/(10^6)

% Equations of Motion
angular_momentum = flywheelInertia * angular_speed

% Iz*omega_dot + C*omega = Torque(t)
rotor_damping = 0;


% Equations for Gimbal

% Iy * theta_ddot + C*theta_dot + k*theta = Torque(t)

syms theta(t)
thetaDot = diff(theta, t);
thetaDDot = diff(theta, t, 2);

Iy = flywheelMass*flywheelRadius^2/4;

gimbal_damping = 0;
gimbal_stiffness = 0;
gimbal_torque = 1; %Nm


gimbalEquation = Iy*thetaDDot + gimbal_damping*thetaDot + ...
    gimbal_stiffness*theta == gimbal_torque;
