%% Single Scissor Pair Gimbal Control System
clear; clc; close;

% Define symbols / conventions
% alpha = Gimbal A angle, beta = gimbal B angle, theta = (alpha + beta)/2
% (common angle)
% angle, gamma = alpha - Beta

% Define angles as a function of time
syms alpha(t) beta(t) theta(t) gamma(t)
assume([alpha(t) beta(t) theta(t) gamma(t)],'real')

theta = (alpha + beta)/2;
gamma = alpha - beta;


% Define inertia
% Ixx = rotation axis of flywheel, Iyy, Izz - Gimbal Axis

m = 0.25; % kg 
radius = 0.025; % meters

Ixx = 1/2*m*radius^2;
Iyy = Ixx/2;
Izz = Iyy;
inertiaMatrix = diag([Ixx, Iyy, Izz]);

% Define angular velocities and transform to common/differential coordinates
alphaDot = diff(alpha, t);
betaDot = diff(beta, t);
thetaDot = (alphaDot + betaDot)/2;
gammaDot = alphaDot - betaDot;

% Define Damping
B = 0.005;  % N m s/rad

J = Izz;
V_dc   = 20;
w_nl   = 513*2*pi/60;          % rad/s (use 540 RPM mid-range to see spread)
I_nl   = 0.07;                 % A
R_LL   = 11;                   % ohm, line-to-line (verify)
R_ph   = R_LL/2;               % ohm, Y connection

Kv_LL  = (V_dc - I_nl*R_LL)/w_nl;   % ~0.358 V·s/rad, line-to-line DC-equivalent
lam    = Kv_LL/sqrt(3);             % ~0.207 V·s/rad, peak phase flux linkage * pole pairs
Kt_q   = 1.5*lam;                   % ~0.31 N·m per A of q-axis current

Km     = Kt_q/R_ph;                 % ~0.056 N·m per volt of Uq
B_bemf = Kt_q*lam/R_ph;             % ~0.0117 N·m·s/rad
B_eff  = B + B_bemf;                % ~0.017

P = Km/(J*s + B_eff) * (1/s);

Kfoc = 1; % FOC amplifier gain
Ke = 1; % encoder gain

%% Building Model

Ti = 0.001;


% Plant
s = tf('s');
P = (Kfoc*Km/(J*s+B))*(1/s)*Ke;
Pz = c2d(P,Ti,'zoh');
bode(Pz);
Cz = pidtune(Pz,'PID',100*2*pi);
Kp = Cz.Kp; Ki = Cz.Ki; Kd = Cz.Kd;
C = Kp + Ki/s + Kd*s;
Cz2 = c2d(C,Ti,'tustin');
%% Closed-Loop Discrete Model
Lz = Pz*Cz2;
Tz = feedback(Lz,1);

margin(Lz)

BW = bandwidth(Tz)/(2*pi)
%% Generate Trajectory

% General acceleration and time step
A = pi/0.05;
D = -A;
% Position
alpha0 = 0; beta0 = 0;
alpha1 = pi; beta1 = pi;

% Velocity
alphaDot0 = 0; betaDot0 = 0;
alphaDotc1 = pi; betaDotc1 = pi;
alphaDot1 = 0; betaDot1 = 0;

% Acceleration
alphaDDot0 = A; betaDDot0 = A;
alphaDDot1 = -A; betaDDot1 = -A;

lalpha = alpha1-alpha0;
lbeta = beta1-beta0;

T2 = lalpha/alphaDotc1 - (1/(2*A) - 1/(2*D))*alphaDotc1;
if(T2 < 0)
    alphaDotcm = sqrt(lalpha/(1/2/A-1/2/D)); % Simplified form 
else
    alphaDotcm = alphaDotc1;
end

T1 = alphaDotcm/A;
N1 = ceil(T1/Ti);
T1_=N1*Ti;
fcn = lalpha/T1_;
A_=fcn/T1_;     D_=-fcn/T1_;


% Make trajectory
for i=1:N1
    k = i;
    L(k) = 0.5*A_*(i*Ti)^2;
    f(k) = A_*i*Ti;
    a(k) = A_;
    t(k) = k*Ti;
end

L1e = L(N1);

for i=1:N1
    k=i+N1;
    L(k) = L1e + fcn*i*Ti + 0.5*D_*(i*Ti)^2;
    f(k) = fcn + D_*i*Ti;
    a(k) = D_;
    t(k) = k*Ti;    
end

t_tot = t;

alpha_tot = alpha0 + L*(alpha1-alpha0)/lalpha;
beta_tot = beta0 + L*(beta1-beta0)/lbeta;
alphaDot_tot = f*(alpha1-alpha0)/lalpha;
betaDot_tot = f*(beta1-beta0)/lbeta;
alphaDDot_tot = a*(alpha1-alpha0)/lalpha;
betaDDot_tot = a*(beta1-beta0)/lbeta;

% Simulate the discrete closed-loop response to the generated trajectory

%Results = sim('Scissor_Pair_Control_Draft');