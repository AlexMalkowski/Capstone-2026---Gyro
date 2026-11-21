%% Scissor-Pair Gimbal: Two-Motor Tracking Simulation (benchtop)
% Conventions (from the original script):
%   alpha = gimbal A angle, beta = gimbal B angle           [rad]
%   theta = (alpha + beta)/2   common angle
%   gamma =  alpha - beta      differential (scissor) angle
%   =>  alpha = theta + gamma/2,   beta = theta - gamma/2
%
% Plant (each motor, independent SISO loop):
%   SimpleFOC voltage-torque mode: command is q-axis voltage Uq [V]
%   J*w' = Km*Uq - B_eff*w ,  angle' = w      (resistive-only winding model)
%   B_eff = B_mech + back-EMF damping.  No gyroscopic coupling between
%   gimbals, no head-motion disturbance, no winding inductance.
%
% Structure:
%   1. HARDWARE + INERTIA INPUTS   (EDIT as parameters are measured)
%   2. TRAJECTORY INPUTS           (EDIT)
%   3. CONTROLLER TUNING           (EDIT)
%   4. ---- no edits below, except the local functions are readable ----
%
% Requires: Control System Toolbox (pidtune, c2d, margin).
clear; clc; close all;

%% ======================================================================
%  1. HARDWARE + INERTIA INPUTS  (EDIT)
%  ======================================================================
% --- 1a. Motor and driver (datasheet no-load row, SimpleFOC Mini / DRV8313)
hw.V_dc         = 20;      % supply voltage [V]
hw.omega_nl_rpm = 513;     % no-load speed at V_dc [rpm]  (datasheet: 513-567)
hw.I_nl         = 0.07;    % no-load current [A]
hw.R_LL         = 11;      % line-to-line resistance [ohm]  (VERIFY with multimeter)
hw.N_gear       = 1;       % gear ratio motor:gimbal (1 = direct drive)
hw.Uq_limit     = 6;       % Uq clamp used in firmware [V] (hard max = V_dc/sqrt(3))
hw.I_limit      = 2.5;     % driver peak current to compare against [A] (VERIFY vs board)
hw.B_mech       = 0.005;   % viscous friction at gimbal axis [N m s/rad] (GUESS: measure by coast-down)

% --- 1b. Digital / sensing
hw.Ts            = 1e-3;   % control period [s]
hw.delay_samples = 1;      % compute + sensor-read delay [samples]
hw.enc_cpr       = 4096;   % encoder counts per rev (12-bit, e.g. AS5600)

% --- 1c. Inertia about the GIMBAL axis [kg m^2]
inertia.m_flywheel   = 0.25;    % flywheel mass [kg]
inertia.R_flywheel   = 0.025;   % flywheel radius [m]
inertia.J_other_load = 2e-5;    % PLACEHOLDER: spin motor + housing + frame about gimbal axis (from CAD)
inertia.J_rotor      = 5e-6;    % PLACEHOLDER: gimbal motor rotor inertia (multiplied by N_gear^2)

% --- 1d. Parameter-uncertainty cases (all simulated with the same controller)
%   fluxX = scale on flux linkage (Kt and back-EMF constant move together;
%           use it to cover the +-5% speed and resistance spread)
%   JX    = scale on total inertia.   X = A or B motor.
cases = struct( ...
    'name', {'Nominal','Low gain / high J','High gain / low J','A-B mismatch'}, ...
    'fluxA',{1,    0.92, 1.08, 1.05}, ...
    'JA',   {1,    1.20, 0.80, 1.10}, ...
    'fluxB',{1,    0.92, 1.08, 0.95}, ...
    'JB',   {1,    1.20, 0.80, 0.90});

%% ======================================================================
%  2. TRAJECTORY INPUTS  (EDIT)
%  Sequence: theta 0 -> theta_target (gamma = 0), dwell,
%            gamma 0 -> gamma_target at gamma_dot, dwell,
%            optionally gamma back to 0, dwell.
%  Each move is a trapezoidal velocity profile (triangular if the move is
%  too short to reach the speed limit), optionally jerk-smoothed.
%  ======================================================================
traj.theta_target = deg2rad(45);    % common-mode target [rad]
traj.theta_vmax   = deg2rad(90);    % max theta speed [rad/s]
traj.theta_amax   = deg2rad(720);   % max theta accel [rad/s^2]

traj.gamma_target = deg2rad(60);    % differential target [rad]
traj.gamma_dot    = deg2rad(120);   % gamma cruise speed [rad/s]
traj.gamma_amax   = deg2rad(1200);  % max gamma accel [rad/s^2]
% Spec check: 60 deg in 50 ms is a triangular move with
%   gamma_amax = 4*0.9996/0.05^2 ~ 1600 rad/s^2 and peak ~ 42 rad/s -- try it
%   once Uq_limit and inertia are real; the feasibility report will flag it.

traj.return_gamma = true;           % also drive gamma back to 0?

traj.dwell_after_theta = 0.5;       % [s]
traj.dwell_after_gamma = 0.5;       % [s]  (after each gamma move)
traj.jerk_time         = 0.02;      % jerk-limiting smoothing time [s]; 0 = pure trapezoid

% --- Accuracy-metric settings
analysis.settle_band_deg = 0.5;     % settling band [deg]
analysis.lag_search_ms   = 100;     % max delay searched for best-fit lag [ms]

%% ======================================================================
%  3. CONTROLLER TUNING  (EDIT -- re-tune whenever section 1 changes)
%  Controller: discrete PID with filtered derivative (on measurement),
%  velocity/accel feedforward from the nominal model, Uq clamp with
%  conditional-integration anti-windup.  Gains are in V/rad.
%  ======================================================================
ctrl.method       = 'pidtune';  % 'pidtune' (auto) or 'manual'
ctrl.bw_Hz        = 15;         % target loop crossover [Hz]          (pidtune)
ctrl.phase_margin = 60;         % target phase margin [deg]           (pidtune)
ctrl.manual = struct('Kp',30,'Ki',100,'Kd',0.1,'Tf',1e-3);  % used if 'manual' (starting guess)
ctrl.ff_vel = 1;                % 0..1 scale on velocity feedforward (B_eff*v/Km)
ctrl.ff_acc = 1;                % 0..1 scale on accel feedforward   (J*a/Km)

%% ======================================================================
%  4. DERIVED PARAMETERS  (no edits)
%  ======================================================================
hw.R_ph      = hw.R_LL/2;                    % equivalent star phase resistance
hw.omega_nl  = hw.omega_nl_rpm*2*pi/60;      % [rad/s]
hw.Kv_LL     = (hw.V_dc - hw.I_nl*hw.R_LL)/hw.omega_nl;  % line-line back-EMF const [V s/rad]
hw.lam       = hw.Kv_LL/sqrt(3);             % flux linkage x pole pairs [V s/rad]
hw.Kt_q      = 1.5*hw.lam;                   % torque per q-current [N m/A]
hw.Uq_hard   = hw.V_dc/sqrt(3);              % SVPWM linear limit [V]
hw.Uq_lim    = min(hw.Uq_limit, hw.Uq_hard); % effective clamp [V]

% Inertia: flywheel is a solid disc; Ixx = spin axis, Iyy = Izz = gimbal axes
Ixx = 0.5*inertia.m_flywheel*inertia.R_flywheel^2;
Iyy = Ixx/2;
Izz = Iyy;
inertiaMatrix = diag([Ixx, Iyy, Izz]); %#ok<NASGU>  kept for later 3-D / coupled model
inertia.J_flywheel = Izz;
J_total = inertia.J_flywheel + inertia.J_other_load + hw.N_gear^2*inertia.J_rotor;

pl_nom = make_plant(hw, J_total, 1, 1);

fprintf('--- Plant (nominal) ---\n');
fprintf('J_total = %.3g kg m^2 | Km = %.4f N m/V | B_eff = %.4f N m s/rad\n', pl_nom.J, pl_nom.Km, pl_nom.B_eff);
fprintf('Mechanical pole B_eff/J = %.0f rad/s | Kt = %.3f N m/A | Uq clamp = %.1f V\n', ...
    pl_nom.B_eff/pl_nom.J, hw.Kt_q, hw.Uq_lim);
fprintf('Max speed at clamp (no load) = %.1f rad/s (%.0f deg/s)\n\n', ...
    hw.Uq_lim/pl_nom.lamN, rad2deg(hw.Uq_lim/pl_nom.lamN));

%% ======================================================================
%  5. TRAJECTORY GENERATION  (independent of the plant)
%  ======================================================================
ref = build_trajectory(traj, hw.Ts);

% Feedforward voltage from the NOMINAL model, per gimbal
ffA = (ctrl.ff_acc*pl_nom.J*ref.alphaDDot + ctrl.ff_vel*pl_nom.B_eff*ref.alphaDot)/pl_nom.Km;
ffB = (ctrl.ff_acc*pl_nom.J*ref.betaDDot  + ctrl.ff_vel*pl_nom.B_eff*ref.betaDot )/pl_nom.Km;

% Feasibility report
peakFF   = max(abs([ffA; ffB]));
peakSpd  = max(abs([ref.alphaDot; ref.betaDot]));
fprintf('--- Trajectory feasibility ---\n');
fprintf('Duration %.2f s | peak gimbal speed %.1f rad/s (%.0f deg/s)\n', ref.t(end), peakSpd, rad2deg(peakSpd));
fprintf('Peak feedforward Uq = %.2f V vs clamp %.2f V (%.0f%%)\n\n', peakFF, hw.Uq_lim, 100*peakFF/hw.Uq_lim);
if peakFF > hw.Uq_lim
    warning('Trajectory needs more voltage than the Uq clamp allows: reduce speed/accel, raise the clamp, or revisit J / gearing.');
end

%% ======================================================================
%  6. TUNING + LINEAR ANALYSIS
%  ======================================================================
s = tf('s');
P_tune = pl_nom.Km/(s*(pl_nom.J*s + pl_nom.B_eff));
P_tune.InputDelay = (hw.delay_samples + 0.5)*hw.Ts;   % compute delay + ZOH (~half sample)

switch lower(ctrl.method)
    case 'pidtune'
        opts  = pidtuneOptions('PhaseMargin', ctrl.phase_margin);
        Cc    = pidtune(P_tune, 'PIDF', 2*pi*ctrl.bw_Hz, opts);
        gains = struct('Kp',Cc.Kp,'Ki',Cc.Ki,'Kd',Cc.Kd,'Tf',Cc.Tf);
    case 'manual'
        gains = ctrl.manual;
    otherwise
        error('ctrl.method must be ''pidtune'' or ''manual''.');
end
fprintf('--- Controller gains ---\n');
fprintf('Kp = %.4g V/rad | Ki = %.4g V/(rad s) | Kd = %.4g V s/rad | Tf = %.3g s\n\n', ...
    gains.Kp, gains.Ki, gains.Kd, gains.Tf);


nC = numel(cases);
for c = 1:nC
    plA(c) = make_plant(hw, J_total, cases(c).fluxA, cases(c).JA); %#ok<SAGROW>
    plB(c) = make_plant(hw, J_total, cases(c).fluxB, cases(c).JB); %#ok<SAGROW>
    la(c)  = linear_analysis(plA(c), gains, hw);                   %#ok<SAGROW>
end

%% ======================================================================
%  7. SIMULATION (both motors, every uncertainty case)
%  ======================================================================
for c = 1:nC
    simA = simulate_axis(ref.alpha, ffA, plA(c), gains, hw);
    simB = simulate_axis(ref.beta,  ffB, plB(c), gains, hw);
    res(c) = compute_metrics(ref, simA, simB, hw, traj, analysis); %#ok<SAGROW>
end

%% ======================================================================
%  8. RESULTS TABLES
%  ======================================================================
names = {cases.name}';
col   = @(f) [res.(f)]';
lin   = @(f) [la.(f)]';

T_loop = table(lin('PM_deg'), lin('GM_dB'), lin('Wcp_Hz'), lin('BW_Hz'), lin('Ms'), ...
    'VariableNames', {'PM_deg','GM_dB','Crossover_Hz','ClosedLoopBW_Hz','PeakSensitivity'}, ...
    'RowNames', names);
disp('--- Linear loop metrics (motor A, discrete, with delay) ---'); disp(T_loop);

T_track = table(col('rmsA_deg'), col('rmsB_deg'), col('rmsT_deg'), col('rmsG_deg'), ...
    col('maxT_deg'), col('maxG_deg'), col('finalT_deg'), col('finalG_deg'), ...
    'VariableNames', {'RMS_alpha','RMS_beta','RMS_theta','RMS_gamma', ...
                      'Max_theta','Max_gamma','Final_theta','Final_gamma'}, ...
    'RowNames', names);
disp('--- Tracking error [deg] ---'); disp(T_track);

T_dyn = table(col('lagT_ms'), col('lagG_ms'), col('lagG_cruise_ms'), ...
    col('overshootT_pct'), col('settleT_s'), col('settleG_s'), ...
    'VariableNames', {'Lag_theta_ms','Lag_gamma_ms','Lag_gamma_cruise_ms', ...
                      'Overshoot_theta_pct','Settle_theta_s','Settle_gamma_s'}, ...
    'RowNames', names);
disp('--- Lag / transient (lag > 0 = actual trails reference) ---'); disp(T_dyn);

T_eff = table(col('peakUq_V'), col('satFrac_pct'), col('peakIq_A'), col('peakSpeed_rpm'), ...
    'VariableNames', {'Peak_Uq_V','Saturated_pct','Peak_Iq_A','Peak_motor_rpm'}, ...
    'RowNames', names);
disp('--- Actuator demand ---'); disp(T_eff);

fprintf('Encoder resolution: %.3f deg/count (quantization floor on error).\n', 360/hw.enc_cpr);
fprintf('Driver current limit used for comparison: %.2f A | Uq clamp: %.2f V\n', hw.I_limit, hw.Uq_lim);

%% ======================================================================
%  9. PLOTS
%  ======================================================================
r2d = 180/pi;

figure('Name','Trajectory');
subplot(3,1,1); plot(ref.t, r2d*[ref.theta ref.gamma ref.alpha ref.beta]); grid on;
ylabel('Angle [deg]'); legend('\theta','\gamma','\alpha','\beta','Location','best');
title('Reference trajectory (plant-independent)');
subplot(3,1,2); plot(ref.t, r2d*[ref.thetaDot ref.gammaDot ref.alphaDot ref.betaDot]); grid on;
ylabel('Rate [deg/s]');
subplot(3,1,3); plot(ref.t, r2d*[ref.thetaDDot ref.gammaDDot ref.alphaDDot ref.betaDDot]); grid on;
ylabel('Accel [deg/s^2]'); xlabel('Time [s]');

m0 = res(1);   % nominal case detail
figure('Name','Nominal tracking');
subplot(4,1,1); plot(ref.t, r2d*ref.alpha, 'k--', ref.t, r2d*m0.yA, ...
    ref.t, r2d*ref.beta, 'k:', ref.t, r2d*m0.yB); grid on;
ylabel('[deg]'); legend('\alpha ref','\alpha','\beta ref','\beta','Location','best');
title(['Nominal case: ' cases(1).name]);
subplot(4,1,2); plot(ref.t, r2d*[m0.eA m0.eB]); grid on;
ylabel('e_\alpha, e_\beta [deg]'); legend('e_\alpha','e_\beta','Location','best');
subplot(4,1,3); plot(ref.t, r2d*[m0.eT m0.eG]); grid on;
ylabel('e_\theta, e_\gamma [deg]'); legend('e_\theta','e_\gamma','Location','best');
subplot(4,1,4); plot(ref.t, [m0.uA m0.uB]); hold on;
yline(hw.Uq_lim,'r--'); yline(-hw.Uq_lim,'r--'); hold off; grid on;
ylabel('U_q [V]'); xlabel('Time [s]'); legend('A','B','clamp','Location','best');

figure('Name','Loop analysis');
margin(la(1).Lz); grid on;

figure('Name','Closed-loop response');
subplot(2,1,1); step(la(1).Tz); grid on; title('Closed-loop step (motor A, nominal)');
subplot(2,1,2); bodemag(la(1).Sz, la(1).Tz); grid on; legend('S','T');

figure('Name','Error across cases');
subplot(2,1,1); hold on; for c = 1:nC, plot(ref.t, r2d*res(c).eT); end; hold off; grid on;
ylabel('e_\theta [deg]'); legend(names,'Location','best'); title('Error across uncertainty cases');
subplot(2,1,2); hold on; for c = 1:nC, plot(ref.t, r2d*res(c).eG); end; hold off; grid on;
ylabel('e_\gamma [deg]'); xlabel('Time [s]');

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function pl = make_plant(hw, J_total, flux_scale, J_scale)
% Voltage (Uq) -> gimbal angle plant parameters, referred to the gimbal axis.
    N         = hw.N_gear;
    pl.J      = J_total*J_scale;
    pl.Km     = N*(hw.Kt_q*flux_scale)/hw.R_ph;               % N m per V
    pl.lamN   = N*hw.lam*flux_scale;                          % back-EMF [V per rad/s of gimbal]
    pl.B_eff  = hw.B_mech + N*(hw.Kt_q*flux_scale)*pl.lamN/hw.R_ph;  % friction + back-EMF damping
end

function la = linear_analysis(pl, g, hw)
% Discrete loop analysis: ZOH plant + delay, PID with filtered derivative.
    Ts   = hw.Ts;
    sysd = c2d(ss([0 1; 0 -pl.B_eff/pl.J], [0; pl.Km/pl.J], [1 0], 0), Ts, 'zoh');
    Pz   = tf(sysd);
    Pz.InputDelay = hw.delay_samples;
    Pz   = absorbDelay(Pz);

    z  = tf('z', Ts);
    Cz = g.Kp + g.Ki*Ts*z/(z-1) + g.Kd/(g.Tf+Ts)*(z-1)/(z - g.Tf/(g.Tf+Ts));

    Lz = Cz*Pz;
    Tz = feedback(Lz, 1);
    Sz = feedback(1, Lz);
    [Gm, Pm, ~, Wcp] = margin(Lz);

    la.Lz = Lz; la.Tz = Tz; la.Sz = Sz;
    la.GM_dB  = 20*log10(Gm);
    la.PM_deg = Pm;
    la.Wcp_Hz = Wcp/(2*pi);
    la.BW_Hz  = bandwidth(Tz)/(2*pi);
    la.Ms     = norm(Sz, Inf);
end

function ref = build_trajectory(traj, Ts)
% Builds theta/gamma references (pos, vel, accel) then maps to alpha/beta.
    thP = 0;  gaP = 0;

    % Move 1: theta, gamma held at 0
    p  = point_to_point(traj.theta_target, traj.theta_vmax, traj.theta_amax, traj.jerk_time, Ts);
    n  = numel(p) - 1;
    thP = [thP; p(2:end)];           gaP = [gaP; zeros(n,1)];
    ev.theta_end = (numel(thP)-1)*Ts;
    [thP, gaP] = hold_both(thP, gaP, traj.dwell_after_theta, Ts);
    ev.gamma_start = (numel(thP)-1)*Ts;

    % Move 2: gamma to target at gamma_dot
    p  = point_to_point(traj.gamma_target, traj.gamma_dot, traj.gamma_amax, traj.jerk_time, Ts);
    n  = numel(p) - 1;
    thP = [thP; thP(end)*ones(n,1)]; gaP = [gaP; gaP(end) + p(2:end)];
    ev.gamma_end_last = (numel(thP)-1)*Ts;
    [thP, gaP] = hold_both(thP, gaP, traj.dwell_after_gamma, Ts);

    % Optional move 3: gamma back to 0
    if traj.return_gamma
        p  = point_to_point(-traj.gamma_target, traj.gamma_dot, traj.gamma_amax, traj.jerk_time, Ts);
        n  = numel(p) - 1;
        thP = [thP; thP(end)*ones(n,1)]; gaP = [gaP; gaP(end) + p(2:end)];
        ev.gamma_end_last = (numel(thP)-1)*Ts;
        [thP, gaP] = hold_both(thP, gaP, traj.dwell_after_gamma, Ts);
    end

    ref.t = (0:numel(thP)-1)'*Ts;
    [ref.theta, ref.thetaDot, ref.thetaDDot] = derive_signals(thP, Ts);
    [ref.gamma, ref.gammaDot, ref.gammaDDot] = derive_signals(gaP, Ts);

    % alpha = theta + gamma/2 ; beta = theta - gamma/2
    ref.alpha = ref.theta + ref.gamma/2;  ref.alphaDot = ref.thetaDot + ref.gammaDot/2;  ref.alphaDDot = ref.thetaDDot + ref.gammaDDot/2;
    ref.beta  = ref.theta - ref.gamma/2;  ref.betaDot  = ref.thetaDot - ref.gammaDot/2;  ref.betaDDot  = ref.thetaDDot - ref.gammaDDot/2;

    ev.t_end = ref.t(end);
    ref.ev = ev;
end

function [thP, gaP] = hold_both(thP, gaP, dwell, Ts)
    nd  = round(dwell/Ts);
    thP = [thP; thP(end)*ones(nd,1)];
    gaP = [gaP; gaP(end)*ones(nd,1)];
end

function [p, v, a] = derive_signals(p, Ts)
% Velocity/accel from the sampled position so all three stay consistent.
    v = gradient(p, Ts);  v(1) = 0;  v(end) = 0;
    a = gradient(v, Ts);
end

function p = point_to_point(dist, vmax, amax, tj, Ts)
% Rest-to-rest move of signed length 'dist' starting at 0 (first sample = 0).
% Trapezoidal velocity (triangular if too short), smoothed with a moving
% average of length tj to limit jerk. Total distance is preserved exactly.
    d = abs(dist);
    if d == 0, p = 0; return; end

    if d >= vmax^2/amax
        vpk = vmax;           t_acc = vmax/amax;  t_cru = (d - vmax^2/amax)/vmax;
    else
        vpk = sqrt(d*amax);   t_acc = vpk/amax;   t_cru = 0;
    end
    T  = 2*t_acc + t_cru;
    N  = ceil(T/Ts);
    tm = ((1:N)' - 0.5)*Ts;                          % sample mid-points

    v = zeros(N,1);
    acc = tm < t_acc;                      v(acc) = amax*tm(acc);
    cru = tm >= t_acc & tm < t_acc+t_cru;  v(cru) = vpk;
    dec = tm >= t_acc+t_cru;               v(dec) = max(vpk - amax*(tm(dec) - t_acc - t_cru), 0);

    Nj = round(tj/Ts);
    if Nj > 1
        v = conv(v, ones(Nj,1)/Nj);        % jerk limiting
    end
    pos = cumsum(v)*Ts;
    pos = pos*(d/pos(end));                % hit the target exactly
    p   = sign(dist)*[0; pos];
end

function out = simulate_axis(pos_ref, u_ff, pl, g, hw)
% One motor + controller, discrete time.
    Ts   = hw.Ts;
    N    = numel(pos_ref);
    sysd = c2d(ss([0 1; 0 -pl.B_eff/pl.J], [0; pl.Km/pl.J], eye(2), zeros(2,1)), Ts, 'zoh');
    Ad   = sysd.A;  Bd = sysd.B;
    q    = 2*pi/hw.enc_cpr;                 % encoder quantum [rad]
    ulim = hw.Uq_lim;

    x = [pos_ref(1); 0];
    y = zeros(N,1); w = zeros(N,1); u = zeros(N,1); sat = false(N,1);
    I = 0; D = 0; yPrev = round(x(1)/q)*q;

    for k = 1:N
        y(k) = x(1);  w(k) = x(2);
        ym = round(x(1)/q)*q;                       % quantized measurement
        e  = pos_ref(k) - ym;

        D     = (g.Tf*D - g.Kd*(ym - yPrev))/(g.Tf + Ts);   % filtered derivative on measurement
        yPrev = ym;

        Itry = I + g.Ki*Ts*e;
        uTry = g.Kp*e + Itry + D + u_ff(k);
        if abs(uTry) > ulim && sign(uTry) == sign(e)        % conditional integration
            uRaw = g.Kp*e + I + D + u_ff(k);
        else
            I = Itry;  uRaw = uTry;
        end
        u(k)   = min(max(uRaw, -ulim), ulim);
        sat(k) = abs(uRaw) > ulim;

        if k > hw.delay_samples, uApp = u(k - hw.delay_samples); else, uApp = 0; end
        x = Ad*x + Bd*uApp;
    end

    out.y = y;  out.w = w;  out.u = u;  out.sat = sat;
    out.iq = (u - pl.lamN*w)/hw.R_ph;       % q-axis current from resistive model [A]
    out.motor_rpm = abs(w)*hw.N_gear*60/(2*pi);
end

function m = compute_metrics(ref, A, B, hw, traj, an)
% Accuracy / lag / actuator metrics for one case.
    rmsf = @(x) sqrt(mean(x.^2));
    d    = 180/pi;
    t    = ref.t;  Ts = hw.Ts;  ev = ref.ev;

    thAct = (A.y + B.y)/2;   gaAct = A.y - B.y;
    eA = ref.alpha - A.y;    eB = ref.beta - B.y;
    eT = ref.theta - thAct;  eG = ref.gamma - gaAct;

    % time series for plotting
    m.yA = A.y; m.yB = B.y; m.eA = eA; m.eB = eB; m.eT = eT; m.eG = eG;
    m.uA = A.u; m.uB = B.u;

    % RMS and peak error [deg]
    m.rmsA_deg = d*rmsf(eA);  m.rmsB_deg = d*rmsf(eB);
    m.rmsT_deg = d*rmsf(eT);  m.rmsG_deg = d*rmsf(eG);
    m.maxT_deg = d*max(abs(eT));  m.maxG_deg = d*max(abs(eG));

    % Final (steady-state) error: mean over last 0.1 s
    nEnd = min(numel(t), max(1, round(0.1/Ts)));
    m.finalT_deg = d*mean(eT(end-nEnd+1:end));
    m.finalG_deg = d*mean(eG(end-nEnd+1:end));

    % Theta overshoot (percent of move) before the gamma move starts
    win = t <= ev.gamma_start;
    tgt = traj.theta_target;
    if abs(tgt) > 0
        m.overshootT_pct = 100*max(0, max(sign(tgt)*(thAct(win) - tgt)))/abs(tgt);
    else
        m.overshootT_pct = NaN;
    end

    % Settling times
    band = an.settle_band_deg/d;
    m.settleT_s = settling_time(t, eT, ev.theta_end,      ev.gamma_start, band);
    m.settleG_s = settling_time(t, eG, ev.gamma_end_last, ev.t_end,       band);

    % Best-fit pure delay between reference and response [ms]
    maxN = round(an.lag_search_ms/1000/Ts);
    m.lagT_ms = best_lag_ms(ref.theta, thAct, Ts, maxN);
    m.lagG_ms = best_lag_ms(ref.gamma, gaAct, Ts, maxN);

    % Cruise lag: steady following error / reference speed [ms]
    m.lagG_cruise_ms = cruise_lag_ms(ref.gammaDot, eG);

    % Actuator demand
    uAll = [A.u; B.u];
    m.peakUq_V      = max(abs(uAll));
    m.satFrac_pct   = 100*mean([A.sat; B.sat]);
    m.peakIq_A      = max(abs([A.iq; B.iq]));
    m.peakSpeed_rpm = max([A.motor_rpm; B.motor_rpm]);
end

function ts = settling_time(t, e, t0, t1, band)
% Time after t0 until |e| stays within band (window ends at t1). NaN if never.
    idx = find(t >= t0 & t <= t1);
    if numel(idx) < 2, ts = NaN; return; end
    bad = find(abs(e(idx)) > band, 1, 'last');
    if isempty(bad)
        ts = 0;
    elseif bad == numel(idx)
        ts = NaN;
    else
        ts = t(idx(bad+1)) - t0;
    end
end

function lag = best_lag_ms(r, y, Ts, maxN)
% Delay (whole samples) applied to r that best matches y in least squares.
    best = inf; lag = 0;
    for L = 0:maxN
        rd   = [r(1)*ones(L,1); r(1:end-L)];
        cost = mean((y - rd).^2);
        if cost < best, best = cost; lag = L; end
    end
    lag = 1000*lag*Ts;
end

function lag = cruise_lag_ms(v, e)
% Mean following error divided by speed over the near-peak-speed samples.
    vpk = max(abs(v));
    if vpk == 0, lag = NaN; return; end
    mask = abs(v) >= 0.95*vpk;
    lag  = 1000*mean(sign(v(mask)).*e(mask))/mean(abs(v(mask)));
end
