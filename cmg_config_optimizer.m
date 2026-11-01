%% CMG_CONFIG_OPTIMIZER
% Compares CMG array configurations (single, scissor pair, two scissor pairs,
% pyramid, double-gimbal) under a common MASS and ENVELOPE-VOLUME budget.
%
% For every configuration the script:
%   1. Finds the largest swept sphere each flywheel may occupy inside the envelope
%   2. Optimizes a spoked (or solid) flywheel for maximum angular momentum H = I*omega
%      subject to mass, swept-volume, rim-stress and rpm limits
%   3. Evaluates how much of H the array can turn into torque about each axis
%      (LP over gimbal rates, sampled over gimbal states)
%   4. Reports torque [N*m], torque/kg, and optional mass sweep + torque clouds
%
% Key physics: torque = H * (J_norm * gimbal_rate), where J_norm depends only on the
% configuration geometry. So the flywheel optimization (max H) and the geometry
% evaluation (kappa) decouple, and configurations differ through: number of
% flywheels (mass/volume split), packing, overhead mass, and kappa.
%
% Requirements: base MATLAB. Optimization Toolbox (fmincon, linprog) is used when
% present; otherwise slower/rougher fallbacks (fminsearch, pinv) are used.
%
% All numbers marked PLACEHOLDER should be replaced with your project values.

close all; clear; clc;
rng(1);                                   % repeatable Monte Carlo / multistart

%% ============================ USER INPUTS ============================
% --- Envelope and budget (whole CMG array, all units together) ---
p.env.type = 'box';             % 'box' | 'cylinder' | 'sphere'
p.env.dims = [0.15 0.10 0.10];  % box [Lx Ly Lz] m | cylinder [D H] m | sphere [D] m   PLACEHOLDER
p.M_total  = 1.5;               % total mass budget, flywheels + motors + gimbals [kg]  PLACEHOLDER

% --- Actuation and evaluation ---
p.rate_max     = 30;            % max gimbal rate per gimbal [rad/s]
p.gim_lim      = deg2rad(180);  % gimbal angle sampling range, +/- [rad] (reduce for limited-travel gimbals)
p.n_states     = 1500;          % gimbal states sampled per configuration/axis
p.scissor_strict = true;        % true: scissor pairs only counter-rotate (gdot2 = -gdot1), i.e. true single-axis scissor operation
p.off_axis_tol = 0.25;          % allowed off-axis torque as a fraction of on-axis torque
p.rank_metric  = 'mean';        % 'mean' | 'peak' | 'p10' (torque available in >=90% of states)
p.focus_axis   = 3;             % axis for sweeps / torque-per-kg: 1=x, 2=y, 3=z
p.req_axes     = [2 3];         % axes that must ALL be covered (e.g. pitch and yaw): min over these is reported

% --- Flywheel geometry rules ---
p.fw.shape     = 'spoked';      % 'spoked' (rim + spokes + hub) | 'solid' (plain disc, for baseline comparison)
p.fw.N_spokes  = 4;
p.fw.r_hub     = 0.006;         % hub radius (motor shaft/rotor interface) [m]          PLACEHOLDER
p.fw.min_rim   = 0.002;         % minimum rim radial thickness [m]
p.fw.min_spoke = 0.0015;        % minimum spoke tangential width [m]
p.fw.min_Ls    = 0.002;         % minimum spoke length (rim inner radius - hub radius) [m]
p.fw.min_w     = 0.004;         % minimum flywheel axial width [m]
p.fw.max_w     = 0.050;         % maximum flywheel axial width [m]
p.fw.rpm_max   = 30000;         % motor/bearing speed limit [rpm]                       PLACEHOLDER
p.fw.clear     = 0.004;         % housing + running clearance added to flywheel radius and half-width [m]

% --- Material (steel-like; density matches your original scripts) ---
p.mat.rho      = 7500;          % [kg/m^3]
p.mat.sigma_y  = 600e6;         % yield strength [Pa]                                   PLACEHOLDER
p.mat.SF       = 3;             % safety factor on yield
p.mat.nu       = 0.3;           % Poisson ratio
p.mat.Kspoke   = 1.5;           % stress concentration / spoke-bending knockdown on rim hoop stress (first-order guess)

% --- Non-flywheel mass ---
p.ovh.frac     = 0.8;           % motor + bearings + housing mass as a fraction of flywheel mass  PLACEHOLDER
p.ovh.fixed    = 0.040;         % fixed mass per GIMBAL (motor, encoder, bracket) [kg]            PLACEHOLDER

% --- Configuration options ---
p.n_pyr    = 4;                 % flywheels in pyramid
p.beta_deg = 54.74;             % pyramid gimbal-axis angle from z [deg] (54.74 = isotropic 4-CMG)

% --- Solver / output switches ---
p.n_starts    = 12;             % multistart count for flywheel optimization
p.do_sweep    = true;           % torque vs. total-mass sweep
p.sweep_mult  = linspace(0.4, 2.5, 7);   % multiples of M_total for the sweep
p.do_envelope = true;           % 3D torque point clouds (like your original script)
p.env_samples = 20000;
p.use_lp      = exist('linprog','file') > 0;
p.use_fmincon = exist('fmincon','file') > 0;
if ~p.use_lp,      warning('linprog not found: using pinv fallback (ignores off_axis_tol, lower bound on torque).'); end
if ~p.use_fmincon, warning('fmincon not found: using penalized fminsearch fallback (rougher optimum).'); end

%% ============================ MAIN ============================
cfgs = make_configs(p);
nC   = numel(cfgs);
axn  = {'x','y','z'};

res = struct('cfg',cell(1,nC),'S',[],'des',[],'kappa',[],'T',[]);
for i = 1:nC
    fprintf('Evaluating %-28s ...\n', cfgs(i).name);
    res(i).cfg   = cfgs(i);
    res(i).S     = torque_stats(cfgs(i), p);
    res(i).kappa = get_metric(res(i).S, p.rank_metric);          % 1x3 normalized torque gain
    res(i).des   = optimize_flywheel(cfgs(i), p, p.M_total);
    if res(i).des.feasible
        res(i).T = res(i).des.H * p.rate_max * res(i).kappa;      % [N*m] about x,y,z
    else
        res(i).T = zeros(1,3);
    end
end

%% ============================ REPORT ============================
fprintf('\n=== Envelope: %s %s m | Mass budget %.2f kg | Flywheel: %s | Rank metric: %s ===\n', ...
    p.env.type, mat2str(p.env.dims), p.M_total, p.fw.shape, p.rank_metric);
fprintf('%-28s %2s %6s %6s %8s %7s %6s %7s %8s | %7s %7s %7s | %7s %9s\n', ...
    'Config','n','R[mm]','w[mm]','rpm','H[Nms]','m_fw','M_tot','Vol use', ...
    'Tx','Ty','Tz','min(req)','Nm/kg(f)');
for i = 1:nC
    d = res(i).des;
    if ~d.feasible
        fprintf('%-28s  INFEASIBLE (mass or envelope too small for this layout)\n', res(i).cfg.name);
        continue;
    end
    Treq = min(res(i).T(p.req_axes));
    fprintf('%-28s %2d %6.1f %6.1f %8.0f %7.3f %6.3f %7.3f %7.0f%% | %7.2f %7.2f %7.2f | %7.2f %9.2f\n', ...
        res(i).cfg.name, d.n, d.g.R*1e3, d.g.w*1e3, d.g.rpm, d.H, d.m, d.M_used, 100*d.vol_frac, ...
        res(i).T(1), res(i).T(2), res(i).T(3), Treq, res(i).T(p.focus_axis)/d.M_used);
end

fprintf('\n--- Peak / mean / p10 of normalized gain kappa (torque = kappa*H*rate) ---\n');
for i = 1:nC
    S = res(i).S;
    fprintf('%-28s peak [%.2f %.2f %.2f]  mean [%.2f %.2f %.2f]  p10 [%.2f %.2f %.2f]\n', ...
        res(i).cfg.name, S.peak, S.mean, S.p10);
end

fprintf('\n--- Optimized flywheel details (spoked: rim t, spoke b x ws) & binding limits ---\n');
for i = 1:nC
    d = res(i).des;  if ~d.feasible, continue; end
    if strcmp(p.fw.shape,'spoked')
        fprintf('%-28s rim t=%.1f mm, spoke %.1f x %.1f mm, tip %.0f m/s, stress %.0f%% of allowable | binding: %s\n', ...
            res(i).cfg.name, d.g.t*1e3, d.g.b*1e3, d.g.ws*1e3, d.vtip, 100*d.sigma/(p.mat.sigma_y/p.mat.SF), d.binding);
    else
        fprintf('%-28s tip %.0f m/s, stress %.0f%% of allowable | binding: %s\n', ...
            res(i).cfg.name, d.vtip, 100*d.sigma/(p.mat.sigma_y/p.mat.SF), d.binding);
    end
end

%% ============================ PLOTS ============================
% 1) Torque comparison
figure('Color','w','Name','Torque comparison');
Tmat = vertcat(res.T);
subplot(1,2,1);
bar(Tmat);
set(gca,'XTickLabel',{cfgs.name},'XTickLabelRotation',20);
ylabel(sprintf('Torque (%s) [N\\cdot m]', p.rank_metric)); legend(axn,'Location','best'); grid on;
title('Torque by axis, same mass + envelope');
subplot(1,2,2);
Mused = ones(nC,1);
for i = 1:nC, if res(i).des.feasible, Mused(i) = res(i).des.M_used; end, end
bar(Tmat(:,p.focus_axis) ./ Mused);
set(gca,'XTickLabel',{cfgs.name},'XTickLabelRotation',20);
ylabel(sprintf('Torque about %s per kg [N\\cdot m/kg]', axn{p.focus_axis})); grid on;
title('Specific torque (focus axis)');

% 2) Flywheel drawings
figure('Color','w','Name','Optimized flywheels');
for i = 1:nC
    subplot(1,nC,i); hold on; axis equal off;
    if res(i).des.feasible, draw_flywheel(res(i).des, p); end
    title(sprintf('%s\n%s', res(i).cfg.name, flywheel_label(res(i).des)), 'FontSize',8);
end

% 3) Mass sweep
if p.do_sweep
    p2 = p;  p2.n_starts = max(4, round(p.n_starts/2));
    Tsw = nan(nC, numel(p.sweep_mult));
    for i = 1:nC
        for j = 1:numel(p.sweep_mult)
            d = optimize_flywheel(cfgs(i), p2, p.sweep_mult(j)*p.M_total);
            if d.feasible, Tsw(i,j) = d.H * p.rate_max * res(i).kappa(p.focus_axis); end
        end
    end
    figure('Color','w','Name','Mass sweep'); hold on; grid on;
    for i = 1:nC, plot(p.sweep_mult*p.M_total, Tsw(i,:), '-o','LineWidth',1.5); end
    xlabel('Total mass budget [kg]');
    ylabel(sprintf('Torque about %s (%s) [N\\cdot m]', axn{p.focus_axis}, p.rank_metric));
    legend({cfgs.name},'Location','best');
    title('Torque vs. mass (envelope fixed): flattening = envelope/stress limited');
end

% 4) 3D torque clouds (same idea as the original cmg_config_calculator.m)
if p.do_envelope
    figure('Color','w','Name','Torque envelopes');
    for i = 1:nC
        subplot(2, ceil(nC/2), i);
        if ~res(i).des.feasible, title([res(i).cfg.name ' (infeasible)']); continue; end
        cloud = zeros(3, p.env_samples);
        for k = 1:p.env_samples
            gam = cfgs(i).sampler(1);                 % respects scissor symmetry / gimbal limits
            gd  = p.rate_max*(2*rand(size(cfgs(i).rate_map,2),1) - 1);
            cloud(:,k) = res(i).des.H * build_jacobian(cfgs(i), gam) * cfgs(i).rate_map * gd;
        end
        scatter3(cloud(1,:), cloud(2,:), cloud(3,:), 3, sqrt(sum(cloud.^2,1)), 'filled');
        axis equal; grid on; colormap(jet); view(45,30);
        xlabel('X'); ylabel('Y'); zlabel('Z'); title(cfgs(i).name);
    end
end

%% ============================ LOCAL FUNCTIONS ============================

% ---------- Configurations ----------
function cfgs = make_configs(p)
    x = [1;0;0]; y = [0;1;0]; z = [0;0;1];
    lim = p.gim_lim;
    uni = @(N,m) lim*(2*rand(N,m) - 1);

    c1 = new_cfg('Single CMG', {mk_unit(z, x, 1)}, 0);
    c1.sampler = @(N) uni(N,1);

    c2 = new_cfg('Scissor pair (1-axis)', {mk_unit(z, x, 1), mk_unit(z, x, 2)}, 0);
    c2.sampler = @(N) sample_scissor(N, lim, 1);          % gamma2 = -gamma1
    if p.scissor_strict, c2.rate_map = [1; -1]; end

    % Second pair: momentum along y so the two pairs output torque about z and y
    c3 = new_cfg('2 scissor pairs (2-axis)', ...
        {mk_unit(z,x,1), mk_unit(z,x,2), mk_unit(y,x,3), mk_unit(y,x,4)}, 0);
    c3.sampler = @(N) sample_scissor(N, lim, 2);
    if p.scissor_strict, c3.rate_map = blkdiag([1;-1], [1;-1]); end

    beta = deg2rad(p.beta_deg);  n = p.n_pyr;  u = cell(1,n);
    for i = 1:n
        psi = (i-1)*2*pi/n;
        g = [sin(beta)*cos(psi); sin(beta)*sin(psi); cos(beta)];
        h = cross([cos(psi); sin(psi); 0], [0;0;1]);       % same construction as your original script
        u{i} = mk_unit(h, g, i);
    end
    c4 = new_cfg(sprintf('%d-CMG pyramid', n), u, 0);
    c4.sampler = @(N) uni(N,n);

    % Double gimbal: outer axis z, inner axis y (carried by outer), spin axis x at zero angle.
    % Larger swept volume (outer ring), so extra clearance is added.
    c5 = new_cfg('Double-gimbal CMG', {mk_unit(x, [z y], [1 2])}, 0.004);
    c5.sampler = @(N) uni(N,2);

    cfgs = [c1 c2 c3 c4 c5];
end

function u = mk_unit(h0, g0, idx)
    u.h0 = h0(:)/norm(h0);
    u.g0 = g0./vecnorm(g0);      % columns = gimbal axes, outermost first (reference frame, zero angles)
    u.idx = idx;                 % indices into the global gimbal-angle vector
end

function c = new_cfg(name, units, extra_clear)
    c.name = name;
    c.units = [units{:}];
    c.n_fw = numel(units);
    c.n_gim = sum(arrayfun(@(u) numel(u.idx), c.units));
    c.gim_per_unit = c.n_gim / c.n_fw;
    c.extra_clear = extra_clear;
    c.sampler = [];
    c.rate_map = eye(c.n_gim);   % gimbal rates = rate_map * independent rate DOFs (each limited to +/- rate_max)
end

function G = sample_scissor(N, lim, npairs)
    G = zeros(N, 2*npairs);
    for k = 1:npairs
        g = lim*(2*rand(N,1) - 1);
        G(:, 2*k-1:2*k) = [g, -g];
    end
end

% ---------- Kinematics ----------
function R = rodrigues(a, th)
    a = a/norm(a);
    K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    R = eye(3) + sin(th)*K + (1-cos(th))*K*K;
end

function J = build_jacobian(cfg, gam)
    % Columns give torque-rate per unit gimbal rate for H = 1:  dh/dt = sum_j gdot_j * (a_j x h)
    % (the torque on the vehicle is the negative of this; magnitudes are identical)
    J = zeros(3, cfg.n_gim);
    for k = 1:numel(cfg.units)
        u = cfg.units(k);
        m = size(u.g0, 2);
        R = eye(3);  a = zeros(3, m);
        for j = 1:m
            a(:,j) = R*u.g0(:,j);                  % gimbal axis in world frame at current angles
            R = rodrigues(a(:,j), gam(u.idx(j)))*R;
        end
        h = R*u.h0;
        for j = 1:m
            J(:, u.idx(j)) = cross(a(:,j), h);
        end
    end
end

% ---------- Torque gain statistics ----------
function S = torque_stats(cfg, p)
    N = p.n_states;
    G = cfg.sampler(N);
    s = zeros(N,3);
    E = eye(3);
    Pp = cell(1,3);
    for a = 1:3, Pp{a} = null(E(:,a)'); end        % 3x2 basis perpendicular to each axis
    lp.opts = [];
    if p.use_lp, lp.opts = optimoptions('linprog','Display','off'); end
    for k = 1:N
        J = build_jacobian(cfg, G(k,:)) * cfg.rate_map;   % restrict to allowed rate combinations
        for a = 1:3
            s(k,a) = max_axis_gain(J, E(:,a), Pp{a}, p.off_axis_tol, p.use_lp, lp.opts);
        end
    end
    ss = sort(s, 1);
    S.s = s;
    S.peak = max(s, [], 1);
    S.mean = mean(s, 1);
    S.p10  = ss(max(1, ceil(0.10*N)), :);
end

function k = get_metric(S, metric)
    switch metric
        case 'peak', k = S.peak;
        case 'mean', k = S.mean;
        case 'p10',  k = S.p10;
        otherwise, error('rank_metric must be mean, peak or p10');
    end
end

function s = max_axis_gain(J, a, Pperp, eps_off, use_lp, opts)
    % Max of a'*J*gd  s.t. |gd_i| <= 1 (normalized) and |Pperp'*J*gd| <= eps*(a'*J*gd)
    n  = size(J,2);
    ta = a'*J;
    if use_lp
        A = [ Pperp'*J - eps_off*ta ; -Pperp'*J - eps_off*ta ];
        [gd,~,flag] = linprog(-ta', A, zeros(4,1), [], [], -ones(n,1), ones(n,1), opts);
        if flag <= 0, s = 0; else, s = max(ta*gd, 0); end
    else
        gd = pinv(J)*a;                               % fallback: exact pure-axis, min-norm direction
        if norm(J*gd - a) > 1e-6, s = 0; else, s = 1/max(abs(gd)); end
    end
end

% ---------- Packing: largest sphere diameter so that n spheres fit in the envelope ----------
function D = max_sphere_diameter(env, n)
    switch env.type
        case 'box'
            D = box_pack(env.dims(:)', n);
        case 'cylinder'
            D = cyl_pack(env.dims(1), env.dims(2), n);
        case 'sphere'
            Dc = env.dims(1);
            if n == 1,     D = Dc;
            elseif n == 2, D = Dc/2;
            else,          s = Dc/sqrt(3);  D = box_pack([s s s], n);   % conservative: inscribed cube
            end
        otherwise
            error('env.type must be box, cylinder or sphere');
    end
end

function D = box_pack(L, n)      % cubic-grid packing (conservative vs. true sphere packing)
    cand = [];
    for i = 1:3, cand = [cand, L(i)./(1:n)]; end %#ok<AGROW>
    cand = sort(unique(cand), 'descend');
    D = 0;
    for c = cand
        if prod(floor(L./c + 1e-9)) >= n, D = c; return; end
    end
end

function D = cyl_pack(Dc, H, n)  % layers x circles-in-circle (exact table for up to 10 per layer)
    ratio = [1 2 2.1547 2.4142 2.7013 3 3 3.3048 3.6131 3.8130];  % container diam / circle diam
    cand = [H./(1:n), Dc./ratio];
    cand = sort(unique(cand), 'descend');
    D = 0;
    for c = cand
        kmax = find(ratio <= Dc/c + 1e-9, 1, 'last');
        if isempty(kmax), kmax = 0; end
        if floor(H/c + 1e-9)*kmax >= n, D = c; return; end
    end
end

function V = env_volume(env)
    switch env.type
        case 'box',      V = prod(env.dims);
        case 'cylinder', V = pi/4*env.dims(1)^2*env.dims(2);
        case 'sphere',   V = pi/6*env.dims(1)^3;
    end
end

% ---------- Flywheel model ----------
function g = unpack(x, p)
    g.R = x(1)*1e-3;  g.w = x(2)*1e-3;
    if strcmp(p.fw.shape,'spoked')
        g.t = x(3)*1e-3;  g.b = x(4)*1e-3;  g.ws = x(5)*g.w;  g.rpm = x(6)*1e3;
    else
        g.rpm = x(3)*1e3;
    end
end

function e = fw_eval(x, p, cl)
    g = unpack(x, p);  m = p.mat;  fw = p.fw;
    om = g.rpm*2*pi/60;
    if strcmp(fw.shape,'spoked')
        Ro = g.R;  Ri = g.R - g.t;  rh = fw.r_hub;  N = fw.N_spokes;
        m_rim = m.rho*pi*(Ro^2 - Ri^2)*g.w;      I_rim = 0.5*m_rim*(Ro^2 + Ri^2);
        Ls    = Ri - rh;
        m_sp  = N*m.rho*g.b*g.ws*Ls;             I_sp  = N*m.rho*g.b*g.ws*(Ri^3 - rh^3)/3;
        m_hub = m.rho*pi*rh^2*g.w;               I_hub = 0.5*m_hub*rh^2;
        e.m = m_rim + m_sp + m_hub;
        e.I = I_rim + I_sp + I_hub;
        % Peak hoop stress at inner radius of a free rotating thick ring, times knockdown
        e.sigma = m.Kspoke * m.rho*om^2/4 * ((3+m.nu)*Ro^2 + (1-m.nu)*Ri^2);
        e.Ls = Ls;
    else
        e.m = m.rho*pi*g.R^2*g.w;
        e.I = 0.5*e.m*g.R^2;
        e.sigma = (3+m.nu)/8 * m.rho*om^2*g.R^2;   % centre of a solid rotating disc
        e.Ls = Inf;
    end
    e.omega = om;  e.H = e.I*om;  e.vtip = om*g.R;
    e.Ds = 2*sqrt((g.R + cl)^2 + (g.w/2 + cl)^2);  % diameter of sphere swept when the unit gimbals through 360 deg
    e.g = g;
end

function [c, ceq] = fw_constraints(e, m_max, Dmax, p)
    c = [e.m/m_max - 1;
         e.Ds/Dmax - 1;
         e.sigma/(p.mat.sigma_y/p.mat.SF) - 1];
    if isfinite(e.Ls), c(end+1,1) = (p.fw.min_Ls - e.Ls)/p.fw.min_Ls; end
    ceq = [];
end

function f = penalized(x, lb, ub, obj, nonl)
    xc = min(max(x, lb), ub);
    c  = nonl(xc);
    f  = obj(xc) + 1e4*sum(max(c,0).^2) + 1e4*sum(((x-xc)./max(ub-lb,eps)).^2);
end

function des = optimize_flywheel(cfg, p, M_total)
    des.feasible = false;
    n = cfg.n_fw;
    m_fixed = p.ovh.fixed*cfg.gim_per_unit;
    m_max = (M_total/n - m_fixed)/(1 + p.ovh.frac);       % flywheel mass allowed per unit
    cl = p.fw.clear + cfg.extra_clear;
    Dmax = max_sphere_diameter(p.env, n);
    fw = p.fw;  rh = fw.r_hub*1e3;
    if m_max <= 0 || Dmax <= 2*cl, return; end

    R_ub = (Dmax/2 - cl)*1e3;
    w_ub = min(fw.max_w, Dmax - 2*cl)*1e3;   w_lb = fw.min_w*1e3;
    if strcmp(fw.shape,'spoked')
        R_lb = rh + (fw.min_Ls + fw.min_rim)*1e3;
        lb = [R_lb, w_lb, fw.min_rim*1e3, fw.min_spoke*1e3, 0.2, 2];
        ub = [R_ub, w_ub, R_ub - rh - fw.min_Ls*1e3, 0.8*2*pi*rh/fw.N_spokes, 1.0, fw.rpm_max/1e3];
    else
        lb = [rh + 1, w_lb, 2];
        ub = [R_ub, w_ub, fw.rpm_max/1e3];
    end
    if any(ub < lb), return; end

    obj  = @(x) -fw_eval(x, p, cl).H;
    nonl = @(x) fw_constraints(fw_eval(x, p, cl), m_max, Dmax, p);
    if p.use_fmincon
        opt = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                           'MaxFunctionEvaluations',4000,'MaxIterations',400);
    end

    bestH = -inf;  xbest = [];
    for s = 1:p.n_starts
        x0 = lb + rand(size(lb)).*(ub - lb);
        try
            if p.use_fmincon
                xs = fmincon(obj, x0, [], [], [], [], lb, ub, nonl, opt);
            else
                xs = fminsearch(@(x) penalized(x, lb, ub, obj, @(q) nonl_c(nonl,q)), x0, ...
                                optimset('MaxFunEvals',6000,'MaxIter',3000,'Display','off'));
                xs = min(max(xs, lb), ub);
            end
        catch
            continue;
        end
        e = fw_eval(xs, p, cl);
        if all(nonl_c(nonl, xs) <= 1e-4) && e.H > bestH, bestH = e.H; xbest = xs; end
    end
    if isempty(xbest), return; end

    e = fw_eval(xbest, p, cl);
    c = nonl_c(nonl, xbest);
    names = {'mass','envelope','stress','min spoke length'};
    bind = names(abs(c(:)') < 2e-3);
    if xbest(end) >= ub(end)*0.999, bind{end+1} = 'rpm limit'; end
    if isempty(bind), bind = {'none'}; end

    des = e;
    des.feasible = true;
    des.n = n;
    des.M_used = n*(e.m*(1 + p.ovh.frac) + m_fixed);
    des.vol_frac = n*pi/6*e.Ds^3 / env_volume(p.env);
    des.binding = strjoin(bind, ', ');
    des.Dmax = Dmax;
    des.m_max = m_max;
end

function c = nonl_c(nonl, x)
    [c, ~] = nonl(x);
end

% ---------- Drawing ----------
function draw_flywheel(des, p)
    g = des.g;  R = g.R*1e3;  th = linspace(0, 2*pi, 120);
    if strcmp(p.fw.shape,'spoked')
        Ri = (g.R - g.t)*1e3;  rh = p.fw.r_hub*1e3;  N = p.fw.N_spokes;
        patch([R*cos(th), Ri*cos(fliplr(th))], [R*sin(th), Ri*sin(fliplr(th))], [0.4 0.5 0.7], 'EdgeColor','none');
        for k = 1:N
            ph = 2*pi*(k-1)/N;  Rm = [cos(ph) -sin(ph); sin(ph) cos(ph)];
            q = Rm*[rh Ri Ri rh; -g.b*1e3/2 -g.b*1e3/2 g.b*1e3/2 g.b*1e3/2];
            patch(q(1,:), q(2,:), [0.4 0.5 0.7], 'EdgeColor','none');
        end
        patch(rh*cos(th), rh*sin(th), [0.25 0.3 0.45], 'EdgeColor','none');
    else
        patch(R*cos(th), R*sin(th), [0.4 0.5 0.7], 'EdgeColor','none');
    end
    plot(R*cos(th), R*sin(th), 'k-');
end

function s = flywheel_label(d)
    if ~d.feasible, s = 'infeasible'; return; end
    s = sprintf('D=%.0f mm, w=%.0f mm, %.0f rpm', 2*d.g.R*1e3, d.g.w*1e3, d.g.rpm);
end
