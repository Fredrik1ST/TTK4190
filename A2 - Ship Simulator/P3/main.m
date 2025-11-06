% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

% Add folder for 3-D visualization files
addpath(genpath('flypath3d_v2'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; clear WP_selector ILOS_guidance; close all; 
T_final = 6400;	        % Final simulation time (s)
h = 0.1;                % Sampling time (s)  (10 Hz)
U_ref   = 9;            % desired surge speed (m/s)

% initial states
eta_0 = [0 0 -110*pi/180]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
Qm_0 = 0;
x = [nu_0' eta_0' delta_0 n_0 Qm_0]'; % The state vector can be extended with additional states here

% === Estimator toggles ===
USE_ESKF = true;   % true -> use INS/ESKF (ins_euler), false -> use old model-based KF / noisy measurements
USE_KF   = ~USE_ESKF; % keep old toggle if used elsewhere

% Reference model initialization
xd = [0; 0; 0] ; % [psi_d, r_d, v_d]

% PID control initialization
e_int = 0;
%wb   = 0.06; zeta = 1.0; % Normal tuning from p2 to p3 4c
wb   = 0.03; zeta = 1.8; % Tuning for task 4d (good one)
%wb   = 0.1; zeta = 2.8; % Tuning for task 4d (ignore)

wn = wb/sqrt(1-2*zeta^2+sqrt(4*zeta^4-4*zeta^2+2));
K_nom = 7.4931e-03;
T_nom = 169.55;
m = T_nom/K_nom;
d = 1/K_nom;
k = 0;

kp = wn^2*m-k;
kd = 2*zeta*wn*m - d;
ki = (wn/10)*kp;

% Actuator limits
delta_max  = deg2rad(40);   % max rudder angle [rad]
Ddelta_max = deg2rad(5);    % max rudder rate [rad/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:h:T_final;                % Time vector
nTimeSteps = length(t);         % Number of time steps

simdata = zeros(nTimeSteps, 17); % Pre-allocate matrix for efficiency

%% -------------------- Shared noise & logs --------------------
rng(1)                             % reproducibility

% For the old model-based KF / noisy meas (kept for A/B testing)
sigma_psi = deg2rad(0.5);          % [rad] std of yaw measurement noise
sigma_r   = deg2rad(0.1);          % [rad/s] std of yaw-rate measurement noise

psi_meas_hist = zeros(nTimeSteps,1);
r_meas_hist   = zeros(nTimeSteps,1);
psi_true_hist = zeros(nTimeSteps,1);
r_true_hist   = zeros(nTimeSteps,1);
psi_hat_hist  = zeros(nTimeSteps,1);
r_hat_hist    = zeros(nTimeSteps,1);
b_hat_hist    = zeros(nTimeSteps,1);

%% -------------------- Old model-based KF matrices (fallback) --------------------
A = [0  1        0;
     0 -1/T_nom  -K_nom/T_nom;
     0  0        0];
B = [0; K_nom/T_nom; 0];
C = [1 0 0];
E = [0 0; 1 0; 0 1];
D = 0;

% First-order discretization
Ad = eye(3) + h*A;
Bd = h*B;
Cd = C;
Ed = h*E;

% KF initialization (fallback)
x_prd_KF = [0;0;0];  % [psi_hat; r_hat; b_hat]
P_prd_KF = diag([(deg2rad(30))^2, (deg2rad(0.01))^2, (deg2rad(1))^2]);
Qd_KF    = diag([1e-11, 1e-7]);    % var{w_r}, var{w_b}
Rd_KF    = sigma_psi^2;            % var of psi measurement

%% -------------------- ESKF/INS initialization (matches YOUR ins_euler.m) --------------------
% Rates for aiding
h_pos  = 0.2;          % 5 Hz GNSS
t_slow = 0;            % next slow measurement time

% IMU & aiding sensor noise (std)
sigma_accel       = 0.001;         % m/s^2
sigma_gyro        = 0.0000175;     % rad/s
sigma_gps_pos     = 0.025;         % m (RTK)
sigma_gps_vel     = 0.02;          % m/s (unused here, but keep for later)
sigma_compass_psi = deg2rad(0.5);  % rad

% Latitude (for gravity(mu) inside ins_euler). Pick your test latitude:
mu = deg2rad(63);                   % Trondheim-ish; adjust if needed

% INS full state vector (15x1): [p(3); v(3); b_acc(3); theta(3); b_ars(3)]
x_ins = [ x(4); x(5); 0;            % p_NED
          0; 0; 0;                  % v_NED
          0; 0; 0;                  % acc bias
          0; 0; x(6);               % [roll; pitch; yaw]
          0; 0; 0 ];                % gyro bias

P_prd = eye(15);

% Process & measurement noise (initial guesses; tune as needed)
Qd = diag([ 0.1 0.1 0,  0.001 0.001 0,  0 0 0,  0.1 0 0,  0 0 0.001 ]);
Rd = diag([0.1 0.1 0.1,   1 1 1,    0.1]);   % [pos xyz, gravity abc, compass psi]

% Gravity in NED for IMU synthesis
gN = [0;0;9.81];

% Velocity for finite-difference acceleration
vNED_prev = [0;0;0];

% ESKF logs
psi_ins_hist = zeros(nTimeSteps,1);
r_ins_hist   = zeros(nTimeSteps,1);
x_ins_hist   = zeros(nTimeSteps,1);
y_ins_hist   = zeros(nTimeSteps,1);

%% -------------------- Guidance setup --------------------
S = load('WP.mat');                 
WP = S.WP;

Delta_h = 600; % Lookahead distance
Rsw = 500;     % Switch radius
Rstop = 100;   % End radius
kappa = 1;     % ILOS design parameter

i_end = nTimeSteps;

%% =============================== MAIN SIM LOOP ===============================
for i = 1:nTimeSteps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disturbances: current (Part 2, 1a) and wind (Part 2, 1c)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Ocean current ---
    Vc = 1;                         % set 0 to disable
    betaVc = 45*pi/180;             % 45 deg CW from N -> NE (rad)
    uc = Vc * cos(betaVc - x(6));   % BODY Surge current (m/s)
    vc = Vc * sin(betaVc - x(6));   % BODY Sway current (m/s)
    nu_c = [ uc vc 0 ]';

    % --- Wind load ---
    Vw = 10;                % wind speed (m/s)
    beta_Vw = deg2rad(135); % wind direction
    rho_a = 1.247;          % air density
    c_y = 0.95;             % sway wind-force coefficient
    c_n = 0.15;             % yaw wind-moment coefficient
    L = 161;                % Length of ship
    A_Lw = 10*L;            % Area of side view above water.
    
    uw = Vw * cos(beta_Vw - x(6));
    vw = Vw * sin(beta_Vw - x(6));
    u_rw = x(1) - uw;
    v_rw = x(2) - vw;
    V_rw = sqrt(u_rw^2 + v_rw^2);
    gamma_rw = -atan2(v_rw,u_rw);
    Ywind = 0.5 * rho_a * V_rw^2 * c_y*sin(gamma_rw) * A_Lw;       % Fossen ch 10.1
    Nwind = 0.5 * rho_a * V_rw^2 * c_n*sin(2*gamma_rw) * A_Lw*L;   % Fossen ch 10.1
    tau_wind = [0 Ywind Nwind]';

    % --- Relative speed ocean currents ---
    u_rc = x(1) - uc;
    v_rc = x(2) - vc;
    U_rc = sqrt(u_rc^2 + v_rc^2);

    % --- Crab & sideslip ---
    beta_c   = atan2(x(2), max(1e-9, x(1)));   % crab angle
    psi      = x(6);
    chi      = ssa(psi + beta_c);              % course (approx if Vc=0)
    beta     = atan2(v_rc, max(1e-9, u_rc));   % sideslip

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guidance (Part 3 - Task 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xk1,yk1,xk,yk,last] = WP_selector(x(4),x(5), WP, Rsw, Rstop);
    [e_y,pi_p] = crossTrackError(xk1,yk1,xk,yk,x(4),x(5));
    chi_d  = LOS_guidance(e_y,pi_p, Delta_h);
    psi_ref = ILOS_guidance(e_y, pi_p, kappa, Delta_h, h);
    xd_dot = ref_model(xd, psi_ref);
    xd     = xd + h * xd_dot;
    psi_d  = xd(1);
    r_d    = xd(2);
    r_d_dot = xd_dot(2); %#ok<NASGU>
    u_d    = U_ref;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurements & Estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- True (for logging/comparison) -----
    psi_true = ssa(x(6));           % x = [u v r x y psi delta n Qm]'
    r_true   = x(3);
    psi_true_hist(i) = psi_true;
    r_true_hist(i)   = r_true;

    if USE_ESKF
        % ----- Compass (yaw) -----
        y_psi  = ssa(psi_true + sigma_compass_psi*randn);

        % Rotation NED<->BODY (yaw)
        c = cos(psi_true); s = sin(psi_true);
        Rnb = [ c -s 0;  s  c 0;  0  0 1];   % NED->BODY
        Rbn = Rnb.';

        % Body -> NED velocity & acceleration
        vNED = Rnb * [x(1); x(2); 0];
        aNED = (vNED - vNED_prev)/h;
        vNED_prev = vNED;

        % Specific force in BODY: f_b = R_bn*(aNED + gN)
        f_true = Rbn*(aNED + gN);
        w_true = [0; 0; r_true];              % p=q=0 for 3-DoF

        % Add IMU noise/bias (bias = 0 initially)
        f_imu = f_true + x_ins(7:9) + sigma_accel*randn(3,1);
        w_imu = w_true + x_ins(13:15) + sigma_gyro*randn(3,1);

        % GNSS aiding at 5 Hz
        do_pos = (t(i) > t_slow);
        if do_pos
            y_pos = [x(4); x(5); 0] + sigma_gps_pos*randn(3,1);
            % Your ins_euler signature:
            % [x_ins, P_prd] = ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos)
            [x_ins, P_prd] = ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi, y_pos);
            t_slow = t_slow + h_pos;
        else
            % [x_ins, P_prd] = ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi)
            [x_ins, P_prd] = ins_euler(x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, y_psi);
        end

        % Feedback signals from INS/ESKF
        psi_fb = ssa(x_ins(12));                 % theta(3) = yaw
        r_fb   = w_imu(3) - x_ins(15);           % gyro_z - b_gz

        % Log INS estimates
        psi_ins_hist(i) = psi_fb;
        r_ins_hist(i)   = r_fb;
        x_ins_hist(i)   = x_ins(1);
        y_ins_hist(i)   = x_ins(2);

    else
        % ----- Noisy measurements (old path) -----
        psi_meas = ssa(psi_true + randn*sigma_psi );
        r_meas   = r_true   + randn*sigma_r;
        psi_meas_hist(i) = psi_meas;
        r_meas_hist(i)   = r_meas;

        % ----- Old model-based KF -----
        [x_pst,P_pst,x_prd_KF,P_prd_KF] = KF(x_prd_KF,P_prd_KF,Ad,Bd,Ed,Cd,Qd_KF,Rd_KF,psi_meas,x(7));
        psi_hat_hist(i) = x_pst(1);
        r_hat_hist(i)   = x_pst(2);
        b_hat_hist(i)   = x_pst(3);

        psi_fb = x_pst(1);
        r_fb   = x_pst(2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Heading controller (PID with anti-windup)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e_psi = ssa(psi_fb - psi_d);
    e_r   = r_fb - r_d;
    e_u   = x(1) - u_d; %#ok<NASGU>
    delta_unsat = -(kp*e_psi + kd*e_r + ki*e_int);

    % Saturation
    delta_c = min(max(delta_unsat, -delta_max), delta_max);

    % Anti-windup (back-calculation)
    e_int_dot = e_psi - (1/ki)*(delta_c - delta_unsat);
    e_int = e_int + h * e_int_dot;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Speed control (open loop or closed loop)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_c = open_loop_speed_control(U_ref);
    % n_c = closed_loop_speed_control(u_d, x(1), h);  % if/when you want to switch

    % ship dynamics
    u_in = [delta_c n_c]';
    [xdot,tau_total] = ship(x,u_in,nu_c,tau_wind); %#ok<ASGLU>

    % store simulation data in a table (for testing)
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u_in(1) u_in(2) u_d psi_d r_d, chi, chi_d, beta_c, beta];     
 
    % Integrate dynamics
    x = rk4(@ship,h,x,u_in,nu_c,tau_wind);

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end

    % --- Stop on final waypoint ---
    if last
        i_end = i;  % last valid index
        % Log one last row for tidy plots
        psi_d = x(6); r_d = 0; u_d = U_ref;
        simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u_in(1) u_in(2) u_d psi_d r_d, chi, chi_d, beta_c, beta];  
        fprintf('  Reached final waypoint at t = %.1f s (step %d)\n', (i-1)*h, i);
        break
    end
end

% Trim outputs to actual length
simdata = simdata(1:i,:);
t = t(1:i);

% Trim estimator logs
psi_meas_hist = psi_meas_hist(1:length(t));
r_meas_hist   = r_meas_hist(1:length(t));
psi_true_hist = psi_true_hist(1:length(t));
r_true_hist   = r_true_hist(1:length(t));
psi_hat_hist  = psi_hat_hist(1:length(t));
r_hat_hist    = r_hat_hist(1:length(t));
b_hat_hist    = b_hat_hist(1:length(t));
psi_ins_hist  = psi_ins_hist(1:length(t));
r_ins_hist    = r_ins_hist(1:length(t));
x_ins_hist    = x_ins_hist(1:length(t));
y_ins_hist    = y_ins_hist(1:length(t));

% Save ESKF results (Part 3b)
simdata_ESKF.t   = t;
simdata_ESKF.psi = psi_ins_hist;
simdata_ESKF.r   = r_ins_hist;
simdata_ESKF.x   = x_ins_hist;
simdata_ESKF.y   = y_ins_hist;
save('simdata_ESKF.mat','simdata_ESKF');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
xN          = simdata(:,4);                 % m (North)
yE          = simdata(:,5);                 % m (East)
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c_out     = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         = simdata(:,13);                % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
chi         = simdata(:,14);                % rad
chi_deg     = (180/pi)* chi;                % deg
chi_d       = simdata(:,15);                % rad
chi_d_deg   = (180/pi)*chi_d;               % deg
beta_c      = simdata(:,16);                % rad
beta_c_deg  = (180/pi)*beta_c;              % deg
beta        = simdata(:,17);                % rad
beta_deg    = (180/pi)*beta;                % deg

figure(3); clf; figure(gcf)
subplot(311)
plot(yE,xN,'linewidth',2); hold on; axis('equal')
if USE_ESKF
    plot(y_ins_hist, x_ins_hist, '--', 'LineWidth', 1.5);
    legend('true track','ESKF track','Location','best')
end
title('North-East positions'); xlabel('East (m)'); ylabel('North (m)'); 

subplot(312)
plot(t,psi_deg,'linewidth',2); hold on
plot(t,psi_d_deg,'linewidth',2);
if USE_ESKF
    plot(t, rad2deg(psi_ins_hist), '--', 'LineWidth', 1.5);
    legend('actual yaw','desired yaw','ESKF \psi','Location','best')
else
    legend('actual yaw','desired yaw','Location','best')
end
title('Yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 

subplot(313)
plot(t,r_deg,'linewidth',2); hold on
plot(t,r_d_deg,'linewidth',2);
if USE_ESKF
    plot(t, rad2deg(r_ins_hist), '--', 'LineWidth', 1.5);
    legend('actual yaw rate','desired yaw rate','ESKF r','Location','best')
else
    legend('actual yaw rate','desired yaw rate','Location','best')
end
title('Yaw rate'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 

figure(2); clf; figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c_out,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')

figure(4); clf; figure(gcf)
plot(t, chi_deg,   'LineWidth', 2); hold on;
plot(t, chi_d_deg, 'LineWidth', 2);
plot(t, psi_deg,   'LineWidth', 2);
plot(t, psi_d_deg, 'LineWidth', 2);
plot(t, beta_c_deg,'--',        'LineWidth', 1.5);
plot(t, beta_deg,  '--',        'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Angle (deg)');
title('Course (χ), Desired Course (χ_d), Heading (ψ), Desired Heading (ψ), Crab (β_c), Sideslip (β)');
legend('\chi','\chi_d','\psi', '\psi_d','\beta_c','\beta','Location','best');

% Old-KF specific comparison plots (only if using old KF)
if ~USE_ESKF
    figure(5); clf; figure(gcf)
    subplot(2,1,1)
    plot(t, rad2deg(psi_meas_hist), 'LineWidth', 1.5); hold on
    plot(t, rad2deg(psi_true_hist), '--', 'LineWidth', 2)
    grid on; xlabel('Time (s)'); ylabel('\psi (deg)')
    title('Yaw angle: true vs noisy measurement')
    legend('noisy \psi (0.5^\circ \sigma)', 'true \psi','Location','best')
    subplot(2,1,2)
    plot(t, rad2deg(r_meas_hist), 'LineWidth', 1.5); hold on
    plot(t, rad2deg(r_true_hist),'--', 'LineWidth', 2)
    grid on; xlabel('Time (s)'); ylabel('r (deg/s)')
    title('Yaw rate: true vs noisy measurement')
    legend('noisy r (0.1^\circ/s \sigma)','true r','Location','best')

    figure(6); clf; figure(gcf)
    subplot(3,1,1)
    plot(t, rad2deg(psi_true_hist), 'LineWidth', 2); hold on
    plot(t, rad2deg(psi_hat_hist),  '--', 'LineWidth', 1.5)
    grid on; ylabel('\psi (deg)'); title('Yaw angle: true vs KF estimate')
    legend('true \psi','estimated $\hat{\psi}$', 'Interpreter','latex','Location','best')
    subplot(3,1,2)
    plot(t, rad2deg(r_true_hist), 'LineWidth', 2); hold on
    plot(t, rad2deg(r_hat_hist),  '--', 'LineWidth', 1.5)
    grid on; ylabel('r (deg/s)'); title('Yaw rate: true vs KF estimate')
    legend('true r','estimated $\hat{r}$','Interpreter','latex', 'Location','best')
    subplot(3,1,3)
    plot(t, rad2deg(b_hat_hist), 'LineWidth', 2)
    grid on; xlabel('Time (s)'); ylabel('b (deg)')
    title('Rudder bias estimate (no true bias available)')
end

%% 3-D visualization objects 
z = zeros(length(xN),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);

new_object('flypath3d_v2/ship1.mat',[xN,yE,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(xN)),max(abs(yE)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);

% Plot the path
pathplotter(xN, yE);
