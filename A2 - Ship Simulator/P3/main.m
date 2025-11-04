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
h = 0.1;                % Sampling time (s)
U_ref   = 9;            % desired surge speed (m/s)

% initial states
eta_0 = [0 0 -110*pi/180]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
Qm_0 = 0;
x = [nu_0' eta_0' delta_0 n_0 Qm_0]'; % The state vector can be extended with addional states here

USE_KF = true;  % toggle: true -> use KF estimates; false -> use noisy measurements


% Reference model initialization
xd = [0; 0; 0];  % [psi_d, r_d, v_d]

% PID control initialization
e_int = 0;
%wb   = 0.06; zeta = 1.0;
wb   = 0.03; zeta = 1.8; % Tuning for task 4d

wn = wb/sqrt(1-2*zeta^2+sqrt(4*zeta^4-4*zeta^2+2));
K_nom = 7.4931e-03;
T_nom = 169.55;
% Nomoto when Uref = 9 m/s
%K_nom = 7.68e-03;
%T_nom = 174.2;
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

% Kalman filter initialization
rng(1)                             % reproducibility
sigma_psi = deg2rad(0.5);          % [rad] std of yaw measurement noise
sigma_r   = deg2rad(0.1);          % [rad/s] std of yaw-rate measurement noise

psi_meas_hist = zeros(nTimeSteps,1);
r_meas_hist   = zeros(nTimeSteps,1);
psi_true_hist = zeros(nTimeSteps,1);
r_true_hist   = zeros(nTimeSteps,1);
psi_hat_hist = zeros(nTimeSteps,1);
r_hat_hist   = zeros(nTimeSteps,1);
b_hat_hist   = zeros(nTimeSteps,1);

% === Continous model ===
A = [0  1     0;
     0 -1/T_nom  -K_nom/T_nom;
     0  0     0];
B = [0; K_nom/T_nom; 0];
C = [1 0 0];
E = [0 0; 1 0; 0 1];
D = 0;

% First-order discretization
Ad = eye(3) + h*A;
Bd = h*B;
Cd = C;
Ed = h*E;

% KF initialization
x_prd = [0;0;0];            % [psi_hat; r_hat; b_hat]
P_prd = diag([(deg2rad(30))^2, (deg2rad(0.01))^2, (deg2rad(1))^2]);

% Tuning (starting point – adjust during tests)
Qd = diag([1e-11, 1e-7]);    % var{w_r}, var{w_b}
Rd = sigma_psi^2;           % var of psi measurement

% Guidance model initialization
% --- Waypoints ---
S = load('WP.mat');                 
WP = S.WP;

Delta_h = 600; % Lookahead distance
Rsw = 500; % Switch radius
Rstop = 100; % End radius

kappa = 1; % Design parameter for ILOS

i_end = nTimeSteps;

for i = 1:nTimeSteps
    % --- Time-varying heading reference ---
    %{
    if t(i) < 500
        psi_ref = 10 * pi/180;     % +10 deg
    else
        psi_ref = -20 * pi/180;    % -20 deg
    end
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1a) 2D irrotational current  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vc = 1;
    %Vc = 0;                         % (m/s)
    betaVc = 45*pi/180;             % 45 deg CW from N -> NE (rad)
    
    uc = Vc * cos(betaVc - x(6));   % BODY Surge current (m/s)
    vc = Vc * sin(betaVc - x(6));   % BODY Sway current (m/s)
    nu_c = [ uc vc 0 ]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1c) Add wind here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vw = 10;                % wind speed (m/s)
    beta_Vw = deg2rad(135); % wind speed direction
    rho_a = 1.247;          % air density
    c_y = 0.95;             % sway wind-force coefficient
    c_n = 0.15;             % Yaw wind-moment coefficient
    L = 161;                % Length of ship
    A_Lw = 10*L;            % Area of side view above water.
    
    % Tilbakemelding fra studass
    uw = Vw * cos(beta_Vw - x(6));
    vw = Vw * sin(beta_Vw - x(6));
    u_rw = x(1) - uw;
    v_rw = x(2) - vw;
    V_rw = sqrt(u_rw^2 + v_rw^2);
    
    gamma_rw = -atan2(v_rw,u_rw);
    Ywind = 0.5 * rho_a * V_rw^2 * c_y*sin(gamma_rw) * A_Lw;       % Equation from Fossen ch 10.1
    Nwind = 0.5 * rho_a * V_rw^2 * c_n*sin(2*gamma_rw) * A_Lw*L;   % Equation from Fossen ch 10.1
    
    tau_wind = [0 Ywind Nwind]';

    % Relative speed ocean currents
    u_rc = x(1) - uc;
    v_rc = x(2) - vc;
    U_rc = sqrt(u_rc^2 + v_rc^2);

    % Defining sidelsip and crab angle
    beta_c   = atan2(x(2), max(1e-9, x(1))); % Crab angle

    % Course: chi = psi + beta_c  (gjelder eksakt når Vc=0)
    psi    = x(6);
    chi    = ssa(psi + beta_c);

    % Sideslip angle
    beta = atan2(v_rc, max(1e-9, u_rc));

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3 - Task 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guidance law
    [xk1,yk1,xk,yk,last] = WP_selector(x(4),x(5), WP, Rsw, Rstop);
    [e_y,pi_p] = crossTrackError(xk1,yk1,xk,yk,x(4),x(5));
    chi_d = LOS_guidance(e_y,pi_p, Delta_h);
    %psi_ref = chi_d - beta_c;
    psi_ref = ILOS_guidance(e_y, pi_p, kappa, Delta_h, h);
    xd_dot = ref_model(xd, psi_ref);
    xd = xd + h * xd_dot;
    psi_d = xd(1);
    r_d   = xd(2);
    r_d_dot = xd_dot(2);
    u_d = U_ref;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add the heading controller here 
    % Define it as a function
    %
    % The result should look like this:
    % delta_c = PID_heading(e_psi,e_r,e_int);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % ----- Noisy measurements (from true states) -----
    psi_true = ssa(x(6));          % x = [u v r x y psi delta n Qm]'
    r_true   = x(3);
    
    psi_meas = ssa(psi_true + randn*sigma_psi );
    r_meas   = r_true   + randn*sigma_r;   % not used by KF, but useful for plots
    
    % Log for plotting
    psi_meas_hist(i) = psi_meas;
    r_meas_hist(i)   = r_meas;
    psi_true_hist(i) = psi_true;
    r_true_hist(i)   = r_true;

    % Run the kalman filter

    [x_pst,P_pst,x_prd,P_prd] = KF(x_prd,P_prd,Ad,Bd,Ed,Cd,Qd,Rd,psi_meas,x(7));

    % Log for plotting
    % Logg estimater
    psi_hat_hist(i) = x_pst(1);
    r_hat_hist(i)   = x_pst(2);
    b_hat_hist(i)   = x_pst(3);


    if USE_KF
        psi_fb = x_pst(1);     % estimated yaw
        r_fb   = x_pst(2);     % estimated yaw rate
    else
        %psi_fb = psi_meas; % noisy yaw
        %r_fb = r_meas;     % noisy yaw rate
        psi_fb = x(6);     % noisy yaw
        r_fb   = x(3);       % noisy yaw rate
    end

    e_psi = ssa(psi_fb - psi_d);
    e_r   = r_fb - r_d;
    e_u   = x(1) - u_d;
    delta_unsat = -(kp*e_psi + kd*e_r + ki*e_int);

    %delta_step = delta_unsat - delta_cmd;
    %delta_step = max(-Ddelta_max*h, min(Ddelta_max*h, delta_step));
    %delta_cmd  = delta_cmd + delta_step;

    % Saturation
    delta_c = min(max(delta_unsat, -delta_max), delta_max);
    %delta_c = max(-delta_max, min(delta_max, delta_cmd));

    % Anti-windup
    e_int_dot = e_psi - (1/ki)*(delta_c - delta_unsat);
    e_int = e_int + h * e_int_dot;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3e) Add open loop speed control here
    % Define it as a function
    %
    % The result should look like this:
    % n_c = open_loop_speed_control(U_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_c = open_loop_speed_control(U_ref);
    %n_c = 10;                   % propeller speed [radians per second (rps)]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3f) Replace the open loop speed controller, 
    % with a closed loop speed controller here 
    % Define it as a function
    %
    % The result should look like this:
    % n_c = closed_loop_speed_control(u_d,e_u,e_int_u);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %n_c = closed_loop_speed_control(u_d, x(1), h);
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,tau_total] = ship(x,u,nu_c,tau_wind);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d, chi, chi_d, beta_c, beta];     
 
    % Euler integration
    % x = euler2(xdot,x,h); 
    % Runge Kutta 4 integration
    x = rk4(@ship,h,x,u,nu_c,tau_wind);

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end

    if last
        i_end = i;  % siste gyldige indeks i denne simuleringen

        % (valgfritt) Logg en siste rad så plott ikke blir tomt helt på slutten
        % Her logger vi med "null" kommandoer for tydelig slutt
        psi_d = x(6); r_d = 0; u_d = U_ref;   % beholder samme referanser
        simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d, chi, chi_d, beta_c, beta];  

        fprintf('  Reached final waypoint at t = %.1f s (step %d)\n', (i-1)*h, i);
        break
    end

end
simdata = simdata(1:i,:);
t = t(1:i);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
x           = simdata(:,4);                 % m
y           = simdata(:,5);                 % m
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c         = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         = simdata(:,13);               % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
chi         = simdata(:,14);                % rad
chi_deg     = (180/pi)* chi;                % deg
chi_d       = simdata(:, 15);               % rad
chi_d_deg   = (180/pi)*chi_d;               % deg
beta_c      = simdata(:, 16);               % rad
beta_c_deg  = (180/pi)*beta_c;               % deg
beta        = simdata(:, 17);               % rad
beta_deg    = (180/pi)*beta;                % deg

psi_meas_hist = psi_meas_hist(1:length(t));
r_meas_hist   = r_meas_hist(1:length(t));
psi_true_hist = psi_true_hist(1:length(t));
r_true_hist   = r_true_hist(1:length(t));
psi_hat_hist = psi_hat_hist(1:length(t));
r_hat_hist   = r_hat_hist(1:length(t));
b_hat_hist   = b_hat_hist(1:length(t));

%%
figure(3)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi_deg,t,psi_d_deg,'linewidth',2);
title('Actual and desired yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 
legend('actual yaw','desired yaw')
subplot(313)
plot(t,r_deg,t,r_d_deg,'linewidth',2);
title('Actual and desired yaw rates'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 
legend('actual yaw rate','desired yaw rate')

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')

figure(4)
figure(gcf)
plot(t, chi_deg,   'LineWidth', 2); hold on;
plot(t, chi_d_deg, 'LineWidth', 2);
plot(t, psi_deg,   'LineWidth', 2);
plot(t, beta_c_deg,'--',        'LineWidth', 1.5);
plot(t, beta_deg,  '--',        'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Angle (deg)');
title('Course (χ), Desired Course (χ_d), Heading (ψ), Crab (β_c), Sideslip (β)');
legend('\chi','\chi_d','\psi','\beta_c','\beta','Location','best');

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


%% Create objects for 3-D visualization 
% Since we only simulate 3-DOF we need to construct zero arrays for the 
% excluded dimensions, including height, roll and pitch
z = zeros(length(x),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);

% create object 1: ship (ship1.mat)
new_object('flypath3d_v2/ship1.mat',[x,y,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(x)),max(abs(y)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);

% Plot trajectories
%{
flypath('flypath3d_v2/ship1.mat',...
'animate','on','step',500,...
'axis','on','axiscolor',[0 0 0],'color',[1 1 1],...
'font','Georgia','fontsize',12,...
'view',[-25 35],'window',[900 900],...
'xlim', [min(y)-0.1*max(abs(y)),max(y)+0.1*max(abs(y))],... 
'ylim', [min(x)-0.1*max(abs(x)),max(x)+0.1*max(abs(x))], ...
'zlim', [-max(max(abs(x)),max(abs(y)))/100,max(max(abs(x)),max(abs(y)))/20]);
%}

% Plot the path
pathplotter(x, y);