% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions using a MAIN LOOP with fixed step. 
% Euler's method, RK4, and exact discetization are investigated.
% 
% The MSS m-files must be on your Matlab path in order to run the script.
% See How to install Matlab for MSS: https://github.com/cybergalactic/MSS/
%
% System:                      
%                            q_dot = T(q)w
%                             
%                            Ib * w_dot - S(Ib * w) * w = tau
% Control law:
%                            tau = constant (3x1)
% 
% Definitions:             
%                            Ib : inertia matrix (3x3)
%                            S(w) : Skew-symmetric matrix (3x3)
%                            T(q) : transformation matrix (4x3)
%                            tau : Control input (3x1)
%                            w : Angular velocity vector (3x1)
%                            q : unit quaternion vector (4x1)

%% USER INPUTS
clc; clear;

T_final = 400;	             % Final simulation time (s)
h = 0.1;                     % Sampling time (s)

% Model parameters
m = 10;                      % Mass
R = 0.1 / sqrt(6);           % Radius of gyration
Ib = m * R^2 * eye(3);       % Inertia matrix in CO
I_inv = inv(Ib);

% Initial states
phi = -deg2rad(10);          % Initial Euler angles
theta = deg2rad(10);
psi = deg2rad(5);

q = euler2q(phi,theta,psi);  % Transform initial Euler angles to q

w = [0 0 0]';                % Initial angular rates

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

%% MAIN LOOP
simdata = zeros(nTimeSteps, 13); % Pre-allocate table for simdata

for i = 1:nTimeSteps

   % Control law
   tau = 1e-4*[0.5 1 -1]';      

   [phi,theta,psi] = q2euler(q); % Transform q to Euler angles
   
   % Store data for presentation
   simdata(i,:) = [q' phi theta psi w' tau'];  % Store data in table
   
   % Numerical integration q(k+1)
   % q_dot = Rquat(q) * w;          % Quaternion kinematics
   % q = q + h * q_dot;	            % Euler's method (not reccomended)
   q = expm( Tquat(w) * h ) * q;    % Exact discretization
   q  = q / norm(q);                % Unit quaternion normalization

   % Numerical integration w(k+1)
   % w_dot = I_inv * (Smtrx(Ib * w) * w + tau); % Rigid-body kinetics
   % w = w + h * w_dot;              % Euler's method (not recommended)
   function_w_dot = @(w, I_inv, Ib, tau) I_inv * (Smtrx(Ib * w) * w + tau); 
   w = rk4(function_w_dot, h, w, I_inv, Ib, tau);

end 

%% PLOTS 
close all;
q       = simdata(:,1:4); 
phi     = rad2deg(simdata(:,5));
theta   = rad2deg(simdata(:,6));
psi     = rad2deg(simdata(:,7));
w       = rad2deg(simdata(:,8:10));  
tau     = simdata(:,11:13);


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('p', 'q', 'r');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)