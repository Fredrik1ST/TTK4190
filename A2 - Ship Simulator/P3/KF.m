function [x_pst,P_pst,x_prd,P_prd] = KF(x_prd,P_prd,Ad,Bd,Ed,Cd,Qd,Rd,psi_meas, delta_c)

% Discrete-time Kalman filter for x = [psi; r; b]
% Inputs:
%   x_prd,P_prd : prior state and covariance
%   Ad,Bd,Ed,Cd : discrete model matrices
%   Qd,Rd       : process/measurement noise covariances
%   psi_meas    : yaw angle measurement (rad)
%   delta       : rudder input (rad)
% Outputs:
%   x_pst,P_pst : posterior state and covariance
%   x_prd,P_prd : predicted state and covariance for next step

% --- Corrector ---
Kk    = P_prd*Cd' / (Cd*P_prd*Cd' + Rd);
innov = ssa( psi_meas - Cd*x_prd );  % wrap angle innovation
x_pst = ssa(x_prd + Kk*innov);
IkC   = eye(size(P_prd)) - Kk*Cd;
P_pst = IkC*P_prd*IkC' + Kk*Rd*Kk';

% --- Predictor ---
u     = delta_c;                        % rudder input to Nomoto model
x_prd = ssa(Ad*x_pst + Bd*u);
P_prd = Ad*P_pst*Ad' + Ed*Qd*Ed';

end