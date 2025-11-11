% Integral Line-of-sight guidance law
function psi_ref = ILOS_guidance(e_y, pi_h, kappa, Delta_h, h)

% Ensure the integral state is updated correctly
persistent e_y_int;    % integral state

if isempty(e_y_int)
    e_y_int = 0;
end

% ILOS guidance law
Kp = 1 / Delta_h;
Ki = kappa*Kp;
%psi_ref = pi_h - atan( Kp * (e_y + kappa * e_y_int) ); From eq 12.108 in Fossen 2021   
psi_ref = pi_h - atan( Kp * e_y + Ki * e_y_int); % From eq 12.108 in Fossen 2021


% Propagation of states to time k+1
e_y_int = e_y_int + h * Delta_h * e_y/(Delta_h^2 + (e_y + kappa * e_y_int)^2); % From eq 12.109 in Fossen 2021

end