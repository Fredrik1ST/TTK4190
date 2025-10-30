function [delta_c] = PID_heading(e_psi, e_r, e_int)
    P = evalin('base','PIDH_PERSIST');  % kan beholdes om du vil hente PID-parametere herfra

    % Compute unsaturated control
    delta_unsat = -(P.kp*e_psi + P.kd*e_r + P.ki*e_int);

    % Saturation
    delta_c = min(max(delta_unsat, -P.delta_max), P.delta_max);

    % Anti-windup
    e_int_dot = e_psi - (1/P.ki)*(delta_c - delta_unsat);
    e_int = e_int + P.h * e_int_dot;
end
