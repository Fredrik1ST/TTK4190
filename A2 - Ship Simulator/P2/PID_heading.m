function [delta_c, e_int] = PID_heading(e_psi, e_r, e_int)
    % PID heading controller with anti-windup
    P = evalin('base','PIDH_PERSIST');

    % Compute the unsaturated control signal
    delta_unsat = -(P.kp*e_psi + P.kd*e_r + P.ki*e_int);

    delta_c = delta_unsat;
    e_int = e_int + P.h*e_psi;

    % Saturated control law and integrator anti-windup
    if delta_unsat > P.delta_max
        delta_c = P.delta_max; % Saturation
        e_int = e_int - (P.h/P.ki) * (delta_c - delta_unsat); % Anti-windup

    elseif delta_unsat < -P.delta_max
        delta_c = -P.delta_max; % Saturation
        e_int = e_int - (P.h/P.ki) * (delta_c - delta_unsat); % Anti-windup
    end

    % Save state
    P.e_int = e_int;
    assignin('base','PIDH_PERSIST', P);
end