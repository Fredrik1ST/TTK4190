function [e_y, pi_p] = crossTrackError(xk1, yk1, xk, yk, xE, yN)
% e_y  : tverr-avvik (pos) til linjesegmentet (positiv til venstre for ruten)
% pi_p : tangentretn. (coursen) til linjesegmentet

    dx   = xk1 - xk; 
    dy   = yk1 - yk;
    L    = hypot(dx, dy) + 1e-9;
    pi_p = atan2(dy, dx);  % [rad], N->E konvensjon

    % Rotér posfeil inn i bane-ramme (t: langs, n: tverrsnitt)
    R = [ cos(pi_p)  sin(pi_p);
         -sin(pi_p)  cos(pi_p)];
    en = R * ([xE; yN] - [xk; yk]);  % [e_t; e_n]
    e_y = en(2);
end
