function n_c = open_loop_speed_control(U_ref)
% Maps desired surge speed U_ref [m/s] -> commanded shaft speed n_c [rpm]
% using cruise-speed relation and Wageningen KT at Ja = 0.

    % --- Surge damping (same calc as in ship.m) ---
    Xudot = -8.9830e5;
    m     = 17.0677e6;  
    T1    = 20;
    Xu    = -(m - Xudot)/T1;      % [N·s/m] (negative)

    % --- Wageningen at J_a = 0 ---
    KT0 = 0.6367;  % KT(0)

    % --- Ship/prop constants (match ship.m) ---
    rho   = 1025;   % [kg/m^3]
    D     = 3.3;    % [m]
    t_thr = 0.05;   % thrust deduction number used in ship.m

    % --- Desired steady-state thrust (cruise condition) ---
    Td = (Xu / (t_thr - 1)) * U_ref;          % [N]

    % --- Output in rpm (ship.m expects rpm) ---
    n_c = sign(Td) * sqrt( max(0,abs(Td)) / (rho * D^4 * KT0) );
end
