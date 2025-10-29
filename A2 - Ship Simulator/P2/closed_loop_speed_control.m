function n_c = closed_loop_speed_control(U_ref, U, dt)
% CLOSED_LOOP_SPEED_CONTROL
% Implements a PI-based closed-loop speed controller for a ship model.
%
%   n_c = closed_loop_speed_control(U_ref, U, dt)
%
% Inputs:
%   U_ref : desired surge speed [m/s]
%   U     : measured surge speed [m/s]
%   dt    : sampling time [s]
%
% Output:
%   n_c   : commanded shaft speed [rpm]
%
% The controller uses an open-loop feedforward term based on the
% Wageningen KT relation and a PI feedback term to correct speed error.
%
% Structure:
%   1. Compute speed error
%   2. Integrate error (PI control)
%   3. Combine open-loop and feedback terms
%   4. Apply output saturation limits

    % --- Persistent integrator for accumulated error ---
    persistent e_int ,
    if isempty(e_int)
        e_int = 0;
    end

    % --- Compute speed error ---
    e = U_ref - U;

    % --- PI controller gains (tuning required for your vessel) ---
    Kp = 30;   % proportional gain
    Ki = 5;    % integral gain

    % --- PI controller output (feedback correction) ---
    n_c_fb = Kp * e + Ki * e_int;

    % --- Feedforward term from open-loop speed control ---
    n_c_ff = open_loop_speed_control(U_ref);

    % --- Combine feedforward and feedback commands ---
    n_c = n_c_ff + n_c_fb;

    % --- Limit commanded RPM to safe range ---
    n_c = max(min(n_c, 200), -200);

    % --- Integrate error over time ---
    e_int = e_int + e * dt;
end