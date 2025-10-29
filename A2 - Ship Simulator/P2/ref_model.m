function xd_dot = ref_model(xd, psi_ref)
% REF_MODEL - 3rd order reference model from Fossen (2021, eq. 15.143)
%
% State vector:
%   xd = [psi_d; r_d; v_d]
%
% Inputs:
%   psi_ref  - reference heading (rad)
%
% Parameters:
%   wb_ref   - natural frequency (rad/s)
%   zeta_ref - damping ratio
%
% Output:
%   xd_dot = [psi_d_dot; r_d_dot; v_d_dot]

    % --- Model parameters ---
    wb_ref   = 0.03;   % natural frequency (rad/s)
    zeta_ref = 1.0;    % damping ratio

    % --- State variables ---
    psi_d = xd(1);
    r_d   = xd(2);
    v_d   = xd(3);

    % --- Coefficients from Fossen (15.143) ---
    a2 = (2*zeta_ref + 1) * wb_ref;
    a1 = (1 + 2*zeta_ref) * wb_ref^2;
    a0 = wb_ref^3;
    g  = wb_ref^3;

    % --- State derivatives ---
    psi_d_dot = r_d;
    r_d_dot   = v_d;
    v_d_dot   = -a2*v_d - a1*r_d - a0*psi_d + g*psi_ref;

    % --- Output ---
    xd_dot = [psi_d_dot; r_d_dot; v_d_dot];
end
