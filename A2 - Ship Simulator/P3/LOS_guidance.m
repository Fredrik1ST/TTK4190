function chi_d = LOS_guidance(e_y, pi_p)

delta_h = 600; % Lookahead distance
    
chi_d = pi_p - atan(e_y / max(delta_h,1e-6));
chi_d = ssa(chi_d);

end