% Line of sight guidance law
function chi_d = LOS_guidance(e_y, pi_p, Delta_h)
    
chi_d = pi_p - atan(e_y / max(Delta_h,1e-6)); % Eq. 12.38 from Fossen 2021
chi_d = ssa(chi_d);

end