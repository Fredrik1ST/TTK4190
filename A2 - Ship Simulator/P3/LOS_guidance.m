function chi_d = LOS_guidance(e_y, pi_p, Delta_h)
    
chi_d = pi_p - atan(e_y / max(Delta_h,1e-6));
chi_d = ssa(chi_d);

end