function J_int=volume_integral_junction(vol_int_J_z,vol_int_J_rt,Junction_conc_IC,particular_soln_J)

J_int=vol_int_J_z*(Junction_conc_IC-particular_soln_J)*vol_int_J_rt + particular_soln_J;

% m=Micro_conc_IC-particular_soln_M;
% minM=min(m(end,:))
end