function [vol_int_S_z, vol_int_S_rt]=create_vol_integration_matrices_S_10082017(r_a,dt,val_De,mu_max_s,beta_max_s,n_max_s,dr_s,dz_s,dtheta,theta_prime,r_s_prime,z_s_prime,L1)
    
    trap_z_sub=dz_s*[0.5 ones(1,length(z_s_prime)-2) 0.5];
    trap_r_sub=dr_s*r_s_prime.*[0.5 ones(1,length(r_s_prime)-2) 0.5];
    trap_theta_sub=dtheta*[0.5 ones(1,length(theta_prime)-2) 0.5];
    trap_r_theta_sub=kron(trap_theta_sub,trap_r_sub);
    
    fn_greens=['./code_matrices/Greens_function_sub-L1_',num2str(L1),'-dt_',num2str(dt),'-De_',num2str(val_De),'-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s)...
    ,'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'.mat'];
    
if exist('GF_sub_rt_rptp','var')
        vol_int_S_z=GF_sub_z_zp*diag(trap_z_sub);
    vol_int_S_rt=(GF_sub_rt_rptp*diag(trap_r_theta_sub))';
else
load(fn_greens)
    
    vol_int_S_z=GF_sub_z_zp*diag(trap_z_sub);
    vol_int_S_rt=(GF_sub_rt_rptp*diag(trap_r_theta_sub))';
end
    % make sure that all primed variables are multiplied together so we
    % integrate out the primed variables

end
