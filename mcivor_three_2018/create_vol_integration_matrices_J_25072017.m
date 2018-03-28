function [vol_int_J_z, vol_int_J_rt]=create_vol_integration_matrices_J_25072017(r_a,dt,val_Dm,mu_max_m,beta_max_m,n_max_m,dr_m,dz_m,dtheta,theta_prime,r_m_prime,z_m_prime)

    
    trap_z_junction=dz_m*[0.5 ones(1,length(z_m_prime)-2) 0.5];
    trap_r_junction=dr_m*r_m_prime.*[0.5 ones(1,length(r_m_prime)-2) 0.5];
    trap_theta_junction=dtheta*[0.5 ones(1,length(theta_prime)-2) 0.5];
    trap_r_theta_junction=kron(trap_theta_junction,trap_r_junction);

fn_greens=['./code_matrices/Greens_function_junction-dt_',num2str(dt),'-Dm_',num2str(val_Dm),'-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'.mat'];
if exist('GF_junction_rt_rptp','var')
        vol_int_J_z=GF_junction_z_zp*diag(trap_z_junction);
    vol_int_J_rt=(GF_junction_rt_rptp*diag(trap_r_theta_junction))';
else
load(fn_greens)
    
    vol_int_J_z=GF_junction_z_zp*diag(trap_z_junction);
    vol_int_J_rt=(GF_junction_rt_rptp*diag(trap_r_theta_junction))';
end

    % make sure that all primed variables are multiplied together so we
    % integrate out the primed variables

end
