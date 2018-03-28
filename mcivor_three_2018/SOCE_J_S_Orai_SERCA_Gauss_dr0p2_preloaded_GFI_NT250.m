
tic
   %%%% T= 1ms %%%%
    T_ms=length_t*dt*1e3
   n_timesteps=20;
 
%% SERCA flux parameters
    %SERCA pump (bidirectional) constants
    Kmr_Bers1998=1700; %micro molar
    Q_Shannon2004=2.6; %Temperature factor
  
%% load SERCA fluxes
  tic   
    tts=0:5:45;
    theta_SERCA_all=1+phi_mult*tts; 
    r_SERCA_all=(r_SERCA/dr_m)*ones(1,length(theta_SERCA_all)); 
    
    flux_SERCA_junction_all=zeros(length(z_m),length(r_m)*length(theta),length(theta_SERCA_all));
    flux_SERCA_subPMER_all=zeros(length(z_s),length(r_s)*length(theta),length(theta_SERCA_all));

    for ii=1:length(theta_SERCA_all)
        r_SERCA_index=r_SERCA_all(ii);
        theta_SERCA_index=theta_SERCA_all(ii);

    fn_particular_soln=['./code_matrices/Boundary_solution_z_L2_junction-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
    ,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'-r_serca_index_',num2str(r_SERCA_index),'-theta_serca_index_',num2str(theta_SERCA_index),'-sr_',num2str(sr),'_st_',num2str(st),'.mat'];
    load(fn_particular_soln,'Boundary_solution_junction_z_L2_rtz')
    Boundary_solution_junction_z_L2_rt_z=reshape(Boundary_solution_junction_z_L2_rtz,length(r_m)*length(theta),length(z_m));
    flux_SERCA_junction_all(:,:,ii)=volume_integral_junction(vol_int_J_z_NT250,vol_int_J_rt_NT250,0,Boundary_solution_junction_z_L2_rt_z')*vol_cyt_nm/vol_cyt_L/Dm;
   
    fn_particular_soln=['./code_matrices/Boundary_solution_z_L2_sub-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s)...
    ,'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'-r_SERCA_index_',num2str(r_SERCA_index),'-theta_SERCA_index_',num2str(theta_SERCA_index),'-sr_',num2str(sr),'_st_',num2str(st),'.mat'];
    load(fn_particular_soln,'Boundary_solution_subPMER_z_L2_rtz')
    Boundary_solution_subPMER_z_L2_rt_z=reshape(Boundary_solution_subPMER_z_L2_rtz,length(r_s)*length(theta),length(z_s));
    flux_SERCA_subPMER_all(:,:,ii)=volume_integral_S_flux_BC(vol_int_S_z_NT250,vol_int_S_rt_NT250,0,Boundary_solution_subPMER_z_L2_rt_z',z_s,r_s,theta)*vol_ER_nm/vol_ER_L/De;
   
    end
    
    disp('loaded all SERCA fluxes')
toc

%% Moving time forwards

    ER_PM_junction_soln=zeros(length(r_m)*length(theta),length(z_m),n_timesteps);
    Sub_PM_ER_soln=zeros(length(r_s)*length(theta),length(z_s),n_timesteps);

    q=repmat(1:n_timesteps,1,round(length_t/n_timesteps));
    P=round(length_t/n_timesteps);  

    
    %create folder to save simulation data in (scratch drive)
    if n_H==n_H_SERCA2b
    cd_data_save=['./simulations/SERCA2b_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/'];
    elseif n_H==n_H_SERCA2a
    cd_data_save=['./simulations/SERCA2a_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/'];
    end    
    if exist(cd_data_save,'dir')
    else
        mkdir(cd_data_save);
    end
    
    SERCA_activity=zeros(length(theta_SERCA_all),length_t);
    
    area_SERCA_grid=r_SERCA*dtheta*dr_m;

    tic
    for p=1:P
    if n_H==n_H_SERCA2b
    fn_data=['./simulations/SERCA2b_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/SOCE-p_',num2str(p),'.mat'];
    elseif n_H==n_H_SERCA2a
    fn_data=['./simulations/SERCA2a_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/SOCE-p_',num2str(p),'.mat'];
    end    

    tic
    for i=(p-1)*n_timesteps+1:p*n_timesteps
        

        % calculate SERCA flux magnitude at each time step dt for each
        % SERCA pump
        flux_SERCA_junction=zeros(length(z_m),length(r_m)*length(theta));
        flux_SERCA_S=zeros(length(z_s),length(r_s)*length(theta));

        for ii=1:length(theta_SERCA_all)
        theta_SERCA_index=theta_SERCA_all(ii);
        r_SERCA_index=r_SERCA_all(ii);

        J_SERCA_num=(Junction_conc(r_SERCA_index,theta_SERCA_index,1)/Kmf)^n_H - (Sub_conc(r_SERCA_index,theta_SERCA_index,end)/Kmr_Bers1998)^n_H;
        J_SERCA_den=1+(Junction_conc(r_SERCA_index,theta_SERCA_index,1)/Kmf)^n_H +(Sub_conc(r_SERCA_index,theta_SERCA_index,end)/Kmr_Bers1998)^n_H;
        J_SERCA_flux_magnitude_shannon=Q_Shannon2004*Vmax*J_SERCA_num/J_SERCA_den;
        J_SERCA_flux_per_unit_area=J_SERCA_flux_magnitude_shannon/area_SERCA_grid;

        flux_SERCA_junction=flux_SERCA_junction + J_SERCA_flux_per_unit_area*flux_SERCA_junction_all(:,:,ii);
        flux_SERCA_S=flux_SERCA_S + J_SERCA_flux_per_unit_area*flux_SERCA_subPMER_all(:,:,ii);

        SERCA_activity(ii,i)=J_SERCA_num/J_SERCA_den;
        end        

        Micro_conc_IC=reshape(Junction_conc,length(r_m)*length(theta),length(z_m))';
        Sub_conc_IC=reshape(Sub_conc,length(r_s)*length(theta),length(z_s))';

        J_int=volume_integral_junction(vol_int_J_z_NT250,vol_int_J_rt_NT250,Micro_conc_IC,J_ra);
        J=J_int+J_CRAC_flux_int+flux_SERCA_junction;
                
        ER_PM_junction_soln(:,:,q(i))=J';
        
        S_int=volume_integral_S(vol_int_S_z_NT250,vol_int_S_rt_NT250,Sub_conc_IC,S_zL1,r_s,theta,z_s);
        Sub_int=S_int+flux_SERCA_S;
        
        Sub_PM_ER_soln(:,:,q(i))=Sub_int';

        % Reset initial conditions
        Junction_conc=reshape(J',length(r_m),length(theta),length(z_m));
        Sub_conc=reshape(Sub_int',length(r_s),length(theta),length(z_s));
    
    end
  
    if p==1
    save(fn_data,'ER_PM_junction_soln','Sub_PM_ER_soln','-v7.3')
    p
    elseif p==25
    save(fn_data,'ER_PM_junction_soln','Sub_PM_ER_soln','-v7.3')
    p
    elseif p==50
    save(fn_data,'ER_PM_junction_soln','Sub_PM_ER_soln','-v7.3')
    p
    end

    if n_H==n_H_SERCA2b
    fn_data_serca=['./simulations/SERCA2b_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/SERCA_activity.mat'];
    elseif n_H==n_H_SERCA2a
    fn_data_serca=['./simulations/SERCA2a_ICRAC_',num2str(I_CRAC),'-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_',num2str(r_CRAC_all(1)),'-',num2str(length(theta_CRAC_all)),'CRACs-r_SERCA_',num2str(r_SERCA_all(1)),'-theta_SERCA_',num2str(theta_SERCA_all(1)),'-',num2str(length(theta_SERCA_all)),'_SERCAs-Vmax_',num2str(Vmax),'_n_max_',num2str(n_max_m),'_Dm',num2str(val_Dm),'_De_',num2str(val_De),'_dt_',num2str(dt),'_dr_',num2str(dr_m),'_dphi_',num2str(dtheta),'_dz_',num2str(dz_m),'r_a_',num2str(r_a),'_T_',num2str(length_t),'-sr_',num2str(sr),'-st_',num2str(st),'/SERCA_activity.mat'];
    end    

save(fn_data_serca,'SERCA_activity')

    toc
    end
toc
