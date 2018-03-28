tic
disp('preload GF check')
    
%% Geometry of domains
    eps=1e-16;
    r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=2000;
    L0=10;
    
    %Volume of microdomain
    vol_junction_nm=pi*r_a^2*(H-L2); %nm^3
    vol_junction_m=vol_junction_nm*1e-27; %m^3
    vol_junction_L=vol_junction_m*1e3; %volume in litres (L)
   
    %Volume of Sub-PM ER
    vol_subPMER_nm=pi*r_a^2*(L2-L1); %nm^3
    vol_subPMER_m=vol_subPMER_nm*1e-27; %m^3
    vol_subPMER_L=vol_subPMER_m*1e3; %volume in litres (L)
    
    %Volume of cytosol (Berlin, Bassani, Bers 1994)
    vol_cyt_L=24.5e-12; %(L) volume of cytosol of ventricular myocyte in litres
    vol_cyt_m=vol_cyt_L*1e-3;    %(m^3) volume of cytosol of ventricular myocyte in metres
    vol_cyt_nm=vol_cyt_m*1e27;
    
    %Volume of Bulk ER
    vol_bulk_ER_nm=pi*r_b^2*(L1-L0); %nm^3
    vol_bulk_ER_m=vol_bulk_ER_nm*1e-27; %m^3
    vol_bulk_ER_L=vol_bulk_ER_m*1e3; %volume in litres (L)

    %Volume of SR in Shannon and Bers (2004) 
    vol_SR_Shannon_Bers_2004=3.5*33e-12/100; %volume in litres
    
    %Volume of ER in Alberst table 12-1: 6 times smaller than cytosolic
    %cell volume
    vol_ER_L=vol_cyt_L/6; %(L) volume of cytosol of ventricular myocyte in litres
    vol_ER_nm=vol_cyt_nm/6;
   
%% estimated area of Orai channel 
area_orai_channel=pi*(.55/2)^2; % max diameter Orai=0.55nm, area = pi*R^2

%% ER-PM micro-domain 

%Geometry specific to micro-domain
%     dr_m=r_a/NR;      %gives 75 mesh points in r array
    dz_m=0.15;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array
    
    NR=r_a/dr_m;
    NZ=(H-L2)/dz_m;
    
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    
    r_m=dr_m:dr_m:r_a;
    r_m_prime=r_m;
    
    z_m=L2:dz_m:H;
    z_m_prime=z_m;
% 
%     [THETA,RHO] = meshgrid(theta,r_m);
%     [Xmicro,Ymicro] = pol2cart(THETA,RHO);
%     [zzmicro,rrmicro]=meshgrid(r_m,z_m);

    mu_max_m=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_m=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_m=NT/2;            % maximal number of Bessel functions (n values)   

    % Initial conditions in ER-PM micro-domain
    rest_conc_J=0.1; %micro moles per litre cytosol
    Junction_conc=rest_conc_J*ones(length(r_m),length(theta),length(z_m));
        
    if exist('vol_int_J_z_NT250')
    disp('GF J integrators preloaded')
    else
    [vol_int_J_z_NT250, vol_int_J_rt_NT250]=create_vol_integration_matrices_J_25072017(r_a,dt,val_Dm,mu_max_m,beta_max_m,n_max_m,dr_m,dz_m,dtheta,theta_prime,r_m_prime,z_m_prime);
        disp('created GF J integrators')
    end
    
%% SubPM ER domain
% Geometry specific to subPM ER
%     dr_s=r_a/NR;      %gives 75 mesh points in r array
    dz_s=4.85;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array
    
    NR=r_a/dr_s;
    NZ=(L2-L1)/dz_s;
    
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    
    r_s=dr_s:dr_s:r_a;
    r_s_prime=r_s;
    
    z_s=L1:dz_s:L2;
    z_s_prime=z_s;

    [THETA,RHO] = meshgrid(theta,r_s);
    [Xsub,Ysub] = pol2cart(THETA,RHO);
    [zzsub,rrsub]=meshgrid(r_s,z_s);    
    
    mu_max_s=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_s=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_s=NT/2;            % maximal number of Bessel functions (n values)   
    
    % Initial conditions in ER-PM micro-domain
    rest_conc_S=150; %micro moles per litre cytosol
    Sub_conc=rest_conc_S*ones(length(r_s),length(theta),length(z_s));
    
    % Boundary conditions in ER-PM micro-domain
    S_zL1=rest_conc_S;
    
    if exist('vol_int_S_z_NT250')
    disp('GF S integrators preloaded')
    else
    [vol_int_S_z_NT250, vol_int_S_rt_NT250]=create_vol_integration_matrices_S_10082017(r_a,dt,val_De,mu_max_s,beta_max_s,n_max_s,dr_s,dz_s,dtheta,theta_prime,r_s_prime,z_s_prime,L1);
        disp('created GF S integrators')
    end
 toc
%% load Orai channel fluxes
    %% Orai flux parameters
    % conversion constant to change SERCA strength from Shannon's model to
    % ours:
    SERCA_strength_conversion=vol_subPMER_L/vol_SR_Shannon_Bers_2004;

    %CRAC channel flux
    F=96485.3365; %Faraday's constant (C/mol)
    I_CRAC=2.1e-15; %Current through CRAC channel (C/s)
    z_val=2; %Valency of calcium ions (dimensionless)
    
    f_CRAC=I_CRAC/(F*z_val); %moles of calcium entering microdomain per second (mol/s) through one Orai channel
    f_CRAC_micro_mol=f_CRAC*1e6;    %micro moles per second per unit area
    
    J_CRAC=f_CRAC_micro_mol*vol_cyt_nm/vol_cyt_L;
    J_ra=rest_conc_J;
       
    J_CRAC_per_unit_area=J_CRAC/area_orai_channel;
    flux_CRAC_magnitude=J_CRAC_per_unit_area/Dm;

tic
     ttc=5:10:45;
    theta_CRAC_all=1+phi_mult*ttc; 
    r_CRAC_all=(r_CRAC/dr_m)*ones(1,length(theta_CRAC_all)); 

   particular_solution_z_H_all=zeros(length(r_m)*length(theta),length(z_m));
   for ii=1:length(theta_CRAC_all)
       r_CRAC_index=r_CRAC_all(ii);
       theta_CRAC_index=theta_CRAC_all(ii);

   fn_particular_soln=['./code_matrices/Boundary_solution_z_H_junction-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
   ,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'-r_crac_index_',num2str(r_CRAC_index),'-theta_crac_index_',num2str(theta_CRAC_index),'-sr_',num2str(sr),'-st_',num2str(st),'.mat'];

    load(fn_particular_soln,'Boundary_solution_junction_z_H_rtz')
    particular_solution_z_H_all=particular_solution_z_H_all+reshape(Boundary_solution_junction_z_H_rtz,length(r_m)*length(theta),length(z_m));
    
    end
    flux_CRAC_M=flux_CRAC_magnitude*particular_solution_z_H_all';

    J_CRAC_flux_int=volume_integral_junction(vol_int_J_z_NT250,vol_int_J_rt_NT250,0,flux_CRAC_M);

    
    disp('loaded all CRAC fluxes')
toc

