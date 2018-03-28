%Green's function for ER-PM microdomain Nov 2016

    r_a=100;
    val_Dm=220;
    dt=1e-6
    
    
    eps=1e-16;
    Dm=val_Dm*1e6;     %Diffusion coefficent of calcium in the cytoplasm (nm^2/s) <--- CORRECT DIFF COEFF!!!

    r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=1000;

    NT=250;
    
    dr_m=0.2;      %gives 75 mesh points in r array
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

        % How many would I like to use for Green's function?
    mu_max_m=NZ       % maximal number of eigenvalues in z - direction in microdomain
    beta_max_m=NR          % maximal number of eigenvalues in r - direction  (j values)
    n_max_m=NT/2            % maximal number of Bessel functions (n values)
    n=0:n_max_m;  

    fn_zeros='./code_matrices/ER_PM_microdomain_radial_zeros-DBC_r_a-J_20000-N_2500.mat';
    
   fn_greens=['./code_matrices/Greens_function_junction-dt_',num2str(dt),'-Dm_',num2str(val_Dm),'-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'.mat'];


    
if exist(fn_greens)
load(fn_greens)
disp('loaded GF')

        %% Make figures of GF
            %Make directory to save figures
    cd_data=['./code_matrices/Figs_Greens_function_junction-dt_',num2str(dt),'-Dm_',num2str(val_Dm),'-NZ_',num2str(NZ),'-NR_',num2str(NR),'-NT_',num2str(NT),'-r_a_',num2str(r_a),'/']; 
    if exist(cd_data,'dir')
        cd(cd_data)
    else
        mkdir(cd_data);
        cd(cd_data)
    end 
        
        figure(1);
        plot(z_m,GF_junction_z_zp(:,50),'b','linewidth',3);
        xlabel('z','fontsize',30,'fontweight','bold');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f1','GF_z_z_p_50','-djpeg')

        GF_r_t_rp_tp=reshape(GF_junction_rt_rptp,length(r_m),length(theta),length(r_m_prime),length(theta_prime));
        [THETA,RHO] = meshgrid(theta,r_m);
        [Xmicro,Ymicro] = pol2cart(THETA,RHO);
        figure(2);
        pcolor(Xmicro,Ymicro,GF_r_t_rp_tp(:,:,30,30)); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('r');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f2','GF_r_t_rp_30_tp_30','-djpeg')
%         
        GF_r_z=GF_junction_z_zp(:,end)*GF_r_t_rp_tp(:,30,60,30)';
        [zz,rr] = meshgrid(z_m,r_m);
        figure(3);
        pcolor(rr,zz,GF_r_z'); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f3','GF_r_z','-djpeg')
        
        GF_t_z=GF_junction_z_zp(:,end)*GF_r_t_rp_tp(60,:,60,30);
        [zzt,ttt] = meshgrid(z_m,theta);
        figure(4);
        pcolor(ttt,zzt,GF_t_z'); shading interp; colormap(jet); colorbar;
        xlabel('\phi');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f4','GF_theta_z','-djpeg')

else
  
%% ER GF z comp
% Load eigenvalues
        
zero_mu=1:mu_max_m;
mu=zero_mu*pi/(H-L2);
ef_z_m=cos((z_m'-L2)*mu);
ef_z_prime_m=cos(mu'*(z_m_prime-L2));
norm_z_m=exp(-mu.^2*Dm*dt);
ef_z_z_prime_m=ef_z_m*diag(norm_z_m)*ef_z_prime_m;

GF_junction_z_zp=(1+2*ef_z_z_prime_m)/(H-L2);


         
%% Polar GF comp
 load(fn_zeros)
alpha=eta_DBC(1:n_max_m+1,1:beta_max_m);

        % Structure of ef_beta_m_r_r_prime_dt:
        % i - index of r - coordinate
        % j - index of r_prime coordinate;

        ef_theta_theta_prime_r_r_prime_n_non0=zeros(length(r_m)*length(theta),length(r_m_prime)*length(theta_prime));

        theta_theta_prime=repmat(theta',1,length(theta_prime))-repmat(theta_prime,length(theta),1);

        for k=n     %Does k=0, then k=1, ...
            % Structure of ef_m_theta_theta_prime(i,j):
            % i - index of theta - coordinate
            % j - index of theta_prime coordinate

            ef_m_theta_theta_prime=cos(k*theta_theta_prime);
            
            norm_r=besselj(k+1,alpha(k+1,:)).^2;
            norm_r_exp=exp(-(alpha(k+1,:)/r_a).^2*Dm*dt)./norm_r;          
            
            ef_alpha_m_r=besselj(k,r_m'*alpha(k+1,:)/r_a)*diag(norm_r_exp);
            
            ef_alpha_m_r_prime=besselj(k,alpha(k+1,:)'*r_m_prime/r_a); %*diag(r_prime); only include this if dtau=0

            ef_alpha_m_r_r_prime_dt=ef_alpha_m_r*ef_alpha_m_r_prime;   %sum over j (radial e-values?)
            
            clear norm_r norm_r_exp ef_alpha_m_r ef_alpha_m_r_prime 
            
            if k==0
                
                ef_theta_theta_prime_r_r_prime_n0= kron(ef_m_theta_theta_prime,ef_alpha_m_r_r_prime_dt);
                clear ef_m_theta_theta_prime ef_alpha_m_r_r_prime_dt
            else
                ef_theta_theta_prime_r_r_prime_n_non0=ef_theta_theta_prime_r_r_prime_n_non0+2*kron(ef_m_theta_theta_prime,ef_alpha_m_r_r_prime_dt); % count k=0 once, but k>0 twice (how k's are incorporated)
                clear ef_m_theta_theta_prime ef_alpha_m_r_r_prime_dt
            end
            k
        end
        
        % Strucure of ef_r_theta_r_prime_theta_prime_IC_inc_r_prime
        % L(r)*L(theta) x L(r_prime)*L(theta_prime)
        ef_theta_theta_prime_r_r_prime_n0=(1/(pi*r_a^2))*ef_theta_theta_prime_r_r_prime_n0;
        ef_theta_theta_prime_r_r_prime_n_non0=(1/(pi*r_a^2))*ef_theta_theta_prime_r_r_prime_n_non0;

        GF_junction_rt_rptp=ef_theta_theta_prime_r_r_prime_n0+ef_theta_theta_prime_r_r_prime_n_non0;

        save(fn_greens,'GF_junction_z_zp','-v7.3');
        save(fn_greens,'GF_junction_rt_rptp','-append');
        
        disp('Greens Function computed')
        %% Make figures of GF
    %Make directory to save figures
    cd_data=['./code_matrices/Figs_Greens_function_junction-dt_',num2str(dt),'-Dm_',num2str(val_Dm),'-NZ_',num2str(NZ),'-NR_',num2str(NR),'-NT_',num2str(NT),'-r_a_',num2str(r_a),'/']; 
    if exist(cd_data,'dir')
        cd(cd_data)
    else
        mkdir(cd_data);
        cd(cd_data)
    end
    
        figure(1);
        plot(z_m,GF_junction_z_zp(:,50),'b','linewidth',3);
        xlabel('z','fontsize',30,'fontweight','bold');
        set(gca,'fontsize',30,'fontweight','bold')
        print('-f1','GF_z_z_p_50_t_1e-6_junction','-djpeg')

        GF_r_t_rp_tp=reshape(GF_junction_rt_rptp,length(r_m),length(theta),length(r_m_prime),length(theta_prime));
        [THETA,RHO] = meshgrid(theta,r_m);
        [Xmicro,Ymicro] = pol2cart(THETA,RHO);
        figure(2);
        pcolor(Xmicro,Ymicro,GF_r_t_rp_tp(:,:,30,30)); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('r');
        set(gca,'fontsize',30,'fontweight','bold')
        print('-f2','GF_r_t_rp_30_tp_30_t_1e-6_junction','-djpeg')
%         
        GF_r_z=GF_junction_z_zp(:,end)*GF_r_t_rp_tp(:,30,60,30)';
        [zz,rr] = meshgrid(z_m,r_m);
        figure(3);
        pcolor(rr,zz,GF_r_z'); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f3','GF_r_z_dt_1e-6_junction','-djpeg')
        
        GF_t_z=GF_junction_z_zp(:,end)*GF_r_t_rp_tp(60,:,60,30);
        [zzt,ttt] = meshgrid(z_m,theta);
        figure(4);
        pcolor(ttt,zzt,GF_t_z'); shading interp; colormap(jet); colorbar;
        xlabel('\phi');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f4','GF_theta_z','-djpeg')
end
%         %% Perform checks
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/Junction/';
%         cd(cd_check)
%         display('Performing GF check (rho)')
%         rho_r_delta_function_check(dt,val_Dm,mu_max_m,beta_max_m,n_max_m,dr_m,dtheta,dz_m,r_a,1)
%         display('GF check (rho) finished')
%         
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/Junction/';
%         cd(cd_check)
%         display('Performing GF check (theta)')
%         theta_delta_function_check(dt,val_Dm,mu_max_m,beta_max_m,n_max_m,dr_m,dtheta,dz_m,r_a,1)
%         display('GF check (theta) finished')      
%         
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/Junction/';
%         cd(cd_check)       
%         display('Performing GF check (h)')
%         h_z_delta_function_check(dt,val_Dm,mu_max_m,beta_max_m,n_max_m,dr_m,dtheta,dz_m,r_a,1)
%         display('GF check (h) finished')        
%         
%         display('All GF checks finished - run Latex to check figures')

