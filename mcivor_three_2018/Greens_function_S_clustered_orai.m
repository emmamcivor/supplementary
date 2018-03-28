%Green's function for ER-PM microdomain Nov 2016
    r_a=100;
    
    val_De=10;
    dt=1e-6
    
    eps=1e-16;
    De=val_De*1e6;     %Diffusion coefficent of calcium in the cytoplasm (nm^2/s) <--- CORRECT DIFF COEFF!!!

    r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=2000;
    
    NT=150;
    
    dr_s=0.2;      %gives 75 mesh points in r array
    dz_s=4.85;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;
    
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    r_s=dr_s:dr_s:r_a;
    r_s_prime=r_s;
    
    NZ=(L2-L1)/dz_s;
    NR=r_a/dr_s;
    
%     dz_s_small=dz_s/1;
%     z_s=[L1:dz_s:L2-5*dz_s_small,L2-4*dz_s_small:dz_s_small:L2];
    z_s=L1:dz_s:L2;
    z_s_prime=z_s;

        % How many would I like to use for Green's function?
    mu_max_s=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_s=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_s=ceil(NT/2);            % maximal number of Bessel functions (n values)
    n=0:n_max_s;  

    fn_zeros='./code_matrices/Sub_PM_ER_radial_zeros-NBC_r_a-beta_max_10000-N_bessel_10000.mat';
    
    
fn_greens=['./code_matrices/Greens_function_sub-L1_',num2str(L1),'-dt_',num2str(dt),'-De_',num2str(val_De),'-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s)...
,'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'.mat'];
        
        %% Make figures of GF
        
if exist(fn_greens)
    load(fn_greens)
    disp('loaded GF')

            %Make directory to save figures
        cd_data=['./code_matrices/Figs/Greens_function_sub-L1_',num2str(L1),'-dt_',num2str(dt),'-De_',num2str(val_De),'-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s),'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'/']; 
        if exist(cd_data,'dir')
            cd(cd_data)
        else
            mkdir(cd_data);
            cd(cd_data)
        end    

    
        figure(1);
        plot(z_s,GF_sub_z_zp(:,50),'b','linewidth',3);
        xlabel('z','fontsize',16,'fontweight','bold');
        set(gca,'fontsize',16,'fontweight','bold')
        print('-f1','GF_z_z_p_50','-djpeg')

        GF_r_t_rp_tp=reshape(GF_sub_rt_rptp,length(r_s),length(theta),length(r_s_prime),length(theta_prime));
        [THETA,RHO] = meshgrid(theta,r_s);
        [Xsub,Ysub] = pol2cart(THETA,RHO);
        figure(2);
        pcolor(Xsub,Ysub,GF_r_t_rp_tp(:,:,60,30)); shading interp; colormap(jet); colorbar;
        xlabel('z','fontsize',16,'fontweight','bold');
        set(gca,'fontsize',16,'fontweight','bold')
        print('-f2','GF_r_t_rp_30_tp_30','-djpeg')
        
       GF_r_z=GF_sub_z_zp(:,end)*GF_r_t_rp_tp(:,30,60,30)';
        [zz,rr] = meshgrid(z_s,r_s);
        figure(3);
        pcolor(rr,zz,GF_r_z'); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f3','GF_r_z','-djpeg')
        
        GF_t_z=GF_sub_z_zp(:,end)*GF_r_t_rp_tp(60,:,60,30);
        [zzt,ttt] = meshgrid(z_s,theta);
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
        
zero_mu=1:mu_max_s;
mu=(2*zero_mu-1)*pi*0.5/(L2-L1);
ef_z_m=cos((L2-z_s')*mu);
ef_z_prime_m=cos(mu'*(L2-z_s_prime));

ef_z_z_prime_m=ef_z_m*diag(exp(-mu.^2*De*dt))*ef_z_prime_m;

GF_sub_z_zp=2*ef_z_z_prime_m/(L2-L1);

         
%% Polar GF comp
load(fn_zeros)
        beta=zero_beta(1:n_max_s+1,1:beta_max_s);
    
        % Structure of ef_beta_m_r_r_prime_dt:
        % i - index of r - coordinate
        % j - index of r_prime coordinate;

        ef_theta_theta_prime_r_r_prime_n_non0=zeros(length(r_s)*length(theta),length(r_s_prime)*length(theta_prime));

        theta_theta_prime=repmat(theta',1,length(theta_prime))-repmat(theta_prime,length(theta),1);

        for k=n    %Does k=0, then k=1, ...

            % Structure of ef_m_theta_theta_prime(i,j):
            % i - index of theta - coordinate
            % j - index of theta_prime coordinate

            ef_m_theta_theta_prime=cos(k*theta_theta_prime);

            norm_r=(1-(k./beta(k+1,:)).^2).*besselj(k,beta(k+1,:)).^2;
            norm_r=1./norm_r;
            if k==0
            norm_r(1)=1;
            end
            
            ef_beta_m_r=besselj(k,r_s'*beta(k+1,:)/r_a)*diag(norm_r.*exp(-(beta(k+1,:)/r_a).^2*De*dt));
            
            
            ef_beta_m_r_prime=besselj(k,beta(k+1,:)'*r_s_prime/r_a); %*diag(r_prime); only include this if dtau=0

            ef_beta_m_r_r_prime_dt=ef_beta_m_r*ef_beta_m_r_prime;   %sum over j (radial e-values?)
            
            if k==0
               ef_theta_theta_prime_r_r_prime_n0= kron(ef_m_theta_theta_prime,ef_beta_m_r_r_prime_dt);
            else
            ef_theta_theta_prime_r_r_prime_n_non0=ef_theta_theta_prime_r_r_prime_n_non0+2*kron(ef_m_theta_theta_prime,ef_beta_m_r_r_prime_dt); % count k=0 once, but k>0 twice (how k's are incorporated)
            end
            k
        end
        
        % Strucure of ef_r_theta_r_prime_theta_prime_IC_inc_r_prime
        % L(r)*L(theta) x L(r_prime)*L(theta_prime)
        ef_theta_theta_prime_r_r_prime_n0=(1/(pi*r_a^2))*ef_theta_theta_prime_r_r_prime_n0;
        ef_theta_theta_prime_r_r_prime_n_non0=(1/(pi*r_a^2))*ef_theta_theta_prime_r_r_prime_n_non0;

        GF_sub_rt_rptp=ef_theta_theta_prime_r_r_prime_n0+ef_theta_theta_prime_r_r_prime_n_non0;
        
%%%% save Greens functions        
        save(fn_greens,'GF_sub_z_zp','-v7.3');
        save(fn_greens,'GF_sub_rt_rptp','-append');
        disp('Greens Function computed')

        %Make directory to save figures
        cd_data=['./code_matrices/Figs/Greens_function_sub-L1_',num2str(L1),'-dt_',num2str(dt),'-De_',num2str(val_De),'-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s),'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'/']; 
        if exist(cd_data,'dir')
            cd(cd_data)
        else
            mkdir(cd_data);
            cd(cd_data)
        end    
        
        figure(1);
        plot(z_s,GF_sub_z_zp(:,50),'b','linewidth',3);
        xlabel('z','fontsize',16,'fontweight','bold');
        set(gca,'fontsize',16,'fontweight','bold')
        print('-f1','GF_z_z_p_50','-djpeg')

        GF_r_t_rp_tp=reshape(GF_sub_rt_rptp,length(r_s),length(theta),length(r_s_prime),length(theta_prime));
        [THETA,RHO] = meshgrid(theta,r_s);
        [Xsub,Ysub] = pol2cart(THETA,RHO);
        figure(2);
        pcolor(Xsub,Ysub,GF_r_t_rp_tp(:,:,60,30)); shading interp; colormap(jet); colorbar;
        xlabel('z','fontsize',16,'fontweight','bold');
        set(gca,'fontsize',16,'fontweight','bold')
        print('-f2','GF_r_t_rp_30_tp_30','-djpeg')
        
        GF_r_z=GF_sub_z_zp(:,end)*GF_r_t_rp_tp(:,30,60,30)';
        [zz,rr] = meshgrid(z_s,r_s);
        figure(3);
        pcolor(rr,zz,GF_r_z'); shading interp; colormap(jet); colorbar;
        xlabel('r');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f3','GF_r_z','-djpeg')
        
        GF_t_z=GF_sub_z_zp(:,end)*GF_r_t_rp_tp(60,:,60,30);
        [zzt,ttt] = meshgrid(z_s,theta);
        figure(4);
        pcolor(ttt,zzt,GF_t_z'); shading interp; colormap(jet); colorbar;
        xlabel('\phi');
        ylabel('z');
        set(gca,'fontsize',30,'fontweight','bold')
        set(gca,'OuterPosition',[0 0 0.93 0.93])       
        print('-f4','GF_theta_z','-djpeg')

end
%         %% Perform checks
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/SubPMER/';
%         cd(cd_check)
%         display('Performing GF check (rho)')
%         rho_r_delta_function_check(dt,val_De,mu_max_s,beta_max_s,n_max_s,dr_s,dtheta,dz_s,r_a,L1)
%         display('GF check (rho) finished')
%         
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/SubPMER/';
%         cd(cd_check)
%         display('Performing GF check (theta)')
%         theta_delta_function_check(dt,val_De,mu_max_s,beta_max_s,n_max_s,dr_s,dtheta,dz_s,r_a,L1)
%         display('GF check (theta) finished')      
%         
%         cd_check='/maths/staff/bhzem/SOCE_model_July_2017/SubPMER/';
%         cd(cd_check)       
%         display('Performing GF check (h)')
%         h_z_delta_function_check(dt,val_De,mu_max_s,beta_max_s,n_max_s,dr_s,dtheta,dz_s,r_a,L1,L2)
%         display('GF check (h) finished')        
%         
%         display('All GF checks finished - run Latex to check figures')

