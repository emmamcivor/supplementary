function Boundary_soln_z_H_gaussian_21092017(r_bar)

t=5:10:45;
 dr_m_all= 0.2;

sr=0.05;

if r_bar==50
st=0.03;
NT=250;
theta_bar_all=1+5*t;
else
st=0.05;
NT=150;
theta_bar_all=1+3*t;
end

for i=1:length(t)
    tic
    theta_bar_index=theta_bar_all(i);
    dr_m=dr_m_all
    r_bar_index=r_bar/dr_m
    
    r_a=100;
    eps=1e-16;

        
    r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=2000;
    
%     dr_m=0.25;      %gives 75 mesh points in r array
    dz_m=0.15;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array

    NR=r_a/dr_m;    
    NZ=(H-L2)/dz_m;
    
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    r_m=dr_m:dr_m:r_a;
    r_m_prime=r_m;
        [r_m(r_bar_index) theta(theta_bar_index)]

%     dz_m_small=dz_m/3;
%     z_m_large_mesh=L2:dz_m:H-5*dz_m_small;
%     z_m_small_mesh=H-5*dz_m_small:dz_m_small:H;
%     z_m=[z_m_large_mesh z_m_small_mesh];
    z_m=L2:dz_m:H;
    z_m_prime=z_m;
    
        mu_max_m=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_m=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_m=NT/2;            % maximal number of Bessel functions (n values)

%     r_bar_index=80; 
%     theta_bar_index=30; 

%     if dr_m==0.5
%         sr=0.1;
%     elseif dr_m==0.25
%         sr=0.25;
%     elseif dr_m==0.2
%         sr=0.05;
%     end


    fn_particular_soln=['./code_matrices/Boundary_solution_z_H_junction-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
    ,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'-r_crac_index_',num2str(r_bar_index),'-theta_crac_index_',num2str(theta_bar_index),'-sr_',num2str(sr),'-st_',num2str(st),'.mat'];

if exist(fn_particular_soln)
    disp('exist')
else

    fn_zeros='./code_matrices/ER_PM_microdomain_radial_zeros-DBC_r_a-J_20000-N_2500.mat';
    load(fn_zeros)

    trap_r=dr_m*r_m_prime.*[0.5 ones(1,length(r_m_prime)-2) 0.5];
    trap_t=dtheta*[0.5 ones(1,length(theta_prime)-2) 0.5];
        
    tt=(theta-theta(theta_bar_index))'/st;
    rr=(r_m-r_m(r_bar_index))'/sr;
    
    tt_exp=exp(-tt.^2);
    rr_exp=exp(-rr.^2);
    
    tt_reduced=tt(tt_exp>1e-10);
    rr_reduced=rr(rr_exp>1e-10);
    
    trap_r_reduced=trap_r(rr_exp>1e-10);
    trap_t_reduced=trap_t(tt_exp>1e-10);
    trap_rt_reduced=kron(trap_t_reduced,trap_r_reduced);
    
    gauss_rt_reduced=exp(-(kron(tt_reduced.^2,ones(length(rr_reduced),1)) + kron(ones(length(tt_reduced),1),rr_reduced.^2)));
    
    r_m_prime_reduced=r_m_prime(rr_exp>1e-10);
    theta_prime_reduced=theta_prime(tt_exp>1e-10);

    
    r_m_all=kron(ones(length(r_m_prime_reduced),1),r_m'); % r.rp x 1
    r_m_prime_all=kron(r_m_prime_reduced',ones(length(r_m),1)); %r.rp x 1
    theta_theta_prime_all=repmat(theta',1,length(theta_prime_reduced))-repmat(theta_prime_reduced,length(theta),1);%t x tp

    
    alpha=eta_DBC(1:n_max_m+1,1:beta_max_m);

    ef_alpha_r_n_0=besselj(0,r_m_all*alpha(1,:)/r_a).*besselj(0,r_m_prime_all*alpha(1,:)/r_a)*diag(1./besselj(1,alpha(1,:)).^2);
    ef_theta_n_0=cos(0*(theta_theta_prime_all));

    ef_alpha_z_n_0_num=r_a*(exp(-alpha(1,:)'*(H-z_m)/r_a)+exp(-alpha(1,:)'*(H+z_m-2*L2)/r_a));
    ef_alpha_z_n_0_den=alpha(1,:)'.*(1-exp(-2*alpha(1,:)'*(H-L2)/r_a));
    ef_alpha_z_n_0=ef_alpha_z_n_0_num./(ef_alpha_z_n_0_den*ones(1,length(z_m)));

    ef_alpha_r_z_n_0 = ef_alpha_r_n_0*ef_alpha_z_n_0;
    ef_r_z_rp=reshape(ef_alpha_r_z_n_0',length(z_m)*length(r_m),length(r_m_prime_reduced));
    
    ef_alpha_r_theta_z_n_0=kron(ef_theta_n_0,ef_r_z_rp);

    ef_alpha_r_theta_z_n_non_0=zeros(length(r_m)*length(theta)*length(z_m),length(theta_prime_reduced)*length(r_m_prime_reduced));

    for k=1:n_max_m

        ef_alpha_r_n_non_0=besselj(k,r_m_all*alpha(k+1,:)/r_a).*besselj(k,r_m_prime_all*alpha(k+1,:)/r_a)*diag(1./besselj(k+1,alpha(k+1,:)).^2);
        ef_alpha_z_n_non_0_num=r_a*(exp(-alpha(k+1,:)'*(H-z_m)/r_a)+exp(-alpha(k+1,:)'*(H+z_m-2*L2)/r_a));
        ef_alpha_z_n_non_0_den=alpha(k+1,:)'.*(1-exp(-2*alpha(k+1,:)'*(H-L2)/r_a));
        ef_alpha_z_n_non_0=ef_alpha_z_n_non_0_num./(ef_alpha_z_n_non_0_den*ones(1,length(z_m)));

        ef_theta_n=cos(k*theta_theta_prime_all);
        ef=ef_alpha_r_n_non_0*ef_alpha_z_n_non_0;
        ef_r_z_rp=reshape(ef',length(z_m)*length(r_m),length(r_m_prime_reduced));

        ef_alpha_r_theta_z_n_non_0=ef_alpha_r_theta_z_n_non_0+2*kron(ef_theta_n,ef_r_z_rp);

    end
    Boundary_solution_junction_z_H=(ef_alpha_r_theta_z_n_0+ef_alpha_r_theta_z_n_non_0)/(pi*r_a^2);

    clear ef_alpha_r_theta_z_n_0 ef_alpha_r_theta_z_n_non_0
    B_zrt=Boundary_solution_junction_z_H*diag(trap_rt_reduced)*gauss_rt_reduced;
    B_z_r_t=reshape(B_zrt,length(z_m),length(r_m),length(theta));
    Boundary_solution_junction_z_H_rtz=permute(B_z_r_t,[2 3 1]);
    
    clear Boundary_solution_junction_z_H B_zrt B_z_r_t 
    
        
    save(fn_particular_soln,'Boundary_solution_junction_z_H_rtz','-v7.3')
end
        toc
%     figure;
%     plot(r_m,Boundary_solution_junction_z_H_rtz(:,30,end));savefig('r_plot')
%     
%     figure;
%     plot(theta,Boundary_solution_junction_z_H_rtz(80,:,end));savefig('theta_plot')
%     
%     d=(Boundary_solution_junction_z_H_rtz(:,:,end)-Boundary_solution_junction_z_H_rtz(:,:,end-1))/dz_m_small;
%     
%     figure;
%     plot(r_m,d(:,30));savefig('r_gradient_plot')
%     
%     figure;
%     plot(theta,d(80,:));savefig('theta_gradient_plot')
end
end