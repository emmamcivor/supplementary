function Boundary_soln_S_z_L2_Gaussian_04102017(r_bar)

%%% We are looking at the particular solution with boundary condition 
%%% De dv/dz (@z=L2) = J_{SERCA}. The shape is always the same as it just
%%% depends on the boundary conditions and geometry of the domain so we
%%% save the matrix solving dv/dz(@z=L2)=1 so we can multiply by J_{SERCA}/De later 
% theta_bar_all=[2 16:15:150 2 16:15:150 2 16:15:150 2 16:15:150];
% r_bar_all=[120*ones(1,10) 150*ones(1,10) 180*ones(1,10) 320*ones(1,10)];
t=0:5:45;

dr_s=0.2;
sr=0.05; 
if r_bar==80
st=0.03;  %%%%0.03 if NT=250. 0.0375 if NT=200. 0.05 if NT=150
NT=250;
theta_bar_all=1+5*t;
else
st=0.05;  %%%%0.03 if NT=250. 0.0375 if NT=200. 0.05 if NT=150
NT=150;
theta_bar_all=1+3*t;
end

r_bar_all=(r_bar/dr_s)*ones(1,length(theta_bar_all));

for i=1:length(theta_bar_all)
    theta_serca_index=theta_bar_all(i);
    r_serca_index=r_bar_all(i);

tic

    r_a=100;
    eps=1e-16;

%     if r_serca_index==120
%         dr_s=0.5;
%         sr=0.1;
%     elseif r_serca_index==150
%         dr_s=0.5;
%         sr=0.1;
%     elseif r_serca_index==180
%         dr_s=0.5;
%         sr=0.1;
%     elseif r_serca_index==320
%         dr_s=0.25;
%         sr=0.275;
%     end

r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=2000;
    
    dz_s=4.85;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array
    NZ=(L2-L1)/dz_s;
    
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    r_s=dr_s:dr_s:r_a;
    r_s_prime=r_s;
    
    z_s=L1:dz_s:L2;
    z_s_prime=z_s;
     
    [r_s(r_serca_index) theta(theta_serca_index)]
    NR=r_a/dr_s;

    % How many would I like to use for Green's function?
    mu_max_s=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_s=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_s=ceil(NT/2);            % maximal number of Bessel functions (n values)

    fn_zeros='./code_matrices/Sub_PM_ER_radial_zeros-NBC_r_a-beta_max_10000-N_bessel_10000.mat';
    load(fn_zeros)
        
    
    trap_r=dr_s*r_s_prime.*[0.5 ones(1,length(r_s_prime)-2) 0.5];
    trap_t=dtheta*[0.5 ones(1,length(theta_prime)-2) 0.5];
        
    tt=(theta-theta(theta_serca_index))'/st;
    rr=(r_s-r_s(r_serca_index))'/sr;
    
    tt_exp=exp(-tt.^2);
    rr_exp=exp(-rr.^2);
    
    tt_reduced=tt(tt_exp>1e-10);
    rr_reduced=rr(rr_exp>1e-10);
    
    trap_r_reduced=trap_r(rr_exp>1e-10);
    trap_t_reduced=trap_t(tt_exp>1e-10);
    trap_rt_reduced=kron(trap_t_reduced,trap_r_reduced);
    
    gauss_rt_reduced=exp(-(kron(tt_reduced.^2,ones(length(rr_reduced),1)) + kron(ones(length(tt_reduced),1),rr_reduced.^2)));
    
    r_s_prime_reduced=r_s_prime(rr_exp>1e-10);
    theta_prime_reduced=theta_prime(tt_exp>1e-10);

    
    r_s_all=kron(ones(length(r_s_prime_reduced),1),r_s'); % r.rp x 1
    r_s_prime_all=kron(r_s_prime_reduced',ones(length(r_s),1)); %r.rp x 1
    theta_theta_prime_all=repmat(theta',1,length(theta_prime_reduced))-repmat(theta_prime_reduced,length(theta),1);%t x tp
    

    beta=zero_beta(1:n_max_s+1,1:beta_max_s);
    
    ef_n_0_m_0=repmat(z_s'-L1,length(r_s)*length(theta),length(theta_prime_reduced)*length(r_s_prime_reduced));
    ef_theta_n_0=cos(0*(theta_theta_prime_all));

    ef_n_0_m_0_r_j=besselj(0,r_s_all*beta(1,2:end)/r_a).*besselj(0,r_s_prime_all*beta(1,2:end)/r_a)*diag(1./(besselj(0,beta(1,2:end)).^2));
    ef_n_0_m_0_j_z_num=exp(beta(1,2:end)'*(z_s-L2)/r_a)-exp(-beta(1,2:end)'*(z_s+L2-2*L1)/r_a);
    ef_n_0_m_0_j_z_den=(1+exp(-2*beta(1,2:end)'*(L2-L1)/r_a)).*beta(1,2:end)'/r_a;
    ef_n_0_m_0_j_z=ef_n_0_m_0_j_z_num'*diag(1./ef_n_0_m_0_j_z_den);

    ef_alpha_r_z_n_0_j2end = ef_n_0_m_0_r_j*ef_n_0_m_0_j_z';
    ef_r_z_rp=reshape(ef_alpha_r_z_n_0_j2end',length(z_s)*length(r_s),length(r_s_prime_reduced));
    
    ef_alpha_r_theta_z_n_0_j2end=kron(ef_theta_n_0,ef_r_z_rp);

    ef_n_non_0_m_non_0_r_theta_z=zeros(length(r_s)*length(theta)*length(z_s),length(theta_prime_reduced)*length(r_s_prime_reduced));

    for k=1:n_max_s
    ef_m_theta_theta_prime=cos(k*theta_theta_prime_all);
    norm=(besselj(k,beta(k+1,:)).^2).*(1-(k./beta(k+1,:)).^2);
    ef_n_non_0_m_non_0_r_j=besselj(k,r_s_all*beta(k+1,:)/r_a).*besselj(k,r_s_prime_all*beta(k+1,:)/r_a)*diag(1./norm);
    ef_n_non_0_m_non_0_j_z_num=exp(beta(k+1,:)'*(z_s-L2)/r_a)-exp(-beta(k+1,:)'*(z_s+L2-2*L1)/r_a);
    ef_n_non_0_m_non_0_j_z_den=(1+exp(-2*beta(k+1,:)'*(L2-L1)/r_a)).*beta(k+1,:)'/r_a;
    ef_n_non_0_m_non_0_j_z=ef_n_non_0_m_non_0_j_z_num'*diag(1./ef_n_non_0_m_non_0_j_z_den);

    ef=ef_n_non_0_m_non_0_r_j*ef_n_non_0_m_non_0_j_z';
    ef_r_z_rp=reshape(ef',length(z_s)*length(r_s),length(r_s_prime_reduced));
    
    
    ef_n_non_0_m_non_0_r_theta_z=ef_n_non_0_m_non_0_r_theta_z+2*kron(ef_m_theta_theta_prime,ef_r_z_rp);
        
    end

    Boundary_solution_subPMER_z_L2=(ef_n_0_m_0+ef_alpha_r_theta_z_n_0_j2end+ef_n_non_0_m_non_0_r_theta_z)/(pi*r_a^2);
    B_zrt=Boundary_solution_subPMER_z_L2*diag(trap_rt_reduced)*gauss_rt_reduced;
    B_z_r_t=reshape(B_zrt,length(z_s),length(r_s),length(theta));
    Boundary_solution_subPMER_z_L2_rtz=permute(B_z_r_t,[2 3 1]);
    
    clear Boundary_solution_subPMER_z_L2 B_zrt B_z_r_t 
    
    fn_particular_soln=['./code_matrices/Boundary_solution_z_L2_sub-mu_max_',num2str(mu_max_s),'-beta_max_',num2str(beta_max_s),'-n_max_',num2str(n_max_s)...
    ,'-dr',num2str(dr_s),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_s),'-r_a_',num2str(r_a),'-r_SERCA_index_',num2str(r_serca_index),'-theta_SERCA_index_',num2str(theta_serca_index),'-sr_',num2str(sr),'_st_',num2str(st),'.mat'];
    
    save(fn_particular_soln,'Boundary_solution_subPMER_z_L2_rtz','-v7.3')
%     display('Particular solution computed for z=L2')
    

    toc
end