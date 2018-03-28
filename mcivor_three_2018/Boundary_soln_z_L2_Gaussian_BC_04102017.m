% function particular_soln_M_z_L2_05052017(r_a)

%%% We are looking at the particular solution with boundary condition 
%%% Dm dv/dz (@z=L2) = J_{SERCA}. The shape is always the same as it just
%%% depends on the boundary conditions and geometry of the domain so we
%%% save the matrix solving dv/dz(@z=L2)=1 so we can multiply by J_{SERCA}/Dm later 
function Boundary_soln_z_L2_Gaussian_BC_04102017(r_bar)

% theta_bar_all=[2 16:15:150 2 16:15:150 2 16:15:150 2 16:15:150];
% r_bar_all=[120*ones(1,10) 150*ones(1,10) 180*ones(1,10) 320*ones(1,10)];
t=0:5:45;

dr_m=0.2;
sr=0.05;

if r_bar==80
st=0.03;
NT=250;
theta_bar_all=1+5*t;
else
st=0.05;
NT=150;
theta_bar_all=1+3*t;
end

r_bar_all=(r_bar/dr_m)*ones(1,length(theta_bar_all));
    

for i=1:length(theta_bar_all)
    theta_serca_index=theta_bar_all(i);
    r_serca_index=r_bar_all(i);

    tic
    r_a=100;
    eps=1e-16;

%     if r_serca_index==120
%         dr_m=0.5;
%         sr=0.1;
%     elseif r_serca_index==150
%         dr_m=0.5;
%         sr=0.1;
%     elseif r_serca_index==180
%         dr_m=0.5;
%         sr=0.1;
%     elseif r_serca_index==320
%         dr_m=0.25;
%         sr=0.25;
%     end
    
    r_b=2000;       %outer radius of annulus (nm)
    H=2500;         %height of annulus (nm)
    L2=2485;         %Height of ER (nm)
    L1=2000;
    
    dz_m=0.15;  %gives 101 meshpoints in z array
    dtheta=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array
       NZ=(H-L2)/dz_m;
 
    theta=0:dtheta:2*pi-eps;
    theta_prime=theta;
    r_m=dr_m:dr_m:r_a;
    r_m_prime=r_m;
    [r_m(r_serca_index) theta(theta_serca_index)]

    z_m=L2:dz_m:H;
    z_m_prime=z_m;
     
    NR=r_a/dr_m;
    % How many would I like to use for Green's function?
    mu_max_m=NZ;         % maximal number of eigenvalues in z - direction in microdomain
    beta_max_m=NR;          % maximal number of eigenvalues in r - direction  (j values)
    n_max_m=NT/2;            % maximal number of Bessel functions (n values)

    fn_zeros='./code_matrices/ER_PM_microdomain_radial_zeros-DBC_r_a-J_5000-N_2500.mat';
    load(fn_zeros)
    
    trap_r=dr_m*r_m_prime.*[0.5 ones(1,length(r_m_prime)-2) 0.5];
    trap_t=dtheta*[0.5 ones(1,length(theta_prime)-2) 0.5];
        
    tt=(theta-theta(theta_serca_index))'/st;
    rr=(r_m-r_m(r_serca_index))'/sr;
    
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

    ef_alpha_z_n_0_num=r_a*(exp(-alpha(1,:)'*(2*H-z_m-L2)/r_a)+exp(-alpha(1,:)'*(z_m-L2)/r_a));
    ef_alpha_z_n_0_den=alpha(1,:)'.*(exp(-2*alpha(1,:)'*(H-L2)/r_a)-1);
    ef_alpha_z_n_0=ef_alpha_z_n_0_num./(ef_alpha_z_n_0_den*ones(1,length(z_m)));

    ef_alpha_r_z_n_0 = ef_alpha_r_n_0*ef_alpha_z_n_0;
    ef_r_z_rp=reshape(ef_alpha_r_z_n_0',length(z_m)*length(r_m),length(r_m_prime_reduced));
    
    ef_alpha_r_theta_z_n_0=kron(ef_theta_n_0,ef_r_z_rp);

    ef_alpha_r_theta_z_n_non_0=zeros(length(r_m)*length(theta)*length(z_m),length(theta_prime_reduced)*length(r_m_prime_reduced));
    for k=1:n_max_m

        ef_alpha_r_n_non_0=besselj(k,r_m_all*alpha(k+1,:)/r_a).*besselj(k,r_m_prime_all*alpha(k+1,:)/r_a)*diag(1./besselj(k+1,alpha(k+1,:)).^2);
        ef_alpha_z_n_non_0_num=r_a*(exp(-alpha(k+1,:)'*(2*H-z_m-L2)/r_a)+exp(-alpha(k+1,:)'*(z_m-L2)/r_a));
        ef_alpha_z_n_non_0_den=alpha(k+1,:)'.*(exp(-2*alpha(k+1,:)'*(H-L2)/r_a)-1);
        ef_alpha_z_n_non_0=ef_alpha_z_n_non_0_num./(ef_alpha_z_n_non_0_den*ones(1,length(z_m)));

        ef_theta_n=cos(k*theta_theta_prime_all);
        ef=ef_alpha_r_n_non_0*ef_alpha_z_n_non_0;
        ef_r_z_rp=reshape(ef',length(z_m)*length(r_m),length(r_m_prime_reduced));

        ef_alpha_r_theta_z_n_non_0=ef_alpha_r_theta_z_n_non_0+2*kron(ef_theta_n,ef_r_z_rp);
     
    end

    Boundary_solution_junction_z_L2=(ef_alpha_r_theta_z_n_0+ef_alpha_r_theta_z_n_non_0)/(pi*r_a^2);
   
    clear ef_alpha_r_theta_z_n_0 ef_alpha_r_theta_z_n_non_0
    B_zrt=Boundary_solution_junction_z_L2*diag(trap_rt_reduced)*gauss_rt_reduced;
    B_z_r_t=reshape(B_zrt,length(z_m),length(r_m),length(theta));
    Boundary_solution_junction_z_L2_rtz=permute(B_z_r_t,[2 3 1]);
    
    clear Boundary_solution_junction_z_L2 B_zrt B_z_r_t 

    
    fn_particular_soln=['./code_matrices/Boundary_solution_z_L2_junction-mu_max_',num2str(mu_max_m),'-beta_max_',num2str(beta_max_m),'-n_max_',num2str(n_max_m)...
    ,'-dr',num2str(dr_m),'-dtheta_',num2str(dtheta),'-dz_',num2str(dz_m),'-r_a_',num2str(r_a),'-r_serca_index_',num2str(r_serca_index),'-theta_serca_index_',num2str(theta_serca_index),'-sr_',num2str(sr),'_st_',num2str(st),'.mat'];
    
    save(fn_particular_soln,'Boundary_solution_junction_z_L2_rtz','-v7.3')
%     display('Particular solution computed for z=L2')
toc
end
