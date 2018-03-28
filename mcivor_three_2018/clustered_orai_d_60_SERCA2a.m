% lengths and r, phi, z grids
r_a=100;
val_Dm=220;
val_De=10;
Dm=val_Dm*1e6;     %Diffusion coefficent of calcium in the cytoplasm (nm^2/s) <--- CORRECT DIFF COEFF!!!
De=val_De*1e6;     %Diffusion coefficent of calcium in the ER (nm^2/s) <--- CORRECT DIFF COEFF!!!
n_timesteps=20;

eps=1e-16;
r_b=2000;       %outer radius of annulus (nm)
H=2500;         %height of annulus (nm)
L2=2485;         %Height of ER (nm)
L1=2000;
L0=10;

dz_m=0.15;  %gives 101 meshpoints in z array
z_m=L2:dz_m:H;
dz_s=4.85;  %gives 101 meshpoints in z array
z_s=L1:dz_s:L2;

NT=150;
dtheta_150=(2*pi-eps)/NT;  %gives 101 meshpoints in theta array
dtheta_250=dtheta_150*30/50;  %gives 101 meshpoints in theta array

theta_150=0:dtheta_150:2*pi-eps;
theta_250=0:dtheta_250:2*pi-eps;

dt=1e-6;


dr=0.2;      %gives 75 mesh points in r array
r_02=dr:dr:r_a;

% build matrices for plotting concentration profiles
[THETA_150,RHO] = meshgrid(theta_150,r_02);
[X_150,Y_150] = pol2cart(THETA_150,RHO);
[Z_150_s,t_150_s]=meshgrid(theta_150,z_s);
[Z_150_j,t_150_j]=meshgrid(theta_150,z_m);

[THETA_250,RHO] = meshgrid(theta_250,r_02);
[X_250,Y_250] = pol2cart(THETA_250,RHO);
[Z_250_s,t_250_s]=meshgrid(theta_250,z_s);
[Z_250_j,t_250_j]=meshgrid(theta_250,z_m);

%% load calcium concentration

% r_orai=30, r_SERCA=90, SERCA2a
fn_data='./simulations/SERCA2a_ICRAC_2.1e-15-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_150-5CRACs-r_SERCA_450-theta_SERCA_1-10_SERCAs-Vmax_1.1956e-16_n_max_75_Dm220_De_10_dt_1e-06_dr_0.2_dphi_0.041888_dz_0.15r_a_100_T_1000-sr_0.05-st_0.05/SOCE-p_50.mat';   
load(fn_data)
JrtzT_o30_s90_2a=reshape(ER_PM_junction_soln(:,:,end),length(r_02),length(theta_150),length(z_m));
JzrtT_o30_s90_2a=reshape(ER_PM_junction_soln(:,:,end)',length(z_m),length(r_02),length(theta_150));
JztrT_o30_s90_2a=permute(JzrtT_o30_s90_2a,[1 3 2]);

SrtzT_o30_s90_2a=reshape(Sub_PM_ER_soln(:,:,end),length(r_02),length(theta_150),length(z_m));
SzrtT_o30_s90_2a=reshape(Sub_PM_ER_soln(:,:,end)',length(z_m),length(r_02),length(theta_150));
SztrT_o30_s90_2a=permute(SzrtT_o30_s90_2a,[1 3 2]);

fn_data='./simulations/SERCA2a_ICRAC_2.1e-15-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_150-5CRACs-r_SERCA_450-theta_SERCA_1-10_SERCAs-Vmax_1.1956e-16_n_max_75_Dm220_De_10_dt_1e-06_dr_0.2_dphi_0.041888_dz_0.15r_a_100_T_1000-sr_0.05-st_0.05/SERCA_activity.mat';   
load(fn_data)
SERCA2a_activity_o30_s90=SERCA_activity;

%% make directory to save figures
cd_data_save='./simulation_figures';

if exist(cd_data_save,'dir')
else
    mkdir(cd_data_save);
end
cd(cd_data_save)

%% Make figures 



q7='../code_matrices/colormap_refill_rO90ab.mat';
load(q7)

figure(108);
pcolor(Z_150_s(60:end,:),t_150_s(60:end,:),SztrT_o30_s90_2a(60:end,:,90/dr))
shading interp
colormap(cmap_refill_rS90ab)
caxis([150 158.5])
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('\theta')
yticks([2290 2485])
yticklabels({'ER_{i}','ERM'})
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .89 1])
savefig('subPMER_rO30_rS90_radius_slice_r90_SERCA2a')
print('-f108','-bestfit','subPMER_rO30_rS90_radius_slice_r90_SERCA2a','-dpdf','-opengl')

