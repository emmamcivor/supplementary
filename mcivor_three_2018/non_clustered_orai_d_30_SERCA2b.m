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

% r_orai=50, r_SERCA=80, SERCA2b
fn_data='./simulations/SERCA2b_ICRAC_2.1e-15-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_250-5CRACs-r_SERCA_400-theta_SERCA_1-10_SERCAs-Vmax_5.9781e-17_n_max_125_Dm220_De_10_dt_1e-06_dr_0.2_dphi_0.025133_dz_0.15r_a_100_T_1000-sr_0.05-st_0.03/SOCE-p_50.mat';   
load(fn_data)
JrtzT_o50_s80_2b=reshape(ER_PM_junction_soln(:,:,end),length(r_02),length(theta_250),length(z_m));
JzrtT_o50_s80_2b=reshape(ER_PM_junction_soln(:,:,end)',length(z_m),length(r_02),length(theta_250));
JztrT_o50_s80_2b=permute(JzrtT_o50_s80_2b,[1 3 2]);

SrtzT_o50_s80_2b=reshape(Sub_PM_ER_soln(:,:,end),length(r_02),length(theta_250),length(z_m));
SzrtT_o50_s80_2b=reshape(Sub_PM_ER_soln(:,:,end)',length(z_m),length(r_02),length(theta_250));
SztrT_o50_s80_2b=permute(SzrtT_o50_s80_2b,[1 3 2]);

fn_data='./simulations/SERCA2b_ICRAC_2.1e-15-flux_per_estimated_area_orai_channel-J_S_Gaussian_BC-r_CRAC_250-5CRACs-r_SERCA_400-theta_SERCA_1-10_SERCAs-Vmax_5.9781e-17_n_max_125_Dm220_De_10_dt_1e-06_dr_0.2_dphi_0.025133_dz_0.15r_a_100_T_1000-sr_0.05-st_0.03/SERCA_activity.mat';   
load(fn_data)
SERCA2b_activity_o50_s80=SERCA_activity;

%% make directory to save figures
cd_data_save='./simulation_figures';

if exist(cd_data_save,'dir')
else
    mkdir(cd_data_save);
end
cd(cd_data_save)

%% Make figures 
% cut along diameter of junction

r=[-fliplr(r_02),r_02];

JrtzT_rORAI_50_PM=[fliplr(JrtzT_o50_s80_2b(:,76+75,end)'),JrtzT_o50_s80_2b(:,76,end)'];
JrtzT_rORAI_50_ER=[fliplr(JrtzT_o50_s80_2b(:,76+75,1)'),JrtzT_o50_s80_2b(:,76,1)'];


figure(1);
plot(r,JrtzT_rORAI_50_PM,'r','linewidth',3)
xlim([-r_a r_a])
ylim([0 65.4])
xlabel('r (nm)')
ylabel('C (\mu M)')
set(gca,'fontsize',30,'fontweight','bold')
savefig('junction_diameter_J_PM_rO30_rS60_rO50_rS80')
print('-f1','-bestfit','junction_diameter_J_PM_rO50_rS80','-dpdf','-opengl')

figure(2);
plot(r,JrtzT_rORAI_50_ER,'r','linewidth',3)
xlim([-r_a r_a])
ylim([0 8])
xlabel('r (nm)')
ylabel('C (\mu M)')
set(gca,'fontsize',30,'fontweight','bold')
savefig('junction_diameter_J_ER_rO30_rS60_rO50_rS80')
print('-f2','-bestfit','junction_diameter_J_ER_rO50_rS80','-dpdf','-opengl')

q1='../code_matrices/colormap_zPM.mat';
load(q1)


figure(4);
pcolor(X_250,Y_250,JrtzT_o50_s80_2b(:,:,end))
shading interp
colormap(cmap_zPM)
caxis([cmin_H cmax_H])
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('r')
ylabel('r')
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])

[x60,y60]=pol2cart(theta_150(46),60);
[x80,y80]=pol2cart(theta_250(76),80);
[x90,y90]=pol2cart(theta_150(46),90);
[xd1,yd1]=pol2cart(theta_150(16),r_02);
[xd2,yd2]=pol2cart(theta_150(16+75),r_02);
[xd3,yd3]=pol2cart(theta_250(26),r_02);
[xd4,yd4]=pol2cart(theta_250(26+125),r_02);

figure(4);
hold on;
plot(xd3,yd3,'k--',xd4,yd4,'k--','linewidth',2)
savefig('junction_rO50_rS80_PM_SERCA2b')
print('-f4','-bestfit','junction_rO50_rS80_PM_SERCA2b','-dpdf','-opengl')


figure(6);
pcolor(X_250,Y_250,JrtzT_o50_s80_2b(:,:,1))
shading interp
colormap(jet)
caxis([cmin_L2 cmax_L2])
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('r')
ylabel('r')
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])

figure(6);
hold on;
plot(x80,y80,'xk','markersize',10,'linewidth',3)
plot(xd3,yd3,'k--',xd4,yd4,'k--','linewidth',3)
savefig('junction_rO50_rS80_ER_SERCA2b')
print('-f6','-bestfit','junction_rO50_rS80_ER_SERCA2b','-dpdf','-opengl')



q6='../code_matrices/colormap_refill_rO3050.mat';
load(q6)

figure(104);
pcolor(Z_250_s(60:end,:),t_250_s(60:end,:),SztrT_o50_s80_2b(60:end,:,80/dr))
shading interp
colormap(cmap_refill_rO3050)
caxis([d1 d2])
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('\theta')
yticks([2290 2485])
yticklabels({'ER_{i}','ERM'})
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .89 1])
savefig('subPMER_rO50_rS80_radius_slice_r80_SERCA2b')
print('-f104','-bestfit','subPMER_rO50_rS80_radius_slice_r80_SERCA2b','-dpdf','-opengl')

