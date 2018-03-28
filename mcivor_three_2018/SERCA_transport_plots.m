

cd_data_save='./simulation_figures';

if exist(cd_data_save,'dir')
else
    mkdir(cd_data_save);
end
cd(cd_data_save)

kb=.27; nb=1.7; vb=36;
ka=.38; na=2.2; va=72;

x=0.01:0.01:3;

sanum=(x/ka).^na - (150/1700)^na;
saden=1 + (x/ka).^na + (150/1700)^na;
sa=va*sanum./saden;

sbnum=(x/kb).^nb - (150/1700)^nb;
sbden=1 + (x/kb).^nb + (150/1700)^nb;
sb=vb*sbnum./sbden;

x3060b=2.658;
sb3060num=(x3060b/kb).^nb - (150/1700)^nb;
sb3060den=1 + (x3060b/kb).^nb + (150/1700)^nb;
sb3060=vb*sb3060num./sb3060den;
sb3060act=sb3060/vb;

x3090b=0.49;
sb3090num=(x3090b/kb).^nb - (150/1700)^nb;
sb3090den=1 + (x3090b/kb).^nb + (150/1700)^nb;
sb3090=vb*sb3090num./sb3090den;
sb3090act=sb3090/vb;

x5080b=1.032;
sb5080num=(x5080b/kb).^nb - (150/1700)^nb;
sb5080den=1 + (x5080b/kb).^nb + (150/1700)^nb;
sb5080=vb*sb5080num./sb5080den;
sb5080act=sb5080/vb;

x3060a=2.18;
sa3060num=(x3060a/ka).^na - (150/1700)^na;
sa3060den=1 + (x3060a/ka).^na + (150/1700)^na;
sa3060=va*sa3060num./sa3060den;
sa3060act=sa3060/va;

x3090a=0.39;
sa3090num=(x3090a/ka).^na - (150/1700)^na;
sa3090den=1 + (x3090a/ka)^na + (150/1700)^na;
sa3090=va*sa3090num./sa3090den;
sa3090act=sa3090/va;
close all

figure(202); hold on;
plot(x,sb,'b','linewidth',3);
plot(x3060b,sb3060,'ro','markersize',10,'linewidth',3)
plot(x5080b,sb5080,'r+','markersize',10,'linewidth',3)
xlabel('C (\mu M)'); 
ylabel('Ca^{2+} ions/s');
ylim([0 40]);
xlim([x(1) x(end)])
set(gca,'fontsize',30,'fontweight','bold');
set(gca,'OuterPosition',[0 0 .9 0.9])
savefig('SERCA_activity_calcium_per_sec_SERCA2b_clustered_vs_non_clustered.fig')
print('-f202','SERCA_activity_calcium_per_sec_SERCA2b_clustered_vs_non_clustered','-bestfit','-dpdf','-opengl')

figure(203); hold on;
plot(x,sa,'k',x,sb,'b','linewidth',3);
plot(x3060a,sa3060,'x',x3090a,sa3090,'o','markersize',15,'markeredgecolor',[0,.9,0],'linewidth',3)
plot(x3060b,sb3060,'rx',x3090b,sb3090,'ro','markersize',15,'linewidth',3)
xlabel('C (\mu M)'); 
ylabel('Ca^{2+} ions/s');
ylim([0 80]);
xlim([x(1) x(end)])
set(gca,'fontsize',30,'fontweight','bold');
set(gca,'OuterPosition',[0 0 .9 0.9])
savefig('SERCA_activity_calcium_per_sec_SERCA2a_SERCA2b.fig')
print('-f203','SERCA_activity_calcium_per_sec_SERCA2a_SERCA2b','-bestfit','-dpdf','-opengl')
