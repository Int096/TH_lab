global mn_P_first N_med p_00 err_00 freq S nu disp num_ch2 experiment_time
% load('new_method_cor4.mat')
load('QD9_autocor');
t=g2_times(2:end)';
p=get_cor(2:end)'-1;
CC=CoinCounts(2:end);

p_00=get_cor(1)-1;
err_00=p_00/sqrt(CoinCounts(1));

raw_data=load('Traject.mat'); 
data=raw_data.Intensity; 

[Y_exp,X_exp]=hist(data,min(data):max(data));
Y_exp=Y_exp/sum(Y_exp);
N_med=Y_exp*X_exp';

diff_t=diff([0 t']);

new_bins=[logspace(-(5-1/2),-4,2) logspace(-(4-1/2),-3,2) logspace(-(3-1/3),-2,3) logspace(-(2-1/3),-1,3) logspace(-(1-1/5),0,5) logspace((0+1/10),1,10) logspace((1+1/15),2,10) logspace((2+1/10),log10(max(t)),10)];
% new_bins=[logspace(-(5-1/2),-4,2) logspace(-(4-1/2),-3,2) logspace(-(3-1/3),-2,3) logspace(-(2-1/3),-1,3) logspace(-(1-1/5),0,5) logspace((0+1/10),1,10) logspace((1+1/15),log10(max(t)),10)];
num_new_bins=numel(new_bins);

new_p=new_bins*0;
new_t=new_p;
new_t2=new_t;
new_p2=new_p;

kk=find((t<=new_bins(1)&t>0));

new_p(1)=mean(p(kk));
new_N(1)=mean(CC(kk));
new_t(1)=mean(t(kk));
new_p2(1)=sum(p(kk)'.*diff([0 t(kk)']))/sum(diff([0 t(kk)']));
new_t2(1)=max(t(kk));

for ii=2:num_new_bins
    kk=find(t<=new_bins(ii)&t>new_bins(ii-1));
    new_p(ii)=mean(p(kk));
    new_N(ii)=mean(CC(kk));
    new_p2(ii)=sum(p(kk)'.*diff([new_bins(ii-1) t(kk)']))/sum(diff([new_bins(ii-1) t(kk)']));
    new_t(ii)=mean(t(kk));
    new_t2(ii)=max(t(kk));
end


new_delta_t=diff([0 new_t]);
n_m=N_med/0.01;
n_m=num_ch2/experiment_time;
new_N_kor=(new_p+1).*new_delta_t*n_m^2*experiment_time;
% new_N_kor=new_N;
% sigma_p=sqrt(new_N_kor)./new_delta_t/n_m^2/max(t)/2;
% sigma_p=new_p./sqrt(new_N);

delta_t=diff([0 t']);
N_p=(p'+1).*delta_t*n_m^2*max(t)*2;
sig_np=sqrt(N_p);
sig_p=sqrt(N_p)./delta_t/n_m^2/max(t)/2;

ff=figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
set(gca,'Xscale','log')
% semilogx(t,p,'LineWidth',1);
hold on

PP=new_p;
tl=new_t;
N_pp=(PP+1).*diff([0 tl])*n_m^2*experiment_time;
sigma_PP=sqrt(N_pp)./tl/n_m^2/experiment_time;

plot(t,p)
hold on
errorbar(tl,PP,(new_N_kor-1/2*chi2inv(0.025,2*new_N_kor))./new_delta_t/n_m^2/experiment_time,(1/2*chi2inv(0.975,2*new_N_kor+2)-new_N_kor)./new_delta_t/n_m^2/experiment_time,'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1.5,'Color','k');
hold on
grid on

% errorbar(new_t2,new_p2,(new_N_kor-1/2*chi2inv(0.025,2*new_N_kor))./new_delta_t/n_m^2/max(t)/2,(1/2*chi2inv(0.975,2*new_N_kor+2)-new_N_kor)./new_delta_t/n_m^2/max(t)/2,'s','MarkerSize',6,'MarkerEdgeColor','r','LineWidth',1.5,'Color','r');

xlabel('Time [s]','Fontsize',24);
ylabel('Autocorrelation  function','Fontsize',24);
set(gca,'YLim',[0.99*min(p+1)-1 max(p)],'FontSize',24,'XLim',[0.8*min(tl) 1.02*max(t)]);
box on
% 
text(0.85,0.85,'QD 2 135','FontSize',24,'Color','black','Units','Normalized');
% saveas(ff,'Autocor_qd2_135','jpeg')

[freq,S,nu,disp]=PSD_Black(g2_times(2:end)',get_cor(2:end)'-1);

gg=figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
PSDplot()
hold on
% plot(freq,4e-2*1./freq,'m','LineWidth',2)
% plot(freq,1e-1*1./freq.^1.3,'r','LineWidth',2)
xlabel('Frequency, [Hz]')
ylabel('PSD, [Hz^{-1}]')
set(gca,'FontSize',24)
% legend('Estimated PSD','1/f','1/f^{1.3}')
% text(0.85,0.85,'QD 2 135','FontSize',24,'Color','black','Units','Normalized');
% saveas(gg,'PSD_qd2_135_3','jpeg')

ggg=figure('Color','white','Units','normalized','position',[0.05 0.05 0.85 0.85]);

subplot(3,4,[1 2 3])
plot(0.01*(1:numel(data)),data,'k','LineWidth',.5)
set(gca,'xlim',[0 0.01*numel(data)],'FontSize',26,'ylim',[min(X_exp) max(X_exp)],'ticklabelinterpreter','latex','XColor','k','Ycolor','k')
xl0=xlabel('Time [s]','units','normalized','interpreter','latex','Color','k');
ylabel('cnts/10 ms','interpreter','latex','Color','k')
annotation('textbox',[.05 .79 .1 .2], ...
    'String','(a)','EdgeColor','none','FontSize',24,'interpreter','latex')
xticks([0 500 1000 2000 2500 fix(0.01*numel(data))])


subplot(3,4,4)
br=bar(X_exp,Y_exp,1);
set(br,'Facecolor','k')
set(gca,'xdir','reverse','xlim',[min(X_exp) max(X_exp)],'XColor','k','Ycolor','k','FontSize',26,'ylim',[0 1.05*max(Y_exp)])
xticks([]);
yticks([])
view([90 90])

annotation('textbox',[.71 .79 .1 .2], ...
    'String','(b)','EdgeColor','none','FontSize',24,'interpreter','latex')

% x_pdf=xlabel('PDF','interpreter','latex','Color','k','units','normalized','rotation',-90);
% ps_dpf=get(x_pdf,'position');
% set(x_pdf,'position',[1.35 ps_dpf(2) ps_dpf(3)]);

subplot(3,4,[5 6 9 10])

annotation('textbox',[.05 .46 .1 .2], ...
    'String','(c)','EdgeColor','none','FontSize',24,'interpreter','latex')

errorbar(tl,PP,(new_N_kor-1/2*chi2inv(0.025,2*new_N_kor))./new_delta_t/n_m^2/experiment_time,(1/2*chi2inv(0.975,2*new_N_kor+2)-new_N_kor)./new_delta_t/n_m^2/experiment_time,'s','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1.5,'Color','k');
set(gca,'xscale','log','Xlim',[1e-5 5e3],'FontSize',26,'ticklabelinterpreter','latex','XColor','k','Ycolor','k')
xl=xlabel('Time [s]','units','normalized','interpreter','latex','Color','k');
ylabel('Autocorrelation Function','units','normalized','interpreter','latex','Color','k');

ps=get(xl,'position');
set(xl,'position',[ps(1) -.12 ps(3)])
grid on
xticks([1e-5 1e-3 1e-1 1e1 1e3])

subplot(3,4,[7 8 11 12])

annotation('textbox',[.51 .46 .1 .2], ...
    'String','(d)','EdgeColor','none','FontSize',24,'interpreter','latex')
yyaxis left
set(gca,'XColor','k','Ycolor','k')
yticks([]);
yyaxis right
PSDplot();
hold on
set(gca,'FontSize',26,'ticklabelinterpreter','latex','XColor','k','Ycolor','k')
xl=xlabel('Frequency [Hz]','units','normalized','interpreter','latex','Color','k');

ps=get(xl,'position');
set(xl,'position',[ps(1) -.12 ps(3)],'Color','k')

yl=ylabel('Power Spectral Density [Hz$^{-1}$]','interpreter','latex','Color','k','rotation',-90,'units','normalized');
ps=get(yl,'position');
set(yl,'position',[1.22 ps(2) ps(3)],'Color','k')

% plot(freq,4e-2*1./freq,'LineStyle','-','Color','m','LineWidth',2)
% 
% plot(freq,1e-1*1./freq.^1.3,'LineStyle','-','Color','r','LineWidth',2)
% text(0.43,0.15,'$PSD(f) \propto 1/f^{1.3}$','Units','Normalized','FontSize',24,'interpreter','latex')
% text(0.02,0.6,'$PSD(f) \propto 1/f$','Units','Normalized','FontSize',24,'interpreter','latex')


xticks([1e-2 1e-1 1e-0 1e1 1e2 1e3 1e4])
yticks([1e-7 1e-5 1e-3 1e-1 1e1 1e3])

grid on

ps=get(xl0,'position');
set(xl0,'position',[ps(1) -.07 ps(3)])

% saveas(ggg,'InitData','epsc')







