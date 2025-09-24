global decay_time decay
raw_data2=load(char('foo2.csv'),'-ascii');
delays=raw_data2(:,4)*1e-12;
arrivals=raw_data2(:,1)*1e-12;
detector_flag=raw_data2(:,2);
first_dec_indexes=find(detector_flag==0);
second_dec_indexes=find(detector_flag==1);
delays(first_dec_indexes)=delays(first_dec_indexes)-2400*1*1e-12;
delays(second_dec_indexes)=delays(second_dec_indexes)-2900*1*1e-12;


arrivals(delays<0)=[];
delays(delays<0)=[];
arrivals(delays>1e-6)=[];
delays(delays>1e-6)=[];

arrivals=arrivals-arrivals(1);
nm_arrivals=numel(arrivals);
bin_size=0.01;

num_ph_bin=1000;
n_max=fix(numel(delays)/num_ph_bin);

delay_matrix=zeros(n_max,num_ph_bin);

for ii=1:n_max
    k=((ii-1)*num_ph_bin+1):ii*num_ph_bin;
    delay_matrix(ii,:)=delays(k)*1e9;
    right_time(ii)=arrivals(k(end));
end

lifetime=ones(1,n_max);
back=lifetime;
dec_time2=lifetime;
coef=lifetime;

for ii=1:n_max
    ii/n_max*100
    data=delay_matrix(ii,:);
    [decay_hist,decay_time]=hist(data,(.5):(ceil(max(data))-.5));
    %     decay_time=decay_time+.5;
    decay=decay_hist/sum(decay_hist);
    x0=[10 1e-3 1 100];
    param=fminsearch(@MaxLikelihoodPower,x0,optimset('MaxFunEvals',10000,'MaxIter',10000));
    lifetime(ii)=param(1);
    back(ii)=param(2);
    coef(ii)=param(3);
    dec_time2(ii)=param(4);
%     fn=1/lifetime(ii)*exp(-decay_time/lifetime(ii))+back(ii)+coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii));
%     perc(ii)=sum(1/lifetime(ii)*exp(-decay_time/lifetime(ii)))/sum(fn);
%     perc_l(ii)=sum(coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii)))/sum(fn);
%     fn=fn/sum(fn);
    if ii>1
        intensity_1000=1000/(right_time(ii)-right_time(ii-1))/100;
    else
        intensity_1000=1000/(right_time(ii))/100;
    end
    inten_short(ii)=intensity_1000*perc(ii);
    inten_long(ii)=intensity_1000*perc_l(ii);
    %     figure()
    %     plot(decay_time,fn)
    %     hold on
    %     plot(decay_time,decay)
end

% save('fit_res','lifetime','back','coef','dec_time2');
load('fit_res')
load('Traject.mat')

for ii=1:n_max
    data=delay_matrix(ii,:);
    [decay_hist,decay_time]=hist(data,(.5):(ceil(max(data))-.5));
    %     decay_time=decay_time+.5;
    decay=decay_hist/sum(decay_hist);
    fn=1/lifetime(ii)*exp(-decay_time/lifetime(ii))+back(ii)+coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii));
    perc(ii)=sum(1/lifetime(ii)*exp(-decay_time/lifetime(ii)))/sum(fn);
    perc_l(ii)=sum(coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii)))/sum(fn);
    fn=fn/sum(fn);
    if ii>1
        intensity_1000=1000/(right_time(ii)-right_time(ii-1))/100;
    else
        intensity_1000=1000/(right_time(ii))/100;
    end
    inten_short(ii)=intensity_1000*perc(ii);
    inten_long(ii)=intensity_1000*perc_l(ii);
end

fg=figure('Color','white','Units','normalized','position',[0.02 0.02 0.95 0.95]);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

annotation('textbox',[.03 .8 .1 .2], ...
    'String','(a)','EdgeColor','none','FontSize',34,'interpreter','latex')

subplot(3,2,1)
plot((0:(numel(Intensity)-1))*1e-2,Intensity,'k');
set(gca,'XLim',[0 (numel(Intensity)-1)*1e-2],'FontSize',34,'YLim',[0 max(Intensity)+5],'Xtick',[],'ticklabelinterpreter','latex','xcolor','k','ycolor','k')
y1=ylabel('cnts/10 ms','interpreter','latex','Color','k','units','normalized');
psy=get(y1,'Position');


annotation('textbox',[.03 .48 .1 .2], ...
    'String','(b)','EdgeColor','none','FontSize',34,'interpreter','latex')
subplot(3,2,3)
plot(right_time,lifetime,'k')
set(gca,'XLim',[0 (numel(Intensity)-1)*1e-2],'FontSize',34,'Xtick',[],'ticklabelinterpreter','latex','xcolor','k','ycolor','k')
t1=ylabel('$\tau_F$ [ns]','interpreter','latex','units','normalized');
ps1=get(t1,'Position');
yticks([0 25 50])

annotation('textbox',[.03 .18 .1 .2], ...
    'String','(c)','EdgeColor','none','FontSize',34,'interpreter','latex')
subplot(3,2,5)
plot(right_time,dec_time2,'k')
set(gca,'XLim',[0 (numel(Intensity)-1)*1e-2],'YLim',[0 500],'FontSize',34,'ticklabelinterpreter','latex','xcolor','k','ycolor','k')
t=ylabel('$\tau_D$ [ns]','interpreter','latex','Color','k','units','normalized');
xlab_1=xlabel('Time [s]','FontSize',34,'interpreter','latex','Color','k','units','normalized');

yticks([0 250 500])
xticks([0 1000 2000 fix(0.01*(numel(Intensity)-1))])
ps=get(t,'Position');
set(t1,'position',[psy(1) ps1(2) ps1(3)])
set(t,'position',[psy(1) ps(2) ps(3)])
annotation('textbox',[.47 .8 .1 .2], ...
    'String','(d)','EdgeColor','none','FontSize',34,'interpreter','latex')


subplot(3,2,[2 4 6])
hist3([ [0 lifetime]' [0 dec_time2]'],[175 175])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','EdgeColor','none');
view(2)

ylabel('Delayed time [ns]','interpreter','latex','Color','k');
set(gca,'FontSize',34,'Color',[0 0 0.49],'ticklabelinterpreter','latex','xcolor','k','ycolor','k','xlim',[0 30],'ylim',[0 300]);
colormap parula
axis([0 30 0 300])
cmap=colormap;
cmap(1,:)=[0 0 0.49];
colormap(cmap)
cb=colorbar('location','east','Color','w');
cb.Label.String='Occurence';
cb.Label.Interpreter='latex';
cb.Label.FontSize=34;
cb.TickLabelInterpreter='latex';

xlab=xlabel('Fast time [ns]','interpreter','latex','units','normalized');
ps_x=get(xlab,'position');
ps_x1=get(xlab_1,'position');
set(xlab_1,'position',[ps_x1(1) -0.27 ps_x1(3)])
set(xlab,'position',[ps_x(1) -0.065  ps_x(3)])

% saveas(gcf,'TrajLifetimes','jpeg')
% saveas(gcf,'TrajLifetimes','epsc')

figure('Color','white','Units','normalized','position',[0.02 0.02 0.95 0.95]);

hist3([ [0 lifetime]' [0 inten_short]'/max(inten_short+inten_long)],[175 175])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','EdgeColor','none');
view(2)

ylabel('Short yield','interpreter','latex','Color','k');
xlab=xlabel('Fast time [ns]','interpreter','latex','units','normalized');
set(gca,'FontSize',34,'Color',[0 0 0.49],'ticklabelinterpreter','latex','xcolor','k','ycolor','k','xlim',[0 25],'ylim',[0 .4]);
colormap parula
% axis([0 30 0 300])
cmap=colormap;
cmap(1,:)=[0 0 0.49];
colormap(cmap)
cb=colorbar('location','east','Color','w');
cb.Label.String='Occurence';
cb.Label.Interpreter='latex';
cb.Label.FontSize=34;
cb.TickLabelInterpreter='latex';

figure('Color','white','Units','normalized','position',[0.02 0.02 0.95 0.95]);

hist3([ [0 dec_time2]' [0 inten_long]'/max(inten_short+inten_long)],[175 175])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','EdgeColor','none');
view(2)

ylabel('Short yield','interpreter','latex','Color','k');
xlab=xlabel('Fast time [ns]','interpreter','latex','units','normalized');
set(gca,'FontSize',34,'Color',[0 0 0.49],'ticklabelinterpreter','latex','xcolor','k','ycolor','k','xlim',[0 300],'ylim',[0 .8]);
colormap parula
% axis([0 30 0 300])
cmap=colormap;
cmap(1,:)=[0 0 0.49];
colormap(cmap)
cb=colorbar('location','east','Color','w');
cb.Label.String='Occurence';
cb.Label.Interpreter='latex';
cb.Label.FontSize=34;
cb.TickLabelInterpreter='latex';

figure('Color','white','Units','normalized','position',[0.02 0.02 0.95 0.95]);

hist3([ [0 inten_long+inten_short]'/max(inten_short+inten_long) [0 inten_short]'/max(inten_short+inten_long)],[175 175])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto','EdgeColor','none');
view(2)

ylabel('Short yield','interpreter','latex','Color','k');
xlab=xlabel('Fast time [ns]','interpreter','latex','units','normalized');
set(gca,'FontSize',34,'Color',[0 0 0.49],'ticklabelinterpreter','latex','xcolor','k','ycolor','k','xlim',[0 1],'ylim',[0 .4]);
colormap parula
% axis([0 30 0 300])
cmap=colormap;
cmap(1,:)=[0 0 0.49];
colormap(cmap)
cb=colorbar('location','east','Color','w');
cb.Label.String='Occurence';
cb.Label.Interpreter='latex';
cb.Label.FontSize=34;
cb.TickLabelInterpreter='latex';

% saveas(gcf,'each_1000_FLID','jpeg')


% % fg0 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% % plot(right_time,dec_time2,'k','LineWidth',1.5)
% % hold on
% % xlabel('Time,[s]')
% % ylabel('Delayed lifetime, [ns]')
% % set(gca,'XLim',[ch1(1) ch1(end)],'FontSize',24);
% % % saveas(fg0,'Delay_vs_time','jpeg')
% % 
% % fg1= figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% % plot(right_time,lifetime,'k','LineWidth',1.5)
% % hold on
% % xlabel('Time,[s]')
% % ylabel('Lifetime,[ns]');
% % set(gca,'XLim',[ch1(1) ch1(end)],'FontSize',24);
% % % saveas(fg1,'Lifetime_vs_time','jpeg')
% % 
% % fg2 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% % % plot(dec_time2,lifetime,'o')
% % hist3([lifetime' dec_time2'],[100 100])
% % set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% % view(2)
% % xlabel('Short lifetime, [ns]');
% % ylabel('Delayed lifetime, [ns]');
% % set(gca,'FontSize',24);
% % colormap parula
% % axis([min(lifetime) max(lifetime) min(dec_time2) max(dec_time2)])
% % cmap=colormap;
% % cmap(1,:)=[0 0 0.49];
% % colormap(cmap)
% % cb=colorbar('location','east','Color','w');
% % % saveas(fg2,'Delay_vs_lifetime','jpeg')

% fg3 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% plot(right_time,inten_short)
% fg4 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% plot(right_time,inten_long)


% for jj=1:9
% fg(jj) = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% ii=jj*1000; %3000 2 1 все кратны 1000 4 5
% data=delay_matrix(ii,:);
% [decay_hist,decay_time]=hist(data,(.5):(ceil(max(data))-.5));
% fn=1/lifetime(ii)*exp(-decay_time/lifetime(ii))+back(ii)+coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii));
% fn=fn/sum(fn);
% 
% cm_sum=1-cumsum(decay_hist)/sum(decay_hist);
% cm_sum_fit=1-cumsum(fn);
% 
% % semilogy(decay_time,cm_sum,'s','MarkerEdgeColor','k','MarkerSize',6)
% % hold on
% % plot(decay_time,cm_sum_fit,'r','LineWidth',2.5);
% % grid on
% 
% semilogy(decay_time,decay_hist/sum(decay_hist),'s','MarkerEdgeColor','k','MarkerSize',6)
% hold on
% plot(decay_time,fn,'r','LineWidth',2.5);
% grid on
% % saveas(fg(jj),['Cumsum_fit_' num2str(jj)],'jpeg')
% end