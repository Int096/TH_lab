global decay decay_time
load('china_622_2_pure_numbers_dot.mat')

flags=M(:,1);
delays=M(:,2)*64*1e-12;
synch=M(:,3);

delays(flags==63)=[];
synch(flags==63)=[];
flags(flags==63)=[];

delays1=delays(flags==0);
delays2=delays(flags==1);
mx_del=max(delays);

[y_del,x_del]=hist(delays1*1e9,(0:1:fix(mx_del*1e9))+.5);
% bar(x_del,y_del);

max_num=find(y_del==max(y_del));
max_x=x_del(max_num)-.5;

% delays=delays-max_x*1e-9;

delays(flags==0)=delays1-max_x*1e-9;

[y_del,x_del]=hist(delays2*1e9,(0:1:fix(mx_del*1e9))+.5);
% bar(x_del,y_del);

max_num=find(y_del==max(y_del));
max_x=x_del(max_num)-.5;
delays(flags==1)=delays2-max_x*1e-9;

synch(delays<0)=[];
flags(delays<0)=[];
delays(delays<0)=[];

arrivals=(synch*2e-7+delays);
arrivals=arrivals-arrivals(1);

nm_arrivals=numel(arrivals);
bin_size=0.001;

num_of_bins=ceil(arrivals(end)/bin_size);
Intensity=zeros(1,num_of_bins);
index_min=Intensity;
index_max=Intensity;

for jj=1:num_of_bins
    jj/num_of_bins*100
    if jj==1
        index_min(jj)=1;
        count=1;
    else
        index_min(jj)=index_max(jj-1)+1;
        count=index_min(jj);
    end
    while arrivals(count)<jj*bin_size
        count=count+1;
        if count+1>nm_arrivals
            break
        end
    end
    if count==index_min(jj)
        Intensity(jj)=0;
        index_max(jj)=index_min(jj)-1;
    else
        index_max(jj)=count-1;
        Intensity(jj)=index_max(jj)-index_min(jj)+1;
    end
end

hhh=figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
plot((1:num_of_bins)*bin_size,Intensity);
xlabel('Time [s]')
ylabel('cnts/ms')
set(gca,'FontSize',24)

NumOfLevels=100; %%100
Intensity_borders=linspace(min(Intensity),max(Intensity),NumOfLevels+1);
av_int_levels=(Intensity_borders(1:end-1)+Intensity_borders(2:end))/2;
av_int_levels2=Intensity_borders(1:end-1);

Arrays_cell=cell(1,NumOfLevels);
for ii=1:NumOfLevels
    if ii==1
        kk=find(Intensity<=Intensity_borders(ii+1));
    else
        kk=find(Intensity<=Intensity_borders(ii+1)&Intensity>Intensity_borders(ii));
    end
    Arrays_cell(1,ii)={kk};
end

Indexes_cell=cell(1,NumOfLevels);
for ii=1:NumOfLevels
    ii
    Intensity_indexes=Arrays_cell{1,ii};
    PhotonIndexes=[];
    for jj=Intensity_indexes
        PhotonIndexes=[PhotonIndexes index_min(jj):1:index_max(jj)];
    end
    Indexes_cell(1,ii)={PhotonIndexes};
end

% save('Indexes_cell_1','Indexes_cell');

% load('Indexes_cell_1')
numOfPhPerFit=1000;
count=0;
abc=0;
for ii=1:NumOfLevels
    ii
    indexes_in_level=Indexes_cell{ii};
    nm_ph=numel(indexes_in_level);
    num_of_fits=fix(nm_ph/numOfPhPerFit);
    for jj=1:num_of_fits
        count=count+1;
        if jj~=num_of_fits
            data_indexes=indexes_in_level(((jj-1)*numOfPhPerFit+1):jj*numOfPhPerFit);
            data=delays(data_indexes)*1e9;
        else
            data_indexes=indexes_in_level(((jj-1)*numOfPhPerFit+1):end);
            data=delays(data_indexes)*1e9;
        end
        [decay_hist,decay_time]=hist(data,0.5:1:ceil(max(data)));
        decay=decay_hist/sum(decay_hist);
%         decay_time=decay_time;
%         decay_time(decay==0)=[];
%         decay(decay==0)=[];
        x0=[10 1e-3 1 100];
        [param,fval,exit_flag]=fminsearch(@MaxLikelihoodPower,x0,optimset('MaxFunEvals',10000,'MaxIter',10000));
        while exit_flag==0
            [param,fval,exit_flag]=fminsearch(@MaxLikelihoodPower,param,optimset('MaxFunEvals',10000,'MaxIter',10000));
        end
        tau=param(1);
        back2=param(2);
        A2=param(3);
        tau2=param(4);
        
        %         if mod(count,20)==0
        %             figure()
        %             cm=1-cumsum(decay)/sum(decay);
        %             plot(decay_time,cm,'o')
        %             hold on
        %             plot(decay_time,1-cumsum(fn)/sum(fn),'r')
        %             pause(0.01)
        %         end
            fn=1/tau*exp(-decay_time/tau)+A2/tau2*exp(-decay_time/tau2)+back2;
            percent_short=sum(1/tau*exp(-decay_time/tau))/sum(fn);
            percent_long=sum(A2/tau2*exp(-decay_time/tau2))/sum(fn);
            lifetime(count)=param(1);
            back(count)=param(2);
            coef(count)=param(3);
            dec_time2(count)=param(4);
            intensityFLID(count)=av_int_levels(ii);
            intensityFLID2(count)=av_int_levels2(ii);
            intensityFLIDshort(count)=av_int_levels2(ii)*percent_short;
            intensityFLIDlong(count)=av_int_levels2(ii)*percent_long;
            perc(count)=percent_short;
            perc_l(count)=percent_long;
            
%             if intensityFLIDlong(count)>20
%             x0=[tau2 back2];
%             [param,fval,exit_flag]=fminsearch(@MaxLikelihoodMono,x0,optimset('MaxFunEvals',10000,'MaxIter',10000));
%             while exit_flag==0
%                 [param,fval,exit_flag]=fminsearch(@MaxLikelihoodMono,param,optimset('MaxFunEvals',10000,'MaxIter',10000));
%             end
%             
%             tau=param(1);
%             back2=param(2);
%             fn=1/tau*exp(-decay_time/tau)+back2;
%             
%             lifetime(count)=param(1);
%             back(count)=param(2);
%             coef(count)=NaN;
%             dec_time2(count)=200;
%             intensityFLIDshort(count)=av_int_levels2(ii);
%             intensityFLIDlong(count)=0;
%             end

%         if intensityFLIDlong(count)>50
%             if mod(abc,50)==0
%             figure()
%             semilogy(decay_time,fn/sum(fn),'r')
%             hold on
%             plot(decay_time,decay/sum(decay),'o')
%             end
%             abc=abc+1;
%         end
    end
end

% for jj=1:numel(lifetime)
%     dc=dec_time2(jj);
%     lf=lifetime(jj);
%     if dc<lf
%         lifetime(jj)=dc;
%         dec_time2(jj)=lf;
%     end
% end

fg1 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
lifetime_add=lifetime(dec_time2<300);
dec_time2_add=dec_time2(dec_time2<300);
intensityFLID_add=intensityFLID2(dec_time2<300);
dec_time2_add(lifetime_add>45)=[];
intensityFLID_add(lifetime_add>45)=[];
lifetime_add(lifetime_add>45)=[];
intensityFLID_add=intensityFLID_add+Intensity_borders(2);

lifetime_add=lifetime(dec_time2<300);
dec_time2_add=dec_time2(dec_time2<300);
intensityFLIDshort_add=intensityFLIDshort(dec_time2<300);
intensityFLIDlong_add=intensityFLIDlong(dec_time2<300);
intensityFLID_add=intensityFLID2(dec_time2<300);
% % % 
dec_time2_add(lifetime_add>45)=[];
intensityFLIDshort_add(lifetime_add>45)=[];
intensityFLIDlong_add(lifetime_add>45)=[];
intensityFLID_add(lifetime_add>45)=[];
lifetime_add(lifetime_add>45)=[];


lifetime=lifetime_add;
dec_time2=dec_time2_add;
intensityFLIDshort=intensityFLIDshort_add;
intensityFLIDlong=intensityFLIDlong_add;
intensityFLID2=intensityFLID_add;
% save('Model_for_fit_v2','lifetime','dec_time2','intensityFLIDshort','intensityFLIDlong','intensityFLID2')

hist3([lifetime_add' dec_time2_add'],[100 100])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
hold on
xlabel('Short lifetime, [ns]');
ylabel('Long lifetime, [ns]');
set(gca,'FontSize',24);

% % colormap parula
% % axis([0 max(lifetime_add) 0 max(dec_time2_add)])
% % % axis([0 45 0 300])
% % cmap=colormap;
% % cmap(1,:)=[0 0 0.49];
% % colormap(cmap)
% % cb=colorbar('location','east','Color','w');
% % hold on
% % xx=[0 50];
% % yy=10*(xx);
% % plot3(xx,yy,[1000 1000],'r','LineWidth',2.5)
% % legend('','\tau_{del} = 10 (\tau-1)')

% plot(lifetime,dec_time2,'s','MarkerEdgeColor','k')
% set(gca,'XLim',[0 50],'YLim',[0 300])

fg2 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% hist3([lifetime_add' intensityFLID_add'],[50 45])
ctrs={.5:44.5 .5:2:79.5};
hist3([lifetime_add' (intensityFLIDshort_add+intensityFLIDlong_add)'],[80 80])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
hold on
xlabel('Short lifetime, [ns]');
ylabel('cts per 1 ms');
set(gca,'FontSize',24);
set(gca,'Xlim',[0 max(lifetime_add)],'Ylim',[0 max(intensityFLIDshort_add+intensityFLIDlong_add)])
plot3(0:45,14*(0:45),1000*ones(1,46),'r','lineWidth',3)
plot3(0:45,28*(0:45),1000*ones(1,46),'r','lineWidth',3)
% % % colormap parula
% % % % axis([0 max(lifetime_add) 0 max(intensityFLIDshort_add)])
% % % axis([0 max(lifetime_add) 0 80])
% % % % axis([0 50 0 200])
cmap=colormap;
cmap(1,:)=[0 0 0.49];
colormap(cmap)
cb=colorbar('location','east','Color','w');
saveas(gcf,'Hist_full_int','jpeg')
% % % xx=[0 25];
% % % yy=8*xx;
% % % hold on
% % % plot3(xx,yy,[1000 1000],'r','LineWidth',2.5)
% % % xx=[0 25];
% % % yy=10*(xx-2);
% % % plot3(xx,yy,[1000 1000],'g','LineWidth',2.5)
% % % legend('','N = 8 \tau','N = 10 (\tau-2)')
% plot(lifetime,intensityFLIDshort,'s','MarkerEdgeColor','k')
% set(gca,'XLim',[0 50],'YLim',[0 40])

fg3 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% % % ctrs={2.5:5:300  1.5:3:185};
% % % hist3([dec_time2_add' intensityFLIDlong_add'],ctrs)
% % % set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% % % view(2)
hold on
xlabel('Long lifetime, [ns]');
ylabel('cts per 10 ms');
% set(gca,'FontSize',24);
% % % colormap parula
% % % % axis([0 max(dec_time2_add) 0 max(intensityFLIDlong_add)])
% % % axis([0 max(dec_time2_add) 0 185])
% % % % axis([0 350 0 200])
% cmap=colormap;
% cmap(1,:)=[0 0 0.49];
% colormap(cmap)
% cb=colorbar('location','east','Color','w');
% % % xx=[0 200];
% % % yy=0.9*xx;
% % % hold on
% % % plot3(xx,yy,[1000 1000],'r','LineWidth',2.5)
% % % legend('','N = 0.9 \tau_{del}')
plot(dec_time2,intensityFLIDlong,'s','MarkerEdgeColor','k')
% set(gca,'XLim',[0 200],'YLim',[0 40])
% 
% saveas(fg1,'Both lifitemes_100ps3_per_135','jpeg')
% saveas(fg2,'FLID_short_100ps3_per_135','jpeg')
% saveas(fg3,'FLID_delayed_100ps3_per_135','jpeg')
% save('Subst_prob')
