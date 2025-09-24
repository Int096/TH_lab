x=input('Number of dot:');
foldername=['china_622_' num2str(x) '_spec'];
mkdir(foldername);
fid=fopen(num2str(x));
hexStr=dec2hex(fread(fid));
trhexstr=hexStr';
StrReshape=reshape(trhexstr,8,[])';
s=size(StrReshape);
rws=s(1);
binary1=dec2bin(hex2dec(StrReshape(1:floor(rws/3),:)));
binary2=dec2bin(hex2dec(StrReshape(floor(rws/3)+1:2*floor(rws/3),:)));
binary3=dec2bin(hex2dec(StrReshape(2*floor(rws/3)+1:end,:)));
binary=[binary1; binary2; binary3];
% binary=dec2bin(hex2dec(StrReshape));
strs=size(binary);
M=[];
T=0;
M(:,2)=bin2dec(binary1(:, 8:22));
M(:,3)=bin2dec(binary1(:, 23:end));
M(:,1)=bin2dec(binary1(:, 2:7));
M(1,:)=[];
M=[M; [bin2dec(binary2(:, 2:7)), bin2dec(binary2(:, 8:22)), bin2dec(binary2(:, 23:end))]];
M=[M; [bin2dec(binary3(:, 2:7)), bin2dec(binary3(:, 8:22)), bin2dec(binary3(:, 23:end))]];
% M(:,2)=[M(:,2); bin2dec(binary2(:, 8:22))];
% M(:,3)=[M(:,3); bin2dec(binary2(:, 23:end))];
% M(:,1)=[M(:,1); bin2dec(binary2(:, 2:7))];
M(end,:)=[];
indexes_reset=(find(M(:,1)==63))';
reset_vector=0:1:length(indexes_reset)-1;
num_of_repeats=[indexes_reset(1)-1 diff(indexes_reset)];
bonus_1023=(repelem(reset_vector,num_of_repeats))';
bonus_1023=[bonus_1023; length(indexes_reset);];
if (length(M(:,3))-indexes_reset(end))~=0
    bonus_1023=[bonus_1023; length(indexes_reset)*ones(length(M(:,3))-indexes_reset(end),1)];
end;
M(indexes_reset,3)=0;
M(:,3)=M(:,3)+bonus_1023*1023;
time_between=max(M(:,2))*64;
ch_number=M(:,3);
M(:,3)=time_between*ch_number+M(:,2);
savename=['china_622_' num2str(x) '_spec/china_622_' num2str(x)];
savename1=[savename '_pure_numbers_dot'];
save(savename1, 'M');

global decay decay_time
load(savename1);

flags=M(:,1);
delays=M(:,2)*64*1e-12;
synch=M(:,3);

delays(flags==63)=[];
synch(flags==63)=[];
flags(flags==63)=[];
mx_del=max(delays);

[y_del,x_del]=hist(delays*1e9,(0:1:fix(mx_del*1e9))+.5);
% bar(x_del,y_del);

max_num=find(y_del==max(y_del));
max_x=x_del(max_num)-.5;

delays=delays-max_x*1e-9;

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

% [time, occurance]=hist(Intensity, num_of_bins);
% figure();
% plot(time, occurance);
% xlabel('Occurrence');
% ylabel('Intensity');
% grid on;

% % hhh=figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% % plot((1:num_of_bins)*bin_size,Intensity);
% % xlabel('Time [s]')
% % ylabel('cnts/ms')
% % set(gca,'FontSize',24)

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
fit_inital_tau=linspace(0, 10, NumOfLevels);
fit_initial_int=(max(av_int_levels2)/10)*fit_inital_tau;
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
%         x0=[10 1e-3 1 100];
%         x0=[10 1e-3 1e-5 100];
%         if fit_initial_int(ii)==av_int_levels2(ii)
        x0=[av_int_levels2(ii)/(max(av_int_levels2)*10) 1e-3 1e-5 100];
        [param,fval,exit_flag]=fminsearch(@China_Ln,x0,optimset('MaxFunEvals',10000,'MaxIter',10000));
        while exit_flag==0
            [param,fval,exit_flag]=fminsearch(@China_Ln,param,optimset('MaxFunEvals',10000,'MaxIter',10000));
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
        if param(1)<param(4)
            fn=1/tau*exp(-decay_time/tau)+back2+A2/tau2*exp(-decay_time/tau2);
%             fn=A2/tau*exp(-decay_time/tau)+A2/tau2*exp(-decay_time/tau2)+back2;
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
        else
            fn=1/tau*exp(-decay_time/tau)+back2+A2/tau2*exp(-decay_time/tau2);
%             fn=A2/tau*exp(-decay_time/tau)+A2/tau2*exp(-decay_time/tau2)+back2;
            percent_long=sum(1/tau*exp(-decay_time/tau))/sum(fn);
            percent_short=sum(A2/tau2*exp(-decay_time/tau2))/sum(fn);
            lifetime(count)=param(4);
            back(count)=param(2);
            coef(count)=param(3);
             dec_time2(count)=param(1);
            intensityFLID(count)=av_int_levels(ii);
            intensityFLID2(count)=av_int_levels2(ii);      
            intensityFLIDshort(count)=av_int_levels2(ii)*percent_short;
            intensityFLIDlong(count)=av_int_levels2(ii)*percent_long;
            perc(count)=percent_short;
            perc_l(count)=percent_long;
        end
    end
end
% save('china_622_97_spec/China_flid_data_97_dot', 'lifetime', 'back', 'coef', 'dec_time2', 'intensityFLID', 'intensityFLID2', 'intensityFLIDshort', 'intensityFLIDlong', 'perc', 'perc_l');   
% for jj=1:numel(lifetime)
%     dc=dec_time2(jj);
%     lf=lifetime(jj);
%     if dc<lf
%         lifetime(jj)=dc;
%         dec_time2(jj)=lf;
%     end
% end

fg1 = figure('Color','white','Units','normalized','position',[0.15 0.15 0.7 0.7]);
% lifetime_add=lifetime(dec_time2<300);
% dec_time2_add=dec_time2(dec_time2<300);
% intensityFLID_add=intensityFLID2(dec_time2<300);
% dec_time2_add(lifetime_add>45)=[];
% intensityFLID_add(lifetime_add>45)=[];
% lifetime_add(lifetime_add>45)=[];
% intensityFLID_add=intensityFLID_add+Intensity_borders(2);

lifetime_add=lifetime(dec_time2<300);
intensityFLIDshort_add=intensityFLIDshort(dec_time2<300);
intensityFLIDlong_add=intensityFLIDlong(dec_time2<300);
intensityFLID_add=intensityFLID2(dec_time2<300);
dec_time2_add=dec_time2(dec_time2<300);


dec_time2_add(lifetime_add>45)=[];
intensityFLIDshort_add(lifetime_add>45)=[];
intensityFLIDlong_add(lifetime_add>45)=[];
lifetime_add(lifetime_add>45)=[];
intensityFLID_add(lifetime_add>45)=[];

dec_time2_add(intensityFLIDshort_add>150)=[];
intensityFLIDlong_add(intensityFLIDshort_add>150)=[];
lifetime_add(intensityFLIDshort_add>150)=[];
intensityFLID_add(intensityFLIDshort_add>150)=[];
intensityFLIDshort_add(intensityFLIDshort_add>150)=[];

lifetime_add=lifetime(intensityFLIDshort<100);
dec_time2_add=dec_time2(intensityFLIDshort<100);
intensityFLIDlong_add=intensityFLIDlong(intensityFLIDshort<100);
intensityFLID_add=intensityFLID2(intensityFLIDshort<100);
intensityFLIDshort_add=intensityFLIDshort(intensityFLIDshort<100);


% lifetime=lifetime_add;
% dec_time2=dec_time2_add;
% intensityFLIDshort=intensityFLIDshort_add;
% intensityFLIDlong=intensityFLIDlong_add;
% intensityFLID2=intensityFLID_add;
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
hist3([lifetime_add' intensityFLIDshort_add'],[100 100]) %intensityFLIDshort_add
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
hold on;
xlabel('Short lifetime, [ns]');
ylabel('cts per 1 ms');
set(gca,'FontSize',24);
set(gca, 'XLim', [0 max(lifetime_add)], 'YLim', [0 max(intensityFLIDshort_add)]);
% saveas(fg2, 'china_622_97_spec/China_flid_LI_97_dot', 'jpeg'); %95 заново сделать 
% % % colormap parula
% % % % axis([0 max(lifetime_add) 0 max(intensityFLIDshort_add)])
% % % axis([0 max(lifetime_add) 0 80])
% % % % axis([0 50 0 200])
% % % cmap=colormap;
% % % cmap(1,:)=[0 0 0.49];
% % % colormap(cmap)
% % % cb=colorbar('location','east','Color','w');
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

lifetime=lifetime(dec_time2<300);
intensityFLIDshort=intensityFLIDshort(dec_time2<300);
intensityFLIDlong=intensityFLIDlong(dec_time2<300);
intensityFLID=intensityFLID2(dec_time2<300);
dec_time2=dec_time2(dec_time2<300);


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
