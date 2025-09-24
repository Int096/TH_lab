function [] = CalculateFLID_new(QD_num, QD_size, dir_name)

load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_parsed_data']);
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_binned_data']);

Intensity = binned_data.intensity;
index_min = binned_data.index_min;
index_max = binned_data.index_max;
delays = parsed_data.delays;

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
plot(dec_time2,intensityFLIDlong,'s','MarkerEdgeColor','k')
end

