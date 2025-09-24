function [] = CalculateFLID(QD_num, QD_size, dir_name)
global decay_time decay

load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_parsed_data']);
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_binned_data']);

delays        = parsed_data.delays * 1e-2;
arrivals      = parsed_data.arrivals * 1e-12;
detector_flag = parsed_data.channel;

first_dec_indexes  = find(detector_flag == 0);
second_dec_indexes = find(detector_flag == 1);

first_dec_delay  = delays(first_dec_indexes);
second_dec_delay = delays(second_dec_indexes);

[decay1, decay1_time] = ...
    hist(first_dec_delay, (.5):(ceil(max(first_dec_delay))-.5));
[decay2, decay2_time] = ...
    hist(second_dec_delay, (.5):(ceil(max(second_dec_delay))-.5));

[trash, decay1_max_idx] = max(decay1);
[trash, decay2_max_idx] = max(decay2);

delays(first_dec_indexes)  = ...
    delays(first_dec_indexes) - decay1_time(decay1_max_idx);
delays(second_dec_indexes) = ...
    delays(second_dec_indexes) - decay2_time(decay2_max_idx);

delays = delays * 1e-10;
time_between = max(delays);

arrivals(delays < 0) = [];
delays(delays < 0) = [];
arrivals(delays > time_between) = [];
delays(delays > time_between) = [];

arrivals = arrivals - arrivals(1);
nm_arrivals = numel(arrivals);

num_ph_bin = 1000;
n_max = fix(numel(delays) / num_ph_bin);

delay_matrix = zeros(n_max, num_ph_bin);

for ii = 1:n_max
    k = ((ii-1)*num_ph_bin+1):ii*num_ph_bin;
    delay_matrix(ii, :) = delays(k) * 1e9;
    right_time(ii) = arrivals(k(end));
end

lifetime = ones(1, n_max);
back = lifetime;
dec_time2 = lifetime;
coef = lifetime;

for ii = 1:n_max
    data = delay_matrix(ii, :);
    [decay_hist, decay_time] = hist(data, (.5):(ceil(max(data))-.5));
    decay = decay_hist / sum(decay_hist);
    x0 = [10 1e-3 1 100];
    param = fminsearch(@MaxLikelihoodPower, x0, ...
        optimset('MaxFunEvals', 10000, 'MaxIter', 10000));

    lifetime(ii) = param(1);
    back(ii) = param(2);
    coef(ii) = param(3);
    dec_time2(ii) = param(4);
    fn = 1/lifetime(ii) * exp(-decay_time/lifetime(ii)) + back(ii) + ...
        coef(ii)/dec_time2(ii) * exp(-decay_time/dec_time2(ii));
    perc(ii) = sum(1/lifetime(ii)*exp(-decay_time/lifetime(ii)))/sum(fn);
    percl(ii) = sum(coef(ii)/dec_time2(ii)*exp(-decay_time/dec_time2(ii)))/sum(fn);
    fn = fn/sum(fn);

    if ii > 1
        intensity_1000 = 1000 / (right_time(ii) - right_time(ii-1)) / 100;
    else
        intensity_1000 = 1000 / (right_time(ii)) / 100;
    end
    inten_short(ii) = intensity_1000 * perc(ii);
    inten_long(ii) = intensity_1000 * percl(ii);
end

FLID_data = struct('lifetime', lifetime, ...
                   'inten_short', inten_short, ...
                   'inten_long', inten_long, ...
                   'back', back, ...
                   'dec_time2', dec_time2, ...
                   'coef', coef, ...
                   'right_time', right_time);

save(['results/' dir_name ...
        'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' ...
        QD_size '_nm_FLID_data'], 'FLID_data');
end