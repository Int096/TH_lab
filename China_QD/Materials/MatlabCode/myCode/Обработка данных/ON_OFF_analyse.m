function [] = ON_OFF_analyse(QD_num, QD_size, dir_name)
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_binned_data']);

intensity = binned_data.intensity;

ON_state = 1;
ON_times = [0];
OFF_times = [0];
time = 0;
level = 0.5 * max(intensity)
for ii = 1:numel(intensity)
    if ON_state == 1 && intensity(ii) >= level
        time = time + 1
    else if ON_state == 1 && intensity(ii) < level
        
end
end

