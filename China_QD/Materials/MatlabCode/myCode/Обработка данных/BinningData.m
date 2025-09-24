function [] = BinningData(QD_num, QD_size, dir_name, bin_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Подпрограмма для бинирования данных. Для использования необходим
% parsed_data/. Ширина бина bin_size задается в секундах.
%
% Сохраняются: 
%       Структура binned_data, состоящая из 3 колонок: 
%           intensity (a.u.)
%           index_min (a.u.)
%           index_max (a.u.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_parsed_data']);

arrivals = parsed_data.arrivals * 1e-12;

% Вспомогательные данные для цикла
nums_of_arrivals = numel(arrivals);
nums_of_bin      = ceil(arrivals(end) / bin_size);

% Выделяем память заранее для ускорения
intensity = zeros(1, nums_of_bin);
index_min = intensity;
index_max = intensity; 

% Алгоритм подсчета фотонов в бине и индекса первого и последнего фотона.
% (используется тот факт, что времена прилетов отсортированы по
% возрастанию).
for ii = 1 : nums_of_bin
    if ii == 1
        index_min(ii) = 1;
        index         = 1;
    else
        index_min(ii) = index_max(ii-1)+1;
        index         = index_min(ii);
    end

    while arrivals(index) < ii*bin_size
        index = index + 1;

        if index+1 > nums_of_arrivals
            break
        end
    end

    if index == index_min(ii)
        intensity(ii) = 0;
        index_max(ii) = index_min(ii) - 1;
    else
        index_max(ii) = index - 1;
        intensity(ii) = index_max(ii) - index_min(ii) + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Формируем структуру сохраняемую
binned_data = struct( ...
    'intensity', intensity, ...
    'index_min', index_min, ...
    'index_max', index_max  ...
    );

save(['results/' dir_name ...
        'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' ...
        QD_size '_nm_binned_data'], 'binned_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

