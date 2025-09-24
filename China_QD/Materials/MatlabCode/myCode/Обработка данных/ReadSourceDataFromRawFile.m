function [] = ...
    ReadSourceDataFromRawFile(QD_num, QD_size, dir_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Подпрограмма для чтения данных из исходных файлов и их последующего
% бинирования. 
%
% Сохраняются:
%       Структура parsed_data, состоящая из 3 колонок: channel (1/0),
%           arrivals ([ps]), delays ([ps]).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Открытие файла с данными. Выкидывает исключение, если произошла ошибка
% открытия файла. 
fileID = fopen(['source_data/' dir_name QD_num]);
assert(fileID ~= -1, 'Ошибка открытия файла');

% Данные считываются из файла. Далее происходит их приведение к матрице для
% дальнейшей работы (decimal_data). Структура матрицы имеет вид:
% control bit || channel || delays (ps/64) || impulse number
source_data  = transpose(dec2hex(fread(fileID)));
reshape_data = transpose(reshape(source_data, 8, []));
binary_data  = dec2bin(hex2dec(reshape_data));
decimal_data = [bin2dec(binary_data(:, 1)), ...
                bin2dec(binary_data(:, 2  : 7)), ...
                bin2dec(binary_data(:, 8  : 22)), ...
                bin2dec(binary_data(:, 23 : end)), ...
                ];

% Поиск и обработка битов переполнения
indexes_overflow_bits = ...           % Индексы битов переполнения 
    find(decimal_data(:, 1) == 1);
nums_bins_between_overflow_bits = ... % Количество элементов между битами
    diff(indexes_overflow_bits)-1;
range_0_to_length = ...               % range от 0 до кол-ва всех фотонов
    0:length(nums_bins_between_overflow_bits)-1;
overflow_array = ...                  % Вектор с количеством переполнений
    transpose(repelem(range_0_to_length, nums_bins_between_overflow_bits));

% Удаляем строки с информацией о переполнении. Далее корректируем номера
% импульсов с учетом информации о количестве переполнений, полученной на
% предыдущем шаге алгоритма. 
decimal_data(indexes_overflow_bits, :) = [];
decimal_data(:, 4) = decimal_data(:, 4) + overflow_array * 1023;

% Корректируем времена задержки (в [ps]). Вычисляем время между импульсами.
decimal_data(:, 3) = decimal_data(:, 3) * 64;
time_between_impulses = max(decimal_data(:, 3));

% Формируем времена прилета фотонов по формуле:
% номер_импульса * время_между_импульсами + задержка_относительно_импульса
decimal_data(:, 1) = ...
    time_between_impulses * decimal_data(:, 4) + decimal_data(:, 3);

% Сортируем массив пузырьком для устранения ситуаций с фотонами, отстающими
% во времени.
while ~issorted(decimal_data(:, 1))
    for ii = 1 : length(decimal_data(:, 1))-1
        if decimal_data(ii, 1) > decimal_data(ii+1, 1)
            tmp = decimal_data(ii, :);
            decimal_data(ii, :)   = decimal_data(ii+1, :);
            decimal_data(ii+1, :) = tmp;
        end
    end
end

% Начало отсчета времен = 0
decimal_data(:, 1) = decimal_data(:, 1) - decimal_data(1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Формируем структуру сохраняемую
parsed_data = struct( ...
    'channel',  decimal_data(:, 2), ...
    'arrivals', decimal_data(:, 1), ...
    'delays',   decimal_data(:, 3) ...
    );

save(['results/' dir_name ...
        'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' ...
        QD_size '_nm_parsed_data'], 'parsed_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end