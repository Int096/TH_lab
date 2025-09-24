% Открытие файла с данными
filename = num2str(input('Number of dot: '));
fid = fopen(filename);
if fid == -1
    fprintf("File %s don`t exist!!!", filename);
    return;
end

% Чтение и подгонка под формат
data_source = transpose(dec2hex(fread(fid)));
data_reshape = transpose(reshape(data_source, 8, []));
data_binary = dec2bin(hex2dec(data_reshape));

% Формируем необходимые данные 
% SPSN bits | CH bits | time bits
% =>
% photon_time | CH number | time_delay
% photon_time = imp_number*time_between + time_delay
data_result = [bin2dec(data_binary(:, 23:end)), bin2dec(data_binary(:, 2:7)), bin2dec(data_binary(:, 8:22))];
data_result = data_result(2:end-1, :);

% Обработка переполнения тактового сигнала (оставил как есть)
idx_reset = transpose(find(data_result(:, 2) == 63));
vec_reset = 0:1:length(idx_reset)-1;
repeats = [idx_reset(1)-1 diff(idx_reset)];
overflow = transpose(repelem(vec_reset, repeats));
overflow = [overflow; length(idx_reset);];
if (length(data_result(:, 1))-idx_reset(end)) ~= 0
    overflow = [overflow; lenght(idx_reset)*ones(length(data_result(:, 1))-idx_reset(end), 1)];
end
data_result(:, 1) = data_result(:, 1) + overflow*1023;
data_result(idx_reset, :) = [];
%%%

time_between = max(data_result(:, 3)) * 64;
data_result(:, 1) = time_between * data_result(:, 1) + data_result(:, 3);

% Сохранение всего непотребства
foldername = ['China_622_' filename '_spec'];
mkdir(foldername);
savename = [foldername '/China_622_' filename '_pure_numbers_dot'];
save(savename, 'data_result');
