function [] = MakeDir(QD_num, QD_size, dir_name)
    % Общий каталог
    mkdir(['results/' dir_name ...
        'QD_' QD_num '_' QD_size 'nm']);

    % Каталог для сохраненных в расчетах данных
    mkdir(['results/' dir_name ...
        'QD_' QD_num '_' QD_size 'nm/data_files']);

    % Каталог для изображений
    mkdir(['results/' dir_name ...
        'QD_' QD_num '_' QD_size 'nm/images']);
end