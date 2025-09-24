function [] = CorrelateFunctionPlotter(QD_num, QD_size, dir_name)
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_corr_data']);

p = corr_data.g2 - 1;
t = corr_data.times;

new_p = corr_data.new_p;
new_t = corr_data.new_t;
new_N_kor = corr_data.new_N_kor;
experiment_time = corr_data.experiment_time;
new_delta_t = corr_data.new_delta_t;
n_m = corr_data.n_m;


corr_figure  = figure('Color', 'white', ...
                      'Units', 'normalized', ...
                      'Name', 'CorrFunction plot', ...
                      'Position', [0.2 0.2 0.7 0.7] ...
                      );
set(gca, 'Xscale', 'log');

hold on;
plot(t, p);

hold on
errorbar(new_t, new_p, ...
    (new_N_kor-1/2*chi2inv(0.025,2*new_N_kor))./new_delta_t/n_m^2/experiment_time, ...
    (1/2*chi2inv(0.975,2*new_N_kor+2)-new_N_kor)./new_delta_t/n_m^2/experiment_time, ...
    's', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Color', 'k');

hold on
grid on
box on

xlabel('Time [s]', 'FontSize', 24);
ylabel('Autocorrelation function', 'FontSize', 24);
text(0.7, 0.92, [QD_num ' QD ' QD_size ' nm'], ...
    'FontSize', 24, ...
    'Color', 'black', ...
    'Units', 'normalized' ...
    );

set(gca, 'YLim', [0.99*min(p) max(p)], ...
         'XLim', [0.8 * min(t) 1.02*max(t)], ...
         'FontSize', 24);

saveas(corr_figure, ['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/images/QD_' QD_num '_' QD_size '_correlate_function_plot.png']);
end