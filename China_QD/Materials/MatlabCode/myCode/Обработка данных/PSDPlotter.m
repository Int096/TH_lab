function [] = PSDPlotter(f, S, nu, disp, QD_num, QD_size, dir_name)
S_up = S*0;
S_dn = S*0;
Sigma_gs = sqrt(disp');

for ii = 1:numel(nu)
    if (nu(ii) < 100)
        S_up(ii) = nu(ii).*S(ii)./chi2inv(0.025, nu(ii));
        S_dn(ii) = nu(ii).*S(ii)./chi2inv(0.975, nu(ii));
    else
        S_up(ii) = S(ii) + 1.96*Sigma_gs(ii);
        S_dn(ii) = S(ii) - 1.96*Sigma_gs(ii);
    end
end

PSD_figure = figure(...
    'Color', 'white', ...
    'Units', 'normalized', ...
    'Name', 'PSD', ...
    'Position', [0.15 0.15 0.7 0.7] ...
    );
hold on
errorbar(f, S, S-S_dn, S_up-S, 's', 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
set(gca, 'YScale', 'log', 'XScale', 'log', ...
    'XLim', [min(f)/2 max(f)*2], ...
    'YLim', [0.5*min(S_dn) 2*max(S_up)]);

xlabel('Frequency, [Hz]');
ylabel('PSD, [Hz^{-1}]');
set(gca, 'FontSize', 24);

text(0.7, 0.92, [QD_num ' QD ' QD_size ' nm'], ...
    'FontSize', 24, ...
    'Color', 'black', ...
    'Units', 'normalized' ...
    );

saveas(gcf, ['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/images/QD_' QD_num '_' QD_size '_PSD_plot.png']);
end