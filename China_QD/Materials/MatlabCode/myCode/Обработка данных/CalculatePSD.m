function [f, S, nu, disp] = CalculatePSD(QD_num, QD_size, dir_name)
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_corr_data']);
load(['results/' dir_name 'QD_' QD_num '_' QD_size ...
        'nm/data_files/QD_' QD_num '_' QD_size '_nm_parsed_data']);

p_00 = corr_data.p_00;
err_00 = corr_data.err_00;
num_ch2 = corr_data.nums_of_ch2;
experiment_time = corr_data.experiment_time;

t = corr_data.times;
p = corr_data.g2 - 1;

new_bins = [ ...
    logspace(-(5-1/2),-4,2) ... 
    logspace(-(4-1/2),-3,2) ...
    logspace(-(3-1/4),-2,4) ...
    logspace(-(2-1/4),-1,4) ...
    logspace(-(1-1/5),0,5) ...
    logspace((0+1/10),1,10) ...
    logspace((1+1/15),2,10) ...
    logspace((2+1/10),log10(max(t)),5)...
    ];

num_new_bins = numel(new_bins);

new_p = new_bins*0;
new_t = new_p;
sigma_p = new_p; 

for ii = 2:num_new_bins
    kk = find(t <= new_bins(ii) & t > new_bins(ii-1));
    new_p(ii) = mean(p(kk));
    new_t(ii) = mean(t(kk));
    sigma_p(ii) = std(p(kk)) / sqrt(numel(kk));
end

kk = find(t <= new_bins(1) & t > 0);
new_p(1) = mean(p(kk));
new_t(1) = mean(t(kk));

PP = new_p(1:end);
tl = new_t(1:end);

f = logspace(log10(10/experiment_time/0.75), 5, 100);

char_time = 1./f;
S = f*0;

nm_P_first = p_00;
sigma_p0 = err_00;
nu = f*0;
n_m = num_ch2 / experiment_time;

for ii = 1:numel(f)
    bin_time = char_time(ii) / 100;
    time_len = char_time(ii) * 20; 
    binning = bin_time : bin_time : time_len/2;

    Ps = interp1([0 tl], [nm_P_first PP], binning, 'pchip');

    N_ps = (Ps+1).*bin_time*n_m^2*experiment_time;
    sigma_interp = sqrt(N_ps)./bin_time/n_m^2/experiment_time;

    a0 = 0.3635819;
    a1 = 0.4891775;
    a2 = 0.1365995;
    a3 = 0.0106411;

    w = a0 - ...
        a1*cos(2*pi*binning/time_len+pi) + ...
        a2*cos(4*pi*binning/time_len) - ...
        a3*cos(6*pi*binning/time_len+pi);
    S(ii) = 2*bin_time*(nm_P_first+2*sum(Ps .* w .* cos(2*pi*f(ii)*binning)));

    Sigma_from_cor(ii)=2^2*bin_time^2*(sigma_p0^2+2^2*sum(sigma_interp.^2.*(w.*cos(2*pi*f(ii)*binning)).^2));
    Sigma_full=sum(2*w.^2).*bin_time/experiment_time*S(ii)^2+Sigma_from_cor(ii);
    
    nu(ii) = 2*experiment_time/sum(2*w.^2)/bin_time;
    disp(ii) = Sigma_full;
end
kk = find(S < 0);
S(kk) = [];
f(kk) = [];
nu(kk) = [];
disp(kk) = [];
PSDPlotter(f, S, nu, disp, QD_num, QD_size, dir_name);
end

