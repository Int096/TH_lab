global N alpha gamma0 Y_exp X_exp N_med param data conf_full lgg noise name N_med_2 Norm_coef Log_norm freq S nu disp

name='QD6_1MHz_580_100nW-8_008_new';
N=10;
alpha=10;

jj_full=1:2^N;
conf_full=dec2bin(jj_full-1,N)-48;

raw_data=load(['Trace/' name '.txt'],'-ascii'); 
data=raw_data(1:end,2);

[Y_exp,X_exp]=hist(data,min(data):max(data));
Norm_coef=sum(Y_exp);
gamma_ln_y=gammaln(Y_exp+1);
gamma_ln_y(gamma_ln_y==inf)=0;
Log_norm=gammaln(Norm_coef+1)-sum(gamma_ln_y);
Y_exp=Y_exp/sum(Y_exp);
N_med=Y_exp*X_exp';

zz=0:max(X_exp);
Z=find(zz<100);
lgg(Z)=log(gamma(zz(Z)+1));
Z=find(zz>=100);
lgg(Z)=gammaln(zz(Z)+1);

% load('FINAL_SPECS_PARAMS\QD8_2.mat')
% gamma1=logspace(log10(gamma0(1)),log10(gamma0(end)),N);
% gamma0=gamma1;

% load('Final_full_fit_params/qd36_with_8_TLS');
gamma0=logspace(-2.4,4.5,N);

% % % gamma0=logspace(-1.2,3.5,N); так было с лучшим фитом через x0 qd82

time_and_cor=load(['Kor/' name '.txt']);
time=time_and_cor(:,1)*1e-9;
cor_fun=time_and_cor(:,2)-1;
[freq,S,nu,disp]=PSD_Black(time,cor_fun); % output PSD

% s_00=param_tot(1);
% p_0=[param_tot(10:17) .5 .5];
% s_0=[param_tot(2:9) 1e-4 1e-4];
% K_0=param_tot(end-1);
% back=param_tot(end);

s_00=0.1*rand(1);
p_0=0.7*rand(1,N);
s_0=0.1*rand(1,N);
K_0=1e7;
back0=0;

x0=[s_00,s_0,p_0,K_0,back0];

% gamma0(9)=[];
% gamma0(2)=[];
% % gamma0(6)=[];
% param(20)=[];
% param(13)=[];
% param(10)=[];
% param(3)=[];
% % % % param(9)=[];
% % % % param(7)=[];

x0=param;

param=fminsearch(@Ln_Gathered,x0,optimset('MaxFunEvals',50000,'MaxIter',50000));

% param(7)=0;
% param(9)=0;
% param(11)=0;
% param(17)=0;
% param(19)=0;
% param(21)=0;

subplot(1,2,1)
PSDplot();
hold on
spec_total=spectra(param(1),param(2:(N+1)),param((N+2):(2*N+1)),param(end-1),param(end),freq);
loglog(freq,spec_total,'r','LineWidth',2.5);

subplot(1,2,2)
hold on
bar(X_exp,Y_exp,1,'hist');
int_tot=int_dist2(param(1),param(2:(N+1)),param((N+2):(2*N+1)),param(end-1),param(end),X_exp);
plot(X_exp,int_tot,'r','LineWidth',3.5);
axis([min(X_exp) max(X_exp) 0 1.1*max(Y_exp)]);

% save('Final_full_fit_params/qd68_without_v3','param','gamma0')

% save('qd36_with_8_TLS','param','gamma0');