function LG = Ln_Gathered(par)
global Y_exp X_exp N Norm_coef Log_norm freq S nu disp

s_0=par(1);
s=par(2:(N+1));
p=par((N+2):(2*N+1));
K=par(end-1);
back=par(end);

% s(6)=0;
% s(8)=0;
% s(10)=0;
% p(6)=0;
% p(8)=0;
% p(10)=0;

Y_t=int_dist2(s_0,s,p,K,back,X_exp);
if Y_t ~= inf
    ls=-Y_exp*log(Y_t)'*Norm_coef - Log_norm;
else
    ls=inf;
end

K_par=par(end-1);
back_par=par(end);

St=spectra(s_0,s,p,K_par,back_par,freq);
K=nu/2;
X=S./St.*K;


if isinf(St)
    LP=inf;
else
    kk=find(nu<100);
    LP=sum(-(K(kk)-1).*log(2*X(kk)) + X(kk) + gammaln(K(kk)) + K(kk)*log(2));
    jj=find(nu>100);
    LN = sum((S(jj)-St(jj)).^2./(2*disp(jj))) + sum(log(sqrt(2*pi*jj./jj)));
    LP=LP + LN;
end
LG=LP+ls
