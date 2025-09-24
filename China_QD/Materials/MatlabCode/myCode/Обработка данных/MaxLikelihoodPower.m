function ls=MaxLikelihoodPower(p)
global decay decay_time
tau=p(1);
back=p(2);
A2=p(3);
tau2=p(4);
% Likelihood=1/tau*exp(-decay_time'/tau)+back;

Likelihood=1/tau*exp(-decay_time/tau)+back+A2/tau2*exp(-decay_time/tau2);
Likelihood=Likelihood/sum(Likelihood);
ls=-sum(decay.*log(Likelihood));

if tau<0 || back<0 || tau2<0 || A2<0 || tau>tau2
    ls=inf;
end