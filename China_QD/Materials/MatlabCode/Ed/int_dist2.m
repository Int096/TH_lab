function ff=int_dist2(s_0,s,p,K,back,n)
global N alpha N_med conf_full lgg gamma0 lnprob coef Left_part Right_part lambda_new lambda_sum n_med1 n_med
% if (~isempty(find([s_0,s,p,back] <0))) || (~isempty(find(p>0.95 | p<0.001 | s<0.00001 | s>0.5))) || (s_0+sum(s))>1
if (~isempty(find([s_0,s,p,back] <0))) || (~isempty(find(p>0.999))) || (~isempty(find(s>0.25)))
    ff=inf;
else

S=(s_0+conf_full*s')';
p_st=prod(conf_full.*(ones(2^N,1)*p)+abs(conf_full-1).*(ones(2^N,1)*(1-p)),2)';
Yield=ones(1,2^N)./(1+K*S.^alpha);

Yield_med=p_st*Yield';
% beta0=(N_med-back)/Yield_med;
% n_med=beta0*Yield+back;
% n_med1=n_med;

lambda_sum=(conf_full*gamma0')';
% Q_right_transpose=cell(1,8);
% Q_left=cell(1,8);
for ii=1:N
    Q_right_transpose{ii}=transpose([[1-p(ii);p(ii)],[-1;1]]);
%     Q_left{ii}=[[1,1];[-p(ii),1-p(ii)]];
%     w_right(:,:,ii)=[[1-p(ii);p(ii)],[-1;1]];
    w_left(:,:,ii)=[[1,1];[-p(ii),1-p(ii)]];
end

% V_right=w_right(:,:,1);
V_left=w_left(:,:,1);
for jj=2:N
%     V_right=kron(V_right,w_right(:,:,jj));
    V_left=kron(V_left,w_left(:,:,jj));
end

lambda_new=1./(0.01*lambda_sum).*(1-exp(-0.01*lambda_sum));
lambda_new(lambda_sum==0)=1;


Left_part=kronm(fliplr(Q_right_transpose),Yield');
% Right_part=kronm(fliplr(Q_left),p_st');
% coef=(Left_part.*Right_part)';
% lambda_new=1./(0.01*lambda_sum).*(1-exp(-0.01*lambda_sum));
% lambda_new(lambda_sum==0)=1;
% Yield_med=sum(coef.*lambda_new);

beta0=(N_med-back)/Yield_med;
% Y_new=((Yield*V_right).*lambda_new)*V_left;
Y_new=(Left_part'.*lambda_new)*V_left;

% Y_new=Yield*(V_left*diag(lambda_new)*V_right);
% Y_new=((V_right*diag(lambda_new)*V_left)*Yield');
n_med=beta0*Y_new+back;

% Left_part=kronm(fliplr(Q_right_transpose),Yield');
% Right_part=kronm(fliplr(Q_left),p_st');
% coef=(Left_part.*Right_part)';
% lambda_new=1./(0.01*lambda_sum).*(1-exp(-0.01*lambda_sum));
% lambda_new(lambda_sum==0)=1;
% Yield_med=sum(coef.*lambda_new);
% beta0=(N_med-back)/Yield_med;
% n_med=beta0*Yield+back;

% lg_nmed=log(n_med);
% % lg_p_st=log(p_st);
% Right_part=kronm(fliplr(Q_left),p_st');
% % part=-n_med+lg_p_st;
% 
% ff=1e-60*ones(1,max(n)+1);
% 
% for kk=1:(max(n)+1)
%     k=kk-1;
%     lnprob1=k*lg_nmed-lgg(kk)-n_med;
%     lnprob=lnprob1';
%     Left_part=kronm(fliplr(Q_right_transpose),exp(lnprob));
%     coef=(Left_part.*Right_part)';
%     lambda_new=1./(0.01*lambda_sum).*(1-exp(-0.01*lambda_sum));
%     lambda_new(lambda_sum==0)=1;
%     ff(kk)=ff(kk)+sum(coef.*lambda_new);
%         
% end
% 
% 
% Left_part=kronm(fliplr(Q_right_transpose),Yield');

lg_nmed=log(n_med);
lg_p_st=log(p_st);
part=-n_med+lg_p_st;

ff=1e-60*ones(1,max(n)+1);

for kk=1:(max(n)+1)
    k=kk-1;
    lnprob=k*lg_nmed-lgg(kk)+part;     
    ff(kk)=ff(kk)+sum(exp(lnprob(lnprob > -30 & ~isnan(lnprob))));
end

end