function sp=spectra(s_0,s,p,K,back,x)
global N alpha N_med gamma0 noise conf_full coef data N_med_2
    S=(s_0+conf_full*s')';
    p_st=prod(conf_full.*(ones(2^N,1)*p)+abs(conf_full-1).*(ones(2^N,1)*(1-p)),2)';
    Yield=ones(1,2^N)./(1+K*S.^alpha);
    Yield_ps=Yield.*p_st;
    Yield_med=p_st*Yield';
%     beta0=(N_med-back)/Yield_med;
%     B=beta0^2/N_med^2;
    B=1/Yield_med;
    B=B*B;
    
    lambda_sum=-(conf_full*gamma0')';
    Q_right_transpose=cell(1,N);
    Q_left=cell(1,N);
    for ii=1:N
        %         w_right(:,:,ii)=sqrt(2)*[[1-p(ii);p(ii)],[-0.5;0.5]];
        %         w_left(:,:,ii)=sqrt(2)*[[0.5,0.5];[-p(ii),1-p(ii)]];
        Q_right_transpose{ii}=transpose([[1-p(ii);p(ii)],[-1;1]]);
        Q_left{ii}=[[1,1];[-p(ii),1-p(ii)]];
    end
    
%     Right_part=kronmult(Q_left,Yield_ps');
%     Left_part=kronmult(Q_right_transpose,Yield');
%     coef=(Left_part.*Right_part)';
    Right_part=kronm(fliplr(Q_left),Yield_ps');
    Left_part=kronm(fliplr(Q_right_transpose),Yield');
    coef=(Left_part.*Right_part)';
    
    %     V_right=w_right(:,:,1);
    %     V_left=w_left(:,:,1);
    %     for jj=2:N
    %         V_right=kron(V_right,w_right(:,:,jj));
    %         V_left=kron(V_left,w_left(:,:,jj));
    %     end
    %
    %     coef=(Yield*V_right).*(V_left*Yield_ps')';
    
    sp=x*0;
    for jj=1:numel(x)
        lambda_new=-2*lambda_sum./(lambda_sum.^2+4*pi^2*x(jj)^2);
        sp(jj) =2*B*sum(coef.*lambda_new);
    end