function PSDplot()
global freq S nu disp

S_up=S*0;
S_dn=S*0;
Sigma_gs=sqrt(disp');

for ii=1:numel(nu)
    if nu(ii)<100
        S_up(ii)=nu(ii).*S(ii)./chi2inv(0.025,nu(ii));
        S_dn(ii)=nu(ii).*S(ii)./chi2inv(0.975,nu(ii));
    else
        S_up(ii)=S(ii)+1.96*Sigma_gs(ii);
        S_dn(ii)=S(ii)-1.96*Sigma_gs(ii);
    end
end

errorbar(freq,S,S-S_dn,S_up-S,'s','Color','k','MarkerSize',8,'LineWidth',1.5);
set(gca,'YScale','log','XScale','log','XLim',[min(freq)/2 max(freq)*2],'YLim',[0.5*min(S_dn),2*max(S_up)])
hold on
grid on

end

