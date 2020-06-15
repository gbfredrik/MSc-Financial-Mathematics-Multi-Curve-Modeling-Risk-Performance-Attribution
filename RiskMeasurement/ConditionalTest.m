load ESBacktestBySimData
%%
VaRInd = 2;
figure();
plot(Dates,Returns,Dates,-VaR(:,VaRInd),Dates,-ES(:,VaRInd))
legend('Returns','VaR','ES')
title(['Test Data, ' num2str(VaRLevel(VaRInd)*100) '% Confidence'])
grid on
%%
rng('default'); % for reproducibility
IDs = ["t(dof) 95%","t(dof) 97.5%","t(dof) 99%"];
IDs = strrep(IDs,"dof",num2str(DoF));
ebts = esbacktestbysim(Returns,VaR,ES,Distribution,...
    'DegreesOfFreedom',DoF,...
    'Location',Mu,...
    'Scale',Sigma,...
    'PortfolioID',"S&P",...
    'VaRID',IDs,...
    'VaRLevel',VaRLevel);
disp(ebts)

disp(ebts.Distribution)

S = summary(ebts);
disp(S)
%%
t = runtests(ebts);
%%
[t,s] = conditional(ebts);


%