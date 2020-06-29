function  [fTermSimulationsAllDays, fBaseSimulationsAllDays, fDemandSimulationsAllDays] = simulateFXCurves(fBase, fTerm, fDemand, N, inSample, outOfSample, curveLength)

% Simulate curves - All arguments

D_Term = fTerm(inSample(1)+1:inSample(end),1:curveLength) - fTerm(inSample(1):inSample(end)-1,1:curveLength);
D_Base = fBase(inSample(1)+1:inSample(end),1:curveLength) - fBase(inSample(1):inSample(end)-1,1:curveLength);
D_Demand = fDemand(inSample(1)+1:inSample(end),1:curveLength) - fDemand(inSample(1):inSample(end)-1,1:curveLength);

k_Term = 3;
k_Base = 3;
k_Demand = 6;

C_Term = cov(D_Term);
C_Base = cov(D_Base);
C_Demand = cov(D_Demand);

[V,D_dom] = eigs(C_Term, k_Term);
[e,ind] = sort(diag(D_dom),1, 'descend');
E = {};
E.Term = V(:,ind);

[V,D_for] = eigs(C_Base, k_Base);
[e,ind] = sort(diag(D_for),1, 'descend');
E.Base = V(:,ind);

[V,D_dem] = eigs(C_Demand, k_Demand);
[e,ind] = sort(diag(D_dem),1, 'descend');
E.Demand = V(:,ind);

% Parameters
[params_Term, rhoHat_Term, df_Copula_Term] = modelRiskFactorsT(D_Term, E.Term);
[params_Base, rhoHat_Base, df_Copula_Base] = modelRiskFactorsT(D_Base, E.Base);
[params_Demand, rhoHat_Demand, df_Copula_Demand] = modelRiskFactorsT(D_Demand, E.Demand);


marginalT = {};
marginalT.Term = 't';
marginalT.Base = 't';
marginalT.Demand = 't';

copulaT = {};
copulaT.Term = 't';
copulaT.Base = 't';
copulaT.Demand = 't';

varRedType = {};
varRedType.Term = 'none';
varRedType.Base = 'none';
varRedType.Demand = 'none';

muT = {};
muT.Term = params_Term(1,:);
muT.Base= params_Base(1,:);
muT.Demand = params_Demand(1,:);

dfC = {};
dfC.Term = [df_Copula_Term, df_Copula_Term, df_Copula_Term, df_Copula_Term, df_Copula_Term, df_Copula_Term];
dfC.Base= [df_Copula_Base, df_Copula_Base, df_Copula_Base, df_Copula_Base, df_Copula_Base, df_Copula_Base];
dfC.Demand = [df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand];
dfM = {};
dfM.Term = params_Term(5,:);
dfM.Base= params_Base(5,:);
dfM.Demand = params_Demand(5,:);

rhoT.Term = rhoHat_Term;
rhoT.Base = rhoHat_Base;
rhoT.Demand = rhoHat_Demand;

omegaT = {};
alphaT = {};
betaT = {};
omegaT.Term = params_Term(2,:);
omegaT.Base = params_Base(2,:);
omegaT.Demand = params_Demand(2,:);
alphaT.Term = params_Term(4,:);
alphaT.Base = params_Base(4,:);
alphaT.Demand = params_Demand(4,:);
betaT.Term = params_Term(3,:);
betaT.Base = params_Base(3,:);
betaT.Demand = params_Demand(3,:);

fRes = {};

fRes.Term = zeros(curveLength, N);
fRes.Base = zeros(curveLength, N);
fRes.Demand = zeros(curveLength, N);

gammaGauss = {};
gammaGauss.Term = zeros(k_Term, 1);
gammaGauss.Base = zeros(k_Base, 1);
gammaGauss.Demand = zeros(k_Demand, 1);

xiHat = {};
xiHat.Term = zeros(k_Term, 1);
xiHat.Base = zeros(k_Base, 1);
xiHat.Demand = zeros(k_Demand, 1);

kappa = zeros(2, 1);

% Execute simulation and plot results

%Argument list:
%   1 - E
%   2 - rho
%   3 - mu
%   4 - omega
%   5 - alpha
%   6 - beta
%   7 - hist
%   8 - marginal     (gaussian or t)
%   9 - copula       (gaussian or t)
%   10 - varRedType  (lhsd or none)
%   11 - d
%   12 - N
%   13 - fRes
%   ----------------- Optional values (not optional in the mex file...)
%   14 - gamma
%   15 - kappa
%   16 - xiHat
%   17 - dfC
%   18 - dfM

d = 1;

fTermSimulationsAllDays = cell(1,size(outOfSample,2)-1);
fBaseSimulationsAllDays = cell(1,size(outOfSample,2)-1);
fDemandSimulationsAllDays = cell(1,size(outOfSample,2)-1);

h = waitbar(0, "Simulating...");
for tradeDate = outOfSample(1):outOfSample(end-1)
    waitbar((tradeDate-inSample(end))/(size(outOfSample,2)-1), h)
    hist = {};
    hist.Term = fTerm(inSample(end)-10:tradeDate,1:curveLength);
    hist.Base = fBase(inSample(end)-10:tradeDate,1:curveLength);
    hist.Demand = fDemand(inSample(end)-10:tradeDate,1:curveLength);

    [fTermSimulationsAllDays{1, tradeDate-outOfSample(1)+1}, fBaseSimulationsAllDays{1, tradeDate-outOfSample(1)+1}, fDemandSimulationsAllDays{1, tradeDate-outOfSample(1)+1}] ... 
        = FXsimMultipleYield(E, rhoT, muT, omegaT, alphaT, betaT, hist, ...
                                marginalT, copulaT, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);

    
%     subplot(3,1,1);
%     plot(fTermOutT)
%     subplot(3,1,2);
%     plot(fBaseOut)
%     subplot(3,1,3)
%     plot(fDemandOutT)

end
close(h);


%mean(test,2)
clear FXsimMultipleYield.mexw64
%pause(0.01);
%end
%mex ../FXsimMultipleYield.cpp ../FXMultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -I../../../Dependencies/boost_1_71_0 -I../../../Dependencies/eigen-3.3.7

end