function  [fDomesticSimulationsAllDays, fForeignSimulationsAllDays, fDemandSimulationsAllDays] = simulateFXCurves(fBase, fTerm, fDemand, N, inSample, outOfSample, curveLength)

% fileNameBase = "FX_SEK_base_";
% fileNameTerm = "FX_SEK_term_";
% fileNameFx = "FX_SEK_fxAvg_";
% fileNameBase = "IRS_USD_ois";
% fileNameTerm = "IRS_USD_tenor";
% fileExt = "";
% 
% fileExt = "FX1-base1000-term100_cb10_jumpY"; % Change
% csvStr = ".csv";
% 
% fForeign = readmatrix(strcat(fileNameBase, fileExt, csvStr));
% fDomestic = readmatrix(strcat(fileNameTerm, fileExt, csvStr));
% fDemand = readmatrix(strcat(fileNameFx, fileExt, csvStr));



% Simulate curves - All arguments

D_Term = fTerm(inSample(1)+1:inSample(end),1:curveLength) - fTerm(inSample(1):inSample(end)-1,1:curveLength);
D_Base = fBase(inSample(1)+1:inSample(end),1:curveLength) - fBase(inSample(1):inSample(end)-1,1:curveLength);
D_Demand = fDemand(inSample(1)+1:inSample(end),1:curveLength) - fDemand(inSample(1):inSample(end)-1,1:curveLength);

% D_Domestic = fDomestic(2444:3845,1:730) - fDomestic(2443:3844,1:730);
% D_Foreign = fForeign(2444:3845,1:730) - fForeign(2443:3844,1:730);
% D_Demand = fDemand(2444:3845,1:730) - fDemand(2443:3844,1:730);

k_Domestic = 3;
k_Foreign = 3;
k_Demand = 6;

C_Domestic = cov(D_Term);
C_Foreign = cov(D_Base);
C_Demand = cov(D_Demand);

[V,D_dom] = eigs(C_Domestic, k_Domestic);
[e,ind] = sort(diag(D_dom),1, 'descend');
E = {};
E.Domestic = V(:,ind);

[V,D_for] = eigs(C_Foreign, k_Foreign);
[e,ind] = sort(diag(D_for),1, 'descend');
E.Foreign = V(:,ind);

[V,D_dem] = eigs(C_Demand, k_Demand);
[e,ind] = sort(diag(D_dem),1, 'descend');
E.Demand = V(:,ind);
%%
[params_Domestic, rhoHat_Domestic, df_Copula_Domestic] = modelRiskFactorsT(D_Term, E.Domestic);
[params_Foreign, rhoHat_Foreign, df_Copula_Foreign] = modelRiskFactorsT(D_Base, E.Foreign);
[params_Demand, rhoHat_Demand, df_Copula_Demand] = modelRiskFactorsT(D_Demand, E.Demand);


marginalT = {};
marginalT.Domestic = 't';
marginalT.Foreign = 't';
marginalT.Demand = 't';

copulaT = {};
copulaT.Domestic = 't';
copulaT.Foreign = 't';
copulaT.Demand = 't';

varRedType = {};
varRedType.Domestic = 'none';
varRedType.Foreign = 'none';
varRedType.Demand = 'none';

muT = {};
muT.Domestic = params_Domestic(1,:);
muT.Foreign= params_Foreign(1,:);
muT.Demand = params_Demand(1,:);

dfC = {};
dfC.Domestic = [df_Copula_Domestic, df_Copula_Domestic, df_Copula_Domestic, df_Copula_Domestic, df_Copula_Domestic, df_Copula_Domestic];
dfC.Foreign= [df_Copula_Foreign, df_Copula_Foreign, df_Copula_Foreign, df_Copula_Foreign, df_Copula_Foreign, df_Copula_Foreign];
dfC.Demand = [df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand, df_Copula_Demand];
dfM = {};
dfM.Domestic = params_Domestic(5,:);
dfM.Foreign= params_Foreign(5,:);
dfM.Demand = params_Demand(5,:);

rhoT.Domestic = rhoHat_Domestic;
rhoT.Foreign = rhoHat_Foreign;
rhoT.Demand = rhoHat_Demand;

omegaT = {};
alphaT = {};
betaT = {};
omegaT.Domestic = params_Domestic(2,:);
omegaT.Foreign = params_Foreign(2,:);
omegaT.Demand = params_Demand(2,:);
alphaT.Domestic = params_Domestic(4,:);
alphaT.Foreign = params_Foreign(4,:);
alphaT.Demand = params_Demand(4,:);
betaT.Domestic = params_Domestic(3,:);
betaT.Foreign = params_Foreign(3,:);
betaT.Demand = params_Demand(3,:);

%hist.Domestic = fAll(end-1000:end-300,:);
%hist.Foreign = fAll(end-1000:end-300,:);
% hist.Demand = piAll(end-1000:end-300,:);

%hist = {};


%hist.Domestic = fDomestic(1001:1356,1:730);
% hist.Domestic = fDomestic(1001:1356,1:730);
% hist.Foreign = fForeign(1001:1356,1:730);
% hist.Demand = fDemand(1001:1356,1:730);
% hist.Domestic = fDomestic(1001:1385,1:730);
% hist.Foreign = fForeign(1001:1385,1:730);
% hist.Demand = fDemand(1001:1385,1:730);
% hist.Domestic = fDomestic(endDay:endSample,1:730);
% hist.Foreign = fForeign(endDay:endSample,1:730);
%hist.Demand = fDemand(endDay:endSample,1:730);

fRes = {};
%fRes.Domestic = zeros(size(fAll, 2), N);
% fRes.Foreign = zeros(size(fAll, 2), N);
% fRes.Demand = zeros(size(fAll, 2), N);

%fRes = {};
fRes.Domestic = zeros(curveLength, N);
fRes.Foreign = zeros(curveLength, N);
fRes.Demand = zeros(curveLength, N);

gammaGauss = {};
gammaGauss.Domestic = zeros(k_Domestic, 1);
gammaGauss.Foreign = zeros(k_Foreign, 1);
gammaGauss.Demand = zeros(k_Demand, 1);

xiHat = {};
xiHat.Domestic = zeros(k_Domestic, 1);
xiHat.Foreign = zeros(k_Foreign, 1);
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

%for i = 1:100
%  
%tic
%[fDomesticOutGauss, fDemandOutGauss, yee] = simMultipleYield(E, rhoGauss, muGauss, omegaGauss, alphaGauss, betaGauss, hist, ...
%    marginalGauss, copulaGauss, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
%tid = toc
%clear simMultipleYield.mexw64
%pause(0.01);
d = 1;

fDomesticSimulationsAllDays = cell(1,size(outOfSample,2));
fForeignSimulationsAllDays = cell(1,size(outOfSample,2));
fDemandSimulationsAllDays = cell(1,size(outOfSample,2));

h = waitbar(0, "Simulating...");
for tradeDate = outOfSample(1):outOfSample(end-1)
    waitbar((tradeDate-inSample(end))/(size(outOfSample,2)-1), h)
    hist = {};
    hist.Domestic = fTerm(inSample(end):tradeDate,1:curveLength);
    hist.Foreign = fBase(inSample(end):tradeDate,1:curveLength);
    hist.Demand = fDemand(inSample(end):tradeDate,1:curveLength);

    [fDomesticOutT, fForeignOut, fDemandOutT] = FXsimMultipleYield(E, rhoT, muT, omegaT, alphaT, betaT, hist, ...
        marginalT, copulaT, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);

    fDomesticSimulationsAllDays{1, tradeDate-outOfSample(1)+1} = fDomesticOutT;
    fForeignSimulationsAllDays{1, tradeDate-outOfSample(1)+1} = fForeignOut;
    fDemandSimulationsAllDays{1, tradeDate-outOfSample(1)+1} = fDemandOutT;
    
    subplot(3,1,1);
    plot(fDomesticOutT)
    subplot(3,1,2);
    plot(fForeignOut)
    subplot(3,1,3)
    plot(fDemandOutT)

end





%mean(test,2)
clear FXsimMultipleYield.mexw64
%pause(0.01);
%end
%mex ../FXsimMultipleYield.cpp ../FXMultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -I../../../Dependencies/boost_1_71_0 -I../../../Dependencies/eigen-3.3.7

end