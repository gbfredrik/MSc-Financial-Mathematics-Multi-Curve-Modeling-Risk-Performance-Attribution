clear
close all
%%

fileNameBase = "FX_SEK_base_";
fileNameTerm = "FX_SEK_term_";
fileNameFx = "FX_SEK_fxAvg_";
% fileNameBase = "IRS_USD_ois";
% fileNameTerm = "IRS_USD_tenor";
% fileExt = "";

fileExt = "FX1-base1000-term100_cb10_jumpY"; % Change
csvStr = ".csv";

fForeign = readmatrix(strcat(fileNameBase, fileExt, csvStr));
fDomestic = readmatrix(strcat(fileNameTerm, fileExt, csvStr));
fDemand = readmatrix(strcat(fileNameFx, fileExt, csvStr));



%% Simulate curves - All arguments
startDay = 1;
endDay = 1500;
endSample = 1783;

D_Domestic = fDomestic(2:endDay,1:730) - fDomestic(1:endDay-1,1:730);
D_Foreign = fForeign(2:endDay,1:730) - fForeign(1:endDay-1,1:730);
D_Demand = fDemand(2:endDay,1:730) - fDemand(1:endDay-1,1:730);

% D_Domestic = fDomestic(2444:3845,1:730) - fDomestic(2443:3844,1:730);
% D_Foreign = fForeign(2444:3845,1:730) - fForeign(2443:3844,1:730);
% D_Demand = fDemand(2444:3845,1:730) - fDemand(2443:3844,1:730);

k_Domestic = 3;
k_Foreign = 3;
k_Demand = 6;

C_Domestic = cov(D_Domestic);
C_Foreign = cov(D_Foreign);
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
[params_Domestic, rhoHat_Domestic, df_Copula_Domestic] = modelRiskFactorsT(D_Domestic, E.Domestic);
[params_Foreign, rhoHat_Foreign, df_Copula_Foreign] = modelRiskFactorsT(D_Foreign, E.Foreign);
[params_Demand, rhoHat_Demand, df_Copula_Demand] = modelRiskFactorsT(D_Demand, E.Demand);

d = 1;
N = 2000;

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

hist = {};
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
hist.Domestic = fDomestic(endDay:endSample,1:730);
hist.Foreign = fForeign(endDay:endSample,1:730);
hist.Demand = fDemand(endDay:endSample,1:730);

fRes = {};
%fRes.Domestic = zeros(size(fAll, 2), N);
% fRes.Foreign = zeros(size(fAll, 2), N);
% fRes.Demand = zeros(size(fAll, 2), N);

%fRes = {};
fRes.Domestic = zeros(730, N);
fRes.Foreign = zeros(730, N);
fRes.Demand = zeros(730, N);

gammaGauss = {};
gammaGauss.Domestic = zeros(k_Domestic, 1);
gammaGauss.Foreign = zeros(k_Foreign, 1);
gammaGauss.Demand = zeros(k_Demand, 1);

xiHat = {};
xiHat.Domestic = zeros(k_Domestic, 1);
xiHat.Foreign = zeros(k_Foreign, 1);
xiHat.Demand = zeros(k_Demand, 1);

kappa = zeros(2, 1);
%% Execute simulation and plot results

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
tic
[fDomesticOutT, fForeignOut, fDemandOutT] = FXsimMultipleYield(E, rhoT, muT, omegaT, alphaT, betaT, hist, ...
    marginalT, copulaT, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
tid = toc


%T = 1:3650;
T = 1:730;

%subplot(1,2,1);
%plot(T, fDomesticOutGauss', T, fDemandOutGauss')
%subplot(1,2,2);
%plot(T, fDomesticOutT')
%plot(T, fForeignOut')
% plot(T, fDomesticOutT')
%pause(0.01);


subplot(3,1,1)
plot(T, fDomesticOutT')
title('Quote OIS')

subplot(3,1,2)
plot(T, fForeignOut')
title('Base OIS')

subplot(3,1,3)
plot(T, fDemandOutT')
title('Demand curve')
sgtitle(fileExt)

figure;
subplot(2,1,1)
plot(var(fDomesticOutT'));
title('Domestic curve')
subplot(2,1,2)
plot(var(fForeignOut'));
title('Foreign')
sgtitle(fileExt)
 
%mean(test,2)
clear FXsimMultipleYield.mexw64
%pause(0.01);
%end
%mex ../FXsimMultipleYield.cpp ../FXMultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -I../../../Dependencies/boost_1_71_0 -I../../../Dependencies/eigen-3.3.7

