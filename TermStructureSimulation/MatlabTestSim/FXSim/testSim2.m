%% Estimate forward curves
load('10YrCurves.mat')
%load('fHist.mat')
%load('piHist.mat')



%% Simulate curves - All arguments
startDay = 1;
endDay = 3400;

DZero = fAll(2:end-1000,:) - fAll(1:end-1001,:);
DTau = piAll(2:end-1000,:) - piAll(1:end-1001,:);

kZero = 6;
kTau = 6;

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
E = {};
E.Zero = V(:,ind);

[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
E.Tau = V(:,ind);



%% Model risk factors

[params_Zero, rhoHat_Zero, df_Copula_Zero] = modelRiskFactorsT(fAll(2:end,:) - fAll(1:end-1,:), E.Zero);
[params_Tau, rhoHat_Tau, df_Copula_Tau] = modelRiskFactorsT(piAll(2:end,:) - piAll(1:end-1,:), E.Tau);

%%

d = 1;
N = 2000;

marginalGauss = {};
marginalGauss.Zero = 'normal';
marginalGauss.Tau = 'normal';

marginalT = {};
marginalT.Zero = 't';
marginalT.Tau = 't';

copulaGauss = {};
copulaGauss.Zero = 'normal';
copulaGauss.Tau = 'normal';

copulaT = {};
copulaT.Zero = 't';
copulaT.Tau = 't';

varRedType = {};
varRedType.Zero = 'none';
varRedType.Tau = 'none';

% muGauss = {};
% muGauss.Zero = params_Zero(1,:);
% muGauss.Tau = params_Tau(1,:);

muT = {};
muT.Zero = params_Zero(1,:)
muT.Tau = params_Tau(1,:)

dfC = {};
dfC.Zero =  [df_Copula_Zero,df_Copula_Zero,df_Copula_Zero,df_Copula_Zero,df_Copula_Zero,df_Copula_Zero];
dfC.Tau = [df_Copula_Tau,df_Copula_Tau,df_Copula_Tau,df_Copula_Tau,df_Copula_Tau,df_Copula_Tau];
dfM = {};
dfM.Zero =  params_Zero(5,:);
dfM.Tau =  params_Tau(5,:);

% rhoGauss = {};
% rhoT = {};
% rhoGauss.Zero = [1 0.00123861276154339 0.146108166779178; 0.00123861276154339 1 0.0595445213255097; ...
%     0.146108166779178 0.0595445213255097 1];
% rhoGauss.Tau = [1 0.137704279260191 -0.000136374806469469; 0.137704279260191 1 0.00883802054182843; ...
%     -0.000136374806469469 0.00883802054182843 1];
rhoT.Zero = rhoHat_Zero;
rhoT.Tau = rhoHat_Tau;

% omegaGauss = {};
% alphaGauss = {};
% betaGauss = {};
% omegaGauss.Zero = params_Zero(2,:);
% omegaGauss.Tau = parms_Tau(2,:);
% alphaGauss.Zero = [0.0706202403513367, 0.0496163281016998, 0.0340474887092155];
% alphaGauss.Tau = [0.153143545756948, 0.12216463832376, 0.11827068800463];
% betaGauss.Zero = [0.929378837402928, 0.950383257355491, 0.965859147308582];
% betaGauss.Tau = [0.846854848083354, 0.877828814698064, 0.881729172332618];

omegaT = {};
alphaT = {};
betaT = {};
omegaT.Zero = params_Zero(2,:);
omegaT.Tau = params_Tau(2,:);
alphaT.Zero = params_Zero(4,:)
alphaT.Tau = params_Tau(4,:);
betaT.Zero = params_Zero(3,:);
betaT.Tau = params_Tau(3,:);

hist = {};
hist.f = fAll(end-1000:end-300,:);
hist.pi = piAll(end-1000:end-300,:);

fRes = {};
fRes.Zero = zeros(size(fAll, 2), N);
fRes.Tau = zeros(size(fAll, 2), N);

gammaGauss = {};
gammaGauss.Zero = zeros(kZero, 1);
gammaGauss.Tau = zeros(kTau, 1);

xiHat = {};
xiHat.Zero = zeros(kZero, 1);
xiHat.Tau = zeros(kTau, 1);

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
%[fZeroOutGauss, fTauOutGauss, yee] = simMultipleYield(E, rhoGauss, muGauss, omegaGauss, alphaGauss, betaGauss, hist, ...
%    marginalGauss, copulaGauss, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
%tid = toc
%clear simMultipleYield.mexw64
%pause(0.01);
tic
[fZeroOutT, fTauOutT] = simMultipleYield(E, rhoT, muT, omegaT, alphaT, betaT, hist, ...
    marginalT, copulaT, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
tid = toc


T = 1:3650;

%subplot(1,2,1);
%plot(T, fZeroOutGauss', T, fTauOutGauss')
%subplot(1,2,2);
plot(T, fTauOutT')
%pause(0.01);


%figure;
%plot(T, fTauOut')
 
%mean(test,2)
clear simMultipleYield.mexw64
%pause(0.01);
%end
%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:/boost_1_72_0 -IX:\exjobb\eigen-eigen-323c052e1731


