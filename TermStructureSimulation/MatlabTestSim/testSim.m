%% Estimate forward curves
load('fHist.mat')
load('piHist.mat')

%% Calculate eigenvector matrices
DZero = fHist(2:end,:) - fHist(1:end-1,:);
DTau = piHist(2:end,:) - piHist(1:end-1,:);

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, 6);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);

[V,D] = eigs(CTau, 6);
[e,ind] = sort(diag(D),1, 'descend');
ETau = V(:,ind);


plot(EZero(:,1:3))



%% Estimate simulation parameters
T = 1:730;
fZero = fHist(end-100,:);
pi = piHist(end,:);

X = randn(6,1000)*0.001;
deltaF = EZero * X;

fTest = fZero' + deltaF;
plot(T,fTest, T, fZero');

kappa = zeros(2,1);
xiHat = zeros(6,2);

d = 1;
fRes = [];



%% Simulate curves
tic

N = 2000;
[fZeroOut, fTauOut] = simMultipleYield(d, N, EZero, ETau, fZero, pi, kappa, xiHat);
tid = toc
%%
x = fZeroOut - fZero';

y = fTauOut - fZeroOut - pi';
%%
T = 1:730;

figure;
plot(T,fTauOut(1:end, 1), T, fZeroOut(1:end, 1), T, fZero)
figure;
plot(T,fTauOut', T, fZeroOut')

figure;
plot(T,fZero)




%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../lhsd.cpp ../unfGenT.cpp ../unfGenGauss.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:/boost_1_72_0 -IX:\exjobb\eigen-eigen-323c052e1731


