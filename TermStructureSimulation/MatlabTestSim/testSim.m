%% Estimate forward curves
load('fHist.mat')
load('piHist.mat')

%% Calculate eigenvector matrices
DZero = fHist(2:end,:) - fHist(1:end-1,:);
DTau = piHist(2:end,:) - piHist(1:end-1,:);

[Q1, R1] = qr(DZero', 0);
[U1, D1, V] = svd(R1', 0);
EZero = Q1 * V(:,1:6);

[Q1, R1] = qr(DTau', 0);
[U1, D1, V] = svd(R1', 0);
ETau = Q1 * V(:,1:6);


%% Estimate simulation parameters

fZero = fHist(end,:);
pi = piHist(end,:);

kappa = zeros(2,1);
xiHat = zeros(6,2);

d = 1;
fRes = [];



%% Simulate curves

[fZeroOut, fTauOut] = simMultipleYield(d, EZero, ETau, fZero, pi, kappa, xiHat);

%%
x = fZeroOut - fZero'

y = fTauOut - fZeroOut - pi'






