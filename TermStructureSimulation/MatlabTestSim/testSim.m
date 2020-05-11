%% Estimate forward curves
load('fHist.mat')
load('piHist.mat')



%% Calculate eigenvector matrices
DZero = fHist(2:end,:) - fHist(1:end-1,:);
DTau = piHist(2:end,:) - piHist(1:end-1,:);

kZero = 3;
kTau = 3;

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);

[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
ETau = V(:,ind);


plot(EZero(:,1:3))



%% Estimate simulation parameters
T = 1:730;
fZero = fHist(end-50,:);
pi = piHist(end,:);


omega = [2.6*10^(-5), 4*10^(-5), 2.5*10^(-6)]';
alpha = [0.041, 0.24, 0.069]';
gamma = [0, 0, 0]';
beta = [0.94, 0.66, 0.93]';

%U = chol(cov(EZero));
c = EZero(:,1:end) -  mean(EZero(:,1:end));
%test = cov(EZero) - (c'*c)/(size(EZero,1) - 1);

U = chol((c'*c)/(size(EZero,1)-1));


sigma = GJR_GARCH(omega, alpha, gamma, beta, EZero, fHist);


X = randn(kZero,2000);
%U'*U = Sigma
%Z = U' * X;

Z = lhsnorm([0, 0, 0], cov(EZero), 2000);

deltaF = EZero * (Z' .* repmat(sigma,1,2000));


fTest = repmat(fZero',1,2000) + deltaF;
plot(T,fTest, T, fZero');

kappa = zeros(2,1);
xiHat = zeros(6,2);
d = 1;



%% Simulate curves - All arguments
startDay = 1;
endDay = 3400;

DZero = fHist(startDay + 1: endDay + 1,:) - fHist(startDay : endDay,:);
DTau = piHist(startDay + 1 : endDay + 1,:) - piHist(startDay : endDay,:);

fZero = fHist(end-50,:);
pi = piHist(end-50,:);

omega = [2.6*10^(-5), 4*10^(-5), 2.5*10^(-6)]';
alpha = [0.041, 0.24, 0.069]';
gamma = [0, 0, 0]';
beta = [0.94, 0.66, 0.93]';
kZero = 3;
kTau = 3;

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);

[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
ETau = V(:,ind);

kappa = zeros(2,1);
xiHat = zeros(6,2);
d = 1;
N = 2000;
marginal = ["normal", "normal"];
copula = ["normal", "normal"];
varRedType = ["lhsd", "lhsd"];
mu = [0 0; 0 0; 0 0];
df = [6.7 6.7; 5.1 5.1; 2.9 2.9];

%% Execute simulation and plot results

%Argument list:
%   d: number of days ahead to simulate
%   N: number of simulations
%   EZero: eigenvectors of zero diffs
%   ETau: eigenvectors of the tau diffs
%   fZero: the curve that the simulation is based on
%   pi: the tau that the simulation is based on
%   kappa: reversion speed
%   xiHat: long term average
%   omega, alpha, beta: garch parameters
%   fHist: term structure history, used in garch
%   piHist: tau history, used in garch
%   marginal: marginal distribution for each tenor
%   copula: copula used for each tenor
%   varRedType: chose type of variance reduction technique
%   mu: mean value for each risk factors
%   volatility of each risk factor
%   df: degrees of freedom for each risk factor

tic
[fZeroOut, fTauOut] = simMultipleYield(d, N, EZero, ETau, fZero, pi, kappa, ...
    xiHat, omega, alpha, beta, fHist, piHist, mu, df, 'normal', 'normal', 't', ...
    't', 'lhsd', 'lhsd');
tid = toc


T = 1:730;

figure;
plot(T,fTauOut', T, fZeroOut')





%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:/boost_1_72_0 -IX:\exjobb\eigen-eigen-323c052e1731


