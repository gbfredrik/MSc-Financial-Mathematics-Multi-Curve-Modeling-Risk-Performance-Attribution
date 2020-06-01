%% Get data
load('Data/fHist.mat')
load('Data/piHist.mat')
load('Data/times.mat')
load('Data/yield.mat')
fHistOOS = fAll(1:2757,:);   %2012-04-02          - 2016-03-03 (friday)
fHistIS = fAll(2758:end,:);     %2016-03-07 (monday) - 2018-12-11
piHistOOS = piAll(1:2757,:); %2012-04-02          - 2016-03-03 (friday)
piHistIS = piAll(2758:end,:);   %2016-03-07 (monday) - 2018-12-11
%% Calculate eigenvector matrices
DZero = fHistOOS(2:end,:) - fHistOOS(1:end-1,:);
DTau =  piHistOOS(2:end,:) - piHistOOS(1:end-1,:);

CZero = cov(DZero);
CTau = cov(DTau);
kZero = size(CZero, 1);
kTau = size(CZero, 1);
[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);
[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
ETau = V(:,ind);
E = {};
E.Zero = EZero;
E.Tau = ETau;

kZero = 9;
kTau = 9;

E_k = {};
EZero_k = EZero(:,1:kZero);
ETau_k = ETau(:,1:kTau);
E_k.Zero = EZero_k;
E_k.Tau = ETau_k;

%% Init parameters
n = size(fHistIS, 2);
A = intMatrix(n);
r = A*fHistIS';
piSpot = A * piHistIS';
y = y(1, 2) / 100;
f = fHistIS';
pi = piHistIS';
fOOS = fHistOOS';
piOOS = piHistOOS';
N = 1000;

%% Excecute performance attribution

% Input to mex-file:
%   0 - N: nominal amount                   - vector
%   1 - y: instrument yield                 - vector
%   2 - E: eigenvector matrices             - struct of matrices
%   3 - k: num. of. risk factors,           - vector of scalars        - scalar
%   4 - floatCashFlows                      - struct of vectors
%   5 - fixCashFlows                        - struct of vectors
%   6 - curveData                           - matrix
%   7 - times, days from origin             - vector
%   8 - startDate
%   9 - endDate

[NPV, carry, sumRiskFactors, epsI, epsA, epsP] = paMex(N, y, E, kZero, kTau, floatCashFlows, fixCashFlows, f, pi, times, startDate, endDate);





