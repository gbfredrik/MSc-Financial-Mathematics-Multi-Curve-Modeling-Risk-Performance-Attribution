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
kPi = 9;

%E_k = {};
%EZero_k = EZero(:,1:kZero);
%ETau_k = ETau(:,1:kTau);
%E_k.Zero = EZero_k;
%E_k.Tau = ETau_k;

%% Init parameters
n = size(fHistIS, 2);
y = y(1, 2) / 100;
f = fHistIS';
pi = piHistIS';
N = 1000;
k = {};
k.zero = kZero;
k.Pi = kPi;
curveData = {};
curveData.zero = fHistIS;
curveData.pi = piHistIS;
times = times(2758:end);
startDate(1) = times(1);
endDate(1) = times(end);
floatCashFlows = {};
fixCashFlows = {};
floatCashFlows.oneY = [2; 94; 186; 277; 367]';
fixCashFlows.oneY  = [367]';

%%
n = size(fHistIS, 2);
A = intMatrix(n);
tic
r = A*fHistIS';
toc
%% Excecute performance attribution

% Input to mex-file:
%   0 - N: nominal amount                   - vector
%   1 - y: instrument yield                 - vector
%   2 - E: eigenvector matrices             - struct of matrices
%   3 - k: num. of. risk factors,           - vector of scalars
%   4 - floatCashFlows                      - struct of vectors
%   5 - fixCashFlows                        - struct of vectors
%   6 - curveData                           - matrix
%   7 - times, days from origin             - vector
%   8 - startDateIS
%   9 - endDateIS

tic
[NPV, carry, sumRiskFactors, epsI, epsA, epsP] = paMex(N, y, E, k, floatCashFlows, fixCashFlows, curveData, times, startDate, endDate);
toc
%mex ../paMex.cpp ../pa.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:/boost_1_72_0 -IX:\exjobb\eigen-eigen-323c052e1731
