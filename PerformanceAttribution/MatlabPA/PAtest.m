%% Get data
load('Data/fHist.mat')
load('Data/piHist.mat')
load('Data/times.mat')
load('Data/yield.mat')
dates  = datevec(times);
fHistOOS = fAll(1777:2757,:);   %2012-04-02 - 2016-03-06
fHistIS = fAll(2758:end,:);    %2016-03-07 - 2018-12-11
piHistOOS = piAll(1777:2757,:); %2012-04-02 - 2016-03-06
piHistIS = piAll(2758:end,:);   %2016-03-07 - 2018-12-11

%% Calculate eigenvector matrices
DZero = fHistOOS(2:end,:) - fHistOOS(1:end-1,:);
DTau = piHistOOS(2:end,:) - piHistOOS(1:end-1,:);

kZero = 6;
kTau = 6;

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);

[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
ETau = V(:,ind);

E = {}
E.Zero = EZero;
E.Tau = ETau;
%plot(EZero(:,1:3))

%% Set parameters
n = size(fHistIS, 2);
A = intMatrix(n);
r = {};
rZero = A*fHistIS';
rTau = A*(fHistIS + piHistIS)';
r.rZero = rZero;
r.rTau = rTau;

deltaTj = 1.01388888888889;
floatCashFlows = [2; 92; 184; 275; 365]';
fixCashFlows = [365]';

r.rZero = rZero(:,1);
r.rTau = rTau(:,1);

fZero = fHistIS(1:2,:)';
fTau = (fHistIS(1:2,:) + piHistIS(1:2,:))';

%% Calculate gradient and hessian
N = 1000;
g = grad(N, y(1,2) / 100, floatCashFlows, fixCashFlows, A, E, deltaTj, r)
H = hes(N, y(1,2) / 100, floatCashFlows, fixCashFlows, A, E, deltaTj, r);

%%
dXi = [EZero' * (fZero(:,2) - fZero(:,1)); ETau' * (fTau(:,2) - fTau(:,1))]; 
dP = g'*dXi + (1/2)*dXi'*H*dXi;






%%






















%% Test price and yield
%P = irsPrice(N, y(1,2) / 100, floatCashFlows, fixCashFlows, deltaTj, r)
%yCalc = irsYield(floatCashFlows, fixCashFlows, deltaTj, r)
%% Numerical comparison

%P = irsPrice(N, y, floatCashFlows, fixCashFlows, deltaTj, r);
%fx = 
%fxph = 
%fxmh = 
%h = 10^(-6);

%gNum = firstOrder(fxph, fxmh, h);
%HNum = secondOrder(fx, fxph, fxmh, h);





