%% Get data
load('fHistOOS.mat')
load('piHistOOS.mat')

load('fHistIS.mat')
load('fHistIS.mat')

%% Calculate eigenvector matrices
DZero = fHist(2:end-1000,:) - fHist(1:end-1001,:);
DTau = piHist(2:end-1000,:) - piHist(1:end-1001,:);

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
%% Calculate IRS price

%P = irsPrice(N, y, )


%% Calculate gradient and hessian

