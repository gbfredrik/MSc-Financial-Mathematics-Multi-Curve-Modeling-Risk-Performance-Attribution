%% Estimate forward curves
%load('10YrCurves.mat')
%load('fHist.mat')
%load('piHist.mat')

clearvars -except fForeign fDomestic fDemand fAll piAll FXSEKcurveRi100
demand = table2array(FXSEKcurveRi100);

%% Simulate N curves one day ahead
outOfSample = [1 80];
inSample = 81:100;
lastDiscPoint = 649;
curve = demand;git 
N = 2000;
fSimulationsAllDays = cell(1,size(inSample,2));

T = 1:lastDiscPoint;

D_Demand = curve(outOfSample(1)+1:outOfSample(end),1:lastDiscPoint) - ...
            curve(outOfSample(1):outOfSample(end)-1,1:lastDiscPoint);

mu = mean(D_Demand);
CDemand = cov(D_Demand);

% [U,e] = eig(CDemand, 'vector');
% e(e<0) = 0;
% [~,L] = qr(diag(sqrt(d))*U');

% L = chol(CDemand+1e-13*eye(size(CDemand)));
% r = normrnd(0,1,730,N);

deltaF = lhsnorm(mu', CDemand, N)';


for i = 1:size(inSample,2)
    fCurrent = curve(i,1:lastDiscPoint);
    fSimulations = (repmat(fCurrent',1,N) + deltaF)';
    fSimulationsAllDays{1,i} = fSimulations;
end



%deltaF = L*r;

%plot(T,fSimulationsAllDays{10});

%%



