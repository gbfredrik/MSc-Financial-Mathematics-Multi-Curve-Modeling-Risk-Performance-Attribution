%% Estimate forward curves
load('fHist.mat')
load('piHist.mat')

%% Simulate curves - All arguments
startDay = 1;
endDay = 3400;

DZero = fHist(2:end-1000,:) - fHist(1:end-1001,:);
DTau = piHist(2:end-1000,:) - piHist(1:end-1001,:);

kZero = 3;
kTau = 3;

CZero = cov(DZero);
CTau = cov(DTau);

[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
E = {};
E.Zero = V(:,ind);

[V,D] = eigs(CTau, kTau);
[e,ind] = sort(diag(D),1, 'descend');
E.Tau = V(:,ind);


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
varRedType.Zero = 'lhsd';
varRedType.Tau = 'lhsd';

muGauss = {};
muGauss.Zero = [-0.000204585228453844, -3.19909066455692*10^(-5), 1.00659791652654*10^(-5)];
muGauss.Tau = [-9.1097165867084*10^(-6), 7.11676370797266*10^(-6), -5.8568881762477*10^(-6)];

muT = {};
muT.Zero = [-0.000190206574936746, -7.36939259476214*10^(-5), -5.04450084368916*10^(-7)];
muT.Tau = [1.16800106617027*10^(-5),  1.01999993427177*10^(-5), -2.6233489649583*10^(-5)];

dfC = {};
dfC.Zero = [4.19891971204304, 4.19891971204304, 4.19891971204304];
dfC.Tau = [4.34478696553142, 4.34478696553142, 4.34478696553142];
dfM = {};
dfM.Zero = [5.28764683561276, 3.68715690722737, 4.35075491326916];
dfM.Tau = [3.91467935574674,  5.26185062542266,  3.25530037044725];

rhoGauss = {};
rhoT = {};
rhoGauss.Zero = [1 0.00123861276154339 0.146108166779178; 0.00123861276154339 1 0.0595445213255097; ...
    0.146108166779178 0.0595445213255097 1];
rhoGauss.Tau = [1 0.137704279260191 -0.000136374806469469; 0.137704279260191 1 0.00883802054182843; ...
    -0.000136374806469469 0.00883802054182843 1];
rhoT.Zero = [1 0.00746556043014991 0.18573702420161; 0.00746556043014991 1 0.0257979549771841; ...
    0.18573702420161 0.0257979549771841 1];
rhoT.Tau = [1 0.114548973189117 0.0221613642424653; 0.114548973189117 1 -0.045675862486691; ...
    0.0221613642424653 -0.045675862486691 1];

omegaGauss = {};
alphaGauss = {};
betaGauss = {};
omegaGauss.Zero = [2.19523356684698*10^(-7), 1.05143756909011*10^(-8), 1.6711102899898*10^(-9)];
omegaGauss.Tau = [8.23935679078583*10^(-7), 1.40669691906538*10^(-7), 5.64017944155194*10^(-8)];
alphaGauss.Zero = [0.0706202403513367, 0.0496163281016998, 0.0340474887092155];
alphaGauss.Tau = [0.153143545756948, 0.12216463832376, 0.11827068800463];
betaGauss.Zero = [0.929378837402928, 0.950383257355491, 0.965859147308582];
betaGauss.Tau = [0.846854848083354, 0.877828814698064, 0.881729172332618];

omegaT = {};
alphaT = {};
betaT = {};
omegaT.Zero = [4.93467318116551*10^(-8), 1.05694526687337*10^(-8), 4.73538895325478*10^-(16)];
omegaT.Tau = [2.75220171510865*10^(-7), 2.33577978973616*10^(-8), 2.52912309240021*10^(-8)];
alphaT.Zero = [0.0436131220019906, 0.0683945983845004, 0.0532838920915436];
alphaT.Tau = [0.126040521253363, 0.0785451521627585, 0.184681025423191];
betaT.Zero = [0.935425526459519, 0.88423290065623, 0.917777792980412];
betaT.Tau = [0.799200766350339, 0.887354806233473, 0.685203018906688];

hist = {};
hist.f = fHist(end-1000:end,:);
hist.pi = piHist(end-1000:end,:);

fRes = {};
fRes.Zero = zeros(size(fHist, 2), N);
fRes.Tau = zeros(size(fHist, 2), N);

gammaGauss = {};
gammaGauss.Zero = zeros(kZero, 1);
gammaGauss.Tau = zeros(kTau, 1);

xiHat = {};
xiHat.Zero = zeros(kZero, 1);
xiHat.Tau = zeros(kTau, 1);

kappa = zeros(2, 1);

%% Sim one day
T = 1:730;

fZero = fHist(end,:);
pi = piHist(end,:);
Z = zeros(kZero, N);

%c = EZero(:,1:end) -  mean(EZero(:,1:end));
%U = chol(rhoGauss(1));
%c = EZero(:,1:end) -  mean(EZero(:,1:end));
%test = cov(EZero) - (c'*c)/(size(EZero,1) - 1);
%U = chol((c'*c)/(size(EZero,1)-1));
sigma = GJR_GARCH(omegaGauss.Zero', alphaGauss.Zero', gammaGauss.Zero, betaGauss.Zero', E.Zero, fHist(end-1000:end,:));
%X = randn(kZero,N);
%U'*U = Sigma
%Z = U' * X;
Z = lhsnorm([0 0 0], rhoGauss.Zero, N)';
%X = mvnrnd([0 0 0], rho.Zero, N);
%U = normcdf(X);
%X = mvtrnd(rhoF, 4.19891971204304, N)';
%U = tcdf(X, 4.19891971204304);
%Z = norminv(U, 0, 1)';
%Z = tinv(U, 6);

deltaF = E.Zero * (muGauss.Zero' + Z.* repmat(sigma,1,N));
fZeroRes = (repmat(fZero',1,N) + deltaF)';

%U = chol(cov(ETau));
%c = EZero(:,1:end) -  mean(EZero(:,1:end));
%test = cov(EZero) - (c'*c)/(size(EZero,1) - 1);
%U = chol((c'*c)/(size(EZero,1)-1));
%sigma = GJR_GARCH(omega, alpha, gamma, beta, ETau, piHist(end-1000:end,:));
%X = randn(kTau,2000);
%U'*U = Sigma
%Z = U' * X;
%Z = lhsnorm(gamma, cov(ETau), 2000)';
%deltaTau = EZero * (Z .* repmat(sigma,1,2000));
%fTauRes = fZeroRes + repmat(pi',1,2000) + deltaTau;
%plot(T,fZeroRes*100, T, fTauRes*100);
figure;
plot(T,fZeroRes*100);
title(['Day: ', num2str(1)]);

%%
d = 10;
plot(T,fZeroRes*100);
title(['Day: ', num2str(1)]);
pause(0.1);
Z = zeros(kZero, 1);
fZeroPrev = fZeroRes; 
for j = 1:N

    sigma = GJR_GARCH_d(omega, alpha, gamma, beta, EZero, fZeroRes(j,:), fZero, sigma);
    %X = randn(kZero,1);
    %Z = U' * X;
    %Z = lhsnorm(mu, cov(EZero), 1)';
    X = mvnrnd([0 0 0], corr(EZero, 'type','Spearman'), 1)';
    U = normcdf(X);
    %X = mvtrnd(corr(EZero, 'type','Kendall'), 3, 1)';
    %U = tcdf(X, 3);
    Z = tinv(U, 3);
    deltaF = (EZero * (Z .* sigma))';
    fZeroRes(j,:) = fZeroRes(j,:) + deltaF;  
    
end
plot(T,fZeroRes*100);
title(['Day: ', num2str(2)]);
pause(0.1);

for i = 2:d
    fZeroPrev = fZeroRes; 
    for j = 1:N

        sigma = GJR_GARCH_d(omega, alpha, gamma, beta, EZero, fZeroRes(j,:), fZeroPrev(j,:), sigma);
        %X = randn(kZero,1);
        %Z = U' * X;
        %Z = lhsnorm(mu, cov(EZero), 1)';
        X = mvnrnd([0 0 0], corr(EZero, 'type','Spearman'), 1)';
        U = normcdf(X);
        %X = mvtrnd(corr(EZero, 'type','Kendall'), 3, 1)';
        %U = tcdf(X, 3);
        Z = tinv(U, 3);
        deltaF = (EZero * (Z .* sigma))';
        fZeroRes(j,:) = fZeroRes(j,:) + deltaF;  
        
    end
    plot(T,fZeroRes*100);
    title(['Day: ', num2str(i)]);
    pause(0.1);
end


