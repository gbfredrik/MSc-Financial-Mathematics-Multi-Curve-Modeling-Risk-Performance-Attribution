%% Estimate forward curves
load('10YrCurves.mat')
%load('fHist.mat')
%load('piHist.mat')

%% Simulate curves - All arguments
startDay = 1;
endDay = 3400;

DZero = fAll(2:end-1000,:) - fAll(1:end-1001,:);
DTau = piAll(2:end-1000,:) - piAll(1:end-1001,:);

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
varRedType.Zero = 'none';
varRedType.Tau = 'none';

muGauss = {};
muGauss.Zero = [-0.000204585228453844, -3.19909066455692*10^(-5), 1.00659791652654*10^(-5)];
muGauss.Tau = [-9.1097165867084*10^(-6), 7.11676370797266*10^(-6), -5.8568881762477*10^(-6)];

muT = {};
muT.Zero = [-0.000913425664668150,-0.000263398104793440,0.000247219225331290,0.000227770223535396,-5.72460180285886*10^(-5),8.06784537339761*10^(-5)];
muT.Tau = [0.000224532985824916,-1.32238887996754*10^(-5),4.67129307165196*10^(-5),-3.44411375152431*10^(-5),-2.55297225541061*10^(-5),1.91624593379172*10^(-5)];

dfC = {};
dfC.Zero = [4.009902956194651, 4.009902956194651, 4.009902956194651, 4.009902956194651, 4.009902956194651, 4.009902956194651];
dfC.Tau = [4.538249855317766, 4.538249855317766, 4.538249855317766, 4.538249855317766, 4.538249855317766, 4.538249855317766];
dfM = {};
dfM.Zero = [5.31609159530994,5.59774459787255,4.03382812087348,5.28258245929443,4.87376341136777,4.58739791424502];
dfM.Tau = [3.96922820601562,5.03678822659782,4.38691286041102,5.22302160698119,4.87285061419161,5.15846764485123];

rhoGauss = {};
rhoT = {};
rhoGauss.Zero = [1 0.00123861276154339 0.146108166779178; 0.00123861276154339 1 0.0595445213255097; ...
    0.146108166779178 0.0595445213255097 1];
rhoGauss.Tau = [1 0.137704279260191 -0.000136374806469469; 0.137704279260191 1 0.00883802054182843; ...
    -0.000136374806469469 0.00883802054182843 1];
rhoT.Zero = [1,-0.389404245331893,0.206167401437967,0.108250511006425,-0.0781904646288599,-0.0707307890908564;
    -0.389404245331893,1,-0.444471508814934,-0.0251155999672263,0.0110085002052923,-0.210493581704619;
    0.206167401437967,-0.444471508814934,1,0.0873380037802093,0.117917490275123,0.193024859550071;
    0.108250511006425,-0.0251155999672263,0.0873380037802093,1,0.0349108950102930,-0.180633890368627;
    -0.0781904646288599,0.0110085002052923,0.117917490275123,0.0349108950102930,1,0.000890674590012762;
    -0.0707307890908564,-0.210493581704619,0.193024859550071,-0.180633890368627,0.000890674590012762,1];
rhoT.Tau = [1,0.00198437626602062,0.166021989607804,-0.0284519082538309,0.138516932276534,0.142909750222074;
    0.00198437626602062,1,0.157206117187295,0.131261834609791,0.0992638373871272,0.00539995560973405;
    0.166021989607804,0.157206117187295,1,0.0256432496103997,0.144482310166173,0.132591672218439;
    -0.0284519082538309,0.131261834609791,0.0256432496103997,1,0.141724064335315,-0.223951910824899;
    0.138516932276534,0.0992638373871272,0.144482310166173,0.141724064335315,1,0.0317193599909348;
    0.142909750222074,0.00539995560973405,0.132591672218439,-0.223951910824899,0.0317193599909348,1];

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
omegaT.Zero = [1.95871766869112*10^(-5),1.50663861908391*10^(-6),1.04169496563947*10^(-6),1.27903019674284*10^(-7),6.65635699439975*10^(-8),2.66215571427713*10^(-8)];
omegaT.Tau = [3.63844318642688*10^(-6),3.12223741678781*10^(-6),1.11388691798752*10^(-5),5.45268740367561*10^(-7),9.18507414853582*10^(-7),4.70724330327256*10^(-8)];
alphaT.Zero = [0.104045217423394,0.0498963341745904,0.0713723919341992,0.0336351389809415,0.0343301092453374,0.0581964257227426];
alphaT.Tau = [0.252637502820338,0.181396309410320,0.314955115122157,0.0650860273837165,0.103175450695571,0.0808177724021258];
betaT.Zero = [0.843980311870861,0.920596231112623,0.856942944198292,0.943883288654432,0.938911438424578,0.903226936583101];
betaT.Tau = [0.648908098567852,0.732681595799868,0.317355768088827,0.878601057895940,0.774085728040299,0.877225162278746];

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
%   8 - marginal     (normal or t)
%   9 - copula       (normal or t)
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

for i = 1:10
%  
tic
[fZeroOutGauss, fTauOutGauss] = simMultipleYield(E, rhoGauss, muGauss, omegaGauss, alphaGauss, betaGauss, hist, ...
    marginalGauss, copulaGauss, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
%tid = toc
%clear simMultipleYield.mexw64
%pause(0.01);
%tic
%[fZeroOutT, fTauOutT] = simMultipleYield(E, rhoT, muT, omegaT, alphaT, betaT, hist, ...
%    marginalT, copulaT, varRedType, d, N, fRes, gammaGauss, kappa, xiHat, dfC, dfM);
tid = toc


T = 1:3650;

figure(1)
%subplot(1,2,1);
%plot(T, fZeroOutGauss', T, fTauOutGauss')
%subplot(1,2,2);
plot(T, fZeroOutGauss')
%pause(0.01);


%figure;
%plot(T, fTauOut')
 
%mean(test,2)
pause(0.01);
end
clear simMultipleYield.mexw64

%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:/boost_1_72_0 -IX:\exjobb\eigen-eigen-323c052e1731

%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -IX:\Documents\MSc Git\MScCurveModeling\Dependencies\boost_1_72_0 -IX:\Documents\MSc Git\MScCurveModeling\Dependencies\eigen-3.3.7
%mex ../simMultipleYield.cpp ../MultipleYieldSim.cpp ../varRed.cpp ../unfGen.cpp ../../Mathlibrary/statisticsOperations.cpp ../../Mathlibrary/matrixOperations.cpp ../../Mathlibrary/rvSim.cpp -I..\..\Dependencies\boost_1_72_0 -I..\..\Dependencies\eigen-3.3.7

