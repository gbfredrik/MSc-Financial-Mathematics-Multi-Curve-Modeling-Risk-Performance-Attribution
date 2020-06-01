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


%%
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
%% Init parameters and calculate first day pricing error term
floatCashFlows = [2; 94; 186; 277; 367]';
fixCashFlows = [367]';
deltaTj = fixCashFlows / 360;
dXi = [EZero_k' * (f(:,1) - fOOS(:,end)); ETau_k' * (pi(:,1) - piOOS(:,end))]; 
XiBar = [EZero' * f(:,1); ETau' * pi(:,1)]; 

epsP = zeros(367, 1);
epsI = zeros(367, 1);
epsA = zeros(367, 1);
shiftZero = zeros(367, 1);
twistZero = zeros(367, 1);
butterflyZero = zeros(367, 1);
fourthZero = zeros(367, 1);
fifthZero = zeros(367, 1);
sixthZero = zeros(367, 1);
shift2Zero = zeros(367, 1);
twist2Zero = zeros(367, 1);
butterfly2Zero = zeros(367, 1);
fourth2Zero = zeros(367, 1);
fifth2Zero = zeros(367, 1);
sixth2Zero = zeros(367, 1);
shiftTau= zeros(367, 1);
twistTau = zeros(367, 1);
butterflyTau = zeros(367, 1);
fourthTau = zeros(367, 1);
fifthTau = zeros(367, 1);
sixthTau = zeros(367, 1);
shift2Tau = zeros(367, 1);
twist2Tau = zeros(367, 1);
butterfly2Tau = zeros(367, 1);
fourth2Tau = zeros(367, 1);
fifth2Tau = zeros(367, 1);
sixth2Tau = zeros(367, 1);
approx = zeros(367, 1);
carry = zeros(367, 1);
NPV = zeros(367, 1);

%------------------- PRICE ERR ---------------------------
epsP(1, 1) = irsPrice(N, y, floatCashFlows, fixCashFlows, deltaTj, r(:,1), piSpot(:,1));
NPV(1, 1) = epsP(1, 1);

%---------------------------------------------------------
%% Performance attribution for 1 year IRS

    postFirstDate = 0;
% Skapa en till for-loop för att hantera flera kontrakt. Initiera den inre
% loopen med kontraktspecifik data från den yttre loopen.
    
for i = 2:252

    % Hantera kassaflödesdatum
    datediff = times(2757 + i) - times(2756 + i)
    floatCashFlowsPrev = floatCashFlows;
    fixCashFlowsPrev = fixCashFlows;
    deltaTjPrev = deltaTj;
    if (postFirstDate == 0)
        floatCashFlows(1:end) = floatCashFlows(1:end) - datediff;
        fixCashFlows(1:end) = fixCashFlows(1:end) - datediff;
        if floatCashFlows(1) < 0 
            floatCashFlows = [0, floatCashFlows(2:end)];
            postFirstDate = 1;
        end
    elseif (postFirstDate == 1)
        floatCashFlows(2:end) = floatCashFlows(2:end) - datediff;
        fixCashFlows(1:end) = fixCashFlows(1:end) - datediff;
        if floatCashFlows(2) == 0 
            floatCashFlows = [0, floatCashFlows(3:end)];
        end
    end

    % Sätt övriga variabler
    deltaTj = fixCashFlows / 360;
    dXiPrev = dXi;
    XiBarPrev = XiBar;
    dXi = [EZero_k' * (f(:,i) - f(:,i-1)); ETau_k' * (pi(:,i) - pi(:,i-1))]; 
    XiBar = [EZero' * f(:,i); ETau' * pi(:,i)];
    
    gPrev = grad(N, y, floatCashFlows, fixCashFlows, A, E_k, deltaTj, r(:,i-1), piSpot(:,i-1));
    HPrev = hes(N, y, floatCashFlows, fixCashFlows, A, E_k, deltaTj, r(:,i-1), piSpot(:,i-1));
    
 
    %------------------- CARRY ----------------------------------------
    carry(i, 1) = irsPriceRiskFactor(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, XiBarPrev) ...
           - irsPriceRiskFactor(N, y, floatCashFlowsPrev, fixCashFlowsPrev, A, E, deltaTjPrev, XiBarPrev);      
    %------------------- TRUNC. ERR. ---------------------------------- 
    epsI(i, 1) = irsPrice(N, y, floatCashFlows, fixCashFlows, deltaTj, r(:,i), piSpot(:,i))...  
        - irsPriceRiskFactor(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, XiBarPrev + [dXi(1:kZero); zeros((length(XiBarPrev)/2 - kZero), 1); dXi(kZero+1:end); zeros((length(XiBarPrev)/2 - kTau), 1)]);
    %------------------- TAYL. APPROX. ERR ---------------------------- 
    epsA(i, 1) = irsPriceRiskFactor(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, XiBarPrev + [dXi(1:kZero); zeros((length(XiBarPrev)/2 - kZero), 1); dXi(kZero+1:end); zeros((length(XiBarPrev)/2 - kTau), 1)]) ...
           - irsPriceRiskFactor(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, XiBarPrev) ...
           - carry(i, 1) - gPrev'*dXi - (1/2)*dXi'*HPrev*dXi;
    %------------------- Shift, twist, butterfly... -------------------    
    dXiMat = diag(dXi);

    shiftZero(i, 1) = gPrev(1)*dXi(1);
    twistZero(i, 1) = gPrev(2)*dXi(2);
    butterflyZero(i, 1) = gPrev(3)*dXi(3);
    fourthZero(i, 1) = gPrev(4)*dXi(4);
    fifthZero(i, 1) = gPrev(5)*dXi(5);
    sixthZero(i, 1) = gPrev(6)*dXi(6);

    shift2Zero(i, 1) = (1/2) * dXiMat(:,1)' * HPrev * dXiMat(:,1);
    twist2Zero(i, 1) = (1/2) * dXiMat(:,2)' * HPrev * dXiMat(:,2);
    butterfly2Zero(i, 1) = (1/2) * dXiMat(:,3)' * HPrev * dXiMat(:,3);
    fourth2Zero(i, 1) = (1/2) * dXiMat(:,4)' * HPrev * dXiMat(:,4);
    fifth2Zero(i, 1) = (1/2) * dXiMat(:,5)' * HPrev * dXiMat(:,5);
    sixth2Zero(i, 1) = (1/2) * dXiMat(:,6)' * HPrev * dXiMat(:,6);

    shiftTau(i, 1) = gPrev(7)*dXi(7);
    twistTau(i, 1) = gPrev(8)*dXi(8);
    butterflyTau(i, 1) = gPrev(9)*dXi(9);
    fourthTau(i, 1) = gPrev(10)*dXi(10);
    fifthTau(i, 1) = gPrev(11)*dXi(11);
    sixthTau(i, 1) = gPrev(12)*dXi(12);

    shift2Tau(i, 1) = (1/2) * dXiMat(:,7)' * HPrev * dXiMat(:,7);
    twist2Tau(i, 1) = (1/2) * dXiMat(:,8)' * HPrev * dXiMat(:,8);
    butterfly2Tau(i, 1) = (1/2) * dXiMat(:,9)' * HPrev * dXiMat(:,9);
    fourth2Tau(i, 1) = (1/2) * dXiMat(:,10)' * HPrev * dXiMat(:,10);
    fifth2Tau(i, 1) = (1/2) * dXiMat(:,11)' * HPrev * dXiMat(:,11);
    sixth2Tau(i, 1) = (1/2) * dXiMat(:,12)' * HPrev * dXiMat(:,12);
    
    approx(i, 1) = gPrev'*dXi + (1/2)*dXi'*HPrev*dXi;
    
    %------------------- NPV ------------------------------------------    
    NPV(i, 1) = irsPrice(N, y, floatCashFlows, fixCashFlows, deltaTj, r(:,i), piSpot(:,i)) ...
        - irsPrice(N, y, floatCashFlowsPrev, fixCashFlowsPrev, deltaTjPrev, r(:,i-1), piSpot(:,i-1));
         
    i
end


%% Plot
T = 1:367;

sumRiskFactorsFirstZero = shiftZero + twistZero + butterflyZero + fourthZero + fifthZero + sixthZero;
sumRiskFactorsSecondZero = shift2Zero + twist2Zero + butterfly2Zero + fourth2Zero + fifth2Zero + sixth2Zero;
sumRiskFactorsFirstTau = shiftTau + twistTau + butterflyTau + fourthTau + fifthTau + sixthTau;
sumRiskFactorsSecondTau = shift2Tau + twist2Tau + butterfly2Tau + fourth2Tau + fifth2Tau + sixth2Tau;
sumRiskFactorsZero = sumRiskFactorsFirstZero + sumRiskFactorsSecondZero;
sumRiskFactorsTau = sumRiskFactorsFirstTau + sumRiskFactorsSecondTau;
SumRiskFactors = sumRiskFactorsZero + sumRiskFactorsTau;

cNPV = cumsum(NPV);
cEpsA = cumsum(epsA);
cEpsI = cumsum(epsI);
cCarry = cumsum(carry);
cApprox = cumsum(approx);
cSumRiskFactors = cumsum(SumRiskFactors);
%%
plot(T(1:252), cNPV(1:252), T(1:252), cApprox(1:252), T(1:252), cCarry(1:252), T(1:252), cEpsA(1:252), T(1:252), cEpsI(1:252), T(1:252), epsP(1:252));

%areaMatrix = [cNPV(1:252), cSumRiskFactors(1:252), cCarry(1:252), cEpsA(1:252), cEpsI(1:252), epsP(1:252)];
%area(T(1:252), areaMatrix)

legend('NPV', 'SumRiskFactors', 'carry', 'epsA', 'epsI', 'epsP');


%%
figure;
plot(T, cSumRiskFactors, T, carry, T, epsI, T, epsA);
legend('SumRiskFactors', 'carry', 'epsI', 'epsA');
figure;
plot(T, epsI, T, epsA);
legend('epsI', 'epsA');
figure;
plot(T, shiftZero, T, twistZero, T, butterflyZero, T, fourthZero, T, fifthZero, T, sixthZero, T, shift2Zero);
legend('shiftZero', 'twistZero', 'butterflyZero', 'fourthZero', 'fifthZero', 'sixthZero', 'shift2Zero');

%% Validering
validation = NPV - approx - epsP - epsI - carry - epsA - carry;

sum(abs(validation))

%cPerformance = cumsum(SumRiskFactors) + cumsum(carry);

%figure;
%plot(T, NPV, T, cPerformance);
%legend('NPV', 'RiskFactors + carry');




