%TODO-list:
%---------------- Kontrollera nuvärdet mot reuters swap-prissättning
%1. Få det att fungera oavsett pay eller receive
%2. Kontrollera dom små hoppen vid rörliga utbetalningar
%3. lägg till accrued interest för fixa benet, kontrollera npv mot SWAPR i
%reuters
%Felet nu är att det fixa benet 


%% Get data
load('Data/fHist.mat')
load('Data/piHist.mat')
load('Data/times.mat')
load('Data/yield.mat')
load('Data/libor.mat')
fHistOOS = fAll(1:2757,:);   %2012-04-02          - 2016-03-03 (friday)
fHistIS = fAll(2758:end,:);     %2016-03-07 (monday) - 2018-12-11
piHistOOS = piAll(1:2757,:); %2012-04-02          - 2016-03-03 (friday)
piHistIS = piAll(2758:end,:);   %2016-03-07 (monday) - 2018-12-11
libor = libor(2826:end,:) / 100;

%% Calculate eigenvector matrices
DZero = fHistOOS(2:end,:) - fHistOOS(1:end-1,:);
DPi =  piHistOOS(2:end,:) - piHistOOS(1:end-1,:);
CZero = cov(DZero);
CPi = cov(DPi);
kZero = size(CZero, 1);
kPi = size(CZero, 1);
[V,D] = eigs(CZero, kZero);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);
[V,D] = eigs(CPi, kPi);
[e,ind] = sort(diag(D),1, 'descend');
EPi = V(:,ind);
kZero = 8;
kPi = 8;
EZero_k = EZero(:,1:kZero);
EPi_k = EPi(:,1:kPi);

%%
n = size(fHistIS, 2);
A = intMatrix(n);
f = fHistIS';
pi = piHistIS';
r = A * f;
piSpot = A * pi;
times = times(2758:end);

%% Cash flows of instruments
floatCashFlows_1Y = [94; 186; 277; 367]';
fixCashFlows_1Y = [367]';
floatCashFlows_2Y = [94 186 277 367 459 553 644 732]';
fixCashFlows_2Y = [367 732]';
floatCashFlows_3Y = [94 186 277 367 459 553 644 732 826 917 1008 1099]';
fixCashFlows_3Y = [367 732 1099]';
floatCashFlows_4Y = [94 186 277 367 459 553 644 732 826 917 1008 1099 1190 1281 1372]';
fixCashFlows_4Y = [367 732 1099 1372]';
floatCashFlows =  {floatCashFlows_1Y, floatCashFlows_2Y, floatCashFlows_3Y, floatCashFlows_4Y};
fixCashFlows = {fixCashFlows_1Y, fixCashFlows_2Y, fixCashFlows_3Y, fixCashFlows_4Y};
aZero = [A*EZero, zeros(n, n)]; 
aPi = [zeros(n, n), A*EPi];    
aZero_k = [A*EZero_k, zeros(n, kPi)]; 
aPi_k = [zeros(n, kZero), A*EPi_k];
N(1) = 1000;
N(2) = 1000;
N(3) = 1000;
N(4) = 1000;
numContracts = 1;
yield(1) = 0.769 / 100;
yield(2) = 0.9395 / 100;
yield(3) = 1.062 / 100;
yield(4) = 1.2006 / 100;
contractDays(1) = 253;
contractDays(2) = 504;
contractDays(3) = 1000;
contractDays(4) = 1000;
ropFix(1) = 'p';
ropFix(2) = 'r';
ropFix(3) = 'r';
ropFix(4) = 'p';
dcIbor = 'act/360';
dcFix = 'act/360';
fixingDate = 2;
%ttValueDate = fixingDate;
%%
epsP = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsI = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsA = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
sumRiskFactors = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
carry = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
NPV = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
Dt = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
cash = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
carryCash = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
%% Performance attribution for 1 year IRS
ttValueDate = fixingDate;   

accruedFloat = 0;
accruedFix = 0;

%Loop over all contracts
for i = 1:numContracts
    
    %Instantiate cash flows
    floatCashFlowsUnknown = floatCashFlows{i};  %All unknown floating cash flows
    floatCashFlowsKnown = floatCashFlows{i}(1); %All fixed floating cash flows
    fixCashFlowsCurr = fixCashFlows{i};         %All fix cash flows
    
    %Loop over the contract's lifetime
    for j = 1:min(contractDays(i), length(times) - 4)       
        
        %Handle the initial contract day when the price is determined
        if j == 1
            
            %The first floating cash flow is fixed according to Blomvall
            dt = handleDaycount(dcIbor, floatCashFlowsKnown(1));
            %floatcf = N(i) * dt * (r(floatCashFlowsKnown(1)+1) + pi(floatCashFlowsKnown(1) + 1));
               
            %The first floating cash flow is fixed
            blomvallFixing = N(i) * dt * (r(floatCashFlowsKnown(1)+1) + pi(floatCashFlowsKnown(1) + 1));
            liborFixing = N(i) * dt * libor(j,2);
            floatcf = liborFixing;
            
            %Dt is the difference in Libor and Blomvall model
            Dt{i}(j, 1) = liborFixing - blomvallFixing;

            %Set fix leg dt
            deltaTj(1) = handleDaycount(dcFix, fixCashFlowsCurr(1)-ttValueDate);
            if length(fixCashFlowsCurr) > 1
                for k = 2:length(fixCashFlowsCurr)
                    deltaTj(k) = handleDaycount(dcFix, fixCashFlowsCurr(k)-fixCashFlowsCurr(k-1));
                end
            end
            deltaTjInit = deltaTj;   
            
            %First fix cash flow
            nextFix = N(i) * yield(i) * deltaTj(1);
            
            %Set risk factors
            XiBar = [EZero' * f(:,j); EPi' * pi(:,j)]; 
            XiBarZero = [XiBar(1:n); zeros(n, 1)];
            XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];
            
            %Set first discount factor
            dcf = exp(ttValueDate/365 * r(ttValueDate+1) - fixCashFlowsCurr(1)/365 * r(fixCashFlowsCurr(1)+1));
            
            %Set price
            [accrualFloat, accrualFix] = accrual(liborFixing, nextFix, dcf, );
            currPrice = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZero, XiBarPi, floatcf, floatCashFlowsKnown, ropFix(i));
            %currPrice = irsPrice(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, r(:,j), piSpot(:,j), floatcf, floatCashFlowsKnown, dcIbor);
            
            %------------------- PRICE ERR ---------------------------
            epsP{i}(j, 1) = currPrice;
            
            %------------------- NPV ---------------------------------
            NPV{i}(j, 1) = epsP{i}(j, 1);
            
        else
                    
            %Change dates and handle cash flows
            [floatcf, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlowsCurr, Dt{i}(j, 1), deltaTj, cash{i}(j, 1), ttValueDate, liborFixing, blomvallFixing, dt] = ...
                handleDates(deltaTj, N(i), yield(i), floatcf, times, j, fixingDate, floatCashFlowsUnknown, ...
                floatCashFlowsKnown, fixCashFlowsCurr, r(:,j), pi(:,j), ropFix(i), libor(:,2), ttValueDate, dcIbor, dcFix, liborFixing, blomvallFixing, dt);           
                                  
            %Set risk factors
            XiBarPrev = XiBar;
            XiBarZeroPrev = XiBarZero;
            XiBarPiPrev = XiBarPi;
            dXi = [EZero_k' * (f(:,j)-f(:,j-1)); EPi_k' * (pi(:,j) - pi(:,j-1))];
            XiBar = [EZero' * f(:,j); EPi' * pi(:,j)]; 
            XiBarZero = [XiBar(1:n); zeros(n, 1)];
            XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];
            XiBardXiZero = XiBarZeroPrev + [dXi(1:kZero); zeros(2*n-kZero, 1)];
            XiBardXiPi = XiBarPiPrev + [zeros(n, 1); dXi(kZero+1:end); zeros(n-kPi, 1)];
            
            %Update price and calculate gradient and hessian
            prevPrice = currPrice;
            %currPrice = irsPrice(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, r(:,j), piSpot(:,j), floatcf, floatCashFlowsKnown, dcIbor);
            currPrice = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZero, XiBarPi, floatcf, floatCashFlowsKnown, ropFix(i));
            currPricePrevRiskFactor =  irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZeroPrev, XiBarPiPrev, floatcf, floatCashFlowsKnown, ropFix(i));
            currPriceDeltaRiskFactor = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBardXiZero, XiBardXiPi, floatcf, floatCashFlowsKnown, ropFix(i));
            gPrev = grad(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatcf, floatCashFlowsKnown, ropFix(i));
            HPrev = hes(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatcf, floatCashFlowsKnown, ropFix(i));

            %------------------- CARRY ----------------------------------------
            carry{i}(j, 1) = currPricePrevRiskFactor - prevPrice + Dt{i}(j, 1);
            dtCurr = handleDaycount(dcIbor, floatCashFlowsKnown(1));
            prevCarryCashFloat = currCarryCashFloat;
            
            
            prevCarryCashFix = currCarryCashFix;
            if ropFix(i) == 'r'
                currCarryCashFloat = (1 - dtCurr/dt(1)) * floatcf(1);
                currCarryCashFix = (1 - deltaTj(1)/deltaTjInit(totFixCashFlows - length(fixCashFlowsCurr)+1)) * N(i) * yield(i);
            elseif ropFix(i) == 'p'
                currCarryCashFloat = (1 - dtCurr/dt(1)) * floatcf(1);
                currCarryCashFix = -(1 - deltaTj(1)/deltaTjInit(totFixCashFlows - length(fixCashFlowsCurr)+1)) * N(i) * yield(i);
            end
            
            
            if abs(currCarryCashFloat) - abs(prevCarryCashFloat) >= 0
                accruedFloat = accruedFloat + currCarryCashFloat - prevCarryCashFloat;
                carryCash{i}(j, 1) = (currCarryCashFix - prevCarryCashFix) + currCarryCashFloat - prevCarryCashFloat;
            else
                accruedFloat = accruedFloat + currCarryCashFloat;
                carryCash{i}(j, 1) = (currCarryCashFix - prevCarryCashFix) + currCarryCashFloat;
            end
            
            %accruedFix = accruedFix + currCarryCashFix - prevCarryCashFix;
            
            %carryCash{i}(j, 1) = (currCarryCashFix - prevCarryCashFix); %+ currCarryCashFloat - prevCarryCashFloat;
     
            %------------------- TRUNC. ERR. ---------------------------------- 
            epsI{i}(j, 1) = currPrice - currPriceDeltaRiskFactor;
            
            %------------------- TAYL. APPROX. ERR ---------------------------- 
            epsA{i}(j, 1) = currPriceDeltaRiskFactor - currPricePrevRiskFactor - gPrev'*dXi - (1/2)*dXi'*HPrev*dXi;
            
            %------------------- Shift, twist, butterfly... -------------------    
            %dXiMat = diag(dXi);
            %shiftZero(i, 1) = gPrev(1)*dXi(1);
            %shift2Zero(i, 1) = (1/2) * dXiMat(:,1)' * HPrev * dXiMat(:,1);
            sumRiskFactors{i}(j, 1) = gPrev'*dXi + (1/2)*dXi'*HPrev*dXi;
           
            %------------------- NPV ------------------------------------------    
            %if j == 252 && i == 2
            %    return
            %end
            NPV{i}(j, 1) = currPrice - prevPrice + cash{i}(j, 1) - carryCash{i}(j, 1); %Dt{i}(j, 1)   %+  
        
            j
            
        end
    end
end
%%

for i = 1:numContracts
    cNPV{i} = cumsum(NPV{i});
    cEpsA{i} = cumsum(epsA{i});
    cEpsI{i} = cumsum(epsI{i});
    cEpsP{i} = cumsum(epsP{i});
    cCarry{i} = cumsum(carry{i});
    cSumRiskFactors{i} = cumsum(sumRiskFactors{i});
    %cCash{i} = cumsum(cash{i});
    cCarryCash{i} = cumsum(carryCash{i});
end

    totNPV = zeros(690, 1);
    totEpsA = zeros(690, 1);
    totEpsI = zeros(690, 1);
    totEpsP = zeros(690, 1);
    totCarry = zeros(690, 1);
    totSumRiskFactors = zeros(690, 1);
    totCarryCash = zeros(690, 1);
    totCarryTest = zeros(690, 1);
    
for i = 1:numContracts
    for j = 1:min(length(cNPV{i}), 690)
        totNPV(j, 1) = totNPV(j) + cNPV{i}(j);
        totEpsA(j, 1) = totEpsA(j) + cEpsA{i}(j);
        totEpsI(j, 1) = totEpsI(j) + cEpsI{i}(j);
        totCarry(j, 1) = totCarry(j) + cCarry{i}(j);
        totSumRiskFactors(j, 1) = totSumRiskFactors(j) + cSumRiskFactors{i}(j);
        totEpsP(j, 1) = totEpsP(j) + cEpsP{i}(j);
        totCarryCash(j, 1) = totCarryCash(j) + cCarryCash{i}(j);
    end
end

totCarryTest = totCarryCash + totCarry;
%
figure(1);
T = 1:253;
plot(T, cNPV{1}, T, cSumRiskFactors{1}, T,  cCarryCash{1} + cCarry{1}, T, cEpsA{1}, T, cEpsI{1}, T, cEpsP{1});%, T, cCash{1});

figure(2);
T = 1:690;
plot(T, totNPV(:, 1), T, totSumRiskFactors(:, 1), T, totCarryTest(:, 1), T, totEpsA(:, 1), T, totEpsI(:, 1), T, totEpsP(:, 1));


%plot(T, cNPV{i}, T, cSumRiskFactors{i}, T, cEpsA{i}, T, epsP{i});
%areaMatrix = [cNPV(1:252), cSumRiskFactors(1:252), cCarry(1:252), cEpsA(1:252), cEpsI(1:252), epsP(1:252)];

%area(T(1:252), areaMatrix)
legend('NPV', 'SumRiskFactors', 'carry', 'epsA', 'epsI', 'epsP');%, 'cash');
%legend('NPV', 'SumRiskFactors', 'epsA', 'epsP');

%validation = NPV{i} - sumRiskFactors{i} - epsP{i} - epsI{i} - carry{i} - epsA{i} - totCarryTest;
validation = totNPV - totSumRiskFactors - totEpsP - totEpsI - totCarryTest - totEpsA;
sum(abs(validation))




