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
numContracts = 4;
yield(1) = 0.769 / 100;
yield(2) = 0.9395 / 100;
yield(3) = 1.062 / 100;
yield(4) = 1.2006 / 100;
contractDays(1) = 253;
contractDays(2) = 504;
contractDays(3) = 1000;
contractDays(4) = 1000;
%%
epsP = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsI = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsA = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
sumRiskFactors = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
carry = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
NPV = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
Dt = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
%% Performance attribution for 1 year IRS

    fixingDate = 2;
    startdate = 2;
for i = 1:numContracts
    
    floatCashFlowsUnknown = floatCashFlows{i};
    floatCashFlowsKnown = floatCashFlows{i}(1);
    fixCashFlowsCurr = fixCashFlows{i};
    startdate = 2;

    for j = 1:min(contractDays(i), length(times) - 4)       
        if j == 1
            floatcf = N(i) * (exp(floatCashFlowsKnown(1)/365 * (r(floatCashFlowsKnown(1)+1,j)+pi(floatCashFlowsKnown(1)+1,j)) - startdate/365 * (r(startdate+1,j)+pi(startdate+1,j)))-1);
            deltaTj(1) = (fixCashFlowsCurr(1)-startdate) / 360;
            if length(fixCashFlowsCurr) > 1
                for k = 2:length(fixCashFlowsCurr)
                    deltaTj(k) = (fixCashFlowsCurr(k)-fixCashFlowsCurr(k-1)) / 360;
                end
            end
            XiBar = [EZero' * f(:,j); EPi' * pi(:,j)]; 
            XiBarZero = [XiBar(1:n); zeros(n, 1)];
            XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];
            %------------------- PRICE ERR ---------------------------
            %currPrice = irsPrice(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, r(:,j), piSpot(:,j), floatcf, floatCashFlowsKnown);
            currPrice = irsPriceRiskFactor(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZero, XiBarPi, floatcf, floatCashFlowsKnown);
            epsP{i}(j, 1) = currPrice;
            NPV{i}(j, 1) = epsP{i}(j, 1);
            
            startdate = startdate - 1;
        else
            
            % change dates
            [floatcf, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlowsCurr, Dt{i}(j, 1), deltaTj] = handleDates(deltaTj, N(i), yield(i), floatcf, times, j, fixingDate, startdate, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlowsCurr, r(:,j), pi(:,j));           
                       
            % Set risk factors
            XiBarPrev = XiBar;
            XiBarZeroPrev = XiBarZero;
            XiBarPiPrev = XiBarPi;
            dXi = [EZero_k' * (f(:,j)-f(:,j-1)); EPi_k' * (pi(:,j) - pi(:,j-1))];
            XiBar = [EZero' * f(:,j); EPi' * pi(:,j)]; 
            XiBarZero = [XiBar(1:n); zeros(n, 1)];
            XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];
            XiBardXiZero = XiBarZeroPrev + [dXi(1:kZero); zeros(2*n-kZero, 1)];
            XiBardXiPi = XiBarPiPrev + [zeros(n, 1); dXi(kZero+1:end); zeros(n-kPi, 1)];
            
            
            % Calc variables
            prevPrice = currPrice;
            %currPrice = irsPrice(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, r(:,j), piSpot(:,j), floatcf, floatCashFlowsKnown);
            currPrice = irsPriceRiskFactor(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZero, XiBarPi, floatcf, floatCashFlowsKnown);
            currPricePrevRiskFactor =  irsPriceRiskFactor(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBarZeroPrev, XiBarPiPrev, floatcf, floatCashFlowsKnown);
            currPriceDeltaRiskFactor = irsPriceRiskFactor(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero, aPi, XiBardXiZero, XiBardXiPi, floatcf, floatCashFlowsKnown);
            gPrev = grad(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatcf, floatCashFlowsKnown);
            HPrev = hes(N(i), yield(i), startdate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatcf, floatCashFlowsKnown);

            %------------------- CARRY ----------------------------------------
            carry{i}(j, 1) = currPricePrevRiskFactor - prevPrice + Dt{i}(j, 1);
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
            %if j == 129
            %    return
            %end
            NPV{i}(j, 1) = currPrice - prevPrice + Dt{i}(j, 1);
        
            j
            
            if (j == 2)
                startdate = startdate - 1;
            end
        end
    end
end
%%

for i = 1:numContracts
    cNPV{i} = cumsum(NPV{i});
    cEpsA{i} = cumsum(epsA{i});
    cEpsI{i} = cumsum(epsI{i});
    cCarry{i} = cumsum(carry{i});
    cSumRiskFactors{i} = cumsum(sumRiskFactors{i});
end

    totNPV = zeros(690, 1);
    totEpsA = zeros(690, 1);
    totEpsI = zeros(690, 1);
    totCarry = zeros(690, 1);
    totSumRiskFactors = zeros(690, 1);
    totEpsP = zeros(690, 1);
    
for i = 1:numContracts
    for j = 1:min(length(cNPV{i}), 690)
        totNPV(j, 1) = totNPV(j) + cNPV{i}(j);
        totEpsA(j, 1) = totEpsA(j) + cEpsA{i}(j);
        totEpsI(j, 1) = totEpsI(j) + cEpsI{i}(j);
        totCarry(j, 1) = totCarry(j) + cCarry{i}(j);
        totSumRiskFactors(j, 1) = totSumRiskFactors(j) + cSumRiskFactors{i}(j);
        totEpsP(j, 1) = totEpsP(j) + epsP{i}(j);
    end
end

%T = 1:253;
%plot(T, cNPV{i}, T, cSumRiskFactors{i}, T, cCarry{i}, T, cEpsA{i}, T, cEpsI{i}, T, epsP{i});

T = 1:690;
plot(T, totNPV(:, 1), T, totSumRiskFactors(:, 1), T, totCarry(:, 1), T, totEpsA(:, 1), T, totEpsI(:, 1), T, totEpsP(:, 1));


%plot(T, cNPV{i}, T, cSumRiskFactors{i}, T, cEpsA{i}, T, epsP{i});
%areaMatrix = [cNPV(1:252), cSumRiskFactors(1:252), cCarry(1:252), cEpsA(1:252), cEpsI(1:252), epsP(1:252)];

%area(T(1:252), areaMatrix)
legend('NPV', 'SumRiskFactors', 'carry', 'epsA', 'epsI', 'epsP');
%legend('NPV', 'SumRiskFactors', 'epsA', 'epsP');

validation = NPV{i} - sumRiskFactors{i} - epsP{i} - epsI{i} - carry{i} - epsA{i};
sum(abs(validation))




