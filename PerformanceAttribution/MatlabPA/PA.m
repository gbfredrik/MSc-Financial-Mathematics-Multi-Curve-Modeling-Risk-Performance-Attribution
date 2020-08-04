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
floatCashFlows = {};
fixCashFlows = {};
yield = [];
N = [];
contractDays = [];
loop = 'A':'D';
for i = 1:length(loop)
    currColumn = strcat(loop(i), ':', loop(i));
    floatCashFlows = [floatCashFlows, xlsread('Data/PortfolioData.xlsx', 'Float Cashflows', currColumn)];
    fixCashFlows = [fixCashFlows, xlsread('Data/PortfolioData.xlsx', 'Fix Cashflows', currColumn)];
    yield = [yield, xlsread('Data/PortfolioData.xlsx', 'Yield', currColumn)];
    N = [yield, xlsread('Data/PortfolioData.xlsx', 'Nominal', currColumn)];
    contractDays = [contractDays, xlsread('Data/PortfolioData.xlsx', 'ContractDays', currColumn)];
end
n = size(fHistIS, 2);
A = intMatrix(n);
f = fHistIS';
pi = piHistIS';
r = A * f;
piSpot = A * pi;
times = times(2758:end);

%Calculate eigenvector matrices
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
aZero = [A*EZero, zeros(n, n)]; 
aPi = [zeros(n, n), A*EPi];    
aZero_k = [A*EZero_k, zeros(n, kPi)]; 
aPi_k = [zeros(n, kZero), A*EPi_k];

%% Set variables
ropFix(1) = 'r';
ropFix(2) = 'p';
ropFix(3) = 'r';
ropFix(4) = 'r';
dcIbor = 'act/360';
dcFix = 'act/360';
fixingDate = 2;
%%



epsP = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsI = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
epsA = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
sumRiskFactors = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
carry = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
NPV = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
Dt = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
accrualCumulative = { zeros(contractDays(1), 1), zeros(contractDays(2), 1), zeros(contractDays(3), 1), zeros(contractDays(4), 1) };
%% Performance attribution for 1 year IRS
numContracts = 4;
%Loop over all contrac+ts
for i = 1:numContracts

    workdaysToFixPayment = 1;
    workdaysToFloatPayment = 1;
    
    ttValueDate = fixingDate;   
    nextFix = 0;
    liborFixingPrev = 0;
    %Instantiate cash flows
    floatCashFlowsCurr = floatCashFlows{i};  %All unknown floating cash flows
    fixCashFlowsCurr = fixCashFlows{i};         %All fix cash flows
    daysToNextFloat = floatCashFlowsCurr(1);
    daysToNextFix = fixCashFlowsCurr(1);
    
    dtFloat = handleDaycount(dcIbor, daysToNextFloat - ttValueDate);
    dtFixNext =  handleDaycount(dcFix, daysToNextFix - ttValueDate);
    
    last = 0;
    
    %Loop over the contract's lifetime
    for j = 1:min(contractDays(i), length(times) - 4)       
        
        %Handle the initial contract day when the price is determined
        if j == 1
            
            %The first floating cash flow is fixed according to Blomvall
            dtFloat = handleDaycount(dcIbor, daysToNextFloat - ttValueDate);
               
            %The first floating cash flow is fixed
            blomvallFixing = N(i) * dtFloat * (r(daysToNextFloat+1) + pi(daysToNextFloat + 1));
            liborFixing = N(i) * dtFloat * libor(j,2);
            %Dt is the difference in Libor and Blomvall model
            if ropFix(i) == 'r'
                 Dt{i}(j, 1) = -(blomvallFixing - liborFixing); %Difference between Blomvall and actual cash flow
            elseif ropFix(i) == 'p'
                 Dt{i}(j, 1) = blomvallFixing - liborFixing; %Difference between Blomvall and actual cash flow
            end
            
            %Set fix leg dt
            dtFix(1) = handleDaycount(dcFix, daysToNextFix - ttValueDate);
            if length(fixCashFlowsCurr) > 1
                for k = 2:length(fixCashFlowsCurr)
                    dtFix(k) = handleDaycount(dcFix, fixCashFlowsCurr(k)-fixCashFlowsCurr(k-1));
                end
            end
            
            %First fix cash flow
            nextFix = N(i) * yield(i) * dtFixNext(1);


            %Set risk factors
            XiBar = [EZero' * f(:,j); EPi' * pi(:,j)]; 
            XiBarZero = [XiBar(1:n); zeros(n, 1)];
            XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];
            
            %Set first discount factor
            dtFloatCurr = handleDaycount(dcIbor, daysToNextFloat);
            dtFixCurr = handleDaycount(dcFix, daysToNextFix);
            
            %Set price
            
            currPrice = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBarZero, XiBarPi, ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);
            accrualFloatPrev = 0;
            accrualFixPrev = 0;

            prevPrice = currPrice;
            
            %------------------- PRICE ERR ---------------------------
            epsP{i}(j, 1) = currPrice;
            
            %------------------- NPV ---------------------------------
            NPV{i}(j, 1) = epsP{i}(j, 1);
            
        else
            
            wdToFloatPrev = workdaysToFloatPayment;
            wdToFixPrev = workdaysToFixPayment;
            %Change dates and handle cash flows
            [floatCashFlowsCurr, fixCashFlowsCurr, daysToNextFloat, daysToNextFix, ttValueDate, dtFixNext, dtFloat, liborFixing, blomvallFixing, Dt{i}(j,1), nextFix, workdaysToFloatPayment, last, workdaysToFixPayment, liborFixingPrev, dtFix] = ...
                handleDates(N(i), yield(i), times, j, fixingDate, floatCashFlowsCurr, fixCashFlowsCurr, daysToNextFloat, daysToNextFix, r(:,j), pi(:,j), ropFix(i), libor(:,2), ...
                ttValueDate, dcIbor, dcFix, liborFixing, blomvallFixing, dtFixNext, dtFloat, nextFix, last, liborFixingPrev, dtFix);
            
            
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
    
            %currPrice = irsPrice(N(i), yield(i), ttValueDate, floatCashFlowsUnknown, fixCashFlowsCurr, deltaTj, r(:,j), piSpot(:,j), floatcf, floatCashFlowsKnown, dcIbor);
            if (workdaysToFloatPayment <= 2) && length(floatCashFlowsCurr) > 1
                currPrice = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBarZero, XiBarPi, ropFix(i), nextFix, liborFixingPrev, daysToNextFix, daysToNextFloat);
                currPricePrevRiskFactor =  irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBarZeroPrev, XiBarPiPrev, ropFix(i), nextFix, liborFixingPrev, daysToNextFix, daysToNextFloat);
                currPriceDeltaRiskFactor = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBardXiZero, XiBardXiPi, ropFix(i), nextFix, liborFixingPrev, daysToNextFix, daysToNextFloat);
                gPrev = grad(N(i), yield(i), ttValueDate, fixCashFlowsCurr, dtFix, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatCashFlowsCurr, ropFix(i), nextFix, liborFixingPrev, daysToNextFix, daysToNextFloat);
                HPrev = hes(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), ropFix(i), nextFix, liborFixingPrev, daysToNextFix, daysToNextFloat);
            else
                currPrice = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBarZero, XiBarPi, ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);   
                currPricePrevRiskFactor =  irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBarZeroPrev, XiBarPiPrev, ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);
                currPriceDeltaRiskFactor = irsPriceRiskFactor(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero, aPi, XiBardXiZero, XiBardXiPi, ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);
                gPrev = grad(N(i), yield(i), ttValueDate, fixCashFlowsCurr, dtFix, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), floatCashFlowsCurr, ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);      
                HPrev = hes(N(i), yield(i), ttValueDate, floatCashFlowsCurr, fixCashFlowsCurr, dtFix, aZero_k, aPi_k, r(:,j-1), piSpot(:,j-1), ropFix(i), nextFix, liborFixing, daysToNextFix, daysToNextFloat);
            end
            
       
            dtFixCurr = handleDaycount(dcFix, daysToNextFix);
            dtFloatCurr = handleDaycount(dcIbor, daysToNextFloat);
            dcfFix = exp(ttValueDate/365 * r(ttValueDate + 1) - (daysToNextFix)/365 * r(daysToNextFix + 1)); 
            dcfFloat = exp(ttValueDate/365 * r(ttValueDate + 1) - daysToNextFloat/365 * r(daysToNextFloat + 1));       
            %utbetalning rörlig
            if (workdaysToFloatPayment == 0) && (workdaysToFixPayment > 0)
               [accrualCumulative{i}(j, 1), accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixingPrev, nextFix, dcfFix, 1, 0, dtFloat, dtFixCurr, dtFixNext, ropFix(i), accrualFloatPrev, accrualFixPrev);
            %utbetalning rörlig och fix
            elseif (workdaysToFloatPayment == 0) && (workdaysToFixPayment == 0) && (last == 0) 
               [accrualCumulative{i}(j, 1), accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixingPrev, nextFix, 1, 1, 0, dtFloat, 0, dtFixNext, ropFix(i), accrualFloatPrev, accrualFixPrev);
            %utbetalning rörlig och fix, sista  
            elseif (workdaysToFloatPayment == 0) && (workdaysToFixPayment == 0) && (last == 1)    
               [accrualCumulative{i}(j, 1), accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixing, nextFix, 1, 1, 0, dtFloat, 0, dtFixNext, ropFix(i), accrualFloatPrev, accrualFixPrev);
            % < fixingdays dagar kvar till rörlig, ej inför sista
            elseif (workdaysToFloatPayment <= 2) && (fixCashFlowsCurr(1) > 10 || length(floatCashFlowsCurr) > 2)
               [accrualCumulative{i}(j, 1), accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixingPrev, nextFix, dcfFix, dcfFloat, dtFloatCurr, dtFloat, dtFixCurr, dtFixNext, ropFix(i), accrualFloatPrev, accrualFixPrev);
            %Före första dagen
            elseif ttValueDate ~= 0 
               accrualCumulative{i}(j, 1) = 0;
               accrualFloat = 0;
               accrualFix = 0; 
               accrualFloatPrev = 0;
               accrualFixPrev = 0;                  
            %vanliga dagar
            else    
               [accrualCumulative{i}(j, 1), accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixing, nextFix, dcfFix, dcfFloat, dtFloatCurr, dtFloat, dtFixCurr, dtFixNext, ropFix(i), accrualFloatPrev, accrualFixPrev);               
            end
            
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
            if j == 381
            %    return
            end
            
            %------------------- CARRY ----------------------------------------
            carry{i}(j, 1) = calcCarry(currPricePrevRiskFactor, prevPrice, wdToFloatPrev, wdToFixPrev, last, ropFix(i), liborFixing, liborFixingPrev, nextFix);
            carry{i}(j, 1) = carry{i}(j, 1);% + Dt{i}(j, 1);
            NPV{i}(j, 1) = calcNPV(wdToFloatPrev, wdToFixPrev, last, ropFix(i), currPrice, prevPrice, liborFixingPrev, liborFixing, nextFix);
            NPV{i}(j, 1) = NPV{i}(j, 1);% + accrualCumulative{i}(j, 1);
            
            j
            prevPrice = currPrice;
           
            
            [floatCashFlowsCurr, fixCashFlowsCurr, daysToNextFloat, daysToNextFix, dtFix, dtFixNext] = removeCashFlows(workdaysToFloatPayment, workdaysToFixPayment, ...
                floatCashFlowsCurr, fixCashFlowsCurr, fixingDate, daysToNextFloat, daysToNextFix, dtFix, dtFixNext);
            
        end
    end
end

%%

    totNPV = zeros(504, 1);
    totEpsA = zeros(504, 1);
    totEpsI = zeros(504, 1);
    totEpsP = zeros(504, 1);
    totCarry = zeros(504, 1);
    totSumRiskFactors = zeros(504, 1);
    
for i = 1:numContracts
    for j = 1:min(length(NPV{i}), 504)
        totNPV(j, 1) = totNPV(j) + NPV{i}(j);
        totEpsA(j, 1) = totEpsA(j) + epsA{i}(j);
        totEpsI(j, 1) = totEpsI(j) + epsI{i}(j);
        totCarry(j, 1) = totCarry(j) + carry{i}(j);
        totSumRiskFactors(j, 1) = totSumRiskFactors(j) + sumRiskFactors{i}(j);
        totEpsP(j, 1) = totEpsP(j) + epsP{i}(j);
    end
end


cNPV = cumsum(totNPV(:, 1));
cEpsA = cumsum(totEpsA(:, 1));
cEpsI = cumsum(totEpsI(:, 1));
cEpsP = cumsum(totEpsP(:, 1));
cCarry = cumsum(totCarry(:, 1));
cSumRiskFactors = cumsum(totSumRiskFactors(:, 1));

T = 1:504;
plot(T, cNPV(:,1), T, cSumRiskFactors(:, 1), T, cEpsA(:, 1), T, cEpsI(:, 1), T, cEpsP(:, 1), T, cCarry(:, 1));

%figure(1);
%T = 1:253;
%plot(T, cNPV{1}, T, cSumRiskFactors{1}, T,  cCarryCash{1} + cCarry{1}, T, cEpsA{1}, T, cEpsI{1}, T, cEpsP{1});%, T, cCash{1});

%figure(2);
%T = 1:690;
%plot(T, totNPV(:, 1), T, totSumRiskFactors(:, 1), T, totCarryTest(:, 1), T, totEpsA(:, 1), T, totEpsI(:, 1), T, totEpsP(:, 1));


%plot(T, cNPV{i}, T, cSumRiskFactors{i}, T, cEpsA{i}, T, epsP{i});
%areaMatrix = [cNPV(1:252), cSumRiskFactors(1:252), cCarry(1:252), cEpsA(1:252), cEpsI(1:252), epsP(1:252)];

%area(T(1:252), areaMatrix)
legend('NPV', 'SumRiskFactors', 'epsA', 'epsI', 'epsP', 'Carry');
%legend('NPV', 'SumRiskFactors', 'epsA', 'epsP');

validation = totNPV(:,1) - totSumRiskFactors(:, 1) - totEpsA(:, 1) - totEpsP(:, 1) - totEpsI(:, 1) - totCarry(:, 1);
%validation = totNPV - totSumRiskFactors - totEpsP - totEpsI - totCarryTest - totEpsA;
sum(abs(validation))




