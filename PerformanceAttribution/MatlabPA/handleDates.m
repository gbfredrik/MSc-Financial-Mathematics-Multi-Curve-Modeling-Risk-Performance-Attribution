function [floatcf, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlows, Dt, deltaTj] = handleDates(deltaTj, N, yield, floatcf, times, j, fixingDate, startdate, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlows, r, pi)

floatCashKnownTemp = floatCashFlowsKnown;
workdaysToFloatPayment = 0;
i = j;
while floatCashKnownTemp >= fixingDate
    dateDiff = times(i) - times(i - 1);
    floatCashKnownTemp = floatCashKnownTemp - dateDiff;
    workdaysToFloatPayment = workdaysToFloatPayment + 1;
    i = i + 1;
    if workdaysToFloatPayment > 2
        break
    end
end

fixCashTemp = fixCashFlows;
workdaysToFixPayment = 0;
i = j;
while fixCashTemp >= fixingDate
    dateDiff = times(i) - times(i - 1);
    fixCashTemp = fixCashTemp - dateDiff;
    workdaysToFixPayment = workdaysToFixPayment + 1;
    i = i + 1;
    if workdaysToFixPayment > 2
        break
    end
end

dateDiff = times(j) - times(j-1);
floatCashFlowsUnknown = floatCashFlowsUnknown - dateDiff;
floatCashFlowsKnown = floatCashFlowsKnown - dateDiff;
fixCashFlows = fixCashFlows - dateDiff;

% case 1: dag 0 för ett kontrakt (settlement - 2)
% hanteras utanför

% case 2: dag då rörlig kupong fixeras men ej fix byts
if ((workdaysToFloatPayment == fixingDate) && (length(floatCashFlowsUnknown) > 1)) %&& (workdaysToFixPayment > fixingDate))   
   
   floatCashFlowsKnown = [floatCashFlowsKnown, floatCashFlowsUnknown(2)];
   floatcf = [floatcf, N * (exp(floatCashFlowsKnown(2)/365 * (r(floatCashFlowsKnown(2)+1)+pi(floatCashFlowsKnown(2)+1)) - startdate/365 * (r(startdate+1)+pi(startdate+1)))-1)];
   Dt = 0;
   floatCashFlowsUnknown = floatCashFlowsUnknown(2:end);
% case 3: gammla rörliga faller bort
elseif (((floatCashFlowsKnown(1) == 0) && (fixCashFlows(1) > 0) && length(floatCashFlowsKnown) > 1))
    
    floatCashFlowsKnown = floatCashFlowsKnown(2);
    Dt = floatcf(1);
    floatcf = floatcf(2);
    %floatCashFlowsUnknown = floatCashFlowsUnknown(2:end);
    
% case 4: gammal fix samt rörlig faller bort
elseif ((floatCashFlowsKnown(1) == 0) && (fixCashFlows(1) == 0))
    
    if length(floatCashFlowsUnknown) > 1
        floatCashFlowsKnown = floatCashFlowsKnown(2);
        fixCashFlows = fixCashFlows(2:end);
        floatcf = floatcf(2);
    end
    
    Dt = floatcf(1) - yield * N;
    
else

    Dt = 0;
    
end


    deltaTj(1) = (fixCashFlows(1) - startdate) / 360;
    if length(fixCashFlows) > 1
        for k = 2:length(fixCashFlows) - 1
            deltaTj(k) = (fixCashFlows(k+1)-fixCashFlows(k)) / 360;
        end
    end
            









end