function [floatCashFlows, fixCashFlows, daysToNextFloat, daysToNextFix, ttValueDate, dtFixNext, dtFloat, liborFixing, blomvallFixing, Dt, nextFix, workdaysToFloatPayment, last, workdaysToFixPayment, liborFixingPrev, dtFix] = ...
            handleDates(N, yield, times, j, fixingDate, floatCashFlows, fixCashFlows, daysToNextFloat, daysToNextFix, ...
            r, pi, ropFix, libor, ttValueDate, dcIbor, dcFix, liborFixing, blomvallFixing, dtFixNext, dtFloat, nextFix, last, liborFixingPrev, dtFix)


%Decrease time to value date with one day
dateDiff = times(j) - times(j - 1);
ttValueDate = max(0, ttValueDate - dateDiff);

%Check if two (fixingDate) days away from float cashflow
daysToNextFloatTemp = daysToNextFloat;
workdaysToFloatPayment = 0;
i = j;
while daysToNextFloatTemp >= 0
    dateDiff = times(i) - times(i - 1);
    daysToNextFloatTemp = daysToNextFloatTemp - dateDiff;
    if daysToNextFloatTemp <= 0
        break
    end
    workdaysToFloatPayment = workdaysToFloatPayment + 1;
    i = i + 1;
    if workdaysToFloatPayment > 4
        break
    end
end

%Check if two (fixingDate) days away from fix cashflow
daysToNextFixTemp = daysToNextFix;
workdaysToFixPayment = 0;
i = j;
while daysToNextFixTemp >= 0
    dateDiff = times(i) - times(i - 1);
    daysToNextFixTemp = daysToNextFixTemp - dateDiff;
    if daysToNextFixTemp <= 0
        break
    end
    workdaysToFixPayment = workdaysToFixPayment + 1;
    i = i + 1;
    if workdaysToFixPayment > 4
        break
    end
end

%Decrease days to cashflows
dateDiff = times(j) - times(j-1);
floatCashFlows = floatCashFlows - dateDiff;
fixCashFlows = fixCashFlows - dateDiff;
daysToNextFloat = daysToNextFloat - dateDiff;
daysToNextFix = daysToNextFix - dateDiff;

%Case 1: day 0 for a contract is handled separately in main function

%Case 2: day when a floating coupon is fixed
if ((workdaysToFloatPayment == fixingDate) && (length(floatCashFlows) > 1)) && (workdaysToFixPayment > fixingDate)   

   dtFloat = handleDaycount(dcIbor, floatCashFlows(2));
   liborFixingPrev = liborFixing;
   liborFixing = N * dtFloat * libor(j);
   %blomvallFixing = N * dtFloat * (r(daysToNextFloat + 1) + pi(daysToNextFloat + 1))
   blomvallFixing = N * dtFloat * (r(floatCashFlows(2) + 1) + pi(floatCashFlows(2) + 1));
   
   if ropFix == 'r'
        Dt = -(blomvallFixing - liborFixing); %Difference between Blomvall and actual cash flow
   elseif ropFix == 'p'
        Dt = blomvallFixing - liborFixing; %Difference between Blomvall and actual cash flow
   end
 
elseif ((workdaysToFloatPayment == 0) && (workdaysToFixPayment == 0) && (length(floatCashFlows) == 1))

    last = 1;
    Dt = 0;

% case 4: date of fix and float cash flow
elseif ((workdaysToFloatPayment == fixingDate) && (workdaysToFixPayment == fixingDate) && (length(floatCashFlows) > 1))
       
   dtFloat = handleDaycount(dcIbor, floatCashFlows(2)); %ny
   liborFixingPrev = liborFixing;
   liborFixing = N * dtFloat * libor(j);

   blomvallFixing = N * dtFloat * (r(floatCashFlows(2) + 1) + pi(floatCashFlows(2) + 1));
   
    nextFix = N * yield * dtFix(2); %ny
    
   if ropFix == 'r'
        Dt = -(blomvallFixing - liborFixing); %Difference between Blomvall and actual cash flow
   elseif ropFix == 'p'
        Dt = blomvallFixing - liborFixing; %Difference between Blomvall and actual cash flow
   end
    
%If regular day
else

    Dt = 0;

end

end