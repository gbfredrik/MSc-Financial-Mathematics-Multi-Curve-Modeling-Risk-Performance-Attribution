function [floatcf, floatCashFlowsUnknown, floatCashFlowsKnown, fixCashFlows, Dt, deltaTj, cash, ttValueDate, liborFixing, blomvallFixing, dt] = ...
    handleDates(deltaTj, N, yield, floatcf, times, j, fixingDate, floatCashFlowsUnknown, ...
    floatCashFlowsKnown, fixCashFlows, r, pi, ropFix, libor, ttValueDate, dcIbor, dcFix, liborFixing, blomvallFixing, dt)

%Decrease time to value date with one day
dateDiff = times(j) - times(j - 1);
ttValueDate = max(0, ttValueDate - dateDiff);

%Check if two (fixingDate) days away from float cashflow
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

%Check if two (fixingDate) days away from fix cashflow
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

%Decrease days to cashflows
dateDiff = times(j) - times(j-1);
floatCashFlowsUnknown = floatCashFlowsUnknown - dateDiff;
floatCashFlowsKnown = floatCashFlowsKnown - dateDiff;
fixCashFlows = fixCashFlows - dateDiff;



%Case 1: day 0 for a contract is handled separately in main function

%Case 2: day when a floating coupon is fixed
if ((workdaysToFloatPayment == fixingDate) && (length(floatCashFlowsUnknown) > 1))  
   
   floatCashFlowsKnown = [floatCashFlowsKnown, floatCashFlowsUnknown(2)]; %New floating cash flow is determined
   dt = [dt, handleDaycount(dcIbor, floatCashFlowsKnown(2))];
   liborFixing = N * dt(2) * libor(j); 
   blomvallFixing = N * dt(2) * (r(floatCashFlowsKnown(2)+1) + pi(floatCashFlowsKnown(2) + 1));
   floatcf = [floatcf, blomvallFixing]; %Set cash flow size
   floatCashFlowsUnknown = floatCashFlowsUnknown(2:end); %Reduce number of unknown cash flows
   
   if ropFix == 'r'
        Dt = -(liborFixing - blomvallFixing) %Difference between Blomvall and actual cash flow
   elseif ropFix == 'p'
        Dt = liborFixing - blomvallFixing; %Difference between Blomvall and actual cash flow
   end
   
   cash = 0; %No payout on these dates
   
%Case 3: date of float cash flow
elseif (((floatCashFlowsKnown(1) == 0) && (fixCashFlows(1) > 0) && length(floatCashFlowsKnown) > 1))
    
    %floatCashFlowsKnown = floatCashFlowsKnown(2); %Known float cash flows has decreased
    %floatcf = floatcf(2); 
    
    %Cash is payed
    %if ropFix == 'r'
    %    cash = -floatcf(1); 
    %elseif ropFix == 'p'
    %    cash = floatcf(1);   
    %end
    cash = 0;
    dt = dt(2);
    Dt = 0; %No fixing difference on these dates
    
elseif (((floatCashFlowsKnown(1) < 0) && (fixCashFlows(1) > 0) && length(floatCashFlowsKnown) > 1))
    
    floatCashFlowsKnown = floatCashFlowsKnown(2); %Known float cash flows has decreased
    floatcf = floatcf(2); 
    cash = 0;
    Dt = 0;
    
% case 4: date of fix and float cash flow
elseif ((floatCashFlowsKnown(1) == 0) && (fixCashFlows(1) == 0) && (length(floatCashFlowsKnown) == 1))
    
    %Do if contract doesn't mature
    %if length(floatCashFlowsUnknown) > 1
        
    %    floatCashFlowsKnown = floatCashFlowsKnown(2); %Known float cash flows has decreased
    %    floatcf = floatcf(2);
    %    fixCashFlows = fixCashFlows(2:end); %Known fix cash flows has decreased

    %end
    
    %Cash payout is the difference in floating and fix
    if ropFix == 'r'
        cash = N * yield - floatcf(1);
    elseif ropFix == 'p'
        cash = -N * yield + floatcf(1);
    end
    
    Dt = 0; %No fixing difference on these dates
    
elseif ((floatCashFlowsKnown(1) < 0) && (fixCashFlows(1) < 0))
    
    if length(floatCashFlowsUnknown) > 1
        
        floatCashFlowsKnown = floatCashFlowsKnown(2); %Known float cash flows has decreased
        floatcf = floatcf(2);
        fixCashFlows = fixCashFlows(2:end); %Known fix cash flows has decreased

    end  
    
    cash = 0;
    Dt = 0;
    
%If regular day
else

    cash = 0;
    Dt = 0;
    
end
%Update fix leg dt
deltaTj(1) = handleDaycount(dcFix, fixCashFlows(1)-ttValueDate);
if length(fixCashFlows) > 1
    for k = 2:length(fixCashFlows) - 1
        deltaTj(k) = handleDaycount(dcFix, fixCashFlows(k+1)-fixCashFlows(k));
    end
end
end