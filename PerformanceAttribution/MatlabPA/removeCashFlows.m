function [floatCashFlows, fixCashFlows, daysToNextFloat, daysToNextFix, dtFix, dtFixNext] = removeCashFlows(workdaysToFloatPayment, workdaysToFixPayment, floatCashFlows, fixCashFlows, fixingDate, daysToNextFloat, ...
    daysToNextFix, dtFix, dtFixNext)

if ((workdaysToFloatPayment == 0) && (length(floatCashFlows) > 1)) && (workdaysToFixPayment > fixingDate)%ny
   floatCashFlows = floatCashFlows(2:end); %ny
   daysToNextFloat = floatCashFlows(1); %ny


elseif (workdaysToFloatPayment == 0) && (workdaysToFixPayment == 0) && (length(floatCashFlows) > 1)

   floatCashFlows = floatCashFlows(2:end); %ny
   daysToNextFloat = floatCashFlows(1); %ny
   
   dtFix = dtFix(2:end);
   dtFixNext = dtFix(1);
   daysToNextFix = fixCashFlows(2);
   if length(fixCashFlows) > 1
       fixCashFlows = fixCashFlows(2:end);
   else
       fixCashFlows = fixCashFlows(2:end);
   end
   
end


end