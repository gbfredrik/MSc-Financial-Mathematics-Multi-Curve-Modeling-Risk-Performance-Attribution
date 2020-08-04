function [NPV] = calcNPV(workdaysToFloatPayment, workdaysToFixPayment, last, ropFix, currPrice, prevPrice, liborFixingPrev, liborFixing, nextFix)

if workdaysToFloatPayment == 0 && workdaysToFixPayment ~= 0 && last == 0 && ropFix == 'p'
    NPV = currPrice - prevPrice + liborFixingPrev; 
    liborFixingPrev;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment ~= 0 && last == 0 && ropFix == 'r'
    NPV = currPrice - prevPrice - liborFixingPrev;% 
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 0 && ropFix == 'p'
    NPV = currPrice - prevPrice + liborFixingPrev - nextFix;  
    liborFixingPrev;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 0 && ropFix == 'r'
    NPV = currPrice - prevPrice - liborFixingPrev + nextFix;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 1 && ropFix == 'r'
    NPV = currPrice - prevPrice;% - nextFix + liborFixing;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 1 && ropFix == 'p'
    NPV = currPrice - prevPrice;% + nextFix - liborFixing;
    liborFixing;
else
    NPV = currPrice - prevPrice;
end



end