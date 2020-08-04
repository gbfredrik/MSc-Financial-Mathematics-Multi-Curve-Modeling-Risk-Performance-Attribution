function [carry] = calcCarry(currPricePrevRiskFactor, prevPrice, workdaysToFloatPayment, workdaysToFixPayment, last, ropFix, liborFixing, liborFixingPrev, nextFix)



if workdaysToFloatPayment == 0 && workdaysToFixPayment ~= 0 && last == 0 && ropFix == 'p'
    carry = currPricePrevRiskFactor - prevPrice + liborFixingPrev; 
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment ~= 0 && last == 0 && ropFix == 'r'
    carry = currPricePrevRiskFactor - prevPrice - liborFixingPrev;% 
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 0 && ropFix == 'p'
    carry = currPricePrevRiskFactor - prevPrice + liborFixingPrev - nextFix;  
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 0 && ropFix == 'r'
    carry = currPricePrevRiskFactor - prevPrice - liborFixingPrev + nextFix;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 1 && ropFix == 'r'
    carry = currPricePrevRiskFactor - prevPrice;% - nextFix + liborFixing;
elseif workdaysToFloatPayment == 0 && workdaysToFixPayment == 0 && last == 1 && ropFix == 'p'
    carry = currPricePrevRiskFactor - prevPrice;% + nextFix - liborFixing;
else
    carry = currPricePrevRiskFactor - prevPrice;
end


end