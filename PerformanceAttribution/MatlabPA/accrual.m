function [accrualCumulative, accrualFloat, accrualFix, accrualFloatPrev, accrualFixPrev] = accrual(liborFixing, nextFix, dcFix, dcFloat, dtFloatCurr, dtFloatInit, dtFixCurr, dtFixInit, ropFix, accrualFloatPrev, accrualFixPrev)

fracFloat = dtFloatCurr/dtFloatInit;
fracFix = dtFixCurr/dtFixInit;

accrualFloat = liborFixing * dcFloat * (1 - fracFloat);
accrualFix = nextFix * dcFix * (1 - fracFix);

if abs(accrualFloat) < abs(accrualFloatPrev)
   accrualFloatPrev = 0; 
end

if abs(accrualFix) < abs(accrualFixPrev)
   accrualFixPrev = 0; 
end

if ropFix == 'p'
    accrualCumulative =  - (accrualFloat - accrualFloatPrev) + (accrualFix - accrualFixPrev);% 
elseif ropFix == 'r'
    accrualCumulative = -(accrualFix - accrualFixPrev) + accrualFloat - accrualFloatPrev;     
end
            
accrualFloatPrev = accrualFloat;
accrualFixPrev = accrualFix;


end