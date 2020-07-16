function [accrualFloat, accrualFix] = accrual(liborFixing, nextFix, dcf, dtFloatCurr, dtFloatInit, dtFixCurr, dtFixInit)

fracFloat = dtFloatCurr/dtFloatInit;
fracFix = dtFixCurr/dtFixInit;

accrualFloat = liborFixing * dcf * (1 - fracFloat);
accrualFix = nextFix * dcf * (1 - fracFix);


end