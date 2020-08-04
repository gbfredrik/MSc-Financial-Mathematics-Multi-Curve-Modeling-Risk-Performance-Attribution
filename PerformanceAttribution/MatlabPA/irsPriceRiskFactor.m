function [P] = irsPriceRiskFactor(N, y, ttValueDate, floatCashFlows, fixCashFlows, dtFix, aZero, aPi, XiZero, XiPi, ropFix, nextFix, blomvallFixing, daysToNextFix, daysToNextFloat)

numFix = length(fixCashFlows);
numFloat = length(floatCashFlows);


%Calc fix leg
fix = nextFix * exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - daysToNextFix/365 * aZero(daysToNextFix + 1,:) * XiZero);
    
for i = 2:numFix
    fix = fix + N * y * dtFix(i) * exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:) * XiZero)';
end


%Calc float leg
float = blomvallFixing * exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - daysToNextFloat/365 * aZero(daysToNextFloat + 1,:) * XiZero);
    
for i = 1:numFloat - 1
    float = float + N * (exp(floatCashFlows(i + 1)/365 * (aPi(floatCashFlows(i + 1) + 1,:) * XiPi) ...
        - floatCashFlows(i)/365 * (aZero(floatCashFlows(i) + 1,:) * XiZero + aPi(floatCashFlows(i) + 1,:) * XiPi)...
        + ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero)' ...
        - exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - floatCashFlows(i + 1)/365 * aZero(floatCashFlows(i + 1) + 1,:) * XiZero)');
end

if ropFix == 'r'
    P = fix - float;
elseif ropFix == 'p'
    P = float - fix;
end

end