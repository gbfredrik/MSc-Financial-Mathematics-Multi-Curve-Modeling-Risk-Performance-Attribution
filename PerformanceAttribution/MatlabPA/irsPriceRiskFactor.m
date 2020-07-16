function [P] = irsPriceRiskFactor(N, y, ttValueDate, floatCashFlowsUnknown, fixCashFlows, ...
    deltaTj, aZero, aPi, XiZero, XiPi, floatcf, floatCashFlowsKnown, ropFix)

numFix = length(fixCashFlows);
numUnknownFloat = length(floatCashFlowsUnknown);
numKnownFloat = length(floatCashFlowsKnown);

%* deltaTj(i)
%Calc fix leg
fix = 0;
for i = 1:numFix
    fix = fix + N * y  * exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:) * XiZero)';
end


%Calc float leg
float = 0;
for i = 1:numKnownFloat
    float = float + floatcf(i) * exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - floatCashFlowsKnown(i)/365 * aZero(floatCashFlowsKnown(i) + 1,:) * XiZero)';
end

for i = 1:numUnknownFloat - 1
    float = float + N * (exp(floatCashFlowsUnknown(i + 1)/365 * (aPi(floatCashFlowsUnknown(i + 1) + 1,:) * XiPi) ...
        - floatCashFlowsUnknown(i)/365 * (aZero(floatCashFlowsUnknown(i) + 1,:) * XiZero + aPi(floatCashFlowsUnknown(i) + 1,:) * XiPi)...
        + ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero)' ...
        - exp(ttValueDate/365 * aZero(ttValueDate + 1,:) * XiZero ...
        - floatCashFlowsUnknown(i + 1)/365 * aZero(floatCashFlowsUnknown(i + 1) + 1,:) * XiZero)');
end

if ropFix == 'r'
    P = float - fix;
elseif ropFix == 'p'
    P = fix - float;
end

end