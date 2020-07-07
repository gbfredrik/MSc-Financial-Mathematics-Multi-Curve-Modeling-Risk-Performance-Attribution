function [P] = irsPriceRiskFactorOLD(N, y, floatCashFlows, fixCashFlows, deltaTj, aZero, aPi, XiZero, XiPi)

numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);

fix = 0;
% Calc fix leg
for i = 1:numCashFix
    fix = fix + deltaTj(i) * exp(floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) * XiZero ...
        - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:) * XiZero)';
end
fix = N * y * fix;

% Calc float leg
float = 0;
for i = 1:numCashFloat - 1
    float = float + exp(floatCashFlows(i + 1)/365 * (aPi(floatCashFlows(i + 1) + 1,:) * XiPi) ...
        - floatCashFlows(i)/365 * (aZero(floatCashFlows(i) + 1,:) * XiZero + aPi(floatCashFlows(i) + 1,:) * XiPi)...
        + floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) * XiZero)' ...
        - exp(floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) * XiZero ...
        - floatCashFlows(i + 1)/365 * aZero(floatCashFlows(i + 1) + 1,:) * XiZero)';
end

float = N * float;

P = float - fix;

end