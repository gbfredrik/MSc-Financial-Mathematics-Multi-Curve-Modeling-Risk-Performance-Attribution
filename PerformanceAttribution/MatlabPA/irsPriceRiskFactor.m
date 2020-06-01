function [P] = irsPriceRiskFactor(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, Xi)

numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);

n = size(E.Zero, 1);
kZero = size(E.Zero, 2);
kTau = size(E.Tau, 2);
aZero = [A*E.Zero, zeros(n, kTau)]; 
aTau = [zeros(n, kZero), A*E.Tau];

XiZero = [Xi(1:kZero); zeros(kTau, 1)];
XiTau = [zeros(kZero, 1); Xi(kZero+1 : end)];

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
    float = float + exp(floatCashFlows(i + 1)/365 * (aTau(floatCashFlows(i + 1) + 1,:) * XiTau) ...
        - floatCashFlows(i)/365 * (aZero(floatCashFlows(i) + 1,:) * XiZero + aTau(floatCashFlows(i) + 1,:) * XiTau)...
        + floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) * XiZero)' ...
        - exp(floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) * XiZero ...
        - floatCashFlows(i + 1)/365 * aZero(floatCashFlows(i + 1) + 1,:) * XiZero)';
end

float = N * float;

P = float - fix;

end