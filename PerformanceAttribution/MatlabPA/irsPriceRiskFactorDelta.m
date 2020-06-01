function [P] = irsPriceRiskFactorDelta(N, y, floatCashFlows, fixCashFlows, A, E, deltaTj, r, pi, dXi)

numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);

n = size(E.Zero, 1);
k = size(E.Zero, 2);
aZero = [A*E.Zero, zeros(n, k)]; 
aTau = [zeros(n, k), A*E.Tau];

dXiZero = [dXi(1:k); zeros(k, 1)];
dXiTau = [zeros(k, 1); dXi(k + 1:2*k)];


fix = 0;
% Calc fix leg
for i = 1:numCashFix
    fix = fix + deltaTj(i) * exp(floatCashFlows(1)/365 * (r(floatCashFlows(1) + 1) + aZero(floatCashFlows(1) + 1,:) * dXiZero) ...
        - fixCashFlows(i)/365 * (r(fixCashFlows(i) + 1) + aZero(fixCashFlows(i) + 1,:) * dXiZero))';
end
fix = N * y * fix;

% Calc float leg
float = 0;
for i = 1:numCashFloat - 1
    float = float + exp(floatCashFlows(i + 1)/365 * (pi(floatCashFlows(i + 1) + 1) + aTau(floatCashFlows(i + 1) + 1,:) * dXiTau) ...
        - floatCashFlows(i)/365 * ((r(floatCashFlows(i) + 1) + aZero(floatCashFlows(i) + 1,:) * dXiZero) + (pi(floatCashFlows(i) + 1) + aTau(floatCashFlows(i) + 1,:) * dXiTau)) ...
        + floatCashFlows(1)/365 * (r(floatCashFlows(1) + 1) + aZero(floatCashFlows(1) + 1,:) * dXiZero)) ...
        - exp(floatCashFlows(1)/365 * (r(floatCashFlows(1) + 1) + aZero(floatCashFlows(1) + 1,:) * dXiZero) ...
        - floatCashFlows(i + 1)/365 * (r(floatCashFlows(i + 1) + 1) + aZero(floatCashFlows(i + 1) + 1,:) * dXiZero))';
end

float = N * float;

P = float - fix;

end