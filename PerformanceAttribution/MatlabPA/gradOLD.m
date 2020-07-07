function [g] = gradOLD(N, y, floatCashFlows, fixCashFlows, deltaTj, aZero, aPi, r, pi)


numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);

fix = 0;
% Calc fix leg
for i = 1:numCashFix
    fix = fix + deltaTj(i) * (floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:))' ...
        * exp(floatCashFlows(1)/365 * r(floatCashFlows(1) + 1) - fixCashFlows(i)/365 * r(fixCashFlows(i) + 1))';
end
fix = N * y * fix;

% Calc float leg
float = 0;
for i = 1:numCashFloat - 1
    float = float + (floatCashFlows(i + 1)/365 * aPi(floatCashFlows(i + 1) + 1,:) - floatCashFlows(i)/365 * (aZero(floatCashFlows(i) + 1,:) + aPi(floatCashFlows(i) + 1,:)) ...
        + floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:))' ...
        * exp(floatCashFlows(i + 1)/365 * pi(floatCashFlows(i + 1) + 1) - floatCashFlows(i)/365 * (r(floatCashFlows(i) + 1) + pi(floatCashFlows(i) + 1)) ...
        + floatCashFlows(1)/365 * r(floatCashFlows(1) + 1))' ...
        - (floatCashFlows(1)/365 * aZero(floatCashFlows(1) + 1,:) - floatCashFlows(i + 1)/365 * aZero(floatCashFlows(i + 1) + 1,:))' ...
        * exp(floatCashFlows(1)/365 * r(floatCashFlows(1) + 1) - floatCashFlows(i + 1)/365 * r(floatCashFlows(i + 1) + 1))';
end

float = N * float;

g = float - fix;

end