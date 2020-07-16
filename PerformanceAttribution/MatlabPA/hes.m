function [H] = hes(N, y, startdate, floatCashFlowsUnknown, fixCashFlows, deltaTj, aZero, aPi, r, pi, floatcf, floatCashFlowsKnown, ropFix)

numFix = length(fixCashFlows);
numUnknownFloat = length(floatCashFlowsUnknown);
numKnownFloat = length(floatCashFlowsKnown);

fix = 0;
% Calc fix leg
for i = 1:numFix
    fix = fix + deltaTj(i) * (startdate/365 * aZero(startdate + 1,:) - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:))' ...
        *(startdate/365 * aZero(startdate + 1,:) - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:)) ...
        * exp(startdate/365 * r(startdate + 1) - fixCashFlows(i)/365 * r(fixCashFlows(i) + 1))';
end
fix = N * y * fix;

% Calc float leg
float = 0;
for i = 1:numKnownFloat
    float = float + floatcf(i) * (startdate/365 * aZero(startdate + 1,:) - floatCashFlowsKnown(i)/365 * aZero(floatCashFlowsKnown(i) + 1,:))' ...
        *(startdate/365 * aZero(startdate + 1,:) - floatCashFlowsKnown(i)/365 * aZero(floatCashFlowsKnown(i) + 1,:)) ...
        * exp(startdate/365 * r(startdate + 1) - floatCashFlowsKnown(i)/365 * r(floatCashFlowsKnown(i) + 1))';
end
for i = 1:numUnknownFloat - 1
    float = float + N*((floatCashFlowsUnknown(i + 1)/365 * aPi(floatCashFlowsUnknown(i + 1) + 1,:) - floatCashFlowsUnknown(i)/365 * (aZero(floatCashFlowsUnknown(i) + 1,:) + aPi(floatCashFlowsUnknown(i) + 1,:))...
        + startdate/365 * aZero(startdate + 1,:))' ...
        * (floatCashFlowsUnknown(i + 1)/365 * aPi(floatCashFlowsUnknown(i + 1) + 1,:) - floatCashFlowsUnknown(i)/365 * (aZero(floatCashFlowsUnknown(i) + 1,:) + aPi(floatCashFlowsUnknown(i) + 1,:))...
        + startdate/365 * aZero(startdate + 1,:)) ...
        * exp(floatCashFlowsUnknown(i + 1)/365 * pi(floatCashFlowsUnknown(i + 1) + 1) - floatCashFlowsUnknown(i)/365 * (r(floatCashFlowsUnknown(i) + 1) + pi(floatCashFlowsUnknown(i) + 1)) ...
        + startdate/365 * r(startdate + 1))' ...
        - (startdate/365 * aZero(startdate + 1,:) - floatCashFlowsUnknown(i + 1)/365 * aZero(floatCashFlowsUnknown(i + 1) + 1,:))' ...
        * (startdate/365 * aZero(startdate + 1,:) - floatCashFlowsUnknown(i + 1)/365 * aZero(floatCashFlowsUnknown(i + 1) + 1,:)) ...
        * exp(startdate/365 * r(startdate + 1) - floatCashFlowsUnknown(i + 1)/365 * r(floatCashFlowsUnknown(i + 1) + 1))');
end

if ropFix == 'r'
    H = fix - float;
elseif ropFix == 'p'
    H = float - fix;
end
end