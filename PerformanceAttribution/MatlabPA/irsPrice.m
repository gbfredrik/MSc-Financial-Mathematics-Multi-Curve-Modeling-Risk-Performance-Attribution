function [P] = irsPrice(N, y, startdate, floatCashFlowsUnknown, fixCashFlows, deltaTj, r, pi, floatcf, floatCashFlowsKnown)

numFix = length(fixCashFlows);
numUnknownFloat = length(floatCashFlowsUnknown);
numKnownFloat = length(floatCashFlowsKnown);

% Calc value of fix leg
fix = 0;
for i = 1:numFix
    fix = fix + N * y * deltaTj(i) * exp(startdate/365 * r(startdate+1) - fixCashFlows(i)/365 * r(fixCashFlows(i)+1));
end


% Calc value of float leg
float = 0;
for i = 1:numKnownFloat
    float = float + floatcf(i) * exp(startdate/365 * r(startdate+1) - floatCashFlowsKnown(i)/365 * r(floatCashFlowsKnown(i)+1));
end
for i = 1:numUnknownFloat-1
    float = float + N * (exp(floatCashFlowsUnknown(i+1)/365 * pi(floatCashFlowsUnknown(i+1) + 1) - floatCashFlowsUnknown(i)/365 * (r(floatCashFlowsUnknown(i) + 1) + pi(floatCashFlowsUnknown(i) + 1)) ...
        + startdate/365 * r(startdate + 1)) ...
        - exp(startdate/365 * r(startdate + 1) - floatCashFlowsUnknown(i+1)/365 * r(floatCashFlowsUnknown(i+1) + 1)))
end

P = float - fix;

end