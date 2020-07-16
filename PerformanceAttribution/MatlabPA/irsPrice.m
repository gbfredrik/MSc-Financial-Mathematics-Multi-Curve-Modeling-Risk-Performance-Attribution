function [P] = irsPrice(N, y, ttValueDate, floatCashFlowsUnknown, fixCashFlows, deltaTj, r, pi, floatcf, floatCashFlowsKnown, ropFix)

numFix = length(fixCashFlows);
numUnknownFloat = length(floatCashFlowsUnknown);
numKnownFloat = length(floatCashFlowsKnown);

% Calc value of fix leg
fix = 0;
for i = 1:numFix
    fix = fix + N * y *  exp(ttValueDate/365 * r(ttValueDate+1) - fixCashFlows(i)/365 * r(fixCashFlows(i)+1));
end
%deltaTj(i) *

% Calc value of float leg
float = 0;
for i = 1:numKnownFloat
    float = float + floatcf(i) * exp(ttValueDate/365 * r(ttValueDate+1) - floatCashFlowsKnown(i)/365 * r(floatCashFlowsKnown(i)+1));
end
for i = 1:numUnknownFloat-1
    float = float + N * (exp(floatCashFlowsUnknown(i+1)/365 * pi(floatCashFlowsUnknown(i+1) + 1) - floatCashFlowsUnknown(i)/365 * (r(floatCashFlowsUnknown(i) + 1) + pi(floatCashFlowsUnknown(i) + 1)) ...
        + ttValueDate/365 * r(ttValueDate + 1)) ...
        - exp(ttValueDate/365 * r(ttValueDate + 1) - floatCashFlowsUnknown(i+1)/365 * r(floatCashFlowsUnknown(i+1) + 1)))
end

if ropFix == 'r'
    P = float - fix;
elseif ropFix == 'p'
    P = float - fix;
end

end