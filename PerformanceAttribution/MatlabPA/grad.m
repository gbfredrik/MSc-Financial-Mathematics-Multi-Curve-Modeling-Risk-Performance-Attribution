function [g] = grad(N, y, ttValueDate, fixCashFlows, deltaTj, aZero, aPi, r, pi, floatCashFlows, ropFix, nextFix, blomvallFixing, daysToNextFix, daysToNextFloat)


numFix = length(fixCashFlows);
numFloat = length(floatCashFlows);


fix = nextFix * (ttValueDate/365 * aZero(ttValueDate + 1,:) - daysToNextFix/365 * aZero(daysToNextFix + 1,:))' ...
        * exp(ttValueDate/365 * r(ttValueDate + 1) - daysToNextFix/365 * r(daysToNextFix + 1))';

% Calc fix leg
for i = 1:numFix
    fix = fix + N * y * deltaTj(i) * (ttValueDate/365 * aZero(ttValueDate + 1,:) - fixCashFlows(i)/365 * aZero(fixCashFlows(i) + 1,:))' ...
        * exp(ttValueDate/365 * r(ttValueDate + 1) - fixCashFlows(i)/365 * r(fixCashFlows(i) + 1))';
end


float = blomvallFixing * (ttValueDate/365 * aZero(ttValueDate + 1,:) - daysToNextFloat/365 * aZero(daysToNextFloat + 1,:))' ...
        * exp(ttValueDate/365 * r(ttValueDate + 1) - daysToNextFloat/365 * r(daysToNextFloat + 1))';

for i = 1:numFloat - 1
    float = float + N * ((floatCashFlows(i + 1)/365 * aPi(floatCashFlows(i + 1) + 1,:) - floatCashFlows(i)/365 * (aZero(floatCashFlows(i) + 1,:) + aPi(floatCashFlows(i) + 1,:)) ...
        + ttValueDate/365 * aZero(ttValueDate + 1,:))' ...
        * exp(floatCashFlows(i + 1)/365 * pi(floatCashFlows(i + 1) + 1) - floatCashFlows(i)/365 * (r(floatCashFlows(i) + 1) + pi(floatCashFlows(i) + 1)) ...
        + ttValueDate/365 * r(ttValueDate + 1))' ...
        - (ttValueDate/365 * aZero(ttValueDate + 1,:) - floatCashFlows(i + 1)/365 * aZero(floatCashFlows(i + 1) + 1,:))' ...
        * exp(ttValueDate/365 * r(ttValueDate + 1) - floatCashFlows(i + 1)/365 * r(floatCashFlows(i + 1) + 1))');
end


if ropFix == 'r'
    g = fix - float;
elseif ropFix == 'p'
    g = float - fix;
end

end