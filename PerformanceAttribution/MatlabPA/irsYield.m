function [y] = irsYield(floatCashFlows, fixCashFlows, deltaTj, r)

numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);

rZero = r.rZero;
rTau = r.rTau;

% Fix leg
fix = 0;
for i = 1:numCashFix
   
    fix = fix + deltaTj(i) * exp(floatCashFlows(1)/365 * rZero(floatCashFlows(1) + 1) - fixCashFlows(i)/365 * rZero(fixCashFlows(i) + 1));
    
end


% Float leg
float = 0;
for i = 1:numCashFloat - 1
    float = float + exp(floatCashFlows(i + 1)/365 * rTau(floatCashFlows(i + 1) + 1) - floatCashFlows(i)/365 * rTau(floatCashFlows(i) + 1) ...
        + floatCashFlows(1)/365 * rZero(floatCashFlows(1) + 1) - floatCashFlows(i + 1)/365 * rZero(floatCashFlows(i + 1) + 1)) ...
        - exp(floatCashFlows(1)/365 * rZero(floatCashFlows(1) + 1) - floatCashFlows(i + 1)/365 * rZero(floatCashFlows(i + 1) + 1));
end

y = float / fix;

end