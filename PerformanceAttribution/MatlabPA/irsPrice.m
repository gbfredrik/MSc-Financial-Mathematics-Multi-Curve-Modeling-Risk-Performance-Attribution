function [P] = irsPrice(N, y, floatCashFlows, fixCashFlows, deltaTj, r, pi)

numCashFix = length(fixCashFlows);
numCashFloat = length(floatCashFlows);


% Fix leg
fix = 0;
for i = 1:numCashFix
   
    fix = fix + deltaTj(i) * exp(floatCashFlows(1)/365 * r(floatCashFlows(1) + 1) - fixCashFlows(i)/365 * r(fixCashFlows(i) + 1));
    
end
fix = N * y * fix;

% Float leg
float = 0;
for i = 1:numCashFloat - 1
    float = float + exp(floatCashFlows(i + 1)/365 * pi(floatCashFlows(i + 1) + 1) - floatCashFlows(i)/365 * (r(floatCashFlows(i) + 1) + pi(floatCashFlows(i) + 1)) ...
        + floatCashFlows(1)/365 * r(floatCashFlows(1) + 1)) ...
        - exp(floatCashFlows(1)/365 * r(floatCashFlows(1) + 1) - floatCashFlows(i + 1)/365 * r(floatCashFlows(i + 1) + 1));
end

float = N * float;

P = float - fix;

end
