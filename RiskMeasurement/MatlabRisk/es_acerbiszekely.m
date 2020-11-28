function [reject, p] = es_acerbiszekely(X, VaRs, ESs, PnLs, phi)
%ES_ACERBISZEKELY Summary of this function goes here
%   Detailed explanation goes here

X(end) = [];
M = size(X, 2);
Z = zeros(M,1);
Z_X_actual = teststat_Z1(VaRs, ESs, PnLs);
n = 0;

for i = 1:M
    Z(i) = teststat_Z1(VaRs, ESs, X{i});
    
    if (Z(i) < Z_X_actual)
        n = n + 1;
    end
end

p = n / M;

reject = (p < phi);
end

function [Z] = teststat_Z1(VaRs, ESs, PnLs)
indicator = var_breaches(VaRs, PnLs);
N_T = sum(indicator);

if (N_T == 0)
    Z = 0;
else
    Z = sum((PnLs .* indicator) ./ ESs) / N_T + 1;
end
end

% Helper
function [indicator] = var_breaches(VaRs, PnLs)
indicator = -VaRs > PnLs;
end
