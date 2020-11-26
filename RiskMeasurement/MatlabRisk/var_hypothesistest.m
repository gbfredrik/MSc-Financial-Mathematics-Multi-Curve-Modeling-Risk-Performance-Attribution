function [reject, Z, tilde_m] = var_hypothesistest(VaRs, PnLs, c, alpha)
%VAR_HYPOTHESISTEST Summary of this function goes here
%   Detailed explanation goes here

X_T = teststatistic_XT(VaRs, PnLs);
T = length(VaRs);
p = 1 - c;
Z = teststatistic_Z(X_T, T, p);
tilde_m = norminv(1 - alpha/2); % Defaults mu = 0, sigma = 1

reject = abs(Z) >= tilde_m;
end

function [X_T] = teststatistic_XT(VaRs, PnLs)
X_T = sum(var_breaches(VaRs, PnLs));
end

function [Z] = teststatistic_Z(X_T, T, p)
Z = (X_T - T * p) / sqrt(T * p * (1 - p));
end

% Helper
function [indicator] = var_breaches(VaRs, PnLs)
indicator = -VaRs > PnLs;
end
