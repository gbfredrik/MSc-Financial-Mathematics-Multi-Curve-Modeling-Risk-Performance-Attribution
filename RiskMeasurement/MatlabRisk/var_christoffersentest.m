function [reject] = var_christoffersentest(VaRs, PnLs, alpha)
%VAR_CHRISTOFFERSENTEST Summary of this function goes here
%   Detailed explanation goes here

[u00, u01, u10, u11] = numberofperiods(VaRs, PnLs);

pi = (u01 + u11) / (u00 + u01 + u10 + u11);
pi01 = (u01) / (u00 + u01);
pi11 = (u11) / (u10 + u11);

c_stat = teststat_christoffersen(pi, pi01, pi11, u00, u01, u10, u11);

reject = (c_stat > chi2inv(1 - alpha, 1));
end

function [u00, u01, u10, u11] = numberofperiods(VaRs, PnLs)
u00 = 0;
u01 = 0;
u10 = 0;
u11 = 0;

indicator = var_breaches(VaRs, PnLs);

for i = 1:length(VaRs) - 1
    exceed_current = (indicator(i) == 1);
    exceed_next = (indicator(i+1) == 1);
    
    if (~exceed_current && ~exceed_next)
        u00 = u00 + 1;
    elseif (~exceed_current && exceed_next)
        u01 = u01 + 1;
    elseif (exceed_current && ~exceed_next)
        u10 = u10 + 1;
    elseif (exceed_current && exceed_next)
        u11 = u11 + 1;
    end
end

end

function [c_stat] = teststat_christoffersen(pi, pi01, pi11, u00, u01, u10, u11)
c_stat = -2.0 * ((u00 + u10) * log(1.0 - pi) + (u01 + u11) * log(pi)) ...
    + 2.0 * (u00 * log(1.0 - pi01) + u01 * log(pi01) + u10 * log(1.0 - pi11) + u11 * log(pi11));
end

% Helper
function [indicator] = var_breaches(VaRs, PnLs)
indicator = -VaRs > PnLs;
end
