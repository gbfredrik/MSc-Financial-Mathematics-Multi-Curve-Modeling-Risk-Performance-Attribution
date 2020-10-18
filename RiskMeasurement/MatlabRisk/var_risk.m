function [var, ind] = var_risk(sims, conf)

sorted_outcomes = sort(-sims);

ind = floor(length(sorted_outcomes) * conf) + 1;
var = sorted_outcomes(ind);
end
