function [es] = es_risk(sims, conf)

sorted_outcomes = sort(-sims);

ind = floor(length(sorted_outcomes) * conf) + 1;
es = mean(sorted_outcomes(ind:end));
end
