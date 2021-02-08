function [f] = negLikelihood(x,r,dt,varianceFunction,likelihoodFunction,useMR)
% Computes the log likelihood value with opposite sign (fmincon solves a minimization problem)

l = likelihoodFunction(x, r, dt, varianceFunction);

if (length(l) ~= length(r)-useMR)
  error('The function likelihoodNormal should compute the log likelyhood for each realization in r minus the first')
end

f = -sum(l); % Switch sign (fmincon solves a minimization problem)
