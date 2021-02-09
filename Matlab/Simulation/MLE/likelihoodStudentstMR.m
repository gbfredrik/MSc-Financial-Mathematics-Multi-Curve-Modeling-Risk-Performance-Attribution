function [l] = likelihoodStudentstMR(x, r, dt, varianceFunction)
% Computes the log likelihood value for each realization of 
% a Student's t-distributed variable

kappa = x(1);
xiHat = mean(r);

df = x(end);

v = varianceFunction(x, r, dt);

l =  gammaln((df+1)/2) - gammaln(df/2) - 1/2*log(pi) - 1/2*log(df) - ...
     log(sqrt(v(2:end-1)*dt)) - ...
     (df+1)/2 * log(1 + 1/df * ((r(2:end) - kappa .* (xiHat - r(1:end-1))*dt).^2) ./ (v(2:end-1)*dt));


