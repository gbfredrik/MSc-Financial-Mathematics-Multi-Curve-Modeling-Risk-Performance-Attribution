function [l]=likelihoodNormal(x,r,dt,varianceFunction)
% Computes the log likelihood value for each realization of a normally distributed variable

mu = x(1);

v = varianceFunction(x,r,dt);

l = -1/2*(log(2*pi*v(1:end-1)*dt)+((r - mu*dt).^2)./(v(1:end-1)*dt));
end

% l = % log likelihood value for a normally distributed variables

