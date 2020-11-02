function [sigma] = GJR_GARCH_d(omega, alpha, beta, E, fPrev, fPrevPrev, sigmaPrev)

dXi = E' * (fPrev - fPrevPrev)';
   
  
sigmaSq = omega + alpha .* dXi.^2 + beta .* sigmaPrev.^2;
sigma = sqrt(sigmaSq);

end