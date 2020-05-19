function [sigma] = GJR_GARCH_d(omega, alpha, gamma, beta, E, fPrev, fPrevPrev, sigmaPrev)

dXi = E' * (fPrev - fPrevPrev)';

for i = 1:size(dXi)
    if dXi(i) < 0
        ind(i) = 1;
    else
        ind(i) = 0;
    end
end
ind = ind';
  
sigmaSq = omega + alpha .* dXi.^2 + gamma .* dXi.^2 .* ind + beta .* sigmaPrev.^2;
sigma = sqrt(sigmaSq);

end