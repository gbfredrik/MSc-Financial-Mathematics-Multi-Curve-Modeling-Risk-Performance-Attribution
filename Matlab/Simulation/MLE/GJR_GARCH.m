function [sigma] = GJR_GARCH(omega, alpha, beta, E, fHist)

dXi = E' * (fHist(2,:) - fHist(1,:))';



sigmaPrev = omega + alpha .* dXi.^2 + beta .* dXi.^2;

for i = 3:size(fHist,1)
    dXi = E' * (fHist(i,:) - fHist(i - 1,:))';
    sigmasq = omega + alpha .* dXi.^2 + beta .* sigmaPrev;
    sigmaPrev = sigmasq;
end

sigma = sqrt(sigmasq);



end