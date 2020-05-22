function [sigma] = GJR_GARCH(omega, alpha, gamma, beta, E, fHist)

dXi = E' * (fHist(2,:) - fHist(1,:))';
for i = 1:size(dXi)
    if dXi(i) < 0
        ind(i) = 1;
    else
        ind(i) = 0;
    end
end
ind = ind';
    

sigmaPrev = omega + alpha .* dXi.^2 + gamma .* dXi.^2 .* ind + beta .* dXi.^2;

for i = 3:size(fHist,1)
    dXi = E' * (fHist(i,:) - fHist(i - 1,:))';
    sigmasq = omega + alpha .* dXi.^2 + gamma .* dXi.^2 .* ind + beta .* sigmaPrev;
    sigmaPrev = sigmasq;
end

sigma = sqrt(sigmasq);



end