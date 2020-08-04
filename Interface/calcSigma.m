function [sigma] = calcSigma(data, M, meanDiff)

    sigmaSq = (1/(M-1)) * sum((data(1:end) - meanDiff).^2);
    
    sigma = sqrt(sigmaSq);

end