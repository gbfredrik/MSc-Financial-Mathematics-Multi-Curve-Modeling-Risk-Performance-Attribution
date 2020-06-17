function [isBetter, probability] = likelihoodRatioTest(di, confidenceLevel)
    N = length(di);
    d = mean(di);
    variance = (1/(N-1))*sum((di-d).^2);
    stdDev = sqrt(variance);
    s = stdDev/sqrt(N);
    
    probability = d/s;
    cdfInv = norminv(confidenceLevel);
    
    isBetter = probability > cdfInv;
    
end