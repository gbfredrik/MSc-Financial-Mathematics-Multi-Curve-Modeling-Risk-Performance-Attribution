function portfolioValue = portfolioValueFxSwap(forwardPrices, K, N, fTerm, discPoints, maturityContracts)    
    maturityContracts = maturityContracts(~isnan(maturityContracts));
 
    K = K(~isnan(forwardPrices));
    forwardPrices = forwardPrices(~isnan(forwardPrices));

    
    expVec = length(forwardPrices);
    %portfolioValue = 0;
    for i = 1:length(forwardPrices)
        if isequal(num2str(maturityContracts(i)), 'NaN')
            error('Input maturity not a member of term structure maturities')
        elseif ismembertol(maturityContracts(i), discPoints, 10^(-4))
            ind = find(abs(discPoints-maturityContracts(i)) < 0.00000000001); 
        else
            error('Input maturity not a member of term structure maturities')
        end
        fTermTmp = fTerm(1:ind-1);
        T0 = discPoints(1:ind-1);
        T = discPoints(2:ind);
        deltaT = T-T0;
        expVec(i) = exp(-fTermTmp * deltaT');
        %sumExp = sumExp + exp(-fTermTmp * deltaT');
        %portfolioValue =  portfolioValue + N(i) * (forwardPrices(i) - K(i)) * exp(-fTermTmp * deltaT');
    end
    portfolioValue = N * sum(((forwardPrices - K) .* expVec));
end