function portfolioValue = portfolioValueFxSwap(forwardPrices, K, N, fTerm, maturityT, maturityContracts)    
    expVec = length(forwardPrices);
    %portfolioValue = 0;
    for i = 1:length(forwardPrices)
        if isequal(num2str(maturityContracts(i)), 'NaN')
            error('Input maturity not a member of term structure maturities')
        elseif ismembertol(maturityContracts(i), maturityT, 10^(-4))
            ind = find(maturityT == maturityContracts(i)); 
        else
            error('Input maturity not a member of term structure maturities')
        end
        fTermTmp = fTerm(1:ind-1);
        T0 = maturityT(1:ind-1);
        T = maturityT(2:ind);
        deltaT = T-T0;
        expVec(i) = exp(-fTermTmp * deltaT');
        %sumExp = sumExp + exp(-fTermTmp * deltaT');
        %portfolioValue =  portfolioValue + N(i) * (forwardPrices(i) - K(i)) * exp(-fTermTmp * deltaT');
    end
    portfolioValue = N * ((forwardPrices - K) .* expVec)';
end