function price = priceFXswapRawInterpol(spotPrice, prices, maturityContracts, maturity)
    maturityContracts = maturityContracts(~isnan(maturityContracts));
    prices = prices(~isnan(prices));
    
    if isequal(num2str(maturity), 'NaN')
        price = "NaN";
        return
    elseif ismembertol(maturity, maturityContracts, 10^(-6))
        ind = find(maturityContracts == maturity);
        useBasis_m = false;
    elseif maturity > max(maturityContracts)
        error('Input maturity not a member of term structure maturities')
    elseif maturity < min(maturityContracts(2))
        ind = [];
        ind_m = 2;
        useBasis_m = true;
    else
        ind = find(maturity >= maturityContracts, 1, 'last');
        ind_m = ind+1;
        useBasis_m = true;
    end
    basis0 = prices(1:ind-1);
    basis = prices(2:ind);
    deltaBasis = basis-basis0;
    T0 = maturityContracts(1:ind-1);
    T = maturityContracts(2:ind);
    deltaT = T-T0;
    
    basisRI = deltaBasis./deltaT;
    
    if useBasis_m
        %deltaTT_m = maturitiesRIC(ind_m)-maturitiesRIC(ind_m-1)-(maturitiesRIC(ind_m)-maturity);
        deltaTT_m = maturity - maturityContracts(ind_m-1);
        deltaT_m = maturityContracts(ind_m)-maturityContracts(ind_m-1);
        deltaBasis_m = prices(ind_m)-prices(ind_m-1);
        basisRI_m = deltaBasis_m/deltaT_m;
    else
        deltaTT_m = 0;
        basisRI_m = 0;
    end
   
    price = spotPrice + basisRI * deltaT' + basisRI_m * deltaTT_m;
    
    
end
