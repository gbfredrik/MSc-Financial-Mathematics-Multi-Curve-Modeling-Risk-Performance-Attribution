function price = priceFXswapRawInterpol(spotPrice, basisRIC, maturitiesRIC, maturity)
    if isequal(num2str(maturity), 'NaN')
        price = "NaN";
        return
    elseif ismembertol(maturity, maturitiesRIC, 10^(-6))
        ind = find(maturitiesRIC == maturity);
        useBasis_m = false;
    elseif maturity > max(maturitiesRIC)
        error('Input maturity not a member of term structure maturities')
    elseif maturity < min(maturitiesRIC(2))
        ind = [];
        ind_m = 2;
        useBasis_m = true;
    else
        ind = find(maturity >= maturitiesRIC, 1, 'last');
        ind_m = ind+1;
        useBasis_m = true;
    end
    basis0 = basisRIC(1:ind-1);
    basis = basisRIC(2:ind);
    deltaBasis = basis-basis0;
    T0 = maturitiesRIC(1:ind-1);
    T = maturitiesRIC(2:ind);
    deltaT = T-T0;
    
    basisRI = deltaBasis./deltaT;
    
    if useBasis_m
        %deltaTT_m = maturitiesRIC(ind_m)-maturitiesRIC(ind_m-1)-(maturitiesRIC(ind_m)-maturity);
        deltaTT_m = maturity - maturitiesRIC(ind_m-1);
        deltaT_m = maturitiesRIC(ind_m)-maturitiesRIC(ind_m-1);
        deltaBasis_m = basisRIC(ind_m)-basisRIC(ind_m-1);
        basisRI_m = deltaBasis_m/deltaT_m;
    else
        deltaTT_m = 0;
        basisRI_m = 0;
    end
   
    price = spotPrice + basisRI * deltaT' + basisRI_m * deltaTT_m;
    
    
end
