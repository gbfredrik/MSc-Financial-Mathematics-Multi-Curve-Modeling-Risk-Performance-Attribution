function price = priceFXswap(spotPrice, fBase, fTerm, pi, T, maturity)
    if isequal(num2str(maturity), 'NaN')
        price = "NaN";
        return
    elseif ismembertol(maturity, T, 10^(-6))
        ind = find(T == maturity);
    else
        error('Input maturity not a member of term structure maturities')
    end
    fBase = fBase(1:ind-1);
    fTerm = fTerm(1:ind-1);
    pi = pi(1:ind-1);
    T0 = T(1:ind-1);
    T = T(2:ind);
    deltaT = T-T0;
    price = spotPrice*exp(fTerm*deltaT'-fBase*deltaT'+pi*deltaT');
end