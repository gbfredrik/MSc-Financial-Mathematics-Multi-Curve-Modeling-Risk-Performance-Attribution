function curveOIS = curveOISRawInterpol(maturitiesRIC, priceRIC, curveLength)

    indEnd = find(isnan(maturitiesRIC), 1, 'first')-1;
       
    T0 = maturitiesRIC(1:indEnd-1);
    T = maturitiesRIC(2:indEnd);
    deltaT = T-T0;
    %curveLength = round(sum(deltaT)*365);
    pos = 1;
    curveOIS = ones(1, curveLength)*inf;
    for nMaturity = 1:indEnd -1
        for nDays = 1:round(deltaT(nMaturity)*365)
            curveOIS(pos) = priceRIC(nMaturity);
            pos = pos + 1;
        end
    end 
end