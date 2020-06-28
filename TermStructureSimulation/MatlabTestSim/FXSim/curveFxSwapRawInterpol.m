function curveRI = curveFxSwapRawInterpol(maturitiesRIC, basisRIC, curveLength)
    
    T0 = maturitiesRIC(1:end-1);
    T = maturitiesRIC(2:end);
    deltaT = T-T0;
    basis0 = basisRIC(1:end-1);
    basis = basisRIC(2:end);
    deltaBasis = basis-basis0; 
    basisRI = deltaBasis./deltaT;
    %curveLength = round(sum(deltaT)*365);
    pos = 1;
    curveRI = ones(1,curveLength)*inf;
    for nMaturity = 1:length(basisRI)
        for nDays = 1:round(deltaT(nMaturity)*365)
            curveRI(pos) = basisRI(nMaturity);
            pos = pos + 1;
        end
    end
    
end