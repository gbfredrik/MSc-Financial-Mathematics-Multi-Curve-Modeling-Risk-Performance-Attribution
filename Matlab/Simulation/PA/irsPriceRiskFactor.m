function  [P, nextEstFloat] = irsPriceRiskFactor(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiZero, XiPi)

    %Init parameters
    numUnknown = length(timeFracFloat) - numKnown;
    numFix = length(timeFracFix); 
    numFloat = numUnknown + numKnown;

    % Split float cash flows into known and unknown
    timeFracFloatKnown = timeFracFloat(1:numKnown+1);
    floatDatesKnownIndexes = floatDatesIndexes(1:numKnown+1);    
        
    
    %Calc fix leg
    fix = 0;   
    for i = 1:numFix-1
        fix = fix + N * y * dtFix(i) * exp(-timeFracFix(i+1) * aZero(fixDatesIndexes(i+1),:) * XiZero);
    end
    
    
    % Calc value of float leg
    float = 0;
    for i = 1:numKnown
        float = float + fixing(i) * dtFloat(i) * N * exp(-timeFracFloatKnown(i+1)*aZero(floatDatesKnownIndexes(i+1),:)* XiZero);
    end
  
    
    i = numKnown+1;
    
    if numKnown+1 <= numFloat-1
        nextEstFloat = N * (exp(timeFracFloat(i+1) * (aPi(floatDatesIndexes(i + 1),:) * XiPi) - timeFracFloat(i) * (aZero(floatDatesIndexes(i)+1,:) * XiZero + aPi(floatDatesIndexes(i)+1,:) * XiPi))...
                - exp(-timeFracFloat(i + 1) * aZero(floatDatesIndexes(i + 1),:) * XiZero)');
    else
        nextEstFloat = 0;
    end
        
    for i = numKnown+1:numFloat-1
        float = float + N * (exp(timeFracFloat(i+1) * (aPi(floatDatesIndexes(i + 1),:) * XiPi) - timeFracFloat(i) * (aZero(floatDatesIndexes(i)+1,:) * XiZero + aPi(floatDatesIndexes(i)+1,:) * XiPi))...
            - exp(-timeFracFloat(i + 1) * aZero(floatDatesIndexes(i + 1),:) * XiZero)');
    end
    
    if RoP == 'r'
        P = fix - float;
    elseif RoP == 'p'
        P = float - fix;
    end

end