function [g] = grad(r, pi, fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi)


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
        fix = fix + N * y * dtFix(i) * (-timeFracFix(i+1) * aZero(fixDatesIndexes(i+1),:))'...
            * exp(-timeFracFix(i+1) * r(fixDatesIndexes(i+1)))'; 
    end

    % Calc value of float leg
    float = 0;
    for i = 1:numKnown
        float = float + fixing(i) * dtFloat(i) * N * (-timeFracFloatKnown(i+1) * aZero(floatDatesKnownIndexes(i+1),:))'...
            * exp(-timeFracFloatKnown(i+1) * r(floatDatesKnownIndexes(i+1)))';
    end
    
    for i = numKnown+1:numFloat-1
        float = float + N * ((timeFracFloat(i + 1) * aPi(floatDatesIndexes(i + 1),:) - timeFracFloat(i) * (aZero(floatDatesIndexes(i),:) + aPi(floatDatesIndexes(i),:)))' ...
        * exp(timeFracFloat(i+1) * pi(floatDatesIndexes(i + 1)) - timeFracFloat(i) * (r(floatDatesIndexes(i)) + pi(floatDatesIndexes(i))))'...
        - (-timeFracFloat(i+1) * aZero(floatDatesIndexes(i+1),:))' ...
        * exp(-timeFracFloat(i + 1) * r(floatDatesIndexes(i + 1)))');
    end
    
    if RoP == 'r'
        g = fix - float;
    elseif RoP == 'p'
        g = float - fix;
    end

end