function [P] = irsPrice(r, pi, fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown)
    
    %Init parameters
    numUnknown = length(timeFracFloat) - numKnown;
    numFix = length(timeFracFix); 
    numFloat = numUnknown + numKnown;
    
    % Split float cash flows into known and unknown
    timeFracFloatKnown = timeFracFloat(1:numKnown+1);
    floatDatesKnownIndexes = floatDatesIndexes(1:numKnown+1);

 
    %timeFracFix(1) * r(fixDatesIndexes(1))
    % Calc value of fix leg
    fix = 0;
    for i = 1:numFix-1
        fix = fix + N * y * dtFix(i) * exp(-timeFracFix(i+1) * r(fixDatesIndexes(i+1))); 
    end

    % Calc value of float leg
    float = 0;
    for i = 1:numKnown
        float = float + fixing(i) * dtFloat(i) * N * exp(-timeFracFloatKnown(i+1) * r(floatDatesKnownIndexes(i+1)));
    end
    %timeFracFloatKnown(1) * r(floatDatesKnownIndexes(1))
    for i = numKnown+1:numFloat-1
        float = float + N * (exp(timeFracFloat(i+1) * pi(floatDatesIndexes(i+1)) - timeFracFloat(i) * (r(floatDatesIndexes(i)) + pi(floatDatesIndexes(i)))) - exp(-timeFracFloat(i+1) * r(floatDatesIndexes(i+1))));
    end
    
    %for i = numKnown+1:numFloat-1
    %   float = float + N * (exp(timeFracFloat(i+1) * (r(floatDatesIndexes(i+1))+pi(floatDatesIndexes(i+1))) - timeFracFloat(i) * (r(floatDatesIndexes(i))+pi(floatDatesIndexes(i))))-1) ...
    %       *exp(-timeFracFloat(i+1) * (r(floatDatesIndexes(i+1)))); 
    %end
    
%timeFracFloat(1) * r(floatDatesIndexes(1))
    % Set value based on receive or pay
    if char(RoP) == 'r'
        P = fix - float;
    elseif char(RoP) == 'p'
        P = float - fix;
    end
end
