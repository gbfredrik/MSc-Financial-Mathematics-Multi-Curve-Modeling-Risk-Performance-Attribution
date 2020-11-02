function floatFixing = getIbor(currDate, Ibor, IborDates, prevNumKnown, numKnown, fixingDates, floatFixing)

    % Cases for handling the Ibor fixings
    % Set the first fixing
    if prevNumKnown == 0 && currDate == fixingDates(1)
        for j = 1:length(IborDates)
            if IborDates(j) == currDate
                iborIndex = j;
            end
        end
        floatFixing = Ibor(iborIndex);
    % Set a fixing for cases when two are known
    elseif numKnown == 2 && currDate == fixingDates(1)
        for j = 1:length(IborDates)
            if IborDates(j) == currDate
                iborIndex = j;
            end
        end            
        floatFixing = [floatFixing; Ibor(iborIndex)];
    % Decrease number of known fixings when previous == 2 and we have passed
    % the end date of a floating period
    elseif prevNumKnown == 2 && numKnown == 1 && currDate == floatDates(1)
        floatFixing = floatFixing(2);
    end

end