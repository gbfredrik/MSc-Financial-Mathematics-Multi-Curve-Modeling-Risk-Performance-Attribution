function [valueParams] = setActiveStatus(valueParams, currDate)
    
simDate = valueParams{1};
floatDates = valueParams{5};
numContracts = length(floatDates);
activeStatus = valueParams{12};
fixingDates = valueParams{7};

% Set active contracts
    for j = 1:numContracts
        if activeStatus(j) == 1 && floatDates{j}(end) < simDate
            activeStatus(j) = 0;
        elseif length(fixingDates{j}) > 0
            if currDate == fixingDates{j}(1) && activeStatus(j) == 0
               activeStatus(j) = 1;
            end
        end
    end

    valueParams{12} = activeStatus;
end

