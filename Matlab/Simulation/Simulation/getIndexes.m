function [floatDatesIndexes, fixDatesIndexes] = getIndexes(simDate, fixDates, floatDates)
        
        
    fixDatesIndexes = zeros(length(fixDates),1);
    floatDatesIndexes = zeros(length(floatDates),1);
    
    for i = 1:length(fixDates)
       fixDatesIndexes(i) = daysdif(simDate, fixDates(i));
    end
    
    for i = 1:length(floatDates)
       floatDatesIndexes(i) = daysdif(simDate, floatDates(i));
    end    
    
end

