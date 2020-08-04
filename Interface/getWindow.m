function [startIndex, endIndex] = getWindow(row, cleanPeriodDates, candidateDatesAdditional, j)

    
    currDate = cleanPeriodDates(row(j));
    
    startIndex = 0;
    endIndex = 0;
    
    for i = 1:30
       tempDate = currDate - 365 + i;
       for j = 1:size(candidateDatesAdditional, 1)
           if (tempDate == candidateDatesAdditional(j))
               startIndex = j;
               break
           end
       end
       if (startIndex ~= 0)
           break
       end
    end
  
    
    for i = 1:30
       tempDate = currDate + 365 - i;
       for j = 1:size(candidateDatesAdditional, 1)
           if (tempDate == candidateDatesAdditional(j))
                endIndex = j;
                break
           end
       end
       if (endIndex ~= 0)
           break
       end
    end
    
    
end