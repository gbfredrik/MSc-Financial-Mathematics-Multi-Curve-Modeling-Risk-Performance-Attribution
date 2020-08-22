function [startIndex, endIndex] = getWindow(row, cleanPeriodDates, candidateDatesAdditional, j)

    
    currDate = cleanPeriodDates(row(j));
    
    startIndex = 0;
    endIndex = 0;
    if j == 1158
        'japp'
    end
    
    %for k = 1:size(candidateDatesAdditional, 1)
    %for k = row(j):row(j)+10
        %if (tempDate == candidateDatesAdditional(k) || (tempDate > candidateDatesAdditional(k) && tempDate < candidateDatesAdditional(k+1)))
            %startIndex = k;
            %break
        %end
    %end
   
    for i = 1:370
        tempDate = currDate - 365 + i;
        startIndex = datefind(tempDate, candidateDatesAdditional, 5);
        if startIndex ~= 0
            break
        end
    end

    for i = 1:370
        tempDate = currDate + 365 - i;
        endIndex = datefind(tempDate, candidateDatesAdditional, 5);
        if endIndex ~= 0
            break
        end
    end
    
    if (endIndex == 0)
        'hehe2' 
    end    
    
end