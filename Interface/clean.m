%% Read data
dirtyData = xlsread('layout_v0.9.xlsm', 'bid', "C:G");
dirtyDataAdditional = xlsread('layout_v0.9.xlsm', 'bidDirty', "C:G");
cleanPeriodDates = table2array(readtable('layout_v0.9.xlsm', 'Sheet', 'bid', 'Range', "B3:B1831"));
%cleanPeriodDates = cleanPeriodDates(2:end);
additionalPeriodDates = table2array(readtable('layout_v0.9.xlsm', 'Sheet', 'bidDirty', 'Range', "B3:B2352"));
%additionalPeriodDates = additionalPeriodDates(2:end);
%% Get candidates
diffs = dirtyData(2:end,:) - dirtyData(1:end-1,:);
diffsAdditional = dirtyDataAdditional(2:end,:) - dirtyDataAdditional(1:end-1,:);
diffs(isnan(diffs))=0;
diffsAdditional(isnan(diffsAdditional))=0;


%%

for i = 1:size(diffs, 2)
    currInstr = diffs(:,i);
    count = 1;
    for j = 1:size(diffs,1)
        if currInstr(j,1) ~= 0
            candidates(count,1) = currInstr(j,1);
            row(count,1) = j;
            candidateDates(count, 1) = cleanPeriodDates(j, 1); 
            count = count + 1;
        end
    end  
    
    count = 1;
    currInstr = diffsAdditional(:,i);
    for j = 1:size(diffsAdditional, 1)
        candidatesAdditional(count,1) = currInstr(j,1);
        rowAdditional(count,1) = j;
        candidateDatesAdditional(count, 1) = additionalPeriodDates(j, 1);
        count = count + 1;
        
    end
    
    
    meanDiff = mean(candidates);
    countTest = 0;
    for j = 1:size(candidates, 1) - 1
        [startIndex, endIndex] = getWindow(row, cleanPeriodDates, candidateDatesAdditional, j);
        M = size(candidatesAdditional(startIndex:endIndex, 1), 1);
        sigma = calcSigma(candidatesAdditional(startIndex:endIndex, 1), M, meanDiff);
       
        if candidates(j) > meanDiff
            prob1 = normcdf(candidates(j), meanDiff, sigma, 'upper');
        else
            prob1 = normcdf(candidates(j), meanDiff, sigma);
        end
        
        if candidates(j+1) > meanDiff
            prob2 = normcdf(candidates(j+1), meanDiff, sigma, 'upper');
        else
            prob2 = normcdf(candidates(j+1), meanDiff, sigma);
        end
                
        if prob1 < 10^(-4) && prob2 < 10^(-4) && ...
                sign(candidates(j)) ~= sign(candidates(j+1))
           countTest = countTest + 1
           
           removeVec(i, j) = row(j);
           
        end
        
    end
    
end




%% Write Data

%xlswrite('layout_v0.9.xlsm', dirtyData, 'bidCleaned', "C3:G500")