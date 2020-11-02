function [valueParams] = getParameters(valueParams, f, pi)  
    
    simDate = valueParams{1};        
    Ibor = valueParams{3};
    IborDates = valueParams{4};
    numContracts = length(valueParams{5});
    activeStatus = valueParams{12};
    cash = 0;
    fixingDay = 0;
    for i = 1:numContracts
        
        if activeStatus(i) == 1
            floatLegDCC = valueParams{2}(i);
            floatDates = valueParams{5}{i};
            fixDates = valueParams{6}{i};
            fixingDates = valueParams{7}{i};
            prevNumKnown = valueParams{8}(i);
            dtFloat = valueParams{10}{i};
            floatFixing = valueParams{11}{i};
            N = valueParams{14}(i);
            RoP = valueParams{13}(i);
            dtFix = valueParams{9}{i};
            y = valueParams{15}(i);
            cash = 0;
            floatNext = floatDates(2);
            fixNext = fixDates(2);
            
            % Set new active dates
            % floatDates
            if simDate == floatDates(2)
                floatDates = [simDate; floatDates(3:end)]; 
            else
                floatDates = [simDate; floatDates(2:end)];
            end

            % fixDates
            if simDate == fixDates(2)
                fixDates = [simDate; fixDates(3:end)];
            else
                fixDates = [simDate; fixDates(2:end)];
            end

            % Set time fracs for all cash flows
            timeFracFloat = zeros(length(floatDates),1);
            for j = 1:length(floatDates)
                timeFracFloat(j) = daysdif(simDate, floatDates(j), 0) / 365;
            end

            % simDate
            timeFracFix = zeros(length(fixDates),1);
            for j = 1:length(fixDates)
                timeFracFix(j) = daysdif(simDate, fixDates(j), 0) / 365; 
            end    

            % Get number of known cash flows
            % Between a fixing date and the period start, except first => 2
            % known
            if length(fixingDates) > 0 && simDate == fixingDates(1) && prevNumKnown ~= 0
                numKnown = 2;    
                for j = 1:length(IborDates)
                    if IborDates(j) == simDate
                        iborIndex = j;
                        break
                    end
                end                                          
                floatFixing = [floatFixing; Ibor(iborIndex)];
                
                
                if floatLegDCC == "act/360"
                    dtFloat = [dtFloat, daysdif(floatDates(2), floatDates(3), 2) / 360];
                elseif floatLegDCC == "30/360"
                    dtFloat = [dtFloat, daysdif(floatDates(2), floatDates(3), 5) / 360];               
                end        
                
                fixingDay = 1;
                
            % Between a started period and the next fixing date => 1 known 
            elseif prevNumKnown == 2 && simDate == floatNext
                numKnown = 1; 
                if RoP == 'p'
                    cash = floatFixing(1) * dtFloat(1) * N;
                elseif RoP == 'r'
                    cash = -floatFixing(1) * dtFloat(1) * N;
                end

                floatFixing = floatFixing(2);
                dtFloat = dtFloat(2);
                
             % Last CF
             elseif prevNumKnown == 1 && simDate == floatNext
                numKnown = 1; 
                if RoP == 'p'
                    cash = floatFixing(1) * dtFloat(1) * N;
                elseif RoP == 'r'
                    cash = -floatFixing(1) * dtFloat(1) * N;
                end
                
            elseif prevNumKnown == 0
                 numKnown = 1;
                 for j = 1:length(IborDates)
                    if IborDates(j) == simDate
                        iborIndex = j;
                        break
                    end
                 end
                
                floatFixing = Ibor(iborIndex);
                                
                if floatLegDCC == "act/360"
                    dtFloat = daysdif(floatDates(2), floatDates(3), 2) / 360;
                elseif floatLegDCC == "30/360"
                    dtFloat = daysdif(floatDates(2), floatDates(3), 5) / 360;                
                end
            else
                numKnown = prevNumKnown;
            end

            if simDate == fixNext
                if RoP == 'p'
                    cash = cash - N * y * dtFix(1);
                elseif RoP == 'r'
                    cash = cash + N * y * dtFix(1);
                end 
                dtFix = dtFix(2:end);
            end

            % fixingDates
            if length(fixingDates) > 0 && simDate == fixingDates(1)
                fixingDates = fixingDates(2:end);
            end 

            valueParams{5}{i} = floatDates;
            valueParams{6}{i} = fixDates;
            valueParams{7}{i} = fixingDates;
            valueParams{8}(i) = numKnown;
            valueParams{10}{i} = dtFloat;
            valueParams{11}{i} = floatFixing;            
            valueParams{16}{i} = timeFracFix;
            valueParams{17}{i} = timeFracFloat; 
            valueParams{19}{i} = cash;
            valueParams{9}{i} = dtFix;
            valueParams{20} = fixingDay;
            
        end
    end
        
    Ibor = valueParams{3}(2:end);
    IborDates = valueParams{4}(2:end);   
     
end