function [valueParams, portfolioValues, cash, realizedPnL] = initValueParams(tradeDatesAll_OOS, fixingDates, floatDates, fixDates, ccy, Ibor, IborDates, RoP, Nom, yield, f, pi)

    % Set number of active contracts and split parameter 
    numContracts = length(floatDates);
    numActive = 0;
    activeStatus = zeros(numContracts, 1);
    for i = 1:numContracts
        if daysdif(tradeDatesAll_OOS(1), fixingDates{i}(1)) == 0
            numActive = numActive + 1;
            activeStatus(i) = 1;  
        end
    end

    % Define data set to store results
    portfolioValues = cell(length(tradeDatesAll_OOS), 1);
    realizedPnL = cell(length(tradeDatesAll_OOS), 1);
    
    % Set day count conventions
    if ccy == "EUR"
        floatLegDCC = "act/360";
        fixLegDCC = "30/360"; 
    elseif ccy == "SEK"
        floatLegDCC = "act/360";
        fixLegDCC = "30/360";    
    elseif ccy == "USD"
        floatLegDCC = "act/360";
        fixLegDCC = "act/360";
    end

    currDate = datestr(tradeDatesAll_OOS(1));
    simDate = currDate;

    % Set year fractions for fix cash flows
    dtFix = cell(numContracts,1);
    for i = 1:numContracts
        if fixLegDCC == "act/360"
            dtFix{i}(1) = daysdif(fixingDates{i}(1), fixDates{i}(1), 2) / 360;
            for j = 1:length(fixDates{i})-1
                dtFix{i}(j+1) = daysdif(fixDates{i}(j), fixDates{i}(j+1), 2) / 360;
            end
        elseif fixLegDCC == "30/360"
            dtFix{i}(1) = daysdif(fixingDates{i}(1), fixDates{i}(1), 5) / 360;
            for j = 1:length(fixDates{i})-1
                dtFix{i}(j+1) = daysdif(fixDates{i}(j), fixDates{i}(j+1), 5) / 360;
            end
        end         
    end

    % Initialize parameters for the first day
    timeFracFix{i} = cell(numContracts, 1);
    timeFracFloat{i} = cell(numContracts, 1);
    floatFixing{i} = cell(numContracts, 1);
    dtFloat = cell(numContracts, 1); 
    cash = cell(numContracts, 1); 
    valueParams{1} = simDate;
    valueParams{3} = Ibor;
    valueParams{4} = IborDates;
    valueParams{19} = cash;
    
    for i = 1:numContracts
        
        valueParams{2}(i) = floatLegDCC;
        valueParams{5}{i} = [currDate; floatDates{i}];
        valueParams{6}{i} = [currDate; fixDates{i}];
        valueParams{7}{i} = fixingDates{i};
        valueParams{8}(i) = 0;
        valueParams{9}{i} = dtFix{i};
        valueParams{10}{i} = dtFloat{i};
        valueParams{11}{i} = floatFixing{i};
        valueParams{12}(i) = activeStatus(i);
        valueParams{13}(i) = RoP{i};
        valueParams{14}(i) = Nom(i);
        valueParams{15}(i) = yield(i);
        valueParams{16}{i} = timeFracFix{i};
        valueParams{17}{i} = timeFracFloat{i};  
        valueParams{18}{i} = fixLegDCC;
        valueParams{19}{i} = cash;
        
    end
    
    valueParams = getParameters(valueParams, f, pi);
    

end

