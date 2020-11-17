function portfolioValues = valuePortfolio(rSimulated, piSpotSimulated, valueParams)

numSim = size(rSimulated, 2);
numContracts = length(valueParams{5});


simDate = valueParams{1};
floatDates = valueParams{5};
fixDates = valueParams{6};
numKnown = valueParams{8};
dtFix = valueParams{9};
dtFloat = valueParams{10};
floatFixing = valueParams{11};
activeStatus = valueParams{12};
RoP = valueParams{13};
N = valueParams{14};
y = valueParams{15};
timeFracFix = valueParams{16};
timeFracFloat = valueParams{17};   


% Get indexes to interest rate dates
floatDatesIndexes = cell(numContracts, 1);
fixDatesIndexes = cell(numContracts, 1);
lastday = zeros(numContracts, 1);
for i = 1:numContracts
    if activeStatus(i) == 1
        [floatDatesIndexes{i}, fixDatesIndexes{i}] = getIndexes(simDate, fixDates{i}, floatDates{i});
         % Set flag for last contract day
        if length(floatDates{i}) == 1
            lastday(i) = 1;
        else
            lastday(i) = 0;
        end
    end
end




% Create data result matrix
contractValues = zeros(numSim, numContracts);
portfolioValues = zeros(numSim, 1);

% Value the portfolio for each curve
for j = 1:numSim
        
    r = rSimulated(:,j);
    piSpot = piSpotSimulated(:,j);

    % Value each contract
    for k = 1:numContracts
        % Calculate the portfolio value
        if activeStatus(k) == 1
            if lastday(k) == 0
                contractValues(j,k) = irsPrice(r, piSpot, fixDatesIndexes{k}, floatDatesIndexes{k}, timeFracFix{k}, timeFracFloat{k}, dtFix{k}, dtFloat{k}, floatFixing{k}, RoP(k), y(k), N(k), numKnown(k));
            else
                contractValues(j,k) = 0;
            end
        end
    end
end

portfolioValues = sum(contractValues, 2);

end