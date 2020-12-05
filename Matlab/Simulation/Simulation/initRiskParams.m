function [portfolioValues, portfolioValuesNext, portfolioValuesSimMC, Risk] = initRiskParams(numDates, numSims)
%INITRISKPARAMS Preallocation of arrays and cell arrays for risk measuring

portfolioValues = cell(numDates, 1);
portfolioValuesNext = cell(numDates, 1);
portfolioValuesNext(1, :) = {0};
portfolioValuesSimMC = cell(1, numDates);
portfolioValuesSimMC(:, :) = {zeros(numSims, 1)};

Risk.VaR_95s = zeros(1, numDates);
Risk.VaR_975s = zeros(1, numDates);
Risk.VaR_99s = zeros(1, numDates);
Risk.ES_975s = zeros(1, numDates);
Risk.PnL = zeros(1, numDates);
Risk.Tail = cell(1, numDates);
end
