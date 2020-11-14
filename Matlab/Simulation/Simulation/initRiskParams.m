function [portfolioValues, portfolioValuesNext, portfolioValuesSimMC, Risk] = initRiskParams(N)
%INITRISKPARAMS Preallocation of arrays and cell arrays for risk measuring

portfolioValues = cell(N, 1);
portfolioValuesNext = cell(N, 1);
portfolioValuesNext(1, :) = {0};
portfolioValuesSimMC = cell(1, N);
portfolioValuesSimMC(:, :) = {zeros(N, 1)};

Risk.VaR_95s = zeros(1, N);
Risk.VaR_99s = zeros(1, N);
Risk.ES_975s = zeros(1, N);
Risk.PnL= zeros(1, N);

end
