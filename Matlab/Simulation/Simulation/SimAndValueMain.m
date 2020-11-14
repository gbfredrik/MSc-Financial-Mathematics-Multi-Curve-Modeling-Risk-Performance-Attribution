clear all
%% Import Forward rates
path_IS = 'Data/Curves/EUR_IS_10YrCurves_Clean_Final2.mat';
path_OOS = 'Data/Curves/EUR_OOS_10YrCurves_Clean_Final.mat';
[fAll_IS, piAll_IS, tradeDatesAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, ccy, fDatesAll] = getData(path_IS, path_OOS);
clearvars path_IS path_OOS
%% Get portfolio data
instruments = 'A':'P';
[floatDates, fixDates, yield, fixingDates, RoP, IborDates, Ibor, Nom] = getPortfolioData(instruments, ccy);

%% Get/set risk factors
kZero = 6;
kPi = 6;
% type 1: Matlab eigenvectors
% type 2: BDCSVD
% type 3: BDCSVD and IRAM
type = 1;
[E, E_k, DZero, DPi] = getRiskFactors(kZero, kPi, fAll_IS, piAll_IS, type);
A = intMatrix(size(E.Zero,1));

%% Initialize sim params
N = 2000; % Number of curves to be simulated
simParams = initSimParams(N, E_k, DZero, DPi, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS);

%% Init PA parameters
[valueParams, ~] = initValueParams(tradeDatesAll_OOS, fixingDates, floatDates, fixDates, ccy, Ibor, IborDates, RoP, Nom, yield, fAll_OOS(1,:)', piAll_OOS(1,:)');
numContracts = length(valueParams{2});
[paParams, paResult, plotPAParams] = initPAParams(numContracts, A, E, E_k, fAll_OOS(1,:)', piAll_OOS(1,:)', fAll_IS(end,:)', piAll_IS(end,:)');

%% Initialize valuation and risk parameters, as well as run main loop
% Init valuation parameters (redo just incase we only run the simulation
% multiple times without full clearing of data)
[valueParams, cash] = initValueParams(tradeDatesAll_OOS, fixingDates, floatDates, fixDates, ccy, Ibor, IborDates, RoP, Nom, yield, fAll_OOS(1,:)', piAll_OOS(1,:)');
fprintf('First day valuation is %.3f.\n', valuePortfolio(A*fAll_OOS(1,:)', A*piAll_OOS(1,:)', valueParams)); % Test the first day valuation

% Init risk parameters
[portfolioValues, portfolioValuesNext, portfolioValuesSimMC, Risk] = initRiskParams(length(tradeDatesAll_OOS));

% Loop over all days
useMR = false; % Use mean-reversion for simulation of curves

for i = 1:1%length(tradeDatesAll_OOS)
    fprintf('Currently on iteration %i.\n', i)
    currDate = datestr(tradeDatesAll_OOS(i));
    
    % Calculate realized portfolio values today and for next wd (backtesting)
    portfolioValues{i} = valuePortfolio(A*fAll_OOS(i,:)', A*piAll_OOS(i,:)', valueParams);
    portfolioValuesNext{i} = valuePortfolio(A*fAll_OOS(i+1,:)', A*piAll_OOS(i+1,:)', valueParams);
    fprintf('portfolioValues{%i} = %.3f, portfolioValuesNext{%i} = %.3f.\n', i, portfolioValues{i}, i, portfolioValuesNext{i})
    
    %Simulate 1d ahead
    [fSimulated, piSimulated, simParams] = TermStructureSim(i+1, simParams, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, useMR);
    
    % Value Portfolio with MC simulation
    if sum(valueParams{12}) > 0
        portfolioValuesSimMC{i} = valuePortfolio(A*fSimulated', A*piSimulated', valueParams);
    end
    
    % Risk Measurement
    Risk.VaR_95s(i) = var_risk(portfolioValuesSimMC{i} - portfolioValues{i}, 0.95);
    Risk.VaR_99s(i) = var_risk(portfolioValuesSimMC{i} - portfolioValues{i}, 0.99);
    Risk.ES_975s(i) = es_risk(portfolioValuesSimMC{i} - portfolioValues{i}, 0.975);
    Risk.PnL(i) = portfolioValuesNext{i} - portfolioValues{i};
    
    % Call Performance attribution
%     if sum(valueParams{12}) > 0
%         paParams = getPAParameters(paParams, A, fAll_OOS(i,:)', piAll_OOS(i,:)');
%         [paParams, paResult] = PA(paParams, paResult, valueParams);
%         plotPAParams = plotPA(paResult, plotPAParams, currDate, valueParams{12});
%     end
    
    % Break on last day (lite onajs, men nu blev det så)
    if i == length(tradeDatesAll_OOS)
        break
    end
    
    nextWD = datestr(tradeDatesAll_OOS(i+1));
    valueParams{1} = nextWD;

    % Update cashflows
    valueParams = setActiveStatus(valueParams, nextWD);
    valueParams = getParameters(valueParams, fAll_OOS(i,:)', piAll_OOS(i,:)');
end

save("Data/Risk/Risk_" + instruments, "Risk")
%%
ttt = 1:252;
figure(101)
hold on
%plot(ttt, Risk.PnL(ttt), ttt, -Risk.VaR_95s(ttt), ttt, -Risk.VaR_99s(ttt), ttt, -Risk.ES_975s(ttt))
scatter(ttt, Risk.PnL(ttt), '*')
plot(ttt, -Risk.VaR_95s(ttt))
plot(ttt, -Risk.VaR_99s(ttt))
plot(ttt, -Risk.ES_975s(ttt))
legend("TheoPnL", "VaR95", "VaR99", "ES97.5")
hold off

clearvars ttt
%%
%load('paResult')
%load('plotPAParams')
%% Plot PA
%plotPAFinal(plotPAParams, tradeDatesAll_OOS)
%diff = paResult{1}{1}' - paResult{2}{1}' - paResult{3}{1}' - paResult{4}{1}' - paResult{5}{1}' - paResult{16}{1}' 

%% Calculate covariance matrix
%[Exp, Cov, Tot_exp] = calcCovar(plotPAParams);

%% 1Y IRS table
%[oneY_IRS_table, allY_first_IRS_table, abs_cont_all, abs_cont_mature] = sum_1Y_results(paResult, tradeDatesAll_OOS);


%%





