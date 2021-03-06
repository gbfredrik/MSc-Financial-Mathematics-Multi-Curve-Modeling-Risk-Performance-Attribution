clear all
%% Import Forward rates
path_IS = 'Data/Curves/EUR_IS_10YrCurves_Clean_Final2.mat';
path_OOS = 'Data/Curves/EUR_OOS_10YrCurves_Clean_Final.mat';
[fAll_IS, piAll_IS, tradeDatesAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, ccy, fDatesAll] = getData(path_IS, path_OOS);
if ccy == "SEK"
    fAll_IS = [fAll_IS(1:1005,:)];
    piAll_IS = [piAll_IS(1:1005,:)];
    fAll_OOS = [fAll_OOS(1006:end,:)];
    piAll_OOS = [piAll_OOS(1006:end,:)];
    tradeDatesAll_OOS = (tradeDatesAll_OOS(1006:end));
end
clearvars path_IS path_OOS
%% Get portfolio data
instruments = 'A':'P';
if ccy == "SEK"
    instruments = 'A':'L';
end
portType = 'sim'; %'sim' or 'PA'
[floatDates, fixDates, yield, fixingDates, RoP, IborDates, Ibor, Nom] = getPortfolioData(instruments, ccy, portType);

%% Get/set risk factors
kZero = 6;
kPi = 6;
% type 1: Matlab eigenvectors
% type 2: BDCSVD
% type 3: BDCSVD and IRAM
type = 2;
[E, E_k, DZero, DPi] = getRiskFactors(kZero, kPi, fAll_IS, piAll_IS, type);
A = intMatrix(size(E.Zero,1));

%% Initialize sim params
N = 2000; % Number of curves to be simulated
useMR = true;
simParams = initSimParams(N, E_k, DZero, DPi, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, useMR);

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
[portfolioValues, portfolioValuesNext, portfolioValuesSimMC, Risk] = initRiskParams(length(tradeDatesAll_OOS), N);

% Loop over all days
simHorizon = 1; % Number of trade days ahead to simulate

for i = 1:100%length(tradeDatesAll_OOS) - simHorizon
    fprintf('Currently on iteration %i.\n', i)
    currDate = datestr(tradeDatesAll_OOS(i));
    
    % Calculate realized portfolio values today and for next wd (backtesting)
    portfolioValues{i} = valuePortfolio(A*fAll_OOS(i,:)', A*piAll_OOS(i,:)', valueParams);
    portfolioValuesNext{i} = valuePortfolio(A*fAll_OOS(i + simHorizon,:)', A*piAll_OOS(i + simHorizon,:)', valueParams);
    fprintf('portfolioValues{%i} = %.3f, portfolioValuesNext{%i} = %.3f.\n', i, portfolioValues{i}, i, portfolioValuesNext{i})
    %Simulate 1d ahead

    [fSimulated, piSimulated, simParams] = TermStructureSim_10d(i+1, simParams, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, useMR, simHorizon);    
    
    % Value Portfolio with MC simulation
    if sum(valueParams{12}) > 0
        portfolioValuesSimMC{i} = valuePortfolio(A*fSimulated', A*piSimulated', valueParams);
    end
    
    % Risk Measurement

    simulatedPnLs = portfolioValuesSimMC{i} - portfolioValues{i};
    Risk.VaR_95s(i) = var_risk(simulatedPnLs, 0.95);
    Risk.VaR_975s(i) = var_risk(simulatedPnLs, 0.975);
    Risk.VaR_99s(i) = var_risk(simulatedPnLs, 0.99);
    Risk.ES_975s(i) = es_risk(simulatedPnLs, 0.975);
    Risk.PnL(i) = portfolioValuesNext{i} - portfolioValues{i};

    Risk.Tail{i} = sort(simulatedPnLs(simulatedPnLs <= -Risk.VaR_975s(i)));
    fprintf('Measured 95-VaR: %.3f. 99-VaR: %.3f. 97.5-ES: %.3f.\n', Risk.VaR_95s(i), Risk.VaR_99s(i), Risk.ES_975s(i));
    % Call Performance attribution
%     if sum(valueParams{12}) > 0
%         paParams = getPAParameters(paParams, A, fAll_OOS(i,:)', piAll_OOS(i,:)');
%         [paParams, paResult] = PA(paParams, paResult, valueParams);
%         plotPAParams = plotPA(paResult, plotPAParams, currDate, valueParams{12});
%     end

    nextWD = datestr(tradeDatesAll_OOS(i+1));
    valueParams{1} = nextWD;

    % Update cashflows
    valueParams = setActiveStatus(valueParams, nextWD);
    valueParams = getParameters(valueParams, fAll_OOS(i,:)', piAll_OOS(i,:)');
end
%%

%save("Data/Risk/Risk_" + ccy + "_" + instruments + "_useMR" + useMR, "Risk")
ccy = "USD";
useMR = "false";
type = "E_100";
%save("Data/Risk/Risk_" + ccy + "_" + type + "_useMR" + useMR, "Risk")
%%

ttt = 1:length(Risk.PnL(Risk.PnL ~= 0));
figure(101)
hold on
%plot(ttt, Risk.PnL(ttt), ttt, -Risk.VaR_95s(ttt), ttt, -Risk.VaR_99s(ttt), ttt, -Risk.ES_975s(ttt))
scatter(ttt, Risk.PnL(ttt), '*')
plot(ttt, -Risk.VaR_975s(ttt))
plot(ttt, -Risk.VaR_99s(ttt))
plot(ttt, -Risk.ES_975s(ttt))
legend("TheoPnL", "VaR97.5", "VaR99", "ES97.5")
hold off

%%
% 
% for i = 1390:1400
%    figure(1)
%    plot(1:3650, fAll_OOS(i,:)) 
%    figure(2)
%    plot(1:3650, fAll_OOS(i,:) + piAll_OOS(i,:))
%    pause(0.1)
%    i
% end


%%
% deleted = 0;
% for i = 1:1884 
%     if (i == 1) || (mod(i,10) == 0)
%         'hej'
%     else
%         Risk.VaR_95s(i-deleted) = [];
%         Risk.VaR_975s(i-deleted) = [];
%         Risk.VaR_99s(i-deleted) = [];
%         Risk.ES_975s(i-deleted) = [];
%         Risk.PnL(i-deleted) = [];
%         Risk.Tail(i-deleted) = [];
%         deleted = deleted + 1;
%     end
% end


%%
%plotRiskFactors(E_k.Zero, "Forward rate risk factors, risk-free")
%plotTenorSpreads(fAll_OOS(1:252,:), " ")
%%
%load('paResult_SEK')
%load('plotPAParams_SEK')


%% Plot PA
%plotPAFinal(plotPAParams, tradeDatesAll_OOS)
%diff = paResult{1}{1}' - paResult{2}{1}' - paResult{3}{1}' - paResult{4}{1}' - paResult{5}{1}' - paResult{16}{1}' 

%% Calculate covariance matrix
%[Exp, Cov, Tot_exp, type] = calcCovar(plotPAParams);

%% 1Y IRS table
%[oneY_IRS_table, allY_first_IRS_table, abs_cont_all, abs_cont_mature] = sum_1Y_results(paResult, tradeDatesAll_OOS);
