%% Risk Data and Backtesting Wrapper
% This is a wrapper for the generation of risk data plots and 
% risk backtesting
clear all
%% Data loading
%risk_path = 'Data/Risk/Risk_ABCD_750_noMR.mat';
%risk_path = 'Data/Risk/Risk_A_750_MR.mat';
risk_path = 'Data/Risk/Risk_SEK_E_IS_01_OOS_01_useMRfalse.mat';
load(risk_path);

%% Risk measurement plots
close all
f_rm = figure(1);
tt = 1:length(Risk.PnL(Risk.PnL ~= 0));

hold on
%plot(tt, Risk.PnL(tt), tt, -Risk.VaR_975s(tt), tt, -Risk.VaR_99s(tt), tt, -Risk.ES_975s(tt))
scatter(tt, Risk.PnL(tt), '*')
plot(tt, -Risk.VaR_975s(tt))
plot(tt, -Risk.VaR_99s(tt))
plot(tt, -Risk.ES_975s(tt))
axis tight
legend("PnL", "VaR975", "VaR99", "ES97.5")
hold off

%% VaR Backtest
fprintf('\n\nValue-at-Risk backtesting for %i trade dates:\n', length(tt));
fprintf('Number of breaches for 975-VaR: %i.\n', sum(-Risk.VaR_975s(tt) > Risk.PnL(tt)));
fprintf('Number of breaches for 99-VaR: %i.\n', sum(-Risk.VaR_99s(tt) > Risk.PnL(tt)));
% Test 1: Hypothesis test
[Risk.Test.Hyp975, Risk.Test.HypZ975, tilde_m] = var_hypothesistest(Risk.VaR_975s(tt), Risk.PnL(tt), 0.975, 0.05);
fprintf('Hypothesis test for 975-VaR is rejected: %i, for Z statistic vs tilde_m "%.3f <---> %.3f".\n', Risk.Test.Hyp975, abs(Risk.Test.HypZ975), tilde_m);

[Risk.Test.Hyp99, Risk.Test.HypZ99, tilde_m] = var_hypothesistest(Risk.VaR_99s(tt), Risk.PnL(tt), 0.99, 0.05);
fprintf('Hypothesis test for 99-VaR is rejected: %i, for Z statistic vs tilde_m "%.3f <---> %.3f".\n', Risk.Test.Hyp99, abs(Risk.Test.HypZ99), tilde_m);

% Test 2: Christoffersen's test
[Risk.Test.Christ975, Risk.Test.ChristC975, chi_val] = var_christoffersentest(Risk.VaR_975s(tt), Risk.PnL(tt), 0.05);
fprintf('Christoffersen´s test for 975-VaR is rejected: %i, for C vs chiInv "%.3f <---> %.3f".\n', Risk.Test.Christ975, Risk.Test.ChristC975, chi_val);

[Risk.Test.Christ99, Risk.Test.ChristC99, chi_val] = var_christoffersentest(Risk.VaR_99s(tt), Risk.PnL(tt), 0.05);
fprintf('Christoffersen´s test for 99-VaR is rejected: %i, for C vs chiInv "%.3f <---> %.3f".\n', Risk.Test.Christ99, Risk.Test.ChristC99, chi_val);


%% ES Backtest
fprintf('\nExpected Shortfall backtesting:\n');

phi = 0.95;
[reject, p] = es_acerbiszekely(Risk.Tail, Risk.VaR_975s(tt), Risk.ES_975s(tt), Risk.PnL(tt), phi);
fprintf('Acerbi & Szekely test for 97.5-ES is rejected: %f vs %f.\n', p, phi);


%%
%f = surf(fAll_OOS(1:252,:));
%shading flat


