%% Risk Data and Backtesting Wrapper
% This is a wrapper for the generation of risk data plots and 
% risk backtesting
clear all
%% Data loading
risk_path = 'Data/Risk/';
file_name = 'Risk.mat';
Risk = load(risk_path + file_name);

%% Risk measurement plots
close all
f_rm = figure(1);
tt = 1:length(Risk.PnL(Risk.PnL ~= 0));

hold on
%plot(tt, Risk.PnL(tt), tt, -Risk.VaR_95s(tt), tt, -Risk.VaR_99s(tt), tt, -Risk.ES_975s(tt))
scatter(tt, Risk.PnL(tt), '*')
plot(tt, -Risk.VaR_95s(tt))
plot(tt, -Risk.VaR_99s(tt))
plot(tt, -Risk.ES_975s(tt))
axis tight
legend("PnL", "VaR95", "VaR99", "ES97.5")
hold off

%% VaR Backtest
% Test 1: Hypothesis test
Risk.Var_hyp = var_hypothesistest(Risk.VaR_99s(tt), Risk.PnL(tt), 0.05, 0.9999);
fprintf('Hypothesis test for VaR 95 is rejected: %b.\n', Risk.Var_hyp);

%%
%%
%f = surf(fAll_OOS(1:252,:));
%shading flat


