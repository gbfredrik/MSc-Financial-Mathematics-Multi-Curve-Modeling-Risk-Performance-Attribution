%% Setup data and parameters
clear
outcomes = randn(100,1); % Simulate random dsitributed outcomes
VaRs = abs(randn(500,1) .* 0.5 - 5);
PnLs = randn(500,1) .* 2.5;
%X = randn(500, 252)
c = 0.95; % Confidence level

% sum(-VaRs > PnLs) % For controlling number of breaches
%% Test VaR
[var, ind] = var_risk(outcomes, c);

%% Test ES
es = es_risk(outcomes, c);

%% Backtest VaR
%  Hypothesis testing: Reject H0 for H1 if true
disp("Reject H0 for H1: " + var_hypothesistest(VaRs, PnLs, c, 0.95));

% Christoffersen testing: Reject if true
disp("Reject model due to clustering: " + var_christoffersentest(VaRs, PnLs, 0.95));

%% Backtest ES
%disp("Reject ES: " + es_acerbiszekely(X, VaRs, ESs, PnLs, 0.95));

