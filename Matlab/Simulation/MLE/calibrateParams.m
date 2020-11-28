function [mu, omega, beta, alpha, dfM, rho, dfC, kappa, like_t, like_garch] = calibrateParams(delta_curves, E)

% OptParamsAll:
    % Rad 1: mu
    % Rad 2: omega: GARCH-parameter
    % Rad 3: beta: GARCH-parameter
    % Rad 4: alpha: GARCH-parameter
    % Rad 5: dfM - frihetsgrader marginaler

[OptParamsAll, rhoHat, nuCopula, kappa, like_t, like_garch] = modelRiskFactorsT(delta_curves, E);

    

mu = OptParamsAll(1,:);
omega = OptParamsAll(2,:);
beta = OptParamsAll(3,:);
alpha = OptParamsAll(4,:);
dfM = OptParamsAll(5,:);
rho = rhoHat;
dfC = nuCopula;

end
