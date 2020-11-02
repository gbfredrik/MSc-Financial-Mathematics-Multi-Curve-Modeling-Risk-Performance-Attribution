%% Data loading
clear
IS = load("Data/Curves/EUR_IS_10YrCurves_Clean_Final2.mat");
OOS = load("Data/Curves/EUR_OOS_10YrCurves_Clean_Final.mat");


%% Cut EUR IS
% IS.TAll = IS.TAll(256:end,:);
% IS.fAll = IS.fAll(256:end,:);
% IS.fDatesAll = IS.fDatesAll(256:end,:);
% IS.piAll = IS.piAll(256:end,:);
% IS.tradeDatesAll = IS.tradeDatesAll(1,256:end);
% IS.zAll = IS.zAll(256:end,:);
%save("Data/Curves/EUR_IS_10YrCurves_Clean_Final2.mat", '-struct', 'IS');

%% Test fitting
delta_curvesIS = diff(IS.fAll(:,1:3650));
cIS = cov(delta_curvesIS);
[EIS, DIS] = eigs(cIS, 6);
RiskFactors = (EIS' * delta_curvesIS');

%%
[X, FVAL] = fmincon(@(kappa) sum(sum((delta_curvesIS - (EIS * (kappa .* (mean(RiskFactors, 2) - RiskFactors)))').^2)), [0.05]);
%%
%(delta_curvesIS - (EIS * (kappa .* (mean(RiskFactors, 2) - RiskFactors)))');
%(((EIS' * delta_curvesIS')' - kappa .* (mean(RiskFactors, 2) - RiskFactors)).^2);
%mean(RiskFactors, 2) - RiskFactors;

%test_kappa = fit_kappa(IS.fAll(2,1:3650), IS.fAll(1,1:3650), EIS, mean(RiskFactors, 2), RiskFactors(:,1))

%kappa = fmincon(@(kappa) sum_squared_err(IS.fAll(:, 1:3650), EIS, mean(RiskFactors, 2), RiskFactors, kappa), 0.5);
%%

function [se] = sum_squared_err(f, E, RF_hat, RF, kappa)
    se = 0.0;
    
    for i = 1 : size(f, 1)-1
        se = se + fit_kappa(f(i+1, :), f(i, :), E, RF_hat, RF(:, i), kappa);
    end
end

function [f] = fit_kappa(f_tplus1, f_t, E, RF_hat, RF_t, kappa_old)
    %[k, f] = fmincon(@(kappa) sum((f_tplus1 - (f_t + (E * (kappa .* (RF_hat - RF_t)))')).^2), kappa_old);
    f = sum((f_tplus1 - (f_t + (E * (kappa_old .* (RF_hat - RF_t)))')).^2);
end
