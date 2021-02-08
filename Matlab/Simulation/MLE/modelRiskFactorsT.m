function [OptParamsAll, rhoHat, nuCopula, like_t, like_garch] = modelRiskFactorsT(delta_curves, E, useMR)

%delta_curves(1:3450,:) = curves(2:end,:) - curves(1:end-1,:);

% Get historic risk factors
RiskFactors = (E'*delta_curves')';

k = size(E,2);
if (useMR)
    U = zeros(size(delta_curves,1)-1, k);
else
    U = zeros(size(delta_curves,1), k);
end
OptParamsAll = zeros(5,k);

for i = 1:k
    % x = [ mu    omega   beta  alpha  df]
    x0 = [0,0.001,0.92,0.07,5];
    x0_garch = [0,0.001,0.92,0.07];
    % Estimate Students t-parameters
    if (useMR)
        [xOpt,ll_t] = mlGARCHStudents(RiskFactors(:,i),1, @likelihoodStudentstMR, x0, useMR);
        [xOpt_garch,ll_garch] = mlGARCH(RiskFactors(:,i),1, @likelihoodNormalMR, x0_garch, useMR);
    else
        [xOpt,ll_t] = mlGARCHStudents(RiskFactors(:,i),1, @likelihoodStudentst, x0, useMR);
        [xOpt_garch,ll_garch] = mlGARCH(RiskFactors(:,i),1, @likelihoodNormal, x0_garch, useMR);
    end
    like_t(i) = sum(ll_t);
    like_garch(i) = sum(ll_garch);
    v = varGARCH(xOpt,RiskFactors(:,i),1);
    if (useMR)
        xi = (RiskFactors(2:end,i) - xOpt(1) .* (mean(RiskFactors(:,i)) - RiskFactors(1:end-1,i)))./(sqrt(v(2:end-1)));
    else
        xi = (RiskFactors(:,i)-xOpt(1))./(sqrt(v(1:end-1)));
    end
    u = tcdf(xi,xOpt(5));
    U(:,i) = u;
    OptParamsAll(:,i) = xOpt;
end


[rhoHat, nuCopula] = copulafit('t',U);
%y_T = copulapdf('t', U, rhohat_T, nuhat);
%t = sum(log(y_T));

end