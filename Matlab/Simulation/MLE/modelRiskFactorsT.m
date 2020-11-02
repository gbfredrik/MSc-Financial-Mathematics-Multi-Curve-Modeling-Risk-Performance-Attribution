function [OptParamsAll, rhoHat, nuCopula, kappa] = modelRiskFactorsT(delta_curves, E)

%delta_curves(1:3450,:) = curves(2:end,:) - curves(1:end-1,:);

% Get historic risk factors
RiskFactors = (E'*delta_curves')';

k = size(E,2);

U = zeros(size(delta_curves,1), k);
OptParamsAll = zeros(5,k);

for i = 1:k
    % x = [ nu    omega   beta  alpha  df]
    x0 = [0,0.001,0.92,0.07,5];

    % Estimate Students t-parameters
    [xOpt,~] = mlGARCHStudents(RiskFactors(:,i),1, @likelihoodStudentst, x0);

    v = varGARCH(xOpt,RiskFactors(:,i),1);
    xi = (RiskFactors(:,i)-xOpt(1))./(sqrt(v(1:end-1)));

    u = tcdf(xi,xOpt(5));
    U(:,i) = u;
    OptParamsAll(:,i) = xOpt;
end


[rhoHat, nuCopula] = copulafit('t',U);
%y_T = copulapdf('t', U, rhohat_T, nuhat);
%t = sum(log(y_T));


muHat = mean(RiskFactors, 2);
lh = E' * delta_curves';
rh = (muHat - RiskFactors)';

[kappa, FVal] = fmincon(@(kappa) sum((lh(:,2:end) - kappa .* rh(:,1:end-1)).^2, 'all'), 0.15);

end