function [simParams] = initSimParams(N, E, DZero, DPi, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS)
% Initialize sim params

T = 1:3650;
simParams = {}; 
for i = 1:8
    for j = 1:2
        simParams{i}{j} = 0;
    end
end

[muZero, omegaZero, betaZero, alphaZero, dfMZero, rhoZero, dfCZero, kappaZero, like_t_zero, like_garch_zero] = calibrateParams(DZero, E.Zero);
[muPi, omegaPi, betaPi, alphaPi, dfMPi, rhoPi, dfCPi, kappaPi, like_t_pi, like_garch_pi] = calibrateParams(DPi, E.Pi);

% Calculate sigma the first day based on the IS data
sigmaZero = GJR_GARCH(omegaZero', alphaZero', betaZero', E.Zero, fAll_IS);
sigmaPi = GJR_GARCH(omegaPi', alphaPi', betaPi', E.Pi, piAll_IS); 



simParams{1}{1} = muZero;
simParams{1}{2} = muPi;
simParams{2}{1} = omegaZero;
simParams{2}{2} = omegaPi;
simParams{3}{1} = betaZero;
simParams{3}{2} = betaPi;
simParams{4}{1} = alphaZero;
simParams{4}{2} = alphaPi;
simParams{5}{1} = dfMZero;
simParams{5}{2} = dfMPi;
simParams{6}{1} = rhoZero;
simParams{6}{2} = rhoPi;
simParams{7}{1} = dfCZero;
simParams{7}{2} = dfCPi;
simParams{8}{1} = sigmaZero; 
simParams{8}{2} = sigmaPi;
simParams{9}{1} = DZero;
simParams{9}{2} = DPi;
simParams{10}{1} = E.Zero;
simParams{10}{2} = E.Pi;
simParams{11} = N;
simParams{12} = T;
simParams{13}{1} = kappaZero;
simParams{13}{2} = kappaPi;
end

