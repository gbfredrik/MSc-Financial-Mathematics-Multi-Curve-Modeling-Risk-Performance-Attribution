function [simParams] = estimateSimParams(DZero, EZero, DPi, EPi)

[muZero, omegaZero, betaZero, alphaZero, dfMZero, rhoZero, dfCZero, kappaZero] = calibrateParams(DZero, EZero);
[muPi, omegaPi, betaPi, alphaPi, dfMPi, rhoPi, dfCPi, kappaPi] = calibrateParams(DPi, EPi);
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
simParams{13}{1} = kappaZero;
simParams{13}{2} = kappaPi;
end
