function [fSimulated, piSimulated, simParams] = TermStructureSim(i, simParams, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, useMR)

% Read parameters from simParams
muZero = simParams{1}{1};
muPi = simParams{1}{2};
omegaZero = simParams{2}{1};
omegaPi = simParams{2}{2};
betaZero = simParams{3}{1};
betaPi = simParams{3}{2};
alphaZero = simParams{4}{1};
alphaPi = simParams{4}{2};
dfMZero = simParams{5}{1};
dfMPi = simParams{5}{2};
rhoZero = simParams{6}{1};
rhoPi = simParams{6}{2};
dfCZero = simParams{7}{1};
dfCPi = simParams{7}{2};
sigmaZero = simParams{8}{1};
sigmaPi = simParams{8}{2};
DZero = simParams{9}{1};
DPi = simParams{9}{2};
EZero = simParams{10}{1};
EPi = simParams{10}{2};
N = simParams{11};
T = simParams{12};
% Mean-reversion parameters
kappaZero = simParams{13}{1};
kappaPi = simParams{13}{2};
% Simulate 1d ahead

% Recalibrate parameters every two years     
if i == 1 || mod(i,504) == 0
   %if i <= 1512 && i ~= 1
   if i <= size(fAll_IS,1) && i ~= 1     
       %start = size(DZero,1) - (1512 - i);
       start = size(DZero,1) - (size(fAll_IS,1) - i);
       DZero = [DZero(start:end,:); fAll_OOS(2:i,:) - fAll_OOS(1:i-1,:)];
       DPi = [DPi(start:end,:); piAll_OOS(2:i,:) - piAll_OOS(1:i-1,:)];
   elseif i > size(fAll_IS,1)
   %elseif i > 1005
       %start = i - 1511;
       start = i - size(fAll_IS,1);
       DZero = fAll_OOS(start+1:i,:) - fAll_OOS(start:i-1,:);
       DPi = piAll_OOS(start+1:i,:) - piAll_OOS(start:i-1,:);
   end
   
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

% Set active curve at day i-1 (the trade date)
fZero = fAll_OOS(i-1,:);
pi = piAll_OOS(i-1,:);

if i == 2
   % Calculate sigma the first day based on the IS data (i.e. trade date 1)
    sigmaZero = GJR_GARCH(omegaZero', alphaZero', betaZero', EZero, fAll_IS);
    sigmaPi = GJR_GARCH(omegaPi', alphaPi', betaPi', EPi, piAll_IS); 
else
    % Calc sigma based on last day's vol
    sigmaZero = GJR_GARCH_d(omegaZero', alphaZero', betaZero', EZero, fZero, fAll_OOS(i-2,:), sigmaZero);
    sigmaPi = GJR_GARCH_d(omegaPi', alphaPi', betaPi', EPi, pi, piAll_OOS(i-2,:), sigmaPi);
end
% Calc random numbers
UZero = lhsstud(zeros(1,length(omegaZero)), rhoZero, dfCZero, N);
for j=1:length(omegaZero)
    ZZero(:,j) = tinv(UZero(:,j), dfMZero(j));
end
UPi = lhsstud(zeros(1,length(omegaPi)), rhoPi, dfCPi, N);
for j=1:length(omegaPi)
    ZPi(:,j) = tinv(UPi(:,j), dfMPi(j));
end

% Calc the 1d ahead simulated curves
if (useMR)
    xiHatZero = mean(EZero' * DZero', 2);
    xiHatPi = mean(EPi' * DPi', 2);
    mrZero = kappaZero .* (xiHatZero - EZero' * DZero(i-1,:)'); % i-1?
    mrPi = kappaPi .* (xiHatPi - EPi' * DPi(i-1,:)'); % i-1?

    deltaF = EZero * (mrZero + ZZero'.* repmat(sigmaZero,1,N));
    deltaPi = EPi * (mrPi + muPi' + ZPi'.* repmat(sigmaPi,1,N));
else
    deltaF = EZero * (muZero' + ZZero'.* repmat(sigmaZero,1,N));
    deltaPi = EPi * (muPi' + ZPi'.* repmat(sigmaPi,1,N));
end
% Compute simulated curves
fSimulated = (repmat(fZero',1,N) + deltaF)';
piSimulated = (repmat(pi',1,N) + deltaPi)';

% Plot results
% Zero
%plotResults(fZero, T, fSimulated, tradeDatesAll_OOS, i, ['Simulated 1d risk-free IR curves, Date: ', num2str(datestr(tradeDatesAll_OOS(i)))], 1);
%  Pi
%plotResults(pi, T, piSimulated, tradeDatesAll_OOS, i, ['Simulated 1d Tenor spreads, Date: ', num2str(datestr(tradeDatesAll_OOS(i)))], 2);
 % Zero + Pi   
%plotResults(fZero + pi, T, fSimulated + piSimulated, tradeDatesAll_OOS, i, ['Simulated 1d 3M IR curves, Date: ', num2str(datestr(tradeDatesAll_OOS(i)))], 3);

% Set sigma to be used in next day simulation
simParams{8}{1} = sigmaZero; 
simParams{8}{2} = sigmaPi;
end
