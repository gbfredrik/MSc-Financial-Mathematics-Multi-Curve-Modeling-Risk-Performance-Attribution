clear all

measurementPath = 'X:\Examensarbete\Data\';
filename = 'FX_';
contract = 'SEK_';
fBase = readmatrix(strcat(measurementPath,filename,contract,'base.csv'),'Delimiter',';'); %, 'NumHeaderLines',1);
fTerm = readmatrix(strcat(measurementPath,filename,contract,'term.csv'),'Delimiter',';');
pi = readmatrix(strcat(measurementPath,filename,contract,'fx.csv'),'Delimiter',';');
piAvg = readmatrix(strcat(measurementPath,filename,contract,'fxAvg.csv'),'Delimiter',';');
basis = readmatrix(strcat(measurementPath,filename,contract,'basis.csv'),'Delimiter',';'); %Basis from Reuters for FX swap contracts
maturitiesFx = readmatrix(strcat(measurementPath,filename,contract,'maturityDateFX.csv'),'Delimiter',';'); %Maturities for FX swap contracts
T = readmatrix(strcat(measurementPath,filename,contract,'T.csv'),'Delimiter',';');
exchangeRate = readmatrix(strcat(measurementPath,filename,contract,'exchangeRateMid.csv'),'Delimiter',';');
midPriceFX = readmatrix(strcat(measurementPath,filename,contract,'midPrices.csv'),'Delimiter',';'); %Calculated price
maturitiesBase = readmatrix(strcat(measurementPath,filename,contract,'maturitiesBase.csv'),'Delimiter',';');
maturitiesTerm = readmatrix(strcat(measurementPath,filename,contract,'maturitiesTerm.csv'),'Delimiter',';');
midPriceTerm = readmatrix(strcat(measurementPath,filename,contract,'midPriceTerm.csv'),'Delimiter',';');
exchangeRate = exchangeRate(2,:);

[nTradeDates, ~] = size(maturitiesFx); 
%basis = [zeros(nTradeDates, 1) basis];
%maturitiesFx = [zeros(nTradeDates, 1) maturitiesFx];

maturitiesBase = [zeros(nTradeDates, 1) maturitiesBase]; %Ta bort och lägg i funktionerna
maturitiesTerm = [zeros(nTradeDates, 1) maturitiesTerm];
midPriceTerm = [zeros(nTradeDates, 1) midPriceTerm];

%% RAW INTERPOLATION

%% Simulate basis curves for in-sample data set
inSample = 1:20;
N = 2000;

nFXContracts = 22;
%curveLength = 730;
midPriceFX = midPriceFX(:,1:nFXContracts);
maturitiesFx = maturitiesFx(:,1:nFXContracts);
basis = basis(:,1:nFXContracts);
%basisCurve = zeros(nTradeDates,curveLength);

[nTradeDates, nRICs] = size(maturitiesFx);
outOfSample = inSample(end)+1:nTradeDates; %1783
outOfSampleSize = size(outOfSample,2);
 
[~, curveLength] = size(fTerm);
basisCurve = ones(nTradeDates, curveLength)*inf;

basisSimulationsAllDays = cell(1,outOfSampleSize);
lastBasisPoint = curveLength-20;

% Calculate the basis curve
for tradeDate = 1:nTradeDates
    basisCurve(tradeDate,:) = curveFxSwapRawInterpol(maturitiesFx(tradeDate,:), basis(tradeDate,:), curveLength);
    plot(1:curveLength, basisCurve(tradeDate,:))
    hold on
end
hold off

basisCurve = basisCurve(:,1:lastBasisPoint); % Truncate NaN

D_Demand = basisCurve(inSample(1)+1:inSample(end),:) - ...
            basisCurve(inSample(1):inSample(end)-1,:);
        
mu = mean(D_Demand);
CDemand = cov(D_Demand);

h = waitbar(0, "Simulating...");
for i = 1:size(outOfSample,2)-1
    waitbar(i/(size(outOfSample,2)-1), h)
    basisCurrent = basisCurve(i,:);
    deltaF = lhsnorm(mu', CDemand, N)';
    basisSimulations = (repmat(basisCurrent', 1, N) + deltaF)';
    basisSimulationsAllDays{1,i+1} = basisSimulations;
end

% Add spike to all curves

close(h)


%% Calculate realised portfolio value changes for all out-of-sample dates

RealisedPortfolioValueChanges = zeros(1,size(outOfSample,2) - 1);



for i = 1:size(outOfSample,2)-1
     forwardPrices_t2 = zeros(1,nFXContracts);  
     for k = 1:nFXContracts
        %Morgondagens priser, sätt in kurva för morgondagen.
        maturityContractNextDay = (maturitiesFx(outOfSample(1)+i,k) - 1/365);
        FXMaturitiesNextDay = [0 (maturitiesFx(outOfSample(1)+i,:) - 1/365)];
        basisCurve = [0 basis(outOfSample(1)+i,:)];
        forwardPrices_t2(k) = priceFXswapRawInterpol(exchangeRate(outOfSample(1)+i), ...
                                                     basisCurve, ...
                                                     FXMaturitiesNextDay, ...
                                                     maturityContractNextDay); 
        %Dagens priser, använd dagens kurvor.
        %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, ...
        %                                             basisRIC, ...
        %                                             maturitiesRIC, ...
        %                                             maturity(k,i)); 
     end

    %indNaN = find(isnan(maturitiesFx(outOfSample(1)+i,:)), 1, 'first');
    quotedPrices_t1 = midPriceFX(outOfSample(1)+i-1,:);
    maturities = (maturitiesFx(outOfSample(1)+i-1,:) - ones(1,nFXContracts)./365);
    %forwardPrices_t2 = forwardPrices_t2(1:indNaN-1);
    
    %Portfolio value is zero at trade date so portfolio change from the day before is portfolio
    %value of tomorrow.
    
    portfolioValueNextDay = portfolioValueFxSwap(forwardPrices_t2, quotedPrices_t1, 1, ...
                                                 fTerm(outOfSample(1)+i,:), ...
                                                 T(outOfSample(1)+i,:), ...
                                                 maturities);
                                             
    RealisedPortfolioValueChanges(i+1) = portfolioValueNextDay;
end


%% Calculate simulated portfolio value changes given all next day curves for each in sample date and current portfolio value

simulatedPortfolioValues_t2 = zeros(size(outOfSample,2) - 1, N);
simulatedPortfolioValues_t1 = zeros(size(outOfSample,2) - 1, N);
simPortfolioValueChanges =  zeros(size(outOfSample,2) - 1, N);


%Iterera över alla in-sample dagar
for i = 1:size(outOfSample,2)-1
        %iterera över alla simulerade kurvor
        for n = 1:N
            forwardPrices_t2 = zeros(1,nFXContracts);    
            forwardPrices_t1 = zeros(1,nFXContracts);    
            %Iterera över alla kontrakt
            for k = 1:nFXContracts
                %Morgondagens priser, sätt in kurva för morgondagen.
                maturityContractNextDay = (maturitiesFx(outOfSample(1)+i,k) - 1/365);
                FXMaturitiesNextDay = [0 (maturitiesFx(outOfSample(1)+i,:) - 1/365)];
                simulatedCurve = [0 basisSimulationsAllDays{1,i+1}(n,:)];
                forwardPrices_t2(k) = priceFXswapRawInterpol(exchangeRate(outOfSample(1)+i), ...
                                                            simulatedCurve, ...
                                                            FXMaturitiesNextDay, ... 
                                                            maturityContractNextDay); 
                %Dagens priser, använd dagens kurvor.
                %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, basisRIC, maturitiesRIC, maturity(k,i)); %Sätta in kurva?
            end
            quotedForwardPrices = midPriceFX(outOfSample(1)+i-1,1:nFXContracts);
            maturitiesNextDay = (maturitiesFx(outOfSample(1)+i,1:nFXContracts) - ones(1,nFXContracts)./365);
            simulatedPortfolioValues_t2(i+1,n) = portfolioValueFxSwap(forwardPrices_t2, ...
                                                                       quotedForwardPrices, 1, ...
                                                                       fTerm(outOfSample(1)+i,:), ...
                                                                       T(outOfSample(1)+i,:), ...
                                                                       maturitiesNextDay);
            %Beräknat Portföljvärde med våra kurvor, ska egentligen bli
            %noll.
            %simulatedPortfolioValues_t1(i, n) = portfolioValueFxSwap(forwardPrices_t1, ...
            %                                                         quotedForwardPrices, ...
            %                                                         N, ...
            %                                                         fTerm, ...
            %                                                         maturitiesFx(i), ...
            %                                                         maturityContracts);
        
        simPortfolioValueChanges(i+1,n) = simulatedPortfolioValues_t2(i+1,n); % - simPortfolioValues(i) if initial value is not zero.
        end

        
        
end


%% Caluculate kernel for each in-saple date

kernelsRawInterpol = kdeMulti(simPortfolioValueChanges, RealisedPortfolioValueChanges);







