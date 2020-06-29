clear all

measurementPath = 'X:\Examensarbete\Data\';
filename = 'FX_';
contract = 'SEK_';
fBase = readmatrix(strcat(measurementPath,filename,contract,'base.csv'),'Delimiter',';'); %, 'NumHeaderLines',1);
fTerm = readmatrix(strcat(measurementPath,filename,contract,'term.csv'),'Delimiter',';');
demand = readmatrix(strcat(measurementPath,filename,contract,'fx.csv'),'Delimiter',';');
demandAvg = readmatrix(strcat(measurementPath,filename,contract,'fxAvg.csv'),'Delimiter',';');
basis = readmatrix(strcat(measurementPath,filename,contract,'basis.csv'),'Delimiter',';'); %Basis from Reuters for FX swap contracts
maturitiesFx = readmatrix(strcat(measurementPath,filename,contract,'maturityDateFX.csv'),'Delimiter',';'); %Maturities for FX swap contracts
T = readmatrix(strcat(measurementPath,filename,contract,'T.csv'),'Delimiter',';');
exchangeRate = readmatrix(strcat(measurementPath,filename,contract,'exchangeRateMid.csv'),'Delimiter',';');
midPriceFX = readmatrix(strcat(measurementPath,filename,contract,'midPrices.csv'),'Delimiter',';'); %Calculated price
maturitiesBase = readmatrix(strcat(measurementPath,filename,contract,'maturitiesBase.csv'),'Delimiter',';');
maturitiesTerm = readmatrix(strcat(measurementPath,filename,contract,'maturitiesTerm.csv'),'Delimiter',';');
midPriceTerm = readmatrix(strcat(measurementPath,filename,contract,'midPriceTerm.csv'),'Delimiter',';');
jumpDates = readmatrix(strcat(measurementPath,filename,contract,'jumpDatesH.csv'),'Delimiter',';');
exchangeRate = exchangeRate(2,:);

[nTradeDates, ~] = size(maturitiesFx); 
%basis = [zeros(nTradeDates, 1) basis];
%maturitiesFx = [zeros(nTradeDates, 1) maturitiesFx];

maturitiesBase = [zeros(nTradeDates, 1) maturitiesBase]; %Ta bort och lägg i funktionerna
maturitiesTerm = [zeros(nTradeDates, 1) maturitiesTerm];
midPriceTerm = [zeros(nTradeDates, 1) midPriceTerm];

%% Setup
inSample = 1:1750;
N = 2000;

nFXContracts = 22;
midPriceFX = midPriceFX(:,1:nFXContracts);
maturitiesFx = maturitiesFx(:,1:nFXContracts);

[nTradeDates, nRICs] = size(maturitiesFx);
outOfSample = inSample(end)+1:nTradeDates; %1783
outOfSampleStart = size(inSample,2) + 1;
outOfSampleEnd = outOfSample(end);
outOfSampleSize = size(outOfSample, 2);

curveLength = size(fTerm,2);

lastBasisPoint = curveLength-20;


%% OPTIMISED CURVES

[fTermSimulationsAllDays, fBaseSimulationsAllDays, fDemandSimulationsAllDays] = simulateFXCurves(fBase, fTerm, demandAvg, N, inSample, outOfSample, lastBasisPoint);

%%

fDemandSimulationsAllDays = addJmp(demand(outOfSampleStart:outOfSampleEnd-1,:), fDemandSimulationsAllDays, jumpDates(outOfSampleStart:outOfSampleEnd-1,:)); 

%% Calculate realised portfolio value changes for all in-date samples

realisedPortfolioValueChanges_O = zeros(1,size(outOfSample,2) - 1);

for i = 1:size(outOfSample,2)-1
    forwardPrices_t2 = zeros(1,nFXContracts);  
    for k = 1:nFXContracts
        %Morgondagens priser, sätt in realiserade kurvor för morgondagen.
        maturityContractNextDay = (maturitiesFx(outOfSample(1)+i-1,k) - 1/365);
        forwardPrices_t2(k) = priceFXswap(exchangeRate(outOfSample(1)+i), ...
                                          fBase(outOfSample(1)+i,:), ...
                                          fTerm(outOfSample(1)+i,:), ...
                                          demand(outOfSample(1)+i,:), ...
                                          T(outOfSample(1)+i-1,:), ...
                                          maturityContractNextDay); 
        %Dagens priser, använd dagens kurvor.
        %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, basisRIC, maturitiesRIC, maturity(k,i)); %Sätta in kurva?
    end

    quotedPrices_t1 = midPriceFX(outOfSample(1)+i-1,1:nFXContracts);
    
    maturities = (maturitiesFx(outOfSample(1)+i-1,1:nFXContracts) - ones(1,nFXContracts)./365);
    
    %Portfolio value is zero at trade date so portfolio change from the day before is portfolio
    %value of tomorrow.
    portfolioValueNextDay = portfolioValueFxSwap(forwardPrices_t2, ...
                                                 quotedPrices_t1, 1, ...
                                                 fTerm(outOfSample(1)+i,:), ...
                                                 T(outOfSample(1)+i-1,:), ...
                                                 maturities);
                                             
    realisedPortfolioValueChanges_O(i+1) = portfolioValueNextDay;
end


%% Calculate simulated portfolio value changes given all next day curves for each in sample date and current portfolio value

simulatedPortfolioValues_O_t2 = zeros(size(outOfSample,2) - 1, N);
simulatedPortfolioValues_O_t1 = zeros(size(outOfSample,2) - 1, N);
simPortfolioValueChanges_O =  zeros(size(outOfSample,2) - 1, N);

h = waitbar(0, "Calculate portfolio values");

%Iterera över alla in-sample dagar
for i = 1:size(outOfSample,2)-1
    waitbar((i)/(size(outOfSample,2)-1), h, (i)/(size(outOfSample,2)-1))
    %iterera över alla simulerade kurvor
    for n = 1:N
        forwardPrices_O_t2 = zeros(1,nFXContracts);    
        forwardPrices_O_t1 = zeros(1,nFXContracts);    
        %Iterera över alla kontrakt
        for k = 1:nFXContracts
            %Morgondagens priser, sätt in kurva för morgondagen.
            maturityContractNextDay = (maturitiesFx(outOfSample(1)+i-1,k) - 1/365);
            forwardPrices_O_t2(k) = priceFXswap(exchangeRate(outOfSample(1)+i), ...
                                                fBaseSimulationsAllDays{1,i+1}(:,n)', ...
                                                fTermSimulationsAllDays{1,i+1}(:,n)', ...
                                                fDemandSimulationsAllDays{1,i+1}(:,n)', ...
                                                T(outOfSample(1)+i-1,:), maturityContractNextDay); 
            %Dagens priser, använd dagens kurvor.
            %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, basisRIC, maturitiesRIC, maturity(k,i)); %Sätta in kurva?
        end
        quotedForwardPrices = midPriceFX(outOfSample(1)+i-1,1:nFXContracts);
        maturitiesNextDay = (maturitiesFx(outOfSample(1)+i,1:nFXContracts) - ones(1,nFXContracts)./365);
        simulatedPortfolioValues_O_t2(i+1,n) = portfolioValueFxSwap(forwardPrices_O_t2, ...
                                                                    quotedForwardPrices, 1, ...
                                                                    fTerm(outOfSample(1) + i,:), ...
                                                                    T(outOfSample(1) + i,:), ...
                                                                    maturitiesNextDay);
        %Beräknat Portföljvärde med våra kurvor, ska egentligen bli
        %noll.
        %simulatedPortfolioValues_t1(i, n) = portfolioValueFxSwap(forwardPrices_t1, quotedForwardPrices, N, fTerm, maturitiesFx(i), maturityContracts);

    simPortfolioValueChanges_O(i+1,n) = simulatedPortfolioValues_O_t2(i+1,n); % - simPortfolioValues(i) if initial value is not zero.
    end
end

close(h)

%% Kernel Density Estimation
kernelO = kdeMulti(simPortfolioValueChanges_O(2:end,:), realisedPortfolioValueChanges_O(2:end));


%% Calculate price differences for out of sample

priceDiffsDaily = zeros(1,size(outOfSample,2)-1);
forwardPrice = zeros(size(outOfSample,2)-1);
priceDiffsContracts = zeros(size(outOfSample,2)-1,k);
for i = 1:size(outOfSample,2)-1
   for k = 1:nFXContracts
        forwardPrice(k) = priceFXswap(exchangeRate(outOfSample(1)+i), ...
                                      fBase(outOfSample(1)+i,:), ...
                                      fTerm(outOfSample(1)+i,:), ...
                                      demand(outOfSample(1)+i,:), ...
                                      T(outOfSample(1)+i-1,:), ...
                                      maturitiesFx(outOfSample(1)+i-1,k)); 
        priceDiffsContracts(i,k) = abs(forwardPrice(k) - midPriceFX(outOfSample(1)+i));
   end
    todaysPriceDiff = priceDiffsContracts(i,~isnan(priceDiffsContracts(i,:)));
    priceDiffsDaily(i) = sum(todaysPriceDiff);
                                  
end



