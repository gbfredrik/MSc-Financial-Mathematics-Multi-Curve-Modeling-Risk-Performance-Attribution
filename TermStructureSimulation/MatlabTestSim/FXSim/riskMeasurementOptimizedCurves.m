%% OPTIMISED CURVES

%Get all realised and simulated curves

%% Calculate realised portfolio value changes for all in-date samples

RealisedPortfolioValueChanges = zeros(1,size(outOfSample,2) - 1);

for i = 1:size(outOfSample,2)-1
     forwardPrices_t2 = zeros(1,nFXContracts);  
     for k = 1:nFXContracts
        %Morgondagens priser, sätt in realiserade kurvor för morgondagen.
        maturityContractNextDay = (maturitiesFx(outOfSample(1)+i,k) - 1/365);
        forwardPrices_t2(k) = priceFXswap(midPrice(outOfSample(1)+i,k), ...
                                          fBase(outOfSample(1)+i,:), ...
                                          fTerm(outOfSample(1)+i,:), ...
                                          pi(outOfSample(1)+i,:), ...
                                          TT(outOfSample(1)+1,:), ...
                                          maturityContractNextDay); 
        %Dagens priser, använd dagens kurvor.
        %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, basisRIC, maturitiesRIC, maturity(k,i)); %Sätta in kurva?
     end

    quotedPrices_t1 = midPrice(outOfSample(1)+i-1,1:nFXContracts);
    
    maturities = (maturitiesFx(outOfSample(1)+i,1:nFXContracts) - ones(1,nFXContracts)./365);
    
    %Portfolio value is zero at trade date so portfolio change from the day before is portfolio
    %value of tomorrow.
    portfolioValueNextDay = portfolioValueFxSwap(forwardPrices_t2, quotedPrices_t1, 1, fTerm(outOfSample(1)+i,:), T(outOfSample(1)+1,:), maturities);
    RealisedPortfolioValueChanges(i+1) = portfolioValueNextDay;
end


%% Calculate simulated portfolio value changes given all next day curves for each in sample date and current portfolio value

simulatedPortfolioValues_O_t2 = zeros(size(outOfSample,2) - 1, N);
simulatedPortfolioValues_O_t1 = zeros(size(outOfSample,2) - 1, N);
simPortfolioValueChanges_O =  zeros(size(outOfSample,2) - 1, N);


%Iterera över alla in-sample dagar
for i = 1:size(outOfSample,2)-1
        %iterera över alla simulerade kurvor
        for n = 1:N
            forwardPrices_O_t2 = zeros(1,nFXContracts);    
            forwardPrices_O_t1 = zeros(1,nFXContracts);    
            %Iterera över alla kontrakt
            for k = 1:nFXContracts
                %Morgondagens priser, sätt in kurva för morgondagen.
                maturityContractNextDay = (maturitiesFx(outOfSample(1)+i,k) - 1/365);
                simulatedCurve = [0 basisSimulationsAllDays{1,i+1}(n,:)];
                forwardPrices_O_t2(k) = priceFXswap(midPrice(outOfSample(1)+i,k), fBase(n,:), fTerm(n,:), pi(n,:), T, maturityContractNextDay); 
                %Dagens priser, använd dagens kurvor.
                %forwardPrices_t1(k) = priceFXswapRawInterpol(midPrice(inSample(1)+i-1, basisRIC, maturitiesRIC, maturity(k,i)); %Sätta in kurva?
            end
            quotedForwardPrices = midPrice(outOfSample(1)+i-1,1:nFXContracts);
            maturitiesNextDay = [0 (maturitiesFx(outOfSample(1)+i,1:nFXContracts) - ones(1,nFXContracts)./365)];
            simulatedPortfolioValues_O_t2(i+1,n) = portfolioValueFxSwap(forwardPrices_O_t2, quotedForwardPrices, 1, fTerm(outOfSample(1)+i,:), T(outOfSample(1)+i,:), maturitiesNextDay);
            %Beräknat Portföljvärde med våra kurvor, ska egentligen bli
            %noll.
            %simulatedPortfolioValues_t1(i, n) = portfolioValueFxSwap(forwardPrices_t1, quotedForwardPrices, N, fTerm, maturitiesFx(i), maturityContracts);
        
        simPortfolioValueChanges_O(i+1,n) = simulatedPortfolioValues_O_t2(i+1,n); % - simPortfolioValues(i) if initial value is not zero.
        end

        
        
end


