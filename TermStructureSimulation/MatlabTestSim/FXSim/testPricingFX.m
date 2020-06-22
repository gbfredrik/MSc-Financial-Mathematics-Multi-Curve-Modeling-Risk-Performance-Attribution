measurementPath = 'X:\Examensarbete\Data\';
filename = 'FX_';
contract = 'SEK_';
fBase = readmatrix(strcat(measurementPath,filename,contract,'base.csv'),'Delimiter',';'); %, 'NumHeaderLines',1);
fTerm = readmatrix(strcat(measurementPath,filename,contract,'term.csv'),'Delimiter',';');
pi = readmatrix(strcat(measurementPath,filename,contract,'fx.csv'),'Delimiter',';');
piAvg = readmatrix(strcat(measurementPath,filename,contract,'fxAvg.csv'),'Delimiter',';');
basis = readmatrix(strcat(measurementPath,filename,contract,'basis.csv'),'Delimiter',';');
maturitiesFx = readmatrix(strcat(measurementPath,filename,contract,'maturityDateFX.csv'),'Delimiter',';');
T = readmatrix(strcat(measurementPath,filename,contract,'T.csv'),'Delimiter',';');
exchangeRate = readmatrix(strcat(measurementPath,filename,contract,'exchangeRateMid.csv'),'Delimiter',';');
midPrice = readmatrix(strcat(measurementPath,filename,contract,'midPrices.csv'),'Delimiter',';');
maturitiesBase = readmatrix(strcat(measurementPath,filename,contract,'maturitiesBase.csv'),'Delimiter',';');
maturitiesTerm = readmatrix(strcat(measurementPath,filename,contract,'maturitiesTerm.csv'),'Delimiter',';');
midPriceTerm = readmatrix(strcat(measurementPath,filename,contract,'midPriceTerm.csv'),'Delimiter',';');

[nTradeDates, ~] = size(maturitiesFx); 
%basis = [zeros(nTradeDates, 1) basis];
maturitiesFx = [zeros(nTradeDates, 1) maturitiesFx];
maturitiesBase = [zeros(nTradeDates, 1) maturitiesBase];
maturitiesTerm = [zeros(nTradeDates, 1) maturitiesTerm];
midPriceTerm = [zeros(nTradeDates, 1) midPriceTerm];

%% Pricing FX swaps - optimization model
[nTradeDates, nContractDates] = size(fBase); 
[~,nRICs] = size(basis);

prices = zeros(nTradeDates, nContractDates);

% Calculate the price on a specified time points
for tradeDate = 1:nTradeDates
    for contractDate = 1:nContractDates
        prices(tradeDate,contractDate) = priceFXswap(exchangeRate(tradeDate), fBase(tradeDate,:), fTerm(tradeDate,:), pi(tradeDate,:), T(tradeDate,:), T(tradeDate,contractDate));
    end
end


%% Pricing Fx swaps - Raw Interpolation model
[nTradeDates, nRICs] = size(basis); 
[~,nDates] = size(T);

pricesRi = zeros(nTradeDates, nContractDates);

% Calculate the price on a specified time points
for tradeDate = 1:nTradeDates
    for contractDate = 1:nDates
        pricesRi(tradeDate,contractDate) = priceFXswapRawInterpol(exchangeRate(tradeDate), basis(tradeDate,:), maturitiesFx(tradeDate,:), T(tradeDate,contractDate));
    end
end

%% Curve Fx swaps - Raw interpolation model
[nTradeDates, nRICs] = size(maturitiesFx);
 
[~, curveLength] = size(fTerm);
curveRi = ones(nTradeDates, curveLength)*inf;

% Calculate the basis curve
for tradeDates = 1:nTradeDates
    curveRi(tradeDates,:) = curveFxSwapRawInterpol(maturitiesFx(tradeDates,:), basis(tradeDates,:), curveLength);
    plot(1:curveLength, curveRi(tradeDates,:))
    hold on
end
hold off

%% Curve OIS term - Raw interpolation model
[nTradeDates, nRICs] = size(midPriceTerm);
 
[~, curveLength] = size(fTerm);
curveRiTerm = ones(nTradeDates, curveLength)*inf;

% Calculate the OIS term curve
for tradeDates = 1:nTradeDates
    curveRiTerm(tradeDates,:) = curveFxSwapRawInterpol(maturitiesTerm(tradeDates,:), midPriceTerm(tradeDates,:), curveLength);
    plot(1:curveLength, curveRiTerm(tradeDates,:))
    hold on
end
hold off

%% Portfolio valuation - optimization model
[nTradeDates, nContractDates] = size(fTerm); 
[~,nRICs] = size(basis);

forwardPricesPV = zeros(nTradeDates, nRICs);
portfolioValue = zeros(nTradeDates-1, 1);

N = ones(nTradeDates-1, nRICs);


for tradeDate = 1:nTradeDates
    for j = 1:nRICs
        forwardPricesPV(tradeDate,j) = priceFXswap(exchangeRate(tradeDate), fBase(tradeDate,:), fTerm(tradeDate,:), pi(tradeDate,:), T(tradeDate,:), maturitiesFx(tradeDate, j));
        %prices(tradeDate,contractDate) = priceFXswap(exchangeRate(tradeDate), fBase(tradeDate,:), fTerm(tradeDate,:), pi(tradeDate,:), T(tradeDate,:), T(tradeDate,contractDate));
    end
end

for tradeDate = 2:nTradeDates
    portfolioValue(tradeDate-1) = portfolioValueFxSwap(forwardPricesPV(tradeDate,:), forwardPricesPV(tradeDate-1,:), N(tradeDate-1,:), fTerm(tradeDate,:), T(tradeDate,:), maturitiesFx(tradeDate,:));  
end

%% Portfolio valuation - Raw interpolation
[nTradeDates, nRICs] = size(basis); 
[~,nContractDates] = size(T);

N = ones(nTradeDates-1, nRICs);

forwardPricesRiPV = zeros(nTradeDates, nRICs);
curveRiTerm = ones(nTradeDates, nContractDates)*inf;
portfolioValueRi = zeros(nTradeDates-1, 1);
git 
% Calculate forward prices and 1-day disc for term interest rate
for tradeDate = 1:nTradeDates
    for j = 1:nRICs
        forwardPricesRiPV(tradeDate, j) = priceFXswapRawInterpol(exchangeRate(tradeDate), basis(tradeDate,:), maturitiesFx(tradeDate,:), maturitiesFx(tradeDate,j));
    end
    curveRiTerm(tradeDates,:) = curveFxSwapRawInterpol(maturitiesTerm(tradeDates,:), midPriceTerm(tradeDates,:), nContractDates);
end
  
% Calculate portfolio value
for tradeDate = 2:nTradeDates
    portfolioValueRi(tradeDate-1) = portfolioValueFxSwap(forwardPricesRiPV(tradeDate,:), forwardPricesRiPV(tradeDate-1,:), N(tradeDate-1,:), curveRiTerm(tradeDate,:), T(tradeDate,:), maturitiesFx(tradeDate,:));  
end



