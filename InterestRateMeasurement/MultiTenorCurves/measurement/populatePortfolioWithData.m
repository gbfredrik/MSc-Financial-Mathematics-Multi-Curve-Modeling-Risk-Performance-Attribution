function [pl, ric, times, assetID, assetType, assetTenor, iborTenor, currencies] = populatePortfolioWithData(fileName, currency, currencyTimeZone, iborTimeZone)

mexPortfolio('clear')
mexPortfolio('init', currency)

% Create portfolio lexicon (pl) for constants in portfolio
pl.atFuture = mexPortfolio('assetTypeCode', 'FutureIR');
pl.atDeposit = mexPortfolio('assetTypeCode', 'Deposit');
pl.atIBOR = mexPortfolio('assetTypeCode', 'IBOR');
pl.atFRAG = mexPortfolio('assetTypeCode', 'FRAG');
pl.atFRA  = mexPortfolio('assetTypeCode', 'FRA');
pl.atIRSG = mexPortfolio('assetTypeCode', 'IRSG');
pl.atIRS  = mexPortfolio('assetTypeCode', 'IRS');
pl.atOISG = mexPortfolio('assetTypeCode', 'OISG');
pl.atOIS  = mexPortfolio('assetTypeCode', 'OIS');
pl.atFxSwapG = mexPortfolio('assetTypeCode', 'FxSwapG');
pl.atFxSwap  = mexPortfolio('assetTypeCode', 'FxSwap');
pl.atFX  = mexPortfolio('assetTypeCode', 'FX');



load([fileName '.mat']); % File created by createMatFromExcel

assetID = cell(1, length(data));
assetType = cell(1, length(data));
for k=1:length(data)
  assetID{k} = zeros(length(ric{k}),1);
  assetType{k} = zeros(length(ric{k}),1);
  for j=1:length(ric{k})
    [assetID{k}(j), assetType{k}(j)] = mexPortfolio('addRIC', ric{k}{j}, dates{k}(1));
  end
end

assetTenor = cell(length(ric{1}), 1);
iborTenor = cell(length(ric{3}), 1);
currencies = ones(length(ric{1}), 2)*inf;

for j=1:length(ric{1})
  [iborID] = mexPortfolio('underlyingIBOR', assetID{1}(j));
  if (~isempty(iborID))
    [assetTenor{j}] = mexPortfolio('tenorIBOR', iborID);
  else
    [assetTenor{j}] = '';    
  end
  [cur] = mexPortfolio('assetCurrency', assetID{1}(j));
  if (length(cur) == 2)
    currencies(j,:) = cur;
  elseif (length(cur) == 1)
    currencies(j,1) = cur;
  else
    error('Missing currency for asset');
  end
  %fprintf('%18s %3d %3d %3d %3s\n', ric{1}{j}, assetID{1}(j), assetType{1}(j), iborID, assetTenor{j});
end
for j=1:length(ric{3})
  [iborID] = mexPortfolio('underlyingIBOR', assetID{3}(j)); 
  [iborTenor{j}] = mexPortfolio('tenorIBOR', iborID);
  %fprintf('%18s %3d %3d %3d %3s\n', ric{3}{j}, assetID{3}(j), assetType{3}(j), iborID, iborTenor{j});
end

% Store IRS/FRA/OIS data
if (sum(abs(dates{1}-dates{2})) > 0)
  % Need to synchronize data
  datesAll = unique([dates{1} ; dates{2}]);
  datesMin = datesAll(1);
  datesMax = datesAll(end);
  datesRangeN = datesMax-datesMin+1;
  ind = zeros(datesRangeN,1);
  ind(datesAll-datesMin+1) = 1:length(datesAll);
  for j=1:2
    tmp = data{j};
    data{j} = ones(length(datesAll), size(tmp,2))*NaN;
    data{j}(ind(dates{j}-datesMin+1),:) = tmp;
    dates{j} = datesAll;
  end
end

tmp = datevec(dates{1});
tmp(:,4) = 22; % Assumes that daily prices are set at 22:00:00 (ensures that IBOR has been set also in JPY)
times = datenum(datetime(tmp));

for j=1:length(ric{1})
  if (assetID{1}(j) ~= assetID{2}(j))
    error('Assumes that BID/ASK exist for both assets');
  end
  if (assetType{1}(j) == pl.atFxSwapG)
    v = [data{1}(:,j) data{2}(:,j)]/10000;
  elseif (assetType{1}(j) == pl.atFX)
    v = [data{1}(:,j) data{2}(:,j)];
  else
    v = [data{1}(:,j) data{2}(:,j)]/100;
  end
  tmp = isnan(v);
  ind = tmp(:,1) & tmp(:,2);
  timesData = times;
  timesData(ind) = [];
  v(ind,:) = [];

  mexPortfolio('setHistory', assetID{1}(j), timesData, currencyTimeZone, {'BID', 'ASK'}, v);
end

% Store IBOR data

tmp = datevec(dates{3});
tmp(:,4) = 11; % Assumes that IBORS are set at 11:00:00
timesIBOR = datenum(datetime(tmp));

for j=1:length(ric{3})
  v = data{3}(:,j)/100;
  ind = isnan(v);
  timesData = timesIBOR;
  timesData(ind) = [];
  v(ind,:) = [];

%   mexPortfolio('setHistory', assetID{3}(j), timesIBOR, iborTimeZone, 'LAST', v);
  mexPortfolio('setHistory', assetID{3}(j), timesData, iborTimeZone, 'LAST', v);
end