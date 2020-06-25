measurementPath = '.\measurement';

addpath(measurementPath)
fxSwapCalcPrices = [];
fxSwapDemand = []; 
 

c = 1;
if (c==1) % USDSEK
  currencyTerm = 'SEK'; calTerm = 'SWE'; currencyTermTimeZone = 'Europe/Stockholm';
  onTermName = 'STISEKTNDFI='; onTermCal = 'SWE'; onTermCalFixing = 'SWE'; tenorTermON = 'TN'; onTermTimeZone = 'Europe/Stockholm'; settlementLagTermON = 0; onTermBDC = 'M'; onTermDCC = 'MMA0';
  oisTenorStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:SWE PDELAY:0 LFLOAT IDX:OSEKSTI'; % OIS_SEKSTI
  currencyBase = 'USD'; calBase = 'USA'; currencyBaseTimeZone = 'America/New_York';
  curencyTimeZone = 'Europe/Paris';  % Bugg gör att tidszonerna inte fungerar korrekt
  onBaseName = 'USONFFE='; onBaseCal = 'USA'; onBaseCalFixing = 'USA'; tenorBaseON = 'ON'; onBaseTimeZone = 'America/New_York'; settlementLagBaseON = 0; onBaseBDC = 'M'; onBaseDCC = 'MMA0';
  oisBaseStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:USA PDELAY:2 LFLOAT IDX:OFFER'; % 
end
fileName = strcat(currencyBase, currencyTerm);

[pl, ric, times, assetID, assetType, assetTenor, iborTenor, currencies] = populatePortfolioWithData(fileName, currencyTerm, currencyTermTimeZone, onTermTimeZone);

currencyBaseID = mexPortfolio('currencyCode', currencyBase);
currencyTermID = mexPortfolio('currencyCode', currencyTerm);

% indOIS = find(strcmp(tenorON, assetTenor));
% indIBOR = find(strcmp(tenorIRS, iborTenor));

if (c==1) % Currently hardcoded
  nY = 2; % Number of years the curves span
  indBaseOIS = (1:16)';
  indTermOIS = (30:42)';
  indFx = 57;
  indFxSwap = (58:79)';
end
% % Small test example
% indOIS = indOIS(3); 
% indFRA = [];
% indIRS = [];

% accountName = 'Cash';
% cashDCC = 'MMA0';
% cashFrq = '1D';
% cashEom = 'S';
% cashBDC = 'F';
% irStartDate = datenum(2010,1,1);
% cashID = mexPortfolio('createCash', currency, accountName, cashDCC, cashFrq, cashEom, cashBDC, cal, irStartDate);

onID = -1;

tnID = -1;
for i=1:length(indFxSwap)
  tmp = ric{1}{indFxSwap(i)};
  j = strfind(tmp, '=');
  tmp = extractBetween(tmp, j-2, j-1);
  if (length(tmp{1}) == 2)
    if (strcmp(tmp, 'ON'))
      onID = mexPortfolio('addRIC', ric{1}{indFxSwap(i)}, floor(times(1)));
    elseif (strcmp(tmp, 'TN'))
      tnID = mexPortfolio('addRIC', ric{1}{indFxSwap(i)}, floor(times(1)));
    end
  end
end
isHolidayBase = mexPortfolio('isHoliday', onBaseCal, floor(times), floor(times)); % Note that projected holidays may change over time (holidays are added and removed, hence the date when the calendar is defined is important)
isHolidayTerm = mexPortfolio('isHoliday', onTermCal, floor(times), floor(times)); % Note that projected holidays may change over time (holidays are added and removed, hence the date when the calendar is defined is important)
times((isHolidayBase==1) | (isHolidayTerm==1)) = []; % Removes dates when there is a holiday in at least one country


indAll = [indBaseOIS ; indTermOIS ; indFxSwap];
ricAll = [ric{1}(indBaseOIS) ric{1}(indTermOIS) ric{1}(indFxSwap)];

fBaseH = ones(length(times), nY*365+20)*inf;
fTermH = ones(length(times), nY*365+20)*inf;
piH = ones(length(times), nY*365+20)*inf;
zH = ones(length(times), length(indAll))*inf;
yieldH = ones(length(times), length(indAll))*inf;
firstDates = ones(length(times), 1)*inf;
lastDates = ones(length(times), 1)*inf;
piAvgH = ones(length(times), nY*365+20)*inf;

h = waitbar(0, 'Calculating.....');

%setappdata(h,'canceling',0);

% for k=1:length(times)
%start = 2547;
start = 1642;
endVar = length(times);
%start = 19;
%endVar = 25;
for k=start:endVar

  tradeDate = floor(times(k));
%   datestr(tradeDate)

  indOIS = [indBaseOIS ; indTermOIS];
  oisID = zeros(size(indOIS));
  oisPrices = ones(3,length(indOIS))*Inf;
  for j=1:length(indOIS)
    % Retrieve data
    [timeData, data] = mexPortfolio('getValues', assetID{1}(indOIS(j)), times(k), currencyTermTimeZone, {'BID', 'ASK'});
    if (floor(timeData) ~= floor(times(k)))
      data = [NaN ; NaN];
    end
    oisPrices(1:2,j) = data;
    if (sum(isnan(data))==0)
      mid = mean(data);
    else
      mid = 0;
    end
    % Create OIS
    [oisID(j)] = mexPortfolio('initFromGenericInterest', assetID{1}(indOIS(j)), 1, mid, tradeDate, 1);
    
    midAll(k-(start-1),j) = mid;
  end

%   % IBOR
%   % Retrieve data
%   [timeData, data] = mexPortfolio('getValues', assetID{3}(indIBOR), times(k), iborTimeZone, 'LAST');
%   if (floor(timeData) ~= floor(times(k)))
%     iborPrices = [NaN ; NaN ; Inf];
%   else
%     iborPrices = [data ; data ; Inf];
%   end
%   % Create IBOR
%   [iborID] = mexPortfolio('initFromGenericInterest', assetID{3}(indIBOR), 1, data, tradeDate, 0);

  
  % FxSwap
  [timeDataFx, dataFx] = mexPortfolio('getValues', assetID{1}(indFx(1)), times(k), currencyTermTimeZone, {'BID', 'ASK'});
  if (sum(isnan(dataFx))==0)
      exchangeRate = mean(dataFx);
  else
    error('Zero exchange rate');
  end
  settlementDate = mexPortfolio('settlementDate', assetID{1}(indFx(1)), tradeDate);

  fxSwapID = zeros(size(indFxSwap));
  fxSwapPrices = ones(3,length(indFxSwap))*Inf;
  for j=1:length(indFxSwap)
    waitbar((k-start) / (endVar-start), h, sprintf('%12.9f', (k-start) / (endVar-start)))
    % Retrieve data
    [timeData, data] = mexPortfolio('getValues', assetID{1}(indFxSwap(j)), times(k), currencyTermTimeZone, {'BID', 'ASK'});
    if (floor(timeData) ~= floor(times(k)))
      data = [NaN ; NaN];
    end
    midBasisH(k-(start-1),j) = mean(data(1:2,1));
    if (assetID{1}(indFxSwap(j)) == onID) % Bid and ask prices switch order for ON
      onSettlementDate = mexPortfolio('settlementDate', assetID{1}(indFxSwap(j)), tradeDate);
      onMaturityDate = mexPortfolio('maturityDate', assetID{1}(indFxSwap(j)), tradeDate);
      if (onMaturityDate < settlementDate) % Need to add rate for TN      
        [timeData, dataTN] = mexPortfolio('getValues', tnID, times(k), currencyTermTimeZone, {'BID', 'ASK'});
        %midBasisH(k-(start-1),j) = mean(dataTN(1:2,1));
        if (floor(timeData) ~= floor(times(k)))
          dataTN = [NaN ; NaN];
        end
      else % Already at settlement date, no TN
        dataTN = [0 ; 0];        
      end
      fxSwapPrices(1:2,j) = dataFx - data(2:-1:1) - dataTN(2:-1:1);      
    elseif (assetID{1}(indFxSwap(j)) == tnID) % Bid and ask prices switch order for TN
      fxSwapPrices(1:2,j) = dataFx - data(2:-1:1);
    else
      fxSwapPrices(1:2,j) = dataFx + data;
    end
    
    if (sum(isnan(data))==0)
      mid = mean(fxSwapPrices(1:2,j));
    else
      mid = 0;
    end
    midFxH(k-(start-1),j) = mid;
    
    % Create FxSwap
    [fxSwapID(j)] = mexPortfolio('initFromGenericFxSwap', assetID{1}(indFxSwap(j)), 1, exchangeRate, mid, tradeDate, 1);
  end
%   instrID = [oisID ; iborID ; fxSwapID];
%   prices = [oisPrices iborPrices fxSwapPrices];
%   ricInstr = [ric{1}(indOIS) ric{3}(indIBOR) ric{1}(indFxSwap)];

  instrID = [oisID ; fxSwapID];
  prices = [oisPrices fxSwapPrices];
  ricInstr = [ric{1}(indOIS) ric{1}(indFxSwap)];
  indInstr = 1:length(instrID);
  
  
  conTransf = []; % If this is not set, then the standard setting will be applied in computeForwardRates
  conType = ones(size(instrID))*6; % Unique mid price
  useEF = zeros(size(instrID)); % No deviation from market prices
  [instr,removeAsset,flg] = createInstrumentsClasses(pl, tradeDate, ricInstr, instrID, prices, conType);
  instrID(removeAsset) = [];
  indInstr(removeAsset) = [];
  removeFX = removeAsset(length(oisID):end);
  fxSwapID(removeAsset(length(oisID):end)) = [];
  oisID(removeAsset(1:length(oisID))) = [];
  if (~isempty(conTransf))
    conTransf(removeAsset) = [];
  end
  if (~isempty(useEF))
    useEF(removeAsset) = [];
  end

  lastDate = mexPortfolio('lastDate', instrID);
  nF = lastDate-tradeDate;
  fDates = tradeDate:lastDate;
  f = ones(nF,1)*0.05;
  piSpread = ones(nF,1)*0.01;

  if (nF >= 5*365+100)
    error('Currently limitation in the time horizon for solver, change indOIS and indIRS which are set manually')
  end  
 
  snapTime = tradeDate + 0.5; % Sets time of snapshot to 12:00 in local time
  mexPortfolio('initMarketState', snapTime, currencyTermTimeZone);
  fxRateID = mexPortfolio('setMarketStateExchangeRate', snapTime, currencyTermTimeZone, currencyBase, currencyTerm, exchangeRate, settlementDate);
  fBaseCurveID = mexPortfolio('setMarketStateTermStructureDiscreteForward', snapTime, currencyTermTimeZone, currencyBase, 'Interbank', 'ON', 'Forward', fDates, f);
  fTermCurveID = mexPortfolio('setMarketStateTermStructureDiscreteForward', snapTime, currencyTermTimeZone, currencyTerm, 'Interbank', 'ON', 'Forward', fDates, f);
  piCurveID = mexPortfolio('setMarketStateCurrencyDiscreteForward', snapTime, currencyTermTimeZone, currencyBase, currencyTerm, 'SpreadRelative', fDates, piSpread);
  fxCurveID = mexPortfolio('setMarketStateCurrencyRelative', snapTime, currencyTermTimeZone, fBaseCurveID, fTermCurveID, piCurveID, currencyBase, currencyTerm, 'Forward');
  curveIDs = [fBaseCurveID ; fTermCurveID ; piCurveID];
  
	mu = 1.0;
	nIterations = 200;
    iterationPrint = 0; %true;
	checkKKT = 1; %true;
	precision = 1E-6;
	precisionEq = 1E-10;
	precisionKKT = 1E-10;
	muDecrease = 0.1;
	maxRelStep = 0.99;
  mexPortfolio('initBlomvallNdengo', mu, nIterations, iterationPrint, checkKKT, precision, precisionEq, precisionKKT, muDecrease, maxRelStep, curveIDs);

  for i=1:length(curveIDs)
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'ul', -1*ones(nF,1));
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'uu', ones(nF,1)*inf);

    T = (fDates-fDates(1))/365;
    knowledgeHorizon = 2; % After knowledge horizon second order derivative have equal weight
    informationDecrease = 2; % How much information decrease in one year
    cb = ones(nF-1,1);
    cb(1:min(round(knowledgeHorizon*365),length(cb))) = exp((T(1:min(round(knowledgeHorizon*365),length(cb)))-knowledgeHorizon)*log(informationDecrease^2));
    cb(1) = 0;
    cb = cb*10; 
   
    if (curveIDs(i) == piCurveID) % Example where a jump for one specific day is allowed
        cpL = zeros(nF-1,1);
        cpL(1) = 0;
        nJumps = 0;
        jumpDates = zeros(1,1);
        cb = cb/10;
        
        % Set jump at specified dates
        for n = 1:nF-1
            [yTrade,mTrade,~] = datevec(tradeDate);
            [y,m,d] = datevec(fDates(n));
           %if d == 1 && (m == 1 || m == 4 || m == 7 || m == 10) %&& y == 2017 % m = 1 % && (m == 1 || m == 4 || m == 7 || m == 10) 
           if d == 1 && m == 1 && y == yTrade+1 %&& mTrade == 12 % m == 1 % && (m == 1 || m == 4 || m == 7 || m == 10) 
                nJumps = nJumps + 1 ;
                if n == 1
                    jumpDate = 1;
                else
                    jumpDate = n-1;
                end
                if fDates(jumpDate) ~= tradeDate
                     cb(jumpDate) = 0; cb(jumpDate-1) = 0; cb(jumpDate+1) = 0; 
                else 
                    cb(jumpDate) = 0; cb(jumpDate+1) = 0;
                end
                cpL(jumpDate) = 100;
                jumpDates(nJumps) = jumpDate;    
            end
        end
      mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'cpL', cpL); % First order derivative, with middle removed (central difference)
    end    
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'cb', cb); % Second order derivative
  end
  
  zType = ones(size(instrID))*2; % Type of z-variable (price error): 0 = None, 1 = Price, 2 = Parallel shift of forward rates or spread from forward rates
  E = ones(size(instrID))*100;
  for i=1:length(instrID)
    at = mexPortfolio('assetType', instrID(i));
    if (at == pl.atOIS)
      curr = mexPortfolio('assetCurrency', instrID(i));
      if (curr == currencyBaseID)
        E(i) = 100;
      elseif (curr == currencyTermID)
        E(i) = 100;
      else
        error('Incorrect currency')
      end
    elseif (at == pl.atFxSwap)      
      E(i) = 1;     
    end
  end
     % Objective function coefficient for z-variables
  F = ones(size(instrID));       % Scaling of z-variable in constraint

  [u, P, z, theoreticalP, marketQuoteP] = mexPortfolio('solveBlomvallNdengo', instrID, zType, E, F);
  
  fBase = u(1:nF);
  fTerm = u(nF+1:2*nF);
  pi = u(2*nF+1:end);
  fBaseH(k-(start-1),1:nF) = fBase;
  fTermH(k-(start-1),1:nF) = fTerm;
  piH(k-(start-1),1:nF) = pi;
  zH(k-(start-1),indInstr) = abs(z);
  for j=1:length(instrID)
    yieldH(k-(start-1),indInstr(j)) = instr.data{j}.price(3);
  end
  firstDates(k) = tradeDate;
  lastDates(k) = lastDate;

  instrT = zeros(size(z));
  instrf = zeros(size(z));
  for i=1:length(instrID)
    at = mexPortfolio('assetType', instrID(i));
    n = instr.maturityDate(i)-tradeDate;
    instrT(i) = T(n);
    if (at == pl.atOIS)
      curr = mexPortfolio('assetCurrency', instrID(i));
      if (curr == currencyBaseID)
        instrf(i) = fBase(n);
      elseif (curr == currencyTermID)
        instrf(i) = fTerm(n);
      else
        error('Incorrect currency')
      end
    elseif (at == pl.atFxSwap)      
      instrf(i) = pi(n);      
    end
  end
  
  plot(T(1:end-1), fBase, T(1:end-1), fTerm, T(1:end-1), pi, instrT, instrf+z, 'o');
  [~, eInd] = sort(abs(z), 1, 'descend');
  for j = 1:min(3,length(eInd))
    text(instrT(eInd(j))+0.1, instrf(eInd(j)) + z(eInd(j)), instr.assetRIC(eInd(j)));
  end   
  title(datestr(times(k)));
  
  for j=1:length(instrID)
%     [timeData, data] = mexPortfolio('getValues', instrID(j), times(k), currencyTermTimeZone, {'BID', 'ASK'});
    fprintf('%3d %18s %12s %9f %9f %9f %9f\n',j,instr.assetRIC{j}, datestr(instr.maturityDate(j)), instr.data{j}.price(3), z(j), theoreticalP(j), marketQuoteP(j));    
  end
  
  % Statistics from pricing error z. Saved for every k
  zAll(k-(start-1),:) = [mean(abs(z)), median(abs(z)), min(abs(z)),  max(abs(z))];
  
  % Calculate fxSwap prices 
  x = 1;
  for i=length(oisID)+1:length(instrID)
      if removeFX(i-(length(oisID)-1))
        x = x + 1;
      end
      n = instr.maturityDate(i)-tradeDate;
      fxSwapCalcPrices(k-(start-1), x) = exchangeRate*exp((fTerm(1:n)./365)'*T(2:n+1)'-(fBase(1:n)./365)'*T(2:n+1)'+(pi(1:n)./365)'*T(2:n+1)');
      fxSwapDemand(k-(start-1), x) = fxSwapCalcPrices(k-(start-1),x)-exchangeRate;
      x = x + 1;
  end
  tradeDates(k-(start-1),:) = datevec(tradeDate);
%  pause;  
%  pause(0.5);
 % stop
  if getappdata(h, 'cancel_callback') == 1
    errordlg('Stopped by user');    
    break;  % Stop the loop: for coubnt_i1
  end
  
 % Modify demand curve for PCA
 piAvgH(k-(start-1),1:nF) = pi;
 for i = 1:length(jumpDates)
     if jumpDates(i) == 1
         piAvgH(k-(start-1),jumpDates(i)) = pi(jumpDates(i)+1);
     else
        piAvgH(k-(start-1),jumpDates(i)) = (pi(jumpDates(i)-1) + pi(jumpDates(i)+1))/2;
     end
 end
end

rmpath(measurementPath)
close(h)
% hej
%%
filename = 'X:\Examensarbete\DataTest\test_14-20_FX1-Base100-Term1000_cb10_';
writematrix(fBaseH, strcat(filename,'fBase.csv'), 'Delimiter', 'semi')
writematrix(fTermH, strcat(filename,'fTerm.csv'), 'Delimiter', 'semi')
writematrix(piH, strcat(filename,'pi.csv'), 'Delimiter', 'semi')
writematrix(piAvgH, strcat(filename,'piAvg.csv'), 'Delimiter', 'semi')


