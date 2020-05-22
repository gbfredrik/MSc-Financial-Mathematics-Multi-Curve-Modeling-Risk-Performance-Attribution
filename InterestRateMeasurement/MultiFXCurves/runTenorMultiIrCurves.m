measurementPath = 'Z:\prog\XiFinPortfolio\matlab\measurement';
% measurementPath = '/home/jorbl45/axel/jorbl45/prog/XiFinPortfolio/matlab/measurement';

addpath(measurementPath)

c = 7;
if (c==1) % CHF
  currency = 'CHF'; cal = 'SWI'; currencyTimeZone = 'Europe/Paris';
  iborName = 'LIBORCHF3M'; iborCal = 'SWI,UKG'; iborCalFixing = 'UKG'; tenorIRS = '6M'; iborTimeZone = 'Europe/London'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA0';
  onName = 'SARON.S'; onCal = 'SWI'; onCalFixing = 'SWI'; tenorON = 'ON'; onTimeZone = 'Europe/Paris'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA0';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:2WD DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:SWI,UKG PDELAY:0 LFIXED FRQ:1 CCM:BB00 LFLOAT FRQ:2 CCM:MMA0'; % CHF_AB6L
  oisStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:SWI PDELAY:0 LFLOAT IDX:OCHFTOIS'; % OIS_CHFTOIS
  timeToMaturityFRA = '6X12'; fraDCC = 'MMA0';
elseif (c==2) % EUR
  currency = 'EUR'; cal = 'EMU'; currencyTimeZone = 'Europe/Paris';
  iborName = 'EURIBOR6M'; iborCal = 'EMU'; iborCalFixing = 'EUR'; tenorIRS = '6M'; iborTimeZone = 'Europe/Paris'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA0';
  onName = 'EONIA='; onCal = 'EMU'; onCalFixing = 'EMU'; tenorON = 'ON'; onTimeZone = 'Europe/Paris'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA0';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:2WD DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:EMU PDELAY:0 LFIXED FRQ:1 CCM:BB00 LFLOAT FRQ:2 CCM:MMA0'; % EUR_AB6E
  oisStructure = 'PAID:FIXED LBOTH SETTLE:0WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:EMU PDELAY:1 LFLOAT IDX:OEONIA'; % OIS_EONIA
  timeToMaturityFRA = '6X12'; fraDCC = 'MMA0';
elseif (c==3) % GBP, Check yearly time frac
  currency = 'GBP'; cal = 'UKG'; currencyTimeZone = 'Europe/London';
  iborName = 'LIBORGBP6M'; iborCal = 'UKG'; iborCalFixing = 'UKG'; tenorIRS = '6M'; iborTimeZone = 'Europe/London'; settlementLagIRS = 0; irsBDC = 'M'; iborDCC = 'MMA5';
  onName = 'SONIAOSR='; onCal = 'UKG'; onCalFixing = 'UKG'; tenorON = 'ON'; onTimeZone = 'Europe/London'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA5';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:0WD FRQ:2 DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:UKG PDELAY:0 LFIXED CCM:BBA5 LFLOAT CCM:MMA5'; % GBP_SB6L (fix)
  oisStructure = 'PAID:FIXED LBOTH SETTLE:0WD FRQ:1 CCM:MMA5 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:EMU PDELAY:0 LFLOAT IDX:OSONIA'; % OIS_SONIA
  timeToMaturityFRA = '6X12'; fraDCC = 'MMA5';
elseif (c==4) % JPY
  currency = 'JPY'; cal = 'JAP'; currencyTimeZone = 'Asia/Tokyo';
  iborName = 'LIBORJPY6M'; iborCal = 'JAP,UKG'; iborCalFixing = 'UKG'; tenorIRS = '6M'; iborTimeZone = 'Europe/London'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA0';
  onName = 'JPONMU=RR'; onCal = 'JAP'; onCalFixing = 'JAP'; tenorON = 'ON'; onTimeZone = 'Asia/Tokyo'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA5';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:2 DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:JAP,UKG PDELAY:0 LFIXED CCM:BBA5 LFLOAT CCM:MMA5'; % JPY_SB6L (fix)
  oisStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA5 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:JAP PDELAY:2 LFLOAT IDX:OJPYONMU'; % OIS_JPYONMU
  timeToMaturityFRA = '6X12'; fraDCC = 'MMA0';
elseif (c==5) % KRW
  currency = 'KRW'; cal = 'KOR'; currencyTimeZone = 'Asia/Seoul';
  iborName = 'KRWCD3M'; iborCal = 'KOR'; iborCalFixing = 'KOR'; tenorIRS = '3M'; iborTimeZone = 'Asia/Seoul'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA5';
  clear onName; % Do hot have OIS
  irsStructure = 'PAID:FIXED LBOTH SETTLE:1WD CCM:MMA5 DMC:MOD EMC:S CFADJ:YES REFDATE:MATURITY PDELAY:0 CLDR:KOR LFIXED FRQ:4 LFLOAT FRQ:4'; % KRW_QMCD
  clear oisStructure; % Do hot have OIS
  timeToMaturityFRA = '3X6'; fraDCC = 'MMA5';
elseif (c==6) % SEK
  currency = 'SEK'; cal = 'SWE'; currencyTimeZone = 'Europe/Stockholm';
  iborName = 'STIBOR3M'; iborCal = 'SWE'; iborCalFixing = 'SWE'; tenorIRS = '3M'; iborTimeZone = 'Europe/Stockholm'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA0';
  onName = 'STISEKTNDFI='; onCal = 'SWE'; onCalFixing = 'SWE'; tenorON = 'TN'; onTimeZone = 'Europe/Stockholm'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA0';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:2WD DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:SWE PDELAY:0 LFIXED FRQ:1 CCM:BB00 LFLOAT FRQ:4 CCM:MMA0'; % SEK_AB3S
  oisStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:SWE PDELAY:0 LFLOAT IDX:OSEKSTI'; % OIS_SEKSTI
  timeToMaturityFRA = '3F1'; fraDCC = 'MMA0';
elseif (c==7) % USD
  currency = 'USD'; cal = 'USA'; currencyTimeZone = 'America/New_York';
% curencyTimeZone = 'Europe/Paris';  % Bugg gör att tidszonerna inte fungerar korrekt
  iborName = 'LIBORUSD3M'; iborCal = 'USA,UKG'; iborCalFixing = 'UKG'; tenorIRS = '3M'; iborTimeZone = 'Europe/London'; settlementLagIRS = 2; irsBDC = 'M'; iborDCC = 'MMA0';
  onName = 'USONFFE='; onCal = 'USA'; onCalFixing = 'USA'; tenorON = 'ON'; onTimeZone = 'America/New_York'; settlementLagON = 0; onBDC = 'M'; onDCC = 'MMA0';
  irsStructure = 'PAID:FIXED LBOTH SETTLE:2WD CCM:MMA0 DMC:M EMC:S CFADJ:YES REFDATE:MATURITY CLDR:USA PDELAY:0 LFIXED FRQ:1 LFLOAT FRQ:4'; % USD_AM3L
  oisStructure = 'PAID:FIXED LBOTH SETTLE:2WD FRQ:1 CCM:MMA0 DMC:FOL EMC:SAME CFADJ:NO REFDATE:MATURITY CLDR:USA PDELAY:2 LFLOAT IDX:OFFER'; % 
  timeToMaturityFRA = '3X6'; fraDCC = 'MMA0';
end
fileName = currency;

[pl, ric, times, assetID, assetType, assetTenor, iborTenor] = populatePortfolioWithData(fileName, currencyTerm, currencyTermTimeZone, onTermTimeZone);


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


indOIS = find(strcmp(tenorON, assetTenor));
indFRA = find(strcmp(tenorIRS, assetTenor) &  pl.atFRAG == assetType{1});
indIRS = find(strcmp(tenorIRS, assetTenor) &  pl.atIRSG == assetType{1});
indIBOR = find(strcmp(tenorIRS, iborTenor));

if (c==6)
  indOIS = indOIS(1:16); % Only keep up to five years (currently to slow solver)
  indIRS = indIRS(1:4);    % Only keep up to five years (currently to slow solver)
elseif (c==7)
  indOIS = indOIS(1:16); % Only keep up to two years (currently to slow solver)
  indIRS = indIRS(1);    % Only keep up to two years (currently to slow solver)
end
% % Small test example
% indOIS = indOIS(3); 
% indFRA = [];
% indIRS = [];

accountName = 'Cash';
cashDCC = 'MMA0';
cashFrq = '1D';
cashEom = 'S';
cashBDC = 'F';
irStartDate = datenum(2010,1,1);
cashID = mexPortfolio('createCash', currency, accountName, cashDCC, cashFrq, cashEom, cashBDC, cal, irStartDate);

% for k=length(times):length(times)
for k=1:length(times)
  tradeDate = floor(times(k));
%   datestr(tradeDate)
  oisID = zeros(size(indOIS));
  oisPrices = ones(3,length(indOIS))*Inf;
  for j=1:length(indOIS)
    % Retrieve data
    [timeData, data] = mexPortfolio('getValues', assetID{1}(indOIS(j)), times(k), currencyTimeZone, {'BID', 'ASK'});
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
  end

  % IBOR
  % Retrieve data
  [timeData, data] = mexPortfolio('getValues', assetID{3}(indIBOR), times(k), iborTimeZone, 'LAST');
  if (floor(timeData) ~= floor(times(k)))
    iborPrices = [NaN ; NaN ; Inf];
  else
    iborPrices = [data ; data ; Inf];
  end
  % Create IBOR
  [iborID] = mexPortfolio('initFromGenericInterest', assetID{3}(indIBOR), 1, data, tradeDate, 0);

  
  % FRA
  fraID = zeros(size(indFRA));
  fraPrices = ones(3,length(indFRA))*Inf;
  for j=1:length(indFRA)
    % Retrieve data
    [timeData, data] = mexPortfolio('getValues', assetID{1}(indFRA(j)), times(k), currencyTimeZone, {'BID', 'ASK'});
    if (floor(timeData) ~= floor(times(k)))
      data = [NaN ; NaN];
    end
    fraPrices(1:2,j) = data;
    if (sum(isnan(data))==0)
      mid = mean(data);
    else
      mid = 0;
    end
    % Create FRA
    [fraID(j)] = mexPortfolio('initFromGenericInterest', assetID{1}(indFRA(j)), 1, mid, tradeDate, 1);
  end
  % IRS
  irsID = zeros(size(indIRS));
  irsPrices = ones(3,length(indIRS))*Inf;
  for j=1:length(indIRS)
    % Retrieve data
    [timeData, data] = mexPortfolio('getValues', assetID{1}(indIRS(j)), times(k), currencyTimeZone, {'BID', 'ASK'});
    if (floor(timeData) ~= floor(times(k)))
      data = [NaN ; NaN];
    end
    irsPrices(1:2,j) = data;
    if (sum(isnan(data))==0)
      mid = mean(data);
    else
      mid = 0;
    end
    % Create IRS
    [irsID(j)] = mexPortfolio('initFromGenericInterest', assetID{1}(indIRS(j)), 1, mid, tradeDate, 1);
  end
  instrID = [oisID ; iborID ; fraID ; irsID];
  prices = [oisPrices iborPrices fraPrices irsPrices];
  ricInstr = [ric{1}(indOIS) ric{3}(indIBOR) ric{1}(indFRA) ric{1}(indIRS)];

  conTransf = []; % If this is not set, then the standard setting will be applied in computeForwardRates
  conType = ones(size(instrID))*6; % Unique mid price
  useEF = zeros(size(instrID)); % No deviation from market prices
  [instr,removeAsset,flg] = createInstrumentsClasses(pl, tradeDate, ricInstr, instrID, prices, conType);
  instrID(removeAsset) = [];
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
  
  if (k~=1)
    mexPortfolio('clearMarketState');
  end
  mexPortfolio('initMarketState', tradeDate, currencyTimeZone);
  fCurveID = mexPortfolio('setMarketStateTermStructureDiscreteForward', tradeDate, currencyTimeZone, currency, 'Interbank', 'ON', 'Forward', fDates, f);
  piCurveID = mexPortfolio('setMarketStateTermStructureDiscreteForward', tradeDate, currencyTimeZone, currency, 'Interbank', tenorIRS, 'SpreadRelative', fDates, piSpread);
  iborCurveID = mexPortfolio('setMarketStateTermStructureRelative', tradeDate, currencyTimeZone, fCurveID, piCurveID, currency, 'Interbank', tenorIRS, 'Forward');
  
	mu = 1.0;
	nIterations = 200;
  iterationPrint = 1; %true;
	checkKKT = 1; %true;
	precision = 1E-6;
	precisionEq = 1E-10;
	precisionKKT = 1E-10;
	muDecrease = 0.1;
	maxRelStep = 0.99;
  curveIDs = [fCurveID ; piCurveID];
  mexPortfolio('initBlomvallNdengo', mu, nIterations, iterationPrint, checkKKT, precision, precisionEq, precisionKKT, muDecrease, maxRelStep, curveIDs);

  for i=1:length(curveIDs)
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'ul', zeros(nF,1));
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'uu', ones(nF,1)*inf);

    T = (fDates-fDates(1))/365;
    knowledgeHorizon = 2; % After knowledge horizon second order derivative have equal weight
    informationDecrease = 2; % How much information decrease in one year
    cb = ones(nF-1,1);
    cb(1:min(round(knowledgeHorizon*365),length(cb))) = exp((T(1:min(round(knowledgeHorizon*365),length(cb)))-knowledgeHorizon)*log(informationDecrease^2));
    cb(1) = 0;
    mexPortfolio('setParamBlomvallNdengo', curveIDs(i), 'cb', cb); % Second order derivative
  end
  
  [u, P, z] = mexPortfolio('solveBlomvallNdengo', instrID);
  
  f = u(1:nF);
  pi = u(nF+1:end);

  plot(T(1:end-1), f, T(1:end-1), pi);
  title(datestr(times(k)));
  pause(0.01);
end

rmpath(measurementPath)



