function [input] = codeInput104(instr, penaltyParam, priceConScale, conTransf)
version = 104;
% Changes from 103:
% Changed order of assetType to: 
% Old: 4 = Bill, 3 = Bond, 1 = FRA, 2 = IRS
% New: 1 = Bill, 2 = Bond, 4 = FRA, 7 = IRS, 104 = TenorFRA, 107 = TenorIRS, 108 = TenorTS
% IRS: Changed order of cfDatesFix and dtFix
% Added TenorFRA, TenorIRS and TS

i = 1;
input(i) = version; i = i+1;
nAssets = 0;

for j=1:length(instr.data) 
  if (strcmp(instr.assetType{j},'Bill'))
    input(i)=1; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.maturityDate-instr.firstDate; i = i+1;   % t1
    input(i) = instr.data{j}.dt; i = i+1;
    price = instr.data{j}.price;        
  elseif (strcmp(instr.assetType{j},'Bond'))
%     input(i)=2; i=i+1; %AssetType
%     input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
%     input(i) = instr.data{j}.c; i = i+1;
%     input(i) = instr.data{j}.frq; i = i+1;
%     input(i) = instr.data{j}.dt; i = i+1;
%     input(i) = instr.data{j}.n; i = i+1;
%     for k = 1:instr.data{j}.n
%       input(i) = instr.data{j}.cfDates(k)-instr.firstDate; i = i+1;    % t
%     end
%     price = instr.data{j}.price;        
    input(i)=9; i=i+1; %AssetType, Switched to one which use cash-flows instead of coupon (CRSP data)
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.c; i = i+1;
    input(i) = instr.data{j}.frq; i = i+1;
    input(i) = instr.data{j}.dt; i = i+1;
    input(i) = instr.data{j}.n; i = i+1;
    for k = 1:instr.data{j}.n
      input(i) = instr.data{j}.cfDates(k)-instr.firstDate; i = i+1;    % t
    end
    for k = 1:instr.data{j}.n
      input(i) = instr.data{j}.cf(k); i = i+1;    % t
    end
    price = instr.data{j}.price;        
  elseif (strcmp(instr.assetType{j},'CD'))
    input(i)=3; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.maturityDate-instr.firstDate; i = i+1;   % t1
    input(i) = instr.data{j}.dt; i = i+1;
    price = instr.data{j}.price;
    if (price(1)~=-inf && price(2)~=inf)
       %price(1:2) = price(2:-1:1);%switch place ask/bid (when bill is used)
    end
  elseif (strcmp(instr.assetType{j},'FRA'))
    input(i)=4; i=i+1; %AssetType
    input(i) = instr.data{j}.valueDate-instr.firstDate; i = i+1;      % t1
    input(i) = instr.data{j}.maturityDate-instr.firstDate; i = i+1;   % t2
    input(i) = instr.data{j}.dt; i = i+1;
    price = instr.data{j}.price;      
  elseif (strcmp(instr.assetType{j},'IBOR'))
    input(i)=3; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.maturityDate-instr.firstDate; i = i+1;   % t1
    input(i) = instr.data{j}.dt; i = i+1;
    price = instr.data{j}.price;
  elseif (strcmp(instr.assetType{j},'IRS'))
    input(i)=7; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.nFix; i = i+1;
    for k = 1:instr.data{j}.nFix
      input(i) = instr.data{j}.cfDatesFix(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFix(k); i = i+1;
    end
    price = instr.data{j}.price;
  elseif (strcmp(instr.assetType{j},'OIS'))
    input(i)=7; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.nFix; i = i+1;
    for k = 1:instr.data{j}.nFix
      input(i) = instr.data{j}.cfDatesFix(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFix(k); i = i+1;
    end
    price = instr.data{j}.price;
  elseif (strcmp(instr.assetType{j},'TenorFRA'))
    input(i)=104; i=i+1; %AssetType
    input(i) = instr.data{j}.valueDate-instr.firstDate; i = i+1;      % t1
    input(i) = instr.data{j}.dValueDate; i = i+1; % d1
    input(i) = instr.data{j}.maturityDate-instr.firstDate; i = i+1;   % t2
    input(i) = instr.data{j}.dMaturityDate; i = i+1; % d2
    input(i) = instr.data{j}.dt; i = i+1;
    price = instr.data{j}.price;      
  elseif (strcmp(instr.assetType{j},'TenorIRS'))
    input(i)=107; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.dSettlementDate; i = i+1; % t0d
    input(i) = instr.data{j}.nFix; i = i+1;
    for k = 1:instr.data{j}.nFix
      input(i) = instr.data{j}.cfDatesFix(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFix(k); i = i+1;
      input(i) = instr.data{j}.dFix(k); i = i+1;    % td
    end
    input(i) = instr.data{j}.nFlt; i = i+1;
    for k = 1:instr.data{j}.nFlt
      input(i) = instr.data{j}.cfDatesFlt(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFlt(k); i = i+1;
      input(i) = instr.data{j}.dFlt(k); i = i+1;    % td
    end
    price = instr.data{j}.price;
  elseif (strcmp(instr.assetType{j},'TS'))
    input(i)=108; i=i+1; %AssetType
    input(i) = instr.data{j}.settlementDate-instr.firstDate; i = i+1; % t0
    input(i) = instr.data{j}.dSettlementDate; i = i+1; % t0d
    input(i) = instr.data{j}.nFlt1; i = i+1;
    for k = 1:instr.data{j}.nFlt1
      input(i) = instr.data{j}.cfDatesFlt1(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFlt1(k); i = i+1;
      input(i) = instr.data{j}.dFlt1(k); i = i+1;    % td
      input(i) = instr.data{j}.LFlt1(k); i = i+1;    % Expected Libor (if constant)
    end
    input(i) = instr.data{j}.nFlt2; i = i+1;
    for k = 1:instr.data{j}.nFlt2
      input(i) = instr.data{j}.cfDatesFlt2(k)-instr.firstDate; i = i+1;    % t
      input(i) = instr.data{j}.dtFlt2(k); i = i+1;
      input(i) = instr.data{j}.dFlt2(k); i = i+1;    % td
      input(i) = instr.data{j}.LFlt2(k); i = i+1;    % Expected Libor (if constant)
    end
    input(i) = instr.data{j}.constFlt1; i = i+1;
    input(i) = instr.data{j}.constFlt2; i = i+1;
    price = instr.data{j}.price;
  else
    error('Asset Type not handled');
  end
  
  
  % First bit (4) = unique, second bit (2) = bid, third bit (1) = ask
    
  if ((price(1) ~= -inf || price(2) ~= inf) && price(3) ~= inf)
    error('Can not have both bid/ask and last')
  elseif (price(1) == -inf && price(2) ~= inf)
    input(i) = 1; i = i+1;
    input(i) = price(2); i = i+1;
  elseif (price(1) ~= -inf && price(2) == inf)
    input(i) = 2; i = i+1;
    input(i) = price(1); i = i+1;
  elseif (price(1) ~= -inf && price(2) ~= inf)
    input(i) = 3; i = i+1;
    input(i) = price(1); i = i+1;
    input(i) = price(2); i = i+1;
  elseif (price(3) ~= inf)
    input(i) = 4; i = i+1;
    input(i) = price(3); i = i+1;
  else
    error('price type error');
  end
  nAssets = nAssets+1;
  input(i) = penaltyParam(nAssets); i = i+1;
  input(i) = priceConScale(nAssets); i = i+1;
  input(i) = conTransf(nAssets); i = i+1;
end