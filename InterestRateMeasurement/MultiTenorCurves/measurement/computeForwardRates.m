function [forw,rSpot,modelSol] = computeForwardRates(instruments,nF,T,model,qStr)

  iterationPrint = 0;
  
  %model.id 0-99 blomvall; 100-199 Interpolation methods;
  %200-299 Splines; 300-399 Least Square methods;
  % contransf: 0 = Linear in price
  %            1 = Linear in yield
  %            2 = Linear in forward rate/spot rate
  
% Model-list
%                       0: Blomvall ...
%			1:step
%			100:Discount
%			101:Spot
%			102:Raw
%			103:Log
%			200:AdamsDeventer
%			201:AdamsDeventerNew
%			202:CubicSpline
%			203:ExpoCubicSpline
%			204:CubicSplineFin
%			205:CubicSplineNatural
%			300:NelsonSiegel


%qStr (om man vill anv?ända fler run/dat/os-filer) kan l?ämnas tom, d?å anv?änds tex AdamsDeventer.run ist?ället f?ör AdamsDeventer'qStr'.run

%amplString = '/sw/modules/bin/cmod csh add cplex'; % Use on MAI
amplString = 'export LD_LIBRARY_PATH= '; % Use on IEI

if nargin<5
qStr='';
end


%------------------------------------


nInstr=length(instruments.data);

%if (isfield(model,'conTransf'))
  conTransf = model.conTransf;
%end

if (isempty(conTransf)) % Not given as input data, use default values
  conTransf = ones(nInstr,1);

  for i=1:nInstr
    if (strcmp(instruments.assetType{i},'Bill'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'Bond')) 
      conTransf(i) = 0; % 0 and 2 are possible choices
    elseif (strcmp(instruments.assetType{i},'CD'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'FRA'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'IBOR'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'IRS'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'OIS'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'TenorFRA'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'TenorIRS'))
      conTransf(i) = 1;
    elseif (strcmp(instruments.assetType{i},'TS'))
      conTransf(i) = 1;
    else
      error('Not able to deal with instrument of type ', assetType{i});
    end
  end
end

f = zeros(nF,1);
xi = ones(nF,1)/365;

modelSol.id = model.id;
%if (isfield(model,'E'))
  penaltyParam = model.E; %E
%end
%if (isfield(model,'E'))
  priceConScale = model.F; %F
%end
if (isempty(penaltyParam))
  penaltyParam = zeros(nInstr,1); %E
end
if (isempty(priceConScale))
  priceConScale = zeros(nInstr,1); %F
end

if (length(conTransf) ~= nInstr || length(penaltyParam) ~= nInstr || length(priceConScale) ~= nInstr)
  error('Incorrect size');
end

modelSol = model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOMVALL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model.id<80
  nOutput2 = 5;
  muS = 1;
  nIterationsS = 500;
  numericalPrecisionS = 1E-6;
  numericalPrecisionEqS = 1E-7;
  gammaS = 0.1; %?
  xiS = 0.99; %?

  if model.id==0 
    modelString='Blomvall';

    gamma = ones(nF-1,1)*0;
    phi = ones(nF-1,1);

    %phi=exp(genPhiGauss(1.4,100,inds)+log(genPhiLn(1,0.8,inds(end))));
    %gamma=genGammaDS(1,0.002,inds(end));
    %gamma(1:90)=1;
%     fl = zeros(size(f)); % Only positive forward rates
    fl = -ones(size(f)); % Allow negative rates (SEK 2015)

  elseif model.id==1
    modelString='StepFunction';
    [tmp,t] = determineKnots(instruments, T);

    [tmp, tind] = min(abs(repmat(T,1,length(t))-repmat(t',length(T),1)));
    tind(1) = [];
    tind(end) = [];

    gamma = ones(nF-1,1); %f?örsta
    gamma(tind-1) = 0;
    phi = ones(nF-1,1)*0; %andra
    fl = zeros(size(f)); % Only positive forward rates

  elseif model.id==2
    modelString='Dynamic';
    
    tempModel.id=1;
    tempModel.E = model.E;
    tempModel.F = model.F;
    tempModel.conTransf = model.conTransf;
    
    [fTemp,~,~]=computeForwardRates(instruments,nF,T,tempModel);
    [tmp,t] = determineKnots(instruments, T);
    [tmp, tind] = min(abs(repmat(T',1,length(t))-repmat(t',length(T),1)));
    tind(1) = 0;
    tind(end) = length(fTemp);
    modelSol.ind=tind;
  
    %[gamma2, phi2, modelSol.f2, modelSol.df, modelSol.d2f]=arctanCurve(tind,fTemp');
    [gamma2, phi2, modelSol.f2, modelSol.df, modelSol.d2f]=erfCurve(tind,fTemp');
    
    lg=length(gamma2);
    lp=length(phi2);
    gamma = ones(nF-1,1)*0;
    phi = ones(nF-1,1);
    if (isfield(model,'param'))
      lambda = model.param(1); %lambda
    else
      lambda = 1; %lambda
    end
    
    gamma(1:min(lg,nF-1))=(gamma2(1:min(lg,nF-1))*lambda)'+(gamma(1:min(lg,nF-1))*(1-lambda));
    phi(1:min(lp,nF-1))=(phi2(1:min(lp,nF-1))*lambda)'+(phi(1:min(lp,nF-1))*(1-lambda));
%gamma = zeros(size(gamma));
%phi = zeros(size(phi));
    scale = max(max(gamma),max(phi));
    gamma = gamma/scale;
    phi = phi/scale;
    fl = zeros(size(f)); % Only positive forward rates
    
  elseif model.id==3
    modelString='generic';
    
    gamma = model.gamma;
    phi = model.phi;
    fl = zeros(size(f)); % Only positive forward rates

  elseif model.id==50
    modelString='BlomvallLS';

    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);

    if (isfield(model,'param'))
      penaltyParam = model.param(1)*ones(nInstr,1); %E
    else
      penaltyParam = 10*ones(nInstr,1); %E
    end
    %penaltyParam(end-4:end) = 0.25;
    priceConScale = ones(nInstr,1); %F
    fl = zeros(size(f)); % Only positive forward rates
  elseif model.id==51
    modelString='BlomvallLSSpec';

    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
phi(1:min(2*365,length(phi))) = phi(1:min(2*365,length(phi)))/100;
    if (isfield(model,'param'))
      penaltyParam = model.param(1)*ones(nInstr,1); %E
    else
      penaltyParam = 10*ones(nInstr,1); %E
    end
    %penaltyParam(end-4:end) = 0.25;
    priceConScale = ones(nInstr,1); %F
    fl = zeros(size(f)); % Only positive forward rates
  elseif model.id==52
    modelString='BlomvallLSSpec';
    
    if (isfield(model,'param'))
      weightE = model.param(1); %E
      knowledgeHorizon = model.param(2); % After knowledge horizon second order derivative have equal weight
      informationDecrease = model.param(3); % How much information decrease in one year
    else
      weightE = 10; %E
      knowledgeHorizon = 5; % After knowledge horizon second order derivative have equal weight
      informationDecrease = 2; % How much information decrease in one year
    end
    
    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
    phi(1:min(round(knowledgeHorizon*365),length(phi))) = exp((T(1:min(round(knowledgeHorizon*365),length(phi)))-knowledgeHorizon)*log(informationDecrease^2));
%     penaltyParam = model.param(1)*penaltyParam; %E
    penaltyParam = weightE*ones(nInstr,1); %E
    %penaltyParam(end-4:end) = 0.25;
    priceConScale = ones(nInstr,1); %F
    fl = -ones(size(f)); % Allow negative forward rates
  elseif model.id==54
    modelString='BlomvallLSSpread';
    
    penaltyParam = model.param(1)*penaltyParam; %E
    gamma = ones(nF-1,1)*model.param(2); 
    phi = ones(nF-1,1)*model.param(3);     
    priceConScale = ones(nInstr,1); %F
    fl = -inf(size(f)); % Allow negative forward rates
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end
  
  modelSol.phi = phi;
  modelSol.gamma = gamma;
  
  if (0) % Old solver (V2)
    input = codeInput102(instruments, penaltyParam, priceConScale);
    [forw, modelSol.dual, modelSol.z, obj, modelSol.objSmooth, modelSol.time, nIter, solverStatus, flg] = mexForwardRateEstimationXiFin(input, f, fl, xi, gamma, phi, nF, nInstr, nOutput2, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);

    if (solverStatus(3) > 0) % Allow negative forward rates
      fl = -10.0*ones(size(f));
      [forw, modelSol.dual, modelSol.z, obj, modelSol.objSmooth, modelSol.time, nIter, solverStatus, flg] = mexForwardRateEstimationXiFin(input, f, fl, xi, gamma, phi, nF, nInstr, nOutput2, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
    end
  else % New solver (V3)
    input = codeInput104(instruments, penaltyParam, priceConScale, conTransf);
% for i=1:length(input)
%   fprintf('%.10f\n',input(i));
% end

    modelType = 0;
    u = f;
    us = fl;
    uh = ones(size(fl))*inf;
    obj = zeros(5,1);
    time = zeros(5,1);
    nIter = zeros(5,1);
    [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, priceDual, obj, time, nIter, solverStatus,g,G] = mexBlomvallIRV3(input, modelType, u, us, uh, gamma, phi, T, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);

%blomvallParamSimple
    if (solverStatus(3) > 0) % Allow negative forward rates
      us = -10.0*ones(size(f));
      u = f;
      [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, priceDual, obj, time, nIter, solverStatus,g,G] = mexBlomvallIRV3(input, modelType, u, us, uh, gamma, phi, T, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
    end
    forw = u;
    modelSol.dual = priceDual;
    modelSol.z = z;
    modelSol.time = time;
    
  end
  modelSol.f = forw;
 
  rSpot=forw2spot(forw,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Blomvall AMPL%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif model.id >= 80 && model.id<100

  if model.id==80 % Discount factor formulation
    modelString='BlomvallD';

    if (isfield(model,'param'))
      weightE = model.param(1); %E
      knowledgeHorizon = model.param(2); % After knowledge horizon second order derivative have equal weight
      informationDecrease = model.param(3); % How much information decrease in one year
    else
      weightE = 10; %E
      knowledgeHorizon = 5; % After knowledge horizon second order derivative have equal weight
      informationDecrease = 2; % How much information decrease in one year
    end
    
    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
    phi(1:min(round(knowledgeHorizon*365),length(phi))) = exp((T(1:min(round(knowledgeHorizon*365),length(phi)))-knowledgeHorizon)*log(informationDecrease^2));
    duration = durationModified(instruments);
    penaltyParam = ones(nInstr,1);
    penaltyParam = weightE*penaltyParam; %E
    priceConScale = duration; %F
    for i=1:length(instruments.data)
      if  (strcmp(instruments.assetType(i),'Bill'))
        y = instruments.data{i}.yieldUnique;
        dt = instruments.data{i}.dt;
        P = priceBillYield(dt,y);
      elseif  (strcmp(instruments.assetType(i),'Bond'))
        P = instruments.data{i}.dirtyUnique;
      elseif  (strcmp(instruments.assetType(i),'IBOR'))
        P = 1;
      elseif  (strcmp(instruments.assetType(i),'FRA'))
        P = 1;
      elseif  (strcmp(instruments.assetType(i),'IRS'))
        P = 1;
      else
        error('unknown');
      end
      priceConScale(i) = priceConScale(i)*P;
    end
  elseif model.id==81 % Discount factor formulation, no scaling of pricing errors
    modelString='BlomvallD';

    if (isfield(model,'param'))
      weightE = model.param(1); %E
      knowledgeHorizon = model.param(2); % After knowledge horizon second order derivative have equal weight
      informationDecrease = model.param(3); % How much information decrease in one year
    else
      weightE = 10; %E
      knowledgeHorizon = 5; % After knowledge horizon second order derivative have equal weight
      informationDecrease = 2; % How much information decrease in one year
    end
    
    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
    phi(1:min(round(knowledgeHorizon*365),length(phi))) = exp((T(1:min(round(knowledgeHorizon*365),length(phi)))-knowledgeHorizon)*log(informationDecrease^2));
    penaltyParam = ones(nInstr,1);
    penaltyParam = weightE*penaltyParam; %E
    priceConScale = ones(nInstr,1); %F
  elseif model.id==82 % Convex discount factor formulation
    modelString='BlomvallDconvex';

    if (isfield(model,'param'))
      weightE = model.param(1); %E
      knowledgeHorizon = model.param(2); % After knowledge horizon second order derivative have equal weight
      informationDecrease = model.param(3); % How much information decrease in one year
      if (length(model.param) <= 3)
        fConst = 0.05;
      else
        fConst = model.param(4);
      end
    else
      weightE = 10; %E
      knowledgeHorizon = 5; % After knowledge horizon second order derivative have equal weight
      informationDecrease = 2; % How much information decrease in one year
      fConst = 0.05;
    end
    
    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
    phi(1:min(round(knowledgeHorizon*365),length(phi))) = exp((T(1:min(round(knowledgeHorizon*365),length(phi)))-knowledgeHorizon)*log(informationDecrease^2));
    duration = durationModified(instruments);
    penaltyParam = ones(nInstr,1);
    penaltyParam = weightE*penaltyParam; %E
    priceConScale = duration; %F
    for i=1:length(instruments.data)
      if  (strcmp(instruments.assetType(i),'Bill'))
        y = instruments.data{i}.yieldUnique;
        dt = instruments.data{i}.dt;
        P = priceBillYield(dt,y);
      elseif  (strcmp(instruments.assetType(i),'Bond'))
        P = instruments.data{i}.dirtyUnique;
      elseif  (strcmp(instruments.assetType(i),'IBOR'))
        P = 1;
      elseif  (strcmp(instruments.assetType(i),'FRA'))
        P = 1;
      elseif  (strcmp(instruments.assetType(i),'IRS'))
        P = 1;
      else
        error('unknown');
      end
      priceConScale(i) = priceConScale(i)*P;
    end
    
    if (~isfield(model,'d_'))
      model.d_ = [1 ; exp(-cumsum(fConst*ones(nF,1).*xi))];
    end
  elseif model.id==83 % Convex discount factor formulation, no scaling of pricing errors
    modelString='BlomvallDconvex';

    if (isfield(model,'param'))
      weightE = model.param(1); %E
      knowledgeHorizon = model.param(2); % After knowledge horizon second order derivative have equal weight
      informationDecrease = model.param(3); % How much information decrease in one year
      if (length(model.param) <= 3)
        fConst = 0.05;
      else
        fConst = model.param(4);
      end
    else
      weightE = 10; %E
      knowledgeHorizon = 5; % After knowledge horizon second order derivative have equal weight
      informationDecrease = 2; % How much information decrease in one year
      fConst = 0.05;
    end
    
    gamma = ones(nF-1,1)*0;
    phi = 1*ones(nF-1,1);
    %phi = T(2:nF);
    phi(1:min(round(knowledgeHorizon*365),length(phi))) = exp((T(1:min(round(knowledgeHorizon*365),length(phi)))-knowledgeHorizon)*log(informationDecrease^2));
    penaltyParam = ones(nInstr,1);
    penaltyParam = weightE*penaltyParam; %E
    priceConScale = ones(nInstr,1); %F

    if (~isfield(model,'d_'))
      model.d_ = [1 ; exp(-cumsum(fConst*ones(nF,1).*xi))];
    end
  elseif model.id==90
    modelString='Blomvall';
    gamma = ones(nF-1,1)*0;
    phi = ones(nF-1,1);
  elseif model.id==95 % Least square in AMPL
    modelString='Blomvall';
    gamma = ones(nF-1,1)*0;
    phi = ones(nF-1,1);
    penaltyParam = ones(nInstr,1); %E
    priceConScale = ones(nInstr,1); %F
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end
  model.gamma = gamma;
  model.phi = phi;

  systemString=[' ; ampl < ' modelString qStr '.run'];
  osString=['os' modelString qStr];
  datString=[modelString qStr '.dat'];
  assetString=['assetData' qStr '.dat'];

  writeModelData(datString,model);
  writeAssetData(instruments, penaltyParam, priceConScale, conTransf, nF, T,assetString);

  system([amplString systemString]);
  fh=str2func(osString);
  fh();

  forw = f;
  modelSol.f = forw;
  modelSol.z = z;
  
  rSpot=forw2spot(forw,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INTERPOLATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif model.id>99 && model.id<150

  if model.id==100
    modelString='InterpolateDiscount';
  elseif model.id==101
    modelString='InterpolateSpot';
  elseif model.id==102
    modelString='InterpolateRaw';
  elseif model.id==103
    modelString='InterpolateLog';
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end

  muS = 1;
  nIterationsS = 100;
  numericalPrecisionS = 1E-9;
  numericalPrecisionEqS = 1E-9;
  gammaS = 0.1; %?
  xiS = 0.99; %?
  input = codeInput104(instruments, penaltyParam, priceConScale, conTransf);

  u = [];
  us = -ones(size(u))*Inf;
  uh = ones(size(u))*Inf;

  B = zeros(0,length(u));
  b = zeros(0,1);
  modelParam = [];
  [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, modelParam, priceDual, obj, time, nIter, flg] = mexForwardRateEstimationXiFinV3(input, model.id, modelParam, T, u, us, uh, B, b, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
  solverStatus = flg;
  
  modelSol.t = modelParam;
  modelSol.u = u;
  if (model.id == 100)
    modelSol.d = u;
    rSpot = spotInterpolateDiscountV3(modelSol.t,modelSol.u,T); 
  elseif (model.id == 101)
    modelSol.r = u;
    rSpot = spotInterpolateSpotV3(modelSol.t,modelSol.u,T);
  elseif (model.id == 102)
    modelSol.d = u;
    rSpot = spotInterpolateRawV3(modelSol.t,modelSol.u,T);
  elseif (model.id==103)
    modelSol.r = u;
    rSpot = spotInterpolateLogV3(modelSol.t,modelSol.u,T);
  end

  forw=spot2forw(rSpot,T);
  
  modelSol.dual = priceDual;
  modelSol.z = z;
  modelSol.time = time;
  
elseif model.id>149 && model.id<200

  if model.id==150
    modelString='InterpolateDiscount';
  elseif model.id==151
    modelString='InterpolateSpot';
  elseif model.id==152
    modelString='InterpolateRaw';
  elseif model.id==153
    modelString='InterpolateLog';
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end
 
  systemString=[' ; ampl < ' modelString qStr '.run'];
  osString=['os' modelString qStr];
  datString=[modelString qStr '.dat'];
  assetString=['assetData' qStr '.dat'];

  writeAssetData(instruments, penaltyParam, priceConScale, conTransf, nF, T,assetString);

  [n,t] = determineKnots(instruments, T);
  if (n < length(instruments.data))
    error('Maturity dates coincide, not sufficient number of instruments to determie curve');
  end
  writeKnotData(n,t,datString,model);
  system([amplString systemString]);
  fh=str2func(osString);
  fh();

  if (model.id == 150)
    rSpot = spotInterpolateDiscount(t,d,T); 
    modelSol.u = d;
    modelSol.d = d;
  elseif (model.id == 151)
    rSpot = spotInterpolateSpot(t,r,T);
    modelSol.u = r;
    modelSol.d = r;
  elseif (model.id == 152)
    rSpot = spotInterpolateRaw(t,f,T);     
    modelSol.u = f;
    modelSol.d = f;
  elseif (model.id==153)
    rSpot = spotInterpolateLog(t,r,T);
    modelSol.u = r;
    modelSol.d = r;
  end

  modelSol.z = z;
  forw=spot2forw(rSpot,T);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPLINES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif model.id>199 && model.id<250

  if model.id==200
    modelString='AdamsDeventer';
  elseif model.id==201
    modelString='McCullochQuadratic';
  elseif model.id==202
    modelString='CubicSpline'; % Do not use
  elseif model.id==203
    modelString='ExpoCubicSpline';
  elseif model.id==204
    modelString='CubicSplineFinR';
  elseif model.id==205
    modelString='CubicSplineNatR';
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end

  muS = 1;
  nIterationsS = 100;
  numericalPrecisionS = 1E-9;
  numericalPrecisionEqS = 1E-9;
  gammaS = 0.1; %?
  xiS = 0.99; %?
  input = codeInput104(instruments, penaltyParam, priceConScale, conTransf);

  u = [];
  us = -ones(size(u))*Inf;
  uh = ones(size(u))*Inf;

  B = zeros(0,length(u));
  b = zeros(0,1);
  modelParam = [];
  [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, modelParam, priceDual, obj, time, nIter, flg] = mexForwardRateEstimationXiFinV3(input, model.id, modelParam, T, u, us, uh, B, b, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
  solverStatus = flg;
  
  modelSol.t = modelParam;
  modelSol.u = u;
  if (model.id == 200)
    modelSol.a = u(5*(0:(length(u)/5-1))+1);
    modelSol.b = u(5*(0:(length(u)/5-1))+2);
    modelSol.c = u(5*(0:(length(u)/5-1))+3);
    modelSol.d = u(5*(0:(length(u)/5-1))+4);
    modelSol.e = u(5*(0:(length(u)/5-1))+5);
    rSpot=spotAdamsDeventerForw(modelSol.t,modelSol.a,modelSol.b,modelSol.c,modelSol.d,modelSol.e,T);    
  elseif (model.id == 204 || model.id == 205)
    modelSol.a = u(4*(0:(length(u)/4-1))+1);
    modelSol.b = u(4*(0:(length(u)/4-1))+2);
    modelSol.c = u(4*(0:(length(u)/4-1))+3);
    modelSol.d = u(4*(0:(length(u)/4-1))+4);
    rSpot=spotCubicSplineRV3(modelSol.t,modelSol.a,modelSol.b,modelSol.c,modelSol.d,T);    
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end

  forw=spot2forw(rSpot,T);
  
  modelSol.dual = priceDual;
  modelSol.z = z;
  modelSol.time = time;
  

elseif model.id>249 && model.id<300

  priceConScale = ones(nInstr,1); %F
  if (sum(penaltyParam) == 0 && (model.id==251 || model.id==256))
    penaltyParam = ones(nInstr,1); %E
  end

  if model.id==250
    modelString='AdamsDeventer';
  elseif model.id==251
    modelString='McCullochQuadratic';
  elseif model.id==252
    modelString='CubicSpline'; % Do not use
  elseif model.id==253
    modelString='ExpoCubicSpline';
  elseif model.id==254
    modelString='CubicSplineFinR';
  elseif model.id==255
    modelString='CubicSplineNatR';
  elseif model.id==256
    modelString='McCullochCubic';
  else
    fprintf('Can not handle model number %d',model.id);
    return
  end

  systemString=[' ; ampl < ' modelString qStr '.run'];
  osString=['os' modelString qStr];
  datString=[modelString qStr '.dat'];
  assetString=['assetData' qStr '.dat'];
  
  writeAssetData(instruments, penaltyParam, priceConScale, conTransf, nF, T,assetString);

  if (model.id==251)
    if (isfield(model,'param'))
      [n,t] = determineKnotsMcCullochQuadratic(instruments, T, model.param(1));
    else
      [n,t] = determineKnotsMcCullochQuadratic(instruments, T);
    end
  elseif (model.id==256)
    if (isfield(model,'param'))
      [n,t] = determineKnotsMcCullochCubic(instruments, T, model.param(1));
    else
      [n,t] = determineKnotsMcCullochCubic(instruments, T);
    end
  else
    [n,t] = determineKnots(instruments, T);
    if (n < length(instruments.data))
      error('Maturity dates coincide, not sufficient number of instruments to determie curve');
    end
  end
  writeKnotData(n,t,datString, model);
  system([amplString systemString]);
pause(1);
  fh=str2func(osString);
  fh();

  if (model.id == 250)
    rSpot = spotAdamsDeventerForw(t, a, b, c, d, e, T);
    modelSol.u = [t' ; a' ; b' ; c' ; d' ; e'];
    modelSol.t = t; modelSol.a = a; modelSol.b = b; modelSol.c = c; 
    modelSol.d = d; modelSol.e = e;
    modelSol.z = z';
  elseif (model.id==251)     
    rSpot = spotMcCullochQuadratic(t, a, T);
    modelSol.u = [t' ; a'];
    modelSol.t = t; modelSol.a = a; 
    modelSol.z = z';
  elseif (model.id == 252)
    rSpot = spotCubicSplineForw(t, a, b, c, d, T);
    % f = evalCubicSplineForw(t, a, b, c, d, T);
  elseif (model.id == 253)
    rSpot = spotExpoCubicSplineForw(t,Alpha, a,b,c,d,T);
    modelSol.t = t; modelSol.a = a; modelSol.b = b; modelSol.c = c; 
    modelSol.d = d; modelSol.Alpha = Alpha;
    modelSol.z = z';
rSpot(200:300)
    %f = evalExpoCubicSplineForw(t,alpha, a,b,c,d,T);
  elseif (model.id==256)     
    rSpot = spotMcCullochCubic(t, a, T);
    modelSol.u = [t' ; a'];
    modelSol.t = t; modelSol.a = a; 
    modelSol.z = z';
    %rSpot(200:300)
    %f = evalExpoCubicSplineForw(t,alpha, a,b,c,d,T);
  elseif (model.id==254 || model.id==255)
    rSpot=spotCubicSplineR(t,a,b,c,d,T);    
    modelSol.u = [t' ; a' ; b' ; c' ; d'];
    modelSol.t = t; modelSol.a = a; modelSol.b = b; modelSol.c = c; 
    modelSol.d = d;
    modelSol.z = z';
  end
  
  forw=spot2forw(rSpot,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEAST SQUARE%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif model.id>299 && model.id<400
  penaltyParam(penaltyParam == 0) = 1; %E
  %penaltyParam = ones(nInstr,1); %E
  priceConScale = ones(nInstr,1); %F

  if model.id==300
    modelString='NelsonSiegel';

    tmpModel.id = 50;
    [forw, rSpot, tmpModel] = computeForwardRates(instruments,nF,T,tmpModel);
    [model.beta0, model.beta1, model.beta2, model.tau] ...
       = initNelsonSiegel(rSpot(1:100:end)',T(1:100:end)); %'
  else
    fprintf('Can not handle model number %d',model.id); 
    return 
  end

  systemString=[' ; ampl < ' modelString qStr '.run'];
  osString=['os' modelString qStr];
  datString=[modelString qStr '.dat'];
  assetString=['assetData' qStr '.dat'];
  
  writeAwssetData(instruments, penaltyParam, priceConScale, conTransf, nF, T,assetString);
  writeModelData(datString,model);
  system([amplString systemString]);
  fh=str2func(osString);
  fh();
  
  modelSol.z = z';
 
  if (model.id == 300)
    rSpot = spotNelsonSiegel(beta0,beta1,beta2,tau,T);
    %rSpotnew = spotNelsonSiegelTest(beta0,beta1,beta2,tau,T);
    modelSol.u = [beta0 ; beta1 ; beta2 ; tau];
    modelSol.beta0 = beta0; modelSol.beta1 = beta1; modelSol.beta2 = beta2; 
    modelSol.tau = tau; 
  end
 
  forw=spot2forw(rSpot,T);

elseif model.id > 399 && model.id<500
  penaltyParam(penaltyParam == 0) = 1; %E
  %penaltyParam = ones(size(penaltyParam)); %E
  priceConScale = ones(size(priceConScale)); %F

  if model.id==400 % interestRatePrimalDualV3 C++
    modelString='NelsonSiegel';
    [forw,rSpot,modelSol] = NelsonSiegelV3(instruments, penaltyParam, priceConScale, conTransf, T);
    %[forw,rSpot,modelSol] = NelsonSiegelV4(instruments, penaltyParam, priceConScale, conTransf, T);
  elseif model.id==401 % interestRatePrimalDualV3 C++
    modelString='Svensson';
    [forw,rSpot,modelSol] = SvenssonV3(instruments, penaltyParam, priceConScale, conTransf, T);
  elseif model.id==402 % interestRatePrimalDualV3 C++
    modelString='Svensson';
    [forw,rSpot,modelSol] = SvenssonAllMinV3(instruments, penaltyParam, priceConScale, conTransf, T);
  elseif model.id==403 % interestRatePrimalDualV3 Matlab
    modelString='NelsonSiegel';

muS = 1;
nIterationsS = 100;
numericalPrecisionS = 1E-9;
numericalPrecisionEqS = 1E-9;
gammaS = 0.1; %?
xiS = 0.99; %?

u = [0.05 ; -0.01 ; 0.02 ; 1.0994];
us = -ones(size(u))*Inf;
us(4) = 0.01;
uh = ones(size(u))*Inf;
nInstr = length(instruments.data);
penaltyParam = 10*ones(nInstr,1); %E
priceConScale = ones(nInstr,1); %F

  B = zeros(0,length(u));
  b = zeros(0,1);
  %B = [0 0 0 1];
  %b = 5513;
  [u, x, ze, zb, sl, su, ss, sh, ye, yb, yB, yl, yu, ys, yh, iter, obj] = interestRatePrimalDualV3(instruments, penaltyParam, priceConScale, conTransf, u, us, uh, B, b, muS, nIterationsS, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);  
  solverStatus = 0;


  else
    fprintf('Can not handle model number %d',model.id); 
    return 
  end

elseif model.id > 499 && model.id<600
  

  if model.id==500 % interestRatePrimalDualV3 C++
    modelString='Blomvall';

    muS = 1;
    nIterationsS = 500;
    numericalPrecisionS = 1E-6;
    numericalPrecisionEqS = 1E-6;
    gammaS = 0.1; %?
    xiS = 0.99; %?

    gamma = ones(nF-1,1)*0;
    phi = ones(nF-1,1);
    u = f;
    us = -ones(size(u))*0;
    uh = ones(size(u))*Inf;
    nInstr = length(instruments.data);
    penaltyParam = zeros(nInstr,1); %E
    priceConScale = zeros(nInstr,1); %F

    modelType = 0;
    input = codeInput103(instruments, penaltyParam, priceConScale, conTransf);
    [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, priceDual, obj, time, nIter, flg, g, G] = mexBlomvallIRV3(input, modelType, u, us, uh, gamma, phi, T, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
    if (flg > 0) % Allow negative forward rates
      us = -10.0*ones(size(u));
      [u, x, z, sl, su, ss, sh, yc, yB, yl, yu, ys, yh, priceDual, obj, time, nIter, flg, g, G] = mexBlomvallIRV3(input, modelType, u, us, uh, gamma, phi, T, muS, nIterationsS, iterationPrint, numericalPrecisionS, numericalPrecisionEqS, gammaS, xiS);
    end
    modelSol.priceDual = priceDual;
    modelSol.solverStatus = flg;
    rSpot=forw2spot(u,T);
    modelSol.u = u;
    modelSol.z = z;
  else
    fprintf('Can not handle model number %d',model.id); 
    return 
  end

  forw=spot2forw(rSpot,T);

else
  fprintf('Can not handle model number %d',model.id);
  return 
end

if (~isfield(modelSol,'obj'))
  modelSol.obj = obj;
end
if (~isfield(modelSol,'nIter'))
  modelSol.nIter = nIter;
end
if (~isfield(modelSol,'solverStatus'))
  modelSol.solverStatus = solverStatus;
end