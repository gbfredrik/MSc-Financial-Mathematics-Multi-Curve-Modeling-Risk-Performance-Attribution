function [paParams, paResult, priceErrorBP] = PA(paParams, paResult, valueParams, priceErrorBP)

% READ PARAMETERS ---------------------------------------------------------
    currDate = valueParams{1};        
    numContracts = length(valueParams{5});
    activeStatus = valueParams{12};
    fixingDay = valueParams{20};

    nextEstFloatPrev = paResult{18};


    
    % IR curves
    f = paParams{14};
    pi = paParams{15};
    fPrev = paParams{16};
    piPrev = paParams{17};
    rPrev = paParams{20};
    piSpotPrev = paParams{21};

    % aZero, aPi, aZero_k, aPi_k
    aZero = paParams{1};
    aPi = paParams{2};
    aZero_k = paParams{3};
    aPi_k = paParams{4};

    EZero = paParams{10};
    EPi = paParams{11};
    EZero_k = paParams{12};
    EPi_k = paParams{13};
    kZero = size(EZero_k,2);
    kPi = size(EPi_k,2);
    
    % XiBarPrev, XiBarZeroPrev, XiBarPiPrev
    XiBarZeroPrev = paParams{6};
    XiBarPiPrev = paParams{7};
    
    n = size(EZero,2);
    
    % Set risk factors     
    dXi = [EZero_k' * (f-fPrev); EPi_k' * (pi - piPrev)];
    XiBar = [EZero' * f; EPi' * pi]; 
    XiBarZero = [XiBar(1:n); zeros(n, 1)];
    XiBarPi = [zeros(n, 1); XiBar(n+1 : end)];   
    XiBardXiZero = XiBarZeroPrev + [dXi(1:kZero); zeros(2*n-kZero, 1)];
    XiBardXiPi = XiBarPiPrev + [zeros(n, 1); dXi(kZero+1:end); zeros(n-kPi, 1)]; 

    % NPV
    npv = paResult{1};
    % Carry
    carry = paResult{2};
    % eps_I
    eps_I = paResult{3};
    % eps_P
    eps_P = paResult{4};
    % eps_A
    eps_A = paResult{5};
    % D_t
    D_t = paResult{6};
    % Shift_1_f
    shift_1_f = paResult{7};
    % Twist_1_f
    twist_1_f = paResult{8};
    % Butterfly_1_f
    butterfly_1_f = paResult{9};    
    % 4th - 6th_1_f
    fourth_sixth_1_f = paResult{10};  
    % sum_second
    sum_second = paResult{11}; 
   
    % Shift_1_pi
    shift_1_pi = paResult{12};
    % Twist_1_pi
    twist_1_pi = paResult{13};
    % Butterfly_1_pi
    butterfly_1_pi = paResult{14};    
    % 4th - 8th_1_pi
    fourth_eight_1_pi = paResult{15};  


    
    % Total risk factor sum
    sumRiskFactors = paResult{16};
    % accrual
    cash = paResult{17}; 
 
% LOOP OVER ALL ACTIVE CONTRACTS ------------------------------------------

    for i = 1:numContracts

        if activeStatus(i) == 1
            
            floatDates = valueParams{5}{i};
            fixDates = valueParams{6}{i};
            numKnown = valueParams{8}(i);
            dtFix = valueParams{9}{i};
            dtFloat = valueParams{10}{i};
            fixing = valueParams{11}{i};
            RoP = valueParams{13}(i);
            N = valueParams{14}(i);
            y = valueParams{15}(i);
            timeFracFix = valueParams{16}{i};
            timeFracFloat = valueParams{17}{i}; 
            cashCurr = valueParams{19}{i};
            [floatDatesIndexes, fixDatesIndexes] = getIndexes(currDate, fixDates, floatDates);
            
            % Get previous price
            prevPrice = paParams{8}(i);

            % IF FIRST CONTRACT DAY
            if isempty(npv{i})
             
                [currPrice, nextEstFloat] = irsPriceRiskFactor(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiBarZero, XiBarPi);
                carry{i}(end+1) = 0;
                eps_I{i}(end+1) = 0;
                eps_P{i}(end+1) = currPrice;
                eps_A{i}(end+1) = 0;
                D_t{i}(end+1) = 0;
                shift_1_f{i}(end+1) = 0;
                twist_1_f{i}(end+1) = 0;
                butterfly_1_f{i}(end+1) = 0;    
                fourth_sixth_1_f{i}(end+1) = 0;  
                sum_second{i}(end+1) = 0; 
  
                shift_1_pi{i}(end+1) = 0;
                twist_1_pi{i}(end+1) = 0;
                butterfly_1_pi{i}(end+1) = 0;    
                fourth_eight_1_pi{i}(end+1) = 0;  
                
                priceErrorBP{i} = pricingErr(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiBarZero, XiBarPi)
                
                sumRiskFactors{i}(end+1) = 0;
                cash{i}(end+1) = cashCurr; %calcAccrual(RoP, dtFloat, dtFix, fixing, N, y, numKnown, currDate, floatLegDCC, fixLegDCC, floatDatesIndexes, fixDatesIndexes); 
                npv{i}(end+1) = currPrice + cash{i}(end);
                
            % IF LAST DAY
            elseif length(floatDates) == 1
                currPrice = 0;
                nextEstFloat = 0;
                eps_I{i}(end+1) = 0;
                eps_P{i}(end+1) = 0;
                eps_A{i}(end+1) = 0;
                D_t{i}(end+1) = 0;
                shift_1_f{i}(end+1) = 0;
                twist_1_f{i}(end+1) = 0;
                butterfly_1_f{i}(end+1) = 0;    
                fourth_sixth_1_f{i}(end+1) = 0;  
                sum_second{i}(end+1) = 0; 
           
                shift_1_pi{i}(end+1) = 0;
                twist_1_pi{i}(end+1) = 0;
                butterfly_1_pi{i}(end+1) = 0;    
                fourth_eight_1_pi{i}(end+1) = 0;  
       
            
                sumRiskFactors{i}(end+1) = 0;
                cash{i}(end+1) = cashCurr; %calcAccrual(RoP, dtFloat, dtFix, fixing, N, y, numKnown, currDate, floatLegDCC, fixLegDCC, floatDatesIndexes, fixDatesIndexes); 
                carry{i}(end+1) = 0 - prevPrice + cash{i}(end);
                npv{i}(end+1) = 0 - prevPrice + cash{i}(end);
                
                
            % DURING THE CONTRACT'S LIFETIME
            else
                
                   
                if fixingDay(i) == 1 && length(valueParams{11}{i}) > 1
                        diffBlomvallLibor = nextEstFloatPrev - valueParams{11}{i}(2) * 1000 * valueParams{10}{i}(2);
                else
                        diffBlomvallLibor = 0;
                end
                
                [currPrice, nextEstFloat]= irsPriceRiskFactor(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiBarZero, XiBarPi);
                [currPricePrevRiskFactor, ~] = irsPriceRiskFactor(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiBarZeroPrev, XiBarPiPrev);
                [currPriceDeltaRiskFactor, ~] = irsPriceRiskFactor(fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero, aPi, XiBardXiZero, XiBardXiPi);
                gPrev = grad(rPrev, piSpotPrev, fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero_k, aPi_k);
                HPrev = hes(rPrev, piSpotPrev, fixDatesIndexes, floatDatesIndexes, timeFracFix, timeFracFloat, dtFix, dtFloat, fixing, RoP, y, N, numKnown, aZero_k, aPi_k);            
                
                %------------------- Pricing Error ----------------------------------
                eps_P{i}(end+1) = 0;
                
                %------------------- TRUNC. ERR. ---------------------------------- 
                eps_I{i}(end+1) = currPrice - currPriceDeltaRiskFactor;

                %------------------- TAYL. APPROX. ERR ---------------------------- 
                eps_A{i}(end+1) = currPriceDeltaRiskFactor - currPricePrevRiskFactor - gPrev'*dXi - (1/2)*dXi'*HPrev*dXi;                

                %------------------- Shift, twist, butterfly... -------------------     
                % Shift_1_f
                shift_1_f{i}(end+1) = gPrev(1) * dXi(1);
                % Twist_1_f
                twist_1_f{i}(end+1) = gPrev(2) * dXi(2);
                % Butterfly_1_f
                butterfly_1_f{i}(end+1) = gPrev(3) * dXi(3);   
                % 4th - 6th_1_f
                fourth_sixth_1_f{i}(end+1) = gPrev(4) * dXi(4) + gPrev(5) * dXi(5) + gPrev(6) * dXi(6);   
                % sum_second
                sum_second{i}(end+1) = (1/2) * (dXi' * HPrev * dXi);
  
                % Shift_1_pi
                shift_1_pi{i}(end+1) = gPrev(7) * dXi(7);
                % Twist_1_pi
                twist_1_pi{i}(end+1) = gPrev(8) * dXi(8);
                % Butterfly_1_pi
                butterfly_1_pi{i}(end+1) = gPrev(9) * dXi(9);   
                % 4th - 8th_1_pi
                fourth_eight_1_pi{i}(end+1) = gPrev(10) * dXi(10) + gPrev(11) * dXi(11) + gPrev(12) * dXi(12) + gPrev(13) * dXi(13) + gPrev(14) * dXi(14);   
           
                sumRiskFactors{i}(end+1) = gPrev' * dXi + (1/2) * dXi' * HPrev * dXi;
                
                
                %------------------- CASH ----------------------------------------
                cash{i}(end+1) = cashCurr;
                
                %------------------- CARRY ----------------------------------------
                if RoP == 'p'
                    carry{i}(end+1) = currPricePrevRiskFactor - prevPrice + cash{i}(end) + diffBlomvallLibor;
                elseif RoP == 'r'
                    carry{i}(end+1) = currPricePrevRiskFactor - prevPrice + cash{i}(end) - diffBlomvallLibor;
                end
                    
                D_t{i}(end+1) = diffBlomvallLibor;
                
                %------------------- NPV ------------------------------------------  
                if RoP == 'p'
                    npv{i}(end+1) = currPrice - prevPrice + cash{i}(end) + diffBlomvallLibor;
                elseif RoP == 'r'
                    npv{i}(end+1) = currPrice - prevPrice + cash{i}(end) - diffBlomvallLibor;
                end

            end
                % prevPrice
                paParams{8}(i) = currPrice;    
        end
         
    end
        % WRITE TO OUTPUT ---------------------------------------------------------
        % Prev risk factors
        paParams{5} = XiBar; 
        paParams{6} = XiBarZero;
        paParams{7} = XiBarPi;   

        % NPV
        paResult{1} = npv;
        % Carry
        paResult{2} = carry;
        % eps_I
        paResult{3} = eps_I;
        % eps_P
        paResult{4} = eps_P;
        % eps_A
        paResult{5} = eps_A;
        % D_t
        paResult{6} = D_t; 
        % Shift_1_f
        paResult{7} = shift_1_f;
        % Twist_1_f
        paResult{8} = twist_1_f;
        % Butterfly_1_f
        paResult{9} = butterfly_1_f;    
        % 4th - 6th_1_f
        paResult{10} = fourth_sixth_1_f;  
        % sum_second
        paResult{11} = sum_second; 
     
        % Shift_1_pi
        paResult{12} = shift_1_pi;
        % Twist_1_pi
        paResult{13} = twist_1_pi;
        % Butterfly_1_pi
        paResult{14} = butterfly_1_pi;    
        % 4th - 8th_1_pi
        paResult{15} = fourth_eight_1_pi;  
 
        
        % Total risk factor sum
        paResult{16} = sumRiskFactors; 
        % cash    
        paResult{17} = cash;
        paResult{18} = nextEstFloat;
    
end

