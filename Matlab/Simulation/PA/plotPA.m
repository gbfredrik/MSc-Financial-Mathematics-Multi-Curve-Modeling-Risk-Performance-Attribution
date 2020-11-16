function plotPAParams = plotPA(paResult, plotPAParams, currDate, activeStatus)

    % READ RESULTS
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
    % 4th - 6th_1_pi
    fourth_eight_1_pi = paResult{15};  
    
    % Total risk factor sum
    sumRiskFactors = paResult{16};

    cNPV = plotPAParams{1};
    totNPV =  plotPAParams{7};
    cSumRiskFactors = plotPAParams{2};
    totSumRiskFactors = plotPAParams{8};
    cEps_A = plotPAParams{3};
    totEps_A =  plotPAParams{9};
    cEps_I = plotPAParams{4};
    totEps_I =  plotPAParams{10};
    cEps_P = plotPAParams{5};
    totEps_P =  plotPAParams{11};
    cCarry = plotPAParams{6};
    totCarry =  plotPAParams{12};
    resultDates = plotPAParams{13};
    
    
    cShift_1_f = plotPAParams{14};
    cTwist_1_f = plotPAParams{15};
    cButterfly_1_f = plotPAParams{16};
    cFourth_sixth_1_f = plotPAParams{17};
    cSum_second = plotPAParams{18};
    
    cShift_1_pi = plotPAParams{19};
    cTwist_1_pi = plotPAParams{20};
    cButterfly_1_pi = plotPAParams{21};
    cFourth_eight_1_pi = plotPAParams{22};    
    
    totShift_1_f = plotPAParams{23};
    totTwist_1_f = plotPAParams{24};
    totButterfly_1_f = plotPAParams{25};
    totFourth_sixth_1_f = plotPAParams{26};
    totSum_second = plotPAParams{27};
    
    totShift_1_pi = plotPAParams{28};
    totTwist_1_pi = plotPAParams{29};
    totButterfly_1_pi = plotPAParams{30};
    totFourth_eight_1_pi = plotPAParams{31};
    
    
    
    
    % SUM UP RESULTS TO BE PLOTTED FOR ALL INSTRUMENTS
    numContracts = length(npv);
    totNPV(end+1) = 0;
    totEps_A(end+1) = 0; 
    totEps_I(end+1) = 0;
    totEps_P(end+1) = 0;
    totSumRiskFactors(end+1) = 0; 
    totCarry(end+1) = 0;
    
    totShift_1_f(end+1) = 0;
    totTwist_1_f(end+1) = 0;
    totButterfly_1_f(end+1) = 0;    
    totFourth_sixth_1_f(end+1) = 0;  
    totSum_second(end+1) = 0; 

    totShift_1_pi(end+1) = 0;
    totTwist_1_pi(end+1) = 0;
    totButterfly_1_pi(end+1) = 0;    
    totFourth_eight_1_pi(end+1) = 0;     
    

    for i = 1:numContracts
        if activeStatus(i) == 1
            totNPV(end) = totNPV(end) + npv{i}(end);
            totSumRiskFactors(end) = totSumRiskFactors(end) + sumRiskFactors{i}(end);
            totEps_A(end) = totEps_A(end) + eps_A{i}(end);
            totEps_I(end) = totEps_I(end) + eps_I{i}(end);
            totEps_P(end) = totEps_P(end) + eps_P{i}(end);
            totCarry(end) = totCarry(end) + carry{i}(end);
            
            totShift_1_f(end) = totShift_1_f(end) + shift_1_f{i}(end);
            totTwist_1_f(end) = totTwist_1_f(end) + twist_1_f{i}(end);
            totButterfly_1_f(end) = totButterfly_1_f(end) + butterfly_1_f{i}(end);
            totFourth_sixth_1_f(end) = totFourth_sixth_1_f(end) + fourth_sixth_1_f{i}(end);
            totSum_second(end) = totSum_second(end) + sum_second{i}(end);
            
            totShift_1_pi(end) = totShift_1_pi(end) + shift_1_pi{i}(end);
            totTwist_1_pi(end) = totTwist_1_pi(end) + twist_1_pi{i}(end);
            totButterfly_1_pi(end) = totButterfly_1_pi(end) + butterfly_1_pi{i}(end);
            totFourth_eight_1_pi(end) = totFourth_eight_1_pi(end) + fourth_eight_1_pi{i}(end);
            
        end
    end
    
    cNPV(end+1) = cNPV(end) + totNPV(end);
    cSumRiskFactors(end+1) = cSumRiskFactors(end) + totSumRiskFactors(end); 
    cEps_A(end+1) = cEps_A(end) + totEps_A(end);
    cEps_I(end+1) = cEps_I(end) + totEps_I(end);
    cEps_P(end+1) = cEps_P(end) + totEps_P(end);
    cCarry(end+1) = cCarry(end) + totCarry(end);
    
    cShift_1_f(end+1) =  cShift_1_f(end) + totShift_1_f(end);
    cTwist_1_f(end+1) = cTwist_1_f(end) + totTwist_1_f(end);
    cButterfly_1_f(end+1) = cButterfly_1_f(end) + totButterfly_1_f(end);
    cFourth_sixth_1_f(end+1) = cFourth_sixth_1_f(end) + totFourth_sixth_1_f(end);
    cSum_second(end+1) = cSum_second(end) + totSum_second(end);
    
    cShift_1_pi(end+1) = cShift_1_pi(end) + totShift_1_pi(end);
    cTwist_1_pi(end+1) = cTwist_1_pi(end) + totTwist_1_pi(end);
    cButterfly_1_pi(end+1) = cButterfly_1_pi(end) + totButterfly_1_pi(end);
    cFourth_eight_1_pi(end+1) = cFourth_eight_1_pi(end) + totFourth_eight_1_pi(end);
    
    %resultDates(end+1) = currDate;
    
    % PLOT RESULTS
    %resultLength  = length(cNPV);
    
    %plot(1:resultLength, cNPV, 1:resultLength, cSumRiskFactors, 1:resultLength, cEps_A, 1:resultLength, cEps_I, 1:resultLength, cEps_P, 1:resultLength, cCarry)
    
    %legend('NPV', 'SumRiskFactors', 'epsA', 'epsI', 'epsP', 'Carry');
    %pause(0.01)
   
    % WRITE TO RESULT DATASET
    plotPAParams{1} = cNPV;
    plotPAParams{7} = totNPV; 
    plotPAParams{2} = cSumRiskFactors;
    plotPAParams{8} = totSumRiskFactors;
    plotPAParams{3} = cEps_A;
    plotPAParams{9} = totEps_A;
    plotPAParams{4} = cEps_I;
    plotPAParams{10} = totEps_I;
    plotPAParams{5} = cEps_P;
    plotPAParams{11} = totEps_P;
    plotPAParams{6} = cCarry;
    plotPAParams{12} = totCarry;
    plotPAParams{13} = resultDates;    
    
    plotPAParams{14} = cShift_1_f;
    plotPAParams{15} = cTwist_1_f; 
    plotPAParams{16} = cButterfly_1_f;
    plotPAParams{17} = cFourth_sixth_1_f;
    plotPAParams{18} = cSum_second;
    plotPAParams{19} = cShift_1_pi;
    plotPAParams{20} = cTwist_1_pi;
    plotPAParams{21} = cButterfly_1_pi;
    plotPAParams{22} = cFourth_eight_1_pi;
    plotPAParams{23} = totShift_1_f;
    plotPAParams{24} = totTwist_1_f;
    plotPAParams{25} = totButterfly_1_f;
    plotPAParams{26} = totFourth_sixth_1_f;     
    
    plotPAParams{27} = totSum_second;
    plotPAParams{28} = totShift_1_pi;
    plotPAParams{29} = totTwist_1_pi;
    plotPAParams{30} = totButterfly_1_pi;
    plotPAParams{31} = totFourth_eight_1_pi;
   
    
end

