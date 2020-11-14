function [paParams, paResult, plotPAParams] = initPAParams(numContracts, A, E, E_k, f, pi, fPrev, piPrev)

    % NPV
    paResult{1} = cell(numContracts, 1);
    % Carry
    paResult{2} = cell(numContracts, 1);
    % eps_I
    paResult{3} = cell(numContracts, 1);
    % eps_P
    paResult{4} = cell(numContracts, 1);
    % eps_A
    paResult{5} = cell(numContracts, 1);
    % D_t
    paResult{6} = cell(numContracts, 1);
    % Shift_1_f
    paResult{7} = cell(numContracts, 1);
    % Twist_1_f
    paResult{8} = cell(numContracts, 1);
    % Butterfly_1_f      
    paResult{9} = cell(numContracts, 1);    
    % 4th - 6th_1_f
    paResult{10} = cell(numContracts, 1);  
    % sum_second
    paResult{11} = cell(numContracts, 1); 
    
    % Shift_1_pi
    paResult{12} = cell(numContracts, 1);
    % Twist_1_pi
    paResult{13} = cell(numContracts, 1);
    % Butterfly_1_pi
    paResult{14} = cell(numContracts, 1);    
    % 4th - 8th_1_pi
    paResult{15} = cell(numContracts, 1);  
    
    % Total risk factor sum
    paResult{16} = cell(numContracts, 1); 
    % cash
    paResult{17} = cell(numContracts, 1); 
    
    paResult{18} = 0; 
    
    % aZero, aPi, aZero_k, aPi_k
    n = size(E.Zero,2);
    m = size(E.Zero,1);
    kZero = size(E_k.Zero,2);
    kPi = size(E_k.Pi,2);
    paParams{1} = [A*E.Zero, zeros(m, n)];
    paParams{2} = [zeros(m, n), A*E.Pi];
    paParams{3} = [A*E_k.Zero, zeros(m, kPi)];
    paParams{4} = [zeros(m, kZero), A*E_k.Pi];
    
    % XiBar, XiBarZero, XiBarPi
    XiBar = [E.Zero' * f; E.Pi' * pi]; 
    XiBarZero = [XiBar(1:n); zeros(n, 1)];
    XiBarPi = [zeros(n, 1); XiBar(n+1 : end)]; 
    
    paParams{5} = XiBar;
    paParams{6} = XiBarZero;
    paParams{7} = XiBarPi;

    % prevPrice
    paParams{8} = zeros(numContracts, 1);
    
    %EZero, EPi, EZero_k, EPi_k
    paParams{10} = E.Zero;
    paParams{11} = E.Pi;
    paParams{12} = E_k.Zero;
    paParams{13} = E_k.Pi;
    
    % f, pi, fPrev, piPrev, r, piSpot, rPrev, piSpotPrev
    paParams{14} = f;
    paParams{15} = pi;
    paParams{16} = fPrev;
    paParams{17} = piPrev;
    paParams{18} = A * f;
    paParams{19} = A * pi;
    paParams{20} = A * fPrev;
    paParams{21} = A * piPrev;    
    
    % cumulative NPV
    plotPAParams{1} = [0];
    % cumulative sumRiskFactors
    plotPAParams{2} = [0];
    % cumulative epsA, epsI, epsP
    plotPAParams{3} = [0];
    plotPAParams{4} = [0];
    plotPAParams{5} = [0];
    % cumulative carry
    plotPAParams{6} = [0];
    % totNPV
    plotPAParams{7} = [];
    % totSumRiskFactors
    plotPAParams{8} = [];
    % totEpsA
    plotPAParams{9} = [];
    % toptEpsI
    plotPAParams{10} = [];
    % totEpsP
    plotPAParams{11} = [];
    % totCarry
    plotPAParams{12} = [];
    % Dates
    plotPAParams{13} = [];
    
    % shift_1_f, twist_1_f, butterfly_1_f
    plotPAParams{14} = [0];
    plotPAParams{15} = [0];
    plotPAParams{16} = [0];
    % fourth_sixth_1_f, sum_second
    plotPAParams{17} = [0];
    plotPAParams{18} = [0];
    % shift_1_pi, twist_1_pi, butterfly_1_pi
    plotPAParams{19} = [0];
    plotPAParams{20} = [0];
    plotPAParams{21} = [0];
    % fourth_eight_1_pi
    plotPAParams{22} = [0];
    
    plotPAParams{23} = [];
    plotPAParams{24} = [];
    plotPAParams{25} = [];
    
    plotPAParams{26} = [];
    plotPAParams{27} = [];
    
    plotPAParams{28} = [];
    plotPAParams{29} = [];
    plotPAParams{30} = [];
    
    plotPAParams{31} = [];
    
    

end

