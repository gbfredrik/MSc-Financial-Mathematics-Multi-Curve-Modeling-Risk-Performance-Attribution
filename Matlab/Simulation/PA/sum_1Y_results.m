function [oneY_IRS_table, allY_first_IRS_table, abs_cont, abs_cont_mature] = sum_1Y_results(paResult, tradeDatesAll_OOS)
    % READ RESULTS
    % NPV
    npv = paResult{1}{1};
    % Carry
    carry = paResult{2}{1};
    % eps_I
    eps_I = paResult{3}{1};
    % eps_P
    eps_P = paResult{4}{1};
    % eps_A
    eps_A = paResult{5}{1};
    % D_t
    D_t = paResult{6}{1};
    % Shift_1_f
    shift_1_f = paResult{7}{1};
    % Twist_1_f
    twist_1_f = paResult{8}{1};
    % Butterfly_1_f
    butterfly_1_f = paResult{9}{1};    
    % 4th - 6th_1_f
    fourth_sixth_1_f = paResult{10}{1};  
    % sum_second
    sum_second = paResult{11}{1}; 
       
    % Shift_1_pi
    shift_1_pi = paResult{12}{1};
    % Twist_1_pi
    twist_1_pi = paResult{13}{1};
    % Butterfly_1_pi
    butterfly_1_pi = paResult{14}{1};    
    % 4th - 6th_1_pi
    fourth_eight_1_pi = paResult{15}{1};  
    
    oneY_IRS_table = [erase(string([num2str(day(datestr(tradeDatesAll_OOS(1:7)))) repelem('/',7)' num2str(month(datestr(tradeDatesAll_OOS(1:7))))])," "), ...
        shift_1_f(1:7)', twist_1_f(1:7)', butterfly_1_f(1:7)', ...
        fourth_sixth_1_f(1:7)', shift_1_pi(1:7)', twist_1_pi(1:7)', butterfly_1_pi(1:7)', fourth_eight_1_pi(1:7)', sum_second(1:7)', ...
        eps_P(1:7)', eps_I(1:7)', eps_A(1:7)', carry(1:7)', D_t(1:7)', npv(1:7)'];

    oneY_IRS_table = [oneY_IRS_table;   erase(string([num2str(day(datestr(tradeDatesAll_OOS(125:131)))) repelem('/',7)' num2str(month(datestr(tradeDatesAll_OOS(125:131))))])," "), ...
        shift_1_f(125:131)', twist_1_f(125:131)', butterfly_1_f(125:131)', ...
        fourth_sixth_1_f(125:131)', shift_1_pi(125:131)', twist_1_pi(125:131)', butterfly_1_pi(125:131)', fourth_eight_1_pi(125:131)', sum_second(125:131)', ...
        eps_P(125:131)', eps_I(125:131)', eps_A(125:131)', carry(125:131)', D_t(125:131)', npv(125:131)'];
    
    oneY_IRS_table = [oneY_IRS_table;   erase(string([num2str(day(datestr(tradeDatesAll_OOS(length(twist_1_pi)-6:length(twist_1_pi))))) repelem('/',7)' num2str(month(datestr(tradeDatesAll_OOS(length(twist_1_pi)-6:length(twist_1_pi)))))])," "), ...
        shift_1_f(end-6:end)', twist_1_f(end-6:end)', butterfly_1_f(end-6:end)', ...
        fourth_sixth_1_f(end-6:end)', shift_1_pi(end-6:end)', twist_1_pi(end-6:end)', butterfly_1_pi(end-6:end)', fourth_eight_1_pi(end-6:end)', sum_second(end-6:end)', ...
        eps_P(end-6:end)', eps_I(end-6:end)', eps_A(end-6:end)', carry(end-6:end)', D_t(end-6:end)', npv(end-6:end)'];
    
    allY_first_IRS_table = [];
    cash = 0;
    
    % READ RESULTS
    for i = 1:7
        % READ RESULTS
        % NPV
        npv = sum(paResult{1}{i});
        % Carry
        carry = sum(paResult{2}{i});
        % eps_I
        eps_I = sum(paResult{3}{i});
        % eps_P
        eps_P = sum(paResult{4}{i});
        % eps_A
        eps_A = sum(paResult{5}{i});
        % D_t
        D_t = sum(paResult{6}{i});
        % Shift_1_f
        shift_1_f = sum(paResult{7}{i});
        % Twist_1_f
        twist_1_f = sum(paResult{8}{i});
        % Butterfly_1_f
        butterfly_1_f = sum(paResult{9}{i});    
        % 4th - 6th_1_f
        fourth_sixth_1_f = sum(paResult{10}{i});  
        % sum_second
        sum_second = sum(paResult{11}{i}); 

        % Shift_1_pi
        shift_1_pi = sum(paResult{12}{i});
        % Twist_1_pi
        twist_1_pi = sum(paResult{13}{i});
        % Butterfly_1_pi
        butterfly_1_pi = sum(paResult{14}{i});    
        % 4th - 6th_1_pi
        fourth_eight_1_pi = sum(paResult{15}{i});
        
        % cash   
        cash = cash + sum(paResult{17}{i});
        
        allY_first_IRS_table = [allY_first_IRS_table; [shift_1_f, twist_1_f, butterfly_1_f, fourth_sixth_1_f, shift_1_pi, twist_1_pi, butterfly_1_pi, fourth_eight_1_pi, sum_second ...
            eps_P, eps_I, eps_A, carry, D_t, npv]];

    end

    allY_first_IRS_table = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cash]; allY_first_IRS_table];
    
    
    
    npv = 0;
    carry = 0; 
    eps_I = 0;
    eps_P = 0;
    eps_A = 0;
    D_t = 0;
    shift_1_f = 0; 
    twist_1_f = 0;
    butterfly_1_f = 0;
    fourth_sixth_1_f = 0; 
    sum_second = 0;
    shift_1_pi  = 0;
    twist_1_pi = 0;
    butterfly_1_pi = 0; 
    fourth_eight_1_pi = 0; 
    
    
    % READ RESULTS
    for i = 1:16
        % READ RESULTS
        % NPV
        npv = npv + sum(paResult{1}{i});
        % Carry
        carry = carry + sum(paResult{2}{i});
        % eps_I
        eps_I = eps_I + sum(paResult{3}{i});
        % eps_P
        eps_P = eps_P + sum(paResult{4}{i});
        % eps_A
        eps_A = eps_A + sum(paResult{5}{i});
        % D_t
        D_t = D_t + sum(paResult{6}{i});
        % Shift_1_f
        shift_1_f = shift_1_f + sum(paResult{7}{i});
        % Twist_1_f
        twist_1_f = twist_1_f + sum(paResult{8}{i});
        % Butterfly_1_f
        butterfly_1_f = butterfly_1_f + sum(paResult{9}{i});    
        % 4th - 6th_1_f
        fourth_sixth_1_f = fourth_sixth_1_f + sum(paResult{10}{i});  
        % sum_second
        sum_second = sum_second + sum(paResult{11}{i}); 

        % Shift_1_pi
        shift_1_pi = shift_1_pi + sum(paResult{12}{i});
        % Twist_1_pi
        twist_1_pi = twist_1_pi + sum(paResult{13}{i});
        % Butterfly_1_pi
        butterfly_1_pi = butterfly_1_pi + sum(paResult{14}{i});    
        % 4th - 6th_1_pi
        fourth_eight_1_pi = fourth_eight_1_pi + sum(paResult{15}{i});
        
    end
          abs_cont = [shift_1_f, twist_1_f, butterfly_1_f, fourth_sixth_1_f, shift_1_pi, twist_1_pi, butterfly_1_pi, fourth_eight_1_pi, sum_second ...
            eps_P, eps_I, eps_A, carry, D_t, npv];  
    
        npv = 0;
        carry = 0; 
        eps_I = 0;
        eps_P = 0;
        eps_A = 0;
        D_t = 0;
        shift_1_f = 0; 
        twist_1_f = 0;
        butterfly_1_f = 0;
        fourth_sixth_1_f = 0; 
        sum_second = 0;
        shift_1_pi  = 0;
        twist_1_pi = 0;
        butterfly_1_pi = 0; 
        fourth_eight_1_pi = 0;         
     % READ RESULTS
    for i = 1:7
        % READ RESULTS
        % NPV
        npv = npv + sum(paResult{1}{i});
        % Carry
        carry = carry + sum(paResult{2}{i});
        % eps_I
        eps_I = eps_I + sum(paResult{3}{i});
        % eps_P
        eps_P = eps_P + sum(paResult{4}{i});
        % eps_A
        eps_A = eps_A + sum(paResult{5}{i});
        % D_t
        D_t = D_t + sum(paResult{6}{i});
        % Shift_1_f
        shift_1_f = shift_1_f + sum(paResult{7}{i});
        % Twist_1_f
        twist_1_f = twist_1_f + sum(paResult{8}{i});
        % Butterfly_1_f
        butterfly_1_f = butterfly_1_f + sum(paResult{9}{i});    
        % 4th - 6th_1_f
        fourth_sixth_1_f = fourth_sixth_1_f + sum(paResult{10}{i});  
        % sum_second
        sum_second = sum_second + sum(paResult{11}{i}); 

        % Shift_1_pi
        shift_1_pi = shift_1_pi + sum(paResult{12}{i});
        % Twist_1_pi
        twist_1_pi = twist_1_pi + sum(paResult{13}{i});
        % Butterfly_1_pi
        butterfly_1_pi = butterfly_1_pi + sum(paResult{14}{i});    
        % 4th - 6th_1_pi
        fourth_eight_1_pi = fourth_eight_1_pi + sum(paResult{15}{i});
        
    end
          abs_cont_mature = [shift_1_f, twist_1_f, butterfly_1_f, fourth_sixth_1_f, shift_1_pi, twist_1_pi, butterfly_1_pi, fourth_eight_1_pi, sum_second ...
            eps_P, eps_I, eps_A, carry, D_t, npv];         
        
        
end

