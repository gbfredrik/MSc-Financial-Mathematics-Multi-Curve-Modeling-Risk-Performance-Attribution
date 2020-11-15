function [Exp, Cov, Tot_exp, type] = calcCovar(plotPAParams)
    
resultMatrix(:,1) = plotPAParams{9};  % eps_A
    resultMatrix(:,2) = plotPAParams{10}; % eps_I 
    resultMatrix(:,3) = plotPAParams{11}; % eps_P
    resultMatrix(:,4) = plotPAParams{12}; % carry
    resultMatrix(:,5) = plotPAParams{23}; % shift_f
    resultMatrix(:,6) = plotPAParams{24}; % twist_f
    resultMatrix(:,7) = plotPAParams{25}; % butterfly_f
    resultMatrix(:,8) = plotPAParams{26}; % fourth_sixth_f
    resultMatrix(:,9) = plotPAParams{27}; % sum_second
    resultMatrix(:,10) = plotPAParams{28}; % shift_pi
    resultMatrix(:,11) = plotPAParams{29}; % twist_pi
    resultMatrix(:,12) = plotPAParams{30}; % butterfly_pi
    resultMatrix(:,13) = plotPAParams{31}; % fourth_eight_pi


    covPortfolio = cov(resultMatrix);
    explVar = diag(covPortfolio);

    for i = 1:13
        for j = 1:13
            eplVar(i,j) = covPortfolio(i,j) / sum(sum(covPortfolio));
        end
    end

    sortedVariances = [];
    sortedIndexes = [];
    tempVar = eplVar;
    for i = 1:13*13
       [maxColumns, rowIndexes] = max(abs(tempVar));
       [sortedVariances(i), columnIndexes] = max(maxColumns);
       sortedIndexes(i,1) = rowIndexes(columnIndexes);
       sortedIndexes(i,2) = columnIndexes;
       tempVar(rowIndexes(columnIndexes), columnIndexes) = 0;
    end

    indexMatrix = ["eps_A", "eps_I", "eps_P", "carry", "shift_f", "twist_f", "butterfly_f", "fourth_sixth_f" ...
        "sum_second", "shift_pi", "twist_pi" "butterfly_pi", "fourth_eight_pi"];



    cum = 0;
    for i = 1:length(sortedVariances)
        Exp(i) = sortedVariances(i)*100';
        type(i,1) = indexMatrix(sortedIndexes(i,1));
        type(i,2) = indexMatrix(sortedIndexes(i,2));    
    end
    Exp = Exp'
    for i = 1:length(sortedVariances)
        Cov(i) = eplVar(sortedIndexes(i,1), sortedIndexes(i,2));  
    end
    Cov = Cov'
    Tot_exp(1) = Cov(1)*100;
    for i = 2:length(sortedVariances)
        Tot_exp(i) = Tot_exp(i-1) + Cov(i)*100;  
    end
    Tot_exp = Tot_exp'
    
    count = 1;
    for i = 1:length(sortedVariances)
        temp1 = type(i,1);
        temp2 = type(i,2);
        for j = i+1:length(sortedVariances)
            if type(j,2) == temp1 && type(j,1) == temp2
                rmIndexes(count) = i
                count = count + 1;
            end
        end
    end
    
    for i = length(sortedVariances):-1:1
        if ismember(i, rmIndexes)
            Exp(i) = []
            Cov(i) = []
            Tot_exp(i) = []
            type(i,:) = []
        end 
    end
    
end

