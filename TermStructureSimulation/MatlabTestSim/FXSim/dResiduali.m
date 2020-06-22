function di = dResiduali(functionValues1, functionValues2, prices1, prices2)
    [rows, columns] = size(functionValues1);
    di = zeros(rows, columns);
    
    di = (prices1-functionValues1).^2 - (prices2-functionValues2).^2;
end