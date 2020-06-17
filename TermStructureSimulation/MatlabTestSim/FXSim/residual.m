function di = residual(functionValues1, functionValues2)
    [rows, columns] = size(functionValues1);
    di = zeros(rows, columns);
    
    di = log(functionValues1) - log(functionValues2);
end