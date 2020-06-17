function f = kdeMulti(xSimulated, xRealized)
    %xSimulated: matrix
    %xRealized: vector
    %f: vector
    
    m = length(xRealized);
    f = zeros(1:m, 1);
    
    for j = 1:m
       f(j) = kde(xSimulated, xRealized(j)); 
       %f(j) = kde(xSimulated(j,:), xRealized(j));
    end
end