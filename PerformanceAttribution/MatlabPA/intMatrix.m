function [A] = intMatrix(n)


for i = 1:n
    for j = 1:n
        if j <= i 
            A(i, j) = 1 / i;
        else
            A(i, j) = 0;
        end
    end
end


end