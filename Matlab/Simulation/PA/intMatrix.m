function [A] = intMatrix(n)
A = zeros(n, n);

for i = 1:n
    A(i, 1:n <= i) = 1 / i;
end

end
