function f = kde(xSimulated, xRealized)
    n = length(xSimulated);
    sigma = sqrt(var(xSimulated));
    h = (4/(3*n))^(1/5)*sigma;
    
    f = (1/(sqrt(2*pi)*h*n))*sum((-(xRealized-xSimulated).^2)/(2*h^2));
end

