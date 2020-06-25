function f = kde(xSimulated, xRealized)
    n = length(xSimulated);
    sigma = sqrt(var(xSimulated));
    h = (4/(3*n))^(1/5)*sigma;
    
    f = (1/(sqrt(2*pi)*h*n))*sum(exp((-(xRealized-xSimulated).^2)/(2*h^2)));
end

% function f = Kernel(x_s, x_r)
%     variance = var(x_s);
%     sigma = sqrt(variance);
%     sum = 0;
%     %h = 0.1;
%     h = ((4)/(3*length(x_s)))^(1/5)*sigma;
% 
%     for i = 1:length(x_s)
%        sum = sum + exp(-(x_r-x_s(i))^2/(2*h^2));
%     end
%     
%     f = sum * (1/(sqrt(2*pi)*h*length(x_s)));
% end

