function [l]=likelihoodGaussianCopula(u, params)
% Computes the log likelihood value for each realization of a normally distributed variable

rho = [1 params(1) params(3); params(1) 1 params(2); params(3) params(2) 1];

norm_inv = norminv(u);
%norm_inv = zeros(size(u));
% 
% for i=1:size(u,1)
%     f1 = norminv(u(i,1));
%     f2 = norminv(u(i,2));
%     f3 = norminv(u(i,3));
%     
%     norm_inv(i,:) = [f1 f2 f3];
% end

l = 0;

for i=1:size(u,1)
    l = l + log(sqrt(det(rho)))+1/2*norm_inv(i,:)*(inv(rho) - eye(3))*norm_inv(i,:)';
end

end

% l = % log likelihood value for a normally distributed variables

