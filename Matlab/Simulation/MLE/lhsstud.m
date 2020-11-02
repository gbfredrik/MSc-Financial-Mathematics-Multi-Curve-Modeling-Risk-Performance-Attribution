function U = lhsstud(mu, sigma, copulaDF,n)
nRiskFactors = size(mu,2);
%Z = mvnrnd(mu,sigma,n);

%xi = chi2rnd(copulaDF(1), 1,n);
%epsilon = sqrt(repmat(copulaDF(1),nRiskFactors,n)./repmat(xi,nRiskFactors,1));

%U = zeros(2000,nRiskFactors);
%for i = 1:nRiskFactors
%    U(:,i) = tcdf(Z(:,i).*epsilon(i,:)',copulaDF(1));
%end

U = copularnd('t', sigma, copulaDF, n);

% Find the ranks of each column
x = zeros(size(U),class(U));
for i=1:nRiskFactors
   x(:,i) = rank(U(:,i));
end

x = x - 0.5;
U = x / n;

% Transform each column back to the desired marginal distribution,
% maintaining the ranks (and therefore rank correlations) from the
% original random sample


U = U;

% -----------------------
function r=rank(x)

% Similar to tiedrank, but no adjustment for ties here
[sx, rowidx] = sort(x);
r(rowidx) = 1:length(x);
r = r(:);
