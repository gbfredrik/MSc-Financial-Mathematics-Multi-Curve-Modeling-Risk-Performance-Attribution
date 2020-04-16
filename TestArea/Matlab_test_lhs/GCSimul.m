function u = GCSimul(Corr, cases, s)
    rng(s);
    n = size(Corr,1);
    X = randn(cases, n); % Independent Gaussian random variables
    U = chol(Corr);; % Cholesky Factorization
    Z = X * U; % Correlated Gaussian random variables
    u = normcdf(Z); % Gaussian copula
    