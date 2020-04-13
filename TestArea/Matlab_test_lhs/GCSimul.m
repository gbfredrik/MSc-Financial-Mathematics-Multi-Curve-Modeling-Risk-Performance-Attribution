function u = GCSimul(Corr, n, s)
    rng(s);
    y = randn(size(Corr,1), n); % Independent Gaussian random variables
    A = chol(Corr); % Cholesky Factorization
    x = (A'*y); % Correlated Gaussian random variables
    u = normcdf(x); % Gaussian copula
    