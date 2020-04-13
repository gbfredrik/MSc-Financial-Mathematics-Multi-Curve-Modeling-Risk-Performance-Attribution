function u = TCSimul(Corr, v, n, s)
    rng(s)
    z = randn(size(Corr,1), n); % Independent Gaussian random variables
    A = chol(Corr); % Cholesky factorization
    y = (A'*z); % Correlated Gaussian random variables
    
    s = chi2rnd(v,[1 n]); % Random number from the chi-square distribution
    
    x = (sqrt(v./s)'*ones(1,size(Corr,1)))'.*y; % Multivariate student-t distribution simulation
    
    u = tcdf(x,v); % Student´s t copula simulation
    