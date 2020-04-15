    function u = TCSimul(Corr, v, cases, s)
    rng(s)
    
    n = size(Corr,1);
    
    X = randn(cases, n); % Independent Gaussian random variables
    U = chol(Corr); % Cholesky factorization
    Z = X * U; % Correlated Gaussian random variables

    x = sqrt(gamrnd(v./2, 2, cases, 1) ./ v);
    Z = Z ./ x(:,ones(n,1));
    
    mean(Z, 1)
    u = tcdf(Z, v); % Return uniform dist
