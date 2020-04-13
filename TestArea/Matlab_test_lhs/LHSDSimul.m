function V = LHSDSimul(corr, n, s, copula, v)
    
    if copula == 'Gauss'
        z = GCSimul(corr, n ,s);
    elseif copula == 't'
        z = TCSimul(corr, v, n, s)
    end
    dim = size(corr,1);
    r = zeros(size(z));
    for i = 1:dim
        r(i,:) = rank(z(i,:))
    end
    
    V = (r - 1/2)/n;
      
    
    
function r = rank(x)
    % Similar to tiedrank, but no adjustment for ties here
    [sx, rowidx] = sort(x);
    r(rowidx) = 1:length(x);
    r = r(:);