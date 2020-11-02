function[corrM] = calcCorr(cov, E)

    corrM(1, 1) = 1;
    corrM(1, 2) = cov(1, 2) / (sqrt(var(E(1,:))) * sqrt(var(E(2,:))));
    corrM(1, 3) = cov(1, 3) / (sqrt(var(E(1,:))) * sqrt(var(E(3,:))));
    
    
    corrM(2, 1) = cov(2, 1) / (sqrt(var(E(2,:))) * sqrt(var(E(1,:))));
    corrM(2, 2) = 1;
    corrM(2, 3) = cov(2, 3) / (sqrt(var(E(2,:))) * sqrt(var(E(3,:))));
    
    corrM(3, 1) = cov(3, 1) / (sqrt(var(E(3,:))) * sqrt(var(E(1,:))));
    corrM(3, 2) = cov(3, 2) / (sqrt(var(E(3,:))) * sqrt(var(E(2,:))));
    corrM(3, 3) = 1;
    
end