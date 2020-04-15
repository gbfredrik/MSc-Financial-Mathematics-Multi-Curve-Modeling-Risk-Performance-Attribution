%test = [5 2 3 1 2 3; 10 5 4 4 5 6; 15 10 5 3 2 1; 20 3 1 4 5 6; 25 2 3 222 1 3; 30 1 2 444 5 123];
test = [4 4 2; 5 6 5; 3 2 3];
corrm = corr(test)

n = size(test,2);
corrm_test = zeros(n,n);
for i = 1:n
    for j = 1:n
        if (i == j)
            corrm_test(i,j) = 1;
        else
            corrm_test(i,j) = corr_egen(test(:,i),test(:,j));
        end
    end
end

corrm_test

U = chol(corrm)
L = U'

function rho = corr_egen(X,Y)
    X_hat = mean(X);
    Y_hat = mean(Y);
    rho = (sum((X - X_hat).*(Y - Y_hat))) / (sqrt(sum((X - X_hat).^2)) * sqrt(sum((Y - Y_hat).^2)));         
end