% Definitions
dim = 3;
n = 10000;
rho = 0.7;

R1 = ones(dim)*rho;
R2 = ones(dim)*1;
R = R1 - diag(diag(R1)) + diag(diag(R2));
%% Gaussian Copula
u = GCSimul(R,n,"default");

subplot(1,2,1)
scatter3(u(:,1),u(:,2),u(:,3));
title('Gaussian copula')

subplot(1,2,2)
scatter3(norminv(u(:,1)),norminv(u(:,2)),norminv(u(:,3)));
title('Gaussian copula and marginals')



%% Student큦 t Copula
v = 10;
u = TCSimul(R,v,n,112234127);

scatter3(u(:,1),u(:,2),u(:,3));

%% Test lim => inf (student큦 t), should converge
n = 10;
V = [2, 3, 4, 7, 15, 20, 30, 50, 100];
uGauss = GCSimul(R,n,s);
diff = zeros(length(V),1);

for i = 1:size(V,2)
    uT = TCSimul(R,V(i),n,s);
    sumGauss = sum(uGauss,2);
    g1 = sumGauss(1);
    g2 = sumGauss(2);
    sumT = sum(uT,2);
    t1 = sumT(1);
    t2 = sumT(2);
    diff(i,1) = abs(g1 - t1 + g2 - t2);
end
%% LHSD
%rho = 0.5;
%R1 = ones(dim)*rho;
%R2 = ones(dim)*1;
%R = R1 - diag(diag(R1)) + diag(diag(R2));
v = 10;
s = rng;
R = randcorr(dim)

V = LHSDSimul(R, n, s, 't', v);
subplot(2,2,1)
scatter3(V(1,:),V(2,:),V(3,:));
title('Student큦 t copula')

V = LHSDSimul(R, n, s, 'Gauss', v);
subplot(2,2,2)
scatter3(V(1,:),V(2,:),V(3,:));
title('Gaussian copula')

V = LHSDSimul(R, n, s, 't', v);
subplot(2,2,3)
scatter3(tinv(V(1,:),v),tinv(V(2,:),v),tinv(V(3,:),v));
title('Student큦 t copula and marginals')

%V = LHSDSimul(R, n, s, 'Gauss', v);
%subplot(2,2,4)
%scatter3(norminv(V(1,:)),norminv(V(2,:)),norminv(V(3,:)));
%title('Gaussian copula and marginals')

V = LHSDSimul(R, n, s, 't', v);
subplot(2,2,4)
scatter3(norminv(V(1,:)),norminv(V(2,:)),norminv(V(3,:)));
title('T copula and Gaussian marginals')

%%

X = lhsdesign(10,2)
scatter(X(:,1),X(:,2))

