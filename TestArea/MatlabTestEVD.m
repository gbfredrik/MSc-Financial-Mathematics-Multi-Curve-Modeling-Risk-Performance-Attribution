A = randn(1500, 10950);
C = cov(A);
%C = A;
%%
[Q,R] = qr(C);

[U, V, D] = svd(R);

error1 = norm(C - Q * U * V * D')


%%
[u, v, d] = svd(C);
error2 = norm(C - u * v * d')

%%
norm(Q * U - u)