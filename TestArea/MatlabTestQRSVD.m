% Från Freddans exempel typ
%m = normrnd(0, 1, 10950, 1500);
%[Q, R, E] = qr(m, 0);
%RR = R' * R;
%[U, S, V] = svd(RR, 0);


%%
clear
H = normrnd(0, 1, 10950, 1500);
%%
tic
[Q1, R1] = qr(H, 0);
[U1, D1, V] = svd(R1, 0);

EVec = Q1 * V(:,1:6);
EVal = diag(D1).^2;
EVal = EVal(1:6,1);
toc

%%
tic
%Cv = H * H';
%[Uc, Dc, Vc] = svd(Cv);
[Uc, Dc, Vc] = svds(H, 6);
EVecc = Uc(:,1:6);
EValc = diag(Dc).^2;
toc

%%
errorEVec = norm(EVecc - EVec) % Near 0, awkward rounding error
ErrorEVal = norm(EValc - EVal) % Approx. 0