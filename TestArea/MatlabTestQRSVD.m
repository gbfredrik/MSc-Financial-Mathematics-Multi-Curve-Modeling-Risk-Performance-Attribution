% Från Freddans exempel typ
%m = normrnd(0, 1, 10950, 1500);
%[Q, R, E] = qr(m, 0);
%RR = R' * R;
%[U, S, V] = svd(RR, 0);


%%
clear
H = normrnd(0, 1, 10950, 1500);
%% QR SVD
tic
[Q1, R1] = qr(H, 0);
[U1, D1, V1] = svd(R1, 0);

EVec1 = Q1 * V1(:,1:6);
EVal1 = diag(D1).^2;
EVal1 = EVal1(1:6,1);
toc

%% Original SVD
tic
%Cv = H * H';
%[U2, D2, V2] = svd(Cv);
[U2, D2, V2] = svds(H, 6);
EVec2 = U2(:,1:6);
EVal2 = diag(D2).^2;
toc

%% QR PCA
% tic
% [Q3, R3] = qr(H, 0);
% Cv = R3' * R3;
% [coeff, score, latent] = pca(Cv);
% trans = Q3(:,1:end-1) * latent;
% 
% toc
% trans = sort(trans, 'descend');

%%
fprintf('Erro between full SVD and qr SVDs: %d', norm(EVec2 - EVec1))
errorEVec = norm(EVec2 - EVec1) % Near 0, awkward rounding error
ErrorEVal = norm(EVal2 - EVal1) % Approx. 0

%%

% A = randn(3,5);
% C = cov(A);
% v = ones(5,1);
% 
% f = C * v



