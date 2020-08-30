% Från Freddans exempel typ
%m = normrnd(0, 1, 10950, 1500);
%[Q, R, E] = qr(m, 0);
%RR = R' * R;
%[U, S, V] = svd(RR, 0);

%%
m_diff = diff(m);
m_diff2 = m(2:end,:) - m(1:end-1,:);

m_centered = m_diff - mean(m_diff);
H = m_centered ./ (size(m_centered, 1) - 1);
C = cov(m_diff);

[U0, D0, V0] = svds(H,6);


%% QR SVD
tic
[Q1, R1] = qr(H', 0);
[U1, D1, V1] = svd(R1', 0);

EVec1 = Q1 * V1(:,1:6);
EVal1 = diag(D1).^2;
EVal1 = EVal1(1:6,1);
toc

%% Original SVD
tic
%Cv = H * H';
%[U2, D2, V2] = svd(Cv);
[U2, D2, V2] = svds(H(:,1:3650)', 6);
EVec2 = U2(:,1:6);
EVal2 = diag(D2).^2;
toc

plot(EVec2(:,1:3));

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
fprintf('Error between full SVD and qr SVDs: %d', norm(EVec2 - EVec1))
errorEVec = norm(EVec2 - EVec1) % Near 0, awkward rounding error
ErrorEVal = norm(EVal2 - EVal1) % Approx. 0

%%
DZero = m_diff;
CZero = cov(DZero);

[V,D] = eigs(CZero, 6);
[e,ind] = sort(diag(D),1, 'descend');
EZero = V(:,ind);
plot(EZero(1:3650,1:3))
