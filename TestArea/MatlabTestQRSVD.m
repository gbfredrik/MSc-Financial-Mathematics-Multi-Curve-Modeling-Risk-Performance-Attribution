% Från Freddans exempel typ
clear all
m = normrnd(0, 1, 1500, 3650);
%[Q, R, E] = qr(m, 0);
%RR = R' * R;
%[U, S, V] = svd(RR, 0);

%%
m_diff = diff(m);

m_centered = m_diff - mean(m_diff);
H = m_centered ./ sqrt(size(m_centered, 1) - 1);
C = cov(m_diff);
C1 = H' * H;

[U0, D0, V0] = svds(H,6);
[E, Lambda] = eigs(H'*H, 6);
%[E, Lambda] = eigs(C, 6);


%% QR SVD
tic
[Q1, R1] = qr(H', 0);
[U1, D1, V1] = svds(R1 * R1', 16);
%[U1, D1, V1] = svds(R1', 6);

V10 = Q1 * V1;
EVal1 = diag(D1);
%EVal1 = EVal1(1:6,1);
U10 = U1;
toc

%RecrR1R1T = U1 * D1 * V1';
%RecrC = V10 * D1.^2 * V10';
RecrC = V10 * D1 * V10';
fprintf('Relative Frobenius error between original cov and approximated cov: %d, as %.3d / %.3d.\n\n', norm(RecrC - C, 'fro') / norm(C, 'fro'), norm(RecrC - C, 'fro'), norm(C, 'fro'))

%%
% fprintf('Error between full SVD and qr SVDs: %d', norm(EVec2 - EVec1))
% errorEVec = norm(EVec2 - EVec1) % Near 0, awkward rounding error
% ErrorEVal = norm(EVal2 - EVal1) % Approx. 0

%%
% DZero = m_diff;
% CZero = cov(DZero);
% 
% [V,D] = eigs(CZero, 6);
% [e,ind] = sort(diag(D),1, 'descend');
% EZero = V(:,ind);
% plot(EZero(1:3650,1:3))
