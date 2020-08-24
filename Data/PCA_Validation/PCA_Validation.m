table_data = readtable("../fHist730.csv");
num_data = table2array(table_data);

%%
clearvars -except num_data
D = diff(num_data); % == num_data(2:end,:)-num_data(1:end-1,:)
C = cov(D);
[U1, Sigma1, V1] = svd(D./sqrt(8));
D1 = Sigma1(1:6,1:6).^2;

[U2, Sigma2, V2] = svd(D);
Sigma2 = (Sigma2./sqrt(8)).^2;
D2 = Sigma2(1:6,1:6).^2;

[E, D3] = eigs(C, 6);

%Sigma.^2
%D