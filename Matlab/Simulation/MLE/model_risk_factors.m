
%clearvars -except fHist EZero fDemand fDomestic fForeign

%%
dt = 1;


%curves = table2array(fAll);
% E = table2array(EZero); 
 %f = table2array(fHist);

delta_curves = demandAvg(2:1500,1:730) - demandAvg(1:1499,1:730);

k = 3;

C = cov(delta_curves);

[V,D] = eigs(C, k);
[e,ind] = sort(diag(D),1, 'descend');
E = V(:,ind);

[Vall,Dall] = eigs(C, 730);
[eall,indall] = sort(diag(Dall),1, 'descend');
Eall = Vall(:,ind);

% Get historic risk factors
RiskFactors = E'*delta_curves';

R1 = RiskFactors(1,:)';
R2 = RiskFactors(2,:)';
R3 = RiskFactors(3,:)';
%R4 = RiskFactors(4,:)';
%R5 = RiskFactors(5,:)';
%R6 = RiskFactors(6,:)';
%%
v(1)=(std(R1))^2/dt;

%% Generate uniform risk factors assuming Gaussian distribution
% x = [ nu    omega   beta  alpha]
%x0 = [0,0.0001,0.92,0.07];
x0 = [0,0.0001,0.5,0.04];

% Estimate parameters
[xOpt_1, ll_1, grad1] = mlGARCH(R1,dt, @likelihoodNormal, x0);
[xOpt_2, ll_2, grad2] = mlGARCH(R2,dt, @likelihoodNormal, x0);
[xOpt_3, ll_3, grad3] = mlGARCH(R3,dt, @likelihoodNormal, x0);
%[xOpt_4, ll_4, grad4] = mlGARCH(R4,dt, @likelihoodNormal, x0);
%[xOpt_5, ll_5, grad5] = mlGARCH(R5,dt, @likelihoodNormal, x0);
%[xOpt_6, ll_6, grad6] = mlGARCH(R6,dt, @likelihoodNormal, x0);

likelihoodGaussian1 = sum(ll_1);
likelihoodGaussian2 = sum(ll_2);
likelihoodGaussian3 = sum(ll_3);
%likelihoodGaussian4 = sum(ll_4);
%likelihoodGaussian5 = sum(ll_5);
%likelihoodGaussian6 = sum(ll_6);

grad_norm1 = norm(grad1);
grad_norm2 = norm(grad2);
grad_norm3 = norm(grad3);
%grad_norm4 = norm(grad4);
%grad_norm5 = norm(grad5);
%grad_norm6 = norm(grad6);

v_1 = varGARCH(xOpt_1,R1,dt);
v_2 = varGARCH(xOpt_2,R2,dt);
v_3 = varGARCH(xOpt_3,R3,dt);
%v_4 = varGARCH(xOpt_4,R4,dt);
%v_5 = varGARCH(xOpt_5,R5,dt);
%v_6 = varGARCH(xOpt_6,R6,dt);
xi_1 = (R1-xOpt_1(1)*dt)./(sqrt(v_1(1:end-1))*sqrt(dt));
xi_2 = (R2-xOpt_2(1)*dt)./(sqrt(v_2(1:end-1))*sqrt(dt));
xi_3 = (R3-xOpt_3(1)*dt)./(sqrt(v_3(1:end-1))*sqrt(dt));
%xi_4 = (R4-xOpt_4(1)*dt)./(sqrt(v_4(1:end-1))*sqrt(dt));
%xi_5 = (R5-xOpt_5(1)*dt)./(sqrt(v_5(1:end-1))*sqrt(dt));
%xi_6 = (R6-xOpt_6(1)*dt)./(sqrt(v_6(1:end-1))*sqrt(dt));

u_1 = normcdf(xi_1);
u_2 = normcdf(xi_2);
u_3 = normcdf(xi_3);
%u_4 = normcdf(xi_4);
%u_5 = normcdf(xi_5);
%u_6 = normcdf(xi_6);

U = [u_1 u_2 u_3];% u_4 u_5 u_6];
%%
%   mu, w, b, a, nu;
x_test = [-0.05,0.004,0.93,0.04,3];
testll = likelihoodStudentst(x_test, R1, 1, @varGARCH);
sum_testll = sum(testll);
var_test = varGARCH(x_test, R1,1);

%% Generate uniform risk factors given Student t-distribution

x0 = [0,0.001,0.92,0.07,5];
%x0 = [0,0.001,0.8,0.008,2];

% Estimate parameters
[xOpt_1, ll_1] = mlGARCHStudents(R1,dt, @likelihoodStudentst, x0);
[xOpt_2, ll_2] = mlGARCHStudents(R2,dt, @likelihoodStudentst, x0);
[xOpt_3, ll_3] = mlGARCHStudents(R3,dt, @likelihoodStudentst, x0);
%[xOpt_4, ll_4] = mlGARCHStudents(R4,dt, @likelihoodStudentst, x0);
%[xOpt_5, ll_5] = mlGARCHStudents(R5,dt, @likelihoodStudentst, x0);
%[xOpt_6, ll_6] = mlGARCHStudents(R6,dt, @likelihoodStudentst, x0);

%xOptAll_T = [xOpt_1; xOpt_2; xOpt_3; xOpt_4; xOpt_5; xOpt_6];
%xOptAll_T = [xOpt_1; xOpt_2; xOpt_3];

likelihoodStudentst1 = sum(ll_1);
likelihoodStudentst2 = sum(ll_2);
likelihoodStudentst3 = sum(ll_3);
%likelihoodStudentst4 = sum(ll_4);
%likelihoodStudentst5 = sum(ll_5);
%likelihoodStudentst6 = sum(ll_6);

v_1 = varGARCH(xOpt_1,R1,dt);
v_2 = varGARCH(xOpt_2,R2,dt);
v_3 = varGARCH(xOpt_3,R3,dt);
%v_4 = varGARCH(xOpt_4,R4,dt);
%v_5 = varGARCH(xOpt_5,R5,dt);
%v_6 = varGARCH(xOpt_6,R6,dt);
xi_1 = (R1-xOpt_1(1)*dt)./(sqrt(v_1(1:end-1))*sqrt(dt));
xi_2 = (R2-xOpt_2(1)*dt)./(sqrt(v_2(1:end-1))*sqrt(dt));
xi_3 = (R3-xOpt_3(1)*dt)./(sqrt(v_3(1:end-1))*sqrt(dt));
%xi_4 = (R4-xOpt_4(1)*dt)./(sqrt(v_4(1:end-1))*sqrt(dt));
%xi_5 = (R5-xOpt_5(1)*dt)./(sqrt(v_5(1:end-1))*sqrt(dt));
%xi_6 = (R6-xOpt_6(1)*dt)./(sqrt(v_6(1:end-1))*sqrt(dt));

u_1 = tcdf(xi_1,xOpt_1(5));
u_2 = tcdf(xi_2,xOpt_2(5));
u_3 = tcdf(xi_3,xOpt_3(5));
%u_4 = tcdf(xi_4,xOpt_4(5));
%u_5 = tcdf(xi_5,xOpt_5(5));
%u_6 = tcdf(xi_6,xOpt_6(5));

%U = [u_1 u_2 u_3 u_4 u_5 u_6];
U = [u_1 u_2 u_3];
%%
norm_inv = norminv(U);
norm_inv(sum(isinf(norm_inv),2) > 0,:) = [];
u_1 = normcdf(norm_inv(:,1));
u_2 = normcdf(norm_inv(:,2));
u_3 = normcdf(norm_inv(:,3));
% u_4 = normcdf(norm_inv(:,4));
% u_5 = normcdf(norm_inv(:,5));
% u_6 = normcdf(norm_inv(:,6));
%U = [u_1 u_2 u_3 u_4 u_5 u_6];
U = [u_1 u_2 u_3];
%% Estimate copula with fmincon

[Opt, ll] = mlCopula(U,dt, @likelihoodGaussianCopula);

OptRho = [1 Opt(1) Opt(3);Opt(1) 1 Opt(2); Opt(3) Opt(2) 1];

%% Estimate Gaussian copula
[rhohat_Gaussian] = copulafit('Gaussian',U);
y_gaussian = copulapdf('Gaussian', U, rhohat_Gaussian);
gaussian = sum(log(y_gaussian));

%% Estimate Students t copula
[rhohat_T, nuhat] = copulafit('t',U);
y_T = copulapdf('t', U, rhohat_T, nuhat);
t = sum(log(y_T));


%%

u_rand = copularnd('t',rhohat_T,nuhat,3451);
%u_rand = copularnd('Gaussian',rhohat_Gaussian,3443);

figure;
scatterhist(u_1, u_2);
title('Historic')

figure;
scatterhist(u_rand(:,1), u_rand(:,2));
title('Copula')

