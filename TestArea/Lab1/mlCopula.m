function [rhoOpt, ll] = mlCopula(u, dt, likelihoodFunction)


optionvec=optimset('MaxFunEvals',2000,'Display','iter','TolX',1e-12,'TolFun',1e-7,'Algorithm','interior-point');

% params0 = [ rho12  rho23   rho13 ]
params0 =   [ 0.1 0.2 0.3];


lb  = [-1 -1 -1];
%lb(1,1) = 1;
%lb(2,2) = 1;
%lb(3,3) = 1;
ub  = [1 1 1];
%ub(1,1) = 1;
%ub(2,2) = 1;
%ub(3,3) = 1;


if (nargin == 4)
  params0 = rh01;
end

% Solve min  f(x)
%       s.t. Ax <= b
%            lb <= x <= ub

[rhoOpt,f,exitflag,output,lambda,grad,hessian] = fmincon(@(params) likelihoodFunction(u, params),params0,[],[],[],[],lb,ub,[],optionvec);

ll = likelihoodFunction(u, rhoOpt);

