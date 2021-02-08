function [xOpt, ll] = mlGARCH(r, dt, likelihoodFunction, x1, useMR)


optionvec=optimset('MaxFunEvals',2000,'TolX',1e-16,'TolFun',1e-10,'Algorithm','interior-point');

% x = [ nu    omega   beta  alpha]
x0  = [ 0.0 ; 0.001 ; 0.80 ; 0.1 ]; %Initial solution
lb  = [-Inf ; 0     ; 0    ; 0    ];
ub  = [ Inf ; Inf   ; 1    ; 1    ];
A   = [ 0     0       1      1    ];
b   = [1];

if (nargin == 4)
  x0 = x1;
end

% Solve min  f(x)
%       s.t. Ax <= b
%            lb <= x <= ub

[xOpt,f,exitflag,output,lambda,grad,hessian] = fmincon(@(x) negLikelihood(x,r,dt,@varGARCH,likelihoodFunction,useMR),x0,A,b,[],[],lb,ub,[],optionvec);

ll = likelihoodFunction(xOpt, r, dt, @varGARCH);

