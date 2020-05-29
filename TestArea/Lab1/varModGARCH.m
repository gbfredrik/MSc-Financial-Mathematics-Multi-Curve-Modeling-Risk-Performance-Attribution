function [v]=varModGARCH(x,r,dt)

beta0  = x(2);
beta1  = x(3);
beta2  = x(4);
alpha0 = x(5); 
alpha1 = x(6);
v=zeros(length(r)+1,1);
v(1)=(std(r))^2/dt;

for i=1:length(r)
    if r(i) < 0
        v(i+1) = beta0 + beta1*v(i) + beta2/dt*(r(i)-alpha0*dt)^2;
    else
        v(i+1) = beta0 + beta1*v(i) + beta2/dt*(r(i)-alpha1*dt)^2;
    end
    
end
end

