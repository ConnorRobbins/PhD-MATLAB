
eps=1E-6;
x_guess=[0,0,0];


jac=zeros(3,3);

vars=x_guess;

for i=1:3
    vars_p=vars;
    vars_p(i)=vars_p(i)+eps;
    
    jac(:,i)=(eqns(vars_p)-eqns(vars))/eps;
end

options=optimset('MaxIter',500);
[Xsolv,resnorm,res,exflag,out,lambda,jac_GN]=lsqnonlin(@eqns,x_guess,[],[],options);
jac_GN=full(jac_GN);



function F= eqns(x)
    F(1)=3*x(1)*(x(2)^2)*x(3) + x(2);
    F(2)=3*x(3)*x(3) +2*x(1);
    F(3)=x(2)*x(2)*x(2) -2*x(1)*x(3);
end