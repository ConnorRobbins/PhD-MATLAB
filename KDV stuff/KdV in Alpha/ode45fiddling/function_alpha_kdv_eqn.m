function [dydx]= function_alpha_kdv_eqn(x,y)
    dydx=zeros(2,1);
    F_hat=0;
    alpha=5;

    
    dydx(1)=y(2);
    dydx(2)=F_hat*y(1) - (y(1))^2 +alpha*exp(-(x^2));
    %dydx(2)=F_hat*y(1) - (y(1))^2 +alpha*(4*x^2 -2)*exp(-(x^2));
end