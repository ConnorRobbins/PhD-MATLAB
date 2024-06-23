%%% u'' +u^2 - F_hat*u = alpha*exp(-x^2)
%
% let  u'=v
% then v'= F_hat*u - v^2 - alpha*exp(-x^2)   
% 
% now let y=[u;v] so dydx=[u',v']
% 
% 
%

x=linspace(0,25,4001);
% initial_u=1;
initial_v=0.0;




%Stolen from Mark for alpha =1.5
u01    = 1.261621093750000; % u
u02    = 1.261645507812500; % d
initial_u = 0.5*(u01+u02)














[x,y] = ode45(@function_alpha_kdv_eqn,x,[initial_u; initial_v]);
figure(1); hold on;
plot(x,y(:,1))
figure(2); hold on;
plot(y(:,1),y(:,2))




function [dydx]= function_alpha_kdv_eqn(x,y)
    dydx=zeros(2,1);
    F_hat=0;
    alpha=1.5;

    
    dydx(1)=y(2);
    %dydx(2)=F_hat*y(1) - (y(1))^2 +alpha*exp(-(x^2));
    dydx(2)=F_hat*y(1) - (y(1))^2 +alpha*(4*x^2 -2)*exp(-(x^2))+alpha*exp(-2*(x^2));
end