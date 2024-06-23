%%% u'' +u^2 - F_hat*u = alpha*exp(-x^2)
%
% let  u'=v
% then v'= F_hat*u - v^2 - alpha*exp(-x^2)   
% 
% now let y=[u;v] so dydx=[u',v']
% 
% 
%



initial_u_guess=2;

U0=fzero(@end_point_finder,initial_u_guess);




function [Yend] = end_point_finder(u_guess)
u_guess
initial_v=0.0;
x_main=[0,3];
options=odeset('RelTol',1e-7,'AbsTol',1e-5);
[x,y] = ode23s(@function_alpha_kdv_eqn,x_main,[u_guess; initial_v],options);
figure(1); hold on;
plot(x,y(:,1))
figure(2); hold on;
plot(y(:,1),y(:,2))

Yend=y(end);
end

