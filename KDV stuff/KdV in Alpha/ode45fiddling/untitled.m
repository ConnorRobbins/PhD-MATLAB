function second_oder_ode
 
% SOLVE  d2x/dt2+5 dx/dt - 4 x = sin(10 t)
% initial conditions: x(0) = 0, x'(0)=0
 
t=0:0.001:3;   % time scale
 
initial_x    = 1;
initial_dxdt = 0;
 
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt] );
 
plot(t,x(:,1));
xlabel('t'); ylabel('x');
 
    function dxdt=rhs(t,x)
        dxdt_1 = x(2);
        dxdt_2 = -5*x(2) + 4*x(1) + sin(10*t);
 
        dxdt=[dxdt_1; dxdt_2];
    end
end