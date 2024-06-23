N=301;
L=30;
F=0.7;

%%%% forcing of form f(x)= a*exp(-((b*(x-c)).^2)) + a*exp(-((b*(x+c)).^2))
a=0.02;
b=1;
c_init=5;

%%%% initial guess for surface
eta_init=zeros(1,N);
%eta_init=2*mu*((sech(x*sqrt(3*mu/2))).^2);


x=linspace(-L,L,N);
del_x=2*L/(N-1);
mu=F-1;
initial_unknowns=[eta_init,c_init];




[eta,c]=kdvsolver_fun(N,a,b,mu,x,del_x,initial_unknowns);
forcing=a*exp(-((b*(x-c)).^2)) +a*exp(-((b*(x+c)).^2));

figure(1); clf; hold on;
plot(x,eta,'-b')
%plot(x,eta_init,'-r')

figure(2); clf; hold on;
plot(x,1+eta,'-b')
plot(x,forcing,'-r')





function [eta_solved,c_solved]=kdvsolver_fun(N,a,b,mu,x,del_x,initial_unknowns);
    unknowns_solved=fsolve(@kdv,initial_unknowns);
    eta_solved=unknowns_solved(1:end-1);
    c_solved=unknowns_solved(end);
    
    function H=kdv(unknowns)
        eta=unknowns(1:end-1);
        c=unknowns(end);
        
        forcing=a*exp(-((b*(x-c)).^2)) +a*exp(-((b*(x+c)).^2));
        
        H(1)=eta(1);
        H(2:N-1)= mu*eta(2:N-1)  - 3*(eta(2:N-1).^2)/4   - ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) /  (6*del_x^2)) - 0.5*forcing(2:N-1);
        
        H(N)=eta(end);
        H(N+1)=eta(end-1);

    end
end
