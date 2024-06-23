N=81;
L=20;
F=0.8;
a=-0.01;
b=1;



x=linspace(-L,L,N);
del_x=2*L/(N-1);





forcing=a*exp(-((b*x).^2));

%forcing=a*exp(-((b*(x-5)).^2)) +a*exp(-((b*(x+5)).^2));

mu=F-1;



eta_init=zeros(1,N);
%eta_init=2*mu*((sech(x*sqrt(3*mu/2))).^2);


eta_solved=kdvsolver_fun(N,forcing,mu,del_x,eta_init);


%%
figure(1); clf; hold on;
plot(x,eta_solved,'-b')
%plot(x,eta_init,'-r')

figure(2); clf; hold on;
plot(x,1+eta_solved,'-b')
plot(x,forcing,'-r')


figure(3); hold on;
plot(x,eta_solved)
legend




function [eta]=kdvsolver_fun(N,forcing,mu,del_x,eta_init)
    eta=fsolve(@kdv,eta_init);
    
    
    function H=kdv(eta)
        H(1)=eta(1);
        H(2:N-1)= mu*eta(2:N-1)  - 3*(eta(2:N-1).^2)/4   - ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) /  (6*del_x^2)) - 0.5*forcing(2:N-1);
        
        if mu<0
            H(N)=eta(2);
        else
            H(N)=eta(end);
        end
    end
end
