N=301;
L=80;
F=1;
a=-0.06;
b=0.3;



x=linspace(-L,L,N);
del_x=2*L/(N-1);





forcing=a*exp(-((b*x).^2));

mu=F-1;



eta_init=zeros(1,N);
eta_init=0.1*(2-x.^2).*exp(-0.5*x.^2);


eta_solved=kdvsolver_fun(N,forcing,mu,del_x,eta_init);

paperalpha=-27*a/(2*b^4)

figure(1); clf; hold on;
plot(x,eta_solved,'-b')
%plot(x,eta_init,'-r')
ylabel('eta')
xlabel('x')
% figure(2); clf; hold on;
% plot(x,9*(b^2)*eta_solved/2,'-b')




function [eta]=kdvsolver_fun(N,forcing,mu,del_x,eta_init)
    eta=fsolve(@kdv,eta_init);
    
    
    function H=kdv(eta)
        H(1)=eta(1);
        H(N)=eta(end);
        H(2:N-1)= mu*eta(2:N-1)  - 3*(eta(2:N-1).^2)/4   - ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) /  (6*del_x^2)) - 0.5*forcing(2:N-1);
        
    end
end
