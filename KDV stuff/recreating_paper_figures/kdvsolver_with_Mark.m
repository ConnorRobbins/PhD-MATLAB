N=301;
L=80;
F=1;



x=linspace(-L,L,N);
del_x=2*L/(N-1);


kdv_forcing=1*exp(-(x.^2));
%kdv_forcing=zeros(1,N);


mu=F-1;
alpha=del_x^2;


%eta_init=zeros(1,N);
eta_init=2*mu*((sech(x*sqrt(3*mu/2))).^2);


eta_solved=kdvsolver_fun(N,kdv_forcing,alpha,mu,eta_init,F);

figure(1); clf; hold on;
plot(x,eta_solved,'-b')
%plot(x,eta_init,'-r')




function [eta]=kdvsolver_fun(N,sigma,alpha,mu,eta_init,F)
    eta=fsolve(@kdv,eta_init);
    
    
    function H=kdv(eta)
        H(1)=eta(1);
        H(2)=eta(end);
        for i=3:N
            H(i)=alpha*sigma(i-1) +(2)*eta(i-1) -(alpha)*(eta(i-1)^2) - eta(i-2) - eta(i);
        end
    end
end
