clc
N=301;
L=80;
F_hat=0;

alpha=100;


x=linspace(-L,L,N);
del_x=2*L/(N-1);





forcing=alpha*exp(-(x.^2));
%forcing=alpha*(4*x.^2 -2).*exp(-(x.^2));




% eta_init=zeros(1,N);
 eta_init=(15-15*x.^2).*exp(-0.5*x.^2);   %Guess for type II
%eta_init=-6*exp(-(x.^2));         %Guess for type I



eta_solved=kdvsolver_fun(N,forcing,F_hat,del_x,eta_init);



figure(1); clf; hold on;
plot(x,eta_solved,'-b')
%plot(x,eta_init,'-r')
ylabel('u')
xlabel('x')
% figure(2); clf; hold on;
% plot(x,9*(b^2)*eta_solved/2,'-b')

% figure(2); clf; hold on; 
% plot(x,1+eta_solved,'-b')
% plot(x,forcing,'-r')




function [eta]=kdvsolver_fun(N,forcing,F_hat,del_x,eta_init)

    options=optimoptions('fsolve','MaxIter',10000,'MaxFunctionEvaluations',1E6);

    eta=fsolve(@kdv,eta_init,options);
    
    
    function H=kdv(eta)
        H(1)=eta(1);
        H(N)=eta(end);
        H(2:N-1)=    ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) / (del_x^2))   + (eta(2:N-1).^2)   -   F_hat*eta(2:N-1)    - forcing(2:N-1);
    end
end
