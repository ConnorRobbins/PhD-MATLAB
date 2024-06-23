k=linspace(0,3,1000);
F=0.95;
A=0.1;



mu=2*(F-1);

figure(1); clf; hold on;
plot(k,c1calc(A,mu,k))
yline(3*(A^2)/4,'--')

if mu<0
    kright=sqrt((9*A/4)-3*mu);
    plot(kright,c1calc(A,mu,kright),'or')
    if 3*A/4 <=abs(mu)
        kleft=sqrt(-((9*A/4)+3*mu));
        plot(kleft,c1calc(A,mu,kleft),'or')
    end
elseif 3*A/4>=mu
    kright=sqrt((9*A/4)-3*mu);
    plot(kright,c1calc(A,mu,kright),'or')
end



function [c1] = c1calc(A,mu,k)
    c1=abs(A*(mu+k.^2/3));
end