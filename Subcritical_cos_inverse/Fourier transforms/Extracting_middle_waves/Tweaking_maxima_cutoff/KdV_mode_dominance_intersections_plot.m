k=linspace(0,3,1000);
F=0.95;
A=0.05;



mu=2*(F-1);

figure(40); clf; hold on;
plot(k,c1calc(A,mu,k),HandleVisibility='off')
yline(3*(A^2)/4,'--',HandleVisibility='off')
xlabel('k')

legendstring={};

if mu<0
    kright=sqrt((9*A/4)-3*mu);
    plot(kright,c1calc(A,mu,kright),'or')
    legendstring{1}=strcat('k=',num2str(kright))
    if 3*A/4 <=abs(mu)
        kleft=sqrt(-((9*A/4)+3*mu));
        plot(kleft,c1calc(A,mu,kleft),'om')
        legendstring{2}=strcat('k=',num2str(kleft))
    end
elseif 3*A/4>=mu
    kright=sqrt((9*A/4)-3*mu);
    plot(kright,c1calc(A,mu,kright),'or')
    legendstring{1}=strcat('k=',num2str(kright))
end

if sum(size(legendstring))~=0
    legend(legendstring,Location='northwest')
end

if mu~=0
    ylim([0,1.5*A*abs(mu)])
end




function [c1] = c1calc(A,mu,k)
    c1=abs(A*(mu+k.^2/3));
end