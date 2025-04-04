F=1.2;
a=0.1;
b=0.3;
N=101;
L=25;
x_0=25;

sigma=10^(-3);




x=linspace(-L,L,N);
eta=a*exp(-(b*x).^2);
topo_no_error=(2*(b^2)*eta - 4*b^4 *(x.^2).*eta - (9/2)*(eta.^2) + 6*(F-1)*eta)/3;


err_variable=ones(1,N)*sigma;

topo_err_term=zeros(1,N);
for j=1:N
   if abs(x(j))<x_0
       topo_err_term(j)=(2*(b^2)*err_variable(j) - 4*(b^4)*(x(j)^2)*err_variable(j) - 9*eta(j)*err_variable(j) - (9*(err_variable(j)^2)/2) + 6*(F-1)*err_variable(j))/3;
   else
       topo_err_term(j)=0;
   end
end
   
figure(1); clf; hold on;
plot(x,topo_no_error,'-b')
plot(x,topo_err_term,'-r')
legend('no error topo','error terms')


figure(2); clf; hold on; 
plot(x,2*(b^2)*err_variable/3)
plot(x,- 4*(b^4)*(x.^2).*err_variable/3)
plot(x,- 9*eta.*err_variable/3)
plot(x,- (9*(err_variable.^2)/2)/3)
plot(x, 6*(F-1)*err_variable/3)
legend('2b^2 err','4b^4 x^2 err','-9 eta err','-9err^2 /2','6(F-1)err')