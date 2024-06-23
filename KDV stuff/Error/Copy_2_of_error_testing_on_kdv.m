F=1.2;
a=0.1;
b=0.3;
N=201;
L=80;
x_0=L;

S=5000;
sigma=10^(-3);




x=linspace(-L,L,N);
eta=a*exp(-(b*x).^2);
topo_no_error=(2*(b^2)*eta - 4*b^4 *(x.^2).*eta - (9/2)*(eta.^2) + 6*(F-1)*eta)/3;




topo_err_term=zeros(1,N);
topo_tracker=zeros(N,S);
for i=1:S
    err_variable=normrnd(0,sigma,1,N);
    for j=1:N
        if abs(x(j))<x_0
            topo_err_term(j)=(2*(b^2)*err_variable(j) - 4*(b^4)*(x(j)^2)*err_variable(j) - 9*eta(j)*err_variable(j) - (9*(err_variable(j)^2)/2) + 6*(F-1)*err_variable(j))/3;
        else
            topo_err_term(j)=0;
        end
    end
    topo_tracker(:,i)=topo_no_error+topo_err_term;
end
   
sample_mean=sum(topo_tracker,2)/S;
sample_std=sqrt(sum((topo_tracker-sample_mean).^2,2)/(S-1));


figure(1); clf; hold on;
plot(x,topo_no_error,'-b')
plot(x,topo_err_term,'-r')

figure(2); clf; hold on;
plot(x,topo_tracker)

figure(3); clf; hold on; 
plot(x,sample_mean)
ylabel('sample mean')

figure(4); clf; hold on;
plot(x,sample_std)
ylabel('sample std')


figure(5); clf; hold on;
plot(x,log(sample_std))
xlim([10,24])



%%

grad_calc_offset=25;
plotting_points=80;
xplot=linspace(log(x(end-1-plotting_points)),log(x(end-1)),plotting_points);

x1=log(x(end-1-grad_calc_offset));
x2=log(x(end-1));
y1=log(sample_std(end-1-grad_calc_offset));
y2=log(sample_std(end-1));
log_log_std_grad=(y2-y1)/(x2-x1);


figure(6); clf; hold on;
plot(log(x(ceil(N/2):end-1)),log(sample_std(ceil(N/2):end-1)))
plot(xplot,log_log_std_grad*(xplot-x1)+y1,'--r')