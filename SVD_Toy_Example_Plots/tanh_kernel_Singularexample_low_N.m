% k(x,y)=0.5(1+tanh( pi(x-y)/2 ))
% h(x)=4exp(x)(Lcosh(L)-sinh(L))/L
%f(y)=sin(y)
%h=\int f k dy
% M f = h


L=pi;
N=25;


x=linspace(-L,L,N)';
y=linspace(-L,L,N);


f_orig=sin(x);

%%%%%% Define M
%M=sin(x-y);
M=(1+tanh(pi*(x-y)/2))/2;

%scale entries as they should be from trapz 
M(:,1)=M(:,1)/2;
M(:,end)=M(:,end)/2;
M=M*(x(2)-x(1));



%%%%% Calculate/define h
h=M*f_orig;




%% SVD 
q=3;







[U,S,V]=svd(M,'vector');
u=U';

truncated_inv=zeros(N,N);
f_rank_k=zeros(N,N);
residuals=zeros(1,N);
f_rank_k_norm=zeros(1,N);
Fourier_coeffs=zeros(1,N);





% Now loop over the ranks 
for k=1:N
    truncated=(1/S(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
    f=truncated_inv*h;
    f_rank_k(:,k)=f;
    f_rank_k_norm(k)=sqrt(trapz(x,f.^2));
    residuals(k)=sqrt(trapz(x,(M*f-h).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*h);
end


[geo_mean_indices,geo_mean] = DPC_geometric_mean(S,Fourier_coeffs,q);




%% plots

example_rank_1=7;
example_rank_2=25;








figure(2);clf;hold on;
plot(log10(residuals),log10(f_rank_k_norm))
xlabel('|Ax-b|')
ylabel('|x|')





figure(7); clf; hold on;
plot(log10(S),DisplayName='singularvals')
plot(log10(Fourier_coeffs),DisplayName='Fourier coeffs |Ui^T b|')
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='moving geometric average')
xlabel('i')
legend



% figure(8); clf; hold on; 
% plot(x,f_orig,LineWidth=2)
% plot(x,f_rank_k(:,example_rank_1),'--')
% %plot(x,L-(x.^2/L),LineWidth=2)


figure(9); clf; hold on;
plot(log10(f_rank_k_norm))




figure(1); clf; hold on;
name1="$\kappa="+num2str(example_rank_1)+"$";
name2="$\kappa="+num2str(example_rank_2)+"$";
plot(x,f_orig,'--r',LineWidth=3,DisplayName='$\sin(x)$')
plot(x,f_rank_k(:,example_rank_1),'-b',LineWidth=1.2,DisplayName=name1)
plot(x,f_rank_k(:,example_rank_2),'-k',LineWidth=1.5,DisplayName=name2)
ylim([-1.5,1.5])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('$f(x)$',Interpreter='latex',FontSize=18)
legend(Location='southeast',Interpreter="latex",FontSize=14)
box on





figure(3); clf; hold on;
plot(log10(S),'-kx',MarkerSize=8)
ylabel('$\log_{10}(\sigma_{i})$',Interpreter='latex',FontSize=18)
xlabel('$i$',Interpreter='latex',FontSize=18)
box on