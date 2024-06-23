% k(x,y)=sin(x-y)
% h(x)=(sin(L)cos(L)-L)cos(x)
%f=sin(x)
%h=\int f k dy
% M f = h


L=pi;
N=25;


x=linspace(-L,L,N)';
y=linspace(-L,L,N);


h=(sin(L)*cos(L)-L)*cos(x); 
% add noise if line below is uncommented
noise_amp_1=0.0002; noise=noise_amp_1*2*(randn(N,1)-0.5); h_1=h+noise;
noise_amp_2=2; noise=noise_amp_2*2*(randn(N,1)-0.5); h_2=h+noise; 


M=sin(x-y);


M(:,1)=M(:,1)/2;
M(:,end)=M(:,end)/2;
M=M*(x(2)-x(1));



%f_inv=M\h;



%% SVD 
q=3;







[U,S,V]=svd(M,'vector');
u=U';

truncated_inv=zeros(N,N);
f_1_rank_k=zeros(N,N);
f_2_rank_k=zeros(N,N);






% Now loop over the ranks 
for k=1:N
    truncated=(1/S(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
    f_1=truncated_inv*h_1;
    f_1_rank_k(:,k)=f_1;
    f_2=truncated_inv*h_2;
    f_2_rank_k(:,k)=f_2;
end





%% plots

noise_rank_1=2;
noise_rank_2=3;


% figure(1); clf; hold on;
% subplot(1,2,1);
% plot(log10(S),'-kx',MarkerSize=8)
% ylabel('$\log_{10}(\sigma_{i})$',Interpreter='latex',FontSize=18)
% xlabel('$i$',Interpreter='latex',FontSize=18)
% box on
% subplot(1,2,2); hold on;
% plot(x,sin(x),'--r',LineWidth=3,DisplayName='$\sin(x)$')
% plot(x,f_rank_k(:,2),'-b',LineWidth=1.5,DisplayName='$\kappa=2$')
% plot(x,f_rank_k(:,3),'-k',LineWidth=1.5,DisplayName='$\kappa=3$')
% ylim([-1.5,1.5])
% xlabel('$x$',Interpreter='latex',FontSize=18)
% ylabel('$f(x)$',Interpreter='latex',FontSize=18)
% legend(Location='southeast',Interpreter="latex",FontSize=14)


figure(1); clf; hold on;
subplot(2,1,1); hold on;
plot(x,sin(x),'--r',LineWidth=3,DisplayName='$\sin(x)$')
plot(x,f_1_rank_k(:,noise_rank_1),'-b',LineWidth=1.2,DisplayName='$\kappa=2$')
%plot(x,f_1_rank_k(:,noise_rank_2),'-k',LineWidth=1.5,DisplayName='$\kappa=3$')
ylim([-1.5,1.5])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('$f(x)$',Interpreter='latex',FontSize=18)
legend(Location='southeast',Interpreter="latex",FontSize=12)
box on
subplot(2,1,2); hold on;
plot(x,sin(x),'--r',LineWidth=3,DisplayName='$\sin(x)$')
plot(x,f_2_rank_k(:,noise_rank_1),'-b',LineWidth=1.2,DisplayName='$\kappa=2$')
%plot(x,f_2_rank_k(:,noise_rank_2),'-k',LineWidth=1.5,DisplayName='$\kappa=3$')
ylim([-1.5,1.5])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('$f(x)$',Interpreter='latex',FontSize=18)
legend(Location='southeast',Interpreter="latex",FontSize=12)
box on






% figure(3); clf; hold on;
% plot(log10(S),'-kx',MarkerSize=8)
% ylabel('$\log_{10}(\sigma_{i})$',Interpreter='latex',FontSize=18)
% xlabel('$i$',Interpreter='latex',FontSize=18)
% box on

% figure(8); clf; hold on; 
% plot(x,f_rank_k(:,2),'-b',LineWidth=1)
% plot(x,f_rank_k(:,3),'-k',LineWidth=1.5)
% plot(x,sin(x),'--r',LineWidth=3)
% ylim([-1.5,1.5])
% xlabel('$x$',Interpreter='latex',FontSize=18)