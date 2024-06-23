% k(x,y)=0.5(1+tanh( pi(x-y)/2 ))
% h(x)=4exp(x)(Lcosh(L)-sinh(L))/L
%f(y)=sin(y)
%h=\int f k dy
% M f = h


savefigures=0; %1 to save figures


L=pi;
N=81;


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


rng(8543)
%%%%% Calculate/define h
h=M*f_orig;
noise_amp=1E-8; noise=noise_amp*2*(randn(N,1)-0.5); h_noise=h+noise;



%% SVD 
q=3;







[U,S,V]=svd(M,'vector');
u=U';

truncated_inv=zeros(N,N);
f_rank_k=zeros(N,N);
residuals=zeros(1,N);
f_rank_k_norm=zeros(1,N);
f_diff_norm=zeros(1,N);
Fourier_coeffs=zeros(1,N);





% Now loop over the ranks 
for k=1:N
    truncated=(1/S(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
    f=truncated_inv*h_noise;
    f_rank_k(:,k)=f;
    f_rank_k_norm(k)=sqrt(trapz(x,f.^2));
    residuals(k)=sqrt(trapz(x,(M*f-h_noise).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*h_noise);
    f_diff_norm(k)=sqrt(trapz(x,(f-f_orig).^2));
end


[geo_mean_indices,geo_mean] = DPC_geometric_mean(S,Fourier_coeffs,q);




%% plots

example_rank_1=27;
example_rank_2=42;




[~,best_rank]=min(f_diff_norm)


figure(43); clf; hold on;
plot(x,h,DisplayName='h')
plot(x,h_noise,DisplayName='noisy h')
legend







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




figure(2);clf;hold on; box on;
plot(log10(residuals),log10(f_rank_k_norm),'-k','LineWidth',2)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{f_{\kappa}}$}$-$\mbox{\boldmath$\underline{g}$}$|)$',Interpreter='latex',FontSize=18)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{f_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
plot(log10(residuals(best_rank)),log10(f_rank_k_norm(best_rank)),'xr',MarkerSize=13,LineWidth=3)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{f_{\kappa}}$}$-$\mbox{\boldmath$\underline{g}$}$_{\epsilon}|)$',Interpreter='latex',FontSize=18)



figure(3); clf; hold on;
plot(log10(S),'-kx',MarkerSize=8)
ylabel('$\log_{10}(\sigma_{i})$',Interpreter='latex',FontSize=18)
xlabel('$i$',Interpreter='latex',FontSize=18)
box on



figure(7); clf; hold on; box on
plot(log10(S),DisplayName='$\log_{10}(\sigma_{i})$',LineWidth=1.5)
%plot(log10(Fourier_coeffs),DisplayName='$\log_{10}($\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$)$',LineWidth=1.5)
plot(log10(Fourier_coeffs),DisplayName='$\log_{10}(|$\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$_{\epsilon}|)$',LineWidth=1.5)
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='$\log_{10}({\rho_i})$',LineWidth=1.5)
yline(log10(noise_amp),'--k',Interpreter='latex',FontSize=18,DisplayName='$\log_{10}(\epsilon)$')
xlabel('$i$',Interpreter='latex',FontSize=18)
legend(Interpreter="latex",FontSize=14,Location='southwest')
%plot(log10(Fourier_coeffs),DisplayName='$\log_{10}($\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$_{\epsilon})$',LineWidth=1.5)





figure(9); clf; hold on; box on;
plot(log10(f_rank_k_norm),'-k',LineWidth=2)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{f_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
xlabel('$\kappa$',Interpreter='latex',FontSize=18)






% if savefigures==1
%     saveas(1,'tanh_kernel_noisy_small_truncation_plots','eps')
%     saveas(1,'tanh_kernel_noisy_small_truncation_plots','fig')
%     %
%     saveas(2,'tanh_kernel_noisy_small_L_curve','eps')
%     saveas(2,'tanh_kernel_noisy_small_L_curve','fig')
%     %
%     saveas(7,'tanh_kernel_noisy_small_DPC_plot','eps')
%     saveas(7,'tanh_kernel_noisy_small_DPC_plot','fig')
%     %
%     saveas(9,'tanh_kernel_noisy_small_norm_rank_plot','eps')
%     saveas(9,'tanh_kernel_noisy_small_norm_rank_plot','fig')
% end


if savefigures==1
    saveas(1,'tanh_kernel_noisy_small_truncation_plots.eps','epsc')
    saveas(1,'tanh_kernel_noisy_small_truncation_plots','fig')
    %
    saveas(2,'tanh_kernel_noisy_small_L_curve.eps','epsc')
    saveas(2,'tanh_kernel_noisy_small_L_curve','fig')
    %
    saveas(7,'tanh_kernel_noisy_small_DPC_plot.eps','epsc')
    saveas(7,'tanh_kernel_noisy_small_DPC_plot','fig')
    %
    saveas(9,'tanh_kernel_noisy_small_norm_rank_plot.eps','epsc')
    saveas(9,'tanh_kernel_noisy_small_norm_rank_plot','fig')
end

