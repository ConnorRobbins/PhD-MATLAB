%% Physical parameters 
L=20;
Froude=1.2;
amplitude=0.1;
b_width=0.3;
topo_centre=10;


N_t=841;
N_s=N_t;




P=zeros(1,N_s);
Phi_s=linspace(-L,L,N_s);

%Ys=1+amplitude*0.5*(tanh(topo_centre+Phi_s)+tanh(topo_centre-Phi_s)).*cos(Phi_s);
Ys=1+amplitude*exp(-b_width*Phi_s.^2);


%[Phi_t,Phi_s,Theta_b_coeff_matrix,singularvals,singularvals_larger,Theta_b,U,S,V,D,maxrank,Theta] = svd_Step_fun(L,Froude,N_s,N_t,P,amplitude,b_width);
[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);

%%

Theta_rank_k_norm=zeros(N_t,1);
Theta_rank_k_eqn_difference_norm=zeros(N_t,1);
Fourier_coeffs=zeros(1,N_t);

for k=1:N_t
    Theta_rank_k=Theta_b(:,k);
    Theta_rank_k_norm(k)=sqrt(Theta_rank_k'*Theta_rank_k);
    Theta_rank_k_eqn_difference=Theta_b_coeff_matrix*Theta_rank_k - D;
    Theta_rank_k_eqn_difference_norm(k)=sqrt(Theta_rank_k_eqn_difference'*Theta_rank_k_eqn_difference);
    Fourier_coeffs(k)=abs(transpose(U(:,k))*D);
end


%[~,mindex]=min(log10(Theta_rank_k_eqn_difference_norm))

q=3;
[geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q);


%%

example_rank=180;

[Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank),Theta);

figure(1); clf; hold on;
plot(log10(Theta_rank_k_norm))
plot(example_rank,log10(Theta_rank_k_norm(example_rank)),'-xr',MarkerSize=8)
ylabel('log10(ThetaK norm)')
xlabel('rank K')

figure(2); clf; hold on;
plot(Theta_rank_k_eqn_difference_norm)
plot(example_rank,Theta_rank_k_eqn_difference_norm(example_rank),'-xr',MarkerSize=8)
ylabel('Ax-b norm')

figure(3); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm))
plot(example_rank,log10(Theta_rank_k_eqn_difference_norm(example_rank)),'-xr',MarkerSize=8)
ylabel('log10(Ax-b norm)')

figure(4); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm),log10(Theta_rank_k_norm))
plot(log10(Theta_rank_k_eqn_difference_norm(example_rank)),log10(Theta_rank_k_norm(example_rank)),'-rx',MarkerSize=8)
xlabel('log10(A thetaK - b norm)')
ylabel('log10(thetaK norm)')




figure(5); clf; hold on;
plot(Phi_t,Yt)
ylabel('Yb')








figure(7); clf; hold on; box on
plot(log10(singularvals),DisplayName='$\log_{10}(\sigma_{i})$',LineWidth=1.5)
%plot(log10(Fourier_coeffs),DisplayName='$\log_{10}($\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$)$',LineWidth=1.5)
plot(log10(Fourier_coeffs),DisplayName='$\log_{10}($\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$_{\epsilon})$',LineWidth=1.5)
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='$\log_{10}({\rho_i})$',LineWidth=1.5)
%yline(log10(noise_amp),'--k',Interpreter='latex',FontSize=18,DisplayName='$\log_{10}(\epsilon)$')
xlabel('$i$',Interpreter='latex',FontSize=18)
legend(Interpreter="latex",FontSize=14,Location='northeast')
%plot(log10(Fourier_coeffs),DisplayName='$\log_{10}($\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{g}$}$_{\epsilon})$',LineWidth=1.5)


