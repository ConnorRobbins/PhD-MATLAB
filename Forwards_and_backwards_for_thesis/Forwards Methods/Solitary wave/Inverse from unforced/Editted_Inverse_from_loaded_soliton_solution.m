clear


%load('solitary_wave_soln.mat')
%load('unforced_solitary_guessN81.mat')
load("unforced_solitary_F1p1_N81.mat")
N_s=N;
N_t=N;






[Phi_t,Phi_s,Theta_inverse,Theta_b_coeff_matrix,singularvals,singularvals_larger,Theta_b_matrix_inverse,U,S,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys_Newton);





u=U';

truncated_inv=zeros(N,N);
f_rank_k=zeros(N,N);
residuals=zeros(1,N);
f_rank_k_norm=zeros(1,N);
Fourier_coeffs=zeros(1,N);




% Now loop over the ranks 
for k=1:N
    truncated=(1/singularvals(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
    f=truncated_inv*D;
    f_rank_k(:,k)=f;
    f_rank_k_norm(k)=sqrt(trapz(Phi_t,f.^2));
    residuals(k)=sqrt(trapz(Phi_t,(Theta_b_coeff_matrix*f-D).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*D);
end


q=3;
[geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q);












%% 
example_rank=81;




Theta_bottom_inverse=Theta_b_matrix_inverse(:,example_rank);





[Yt_inverse,Ys_inverse,Xt_inverse,Xs_inverse] = variables_after_SVD(L,N_t,N_s,Theta_bottom_inverse,Theta_inverse);



%
figure(1); clf; hold on;
plot(Phi_t,Yt_inverse)
plot(Phi_s,Ys_Newton)



%%%%%%%%%%%%%%%% Theta_b
% figure(1); clf; hold on;
% plot(Phi_t,Theta_bottom_inverse,'-r')
% plot(Phi,Theta_bottom_Newton,'--b')
% xlabel('Phi')
% ylabel('Theta bottom')
% legend('Inverse','Newton')
% 
% figure(2); clf; hold on;
% plot(Xt_inverse,Theta_bottom_inverse,'-r')
% plot(Xt_Newton,Theta_bottom_Newton,'--b')
% xlabel('Xt')
% ylabel('Theta bottom')
% legend('Inverse','Newton')



%%%%%%%%%% Yb
figure(3); clf; hold on;
plot(Phi_t,Yt_inverse,'-r')
plot(Phi,Yt_Newton,'--b')
xlabel('Phi')
ylabel('Yb')
legend('Inverse','Newton')

figure(4); clf; hold on;
plot(Xt_inverse,Yt_inverse,'-r')
plot(Xt_Newton,Yt_Newton,'--b')
xlabel('Xt')
ylabel('Yb')
legend('Inverse','Newton')



%%%%%%%%% Ys

figure(5); clf; hold on;
plot(Xs_Newton,Ys_Newton,'-k')
plot(Phi_s,Ys_inverse,'-r')
plot(Phi,Ys_Newton,'--b')
ylabel('Yf')
legend('Xf from Newton','Phi from inverse','Phi from Newton')




figure(7); clf; hold on;
plot(log10(singularvals),DisplayName='singularvals')
plot(log10(Fourier_coeffs),DisplayName='Fourier coeffs |Ui^T b|')
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='moving geometric average')
xlabel('i')
legend

figure(9); clf; hold on;
plot(log10(f_rank_k_norm))
xlabel('rank')
ylabel('theta rank k norm')

figure(10);clf;hold on; box on;
plot(log10(residuals),log10(f_rank_k_norm),'-k','LineWidth',2)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{f_{\kappa}}$}$-$\mbox{\boldmath$\underline{g}$}$|)$',Interpreter='latex',FontSize=18)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{f_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
plot(log10(residuals(example_rank)),log10(f_rank_k_norm(example_rank)),'xr',MarkerSize=13,LineWidth=3)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{f_{\kappa}}$}$-$\mbox{\boldmath$\underline{g}$}$_{\epsilon}|)$',Interpreter='latex',FontSize=18)




