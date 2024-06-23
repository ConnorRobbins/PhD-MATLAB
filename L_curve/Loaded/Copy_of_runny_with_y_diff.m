%% Physical parameters 
clear

%s=load('fall_soln_a0.1_N721.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_Newton; Yt_orig=s.Yt_Newton; N_t=s.N; N_s=s.N; 
%s=load('trapped_wave_401points.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_Newton; Yt_orig=s.Yt_Newton; N_t=s.N; N_s=s.N;
s=load('N1241trappedc2.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_Newton; Yt_orig=s.Yt_Newton; N_t=s.N; N_s=s.N;
%s=load('twice_as_noisy2.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_measurement; Yt_orig=s.Yt; N_t=s.N_t; N_s=s.N_s;
%s=load('forward_test_sub3.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_Newton; Yt_orig=s.Yt_Newton; N_t=s.N; N_s=s.N;
%s=load('one_tenth_noise_N1051_set3.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_measurement; Yt_orig=s.Yt; N_t=s.N_t; N_s=s.N_s;

s=load('subcrittwobump.mat'); L=s.L; Froude=s.Froude; P=s.P; Ys=s.Ys_Newton; N_t=s.N; N_s=s.N; Yt_orig=reshape(s.Yt,[N_t,1]);


%noise_amp_1=0.000002; noise=noise_amp_1*2*(randn(N_s,1)-0.5)'; Ys=Ys+noise;

%[Phi_t,Phi_s,Theta_b_coeff_matrix,singularvals,singularvals_larger,Theta_b,U,S,V,D,maxrank,Theta] = svd_Step_fun(L,Froude,N_s,N_t,P,amplitude,b_width);
[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);



%%

Theta_rank_k_norm=zeros(maxrank,1);
Theta_rank_k_eqn_difference_norm=zeros(maxrank,1);
Yt_minus_Yb_rank_k_difference_norm=zeros(maxrank,1);
Fourier_coeffs=zeros(maxrank,1);



for k=1:maxrank
    Theta_rank_k=Theta_b(:,k);
    %Theta_rank_k_norm(k)=sqrt(Theta_rank_k'*Theta_rank_k);
    Theta_rank_k_norm(k)=sqrt(trapz(Phi_t,((Theta_rank_k).*(conj(Theta_rank_k))).^2));
    Theta_rank_k_eqn_difference=Theta_b_coeff_matrix*Theta_rank_k - D;
    Theta_rank_k_eqn_difference_norm(k)=sqrt(Theta_rank_k_eqn_difference'*Theta_rank_k_eqn_difference);
    [Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_rank_k,Theta);
    %Yt_minus_Yb_rank_k_difference_norm(k)=sqrt((Yt-Yt_orig)'*(Yt-Yt_orig));
    Yt_minus_Yb_rank_k_difference_norm(k)=sqrt(trapz(Phi_t,((Yt-Yt_orig).*(conj(Yt-Yt_orig))).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*D);
end


[~,Theta_mindex]=min(log10(Theta_rank_k_eqn_difference_norm))
[~,Ydiff_mindex]=min(log10(Yt_minus_Yb_rank_k_difference_norm)) 





%%

example_rank=32;
q=3;



[geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q);






[Yt_axminusb_mindex,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,Theta_mindex),Theta);
[Yt_ytminusyb_mindex,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,Ydiff_mindex),Theta);
[Yt_example,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank),Theta);

figure(1); clf; hold on;
plot(log10(Theta_rank_k_norm))
plot(example_rank,log10(Theta_rank_k_norm(example_rank)),'-*r',MarkerSize=10)
plot(Ydiff_mindex,log10(Theta_rank_k_norm(Ydiff_mindex)),'-xm',MarkerSize=10)
plot(Theta_mindex,log10(Theta_rank_k_norm(Theta_mindex)),'-xk',MarkerSize=10)
ylabel('log10(ThetaK norm)')
xlabel('rank K')



%figure(2); clf; hold on;
% plot(Theta_rank_k_eqn_difference_norm)
% plot(example_rank,Theta_rank_k_eqn_difference_norm(example_rank),'-xr',MarkerSize=10)
% ylabel('Ax-b norm')

figure(3); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm))
plot(example_rank,log10(Theta_rank_k_eqn_difference_norm(example_rank)),'-*r',MarkerSize=10)
plot(Ydiff_mindex,log10(Theta_rank_k_eqn_difference_norm(Ydiff_mindex)),'-xm',MarkerSize=10)
plot(Theta_mindex,log10(Theta_rank_k_eqn_difference_norm(Theta_mindex)),'-xk',MarkerSize=10)
ylabel('log10(Ax-b norm)')
xlabel('rank')

figure(4); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm),log10(Theta_rank_k_norm))
plot(log10(Theta_rank_k_eqn_difference_norm(example_rank)),log10(Theta_rank_k_norm(example_rank)),'-*r',MarkerSize=8)
plot(log10(Theta_rank_k_eqn_difference_norm(Ydiff_mindex)),log10(Theta_rank_k_norm(Ydiff_mindex)),'-xm',MarkerSize=10)
plot(log10(Theta_rank_k_eqn_difference_norm(Theta_mindex)),log10(Theta_rank_k_norm(Theta_mindex)),'-xk',MarkerSize=10)
xlabel('log10(A thetaK - b norm)')
ylabel('log10(thetaK norm)')




figure(5); clf; hold on;
plot(Phi_s,Ys,DisplayName='surface')
plot(Phi_t,Yt_orig,'-b',DisplayName='original topo')
plot(Phi_t,Yt_example,'--r',DisplayName='example rank')
ylabel('Yb')
legend

figure(6); clf; hold on;
plot(Phi_t,Yt_orig,'-b',DisplayName='True',LineWidth=2)
%plot(Phi_t,Yt_axminusb_mindex,'-k',DisplayName='Ax-b minrank',LineWidth=2)
plot(Phi_t,Yt_ytminusyb_mindex,'-m',DisplayName='Yt-yb minrank',LineWidth=2)
plot(Phi_t,Yt_example,'-r',DisplayName='example rank',LineWidth=2)
ylabel('Topos')
legend


figure(7); clf; hold on;
plot(log10(singularvals),DisplayName='singularvals')
plot(log10(Fourier_coeffs),DisplayName='Fourier coeffs |Ui^T D|')
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='moving geometric average')
plot(example_rank,log10(singularvals(example_rank)),'-*r',MarkerSize=10,HandleVisibility='off')
plot(Ydiff_mindex,log10(singularvals(Ydiff_mindex)),'-xm',MarkerSize=10,HandleVisibility='off')
plot(Theta_mindex,log10(singularvals(Theta_mindex)),'-xk',MarkerSize=10,HandleVisibility='off')
try 
    yline(log10(noise_amp_1))
end
xlabel('i')
legend


figure(12); clf; hold on;
plot(log10(Yt_minus_Yb_rank_k_difference_norm))
plot(example_rank,log10(Yt_minus_Yb_rank_k_difference_norm(example_rank)),'-*r',MarkerSize=10)
plot(Ydiff_mindex,log10(Yt_minus_Yb_rank_k_difference_norm(Ydiff_mindex)),'-xm',MarkerSize=10)
plot(Theta_mindex,log10(Yt_minus_Yb_rank_k_difference_norm(Theta_mindex)),'-xk',MarkerSize=10)
ylabel('log10(Yt-Yb)')
xlabel('rank')





figure(13); clf; hold on;
plot(log10(Fourier_coeffs./singularvals),'-k',DisplayName='fourier coeffs / singularvals')
plot(smoothdata(log10(Fourier_coeffs./singularvals),'rlowess',40),'-r',DisplayName='smoothed curve')
ylabel('log ( fourier coeffs / sing vals)')
legend