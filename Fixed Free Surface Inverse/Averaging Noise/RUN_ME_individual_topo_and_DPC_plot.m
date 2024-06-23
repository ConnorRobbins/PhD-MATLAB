%% Physical parameters 
clear
saveplots=0;

%loadname="supercritical_gaussian_bump_a2E-1_b1_F1p5"; s=load(loadname,'N','L','b_width','Froude','amplitude','Ys_Newton','Yt');
loadname="hydraulicfall_soln_a4E-1_N721"; s=load(loadname,'N','L','b_width','Froude','amplitude','Ys_Newton','Yt_Newton'); s.Yt=s.Yt_Newton;
%loadname="forced_solitary_a2E-2_F1p1_N641"; s=load(loadname,'N','L','b_width','Froude','amplitude','Ys_Newton','Yt_Newton'); s.Yt=s.Yt_Newton;
%loadname="N641trappedc4"; s=load(loadname,'N','L','b_width','Froude','amplitude','Ys_Newton','Yt_Newton'); s.Yt=s.Yt_Newton;
%loadname="flat_free_Stream"; s=load(loadname,'N','L','b_width','Froude','amplitude','Ys_Newton','Yt_Newton'); s.Yt=s.Yt_Newton;



[N,N_s,N_t]=deal(s.N);
L=s.L;
b_width=s.b_width; 
Froude=s.Froude;
amplitude=s.amplitude; 
Ys=s.Ys_Newton;
Yt_true=reshape(s.Yt,N,1);


Phi_s=linspace(-L,L,N_s);
Phi_t=linspace(-L,L,N_t);

P=zeros(1,N_s);

Ys_unperturbed=Ys;
noise_amp=1E-3; noise=noise_amp*2*(randn(N_s,1))'; Ys=Ys+noise;


[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);



%%

Theta_rank_k_norm=zeros(maxrank,1);
residuals=zeros(maxrank,1);
Yt_minus_Yb_rank_k_difference_norm=zeros(maxrank,1);
Fourier_coeffs=zeros(maxrank,1);



for k=1:maxrank
    Theta_rank_k=Theta_b(:,k);
    %Theta_rank_k_norm(k)=sqrt(Theta_rank_k'*Theta_rank_k);
    Theta_rank_k_norm(k)=sqrt(trapz(Phi_t,((Theta_rank_k).*(conj(Theta_rank_k))).^2));
    Theta_rank_k_eqn_difference=Theta_b_coeff_matrix*Theta_rank_k - D;
    residuals(k)=sqrt(Theta_rank_k_eqn_difference'*Theta_rank_k_eqn_difference);
    [Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_rank_k,Theta);
    %Yt_minus_Yb_rank_k_difference_norm(k)=sqrt((Yt-Yt_orig)'*(Yt-Yt_orig));
    %Yt_minus_Yb_rank_k_difference_norm(k)=sqrt(trapz(Phi_t,((Yt-Yt_orig).*(conj(Yt-Yt_orig))).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*D);
end



[~,Ydiff_mindex]=min(log10(Yt_minus_Yb_rank_k_difference_norm)) 

q=3;
[geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q);

clear vars U V Theta_b_coeff_matrix 


%%



example_rank_1=50;
%example_rank_2=125;    

if example_rank_1>maxrank
    example_rank_1=maxrank;
    disp('Example truncation rank larger than system so it has been reduced')
end 


%[Yt_at_Ydiff_mindex,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,Ydiff_mindex),Theta);
[Yt_example,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank_1),Theta);


figure(11); clf;
plot(Phi_t,Theta_b(:,example_rank_1))

% Topography Profiles
figure(1); clf; hold on;
rankname_1="$y_b$ : $\kappa="+num2str(example_rank_1)+"$";
plot(Phi_t,Yt_true,'-r',LineWidth=2,DisplayName='$y_T$')
plot(Phi_t,Yt_example,'-k',LineWidth=2,DisplayName=rankname_1)
try
    rankname_2="$\kappa="+num2str(example_rank_2)+"$";
    [Yt_example_2,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank_2),Theta);
    plot(Phi_t,Yt_example_2,'-b',LineWidth=2,DisplayName=rankname_2)
end
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$y_b$',Interpreter='latex',FontSize=18)
legend(Location='east',Interpreter="latex",FontSize=14)
box on
xlim([Phi_t(1),Phi_t(end)])


% L curve
figure(2);clf;hold on; box on;
plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
plot(log10(residuals(example_rank_1)),log10(Theta_rank_k_norm(example_rank_1)),'xr',MarkerSize=12,LineWidth=1)




% DPC plot
figure(3); clf; hold on; box on
plot(log10(singularvals),DisplayName='$\log_{10}(\sigma_{i})$',LineWidth=1.5)
plot(log10(Fourier_coeffs),DisplayName='$\log_{10}(|$\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{b}$}$|)$',LineWidth=1.5)
plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='$\log_{10}({\rho_i})$',LineWidth=1.5)
try
    yline(log10(noise_amp),'--k',Interpreter='latex',FontSize=18,DisplayName='$\log_{10}(\epsilon)$')
end
xlabel('$i$',Interpreter='latex',FontSize=18)
legend(Interpreter="latex",FontSize=14,Location='northeast')
xlim([0,200])



% norm rank plot
figure(4); clf; hold on; box on;
plot(log10(Theta_rank_k_norm),'-k',LineWidth=2)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
xlabel('$\kappa$',Interpreter='latex',FontSize=18)









% 
% figure(80); clf; hold on;
% plot(log10(residuals))
% plot(example_rank_1,log10(residuals(example_rank_1)),'-*r',MarkerSize=10)
% plot(Ydiff_mindex,log10(residuals(Ydiff_mindex)),'-xm',MarkerSize=10)
% %plot(Theta_mindex,log10(Theta_rank_k_eqn_difference_norm(Theta_mindex)),'-xk',MarkerSize=10)
% ylabel('log10(Ax-b norm)')
% xlabel('rank')






figure(82); clf; hold on; box on
plot(Phi_s,Ys,'-k',LineWidth=2,DisplayName="Perturbed")
plot(Phi_s,Ys_unperturbed,'-b',LineWidth=2,DisplayName="Unperturbed")
%plot(Phi_t,Yt_example,'-k',DisplayName='$y_b(\phi)$',LineWidth=2)
%plot(Phi_t,Yt_true,'-r',DisplayName='$y_T(\phi)$',LineWidth=2)
ylabel('$y_{f}(\phi)$',Interpreter='latex',FontSize=18)
xlabel('$\phi$',Interpreter='latex',FontSize=18)
legend(Interpreter="latex",FontSize=14,Location='east')
xlim([Phi_t(1),Phi_t(end)])
%ylim([(floor(10*(min(Yt_example)))-1)/10,(floor(10*(max(Ys)))+1)/10])










% L curve
figure(2);clf;hold on; box on;
plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
plot(log10(residuals(example_rank_1)),log10(Theta_rank_k_norm(example_rank_1)),'xr',MarkerSize=12,LineWidth=1)




% % DPC plot and norm rank plot
% figure(903); clf;  box on
% subplot(2,1,1);  box on; hold on;
% plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
% xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
% ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% plot(log10(residuals(example_rank)),log10(Theta_rank_k_norm(example_rank)),'xr',MarkerSize=12,LineWidth=1)
% subplot(2,1,2); box on; hold on;
% plot(log10(Theta_rank_k_norm),'-k',LineWidth=2)
% plot(example_rank,log10(Theta_rank_k_norm(example_rank)),'xr',MarkerSize=12,LineWidth=2)
% ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% xlabel('$\kappa$',Interpreter='latex',FontSize=18)





filename="noisy_"+loadname+"_epsilon1E-3";


if saveplots==1
    saveas(82,filename+"_surfaces.eps",'epsc')
    saveas(82,filename+"_surfaces",'fig')
    %
    saveas(1,filename+"_topography_comparison.eps",'epsc')
    saveas(1,filename+"_topography_comparison",'fig')
    %
    saveas(2,filename+"_Lcurve.eps",'epsc')
    saveas(2,filename+"_Lcurve",'fig')
    %
    saveas(3,filename+"_DPC.eps",'epsc')
    saveas(3,filename+"_DPC",'fig')
    %
    saveplots=saveplots+1
    disp('The figures have been saved, saveplots has been changed to prevent accidentally overwriting.')
end
