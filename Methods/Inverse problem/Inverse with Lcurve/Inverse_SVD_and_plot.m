%% Physical parameters 
clear

%filename="forward_backward_"+"hydraulicfall_soln_a1E-1_N721";

saveplots=0;
SolveForwardAtEnd=1;

ijac=1;

L=20; %domain truncation
Froude=1.4;
amplitude=0.1;  %surface parameter, see form of Ys
b_width=0.3;    %surface parameter, see form of Ys
topo_centre=0;  %surface parameter, see form of Ys

N_t=181;         %Number of points on topography
N_s=N_t;        %Number of points on surface

Phi_t=linspace(-L,L,N_t); %Creates topography mesh
Phi_s=linspace(-L,L,N_s); %Creates surface mesh

P=zeros(1,N_s);     %Pressure on surface
%P=amplitude*exp(-(b_width*(Phi_s-topo_centre)).^2);


%Ys=1+amplitude*exp(-(b_width*(Phi_s-topo_centre)).^2);
%Ys=1+amplitude*(1-tanh(b_width*Phi_s)).*cos(Phi_s)/2;
%Ys=1+amplitude*(sech(b_width*Phi_s)).^2;
%Ys=1+amplitude*(exp(-(b_width*Phi_s).^2)).*cos(Phi_s)/2;

y0=( ((Froude^2)/2)  +   sqrt( ( (Froude^4)/4)+2*Froude^2))/2; c1=0.5*(1+y0); c2=c1-1; Ys=c1+c2*tanh(-Phi_s);


%noise_amp=1E-6; noise=noise_amp*2*(randn(N_s,1)-0.5)'; Ys=Ys+noise;


[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);



%%

Theta_rank_k_norm=zeros(maxrank,1);
residuals=zeros(maxrank,1);
Yt_minus_Yb_rank_k_difference_norm=zeros(maxrank,1);
Fourier_coeffs=zeros(maxrank,1);



for k=1:maxrank
    Theta_rank_k=Theta_b(:,k);
    Theta_rank_k_norm(k)=sqrt(trapz(Phi_t,((Theta_rank_k).*(conj(Theta_rank_k))).^2));
    Theta_rank_k_eqn_difference=Theta_b_coeff_matrix*Theta_rank_k - D;
    residuals(k)=sqrt(Theta_rank_k_eqn_difference'*Theta_rank_k_eqn_difference);
    [Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_rank_k,Theta);
    %Yt_minus_Yb_rank_k_difference_norm(k)=sqrt(trapz(Phi_t,((Yt-Yt_orig).*(conj(Yt-Yt_orig))).^2));
    Fourier_coeffs(k)=abs(transpose(U(:,k))*D);
end



%[Y_diff_min_value_log,Ydiff_mindex]=min(log10(Yt_minus_Yb_rank_k_difference_norm)) 

q=3;
[geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q);

clear vars U V Theta_b_coeff_matrix


%%



example_rank_1=160;
%example_rank_2=166;

%[Yt_at_Ydiff_mindex,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,Ydiff_mindex),Theta);
[Yt_example,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank_1),Theta);

%% Solve forward if asked for 

if SolveForwardAtEnd==1
    [Phi,Xs_forward,Ys_forward,Xt_forward,Yt_forward,Theta_forward,Theta_bottom_forward,~] = Newton_forward_fun(L,N_t,Froude,ijac,P,Theta,Theta_b(:,example_rank_1),Yt_example);
end

%%







% Topography Profiles
figure(1); clf; hold on;
rankname_1="$\kappa="+num2str(example_rank_1)+"$";
%plot(Phi_t,Yt_orig,'-r',LineWidth=2,DisplayName='$y_T$')
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
legend(Interpreter="latex",FontSize=14,Location='southwest')



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
plot(Phi_s,Ys,'-b',LineWidth=2,DisplayName='$y_f(\phi)$')
%plot(Phi_t,Yt_example,'-k',DisplayName='$y_b(\phi)$',LineWidth=2)
%plot(Phi_t,Yt_orig,'-r',DisplayName='$y_T(\phi)$',LineWidth=2)
ylabel('$y(\phi)$',Interpreter='latex',FontSize=18)
xlabel('$\phi$',Interpreter='latex',FontSize=18)
legend(Interpreter="latex",FontSize=14,Location='east')
xlim([Phi_t(1),Phi_t(end)])
%ylim([(floor(10*(min(Yt_example)))-1)/10,(floor(10*(max(Ys)))+1)/10])
plot(Phi,Ys_forward,'--r')









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




figure(6); clf; hold on;
plot(Phi_t,Yt_example)
plot(Phi_s,Ys_forward)







if saveplots==1
    saveas(82,filename+"_surf_and_topo_profiles.eps",'epsc')
    saveas(82,filename+"_surf_and_topo_profiles",'fig')
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
