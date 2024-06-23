clear
saveplots=0;

%% Physical parameters 

%N=641; L=20; amplitude=0.2; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[187,152,140,121,119,110];
%N=641; L=20; amplitude=-0.2; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[114,99,101,102,99,96];

%N=641; L=20; amplitude=0.15; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[153,131,124,114,113,104];
%N=641; L=20; amplitude=-0.15; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[105,102,101,99,98,97];

%N=641; L=20; amplitude=0.1; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[127,115,108,102,101,95];
N=641; L=20; amplitude=-0.1; b_width=0.3; Froudes=[0.8,0.9,1,1.1,1.2,1.3]; truncation_ranks=[99,97,96,96,93,92];



truncation_ranks=[90];


[N_t,N_s]=deal(N);
Phi_s=linspace(-L,L,N_s); Phi_t=linspace(-L,L,N_t);

P=zeros(1,N_s);
Ys=1+amplitude*exp(-(b_width*Phi_s).^2);

Theta_approx_for_pressure=atan(-2*amplitude*(b_width^2)*Phi_s.*exp(-(b_width*Phi_s).^2))';



loopnumber=numel(Froudes);
steps=num2str(loopnumber);


Yt_store=zeros(N_t,loopnumber);
P_store=zeros(N_s,loopnumber);


if numel(truncation_ranks)==1
    truncation_ranks=ones(1,loopnumber)*truncation_ranks;
end

for i=1:loopnumber
    %
    disp("Starting step "+num2str(i)+"/"+steps)
    Froude=Froudes(i);
    truncation_rank=truncation_ranks(i);

    %% %% Perform the topography calculations
    %Form the truncated pseudoinverse
    [truncated_psuedoinv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_rank);
    % form the RHS of the matrix eqn and output Theta_f
    [D,Theta] = FUN_make_SVD_RHS(Ys,P,Froude,L,N_s,N_t);
    %Solve the truncated system
    Theta_b=truncated_psuedoinv*D;
    [Yt_store(:,i),~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b,Theta);


    %% %% Perform the pressure calculations
    [~,~,~,~,~,~,~,inverse_Pressure] = FUN_Inverse_Pressure_NO_TOPOGRAPHY(L,N_s,N_t,Froude,Theta_approx_for_pressure);
    P_store(:,i)=inverse_Pressure;

end






%%





legendstring="$F="+num2str(Froudes')+"$";
% yaxismax=max(max(max(P_store)),max(max(Yt_store)));
% yaxismin=min(min(min(P_store)),min(min(Yt_store)));


figure(436); clf; hold on; box on;
plot(Phi_t,Yt_store,LineWidth=2)
legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$y_b$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])



figure(437); clf; hold on; box on;
plot(Phi_s,P_store,LineWidth=2)
legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$P$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])


% 
% 
% 
% 
% % Topography Profiles
% figure(1); clf; hold on;
% rankname_1="$\kappa="+num2str(example_rank_1)+"$";
% plot(Phi_t,Yt_orig,'-r',LineWidth=2,DisplayName='$y_T$')
% plot(Phi_t,Yt_example,'-k',LineWidth=2,DisplayName=rankname_1)
% try
%     rankname_2="$\kappa="+num2str(example_rank_2)+"$";
%     [Yt_example_2,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank_2),Theta);
%     plot(Phi_t,Yt_example_2,'-b',LineWidth=2,DisplayName=rankname_2)
% end
% xlabel('$\phi$',Interpreter='latex',FontSize=18)
% ylabel('$y_b$',Interpreter='latex',FontSize=18)
% legend(Location='east',Interpreter="latex",FontSize=14)
% box on
% xlim([Phi_t(1),Phi_t(end)])
% 
% 
% % L curve
% figure(2);clf;hold on; box on;
% plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
% xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
% ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% plot(log10(residuals(example_rank_1)),log10(Theta_rank_k_norm(example_rank_1)),'xr',MarkerSize=12,LineWidth=1)
% 
% 
% 
% 
% % DPC plot
% figure(3); clf; hold on; box on
% plot(log10(singularvals),DisplayName='$\log_{10}(\sigma_{i})$',LineWidth=1.5)
% plot(log10(Fourier_coeffs),DisplayName='$\log_{10}(|$\mbox{\boldmath$\underline{u}$}$_{i}^{T}$\mbox{\boldmath$\underline{b}$}$|)$',LineWidth=1.5)
% plot(geo_mean_indices,log10(geo_mean(geo_mean_indices)),DisplayName='$\log_{10}({\rho_i})$',LineWidth=1.5)
% try
%     yline(log10(noise_amp),'--k',Interpreter='latex',FontSize=18,DisplayName='$\log_{10}(\epsilon)$')
% end
% xlabel('$i$',Interpreter='latex',FontSize=18)
% legend(Interpreter="latex",FontSize=14,Location='southwest')
% 
% 
% 
% % norm rank plot
% figure(4); clf; hold on; box on;
% plot(log10(Theta_rank_k_norm),'-k',LineWidth=2)
% ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% xlabel('$\kappa$',Interpreter='latex',FontSize=18)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% % figure(80); clf; hold on;
% % plot(log10(residuals))
% % plot(example_rank_1,log10(residuals(example_rank_1)),'-*r',MarkerSize=10)
% % plot(Ydiff_mindex,log10(residuals(Ydiff_mindex)),'-xm',MarkerSize=10)
% % %plot(Theta_mindex,log10(Theta_rank_k_eqn_difference_norm(Theta_mindex)),'-xk',MarkerSize=10)
% % ylabel('log10(Ax-b norm)')
% % xlabel('rank')
% 
% 
% 
% 
% 
% 
% figure(82); clf; hold on; box on
% plot(Phi_s,Ys,'-b',LineWidth=2,DisplayName='$y_f(\phi)$')
% %plot(Phi_t,Yt_example,'-k',DisplayName='$y_b(\phi)$',LineWidth=2)
% plot(Phi_t,Yt_orig,'-r',DisplayName='$y_T(\phi)$',LineWidth=2)
% ylabel('$y(\phi)$',Interpreter='latex',FontSize=18)
% xlabel('$\phi$',Interpreter='latex',FontSize=18)
% legend(Interpreter="latex",FontSize=14,Location='east')
% xlim([Phi_t(1),Phi_t(end)])
% ylim([(floor(10*(min(Yt_example)))-1)/10,(floor(10*(max(Ys)))+1)/10])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % L curve
% figure(2);clf;hold on; box on;
% plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
% xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
% ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% plot(log10(residuals(example_rank_1)),log10(Theta_rank_k_norm(example_rank_1)),'xr',MarkerSize=12,LineWidth=1)
% 
% 
% 
% 
% % % DPC plot and norm rank plot
% % figure(903); clf;  box on
% % subplot(2,1,1);  box on; hold on;
% % plot(log10(residuals),log10(Theta_rank_k_norm),'-k','LineWidth',2)
% % xlabel('$\log_{10}(|$\mbox{\boldmath$M\underline{\theta_{\kappa}}$}$-$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
% % ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% % plot(log10(residuals(example_rank)),log10(Theta_rank_k_norm(example_rank)),'xr',MarkerSize=12,LineWidth=1)
% % subplot(2,1,2); box on; hold on;
% % plot(log10(Theta_rank_k_norm),'-k',LineWidth=2)
% % plot(example_rank,log10(Theta_rank_k_norm(example_rank)),'xr',MarkerSize=12,LineWidth=2)
% % ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{\theta_{\kappa}}$}$|)$',Interpreter='latex',FontSize=18)
% % xlabel('$\kappa$',Interpreter='latex',FontSize=18)
% 
% 
% 
% 
% 
% 
% 
% 
% if saveplots==1
%     saveas(82,filename+"_surf_and_topo_profiles.eps",'epsc')
%     %saveas(82,filename+"_surf_and_topo_profiles",'fig')
%     saveas(82,filename+"_surf_and_topo_profiles"+".fig")
%     %
%     saveas(1,filename+"_topography_comparison.eps",'epsc')
%     %saveas(1,filename+"_topography_comparison",'fig')
%     saveas(1,filename+"_topography_comparison"+".fig")
%     %
%     saveas(2,filename+"_Lcurve.eps",'epsc')
%     %saveas(2,filename+"_Lcurve",'fig')
%     saveas(2,filename+"_Lcurve"+".fig")
%     %
%     saveas(3,filename+"_DPC.eps",'epsc')
%     %saveas(3,filename+"_DPC",'fig')
%     saveas(3,filename+"_DPC"+".fig")
%     %
%     saveplots=saveplots+1
%     disp('The figures have been saved, saveplots has been changed to prevent accidentally overwriting.')
% end
