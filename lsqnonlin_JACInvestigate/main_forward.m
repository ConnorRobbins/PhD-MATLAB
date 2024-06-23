% clc
% clear all
% close all

makeplots='no';


%%%%%%
% iload: 1= load saved data "saved.mat",
%        2= load only Theta value to look for unforced solitary wave
%                   %% might have to restore from back-up copy first
%                   %% "Forunforced_thetaguess_F1,1_backup.mat"
% isave: 1=save data at end
% ijac:  1= calculate jacobian every loop
% Pressure: 1= include non zero pressure pressure
% Topography: 1= include non zero topography
%%%%%




[Phi,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom,P,jac_N,jac_sol]=eqsolve_forward(L_width,N_points,Froude,iload,isave,ijac,Pressure,Topography,amplitude,b_width);

%%
%figure


if strcmp(makeplots,'yes') == 1;

    if Pressure == 0;

        subplot(2,1,1)
        %plot(Xs,Ys-1,'-ob')
        plot(Phi,Ys-1,'-ob')
        %ylim([-0.1,0])
        hold on
        %plot(Xt,Yt,'-or')
        plot(Phi,Yt,'-or')
        legend('Surface-1','Topography')
        %xlabel('X')
        ylabel('$y$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(2,1,2)
        %plot(Xs,Theta_surface,'-xb')
        plot(Phi,Theta_surface,'-xb')
        hold on
        %plot(Xt,Theta_bottom,'-xr')
        plot(Phi,Theta_bottom,'-xr')
        legend('Surface','Topography')
        %title('Forward Problem Solver')
        %xlabel('X')
        xlabel('$\phi$','interpreter','latex','FontSize',16)
        ylabel('$\theta$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

    end

    if Pressure == 1;

        subplot(3,1,1)
        %plot(Xs,Ys-1,'-ob')
        plot(Phi,Ys-1,'-ob')
        hold on
        %plot(Xt,Yt,'-or')
        plot(Phi,Yt,'-or')
        legend('Surface-1','Topography')
        %xlabel('X')
        ylabel('$y$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(3,1,2)
        %plot(Xs,Theta_surface,'-xb')
        plot(Phi,Theta_surface,'-xb')
        hold on
        %plot(Xt,Theta_bottom,'-xr')
        plot(Phi,Theta_bottom,'-xr')
        legend('Surface','Topography')
        %title('Forward Problem Solver')
        %xlabel('X')
        ylabel('$\theta$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(3,1,3)
        plot(Phi,P)
        ylabel('$P$','interpreter','latex','FontSize',16)
        xlabel('$\phi$','interpreter','latex','FontSize',16)

    end

end




