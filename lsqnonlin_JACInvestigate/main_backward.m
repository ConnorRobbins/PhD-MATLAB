% 
% if FreeSurfaceType==1
%     forcingtype='Gaussian';
% elseif FreeSurfaceType==2
%     forcingtype='Sech^2';
% else
%     forcingtype='';
% end

% titlestring=strcat(forcingtype, ' amplitude=',num2str(amplitude),', b=',num2str(b_width),', N=',num2str(N_points),', L=',num2str(L_width),' F=',num2str(Froude))


[Phi,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom]=eqsolve_backward(L_width,N_points,Froude,iload,isave,ijac,Pressure,FreeSurfaceType,amplitude,b_width);


% figure(1)
% subplot(2,1,1)
% plot(Phi,Ys-1,'-ob')
% hold on
% plot(Phi,Yt,'-or')
% legend('Surface-1','Topography')
% xlabel('X')
% ylabel('y')
% xlim([Xt(1) Xt(end)])
% 
% subplot(2,1,2)
% plot(Xs,Theta_surface,'-xb')
% hold on
% plot(Xt,Theta_bottom,'-xr')
% legend('Surface','Topography')
% title('Backward Problem Solver')
% xlabel('X')
% ylabel('theta')
% xlim([Xt(1) Xt(end)])


% %%
% figure
% subplot(2,1,1)
% plot(Phi,Ys-1,'-ob')
% hold on
% plot(Phi,Yt,'-or')
% legend('Surface-1','Topography')
% xlabel('$\phi$','Interpreter','Latex','FontSize',16)
% ylabel('$y$','Interpreter','Latex','FontSize',18,'Rotation',0)
% %title(titlestring)
% 
% 
% subplot(2,1,2)
% plot(Phi,Theta_surface,'-xb')
% hold on
% plot(Phi,Theta_bottom,'-xr')
% legend('Surface','Topography')
% %title('Inverse Problem Solver')
% xlabel('$\phi$','Interpreter','Latex','FontSize',16)
% ylabel('$\theta$','Interpreter','Latex','FontSize',18,'Rotation',0)


