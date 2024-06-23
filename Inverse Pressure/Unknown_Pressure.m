clc
%close
clear all
%cd ('/home/stusci1/rgq13jzu/MATLAB')
%cd('/local/scratch/rgq13jzu/Dropbox/Connor Notes 2018/matlab/mark_results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 8; %Phi range on either side of the origin
N = 81;
iload = 0;
isave = 0;
F=1.2;
Pzero_origin=0;   %%1=sets P=0 at origin


[Phi,X,Y,P,Theta] = pressurefinder_fixedFTheta(L,N,F,iload,isave,Pzero_origin);

%%  Plotting

%%%%%%%%%%%%%%%%%%%%%%%%%   my plotting
% subplot(3,1,2)
% plot(Phiold,Pold,'-ob')
% hold on
% if Ntemp ~= Nold
%     plot(Phi,Poldspline,'-xr')
%     legend('Exact input','Interpolated input')
% else
%     legend('Exact input')
% end
% ylabel('Pressure input')
% subplot(3,1,1)
% plot(Phi,Pinit)
% ylabel('P initial guess')
% subplot(3,1,3)
% plot(Phi,P,'-xr')
% xlabel('Phi')
% ylabel('Converged pressure')
%  
% figure
% plot(Phiold,Pold,'--rx')
% hold on
% plot(Phi,P,'k')
% if Ntemp ~= Nold
%     plot(Phi,Poldspline,'bx')
%     legend('actual pressure','converged pressure','interpolated old pressure')
% else
%     legend('actual pressure','converged pressure')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Mark's plotting
figure
subplot(3,1,2)
plot(Phi,Y,'b')
hold on
ylabel('Wave shape')
subplot(3,1,1)
plot(Phi,Theta,'b')
ylabel('Theta')
subplot(3,1,3)
plot(Phi,P,'-b')
xlabel('Phi')
ylabel('Converged pressure')
%%%%%%%%%%%%%%%%%%%%%%%%%%
