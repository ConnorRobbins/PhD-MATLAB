clc
clear all
close all


L_width = 20; %Phi range on either side of the origin
N_points = 41;
Froude = 1.2;
iload = 0;
isave = 0;
ijac = 0;
Pressure = 0;
Topography = 1;

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




[Phi,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom]=surfacefinder_topography_fixedFP(L_width,N_points,Froude,iload,isave,ijac,Pressure,Topography);



subplot(2,1,1)
plot(Xs,Ys-1,'-ob')
hold on
plot(Xt,Yt,'-or')
legend('Surface-1','Topography')
xlabel('X')
ylabel('y')
xlim([Xt(1) Xt(end)])

subplot(2,1,2)
plot(Xs,Theta_surface,'-xb')
hold on
plot(Xt,Theta_bottom,'-xr')
legend('Surface','Topography')
xlabel('X')
ylabel('theta')
xlim([Xt(1) Xt(end)])




