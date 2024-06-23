clear;
clc;



%parameters
amplitude=0.001;
b_width=1;
L=20;
N=641;

ijac=1;
P=zeros(N,1);



%%%%%%%% Lower branch: Loop over Froude numbers to find max Y
N_lower=40;
F_start=1.2;
F_end=1.1;
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Upper branch: Loop over max Y to find Froude numbers
N_upper =80;  % Number of points on upper branch
%max_y_start will be set to be final point of lower branch
max_y_end = 1.25;  %End maxY value number
%%%%%%%%%%%%%%




Phi = linspace(-L,L,N);
Yt_true=amplitude*exp(-(b_width*Phi).^2);




%% LOWER BRANCH 


if N_lower ==1
    dF = 0;
else
    dF = (F_end-F_start)/(N_lower-1);
end

max_y_lower=zeros(1,N_lower);
F_lower=zeros(1,N_lower);
Theta_Newton=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%initial guess of roughly gaussian hump derivative
Theta_bottom_Newton=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%initial guess of roughly gaussian hump derivative


for iF=1:N_lower
    disp("Lower branch step "+num2str(iF)+"/"+num2str(N_lower))
    Froude = F_start + (iF-1)*dF;

    [~,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,~] = Newton_forward_fun(L,N,Froude,ijac,P,Theta_Newton,Theta_bottom_Newton,Yt_true);

    %[Phi,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom]=surfacefinder_topography_fixedFP(L,N,Froude,iload,isave,ijac,Pressure,Topography);
    max_y_lower(iF) = max(Ys_Newton);
    F_lower(iF)  = Froude;
end








%% UPPER BRANCH

max_y_start=max(Ys_Newton);
F=Froude;
dmx = (max_y_end-max_y_start)/(N_upper-1);
if N_upper==1 dmx=0; end;



max_y_upper=zeros(1,N_upper);
F_upper=zeros(1,N_upper);

iF=1;



for iF=1:N_upper
    disp("Upper branch step "+num2str(iF)+"/"+num2str(N_upper))
    max_y = max_y_start + (iF-1)*dmx;
    [~,Xs_Newton,Ys_newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,F,~] = Newton_forward_fun_fixed_max_y(L,N,max_y,Theta_Newton,Theta_bottom_Newton,F,Yt_true,P,ijac);
    max_y_upper(iF)=max(Ys_newton);
    F_upper(iF)=F;
end





%%


figure(1); %clf; hold on;    
plot(F_lower,max_y_lower,'-xr')
plot(F_upper,max_y_upper,'-*','markersize',1)
xlabel('Froude number','FontSize',16,'Interpreter','Latex')
ylabel('Wave peak height','FontSize',16,'Interpreter','Latex')

