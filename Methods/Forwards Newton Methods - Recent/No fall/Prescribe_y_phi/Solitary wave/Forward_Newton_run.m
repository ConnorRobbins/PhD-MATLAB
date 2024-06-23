clear 

L=20;
N=641;


Phi = linspace(-L,L,N);

%parameters
amplitude=0.0;
b_width=1;
Froude=1.1;
ijac=1;




%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);



%initial guesses

%Theta_init_guess=zeros(N,1);
%Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
%Theta_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
%Theta_init_guess=load('Forunforced_Thetaguess_N41_F1,1_backup.mat').Theta;
% Theta_init_guess=load('unforced_solitary_guess2.mat').Theta;
% Phi_old=load('unforced_solitary_guess2.mat').Phi;

Theta_init_guess=load('unforced_solitary_guessN1641.mat').Theta_Newton;
Phi_old=load('unforced_solitary_guessN1641.mat').Phi;
Theta_init_guess=interp1(Phi_old,Theta_init_guess,Phi);

%Theta_b_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
Theta_b_init_guess=zeros(N,1);




[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess);


%%

figure(11); clf; hold on;
plot(Phi,Ys_Newton,'-b')
plot(Phi,Yt_Newton,'-k')
% 
% %figure(13); clf;hold on;
% plot(Xs_Newton,Ys_Newton,'--b')
% plot(Xt_Newton,Yt_Newton,'--k')
% 
% 
% 
% 
% figure(12); clf; hold on;
% plot(Phi,Theta_bottom_Newton)
% ylabel('Theta bottom')
% xlabel('Phi')
% 
% 
% figure(13); clf; hold on;
% plot(Phi,Yt_Newton)
% plot(Phi, amplitude*exp(-(b_width*Phi).^2),'--r')


figure(1); clf;
plot(Phi,Ys_Newton)