clear 


%parameters
L=20;
N=241;

amplitude=0;
b_width=1;
ijac=1;




guessfile="unforced_solitary_guessN1641"


max_Froude=1.2;
Froude_N=10;
Froudes=linspace(1.1,max_Froude,Froude_N);



Phi = linspace(-L,L,N);




%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


Theta_Newton=load(guessfile).Theta_Newton;
Theta_bottom_Newton=load(guessfile).Theta_bottom_Newton;
Phi_old=load(guessfile).Phi;
Theta_Newton=interp1(Phi_old,Theta_Newton,Phi);
Theta_bottom_Newton=interp1(Phi_old,Theta_bottom_Newton,Phi);









for F_i=1:Froude_N
    F_i
    Froude=Froudes(F_i);
    [Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_Newton,Theta_bottom_Newton);
end

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