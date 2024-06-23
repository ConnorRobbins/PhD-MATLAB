clear 


%parameters
N=1241;

ijac=1;




%guessfile="unforced_solitary_guessN1641";
guessfile="a2E-2_F1p3_N841";

s=load(guessfile);

max_Froude=1.3;
Froude_N=1;
Froudes=linspace(s.Froude,max_Froude,Froude_N);


L=s.L;
Phi = linspace(-L,L,N);




%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


amplitude=s.amplitude;
b_width=s.b_width;
Theta_Newton=s.Theta_Newton;
Theta_bottom_Newton=s.Theta_bottom_Newton;
Phi_old=s.Phi;
clear s
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