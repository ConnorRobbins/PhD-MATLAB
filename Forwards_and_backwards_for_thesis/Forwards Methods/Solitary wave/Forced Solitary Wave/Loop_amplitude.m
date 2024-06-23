clear 


%parameters

N=241;


ijac=1;




%guessfile="unforced_solitary_guessN1641"
%guessfile="forced_amp_1E-2_N181"
%guessfile="forced_solitary_a2E-2_F1p1_N181"
guessfile="a2E-2"


max_amplitude=0.03;
amplitude_N=20;







s=load(guessfile,'Theta_Newton','Theta_bottom_Newton','Phi','amplitude','Froude','b_width','L');
amplitudes=linspace(s.amplitude,max_amplitude,amplitude_N);
%amplitudes=[0.01]; amplitude_N=1;

L=s.L;
Phi = linspace(-L,L,N);




%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


Theta_Newton=s.Theta_Newton;
Theta_bottom_Newton=s.Theta_bottom_Newton;
Phi_old=s.Phi;
Theta_Newton=interp1(Phi_old,Theta_Newton,Phi);
Theta_bottom_Newton=interp1(Phi_old,Theta_bottom_Newton,Phi);


Froude=s.Froude;
b_width=s.b_width;






for amp_i=1:amplitude_N
    amp_i
    amplitude=amplitudes(amp_i);
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