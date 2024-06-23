clear 

N=1081;
filename='subcrit_guess';



s=load(filename);
L=s.L;
amplitude=s.amplitude;
b_width=s.b_width;
Froude=s.Froude;
Phi_old=s.Phi;




Phi = linspace(-L,L,N);
ijac=1;




%forcings
Yt=amplitude*exp(-(b_width*Phi).^2);
%Yt=amplitude*exp(-(b_width*(Phi-10)).^2) + amplitude*exp(-(b_width*(Phi+10)).^2);

P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses

Theta_init_guess=interp1(Phi_old,s.Theta_Newton,Phi);

Theta_b_init_guess=interp1(Phi_old,s.Theta_bottom_Newton,Phi);





[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,ijac,P,Theta_init_guess,Theta_b_init_guess,Yt);


%%

figure(11); clf; hold on;
plot(Phi,Ys_Newton,'-b')
plot(Phi,Yt_Newton,'-k')

%figure(13); clf;hold on;
plot(Xs_Newton,Ys_Newton,'--b')
plot(Xt_Newton,Yt_Newton,'--k')




figure(12); clf; hold on;
plot(Phi,Theta_bottom_Newton)
ylabel('Theta bottom')
xlabel('Phi')


figure(13); clf; hold on;
plot(Phi,Ys_Newton)

