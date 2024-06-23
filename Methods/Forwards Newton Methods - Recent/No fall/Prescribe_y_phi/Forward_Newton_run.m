clear 

L=20;
N=641;


Phi = linspace(-L,L,N);

%parameters
amplitude=0.0;
b_width=0.3;
Froude=1.2;
ijac=1;




%forcings
%Yt=amplitude*exp(-(b_width*Phi).^2);
Yt=amplitude*exp(-(b_width*Phi).^2).*cos(Phi);
%Yt=amplitude*exp(-(b_width*(Phi-10)).^2) + amplitude*exp(-(b_width*(Phi+10)).^2);

P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses

Theta_init_guess=zeros(N,1);
%Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
%Theta_init_guess=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative

%Theta_b_init_guess=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
Theta_b_init_guess=atan([(Yt(2:end)-Yt(1:end-1))/(Phi(2)-Phi(1)),0]);




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

