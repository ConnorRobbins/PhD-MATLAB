clear 

L=41;
N=161;


Phi = linspace(-L,L,N);

%parameters
amplitude=-0.05;
b_width=0.3;
Froude=0.8;
ijac=1;




%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses

Theta_init_guess=zeros(N,1);
%Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
%Theta_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative

Theta_b_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative





[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess);


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

