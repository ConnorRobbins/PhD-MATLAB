clear 

L=30;
%N=641;

Phi = linspace(-L,L,N);

%parameters
Froude=0.8;
ijac=1;
%%%% forcing of form f(x)= a*exp(-((b*(x-c)).^2)) + a*exp(-((b*(x+c)).^2))
amplitude=-0.04;
b_width=1;
c_init_guess=4;
%c_init_guess=3.268632740853898;
%c_init_guess=2*3.133387698844217
%c_init_guess=3*3.268632740853898;



%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses

Theta_init_guess=zeros(N,1);
%Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
%Theta_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative

Theta_b_init_guess=atan(-2*amplitude*(b_width^2)*((Phi-c_init_guess).*exp(-(b_width*(Phi-c_init_guess)).^2)   + (Phi+c_init_guess).*exp(-(b_width*(Phi+c_init_guess)).^2))   ) ;





[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude,c] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess,c_init_guess);


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
plot(Phi,Yt_Newton)
plot(Phi, amplitude*exp(-(b_width*Phi).^2),'--r')
