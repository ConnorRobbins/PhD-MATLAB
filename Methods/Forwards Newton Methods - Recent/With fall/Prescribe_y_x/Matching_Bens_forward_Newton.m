clear 

L=20;
N=161;


Phi = linspace(-L,L,N);

%parameters
bens_h=0.3;
b_width=0.3;
ijac=1;


amplitude=sign(bens_h)*(bens_h^2)*((2*(b_width^2)/pi)^0.25);



%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses

%Theta_init_guess=zeros(N,1);
Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
Froude_init_guess=1.3;
Theta_b_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative





[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun_hydraulic_fall(L,N, Froude_init_guess,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess);


%%
%Check areas

myarea=sqrt(trapz(Xt_Newton,Yt_Newton.^2))
testarea=sqrt(trapz(Xt_Newton,((amplitude*exp(-(b_width*Xt_Newton).^2)).^2)))
bensarea=bens_h^2


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
plot(Xt_Newton,Yt_Newton)
plot(Xt_Newton,amplitude*exp(-(b_width*Xt_Newton).^2))
