clear 

L=20;
N=561;


Phi = linspace(-L,L,N);

%parameters
amplitude=0.086;
b_width=1;
Froude=1.2;
ijac=1;

kdv_solve_for_init_guess=0; %1=yes do this
load_for_init_guess=1;    guessfilename='testloading';   %



%forcings
P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%initial guesses
if kdv_solve_for_init_guess==1
    %%%% KDV initial guess
    del_Phi=Phi(2)-Phi(1);
    [eta_solved] = FUNCTION_initial_guess_from_kdv(N,amplitude*exp(-(b_width*Phi).^2),Froude,del_Phi,ijac);
    Ys_kdv=1+eta_solved;
    dYs_kdv_dPhi=[0;((Ys_kdv(2:end)-Ys_kdv(1:end-1))/del_Phi)];
    Theta_init_guess=atan(dYs_kdv_dPhi);
    Theta_b_init_guess=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
elseif load_for_init_guess==1
    s=load(guessfilename,'Theta_Newton','Theta_bottom_Newton');
    Theta_init_guess=s.Theta_Newton;
    Theta_b_init_guess=s.Theta_bottom_Newton;
    clearvars s
else
    Theta_init_guess=zeros(N,1);
    %Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
    %Theta_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative

    Theta_b_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative

end




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


figure(14); %clf; hold on;
plot(Phi,Ys_Newton)

