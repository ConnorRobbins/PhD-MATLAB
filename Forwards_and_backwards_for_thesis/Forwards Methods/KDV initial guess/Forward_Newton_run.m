clear;

L=20;
%L=18.6429
N=141;


Phi = linspace(-L,L,N);
del_Phi=2*L/(N-1);

%parameters
amplitude_left=-0.04;
amplitude_right=-0.04;
b_width=1;
c_sep=5;
Froude=0.8;
ijac=1;
kdv_solve_for_init_guess=1; %1=yes do this




%forcings

%Yt=amplitude*exp(-(b_width*Phi).^2);
Yt_true=amplitude_right*exp(-(b_width*(Phi-c_sep)).^2) +amplitude_left*exp(-(b_width*(Phi+c_sep)).^2);

P=zeros(N,1);
%P=amplitude*exp(-(b_width*Phi).^2);


%% initial guesses

if kdv_solve_for_init_guess==1
    %%%% KDV initial guess
    [eta_solved] = FUNCTION_initial_guess_from_kdv(N,Yt_true,Froude,del_Phi,ijac);
    Ys_kdv=1+eta_solved;
    dYs_kdv_dPhi=[0;((Ys_kdv(2:end)-Ys_kdv(1:end-1))/del_Phi)];
    Theta_init_guess=atan(dYs_kdv_dPhi);
else
    Theta_init_guess=zeros(N,1);
    %Theta_init_guess=-0.25*((sech(Phi)).^2); %like an init guess for y of(0.25*tanh(-Phi)+1.25)
    %Theta_init_guess=atan(-2*amplitude*b_width*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
    
    Ys_kdv=NaN*ones(N,1);
end





%Theta_b_init_guess=atan(-2*amplitude*(b_width^2)*Phi.*exp(-(b_width*Phi).^2));%init guess of roughly gaussian hump derivative
Theta_b_init_guess=atan([(Yt_true(2:end)-Yt_true(1:end-1))/(Phi(2)-Phi(1)),0]);


%%

[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,ijac,P,Theta_init_guess,Theta_b_init_guess,Yt_true);
%%

figure(11); clf; hold on;
plot(Phi,Ys_Newton,'-b')
plot(Phi,Yt_Newton,'-k')
plot(Phi, Ys_kdv,'--r')

%figure(13); clf;hold on;
plot(Xs_Newton,Ys_Newton,'--b')
plot(Xt_Newton,Yt_Newton,'--k')




figure(12); clf; hold on;
plot(Phi,Theta_bottom_Newton)
ylabel('Theta bottom')
xlabel('Phi')


figure(13); clf; hold on;
plot(Phi,Yt_Newton,DisplayName="Yt Newton")
plot(Phi,Yt_true,'--r',DisplayName="True topo")
legend
