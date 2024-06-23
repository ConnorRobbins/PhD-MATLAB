guessfile=load("unforced_solitary_guessN1641.mat");

N=641;


F_max=1.2;
F_min=1.01;
FN_up=10;
FN_down=60;

L=guessfile.L;
Phi=linspace(-L,L,N);
P=zeros(N,1);
amplitude=0; b_width=1;
ijac=1;

initial_Theta=interp1(guessfile.Phi,guessfile.Theta_Newton,Phi);
initial_Theta_b=interp1(guessfile.Phi,guessfile.Theta_bottom_Newton,Phi);
initial_F=guessfile.Froude;


Theta_Newton=initial_Theta;
Theta_bottom_Newton=initial_Theta_b;


Fup=linspace(initial_F,F_max,FN_up);
max_y_upper=zeros(1,FN_up);
Fdown=linspace(initial_F,F_min,FN_down);
max_y_down=zeros(1,FN_down);

for iF=1:FN_up
    F=Fup(iF);
    disp("Now beginning step "+num2str(iF)+"/"+num2str(FN_up)+" of moving upwards")
   [~,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,~] = Newton_forward_fun(L,N,F,amplitude,b_width,ijac,P,Theta_Newton,Theta_bottom_Newton);
    max_y_upper(iF)=max(Ys_Newton);
end


Theta_Newton=initial_Theta;
Theta_bottom_Newton=initial_Theta_b;

for iF=1:FN_down
    F=Fdown(iF)
    disp("Now beginning step "+num2str(iF)+"/"+num2str(FN_down)+" of moving downwards")
   [~,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,~] = Newton_forward_fun(L,N,F,amplitude,b_width,ijac,P,Theta_Newton,Theta_bottom_Newton);
    max_y_down(iF)=max(Ys_Newton);
end




figure(1); clf; hold on; box on;
plot(Fup,max_y_upper,'-xb')
plot(Fdown,max_y_down,'-xr')


%[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess);

