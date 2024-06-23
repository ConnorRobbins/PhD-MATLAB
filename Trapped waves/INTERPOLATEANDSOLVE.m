clear 

N=141;

saveresult=0


%for ind=1:8
for ind=4

    s=load(strcat('N1241trappedc',num2str(ind),'.mat'));
    L=s.L;
    Froude=s.Froude;
    amplitude=s.amplitude;
    b_width=s.b_width;
    Phi_old=s.Phi;
    c_init_guess=s.c;
    Phi=linspace(-L,L,N);
    Theta_init_guess=interp1(Phi_old,s.Theta_Newton,Phi);
    Theta_b_init_guess=interp1(Phi_old,s.Theta_bottom_Newton,Phi);




    %parameters
    ijac=1;
    %forcings
    P=zeros(N,1); %havent loaded pressure so just make new vector of zeroes
    %P=amplitude*exp(-(b_width*Phi).^2);


    [Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude,c] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess,c_init_guess);

    clear("Phi_old","s")
    if saveresult==1
        save(strcat("N",num2str(N),'trappedc',num2str(ind),'.mat'))
    end
end

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
