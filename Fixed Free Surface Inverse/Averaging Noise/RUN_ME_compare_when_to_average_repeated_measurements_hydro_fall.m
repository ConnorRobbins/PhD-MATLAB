clear
saveplots=0;
gray=[0.7,0.7,0.7];
rng('default')


%% Physical parameters 

s=load('hydraulicfall_soln_a4E-1_N721.mat','N','L','b_width','Froude','amplitude','Ys_Newton'); Phi=linspace(-s.L,s.L,s.N); s.Yt=s.amplitude*exp(-(s.b_width*Phi).^2);
%s=load('hydraulicfall_soln_a1E-1_N721.mat','N','L','b_width','Froude','amplitude','Ys_Newton'); Phi=linspace(-s.L,s.L,s.N); s.Yt=s.amplitude*exp(-(s.b_width*Phi).^2);

N=s.N;
L=s.L;
b_width=s.b_width; 
Froude=s.Froude;
amplitude=s.amplitude; 
Ys=s.Ys_Newton;
Yt_true=reshape(s.Yt,N,1);

perturbation_N=1000;  %how many different noise stds#
noise_amp=1E-3;

truncation_ranks=[40];


[N_t,N_s]=deal(N);
Phi_s=linspace(-L,L,N_s); Phi_t=linspace(-L,L,N_t);
P=zeros(1,N_s);




loopnumber=perturbation_N;
steps=num2str(loopnumber);

Yt_store=NaN*ones(N_t,loopnumber);
imaginary_counter=0;





if numel(truncation_ranks)==1
    truncation_ranks=ones(1,loopnumber)*truncation_ranks; 
    [truncated_psuedoinv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_ranks(1)); trunc_uniformly_check=1;
else
    trunc_uniformly_check=0;
end



Ys_noisy_store=NaN*ones(N_s,loopnumber);

for i=1:loopnumber
        disp("Starting step "+num2str(i)+"/"+steps)
        Ys_noisy=Ys; noise=noise_amp*2*(randn(N_s-2,1))'; Ys_noisy(2:end-1)=Ys_noisy(2:end-1)+noise;
        Ys_noisy_store(:,i)=Ys_noisy;
        truncation_rank=truncation_ranks(i);
    if trunc_uniformly_check==0
        [truncated_psuedoinv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_rank);
    end



        %% %% Perform the topography calculations
        %Form the truncated pseudoinverse
        % form the RHS of the matrix eqn and output Theta_f
        [D,Theta] = FUN_make_SVD_RHS(Ys_noisy,P,Froude,L,N_s,N_t);
        %Solve the truncated system
        Theta_b=truncated_psuedoinv*D;
        [Yt_one_run,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b,Theta);
        if isreal(Yt_one_run)
            Yt_store(:,i)=Yt_one_run;
        else
            imaginary_counter=imaginary_counter+1;
        end
        

end

Yt_averaged=mean(Yt_store,2,'omitnan');

if imaginary_counter~=0
    if imaginary_counter==1
        disp(num2str(imaginary_counter)+" solution has been rejected for being complex")
    else
        disp(num2str(imaginary_counter)+" solutions have been rejected for being complex")
    end
end





%% Average surfaces and then solve instead
Ys_averaged=mean(Ys_noisy_store,2);
% form the RHS of the matrix eqn and output Theta_f
[D,Theta] = FUN_make_SVD_RHS(Ys_averaged,P,Froude,L,N_s,N_t);
%Solve the truncated system
Theta_b=truncated_psuedoinv*D;
[Yt_from_averaged_Ys,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b,Theta);








%%


figure(1); clf; hold on; box on; 
plot(Phi_s,Ys)
plot(Phi_s,Ys_noisy,'-r')
plot(Phi_s,Ys_averaged,'-b')




figure(436); clf; hold on; box on;
plot(Phi_t,Yt_store,LineWidth=2,Color=gray,HandleVisibility='off')
plot(Phi_t,Yt_true,'-r',LineWidth=2)%,DisplayName="$y_T$")
plot(Phi_t,Yt_averaged,'-k',LineWidth=2)%,DisplayName="$y_b$: averaging $y_b$")
plot(Phi_t,Yt_from_averaged_Ys,'-b',LineWidth=2)%,DisplayName="$y_b$: averaging $y_f$")
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$y_b$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])
%legend('True','averaged Yt','Yt from averaged Ys')
%legend(Interpreter="latex",FontSize=14)

%plot the inset
xstart=0.2; ystart=0.7; xend=0.4; yend=0.9;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on; hold on;
plot(Phi_t,Yt_store,LineWidth=2,Color=gray,HandleVisibility='off')
plot(Phi_t,Yt_true,'-r',LineWidth=3)
plot(Phi_t,Yt_averaged,'-k',LineWidth=2)
plot(Phi_t,Yt_from_averaged_Ys,'-b',LineWidth=2)
xlim([-1.5,1.5])
ylim([0.3,0.45])


% 
% legendstring='';
% %legendstring="$\alpha="+num2str(amplitudes')+"$";
% legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')



