clear
saveplots=0;


%% Physical parameters 

s=load('supercritical_gaussian_bump_a2E-1_b1_F1p5.mat','N','L','b_width','Froude','amplitude','Ys_Newton','Yt');
N=s.N;
[N_t,N_s]=deal(N);
L=s.L;
b_width=s.b_width; 
Froude=s.Froude;
amplitude=s.amplitude; 
Ys=s.Ys_Newton;
Yt_true=reshape(s.Yt,N,1);


measurement_noise_std=1E-2;
Ys_measurement=Ys; noise=measurement_noise_std*2*(randn(N_s-2,1))'; Ys_measurement(2:end-1)=Ys_measurement(2:end-1)+noise;

perturbation_N=100;  %how many different noise stds#
perturbation_amp_std=1E-2;


truncation_ranks=[40];



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





for i=1:loopnumber
        disp("Starting step "+num2str(i)+"/"+steps)
        Ys_noisy=Ys_measurement; noise=perturbation_amp_std*2*(randn(N_s-2,1))'; Ys_noisy(2:end-1)=Ys_noisy(2:end-1)+noise;
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



%%


figure(1); clf; hold on; box on; 
plot(Phi_s,Ys)
plot(Phi_s,Ys_measurement,'-r')
plot(Phi_s,Ys_noisy,'-b')




figure(436); clf; hold on; box on;
plot(Phi_t,Yt_store,LineWidth=2,Color=[0.6,0.6,0.6])
plot(Phi_t,Yt_true,'-r',LineWidth=2)
plot(Phi_t,Yt_averaged,'--k',LineWidth=2)
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$y_b$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])


legendstring='';
%legendstring="$\alpha="+num2str(amplitudes')+"$";
legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')



