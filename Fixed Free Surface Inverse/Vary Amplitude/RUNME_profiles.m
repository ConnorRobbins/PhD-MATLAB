clear
saveplots=0;

%% Physical parameters 


N=641; L=20; b_width=0.3; Froude=0.9; amplitudes=-[0.05,0.1,0.125,0.15,0.2]; %truncation_ranks=[];


truncation_ranks=[90];


[N_t,N_s]=deal(N);
Phi_s=linspace(-L,L,N_s); Phi_t=linspace(-L,L,N_t);
P=zeros(1,N_s);






loopnumber=numel(amplitudes);
steps=num2str(loopnumber);

Yt_store=zeros(N_t,loopnumber);
P_store=zeros(N_s,loopnumber);




if numel(truncation_ranks)==1
    truncation_ranks=ones(1,loopnumber)*truncation_ranks;
end




for i=1:loopnumber
        disp("Starting step "+num2str(i)+"/"+steps)
        amplitude=amplitudes(i);
        Ys=1+amplitude*exp(-(b_width*Phi_s).^2);
        Theta_approx_for_pressure=atan(-2*amplitude*(b_width^2)*Phi_s.*exp(-(b_width*Phi_s).^2))';
        truncation_rank=truncation_ranks(i);



        %% %% Perform the topography calculations
        %Form the truncated pseudoinverse
        [truncated_psuedoinv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_rank);
        % form the RHS of the matrix eqn and output Theta_f
        [D,Theta] = FUN_make_SVD_RHS(Ys,P,Froude,L,N_s,N_t);
        %Solve the truncated system
        Theta_b=truncated_psuedoinv*D;
        [Yt_store(:,i),~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b,Theta);


        %% %% Perform the pressure calculations
        [~,~,~,~,~,~,~,inverse_Pressure] = FUN_Inverse_Pressure_NO_TOPOGRAPHY(L,N_s,N_t,Froude,Theta_approx_for_pressure);
        P_store(:,i)=inverse_Pressure;

end






%%





legendstring="$\alpha="+num2str(amplitudes')+"$";
% yaxismax=max(max(max(P_store)),max(max(Yt_store)));
% yaxismin=min(min(min(P_store)),min(min(Yt_store)));


figure(436); clf; hold on; box on;
plot(Phi_t,Yt_store,LineWidth=2)
legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$y_b$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])



figure(437); clf; hold on; box on;
plot(Phi_s,P_store,LineWidth=2)
legend(legendstring,Interpreter="latex",FontSize=14,Location='southeast')
xlabel('$\phi$',Interpreter='latex',FontSize=18)
ylabel('$P$',Interpreter='latex',FontSize=18)
%ylim([yaxismin,yaxismax])


