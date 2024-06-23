clear
%Calculations at 7 points, the intersections and x intercept and
%equidistant around those


N=1624; 
amplitude=0.05;
Froude=0.95;

mu=Froude-1;
kl=sqrt(-(((9*abs(amplitude)/4) +6*mu)));
km=sqrt(-6*mu);
kr=sqrt((((9*abs(amplitude)/4) -6*mu)));
klm=(km-kl)/7;
kmr=(kr-km)/7;
wave_ks=[kl-6*klm,kl-5*klm,kl-4*klm,kl-3*klm,kl-2*klm,kl-klm,kl,kl+klm,kl+2*klm,kl+3*klm,kl+4*klm,kl+5*klm,kl+6*klm,km,km+kmr,km+2*kmr,km+3*kmr,km+4*kmr,km+5*kmr,km+6*kmr,kr,kr+kmr,kr+2*kmr,kr+3*kmr,kr+4*kmr,kr+5*kmr,kr+6*kmr];



boo_auto_pick_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 
rankSVD=900; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem

grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank

n_turningpoint_cutoff_kdv=2;
n_turningpoint_cutoff_nonlinear=2;


support_gradient=1.2;
length_support_multiplier=2.;



[N_s,N_t]=deal(N);
P=zeros(1,N_s);     %Pressure on surfacgrad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"



for wave_ki=1:numel(wave_ks)




wave_k=wave_ks(wave_ki);


surface_support=13*pi/(2*wave_k);
L=length_support_multiplier*surface_support;
Phi_s=linspace(-L,L,N);
tanh_support=0.5*(tanh(support_gradient*(surface_support-Phi_s))+tanh(support_gradient*(surface_support+Phi_s)));
Ys=1+amplitude*cos(wave_k*Phi_s).*tanh_support;



%% Calculate nonlinear topography, extract centre and Fourier transform

[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt_nonlinear,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol);


[frequencies_nonlinear,yshift_nonlinear] = FUN_FourierTransform_topography(Yt_nonlinear,L,N);
[pos_frequencies_nonlinear,pos_yshift_nonlinear]=FUN_FourierTransform_post_clean(frequencies_nonlinear,yshift_nonlinear);


%middle extraction
[left_index_nonlinear,right_index_nonlinear]=FUN_middle_wave_extraction_n_turningpoint_cutoff(Yt_nonlinear,Phi_t,surface_support,n_turningpoint_cutoff_nonlinear);
middle_L=(Phi_t(right_index_nonlinear)-Phi_t(left_index_nonlinear))/2;
middle_N=right_index_nonlinear-left_index_nonlinear+1;
[frequencies_nonlinear_middle_only,yshift_nonlinear_middle_only] = FUN_FourierTransform_topography(Yt_nonlinear(left_index_nonlinear:right_index_nonlinear),middle_L,middle_N);
[pos_frequencies_nonlinear_middle_only,pos_yshift_nonlinear_middle_only] = FUN_FourierTransform_post_clean(frequencies_nonlinear_middle_only,yshift_nonlinear_middle_only);


[~,maxf]=max(pos_yshift_nonlinear_middle_only);
max_amp_freq_nonlinear_middle_only=pos_frequencies_nonlinear_middle_only(maxf);





%% Calculate KdV topography





[Yt_kdv] = FUN_KdV_InverseForcing_UniformFarField(Ys-1,Phi_s,Froude);

[frequencies_kdv,yshift_kdv] = FUN_FourierTransform_topography(Yt_kdv,L,N);
[pos_frequencies_kdv,pos_yshift_kdv]=FUN_FourierTransform_post_clean(frequencies_kdv,yshift_kdv);

%middle extraction
[left_index_kdv,right_index_kdv]=FUN_middle_wave_extraction_n_turningpoint_cutoff(Yt_kdv,Phi_t,surface_support,n_turningpoint_cutoff_kdv);
middle_L=(Phi_t(right_index_kdv)-Phi_t(left_index_kdv))/2;
middle_N=right_index_kdv-left_index_kdv+1;
[frequencies_kdv_middle_only,yshift_kdv_middle_only] = FUN_FourierTransform_topography(Yt_kdv(left_index_kdv:right_index_kdv),middle_L,middle_N);
[pos_frequencies_kdv_middle_only,pos_yshift_kdv_middle_only] = FUN_FourierTransform_post_clean(frequencies_kdv_middle_only,yshift_kdv_middle_only);



[~,maxf]=max(pos_yshift_kdv_middle_only);
max_amp_freq_kdv_middle_only=pos_frequencies_kdv_middle_only(maxf);


save(strcat(num2str(numel(wave_ks)),'K',num2str(wave_ki),'.mat'))

end






%% Plotting

figure(2); clf; hold on;
plot(Phi_s,Ys,'-b')
plot(Phi_t,Yt_nonlinear,'-r')
plot(Phi_t,Yt_kdv,'-k')
legend('Surface','Nonlinear','KdV',Location='east')
xlim([-L,L])



figure(3); clf; hold on;
plot(Phi_t,Yt_nonlinear,'-r')
plot(Phi_t,Yt_kdv,'-k')
legend('Nonlinear','KdV')
xlim([-L,L])
ylabel('Yt')



figure(4); clf; 
%
subplot(3,1,1)
stem(pos_frequencies_nonlinear,pos_yshift_nonlinear,'-r','MarkerSize',4)
title('nonlinear')
xlim([0,5*wave_k])
%
subplot(3,1,2)
stem(pos_frequencies_kdv,pos_yshift_kdv,'-k','MarkerSize',4)
title('kdv')
xlim([0,5*wave_k])
%
subplot(3,1,3); hold on;
stem(pos_frequencies_nonlinear,pos_yshift_nonlinear,'-r','MarkerSize',3)
stem(pos_frequencies_kdv,pos_yshift_kdv,'-k','MarkerSize',3)
title('both')
xlim([0,5*wave_k])






figure(5); clf; hold on;
plot(Phi_t,Yt_nonlinear,'--r')
plot(Phi_t,Yt_kdv,'--k')
plot(Phi_t(left_index_nonlinear:right_index_nonlinear),Yt_nonlinear(left_index_nonlinear:right_index_nonlinear),'-r','LineWidth',2)
plot(Phi_t(left_index_kdv:right_index_kdv),Yt_kdv(left_index_kdv:right_index_kdv),'-k','LineWidth',2)
xlim([-L,L])
legend('Nonlinear','KdV','Nonlinear middle','KdV middle')
ylabel('Yt')






figure(6); clf; 
%
subplot(3,1,1)
stem(pos_frequencies_nonlinear_middle_only,pos_yshift_nonlinear_middle_only,'-r','MarkerSize',4)
title('nonlinear')
xlim([0,5*wave_k])
%
subplot(3,1,2)
stem(pos_frequencies_kdv_middle_only,pos_yshift_kdv_middle_only,'-k','MarkerSize',4)
title('kdv')
xlim([0,5*wave_k])
%
subplot(3,1,3); hold on;
stem(pos_frequencies_nonlinear_middle_only,pos_yshift_nonlinear_middle_only,'-r','MarkerSize',3)
stem(pos_frequencies_kdv_middle_only,pos_yshift_kdv_middle_only,'-k','MarkerSize',3)
title('both')
xlim([0,5*wave_k])



figure(7); clf; hold on;
plot(log10(Theta_b_norm))
xline(rankSVD,'--r')
legend('log10 norm','chosen rank',Location='southeast')