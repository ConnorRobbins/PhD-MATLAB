clear

N=871;
amplitude=0.01;

Froude=0.8;

support_gradient=1;
length_support_multiplier=2.0;
wave_ks=[0.5,1,1.5,2];
wave_ks=[0.5,0.75,1,1.25];


[N_s,N_t]=deal(N);
wave_k_n=numel(wave_ks);



P=zeros(1,N_s);     %Pressure on surfacgrad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"


grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank

boo_auto_pick_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 
rankSVD=200; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem


Ys_store=zeros(N,wave_k_n);
Yt_store=zeros(N,wave_k_n);
Phi_t_store=zeros(N,wave_k_n);



for i =1:wave_k_n
    wave_k=wave_ks(i);
    surface_support=11*pi/(2*wave_k);
    L=length_support_multiplier*surface_support;
    Ls(i)=L;
    Phi_s=linspace(-L,L,N);
    tanh_support=0.5*(tanh(support_gradient*(surface_support-Phi_s))+tanh(support_gradient*(surface_support+Phi_s)));
    Ys=1+amplitude*cos(wave_k*Phi_s).*tanh_support;


    [Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol);
    Ys_store(:,i)=Ys_after;
    Yt_store(:,i)=Yt;
    Phi_t_store(:,1)=Phi_t;

    [frequencies,yshift] = fun_FT_topography(Yt,L,N);
    [pos_frequencies,pos_yshift]=fun_FT_post_clean(frequencies,yshift);


    figure(2*i);  clf; hold on;
    plot(Phi_s,Yt)
    xline(-surface_support,'--r')
    xline(surface_support,'--r')
    title(num2str(wave_k))



    figure(2*i +1); clf;
    %stem(frequencies,yshift)
    %plot(frequencies,log10(yshift))
    %plot(frequencies,yshift)
    plot(pos_frequencies,pos_yshift)
    title(num2str(wave_k))
end

% 
% figure(22); clf;
% plot(Phi_s,tanh_support)
% 
% 
% figure(23); clf; hold on;
% plot(Phi_s,Ys)
% plot(Phi_t,Yt)