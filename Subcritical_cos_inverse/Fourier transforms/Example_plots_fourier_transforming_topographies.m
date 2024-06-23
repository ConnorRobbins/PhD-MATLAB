N=851;
amplitude=0.01;

Froude=0.8;

support_gradient=1;
length_support_multiplier=2.5;
%wave_ks=[0.5,1,1.5,2];
wave_ks=linspace(0.5,2,30);


plot_every_k=0;


[N_s,N_t]=deal(N);
wave_k_n=numel(wave_ks);



P=zeros(1,N_s);     %Pressure on surfacgrad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"


grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank

boo_auto_pick_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 
rankSVD=300; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem

figure(1); clf; 
peak_freqs=zeros(1,wave_k_n);

for i =1:wave_k_n
    wave_k=wave_ks(i);
    surface_support=9*pi/(2*wave_k);
    L=length_support_multiplier*surface_support;
    Ls(i)=L;
    Phi_s=linspace(-L,L,N);
    tanh_support=0.5*(tanh(support_gradient*(surface_support-Phi_s))+tanh(support_gradient*(surface_support+Phi_s)));
    Ys=1+amplitude*cos(wave_k*Phi_s).*tanh_support;


    [Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol);
    Ys_store(:,i)=Ys_after;
    Yt_store(:,i)=Yt;

    [frequencies,yshift] = fun_FT_topography(Yt,L,N);
    [pos_frequencies,pos_yshift]=fun_FT_post_clean(frequencies,yshift);


    [peak_amp,peak_freq_i]=max(pos_yshift);
    peak_freqs(i)=pos_frequencies(peak_freq_i);
    
    
    
    if plot_every_k==1
        subplot(4,2,2*i -1)
        plot(Phi_s,Yt)
        title(strcat('k=',num2str(wave_k)))
        xlabel('Phi')
        ylabel('Yt')
         
        subplot(4,2,2*i)
        %stem(frequencies,yshift)
        %plot(frequencies,log10(yshift))
        %plot(frequencies,yshift)
        plot(pos_frequencies,pos_yshift)
        title(strcat('k=',num2str(wave_k)))
        xlabel('freq')
        xlim([0,10])
    end
end

figure(2); clf; hold on; 
plot(wave_ks,peak_freqs,'-xr')
plot(wave_ks,wave_ks,'--b')
legend('Data','y=x')