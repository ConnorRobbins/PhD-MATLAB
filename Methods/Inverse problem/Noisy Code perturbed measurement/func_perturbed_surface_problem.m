function [Yt_std,Ys_std,Ys_noisy_store,Yt_noisy_store,Theta_b_store,truncated_inv,Yt_noisy_average]=func_perturbed_surface_problem(L,Froude,N_s,N_t,P,Ys,Phi_s,Phi_t,rankSVD_noisy,U,singularvals,V,unperturbed_endpoints,Noise_range,Random_error_Number,random_noise_std)
%[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol)


%% Solve the perturbed problems

% Initialise storage matrices 

Ys_noisy_store=NaN*ones(N_s,Random_error_Number);   %stores perturbed surfaces
Yt_noisy_store=NaN*ones(N_t,Random_error_Number);   %stores all outputs for topography Yt
Theta_b_store=NaN*ones(N_t,Random_error_Number); %stores all outputs for topography theta_b





%Identify the indices of all surface points to which noise will be added
[closestval,neg_range_index]=min(abs(Phi_s+Noise_range));
[closestval,pos_range_index]=min(abs(Phi_s-Noise_range));

if ~unperturbed_endpoints==0
    if neg_range_index<=unperturbed_endpoints
        neg_range_index=unperturbed_endpoints+1; %shift to exclude end points if needed
    end
    if pos_range_index>=N_s-unperturbed_endpoints
        pos_range_index=N_s-unperturbed_endpoints; %shift to exclude end points if needed
    end
end
noise_range_entries=pos_range_index-neg_range_index+1;


u=U';
truncated_inv=zeros(N_t,N_s);
% Now loop over the ranks to build the TSVD inverse
for k=1:rankSVD_noisy
    truncated=(1/singularvals(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
end




%Loop to repeatedly perturb surface and solve inverse problem by same
%method as previously with no noise
for rn=1:Random_error_Number
    %Generate noise and add to surface 
    random_noise=random_noise_std*randn(1,noise_range_entries);
    Ys_noisy=Ys;
    Ys_noisy(neg_range_index:pos_range_index)=Ys(neg_range_index:pos_range_index)+random_noise;
    Ys_noisy_store(:,rn)=Ys_noisy; %store the noisy surface as a column
    
    % Call svd function to compute theta_b_noisy at all truncation ranks,
    % slightly different function to previously used as SVD already
    % calculated so no need to calculate again
    [Phi_t,Phi_s,Theta_noisy,Theta_b_noisy,D_noisy] = svd_Step_fun_given_SVD(L,Froude,N_s,N_t,P,Ys_noisy,truncated_inv);
    
    
  
    % Calculate Yt_noisy from theta_b at selected rank
    [Yt_noisy,Ys_after_noisy] = Y_after_SVD_fun(L,N_t,N_s,Theta_b_noisy,Theta_noisy);

    %store variables at chosen rank
    Theta_b_store(:,rn)=Theta_b_noisy;
    Yt_noisy_store(:,rn)=Yt_noisy;

end



%Average the  topographies at each meshpoint
Yt_noisy_average=mean(Yt_noisy_store,2,'omitnan');



% Calculate the standard devs
Yt_std=std(Yt_noisy_store,[],2,'omitnan');
Ys_std=std(Ys_noisy_store,[],2,'omitnan');

