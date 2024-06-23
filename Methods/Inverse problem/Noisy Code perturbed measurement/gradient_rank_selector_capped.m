function [maxrank] = gradient_rank_selector_capped(Theta_b_norm,grad_tol,log_normmax)
%GRADIENT_RANK_SELECTOR Summary of this function goes here
%   Detailed explanation goes here

logvals=log10(Theta_b_norm);
logvals_grad=logvals(2:end)-logvals(1:end-1);
% 
norm_val=8; %initialise for a cap on norm
%log_normmax=2; %set a cap on log(normval)

%Find longest region where gradient is small enough
grad_big=logvals_grad>grad_tol; % 1 where gradient larger than tolerance
grad_cum=cumsum(grad_big); % sum truth values, increases when gradient too large

while norm_val>log_normmax
    flat_indices=find(grad_cum==mode(grad_cum)); % find indices of largest non-increasing region
    % Choose maxrank to be ~80% along flat region
    maxrank=flat_indices(floor(0.8*length(flat_indices)));
    norm_val=logvals(maxrank);
    
    if norm_val>log_normmax
        grad_cum(flat_indices)=NaN;
    end
end

end

