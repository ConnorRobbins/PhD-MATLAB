function [first_turning,last_turning] = FUN_middle_wave_extraction_n_turningpoint_cutoff(Y,Phi,surface_support,n_turningpoint_cutoff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

truth_indices=abs(Phi)<surface_support; %extract indices of positions within support 
tfmax=islocalmax(Y); %extract indices of local maxima
tfmin=islocalmin(Y);
tf=tfmax|tfmin;
tf=reshape(tf,size(truth_indices)); %make sure the two truth vectors are same shape
truth_indices=truth_indices.*tf; %obtain truth values for both in support and local maxima

%find and cut off first and last local turningpoint to help remove edge effects
for i=1:n_turningpoint_cutoff
    first_turning=find(truth_indices,1,"first");
    last_turning=find(truth_indices,1,"last");
    truth_indices(first_turning)=0;
    truth_indices(last_turning)=0;
end


%find and return the first and last maxima now the edge cases have been
%removed
first_turning=find(truth_indices,1,"first");
last_turning=find(truth_indices,1,"last");



end