function [first_max,last_max] = FUN_middle_wave_extraction(Y,Phi,surface_support)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

truth_indices=abs(Phi)<surface_support; %extract indices of positions within support 
tf=islocalmax(Y); %extract indices of local maxima
tf=reshape(tf,size(truth_indices)); %make sure the two truth vectors are same shape
truth_indices=truth_indices.*tf; %obtain truth values for both in support and local maxima

%find and cut off first and last local maxima to help remove edge effects
first_max=find(truth_indices,1,"first");
last_max=find(truth_indices,1,"last");
truth_indices(first_max)=0;
truth_indices(last_max)=0;


%find and return the first and last maxima now the edge cases have been
%removed
first_max=find(truth_indices,1,"first");
last_max=find(truth_indices,1,"last");



end