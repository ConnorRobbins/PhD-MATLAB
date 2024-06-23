function [pos_frequencies,pos_yshift] = fun_FT_post_clean(frequencies,yshift)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    N=numel(frequencies);
    mid=ceil((N+1)/2);
    pos_frequencies=frequencies(mid:end);
    pos_yshift=yshift(mid:end);
    pos_yshift(2:end)=2*pos_yshift(2:end);
end