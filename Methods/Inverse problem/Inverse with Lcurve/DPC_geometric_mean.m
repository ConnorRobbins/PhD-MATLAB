function [geo_mean_indices,geo_mean] = DPC_geometric_mean(singularvals,Fourier_coeffs,q)
%UNTITLED2 Summary of this function goes here
% see the discrerte picard condition for discrete ill-posed problems -
% hansen 1990
geo_mean_indices=q+1:1:numel(singularvals)-q;
geo_mean=zeros(numel(singularvals),1)*NaN;

for i=q+1:1:numel(singularvals)-q
    geo_mean(i)=(((prod(Fourier_coeffs(i-q:i+q)))^(1/(2*q -1))))/singularvals(i);
end




end
