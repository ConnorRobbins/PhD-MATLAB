function [frequencies,yshift] = FUN_FourierTransform_topography(Yt,L,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


dx=2*L/(N-1);
fs=1/dx;
% fshift=(-N/2:N/2-1)*(fs/N);
% yshift=abs(fftshift(fft(Yt)))/N;
% frequencies=2*pi*fshift;


fshift=(0:(N-1))*fs/N;
if mod(N,2)==1 %odd
    fshift((N+1)/2+1:end) = fshift((N+1)/2+1:end)-fs;
else % even
    fshift(N/2+1:end) = fshift(N/2+1:end)-fs;
end


frequencies=2*pi*fftshift(fshift);
yshift=abs(fftshift(fft(Yt)))/N;





% figure(2); clf; hold on;
% plot(frequencies,yshift)
% xlabel('frequency')
% 
% 
% figure(3); clf; hold on;
% stem(frequencies,yshift)
% xlabel('frequency')




end