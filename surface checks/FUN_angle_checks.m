function [outputArg1,outputArg2] = FUN_angle_checks(x,Y)
%ANGLE_CHECKS Summary of this function goes here
%   Detailed explanation goes here
theta=[atand((Y(2:end)-Y(1:end-1))./(x(2:end)-x(1:end-1)))];
theta_changes=abs(theta(2:end)-theta(1:end-1));
max_theta_change=max(theta_changes);
if max_theta_change<120
    outputArg1=strcat("Surface is okay, the greatest angle is ",num2str(max_theta_change)," degrees");
    outputArg2=true;
else
    outputArg1=strcat("Surface is too sharp, the greatest angle is ",num2str(max_theta_change)," degrees");
    outputArg2=false;
end

