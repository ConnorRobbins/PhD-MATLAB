function [straight_y,m] = FUN_simple_straight_line_fit(x,y)
%SIMPLE_STRAIGHT_LINE_FIT Fits a straight line on first and last points
%   Detailed explanation goes here
x1=x(1);
x2=x(end);
y1=y(1);
y2=y(end);
m=(y2-y1)/(x2-x1);

straight_y=y1+ (x-x1)*m;
end

