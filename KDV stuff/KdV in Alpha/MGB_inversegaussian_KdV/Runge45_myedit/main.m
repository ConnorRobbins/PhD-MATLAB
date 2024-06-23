%close all
clear all
clc

format long

global alpha 

alpha = 1.5;

u01    = 1.261621093750000; % u
u02    = 1.261645507812500; % d


u0 = 0.5*(u01+u02)

KK = 20000;

tol = 1.d-7;
dx  = 1.d-3;
h   = dx;

x = 0;

y(1) = u0;
y(2) = 0;

u(1,1) = y(1);
u(2,1) = y(2);
xx(1)  = x;

k  = 1;


while k<KK
    k = k+1;

    
    [x,y]  = runge4(y,x,h,@eqsystem2);
    
    u(1,k) = y(1);
    u(2,k) = y(2);
    xx(k)  = x;



end

subplot(1,2,1)
uu = u(1,:);
ud = u(2,:);
plot(uu,ud,'b','linewidth',2.5)
grid on
xlabel('u')
ylabel('du/dx')


subplot(1,2,2)
plot(xx,uu,'k','linewidth',1.5)
hold on
plot(-xx,uu,'k','linewidth',1.5)
xlabel('x')
ylabel('u')
