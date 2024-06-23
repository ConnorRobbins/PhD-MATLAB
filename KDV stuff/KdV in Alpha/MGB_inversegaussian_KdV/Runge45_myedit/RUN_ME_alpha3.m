%close all
clear all
clc

format long

global alpha 

alpha = 3;

u01    =  1.705355362892151; % u
u02    =  1.705355424880981; % d


u0 = 0.5*(u01+u02)

KK = 200000;

tol = 1.d-7;
dx  = 5.d-4;
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

algebraic_x=linspace(-1,0,100);
algebraic_y=sqrt((-2/3)*algebraic_x.^3);


figure(1); clf; hold on; box on;
uu = u(1,:);
ud = u(2,:);
plot(uu,ud,'b','linewidth',2.5)
plot(algebraic_x,algebraic_y,'--r')
plot(algebraic_x,-algebraic_y,'--r')
grid on
xlabel('u')
ylabel('du/dx')
ylim([-1.5,0.5])
xlim([-1,1.5])


figure(2); clf; hold on; box on;
plot(xx,uu,'k','linewidth',1.5)
hold on
plot(-xx,uu,'k','linewidth',1.5)
xlabel('x')
ylabel('u')
