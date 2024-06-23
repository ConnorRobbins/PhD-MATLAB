
close all
clear all
clc

format long

global alpha n h

n = 2;

alpha = 1.5;

u01    = 1.261621093750000; % u
u02    = 1.261645507812500; % d

u0 = 0.5*(u01+u02)

KK = 30000;

tol = 1.d-7;
dx  = 1.d-3;
h   = dx;

x = 0;

y(1) = u0;
y(2) = 0;

u(1,1) = y(1);
u(2,1) = y(2);
xx(1)  = x;
qq(1)  = 99.0; sqrt(6/abs(u(1,1))) - x;

dq = 99.0;
k  = 1;
im = 0;
%while abs(y(1))<blowup
%while dq>tol
while k<KK
    k = k+1;

    [dydx] = eqsystem(x,y);
    [x,y]  = runge4(y,dydx,x);
    
    u(1,k) = y(1);
    u(2,k) = y(2);
    xx(k)  = x;

    if y(1)<0
        im = im+1;
        xx2(im) = x;
        uu2(im) = y(1);
        pp(im)  = -sqrt(-2*y(1).^3/3);
        qq(im)  = x + sqrt(6/abs(y(1)));
        rr(im)  = -(3/2)*y(2)^2/y(1)^3;
    end
    
%    disp([x, y(1), rr(k), x+sqrt(6/abs(y(1)))])
    
    if im>3
%        dq = abs( (qq(im)-qq(im-1))/qq(im-1) )
%        dq = 0.5*abs(qq(im)-qq(im-2))/h;
    end
    
%    if k==100 stop; end

end

subplot(1,2,1)
uu = u(1,:);
ud = u(2,:);
plot(uu,ud,'b','linewidth',2.5)
hold on
plot(uu2,pp,':r','linewidth',1.5)
plot(uu2,-pp,':r','linewidth',1.5)
xlim([-0.5 1.5])
ylim([-1.5 0.25])
grid on
xlabel('u')
ylabel('du/dx')


subplot(1,2,2)
plot(xx,uu,'k','linewidth',1.5)
hold on
plot(-xx,uu,'k','linewidth',1.5)
xlabel('x')
ylabel('u')

%save('output.mat','xx','uu','ud')