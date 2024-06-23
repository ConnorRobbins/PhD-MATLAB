N=1624;
L=30;
amplitude=0.1;

x=linspace(-L,L,N);


ks=0.1:0.1:20;
ks=0.3;

smoothenough=true;

for i=1:numel(ks)
    if smoothenough
        k=ks(i);
        Y=amplitude*cos(k*x);
        [message1,smoothenough]=FUN_angle_checks(x,Y);
    end
end
k
message1