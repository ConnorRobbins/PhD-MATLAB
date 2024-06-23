k=linspace(0,3,1000);
A=0.2;
mus=[0.55,0,-0.55];


figure(17); clf; hold on;

for i=1:numel(mus)
    mu=mus(i);
    plot(k,abs(A*(mu+(k.^2)/3)))
end

yline(A*0.8+0.1,'--k')
yline(A*0.55-0.02,'--k')
ylim([0,0.4])

legend('$\mu>0$: Supercritical',' $\mu=0$: Critical','$\mu<0$: Subcritical','Location','northwest','Interpreter','latex','FontSize',12)
xticks([])
yticks([])
xlab=xlabel('$k$','Interpreter','latex','FontSize',16);
ylab=ylabel('$c_i$','Interpreter','latex','FontSize',17,'Rotation',0);
xlab.Position(1:2)=[1.5,-0.014];
ylab.Position(1:2)=[-0.2,0.2];

%label x intercept
%text(sqrt(-3*mus(end))+0.1,0.01,'$\leftarrow k=\sqrt{-3\mu}$','Interpreter','latex','FontSize',12)
arrowposx=[0.56,0.46];
arrowposy=[0.15,0.11];
annotation('textarrow',arrowposx,arrowposy,'String','$k=\sqrt{-3\mu}$','Interpreter','latex','FontSize',12);



%label y intercept
arrowposx=[0.181,0.13];
arrowposy=[0.45,0.334];
annotation('textarrow',arrowposx,arrowposy,'String','$c_1(0)=A|\mu|$','Interpreter','latex','FontSize',12);




%text(0+0.1,abs(A*mus(1))+0.01,'$\leftarrow k=\sqrt{-3\mu}$','Interpreter','latex','FontSize',12)



% annotation('textarrow',[0+0.1,abs(A*mus(1))+0.01])
% annotation('textarrow',0+0.1,abs(A*mus(1))+0.01,'String','$\leftarrow k=\sqrt{-3\mu}$');%,'Interpreter','latex','FontSize',12)

