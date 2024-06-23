clear

%filename='Apoint01_N1624_k_half_to_2.mat'; xlims=[0.5,1.5]; ylims=[-4.3,-0.9];
%filename='Apoint05_N1624_k_half_to_2.mat'; xlims=[0.5,1.5]; ylims=[-3,-0.5];
%filename='Apoint1_N1624_k_half_to_1p5.mat'; xlims=[0.5,1.5]; ylims=[-2.5,-0.5];
filename='Apoint2_F1p05_p3_to_p9.mat'; xlims=[0.3,0.9]; ylims=[-2.,-1];


load(filename)

for i=1:ki
    kdv_freq=kdv_freqs{i};
    kdv_y=kdv_yshift{i};
    nonlinear_freq=nonlinear_freqs{i};
    nonlinear_y=nonlinear_yshift{i};

    k=ks(i);
    
    % Extract c1 mode
    [~,kdv_c1_index]=min(abs(kdv_freq-k));
    [~,nonlinear_c1_index]=min(abs(nonlinear_freq-k));
    kdv_c1_power(i)=kdv_y(kdv_c1_index);
    nonlinear_c1_power(i)=nonlinear_y(nonlinear_c1_index);


     % Extract c2 mode
    [~,kdv_c2_index]=min(abs(kdv_freq-2*k));
    [~,nonlinear_c2_index]=min(abs(nonlinear_freq-2*k));
    kdv_c2_power(i)=kdv_y(kdv_c2_index);
    nonlinear_c2_power(i)=nonlinear_y(nonlinear_c2_index);


    [~,nonlinear_c3_index]=min(abs(nonlinear_freq-3*k));
    nonlinear_c3_power(i)=nonlinear_y(nonlinear_c3_index);

end

% 
% figure(61); clf; hold on;
% plot(ks,nonlinear_c1_power,'-xk')
% plot(ks,nonlinear_c2_power,'--xk')
% plot(ks,kdv_c1_power,'-xr')
% plot(ks,kdv_c2_power,'--xr')


% plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b')
% plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb')

%legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')
%plot(ks,nonlinear_c3_power,':m')





figure(62); clf; hold on;
plot(ks,log10(nonlinear_c1_power),'-k',LineWidth=1,DisplayName="Nonlinear $c_1$")
plot(ks,log10(abs(amplitude*(2*(Froude-1) +(ks.^2)/3))),':k',LineWidth=1.5,DisplayName='fKdV $c_1$')
plot(ks,log10(nonlinear_c2_power),'-r',LineWidth=1,DisplayName="Nonlinear $c_2$")
plot(ks,log10(ones(size(ks))*3*(amplitude^2)/4),'--r',LineWidth=1,DisplayName='fKdV $c_2$')

legend(Interpreter='latex',FontSize=14,Location='northwest')
xlab=xlabel('$\nu$','Interpreter','latex','FontSize',16);
ylab=ylabel('$\log_{10}(c_i)$','Interpreter','latex','FontSize',17,'Rotation',90);
xlim(xlims)
ylim(ylims)
box on


%%%% Intersection points of nonlinear
[~,intersection_index]=min(abs(nonlinear_c1_power-nonlinear_c2_power));
intersection_nu=ks(intersection_index)
intersection_energy1=nonlinear_c1_power(intersection_index)
intersection_energy2=nonlinear_c2_power(intersection_index)