clear


file1='Apoint01_N1624_k_half_to_2.mat';

%file2="point9_point95";
% file2="60_ks_N1624_k_1_to_2";
%file2="datacleaned";
file2='Apoint05_N1624_k_half_to_2.mat';




%% Do it for first data set
load(file1)




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


figure(51); clf; hold on;
title(strcat("amp=",num2str(amplitude)))
plot(ks,nonlinear_c1_power,'-b')
plot(ks,nonlinear_c2_power,'--b')
plot(ks,kdv_c1_power,'-r')
plot(ks,kdv_c2_power,'--r')
% plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b')
% plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb')
legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')



figure(52); clf; hold on;
title(strcat("amp=",num2str(amplitude)))
[straight_c1,m1]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c1_power));
[straight_c2,m2]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c2_power));
plot(ks,log10(nonlinear_c1_power),'-b',DisplayName=strcat("log(c1)"))
plot(ks,straight_c1,'--m',DisplayName=strcat("ref line gradient=",num2str(m1)))
P1=polyfit(ks,log10(nonlinear_c1_power),1);
plot(ks,polyval(P1,ks),DisplayName=strcat("best fit, gradient=",num2str(P1(1))))
plot(ks,log10(nonlinear_c2_power),'--b',DisplayName=strcat("log(c2)"))
plot(ks,straight_c2,'--',Color=[0.8500 0.3250 0.0980],DisplayName=strcat("ref line gradient=",num2str(m2)))
P2=polyfit(ks,log10(nonlinear_c2_power),1);
plot(ks,polyval(P2,ks),DisplayName=strcat("best fit, gradient=",num2str(P2(1))))
legend




figure(55); clf; hold on;
title("both amps log scale")
plot(ks,log10(nonlinear_c1_power),'-b',DisplayName=strcat("log(c1): amp=",num2str(amplitude)))
plot(ks,log10(nonlinear_c2_power),'--b',DisplayName=strcat("log(c2): amp=",num2str(amplitude)))
%legend('log(c1)','log(c2)')




%% Do it for second data set
clearvars -except file2
load(file2)




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




figure(53); clf; hold on;
title(strcat("amp=",num2str(amplitude)))
plot(ks,nonlinear_c1_power,'-k')
plot(ks,nonlinear_c2_power,'--k')
plot(ks,kdv_c1_power,'-r')
plot(ks,kdv_c2_power,'--r')
% plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b')
% plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb')
legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')


figure(54); clf; hold on;
title(strcat("amp=",num2str(amplitude)))
[straight_c1,m1]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c1_power));
[straight_c2,m2]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c2_power));
plot(ks,log10(nonlinear_c1_power),'-b',DisplayName=strcat("log(c1)"))
plot(ks,straight_c1,'--m',DisplayName=strcat("ref line gradient=",num2str(m1)))
P1=polyfit(ks,log10(nonlinear_c1_power),1);
plot(ks,polyval(P1,ks),DisplayName=strcat("best fit, gradient=",num2str(P1(1))))
plot(ks,log10(nonlinear_c2_power),'--b',DisplayName=strcat("log(c2)"))
plot(ks,straight_c2,'--',Color=[0.8500 0.3250 0.0980],DisplayName=strcat("ref line gradient=",num2str(m2)))
P2=polyfit(ks,log10(nonlinear_c2_power),1);
plot(ks,polyval(P2,ks),DisplayName=strcat("best fit, gradient=",num2str(P2(1))))

legend




figure(55); 
plot(ks,log10(nonlinear_c1_power),'-k',DisplayName=strcat("log(c1): amp=",num2str(amplitude)))
plot(ks,log10(nonlinear_c2_power),'--k',DisplayName=strcat("log(c2): amp=",num2str(amplitude)))
legend

