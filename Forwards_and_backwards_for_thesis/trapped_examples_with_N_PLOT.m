clear
load("trapped_examples_with_N_DATA_c4.mat")

figure(1); clf; hold on; box on;
for i=[1,3,5]
    plot(Phi_store{i},Yt_store{i},DisplayName="$N=$ "+num2str(N_store{i}),LineWidth=1.5)
    xlabel('$\phi$',Interpreter='latex',FontSize=18)
    ylabel('$y_b$',Interpreter='latex',FontSize=18)
end

legend(Interpreter="latex",FontSize=14,Location='southwest')
