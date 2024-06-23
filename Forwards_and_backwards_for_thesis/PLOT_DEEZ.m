Ns=[81,141,241,341,441,841,1241,1641];


try 
    load('unforced_solitary_DEEZ.mat')
    norm_store=zeros(1,numel(Ns));
    figure(1); clf; hold on; box on;
    for i=1:numel(Ns)
        plot(Phi_t_store{i},D_store{i})
        norm_store(i)=sqrt(trapz(Phi_t_store{i},D_store{i}.^2));
    end

    refx=linspace(2,3.2,100);

    legend("N="+string(Ns'))
    figure(2); clf; hold on; box on;
    plot(log10(Ns),log10(norm_store),LineWidth=2)
    xlabel('$\log_{10}(N)$',Interpreter='latex',FontSize=18)
    ylabel('$\log_{10}(|$\mbox{\boldmath$\underline{b}$}$|)$',Interpreter='latex',FontSize=18)
    plot(refx,-2*refx+2.0,'--k')


catch ME
    figure(1); clf; hold on; box on; 
    for i=1:numel(Ns)
        N=Ns(i);
        load("unforced_solitary_F1p1_N"+num2str(N)+".mat")
        [Phi_t_store{i},~,~,~,~,~,~,~,D_store{i},~] = svd_Step_fun(L,Froude,N,N,P,Ys_Newton);
        plot(Phi_t_store{i},D_store{i})
    end
end
        