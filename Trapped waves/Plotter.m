figure(1); clf; hold on;


%loadstring1='trappedc';
%loadstring1='N641trappedc';
loadstring1='N1241trappedc';



for i=1:8
    loadstring=strcat(loadstring1,num2str(i),'.mat');
    load(loadstring);
    plot(Phi,Ys_Newton)
    legend
    cs(i)=c;
end
legend

figure(2); clf; hold on;
plot(cs,'-x')