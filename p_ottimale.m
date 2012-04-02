intervalli=20;
pmed=zeros(1,intervalli);
beta=linspace(0.4,1,intervalli);
samples=ceil(1000/intervalli*10)
fit = [    2.1417   -5.6286   -7.7126];
for i=1:intervalli
    if(beta(i)<0.4)
       prob=0.5;    
    else
       prob=exp(polyval(fit,log(beta(i))));
    end
    
    L = 2*round(double(rand(200)<prob))-1;
    [L,p]=MetroIsing2d(200,beta(i),samples,L);
    pmed(i)=mean(p((end/2):end));
    disp(['beta ' num2str(beta(i)) ', p: ' num2str(pmed(i)) ', stimata: ' num2str(prob)])
    plot(p);
    title(['beta ' num2str(beta(i))]);
    pause(0.1);
end
