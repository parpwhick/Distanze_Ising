x=linspace(0,1,100);
r=1:100;

% for beta=linspace(0.1,4,
% 
% corr1=tanh(beta).^r;
% corr2=corr1;
% 
% p=polyfit(r,log(corr1),1);
% xi1=-1/p(1);
% 
% for i=1:100
% corr2(i)=quad(@(x) x.^i .* tanh(x * beta).^i,0,1,10e-12);
% end
% 
% semilogy(r,corr1,r,corr2)
% p=polyfit(r,log(corr2),1);
% xi2=-1/p(1);
% 
% plot(beta,xi1,'sr',beta,xi2,'db');
% 
% end
beta=linspace(1,10,100);
xi1=zeros(size(beta));
xi2=xi1;

%p=@(x) unifpdf(x,0.5,1.5);
%p=@(x) unifpdf(x,0.5,1);
p=@(x) unifpdf(x,0,1);
%p=@(x) unifpdf(x,1.5,2.5);

for i=1:length(beta)

xi1=-log(tanh(beta(i)));

corr2=quadv(@(x) p(x) .* (log(x) + log(tanh(x*beta(i)))),0,1,10e-8);
xi2(i)=abs(1/corr2);
%plot(beta,xi1,'sr');
end

plot(beta,xi2,'.-');