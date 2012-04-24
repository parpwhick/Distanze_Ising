function stima_media_j(a,b,betamax)

if (nargin < 2)
    a=0.1;
    b=1;
end

if (nargin < 3)
    betamax=10;
end

beta=linspace(1,betamax,100);
xi1=zeros(size(beta));
xi2=xi1;

%f = @(p,x) (p(1)*x+p(2))./(x+p(3));
%pdf=@(x) unifpdf(x,0.1,0.55);
%r=200:20:500;

% alt=(-1).^(1:30);
% alt(end)=0;
% mytan = @(x) last(filter(1,[1 -x],alt));



for i=1:length(beta)
    
xi1(i)=-1/log(tanh(beta(i)));

j=10000000;

%chi=@(p) sum((f(p,r(dove))-and(dove)).^2);
%p=fminsearch(@(p) chi(p),[1,1,1]');


corr2=quadgk(@(x) 1/beta(i) * (tanh(x)).^j,a*beta(i),b*beta(i));
%corr2=quadgk(@(x) tanh(beta(i)* x).^j,a,b);
xi2(i)=-j/log(corr2);
%plot(beta,xi1,'sr');
end

plot(beta,xi2,'+-');

end