function [xtarget] = fitta_andamento(dati,p)
x=dati.x;
y=dati.y;
target=dati.target;
beta=dati.beta;
% x=matrice(:,1);
% y=matrice(:,2);

fun = @(p,x) (p(4)*erf((log(x)-p(1))/p(2))+p(3)).*(log(x).^p(5));
chi = @(p) sum((y - fun(p,x)).^2);

if(nargin<4)
p=[6.0351    1.9022    0.4361    0.3994    0.2328];
end


options= optimset('Display','off','MaxFunEvals',1000, 'MaxIter',1000);
p1=fminsearch(chi,p,options);

% guess circa lo zero, e' verso exp(p(1))
guess=2*exp(p1(1));
if(guess<10)
    guess=10;
end
[xtarget,~,exitflag]=fzero(@(x11) fun(p1,x11)-target, guess,options);
while(exitflag<0)
    guess=guess*2;
    [xtarget,~,exitflag]=fzero(@(x11) fun(p1,x11)-target, guess,options);
end


xsteps=1:6000;
ybuffer=fun(p1,xsteps);
ysteps=linspace(0,y(end)-0.01,50);
xeffettivo=ysteps;

for i=1:length(ysteps)
    xeffettivo(i)=min(find(ybuffer>=ysteps(i)));
    
    %disp(sprintf('%d %f\n',ceil(xeffettivo),j))
end


semilogx(xeffettivo,ysteps,'o',xsteps,ybuffer);
hold on;
semilogx(xtarget,target,'^','MarkerSize',10,'MarkerFaceColor','red');
hold off;
title(['\beta =',num2str(beta)]);
pause(0.4)
end
