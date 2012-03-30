function [xtarget] = fitta_andamento(dati,p)
x=dati.x;
y=dati.y;
target=dati.target;
beta=dati.beta;
% x=matrice(:,1);
% y=matrice(:,2);

fun = @(p,x) (p(4)*erf(real((log(x)-p(1))/p(2)))+p(3)).*(log(x).^p(5));
chi = @(p) sum((y - fun(p,x)).^2);
options= optimset('Display','off','MaxFunEvals',1000, 'MaxIter',1000);


if(nargin<2)
    p01=[1.4223    1.8766    0.4851    0.6500    0.0000];
    p02=[6.0351    1.9022    0.4361    0.3994    0.2328];
    [p1,chi1]=fminsearch(chi,p01,options);
    [p2,chi2]=fminsearch(chi,p02,options);
    if(chi2<chi1)
        p1=p2;
    end
else
    p1=fminsearch(chi,p,options);
end


xsteps=1:max(x);
ybuffer=fun(p1,xsteps);

semilogx(x,y,'o',xsteps,ybuffer);
hold on;

% guess circa lo zero, e' verso exp(p(1))
% guess=2*exp(p1(1));
% if(guess<5)
%     guess=5;
% end
% [xtarget,~,exitflag]=fzero(@(x11) fun(p1,x11)-target, guess,options);
% while(exitflag<0)
%     guess=guess*2;
%     [xtarget,~,exitflag]=fzero(@(x11) fun(p1,x11)-target, guess,options);
% end
xtarget=fminbnd(@(x11) (fun(p1,x11)-target).^2 ,2,max(x)*2/3);


% ysteps=linspace(0,y(end)-0.01,50);
% xeffettivo=ysteps;
% 
% for i=1:length(ysteps)
%     xeffettivo(i)=min(find(ybuffer>=ysteps(i)));
%     
%     disp(sprintf('%d %f\n',ceil(xeffettivo),j))
% end

semilogx(xtarget,target,'^','MarkerSize',10,'MarkerFaceColor','red');
hold off;
title(['\beta =',num2str(beta)]);
pause(0.4)
end
