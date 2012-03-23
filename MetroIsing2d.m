function [rho, L] = MetroIsing2d(N,bt,T,L0)
%function rho = MetroIsing2d(N,bt,T,L0)

rho = zeros(T,1);
% if nargin < 4
%     cicle=30;
% end
    
fid=fopen('ising.bin','w');

if (nargin ==4 && length(L0) == N)
    L=L0;
else
    L = 2*round(rand(N))-1;
end;

nnr = reshape(1:N^2,N,N);
nnd = circshift(nnr,-1);
nnu = circshift(nnr,1);
nnl = (circshift(nnr',1))';
nnr = (circshift(nnr',-1))';


elenco=reshape(ceil((1:N^2)./N),N,N);
elenco=elenco+elenco';
dispari=find(mod(elenco,2));
pari=find(~mod(elenco,2));
beta=bt;

t=1;
for t = 1:T,
%    beta=bt+0.4*sin(2*pi/200*t);
    exp_beta=exp(-beta*[-8 -4 0 4 8]');


%     N2 = randperm(N^2);
%         for j = N2,
%             dH = 2*L(j)*(L(nnr(j))+L(nnl(j))+L(nnu(j))+L(nnd(j)));
%             if (dH < 0) || (rand < exp_beta(dH/4+3))
%                 L(j) = - L(j);
%             end;
%         end;


    j=pari;
    dH = 2*L(j).*(L(nnr(j))+L(nnl(j))+L(nnu(j))+L(nnd(j)));
    accettati= (rand(1,length(j))'<exp_beta(dH/4+3));
    accettati=j(accettati>0);
    L(accettati) = - L(accettati);

    j=dispari;
    dH = 2*L(j).*(L(nnr(j))+L(nnl(j))+L(nnu(j))+L(nnd(j)));
    accettati= (rand(1,length(j))'<exp_beta(dH/4+3));
    accettati=j(accettati>0);
    L(accettati) = - L(accettati);

    if( mod(t,5)==0 )
        fwrite(fid,int32(L),'int32');
        if(nargout==0)
            imagesc(L);
            pause(0.5); 
        end
    end

    rho(t) = mean(mean(L));
end;
